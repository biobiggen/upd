#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (C) 2024 biobiggen
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
Core UPD computation

Method chain used for a UPD calculation:
    get_probe_type
      -> get_background_gt
      -> get_mother_genotype
      -> get_fetal_fraction
      -> predict_fetal_genotype_hypergeom
      -> predict_upd_hmm
"""

import os

# Pin BLAS to a single thread to avoid thread contention during batch runs
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
os.environ.setdefault('VECLIB_NUM_THREADS', '1')
os.environ.setdefault('NUMEXPR_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

import logging
from collections import Counter
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import scipy.stats as ss

from .regions import UPD_TARGETS_HMM

logger = logging.getLogger(__name__)


# ============================================================================
# Probe version configuration
# ============================================================================
PROBE_VERSIONS: Dict[str, Dict] = {
    'NIPT3V3': {
        'expected_rows': 7749,
        'ff_threshold': 0.03,
        'depth_threshold': 650,
        'ignore_snps': ['HBA', 'RHD', 'SMN', 'HAP', 'other'],
    },
    'NIPT3V4': {
        'expected_rows': 13833,
        'ff_threshold': 0.03,
        'depth_threshold': 650,
        'ignore_snps': ['HBA', 'RHD', 'SMN', 'HAP', 'other'],
    },
}


# ============================================================================
# HMM parameters
# ============================================================================
OBS_MAPPING: Dict[str, int] = {
    'BBBB': 0, 'BBBA': 1, 'BBAA': 2,
    'AABB': 3, 'AABA': 4, 'AAAA': 5,
}

STATE_NAMES: List[str] = ['Normal', 'UPDM', 'UPDPI', 'UPDPII']

STARTPROB = np.array([0.9, 0.04, 0.03, 0.03])

TRANSMAT = np.array([
    [0.9, 0.04, 0.03, 0.03],
    [0.05, 0.95, 0.00, 0.00],
    [0.05, 0.00, 0.65, 0.30],
    [0.05, 0.00, 0.30, 0.65],
])

EMISSIONPROB = np.array([
    [0.25, 0.25, 0.0, 0.0, 0.25, 0.25],
    [0.50, 0.00, 0.00, 0.00, 0.00, 0.50],
    [0.125, 0.25, 0.125, 0.125, 0.25, 0.125],
    [0.25, 0.00, 0.25, 0.25, 0.00, 0.25],
])

# Decision thresholds
MIN_REGION_SNPS = 50           # minimum SNPs per region
MIN_HOMOZYGOUS_SNPS = 20       # minimum maternal homozygous sites
MIN_OBSERVATIONS = 20          # minimum observations
MIN_REGION_LENGTH = 20         # minimum segment length
DEPTH_FILTER = 400             # depth filter threshold
NORMAL_RATIO_THRESHOLD = 0.2   # Normal state decision threshold
SIGNIFICANT_UPD_RATIO = 0.1    # significant UPD ratio threshold

# ---- Dual-track fetal fraction estimation parameters ----
# The ffAB estimate from maternal homozygous sites implicitly assumes the fetus
# is heterozygous there (AF = ff/2). If the fetus is homozygous genome-wide it
# is homozygous at those sites too (AF = ff), so this path overestimates FF by
# about 2x, which makes BBAA/AABB observations be misread as BBBA/AABA --
# exactly erasing the only signature of paternal isodisomy.
# At maternal heterozygous (BA) sites |AF - 0.5| = (ff/2) * |fetal alt dosage - 1|,
# which is always ff/2 for a homozygous fetus regardless of whether the fetal
# genome comes from both parents or a single parent, hence robust to
# genome-wide homozygosity.
MIN_HET_SNPS_FOR_FF = 30       # minimum sites required by the BA path
FF_RATIO_HOM_UPPER = 1.5       # ff_hom/ff_het above this suggests genome-wide homozygosity
HET_OFFSET_MIN = 0.02          # |AF-0.5| floor, filters out noise-driven zero offsets
HET_BAND = 0.25                # AF window for maternal BA sites (0.5 +/- HET_BAND)

# Depth threshold for the fetal-fraction homozygous-sites path (separate from
# DEPTH_FILTER, which gates the HMM region filtering).
FF_DEPTH_THRESHOLD = 450


class UPDCalculator:
    """UPD calculator.

    Args:
        probe_version: probe version ('NIPT3V3' / 'NIPT3V4')
        probe_file: path to the probe file (optional, used for SNP_Tag mapping)
        test_code: test code type (affects fetal fraction thresholds)
        cli: sample identifier (derived from the readlist file name by default)
    """

    def __init__(
        self,
        probe_version: str = 'NIPT3V4',
        probe_file: Optional[str] = None,
        test_code: int = 1,
        cli: Optional[str] = None,
    ):
        if probe_version not in PROBE_VERSIONS:
            raise ValueError(
                f'Unsupported probe version: {probe_version}. '
                f'Available: {list(PROBE_VERSIONS.keys())}'
            )

        self.probe_version = probe_version
        self.probe_file = str(probe_file) if probe_file is not None else None
        self.test_code = test_code
        self.cli = cli

        config = PROBE_VERSIONS[probe_version]
        self.expected_rows = config['expected_rows']
        self.ff_threshold = config['ff_threshold']
        self.depth_threshold = config['depth_threshold']
        self.ignore_snps = list(config['ignore_snps'])

        self.readlist_df: Optional[pd.DataFrame] = None
        self.readlist_file: Optional[str] = None

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------
    def load_readlist(self, path: str) -> None:
        """Load a snp.readslist.txt file."""
        try:
            self.readlist_df = pd.read_csv(path, sep='\t')
        except Exception as e:
            logger.error(f'Error reading readlist file {path}: {e}')
            raise

        if not self.readlist_df.shape[0]:
            raise ValueError('The input file is empty.')

        self.readlist_file = path
        if self.cli is None:
            self.cli = '_'.join(os.path.basename(path).split('_')[:-1])

    def validate_row_count(self) -> None:
        """Check that the readlist row count matches the probe version."""
        actual = self.readlist_df.shape[0]
        if actual != self.expected_rows:
            raise ValueError(
                f'Readlist file is not {self.probe_version}: '
                f'expected {self.expected_rows} rows, got {actual}'
            )

    # ------------------------------------------------------------------
    # Probe type mapping
    # ------------------------------------------------------------------
    def get_probe_type(self, probe_file: Optional[str] = None) -> None:
        """Annotate each SNP with a SNP_Tag based on the probe file."""
        probe_file = probe_file if probe_file is not None else self.probe_file

        if probe_file is None or not os.path.exists(probe_file):
            logger.warning(f'Probe file not found: {probe_file}')
            if 'SNP_Tag' not in self.readlist_df.columns:
                self.readlist_df['SNP_Tag'] = 'chrRef'
            return

        try:
            if str(probe_file).endswith('.csv'):
                df_probe = pd.read_csv(probe_file, sep=',')
            else:
                df_probe = pd.read_csv(probe_file, sep='\t')
        except Exception as e:
            logger.error(f'Error reading probe file {probe_file}: {e}')
            if 'SNP_Tag' not in self.readlist_df.columns:
                self.readlist_df['SNP_Tag'] = 'chrRef'
            return

        required_columns = ['#Chr', 'Pos']
        missing = [col for col in required_columns if col not in df_probe.columns]
        if missing:
            logger.warning(f'Required columns missing in probe file: {missing}')
            if 'SNP_Tag' not in self.readlist_df.columns:
                self.readlist_df['SNP_Tag'] = 'chrRef'
            return

        type_column = None
        candidates = ['Type', 'type', 'TYPE', 'SNP_Tag', 'SNP_TAG',
                      'tag', 'TAG', 'Probe_Type', 'probe_type']
        for col in candidates:
            if col in df_probe.columns:
                type_column = col
                break

        if type_column is None:
            logger.warning(
                f'No type column found in probe file. '
                f'Available columns: {list(df_probe.columns)}'
            )
            if 'SNP_Tag' not in self.readlist_df.columns:
                self.readlist_df['SNP_Tag'] = 'chrRef'
            return

        try:
            probe_map = df_probe.set_index(['#Chr', 'Pos'])[type_column].to_dict()
            self.readlist_df['SNP_Tag'] = self.readlist_df.apply(
                lambda row: probe_map.get((row['#Chr'], row['Pos']), 'chrRef'),
                axis=1,
            )
        except Exception as e:
            logger.error(f'Error processing probe file {probe_file}: {e}')
            if 'SNP_Tag' not in self.readlist_df.columns:
                self.readlist_df['SNP_Tag'] = 'chrRef'

    # ------------------------------------------------------------------
    # Genotype inference
    # ------------------------------------------------------------------
    def get_background_gt(self) -> None:
        """Infer the background genotype from pA_Ratio."""
        def get_gt(x) -> str:
            if x < 0.25:
                return 'BB'
            elif x <= 0.75:
                return 'BA'
            return 'AA'

        self.background_gt = self.readlist_df['pA_Ratio'].apply(get_gt)

    def get_mother_genotype(self) -> pd.Series:
        """Return the maternal genotype (identical to the background genotype)."""
        if not hasattr(self, 'mother_genotype'):
            if not hasattr(self, 'background_gt'):
                self.get_background_gt()
            self.mother_genotype = self.background_gt
        return self.mother_genotype

    # ------------------------------------------------------------------
    # Fetal fraction
    # ------------------------------------------------------------------
    def _estimate_ff_het(self, df: pd.DataFrame, dep_threshold: float) -> Dict:
        """Estimate the fetal fraction from maternal heterozygous (BA) sites.

        At maternal BA sites the expected plasma alt ratio is::

            AF = (1 - ff) * 0.5 + ff * fetal_alt_dosage / 2
               = 0.5 + (ff / 2) * (fetal_alt_dosage - 1)

        that is, ``|AF - 0.5| = (ff / 2) * |fetal_alt_dosage - 1|``. When the
        fetus is homozygous (dosage 0 or 2) the offset is always ``ff / 2``.
        This relation holds regardless of whether the fetal genome comes from
        both parents or entirely from one parent, hence it is robust to
        genome-wide homozygosity.

        FF is twice the median offset over sites whose offset exceeds
        ``HET_OFFSET_MIN`` (i.e. sites where the fetus is homozygous).

        Args:
            df: readlist with ignore_snps already removed
            dep_threshold: depth floor

        Returns:
            ``{'ff': float, 'n_het': int, 'n_shifted': int}``; ``ff`` is
            ``None`` when there are too few sites
        """
        chrom_col = df.columns[0]
        index_auto = ~df[chrom_col].isin(['chrX', 'chrY'])
        index_snp = df.Alt != df.Ref
        index_dep = df.Depth >= dep_threshold
        index_het = df.pA_Ratio.between(0.5 - HET_BAND, 0.5 + HET_BAND)

        het_af = df.loc[index_auto & index_snp & index_dep & index_het, 'pA_Ratio']
        n_het = int(len(het_af))
        if n_het < MIN_HET_SNPS_FOR_FF:
            return {'ff': None, 'n_het': n_het, 'n_shifted': 0}

        offset = (het_af - 0.5).abs()
        shifted = offset[offset > HET_OFFSET_MIN]
        n_shifted = int(len(shifted))
        if n_shifted < MIN_HET_SNPS_FOR_FF:
            return {'ff': None, 'n_het': n_het, 'n_shifted': n_shifted}

        ff_het = 2 * float(np.median(shifted))
        return {'ff': ff_het, 'n_het': n_het, 'n_shifted': n_shifted}

    def get_fetal_fraction(self) -> float:
        """Compute the fetal fraction."""
        if 'SNP_Tag' not in self.readlist_df.columns:
            self.get_probe_type(self.probe_file)

        df = self.readlist_df.loc[
            ~self.readlist_df['SNP_Tag'].isin(self.ignore_snps)
        ].copy()
        df.reset_index(drop=True, inplace=True)

        index_snp = df.Alt != df.Ref
        index_chrx = df[df.columns[0]] == 'chrX'
        index_chry = df[df.columns[0]] == 'chrY'

        min_threshold = 0.01
        if self.test_code == 2:
            min_threshold = 0.025

        index_BB = df.pA_Ratio.between(min_threshold, 0.25)
        index_AA = df.pA_Ratio.between(0.75, 1 - min_threshold)
        query_index = index_snp & (index_BB | index_AA) & ~index_chrx & ~index_chry

        base_condition = index_snp & ~index_chrx & ~index_chry
        query_index2 = base_condition & (df.Depth >= FF_DEPTH_THRESHOLD)

        fbref = df.loc[query_index2 & index_BB, 'pA_Ratio']
        faref = df.loc[query_index2 & index_AA, 'pA_Ratio']

        if len(fbref) > 10 and len(faref) > 10:
            ffAA = 2 * (1.0 - np.median(faref))
            ffBB = 2 * np.median(fbref)
            ffAB = (ffAA + ffBB) / 2.0
        else:
            ffAB = 1e-6

        af_cutoff = max(min_threshold, ffAB / 5)
        fbref = df.loc[
            query_index2
            & df.pA_Ratio.gt(af_cutoff)
            & df.pA_Ratio.lt(0.25),
            'pA_Ratio',
        ]
        faref = df.loc[
            query_index2
            & df.pA_Ratio.gt(0.75)
            & df.pA_Ratio.lt(1 - af_cutoff),
            'pA_Ratio',
        ]

        if len(fbref) > 10 and len(faref) > 10:
            ffAA = 2 * (1.0 - np.median(faref))
            ffBB = 2 * np.median(fbref)
            ffAB = (ffAA + ffBB) / 2.0
        else:
            ffAA = ffBB = ffAB = 1e-6

        # ---- Dual-track estimate: maternal homozygous vs maternal BA path ----
        # ffAB comes from homozygous sites and implicitly assumes the fetus is
        # heterozygous there; a genome-wide homozygous fetus makes it about 2x
        # too high. The independent BA estimate is used as a cross-check, and
        # when the ratio is clearly high the BA estimate is used instead so that
        # downstream does not misread BBAA/AABB as BBBA/AABA.
        ff_hom = ffAB
        het_res = self._estimate_ff_het(df, FF_DEPTH_THRESHOLD)
        ff_het = het_res['ff']

        ff_ratio = None
        ff_method = 'homozygous'
        if ff_het is not None and ff_het > 0:
            if ff_hom <= 0 or ff_hom <= 1e-5:
                # The homozygous path failed (set to the 1e-6 sentinel when
                # there are too few sites); use the BA path when available
                ff_method = 'heterozygous'
                ffAB = ff_het
                logger.info(
                    f'FF from homozygous sites unavailable ({ff_hom:.2e}), '
                    f'using maternal heterozygous estimate: {ff_het:.4f}'
                )
            else:
                ff_ratio = ff_hom / ff_het
                if ff_ratio > FF_RATIO_HOM_UPPER:
                    ff_method = 'heterozygous'
                    ffAB = ff_het
                    logger.info(
                        f'FF corrected by maternal heterozygous sites: '
                        f'{ff_hom:.4f} -> {ff_het:.4f} (ratio={ff_ratio:.2f}), '
                        f'genome-wide homozygosity suspected'
                    )

        self.ff_hom = ff_hom
        self.ff_het = ff_het
        self.ff_ratio = ff_ratio
        self.ff_method = ff_method
        self.ff_het_sites = het_res['n_het']
        self.ff_het_shifted_sites = het_res['n_shifted']

        self.fetal_fraction = ffAB
        self.ff = self.fetal_fraction
        return float(self.fetal_fraction)

    # ------------------------------------------------------------------
    # Fetal genotype prediction
    # ------------------------------------------------------------------
    @staticmethod
    def _calculate_expected_alt_ratio(
        mother_genotype: str,
        fetus_genotype: str,
        fetal_fraction: float,
    ) -> float:
        """Compute the expected alt allele ratio from the maternal/fetal
        genotypes and the fetal fraction."""
        mother_alt_ratio = {'BB': 0, 'BA': 0.5, 'AA': 1}[mother_genotype]
        fetus_alt_ratio = {'BB': 0, 'BA': 0.5, 'AA': 1}[fetus_genotype]
        return (1 - fetal_fraction) * mother_alt_ratio + fetal_fraction * fetus_alt_ratio

    def predict_fetal_genotype_hypergeom(self) -> Dict:
        """Predict the fetal genotype from binomial likelihoods."""
        if 'pA_Ratio' not in self.readlist_df.columns:
            raise ValueError('pA_Ratio column not found in data')

        if 'Depth' not in self.readlist_df.columns:
            raise ValueError('Depth column not found in data')

        if not hasattr(self, 'ff'):
            self.get_fetal_fraction()
        fetal_fraction = self.ff

        mother_gt = self.get_mother_genotype()

        genotypes: Dict[str, List] = {
            'mother': [],
            'fetus_predicted': [],
            'fetus_likelihoods': [],
            'chr': [],
            'pos': [],
            'pA_Ratio': [],
            'max_likelihood': [],
            'depth': [],
        }

        possible_genotypes = ['BB', 'BA', 'AA']

        for index, row in self.readlist_df.iterrows():
            chr_name = row['#Chr']
            position = row['Pos']
            pA_ratio = row['pA_Ratio']
            depth = row['Depth']
            mother_genotype = mother_gt.loc[index]

            if pA_ratio < max(0.01, fetal_fraction / 5):
                simple_prediction = 'BB'
            elif pA_ratio > min(0.99, 1 - fetal_fraction / 5):
                simple_prediction = 'AA'
            else:
                simple_prediction = None

            likelihoods = {}
            total_reads = depth
            alt_reads = int(pA_ratio * total_reads)

            for fetus_genotype in possible_genotypes:
                expected_alt_ratio = self._calculate_expected_alt_ratio(
                    mother_genotype, fetus_genotype, fetal_fraction
                )

                if expected_alt_ratio <= 0:
                    expected_alt_ratio = 1e-2
                elif expected_alt_ratio >= 1:
                    expected_alt_ratio = 1 - 1e-2

                if alt_reads == 0:
                    log_likelihood = total_reads * np.log(1 - expected_alt_ratio)
                elif alt_reads == total_reads:
                    log_likelihood = total_reads * np.log(expected_alt_ratio)
                else:
                    log_likelihood = ss.binom.logpmf(
                        alt_reads, total_reads, expected_alt_ratio
                    )

                likelihoods[fetus_genotype] = log_likelihood

            max_likelihood_genotype = max(likelihoods, key=likelihoods.get)

            if simple_prediction is not None:
                predicted_genotype = simple_prediction
            else:
                predicted_genotype = max_likelihood_genotype

            genotypes['mother'].append(mother_genotype)
            genotypes['fetus_predicted'].append(predicted_genotype)
            genotypes['fetus_likelihoods'].append(likelihoods)
            genotypes['chr'].append(chr_name)
            genotypes['pos'].append(position)
            genotypes['pA_Ratio'].append(pA_ratio)
            genotypes['max_likelihood'].append(likelihoods[max_likelihood_genotype])
            genotypes['depth'].append(depth)

        return genotypes

    # ------------------------------------------------------------------
    # UPD HMM prediction
    # ------------------------------------------------------------------
    def predict_upd_hmm(self) -> Dict:
        """Predict fetal UPD with a hidden Markov model.

        Analyses the fetal genotypes at maternal homozygous sites and uses an
        HMM to detect the four states UPDM / UPDPI / UPDPII / Normal.

        Returns:
            Dictionary of per-region UPD prediction results
        """
        required_columns = ['#Chr', 'Pos', 'pA_Ratio', 'Depth']
        for col in required_columns:
            if col not in self.readlist_df.columns:
                raise ValueError(f'{col} column not found in data')

        if not hasattr(self, 'ff'):
            self.get_fetal_fraction()

        if 'SNP_Tag' not in self.readlist_df.columns:
            self.get_probe_type(self.probe_file)

        fetal_genotypes = self.predict_fetal_genotype_hypergeom()

        upd_results: Dict[str, Dict] = {}

        genotypes_df = pd.DataFrame({
            'chr': fetal_genotypes['chr'],
            'pos': fetal_genotypes['pos'],
            'mother_gt': fetal_genotypes['mother'],
            'fetus_gt': fetal_genotypes['fetus_predicted'],
            'pA_Ratio': fetal_genotypes['pA_Ratio'],
            'depth': fetal_genotypes['depth'],
        })

        self.fetal_genotypes = genotypes_df
        genotypes_df.index = self.readlist_df.index

        for region_name, (chrom, start_pos, end_pos) in UPD_TARGETS_HMM.items():
            if region_name == 'chrRef':
                chr_ref_mask = self.readlist_df['SNP_Tag'] == 'chrRef'
                region_indices = self.readlist_df[chr_ref_mask].index
                region_data = genotypes_df.loc[region_indices].copy()
            elif chrom is not None and chrom != 'chrRef':
                region_mask = (
                    (genotypes_df['chr'] == chrom)
                    & (genotypes_df['pos'] >= start_pos)
                    & (genotypes_df['pos'] <= end_pos)
                )
                region_data = genotypes_df[region_mask].copy()
            else:
                region_data = genotypes_df[genotypes_df['chr'] == chrom].copy()

            if len(region_data) < MIN_REGION_SNPS:
                upd_results[region_name] = {
                    'status': 'insufficient_data',
                    'regions': [],
                }
                continue

            depth_filter = region_data['depth'] >= DEPTH_FILTER
            hom_filter = (
                (region_data['mother_gt'] == 'AA')
                | (region_data['mother_gt'] == 'BB')
            )
            filtered_data = region_data[depth_filter & hom_filter]

            if len(filtered_data) < MIN_HOMOZYGOUS_SNPS:
                upd_results[region_name] = {
                    'status': 'insufficient_homozygous',
                    'regions': [],
                }
                continue

            obs_sequence = []
            positions = []

            for _, row in filtered_data.iterrows():
                pos = row['pos']
                mat_gt = row['mother_gt']
                fetal_gt = row['fetus_gt']

                if mat_gt in ['AA', 'BB']:
                    obs_sequence.append(mat_gt + fetal_gt)
                    positions.append(pos)

            if len(obs_sequence) < MIN_OBSERVATIONS:
                upd_results[region_name] = {
                    'status': 'insufficient_observations',
                    'regions': [],
                }
                continue

            try:
                # All mother genotype is AA/BB here (filtered_data), and fetal
                # genotype is BB/BA/AA, so every combo is present in OBS_MAPPING.
                # Map directly so obs_indices stays aligned with positions and
                # filtered_data; an unknown observation raises KeyError instead
                # of being silently dropped (which would misalign downstream).
                obs_indices = [OBS_MAPPING[obs] for obs in obs_sequence]
                if len(obs_indices) == 0:
                    upd_results[region_name] = {
                        'status': 'no_valid_observations',
                        'regions': [],
                    }
                    continue

                obs_array = np.array(obs_indices).reshape(-1, 1)
                obs_array = np.clip(obs_array, 0, 5)
            except Exception as e:
                upd_results[region_name] = {
                    'status': 'observation_conversion_error',
                    'regions': [],
                    'error': str(e),
                }
                continue

            used_obs = sorted(set(obs_indices))
            if len(used_obs) < 2:
                upd_results[region_name] = {
                    'status': 'insufficient_observation_types',
                    'regions': [],
                }
                continue

            try:
                from hmmlearn import hmm

                model = hmm.CategoricalHMM(
                    n_components=4,
                    n_features=6,
                    tol=1e-3,
                    n_iter=100,
                )
                model.startprob_ = STARTPROB
                model.transmat_ = TRANSMAT
                model.emissionprob_ = EMISSIONPROB

                states = model.predict(obs_array)

                state_counts = Counter(states)
                total_states = len(states)

                state_ratios = {
                    STATE_NAMES[i]: count / total_states
                    for i, count in state_counts.items()
                }

                regions = []
                current_state = None
                region_start = None
                region_positions = []

                site_details = []
                for i, (state, pos) in enumerate(zip(states, positions)):
                    site_details.append({
                        'position': pos,
                        'state': STATE_NAMES[state],
                        'state_index': state,
                        'fetal_genotype': (
                            filtered_data.iloc[i]['fetus_gt']
                            if i < len(filtered_data) else None
                        ),
                        'maternal_genotype': (
                            filtered_data.iloc[i]['mother_gt']
                            if i < len(filtered_data) else None
                        ),
                        'pa_ratio': (
                            filtered_data.iloc[i]['pA_Ratio']
                            if i < len(filtered_data) else None
                        ),
                    })
                    if current_state is None:
                        current_state = state
                        region_start = pos
                        region_positions = [pos]
                    elif state == current_state:
                        region_positions.append(pos)
                    else:
                        if len(region_positions) >= MIN_REGION_LENGTH:
                            regions.append({
                                'start': region_start,
                                'end': region_positions[-1],
                                'state': STATE_NAMES[current_state],
                                'length': len(region_positions),
                                'ratio': len(region_positions) / total_states,
                            })

                        current_state = state
                        region_start = pos
                        region_positions = [pos]

                if len(region_positions) >= MIN_REGION_LENGTH:
                    regions.append({
                        'start': region_start,
                        'end': region_positions[-1],
                        'state': STATE_NAMES[current_state],
                        'length': len(region_positions),
                        'ratio': len(region_positions) / total_states,
                    })

                significant_upds = [
                    r for r in regions
                    if r['state'] != 'Normal' and r['ratio'] > SIGNIFICANT_UPD_RATIO
                ]

                if ('Normal' in state_ratios
                        and state_ratios['Normal'] > NORMAL_RATIO_THRESHOLD):
                    final_state = 'Normal'
                    final_ratio = state_ratios['Normal']
                else:
                    final_state = max(state_ratios, key=state_ratios.get)
                    final_ratio = state_ratios[final_state]

                upd_results[region_name] = {
                    'status': 'success',
                    'state_ratios': state_ratios,
                    'regions': regions,
                    'significant_upds': significant_upds,
                    'final_state': final_state,
                    'final_ratio': final_ratio,
                    'total_observations': len(obs_sequence),
                    'chromosome': chrom,
                    'start_pos': (
                        start_pos if start_pos is not None
                        else region_data['pos'].min() if len(region_data) > 0 else 0
                    ),
                    'end_pos': (
                        end_pos if end_pos is not None
                        else region_data['pos'].max() if len(region_data) > 0 else 0
                    ),
                    'site_details': site_details,
                }

            except Exception as e:
                upd_results[region_name] = {
                    'status': 'hmm_error',
                    'error': str(e),
                    'regions': [],
                }
                logger.error(f'Error calculating UPD: {e}')

        self.upd_results = upd_results

        # Fetal fraction diagnostics: a clearly high ff_ratio itself suggests
        # genome-wide homozygosity, so both estimates are kept for review
        self.ff_info = {
            'ff_used': self.ff,
            'ff_homozygous': self.ff_hom,
            'ff_heterozygous': self.ff_het,
            'ff_ratio': self.ff_ratio,
            'ff_method': self.ff_method,
            'het_sites': self.ff_het_sites,
            'het_shifted_sites': self.ff_het_shifted_sites,
        }

        return upd_results


# ============================================================================
# Helper functions
# ============================================================================
def convert_numpy_types(obj):
    """Recursively convert NumPy types to native Python types for JSON."""
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {key: convert_numpy_types(value) for key, value in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [convert_numpy_types(item) for item in obj]
    return obj


def calculate_upd(
    readlist_file: str,
    probe_file: Optional[str] = None,
    probe_version: str = 'NIPT3V4',
    test_code: int = 1,
    with_ff_info: bool = False,
) -> Dict:
    """Convenience function: run a full single-sample UPD calculation.

    Args:
        readlist_file: path to the snp.readslist.txt file
        probe_file: path to the probe file
        probe_version: probe version
        test_code: test code type
        with_ff_info: when True, also return the ``_ff_info`` key (fetal
            fraction diagnostics)

    Returns:
        UPD result dictionary (already converted to JSON-serialisable types)
    """
    calculator = UPDCalculator(
        probe_version=probe_version,
        probe_file=probe_file,
        test_code=test_code,
    )
    calculator.load_readlist(readlist_file)
    results = convert_numpy_types(calculator.predict_upd_hmm())
    if with_ff_info:
        results['_ff_info'] = convert_numpy_types(calculator.ff_info)
    return results
