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
Small simulated demo dataset generator

Real plasma data cannot be distributed publicly, so this script generates a
dataset from the forward model of the UPD computation that runs the whole
:mod:`upd` pipeline end to end:

- probe file ``demo_probe.bed`` (three columns: ``#Chr`` / ``Pos`` / ``Type``)
- one ``*_consensus.mapped.clipped.snp.reads.list`` file per sample

The expected plasma alt allele ratio matches
:meth:`UPDCalculator._calculate_expected_alt_ratio` exactly::

    pA_Ratio = (1 - ff) * maternal_alt_dosage/2 + ff * fetal_alt_dosage/2

Origin of the fetal genotype at maternal homozygous sites per state:

===========  ==========================================================
State        Fetal genotype
===========  ==========================================================
``Normal``   one maternal + one paternal allele
``UPDM``     both from the mother (maternal UPD), no paternal signal at
             homozygous sites
``UPDPI``    both from the father, two different paternal chromosomes
             (heterodisomy)
``UPDPII``   both copies of the same paternal chromosome (isodisomy),
             always homozygous
===========  ==========================================================

Usage::

    python -m upd.simulate_demo_data -o demo_data
"""

import argparse
import os
from typing import Dict, List, Optional, Tuple

import numpy as np

from .regions import UPD_TARGETS_HMM

# readlist column order: the first column must be #Chr
# (core.get_fetal_fraction relies on columns[0])
READLIST_COLUMNS: List[str] = [
    '#Chr', 'Pos', 'Ref', 'Alt', 'Depth', 'Ref_Dep', 'Alt_Dep',
    'pA_Ratio', 'GC', 'DepRegion', 'MapQ',
]

READLIST_SUFFIX = '_consensus.mapped.clipped.snp.reads.list'

BASES: Tuple[str, ...] = ('A', 'C', 'G', 'T')

# chrRef background windows: uniformly distributed on autosomes outside
# any UPD target region
REF_WINDOWS: List[Tuple[str, int, int]] = [
    ('chr1', 20000000, 240000000),
    ('chr3', 20000000, 190000000),
    ('chr4', 20000000, 180000000),
    ('chr5', 20000000, 175000000),
    ('chr8', 20000000, 140000000),
    ('chr10', 10000000, 125000000),
    ('chr12', 10000000, 125000000),
    ('chr19', 5000000, 55000000),
]

# Probe types filtered out by PROBE_VERSIONS['ignore_snps'], used to cover
# the filtering branch
IGNORED_PROBE_TYPES: List[Tuple[str, str, int, int]] = [
    ('HBA', 'chr16', 150000, 230000),
    ('SMN', 'chr5', 70900000, 71000000),
]

# Demo samples: (name, UPD state, UPD region, male)
DEMO_SAMPLES: List[Tuple[str, str, Optional[str], bool]] = [
    ('demoN01', 'Normal', None, True),
    ('demoN02', 'Normal', None, False),
    ('demoN03', 'Normal', None, True),
    ('demoN04', 'Normal', None, False),
    ('demoP01', 'UPDM', '15q11q13', True),
    ('demoP02', 'UPDPI', '11p15', False),
    ('demoP03', 'UPDPII', '11p15', True),
]


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Generate small simulated demo dataset for upd')
    parser.add_argument('-o', '--output', default='demo_data',
                        help='output directory (default: demo_data)')
    parser.add_argument('--n-probe', type=int, default=120,
                        help='number of SNPs per UPD target region (default: 120)')
    parser.add_argument('--n-ref', type=int, default=400,
                        help='total chrRef background SNPs (default: 400)')
    parser.add_argument('--depth', type=int, default=1200,
                        help='average sequencing depth (default: 1200)')
    parser.add_argument('--ff', type=float, default=0.10,
                        help='fetal fraction (default: 0.10)')
    parser.add_argument('--maf', type=float, default=0.5,
                        help='minor allele frequency at simulated sites (default: 0.5)')
    parser.add_argument('--seed', type=int, default=42,
                        help='random seed (default: 42)')
    return parser.parse_args()


def _sample_pos(rng, start: int, end: int, size: int) -> np.ndarray:
    """Sample positions without replacement within [start, end), avoiding
    materialising the entire coordinate range."""
    pos = set()
    while len(pos) < size:
        pos.update(int(x) for x in rng.integers(start, end, size=size * 2))
    return np.array(sorted(pos)[:size])


def make_probes(n_probe: int, n_ref: int, rng) -> List[Tuple[str, int, str]]:
    """Generate probe list ``[(chrom, pos, type), ...]`` sorted by chromosome
    and position."""
    probes: List[Tuple[str, int, str]] = []

    for region, (chrom, start, end) in UPD_TARGETS_HMM.items():
        if start is None or end is None:      # chrRef is generated from REF_WINDOWS
            continue
        for pos in _sample_pos(rng, start + 1, end, n_probe):
            probes.append((chrom, int(pos), region))

    per_window = max(1, n_ref // len(REF_WINDOWS))
    for chrom, start, end in REF_WINDOWS:
        for pos in _sample_pos(rng, start, end, per_window):
            probes.append((chrom, int(pos), 'chrRef'))

    for probe_type, chrom, start, end in IGNORED_PROBE_TYPES:
        for pos in _sample_pos(rng, start, end, 10):
            probes.append((chrom, int(pos), probe_type))

    # Sex chromosome sites: labelled separately to avoid being mixed into
    # chrRef HMM calculation
    for pos in _sample_pos(rng, 3100000, 155000000, 40):
        probes.append(('chrX', int(pos), 'chrX'))
    for pos in _sample_pos(rng, 2900000, 22000000, 20):
        probes.append(('chrY', int(pos), 'chrY'))

    chrom_order = {f'chr{i}': i for i in range(1, 23)}
    chrom_order['chrX'] = 23
    chrom_order['chrY'] = 24
    probes.sort(key=lambda x: (chrom_order[x[0]], x[1]))
    return probes


def write_probe_file(path: str, probes: List[Tuple[str, int, str]]) -> None:
    """Write probe file (tab-separated, used by ``-p/--probe-file``)."""
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('\t'.join(['#Chr', 'Pos', 'Type']) + '\n')
        for chrom, pos, probe_type in probes:
            fh.write(f'{chrom}\t{pos}\t{probe_type}\n')


def _fetal_dosage(state: str, mother_dosage: int, rng, maf: float) -> int:
    """Generate fetal alt allele dosage (0/1/2) according to UPD state."""
    # Allele transmitted from mother to fetus
    if mother_dosage == 2:
        maternal_allele = 1
    elif mother_dosage == 0:
        maternal_allele = 0
    else:
        maternal_allele = int(rng.binomial(1, 0.5))

    if state == 'UPDM':
        # Maternal UPD: both homologous chromosomes come from the mother
        return int(mother_dosage)
    if state == 'UPDPI':
        # Paternal heterodisomy: two different paternal homologous chromosomes
        return int(rng.binomial(2, maf))
    if state == 'UPDPII':
        # Paternal isodisomy: same paternal chromosome duplicated, always homozygous
        return 2 * int(rng.binomial(1, maf))
    # Normal: one maternal + one paternal
    return maternal_allele + int(rng.binomial(1, maf))


def simulate_sample(
    probes: List[Tuple[str, int, str]],
    ff: float,
    depth: int,
    rng,
    state: str = 'Normal',
    upd_region: Optional[str] = None,
    male: bool = True,
    maf: float = 0.5,
) -> List[List]:
    """Simulate readlist records for a single sample."""
    upd_regions = set()
    if state != 'Normal' and upd_region is not None:
        upd_regions = {upd_region}

    rows: List[List] = []
    for chrom, pos, probe_type in probes:
        site_state = state if probe_type in upd_regions else 'Normal'

        mother_dosage = int(rng.binomial(2, maf))
        fetal_dosage = _fetal_dosage(site_state, mother_dosage, rng, maf)

        if chrom == 'chrY':
            af = ff / 2 if male else 0.0
        elif chrom == 'chrX' and male:
            # Male fetus has only one X from the mother
            maternal_allele = (1 if mother_dosage == 2 else
                               0 if mother_dosage == 0 else
                               int(rng.binomial(1, 0.5)))
            af = (1 - ff) * mother_dosage / 2 + ff * maternal_allele
        else:
            af = (1 - ff) * mother_dosage / 2 + ff * fetal_dosage / 2
        af = float(min(max(af, 0.0), 1.0))

        dep = max(50, int(rng.normal(depth, depth * 0.12)))
        if chrom == 'chrY':
            dep = max(5, int(dep * (ff / 2 if male else 0.005)))

        alt_dep = int(rng.binomial(dep, af))
        ref_base, alt_base = rng.choice(BASES, size=2, replace=False)
        rows.append([
            chrom, pos, ref_base, alt_base, dep, dep - alt_dep, alt_dep,
            '%0.4f' % (alt_dep / dep),
            '%0.3f' % float(rng.uniform(0.35, 0.60)),
            int(dep * rng.uniform(0.95, 1.05)),
            60,
        ])
    return rows


def write_readlist(path: str, rows: List[List]) -> None:
    """Write readlist file (tab-separated)."""
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('\t'.join(READLIST_COLUMNS) + '\n')
        for row in rows:
            fh.write('\t'.join(map(str, row)) + '\n')


def generate(
    output_dir: str,
    n_probe: int = 120,
    n_ref: int = 400,
    depth: int = 1200,
    ff: float = 0.10,
    maf: float = 0.5,
    seed: int = 42,
) -> Dict[str, str]:
    """Generate the entire demo dataset, return ``{sample name or 'probe': file path}``."""
    rng = np.random.default_rng(seed)
    os.makedirs(output_dir, exist_ok=True)

    probes = make_probes(n_probe, n_ref, rng)
    probe_path = os.path.join(output_dir, 'demo_probe.bed')
    write_probe_file(probe_path, probes)

    written = {'probe': probe_path}
    for name, state, region, male in DEMO_SAMPLES:
        rows = simulate_sample(probes, ff, depth, rng, state=state,
                               upd_region=region, male=male, maf=maf)
        path = os.path.join(output_dir, f'{name}_{state.lower()}{READLIST_SUFFIX}')
        write_readlist(path, rows)
        written[name] = path
    return written


def main() -> None:
    args = get_args()
    written = generate(
        output_dir=args.output,
        n_probe=args.n_probe,
        n_ref=args.n_ref,
        depth=args.depth,
        ff=args.ff,
        maf=args.maf,
        seed=args.seed,
    )
    probe_path = written.pop('probe')
    print(f'Probe file: {probe_path}')
    print(f'Samples: {len(written)}  (fetal fraction {args.ff}, avg depth {args.depth})')
    for name, path in written.items():
        print(f'  {name}: {os.path.basename(path)}')


if __name__ == '__main__':
    main()
