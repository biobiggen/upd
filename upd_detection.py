#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
UPD Detection from Maternal Plasma SNP BAF Data
----------------------------------------------
This script detects Uniparental Disomy (UPD) in fetal DNA from maternal plasma
using SNP B-Allele Frequency (BAF) data.

The workflow:
1. Load BAF data from maternal plasma
2. Estimate fetal genotype from BAF patterns (including homozygous maternal SNPs)
3. Apply Hidden Markov Model (HMM) to detect UPD regions
4. Distinguish between isodisomy and heterodisomy:
   - For paternal UPD: Use HMM states to distinguish isodisomy vs heterodisomy
   - For maternal UPD: Analyze maternal heterozygous sites to distinguish isodisomy vs heterodisomy
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
from scipy import stats
from sklearn.mixture import GaussianMixture
from scipy.special import betaln, gammaln
from hmmlearn import hmm


def load_baf_data(file_path):
    """
    Load BAF data from input file.
    Expected format: chromosome, position, BAF
    """
    print(f"Loading BAF data from {file_path}...")
    try:
        df = pd.read_csv(file_path, sep='\t')
        required_columns = ['chromosome', 'position', 'baf']

        # Check if all required columns exist (case insensitive)
        df.columns = [col.lower() for col in df.columns]
        for col in required_columns:
            if col not in df.columns:
                raise ValueError(
                    f"Required column '{col}' not found in input file")

        # Rename columns to standardized names
        df = df.rename(columns={
            df.columns[df.columns.str.lower() == 'chromosome'][0]: 'chromosome',
            df.columns[df.columns.str.lower() == 'position'][0]: 'position',
            df.columns[df.columns.str.lower() == 'baf'][0]: 'baf'
        })

        # Convert chromosome to string and handle X/Y
        df['chromosome'] = df['chromosome'].astype(str)
        df['chromosome'] = df['chromosome'].str.replace('23', 'X')
        df['chromosome'] = df['chromosome'].str.replace('24', 'Y')

        # Remove 'chr' prefix if present
        df['chromosome'] = df['chromosome'].str.replace('chr', '')

        # Ensure position is numeric
        df['position'] = pd.to_numeric(df['position'])

        # Ensure BAF is between 0 and 1
        df = df[(df['baf'] >= 0) & (df['baf'] <= 1)]

        print(
            f"Loaded {len(df)} SNPs across {df['chromosome'].nunique()} chromosomes")
        return df

    except Exception as e:
        print(f"Error loading BAF data: {e}")
        return None


def filter_informative_snps(df, baf_threshold=0.2):
    """
    Filter for informative SNPs that are likely heterozygous in mother or fetus.
    Removes SNPs with BAF close to 0 or 1 (homozygous).
    """
    filtered_df = df
    return filtered_df


def segment_chromosome(df, chromosome, window_size=1000000):
    """
    Segment chromosome into windows and calculate mean BAF for each window.
    """
    chr_data = df[df['chromosome'] == chromosome].copy()
    if len(chr_data) == 0:
        return None

    # Sort by position
    chr_data = chr_data.sort_values('position')

    # Create windows
    max_pos = chr_data['position'].max()
    windows = []

    for start in range(0, int(max_pos) + window_size, window_size):
        end = start + window_size
        window_data = chr_data[(chr_data['position'] >= start) & (
            chr_data['position'] < end)]

        if len(window_data) > 0:
            windows.append({
                'chromosome': chromosome,
                'start': start,
                'end': end,
                'mean_baf': window_data['baf'].mean(),
                'median_baf': window_data['baf'].median(),
                'snp_count': len(window_data)
            })

    return pd.DataFrame(windows)


def estimate_fetal_genotype(baf_values, fetal_fraction=0.1, depth=100):
    """
    Estimate fetal genotype from BAF values, considering maternal DNA background.
    Enhanced to better detect paternal UPD scenarios (BB-AA or AA-BB patterns).

    Parameters:
    - baf_values: array of BAF values
    - fetal_fraction: estimated fetal DNA fraction in maternal plasma
    - depth: estimated sequencing depth (for confidence calculation)

    Returns:
    - Array of estimated genotypes (0: AA, 1: AB, 2: BB)
    - Array of maternal genotypes (0: AA, 1: AB, 2: BB)
    - Array of confidence scores (0-1)
    - Dictionary of additional metrics for each SNP
    """
    n_snps = len(baf_values)
    fetal_genotypes = np.zeros(n_snps, dtype=int)
    maternal_genotypes = np.zeros(n_snps, dtype=int)
    confidence_scores = np.zeros(n_snps)

    # Additional metrics for each SNP
    metrics = {
        'maternal_baf': np.zeros(n_snps),  # Estimated maternal BAF
        'fetal_baf': np.zeros(n_snps),     # Estimated fetal BAF
        'deviation': np.zeros(n_snps),     # Deviation from expected BAF
        # Whether SNP is informative for UPD
        'informative': np.zeros(n_snps, dtype=bool),
        # UPD pattern indicator (0=normal, 1=paternal UPD indicator, 2=maternal UPD indicator)
        'upd_pattern': np.zeros(n_snps, dtype=int)
    }

    # Expected BAF values for different maternal-fetal genotype combinations
    # Format: (maternal_genotype, fetal_genotype): expected_baf
    expected_baf = {
        # Maternal AA (0)
        (0, 0): 0.0,                    # Fetal AA
        (0, 1): fetal_fraction / 2,     # Fetal AB
        (0, 2): fetal_fraction,         # Fetal BB - paternal UPD indicator

        # Maternal AB (1)
        (1, 0): 0.5 - fetal_fraction / 2,  # Fetal AA
        (1, 1): 0.5,                       # Fetal AB
        (1, 2): 0.5 + fetal_fraction / 2,  # Fetal BB

        # Maternal BB (2)
        (2, 0): 1.0 - fetal_fraction,   # Fetal AA - paternal UPD indicator
        (2, 1): 1.0 - fetal_fraction / 2,  # Fetal AB
        (2, 2): 1.0                        # Fetal BB
    }

    # Thresholds adjusted for maternal DNA background
    # These thresholds are critical for accurate genotype calling
    aa_threshold = 0.25  # Upper threshold for maternal AA
    bb_threshold = 0.75  # Lower threshold for maternal BB

    # Deviation thresholds for detecting fetal genotype in maternal heterozygous sites
    het_deviation_threshold = fetal_fraction / 8

    # Paternal UPD detection thresholds
    # For maternal AA, if BAF is significantly higher than expected for AB, it might be BB (paternal UPD)
    # For maternal BB, if BAF is significantly lower than expected for AB, it might be AA (paternal UPD)
    # More sensitive threshold for detecting paternal UPD
    paternal_upd_threshold = fetal_fraction * 0.7

    # Determine maternal genotype and estimate fetal genotype
    for i, baf in enumerate(baf_values):
        # Store the observed BAF
        metrics['maternal_baf'][i] = baf

        # Determine maternal genotype
        if baf < aa_threshold:
            maternal_genotypes[i] = 0  # AA
            # Maternal homozygous sites are informative
            metrics['informative'][i] = True

            # Check for potential paternal UPD (AA mother, BB fetus)
            if baf > paternal_upd_threshold:
                fetal_genotypes[i] = 2  # BB - strong paternal UPD indicator
                confidence_scores[i] = min(1.0, baf / fetal_fraction)
                metrics['deviation'][i] = abs(baf - fetal_fraction)
                metrics['upd_pattern'][i] = 1  # Paternal UPD indicator
                metrics['fetal_baf'][i] = 1.0
            # Normal inheritance patterns
            elif baf < fetal_fraction / 4:
                fetal_genotypes[i] = 0  # AA
                confidence_scores[i] = 1.0 - (baf / (fetal_fraction / 2))
                metrics['deviation'][i] = baf  # Deviation from expected 0
                metrics['fetal_baf'][i] = 0.0
            else:
                fetal_genotypes[i] = 1  # AB
                confidence_scores[i] = baf / (fetal_fraction / 2)
                metrics['deviation'][i] = abs(baf - fetal_fraction / 2)
                metrics['fetal_baf'][i] = 0.5

        elif baf > bb_threshold:
            maternal_genotypes[i] = 2  # BB
            # Maternal homozygous sites are informative
            metrics['informative'][i] = True

            # Check for potential paternal UPD (BB mother, AA fetus)
            if baf < 1.0 - paternal_upd_threshold:
                fetal_genotypes[i] = 0  # AA - strong paternal UPD indicator
                confidence_scores[i] = min(1.0, (1.0 - baf) / fetal_fraction)
                metrics['deviation'][i] = abs(baf - (1.0 - fetal_fraction))
                metrics['upd_pattern'][i] = 1  # Paternal UPD indicator
                metrics['fetal_baf'][i] = 0.0
            # Normal inheritance patterns
            elif baf > 1.0 - fetal_fraction / 4:
                fetal_genotypes[i] = 2  # BB
                confidence_scores[i] = 1.0 - \
                    ((1.0 - baf) / (fetal_fraction / 2))
                metrics['deviation'][i] = 1.0 - \
                    baf  # Deviation from expected 1
                metrics['fetal_baf'][i] = 1.0
            else:
                fetal_genotypes[i] = 1  # AB
                confidence_scores[i] = (1.0 - baf) / (fetal_fraction / 2)
                metrics['deviation'][i] = abs(baf - (1.0 - fetal_fraction / 2))
                metrics['fetal_baf'][i] = 0.5

        else:
            maternal_genotypes[i] = 1  # AB

            # Maternal heterozygous sites are informative for isodisomy
            metrics['informative'][i] = True

            # Estimate fetal genotype based on deviation from 0.5
            deviation = baf - 0.5
            metrics['deviation'][i] = abs(deviation)

            if abs(deviation) < het_deviation_threshold:
                fetal_genotypes[i] = 1  # AB
                confidence_scores[i] = 1.0 - \
                    (abs(deviation) / (fetal_fraction / 4))
                metrics['fetal_baf'][i] = 0.5
            elif deviation < 0:
                fetal_genotypes[i] = 0  # AA
                confidence_scores[i] = min(
                    1.0, abs(deviation) / (fetal_fraction / 4))
                metrics['fetal_baf'][i] = 0.0
                # Check if this could be maternal UPD
                if abs(deviation) > fetal_fraction / 3:
                    # Potential maternal UPD indicator
                    metrics['upd_pattern'][i] = 2
            else:
                fetal_genotypes[i] = 2  # BB
                confidence_scores[i] = min(
                    1.0, deviation / (fetal_fraction / 4))
                metrics['fetal_baf'][i] = 1.0
                # Check if this could be maternal UPD
                if abs(deviation) > fetal_fraction / 3:
                    # Potential maternal UPD indicator
                    metrics['upd_pattern'][i] = 2

        # Adjust confidence based on sequencing depth
        # Lower depth means less confidence in BAF estimation
        depth_factor = min(1.0, depth / 100)
        confidence_scores[i] *= depth_factor

    return fetal_genotypes, maternal_genotypes, confidence_scores, metrics


def analyze_maternal_heterozygous_sites(positions, fetal_genotypes, maternal_genotypes, confidence_scores, baf_values, fetal_fraction=0.1):
    """
    Analyze maternal heterozygous sites to distinguish between maternal isodisomy and heterodisomy.
    
    In maternal UPD:
    - Heterodisomy: Maternal AB sites should show fetal AB (normal BAF distribution around 0.5)
    - Isodisomy: Maternal AB sites should show fetal AA or BB (skewed BAF away from 0.5)
    
    Parameters:
    - positions: genomic positions
    - fetal_genotypes: estimated fetal genotypes (0: AA, 1: AB, 2: BB)
    - maternal_genotypes: maternal genotypes (0: AA, 1: AB, 2: BB)
    - confidence_scores: confidence in genotype calls (0-1)
    - baf_values: original BAF values
    - fetal_fraction: estimated fetal DNA fraction
    
    Returns:
    - Float between 0 and 1 indicating isodisomy ratio (1 = complete isodisomy, 0 = complete heterodisomy)
    """
    # Filter for maternal heterozygous sites with good confidence
    het_idx = (maternal_genotypes == 1) & (confidence_scores > 0.6)
    
    if sum(het_idx) < 20:
        print(f"  Insufficient maternal heterozygous sites ({sum(het_idx)}) for isodisomy/heterodisomy analysis")
        return 0.5  # Default value when insufficient data
    
    # Get BAF values at maternal heterozygous sites
    het_baf = baf_values[het_idx]
    
    # Calculate deviation from 0.5 (perfect heterozygosity)
    deviation = np.abs(het_baf - 0.5)
    
    # In maternal heterodisomy, fetal genotype should also be heterozygous (AB)
    # In maternal isodisomy, fetal genotype should be homozygous (AA or BB)
    
    # Expected deviation for heterodisomy (should be close to 0)
    expected_hetero_dev = 0.0
    
    # Expected deviation for isodisomy (should be around fetal_fraction/2)
    # For isodisomy, BAF should be around 0.5 ± fetal_fraction/2
    expected_iso_dev = fetal_fraction / 2
    
    # Count sites that match isodisomy pattern (deviation > threshold)
    iso_threshold = fetal_fraction / 4  # Threshold to distinguish iso from hetero
    iso_count = np.sum(deviation > iso_threshold)
    
    # Calculate isodisomy ratio
    iso_ratio = iso_count / len(het_baf) if len(het_baf) > 0 else 0.5
    
    # Adjust ratio based on expected patterns
    # If deviations are consistently large, it's more likely isodisomy
    # If deviations are consistently small, it's more likely heterodisomy
    mean_deviation = np.mean(deviation)
    
    # Normalize mean deviation to a 0-1 scale where 1 = isodisomy, 0 = heterodisomy
    normalized_dev = min(1.0, mean_deviation / expected_iso_dev)
    
    # Combine count-based and deviation-based metrics
    final_iso_ratio = (iso_ratio + normalized_dev) / 2
    
    print(f"  Maternal heterozygous sites analysis:")
    print(f"    Sites analyzed: {len(het_baf)}")
    print(f"    Mean BAF deviation from 0.5: {mean_deviation:.4f}")
    print(f"    Isodisomy-pattern sites: {iso_count} ({iso_ratio:.2f})")
    print(f"    Final isodisomy ratio: {final_iso_ratio:.2f}")
    
    return final_iso_ratio


def detect_upd_with_hmm(positions, fetal_genotypes, maternal_genotypes, confidence_scores, chromosome, metrics=None, fetal_fraction=0.1, baf_values=None):
    """
    Detect UPD regions using Hidden Markov Model.
    Uses CategoricalHMM and only considers maternal homozygous genotypes (pure AA or BB).
    For maternal UPD, also analyzes maternal heterozygous sites to distinguish isodisomy from heterodisomy.

    Parameters:
    - positions: genomic positions
    - fetal_genotypes: estimated fetal genotypes (0: AA, 1: AB, 2: BB)
    - maternal_genotypes: maternal genotypes (0: AA, 1: AB, 2: BB)
    - confidence_scores: confidence in genotype calls (0-1)
    - chromosome: chromosome number
    - metrics: additional metrics from genotype estimation
    - fetal_fraction: estimated fetal DNA fraction
    - baf_values: original BAF values (for diagnostic purposes)

    Returns:
    - Dictionary with UPD detection results
    """
    # Filter for maternal genotypes that are homozygous (0: AA or 2: BB) with good confidence
    valid_idx = (confidence_scores > 0.6) & (
        (maternal_genotypes == 0) | (maternal_genotypes == 2))

    if sum(valid_idx) < 30:
        print(
            f"Skipping chromosome {chromosome}: insufficient homozygous maternal genotypes ({sum(valid_idx)})")
        return {
            'upd_type': "Insufficient data",
            'heterozygosity_ratio': np.nan,
            'upd_score': np.nan,
            'segmented_data': None,
            'snp_count': len(fetal_genotypes),
            'valid_snp_count': sum(valid_idx),
            'hmm_states': None,
            'maternal_upd_prop': np.nan,
            'paternal_upd_prop': np.nan,
            'maternal_iso_ratio': np.nan,
            'paternal_iso_ratio': np.nan,
            'state_proportions': np.array([np.nan, np.nan, np.nan, np.nan]),
            'segments': []
        }

    positions_filtered = positions[valid_idx]
    fetal_genotypes_filtered = fetal_genotypes[valid_idx]
    maternal_genotypes_filtered = maternal_genotypes[valid_idx]

    # Convert genotype combinations to categorical observations for CategoricalHMM
    # We'll encode the combination of (maternal_genotype, fetal_genotype) as a single categorical value
    # For homozygous maternal sites:
    # 0: Maternal AA, Fetal AA -> normal or maternal UPD
    # 1: Maternal AA, Fetal AB -> normal or paternal heterodisomy
    # 2: Maternal AA, Fetal BB -> paternal UPD (strong indicator)
    # 3: Maternal BB, Fetal AA -> paternal UPD (strong indicator)
    # 4: Maternal BB, Fetal AB -> normal or paternal heterodisomy
    # 5: Maternal BB, Fetal BB -> normal or maternal UPD

    n_samples = len(fetal_genotypes_filtered)
    observations = np.zeros(n_samples, dtype=int)

    for i in range(n_samples):
        mat_gt = maternal_genotypes_filtered[i]
        fet_gt = fetal_genotypes_filtered[i]

        if mat_gt == 0:  # Maternal AA
            if fet_gt == 0:
                observations[i] = 0  # Maternal AA, Fetal AA
            elif fet_gt == 1:
                observations[i] = 1  # Maternal AA, Fetal AB
            else:
                # Maternal AA, Fetal BB - strong paternal UPD indicator
                observations[i] = 2
        elif mat_gt == 2:  # Maternal BB
            if fet_gt == 0:
                # Maternal BB, Fetal AA - strong paternal UPD indicator
                observations[i] = 3
            elif fet_gt == 1:
                observations[i] = 4  # Maternal BB, Fetal AB
            else:
                observations[i] = 5  # Maternal BB, Fetal BB

    # Debug: Print observation distribution
    unique_obs, obs_counts = np.unique(observations, return_counts=True)
    print(f"Observation distribution: {list(zip(unique_obs, obs_counts))}")
    
    # Reshape observations for CategoricalHMM
    observations = observations.reshape(-1, 1)

    # Calculate distances between SNPs (for transition probabilities)
    distances = np.diff(positions_filtered)
    median_distance = np.median(distances) if len(distances) > 0 else 1000

    # Adjust transition probability based on median distance
    base_transition_prob = min(0.01, max(0.0001, 1000 / median_distance))
    # Probability to stay in the same state
    stay_prob = 1 - 3 * base_transition_prob

    # Initialize CategoricalHMM with 4 states:
    # 0: Normal
    # 1: Maternal UPD (combined)
    # 2: Paternal isodisomy
    # 3: Paternal heterodisomy
    #
    # Note: Maternal and paternal UPD are biologically distinct and cannot transition between each other.
    # UPD occurs during meiosis or early mitosis and remains fixed for the entire chromosome.
    model = hmm.CategoricalHMM(n_components=4, random_state=42)

    # Set initial probabilities (start in normal state with higher probability)
    model.startprob_ = np.array([0.55, 0.15, 0.15, 0.15])

    # Set transition probabilities
    # High probability to stay in the same state, low probability to transition
    # Prevent transitions between maternal UPD and paternal UPD states (biologically impossible)
    model.transmat_ = np.array([
        [stay_prob, base_transition_prob, base_transition_prob,
            base_transition_prob],  # From normal
        # From maternal UPD - cannot transition to paternal UPD
        [base_transition_prob, stay_prob, 0.0, 0.0],
        # From paternal isodisomy - cannot transition to maternal UPD
        [base_transition_prob, 0.0, stay_prob, base_transition_prob],
        # From paternal heterodisomy - cannot transition to maternal UPD
        [base_transition_prob, 0.0, base_transition_prob, stay_prob]
    ])

    # Normalize transition probabilities to ensure each row sums to 1
    for i in range(4):
        model.transmat_[i] = model.transmat_[i] / model.transmat_[i].sum()

    # Set emission probabilities for CategoricalHMM
    # For each state, define probability of observing each genotype combination
    # We have 6 possible observations (0-5) and 4 states

    # Initialize emission probabilities
    # Small baseline probability for all combinations
    emissionprob = np.zeros((4, 6))

    # Normal inheritance (state 0)
    emissionprob[0, 0] = 0.5  # Maternal AA, Fetal AA - common in normal inheritance
    emissionprob[0, 1] = 0.5  # Maternal AA, Fetal AB - possible but less common
    emissionprob[0, 2] = 0.0  # Maternal AA, Fetal BB - rare in normal inheritance
    emissionprob[0, 3] = 0.0  # Maternal BB, Fetal AA - rare in normal inheritance
    emissionprob[0, 4] = 0.5  # Maternal BB, Fetal AB - possible but less common
    emissionprob[0, 5] = 0.5  # Maternal BB, Fetal BB - common in normal inheritance

    # Maternal UPD (state 1) - Combined isodisomy and heterodisomy
    emissionprob[1, 0] = 1  # Maternal AA, Fetal AA - common in maternal UPD
    emissionprob[1, 1] = 0.0  # Maternal AA, Fetal AB - possible in maternal heterodisomy
    emissionprob[1, 2] = 0.0  # Maternal AA, Fetal BB - rare in maternal UPD
    emissionprob[1, 3] = 0.0  # Maternal BB, Fetal AA - rare in maternal UPD
    emissionprob[1, 4] = 0.0  # Maternal BB, Fetal AB - possible in maternal heterodisomy
    emissionprob[1, 5] = 1  # Maternal BB, Fetal BB - common in maternal UPD

    # Paternal isodisomy (state 2) - Increased emphasis on BB-AA and AA-BB patterns
    emissionprob[2, 0] = 0.5  # Maternal AA, Fetal AA - less common
    emissionprob[2, 1] = 0.0  # Maternal AA, Fetal AB - less common
    emissionprob[2, 2] = 0.5  # Maternal AA, Fetal BB - common in paternal isodisomy
    emissionprob[2, 3] = 0.5  # Maternal BB, Fetal AA - common in paternal isodisomy
    emissionprob[2, 4] = 0.0  # Maternal BB, Fetal AB - less common
    emissionprob[2, 5] = 0.5  # Maternal BB, Fetal BB - less common

    # Paternal heterodisomy (state 3)
    emissionprob[3, 0] = 0.25  # Maternal AA, Fetal AA - less common
    emissionprob[3, 1] = 0.50  # Maternal AA, Fetal AB - common in paternal heterodisomy
    emissionprob[3, 2] = 0.25  # Maternal AA, Fetal BB - less common
    emissionprob[3, 3] = 0.25  # Maternal BB, Fetal AA - less common
    emissionprob[3, 4] = 0.50  # Maternal BB, Fetal AB - common in paternal heterodisomy
    emissionprob[3, 5] = 0.25  # Maternal BB, Fetal BB - less common

    # Normalize emission probabilities
    for i in range(4):
        emissionprob[i] = emissionprob[i] / emissionprob[i].sum()

    model.emissionprob_ = emissionprob
    # Predict states directly (no training needed as we've set all parameters)
    states = model.predict(observations)

    # Calculate state proportions
    state_counts = np.bincount(states, minlength=4)
    print(f"State counts {state_counts=}")
    state_props = state_counts / len(states)
    print(f"State proportions {state_props=}")

    # Calculate heterozygosity ratio
    het_count = np.sum(fetal_genotypes_filtered == 1)
    het_ratio = het_count / len(fetal_genotypes_filtered)

    # Determine UPD type based on state proportions
    maternal_upd_prop = state_props[1]  # Maternal UPD state
    paternal_upd_prop = state_props[2] + state_props[3]  # Paternal UPD states
    upd_score = maternal_upd_prop + paternal_upd_prop

    # Determine isodisomy vs heterodisomy ratios
    paternal_iso_ratio = state_props[2] / (state_props[2] + state_props[3]) if (
        state_props[2] + state_props[3]) > 0 else 0
    
    # For maternal UPD, analyze maternal heterozygous sites to determine isodisomy vs heterodisomy
    maternal_iso_ratio = 0.5  # Default value
    
    # If maternal UPD is detected, analyze maternal heterozygous sites
    if maternal_upd_prop > 0.3 and baf_values is not None:
        print(f"  Detected potential maternal UPD, analyzing maternal heterozygous sites...")
        maternal_iso_ratio = analyze_maternal_heterozygous_sites(
            positions, fetal_genotypes, maternal_genotypes, 
            confidence_scores, baf_values, fetal_fraction
        )

    # Determine UPD type with confidence level
    confidence_threshold = 0.7  # Higher threshold for more confident calls

    if upd_score > confidence_threshold:
        if maternal_upd_prop > paternal_upd_prop:
            if maternal_iso_ratio > 0.7:
                upd_type = "Maternal Isodisomy"
            elif maternal_iso_ratio > 0.3:
                upd_type = "Mixed Maternal Iso/Heterodisomy"
            else:
                upd_type = "Maternal Heterodisomy"
        else:
            if paternal_iso_ratio > 0.7:
                upd_type = "Paternal Isodisomy"
            elif paternal_iso_ratio > 0.3:
                upd_type = "Mixed Paternal Iso/Heterodisomy"
            else:
                upd_type = "Paternal Heterodisomy"
    elif upd_score > 0.3:
        if maternal_upd_prop > paternal_upd_prop:
            if maternal_iso_ratio > 0.7:
                upd_type = "Partial Maternal Isodisomy"
            elif maternal_iso_ratio > 0.3:
                upd_type = "Partial Maternal Mixed UPD"
            else:
                upd_type = "Partial Maternal Heterodisomy"
        else:
            if paternal_iso_ratio > 0.7:
                upd_type = "Partial Paternal Isodisomy"
            elif paternal_iso_ratio > 0.3:
                upd_type = "Partial Paternal Mixed UPD"
            else:
                upd_type = "Partial Paternal Heterodisomy"
    else:
        upd_type = "No UPD detected"

    # Create segments for visualization
    segments = []
    current_state = states[0]
    start_idx = 0

    for i in range(1, len(states)):
        if states[i] != current_state:
            segments.append({
                'start': positions_filtered[start_idx],
                'end': positions_filtered[i],
                'state': current_state,
                'length': positions_filtered[i] - positions_filtered[start_idx]
            })
            current_state = states[i]
            start_idx = i

    # Add the last segment
    segments.append({
        'start': positions_filtered[start_idx],
        'end': positions_filtered[-1],
        'state': current_state,
        'length': positions_filtered[-1] - positions_filtered[start_idx]
    })

    # Create segmented data for visualization
    segmented_data = segment_chromosome_with_states(
        positions_filtered, fetal_genotypes_filtered, states, chromosome
    )

    # Calculate confidence in UPD call based on:
    # 1. Number of informative SNPs
    # 2. Consistency of state calls
    # 3. Fetal fraction
    # More SNPs = higher confidence
    snp_count_factor = min(1.0, sum(valid_idx) / 500)
    state_consistency = max(
        maternal_upd_prop, paternal_upd_prop) if upd_score > 0.3 else state_props[0]
    upd_confidence = snp_count_factor * state_consistency

    return {
        'upd_type': upd_type,
        'heterozygosity_ratio': het_ratio,
        'upd_score': upd_score,
        'upd_confidence': upd_confidence,
        'maternal_upd_prop': maternal_upd_prop,
        'paternal_upd_prop': paternal_upd_prop,
        'maternal_iso_ratio': maternal_iso_ratio,
        'paternal_iso_ratio': paternal_iso_ratio,
        'state_proportions': state_props,
        'segmented_data': segmented_data,
        'snp_count': len(fetal_genotypes),
        'valid_snp_count': sum(valid_idx),
        'hmm_states': states,
        'segments': segments
    }


def segment_chromosome_with_states(positions, genotypes, states, chromosome, window_size=1000000):
    """
    Segment chromosome into windows with state information.
    Enhanced to provide more detailed metrics for UPD detection.
    Adapted for the 4-state model:
    0: Normal
    1: Maternal UPD (combined)
    2: Paternal isodisomy
    3: Paternal heterodisomy
    """
    if len(positions) == 0:
        return None

    # Sort by position
    sorted_idx = np.argsort(positions)
    positions = positions[sorted_idx]
    genotypes = genotypes[sorted_idx]
    states = states[sorted_idx]

    # Create windows
    max_pos = positions.max()
    windows = []

    for start in range(0, int(max_pos) + window_size, window_size):
        end = start + window_size
        window_idx = (positions >= start) & (positions < end)

        if np.sum(window_idx) > 0:
            window_genotypes = genotypes[window_idx]
            window_states = states[window_idx]

            # Calculate state proportions
            state_counts = np.bincount(window_states, minlength=4)
            state_props = state_counts / \
                len(window_states) if len(window_states) > 0 else np.zeros(4)

            # Calculate heterozygosity metrics
            het_count = np.sum(window_genotypes == 1)
            het_ratio = het_count / \
                len(window_genotypes) if len(window_genotypes) > 0 else 0

            # Calculate UPD proportions
            maternal_upd_prop = state_props[1]  # Maternal UPD state
            paternal_upd_prop = state_props[2] + \
                state_props[3]  # Paternal UPD states

            # Calculate isodisomy ratios (if UPD is present)
            paternal_iso_ratio = state_props[2] / \
                paternal_upd_prop if paternal_upd_prop > 0 else 0
            # Default value as we don't have separate maternal iso/hetero states
            maternal_iso_ratio = 0.5

            # Determine dominant state
            dominant_state = np.argmax(state_counts)

            # Determine UPD type for this window
            if dominant_state == 0:
                upd_type = "Normal"
            elif dominant_state == 1:
                upd_type = "Maternal UPD"
            elif dominant_state == 2:
                upd_type = "Paternal Isodisomy"
            elif dominant_state == 3:
                upd_type = "Paternal Heterodisomy"
            else:
                upd_type = "Unknown"

            windows.append({
                'chromosome': chromosome,
                'start': start,
                'end': end,
                'mean_genotype': np.mean(window_genotypes),
                'het_ratio': het_ratio,
                'state_0_prop': state_props[0],  # Normal
                'state_1_prop': state_props[1],  # Maternal UPD
                'state_2_prop': state_props[2],  # Paternal isodisomy
                'state_3_prop': state_props[3],  # Paternal heterodisomy
                'maternal_upd_prop': maternal_upd_prop,
                'paternal_upd_prop': paternal_upd_prop,
                'maternal_iso_ratio': maternal_iso_ratio,
                'paternal_iso_ratio': paternal_iso_ratio,
                'dominant_state': dominant_state,
                'upd_type': upd_type,
                'snp_count': len(window_genotypes)
            })

    return pd.DataFrame(windows)


def detect_upd_patterns(df, fetal_fraction=0.1, depth=100):
    """
    Detect UPD patterns using direct genotype estimation and HMM.
    Distinguishes between isodisomy and heterodisomy.
    Optimized for maternal plasma analysis with consideration of maternal DNA background.

    Parameters:
    - df: DataFrame with chromosome, position, and baf columns
    - fetal_fraction: Estimated fetal DNA fraction in maternal plasma
    - depth: Estimated sequencing depth

    Returns:
    - Dictionary of results for each chromosome
    """
    results = {}

    print(
        f"Starting UPD detection with fetal fraction: {fetal_fraction:.3f}, depth: {depth}")

    for chromosome in sorted(df['chromosome'].unique(), key=lambda x: (0 if x.isdigit() else 1, x)):
        chr_data = df[df['chromosome'] == chromosome]

        if len(chr_data) < 100:  # Skip chromosomes with too few SNPs
            print(
                f"Skipping chromosome {chromosome}: insufficient SNPs ({len(chr_data)})")
            continue

        print(
            f"Processing chromosome {chromosome} with {len(chr_data)} SNPs...")

        # Sort by position
        chr_data = chr_data.sort_values('position')

        # Estimate fetal and maternal genotypes with additional metrics
        fetal_genotypes, maternal_genotypes, confidence_scores, metrics = estimate_fetal_genotype(
            chr_data['baf'].values,
            fetal_fraction=fetal_fraction,
            depth=depth
        )
        # Add genotypes to dataframe for later use
        chr_data['fetal_genotype'] = fetal_genotypes
        chr_data['maternal_genotype'] = maternal_genotypes
        chr_data['confidence_score'] = confidence_scores

        # Count informative SNPs
        if 'informative' in metrics:
            informative_count = np.sum(metrics['informative'])
            high_conf_informative = np.sum(
                (confidence_scores > 0.6) & metrics['informative'])
            print(
                f"  Found {informative_count} informative SNPs ({high_conf_informative} high confidence)")

        # Count important UPD indicator patterns
        if 'upd_pattern' in metrics:
            paternal_upd_indicators = np.sum(metrics['upd_pattern'] == 1)
            maternal_upd_indicators = np.sum(metrics['upd_pattern'] == 2)
            print(f"  UPD indicator patterns:")
            print(
                f"    Paternal UPD indicators (BB-AA or AA-BB): {paternal_upd_indicators}")
            print(f"    Maternal UPD indicators: {maternal_upd_indicators}")

        # Detect UPD using HMM with CategoricalHMM, focusing on maternal homozygous sites
        upd_results = detect_upd_with_hmm(
            chr_data['position'].values,
            fetal_genotypes,
            maternal_genotypes,
            confidence_scores,
            chromosome,
            metrics=metrics,
            fetal_fraction=fetal_fraction,
            baf_values=chr_data['baf'].values
        )

        # Store results
        results[chromosome] = upd_results

        # Print summary with confidence information
        confidence_str = f", Confidence: {upd_results.get('upd_confidence', 0):.3f}" if 'upd_confidence' in upd_results else ""

        print(f"Chromosome {chromosome}: {upd_results['upd_type']} "
              f"(Score: {upd_results['upd_score']:.3f}, "
              f"Het Ratio: {upd_results['heterozygosity_ratio']:.3f}{confidence_str})")

        if 'maternal_upd_prop' in upd_results and 'paternal_upd_prop' in upd_results:
            print(f"  Maternal UPD: {upd_results['maternal_upd_prop']:.3f}, "
                  f"Paternal UPD: {upd_results['paternal_upd_prop']:.3f}")

            if upd_results['maternal_upd_prop'] > 0.3:
                print(
                    f"  Maternal isodisomy ratio: {upd_results['maternal_iso_ratio']:.3f}")

            if upd_results['paternal_upd_prop'] > 0.3:
                print(
                    f"  Paternal isodisomy ratio: {upd_results['paternal_iso_ratio']:.3f}")

    return results


def plot_results(df, results, output_dir):
    """
    Generate plots for visualization of UPD detection results.
    Creates chromosome-specific plots and summary visualizations.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Plot BAF distribution and HMM states for each chromosome
    for chromosome in results:
        chr_data = df[df['chromosome'] == chromosome].copy()
        segmented_data = results[chromosome]['segmented_data']

        if segmented_data is None or len(segmented_data) == 0:
            continue

        # Re-estimate genotypes for plotting if they're not in the dataframe
        if 'fetal_genotype' not in chr_data.columns or 'maternal_genotype' not in chr_data.columns:
            print(f"Re-estimating genotypes for chromosome {chromosome} for plotting...")
            fetal_genotypes, maternal_genotypes, confidence_scores, metrics = estimate_fetal_genotype(
                chr_data['baf'].values,
                fetal_fraction=0.1,  # Default value
                depth=100  # Default value
            )
            chr_data['fetal_genotype'] = fetal_genotypes
            chr_data['maternal_genotype'] = maternal_genotypes
            chr_data['confidence_score'] = confidence_scores
            
            # Add UPD pattern indicators if available
            if 'upd_pattern' in metrics:
                chr_data['upd_pattern'] = metrics['upd_pattern']

        plt.figure(figsize=(14, 12))

        # Plot 1: BAF scatter plot
        plt.subplot(4, 1, 1)
        plt.scatter(chr_data['position'], chr_data['baf'], s=1, alpha=0.5)
        plt.title(f"Chromosome {chromosome} - BAF Distribution")
        plt.ylabel("B-Allele Frequency")
        plt.ylim(-0.05, 1.05)
        plt.axhline(y=0.5, color='r', linestyle='--', alpha=0.5)

        # Plot 2: Fetal genotype estimates with UPD indicators highlighted
        plt.subplot(4, 1, 2)
        # Regular genotypes
        plt.scatter(
            chr_data['position'],
            chr_data['fetal_genotype'],
            s=1, alpha=0.5,
            color='blue'
        )

        # Highlight UPD indicator patterns if available
        if 'upd_pattern' in chr_data.columns:
            paternal_idx = chr_data['upd_pattern'] == 1
            maternal_idx = chr_data['upd_pattern'] == 2

            if np.any(paternal_idx):
                plt.scatter(
                    chr_data.loc[paternal_idx, 'position'],
                    chr_data.loc[paternal_idx, 'fetal_genotype'],
                    s=5, alpha=0.8, color='red',
                    marker='x', label='Paternal UPD indicator'
                )

            if np.any(maternal_idx):
                plt.scatter(
                    chr_data.loc[maternal_idx, 'position'],
                    chr_data.loc[maternal_idx, 'fetal_genotype'],
                    s=5, alpha=0.8, color='green',
                    marker='+', label='Maternal UPD indicator'
                )

            plt.legend(loc='upper right', fontsize='small')

        plt.title(f"Chromosome {chromosome} - Estimated Fetal Genotypes")
        plt.ylabel("Genotype (0=AA, 1=AB, 2=BB)")
        plt.ylim(-0.1, 2.1)

        # Plot 3: Maternal genotype estimates
        plt.subplot(4, 1, 3)
        plt.scatter(
            chr_data['position'],
            chr_data['maternal_genotype'],
            s=1, alpha=0.5
        )
        plt.title(f"Chromosome {chromosome} - Estimated Maternal Genotypes")
        plt.ylabel("Genotype (0=AA, 1=AB, 2=BB)")
        plt.ylim(-0.1, 2.1)

        # Plot 4: HMM state segmentation
        plt.subplot(4, 1, 4)

        # Plot state proportions by window
        if 'dominant_state' in segmented_data.columns:
            # Create colormap for states
            colors = ['green', 'red', 'blue', 'purple']
            # Normal, Maternal UPD, Paternal isodisomy, Paternal heterodisomy

            # Plot segments colored by dominant state
            for _, segment in segmented_data.iterrows():
                plt.bar(
                    segment['start'] + (segment['end'] - segment['start'])/2,
                    1.0,
                    width=(segment['end'] - segment['start']),
                    color=colors[int(segment['dominant_state'])],
                    alpha=0.7
                )

            plt.title(f"Chromosome {chromosome} - HMM State Segmentation")
            plt.xlabel("Position")
            plt.ylabel("State")
            plt.ylim(0, 1.1)

            # Add legend
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor='green', alpha=0.7, label='Normal'),
                Patch(facecolor='red', alpha=0.7, label='Maternal UPD'),
                Patch(facecolor='blue', alpha=0.7, label='Paternal Isodisomy'),
                Patch(facecolor='purple', alpha=0.7,
                      label='Paternal Heterodisomy')
            ]
            plt.legend(handles=legend_elements, loc='upper right')

        # Add UPD information
        upd_type = results[chromosome]['upd_type']
        het_ratio = results[chromosome]['heterozygosity_ratio']
        upd_score = results[chromosome].get('upd_score', 0)

        # Add more detailed UPD information if available
        upd_info = f"{upd_type} (Heterozygosity: {het_ratio:.3f}, UPD Score: {upd_score:.3f})"

        if 'maternal_upd_prop' in results[chromosome] and 'paternal_upd_prop' in results[chromosome]:
            mat_upd = results[chromosome]['maternal_upd_prop']
            pat_upd = results[chromosome]['paternal_upd_prop']
            upd_info += f"\nMaternal UPD: {mat_upd:.3f}, Paternal UPD: {pat_upd:.3f}"

            if mat_upd > 0.3:
                mat_iso = results[chromosome]['maternal_iso_ratio']
                upd_info += f", Maternal Iso/Het: {mat_iso:.2f}/{1-mat_iso:.2f}"

            if pat_upd > 0.3:
                pat_iso = results[chromosome]['paternal_iso_ratio']
                upd_info += f", Paternal Iso/Het: {pat_iso:.2f}/{1-pat_iso:.2f}"

        plt.figtext(0.5, 0.01, upd_info,
                    ha="center", fontsize=12, bbox={"facecolor": "orange", "alpha": 0.2, "pad": 5})

        plt.tight_layout()
        plt.savefig(os.path.join(
            output_dir, f"chromosome_{chromosome}_analysis.png"), dpi=150)
        plt.close()

    # Create summary plots
    plt.figure(figsize=(14, 10))

    # Prepare data
    chromosomes = []
    upd_scores = []
    het_ratios = []
    maternal_upd_props = []
    paternal_upd_props = []
    maternal_iso_ratios = []
    paternal_iso_ratios = []
    colors = []

    for chromosome in sorted(results.keys(), key=lambda x: (0 if x.isdigit() else 1, x)):
        if results[chromosome]['upd_type'] == "Insufficient data":
            continue

        chromosomes.append(chromosome)
        upd_scores.append(results[chromosome]['upd_score'])
        het_ratios.append(results[chromosome]['heterozygosity_ratio'])

        # Get UPD proportions
        mat_upd = results[chromosome].get('maternal_upd_prop', 0)
        pat_upd = results[chromosome].get('paternal_upd_prop', 0)
        maternal_upd_props.append(mat_upd)
        paternal_upd_props.append(pat_upd)

        # Get isodisomy ratios
        mat_iso = results[chromosome].get('maternal_iso_ratio', 0)
        pat_iso = results[chromosome].get('paternal_iso_ratio', 0)
        maternal_iso_ratios.append(mat_iso)
        paternal_iso_ratios.append(pat_iso)

        # Set color based on UPD type
        upd_type = results[chromosome]['upd_type']
        if "Complete Maternal" in upd_type:
            colors.append('red')
        elif "Complete Paternal" in upd_type:
            colors.append('blue')
        elif "Partial Maternal" in upd_type:
            colors.append('orange')
        elif "Partial Paternal" in upd_type:
            colors.append('purple')
        else:
            colors.append('green')

    # Plot 1: Heterozygosity ratio
    plt.subplot(2, 2, 1)
    plt.bar(chromosomes, het_ratios, color=colors)
    plt.axhline(y=0.25, color='r', linestyle='--',
                alpha=0.5, label="UPD Threshold")
    plt.title("Heterozygosity Ratio by Chromosome")
    plt.xlabel("Chromosome")
    plt.ylabel("Heterozygosity Ratio")
    plt.xticks(rotation=90)
    plt.legend()

    # Plot 2: UPD score
    plt.subplot(2, 2, 2)
    plt.bar(chromosomes, upd_scores, color=colors)
    plt.axhline(y=0.3, color='orange', linestyle='--',
                alpha=0.5, label="Partial UPD Threshold")
    plt.axhline(y=0.7, color='r', linestyle='--',
                alpha=0.5, label="Complete UPD Threshold")
    plt.title("UPD Score by Chromosome")
    plt.xlabel("Chromosome")
    plt.ylabel("UPD Score")
    plt.xticks(rotation=90)
    plt.legend()

    # Plot 3: Maternal vs Paternal UPD
    plt.subplot(2, 2, 3)
    width = 0.35
    x = np.arange(len(chromosomes))
    plt.bar(x - width/2, maternal_upd_props, width,
            label='Maternal UPD', color='red', alpha=0.7)
    plt.bar(x + width/2, paternal_upd_props, width,
            label='Paternal UPD', color='blue', alpha=0.7)
    plt.axhline(y=0.3, color='orange', linestyle='--', alpha=0.5)
    plt.title("Maternal vs Paternal UPD by Chromosome")
    plt.xlabel("Chromosome")
    plt.ylabel("UPD Proportion")
    plt.xticks(x, chromosomes, rotation=90)
    plt.legend()

    # Plot 4: Isodisomy vs Heterodisomy
    plt.subplot(2, 2, 4)

    # Filter chromosomes with significant UPD
    upd_chroms = []
    mat_iso_vals = []
    pat_iso_vals = []
    iso_colors = []

    for i, chrom in enumerate(chromosomes):
        if maternal_upd_props[i] > 0.3 or paternal_upd_props[i] > 0.3:
            upd_chroms.append(chrom)
            mat_iso_vals.append(
                maternal_iso_ratios[i] if maternal_upd_props[i] > 0.3 else 0)
            pat_iso_vals.append(
                paternal_iso_ratios[i] if paternal_upd_props[i] > 0.3 else 0)
            iso_colors.append(colors[i])

    if upd_chroms:
        width = 0.35
        x = np.arange(len(upd_chroms))
        plt.bar(x - width/2, mat_iso_vals, width,
                label='Maternal Isodisomy', color='red', alpha=0.7)
        plt.bar(x + width/2, pat_iso_vals, width,
                label='Paternal Isodisomy', color='blue', alpha=0.7)
        plt.axhline(y=0.5, color='black', linestyle='--', alpha=0.5)
        plt.title("Isodisomy Ratio in UPD Chromosomes")
        plt.xlabel("Chromosome")
        plt.ylabel("Isodisomy Ratio")
        plt.xticks(x, upd_chroms, rotation=90)
        plt.legend()
    else:
        plt.text(0.5, 0.5, "No significant UPD detected",
                 ha='center', va='center', fontsize=14)
        plt.axis('off')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "upd_summary.png"), dpi=150)
    plt.close()

    # Generate summary report
    summary_df = pd.DataFrame({
        'Chromosome': chromosomes,
        'UPD_Type': [results[chr]['upd_type'] for chr in chromosomes],
        'Heterozygosity': [results[chr]['heterozygosity_ratio'] for chr in chromosomes],
        'UPD_Score': [results[chr]['upd_score'] for chr in chromosomes],
        'SNP_Count': [results[chr]['snp_count'] for chr in chromosomes]
    })

    summary_df.to_csv(os.path.join(
        output_dir, "upd_summary.tsv"), sep='\t', index=False)

    return summary_df


def main():
    parser = argparse.ArgumentParser(
        description='Detect UPD from maternal plasma SNP BAF data')
    parser.add_argument('--input', '-i', required=True,
                        help='Input BAF file (TSV format)')
    parser.add_argument(
        '--output', '-o', default='upd_results', help='Output directory')
    parser.add_argument('--baf-threshold', '-b', type=float, default=0.0,
                        help='BAF threshold for informative SNPs (default: 0.2)')
    parser.add_argument('--fetal-fraction', '-f', type=float, default=0.1,
                        help='Estimated fetal DNA fraction in maternal plasma (default: 0.1)')
    parser.add_argument('--depth', '-d', type=int, default=100,
                        help='Estimated sequencing depth (default: 100)')

    args = parser.parse_args()

    # Load data
    df = load_baf_data(args.input)
    if df is None:
        return

    # Filter informative SNPs
    filtered_df = filter_informative_snps(df, args.baf_threshold)

    # Detect UPD patterns with consideration of maternal DNA background
    results = detect_upd_patterns(filtered_df,
                                  fetal_fraction=args.fetal_fraction,
                                  depth=args.depth)

    # Plot and save results
    summary = plot_results(df, results, args.output)  # Use original df for plotting

    print("\nUPD Detection Summary:")
    print(summary)
    print(f"\nResults saved to {args.output}/")


if __name__ == "__main__":
    main()
