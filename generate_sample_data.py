#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Sample BAF Data Generator for UPD Detection Testing
--------------------------------------------------
This script generates synthetic B-Allele Frequency (BAF) data for testing
the UPD detection workflow, including normal chromosomes and chromosomes with UPD.
"""

import numpy as np
import pandas as pd
import argparse
import os

NOISE_LEVEL = 0.005

def generate_normal_baf(n_snps, fetal_fraction=0.1, noise_level=0.005):
    """
    Generate BAF values for a normal chromosome (no UPD).
    Expected pattern: peaks at 0, 0.5, and 1
    """
    # Generate maternal genotypes (0: AA, 1: AB, 2: BB)
    maternal_ratio = [0.25, 0.5, 0.25]  # 25% AA, 50% AB, 25% BB
    maternal_genotypes = np.random.choice([0, 1, 2], size=n_snps, p=maternal_ratio)
    
    # Generate fetal genotypes based on maternal genotypes
    fetal_genotypes = np.zeros(n_snps, dtype=int)
    
    for i, mat_gt in enumerate(maternal_genotypes):
        if mat_gt == 0:  # Maternal AA
            # Fetal can be AA (if paternal is A) or AB (if paternal is B)
            fetal_genotypes[i] = np.random.choice([0, 1], p=[0.5, 0.5])
        elif mat_gt == 1:  # Maternal AB
            # Fetal can be AA, AB, or BB
            fetal_genotypes[i] = np.random.choice([0, 1, 2], p=[0.25, 0.5, 0.25])
        else:  # Maternal BB
            # Fetal can be AB (if paternal is A) or BB (if paternal is B)
            fetal_genotypes[i] = np.random.choice([1, 2], p=[0.5, 0.5])
    
    # Generate BAF values based on maternal and fetal genotypes
    baf_values = np.zeros(n_snps)
    
    for i in range(n_snps):
        mat_gt = maternal_genotypes[i]
        fet_gt = fetal_genotypes[i]
        
        # Calculate expected BAF based on maternal and fetal genotypes
        if mat_gt == 0:  # Maternal AA
            if fet_gt == 0:  # Fetal AA
                expected_baf = 0.0
            else:  # Fetal AB
                expected_baf = fetal_fraction / 2
        elif mat_gt == 1:  # Maternal AB
            if fet_gt == 0:  # Fetal AA
                expected_baf = 0.5 - fetal_fraction / 2
            elif fet_gt == 1:  # Fetal AB
                expected_baf = 0.5
            else:  # Fetal BB
                expected_baf = 0.5 + fetal_fraction / 2
        else:  # Maternal BB
            if fet_gt == 1:  # Fetal AB
                expected_baf = 1.0 - fetal_fraction / 2
            else:  # Fetal BB
                expected_baf = 1.0
        
        # Add noise
        baf_values[i] = np.clip(expected_baf + np.random.normal(0, noise_level), 0, 1)
    
    return baf_values, maternal_genotypes, fetal_genotypes


def generate_maternal_isodisomy_baf(n_snps, fetal_fraction=0.1, noise_level=0.005):
    """
    Generate BAF values for a chromosome with maternal isodisomy.
    Both chromosomes from mother are identical.
    """
    # Generate maternal genotypes (0: AA, 1: AB, 2: BB)
    maternal_ratio = [0.25, 0.5, 0.25]  # 25% AA, 50% AB, 25% BB
    maternal_genotypes = np.random.choice([0, 1, 2], size=n_snps, p=maternal_ratio)
    
    # For maternal isodisomy, fetal genotype matches maternal homozygous genotypes
    # and is homozygous for maternal heterozygous sites
    fetal_genotypes = np.zeros(n_snps, dtype=int)
    
    for i, mat_gt in enumerate(maternal_genotypes):
        if mat_gt == 0:  # Maternal AA
            fetal_genotypes[i] = 0  # Fetal AA
        elif mat_gt == 1:  # Maternal AB
            # For isodisomy, heterozygous maternal sites become homozygous in fetus
            fetal_genotypes[i] = np.random.choice([0, 2], p=[0.5, 0.5])
        else:  # Maternal BB
            fetal_genotypes[i] = 2  # Fetal BB
    
    # Generate BAF values based on maternal and fetal genotypes
    baf_values = np.zeros(n_snps)
    
    for i in range(n_snps):
        mat_gt = maternal_genotypes[i]
        fet_gt = fetal_genotypes[i]
        
        # Calculate expected BAF based on maternal and fetal genotypes
        if mat_gt == 0:  # Maternal AA
            expected_baf = 0.0  # Fetal must be AA
        elif mat_gt == 1:  # Maternal AB
            if fet_gt == 0:  # Fetal AA
                expected_baf = 0.5 - fetal_fraction / 2
            else:  # Fetal BB
                expected_baf = 0.5 + fetal_fraction / 2
        else:  # Maternal BB
            expected_baf = 1.0  # Fetal must be BB
        
        # Add noise
        baf_values[i] = np.clip(expected_baf + np.random.normal(0, noise_level), 0, 1)
    
    return baf_values, maternal_genotypes, fetal_genotypes


def generate_maternal_heterodisomy_baf(n_snps, fetal_fraction=0.1, noise_level=0.005):
    """
    Generate BAF values for a chromosome with maternal heterodisomy.
    Both chromosomes from mother but different homologs.
    """
    # Generate maternal genotypes (0: AA, 1: AB, 2: BB)
    maternal_ratio = [0.25, 0.5, 0.25]  # 25% AA, 50% AB, 25% BB
    maternal_genotypes = np.random.choice([0, 1, 2], size=n_snps, p=maternal_ratio)
    
    # For maternal heterodisomy, fetal genotype matches maternal homozygous genotypes
    # and remains heterozygous for maternal heterozygous sites
    fetal_genotypes = np.zeros(n_snps, dtype=int)
    
    for i, mat_gt in enumerate(maternal_genotypes):
        if mat_gt == 0:  # Maternal AA
            fetal_genotypes[i] = 0  # Fetal AA
        elif mat_gt == 1:  # Maternal AB
            # For heterodisomy, heterozygous maternal sites remain heterozygous in fetus
            fetal_genotypes[i] = 1  # Fetal AB
        else:  # Maternal BB
            fetal_genotypes[i] = 2  # Fetal BB
    
    # Generate BAF values based on maternal and fetal genotypes
    baf_values = np.zeros(n_snps)
    
    for i in range(n_snps):
        mat_gt = maternal_genotypes[i]
        
        # Calculate expected BAF based on maternal and fetal genotypes
        if mat_gt == 0:  # Maternal AA
            expected_baf = 0.0  # Fetal must be AA
        elif mat_gt == 1:  # Maternal AB
            expected_baf = 0.5  # Fetal must be AB
        else:  # Maternal BB
            expected_baf = 1.0  # Fetal must be BB
        
        # Add noise
        baf_values[i] = np.clip(expected_baf + np.random.normal(0, noise_level), 0, 1)
    
    return baf_values, maternal_genotypes, fetal_genotypes


def generate_paternal_isodisomy_baf(n_snps, fetal_fraction=0.1, noise_level=0.005):
    """
    Generate BAF values for a chromosome with paternal isodisomy.
    Both chromosomes from father are identical.
    """
    # Generate maternal genotypes (0: AA, 1: AB, 2: BB)
    maternal_ratio = [0.25, 0.5, 0.25]  # 25% AA, 50% AB, 25% BB
    maternal_genotypes = np.random.choice([0, 1, 2], size=n_snps, p=maternal_ratio)
    
    # For paternal isodisomy, fetal genotype is homozygous and may differ from maternal
    fetal_genotypes = np.zeros(n_snps, dtype=int)
    
    for i, mat_gt in enumerate(maternal_genotypes):
        # For paternal isodisomy, fetal is either AA or BB with equal probability
        # This creates discordance with maternal genotype
        fetal_genotypes[i] = np.random.choice([0, 2], p=[0.5, 0.5])
    
    # Generate BAF values based on maternal and fetal genotypes
    baf_values = np.zeros(n_snps)
    
    for i in range(n_snps):
        mat_gt = maternal_genotypes[i]
        fet_gt = fetal_genotypes[i]
        
        # Calculate expected BAF based on maternal and fetal genotypes
        if mat_gt == 0:  # Maternal AA
            if fet_gt == 0:  # Fetal AA
                expected_baf = 0.0
            else:  # Fetal BB
                expected_baf = fetal_fraction  # Strong deviation
        elif mat_gt == 1:  # Maternal AB
            if fet_gt == 0:  # Fetal AA
                expected_baf = 0.5 - fetal_fraction / 2
            else:  # Fetal BB
                expected_baf = 0.5 + fetal_fraction / 2
        else:  # Maternal BB
            if fet_gt == 0:  # Fetal AA
                expected_baf = 1.0 - fetal_fraction  # Strong deviation
            else:  # Fetal BB
                expected_baf = 1.0
        
        # Add noise
        baf_values[i] = np.clip(expected_baf + np.random.normal(0, noise_level), 0, 1)
    
    return baf_values, maternal_genotypes, fetal_genotypes


def generate_paternal_heterodisomy_baf(n_snps, fetal_fraction=0.1, noise_level=0.005):
    """
    Generate BAF values for a chromosome with paternal heterodisomy.
    Both chromosomes from father but different homologs.
    """
    # Generate maternal genotypes (0: AA, 1: AB, 2: BB)
    maternal_ratio = [0.25, 0.5, 0.25]  # 25% AA, 50% AB, 25% BB
    maternal_genotypes = np.random.choice([0, 1, 2], size=n_snps, p=maternal_ratio)
    
    # For paternal heterodisomy, fetal genotype is often heterozygous
    fetal_genotypes = np.zeros(n_snps, dtype=int)
    
    for i, mat_gt in enumerate(maternal_genotypes):
        # For paternal heterodisomy, fetal is heterozygous with high probability
        fetal_genotypes[i] = np.random.choice([0, 1, 2], p=[0.25, 0.5, 0.25])
    
    # Generate BAF values based on maternal and fetal genotypes
    baf_values = np.zeros(n_snps)
    
    for i in range(n_snps):
        mat_gt = maternal_genotypes[i]
        fet_gt = fetal_genotypes[i]
        
        # Calculate expected BAF based on maternal and fetal genotypes
        if mat_gt == 0:  # Maternal AA
            if fet_gt == 0:  # Fetal AA
                expected_baf = 0.0
            elif fet_gt == 1:  # Fetal AB
                expected_baf = fetal_fraction / 2
            else:  # Fetal BB
                expected_baf = fetal_fraction
        elif mat_gt == 1:  # Maternal AB
            if fet_gt == 0:  # Fetal AA
                expected_baf = 0.5 - fetal_fraction / 2
            elif fet_gt == 1:  # Fetal AB
                expected_baf = 0.5
            else:  # Fetal BB
                expected_baf = 0.5 + fetal_fraction / 2
        else:  # Maternal BB
            if fet_gt == 0:  # Fetal AA
                expected_baf = 1.0 - fetal_fraction
            elif fet_gt == 1:  # Fetal AB
                expected_baf = 1.0 - fetal_fraction / 2
            else:  # Fetal BB
                expected_baf = 1.0
        
        # Add noise
        baf_values[i] = np.clip(expected_baf + np.random.normal(0, noise_level), 0, 1)
    
    return baf_values, maternal_genotypes, fetal_genotypes


def generate_mixed_upd_baf(n_snps, upd_type, upd_fraction=0.5, fetal_fraction=0.1, noise_level=0.005):
    """
    Generate BAF values for a chromosome with partial UPD.
    
    Parameters:
    - n_snps: Number of SNPs to generate
    - upd_type: Type of UPD ('maternal_iso', 'maternal_hetero', 'paternal_iso', 'paternal_hetero')
    - upd_fraction: Fraction of chromosome with UPD
    - fetal_fraction: Fetal DNA fraction in maternal plasma
    - noise_level: Standard deviation of noise to add
    
    Returns:
    - BAF values
    - Maternal genotypes
    - Fetal genotypes
    """
    # Determine number of SNPs with UPD
    n_upd_snps = int(n_snps * upd_fraction)
    n_normal_snps = n_snps - n_upd_snps
    
    # Generate normal and UPD BAF values
    normal_baf, normal_mat_gt, normal_fet_gt = generate_normal_baf(
        n_normal_snps, fetal_fraction, noise_level
    )
    
    if upd_type == 'maternal_iso':
        upd_baf, upd_mat_gt, upd_fet_gt = generate_maternal_isodisomy_baf(
            n_upd_snps, fetal_fraction, noise_level
        )
    elif upd_type == 'maternal_hetero':
        upd_baf, upd_mat_gt, upd_fet_gt = generate_maternal_heterodisomy_baf(
            n_upd_snps, fetal_fraction, noise_level
        )
    elif upd_type == 'paternal_iso':
        upd_baf, upd_mat_gt, upd_fet_gt = generate_paternal_isodisomy_baf(
            n_upd_snps, fetal_fraction, noise_level
        )
    elif upd_type == 'paternal_hetero':
        upd_baf, upd_mat_gt, upd_fet_gt = generate_paternal_heterodisomy_baf(
            n_upd_snps, fetal_fraction, noise_level
        )
    else:
        raise ValueError(f"Unknown UPD type: {upd_type}")
    
    # Combine normal and UPD values
    baf_values = np.concatenate([normal_baf, upd_baf])
    maternal_genotypes = np.concatenate([normal_mat_gt, upd_mat_gt])
    fetal_genotypes = np.concatenate([normal_fet_gt, upd_fet_gt])
    
    # Shuffle to mix normal and UPD regions
    indices = np.arange(n_snps)
    np.random.shuffle(indices)
    
    return baf_values[indices], maternal_genotypes[indices], fetal_genotypes[indices]


def generate_sample_dataset(output_file, n_snps_per_chr=500, upd_chromosomes=None, fetal_fraction=0.1):
    """
    Generate a complete sample dataset with specified UPD chromosomes.
    
    Parameters:
    - output_file: Path to save the generated data
    - n_snps_per_chr: Number of SNPs per chromosome
    - upd_chromosomes: Dictionary mapping chromosome to UPD type and parameters
    - fetal_fraction: Fetal DNA fraction in maternal plasma
    """
    if upd_chromosomes is None:
        upd_chromosomes = {}
    
    # Define chromosomes
    autosomes = [str(i) for i in range(1, 23)]
    sex_chromosomes = ['X']
    all_chromosomes = autosomes + sex_chromosomes
    
    data = []
    
    for chromosome in all_chromosomes:
        # Determine if this chromosome has UPD
        upd_info = upd_chromosomes.get(chromosome, None)
        
        # Number of SNPs for this chromosome (fewer for sex chromosomes)
        if chromosome in sex_chromosomes:
            chr_n_snps = n_snps_per_chr // 2
        else:
            chr_n_snps = n_snps_per_chr
        
        # Generate positions (evenly spaced across chromosome)
        if chromosome == 'Y':
            max_pos = 57_227_415  # Approximate Y chromosome length
        elif chromosome == 'X':
            max_pos = 156_040_895  # Approximate X chromosome length
        else:
            # Approximate autosome lengths (very rough estimates)
            max_pos = 250_000_000
        
        positions = np.linspace(1, max_pos, chr_n_snps).astype(int)
        
        # Generate BAF values based on UPD type
        if upd_info:
            upd_type = upd_info.get('type', 'normal')
            upd_fraction = upd_info.get('fraction', 1.0)
            
            if upd_type == 'normal' or upd_fraction <= 0:
                baf_values, _, _ = generate_normal_baf(chr_n_snps, fetal_fraction,noise_level=NOISE_LEVEL)
            elif upd_type == 'maternal_iso':
                if upd_fraction >= 1.0:
                    baf_values, _, _ = generate_maternal_isodisomy_baf(chr_n_snps, fetal_fraction,noise_level=NOISE_LEVEL)
                else:
                    baf_values, _, _ = generate_mixed_upd_baf(
                        chr_n_snps, 'maternal_iso', upd_fraction, fetal_fraction,noise_level=NOISE_LEVEL
                    )
            elif upd_type == 'maternal_hetero':
                if upd_fraction >= 1.0:
                    baf_values, _, _ = generate_maternal_heterodisomy_baf(chr_n_snps, fetal_fraction,noise_level=NOISE_LEVEL)
                else:
                    baf_values, _, _ = generate_mixed_upd_baf(
                        chr_n_snps, 'maternal_hetero', upd_fraction, fetal_fraction,noise_level=NOISE_LEVEL
                    )
            elif upd_type == 'paternal_iso':
                if upd_fraction >= 1.0:
                    baf_values, _, _ = generate_paternal_isodisomy_baf(chr_n_snps, fetal_fraction,noise_level=NOISE_LEVEL)
                else:
                    baf_values, _, _ = generate_mixed_upd_baf(
                        chr_n_snps, 'paternal_iso', upd_fraction, fetal_fraction,noise_level=NOISE_LEVEL
                    )
            elif upd_type == 'paternal_hetero':
                if upd_fraction >= 1.0:
                    baf_values, _, _ = generate_paternal_heterodisomy_baf(chr_n_snps, fetal_fraction,noise_level=NOISE_LEVEL)
                else:
                    baf_values, _, _ = generate_mixed_upd_baf(
                        chr_n_snps, 'paternal_hetero', upd_fraction, fetal_fraction,noise_level=NOISE_LEVEL
                    )
            else:
                raise ValueError(f"Unknown UPD type: {upd_type}")
        else:
            baf_values, _, _ = generate_normal_baf(chr_n_snps, fetal_fraction,noise_level=NOISE_LEVEL)
        
        # Add to dataset
        for pos, baf in zip(positions, baf_values):
            data.append({
                'chromosome': chromosome,
                'position': pos,
                'baf': baf
            })
    
    # Create DataFrame and save
    df = pd.DataFrame(data)
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Generated sample dataset with {len(df)} SNPs")
    print(f"UPD chromosomes: {upd_chromosomes}")
    print(f"Saved to {output_file}")


def parse_upd_chromosomes(complete_upd, partial_upd, upd_type='maternal_iso'):
    """
    Parse UPD chromosome specifications.
    
    Parameters:
    - complete_upd: Comma-separated list of chromosomes with complete UPD
    - partial_upd: Comma-separated list of chromosomes with partial UPD
    - upd_type: Type of UPD to use
    
    Returns:
    - Dictionary mapping chromosome to UPD type and parameters
    """
    upd_chromosomes = {}
    
    # Parse complete UPD chromosomes
    for chr_spec in complete_upd.split(','):
        if not chr_spec.strip():
            continue
            
        parts = chr_spec.strip().split(':')
        chromosome = parts[0]
        
        if len(parts) > 1:
            upd_type = parts[1]
        
        upd_chromosomes[chromosome] = {
            'type': upd_type,
            'fraction': 1.0
        }
    
    # Parse partial UPD chromosomes
    for chr_spec in partial_upd.split(','):
        if not chr_spec.strip():
            continue
            
        parts = chr_spec.strip().split(':')
        chromosome = parts[0]
        
        if len(parts) > 1:
            upd_type = parts[1]
            
        if len(parts) > 2:
            try:
                fraction = float(parts[2])
            except ValueError:
                fraction = 0.5
        else:
            fraction = 0.5
        
        upd_chromosomes[chromosome] = {
            'type': upd_type,
            'fraction': fraction
        }
    
    return upd_chromosomes


def main():
    parser = argparse.ArgumentParser(description='Generate sample BAF data for UPD detection testing')
    parser.add_argument('--output', '-o', default='sample_data.tsv', help='Output file path')
    parser.add_argument('--snps', '-n', type=int, default=1000, help='Number of SNPs per chromosome')
    parser.add_argument('--complete-upd', type=str, default='7:maternal_iso,15:paternal_iso', 
                        help='Comma-separated list of chromosomes with complete UPD (format: chr:type)')
    parser.add_argument('--partial-upd', type=str, default='9:maternal_hetero:0.5,11:paternal_hetero:0.7', 
                        help='Comma-separated list of chromosomes with partial UPD (format: chr:type:fraction)')
    parser.add_argument('--fetal-fraction', type=float, default=0.1, help='Fetal DNA fraction')
    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')
    
    args = parser.parse_args()
    
    # Set random seed
    np.random.seed(args.seed)
    
    # Parse UPD chromosomes
    upd_chromosomes = parse_upd_chromosomes(args.complete_upd, args.partial_upd)
    
    # Generate dataset
    generate_sample_dataset(args.output, args.snps, upd_chromosomes, args.fetal_fraction)


if __name__ == "__main__":
    main()