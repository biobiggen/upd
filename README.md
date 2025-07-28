# Maternal Plasma Fetal UPD Detection

This tool detects Uniparental Disomy (UPD) in fetal DNA from maternal plasma using SNP B-Allele Frequency (BAF) data.

## Background

Uniparental Disomy (UPD) is a genetic condition where a child inherits both copies of a chromosome from one parent and none from the other parent. This can lead to various genetic disorders depending on the chromosome involved.

This tool analyzes B-Allele Frequency (BAF) patterns from SNP data to detect potential UPD regions in the fetal genome using maternal plasma samples. It uses a beta-binomial model to estimate fetal genotypes and a Hidden Markov Model (HMM) to detect UPD regions.

## Requirements

- Python 3.9+
- Required Python packages:
  - numpy
  - pandas
  - matplotlib
  - scipy
  - scikit-learn
  - hmmlearn
  - seaborn
  - plotly (for interactive visualization)
  - circlify (for circular plots)

Install required packages:

```bash
pip install -r requirements.txt
```

## Input Data Format

The input file should be a tab-separated (TSV) file with the following columns:
- chromosome: Chromosome number (1-22, X, Y)
- position: Genomic position
- baf: B-Allele Frequency (value between 0 and 1)

Example:

```
chromosome	position	baf
1	1000000	0.5
1	1050000	0.0
1	1100000	1.0
...
```

## Usage

```bash
python upd_detection.py --input your_baf_data.tsv --output results_directory
```

### Parameters

- `--input`, `-i`: Input BAF file (TSV format) [required]
- `--output`, `-o`: Output directory (default: 'upd_results')
- `--baf-threshold`, `-b`: BAF threshold for informative SNPs (default: 0.2)
- `--fetal-fraction`, `-f`: Estimated fetal DNA fraction in maternal plasma (default: 0.1)

## Output

The tool generates:

1. BAF distribution plots for each chromosome
2. A summary plot showing heterozygosity ratios across all chromosomes
3. A summary TSV file with UPD detection results

### Interpretation

- **Complete UPD**: Heterozygosity ratio < 0.15
- **Partial UPD**: Heterozygosity ratio between 0.15 and 0.25
- **No UPD**: Heterozygosity ratio > 0.25

## Workflow Description

1. **Data Loading**: Loads SNP BAF data from the input file
2. **Filtering**: Identifies informative SNPs (BAF between threshold and 1-threshold)
3. **Fetal Genotype Estimation**: Uses a beta-binomial model to estimate fetal genotypes from BAF values
4. **UPD Detection with HMM**: Applies a Hidden Markov Model to detect UPD regions based on estimated genotypes
5. **Segmentation**: Divides chromosomes into windows and calculates state probabilities
6. **UPD Scoring**: Calculates heterozygosity ratio and UPD score for each chromosome
7. **Visualization**: Generates plots and summary reports showing BAF values, estimated genotypes, and HMM states

## Example

```bash
python upd_detection.py --input sample_data.tsv --output upd_results
```

## Notes

- The accuracy of UPD detection depends on the quality and coverage of the input BAF data
- Higher coverage leads to more reliable results
- The tool works best with data from high-quality SNP arrays or NGS data