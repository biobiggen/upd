# UPD Calculation Tool

A standalone UPD (Uniparental Disomy) analysis tool for NIPT SNP readlists.

[中文文档](README.zh.md)

## Features

- **Self-contained**: only 5 third-party libraries required
- **Robust fetal fraction**: dual-track FF estimation prevents missed calls in
  genome-wide homozygous samples caused by FF overestimation
- **Batch parallelism**: whole batches processed in parallel across processes

## System Requirements

### Software dependencies

| Component | Version | Notes |
|---|---|---|
| Python | 3.9 – 3.11 | Tested on 3.9 / 3.10 / 3.11 |
| numpy | >=1.20 | Numerical computation |
| pandas | >=1.3 | readlist parsing and group statistics |
| scipy | >=1.7 | `scipy.stats.binom` likelihood |
| hmmlearn | >=0.3.0 | Must provide `CategoricalHMM` (absent in 0.2.x) |
| matplotlib | >=3.4 | Required only for the plotting feature |

Pinned versions are listed in `requirements.txt`.

### Operating systems

- Linux (CentOS 7 / Ubuntu 20.04, 22.04): the production platform
- Windows 10/11, macOS 12+: can run the demo in this document; plotting forces
  the `matplotlib` `Agg` backend, so no display is needed

### Hardware

No special hardware and no GPU. BLAS thread count is pinned to 1 internally;
the `batch` subcommand parallelises per sample across processes. Peak memory
per sample is a few hundred MB (at the 13833-probe scale).

## Installation

```bash
git clone https://github.com/biobiggen/upd
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r upd/requirements.txt
```

The repository root **is** the Python package `upd`, so all commands below must
be run from the directory that *contains* the cloned `upd/` directory (the
modules use relative imports). Do not `cd upd` first.

Pure Python, no compilation needed. Installation time on a normal desktop is
dominated by downloading the numpy/pandas/scipy wheels, usually a few minutes.

## Demo

Real plasma data cannot be distributed publicly, so a simulated data generator
is provided; its output format is identical to production.

### 1. Generate simulated demo data

```bash
python -m upd.simulate_demo_data -o demo_data
```

Generated content (defaults: fetal fraction 0.10, mean depth 1200x, 120 SNPs
per target region, roughly 1200 sites in total):

| File | Description |
|---|---|
| `demo_data/demo_probe.bed` | Probe file (`#Chr` / `Pos` / `Type`) |
| `demo_data/demoN01..N04_normal_*.reads.list` | 4 normal biparental disomy samples |
| `demo_data/demoP01_updm_*.reads.list` | `15q11q13` maternal UPD |
| `demo_data/demoP02_updpi_*.reads.list` | `11p15` paternal heterodisomy |
| `demo_data/demoP03_updpii_*.reads.list` | `11p15` paternal isodisomy |

readlist column format:

```
#Chr  Pos  Ref  Alt  Depth  Ref_Dep  Alt_Dep  pA_Ratio  GC  DepRegion  MapQ
```

### 2. Run a single sample

```bash
python -m upd.cli single \
    -i demo_data/demoP01_updm_consensus.mapped.clipped.snp.reads.list \
    -o demo_out/demoP01_upd_results.json \
    -p demo_data/demo_probe.bed \
    --plot demo_out/demoP01_upd_regions_pa_ratio.png
```

Expected output:

- stdout prints the `final_state` and ratio of each region plus fetal fraction
  diagnostics
- `demo_out/demoP01_upd_results.json`: full per-region results (see
  "Output Format")
- `demo_out/demoP01_upd_regions_pa_ratio.png`: pA_Ratio scatter plot

For this simulated sample the `final_state` of `15q11q13` is expected to be
`UPDM` with all other regions `Normal`, and `FF_Method` to be `homozygous`
(ratio ~ 1).

> Note: the simulated data has far fewer sites than the real probe panel
> (13833 rows), so do not call `validate_row_count()`. `--probe-version` only
> affects `ignore_snps` and thresholds; the default `NIPT3V4` is fine for the
> demo.

### 3. Batch processing and summary report

```bash
python -m upd.cli batch \
    -i demo_data/ -o demo_out/ -p demo_data/demo_probe.bed \
    --threads 4 --plot-dir demo_out/UPD_image/ -v

python -m upd.cli report --results demo_out/ -o demo_out/upd_report.csv
```

Expected output: one `*_upd_results.json` per sample under `demo_out/`, one PNG
per sample under `demo_out/UPD_image/`, and the aggregated
`demo_out/upd_report.csv` (one row per sample, two columns per region: state
and ratio).

### Expected runtime

On a normal desktop (4-core CPU, 16 GB RAM) the whole demo is expected to
finish within a few minutes: data generation takes seconds; single-sample
runtime is dominated by the per-site loop in
`predict_fetal_genotype_hypergeom`, a few seconds at demo scale and about 1–2
minutes for the real 13833-row panel; `batch` scales linearly with
`--threads`.

## Usage

### Single sample

```bash
python -m upd.cli single \
    -i sample_consensus.mapped.clipped.snp.reads.list \
    -o result.json \
    -p /path/to/NIPT3V4_CNV-HAP_targets_hg38.bed \
    --probe-version NIPT3V4
```

Also generate the scatter plot:

```bash
python -m upd.cli single -i sample.readslist -o result.json \
    -p probe.bed --plot sample_upd_regions_pa_ratio.png
```

### Batch processing

```bash
python -m upd.cli batch \
    -i ./readslist_dir/ \
    -o ./results/ \
    -p /path/to/probe.bed \
    --threads 8 \
    --plot-dir ./UPD_image/ \
    -v
```

Arguments:

| Argument | Description |
|---|---|
| `-i, --input` | Directory containing `*.snp.reads.list` files |
| `-o, --output` | Output directory for JSON results |
| `-p, --probe-file` | Path to the probe file |
| `--probe-version` | `NIPT3V3` / `NIPT3V4` (default `NIPT3V4`) |
| `--threads` | Number of parallel processes (default 4) |
| `--plot-dir` | Optional, output directory for scatter plots |
| `--no-recursive` | Do not scan subdirectories recursively |
| `-v, --verbose` | Print the result of every sample |

### Generating the report

```bash
# Scan a directory
python -m upd.cli report --results ./results/ -o upd_report.csv

# Or pass an explicit file list
python -m upd.cli report --upd-jsons a.json b.json -o upd_report.csv
```

## Running on Your Own Data

1. **Prepare readlists**: the upstream alignment pipeline produces one
   `*.snp.reads.list` (tab-separated) per sample. Required columns are `#Chr`,
   `Pos`, `Ref`, `Alt`, `Depth`, `pA_Ratio`. If a `DepRegionGC` column exists it
   is preferred for the chrY depth ratio, otherwise `Depth` is used. The first
   column must be `#Chr`.
2. **Prepare the probe file**: `.bed` / `.xls` (tab-separated) or `.csv`, with
   required columns `#Chr` and `Pos` plus a type column (common names such as
   `Type` / `SNP_Tag` / `Probe_Type` are recognised). Values in the type column
   must match the region names in `regions.py` (`6q24`, `7q32`, `11p15`,
   `14q32`, `15q11q13`, `20q13`, `chrRef`). Sites absent from the probe file are
   labelled `chrRef`.
3. **Choose the probe version**: `--probe-version NIPT3V3` or `NIPT3V4`
   determines `expected_rows`, `ff_threshold`, `depth_threshold` and the ignored
   probe types (`HBA` / `RHD` / `SMN` / `HAP` / `other`). For a different
   upstream panel, add an entry to `core.PROBE_VERSIONS`.
4. **Change the target regions**: if the imprinted regions of interest differ
   from the defaults, edit `UPD_TARGETS_HMM` (narrow intervals for the HMM) and
   `UPD_TARGETS_PLOT` (wide intervals for plotting) in `regions.py`, and keep
   `REPORT_REGIONS` and `PLOT_ORDER` in sync.
5. **Run in batch and aggregate**: run `batch` first, then `report` (see the
   commands above). Keep each batch of samples in one directory and do not set
   `--threads` above the physical core count.
6. **Tunable parameters**: decision thresholds are collected as constants at the
   top of `core.py` (`MIN_REGION_SNPS`, `MIN_HOMOZYGOUS_SNPS`,
   `MIN_OBSERVATIONS`, `MIN_REGION_LENGTH`, `DEPTH_FILTER`,
   `NORMAL_RATIO_THRESHOLD`, `SIGNIFICANT_UPD_RATIO`); HMM parameters are
   `STARTPROB` / `TRANSMAT` / `EMISSIONPROB`. `DEPTH_FILTER` (default 400)
   gates HMM region filtering and may be lowered for low-depth samples; the
   fetal-fraction estimate uses a separate `FF_DEPTH_THRESHOLD` (default 450).
7. **Quality control**: UPD calls are unreliable when the fetal fraction is
   below `ff_threshold` (0.03); the `status` field of each region in the JSON
   explains why a region could not be computed (see "Output Format").

## Programming Interface

```python
from upd import UPDCalculator

calc = UPDCalculator(probe_version='NIPT3V4', probe_file='probe.bed')
calc.load_readlist('sample.snp.reads.list')

ff = calc.get_fetal_fraction()      # fetal fraction (dual-track corrected)
results = calc.predict_upd_hmm()    # UPD prediction
print(calc.ff_info)                 # FF diagnostics (both estimates and ratio)
```

Convenience function:

```python
from upd.core import calculate_upd

results = calculate_upd('sample.readslist', probe_file='probe.bed')

# When FF diagnostics are needed
results = calculate_upd('sample.readslist', probe_file='probe.bed',
                        with_ff_info=True)
print(results['_ff_info'])
```

## Module Layout

The repository root is the package itself; the files below sit directly in it.

| File | Description |
|---|---|
| `core.py` | `UPDCalculator` class, core UPD computation |
| `regions.py` | UPD target region coordinates |
| `plotting.py` | pA_Ratio scatter plotting |
| `cli.py` | Command-line entry point |
| `__main__.py` | Makes `python -m upd` equivalent to `python -m upd.cli` |
| `simulate_demo_data.py` | Small simulated demo dataset generator |
| `requirements.txt` | Dependency list |

## Algorithm

### Computation flow

```
load readlist
  -> get_probe_type       annotate SNP_Tag
  -> get_background_gt    infer maternal genotype
  -> get_fetal_fraction   fetal fraction (dual-track estimate + auto switch)
  -> predict_fetal_genotype_hypergeom   binomial-likelihood fetal genotype
  -> predict_upd_hmm      HMM state prediction
```

### Dual-track fetal fraction estimation

The fetal fraction (FF) is estimated along two mutually independent paths and
cross-checked:

**Path 1 — maternal homozygous sites** (`ffAB`)

At maternal BB sites the expected plasma alt ratio is
`AF = ff * fetal_alt_dosage / 2`. Taking `ff = 2 * median(AF)` implicitly
assumes **the fetus is heterozygous at that site** (i.e. it inherited one alt
allele from the father, dosage 1).

**Path 2 — maternal heterozygous (BA) sites**

At maternal BA sites:

```
AF = (1 - ff) * 0.5 + ff * fetal_alt_dosage / 2
   = 0.5 + (ff / 2) * (fetal_alt_dosage - 1)
```

that is, `|AF - 0.5| = (ff / 2) * |fetal_alt_dosage - 1|`. When the fetus is
homozygous (dosage 0 or 2) the offset is always `ff / 2`, so FF equals twice
the median offset. **This relation holds regardless of whether the fetal genome
comes from both parents or entirely from one parent**, hence it is robust to
genome-wide homozygosity.

**Why genome-wide homozygous samples require path 2**

If the fetus is homozygous genome-wide, it is also homozygous at maternal
homozygous sites (`AF = ff` instead of `ff/2`), so path 1 **overestimates FF by
about 2x**. With FF doubled, `expected_alt_ratio(BB->BA) = ff_est/2` matches the
observed AF exactly, so a fetal **AA is misclassified as BA** and the
observation `BBAA` degenerates into `BBBA`. Since `BBAA`/`AABB` is the only
signature of paternal isodisomy in the emission matrix (probability 0 under
`Normal`), losing it makes the region be called `Normal` and paternal
isodisomy undetectable.

Therefore, when `ff_hom / ff_het > 1.5` (`FF_RATIO_HOM_UPPER`), genome-wide
homozygosity is suspected and the path 2 estimate is used automatically. Normal
samples agree between both paths (ratio ~ 1) and are unaffected; local UPD
(homozygosity in a single region only) does not trigger the switch either,
because FF is a genome-wide statistic.

Both estimates and their ratio are written to the `_ff_info` field of the
results for review.

### HMM states

| State | Meaning |
|---|---|
| `Normal` | Normal biparental disomy |
| `UPDM` | Maternal uniparental disomy |
| `UPDPI` | Paternal uniparental disomy (heterodisomy) |
| `UPDPII` | Paternal uniparental disomy (isodisomy) |

### Target regions

Two sets of region coordinates are used:

- **HMM computation** (narrow intervals): `6q24`, `7q32`, `11p15`, `14q32`,
  `15q11q13`, `20q13`, `chrRef`
- **Plotting** (wide intervals): the first 6 regions (excluding `chrRef`)

### Decision thresholds

| Threshold | Value | Description |
|---|---|---|
| Minimum SNPs per region | 50 | Below this yields `insufficient_data` |
| Minimum maternal homozygous sites | 20 | Below this yields `insufficient_homozygous` |
| Depth filter | 400 | Only sites with depth >=400 are used |
| Minimum segment length | 20 | Shorter segments are discarded |
| Normal decision threshold | 0.2 | Normal ratio >0.2 is called normal |
| Significant UPD ratio | 0.1 | Non-Normal segments with ratio >0.1 enter `significant_upds` |
| FF homozygous-path depth | 450 | Depth floor (`FF_DEPTH_THRESHOLD`) for the maternal homozygous FF estimate; separate from the 400 `DEPTH_FILTER` used for HMM region filtering |
| Minimum sites for the BA path | 30 | Below this no FF cross-check is done |
| FF ratio switch threshold | 1.5 | `ff_hom/ff_het` above this switches to the BA estimate |
| BA offset floor | 0.02 | `\|AF-0.5\|` below this is treated as a heterozygous fetus and excluded from FF estimation |
| BA AF window | 0.25 | `HET_BAND`: maternal BA sites must have `pA_Ratio` within `0.5 +/- 0.25` |

## Output Format

### JSON result

```json
{
  "6q24": {
    "status": "success",
    "state_ratios": {"Normal": 0.95, "UPDM": 0.05},
    "regions": [
      {"start": 142448249, "end": 145502506, "state": "Normal",
       "length": 120, "ratio": 0.95}
    ],
    "significant_upds": [],
    "final_state": "Normal",
    "final_ratio": 0.95,
    "total_observations": 126,
    "chromosome": "chr6",
    "start_pos": 142448249,
    "end_pos": 145502506,
    "site_details": [...]
  },
  "_ff_info": {
    "ff_used": 0.1012,
    "ff_homozygous": 0.1008,
    "ff_heterozygous": 0.1015,
    "ff_ratio": 0.99,
    "ff_method": "homozygous",
    "het_sites": 412,
    "het_shifted_sites": 208
  }
}
```

Besides the regions, the top level contains an `_ff_info` key (underscore
prefix to distinguish it from region names) holding fetal fraction
diagnostics:

| Field | Description |
|---|---|
| `ff_used` | FF actually used downstream |
| `ff_homozygous` | Estimate from the maternal homozygous path |
| `ff_heterozygous` | Estimate from the maternal heterozygous (BA) path; `null` when sites are insufficient |
| `ff_ratio` | `ff_homozygous / ff_heterozygous`; close to 2 suggests genome-wide homozygosity |
| `ff_method` | `homozygous` (default) or `heterozygous` (correction applied) |
| `het_sites` | Number of maternal BA sites passing the depth filter |
| `het_shifted_sites` | Of those, the number above the offset floor used for FF estimation |

Possible `status` values:

| Value | Description |
|---|---|
| `success` | Computation succeeded |
| `insufficient_data` | Fewer than 50 SNPs in the region |
| `insufficient_homozygous` | Fewer than 20 maternal homozygous sites |
| `insufficient_observations` | Fewer than 20 valid observations |
| `no_valid_observations` | No observation could be mapped |
| `insufficient_observation_types` | Fewer than 2 observation types |
| `observation_conversion_error` | Observation conversion failed |
| `hmm_error` | HMM prediction failed |

### CSV report

After the two columns per region (state, ratio), three FF diagnostic columns
are appended:

```csv
Sample_Name,6q24_UPD,6q24_Ratio,...,chrRef_Ratio,FF_Used,FF_Method,FF_Ratio
A346_cf1439_2516654P1DE,Normal,0.9500,...,0.9600,0.1012,homozygous,0.99
demoP03,UPDPII,0.9800,...,0.9700,0.1015,homozygous,1.01
```

Samples with `FF_Method` equal to `heterozygous` had the FF correction
triggered; together with `FF_Ratio` (close to 2) they allow quick screening of
samples suspected of genome-wide homozygosity.

## License

This software is released under the **GNU General Public License v3.0** (an
OSI-approved open source licence); the full text is in the `COPYING` file at the
repository root.

Copyright (C) 2024 biobiggen

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but **WITHOUT
ANY WARRANTY**; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>.

### Commercial licensing

GPL-3.0 is a copyleft licence: derivative works based on this code must also be
distributed under GPL-3.0. To integrate this software into a closed-source
product, or for licensing terms not bound by copyleft, please contact the author
to discuss a separate commercial licence.

A complete description of the algorithm (pseudocode) is given in the
**Methods** section of the paper.
