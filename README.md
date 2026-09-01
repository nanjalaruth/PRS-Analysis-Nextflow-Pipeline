# PRS Analysis — a Nextflow pipeline for computing polygenic scores from the PGS Catalog

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.10.0-brightgreen.svg)](https://www.nextflow.io/)

A [Nextflow](https://www.nextflow.io) pipeline that computes **polygenic scores (PGS)** for many traits and many PGS Catalog score IDs at once. Give it a list of traits, the PGS IDs you want for each trait, and your genotype data in PLINK format; it downloads the scoring files from the [PGS Catalog](https://www.pgscatalog.org/), reformats them, runs `plink --score`, and returns one tidy table of scores per trait, ready for downstream modelling (e.g. combining scores with ElasticNet).

Computing scores by hand for dozens of trait × PGS-ID combinations is slow and error-prone. Running it through Nextflow means every combination is processed in parallel, failures are retried, and a run can be resumed with `-resume` instead of starting again.

---

## Contents

- [What the pipeline does](#what-the-pipeline-does)
- [Requirements](#requirements)
- [Installation](#installation)
- [Preparing your input](#preparing-your-input)
- [Running the pipeline](#running-the-pipeline)
- [Parameters](#parameters)
- [Outputs](#outputs)
- [Repository layout](#repository-layout)
- [Support](#support)
- [Citation](#citation)

---

## What the pipeline does

For every **trait × PGS ID** pair:

1. **Downloads the scoring file** from the PGS Catalog FTP site. The pipeline fetches the *harmonized* file with GRCh37 positions (`<PGS_ID>_hmPOS_GRCh37.txt.gz`).
2. **Downloads the score metadata** (`<PGS_ID>_metadata_scores.csv`) and extracts the reported trait, PGS ID and original genome build.
3. **Strips the comment header** (`#` lines) from the scoring file.
4. **Builds a clean three-column scoring file** (`SNP`, `A1`, `BETA`) in R: variants are identified as `chr:pos` from the harmonized `hm_chr`/`hm_pos` columns, duplicates and missing values are dropped.
5. **Computes scores with PLINK 1.9** (`plink --score <file> 1 2 3 sum`) against your genotype data.
6. **Extracts `FID` and `SCORE`** from the PLINK output and renames the score column to `<PGS_ID>_SCORE`.

Then, for every **dataset × trait**:

7. **Merges all PGS-ID scores** for that trait into a single table (one row per individual, one column per PGS ID).

Two further steps — score integration with phenotype and covariate data, and association analysis with forest plots — are implemented in `main.nf` and `templates/` but currently commented out of the workflow. A liftover branch (hg19 ↔ hg38, using UCSC `liftOver` and SNP-info reference files) is likewise available but disabled, since the harmonized GRCh37 files make it unnecessary in most cases. See `extra_README.md` for the inputs those steps need.

![pipeline](https://github.com/nanjalaruth/PRS-Analysis-Nextflow-Pipeline/blob/main/output/pipeline_info/output.png)

---

## Requirements

| Tool | Notes |
|------|-------|
| [Nextflow](https://www.nextflow.io) ≥ 20.10.0 | Requires Java 11+ |
| [PLINK 1.9](https://www.cog-genomics.org/plink/) | Must be on your `PATH` as `plink` |
| R (≥ 4.x) | With `pacman`, `readr`, `dplyr`, `data.table` (installed automatically on first run via `pacman`) |
| `curl`, `gunzip` | For fetching PGS Catalog files — the compute nodes need internet access |

---

## Installation

### 1. Nextflow

```bash
wget -qO- https://get.nextflow.io | bash
# move the `nextflow` binary somewhere on your PATH
```

### 2. PLINK 1.9

Download a binary from https://www.cog-genomics.org/plink/ and put it on your `PATH`, or:

```bash
conda install -c bioconda plink
```

### 3. R

If you do not already have R, the easiest route is conda:

```bash
# install Miniconda if needed: https://docs.conda.io/en/latest/miniconda.html
conda install -c conda-forge r-base
```

The R packages the templates use are installed automatically on first run.

### 4. The pipeline

Nothing to install: Nextflow fetches the pipeline from GitHub the first time you run it. To update to the latest version later:

```bash
nextflow pull nanjalaruth/PRS-Analysis-Nextflow-Pipeline
```

---

## Preparing your input

### 1. PGS ID lists (one file per trait)

Create a directory (e.g. `all_blood_traits_prs_scores/`) containing one text file per trait, named `<trait>_PGS_score_ids.txt`. Each file lists the PGS Catalog IDs to compute for that trait — **one ID per line, no header**:

```
PGS003940
PGS004727
PGS004728
```

The `<trait>` prefix must match the names you list in `params.bloodCells` (see below). For example, for basophil count you might use `baso`, giving `baso_PGS_score_ids.txt`.

You can find PGS IDs by searching the [PGS Catalog](https://www.pgscatalog.org/) for your trait.

### 2. Genotype data

Your genotypes in PLINK binary format (`.bed`, `.bim`, `.fam`). Two things to check:

- **Genome build:** the scoring files are downloaded in **GRCh37/hg19** coordinates, so your genotypes should be on GRCh37 as well.
- **Variant IDs:** the pipeline matches variants by **`chr:pos`** (e.g. `10:76684698`, no `chr` prefix). The second column of your `.bim` file must use this format, otherwise PLINK will find no overlapping variants. If your `.bim` uses rsIDs, re-label it first, e.g.:

  ```bash
  awk 'BEGIN{OFS="\t"} {$2=$1":"$4; print}' data.bim > data.chrpos.bim
  ```

---

## Running the pipeline

The pipeline is configured through a config file rather than command-line flags. Copy `conf/test.config` and edit these lines:

```groovy
// Trait names — must match the <trait>_PGS_score_ids.txt file names
params.bloodCells = ["baso", "rbc", "wbc"]

// Directory holding the <trait>_PGS_score_ids.txt files
params.basePath = "/path/to/all_blood_traits_prs_scores"

params {
    // Genotype data: [dataset label, bed, bim, fam]
    plink_file = [
        ['UGRC', '/path/to/uganda.bed', '/path/to/uganda.bim', '/path/to/uganda.fam']
    ]
}
```

You can list several datasets in `plink_file`; every trait × PGS ID will be scored in each of them.

Then run:

```bash
nextflow run nanjalaruth/PRS-Analysis-Nextflow-Pipeline \
    -c /path/to/my.config \
    -resume
```

Or, if you have cloned the repository and edited `conf/test.config` in place:

```bash
nextflow run main.nf -c conf/test.config -resume
```

### Running on a cluster

Add `-profile slurm` to submit each task as a SLURM job:

```bash
nextflow run nanjalaruth/PRS-Analysis-Nextflow-Pipeline -profile slurm -c my.config -resume
```

For other schedulers (PBS, SGE, …) add a profile to `nextflow.config` with the appropriate `process.executor`.

---

## Parameters

### Pipeline parameters (set in your config file)

| Parameter | Example | Description |
|-----------|---------|-------------|
| `bloodCells` | `["baso", "rbc", "wbc"]` | List of trait names. Each needs a `<trait>_PGS_score_ids.txt` file in `basePath`. **Required.** |
| `basePath` | `/path/to/all_blood_traits_prs_scores` | Directory containing the PGS ID list files. **Required.** |
| `plink_file` | `[['UGRC', 'x.bed', 'x.bim', 'x.fam']]` | One or more genotype datasets as `[label, bed, bim, fam]`. **Required.** |
| `outdir` | `./output` | Where results are written. Default: `./output`. |

### Nextflow options (`-`)

| Option | Values | Description |
|--------|--------|-------------|
| `-c` | `path/to/file.config` | Your edited configuration file. |
| `-profile` | `standard`, `slurm`, `debug` | Execution profile. `standard` (default) runs locally. |
| `-resume` | — | Reuse cached results from a previous run. |

Resource limits (`max_memory`, `max_cpus`, `max_time`) can be adjusted in `conf/base.config`.

---

## Outputs

All results go under `outdir` (default `./output`):

```
output/
├── <dataset>/
│   └── <trait>/
│       ├── <dataset>_<trait>_<PGS_ID>_prsval.profile   # raw PLINK --score output
│       ├── <dataset>_<trait>_<PGS_ID>_prsval.log
│       └── <dataset>_<trait>_<PGS_ID>_prsval.nosex
├── <trait>/
│   └── <dataset>_<trait>_pgs_scores.txt                # merged scores, one column per PGS ID
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    └── pipeline_dag.dot
```

The main deliverable is **`<trait>/<dataset>_<trait>_pgs_scores.txt`**, a tab-separated table with `FID` followed by one `<PGS_ID>_SCORE` column per score ID:

```
FID        PGS003940_SCORE   PGS004727_SCORE   PGS004728_SCORE
SAMPLE001  0.1234            -0.0456           0.0789
...
```

Examples from a completed run: [analysis output](https://nanje.quarto.pub/intergrated_prs/) and [Nextflow execution report](https://github.com/nanjalaruth/PRS-Analysis-Nextflow-Pipeline/blob/main/output/pipeline_info/execution_report.html).

---

## Repository layout

```
.
├── main.nf                     # pipeline definition (DSL2)
├── nextflow.config             # profiles, resources, reports
├── conf/
│   ├── base.config             # default resource limits
│   └── test.config             # example configuration — copy and edit
├── templates/
│   ├── modify_meta_file.R      # extracts trait / PGS ID / genome build from metadata
│   ├── modify_score_file.R     # builds the SNP / A1 / BETA scoring file
│   ├── merge_scores.R          # merges per-PGS-ID scores into one table per trait
│   ├── intergrate_scores.R     # (optional step) integrates scores with phenotypes/covariates
│   └── prediction.R            # (optional step) association analysis and forest plots
├── extra_README.md             # inputs for the optional liftover branch
└── output/pipeline_info/       # example run reports
```

---

## Support

Questions, bugs and feature requests are tracked in the [GitHub issues](https://github.com/nanjalaruth/PRS-Analysis-Nextflow-Pipeline/issues).

---

## Citation

If you use this pipeline, please cite the PGS Catalog and PLINK:

> Lambert SA, et al. (2021). The Polygenic Score Catalog as an open database for reproducibility and systematic evaluation. *Nature Genetics*, 53, 420–425. https://doi.org/10.1038/s41588-021-00783-5

> Chang CC, et al. (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience*, 4, 7. https://doi.org/10.1186/s13742-015-0047-8
