# Targeted DNA Methylation Sequencing Pipeline

Code accompanying:

> **A Standardized End-to-End Pipeline for Processing Twist Methylation Sequencing Data**
> Giannini H, Longtin A, Cao E, Watowich MM, Pershad Y, Hodges E, Lea A, Bick AG.
> *Current Protocols* (in preparation)

This repository contains the full computational pipeline described in the paper, from raw FASTQ files through cytosine reports, quality-controlled methylation matrices, cell-type deconvolution, epigenetic clock estimation, and methylation risk score (MRS) generation.

## Pipeline overview

The pipeline mirrors the six Basic Protocols (BP1–BP6) in the paper. Each protocol has its own folder so a given section of the paper maps directly to a given folder/script.

| Protocol | What it does | Folder | Environment |
|---|---|---|---|
| BP1 | FASTQ → adapter trimming → Bismark alignment → cytosine (CX) reports | `bp1_fastq_to_cx_report/` | Bash + SLURM |
| BP2 | CX reports → filtered BSseq object | `bp2_bsseq_generation/` | R + SLURM Submission |
| BP3 | CpG-level QC and generation of methylation-level / coverage matrices | `bp3_qc_and_matrices/` | R |
| BP4 | Cell-type deconvolution (EpiDISH) | `bp4_cell_deconvolution/` | R |
| BP5 | Epigenetic clock estimation (methylCIPHER) | `bp5_epigenetic_clocks/` | R |
| BP6 | Methylation risk score (MRS) generation | `bp6_methylation_risk_scores/` | R |

BP1 is designed to run on a SLURM-managed HPC cluster. BP2–BP6 are R scripts that can be run in RStudio, on Terra, or on any machine with the required Bioconductor packages installed — they do not depend on SLURM.

## Repository structure

```
methylome-pipeline/
├── README.md
├── environment/
│   ├── environment.yml        # or renv.lock — R package versions
│   └── slurm_modules.txt      # exact module versions loaded in #SBATCH scripts
├── bp1_fastq_to_cx_report/
│   ├── 00_genome_preparation.sh
│   ├── 01_trimming.sh
│   └── 02_bismark_suite.sh
├── bp2_bsseq_generation/
│   └── bsseq.R
├── bp3_qc_and_matrices/
│   └── qc_matrices.R
├── bp4_cell_deconvolution/
│   └── deconvolution.R
├── bp5_epigenetic_clocks/
│   └── epigenetic_clocks.R
├── bp6_methylation_risk_scores/
│   └── mrs_scoring.R
├── config/
│   └── sample_manifest.tsv    # sample IDs, FASTQ paths, batch info
└── results/                   # pipeline outputs (gitignored if large)
```

## Setup

### HPC / SLURM environment (BP1)

BP1 assumes access to a SLURM-managed cluster with the following modules available:

- Trimmomatic
- Bismark suite
- Bowtie2
- SAMtools

Exact module versions used for this pipeline are listed in `environment/slurm_modules.txt`. **SLURM configuration (partition names, account strings, module names) is site-specific** — adjust the `#SBATCH` headers and `module load` lines in `bp1_fastq_to_cx_report/*.sh` to match your institution's cluster before running.

### R environment (BP2–BP6)

Install R (version X.X used here) and the following Bioconductor/CRAN packages:

- Biostrings
- BSseq
- GenomicRanges
- BSgenome.Hsapiens.UCSC.hg38
- EpiDISH
- methylCIPHER

A full package list with versions is in `environment/environment.yml`. If using Terra, these scripts were originally run in [Terra workspace name / Docker image — fill in].

## Running the pipeline

Steps must be run in order; each protocol's output is the next protocol's input.

1. **BP1 (one-time step):** `sbatch bp1_fastq_to_cx_report/00_genome_preparation.sh` — generates a bisulfite-converted reference genome. Only needs to be run once per reference genome.
2. **BP1:** `sbatch bp1_fastq_to_cx_report/01_trimming.sh` — trims adapters from raw FASTQ files.
3. **BP1:** `sbatch bp1_fastq_to_cx_report/02_bismark_suite.sh` — aligns trimmed reads and generates per-sample CX reports.
4. **BP2:** Run `bp2_bsseq_generation/bsseq.R` in R — reads CX reports and produces a filtered `BSseq` object.
5. **BP3:** Run `bp3_qc_and_matrices/qc_matrices.R` in R — performs CpG-level QC and writes methylation-level and coverage matrices.
6. **BP4:** Run `bp4_cell_deconvolution/deconvolution.R` in R — estimates cell-type proportions with EpiDISH.
7. **BP5:** Run `bp5_epigenetic_clocks/epigenetic_clocks.R` in R — estimates epigenetic age with methylCIPHER.
8. **BP6:** Run `bp6_methylation_risk_scores/mrs_scoring.R` in R — generates methylation risk scores for exposures/diseases.

## Data

- 1 raw FASTQ file (2 paired files) are included in this repository. [Zenodo]
- 8 CX reports are included in this repository. [Zenodo]
- The Twist Methylome manifest (BED file) used for targeted panel filtering is available from [Twist Bioscience](https://www.twistbioscience.com/products/ngs/fixed-panels/human-methylome-panel?tab=resources#accordion-30b78b384b-item-dfa177bd30).
- The Illumina manifest (csv) used for EPIC array filtering is available from [Illumina] (
