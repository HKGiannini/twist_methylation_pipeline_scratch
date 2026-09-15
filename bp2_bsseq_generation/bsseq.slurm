#!/bin/bash
#SBATCH --mem=100GB
#SBATCH -t 0-24
#SBATCH --cpus-per-task=8

set -euo pipefail

# Resource guidance: sufficient for <200 samples. For thousands of samples,
# memory needs can scale up to ~300GB (see BP2 introduction).

module load GCC/10.2.0
module load OpenMPI/4.0.5
module load R/4.0.5

Rscript bsseq.R

# SUBMIT AS: sbatch bsseq.slurm
