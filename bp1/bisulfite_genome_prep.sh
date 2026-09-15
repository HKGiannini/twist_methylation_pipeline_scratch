# One-time step: only needs to be run once per reference genome. If you've already prepared this genome before, skip this script entirely.

#!/bin/bash
#SBATCH --mem=50GB # assign your own requirements
#SBATCH -t 0-4 # assign your own requirements
#SBATCH --cpus-per-task=8 # assign your own requirements

# Prevents code from silently failing
set -euo pipefail

# Loading bismark 
module load Bismark/0.25.1

# Reference genomes can be downloaded on ENCODE [https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz]. Store the genome builds in your genome_folder. 
genome_folder=/path/to/genome/directory/

# Loading bowtie
module load bowtie/v2.5.4

# Creates Bisulfite_Genome/CT_conversion/ and Bisulfite_Genome/GA_conversion/ from the reference genome in genome_folder. It will be indexed with Bowtie2. You must generate this file before running bismark_suite.sh. It could take up to 4hrs to run.
bismark_genome_preparation --bowtie2 --verbose "${genome_folder}"

echo "Finished preparing genome at ${genome_folder}"

# SUBMIT AS: sbatch genome_preparation.sh
