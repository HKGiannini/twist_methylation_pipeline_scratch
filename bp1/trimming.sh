#!/bin/bash
#SBATCH --mem=50GB
#SBATCH -t 0-4
#SBATCH --array=X-XX          # replace with your actual sample count
#SBATCH --cpus-per-task=8

# Exit immediately on any error, undefined variable, or failed pipe step — to prevent a broken step from silently continuing into downstream files.
set -euo pipefail

threads=8
sample_list=/path/to/directory/sampleID_list.txt
input_path=/path/to/fastq/files
output_path=/path/to/output/directory
adapter_file=/path/to/trimmomatic/adapter/file/adapters.fa
trimmomatic_jar=/path/to/trimmomatic/trimmomatic-0.39.jar

sampleID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${sample_list}")
echo "SampleID: ${sampleID}"

java -jar "${trimmomatic_jar}" PE \
    -threads "${threads}" \
    "${input_path}/${sampleID}_R1.fastq.gz" \
    "${input_path}/${sampleID}_R2.fastq.gz" \
    "${output_path}/${sampleID}_paired_R1.fastq.gz" \
    "${output_path}/${sampleID}_unpaired_R1.fastq.gz" \
    "${output_path}/${sampleID}_paired_R2.fastq.gz" \
    "${output_path}/${sampleID}_unpaired_R2.fastq.gz" \
    ILLUMINACLIP:"${adapter_file}":2:30:10:8:true \
    HEADCROP:3 TRAILING:10 MINLEN:25

echo "Finished trimming ${sampleID}"

# SUBMIT AS: sbatch trimming.sh
# Typically ~1 hour per array task
