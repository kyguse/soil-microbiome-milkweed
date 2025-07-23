#!/bin/bash
#SBATCH --job-name=fastqc_ITS
#SBATCH --output=fastqc_ITS_%j.out
#SBATCH --error=fastqc_ITS_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# Load FastQC module
module load biocontainers
module load fastqc/0.12.1
# Define input and output directories
INPUT_DIR="/anvil/projects/x-bio250086/soil_microbiome/ITS1"
OUTPUT_DIR="/anvil/projects/x-bio250086/ITS1/ITS_analysis/fastqc_results"
# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
# Run FastQC on all .fastq.gz files in the input directory
fastqc -t 4 -o "$OUTPUT_DIR" "${INPUT_DIR}"/*.fastq.gz

