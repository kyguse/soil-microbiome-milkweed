#!/bin/bash
#SBATCH --job-name=fastqc_16S
#SBATCH --output=fastqc_16S_%j.out
#SBATCH --error=fastqc_16S_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# Load FastQC module
module load fastqc/0.12.1

# Define input and output directories
INPUT_DIR="/anvil/projects/x-bio250086/soil_microbiome/16S"
OUTPUT_DIR="/anvil/scratch/x-kguse/16S_analysis/fastqc_results"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run FastQC on all .fastq.gz files in the input directory
fastqc -t 4 -o "$OUTPUT_DIR" "${INPUT_DIR}"/*.fastq.gz

