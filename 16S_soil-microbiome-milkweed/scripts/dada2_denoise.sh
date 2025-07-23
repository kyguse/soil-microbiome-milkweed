#!/bin/bash
#SBATCH --job-name=dada2_denoise
#SBATCH --output=dada2_denoise_%j.out
#SBATCH --error=dada2_denoise_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# Load QIIME 2 module
module load biocontainers
module load qiime2/2024.2

# change directory
cd /anvil/scratch/x-kguse/16S_analysis 

# Run DADA2 denoise-paired
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-paired-end.qza \
  --p-trim-left-f 19 \
  --p-trim-left-r 20 \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 200 \
  --p-n-threads 8 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

