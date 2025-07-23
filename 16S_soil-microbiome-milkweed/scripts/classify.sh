#!/bin/bash
#SBATCH --job-name=classify_taxonomy
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=classify.out
#SBATCH --error=classify.err

cd anvil/scratch/x-kguse/16S_analysis
module load biocontainers
module load qiime2/2024.2

qiime feature-classifier classify-sklearn \
  --i-classifier taxonomy.qza \
  --i-reads rep-seqs.qza \
  --o-classification classified_taxonomy.qza \
  --p-n-jobs 8

