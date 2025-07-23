#!/bin/bash
#SBATCH --job-name=ITS_tree
#SBATCH --output=ITS_tree_%j.out
#SBATCH --error=ITS_tree_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=standard

# Change to correct working directory
cd /anvil/projects/x-bio250086/soil_microbiome/ITS1/ITS_analysis

# Load QIIME2 (adjust based on your system module name)
module load qiime2

# Run cutadapt to trim primers
#qiime cutadapt trim-paired \
 # --i-demultiplexed-sequences ITS1_paired_end.qza \
  #--p-front-f CTTGGTCATTTAGAGGAAGTAA \
  #--p-front-r GCTGCGTTCTTCATCGATGC \
  #--p-error-rate 0.1 \
  #--o-trimmed-sequences ITS1_trimmed.qza \
  #--verbose


# Run DADA2 denoising
#qiime dada2 denoise-paired \
 # --i-demultiplexed-seqs ITS1_trimmed.qza \
  #--p-trunc-len-f 227 \
  #--p-trunc-len-r 180 \
  #--p-n-threads 8 \
  #--o-table ITS1_table.qza \
  #--o-representative-sequences ITS1_rep_seqs.qza \
  #--o-denoising-stats ITS1_denoising_stats.qza \
  #--verbose

###Make trees
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ITS1_rep_seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

##qiime diversity metrics
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table ITS1_table.qza \
  --p-sampling-depth 20000 \
  --m-metadata-file soil_metadata.tsv \
  --output-dir core-metrics-results

###Extract primers (515F/806R)
qiime feature-classifier extract-reads \
  --i-sequences ITS1_rep_seqs.qza \
  --p-f-primer CTTGGTCATTTAGAGGAAGTAA \
  --p-r-primer GCTGCGTTCTTCATCGATGC \
  --p-trunc-len 230 \
  --o-reads ref-seqs-trimmed-ITS1.qza


###Train the classifer
#qiime feature-classifier fit-classifier-naive-bayes \
 # --i-reference-reads ref-seqs-trimmed.qza \
  #--i-reference-taxonomy silva-138-99-tax.qza \
  #--o-classifier taxonomy.qza

