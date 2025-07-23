#!/bin/bash
#SBATCH --job-name=unite_classifier
#SBATCH --output=unite_classifier_%j.out
#SBATCH --error=unite_classifier_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=standard

#change directory
cd /anvil/projects/x-bio250086/soil_microbiome/ITS1/ITS_analysis

# Load QIIME2 module (adjust version if needed)
module load biocontainers
module load qiime2

#Train the Classifier
qiime feature-classifier extract-reads \
  --i-sequences unite-refs-99-10.0.qza \
  --p-f-primer CTTGGTCATTTAGAGGAAGTAA \
  --p-r-primer GCTGCGTTCTTCATCGATGC \
  --p-trunc-len 230 \
  --o-reads unite-refs-ITS1.qza

#Train the Classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads unite-refs-ITS1.qza \
  --i-reference-taxonomy unite-tax-99-10.0.qza \
  --o-classifier unite-ITS1-classifier.qza

##Assign Taxonomy to representative sequences
qiime feature-classifier classify-sklearn \
  --i-classifier unite-ITS1-classifier.qza \
  --i-reads ITS1_rep_seqs.qza \
  --o-classification classified_taxonomy.qza

###Visualize taxonomy
qiime metadata tabulate \
  --m-input-file classified_taxonomy.qza \
  --o-visualization classified_taxonomy.qzv

###Create qiime barplot
qiime taxa barplot \
  --i-table ITS1_table.qza \
  --i-taxonomy classified_taxonomy.qza \
  --m-metadata-file soil_metadata.tsv \
  --o-visualization taxa-bar-plots.qzv

###Create taxonomy table level 6 (genus)
qiime taxa collapse --i-table ITS1_table.qza  --i-taxonomy classified_taxonomy.qza --p-level 6 --o-collapsed-table ITS1_table6.qza
