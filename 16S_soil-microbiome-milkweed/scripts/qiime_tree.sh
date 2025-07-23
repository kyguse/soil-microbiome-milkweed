#!/bin/bash
#SBATCH --job-name=qiime_tree
#SBATCH --output=qiime_tree_%j.out
#SBATCH --error=qiime_tree_%j.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

cd anvil/scratch/x-kguse/16S_analysis

module load biocontainers
module load qiime2


###Make trees
#qiime phylogeny align-to-tree-mafft-fasttree \
 # --i-sequences rep-seqs.qza \
  #--o-alignment aligned-rep-seqs.qza \
  #--o-masked-alignment masked-aligned-rep-seqs.qza \
 # --o-tree unrooted-tree.qza \
 # --o-rooted-tree rooted-tree.qza


##qiime diversity metrics
#qiime diversity core-metrics-phylogenetic \
 # --i-phylogeny rooted-tree.qza \
  #--i-table table.qza \
  #--p-sampling-depth 65000 \
  #--m-metadata-file soil_metadata.tsv \
  #--output-dir core-metrics-results

###Extract primers (515F/806R)
#qiime feature-classifier extract-reads \
 # --i-sequences silva-138-99-seqs.qza \
  #--p-f-primer GTGCCAGCMGCCGCGGTAA \
 # --p-r-primer GGACTACHVGGGTWTCTAAT \
 # --p-trunc-len 250 \
  #--o-reads ref-seqs-trimmed.qza

###Train the classifer
#qiime feature-classifier fit-classifier-naive-bayes \
 # --i-reference-reads ref-seqs-trimmed.qza \
 # --i-reference-taxonomy silva-138-99-tax.qza \
 # --o-classifier taxonomy.qza

##Assign taxonomy
qiime feature-classifier classify-sklearn \
  --i-classifier taxonomy.qza \
  --i-reads rep-seqs.qza \
  --o-classification classified_taxonomy.qza


###Visualize taxonomy
qiime metadata tabulate \
  --m-input-file classified_taxonomy.qza \
  --o-visualization classified_taxonomy.qzv

###Create qiime barplot
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy classified_taxonomy.qza \
  --m-metadata-file soil_metadata.tsv \
  --o-visualization taxa-bar-plots.qzv

###Create taxonomy table level 6 (genus)
qiime taxa collapse --i-table table.qza  --i-taxonomy classified_taxonomy.qza --p-level 6 --o-collapsed-table table6.qza











