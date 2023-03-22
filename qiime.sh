#!/bin/sh
#SBATCH -J Sush
#SBATCH --time=00-20:00:00
#SBATCH -p batch 
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -c 1 
#SBATCH --mem=16g
#SBATCH --output=MyJob.%j.%N.out
#SBATCH --error=MyJob.%j.%N.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sushrut.jangi@tufts.edu

#variables
module load qiime2/2021.11
source activate qiime2-2021.11
BASE_PWD=/cluster/tufts/jangilab/shared/its2
MANIFEST=$BASE_PWD/manifest3.txt
METADATA=$BASE_PWD/Fungal_Metadata_Final830.tsv
OUTPUT=$BASE_PWD
TRUNC_FORWARD=280
TRUNC_REVERSE=261
SAMPLE_DEPTH=416
CORE=core-metric-results-5

qiime tools import \--type 'SampleData[PairedEndSequencesWithQuality]' \--input-path $MANIFEST \--output-path $OUTPUT/paired-end-demux2.qza \--input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \--i-data $OUTPUT/paired-end-demux2.qza \--o-visualization $OUTPUT/paired-end-demux2.qzv

qiime dada2 denoise-paired \--i-demultiplexed-seqs $OUTPUT/paired-end-demux2.qza \--p-trunc-len-f 280 \--p-trunc-len-r 261 \--o-representative-sequences $OUTPUT/rep-seqs-dada2.qza \--o-table $OUTPUT/table-dada2.qza \--o-denoising-stats $OUTPUT/stats-dada2.qza

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-dada2.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
  
qiime feature-table filter-samples \
  --i-table table-dada2.qza \
  --m-metadata-file $METADATA \
  --o-filtered-table filtered-table.qza \
  
  qiime feature-table summarize \
  --i-table filtered-table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file $METADATA \
  
  
qiime diversity core-metrics-phylogenetic \--i-phylogeny $OUTPUT/rooted-tree.qza \--i-table $OUTPUT/filtered-table.qza \--p-sampling-depth $SAMPLE_DEPTH \--m-metadata-file $METADATA \--output-dir $BASE_PWD/$CORE \--verbose


# Fungal Classifier training 

qiime tools import \
 --type FeatureData[Sequence] \
 --input-path $INPUT/sh_qiime_release_10.05.2021/sh_refs_qiime_ver8_dynamic_10.05.2021.fasta \
 --output-path $OUTPUT/qiime_ver8_dynamic_10.05.2021.qza 
 
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path $INPUT/sh_qiime_release_10.05.2021/sh_taxonomy_qiime_ver8_dynamic_10.05.2021.txt \
--output-path $OUTPUT/ref-taxonomy.qza

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads $INPUT/qiime_ver8_dynamic_10.05.2021.qza  \
--i-reference-taxonomy $INPUT/ref-taxonomy.qza \
--o-classifier $OUTPUT/classifier.qza

# Creation of ASV Table

qiime feature-classifier classify-sklearn \
  --i-classifier $INPUT/classifier.qza \
  --i-reads $OUTPUT/rep-seqs-dada2.qza \
  --o-classification $OUTPUT/taxonomy.qza

qiime metadata tabulate \
  --m-input-file $OUTPUT/taxonomy.qza \
  --o-visualization $OUTPUT/taxonomy.qzv

qiime taxa barplot \
  --i-table $OUTPUT/filtered-table.qza \
  --i-taxonomy $OUTPUT/taxonomy.qza \
  --m-metadata-file $METADATA \
  --o-visualization $OUTPUT/taxa-bar-plots.qzv


qiime tools export \
--input-path $INPUT/filtered-table.qza \
--output-path $OUTPUT/filtered-table

biom convert \
  -i $OUTPUT/filtered-table/feature-table.biom \
  -o $OUTPUT/otu_table.txt \
  --to-tsv
  
qiime tools export \
--input-path $INPUT/taxonomy.qza \
--output-path $OUTPUT/taxonomy
  
qiime tools export \
--input-path $INPUT/unrooted-tree.qza \
--output-path $OUTPUT/unrooted-tree


source deactivate
echo End: $(date) >> $FOLDER/run.txtsource deactivate
