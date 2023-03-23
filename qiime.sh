#variables
<load qiime2/2021.11>
BASE_PWD=path/to/folder
MANIFEST=path/to/manifest
METADATA=path/to/metadata
OUTPUT=path/to/output/folder
TRUNC_FORWARD=<truncation of forward reads, we used 280>
TRUNC_REVERSE= <truncation of reversereads, we used 261>
SAMPLE_DEPTH=<depth to which to sample, 416>

qiime tools import \--type 'SampleData[PairedEndSequencesWithQuality]' \--input-path $MANIFEST \--output-path $OUTPUT/paired-end-demux2.qza \--input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \--i-data $OUTPUT/paired-end-demux2.qza \--o-visualization $OUTPUT/paired-end-demux2.qzv

qiime dada2 denoise-paired \--i-demultiplexed-seqs $OUTPUT/paired-end-demux2.qza \--p-trunc-len-f $TRUNC_FORWARD \--p-trunc-len-r $TRUNC_REVERSE \--o-representative-sequences $OUTPUT/rep-seqs-dada2.qza \--o-table $OUTPUT/table-dada2.qza \--o-denoising-stats $OUTPUT/stats-dada2.qza

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
  
  
qiime diversity core-metrics-phylogenetic \--i-phylogeny $OUTPUT/rooted-tree.qza \--i-table $OUTPUT/filtered-table.qza \--p-sampling-depth $SAMPLE_DEPTH \--m-metadata-file $METADATA \--output-dir $BASE_PWD/core-metrics/ \--verbose


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
