#!/bin/bash

#########################################################################################################
# Sokolovska_et_al_2022_script1.sh
# Paper: "Dietary vitamin B1, B2, and B6 intake influence the microbial composition and functional potential of the gut microbiome in Parkinson’s disease"
# Authors: Helena Sokolovska, Yixuan Zhang, Ayda Fathi, and Yoyo Lee
# Date: Sep 5, 2022
# Purpose:
  # QIIME2 analysis - generating feature table, alpha/beta diversity analysis, taxonomic classification
  # PICRUSt2 analysis - differential functional potential analysis
#########################################################################################################

mkdir ./data
mkdir ./data/exported

######################################## QIIME2 analysis: Generating Feature Table ########################################

# samples with N/A values are manually removed from manifest

# import demultiplexed sequences from server while in data directory
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path ./data/na_rm_parkinsons_manifest.txt \
  --output-path ./data/demux_seqs.qza

# visualize metadata file
qiime metadata tabulate \
--m-input-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
--o-visualization ./data/colrm_stratified_PD_Disease_Combined_metadata.qzv

# visualize the demultiplexed sequences 
qiime demux summarize \
  --i-data ./data/demux_seqs.qza \
  --o-visualization ./data/demux_seqs.qzv

# denoise with DADA2
qiime dada2 denoise-single \
  --i-demultiplexed-seqs ./data/demux_seqs.qza \
  --p-trunc-len 251 \
  --o-representative-sequences ./data/rep-seqs-dada2.qza \
  --o-table ./data/table-dada2.qza \
  --o-denoising-stats ./data/stats-dada2.qza

# visualize the denoising statistics 
qiime metadata tabulate \
--m-input-file ./data/stats-dada2.qza  \
--o-visualization ./data/stats-dada2.qzv

# visualize feature table
qiime feature-table summarize \
  --i-table ./data/table-dada2.qza \
  --m-sample-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --o-visualization ./data/table-dada2.qzv

# download SILVA classifier
wget \
  -O "gg-13-8-99-515-806-nb-classifier.qza" \
  "https://data.qiime2.org/2022.2/common/silva-138-99-515-806-nb-classifier.qza"

# taxonomy classification using SILVA classifier
qiime feature-classifier classify-sklearn \
  --i-reads ./data/rep-seqs-dada2.qza \
  --i-classifier ./silva-138-99-515-806-nb-classifier.qza \
  --o-classification ./data/table-with-taxonomy.qza

# filter feature table and sequences:

# filter mitochondrial and chloroplastic ASVs from the feature table
qiime taxa filter-table \
  --i-table ./data/table-dada2.qza \
  --i-taxonomy ./data/taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table ./data/feature-table-no-mitochondria-no-chloroplast.qza

# visualize the filtered feature table
qiime feature-table summarize \
  --i-table ./data/feature-table-no-mitochondria-no-chloroplast.qza \
  --m-sample-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --o-visualization ./data/feature-table-no-mitochondria-no-chloroplast.qzv

# filter mitochondrial and chloroplastic ASVs from the representative sequences
qiime taxa filter-seqs \
  --i-sequences ./data/rep-seqs-dada2.qza \
  --i-taxonomy ./data/taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences ./data/sequences-no-mitochondria-no-chloroplast.qza

######################################## QIIME2 analysis: Alpha/Beta Diversity Analysis ########################################

# create unrooted and rooted trees using fasttree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ./data/sequences-no-mitochondria-no-chloroplast.qza \
  --o-alignment ./data/aligned-rep-seqs.qza \
  --o-masked-alignment ./data/masked-aligned-rep-seqs.qza \
  --o-tree ./data/unrooted-tree.qza \
  --o-rooted-tree ./data/rooted-tree.qza

# generate a rarefaction curve
qiime diversity alpha-rarefaction \
  --i-table ./data/table-dada2.qza \
  --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --o-visualization ./data/alpha_rarefaction.qzv \
  --p-max-depth 13000

# perform alpha and beta diversity analyses for each nutrient 
qiime diversity core-metrics-phylogenetic \
  --i-table ./data/feature-table-no-mitochondria-no-chloroplast.qzv \
  --i-phylogeny ./data/rooted-tree.qza \
  --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --p-sampling-depth 7000 \
  --output-dir ./data/core-metrics-results_filtered_no_mitochondria

# visualize the statistical results of alpha diversity analysis (contains Kruskal-Wallis statistics)

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./data/core-metrics-results_filtered_no_mitochondria/faith_pd_vector.qza \
  --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv\
  --o-visualization ./data/core-metrics-results_filtered_no_mitochondria/faiths_pd_statistics.qzv

qiime diversity alpha-group-significance \
 --i-alpha-diversity ./data/core-metrics-results_filtered_no_mitochondria/evenness_vector.qza \
 --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
 --o-visualization ./data/core-metrics-results_filtered_no_mitochondria/evenness_statistics.qzv

qiime diversity alpha-group-significance \
 --i-alpha-diversity ./data/core-metrics-results_filtered_no_mitochondria/shannon_vector.qza \
 --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
 --o-visualization ./data/core-metrics-results_filtered_no_mitochondria/shannon_statistics.qzv

qiime diversity alpha-group-significance \
 --i-alpha-diversity ./data/core-metrics-results_filtered_no_mitochondria/observed_features_vector.qza \
 --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
 --o-visualization ./data/core-metrics-results_filtered_no_mitochondria/observed_features_statistics.qzv

# beta diversity analysis: adonis test and pairwise PERMANOVA

# adonis testing to identify variables which should be analyzed in pairwise PERMANOVA
# weighted UniFrac: PUFAs, vitamins B2/B12/D, disease status
# unweighted UniFrac: PUFAs, vitamins A/B1/B2/B3/B6/C/E, disease status

# pairwise PERMANOVA

# vitamin A (unweighted UniFrac)
qiime diversity beta-group-significance \
  --i-distance-matrix ./data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_A_and_Disease \
  --o-visualization ./data/core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-Total_Vitamin_A_and_Disease-significance.qzv \
  --p-pairwise

# vitamin B1 (unweighted UniFrac)
qiime diversity beta-group-significance \
  --i-distance-matrix ./data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_B1_and_Disease \
  --o-visualization ./data/core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-Total_Vitamin_B1_and_Disease-significance.qzv \
  --p-pairwise

# vitamin B2 (weighted and unweighted UniFrac)
qiime diversity beta-group-significance \
  --i-distance-matrix ./data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_B2_and_Disease \
  --o-visualization ./data/core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-Total_Vitamin_B2_and_Disease-significance.qzv \
  --p-pairwise
qiime diversity beta-group-significance \
  --i-distance-matrix ./data/core-metrics-results_filtered_no_mitochondria/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_B2_and_Disease \
  --o-visualization ./data/core-metrics-results_filtered_no_mitochondria/weighted-unifrac-Total_Vitamin_B2_and_Disease-significance.qzv \
  --p-pairwise

# vitamin B3 (unweighted UniFrac)
qiime diversity beta-group-significance \
  --i-distance-matrix ./data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_B3_and_Disease \
  --o-visualization ./data/core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-Total_Vitamin_B3_and_Disease-significance.qzv \
  --p-pairwise

# vitamin B6 (unweighted UniFrac)
qiime diversity beta-group-significance \
  --i-distance-matrix ./data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_B6_and_Disease \
  --o-visualization ./data/core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-Total_Vitamin_B6_and_Disease-significance.qzv \
  --p-pairwise

# vitamin B12 (weighted UniFrac)
qiime diversity beta-group-significance \
  --i-distance-matrix ./data/core-metrics-results_filtered_no_mitochondria/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_B12_and_Disease \
  --o-visualization ./data/core-metrics-results_filtered_no_mitochondria/weighted-unifrac-Total_Vitamin_B12_and_Disease-significance.qzv \
  --p-pairwise

# vitamin C (unweighted UniFrac)
qiime diversity beta-group-significance \
  --i-distance-matrix ./data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_C_and_Disease \
  --o-visualization ./data/core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-Total_Vitamin_C_and_Disease-significance.qzv \
  --p-pairwise

# vitamin D (weighted UniFrac)
qiime diversity beta-group-significance \
  --i-distance-matrix ./data/core-metrics-results_filtered_no_mitochondria/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_C_and_Disease \
  --o-visualization ./data/core-metrics-results_filtered_no_mitochondria/weighted-unifrac-Total_Vitamin_D_and_Disease-significance.qzv \
  --p-pairwise

# vitamin E (unweighted UniFrac)
qiime diversity beta-group-significance \
  --i-distance-matrix ./data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_E_and_Disease \
  --o-visualization ./data/core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-Total_Vitamin_E_and_Disease-significance.qzv \
  --p-pairwise

# SFA (none)

# PUFA (weighted and unweighted UniFrac)
qiime diversity beta-group-significance \
  --i-distance-matrix ./data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column PUFA_intake_stratified_and_Disease \
  --o-visualization ./data/core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-PUFA_and_Disease-significance.qzv \
  --p-pairwise
qiime diversity beta-group-significance \
  --i-distance-matrix ./data/core-metrics-results_filtered_no_mitochondria/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column PUFA_intake_stratified_and_Disease \
  --o-visualization ./data/core-metrics-results_filtered_no_mitochondria/weighted-unifrac-PUFA_and_Disease-significance.qzv \
  --p-pairwise

# disease status (weighted and unweighted UniFrac)
qiime diversity beta-group-significance \
  --i-distance-matrix ./data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Disease \
  --o-visualization ./data/core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-Disease-significance.qzv \
  --p-pairwise
qiime diversity beta-group-significance \
  --i-distance-matrix ./data/core-metrics-results_filtered_no_mitochondria/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Disease \
  --o-visualization ./data/core-metrics-results_filtered_no_mitochondria/weighted-unifrac-Disease-significance.qzv \
  --p-pairwise

# export files to do analysis in R and PICRUSt2
qiime tools export \
  --input-path ./data/rooted-tree.qza \
  --output-path ./data/exported
qiime tools export \
  --input-path ./data/feature-table-no-mitochondria-no-chloroplast.qzv \
  --output-path ./data
qiime tools export \
  --input-path ./data/taxonomy.qza \
  --output-path ./data
qiime tools export \
  --input-path ./data/sequences-no-mitochondria-no-chloroplast.qza \
  --output-path ./data
qiime tools export \
  --input-path ./data/core-metrics-results_filtered_no_mitochondria/weighted_unifrac_distance_matrix.qza \
  --output-path ./data/exported
qiime tools export \
  --input-path ./data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --output-path ./data/exported

# we’ll need to modify the exported taxonomy file’s header before using it with BIOM software in R. Before modifying that file, make a copy
cp ./data/taxonomy.tsv ./data/biom-taxonomy.tsv

# change the column names: “Feature ID” to “#OTUID”, “Taxon to “taxonomy”, “Confidence” to “confidence”
nano ./data/biom-taxonomy.tsv

# convert feature table to biom file to use in R and PICRUSt2
biom convert \
  -i ./data/feature-table-no-mitochondria-no-chloroplast.qza \
  -o ./data/feature-table-no-mitochondria-no-chloroplast.biom \
  --table-type="OTU table" \
  --to-hdf5

# convert representative sequences to fna file to use in PICRUSt2
mv ./data/sequences-no-mitochondria-no-chloroplast.fasta ./data/sequences-no-mitochondria-no-chloroplast.fna

# create biom file of feature table with taxonomy information to use in R
biom add-metadata \
  -i ./data/feature-table-no-mitochondria-no-chloroplast.biom \
  -o ./data/exported/table-with-taxonomy.biom \
  --observation-metadata-fp ./data/biom-taxonomy.tsv\
  --sc-separated taxonomy

# export biom files to local computer
scp root@10.19.139.102:./data/exported/*.biom .
scp root@10.19.139.102:./data/exported/*.nwk .
scp root@10.19.139.102:./data/exported/*.txt .

######################################## QIIME2 analysis: Taxonomic Classification ########################################
  
# visualize taxonomy associated with the sequences
qiime metadata tabulate \
  --m-input-file ./data/taxonomy.qza \
  --o-visualization ./data/taxonomy.qzv
  
# visualize the representative sequences
qiime feature-table tabulate-seqs \
  --i-data ./data/sequences-no-mitochondria-no-chloroplast.qza \
  --o-visualization ./data/sequences-no-mitochondria-no-chloroplast.qzv

# taxonomy barchart:

# filter out any samples with fewer features than our 
# rarefaction threshold (7000). 
qiime feature-table filter-samples \
  --i-table ./data/feature-table-no-mitochondria-no-chloroplast.qza \
  --p-min-frequency 7000 \
  --o-filtered-table ./data/table_7k_filtered.qza

# use the filtered table to build an interactive barplot of the taxonomy in each sample
qiime taxa barplot \
  --i-table ./data/table_7k_filtered.qza \
  --i-taxonomy ./data/taxonomy.qza \
  --m-metadata-file ./data/colrm_stratified_PD_Disease_Combined_metadata.tsv\
  --o-visualization ./data/taxa_barplot.qzv

######################################## PICRUSt2 analysis: Functional Potential Analysis ########################################

set -ex

# [1] pipe sequence placement, hidden-state prediction of genomes, metagenome prediction and pathway-level predictions
# input representative sequences and feature table as arguments
picrust2_pipeline.py -s ./data/sequences-no-mitochondria-no-chloroplast.fna -i ./data/feature-table-no-mitochondria-no-chloroplast.biom -o ./data/picrust2_out_pipeline -p 1

# [2] generate shuffled predictions
shuffle_predictions.py -i EC_predicted.tsv.gz \
-o EC_predicted_shuffled \
-r 5 \
-s 131

# make folders for shuffled output
mkdir EC_metagenome_out_shuffled
mkdir pathways_out_shuffled

for i in {1..5}; do
    # Define in and out file paths.
    EC_SHUFFLED="EC_predicted_shuffled/EC_predicted_shuf"$i".tsv.gz"
    OUT_META="EC_metagenome_out_shuffled/rep"$i
    OUT_PATHWAYS="pathways_out_shuffled/rep"$i
    # PICRUSt2 scripts to get prediction abundance tables for gene and pathway levels, respectively.
    metagenome_pipeline.py -i ./data/feature-table-no-mitochondria-no-chloroplast.biom -m marker_predicted_and_nsti.tsv.gz -f $EC_SHUFFLED \
    -o $OUT_META \
    --strat_out
    pathway_pipeline.py -i $OUT_META/pred_metagenome_contrib.tsv.gz \
    -o $OUT_PATHWAYS \
    -p 1
done
