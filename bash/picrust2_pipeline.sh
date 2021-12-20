#samples with N/A values are manually removed from manifest.

#import demultiplexed sequences while in data directory
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /data/na_rm_parkinsons_manifest.txt \
  --output-path /data/demux_seqs.qza

#visualize metadata file
qiime metadata tabulate \
--m-input-file /data/na-rm/rm_+B3_stratified_na.rm_PD_metadata.tsv \
--o-visualization /data/na-rm/na_rm_metadata.qzv

#visualize the demultiplexed sequences 
qiime demux summarize \
  --i-data /data/demux_seqs.qza \
  --o-visualization /data/demux_seqs.qzv

#denoise with DADA2
qiime dada2 denoise-single \
  --i-demultiplexed-seqs /data/demux_seqs.qza \
  --p-trunc-len 251 \
  --o-representative-sequences /data/rep-seqs-dada2.qza \
  --o-table /data/table-dada2.qza \
  --o-denoising-stats /data/stats-dada2.qza

#visualize the denoising statistics 
qiime metadata tabulate \
--m-input-file /data/stats-dada2.qza  \
--o-visualization /data/stats-dada2.qzv

#create a feature table
qiime feature-table summarize \
  --i-table /data/table-dada2.qza \
  --m-sample-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --o-visualization /data/table-dada2.qzv

#filter feature table
#filter mitochondrial and chloroplastic ASVs from the feature table
qiime taxa filter-table \
  --i-table /data/table-dada2.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza

#visualize the filtered feature table
qiime feature-table summarize \
  --i-table /data/table-no-mitochondria-no-chloroplast.qza \
  --m-sample-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --o-visualization /data/table-no-mitochondria-no-chloroplast.qzv

#create unrooted and rooted trees using fasttree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /data/rep-seqs-dada2.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

#generate a rarefaction curve
qiime diversity alpha-rarefaction \
  --i-table /data/table-dada2.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --o-visualization /data/alpha_rarefaction.qzv \
  --p-max-depth 13000

#perform alpha and beta diversity analyses for each nutrient 
qiime diversity core-metrics-phylogenetic \
  --i-table /data/table-no-mitochondria-no-chloroplast.qza \
  --i-phylogeny /data/rooted-tree.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --p-sampling-depth 7000 \
  --output-dir /data/core-metrics-results_filtered_no_mitochondria

#visualize the statistical results of alpha diversity analysis
qiime diversity alpha-group-significance \
  --i-alpha-diversity /data/core-metrics-results_filtered_no_mitochondria/faith_pd_vector.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv\
  --o-visualization /data/core-metrics-results_filtered_no_mitochondria/faiths_pd_statistics.qzv

qiime diversity alpha-group-significance \
 --i-alpha-diversity /data/core-metrics-results_filtered_no_mitochondria/evenness_vector.qza \
 --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
 --o-visualization /data/core-metrics-results_filtered_no_mitochondria/evenness_statistics.qzv

qiime diversity alpha-group-significance \
 --i-alpha-diversity /data/core-metrics-results_filtered_no_mitochondria/shannon_vector.qza \
 --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
 --o-visualization /data/core-metrics-results_filtered_no_mitochondria/shannon_statistics.qzv

qiime diversity alpha-group-significance \
 --i-alpha-diversity /data/core-metrics-results_filtered_no_mitochondria/observed_features_vector.qza \
 --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
 --o-visualization /data/core-metrics-results_filtered_no_mitochondria/observed_features_statistics.qzv

#perform multivariant ANOVA for each nutrient
qiime longitudinal anova \
  --m-metadata-file core-metrics-results_filtered_no_mitochondria/faith_pd_vector.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --p-formula 'faith_pd ~ Disease * PUFA_intake_stratified' \
  --o-visualization core-metrics-results_filtered_no_mitochondria/faiths_pd_PUFA_anova.qzv

qiime longitudinal anova \
  --m-metadata-file core-metrics-results_filtered_no_mitochondria/faith_pd_vector.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --p-formula 'faith_pd ~ Disease * SFA_intake_stratified' \
  --o-visualization core-metrics-results_filtered_no_mitochondria/faiths_pd_SFA_anova.qzv

qiime longitudinal anova \
  --m-metadata-file core-metrics-results_filtered_no_mitochondria/faith_pd_vector.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --p-formula 'faith_pd ~ Disease * Total_Vitamin_A_intake_stratified' \
  --o-visualization core-metrics-results_filtered_no_mitochondria/faiths_pd_Total_Vitamin_A_anova.qzv

qiime longitudinal anova \
  --m-metadata-file core-metrics-results_filtered_no_mitochondria/faith_pd_vector.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --p-formula 'faith_pd ~ Disease * Vitamin_B1_intake_stratified' \
  --o-visualization core-metrics-results_filtered_no_mitochondria/faiths_pd_Vitamin_B1_anova.qzv

qiime longitudinal anova \
  --m-metadata-file core-metrics-results_filtered_no_mitochondria/faith_pd_vector.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --p-formula 'faith_pd ~ Disease * Vitamin_B2_intake_stratified' \
  --o-visualization core-metrics-results_filtered_no_mitochondria/faiths_pd_Vitamin_B2_anova.qzv

qiime longitudinal anova \
  --m-metadata-file core-metrics-results_filtered_no_mitochondria/faith_pd_vector.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --p-formula 'faith_pd ~ Disease * Vitamin_B3_intake_stratified' \
  --o-visualization core-metrics-results_filtered_no_mitochondria/faiths_pd_Vitamin_B3_anova.qzv

qiime longitudinal anova \
  --m-metadata-file core-metrics-results_filtered_no_mitochondria/faith_pd_vector.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --p-formula 'faith_pd ~ Disease * Vitamin_B6_intake_stratified' \
  --o-visualization core-metrics-results_filtered_no_mitochondria/faiths_pd_Vitamin_B6_anova.qzv

qiime longitudinal anova \
  --m-metadata-file core-metrics-results_filtered_no_mitochondria/faith_pd_vector.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --p-formula 'faith_pd ~ Disease * Vitamin_B12_intake_stratified' \
  --o-visualization core-metrics-results_filtered_no_mitochondria/faiths_pd_Vitamin_B12_anova.qzv

qiime longitudinal anova \
  --m-metadata-file core-metrics-results_filtered_no_mitochondria/faith_pd_vector.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --p-formula 'faith_pd ~ Disease * Vitamin_C_intake_stratified' \
  --o-visualization core-metrics-results_filtered_no_mitochondria/faiths_pd_Vitamin_C_anova.qzv

qiime longitudinal anova \
  --m-metadata-file core-metrics-results_filtered_no_mitochondria/faith_pd_vector.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --p-formula 'faith_pd ~ Disease * Vitamin_D_intake_stratified' \
  --o-visualization core-metrics-results_filtered_no_mitochondria/faiths_pd_Vitamin_D_anova.qzv

qiime longitudinal anova \
  --m-metadata-file core-metrics-results_filtered_no_mitochondria/faith_pd_vector.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --p-formula 'faith_pd ~ Disease * Vitamin_E_intake_stratified' \
  --o-visualization core-metrics-results_filtered_no_mitochondria/faiths_pd_Vitamin_E_anova.qzv

#for each nutrient, perform pairwise PERMANOVA using weighted & unweighted UniFrac
qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_A_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-Total_Vitamin_A_and_Disease-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_A_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/weighted-unifrac-Total_Vitamin_A_and_Disease-significance.qzv \
  --p-pairwise

#vitamin B1
qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_B1_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-Total_Vitamin_B1_and_Disease-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_B1_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/weighted-unifrac-Total_Vitamin_B1_and_Disease-significance.qzv \
  --p-pairwise

#vitamin B2
qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_B2_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-Total_Vitamin_B2_and_Disease-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_B2_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/weighted-unifrac-Total_Vitamin_B2_and_Disease-significance.qzv \
  --p-pairwise

#vitamin B3
qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_B3_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-Total_Vitamin_B3_and_Disease-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_B3_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/weighted-unifrac-Total_Vitamin_B3_and_Disease-significance.qzv \
  --p-pairwise

#vitamin B6
qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_B6_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-Total_Vitamin_B6_and_Disease-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_B6_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/weighted-unifrac-Total_Vitamin_B6_and_Disease-significance.qzv \
  --p-pairwise

#vitamin B12
qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_B12_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-Total_Vitamin_B12_and_Disease-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_B12_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/weighted-unifrac-Total_Vitamin_B12_and_Disease-significance.qzv \
  --p-pairwise

#vitamin C
qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_C_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-Total_Vitamin_C_and_Disease-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_C_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/weighted-unifrac-Total_Vitamin_C_and_Disease-significance.qzv \
  --p-pairwise

#vitamin D
qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_D_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-Total_Vitamin_D_and_Disease-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_C_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/weighted-unifrac-Total_Vitamin_D_and_Disease-significance.qzv \
  --p-pairwise

#vitamin E
qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_E_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-Total_Vitamin_E_and_Disease-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Total_Vitamin_C_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/weighted-unifrac-Total_Vitamin_E_and_Disease-significance.qzv \
  --p-pairwise

#SFA
qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column SFA_intake_stratified_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-SFA_and_Disease-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column SFA_intake_stratified_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/weighted-unifrac-SFA_and_Disease-significance.qzv \
  --p-pairwise

#PUFA
qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column PUFA_intake_stratified_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-PUFA_and_Disease-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column PUFA_intake_stratified_and_Disease \
  --o-visualization core-metrics-results_filtered_no_mitochondria/weighted-unifrac-PUFA_and_Disease-significance.qzv \
  --p-pairwise

#performed unweighted unifrac on disease status (SILVA)
qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Disease \
  --o-visualization /data/core-metrics-results_filtered_no_mitochondria/unweighted-unifrac-Disease-significance.qzv \
  --p-pairwise

#performed weighted unifrac on disease status (SILVA)
qiime diversity beta-group-significance \
  --i-distance-matrix /data/core-metrics-results_filtered_no_mitochondria/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --m-metadata-column Disease \
  --o-visualization /data/core-metrics-results_filtered_no_mitochondria/weighted-unifrac-Disease-significance.qzv \
  --p-pairwise

#export files for R

#export files to do analysis in R
qiime tools export \
  --input-path rooted-tree.qza \
  --output-path exported

qiime tools export \
  --input-path table-no-mitochondria-no-chloroplast.qza \
  --output-path exported

qiime tools export \
  --input-path taxonomy.qza \
  --output-path exported

#we’ll need to modify the exported taxonomy file’s header before using it with BIOM software in R. Before modifying that file, make a copy
cp exported/taxonomy.tsv biom-taxonomy.tsv

#change the column names: “Feature ID” to “#OTUID”, “Taxon to “taxonomy”, “Confidence” to “confidence”
nano exported/biom-taxonomy.tsv

#convert to biom file
biom convert -i table-no-mitochondria-no-chloroplast.qza -o feature-table-no-mitochondria-no-chloroplast..biom --table-type="OTU table" --to-hdf5

#create biom file with taxonomy information to use in R
biom add-metadata \
  -i exported/feature-table.biom \
  -o exported/table-with-taxonomy.biom \
  --observation-metadata-fp exported/biom-taxonomy.tsv\
  --sc-separated taxonomy

#export biom files to local computer
scp root@10.19.139.102:/data/exported/*.biom .
scp root@10.19.139.102:/data/exported/*.nwk .

#taxonomic classification 
#Filter metadata columns to include only high vitamin B6 intake
qiime feature-table filter-samples \
  --i-table  table-no-mitochondria-no-chloroplast.qza \
  --m-metadata-file colrm_stratified_PD_Disease_Combined_metadata.tsv \
  --p-where "[Total_Vitamin_B6_and_Disease] IN ('high and PD', 'high and Control')" \
  --o-filtered-table high_B6-filtered-table.qza

 #taxonomy classification using machine-learning classifier SILVA
 qiime feature-classifier classify-sklearn \
  --i-reads ./rep-seqs-filtered.qza\
  --i-classifier /mnt/datasets/classifiers/silva-138-99-515-806-nb-classifier.qza \
  --o-classification ./silva-138-99-515-806-nb-taxonomy.qza
  
  #visualize taxonomy associated with the sequences
  qiime metadata tabulate \
  --m-input-file ./silva-138-99-515-806-nb-taxonomy.qza \
  --o-visualization ./silva-138-99-515-806-nb-taxonomy.qzv
  
  #tabulate the representative sequences
  qiime feature-table tabulate-seqs \
  --i-data ./rep-seqs-filtered.qza \
  --o-visualization ./rep-seqs-filtered.qzv

 #taxonomy barchart
 #first filter out any samples with fewer features than our 
 #rarefaction threshold (7000). 

  qiime feature-table filter-samples \
  --i-table ./high_B6-filtered-table.qza\
  --p-min-frequency 7000 \
  --o-filtered-table ./table_7k_filtered.qza

  #use the filtered table to build an interactive barplot of the taxonomy in each sample.
  qiime taxa barplot \
  --i-table ./table_7k_filtered.qza \
  --i-taxonomy ./silva-138-99-515-806-nb-taxonomy.qza \
  --m-metadata-file ./colrm_stratified_PD_Disease_Combined_metadata.tsv\
  --o-visualization ./taxa_barplot.qzv
  
scp root@10.19.139.56:/data/project_2/taxa_barplot.qzv .

#PICRUSt Analysis
#!/usr/bin/env bash
set -ex

# [1] scp study_seqs.fna and study_seqs.biom into the picrust2_out_pipeline directory

# [2] pipe sequence placement, hidden-state prediction of genomes, metagenome prediction and pathway-level predictions
picrust2_pipeline.py -s /data/study_seqs.fna -i /data/study_seqs.biom -o /data/picrust2_out_pipeline -p 1

# [3] add description of each functional category to table of genes or pathways

add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
-o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
-o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
-o pathways_out/path_abun_unstrat_descrip.tsv.gz

# [4] generate shuffled predictions
shuffle_predictions.py -i EC_predicted.tsv.gz \
-o EC_predicted_shuffled \
-r 5 \
-s 131

# [5] Make folders for shuffled output
mkdir EC_metagenome_out_shuffled
mkdir pathways_out_shuffled

for i in {1..5}; do
    # Define in and out file paths.
    EC_SHUFFLED="EC_predicted_shuffled/EC_predicted_shuf"$i".tsv.gz"
    OUT_META="EC_metagenome_out_shuffled/rep"$i
    OUT_PATHWAYS="pathways_out_shuffled/rep"$i
    # [6] PICRUSt2 scripts to get prediction abundance tables for gene and pathway levels, respectively.
    metagenome_pipeline.py -i ../study_seqs.biom -m marker_predicted_and_nsti.tsv.gz -f $EC_SHUFFLED \
    -o $OUT_META \
    --strat_out
    pathway_pipeline.py -i $OUT_META/pred_metagenome_contrib.tsv.gz \
    -o $OUT_PATHWAYS \
    -p 1
done

