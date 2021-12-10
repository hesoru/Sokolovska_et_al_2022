#!/bin/bash
set -ex

# scp study_seqs.fna and study_seqs.biom into the picrust2_out_pipeline directory

# [2] pipe sequence placement, hidden-state prediction of genomes, metagenome prediction and pathway-level predictions
# picrust2_pipeline.py -s /data/study_seqs.fna -i /data/study_seqs.biom -o /data/picrust2_out_pipeline -p 1

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
