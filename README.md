## Dietary vitamin B1, B2, and B6 intake influence the microbial composition and functional potential of the gut microbiome in Parkinson’s disease

**By Helena Sokolovska, Yixuan Zhang, Ayda Fathi, and Yoyo Lee**

This is the full code repository for our study. Processing of 16S rRNA sequencing data; feature table/phylogenetic tree generation; taxonomic classification; alpha/beta diversity analysis (excluding adonis/Dunn's testing); and functional potential analysis were done in bash. R was used for metadata filtering/stratifying nutrient intake data; adonis/Dunn's testing; and differential abundance analysis. 

### The following dependencies are required:
- Ubuntu 16.04.5
- Python 3.9.7 
- QIIME2 2021.4
- RStudio 2021.09.0
  - dplyr
  - readxl
  - xlsx
  - ggplot2
  - tidyverse
  - ape
  - phyloseq
  - DESeq2
  - vegan
  - FSA
- PICRUSt2 2.4.1
- STAMP

## For the full pipeline of analysis, run numbered scripts in order:
- Sokolovska_et_al_2022_script1.sh
- Sokolovska_et_al_2022_script2.R
