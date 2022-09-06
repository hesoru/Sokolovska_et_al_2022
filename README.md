## Dietary vitamin B1, B2, and B6 intake influence the microbial composition and functional potential of the gut microbiome in Parkinson’s disease

**By Helena Sokolovska, Yixuan Zhang, Ayda Fathi, and Yoyo Lee**

This is the full code repository for our study. Processing of 16S rRNA sequencing data; feature table/phylogenetic tree generation; taxonomic classification; alpha/beta diversity analysis (excluding adonis/Dunn's testing); and functional potential analysis were done in bash. R was used for metadata filtering/stratifying nutrient intake data; adonis/Dunn's testing; and differential abundance analysis. 

### The following dependencies are required:
- Ubuntu 16.04.5
- Python 3.9.7 
- QIIME2 2021.4
- RStudio 2021.09.0
  - readxl 1.3.1  
  - xlsx 0.6.5 
  - dplyr 1.0.8
  - tidyverse 1.3.1 
  - ape 5.6-2  
  - phyloseq 1.34.0 
  - DESeq2 1.30.1 
  - vegan 2.5-7
  - FSA 0.9.3
  - ggplot2 3.3.5 
- PICRUSt2 2.4.1
- STAMP 2.1.3

### For the full pipeline of analysis, run numbered scripts in order:
- Sokolovska_et_al_2022_script1.sh
- Sokolovska_et_al_2022_script2.R
