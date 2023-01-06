#########################################################################################################
# Sokolovska_et_al_2022_script2.R
# Paper: "Dietary vitamin B1, B2, and B6 intake influence the microbial composition and functional potential of the gut microbiome in Parkinson’s disease"
# Authors: Helena Sokolovska, Yixuan Zhang, Ayda Fathi, and Yoyo Lee
# Date: Sep 5, 2022
# Purpose: R analysis - nutrient stratification, differential/relative abundance analysis
#########################################################################################################



######################################## Nutrient Stratification ########################################

# load packages
library(ggplot2)
library(dplyr)
library(readxl)
library(xlsx)
library(vegan)
library(tidyverse)

# select a specific set of random numbers and makes the analysis reproducible
set.seed(711)

# import metadata Excel file
parkinsons_metadata <- read_excel("~/4th Year Uni/TERM 2/MICB 447/Datasets/parkinsons_metadata.xlsx")
View(parkinsons_metadata)

# create data frame and clean up data:

# create data frame of metadata
PD_data_frame <- data.frame(parkinsons_metadata)

# change first column name to remove hashtag (interferes with coding)
PD_data_frame <- parkinsons_metadata %>%
  rename("SampleID" = "#SampleID")

# add new column to data frame that sums up vitamin A + equivalents
PD_data_frame$Total_Vitamin_A <- (PD_data_frame$Vitamin_A + PD_data_frame$Vitamin_A_equiv)

# remove N/A values in nutrition data prior to summary stats
na.rm_PD_data_frame <- PD_data_frame %>%
  filter(!is.na(PUFA), !is.na(SFA), !is.na(Total_Vitamin_A), !is.na(Vitamin_B1), !is.na(Vitamin_B2), !is.na(Niacin), !is.na(Vitamin_B6), !is.na(Vitamin_B12), !is.na(Vitamin_C), !is.na(Vitamin_D), !is.na(Vitamin_E))

# determining quartiles for stratification:

# split metadata into 2 data frames: control-only and PD-only
na.rm_control_only_data_frame <- filter(na.rm_PD_data_frame, Disease == "Control")
na.rm_PD_only_data_frame <- filter(na.rm_PD_data_frame, Disease == "PD")

# summary stats for each vitamin intake distribution, using control-only data frame
# (to determine quartiles for stratification)

summary(na.rm_control_only_data_frame$PUFA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.90   10.40   12.90   14.22   17.46   36.20 
summary(na.rm_control_only_data_frame$SFA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.594  17.936  22.600  25.061  29.900  61.900 
summary(na.rm_control_only_data_frame$Total_Vitamin_A)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 374.3  1036.3  1380.3  1622.9  1710.5 12601.6 
summary(na.rm_control_only_data_frame$Vitamin_B1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4000  0.9726  1.2000  1.3078  1.4262  6.8563 
summary(na.rm_control_only_data_frame$Vitamin_B2)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.600   1.200   1.600   1.634   1.990   4.500 
summary(na.rm_control_only_data_frame$Niacin)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8.959  15.183  19.901  20.457  23.400  50.400
summary(na.rm_control_only_data_frame$Vitamin_B6)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9535  1.5805  1.9193  2.0082  2.3000  4.8000 
summary(na.rm_control_only_data_frame$Vitamin_B12)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.02695  2.98709  4.80000  5.44874  7.00000 27.50000 
summary(na.rm_control_only_data_frame$Vitamin_C)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 18.55   80.50  112.80  127.71  157.99  334.02
summary(na.rm_control_only_data_frame$Vitamin_D)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   1.300   2.025   2.387   2.824   8.815 
summary(na.rm_control_only_data_frame$Vitamin_E)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.771   8.800  11.100  12.172  14.878  31.700 

# make boxplots of nutrient intake distributions, confirm intake Q1/Q3 for control match summary stats:

na.rm_PD_data_frame %>%
  ggplot(aes(x = Disease, y = PUFA)) +
  geom_boxplot() +
  labs(title = "PUFA Intake",
       x = "Disease status",
       y = "Intake (g)") +
  theme_bw()

na.rm_PD_data_frame %>%
  ggplot(aes(x = Disease, y = SFA)) +
  geom_boxplot() +
  labs(title = "SFA Intake",
       x = "Disease status",
       y = "Intake (g)") +
  theme_bw()

na.rm_PD_data_frame %>%
  ggplot(aes(x = Disease, y = Total_Vitamin_A)) +
  geom_boxplot() +
  labs(title = "(Vitamin A + Equivalents) Intake",
       x = "Disease status",
       y = "Intake (mcg)") +
  theme_bw()

na.rm_PD_data_frame %>%
  ggplot(aes(x = Disease, y = Vitamin_B1)) +
  geom_boxplot() +
  labs(title = "Vitamin B1 Intake",
       x = "Disease status",
       y = "Intake (mg)") +
  theme_bw()

na.rm_PD_data_frame %>%
  ggplot(aes(x = Disease, y = Vitamin_B2)) +
  geom_boxplot() +
  labs(title = "Vitamin B2 Intake",
       x = "Disease status",
       y = "Intake (mg)") +
  theme_bw()

na.rm_PD_data_frame %>%
  ggplot(aes(x = Disease, y = Niacin)) +
  geom_boxplot() +
  labs(title = "Vitamin B3 Intake",
       x = "Disease status",
       y = "Intake (mg)") +
  theme_bw()

na.rm_PD_data_frame %>%
  ggplot(aes(x = Disease, y = Vitamin_B6)) +
  geom_boxplot() +
  labs(title = "Vitamin B6 Intake",
       x = "Disease status",
       y = "Intake (mg)") +
  theme_bw()

na.rm_PD_data_frame %>%
  ggplot(aes(x = Disease, y = Vitamin_B12)) +
  geom_boxplot() +
  labs(title = "Vitamin B12 Intake",
       x = "Disease status",
       y = "Intake (mcg)") +
  theme_bw()

na.rm_PD_data_frame %>%
  ggplot(aes(x = Disease, y = Vitamin_C)) +
  geom_boxplot() +
  labs(title = "Vitamin C Intake",
       x = "Disease status",
       y = "Intake (mg)") +
  theme_bw()

na.rm_PD_data_frame %>%
  ggplot(aes(x = Disease, y = Vitamin_D)) +
  geom_boxplot() +
  labs(title = "Vitamin D Intake",
       x = "Disease status",
       y = "Intake (mcg)") +
  theme_bw()

na.rm_PD_data_frame %>%
  ggplot(aes(x = Disease, y = Vitamin_E)) +
  geom_boxplot() +
  labs(title = "Vitamin E Intake",
       x = "Disease status",
       y = "Intake (mg)") +
  theme_bw()

# stratify nutrient intake data based on intake Q1/Q3 for control:

stratified_PD_data_frame <- mutate(PD_data_frame, PUFA_intake_stratified = cut(PUFA,
                                                                               breaks = c(0, 10.40, 17.46, Inf),
                                                                               labels = c("low",
                                                                                          "moderate",
                                                                                          "high")))

stratified_PD_data_frame <- mutate(stratified_PD_data_frame, SFA_intake_stratified = cut(SFA,
                                                                                           breaks = c(0, 17.936, 29.900, Inf),
                                                                                           labels = c("low",
                                                                                                      "moderate",
                                                                                                      "high")))

stratified_PD_data_frame <- mutate(stratified_PD_data_frame, Total_Vitamin_A_intake_stratified = cut(Total_Vitamin_A,
                                                                                                     breaks = c(0, 1036.3, 1710.5, Inf),
                                                                                                     labels = c("low",
                                                                                                                "moderate",
                                                                                                                "high")))

stratified_PD_data_frame <- mutate(stratified_PD_data_frame, Vitamin_B1_intake_stratified = cut(Vitamin_B1,
                                                                                                                 breaks = c(0, 0.9726, 1.4262, Inf),
                                                                                                                 labels = c("low",
                                                                                                                            "moderate",
                                                                                                                            "high")))

stratified_PD_data_frame <- mutate(stratified_PD_data_frame, Vitamin_B2_intake_stratified = cut(Vitamin_B2,
                                                                                                            breaks = c(0, 1.200, 1.990, Inf),
                                                                                                            labels = c("low",
                                                                                                                       "moderate",
                                                                                                                       "high")))

stratified_PD_data_frame <- mutate(stratified_PD_data_frame, Vitamin_B3_intake_stratified = cut(Niacin,
                                                                                                            breaks = c(0, 15.183, 23.400, Inf),
                                                                                                            labels = c("low",
                                                                                                                       "moderate",
                                                                                                                       "high")))

stratified_PD_data_frame <- mutate(stratified_PD_data_frame, Vitamin_B6_intake_stratified = cut(Vitamin_B6,
                                                                                                            breaks = c(0, 1.5805, 2.3000, Inf),
                                                                                                            labels = c("low",
                                                                                                                       "moderate",
                                                                                                                       "high")))

stratified_PD_data_frame <- mutate(stratified_PD_data_frame, Vitamin_B12_intake_stratified = cut(Vitamin_B12,
                                                                                                            breaks = c(0, 2.98709, 7.00000, Inf),
                                                                                                            labels = c("low",
                                                                                                                       "moderate",
                                                                                                                       "high")))

stratified_PD_data_frame <- mutate(stratified_PD_data_frame, Vitamin_C_intake_stratified = cut(Vitamin_C,
                                                                                                             breaks = c(0, 80.50, 157.99, Inf),
                                                                                                             labels = c("low",
                                                                                                                        "moderate",
                                                                                                                        "high")))

stratified_PD_data_frame <- mutate(stratified_PD_data_frame, Vitamin_D_intake_stratified = cut(Vitamin_D,
                                                                                                           breaks = c(-Inf, 1.300, 2.824, Inf),
                                                                                                           labels = c("low",
                                                                                                                      "moderate",
                                                                                                                      "high")))

stratified_PD_data_frame <- mutate(stratified_PD_data_frame, Vitamin_E_intake_stratified = cut(Vitamin_E,
                                                                                                           breaks = c(0, 8.800, 14.878, Inf),
                                                                                                           labels = c("low",
                                                                                                                      "moderate",
                                                                                                                      "high")))
# final output = stratified_PD_data_frame (111 columns, 300 rows)

# startifying based on disease and dose intake:

# creating a new variable combining the disease to nutrient intake:
parkinsons_metadata_Stratified <- read_excel("+B3_stratified_na.rm_PD_metadata.xlsx")
View(parkinsons_metadata_Stratified)

# create data frame of metadata
PD_data_frame <- data.frame(parkinsons_metadata_Stratified)

# add new column to data frame that sums up nutrient intake and disease status
PD_data_frame$Total_Vitamin_A_and_Disease <- paste(PD_data_frame$Total_Vitamin_A_intake_stratified, "and", PD_data_frame$Disease)
PD_data_frame$Total_Vitamin_B1_and_Disease <- paste(PD_data_frame$Vitamin_B1_intake_stratified, "and", PD_data_frame$Disease)
PD_data_frame$Total_Vitamin_B2_and_Disease <- paste(PD_data_frame$Vitamin_B2_intake_stratified, "and", PD_data_frame$Disease)
PD_data_frame$Total_Vitamin_B3_and_Disease <- paste(PD_data_frame$Vitamin_B3_intake_stratified, "and", PD_data_frame$Disease)
PD_data_frame$Total_Vitamin_B6_and_Disease <- paste(PD_data_frame$Vitamin_B6_intake_stratified, "and", PD_data_frame$Disease)
PD_data_frame$Total_Vitamin_B12_and_Disease <- paste(PD_data_frame$Vitamin_B12_intake_stratified, "and", PD_data_frame$Disease)
PD_data_frame$Total_Vitamin_C_and_Disease <- paste(PD_data_frame$Vitamin_C_intake_stratified, "and", PD_data_frame$Disease) 
PD_data_frame$Total_Vitamin_D_and_Disease <- paste(PD_data_frame$Vitamin_D_intake_stratified, "and", PD_data_frame$Disease)
PD_data_frame$Total_Vitamin_E_and_Disease <- paste(PD_data_frame$Vitamin_E_intake_stratified, "and", PD_data_frame$Disease)
PD_data_frame$SFA_intake_stratified_and_Disease<- paste(PD_data_frame$SFA_intake_stratified, "and", PD_data_frame$Disease)
PD_data_frame$PUFA_intake_stratified_and_Disease<- paste(PD_data_frame$PUFA_intake_stratified, "and", PD_data_frame$Disease)

# create box plots of nutrients combined with disease

PD_data_frame %>%
  ggplot(aes(x = PUFA_intake_stratified_and_Disease, y = PUFA)) +
  geom_boxplot() +
  labs(title = "PUFA Intake and Disease",
       x = "Intake dose and disease status",
       y = "PUFA (g)") +
  theme_bw()

PD_data_frame %>%
  ggplot(aes(x = SFA_intake_stratified_and_Disease, y = SFA)) +
  geom_boxplot() +
  labs(title = "SFA Intake and Disease",
       x = "Intake dose and disease status",
       y = "SFA (g)") +
  theme_bw()

PD_data_frame %>%
  ggplot(aes(x = Total_Vitamin_A_and_Disease, y = Total_Vitamin_A)) +
  geom_boxplot() +
  labs(title = "(Vitamin A + Equivalents) Intake and Disease",
       x = "Intake dose and disease status",
       y = "Intake (mcg)") +
  theme_bw()

PD_data_frame %>%
  ggplot(aes(x = Total_Vitamin_B1_and_Disease, y = Vitamin_B1)) +
  geom_boxplot() +
  labs(title = "Vitamin B1 Intake and Disease",
       x = "Intake dose and disease status",
       y = "Intake (mg)") +
  theme_bw()

PD_data_frame %>%
  ggplot(aes(x = Total_Vitamin_B2_and_Disease, y = Vitamin_B2)) +
  geom_boxplot() +
  labs(title = "Vitamin B2 Intake and Disease",
       x = "Intake dose and disease status",
       y = "Intake (mg)") +
  theme_bw()

PD_data_frame %>%
  ggplot(aes(x = Total_Vitamin_B3_and_Disease, y = Niacin)) +
  geom_boxplot() +
  labs(title = "Vitamin B3 Intake and Disease",
       x = "Intake dose and disease status",
       y = "Intake (mg)") +
  theme_bw()

PD_data_frame %>%
  ggplot(aes(x = Total_Vitamin_B6_and_Disease, y = Vitamin_B6)) +
  geom_boxplot() +
  labs(title = "Vitamin B6 Intake and Disease",
       x = "Intake dose and disease status",
       y = "Intake (mg)") +
  theme_bw()

PD_data_frame %>%
  ggplot(aes(x = Total_Vitamin_B12_and_Disease, y = Vitamin_B12)) +
  geom_boxplot() +
  labs(title = "Vitamin B12 Intake and Disease",
       x = "Intake dose and disease status",
       y = "Intake (mcg)") +
  theme_bw()

PD_data_frame %>%
  ggplot(aes(x = Total_Vitamin_C_and_Disease, y = Vitamin_C)) +
  geom_boxplot() +
  labs(title = "Vitamin C Intake and Disease",
       x = "Intake dose and disease status",
       y = "Intake (mg)") +
  theme_bw()

PD_data_frame %>%
  ggplot(aes(x = Total_Vitamin_D_and_Disease, y = Vitamin_D)) +
  geom_boxplot() +
  labs(title = "Vitamin D Intake and Disease",
       x = "Intake dose and disease status",
       y = "Intake (mcg)") +
  theme_bw()

PD_data_frame %>%
  ggplot(aes(x = Total_Vitamin_E_and_Disease, y = Vitamin_E)) +
  geom_boxplot() +
  labs(title = "Vitamin E Intake and Disease",
       x = "Intake dose and disease status",
       y = "Intake (mg)") +
  theme_bw()

# export data frame as xlsx
write.xlsx(PD_data_frame, "/Users/Ayda/Desktop/stratified_PD_Disease_Combined_metadata.xlsx")

######################################## Dunn's Test: significant variables from Kruskal-Wallis test ########################################

# download Faith's PD statistics for vitamin B2 (fpd_B2.tsv)

# import file
fpd_B2 <- read.table(file = "fpd_B2.tsv", sep = '\t', header = TRUE)

# confirm Kruskal-Wallis result is significant
kruskal.test(faith_pd ~ Total_Vitamin_B2_and_Disease, data = fpd_B2)

# Dunn's test
dunn <- dunnTest(faith_pd ~ Total_Vitamin_B2_and_Disease,
         data=fpd_B2,
         method="bh")
print(dunn)

######################################## Adonis Analysis ########################################

# import metadata
meta <- read_tsv("colrm_stratified_PD_Disease_Combined_metadata.tsv", col_names = TRUE)

# import weighted/unweighted UniFrac distance matrices
unweighted_csv <- read.csv("unweighted_unifrac_distance_matrix.txt", row.names = 1)
unweighted <- as.matrix(unweighted_csv)

weighted_csv <- read.csv("weighted_unifrac_distance_matrix.txt", row.names = 1)
weighted <- as.matrix(weighted_csv)

unweighted_csv_samples <- cbind(Samples = rownames(unweighted_csv), unweighted_csv)
meta_subset <- meta %>%
  filter(SampleID %in% unweighted_csv$X)

# weighted/unweighted adonis analysis

set.seed(42)

vitamin_A_weighted_adonis <- adonis(weighted ~ Disease*Total_Vitamin_A_intake_stratified, meta_subset, permutations = 999, parallel = 1)
vitamin_A_unweighted_adonis <- adonis(unweighted ~ Disease*Total_Vitamin_A_intake_stratified, meta_subset, permutations = 999, parallel = 1)

vitamin_B1_weighted_adonis <- adonis(weighted ~ Disease*Vitamin_B1_intake_stratified, meta_subset, permutations = 999, parallel = 1)
vitamin_B1_unweighted_adonis <- adonis(unweighted ~ Disease*Vitamin_B1_intake_stratified, meta_subset, permutations = 999, parallel = 1)

vitamin_B2_weighted_adonis <- adonis(weighted ~ Disease*Vitamin_B2_intake_stratified, meta_subset, permutations = 999, parallel = 1)
vitamin_B2_unweighted_adonis <- adonis(unweighted ~ Disease*Vitamin_B2_intake_stratified, meta_subset, permutations = 999, parallel = 1)

vitamin_B3_weighted_adonis <- adonis(weighted ~ Disease*Vitamin_B3_intake_stratified, meta_subset, permutations = 999, parallel = 1)
vitamin_B3_unweighted_adonis <- adonis(unweighted ~ Disease*Vitamin_B3_intake_stratified, meta_subset, permutations = 999, parallel = 1)

vitamin_B6_weighted_adonis <- adonis(weighted ~ Disease*Vitamin_B6_intake_stratified, meta_subset, permutations = 999, parallel = 1)
vitamin_B6_unweighted_adonis <- adonis(unweighted ~ Disease*Vitamin_B6_intake_stratified, meta_subset, permutations = 999, parallel = 1)

vitamin_B12_weighted_adonis <- adonis(weighted ~ Disease*Vitamin_B12_intake_stratified, meta_subset, permutations = 999, parallel = 1)
vitamin_B12_unweighted_adonis <- adonis(unweighted ~ Disease*Vitamin_B12_intake_stratified, meta_subset, permutations = 999, parallel = 1)

vitamin_C_weighted_adonis <- adonis(weighted ~ Disease*Vitamin_C_intake_stratified, meta_subset, permutations = 999, parallel = 1)
vitamin_C_unweighted_adonis <- adonis(unweighted ~ Disease*Vitamin_C_intake_stratified, meta_subset, permutations = 999, parallel = 1)

vitamin_D_weighted_adonis <- adonis(weighted ~ Disease*Vitamin_D_intake_stratified, meta_subset, permutations = 999, parallel = 1)
vitamin_D_unweighted_adonis <- adonis(unweighted ~ Disease*Vitamin_D_intake_stratified, meta_subset, permutations = 999, parallel = 1)

vitamin_E_weighted_adonis <- adonis(weighted ~ Disease*Vitamin_E_intake_stratified, meta_subset, permutations = 999, parallel = 1)
vitamin_E_unweighted_adonis <- adonis(unweighted ~ Disease*Vitamin_E_intake_stratified, meta_subset, permutations = 999, parallel = 1)

SFA_weighted_adonis <- adonis(weighted ~ Disease*SFA_intake_stratified, meta_subset, permutations = 999, parallel = 1)
SFA_unweighted_adonis <- adonis(unweighted ~ Disease*SFA_intake_stratified, meta_subset, permutations = 999, parallel = 1)

PUFA_weighted_adonis <- adonis(weighted ~ Disease*PUFA_intake_stratified, meta_subset, permutations = 999, parallel = 1)
PUFA_unweighted_adonis <- adonis(unweighted ~ Disease*PUFA_intake_stratified, meta_subset, permutations = 999, parallel = 1)

# print results

print(vitamin_A_weighted_adonis)
print(vitamin_A_unweighted_adonis)

print(vitamin_B1_weighted_adonis)
print(vitamin_B1_unweighted_adonis)

print(vitamin_B2_weighted_adonis)
print(vitamin_B2_unweighted_adonis)

print(vitamin_B3_weighted_adonis)
print(vitamin_B3_unweighted_adonis)

print(vitamin_B6_weighted_adonis)
print(vitamin_B6_unweighted_adonis)

print(vitamin_B12_weighted_adonis)
print(vitamin_B12_unweighted_adonis)

print(vitamin_C_weighted_adonis)
print(vitamin_C_unweighted_adonis)

print(vitamin_D_weighted_adonis)
print(vitamin_D_unweighted_adonis)

print(vitamin_E_weighted_adonis)
print(vitamin_E_unweighted_adonis)

print(PUFA_weighted_adonis)
print(PUFA_unweighted_adonis)

print(SFA_weighted_adonis)
print(SFA_unweighted_adonis)

######################################## Differential Abundance Analysis: PD vs. Non-PD ########################################

# loading CRAN packages
library(tidyverse)
library(vegan)
library(ape)

# load Bioconductor packages
library(phyloseq)
library(DESeq2)

# provide custom functions: calculating relative abundance and geometric mean
calculate_relative_abundance <- function(x) x / sum(x)
calculate_gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

# setting random numbers
set.seed(711)

# importing Qiime2 data into R
biom_file <- import_biom("table-with-taxonomy.biom")
metadata  <- import_qiime_sample_data("colrm_stratified_PD_Disease_Combined_metadata.tsv")
tree      <- read_tree_greengenes ("tree.nwk")

# convert tree from multichotomous to dichotomous
tree <- multi2di(tree)

# combine all information into a single phyloseq object
physeq <- merge_phyloseq(biom_file, metadata, tree)

# assign new taxonomic column name
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class","Order", "Family",
                                 "Genus", "Species")
head(tax_table(physeq))

# setting sampling depth
sample_sums(physeq) >= 7000

# filter samples based on sampling depth
at_least_7000 <- prune_samples(sample_sums(physeq) >= 7000, physeq)
sample_sums(at_least_7000)

# getting most abundant taxa at the genus level
# remove low abundance features
total_counts <- taxa_sums(at_least_7000)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, at_least_7000)
abundant_taxa

# set taxonomic level to genus
abundant_genera <- tax_glom(abundant_taxa, taxrank = "Genus")
abundant_genera

# DESeq2 analysis by genus
deseq_genera <- phyloseq_to_deseq2(abundant_genera, ~ Disease)
geo_means_genera<- apply(counts(deseq_genera), 1, calculate_gm_mean)
deseq_genera <- estimateSizeFactors(deseq_genera, geoMeans = geo_means_genera)
deseq_genera <- DESeq(deseq_genera, fitType = "local")

diff_abund_genera<- results(deseq_genera)

# filter for data with p-value below alpha (0.05)
alpha <- 0.05
significant_genera <- as.data.frame(diff_abund_genera)
significant_genera <- filter(significant_genera, padj < alpha)

# merge tables with significant results with table of taxonomic information
genera_df <- as.data.frame(tax_table(abundant_genera))
significant_genera <- merge(significant_genera, genera_df, by = "row.names")
significant_genera <- arrange(significant_genera, log2FoldChange)

dim(significant_genera)
significant_genera

# create differential abundance genera plot
significant_genera <- filter(significant_genera, Genus != "g__")
significant_genera <- mutate(significant_genera,
                             Genus = factor(Genus, levels = Genus))

ggplot(significant_genera, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera control vs PD",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()

######################################## Relative Abundance Analysis: PD vs. Non-PD ########################################

# calculate relative abundance
at_least_7000_RA <- transform_sample_counts(at_least_7000, calculate_relative_abundance)

# remove low abundance features
total_counts <- taxa_sums(at_least_7000)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001
abundant_RA_taxa <- prune_taxa(abundant, at_least_7000_RA)

# set taxonomic rank to genus: g__Bifidobacterium  
abundant_RA_genera <- tax_glom(abundant_RA_taxa, taxrank = "Genus")

bifidobacterium <- subset_taxa(abundant_RA_genera, Genus == "g__Bifidobacterium")
otu_table(bifidobacterium)

bifidobacterium_long <- psmelt(bifidobacterium)
bifidobacterium_long

# plot with ggplot
ggplot(bifidobacterium_long, aes(x = Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Bifidobacterium",
       x     = "Disease",
       y     = "Relative abundance")

# set taxonomic rank to genus: g__Akkermansia 
abundant_RA_genera <- tax_glom(abundant_RA_taxa, taxrank = "Genus")

akkermansia <- subset_taxa(abundant_RA_genera, Genus == "g__Akkermansia")
otu_table(akkermansia)

akkermansia_long <- psmelt(akkermansia)
akkermansia_long

# plot with ggplot
ggplot(akkermansia_long, aes(x = Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Akkermansia",
       x     = "Disease",
       y     = "Relative abundance")

# set taxonomic rank to genus: g__Collinsella
abundant_RA_genera <- tax_glom(abundant_RA_taxa, taxrank = "Genus")

collinsella <- subset_taxa(abundant_RA_genera, Genus == "g__Collinsella")
otu_table(collinsella)

collinsella_long <- psmelt(collinsella)
collinsella_long

# plot with ggplot
ggplot(collinsella_long, aes(x = Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Collinsella",
       x     = "Disease",
       y     = "Relative abundance")

# set taxonomic rank to genus: g__Faecalibacterium 
abundant_RA_genera <- tax_glom(abundant_RA_taxa, taxrank = "Genus")

faecalibacterium <- subset_taxa(abundant_RA_genera, Genus == "g__Faecalibacterium")
otu_table(faecalibacterium)

faecalibacterium_long <- psmelt(faecalibacterium)
faecalibacterium_long

# plot with ggplot
ggplot(faecalibacterium_long, aes(x = Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Faecalibacterium",
       x     = "Disease",
       y     = "Relative abundance")

# set taxonomic rank to genus: g__Roseburia
abundant_RA_genera <- tax_glom(abundant_RA_taxa, taxrank = "Genus")

roseburia <- subset_taxa(abundant_RA_genera, Genus == "g__Roseburia")
otu_table(roseburia)

roseburia_long <- psmelt(roseburia)
roseburia_long

# plot with ggplot
ggplot(roseburia_long, aes(x = Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Roseburia",
       x     = "Disease",
       y     = "Relative abundance")

# set taxonomic rank to genus: g__Oscillibacter 
abundant_RA_genera <- tax_glom(abundant_RA_taxa, taxrank = "Genus")

Oscillibacter <- subset_taxa(abundant_RA_genera, Genus == "g__Oscillibacter")
otu_table(Oscillibacter)

Oscillibacter_long <- psmelt(Oscillibacter)
Oscillibacter_long

# plot with ggplot
ggplot(Oscillibacter_long, aes(x = Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Oscillibacter",
       x     = "Disease",
       y     = "Relative abundance")

# set taxonomic level to genus: g__[Eubacterium]_coprostanoligenes_group
abundant_RA_genera <- tax_glom(abundant_taxa, taxrank = "Genus")

Eubacterium_coprostanoligenes_group <- subset_taxa(abundant_RA_genera, Genus == "g__[Eubacterium]_coprostanoligenes_group")
otu_table(Eubacterium_coprostanoligenes_group)

Eubacterium_coprostanoligenes_group_long <- psmelt(Eubacterium_coprostanoligenes_group)
Eubacterium_coprostanoligenes_group_long

# plot with ggplot
ggplot(Eubacterium_coprostanoligenes_group_long, aes(x = Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Eubacterium coprostanoligenes group",
       x     = "Disease",
       y     = "Relative abundance")

######################################## Differential Abundance Analysis: Non-PD, Low vs. High Vitamin B1 ######################################## 

# filter based on metadata
vitb1 <- subset_samples(at_least_7000, Total_Vitamin_B1_and_Disease %in% c("low and Control", "high and Control"))

# remove low abundance features
total_counts <- taxa_sums(vitb1)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, vitb1)
abundant_taxa

# set taxonomic level to genus
abundant_genera <- tax_glom(abundant_taxa, taxrank = "Genus")
abundant_genera

# convert variable to 2 categories: low and high intake in Control
sample_data(abundant_genera)$Total_Vitamin_B1_and_Disease <-
  factor(sample_data(abundant_genera)$Total_Vitamin_B1_and_Disease,
         levels = c("low and Control", "high and Control"))

# DESeq2 analysis by genus
deseq <- phyloseq_to_deseq2(abundant_genera, ~ Total_Vitamin_B1_and_Disease)
geo_means <- apply(counts(deseq), 1, calculate_gm_mean)
deseq <- estimateSizeFactors(deseq, geoMeans = geo_means)
deseq <- DESeq(deseq, fitType = "local")

diff_abund <- results(deseq)

# filter for data with p-value below alpha (0.05)
alpha <- 0.05
significant <- as.data.frame(diff_abund)
significant <- filter(significant, padj < alpha)

# merge tables with significant results with table of taxonomic information
genera_df <- as.data.frame(tax_table(abundant_genera))
significant <- merge(significant, genera_df, by = "row.names")
significant <- arrange(significant, log2FoldChange)

dim(significant)
significant

# create differential abundance genera plot
significant <- filter(significant, Genus != "g__")
significant <- mutate(significant,
                      Genus = factor(Genus, levels = Genus))

ggplot(significant, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera non-PD low Vitamin B1 vs non-PD high Vitamin B1",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()

######################################## Relative Abundance Analysis: Non-PD, Low vs. High Vitamin B1 ######################################## 

# relative abundance in control gut
control_gut <- subset_samples(at_least_7000, Disease == "Control")
control_gut_RA <- transform_sample_counts(control_gut, calculate_relative_abundance)

# remove low abundance features
total_counts <- taxa_sums(control_gut)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, control_gut_RA)

# set taxonomic level to genus and select genera of interest
abundant_RA_genera <- tax_glom(abundant_taxa, taxrank = "Genus")

Muribaculaceae <- subset_taxa(abundant_RA_genera, Genus == "g__Muribaculaceae")
otu_table(Muribaculaceae)

Muribaculaceae_long <- psmelt(Muribaculaceae)
Muribaculaceae_long

Prevotellaceae_NK3B31_group <- subset_taxa(abundant_RA_genera, Genus == "g__Prevotellaceae_NK3B31_group")
otu_table(Prevotellaceae_NK3B31_group)

Prevotellaceae_NK3B31_group_long <- psmelt(Prevotellaceae_NK3B31_group)
Prevotellaceae_NK3B31_group_long

Bacteroides <- subset_taxa(abundant_RA_genera, Genus == "g__Bacteroides")
otu_table(Bacteroides)

Bacteroides_long <- psmelt(Bacteroides)
Bacteroides_long

Lachnoclostridium <- subset_taxa(abundant_RA_genera, Genus == "g__Lachnoclostridium")
otu_table(Lachnoclostridium)

Lachnoclostridium_long <- psmelt(Lachnoclostridium)
Lachnoclostridium_long

Gastranaerophilales <- subset_taxa(abundant_RA_genera, Genus == "g__Gastranaerophilales")
otu_table(Gastranaerophilales)

Gastranaerophilales_long <- psmelt(Gastranaerophilales)
Gastranaerophilales_long

# reorder factor levels
Muribaculaceae_long$Vitamin_B1_intake_stratified <- factor(Muribaculaceae_long$Vitamin_B1_intake_stratified, levels = c("low", "moderate", "high"), ordered = TRUE)
Prevotellaceae_NK3B31_group_long$Vitamin_B1_intake_stratified <- factor(Prevotellaceae_NK3B31_group_long$Vitamin_B1_intake_stratified, levels = c("low", "moderate", "high"), ordered = TRUE)
Bacteroides_long$Vitamin_B1_intake_stratified <- factor(Bacteroides_long$Vitamin_B1_intake_stratified, levels = c("low", "moderate", "high"), ordered = TRUE)
Lachnoclostridium_long$Vitamin_B1_intake_stratified <- factor(Lachnoclostridium_long$Vitamin_B1_intake_stratified, levels = c("low", "moderate", "high"), ordered = TRUE)
Gastranaerophilales_long$Vitamin_B1_intake_stratified <- factor(Gastranaerophilales_long$Vitamin_B1_intake_stratified, levels = c("low", "moderate", "high"), ordered = TRUE)

# relative abundance plots

ggplot(Muribaculaceae_long, aes(x = Vitamin_B1_intake_stratified, y = Abundance)) +
  geom_boxplot() +
  theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25), plot.title = element_text(size = 25)) +
  labs(title = "Relative abundance of Muribaculaceae",
       x     = "Vitamin B1 intake in non-PD",
       y     = "Relative abundance")

ggplot(Prevotellaceae_NK3B31_group_long, aes(x = Vitamin_B1_intake_stratified, y = Abundance)) +
  geom_boxplot() +
  theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25), plot.title = element_text(size = 25)) +
  labs(title = "Relative abundance of Prevotellaceae_NK3B31_group",
       x     = "Vitamin B1 intake in non-PD",
       y     = "Relative abundance")

ggplot(Bacteroides_long, aes(x = Vitamin_B1_intake_stratified, y = Abundance)) +
  geom_boxplot() +
  theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25), plot.title = element_text(size = 25)) +
  labs(title = "Relative abundance of Bacteroides",
       x     = "Vitamin B1 intake in non-PD",
       y     = "Relative abundance")

ggplot(Lachnoclostridium_long, aes(x = Vitamin_B1_intake_stratified, y = Abundance)) +
  geom_boxplot() +
  theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25), plot.title = element_text(size = 25)) +
  labs(title = "Relative abundance of Lachnoclostridium",
       x     = "Vitamin B1 intake in non-PD",
       y     = "Relative abundance")

ggplot(Gastranaerophilales_long, aes(x = Vitamin_B1_intake_stratified, y = Abundance)) +
  geom_boxplot() +
  theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25), plot.title = element_text(size = 25)) +
  labs(title = "Relative abundance of Gastranaerophilales",
       x     = "Vitamin B1 intake in non-PD",
       y     = "Relative abundance")

######################################## Differential Abundance Analysis: PD, Low vs. High Vitamin B1 ######################################## 

# filter based on metadata
vitb1 <- subset_samples(at_least_7000, Total_Vitamin_B1_and_Disease %in% c("low and PD", "high and PD"))

# remove low abundance features
total_counts <- taxa_sums(vitb1)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, vitb1)
abundant_taxa

# set taxonomic level to genus
abundant_genera <- tax_glom(abundant_taxa, taxrank = "Genus")
abundant_genera

# convert variable to 2 categories: low and high intake in PD
sample_data(abundant_genera)$Total_Vitamin_B1_and_Disease <-
  factor(sample_data(abundant_genera)$Total_Vitamin_B1_and_Disease,
         levels = c("low and PD", "high and PD"))

# DESeq2 analysis by genus
deseq <- phyloseq_to_deseq2(abundant_genera, ~ Total_Vitamin_B1_and_Disease)
geo_means <- apply(counts(deseq), 1, calculate_gm_mean)
deseq <- estimateSizeFactors(deseq, geoMeans = geo_means)
deseq <- DESeq(deseq, fitType = "local")

diff_abund <- results(deseq)

# filter for data with p-value below alpha (0.05)
alpha <- 0.05
significant <- as.data.frame(diff_abund)
significant <- filter(significant, padj < alpha)

# merge tables with significant results with table of taxonomic information
genera_df <- as.data.frame(tax_table(abundant_genera))
significant <- merge(significant, genera_df, by = "row.names")
significant <- arrange(significant, log2FoldChange)

dim(significant)
significant

# create differential abundance genera plot
significant <- filter(significant, Genus != "g__")
significant <- mutate(significant,
                      Genus = factor(Genus, levels = Genus))

ggplot(significant, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera PD low Vitamin B1 vs PD high Vitamin B1",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()

######################################## Differential Abundance Analysis: Non-PD, Low vs. High Vitamin B2 ######################################## 

# filter based on metadata
vitb2 <- subset_samples(at_least_7000, Total_Vitamin_B2_and_Disease %in% c("low and Control", "high and Control"))

# remove low abundance features
total_counts <- taxa_sums(vitb2)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, vitb2)
abundant_taxa

# set taxonomic level to genus
abundant_genera <- tax_glom(abundant_taxa, taxrank = "Genus")
abundant_genera

# convert variable to 2 categories: low and high intake in Control
sample_data(abundant_genera)$Total_Vitamin_B2_and_Disease <-
  factor(sample_data(abundant_genera)$Total_Vitamin_B2_and_Disease,
         levels = c("low and Control", "high and Control"))

# DESeq2 analysis by genus
deseq <- phyloseq_to_deseq2(abundant_genera, ~ Total_Vitamin_B1_and_Disease)
geo_means <- apply(counts(deseq), 1, calculate_gm_mean)
deseq <- estimateSizeFactors(deseq, geoMeans = geo_means)
deseq <- DESeq(deseq, fitType = "local")

diff_abund <- results(deseq)

# filter for data with p-value below alpha (0.05)
alpha <- 0.05
significant <- as.data.frame(diff_abund)
significant <- filter(significant, padj < alpha)

# merge tables with significant results with table of taxonomic information
genera_df <- as.data.frame(tax_table(abundant_genera))
significant <- merge(significant, genera_df, by = "row.names")
significant <- arrange(significant, log2FoldChange)

dim(significant)
significant

# create differential abundance genera plot
significant <- filter(significant, Genus != "g__")
significant <- mutate(significant,
                      Genus = factor(Genus, levels = Genus))

ggplot(significant, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera non-PD low Vitamin B2 vs non-PD high Vitamin B2",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()

######################################## Differential Abundance Analysis: PD, Low vs. High Vitamin B2 ######################################## 

# filter based on metadata
vitB2 <- subset_samples(at_least_7000, Total_Vitamin_B2_and_Disease %in% c("low and PD", "high and PD"))

# remove low abundance features
total_counts <- taxa_sums(vitB2)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, vitB2)
abundant_taxa

# set taxonomic level to genus
abundant_genera <- tax_glom(abundant_taxa, taxrank = "Genus")
abundant_genera

# convert variable to 2 categories: low and high intake in PD
sample_data(abundant_genera)$Total_Vitamin_B2_and_Disease <-
  factor(sample_data(abundant_genera)$Total_Vitamin_B2_and_Disease,
         levels = c("low and PD", "high and PD"))

# DESeq2 analysis genus
deseq <- phyloseq_to_deseq2(abundant_genera, ~ Total_Vitamin_B2_and_Disease)
geo_means <- apply(counts(deseq), 1, calculate_gm_mean)
deseq <- estimateSizeFactors(deseq, geoMeans = geo_means)
deseq <- DESeq(deseq, fitType = "local")

diff_abund <- results(deseq)

# filter for data with p-value below alpha (0.05)
alpha <- 0.05
significant <- as.data.frame(diff_abund)
significant <- filter(significant, padj < alpha)

# merge tables with significant results with table of taxonomic information
genera_df <- as.data.frame(tax_table(abundant_genera))
significant <- merge(significant, genera_df, by = "row.names")
significant <- arrange(significant, log2FoldChange)

dim(significant)
significant

# create differential abundance genera plot
significant <- filter(significant, Genus != "g__")
significant <- mutate(significant,
                      Genus = factor(Genus, levels = Genus))

ggplot(significant, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera PD low Vitamin B2 vs PD high Vitamin B2",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()

######################################################### Relative abundance: vitamin B2 in non-PD #########################################################

# relative abundance of Prevotellaceae_NK3B31_group (non PD Vitamin B6):

# calculate relative abundance
control_gut <- subset_samples(at_least_7000, Disease == "Control")
control_gut_RA <- transform_sample_counts(control_gut, calculate_relative_abundance)

# remove low abundance features
total_counts <- taxa_sums(control_gut)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, control_gut_RA)

# set taxonomic level to genus: g__Butyrivibrio
abundant_RA_genera <- tax_glom(abundant_taxa, taxrank = "Genus")
Butyrivibrio <- subset_taxa(abundant_RA_genera, Genus == "g__Butyrivibrio")
otu_table(Butyrivibrio)

Butyrivibrio_long <- psmelt(Butyrivibrio)
Butyrivibrio_long

Butyrivibrio_long$Vitamin_B2_intake_stratified <- factor(Butyrivibrio_long$Vitamin_B2_intake_stratified, levels = c("low", "moderate", "high"), ordered = TRUE)

# plot with ggplot
ggplot(Butyrivibrio_long, aes(x = Vitamin_B2_intake_stratified, y = Abundance)) +
  geom_boxplot() +
  theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25), plot.title = element_text(size = 25)) +
  labs(title = "Relative abundance of Butyrivibrio",
       x     = "Vitamin B2 intake in non-PD",
       y     = "Relative abundance")

# set taxonomic level to genus: g__Prevotellaceae_NK3B31_group
abundant_RA_genera <- tax_glom(abundant_taxa, taxrank = "Genus")
Prevotellaceae_NK3B31_group <- subset_taxa(abundant_RA_genera, Genus == "g__Prevotellaceae_NK3B31_group")
otu_table(Prevotellaceae_NK3B31_group)

Prevotellaceae_NK3B31_group_long <- psmelt(Prevotellaceae_NK3B31_group)
Prevotellaceae_NK3B31_group_long

Prevotellaceae_NK3B31_group_long$Vitamin_B2_intake_stratified <- factor(Prevotellaceae_NK3B31_group_long$Vitamin_B2_intake_stratified, levels = c("low", "moderate", "high"), ordered = TRUE)

# plot with ggplot
ggplot(Prevotellaceae_NK3B31_group_long, aes(x = Vitamin_B2_intake_stratified, y = Abundance)) +
  geom_boxplot() +
  theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25), plot.title = element_text(size = 25)) +
  labs(title = "Relative abundance of Prevotellaceae_NK3B31_group",
       x     = "Vitamin B2 intake in non-PD",
       y     = "Relative abundance")

######################################################### Relative abundance: vitamin B2 in PD #########################################################

# relative abundance of Prevotellaceae_NK3B31_group (non PD Vitamin B6)
# calculate relative abundance
control_gut <- subset_samples(at_least_7000, Disease == "PD")
control_gut_RA <- transform_sample_counts(control_gut, calculate_relative_abundance)

# remove low abundance features
total_counts <- taxa_sums(control_gut)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, control_gut_RA)

# set taxonomic level to genus: g__Anaeroplasma
abundant_RA_genera <- tax_glom(abundant_taxa, taxrank = "Genus")
Anaeroplasma <- subset_taxa(abundant_RA_genera, Genus == "g__Anaeroplasma")
otu_table(Anaeroplasma)

Anaeroplasma_long <- psmelt(Anaeroplasma)
Anaeroplasma_long

Anaeroplasma_long$Vitamin_B2_intake_stratified <- factor(Anaeroplasma_long$Vitamin_B2_intake_stratified, levels = c("low", "moderate", "high"), ordered = TRUE)

# plot with ggplot
ggplot(Anaeroplasma_long, aes(x = Vitamin_B2_intake_stratified, y = Abundance)) +
  geom_boxplot() +
  theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25), plot.title = element_text(size = 25)) +
  labs(title = "Relative abundance of Anaeroplasma",
       x     = "Vitamin B2 intake in PD",
       y     = "Relative abundance")

######################################## Differential Abundance Analysis: non-PD, Low vs. High Vitamin B6 #########################################

# filter based on metadata
vitb6 <- subset_samples(at_least_7000, Total_Vitamin_B6_and_Disease %in% c("low and Control", "high and Control"))

# remove low abundance features
total_counts <- taxa_sums(vitb6)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, vitb6)
abundant_taxa

# set taxonomic level to genus
abundant_genera <- tax_glom(abundant_taxa, taxrank = "Genus")
abundant_genera

# convert variable to 2 categories: low and high intake in Control
sample_data(abundant_genera)$Total_Vitamin_B6_and_Disease <-
  factor(sample_data(abundant_genera)$Total_Vitamin_B6_and_Disease,
         levels = c("low and Control", "high and Control"))

# DESeq2 analysis genus
deseq <- phyloseq_to_deseq2(abundant_genera, ~ Total_Vitamin_B6_and_Disease)
geo_means <- apply(counts(deseq), 1, calculate_gm_mean)
deseq <- estimateSizeFactors(deseq, geoMeans = geo_means)
deseq <- DESeq(deseq, fitType = "local")

diff_abund <- results(deseq)

# filter for data with p-value below alpha (0.05)
alpha <- 0.05
significant <- as.data.frame(diff_abund)
significant <- filter(significant, padj < alpha)

# merge tables with significant results with table of taxonomic information
genera_df <- as.data.frame(tax_table(abundant_genera))
significant <- merge(significant, genera_df, by = "row.names")
significant <- arrange(significant, log2FoldChange)

dim(significant)
significant

# create differential abundance genera plot
significant <- filter(significant, Genus != "g__")
significant <- mutate(significant,
                      Genus = factor(Genus, levels = Genus))

ggplot(significant, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera non-PD low Vitamin B6 vs non-PD high Vitamin B6",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()

######################################## Relative Abundance Analysis: non-PD, Low vs. High Vitamin B6 #########################################

# relative abundance in control gut
control_gut <- subset_samples(at_least_7000, Disease == "Control")
control_gut_RA <- transform_sample_counts(control_gut, calculate_relative_abundance)

# remove low abundance features
total_counts <- taxa_sums(control_gut)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, control_gut_RA)

# set taxonomic level to genus: g__Prevotellaceae_NK3B31_group
abundant_RA_genera <- tax_glom(abundant_taxa, taxrank = "Genus")
Prevotellaceae_NK3B31_group <- subset_taxa(abundant_RA_genera, Genus == "g__Prevotellaceae_NK3B31_group")
otu_table(Prevotellaceae_NK3B31_group)

Prevotellaceae_NK3B31_group_long <- psmelt(Prevotellaceae_NK3B31_group)
Prevotellaceae_NK3B31_group_long

# plot with ggplot
ggplot(Prevotellaceae_NK3B31_group_long, aes(x = Total_Vitamin_B6_and_Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Prevotellaceae_NK3B31_group",
       x     = "Vitamin B6 intake and disease",
       y     = "Relative abundance")

######################################## Differential Abundance Analysis: PD, Low vs. High vitamin B6 ########################################

# filter based on metadata
vitb6 <- subset_samples(at_least_7000, Total_Vitamin_B6_and_Disease %in% c("low and PD", "high and PD"))

# remove low abundance features
total_counts <- taxa_sums(vitb6)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, vitb6)
abundant_taxa

# set taxonomic level to genus
abundant_genera <- tax_glom(abundant_taxa, taxrank = "Genus")
abundant_genera

# convert variable to 2 categories: low and high intake in PD
sample_data(abundant_genera)$Total_Vitamin_B6_and_Disease <-
  factor(sample_data(abundant_genera)$Total_Vitamin_B6_and_Disease,
         levels = c("low and PD", "high and PD"))

# DESeq2 analysis by genus
deseq <- phyloseq_to_deseq2(abundant_genera, ~ Total_Vitamin_B6_and_Disease)
geo_means <- apply(counts(deseq), 1, calculate_gm_mean)
deseq <- estimateSizeFactors(deseq, geoMeans = geo_means)
deseq <- DESeq(deseq, fitType = "local")

diff_abund <- results(deseq)

# filter for data with p-value below alpha (0.05)
alpha <- 0.05
significant <- as.data.frame(diff_abund)
significant <- filter(significant, padj < alpha)

# merge tables with significant results with table of taxonomic information
genera_df <- as.data.frame(tax_table(abundant_genera))
significant <- merge(significant, genera_df, by = "row.names")
significant <- arrange(significant, log2FoldChange)

dim(significant)
significant

# create differential abundance genera plot
significant <- filter(significant, Genus != "g__")
significant <- mutate(significant,
                      Genus = factor(Genus, levels = Genus))

ggplot(significant, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera PD low Vitamin B6 vs PD high Vitamin B6",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()
