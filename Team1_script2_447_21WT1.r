#loading CRAN packages
library(tidyverse)
library(vegan)
library(ape)

#load Bioconductor packages
library(phyloseq)
library(DESeq2)

#load packages
library(ggplot2)
library(dplyr)
library(readxl)
library(xlsx)

#provide custom function
calculate_relative_abundance <- function(x) x / sum(x)
calculate_gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

#setting random numbers
set.seed(711)

#import metadata Excel file
parkinsons_metadata <- read_excel("~/4th Year Uni/TERM 2/MICB 447/Datasets/parkinsons_metadata.xlsx")
View(parkinsons_metadata)

#create data frame and clean up data:

#create data frame of metadata
PD_data_frame <- data.frame(parkinsons_metadata)

#change first column name to remove hashtag (interferes with coding)
PD_data_frame <- parkinsons_metadata %>%
  rename("SampleID" = "#SampleID")

#add new column to data frame that sums up vitamin A + equivalents
PD_data_frame$Total_Vitamin_A <- (PD_data_frame$Vitamin_A + PD_data_frame$Vitamin_A_equiv)

#remove N/A values in nutrition data prior to summary stats
na.rm_PD_data_frame <- PD_data_frame %>%
  filter(!is.na(PUFA), !is.na(SFA), !is.na(Total_Vitamin_A), !is.na(Vitamin_B1), !is.na(Vitamin_B2), !is.na(Niacin), !is.na(Vitamin_B6), !is.na(Vitamin_B12), !is.na(Vitamin_C), !is.na(Vitamin_D), !is.na(Vitamin_E))

#determining quartiles for stratification:

#split metadata into 2 data frames: control-only and PD-only
na.rm_control_only_data_frame <- filter(na.rm_PD_data_frame, Disease == "Control")
na.rm_PD_only_data_frame <- filter(na.rm_PD_data_frame, Disease == "PD")

#summary stats for each vitamin intake distribution, using control-only data frame
#(to determine quartiles for stratification)

summary(na.rm_control_only_data_frame$PUFA)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#3.90   10.40   12.90   14.22   17.46   36.20 
summary(na.rm_control_only_data_frame$SFA)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#6.594  17.936  22.600  25.061  29.900  61.900 
summary(na.rm_control_only_data_frame$Total_Vitamin_A)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#374.3  1036.3  1380.3  1622.9  1710.5 12601.6 
summary(na.rm_control_only_data_frame$Vitamin_B1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.4000  0.9726  1.2000  1.3078  1.4262  6.8563 
summary(na.rm_control_only_data_frame$Vitamin_B2)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.600   1.200   1.600   1.634   1.990   4.500 
summary(na.rm_control_only_data_frame$Niacin)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#8.959  15.183  19.901  20.457  23.400  50.400
summary(na.rm_control_only_data_frame$Vitamin_B6)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9535  1.5805  1.9193  2.0082  2.3000  4.8000 
summary(na.rm_control_only_data_frame$Vitamin_B12)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.02695  2.98709  4.80000  5.44874  7.00000 27.50000 
summary(na.rm_control_only_data_frame$Vitamin_C)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#18.55   80.50  112.80  127.71  157.99  334.02
summary(na.rm_control_only_data_frame$Vitamin_D)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   1.300   2.025   2.387   2.824   8.815 
summary(na.rm_control_only_data_frame$Vitamin_E)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#3.771   8.800  11.100  12.172  14.878  31.700 

#make boxplots of nutrient intake distributions, confirm intake Q1/Q3 for control match summary stats:

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

#stratify nutrient intake data based on intake Q1/Q3 for control:

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
#final output = stratified_PD_data_frame (111 columns, 300 rows)

#startifying based on disease and dose intake

#creating a new variable combining the disease to nutrient intake:
parkinsons_metadata_Stratified <- read_excel("+B3_stratified_na.rm_PD_metadata.xlsx")
View(parkinsons_metadata_Stratified)

#create data frame of metadata
PD_data_frame <- data.frame(parkinsons_metadata_Stratified)

#add new column to data frame that sums up nutrient intake and disease status
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

#create box plots of nutrients combined with disease

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

#export data frame as xlsx
write.xlsx(PD_data_frame, "/Users/Ayda/Desktop/stratified_PD_Disease_Combined_metadata.xlsx")


#differential abundance for PD vs non-PD

#importing QIIME2 data into R
biom_file <- import_biom("table-with-taxonomy.biom")
metadata  <- import_qiime_sample_data("colrm_stratified_PD_Disease_Combined_metadata.tsv")
tree      <- read_tree_greengenes ("tree.nwk")

#convert tree from multichotomous to dichotomous
tree <- multi2di(tree)

#combine all information into a single phyloseq object
physeq <- merge_phyloseq(biom_file, metadata, tree)

#assign new taxonomic column name
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class","Order", "Family",
                                 "Genus", "Species")
head(tax_table(physeq))

#Setting sampling depth
sample_sums(physeq) >= 7000

#filter samples based on sampling depth
at_least_7000 <- prune_samples(sample_sums(physeq) >= 7000, physeq)
sample_sums(at_least_7000)

#remove low abundance features
total_counts <- taxa_sums(at_least_7000)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, at_least_7000)
abundant_taxa

#set taxonomic level to species
abundant_species <- tax_glom(abundant_taxa, taxrank = "Species")
abundant_species

#DESeq2 analysis species
deseq_genus_species <- phyloseq_to_deseq2(abundant_species, ~ Disease)
geo_means_species <- apply(counts(deseq), 1, calculate_gm_mean)
deseq_species <- estimateSizeFactors(deseq, geoMeans = geo_means_species)
deseq_species <- DESeq(deseq_species, fitType = "local")

diff_abund_species <- results(deseq_species)

#filter for data with p-value below alpha (0.05)
alpha <- 0.05
significant_species <- as.data.frame(diff_abund_species)
significant_species <- filter(significant_species, padj < alpha)

#merge tables with significant results with table of taxonomic information
species_df <- as.data.frame(tax_table(abundant_species))
significant_species <- merge(significant_species, species_df, by = "row.names")
significant_species <- arrange(significant_species, log2FoldChange)

dim(significant_species)
significant_species

#set taxonomic level to family
abundant_family <- tax_glom(abundant_taxa, taxrank = "Family")
abundant_family

#DESeq2 analysis family
deseq_family <- phyloseq_to_deseq2(abundant_family, ~ Disease)
geo_means_family <- apply(counts(deseq_family), 1, calculate_gm_mean)
deseq_family <- estimateSizeFactors(deseq_family, geoMeans = geo_means_family)
deseq_family <- DESeq(deseq_family, fitType = "local")

diff_abund_family <- results(deseq_family)

#filter for data with p-value below alpha (0.05)
alpha <- 0.05
significant_family <- as.data.frame(diff_abund_family)
significant_family <- filter(significant_family, padj < alpha)

#merge tables with significant results with table of taxonomic information
family_df <- as.data.frame(tax_table(abundant_family))
significant_family <- merge(significant_family, family_df, by = "row.names")
significant_family <- arrange(significant_family, log2FoldChange)

dim(significant_family)
significant_family

#create differential abundance family plot
significant_family <- mutate(significant_family,
                             Family = factor(Family, levels = unique(Family)))

ggplot(significant_family, aes(x = log2FoldChange, y = Family)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant family control vs PD",
       x = expression(log[2]~fold~change),
       y = "Family") +
  theme_bw()

#getting most abundant taxa at the genus level
#set taxonomic level to genus
abundant_genera <- tax_glom(abundant_taxa, taxrank = "Genus")
abundant_genera

#DESeq2 analysis genus
deseq_genera <- phyloseq_to_deseq2(abundant_genera, ~ Disease)
geo_means_genera<- apply(counts(deseq_genera), 1, calculate_gm_mean)
deseq_genera <- estimateSizeFactors(deseq_genera, geoMeans = geo_means_genera)
deseq_genera <- DESeq(deseq_genera, fitType = "local")

diff_abund_genera<- results(deseq_genera)

#filter for data with p-value below alpha (0.05)
alpha <- 0.05
significant_genera <- as.data.frame(diff_abund_genera)
significant_genera <- filter(significant_genera, padj < alpha)

#merge tables with significant results with table of taxonomic information
genera_df <- as.data.frame(tax_table(abundant_genera))
significant_genera <- merge(significant_genera, genera_df, by = "row.names")
significant_genera <- arrange(significant_genera, log2FoldChange)

dim(significant_genera)
significant_genera

#create differential abundance genera plot
significant_genera <- filter(significant_genera, Genus != "g__")
significant_genera <- mutate(significant_genera,
                             Genus = factor(Genus, levels = Genus))

ggplot(significant_genera, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera control vs PD",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()

#set taxonomic rank to genus 
abundant_RA_genera <- tax_glom(abundant_RA_taxa, taxrank = "Genus")

roseburia <- subset_taxa(abundant_RA_genera, Genus == "g__Roseburia")
otu_table(roseburia)

roseburia_long <- psmelt(roseburia)
roseburia_long

#plot with ggplot
ggplot(roseburia_long, aes(x = Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Roseburia",
       x     = "Disease",
       y     = "Relative abundance")


#Oscillibacter relative abundance

#set taxonomic rank to genus 
abundant_RA_genera <- tax_glom(abundant_RA_taxa, taxrank = "Genus")

Oscillibacter <- subset_taxa(abundant_RA_genera, Genus == "g__Oscillibacter")
otu_table(Oscillibacter)

Oscillibacter_long <- psmelt(Oscillibacter)
Oscillibacter_long

#plot with ggplot
ggplot(Oscillibacter_long, aes(x = Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Oscillibacter",
       x     = "Disease",
       y     = "Relative abundance")

#Eubacterium_coprostanoligenes_group relative abundance

#set taxonomic level to genus
abundant_RA_genera <- tax_glom(abundant_taxa, taxrank = "Genus")

Eubacterium_coprostanoligenes_group <- subset_taxa(abundant_RA_genera, Genus == "g__[Eubacterium]_coprostanoligenes_group")
otu_table(Eubacterium_coprostanoligenes_group)

Eubacterium_coprostanoligenes_group_long <- psmelt(Eubacterium_coprostanoligenes_group)
Eubacterium_coprostanoligenes_group_long

#plot with ggplot
ggplot(Eubacterium_coprostanoligenes_group_long, aes(x = Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Eubacterium coprostanoligenes group",
       x     = "Disease",
       y     = "Relative abundance")



#differential abundance for non-PD low Vitamin B1 vs non-PD high Vitamin B1

#filter based on metadata
vitb1 <- subset_samples(at_least_7000, Total_Vitamin_B1_and_Disease %in% c("low and Control", "high and Control"))

#remove low abundance features
total_counts <- taxa_sums(vitb1)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, vitb1)
abundant_taxa

#set taxonomic level to genus
abundant_genera <- tax_glom(abundant_taxa, taxrank = "Genus")
abundant_genera

#convert variable to 2 categories
sample_data(abundant_genera)$Total_Vitamin_B1_and_Disease <-
  factor(sample_data(abundant_genera)$Total_Vitamin_B1_and_Disease,
         levels = c("low and Control", "high and Control"))

#DESeq2 analysis genus
deseq <- phyloseq_to_deseq2(abundant_genera, ~ Total_Vitamin_B1_and_Disease)
geo_means <- apply(counts(deseq), 1, calculate_gm_mean)
deseq <- estimateSizeFactors(deseq, geoMeans = geo_means)
deseq <- DESeq(deseq, fitType = "local")

diff_abund <- results(deseq)

#filter for data with p-value below alpha (0.05)
alpha <- 0.05
significant <- as.data.frame(diff_abund)
significant <- filter(significant, padj < alpha)

#merge tables with significant results with table of taxonomic information
genera_df <- as.data.frame(tax_table(abundant_genera))
significant <- merge(significant, genera_df, by = "row.names")
significant <- arrange(significant, log2FoldChange)

dim(significant)
significant

#create differential abundance genera plot
significant <- filter(significant, Genus != "g__")
significant <- mutate(significant,
                      Genus = factor(Genus, levels = Genus))

ggplot(significant, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera non-PD low Vitamin B1 vs non-PD high Vitamin B1",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()

#relative abundance for non-PD low Vitamin B1 vs non-PD high Vitamin B1
#calculate relative abundance
vitb1_RA <- transform_sample_counts(vitb1, calculate_relative_abundance)

#remove low abundance features
total_counts <- taxa_sums(vitb1)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, vitb1_RA)

#set taxonomic level to genus
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

#relative abundance plots
ggplot(Muribaculaceae_long, aes(x = Total_Vitamin_B1_and_Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Muribaculaceae",
       x     = "Vitamin B1 intake and disease",
       y     = "Relative abundance")

ggplot(Prevotellaceae_NK3B31_group_long, aes(x = Total_Vitamin_B1_and_Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Prevotellaceae_NK3B31_group",
       x     = "Vitamin B1 intake and disease",
       y     = "Relative abundance")

ggplot(Bacteroides_long, aes(x = Total_Vitamin_B1_and_Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Bacteroides",
       x     = "Vitamin B1 intake and disease",
       y     = "Relative abundance")

ggplot(Lachnoclostridium_long, aes(x = Total_Vitamin_B1_and_Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Lachnoclostridium",
       x     = "Vitamin B1 intake and disease",
       y     = "Relative abundance")

ggplot(Gastranaerophilales_long, aes(x = Total_Vitamin_B1_and_Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Gastranaerophilales",
       x     = "Vitamin B1 intake and disease",
       y     = "Relative abundance")

#differential abundance for PD low Vitamin B1 vs PD high Vitamin B1

#filter based on metadata
vitb1 <- subset_samples(at_least_7000, Total_Vitamin_B1_and_Disease %in% c("low and PD", "high and PD"))

#remove low abundance features
total_counts <- taxa_sums(vitb1)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, vitb1)
abundant_taxa

#set taxonomic level to genus
abundant_genera <- tax_glom(abundant_taxa, taxrank = "Genus")
abundant_genera

#convert variable to 2 categories
sample_data(abundant_genera)$Total_Vitamin_B1_and_Disease <-
  factor(sample_data(abundant_genera)$Total_Vitamin_B1_and_Disease,
         levels = c("low and PD", "high and PD"))

#DESeq2 analysis genus
deseq <- phyloseq_to_deseq2(abundant_genera, ~ Total_Vitamin_B1_and_Disease)
geo_means <- apply(counts(deseq), 1, calculate_gm_mean)
deseq <- estimateSizeFactors(deseq, geoMeans = geo_means)
deseq <- DESeq(deseq, fitType = "local")

diff_abund <- results(deseq)

#filter for data with p-value below alpha (0.05)
alpha <- 0.05
significant <- as.data.frame(diff_abund)
significant <- filter(significant, padj < alpha)

#merge tables with significant results with table of taxonomic information
genera_df <- as.data.frame(tax_table(abundant_genera))
significant <- merge(significant, genera_df, by = "row.names")
significant <- arrange(significant, log2FoldChange)

dim(significant)
significant

#create differential abundance genera plot
significant <- filter(significant, Genus != "g__")
significant <- mutate(significant,
                      Genus = factor(Genus, levels = Genus))

ggplot(significant, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera PD low Vitamin B1 vs PD high Vitamin B1",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()

#differential abundance for non-PD low Vitamin B6 vs non-PD high Vitamin B6

#filter based on metadata
vitb6 <- subset_samples(at_least_7000, Total_Vitamin_B6_and_Disease %in% c("low and Control", "high and Control"))

#remove low abundance features
total_counts <- taxa_sums(vitb6)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, vitb6)
abundant_taxa

#set taxonomic level to genus
abundant_genera <- tax_glom(abundant_taxa, taxrank = "Genus")
abundant_genera

#convert variable to 2 categories
sample_data(abundant_genera)$Total_Vitamin_B6_and_Disease <-
  factor(sample_data(abundant_genera)$Total_Vitamin_B6_and_Disease,
         levels = c("low and Control", "high and Control"))

#DESeq2 analysis genus
deseq <- phyloseq_to_deseq2(abundant_genera, ~ Total_Vitamin_B6_and_Disease)
geo_means <- apply(counts(deseq), 1, calculate_gm_mean)
deseq <- estimateSizeFactors(deseq, geoMeans = geo_means)
deseq <- DESeq(deseq, fitType = "local")

diff_abund <- results(deseq)

#filter for data with p-value below alpha (0.05)
alpha <- 0.05
significant <- as.data.frame(diff_abund)
significant <- filter(significant, padj < alpha)

#merge tables with significant results with table of taxonomic information
genera_df <- as.data.frame(tax_table(abundant_genera))
significant <- merge(significant, genera_df, by = "row.names")
significant <- arrange(significant, log2FoldChange)

dim(significant)
significant

#create differential abundance genera plot
significant <- filter(significant, Genus != "g__")
significant <- mutate(significant,
                      Genus = factor(Genus, levels = Genus))

ggplot(significant, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera non-PD low Vitamin B6 vs non-PD high Vitamin B6",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()

#relative abundance of Prevotellaceae_NK3B31_group on non-PD low Vitamin B6 vs non-PD high Vitamin B6
#calculate relative abundance
vitb6_RA <- transform_sample_counts(vitb6, calculate_relative_abundance)

#remove low abundance features
total_counts <- taxa_sums(vitb6)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, vitb6_RA)

#set taxonomic level to genus
abundant_RA_genera <- tax_glom(abundant_taxa, taxrank = "Genus")
Prevotellaceae_NK3B31_group <- subset_taxa(abundant_RA_genera, Genus == "g__Prevotellaceae_NK3B31_group")
otu_table(Prevotellaceae_NK3B31_group)

Prevotellaceae_NK3B31_group_long <- psmelt(Prevotellaceae_NK3B31_group)
Prevotellaceae_NK3B31_group_long

#plot with ggplot
ggplot(Prevotellaceae_NK3B31_group_long, aes(x = Total_Vitamin_B6_and_Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Prevotellaceae_NK3B31_group",
       x     = "Vitamin B6 intake and disease",
       y     = "Relative abundance")

#differential abundance for PD low Vitamin B6 vs PD high Vitamin B6

#filter based on metadata
vitb6 <- subset_samples(at_least_7000, Total_Vitamin_B6_and_Disease %in% c("low and PD", "high and PD"))

#remove low abundance features
total_counts <- taxa_sums(vitb6)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, vitb6)
abundant_taxa

#set taxonomic level to genus
abundant_genera <- tax_glom(abundant_taxa, taxrank = "Genus")
abundant_genera

#convert variable to 2 categories
sample_data(abundant_genera)$Total_Vitamin_B6_and_Disease <-
  factor(sample_data(abundant_genera)$Total_Vitamin_B6_and_Disease,
         levels = c("low and PD", "high and PD"))

#DESeq2 analysis genus
deseq <- phyloseq_to_deseq2(abundant_genera, ~ Total_Vitamin_B6_and_Disease)
geo_means <- apply(counts(deseq), 1, calculate_gm_mean)
deseq <- estimateSizeFactors(deseq, geoMeans = geo_means)
deseq <- DESeq(deseq, fitType = "local")

diff_abund <- results(deseq)

#filter for data with p-value below alpha (0.05)
alpha <- 0.05
significant <- as.data.frame(diff_abund)
significant <- filter(significant, padj < alpha)

#merge tables with significant results with table of taxonomic information
genera_df <- as.data.frame(tax_table(abundant_genera))
significant <- merge(significant, genera_df, by = "row.names")
significant <- arrange(significant, log2FoldChange)

dim(significant)
significant

#create differential abundance genera plot
significant <- filter(significant, Genus != "g__")
significant <- mutate(significant,
                      Genus = factor(Genus, levels = Genus))

ggplot(significant, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera PD low Vitamin B6 vs PD high Vitamin B6",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()
