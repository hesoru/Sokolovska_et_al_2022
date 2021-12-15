#differential abundance PD vs Non-PD
#loading CRAN packages
library(tidyverse)
library(vegan)
library(ape)

#load Bioconductor packages
library(phyloseq)
library(DESeq2)

#provide custom function
calculate_relative_abundance <- function(x) x / sum(x)
calculate_gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

#setting random numbers
set.seed(711)

#importing Qiime2 data into R
biom_file <- import_biom("table-with-taxonomy.biom")
metadata  <- import_qiime_sample_data("colrm_stratified_PD_Disease_Combined_metadata.tsv")
tree      <- read_tree_greengenes ("tree.nwk")

#convert tree from multichotomous to dichotomous
tree <- multi2di(tree)


# Combine all information into a single phyloseq object
physeq <- merge_phyloseq(biom_file, metadata, tree)

#assign new taxonomic column name
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class","Order", "Family",
                                 "Genus", "Species")
head(tax_table(physeq))

#Setting sampling depth
sample_sums(physeq) >= 7000

#Filter samples based on sampling depth
at_least_7000 <- prune_samples(sample_sums(physeq) >= 7000, physeq)
sample_sums(at_least_7000)

#Getting most abundant taxa at the genus level
#remove low abundance features
total_counts <- taxa_sums(at_least_7000)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, at_least_7000)
abundant_taxa

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

#relative abundance
#calculate relative abundance
at_least_7000_RA <- transform_sample_counts(at_least_7000, calculate_relative_abundance)

#remove low abundant data
total_counts <- taxa_sums(at_least_7000)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.001
abundant_RA_taxa <- prune_taxa(abundant, at_least_7000_RA)

#set taxonomic rank to genus 
abundant_RA_genera <- tax_glom(abundant_RA_taxa, taxrank = "Genus")

bifidobacterium <- subset_taxa(abundant_RA_genera, Genus == "g__Bifidobacterium")
otu_table(bifidobacterium)

bifidobacterium_long <- psmelt(bifidobacterium)
bifidobacterium_long

#plot with ggplot
ggplot(bifidobacterium_long, aes(x = Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Bifidobacterium",
       x     = "Disease",
       y     = "Relative abundance")

#set taxonomic rank to genus 
abundant_RA_genera <- tax_glom(abundant_RA_taxa, taxrank = "Genus")

akkermansia <- subset_taxa(abundant_RA_genera, Genus == "g__Akkermansia")
otu_table(akkermansia)

akkermansia_long <- psmelt(akkermansia)
akkermansia_long

#plot with ggplot
ggplot(akkermansia_long, aes(x = Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Akkermansia",
       x     = "Disease",
       y     = "Relative abundance")

#set taxonomic rank to genus
abundant_RA_genera <- tax_glom(abundant_RA_taxa, taxrank = "Genus")

collinsella <- subset_taxa(abundant_RA_genera, Genus == "g__Collinsella")
otu_table(collinsella)

collinsella_long <- psmelt(collinsella)
collinsella_long

#plot with ggplot
ggplot(collinsella_long, aes(x = Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Collinsella",
       x     = "Disease",
       y     = "Relative abundance")

#set taxonomic rank to genus 
abundant_RA_genera <- tax_glom(abundant_RA_taxa, taxrank = "Genus")
faecalibacterium <- subset_taxa(abundant_RA_genera, Genus == "g__Faecalibacterium")
otu_table(faecalibacterium)

faecalibacterium_long <- psmelt(faecalibacterium)
faecalibacterium_long

#plot with ggplot
ggplot(faecalibacterium_long, aes(x = Disease, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Faecalibacterium",
       x     = "Disease",
       y     = "Relative abundance")

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
#loading CRAN packages
library(tidyverse)
library(vegan)
library(ape)

#load Bioconductor packages
library(phyloseq)
library(DESeq2)

#provide custom function
calculate_relative_abundance <- function(x) x / sum(x)
calculate_gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

#setting random numbers
set.seed(711)

#importing Qiime2 data into R
biom_file <- import_biom("table-with-taxonomy.biom")
metadata  <- import_qiime_sample_data("colrm_stratified_PD_Disease_Combined_metadata.tsv")
tree      <- read_tree_greengenes("tree.nwk")

#convert tree from multichotomous to dichotomous
tree <- multi2di(tree)

# Combine all information into a single phyloseq object
physeq <- merge_phyloseq(biom_file, metadata, tree)

#assign new taxonomic column name
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class","Order", "Family",
                                 "Genus", "Species")
head(tax_table(physeq))

#Setting sampling depth
sample_sums(physeq) >= 7000

#Filter samples based on sampling depth
at_least_7000 <- prune_samples(sample_sums(physeq) >= 7000, physeq)
sample_sums(at_least_7000)

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

#relative abundance (non PD Vitamin B1)
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
#loading CRAN packages
library(tidyverse)
library(vegan)
library(ape)

#load Bioconductor packages
library(phyloseq)
library(DESeq2)

#provide custom function
calculate_relative_abundance <- function(x) x / sum(x)
calculate_gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

#setting random numbers
set.seed(711)

#importing Qiime2 data into R
biom_file <- import_biom("table-with-taxonomy.biom")
metadata  <- import_qiime_sample_data("colrm_stratified_PD_Disease_Combined_metadata.tsv")
tree      <- read_tree_greengenes("tree.nwk")

#convert tree from multichotomous to dichotomous
tree <- multi2di(tree)

# Combine all information into a single phyloseq object
physeq <- merge_phyloseq(biom_file, metadata, tree)

#assign new taxonomic column name
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class","Order", "Family",
                                 "Genus", "Species")
head(tax_table(physeq))

#Setting sampling depth
sample_sums(physeq) >= 7000

#Filter samples based on sampling depth
at_least_7000 <- prune_samples(sample_sums(physeq) >= 7000, physeq)
sample_sums(at_least_7000)

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
#loading CRAN packages
library(tidyverse)
library(vegan)
library(ape)

#load Bioconductor packages
library(phyloseq)
library(DESeq2)

#provide custom function
calculate_relative_abundance <- function(x) x / sum(x)
calculate_gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

#setting random numbers
set.seed(711)

#importing Qiime2 data into R
biom_file <- import_biom("table-with-taxonomy.biom")
metadata  <- import_qiime_sample_data("colrm_stratified_PD_Disease_Combined_metadata.tsv")
tree      <- read_tree_greengenes("tree.nwk")

#convert tree from multichotomous to dichotomous
tree <- multi2di(tree)

# Combine all information into a single phyloseq object
physeq <- merge_phyloseq(biom_file, metadata, tree)

#assign new taxonomic column name
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class","Order", "Family",
                                 "Genus", "Species")
head(tax_table(physeq))

#Setting sampling depth
sample_sums(physeq) >= 7000

#Filter samples based on sampling depth
at_least_7000 <- prune_samples(sample_sums(physeq) >= 7000, physeq)
sample_sums(at_least_7000)

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

#relative abundance of Prevotellaceae_NK3B31_group (non PD Vitamin B6)
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
#loading CRAN packages
library(tidyverse)
library(vegan)
library(ape)

#load Bioconductor packages
library(phyloseq)
library(DESeq2)

#provide custom function
calculate_relative_abundance <- function(x) x / sum(x)
calculate_gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

#setting random numbers
set.seed(711)

#importing Qiime2 data into R
biom_file <- import_biom("table-with-taxonomy.biom")
metadata  <- import_qiime_sample_data("colrm_stratified_PD_Disease_Combined_metadata.tsv")
tree      <- read_tree_greengenes("tree.nwk")

#convert tree from multichotomous to dichotomous
tree <- multi2di(tree)

# Combine all information into a single phyloseq object
physeq <- merge_phyloseq(biom_file, metadata, tree)

#assign new taxonomic column name
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class","Order", "Family",
                                 "Genus", "Species")
head(tax_table(physeq))

#Setting sampling depth
sample_sums(physeq) >= 7000

#Filter samples based on sampling depth
at_least_7000 <- prune_samples(sample_sums(physeq) >= 7000, physeq)
sample_sums(at_least_7000)

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
