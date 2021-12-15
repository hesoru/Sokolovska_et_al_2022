#load packages
library(ggplot2)
library(dplyr)
library(readxl)
library(xlsx)

#select a specific set of random numbers and makes the analysis reproducible
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



