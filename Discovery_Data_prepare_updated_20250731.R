# This script cleans data for the pipeline
# Outputs include clinical_info.csv and proteinGroups.txt

# Load necessary libraries
library(readxl)
library(tidyverse)
library(dplyr)

# Load data
protein_data <- readr::read_tsv("Input/proteinGroups_CSF_08.2023.txt")
library(readxl)
#clinical_data_ALS <- read_excel('/Users/xliu2942/Documents/Projects/MAXOMOD/data/Updated_list_for_MAXOMOD_CSF_all_omics_and_clinic_2024.12.17 (clean version) - sent to Xuan.xlsx', sheet = 1)
#clinical_data_Control <- read_excel('/Users/xliu2942/Documents/Projects/MAXOMOD/data/Updated_list_for_MAXOMOD_CSF_all_omics_and_clinic_2024.12.17 (clean version) - sent to Xuan.xlsx', sheet = 2)
clinical_data <- read_excel('Input/1.MAXOMOD CSF discovery cohort clinical data - 270514  Xuan (in use final) - clean.xlsx')

# rename columns names in clinical_data
clinical_data <- clinical_data %>% dplyr::rename(
    tube_id = "Tube ID",
    patid = "CSF ID",
    disease = "Group",
    sex = "Sex",
    age = "Age at collection",
    genetics = "Genetic",
    Nfl = "NfL (pg/ml)",
    pNFh = "pNFh (pg/ml)",
    onset = "Disease onset (location of first wekness)",
    limb = "Limb Onset",
    disease_duration = "Disease duration (until collection in days)",
    age_at_onset = "Age at onset",
    progression_rate = "ALS progresion rate per month (delta ALSFRS-R / days *30)",
    progression_group = "Progression group",
    slow_vital_capacity = "slow vital capacity (%)")

# manually clean NA value
clinical_data$Nfl[clinical_data$Nfl == "na"] = NA
clinical_data$pNFh[clinical_data$pNFh == "na"] = NA
clinical_data$limb[clinical_data$limb == "na"] = NA
clinical_data$progression_rate[clinical_data$progression_rate == "na"] = NA
clinical_data$progression_group[clinical_data$progression_group == "na"] = NA
clinical_data$slow_vital_capacity[clinical_data$slow_vital_capacity == "na"] = NA
#add ECAS
clinical_data$ECAS[clinical_data$ECAS == "na"] = NA
clinical_data$ECAS = as.integer(clinical_data$ECAS)

# combind the data
# Load dplyr package
library(dplyr)
#renaming the different values for negative to just "negative"
clinical_data$genetics = as.factor(clinical_data$genetics)
levels(clinical_data$genetics) = c("C9orf72", "negative", "not_performed", "negative","VUS")

#add center variable
clinical_data$batchid = rep("munich", nrow(clinical_data))
clinical_data$batchid[grep("Ctr", clinical_data$'MAXOMOD ID')] = "goettingen" #the maxomod id's starting with "Ctr" are from goettingen
clinical_data$batchid = as.factor(clinical_data$batchid)

# change slow vital capacity from absolute value to percentage value and add limb onset info
#clinical_data = clinical_data %>% dplyr::select(-slow_vital_capacity)

#library(openxlsx)
#clinical = read.xlsx('/Users/xliu2942/Documents/Projects/MAXOMOD/data/Updated list for MAXOMOD CSF_all omics and clinic_2023.11.29 (in use) - sent to clara.xlsx', sheet = 1)

#clinical = clinical[,c("MAXOMOD.ID","Limb.Onset","slow.vital.capacity.in.%.(month.of.CSF.sampling)")]

#clinical = clinical %>% dplyr::rename("limb" = "Limb.Onset",slow_vital_capacity = "slow.vital.capacity.in.%.(month.of.CSF.sampling)")
#clinical$limb[clinical$limb == "na"] = NA

#clinical_data = clinical_data %>% dplyr::left_join(clinical, by = "MAXOMOD.ID")
#clinical_data$slow_vital_capacity[clinical_data$slow_vital_capacity == "na"] = NA
#clinical_data$pNFh[clinical_data$pNFh == "na"] = NA
#clinical_data$pNFh[clinical_data$pNFh == ""] = NA
#clinical_data$neurofilaments[clinical_data$neurofilaments == "na"] = NA


if (!dir.exists("Discovery")) {
  dir.create("Discovery")
}

output_dir <- "Discovery/data"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

write.csv(clinical_data, row.names = FALSE, file.path(output_dir, "clinical_info.csv"))
write.csv(protein_data, row.names = FALSE, file.path(output_dir, "proteinGroups.csv"))
write_delim(protein_data, file.path(output_dir, "proteinGroups.txt"), delim = "\t")
saveRDS(protein_data, file.path(output_dir, "proteinGroups.rds"))