# This script cleans data for the pipeline
# Outputs include clinical_info.csv and proteinGroups.txt

# Load necessary libraries
library(readxl)
library(tidyverse)
library(dplyr)

# Load data
protein_data <- readr::read_tsv("Input/proteinGroups.txt")

library(openxlsx)
#demographic_data = read.xlsx('/Input/validataion/MAXOMOD validation cohort clinical data - Xuan_20241203.xlsx', sheet = 1)
demographic_data = read.xlsx('Input/1.MAXOMOD CSF validation cohort clinical data - 250414 Xuan (in use final).xlsx', sheet = 1)
# select the columns that are needed
demographic_data = demographic_data %>% 
  dplyr::rename("Patient.ID" = "ID",
                "batchid" = "center",
                "disease" = "Group",
                "sex" = "Sex",
                "age" = "Age.at.collection",
                "Nfl" = "NfL.(pg/ml)",
                "pNFh" = "pNFh.(pg/ml)",
                "genetics" = "Genetic",
                "age_at_onset" = "Age.at.onset",
                "onset" = "Disease.onset.(location.of.first.wekness)",
                "limb" = "Limb.onset",
                "progression_rate" = "ALS.progresion.rate.per.month.(delta.ALSFRS-R./.days.*30)",
                "slow_vital_capacity" = "VC%",
                "disease_duration" = "Disease.duration.(until.collection.in.days)",
                "progression_group" = "Progression.group")

demographic_data = demographic_data %>% dplyr::select(`Patient.ID`,disease,sex,age,Nfl,pNFh,genetics,progression_rate,progression_group,slow_vital_capacity,onset,limb,age_at_onset,disease_duration,batchid,ECAS)

# if it's "na" in clinical$limb, replace it with NA
demographic_data$limb[demographic_data$limb == "na"] = NA
demographic_data$batchid = tolower(demographic_data$batchid)
demographic_data$ECAS[demographic_data$ECAS == "na"] = NA

clinical_data <- read.csv("Input/SampleTracking.csv")
#demographic_data <- read_excel("/Input/validataion/validation cohort sex and age.xlsx")

# Join clinical data with demographic data by Patient.ID
clinical_data <- inner_join(clinical_data, demographic_data, by = c("Patient.ID" = "Patient.ID"))

# Remove rows with missing MaxQuant file names
clinical_data <- clinical_data[clinical_data$MaxQuant.File.Name != "Missing", ]
row.names(clinical_data) <- NULL

# Remove Lsmbo.ID column
clinical_data <- clinical_data[ , !colnames(clinical_data) %in% "Lsmbo.ID"]

# remove duplicated samples,sort by Patient.ID
clinical_data <- clinical_data[!clinical_data$`Patient.ID` %in% c(1001101, 2721527, 2721900, 2798778, 2805592, 2814942, 2844960,2821491), ]
row.names(clinical_data) <- NULL


# Add formatted MaxQuant match ID with leading zeroes and create tube_id for matching
clinical_data$MaxQuant_match <- sprintf("%03d", as.integer(clinical_data$MaxQuant.File.Name))
clinical_data$tube_id <- paste0(clinical_data$Group, clinical_data$MaxQuant_match)

# Create unique patient ID and format disease, sex, and age columns
clinical_data$patid <- paste0("CSF", rownames(clinical_data))
clinical_data$disease <- tolower(clinical_data$Group)
clinical_data$genetics = as.factor(clinical_data$genetics)
#levels(clinical_data$genetics) = c("C9orf72", "negative", "not_performed", "negative", "ROCK", "SOD1", "SOD1_FIG4", "VUS")
levels(clinical_data$genetics) = c("C9orf72", "negative", "not_performed", "negative","SOD1", "SOD1_FIG4", "VUS")

########################### Match clinical data with protein data ###########################

# Extract columns with numeric suffixes from protein_data and isolate their numeric parts
numeric_columns <- grep("\\d+", colnames(protein_data), value = TRUE)
col_map <- data.frame(numeric_columns)
col_map$number <- sprintf("%03d", as.integer(gsub("\\D", "", col_map$numeric_columns)))
col_map$prefix <- gsub("\\d+", "", col_map$numeric_columns)

# Match MaxQuant file IDs to tube IDs in clinical_data
col_map <- col_map %>%
  left_join(clinical_data %>% dplyr::select(MaxQuant_match, tube_id), by = c("number" = "MaxQuant_match"))

# Identify columns without matching tube IDs and remove them from protein_data
unmatched_cols <- col_map[is.na(col_map$tube_id), ]
col_map <- col_map[!is.na(col_map$tube_id), ]
protein_data <- protein_data[ , !colnames(protein_data) %in% unmatched_cols$numeric_columns]

# Rename columns in protein_data using tube_id values
col_map <- col_map %>%
  rowwise() %>%
  mutate(mapped_column_name = gsub(number, tube_id, numeric_columns, fixed = TRUE)) %>%
  ungroup()
colnames(protein_data)[colnames(protein_data) %in% col_map$numeric_columns] <- col_map$mapped_column_name

# change into integer
clinical_data$MaxQuant.File.Name = as.integer(clinical_data$MaxQuant.File.Name)
clinical_data$MaxQuant_match = as.integer(clinical_data$MaxQuant_match)
########################### Save the cleaned data ###########################

if (!dir.exists("Validation")) {
  dir.create("Validation")
}

output_dir <- "Validation/data"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

write.csv(clinical_data, row.names = FALSE, file.path(output_dir, "clinical_info.csv"))
write.csv(protein_data, row.names = FALSE, file.path(output_dir, "proteinGroups.csv"))
write_delim(protein_data, file.path(output_dir, "proteinGroups.txt"), delim = "\t")
saveRDS(protein_data, file.path(output_dir, "proteinGroups.rds"))