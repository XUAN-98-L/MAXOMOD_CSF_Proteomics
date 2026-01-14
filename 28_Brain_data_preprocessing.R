library(RColorBrewer)
input = "Input/intensity_mat_filtered_imputed_log2transf_norm.csv"
output_dir = "28_Brain_data/1_pre_processing"
if (!file.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}
library(openxlsx)
clinical = read.xlsx('Input/Supplemental Data 1 - Individual description of human brain samples.xlsx')
set.seed(9)
data = read.csv(input, row.names = 1)
#> setdiff(rownames(data),clinical$Sample.ID)
#[1] "H-DC110" "H-DC112"
# need to rename the clinical data to match the data
clinical$Sample.ID[clinical$Sample.ID == "HDC110"] = "H-DC110"
clinical$Sample.ID[clinical$Sample.ID == "HDC112"] = "H-DC112"

# as the data is log transformed, need to convert it back to normal scale
original_data <- 2^data

original_data = t(original_data)
original_data = as.data.frame(original_data)
library(tidyverse)
original_data = original_data %>% rownames_to_column(var = "UniProtName")
# not all proteins in uniprot_to_genename
#uniprot_to_genename = readRDS("/Users/xliu2942/Documents/Projects/MAXOMOD/Script/Validation/1_pre_processing/uniprot_to_genename.rds")

message("Reading protein to gene mapping...")
database = "Input/HUMAN_9606_idmapping_selected.tab.gz"
prot_to = readr::read_tsv(database, col_names = FALSE)
prot_to_gene <- prot_to[,c("X1","X2", "X19")]
colnames(prot_to_gene) <- c("UniProtAccession","UniProtName","gene_id")
#> length(original_data$UniProtName)
#[1] 2363
#> length(intersect(original_data$UniProtName,prot_to_gene$UniProtName))
#[1] 2317

# Find intersecting UniProt accessions between 'abu_data' and 'prot_to_gene'
inter = intersect(original_data$UniProtName, prot_to_gene$UniProtName)
# Filter 'prot_to_gene' to include only the matching UniProt accessions found in 'abu_data'
match_ids <- prot_to_gene[which(prot_to_gene$UniProtName %in% inter),]
# Sort 'match_ids' by 'UniProtAccession' to ensure it matches the order in 'abu_data'
match_ids = match_ids[match(original_data$UniProtName,match_ids$UniProtName),]
# Remove the "_HUMAN" suffix from 'UniProtName' to obtain a clean 'gene_name'
match_ids$gene_name = gsub("_HUMAN","",match_ids$UniProtName)

# define gene ID use for search
gene_ids <- match_ids$UniProtAccession

library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_info = getBM(
  attributes=c('ensembl_gene_id','uniprotswissprot','hgnc_symbol'), 
  mart = mart) 

# rename specific columns
colnames(gene_info)[colnames(gene_info) == "hgnc_symbol"] <- "SYMBOL"
colnames(gene_info)[colnames(gene_info) == "uniprotswissprot"] <- "UniProtAccession"

# keep only unique uniprot accession numbers
gene_info = gene_info %>% distinct(gene_info$UniProtAccession, .keep_all = TRUE)

df = left_join(match_ids, gene_info, by = "UniProtAccession" )
# if df$SYMBOL is "", replace it into NA
df$SYMBOL[df$SYMBOL == ""] = NA
# if df$SYMBOL is not NA, use SYMBOL instead of gene_name
df$gene_name[!is.na(df$SYMBOL)] = df$SYMBOL[!is.na(df$SYMBOL)]

#rename SYMBOL in df into better_gene_name
uniprot_to_genename <- df %>%
  dplyr::rename(better_gene_name = SYMBOL) %>%
  dplyr::select(ensembl_gene_id,gene_name,UniProtAccession,better_gene_name,UniProtName)

colnames(uniprot_to_genename)[colnames(uniprot_to_genename) == "ensembl_gene_id"] <- "gene_id"

# make rows in uniprot_to_genename the same order as in abu_data
uniprot_to_genename = uniprot_to_genename[match(original_data$UniProtName,uniprot_to_genename$UniProtName),]

abu_data = original_data
# add gene_name and uniprot_accession to abu_data
abu_data$name = uniprot_to_genename$gene_name
abu_data$name = make.unique(abu_data$name)
abu_data$ID = uniprot_to_genename$UniProtAccession
abu_data$ID = make.unique(abu_data$ID)

# if In the setdiff(original_data$UniProtName,prot_to_gene$UniProtName), need to remove these rows
abu_data = abu_data[-which(abu_data$UniProtName %in% setdiff(original_data$UniProtName,prot_to_gene$UniProtName)),]

# remove rows that has NA in uniprot_to_genename
uniprot_to_genename = uniprot_to_genename[!is.na(uniprot_to_genename$UniProtAccession),]

message("Update uniprot to genename")
saveRDS(uniprot_to_genename, file.path(output_dir, "uniprot_to_genename.rds"))
write_delim(uniprot_to_genename, file.path(output_dir, "uniprot_to_genename.txt"), delim = "\t")

saveRDS(abu_data, file.path(output_dir, "abu_data_cleaned.rds"))

message("Protein name cleaning is done!")
#===========================Summarized Experiment Object=============================
message("Making Summarized Experiment Object...")
#make summarized experiment
library(DEP)
library(SummarizedExperiment)
#selection of the clinical variables to include in the summarized experiment
experimental.design = clinical
colnames(experimental.design)[which(colnames(experimental.design) == "Sample.ID")] = "tube_id"
experimental.design$patid = experimental.design$tube_id
colnames(experimental.design)[which(colnames(experimental.design) == "Condition")] = "disease"
experimental.design[which(experimental.design$disease == "control"),"disease"] = "ctrl"
experimental.design$disease = tolower(experimental.design$disease)

# based on protein expression data, need to remove the rows that are not in the abu_data
experimental.design = experimental.design[which(experimental.design$tube_id %in% colnames(abu_data)),]
# remove rownames
rownames(experimental.design) = NULL

colnames(experimental.design)[which(colnames(experimental.design) == "Gender")] = "sex"
experimental.design[which(experimental.design$sex == "M"),"sex"] = "Male"
experimental.design[which(experimental.design$sex == "F"),"sex"] = "Female"

colnames(experimental.design)[which(colnames(experimental.design) == "Age.at.death.(years)")] = "age"

colnames(experimental.design)[which(colnames(experimental.design) == "Post.mortem.interval.(hrs)")] = "Post_mortem_interval"
experimental.design$Post_mortem_interval[which(experimental.design$Post_mortem_interval =="-")] = NA
experimental.design$Post_mortem_interval = as.numeric(experimental.design$Post_mortem_interval)

colnames(experimental.design)[which(colnames(experimental.design) == "Stratum.of.onset")] = "onset"
experimental.design$onset[which(experimental.design$onset =="-")] = "Missing"
experimental.design$onset[which(experimental.design$onset =="n/a - CTR")] = NA

colnames(experimental.design)[which(colnames(experimental.design) == "Limb.of.onset")] = "limb"
experimental.design$limb[which(experimental.design$limb =="-")] = "missing"
experimental.design$limb[which(experimental.design$limb =="n/a - CTR")] = NA

#create extra categorical variable for age based on median
m = median(experimental.design$age)
experimental.design$age_cat = rep(NA, length(experimental.design$age))
experimental.design$age_cat[experimental.design$age>=m] = paste0("over_",m)
experimental.design$age_cat[experimental.design$age<m] = paste0("under_",m)
experimental.design$age_cat = as.factor(experimental.design$age_cat)

#rename the "patid" and "disease" to label and condition, which is needed to make it work
patid_col = which(colnames(experimental.design) == "patid")
colnames(experimental.design)[patid_col] = "label"

groupby_col = which(colnames(experimental.design) == "disease") 
colnames(experimental.design)[groupby_col] = "condition"

experimental.design <- experimental.design %>%
  group_by(condition) %>%  # Group by 'condition'
  mutate(replicate = row_number()) %>%  # Number each occurrence within group
  ungroup() 

abu_data = abu_data[,c("UniProtName",experimental.design$label,"name","ID")]
# get abundance column numbers
abundance.columns <- which(colnames(abu_data) %in% experimental.design$label)

#construct the SE (summarized experiment)
se_abu_data <- make_se(abu_data, abundance.columns, experimental.design) 
saveRDS(se_abu_data, file.path(output_dir, "se_abu_data.rds"))

output_dir = "28_Brain_data/2_Missing_Inspection"
if (!file.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}
# the data is already normalize and imputed
norm_imp_MinProb = se_abu_data
saveRDS(norm_imp_MinProb, file.path(output_dir, "norm_imp_MinProb.rds"))


#===========================Summarized Experiment Object for ALS only=============================
subset = "als"
experimental.design = experimental.design[which(experimental.design$condition == subset),]
experimental.design = experimental.design[,-which(colnames(experimental.design) == "condition")]
colnames(experimental.design)[which(colnames(experimental.design) == "onset")] = "condition"

#create extra categorical variable for age based on median
m = median(experimental.design$age)
experimental.design$age_cat = rep(NA, length(experimental.design$age))
experimental.design$age_cat[experimental.design$age>=m] = paste0("over_",m)
experimental.design$age_cat[experimental.design$age<m] = paste0("under_",m)
experimental.design$age_cat = as.factor(experimental.design$age_cat)

experimental.design <- experimental.design %>%
  group_by(condition) %>%  # Group by 'condition'
  mutate(replicate = row_number()) %>%  # Number each occurrence within group
  ungroup()  # Ungroup to return a regular data frame


abu_data = abu_data[,c("UniProtName",experimental.design$label,"name","ID")]
# get abundance column numbers
abundance.columns <- which(!colnames(abu_data) %in% c("UniProtName", "name", "ID"))
#construct the SE (summarized experiment)
se_abu_data <- make_se(abu_data, abundance.columns, experimental.design) 

output_dir = "28_Brain_data/1_pre_processing_als"
if (!file.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}
saveRDS(se_abu_data, file.path(output_dir, "se_abu_data.rds"))

message("Update uniprot to genename")
saveRDS(uniprot_to_genename, file.path(output_dir, "uniprot_to_genename.rds"))
write_delim(uniprot_to_genename, file.path(output_dir, "uniprot_to_genename.txt"), delim = "\t")

output_dir = "28_Brain_data/2_Missing_Inspection_als"
if (!file.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}
# the data is already normalize and imputed
norm_imp_MinProb = se_abu_data
saveRDS(norm_imp_MinProb, file.path(output_dir, "norm_imp_MinProb.rds"))

