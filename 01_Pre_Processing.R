#=========================Script Description=================================
# This script is used for make pre-processing of the data, including protein name cleaning, make summarized experiment object
# Rscript 01_Pre_Processing.R -i Discovery/data/proteinGroups.txt -c Discovery/data/clinical_info.csv -o Discovery/01_Pre_Processing -e 9 -d Input/HUMAN_9606_idmapping_selected.tab.gz >output.log
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("readr"))
suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library("AnnotationDbi"))
suppressMessages(library("biomaRt"))
suppressMessages(library("DEP"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("visdat"))
suppressMessages(library("tidyr"))
suppressMessages(library("ggplot2"))
#===========================Function Definition=============================

# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "data/proteinGroups.txt",
              help = "cleaned proteinGroups (protein abundance table) file path."
  ),make_option(c("--clinical", "-c"),
              type = "character", default = "data/clinical_info.csv",
              help = "cleaned clinical data file path."
  ),make_option(c("--output", "-o"),
              type = "character", default = "01_Pre_Processing",
              help = "output directory path."
  ),make_option(c("--subset", "-s"),
                type = "character", default = NULL,
                help = "use subset samples only, subset can be ctrl or als"
  ),make_option(c("--seed", "-e"),
              type = "integer", default = 9,
              help = "set.seed"
  ),make_option(c("--database", "-d"),
              type = "character", default = "/Input/HUMAN_9606_idmapping_selected.tab.gz",
              help = "Protein name database"
  ),make_option(c("--threshold", "-t"),
              type = "numeric", default = 0.5,
              help = "filter proteins threshold"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
#============================================================================
if (is.null(opt$output)) {
  print("NO OUTPUT PATH SUPPLIED,current directory will be used!")
  output_dir <- getwd()
} else {
  output_dir <- opt$output
  if (!file.exists(output_dir)) {
    dir.create(output_dir, recursive = T)
  }
}

if (is.null(opt$input)) {
  stop("Please provide the cleaned proteinGroups file path!")
}else if (!file.exists(opt$input)) {
  stop("The cleaned proteinGroups file does not exist!")
}else{
  input = opt$input
  #load expression data
  protein_expression <- readr::read_tsv(input) #the txt file with the proteomic data
}

if (is.null(opt$clinical)) {
  stop("Please provide the cleaned clinical data file path!"
  )}else if (!file.exists(opt$clinical)) {
  stop("The cleaned clinical data file does not exist!")
}else{
  clinical = opt$clinical
  clinical = read.csv(clinical)
}

if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
}else{
  seed = opt$seed
  set.seed(seed)
}

if (is.null(opt$database)) {
  stop("Please provide the protein name database!")
}else if (!file.exists(opt$database)) {
  stop("The protein name database does not exist!")
}else{
  database = opt$database
}
#===========================Data Pre-processing=============================
#pre-filtering
protein_expression = protein_expression %>%
  dplyr::filter(is.na(`Potential contaminant`), 
                #removing possible contaminants -> '+' if possible contaminant, otherwise NA
                is.na(Reverse),
                #removing reverse proteins -> '+' if reverse, otherwise NA
                is.na(`Only identified by site`)
                #removing proteins that have only been identified by peptides containing variable PTMs, these are often seen as less reliable matches
  ) 

#LFQ intensity columns are MaxQuant has derived using the MaxLFQ algorithm  rather than raw intensities
Intensity_names <- colnames(protein_expression)[grep("LFQ intensity ", colnames(protein_expression))]

# split protein Id into SwissProt_Tremble, UniProtAccession, UniProtName
protein_expression <- protein_expression %>%
  mutate(id_first = gsub(";.*", "", `Protein IDs`)) %>%
  tidyr::separate(id_first, into = c("SwissProt_Tremble", "UniProtAccession", "UniProtName"), sep = "\\|")

# remove Biognosys_iRT peptides, which are used as standards in the MS analysis to measure peptide retention times
# remain UniProtAccession
abu_data <- protein_expression %>%
  dplyr::filter(SwissProt_Tremble == "sp", UniProtName != "Biognosys_iRT") %>%
  dplyr::select(UniProtName,UniProtAccession, Intensity_names) %>%
  dplyr::rename_with(~ gsub("LFQ intensity ", "Intensity ", .),
                     starts_with("LFQ intensity "))

tube_id_abu_data <- gsub("Intensity ","",colnames(protein_expression)[grep("Intensity ", colnames(protein_expression))])

####Loading Clinical data
# Define the required and optional columns
required_columns <- c("tube_id", "patid", "disease", "sex", "age")
optional_columns <- c("Nfl", "genetics","progression_rate","slow_vital_capacity","onset","batchid","age_at_onset","disease_duration","limb","pNFh","progression_group","ECAS")

# Check if optional columns exist in clinical and add them to the selection
columns_to_select <- c(required_columns, optional_columns[optional_columns %in% colnames(clinical)])

# Select the columns from clinical
clin <- clinical %>% dplyr::select(all_of(columns_to_select))

#make variables factor or numeric
if ("disease" %in% colnames(clin)) {clin$disease <- factor(clin$disease)}
if ("sex" %in% colnames(clin)) {clin$sex <- factor(clin$sex)}
if ("age" %in% colnames(clin)) {clin$age <- as.integer(clin$age)}

if ("Nfl" %in% colnames(clin)) {clin$Nfl <- as.numeric(clin$Nfl)}
if ("genetics" %in% colnames(clin)) {clin$genetics <- factor(clin$genetics)}
if ("progression_rate" %in% colnames(clin)) {clin$progression_rate <- as.numeric(clin$progression_rate)}
if ("slow_vital_capacity" %in% colnames(clin)) {clin$slow_vital_capacity <- as.integer(clin$slow_vital_capacity)}
if ("onset" %in% colnames(clin)) {clin$onset <- factor(clin$onset)}
if ("pNFh" %in% colnames(clin)) {clin$pNFh <- as.integer(clin$pNFh)}
if ("batchid" %in% colnames(clin)) {clin$batchid <- factor(clin$batchid)}
if ("age_at_onset" %in% colnames(clin)) {clin$age_at_onset <- as.integer(clin$age_at_onset)}
if ("disease_duration" %in% colnames(clin)) {clin$disease_duration <- as.numeric(clin$disease_duration)}
if ("limb" %in% colnames(clin)) {clin$limb <- factor(clin$limb)}
if ("progression_group" %in% colnames(clin)) {clin$progression_group <- factor(clin$progression_group)}
if ("ECAS" %in% colnames(clin)) {clin$ECAS <- as.integer(clin$ECAS)}

#create extra categorical variable for age based on median
m = median(clin$age)
clin$age_cat = rep(NA, length(clin$age))
clin$age_cat[clin$age>=m] = paste0("over_",m)
clin$age_cat[clin$age<m] = paste0("under_",m)
clin$age_cat = as.factor(clin$age_cat)

#check which patient from the clinical table is not in the abundancy table
missing_ids <- clin$tube_id[!clin$tube_id %in% tube_id_abu_data]
if (length(missing_ids) > 0) {
  message(paste0("The following patients are not in the abundancy table: ", paste(missing_ids, collapse = ", ")))
}

#reorder the clinical table so that clinical table and abundancy table are the same order and have the same patids
rownames(clin) = clin$tube_id #add the tube id as rownames for the clinical variables
clin = clin[tube_id_abu_data ,]#align clinical variables with proteomics data order
colnames(abu_data)[grep("Intensity", colnames(abu_data))] = clin$patid

#save the cleaned data
write.csv(clin,row.names = F, file.path(output_dir, "clinical_cleaned.csv"))
saveRDS(abu_data, file.path(output_dir, "abu_data_cleaned.rds"))
#===========================Protein Name Cleaning=============================
message("Reading protein to gene mapping...")
prot_to = readr::read_tsv(database, col_names = FALSE)

# README file as bellow:
# https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
# 	HUMAN_9606_idmapping_selected.tab.gz	2024-10-02 15:00	118M	 
# select UniProtKB-AC, UniProtKB-ID and Ensembl columns
prot_to_gene <- prot_to[,c("X1","X2", "X19")]
colnames(prot_to_gene) <- c("UniProtAccession","UniProtName","gene_id")

# Find intersecting UniProt accessions between 'abu_data' and 'prot_to_gene'
inter = intersect(abu_data$UniProtAccession, prot_to_gene$UniProtAccession)
# Filter 'prot_to_gene' to include only the matching UniProt accessions found in 'abu_data'
match_ids <- prot_to_gene[which(prot_to_gene$UniProtAccession %in% inter),]
# Sort 'match_ids' by 'UniProtAccession' to ensure it matches the order in 'abu_data'
match_ids = match_ids[match(abu_data$UniProtAccession,match_ids$UniProtAccession),]
# Remove the "_HUMAN" suffix from 'UniProtName' to obtain a clean 'gene_name'
match_ids$gene_name = gsub("_HUMAN","",match_ids$UniProtName)

# define gene ID use for search
gene_ids <- match_ids$UniProtAccession

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
uniprot_to_genename = uniprot_to_genename[match(abu_data$UniProtAccession,uniprot_to_genename$UniProtAccession),]

# add gene_name and uniprot_accession to abu_data
abu_data$name = uniprot_to_genename$gene_name
abu_data$name = make.unique(abu_data$name)
abu_data$ID = uniprot_to_genename$UniProtAccession
abu_data$ID = make.unique(abu_data$ID)

message("Update uniprot to genename")
saveRDS(uniprot_to_genename, file.path(output_dir, "uniprot_to_genename.rds"))
write_delim(uniprot_to_genename, file.path(output_dir, "uniprot_to_genename.txt"), delim = "\t")

saveRDS(abu_data, file.path(output_dir, "abu_data_cleaned.rds"))

message("Protein name cleaning is done!")

#===========================Summarized Experiment Object=============================
message("Making Summarized Experiment Object...")
#make summarized experiment
#selection of the clinical variables to include in the summarized experiment
experimental.design = clin
#experimental.design[which(is.na(experimental.design$onset)),"onset"] = "ctrl"

#rename the "patid" and "disease" to label and condition, which is needed to make it work
patid_col = which(colnames(experimental.design) == "patid")
colnames(experimental.design)[patid_col] = "label"

groupby_col = which(colnames(experimental.design) == "disease") 
colnames(experimental.design)[groupby_col] = "condition"

# if there is no subset column name supplied, all samples will be used
# if there is a subset column name supplied, only the subset samples will be used
# if subset is ctrl, use sex as condition
# if subset is als, use onset as condition
if (is.null(opt$subset)) {
  message("No subset column name supplied, all samples will be used!")
  # Add a 'replicate' column to the experimental.design data frame
  experimental.design <- experimental.design %>%
    group_by(condition) %>%  # Group by 'condition'
    mutate(replicate = row_number()) %>%  # Number each occurrence within group
    ungroup()  # Ungroup to return a regular data frame
}else{
  subset = opt$subset
  message(paste(subset,"samples ONLY will be used!"))
  experimental.design = experimental.design[which(experimental.design$condition == subset),]
  
  if(subset == "als"){
    # rename onset column into condition
    experimental.design = experimental.design[,-which(colnames(experimental.design) == "condition")]
    colnames(experimental.design)[which(colnames(experimental.design) == "onset")] = "condition"
    
    experimental.design <- experimental.design %>%
      group_by(condition) %>%  # Group by 'condition'
      mutate(replicate = row_number()) %>%  # Number each occurrence within group
      ungroup()  # Ungroup to return a regular data frame
    
  }else{
    experimental.design = experimental.design[,-which(colnames(experimental.design) == "condition")]
    colnames(experimental.design)[which(colnames(experimental.design) == "sex")] = "condition"
    
    experimental.design <- experimental.design %>%
      group_by(condition) %>%  # Group by 'condition'
      mutate(replicate = row_number()) %>%  # Number each occurrence within group
      ungroup()  # Ungroup to return a regular data frame
  }
}

abu_data = abu_data[,c("UniProtName",experimental.design$label,"name","ID")]
# get abundance column numbers
abundance.columns <- grep("CSF", colnames(abu_data)) 
#construct the SE (summarized experiment)
se_abu_data <- make_se(abu_data, abundance.columns, experimental.design) 
saveRDS(se_abu_data, file.path(output_dir, "se_abu_data.rds"))

#save the missing heatmap
vis_miss(as.data.frame(assay(se_abu_data)) ,show_perc = TRUE, show_perc_col = TRUE, cluster = T) 
ggsave(file.path(output_dir, "missing_vis_miss_heatmap_before.pdf"), width = 11, height = 8, units = "in") 


if (is.null(opt$threshold)) {
  message("Do not provide the threshold for filtering proteins!")
  se_abu_data_filtered = se_abu_data
}else{
  se_abu_data_filtered = DEP::filter_proteins(se_abu_data,type = "fraction", min = opt$threshold)
  saveRDS(se_abu_data_filtered, file.path(output_dir, "se_abu_data_filtered.rds"))
  vis_miss(as.data.frame(assay(se_abu_data_filtered)),show_perc = TRUE, show_perc_col = TRUE, cluster = T)
  ggsave(file.path(output_dir, "missing_vis_miss_heatmap_after.pdf"), width = 11, height = 8, units = "in")
}


# Calculate the number of missing values (0 or NA) per sample (column)
missing_per_sample <- abu_data %>%
  dplyr::select(-UniProtName) %>%  # Exclude 'UniProtName' from calculations
  summarise_all(~ sum(. == 0 | is.na(.))) %>%  # Count 0s and NAs in each column
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Missing_Count") %>%  # Convert to long format
  mutate(
    Missing_Percentage = (Missing_Count / nrow(abu_data)) * 100  # Calculate the percentage of missing values
  )

# Find samples with more than 50% missing values
#missing_per_sample[which(missing_per_sample$Missing_Percentage>50),]
write.csv(missing_per_sample, file.path(output_dir, "missing_per_sample.csv"),row.names = FALSE)

# Calculate the number and percentage of missing values per protein (row)
num_numeric_cols <- abu_data %>%
  dplyr::select(where(is.numeric)) %>%
  ncol()

# Calculate the number and percentage of missing values per protein (row)
missing_summary <- abu_data %>%
  rowwise() %>%  # Row-wise operations
  mutate(
    Missing_Count = sum(c_across(where(is.numeric)) == 0 | is.na(c_across(where(is.numeric)))),  # Count 0s and NAs in each row for numeric columns
    Missing_Percentage = (Missing_Count / num_numeric_cols) * 100  # Use pre-calculated number of numeric columns
  ) %>%
  ungroup() %>%  # Ungroup to return a regular tibble
  dplyr::select(UniProtName,name,ID, Missing_Count, Missing_Percentage)

write.csv(missing_summary, file.path(output_dir, "missing_per_protein.csv"),row.names = FALSE)


# 1. Reshape `abu_data` into long format, excluding `UniProtName`, `name`, and `ID` columns
abu_data_long <- abu_data %>%
  pivot_longer(cols = -c(UniProtName, name, ID), names_to = "label", values_to = "value")

# 2. Join with `experimental.design` to add condition information
abu_data_with_conditions <- abu_data_long %>%
  left_join(experimental.design, by = "label")

# 3. Group by both `UniProtName` and `condition`, and count zeros in each group
zero_counts_all <- abu_data_with_conditions %>%
  group_by(UniProtName, condition) %>%
  summarise(zero_count = sum(value == 0), .groups = "drop")  # Count zeros in each group

# Merge zero_counts_all into missing_summary
missing_summary_combined <- missing_summary %>%
  left_join(zero_counts_all, by = "UniProtName")

write.csv(missing_summary_combined, file.path(output_dir, "missing_per_protein_with_conditions.csv"),row.names = FALSE)