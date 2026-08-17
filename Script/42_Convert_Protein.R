#=========================Script Description=================================
# Reads the Spectronaut protein report for the ELA-LS cohort and converts the
# "PG.ProteinGroups" column (UniProt accession, ";"-separated when a protein
# group lists more than one accession) into a protein/gene name, using the
# same two-step mapping as the "Protein Name Cleaning" section of
# 01_Pre_Processing.R:
#   1) local lookup in HUMAN_9606_idmapping_selected.tab.gz
#      (UniProtName with the "_HUMAN" suffix stripped -> gene_name)
#   2) biomaRt / Ensembl query for the official HGNC symbol, which overrides
#      gene_name when available (better_gene_name)
# and, like 01_Pre_Processing.R, also writes out the resulting
# uniprot_to_genename mapping table (.rds/.txt).
#
# Rscript Script/42_Convert_Protein.R -i Input/20220614_130757_ELA-LS_Report.xls -d Input/HUMAN_9606_idmapping_selected.tab.gz -o 42_Convert_Protein
#===========================Loading Packages=============================
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(biomaRt))

# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character", default = "Input/20220614_130757_ELA-LS_Report.xls",
              help = "Spectronaut protein report (tab-delimited text, despite the .xls extension). [default: %default]"
  ), make_option(c("-d", "--database"),
              type = "character", default = "Input/HUMAN_9606_idmapping_selected.tab.gz",
              help = "UniProt id-mapping file used to look up protein/gene names. [default: %default]"
  ), make_option(c("-o", "--output"),
              type = "character", default = "42_Convert_Protein",
              help = "Output directory. [default: %default]"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
#============================================================================
if (!file.exists(opt$input))    stop("Input report not found: ", opt$input)
if (!file.exists(opt$database)) stop("ID-mapping database not found: ", opt$database)

output_dir <- opt$output
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

#===========================Read Report=============================
message("Reading protein report: ", opt$input)
report <- readr::read_tsv(opt$input, show_col_types = FALSE)

if (!"PG.ProteinGroups" %in% colnames(report)) {
  stop("Column 'PG.ProteinGroups' not found in the report.")
}

# A protein group can list several UniProt accessions separated by ";"
# (e.g. "A0A075B6S2;A2NJV5"); keep the first accession for name lookup, same
# convention as 01_Pre_Processing.R's handling of the MaxQuant "Protein IDs"
# column (id_first = gsub(";.*", "", ...)).
report <- report %>%
  mutate(UniProtAccession = gsub(";.*", "", `PG.ProteinGroups`), .after = `PG.ProteinGroups`)

#===========================Protein Name Cleaning (same as 01_Pre_Processing.R)=============================
message("Reading protein to gene mapping...")
prot_to <- readr::read_tsv(opt$database, col_names = FALSE, show_col_types = FALSE)

# README file as below:
# https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
# select UniProtKB-AC, UniProtKB-ID and Ensembl columns
prot_to_gene <- prot_to[, c("X1", "X2", "X19")]
colnames(prot_to_gene) <- c("UniProtAccession", "UniProtName", "gene_id")

# Find intersecting UniProt accessions between 'report' and 'prot_to_gene'
inter <- intersect(report$UniProtAccession, prot_to_gene$UniProtAccession)
# Filter 'prot_to_gene' to include only the matching UniProt accessions found in 'report'
match_ids <- prot_to_gene[which(prot_to_gene$UniProtAccession %in% inter), ]
# Sort 'match_ids' by 'UniProtAccession' to ensure it matches the order in 'report'
match_ids <- match_ids[match(report$UniProtAccession, match_ids$UniProtAccession), ]
# Remove the "_HUMAN" suffix from 'UniProtName' to obtain a clean 'gene_name'
match_ids$gene_name <- gsub("_HUMAN", "", match_ids$UniProtName)

# query Ensembl BioMart for the official HGNC symbol
# NOTE: Ensembl's BioMart backend was migrated in 2023/2024 and the old
# useMart("ensembl", ...) call now 404s on some biomaRt versions. useEnsembl()
# with biomart = "genes" is the current replacement; we also try the regional
# mirrors in case the default host is unreachable. If none of them work we
# fall back to the local id-mapping only (gene_name), so the script still
# completes instead of hard-failing on a network issue.
message("Querying biomaRt for HGNC symbols...")
get_ensembl_mart <- function(dataset = "hsapiens_gene_ensembl") {
  attempts <- list(
    list(biomart = "genes", mirror = NULL),
    list(biomart = "genes", mirror = "useast"),
    list(biomart = "genes", mirror = "asia"),
    list(biomart = "ensembl", mirror = NULL)
  )
  for (a in attempts) {
    mart <- tryCatch({
      if (is.null(a$mirror)) {
        useEnsembl(biomart = a$biomart, dataset = dataset)
      } else {
        useEnsembl(biomart = a$biomart, dataset = dataset, mirror = a$mirror)
      }
    }, error = function(e) NULL)
    if (!is.null(mart)) return(mart)
  }
  NULL
}
mart <- get_ensembl_mart()

if (is.null(mart)) {
  warning("Could not reach Ensembl BioMart (tried the default host and useast/asia mirrors). ",
          "Skipping the HGNC symbol refinement step; gene_name will be based only on the local ",
          "UniProt id-mapping file (better_gene_name will be all NA).")
  gene_info <- data.frame(ensembl_gene_id = character(0), UniProtAccession = character(0),
                          SYMBOL = character(0), stringsAsFactors = FALSE)
} else {
  gene_info <- getBM(
    attributes = c('ensembl_gene_id', 'uniprotswissprot', 'hgnc_symbol'),
    mart = mart)

  # rename specific columns
  colnames(gene_info)[colnames(gene_info) == "hgnc_symbol"] <- "SYMBOL"
  colnames(gene_info)[colnames(gene_info) == "uniprotswissprot"] <- "UniProtAccession"

  # keep only unique uniprot accession numbers
  gene_info <- gene_info %>% distinct(gene_info$UniProtAccession, .keep_all = TRUE)
}

df <- left_join(match_ids, gene_info, by = "UniProtAccession")

# if df$SYMBOL is "", replace it into NA
df$SYMBOL[df$SYMBOL == ""] <- NA
# if df$SYMBOL is not NA, use SYMBOL instead of gene_name
df$gene_name[!is.na(df$SYMBOL)] <- df$SYMBOL[!is.na(df$SYMBOL)]

# rename SYMBOL in df into better_gene_name
uniprot_to_genename <- df %>%
  dplyr::rename(better_gene_name = SYMBOL) %>%
  dplyr::select(ensembl_gene_id, gene_name, UniProtAccession, better_gene_name, UniProtName)

colnames(uniprot_to_genename)[colnames(uniprot_to_genename) == "ensembl_gene_id"] <- "gene_id"

# make rows in uniprot_to_genename the same order as in report
uniprot_to_genename <- uniprot_to_genename[match(report$UniProtAccession, uniprot_to_genename$UniProtAccession), ]

# add gene_name and uniprot_accession to report
report$gene_name <- uniprot_to_genename$gene_name
report$gene_name <- make.unique(report$gene_name)
report <- relocate(report, gene_name, .after = UniProtAccession)

n_missing <- sum(is.na(uniprot_to_genename$gene_name))
if (n_missing > 0) {
  message(n_missing, " of ", nrow(report), " protein groups had no match in the id-mapping database (gene_name left NA before make.unique).")
}

message("Update uniprot to genename")
saveRDS(uniprot_to_genename, file.path(output_dir, "uniprot_to_genename.rds"))
write_delim(uniprot_to_genename, file.path(output_dir, "uniprot_to_genename.txt"), delim = "\t")

#===========================Save Output=============================
write.csv(report, file.path(output_dir, "ELA-LS_Report_with_gene_names.csv"), row.names = FALSE)
saveRDS(report, file.path(output_dir, "ELA-LS_Report_with_gene_names.rds"))

message("Done. Protein names resolved for ", nrow(report) - n_missing, "/", nrow(report), " protein groups.")
message("Output written to: ", file.path(output_dir, "ELA-LS_Report_with_gene_names.csv"))
