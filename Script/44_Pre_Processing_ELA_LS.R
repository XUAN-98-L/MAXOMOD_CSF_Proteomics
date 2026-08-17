#=========================Script Description=================================
# Extracts the protein quantity ("...htrms.PG.Quantity") columns from the
# Spectronaut report produced by Script/42_Convert_Protein.R
# (42_Convert_Protein/ELA-LS_Report_with_gene_names.csv), builds a clean
# protein x sample expression table, and runs the same pre-processing steps
# used elsewhere in this repo (Script/01_Pre_Processing.R and
# Script/02_Missing_Inspection.R):
#   1) identify the "[N] <sample>.htrms.PG.Quantity" columns and rename them
#      to plain sample names (e.g. "[16] ELA-LS_01.htrms.PG.Quantity" -> "ELA-LS_01")
#   2) assemble a clean protein x sample quantity table (+ QC summaries of
#      missingness per sample / per protein)
#   3) build a DEP SummarizedExperiment and filter proteins on the fraction of
#      samples with a valid (non-missing) value (DEP::filter_proteins, same as
#      01_Pre_Processing.R)
#   4) VSN-normalize and impute missing values with DEP::impute(fun="MinProb")
#      (same as 02_Missing_Inspection.R)
#
# Sample groups are inferred from the sample name prefix: "Ctl_*" -> "Control",
# "ELA-LS_*" -> "ELA-LS", "ELA-SS_*" -> "ELA-SS" (same convention used in
# Script/43_Predict_ELA_LS_Subtype.R).
#
# Rscript Script/44_Pre_Processing_ELA_LS.R -i 42_Convert_Protein/ELA-LS_Report_with_gene_names.csv -o 44_Pre_Processing_ELA_LS -e 9 -t 0.5 -q 0.01 >output.log
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("DEP"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("visdat"))
suppressMessages(library("ggplot2"))

# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "42_Convert_Protein/ELA-LS_Report_with_gene_names.csv",
              help = "Annotated Spectronaut report from Script/42_Convert_Protein.R. [default: %default]"
  ), make_option(c("--output", "-o"),
              type = "character", default = "44_Pre_Processing_ELA_LS",
              help = "Output directory. [default: %default]"
  ), make_option(c("--seed", "-e"),
              type = "integer", default = 9,
              help = "set.seed [default: %default]"
  ), make_option(c("--threshold", "-t"),
              type = "numeric", default = 0.5,
              help = "Filter proteins: minimum fraction of samples with a valid (non-missing) value to keep a protein, passed to DEP::filter_proteins(type = 'fraction', min = ...). Set to 0 to skip filtering. [default: %default]"
  ), make_option(c("--scalar", "-q"),
              type = "numeric", default = 0.01,
              help = "MinProb imputation quantile (0 < q < 1). [default: %default]"
  ), make_option(c("--normalization", "-n"),
              type = "logical", default = TRUE,
              help = "Whether to VSN-normalize before imputation. [default: %default]"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
#============================================================================
if (is.null(opt$output)) {
  message("NO OUTPUT PATH SUPPLIED, current directory will be used!")
  output_dir <- getwd()
} else {
  output_dir <- opt$output
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
}

if (is.null(opt$input)) {
  stop("Please provide the ELA-LS report file path!")
} else if (!file.exists(opt$input)) {
  stop("The ELA-LS report file does not exist: ", opt$input)
}

set.seed(opt$seed)

#===========================Read Report & Extract PG.Quantity Columns=============================
message("Reading report: ", opt$input)
report <- readr::read_csv(opt$input, show_col_types = FALSE, na = c("", "NA", "NaN"))
report <- as.data.frame(report, stringsAsFactors = FALSE)

if (!all(c("gene_name", "UniProtAccession") %in% colnames(report))) {
  stop("Report is missing 'gene_name'/'UniProtAccession' columns -- run Script/42_Convert_Protein.R first.")
}

# Remove the Biognosys iRT standard row (retention-time calibration peptides,
# not a real protein) -- same exclusion 01_Pre_Processing.R applies to
# "Biognosys_iRT". In this Spectronaut report it shows up as
# PG.ProteinGroups == "iRT-Kit_WR_fusion" with gene_name left NA by
# 42_Convert_Protein.R. Also drops any other row with NA gene_name as a
# safety net, since DEP::make_se() cannot have NA in the "name" column.
irt_rows <- grepl("iRT", report$PG.ProteinGroups, ignore.case = TRUE) | is.na(report$gene_name)
if (any(irt_rows)) {
  message(sum(irt_rows), " iRT standard / unannotated row(s) removed: ",
          paste(unique(report$PG.ProteinGroups[irt_rows]), collapse = ", "))
  report <- report[!irt_rows, , drop = FALSE]
}

# columns look like "[16] ELA-LS_01.htrms.PG.Quantity" -> sample "ELA-LS_01"
quant_cols <- grep("\\.htrms\\.PG\\.Quantity$", colnames(report), value = TRUE)
if (length(quant_cols) == 0) stop("No '...htrms.PG.Quantity' columns found in the report.")

sample_names <- gsub("^\\[[0-9]+\\]\\s*", "", quant_cols)
sample_names <- gsub("\\.htrms\\.PG\\.Quantity$", "", sample_names)

message(length(quant_cols), " PG.Quantity (expression) columns found, e.g. '",
        quant_cols[1], "' -> sample '", sample_names[1], "'")

# infer group from sample name prefix: Ctl_* -> Control, ELA-LS_* / ELA-SS_* kept as-is
group <- gsub("_[0-9]+$", "", sample_names)
group[group == "Ctl"] <- "Control"

#===========================Build Clean Protein x Sample Table=============================
quant_mat <- report[, quant_cols, drop = FALSE]
colnames(quant_mat) <- sample_names
quant_mat <- as.data.frame(lapply(quant_mat, as.numeric), check.names = FALSE)

clean_data <- data.frame(
  PG.ProteinGroups   = report$PG.ProteinGroups,
  UniProtAccession   = report$UniProtAccession,
  gene_name          = report$gene_name,
  quant_mat,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.csv(clean_data, file.path(output_dir, "ELA-LS_protein_quantity_clean.csv"), row.names = FALSE)
saveRDS(clean_data, file.path(output_dir, "ELA-LS_protein_quantity_clean.rds"))
message("Clean protein x sample table written: ",
        nrow(clean_data), " proteins x ", length(sample_names), " samples.")

#===========================Sample Metadata / Experimental Design=============================
experimental.design <- data.frame(label = sample_names, condition = group, stringsAsFactors = FALSE) %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(replicate = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  as.data.frame()

write.csv(experimental.design, file.path(output_dir, "sample_info.csv"), row.names = FALSE)
message("Sample groups: ", paste(names(table(group)), table(group), sep = "=", collapse = ", "))

#===========================Build DEP SummarizedExperiment=============================
# DEP::make_se needs a data frame with unique, non-NA "name"/"ID" columns
# plus the abundance columns (same convention as 01_Pre_Processing.R /
# 43_Predict_ELA_LS_Subtype.R). Fall back to UniProtAccession for any
# stray NA gene_name (defensive -- the iRT row is already dropped above).
name_raw <- as.character(report$gene_name)
name_raw[is.na(name_raw)] <- as.character(report$UniProtAccession)[is.na(name_raw)]

abu_data <- data.frame(
  quant_mat,
  name = make.unique(name_raw),
  ID   = make.unique(as.character(report$UniProtAccession)),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

abundance.columns <- match(experimental.design$label, colnames(abu_data))
if (anyNA(abundance.columns)) stop("Could not align sample columns with experimental design labels.")

message("Building SummarizedExperiment (", nrow(abu_data), " proteins x ", length(abundance.columns), " samples)...")
se_abu_data <- DEP::make_se(abu_data, abundance.columns, experimental.design)
saveRDS(se_abu_data, file.path(output_dir, "se_abu_data.rds"))

#===========================Missingness QC (before filtering)=============================
vis_miss(as.data.frame(assay(se_abu_data)), show_perc = TRUE, show_perc_col = TRUE, cluster = TRUE)
ggsave(file.path(output_dir, "missing_vis_miss_heatmap_before.pdf"), width = 11, height = 8, units = "in")

missing_per_sample <- clean_data %>%
  dplyr::select(all_of(sample_names)) %>%
  summarise_all(~ sum(is.na(.))) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Missing_Count") %>%
  mutate(Missing_Percentage = (Missing_Count / nrow(clean_data)) * 100)
write.csv(missing_per_sample, file.path(output_dir, "missing_per_sample.csv"), row.names = FALSE)

missing_per_protein <- clean_data %>%
  dplyr::select(PG.ProteinGroups, UniProtAccession, gene_name, all_of(sample_names)) %>%
  rowwise() %>%
  mutate(
    Missing_Count = sum(is.na(c_across(all_of(sample_names)))),
    Missing_Percentage = (Missing_Count / length(sample_names)) * 100
  ) %>%
  ungroup() %>%
  dplyr::select(PG.ProteinGroups, UniProtAccession, gene_name, Missing_Count, Missing_Percentage)
write.csv(missing_per_protein, file.path(output_dir, "missing_per_protein.csv"), row.names = FALSE)

#===========================Filter Proteins=============================
if (is.null(opt$threshold) || opt$threshold <= 0) {
  message("No (or zero) threshold provided, skipping protein filtering!")
  se_abu_data_filtered <- se_abu_data
} else {
  se_abu_data_filtered <- DEP::filter_proteins(se_abu_data, type = "fraction", min = opt$threshold)
  message(nrow(se_abu_data_filtered), " of ", nrow(se_abu_data),
          " proteins kept after filtering (>= ", opt$threshold * 100, "% valid values per protein).")
  vis_miss(as.data.frame(assay(se_abu_data_filtered)), show_perc = TRUE, show_perc_col = TRUE, cluster = TRUE)
  ggsave(file.path(output_dir, "missing_vis_miss_heatmap_after.pdf"), width = 11, height = 8, units = "in")
}
saveRDS(se_abu_data_filtered, file.path(output_dir, "se_abu_data_filtered.rds"))

#===========================Normalize (VSN) + Impute (MinProb)=============================
if (isTRUE(opt$normalization)) {
  message("VSN-normalizing...")
  norm <- DEP::normalize_vsn(se_abu_data_filtered)
} else {
  norm <- se_abu_data_filtered
}

pdf(file.path(output_dir, "meanSdPlot.pdf"))
meanSdPlot(norm)
dev.off()

pdf(file.path(output_dir, "Normalization_Boxplot.pdf"), height = 12)
plot_normalization(se_abu_data_filtered, norm)
dev.off()

message("MinProb imputation (q = ", opt$scalar, ")...")
norm_imp_MinProb <- DEP::impute(norm, fun = "MinProb", q = opt$scalar)

saveRDS(norm, file.path(output_dir, "norm.rds"))
saveRDS(norm_imp_MinProb, file.path(output_dir, "norm_imp_MinProb.rds"))

message("Done. Outputs written to: ", output_dir)
