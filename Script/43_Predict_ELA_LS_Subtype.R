#=========================Script Description=================================
# Applies the alpha/beta ALS molecular-subtype classifiers trained in
# 15_ML_multi_models.R to the new ELA-LS/ELA-SS cohort (Script/42_Convert_Protein.R
# output), using the 5-protein panel fixed by the user from
# 15_ML_multi_models_500_04_noAge_3sd (PARK7, PTPRS, ATRN, CNTN1, PCSK1N).
#
# Why this needs its own pipeline and can't just call predict() on the raw
# report: the models were trained on Discovery/Validation data that went
# through DEP normalize_vsn() + impute(fun = "MinProb") (see
# 02_Missing_Inspection.R), then 2^-back-transformed and min-max scaled
# per-cohort (scale_manual(), see 15_ML_multi_models.R). The ELA-LS report is
# raw, un-normalised, un-imputed Spectronaut PG.Quantity. This script
# reproduces the same steps on the ELA-LS cohort so the input the models see
# is on the same scale they were trained on:
#   1. Read the annotated report from 42_Convert_Protein.R.
#   2. Keep only ELA-LS_*/ELA-SS_* samples (ALS patients). Control samples are
#      dropped by default -- alpha/beta subtype is only defined for ALS
#      patients (--include_controls keeps them, but their predictions are not
#      biologically meaningful).
#   3. Build a DEP SummarizedExperiment from all 1147 quantified proteins
#      (needed so normalize_vsn()/impute() have enough data to fit on -- the
#      5-protein panel is only selected afterwards), normalize + impute
#      exactly like 02_Missing_Inspection.R.
#   4. 2^-back-transform, subset to the fixed 5-protein panel, min-max scale
#      using this cohort's own range (scale_manual(), same convention
#      15_ML_multi_models.R used for the Validation cohort).
#   5. Load ML_results.rds (all bootstrap-trained caret models for that
#      panel) and predict(..., type = "prob") for every bootstrap x model,
#      then report the mean probability + majority-vote class per sample.
#
# Rscript Script/43_Predict_ELA_LS_Subtype.R
#===========================Loading Packages=============================
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(DEP))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(ggplot2))

#===========================Function Definition=============================
# same min-max scaling convention as 15_ML_multi_models.R::scale_manual(),
# with na.rm added defensively (no NA should remain after MinProb imputation)
scale_manual <- function(df) {
  as.data.frame(apply(df, 2, function(x) (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE))))
}

# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("-r", "--report"), type = "character",
              default = "42_Convert_Protein/ELA-LS_Report_with_gene_names.rds",
              help = "Annotated ELA-LS report from 42_Convert_Protein.R (.rds preferred, .csv also works). [default: %default]"
  ), make_option(c("-m", "--ml_dir"), type = "character",
              default = "15_ML_multi_models_500_04_noAge_3sd",
              help = "Folder with ML_results.rds + selected_features.txt (feature panel and models must come from the same run). [default: %default]"
  ), make_option(c("-o", "--output"), type = "character",
              default = "43_Predict_ELA_LS_Subtype",
              help = "Output directory. [default: %default]"
  ), make_option(c("-q", "--minprob_q"), type = "numeric", default = 0.01,
              help = "MinProb imputation quantile, same default as 02_Missing_Inspection.R. [default: %default]"
  ), make_option(c("--include_controls"), action = "store_true", default = FALSE,
              help = "Also score Control (Ctl_*) samples. Off by default: alpha/beta subtype is only defined for ALS patients."
  ), make_option(c("-e", "--seed"), type = "integer", default = 9,
              help = "set.seed [default: %default]"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
set.seed(opt$seed)
#============================================================================
if (!file.exists(opt$report)) stop("Report not found: ", opt$report)
ml_results_path <- file.path(opt$ml_dir, "ML_results.rds")
features_path    <- file.path(opt$ml_dir, "selected_features.txt")
if (!file.exists(ml_results_path)) stop("ML_results.rds not found in: ", opt$ml_dir)
if (!file.exists(features_path))   stop("selected_features.txt not found in: ", opt$ml_dir)

output_dir <- opt$output
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

#===========================Load Trained Models + Fixed Feature Panel=============================
message("Loading trained models from: ", opt$ml_dir)
ml_results <- readRDS(ml_results_path)
selected_final <- scan(features_path, what = character(), quiet = TRUE)
selected_final <- gsub('"', '', selected_final)
message("Fixed protein panel (", length(selected_final), "): ", paste(selected_final, collapse = ", "))

#===========================Read & Reshape the ELA-LS Report=============================
message("Reading ELA-LS report: ", opt$report)
report <- if (grepl("\\.rds$", opt$report)) {
  readRDS(opt$report)
} else {
  read.csv(opt$report, check.names = FALSE)
}
report <- as.data.frame(report)

if (!all(c("gene_name", "UniProtAccession") %in% colnames(report))) {
  stop("Report is missing 'gene_name'/'UniProtAccession' columns -- run Script/42_Convert_Protein.R first.")
}

quant_cols <- grep("\\.htrms\\.PG\\.Quantity$", colnames(report), value = TRUE)
if (length(quant_cols) == 0) stop("No '...PG.Quantity' columns found in the report.")

sample_names <- gsub("^\\[[0-9]+\\]\\s*", "", quant_cols)
sample_names <- gsub("\\.htrms\\.PG\\.Quantity$", "", sample_names)
group <- gsub("_[0-9]+$", "", sample_names)      # "Ctl" / "ELA-LS" / "ELA-SS"
group[group == "Ctl"] <- "Control"

keep <- if (isTRUE(opt$include_controls)) rep(TRUE, length(sample_names)) else group != "Control"
if (!any(keep)) stop("No ALS (ELA-LS/ELA-SS) samples left after filtering.")

message(sum(keep), " of ", length(sample_names), " samples kept: ",
        paste(names(table(group[keep])), table(group[keep]), sep = "=", collapse = ", "),
        if (!isTRUE(opt$include_controls)) ". Control samples excluded (alpha/beta is only defined for ALS patients; use --include_controls to override)." else "")

use_cols  <- quant_cols[keep]
use_names <- sample_names[keep]
use_group <- group[keep]

mat <- report[, use_cols, drop = FALSE]
colnames(mat) <- use_names
mat <- as.data.frame(lapply(mat, as.numeric))   # "NaN" text (if read from csv) -> proper NA/NaN doubles

mat_df <- data.frame(name = report$gene_name, ID = make.unique(report$UniProtAccession),
                     mat, check.names = FALSE, stringsAsFactors = FALSE)

experimental.design <- data.frame(label = use_names, condition = use_group, stringsAsFactors = FALSE) %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(replicate = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  as.data.frame()

abundance.columns <- match(experimental.design$label, colnames(mat_df))
if (anyNA(abundance.columns)) stop("Could not align sample columns with experimental design labels.")

#===========================Normalize + Impute (same steps as 02_Missing_Inspection.R)=============================
# Done on the FULL protein set (1147 proteins) so normalize_vsn()/impute()
# have enough data to fit on; the panel is subset afterwards.
message("Building SummarizedExperiment (", nrow(mat_df), " proteins x ", ncol(mat), " samples), VSN-normalizing, MinProb-imputing...")
se   <- DEP::make_se(mat_df, abundance.columns, experimental.design)
norm <- DEP::normalize_vsn(se)
imp  <- DEP::impute(norm, fun = "MinProb", q = opt$minprob_q)

saveRDS(imp, file.path(output_dir, "ELA_LS_norm_imp_MinProb.rds"))

#===========================Subset to the Fixed Panel + Min-Max Scale=============================
missing_panel <- setdiff(selected_final, rowData(imp)$name)
if (length(missing_panel) > 0) {
  stop("These panel proteins are missing from the ELA-LS report: ", paste(missing_panel, collapse = ", "))
}

expr_lin <- as.data.frame(2^assay(imp))     # back to linear scale, same convention as 15_ML_multi_models.R
rownames(expr_lin) <- rowData(imp)$name
expr_lin <- expr_lin[selected_final, , drop = FALSE]

X_new <- as.data.frame(t(expr_lin))          # samples x proteins
X_new <- X_new[, selected_final, drop = FALSE]

X_new_scaled <- scale_manual(X_new)
rownames(X_new_scaled) <- rownames(X_new)

write.csv(data.frame(sample = rownames(X_new), group = use_group, X_new, check.names = FALSE),
         file.path(output_dir, "ELA_LS_panel_expression_raw.csv"), row.names = FALSE)
write.csv(data.frame(sample = rownames(X_new_scaled), group = use_group, X_new_scaled, check.names = FALSE),
         file.path(output_dir, "ELA_LS_panel_expression_scaled.csv"), row.names = FALSE)

#===========================Predict with Every Bootstrap Model=============================
# Same "positive class" convention as 15_ML_multi_models.R::get_probs(): the
# training labels are factor levels c("alpha","beta"), so levels(y)[2] ==
# "beta" is the probability that gets reported/compared to 0.5 throughout
# that script's ROC/AUC calculations. We mirror that here.
models_available <- names(ml_results$all_models[[1]])
message("Bootstrap models available: ", paste(models_available, collapse = ", "))

predict_one_model <- function(model_key) {
  probs <- sapply(ml_results$all_models, function(models) {
    fit <- models[[model_key]]
    p <- predict(fit, X_new_scaled, type = "prob")
    pos <- "beta"
    if (!pos %in% colnames(p)) pos <- colnames(p)[2]
    p[[pos]]
  })
  if (is.null(dim(probs))) probs <- matrix(probs, ncol = 1)  # single bootstrap edge case
  rownames(probs) <- rownames(X_new_scaled)
  probs
}

pred_long <- dplyr::bind_rows(lapply(models_available, function(mk) {
  probs <- predict_one_model(mk)
  data.frame(
    sample = rownames(probs),
    group  = use_group,
    model  = mk,
    mean_prob_beta = rowMeans(probs, na.rm = TRUE),
    sd_prob_beta   = apply(probs, 1, sd, na.rm = TRUE),
    predicted_class = ifelse(rowMeans(probs, na.rm = TRUE) > 0.5, "beta", "alpha"),
    stringsAsFactors = FALSE
  )
}))

write.csv(pred_long, file.path(output_dir, "ELA_LS_subtype_predictions_long.csv"), row.names = FALSE)

# wide format: one row per sample, one column per model + a consensus (majority vote across models)
pred_wide <- pred_long %>%
  dplyr::select(sample, group, model, predicted_class) %>%
  tidyr::pivot_wider(names_from = model, values_from = predicted_class)

pred_wide <- as.data.frame(pred_wide, stringsAsFactors = FALSE)
model_cols <- setdiff(colnames(pred_wide), c("sample", "group"))
vote_mat <- as.matrix(pred_wide[, model_cols, drop = FALSE])
pred_wide$n_alpha_votes <- rowSums(vote_mat == "alpha")
pred_wide$n_beta_votes  <- rowSums(vote_mat == "beta")
pred_wide$consensus_class <- ifelse(pred_wide$n_beta_votes > pred_wide$n_alpha_votes, "beta",
                                    ifelse(pred_wide$n_alpha_votes > pred_wide$n_beta_votes, "alpha", "tie"))

write.csv(pred_wide, file.path(output_dir, "ELA_LS_subtype_predictions_wide.csv"), row.names = FALSE)

message("Consensus class counts:")
print(table(pred_wide$group, pred_wide$consensus_class))

#===========================Sanity-check Plot=============================
# predicted P(beta) by model, split by ELA-LS (long survival) vs ELA-SS
# (short survival) -- a quick look at whether predicted subtype tracks with
# the survival grouping already present in this cohort.
p <- ggplot(pred_long, aes(x = group, y = mean_prob_beta, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  facet_wrap(~model) +
  labs(title = "Predicted P(beta) in the ELA-LS/ELA-SS cohort",
       subtitle = paste0("Panel: ", paste(selected_final, collapse = ", ")),
       x = NULL, y = "Mean predicted P(beta) across bootstraps") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")
ggsave(file.path(output_dir, "predicted_prob_beta_by_group.pdf"), p, width = 9, height = 6)
ggsave(file.path(output_dir, "predicted_prob_beta_by_group.png"), p, width = 9, height = 6, dpi = 300)

message("Done. Predictions written to: ", file.path(output_dir, "ELA_LS_subtype_predictions_wide.csv"))
