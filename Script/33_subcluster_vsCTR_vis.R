#=========================Script Description=================================
# This script is used for visualize the subcluster vs control differential
# expression results (volcano plots), from 32_subcluster_vsCTR.R.
# Only the age+sex covariate, all-patients model is shown
# (norm_imp_MinProb_subcluster_age_sex_cov_all_patients); the female-only
# and male-only stratified models are not plotted here.
# Rscript 33_subcluster_vsCTR_vis.R -i Discovery/32_subcluster_vsCTR/Differential_Expression_Results.rds -o Discovery/33_subcluster_vsCTR_vis -e 9
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggrepel"))
suppressMessages(library("ggpubr"))
suppressMessages(library("stringr"))
suppressMessages(library("purrr"))
#===========================Function Definition=============================
volcano_plot <- function(df, alpha_sig, name_title, labels, case_color, ctrl_color) {
  df <- df %>%
    mutate(omic_type = case_when(
      x >= 0 & y >= -log10(alpha_sig) ~ "up",
      x <= 0 & y >= -log10(alpha_sig) ~ "down",
      TRUE ~ "ns"
    ))

  cols <- c("up" = case_color, "down" = ctrl_color, "ns" = "grey")
  ggplot(data = df, aes(x, y)) +
    geom_point(aes(colour = omic_type), alpha = 0.5, shape = 16, size = 3) +
    geom_hline(yintercept = -log10(alpha_sig), linetype = "dashed") +
    geom_text_repel(data = filter(df, name %in% labels),
                    aes(label = name),
                    force = 1,
                    nudge_x = -0.3,
                    nudge_y = 1.5,
                    direction = "both",
                    max.overlaps = 20,
                    size = 4) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_colour_manual(values = cols) +
    labs(title = name_title,
         x = "log2(fold change)",
         y = expression(-log[10] ~ "(adjusted p-value)"),
         colour = "Differential \nExpression") +
    theme_classic() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 15, hjust = 0.5)) +
    annotate("text", x = 1, y = 0.5,
             label = paste0(sum(df$omic_type == "up"), " more abundant\n", sum(df$omic_type == "down"), " less abundant"))
}

# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "32_subcluster_vsCTR/Differential_Expression_Results.rds",
              help = "32_subcluster_vsCTR/Differential_Expression_Results.rds"
  ),make_option(c("--output", "-o"),
                type = "character", default = "33_subcluster_vsCTR_vis",
                help = "output directory path."
  ),make_option(c("--seed", "-e"),
                type = "integer", default = 9,
                help = "set.seed"
  ),make_option(c("--label", "-l"),
                type = "character", default = NULL,
                help = "label proteins"
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
  stop("Please provide the cleaned SummarizedExperiment object after normalize and imputation file path!")
}else if (!file.exists(opt$input)) {
  stop("SummarizedExperiment object after normalize and imputation does not exist!")
}else{
  input = opt$input
  res = readRDS(input)
}


if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
}else{
  seed = opt$seed
  set.seed(seed)
}

if (is.null(opt$label)) {
  message("NO label proteins provided, use default setting!")
  label = opt$label
}else{
  label = opt$label
  label <- unlist(strsplit(label, ","))
}

#===========================Main Script======================================
# only the age+sex covariate, all-patients model (matches
# 32_subcluster_vsCTR.R's title = paste0("norm_imp_MinProb_subcluster_",
# covariates[i], "_", patients[j]))
name_df <- "norm_imp_MinProb_subcluster_age_sex_cov_all_patients"
data_res <- res[[name_df]]

if (is.null(data_res)){
  stop("'", name_df, "' not found in ", input,
      ". Available comparisons: ", paste(names(res), collapse = ", "))
}

# alpha_vs_ctrl and beta_vs_ctrl each have their own diff/fdr columns
comparisons <- c("alpha_vs_ctrl", "beta_vs_ctrl")

plots_FDR0.05 = plots_FDR0.1 = list()

for (i in seq_along(comparisons)) {
  comp <- comparisons[i]
  diff_col <- paste0(comp, "_diff")
  fdr_col <- paste0(comp, "_fdr")

  if (!(diff_col %in% colnames(data_res)) || !(fdr_col %in% colnames(data_res))){
    stop("Expected columns '", diff_col, "' and '", fdr_col, "' not found in '", name_df, "'.")
  }

  df <- data.frame(
    x = data_res[[diff_col]],
    y = -log10(data_res[[fdr_col]]),
    name = data_res$name
  )

  top10_label_up <- df %>%
    filter(x >= 0 & y >= -log10(0.05)) %>%
    arrange(desc(x)) %>%
    slice_head(n = 10)

  top10_label_down <- df %>%
    filter(x <= 0 & y >= -log10(0.05)) %>%
    arrange(x) %>%
    slice_head(n = 10)

  case_color <- "#D73027"
  ctrl_color <- "#4575B4"

  plots_FDR0.05[[i]] <- volcano_plot(df, 0.05,
                                     paste0("Volcano plot FDR 0.05\n", comp),
                                     c(top10_label_up$name, top10_label_down$name, label),
                                     case_color, ctrl_color)

  top10_label_up <- df %>%
    filter(x >= 0 & y >= -log10(0.1)) %>%
    arrange(desc(x)) %>%
    slice_head(n = 10)

  top10_label_down <- df %>%
    filter(x <= 0 & y >= -log10(0.1)) %>%
    arrange(x) %>%
    slice_head(n = 10)

  plots_FDR0.1[[i]] <- volcano_plot(df, 0.1,
                                    paste0("Volcano plot FDR 0.1\n", comp),
                                    c(top10_label_up$name, top10_label_down$name, label),
                                    case_color, ctrl_color)
}

#===========================Save Plots======================================
ggarrange(plotlist = plots_FDR0.05, ncol = 2, nrow = 1)
ggsave(filename = file.path(output_dir, "FDR005.pdf"), width = 12, height = 5, units = "in")

ggarrange(plotlist = plots_FDR0.1, ncol = 2, nrow = 1)
ggsave(filename = file.path(output_dir, "FDR01.pdf"), width = 12, height = 5, units = "in")
