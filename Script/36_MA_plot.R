#=========================Script Description=================================
# This script tests whether differentially expressed proteins (ALS vs CTR)
# are predominantly found among the lower-intensity proteins, using an
# MA-style plot (M = log2 fold-change, A = mean log2 intensity) for the
# Discovery and Validation cohorts.
#
# The statistic used to call a protein "significant" is configurable via
# --stat (either "fdr" or "p.val") and --cutoff (default 0.05), e.g.:
#   Rscript 36_MA_plot.R -s fdr -f 0.05 -o 36_MA_plot     # FDR < 0.05 (default)
#   Rscript 36_MA_plot.R -s p.val -f 0.05 -o 36_MA_plot   # raw p-value < 0.05
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggrepel"))
suppressMessages(library("ggpubr"))
suppressMessages(library("SummarizedExperiment"))
#===========================Function Definition=============================
# Load one cohort's mean intensity (A) and merge with its ALS vs CTR DE
# result (M = als_vs_ctrl_diff, fdr). Returns one row per protein present in
# both the intensity matrix and the DE table.
load_cohort_MA_data = function(cohort, intensity_path, de_path){
  se = readRDS(file.path(cohort, intensity_path))
  A = rowMeans(assay(se), na.rm = TRUE)
  A = data.table(name = names(A), A = as.numeric(A))

  de = fread(file.path(cohort, de_path))
  de = de[, .(name, diff = als_vs_ctrl_diff, p.val = als_vs_ctrl_p.val, fdr)]

  data = merge(A, de, by = "name")
  data$cohort = cohort
  data = as.data.frame(data)
  return(data)
}

# Add sig/direction flags, rank (1 = lowest intensity) and intensity bin
# (ntile of A, bin 1 = lowest intensity ... bin n_bins = highest intensity).
# stat_col selects which column to threshold on (e.g. "fdr" or "p.val").
annotate_MA_data = function(data, stat_col, cutoff, n_bins = 5){
  data$stat_value = data[[stat_col]]
  data = data %>%
    mutate(sig = stat_value < cutoff,
           direction = case_when(sig & diff > 0 ~ "up in ALS",
                                 sig & diff < 0 ~ "down in ALS",
                                 TRUE ~ "ns"),
           rank = rank(A, ties.method = "average"),
           bin = ntile(A, n_bins))
  return(data)}

# MA plot: x = mean log2 intensity (A), y = log2 fold-change (M), coloured
# by DE direction. The n lowest-intensity significant hits in each direction
# are labelled, to make it easy to see whether they cluster at low A.
ma_plot = function(data, cutoff, stat_label, main_title, label_n = 10,
                   case_color = "#D73027", ctrl_color = "#4575B4"){
  cols = c("up in ALS" = case_color, "down in ALS" = ctrl_color, "ns" = "grey70")

  label_df = data %>%
    filter(sig) %>%
    arrange(A) %>%
    group_by(direction) %>%
    slice_head(n = label_n) %>%
    ungroup()

  n_up = sum(data$direction == "up in ALS")
  n_down = sum(data$direction == "down in ALS")

  ggplot(data, aes(A, diff)) +
    geom_point(aes(colour = direction), alpha = 0.5, shape = 16, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_text_repel(data = label_df,
                    aes(label = name),
                    force = 1,
                    max.overlaps = 30,
                    segment.size = 0.2,
                    min.segment.length = 0,
                    size = 2.5) +
    scale_colour_manual(values = cols) +
    labs(title = main_title,
         x = "mean log2 intensity",
         y = "log2 fold-change",
         colour = paste0(stat_label, " < ", cutoff)) +
    theme_classic() +
    theme(axis.title = element_text(size = 13),
          axis.text = element_text(size = 11),
          plot.title = element_text(size = 14, hjust = 0.5)) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5,
             label = paste0(n_up, " up in ALS\n", n_down, " down in ALS"),
             size = 3)}

# MA plot with no significance cutoff: every protein is coloured by a
# continuous -log10(stat) gradient (pale yellow = high stat/not significant,
# dark blue = low stat/highly significant), so the intensity-vs-significance
# trend is visible without picking a threshold.
ma_plot_gradient = function(data, stat_col, stat_label, main_title, scale_max = NULL,
                            sig_size = 3.5, other_size = 1.8){
  data$neg_log10_stat = -log10(data[[stat_col]])
  # draw least-significant points first so the most significant (darkest)
  # points are plotted on top and stay visible
  data = data %>% arrange(neg_log10_stat)

  # colour scale fixed from 0 to floor(max(-log10(stat))), so the legend
  # always starts at 0 and ends on a whole number. Pass a shared scale_max
  # (e.g. from all cohorts combined) so multiple plots share identical
  # colour limits/breaks -- required for ggarrange(common.legend = TRUE)
  # to actually be correct for every panel, not just the first one.
  if (is.null(scale_max)) {
    scale_max = floor(max(data$neg_log10_stat, na.rm = TRUE))
  }

  ggplot(data, aes(A, diff, colour = neg_log10_stat, size = sig)) +
    geom_point(alpha = 0.7, shape = 16) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_colour_distiller(palette = "YlGnBu", direction = 1,
                           limits = c(0, scale_max), breaks = 0:scale_max,
                           oob = scales::squish) +
    scale_size_manual(values = c(`TRUE` = sig_size, `FALSE` = other_size), guide = "none") +
    labs(title = main_title,
         x = "mean log2 intensity",
         y = "log2 fold-change",
         colour = paste0("-log10(", stat_label, ")")) +
    theme_classic() +
    theme(axis.title = element_text(size = 13),
          axis.text = element_text(size = 11),
          plot.title = element_text(size = 14, hjust = 0.5))}

# Barplot: proportion of proteins that are DE-significant within each
# intensity bin (bin 1 = lowest intensity). If DE proteins are predominantly
# low-intensity, the proportion should decrease from bin 1 to bin n_bins.
intensity_bin_plot = function(data, cutoff, stat_label, main_title){
  bin_summary = data %>%
    group_by(bin) %>%
    summarise(n = n(), n_sig = sum(sig), prop_sig = n_sig / n,
             mean_A = mean(A), .groups = "drop") %>%
    mutate(bin_label = paste0("Q", bin))

  bins_sorted = sort(unique(bin_summary$bin))
  qualifier = rep("", length(bins_sorted))
  qualifier[1] = "\n(lowest)"
  qualifier[length(qualifier)] = "\n(highest)"
  x_labels = paste0("Q", bins_sorted, qualifier)

  ggplot(bin_summary, aes(x = factor(bin), y = prop_sig)) +
    geom_col(fill = "#4575B4", alpha = 0.8, width = 0.65) +
    geom_text(aes(label = paste0(n_sig, "/", n)), vjust = -0.4, size = 3) +
    scale_x_discrete(labels = x_labels) +
    labs(title = main_title,
         x = "Intensity bin (equal-sized, ranked low → high)",
         y = paste0("Proportion ", stat_label, " < ", cutoff)) +
    theme_classic() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          plot.title = element_text(size = 14, hjust = 0.5))}

# Formal statistical tests for "DE proteins skew toward lower intensity":
#  - Cochran-Armitage trend test on proportion significant across ordered
#    intensity bins (stats::prop.trend.test), score = bin index.
#  - One-sided Wilcoxon rank-sum test: is A lower in the significant group
#    than in the non-significant group?
run_intensity_tests = function(data, n_bins){
  bin_summary = data %>%
    group_by(bin) %>%
    summarise(n = n(), n_sig = sum(sig), .groups = "drop") %>%
    arrange(bin)

  trend_test = stats::prop.trend.test(x = bin_summary$n_sig, n = bin_summary$n,
                                      score = bin_summary$bin)

  wilcox_test = wilcox.test(data$A[data$sig], data$A[!data$sig], alternative = "less")

  list(bin_summary = bin_summary, trend_test = trend_test, wilcox_test = wilcox_test)}
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--intensity", "-i"),
              type = "character", default = "02_Missing_Inspection/norm_imp_MinProb.rds",
              help = "VSN-normalized, imputed SummarizedExperiment, relative to each cohort folder."
  ),make_option(c("--de", "-d"),
                type = "character", default = "03_Differential_expression_analysis/norm_imp_MinProb_age_sex_cov_all_patients.csv",
                help = "ALS vs CTR DE result CSV (age+sex-adjusted, all patients), relative to each cohort folder."
  ),make_option(c("--stat", "-s"),
                type = "character", default = "fdr",
                help = "Which statistic column to threshold on: 'fdr' or 'p.val'. Default 'fdr'."
  ),make_option(c("--cutoff", "-f"),
                type = "double", default = 0.05,
                help = "Cutoff applied to --stat to call a protein differentially expressed. Default 0.05."
  ),make_option(c("--n_bins", "-n"),
                type = "integer", default = 5,
                help = "Number of equal-sized intensity bins (ranked low to high). Default 5 (quintiles)."
  ),make_option(c("--output", "-o"),
                type = "character", default = "36_MA_plot",
                help = "output directory path."
  ),make_option(c("--seed", "-e"),
                type = "integer", default = 9,
                help = "set.seed"
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

if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
}else{
  seed = opt$seed
  set.seed(seed)
}

stat_col = opt$stat
if (!stat_col %in% c("fdr", "p.val")) {
  stop("--stat must be either 'fdr' or 'p.val' (got '", stat_col, "')")
}
stat_label = if (stat_col == "fdr") "FDR" else "p-value"
cutoff = opt$cutoff
n_bins = opt$n_bins
cohorts = c("Discovery", "Validation")

################################################################################################
#LOAD + ANNOTATE DATA PER COHORT
data_all = list()
for (cohort in cohorts){
  d = load_cohort_MA_data(cohort, opt$intensity, opt$de)
  d = annotate_MA_data(d, stat_col = stat_col, cutoff = cutoff, n_bins = n_bins)
  data_all[[cohort]] = d
  write.csv(d, paste0(output_dir, "/", cohort, "_MA_data.csv"), row.names = FALSE)
}

################################################################################################
#MA PLOTS (per cohort + combined)
ma_plots = list()
for (cohort in cohorts){
  ma_plots[[cohort]] = ma_plot(data_all[[cohort]], cutoff, stat_label,
                               main_title = paste0(cohort, ": ALS vs CTR"))
  ggsave(paste0(output_dir, "/", cohort, "_MA_plot.pdf"), ma_plots[[cohort]],
        width = 6, height = 4)
  ggsave(paste0(output_dir, "/", cohort, "_MA_plot.svg"), ma_plots[[cohort]],
        width = 6, height = 4)
}

ggarrange(plotlist = ma_plots, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
ggsave(paste0(output_dir, "/Combined_MA_plot.pdf"), width = 12, height = 4)
ggsave(paste0(output_dir, "/Combined_MA_plot.svg"), width = 12, height = 4)

################################################################################################
#MA PLOTS - CONTINUOUS GRADIENT, NO CUTOFF (per cohort + combined)
# shared colour scale across cohorts, so common.legend = TRUE below is
# actually valid (both panels use identical limits/breaks)
neg_log10_all = -log10(unlist(lapply(data_all, function(d) d[[stat_col]])))
shared_scale_max = floor(max(neg_log10_all, na.rm = TRUE))

ma_plots_gradient = list()
for (cohort in cohorts){
  ma_plots_gradient[[cohort]] = ma_plot_gradient(data_all[[cohort]], stat_col, stat_label,
                                                 main_title = paste0(cohort, ": ALS vs CTR"),
                                                 scale_max = shared_scale_max)
  ggsave(paste0(output_dir, "/", cohort, "_MA_plot_gradient.pdf"), ma_plots_gradient[[cohort]],
        width = 6, height = 4)
  ggsave(paste0(output_dir, "/", cohort, "_MA_plot_gradient.svg"), ma_plots_gradient[[cohort]],
        width = 6, height = 4)
}

ggarrange(plotlist = ma_plots_gradient, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
ggsave(paste0(output_dir, "/Combined_MA_plot_gradient.pdf"), width = 12, height = 4)
ggsave(paste0(output_dir, "/Combined_MA_plot_gradient.svg"), width = 12, height = 4)
