#=========================Script Description=================================
# Summary of protein detection before and after completeness filtering
# (DEP::filter_proteins, >50% missing values across all samples removed),
# in the Discovery and Validation cohorts.
#
# For each cohort, the number of proteins quantified (non-missing) is
# counted per sample, before filtering (01_Pre_Processing/se_abu_data.rds)
# and after filtering (01_Pre_Processing/se_abu_data_filtered.rds).
# Boxplots show the median and interquartile range (IQR). In addition, one
# per-sample bar plot is produced for each cohort x stage combination
# (disease-coloured, samples ordered ascending within group).
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
#===========================Function Definition=============================
# Per-sample count of quantified (non-missing) proteins for one SE object.
# "sample" is taken from colData(se)$label (the patient/sample ID, e.g.
# "CSF83"), not from colnames(se) (which is "condition_replicate", e.g.
# "als_1"). Also carries colData(se)$condition ("als"/"ctrl") through,
# needed for the per-sample bar plots below.
count_per_sample = function(se, cohort, stage){
  n_quantified = colSums(!is.na(assay(se)))
  data.frame(cohort = cohort,
            stage = stage,
            sample = as.character(colData(se)$label),
            condition = as.character(colData(se)$condition),
            n_proteins = as.numeric(n_quantified),
            stringsAsFactors = FALSE)}

# Boxplot of per-sample protein counts, grouped by cohort and coloured by
# filtering stage (median + IQR box, default Tukey whiskers/outliers).
plot_protein_numbers = function(data, stage_colors){
  ggplot(data, aes(x = cohort, y = n_proteins, fill = stage)) +
    geom_boxplot(outlier.size = 1, width = 0.6) +
    scale_fill_manual(values = stage_colors) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(title = NULL,
         x = NULL,
         y = "Proteins quantified per sample",
         fill = NULL) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_blank(),
         panel.grid.minor = element_blank())}

# Per-sample bar plot of quantified-protein counts for a single cohort/stage
# combination, bars coloured by disease group ("als"/"ctrl") and ordered
# ascending within each group (als block first, then ctrl block).
condition_colors = c(als = "#F8766D", ctrl = "#00BFC4")

plot_protein_numbers_per_sample = function(data, main_title){
  data$condition = factor(data$condition, levels = c("als", "ctrl"))
  data = data[order(data$condition, data$n_proteins), ]
  data$sample = factor(data$sample, levels = data$sample)
  n_als = sum(data$condition == "als")
  n_total = nrow(data)
  x_als_mid = (1 + n_als) / 2
  x_ctrl_mid = (n_als + 1 + n_total) / 2

  ggplot(data, aes(x = sample, y = n_proteins, fill = condition)) +
    geom_col(width = 0.8) +
    geom_vline(xintercept = n_als + 0.5, linetype = "dashed", colour = "grey30") +
    annotate("text", x = x_als_mid, y = Inf, label = "ALS", vjust = 1.4, fontface = "bold", size = 4) +
    annotate("text", x = x_ctrl_mid, y = Inf, label = "Ctrl", vjust = 1.4, fontface = "bold", size = 4) +
    scale_fill_manual(values = condition_colors, name = "disease") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    labs(title = main_title, x = NULL, y = "Number of proteins") +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0, size = 14),
         panel.grid.major.x = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text.x = element_text(size = 3, angle = 90, hjust = 1, vjust = 0.5))}

# Per-sample bar plot for a single cohort, before- and after-filtering counts
# shown as two dodged bars per sample (coloured by stage). Samples are
# ordered ascending by their after-filtering count, within condition group.
plot_protein_numbers_per_sample_stage = function(data, stage_colors, main_title){
  data$condition = factor(data$condition, levels = c("als", "ctrl"))
  order_df = data[data$stage == "After filtering", ]
  order_df = order_df[order(order_df$condition, order_df$n_proteins), ]
  data$sample = factor(data$sample, levels = order_df$sample)
  n_als = sum(order_df$condition == "als")
  n_total = nrow(order_df)
  x_als_mid = (1 + n_als) / 2
  x_ctrl_mid = (n_als + 1 + n_total) / 2

  ggplot(data, aes(x = sample, y = n_proteins, fill = stage)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_vline(xintercept = n_als + 0.5, linetype = "dashed", colour = "grey30") +
    annotate("text", x = x_als_mid, y = Inf, label = "ALS", vjust = 1.4, fontface = "bold", size = 4) +
    annotate("text", x = x_ctrl_mid, y = Inf, label = "Ctrl", vjust = 1.4, fontface = "bold", size = 4) +
    scale_fill_manual(values = stage_colors) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    labs(title = main_title, x = NULL, y = "Number of proteins", fill = NULL) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0, size = 14),
         panel.grid.major.x = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text.x = element_text(size = 3, angle = 90, hjust = 1, vjust = 0.5))}
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--before", "-b"),
              type = "character", default = "01_Pre_Processing/se_abu_data.rds",
              help = "unfiltered SummarizedExperiment object, before completeness filter."
  ),make_option(c("--after", "-a"),
                type = "character", default = "01_Pre_Processing/se_abu_data_filtered.rds",
                help = "filtered SummarizedExperiment object, after completeness filter."
  ),make_option(c("--output", "-o"),
                type = "character", default = "28_Protein_numbers",
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

if (is.null(opt$before) | is.null(opt$after)) {
  stop("Please provide both the before- and after-completeness-filter SummarizedExperiment object file paths!")
}else{
  res_Discovery_before = readRDS(paste0("Discovery/",opt$before))
  res_Validation_before = readRDS(paste0("Validation/",opt$before))
  res_Discovery_after = readRDS(paste0("Discovery/",opt$after))
  res_Validation_after = readRDS(paste0("Validation/",opt$after))
}

if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
}else{
  seed = opt$seed
  set.seed(seed)
}

################################################################################################
#PER-SAMPLE PROTEIN COUNTS
data = rbind(
  count_per_sample(res_Discovery_before, "Discovery", "Before filtering"),
  count_per_sample(res_Discovery_after, "Discovery", "After filtering"),
  count_per_sample(res_Validation_before, "Validation", "Before filtering"),
  count_per_sample(res_Validation_after, "Validation", "After filtering")
)
data$cohort = factor(data$cohort, levels = c("Discovery", "Validation"))
data$stage = factor(data$stage, levels = c("Before filtering", "After filtering"))

write.csv(data, paste0(output_dir, "/protein_numbers_per_sample.csv"), row.names = FALSE)

summary_table = as.data.table(data)[, .(n_samples = .N,
                                        median = median(n_proteins),
                                        Q1 = quantile(n_proteins, 0.25),
                                        Q3 = quantile(n_proteins, 0.75),
                                        min = min(n_proteins),
                                        max = max(n_proteins)),
                                    by = .(cohort, stage)]
write.csv(summary_table, paste0(output_dir, "/protein_numbers_summary.csv"), row.names = FALSE)

################################################################################################
#PLOT
stage_colors = c("Before filtering" = "#F8766D", "After filtering" = "#00BFC4")

p = plot_protein_numbers(data, stage_colors)

ggsave(paste0(output_dir, "/Protein_numbers_before_after_filter.pdf"), p, width = 5, height = 4, dpi = 300)

################################################################################################
#PER-SAMPLE BAR PLOTS, ONE PER COHORT x STAGE, COLOURED BY DISEASE GROUP
per_sample_specs = list(
  list(cohort = "Discovery",  stage = "before", label = "Discovery cohort, before filtering"),
  list(cohort = "Discovery",  stage = "after",  label = "Discovery cohort, after filtering"),
  list(cohort = "Validation", stage = "before", label = "Validation cohort, before filtering"),
  list(cohort = "Validation", stage = "after",  label = "Validation cohort, after filtering")
)
stage_label_map = c(before = "Before filtering", after = "After filtering")

for(spec in per_sample_specs){
  d = data[data$cohort == spec$cohort & data$stage == stage_label_map[[spec$stage]], ]
  p_sample = plot_protein_numbers_per_sample(d, paste0("Proteins per sample (", spec$label, ")"))
  ggsave(paste0(output_dir, "/proteins_per_sample_", spec$cohort, "_", spec$stage, ".pdf"),
        p_sample, width = 9, height = 4.5, dpi = 300)}

################################################################################################
#PER-SAMPLE BAR PLOTS, ONE PER COHORT, BEFORE VS AFTER FILTERING AS TWO DODGED BARS
for(coh in c("Discovery", "Validation")){
  d = data[data$cohort == coh, ]
  p_stage = plot_protein_numbers_per_sample_stage(d, stage_colors,
                                                  paste0("Proteins per sample (", coh, " cohort)"))
  ggsave(paste0(output_dir, "/proteins_per_sample_before_after_", coh, ".pdf"),
        p_stage, width = 9, height = 4.5, dpi = 300)}
