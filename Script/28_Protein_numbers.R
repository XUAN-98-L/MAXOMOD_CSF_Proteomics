#=========================Script Description=================================
# Summary of protein detection before and after completeness filtering
# (DEP::filter_proteins, >50% missing values across all samples removed),
# in the Discovery and Validation cohorts.
#
# For each cohort, the number of proteins quantified (non-missing) is
# counted per sample, before filtering (01_Pre_Processing/se_abu_data.rds)
# and after filtering (01_Pre_Processing/se_abu_data_filtered.rds).
# Boxplots show the median and interquartile range (IQR).
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
#===========================Function Definition=============================
# Per-sample count of quantified (non-missing) proteins for one SE object.
count_per_sample = function(se, cohort, stage){
  n_quantified = colSums(!is.na(assay(se)))
  data.frame(cohort = cohort,
            stage = stage,
            sample = names(n_quantified),
            n_proteins = as.numeric(n_quantified),
            stringsAsFactors = FALSE)}

# Boxplot of per-sample protein counts, grouped by cohort and coloured by
# filtering stage (median + IQR box, default Tukey whiskers/outliers).
plot_protein_numbers = function(data, stage_colors,
                                main_title = "Per-sample protein counts before vs after filter"){
  ggplot(data, aes(x = cohort, y = n_proteins, fill = stage)) +
    geom_boxplot(outlier.size = 1, width = 0.6) +
    scale_fill_manual(values = stage_colors) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(title = main_title,
         x = NULL,
         y = "Proteins quantified per sample",
         fill = NULL) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
         panel.grid.minor = element_blank())}
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

ggsave(paste0(output_dir, "/Protein_numbers_before_after_filter.pdf"), p, width = 6, height = 5, dpi = 300)
