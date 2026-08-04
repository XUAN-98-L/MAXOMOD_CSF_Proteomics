#=========================Script Description=================================
# This script is used for visualizing cross-cohort agreement of the ALS
# subtype vs control differential expression results (alpha vs ctrl, beta
# vs ctrl), the same way 12_Scatterplot_FDR.R does for als vs ctrl: for
# each protein, signed -log10(FDR) [and signed -log10(p-value)] in Discovery
# (x) is plotted against Validation (y), sign taken from the direction of
# change (up/down) in that cohort.
#
# Input: 32_subcluster_vsCTR/Differential_Expression_Results.rds (alpha vs
# ctrl and beta vs ctrl, age+sex-adjusted, all patients), produced by
# 32_subcluster_vsCTR.R. alpha and beta are handled as two SEPARATE
# comparisons, each producing its own FDR and p-value scatter plot.
#
# Rscript 37_Scatterplot_subtype_Ctrl.R -i 32_subcluster_vsCTR/Differential_Expression_Results.rds -o 37_Scatterplot_subtype_Ctrl
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggrepel"))
#===========================Function Definition=============================
# Discovery (x) vs Validation (y) scatter of a signed statistic, same
# function used by 12_Scatterplot_FDR.R / 12_Scatterplot_FDR_subclusters.R.
scatterplot_FDR_male_female = function(data,
                                       cut_off = -log10(0.05),
                                       q = 0.95,
                                       main_title,
                                       max.overlaps = 10,
                                       lab_x = "signed -log10(FDR) for Discovery",
                                       lab_y = "signed -log10(FDR) for Validation",
                                       labels_T_F = T,
                                       annotate_YN = T,
                                       text_y = "significant in Validation",
                                       text_x = "significant in Discovery"){
  data$omic_type = rep("ns", nrow(data))
  data$omic_type[abs(data$y) >= cut_off] = text_y
  data$omic_type[abs(data$x) >= cut_off] = text_x
  data$omic_type[(abs(data$x) >= cut_off) & (abs(data$y) >= cut_off)] = "significant in both"
  cols <- c("x" = "salmon", "y" = "#26b3ff", "ns" = "grey", "significant in both" = "mediumpurple1")
  attributes(cols)$names[1] = text_x
  attributes(cols)$names[2] = text_y

  quantile_y = quantile(abs(data$y),na.rm = T, probs = q)
  quantile_x = quantile(abs(data$x),na.rm = T, probs = q)

  plot = ggplot(data, aes(x,y)) +
    geom_point(aes(colour = omic_type),
               alpha = 0.5,
               shape = 16,
               size = 2) +
    geom_point(data = filter(data, abs(y) >= cut_off | abs(x) >= cut_off),
               aes(colour = omic_type),
               alpha = 0.5,
               shape = 16,
               size = 3) +
    geom_smooth(method = "lm", color = "#2C3E50", se = TRUE) +
    geom_hline(yintercept = cut_off, linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -cut_off, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = cut_off, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = -cut_off, linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey80") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey80") +
    geom_text_repel(data = filter(data, abs(y) >= quantile_y | abs(x) >= quantile_x),
                    aes(label = name),
                    force = 1,
                    hjust = 1,
                    max.overlaps = max.overlaps,
                    segment.size = 0.2,
                    min.segment.length = 0,
                    size = 2)  +
    scale_colour_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(title = main_title,
         x = lab_x,
         y = lab_y,
         colour = "Differential \nExpression") +
    theme_classic() + # Select theme with a white background
    theme(axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 15, hjust = 0.5),
          text = element_text(size = 14)) +
    if(annotate_YN){
      annotate("text", x = -0.75, y = 1.75, label =
                 paste0(sum(data$omic_type==text_y), " ", text_y,"\n",
                        sum(data$omic_type==text_x), " ", text_x,"\n",
                        sum(data$omic_type=="significant in both"), " significant in both"), size = 8/.pt)
    }

  return(plot)}

# Build the signed -log10(stat_col) data.table (x = Discovery, y =
# Validation) for one subtype-vs-ctrl contrast (e.g. "alpha_vs_ctrl"), sign
# taken from that contrast's own "_diff" column in each cohort.
build_signed_data = function(res_Discovery_inter, res_Validation_inter, contrast, stat_col){
  diff_col = paste0(contrast, "_diff")
  value_col = paste0(contrast, "_", stat_col)

  data = as.data.table(cbind(res_Discovery_inter$name,
                             res_Discovery_inter[[value_col]],
                             res_Validation_inter[[value_col]]))
  colnames(data) = c("name", "x", "y")
  data$x = as.numeric(data$x)
  data$y = as.numeric(data$y)
  data$x = -log10(data$x)
  data$y = -log10(data$y)

  data$x[res_Discovery_inter[[diff_col]] < 0] = -data$x[res_Discovery_inter[[diff_col]] < 0]
  data$y[res_Validation_inter[[diff_col]] < 0] = -data$y[res_Validation_inter[[diff_col]] < 0]

  return(data)}

# Run both the FDR-based and p-value-based Discovery-vs-Validation scatter
# plot for one subtype-vs-ctrl contrast, writing plots/CSVs/correlation
# tests to output_dir, filenames prefixed with the contrast name.
run_contrast_scatter = function(contrast, contrast_label, res_Discovery_inter, res_Validation_inter, output_dir){

  #### FDR ####
  data = build_signed_data(res_Discovery_inter, res_Validation_inter, contrast, "fdr")

  p = scatterplot_FDR_male_female(data, q = 0.95, main_title = paste0(contrast_label, " (FDR)"), max.overlaps = Inf)
  ggsave(paste0(output_dir, "/", contrast, "_FDR.pdf"), p, width = 6, height = 4, dpi = 300)

  sink(paste0(output_dir, "/", contrast, "_FDR_cor.txt"))
  print(cor.test(data$x, data$y, method = "pearson"))
  print(cor.test(data$x, data$y, method = "spearman"))
  sink()

  data$signif_dis <- 10^(-abs(data$x)) < 0.05
  data$signif_vali <- 10^(-abs(data$y)) < 0.05
  data$signif_both <- data$signif_dis & data$signif_vali
  write.csv(data, paste0(output_dir, "/", contrast, "_FDR.csv"), row.names = F)

  #### p-value ####
  data = build_signed_data(res_Discovery_inter, res_Validation_inter, contrast, "p.val")

  p = scatterplot_FDR_male_female(data, q = 0.95, main_title = paste0(contrast_label, " (p-value)"), max.overlaps = Inf,
                                  lab_x = "signed -log10(pvalue) for Discovery",
                                  lab_y = "signed -log10(pvalue) for Validation")
  ggsave(paste0(output_dir, "/", contrast, "_pvalue.pdf"), p, width = 6, height = 4, dpi = 300)

  sink(paste0(output_dir, "/", contrast, "_pvalue_cor.txt"))
  print(cor.test(data$x, data$y, method = "pearson"))
  print(cor.test(data$x, data$y, method = "spearman"))
  sink()

  data$signif_dis <- 10^(-abs(data$x)) < 0.05
  data$signif_vali <- 10^(-abs(data$y)) < 0.05
  data$signif_both <- data$signif_dis & data$signif_vali
  write.csv(data, paste0(output_dir, "/", contrast, "_pvalue.csv"), row.names = F)

  invisible(NULL)}
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "32_subcluster_vsCTR/Differential_Expression_Results.rds",
              help = "32_subcluster_vsCTR/Differential_Expression_Results.rds"
  ),make_option(c("--comparison", "-c"),
                type = "character", default = "norm_imp_MinProb_subcluster_age_sex_cov_all_patients",
                help = "which element of Differential_Expression_Results.rds to use (age+sex-adjusted, all patients, by default)."
  ),make_option(c("--output", "-o"),
                type = "character", default = "37_Scatterplot_subtype_Ctrl",
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

if (is.null(opt$input)) {
  stop("Please provide the Differential_Expression_Results.rds file path!")
}else{
  res_Discovery = readRDS(paste0("Discovery/",opt$input))
  res_Validation = readRDS(paste0("Validation/",opt$input))
}

if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
}else{
  seed = opt$seed
  set.seed(seed)
}

comparison = opt$comparison

################################################################################################
#INTERSECT PROTEINS PRESENT IN BOTH COHORTS, SAME ORDER
inter = intersect(res_Discovery[[comparison]]$name,
                  res_Validation[[comparison]]$name)

res_Discovery_inter = res_Discovery[[comparison]][which(res_Discovery[[comparison]]$name %in% inter),]
res_Discovery_inter = res_Discovery_inter[match(inter, res_Discovery_inter$name),]

res_Validation_inter = res_Validation[[comparison]][which(res_Validation[[comparison]]$name %in% inter),]
res_Validation_inter = res_Validation_inter[match(inter, res_Validation_inter$name),]

################################################################################################
#RUN alpha_vs_ctrl AND beta_vs_ctrl SEPARATELY
run_contrast_scatter("alpha_vs_ctrl", "Alpha vs Control", res_Discovery_inter, res_Validation_inter, output_dir)
run_contrast_scatter("beta_vs_ctrl", "Beta vs Control", res_Discovery_inter, res_Validation_inter, output_dir)
