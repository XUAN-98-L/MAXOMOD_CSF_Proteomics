#=========================Script Description=================================
# This script visualizes the directional concordance of differential protein
# expression between the discovery and validation cohorts. It produces a
# scatter plot of log2 fold changes (logFC) for the ALS vs Control comparison
# in the discovery (x-axis) and validation (y-axis) cohorts. Positive values
# indicate proteins upregulated in ALS, negative values indicate proteins
# downregulated in ALS. Proteins reaching statistical significance
# (FDR < 0.05) in either cohort are highlighted.
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggrepel"))
#===========================Function Definition=============================
# plot discovery vs validation logFC in one scatterplot
# takes datatable with "x" for Discovery logFC, "y" for Validation logFC,
# "fdr_x"/"fdr_y" for the corresponding FDR values, and "name" variable
scatterplot_logFC_discovery_validation = function(data,
                                       fdr_cutoff = 0.05,
                                       q = 0.95,
                                       main_title,
                                       max.overlaps = 10,
                                       lab_x = "log2FC (ALS vs Control) for Discovery",
                                       lab_y = "log2FC (ALS vs Control) for Validation",
                                       labels_T_F = T,
                                       annotate_YN = T,
                                       text_y = "significant in Validation",
                                       text_x = "significant in Discovery"){
  data$omic_type = rep("ns", nrow(data))
  data$omic_type[data$fdr_y < fdr_cutoff] = text_y
  data$omic_type[data$fdr_x < fdr_cutoff] = text_x
  data$omic_type[(data$fdr_x < fdr_cutoff) & (data$fdr_y < fdr_cutoff)] = "significant in both"
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
    geom_point(data = filter(data, fdr_y < fdr_cutoff | fdr_x < fdr_cutoff),
               aes(colour = omic_type),
               alpha = 0.5,
               shape = 16,
               size = 3) +
    geom_smooth(method = "lm", color = "#2C3E50", se = TRUE) +
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
      annotate("text", x = -0.75, y = 1.2, label =
                 paste0(sum(data$omic_type==text_y), " ", text_y,"\n",
                        sum(data$omic_type==text_x), " ", text_x,"\n",
                        sum(data$omic_type=="significant in both"), " significant in both"), size = 8/.pt)
    }

  return(plot)}
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "03_Differential_expression_analysis/Differential_Expression_Results.rds",
              help = "03_Differential_expression_analysis/Differential_Expression_Results.rds"
  ),make_option(c("--output", "-o"),
                type = "character", default = "20_scatterplot_logFC",
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
  stop("Please provide the cleaned SummarizedExperiment object after normalize and imputation file path!")
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

################################################################################################
#MAKE PLOT WITH logFC
inter = intersect(res_Discovery$norm_imp_MinProb_age_sex_cov_all_patients$name,
                  res_Validation$norm_imp_MinProb_age_sex_cov_all_patients$name)

res_Discovery_inter = res_Discovery$norm_imp_MinProb_age_sex_cov_all_patients[which(res_Discovery$norm_imp_MinProb_age_sex_cov_all_patients$name %in% inter),]
# sort res_Discovery by the order in inter
res_Discovery_inter = res_Discovery_inter[match(inter, res_Discovery_inter$name),]

res_Validation_inter = res_Validation$norm_imp_MinProb_age_sex_cov_all_patients[which(res_Validation$norm_imp_MinProb_age_sex_cov_all_patients$name %in% inter),]
# sort res_Validation by the order in inter
res_Validation_inter = res_Validation_inter[match(inter, res_Validation_inter$name),]


data = as.data.table(cbind(res_Discovery_inter$name,
                           res_Discovery_inter$als_vs_ctrl_diff,
                           res_Validation_inter$als_vs_ctrl_diff,
                           res_Discovery_inter$fdr,
                           res_Validation_inter$fdr))

# x is Discovery logFC, y is Validation logFC
colnames(data) = c("name", "x", "y", "fdr_x", "fdr_y")
data$x = as.numeric(data$x)
data$y = as.numeric(data$y)
data$fdr_x = as.numeric(data$fdr_x)
data$fdr_y = as.numeric(data$fdr_y)

p = scatterplot_logFC_discovery_validation(data, q = 0.95, main_title = "ALS vs Control (log2FC)", max.overlaps = Inf)

ggsave(paste0(output_dir, "/ALS_vs_Control_logFC.pdf"), p, width = 7, height = 5)
ggsave(paste0(output_dir, "/ALS_vs_Control_logFC.svg"), p, width = 7, height = 5)

sink(paste0(output_dir, "/logFC_cor.txt"))
cor.test(data$x, data$y, method = "pearson")
cor.test(data$x, data$y, method = "spearman")
sink()
