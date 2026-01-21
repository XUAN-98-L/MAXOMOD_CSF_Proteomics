#=========================Script Description=================================
# This script is used for visualize differential expression analysis
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggrepel"))
#===========================Function Definition=============================
#plot male vs female in one scatterplot
# takes datatable with "x" for women, "y" for men, and "name" variables
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
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "03_Differential_expression_analysis/Differential_Expression_Results.rds",
              help = "03_Differential_expression_analysis/Differential_Expression_Results.rds"
  ),make_option(c("--output", "-o"),
                type = "character", default = "20_scatterplot_FDR_signedp",
                help = "output directory path."
  ),make_option(c("--seed", "-e"),
                type = "integer", default = 9,
                help = "set.seed"
  ),make_option(c("--FDR_Cutoff", "-f"),
                type = "numeric", default = 0.05,
                help = "FDR cutoff value (default: 0.05)"
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

if (is.null(opt$FDR_Cutoff)) {
  FDR_Cutoff = 0.05
}else{
  FDR_Cutoff = opt$FDR_Cutoff
}
cut_off = -log10(FDR_Cutoff)

################################################################################################
#MAKE PLOT WITH FDR
inter = intersect(res_Discovery$norm_imp_MinProb_age_sex_cov_all_patients$name,
                  res_Validation$norm_imp_MinProb_age_sex_cov_all_patients$name)

res_Discovery_inter = res_Discovery$norm_imp_MinProb_age_sex_cov_all_patients[which(res_Discovery$norm_imp_MinProb_age_sex_cov_all_patients$name %in% inter),]
# sort res_Discovery by the order in inter
res_Discovery_inter = res_Discovery_inter[match(inter, res_Discovery_inter$name),]

res_Validation_inter = res_Validation$norm_imp_MinProb_age_sex_cov_all_patients[which(res_Validation$norm_imp_MinProb_age_sex_cov_all_patients$name %in% inter),]
# sort res_Validation by the order in inter
res_Validation_inter = res_Validation_inter[match(inter, res_Validation_inter$name),]


data = as.data.table(cbind(res_Discovery_inter$name, 
                           res_Discovery_inter$fdr, 
                           res_Validation_inter$fdr))

# x is Discovery
# y is Validation
colnames(data) = c("name", "x", "y")
data$x = as.numeric(data$x)
data$y = as.numeric(data$y)
data$x = -log10(data$x)
data$y = -log10(data$y)

data$x[res_Discovery_inter$als_vs_ctrl_diff <0] = -data$x[res_Discovery_inter$als_vs_ctrl_diff <0]

data$y[res_Validation_inter$als_vs_ctrl_diff <0] = -data$y[res_Validation_inter$als_vs_ctrl_diff <0]

p = scatterplot_FDR_male_female(data, cut_off = cut_off, q = 0.95, main_title = "ALS vs Control(FDR)",max.overlaps = Inf)

ggsave(paste0(output_dir, "/ALS_vs_Control.pdf"), p, width = 6, height = 4, dpi = 300)

sink(paste0(output_dir, "/FDR_cor.txt"))
cor.test(data$x, data$y, method = "pearson")
sink()

# Create summary dataframe with labels, pvalues, FDRs, and fold changes
summary_df = data.frame(
  name = res_Discovery_inter$name,
  Discovery_pvalue = res_Discovery_inter$als_vs_ctrl_p.val,
  Validation_pvalue = res_Validation_inter$als_vs_ctrl_p.val,
  Discovery_FDR = res_Discovery_inter$fdr,
  Validation_FDR = res_Validation_inter$fdr,
  Discovery_foldchange = res_Discovery_inter$als_vs_ctrl_diff,
  Validation_foldchange = res_Validation_inter$als_vs_ctrl_diff,
  stringsAsFactors = FALSE
)

# Determine labels based on FDR cutoff
summary_df$label = "ns"
summary_df$label[summary_df$Validation_FDR < FDR_Cutoff] = "significant in Validation"
summary_df$label[summary_df$Discovery_FDR < FDR_Cutoff] = "significant in Discovery"
summary_df$label[(summary_df$Discovery_FDR < FDR_Cutoff) & (summary_df$Validation_FDR < FDR_Cutoff)] = "significant in both"

# Reorder columns to put label first
summary_df = summary_df[, c("name", "label", "Discovery_pvalue", "Validation_pvalue", 
                            "Discovery_FDR", "Validation_FDR", 
                            "Discovery_foldchange", "Validation_foldchange")]

# Save to CSV
write.csv(summary_df, file = paste0(output_dir, "/Protein_summary_with_labels.csv"), row.names = FALSE)



#########pvalue
data = as.data.table(cbind(res_Discovery_inter$name, 
                           res_Discovery_inter$als_vs_ctrl_p.val, 
                           res_Validation_inter$als_vs_ctrl_p.val))

# x is Discovery
# y is Validation
colnames(data) = c("name", "x", "y")
data$x = as.numeric(data$x)
data$y = as.numeric(data$y)
data$x = -log10(data$x)
data$y = -log10(data$y)

data$x[res_Discovery_inter$als_vs_ctrl_diff <0] = -data$x[res_Discovery_inter$als_vs_ctrl_diff <0]

data$y[res_Validation_inter$als_vs_ctrl_diff <0] = -data$y[res_Validation_inter$als_vs_ctrl_diff <0]

p = scatterplot_FDR_male_female(data, cut_off = cut_off, q = 0.95, main_title = "ALS vs Control(pvalue)",max.overlaps = Inf,lab_x = "signed -log10(pvalue) for Discovery",
                                lab_y = "signed -log10(pvalue) for Validation")

ggsave(paste0(output_dir, "/ALS_vs_Control_pvalue.pdf"), p, width = 6, height = 4, dpi = 300)

sink(paste0(output_dir, "/pvalue_cor.txt"))
cor.test(data$x, data$y, method = "pearson")
sink()
