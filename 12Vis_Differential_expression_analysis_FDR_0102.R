#=========================Script Description=================================
# This script is used for visualize differential expression analysis between als vs ctrl result (volcano plot)
# cd /Users/xliu2942/Documents/Projects/MAXOMOD/TESTing/Pipeline
# Rscript 12Vis_Differential_expression_analysis.R
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
              type = "character", default = "3_Differential_Expression_Analysis/Differential_Expression_Results.rds",
              help = "3_Differential_Expression_Analysis/Differential_Expression_Results.rds"
  ),make_option(c("--output", "-o"),
                type = "character", default = "12_Vis_Differential_Expression_Analysis",
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
names_dfs <- c("norm_imp_MinProb_age_sex_cov_all_patients", 
               "norm_imp_MinProb_age_cov_only_female",
               "norm_imp_MinProb_age_cov_only_male")

plots_FDR0.2 = plots_FDR0.1 = list()

for (i in seq_along(names_dfs)) {
  data_res <- res[[names_dfs[i]]]
  
  df <- data.frame(
    x = data_res[[grep("diff", colnames(data_res))]],
    y = -log10(data_res[[grep("fdr", colnames(data_res))]]),
    name = data_res$name
  )
  
  top10_label_up <- df %>%
    filter(x >= 0 & y >= -log10(0.2)) %>%
    arrange(desc(x)) %>%
    slice_head(n = 10)
  
  top10_label_down <- df %>%
    filter(x <= 0 & y >= -log10(0.2)) %>%
    arrange(x) %>%
    slice_head(n = 10)

  case_color <- "#D73027"
  ctrl_color <- "#4575B4"
  
  plots_FDR0.2[[i]] <- volcano_plot(df, 0.2, 
                                     paste0("Volcano plot FDR 0.2\n", names_dfs[i]), 
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
                                    paste0("Volcano plot FDR 0.1\n", names_dfs[i]), 
                                    c(top10_label_up$name, top10_label_down$name, label),
                                    case_color, ctrl_color)
}

#===========================Save Plots======================================
ggarrange(plotlist = plots_FDR0.2, ncol = 2, nrow = 2)
ggsave(filename = file.path(output_dir, "FDR02.pdf"), width = 12, height = 10, units = "in")

ggarrange(plotlist = plots_FDR0.1, ncol = 2, nrow = 2)
ggsave(filename = file.path(output_dir, "FDR01.pdf"), width = 12, height = 10, units = "in")