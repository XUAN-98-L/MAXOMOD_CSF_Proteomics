#=========================Script Description=================================
# This script is used for AUC comparison analysis for ML models
# Rscript 17_AUC_comparision.R -o 17_AUC_comparision

#===========================Loading Packages=============================
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))
suppressMessages(library("gridExtra"))
suppressMessages(library("optparse"))

#============================================================================
# Set up command line argument parsing
option_list <- list(
  make_option(c("--sd"), type = "integer", default = 3, 
              help = "Standard deviation parameter (default: 3)"),
  make_option(c("--input"), type = "character", default = NULL,
              help = "Input path for ML results folders (default: current directory)"),
  make_option(c("--output_dir"), type = "character", default = NULL,
              help = "Output directory (default: auto-generated based on sd)")
)

# Parse command line arguments
parser <- OptionParser(option_list = option_list, 
                      description = "AUC Comparison Analysis for ML Models")
args <- parse_args(parser)

# Use parsed arguments
sd <- args$sd
input_path <- if (is.null(args$input)) getwd() else args$input

# Set output directory
if (is.null(args$output_dir)) {
  output_dir <- paste0("17_AUC_comparision_noAge_", sd, "sd")
} else {
  output_dir <- args$output_dir
}
# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

ML_results_10_input_folder = paste0(input_path, "/15_ML_multi_models_500_01_noAge_",sd,"sd")
ML_results_10 = readRDS(paste0(ML_results_10_input_folder, "/ML_results.rds"))
ML_results_10_proteins = read.delim(paste0(ML_results_10_input_folder, "/selected_features.txt"),header = FALSE)

ML_results_20_input_folder = paste0(input_path, "/15_ML_multi_models_500_02_noAge_",sd,"sd")
ML_results_20 = readRDS(paste0(ML_results_20_input_folder, "/ML_results.rds"))
ML_results_20_proteins = read.delim(paste0(ML_results_20_input_folder, "/selected_features.txt"),header = FALSE)

ML_results_30_input_folder = paste0(input_path, "/15_ML_multi_models_500_03_noAge_",sd,"sd")
ML_results_30 = readRDS(paste0(ML_results_30_input_folder, "/ML_results.rds"))
ML_results_30_proteins = read.delim(paste0(ML_results_30_input_folder, "/selected_features.txt"),header = FALSE)

ML_results_40_input_folder = paste0(input_path, "/15_ML_multi_models_500_04_noAge_",sd,"sd")
ML_results_40 = readRDS(paste0(ML_results_40_input_folder, "/ML_results.rds"))
ML_results_40_proteins = read.delim(paste0(ML_results_40_input_folder, "/selected_features.txt"),header = FALSE)

ML_results_50_input_folder = paste0(input_path, "/15_ML_multi_models_500_05_noAge_",sd,"sd")
ML_results_50 = readRDS(paste0(ML_results_50_input_folder, "/ML_results.rds"))
ML_results_50_proteins = read.delim(paste0(ML_results_50_input_folder, "/selected_features.txt"),header = FALSE)

ML_results_60_input_folder = paste0(input_path, "/15_ML_multi_models_500_06_noAge_",sd,"sd")
ML_results_60 = readRDS(paste0(ML_results_60_input_folder, "/ML_results.rds"))
ML_results_60_proteins = read.delim(paste0(ML_results_60_input_folder, "/selected_features.txt"),header = FALSE)

# Extract ROC summary data from both models
roc_10 = ML_results_10$roc_summary
roc_20 = ML_results_20$roc_summary
roc_30 = ML_results_30$roc_summary
roc_40 = ML_results_40$roc_summary
roc_50 = ML_results_50$roc_summary
roc_60 = ML_results_60$roc_summary

# Combine the data into a single dataframe
roc_combined = bind_rows(roc_10, roc_20, roc_30, roc_40, roc_50, roc_60)

roc_10 <- roc_10 %>% dplyr::mutate(method = factor(method, levels = c("lasso","rf","ridge","svm_linear","svm_radial","xgb")))
roc_20 <- roc_20 %>% dplyr::mutate(method = factor(method, levels = c("lasso","rf","ridge","svm_linear","svm_radial","xgb")))
roc_30 <- roc_30 %>% dplyr::mutate(method = factor(method, levels = c("lasso","rf","ridge","svm_linear","svm_radial","xgb")))
roc_40 <- roc_40 %>% dplyr::mutate(method = factor(method, levels = c("lasso","rf","ridge","svm_linear","svm_radial","xgb")))
roc_50 <- roc_50 %>% dplyr::mutate(method = factor(method, levels = c("lasso","rf","ridge","svm_linear","svm_radial","xgb")))
roc_60 <- roc_60 %>% dplyr::mutate(method = factor(method, levels = c("lasso","rf","ridge","svm_linear","svm_radial","xgb")))
# ... existing code ...

# Extract AUC data for each method and percentage
auc_data <- data.frame(
  percentage = rep(c(10, 20, 30, 40, 50, 60), each = length(c("lasso","rf","ridge","svm_linear","svm_radial","xgb"))),
  method = rep(c("lasso","rf","ridge","svm_linear","svm_radial","xgb"), times = 6),
  mean_auc = c(
    roc_10 %>% group_by(method) %>% summarise(mean_auc = mean(mean_auc)) %>% pull(mean_auc),
    roc_20 %>% group_by(method) %>% summarise(mean_auc = mean(mean_auc)) %>% pull(mean_auc),
    roc_30 %>% group_by(method) %>% summarise(mean_auc = mean(mean_auc)) %>% pull(mean_auc),
    roc_40 %>% group_by(method) %>% summarise(mean_auc = mean(mean_auc)) %>% pull(mean_auc),
    roc_50 %>% group_by(method) %>% summarise(mean_auc = mean(mean_auc)) %>% pull(mean_auc),
    roc_60 %>% group_by(method) %>% summarise(mean_auc = mean(mean_auc)) %>% pull(mean_auc)
  ),
  sd_auc = c(
    roc_10 %>% group_by(method) %>% summarise(sd_auc = mean(sd_auc)) %>% pull(sd_auc),
    roc_20 %>% group_by(method) %>% summarise(sd_auc = mean(sd_auc)) %>% pull(sd_auc),
    roc_30 %>% group_by(method) %>% summarise(sd_auc = mean(sd_auc)) %>% pull(sd_auc),
    roc_40 %>% group_by(method) %>% summarise(sd_auc = mean(sd_auc)) %>% pull(sd_auc),
    roc_50 %>% group_by(method) %>% summarise(sd_auc = mean(sd_auc)) %>% pull(sd_auc),
    roc_60 %>% group_by(method) %>% summarise(sd_auc = mean(sd_auc)) %>% pull(sd_auc)
  )
)


# Create the AUC comparison plot with one line per method
ggplot(auc_data, aes(x = factor(percentage), y = mean_auc, color = method, group = method)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc), 
                width = 0.2) +
  geom_line(size = 1) +
  labs(title = "AUC Comparison Across Different LASSO cutoffs",
       x = "LASSO cutoffs (%)",
       y = "Mean AUC",
       color = "Method",
       subtitle = "Error bars show ±1 standard deviation") +
  scale_color_brewer(palette = "Set2") +  # Use a refined color palette
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        legend.position = "right",  # Move legend to the right side
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11))




# Save the plot
ggsave(file.path(output_dir, "AUC_comparison_by_method.pdf"), 
       width = 6, height = 4, dpi = 300)

# Create a barplot showing the number of proteins for each lasso cutoff
protein_counts <- data.frame(
  percentage = c(10, 20, 30, 40, 50, 60),
  protein_count = c(
    length(ML_results_10_proteins$V1),
    length(ML_results_20_proteins$V1),
    length(ML_results_30_proteins$V1),
    length(ML_results_40_proteins$V1),
    length(ML_results_50_proteins$V1),
    length(ML_results_60_proteins$V1)
  )
)

# Create the barplot
protein_barplot <- ggplot(protein_counts, aes(x = factor(percentage), y = protein_count)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_text(aes(label = protein_count), vjust = -0.5, size = 4, fontface = "bold") +
  labs(title = "Number of Selected Proteins by LASSO Cutoff",
       x = "LASSO Cutoffs (%)",
       y = "Number of Proteins") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5))

# Display the plot
print(protein_barplot)

# Save the protein count barplot
ggsave(file.path(output_dir, "protein_count_by_lasso_cutoff.pdf"), 
       width = 6, height = 4, dpi = 300)

# # Also create a combined plot with both AUC and protein count
# combined_plot <- grid.arrange(
#   ggplot(auc_data, aes(x = factor(percentage), y = mean_auc, color = method, group = method)) +
#     geom_point(size = 3) +
#     geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc), 
#                   width = 0.2) +
#     geom_line(size = 1) +
#     labs(title = "AUC Comparison Across Different LASSO cutoffs",
#          x = "LASSO cutoffs (%)",
#          y = "Mean AUC",
#          color = "Method",
#          subtitle = "Error bars show ±1 standard deviation") +
#     scale_color_brewer(palette = "Set2") +
#     theme_minimal() +
#     theme(axis.text.x = element_text(size = 12),
#           axis.text.y = element_text(size = 12),
#           axis.title = element_text(size = 14),
#           plot.title = element_text(size = 16, hjust = 0.5),
#           plot.subtitle = element_text(size = 12, hjust = 0.5),
#           legend.position = "right",
#           legend.title = element_text(size = 12),
#           legend.text = element_text(size = 11)),
  
#   protein_barplot,
#   ncol = 1, nrow = 2
# )


auc_plot =  ggplot(auc_data, aes(x = factor(percentage), y = mean_auc, color = method, group = method)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc), 
                  width = 0.2) +
    geom_line(size = 1) +
    labs(title = "AUC Comparison Across Different LASSO cutoffs",
         x = "LASSO cutoffs (%)",
         y = "Mean AUC",
         color = "Method",
         subtitle = "Error bars show ±1 standard deviation") +
    scale_color_brewer(palette = "Set2") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 16, hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5),
          legend.position = "right",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11))

library(cowplot)
combined_plot <- plot_grid(
  auc_plot + theme(legend.position = "bottom"),  
  protein_barplot,                               
  ncol = 1, align = "v", axis = "lr"             
)


# Save the combined plot
ggsave(file.path(output_dir, "AUC_and_protein_count_combined.pdf"), 
       plot = combined_plot, width = 8, height = 10, dpi = 300)