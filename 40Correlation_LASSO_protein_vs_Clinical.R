#=========================Script Description=================================
# This script is used for correlation analysis between LASSO selected proteins and clinical variables
# cd /Users/xliu2942/Documents/Projects/MAXOMOD/Script
# Rscript 40Correlation_LASSO_protein_vs_Clinical.r --input Discovery/2_Missing_Inspection_als/norm_imp_MinProb.rds -s 31_ML_multi_models_20250819_500_04_test/selected_features.txt -c Discovery/5_Clustering_als/cluster_assignments_2.csv -g Discovery/8_Differential_expression_analysis_subclusters/Differential_expression_analysis_for_k2.csv -o 40_Correlation_LASSO_protein_vs_Clinical_Discovery 
# For Validation
# Rscript 40Correlation_LASSO_protein_vs_Clinical.r --input Validation/2_Missing_Inspection_als/norm_imp_MinProb.rds -s 31_ML_multi_models_20250819_500_04_test/selected_features.txt -c Validation/5_Clustering_als/cluster_assignments_2.csv -g Validation/8_Differential_expression_analysis_subclusters/Differential_expression_analysis_for_k2.csv -o 40_Correlation_LASSO_protein_vs_Clinical_Validation 
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("ggplot2"))
suppressMessages(library("gridExtra"))
suppressMessages(library("ggpubr"))
#===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "2_Missing_Inspection_als/norm_imp_MinProb.rds",
              help = "Input data RDS file path."
  ),
  make_option(c("--selected_features", "-s"),
              type = "character", default = "31_ML_multi_models_20250819_500_04_test/selected_features.txt",
              help = "Selected features file path."
  ),
  make_option(c("--cluster", "-c"),
              type = "character", default = "5_Clustering_als/cluster_assignments_2.csv",
              help = "Cluster assignments file path."
  ),
  make_option(c("--deg_results", "-g"),
              type = "character", default = "8_Differential_expression_analysis_subclusters/Differential_expression_analysis_for_k2.csv",
              help = "DEG results file path."
  ),
  make_option(c("--output", "-o"),
              type = "character", default = "40_Correlation_LASSO_protein_vs_Clinical",
              help = "Output directory path."
  ),
  make_option(c("--p_threshold", "-p"),
              type = "numeric", default = 0.1,
              help = "P-value threshold for significant correlations."
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
#============================================================================
# Parameter validation
if (is.null(opt$output)) {
  print("NO OUTPUT PATH SUPPLIED, current directory will be used!")
  output_dir <- getwd()
} else {
  output_dir <- opt$output
  if (!file.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
}

if (!file.exists(opt$input)) {
  stop("The input data file does not exist!")
}

if (!file.exists(opt$selected_features)) {
  stop("The selected features file does not exist!")
}

if (!file.exists(opt$cluster)) {
  stop("The cluster assignments file does not exist!")
}

if (!file.exists(opt$deg_results)) {
  stop("The DEG results file does not exist!")
}

#===========================Data Loading=============================
Input = readRDS(opt$input)
Input_df = as.data.frame(colData(Input))
Selected_features = read.delim(opt$selected_features, header = FALSE)
Cluster = read.csv(opt$cluster)
DEG_results = read.csv(opt$deg_results)

Input_expression = assay(Input)[Selected_features$V1,]

# Transpose expression data to have samples as rows and proteins as columns
Input_expression_t = t(Input_expression)
colnames(Input_expression_t) = Selected_features$V1

# Add expression data to clinical information
Input_df = cbind(Input_df, Input_expression_t)

# Create correlation plots between selected features and clinical variables
library(ggplot2)
library(gridExtra)
library(ggpubr)

# Define clinical variables of interest
clinical_vars <- c("Nfl", "pNFh", "progression_rate", "age", "slow_vital_capacity", "age_at_onset", "disease_duration")

# Get selected protein names and their differential expression direction
protein_names <- Selected_features$V1

cat("Number of selected proteins:", length(protein_names), "\n")
cat("Selected proteins:", paste(protein_names, collapse = ", "), "\n\n")

# Check which proteins exist in the expression data
protein_names <- protein_names[protein_names %in% rownames(Input_expression)]
cat("Proteins found in expression data:", length(protein_names), "\n")
cat("Available proteins:", paste(protein_names, collapse = ", "), "\n\n")

# Create a mapping of protein names to their colors based on DEG results
protein_colors <- sapply(protein_names, function(protein) {
  # Find the protein in DEG_results
  deg_row <- DEG_results[DEG_results$name == protein, ]
  if (nrow(deg_row) > 0) {
    diff_value <- deg_row$alpha_vs_beta_diff[1]
    if (diff_value > 0) {
      return("#FC8D62")  # Blue for higher expression in beta
    } else if (diff_value < 0) {
      return("#8DA0CB")  # Orange for higher expression in alpha
    } else {
      return("gray")     # Gray for no difference
    }
  } else {
    return("gray")       # Gray if protein not found in DEG results
  }
})

# Create plot list
plot_list <- list()
plot_count <- 0

for (protein in protein_names) {
  for (trait in clinical_vars) {
    # Check if both variables exist in the data
    if (protein %in% colnames(Input_df) && trait %in% colnames(Input_df)) {
      df <- Input_df[, c(protein, trait)]
      df <- df[complete.cases(df), ]
      
      cat(paste("Testing", protein, "vs", trait, "-", nrow(df), "complete cases\n"))
      
      # Only create plot if we have data
      if (nrow(df) > 0) {
        # Pearson correlation test
        cor_test <- cor.test(df[[1]], df[[2]], method = "pearson")
        
        cat(paste("  Correlation:", round(cor_test$estimate, 3), "p-value:", round(cor_test$p.value, 4), "\n"))
        
        # Create plot if p-value is significant
        if (!is.na(cor_test$p.value) && cor_test$p.value < opt$p_threshold) {
          plot_count <- plot_count + 1
          cat(paste("  *** SIGNIFICANT! Adding plot", plot_count, "***\n"))
          
          # Get the color for this protein
          dot_color <- protein_colors[protein]
          
          p <- ggplot(df, aes_string(x = trait, y = protein)) +
            geom_point(color = dot_color, size = 2.5, alpha = 0.7) +
            geom_smooth(method = "lm", se = TRUE, color = "black") +
            stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 3.5) +
            theme_classic(base_size = 12) +
            labs(
              x = trait,
              y = protein,
              title = paste(protein, "vs", trait)
            )
          plot_list[[paste(protein, trait, sep = "_")]] <- p
        }
      }
    }
  }
}

# Create the combined plot
if (length(plot_list) > 0) {
  n_show <- min(12, length(plot_list))
  protein_plot = do.call(gridExtra::grid.arrange, c(plot_list[1:n_show], ncol = 4))
  
  # Calculate height based on number of rows
  n_cols <- 4  # Fixed number of columns
  n_rows <- ceiling(n_show / n_cols)  # Calculate number of rows needed
  plot_height <- n_rows * 4  # 4 inches per row (you can adjust this multiplier)
  
  # Ensure minimum and maximum heights
  plot_height <- max(8, min(20, plot_height))  # Between 8 and 20 inches
  
  # Save the plot
  ggsave(filename = file.path(output_dir, "Selected_Proteins_vs_Clinical_Variables.pdf"), 
         plot = protein_plot, 
         width = 16, 
         height = plot_height)
  
  cat(paste("Created plot with", n_rows, "rows and", n_cols, "columns. Height:", plot_height, "inches\n"))
} else {
  cat(paste("No significant correlations found (p <", opt$p_threshold, ")\n"))
}


