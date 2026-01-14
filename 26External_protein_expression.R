# Load required libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(SummarizedExperiment))

# Set up command line options
option_list <- list(
  make_option(c("-i", "--input"), 
              type="character", 
              default="26External_data/norm_imp_MinProb_als.rds",
              help="Path to input RDS file [default: %default]"),
  make_option(c("-f", "--features"), 
              type="character", 
              default="31_ML_multi_models_20250908_500_04_noR2Age_updated/selected_features.txt",
              help="Path to selected features file [default: %default]"),
  make_option(c("-o", "--output"), 
              type="character", 
              default="26External_data/26External_protein_expression/",
              help="Output directory [default: %default]"),
  make_option(c("--width"), 
              type="numeric", 
              default=8,
              help="Plot width in inches [default: %default]"),
  make_option(c("--height"), 
              type="numeric", 
              default=3,
              help="Plot height in inches [default: %default]"),
  make_option(c("--title"), 
              type="character", 
              default="Protein Expression in External Cohort",
              help="Plot title [default: %default]"),
  make_option(c("-v", "--verbose"), 
              action="store_true", 
              default=FALSE,
              help="Print extra output [default: %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Print options if verbose
if (opt$verbose) {
  cat("Running with options:\n")
  print(opt)
}

######Plot
input = readRDS(opt$input)
# Generate Boxplots for Feature Expression in Discovery Cohort
plot_discovery <- assay(input)
plot_discovery <- as.data.frame(plot_discovery)
plot_discovery = t(plot_discovery)
selected_features = read.delim(opt$features, header = F)
plot_discovery = plot_discovery[,colnames(plot_discovery) %in% selected_features$V1]
plot_discovery = as.data.frame(plot_discovery)
plot_discovery$outcome <- input$condition

output_dir = opt$output
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  if (opt$verbose) {
    cat("Created output directory:", output_dir, "\n")
  }
}

# Convert to long format for ggplot
plot_long_discovery <- data.table::melt(plot_discovery, id.vars = "outcome", variable.name = "Protein", value.name = "Expression")

# Create the boxplot for Discovery Cohort
boxplot_discovery <- ggplot(plot_long_discovery, aes(x = as.factor(outcome), y = Expression, fill = as.factor(outcome))) +
  geom_boxplot(alpha = 1, outlier.shape = NA) +  
  geom_jitter(width = 0.2, alpha = 1, size = 1) +  
  facet_wrap(~Protein, scales = "free_y") +  
  stat_compare_means(aes(label = ..p.format..), method = "wilcox.test", label = "p.format") +  
  theme_classic() +
  labs(title = opt$title, 
       x = "ALS Subtype", y = "log2 Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = brewer.pal(8, "Set2")[c(2, 3)]) 

# Save External Cohort boxplot
output_file = file.path(output_dir, "boxplot_external.pdf")
ggsave(filename = output_file, plot = boxplot_discovery, width = opt$width, height = opt$height)

if (opt$verbose) {
  cat("Plot saved to:", output_file, "\n")
}
