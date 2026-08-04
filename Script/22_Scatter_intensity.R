#=========================Script Description=================================
# Scatter plot comparing the log2-transformed mean protein intensity between
# the Discovery and Validation cohorts.
#
# Protein intensities were derived from the VSN-normalized, imputed LFQ
# intensity matrix (norm_imp_MinProb.rds). For each protein, the mean
# intensity across all samples within each cohort was calculated (on the
# log2 scale produced by VSN normalization). Only proteins detected in both
# cohorts were included. Cross-cohort agreement was assessed using Pearson's
# and Spearman's correlation coefficients.
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("data.table"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggrepel"))
#===========================Function Definition=============================
# format a p-value the way it is usually reported in figure legends
format_pval = function(p){
  if (p < 2.2e-16){
    return("< 2.2 × 10^-16")
  } else {
    return(paste0("= ", format(p, scientific = TRUE, digits = 3)))
  }
}

# scatter plot of mean log2 intensity, Discovery (x) vs Validation (y)
scatterplot_intensity = function(data,
                                 q = 0.97,
                                 #main_title = "Protein Intensity: Discovery vs Validation",
                                 max.overlaps = 10,
                                 lab_x = "log2(mean LFQ intensity) for Discovery",
                                 lab_y = "log2(mean LFQ intensity) for Validation",
                                 annotate_YN = T){

  pearson = cor.test(data$x, data$y, method = "pearson")
  spearman = cor.test(data$x, data$y, method = "spearman")

  quantile_x = quantile(data$x, na.rm = T, probs = q)
  quantile_y = quantile(data$y, na.rm = T, probs = q)

  plot = ggplot(data, aes(x,y)) +
    geom_point(colour = "black",
               alpha = 0.5,
               shape = 16,
               size = 2) +
    geom_smooth(method = "lm", color = "#2C3E50", se = TRUE) +
    geom_text_repel(data = filter(data, x >= quantile_x | y >= quantile_y),
                    aes(label = name),
                    force = 1,
                    hjust = 1,
                    max.overlaps = max.overlaps,
                    segment.size = 0.2,
                    min.segment.length = 0,
                    size = 2.5)  +
    labs(title = NULL,
         x = lab_x,
         y = lab_y) +
    theme_classic() + # Select theme with a white background
    theme(axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text = element_text(size = 12),
          plot.title = element_blank(),
          text = element_text(size = 14)) +
    if(annotate_YN){
      annotate("text",
                x = min(data$x, na.rm = T) + 0.05 * diff(range(data$x, na.rm = T)),
                y = max(data$y, na.rm = T) - 0.05 * diff(range(data$y, na.rm = T)),
                hjust = 0,
                vjust = 1,
                label = paste0("Pearson's r = ", round(pearson$estimate, 4), ", P ", format_pval(pearson$p.value), "\n",
                               "Spearman's ρ = ", round(spearman$estimate, 4), ", P ", format_pval(spearman$p.value)),
                size = 8/.pt)
    }

  return(list(plot = plot, pearson = pearson, spearman = spearman))}
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "02_Missing_Inspection/norm_imp_MinProb.rds",
              help = "02_Missing_Inspection/norm_imp_MinProb.rds"
  ),make_option(c("--output", "-o"),
                type = "character", default = "22_Scatter_intensity",
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
  stop("Please provide the VSN-normalized, imputed SummarizedExperiment object file path!")
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
#MAKE PLOT WITH MEAN INTENSITY
mean_Discovery = rowMeans(assay(res_Discovery), na.rm = T)
mean_Validation = rowMeans(assay(res_Validation), na.rm = T)

inter = intersect(names(mean_Discovery), names(mean_Validation))

data = data.table(name = inter,
                  x = mean_Discovery[inter],
                  y = mean_Validation[inter])

res = scatterplot_intensity(data, q = 0.97, #main_title = "Protein Intensity: Discovery vs Validation",
                            max.overlaps = Inf)

ggsave(paste0(output_dir, "/Scatter_intensity_Discovery_vs_Validation.pdf"), res$plot,
       width = 6, height = 5, device = cairo_pdf)

sink(paste0(output_dir, "/intensity_cor.txt"))
print(res$pearson)
print(res$spearman)
sink()
