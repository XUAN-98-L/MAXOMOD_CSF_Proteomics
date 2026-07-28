#=========================Script Description=================================
# Nfl (neurofilament light chain) by subgroup (control, alpha, beta),
# pooling samples from the Discovery and Validation cohorts together. Nfl
# here is the clinically measured biomarker stored in colData(se)$Nfl
# (collected directly from patients), NOT a proteomics assay row.
#   - one boxplot of Nfl (raw values) across control/alpha/beta, all
#     samples from both cohorts pooled (jitter point shape marks cohort)
#   - one single-row heatmap of Nfl, z-scored WITHIN each cohort (to
#     remove cohort-level scale differences before pooling), samples
#     ordered from lowest to highest z-score, with top annotations for
#     Subgroup and Cohort
#
# Control samples come from the full (all-patients) normalised/imputed
# SummarizedExperiment; ALS patients are split into alpha/beta using the
# k=2 clustering solution (kmeans_k=2), the same approach used in
# 26_vsCtr_proteins.R. Theta and unlabeled samples are excluded.
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggpubr"))
suppressMessages(library("ComplexHeatmap"))
suppressMessages(library("circlize"))
#===========================Function Definition=============================
# Per-sample value of one clinical colData variable + control/alpha/beta
# group, for one cohort. Returns a data.frame(cohort, sample, group, value).
get_clinical_data = function(cohort, intensity_path, clusters_path, variable){
  se = readRDS(file.path(cohort, intensity_path))
  meta = as.data.frame(colData(se))

  if (!(variable %in% colnames(meta))){
    stop(cohort, ": clinical variable '", variable, "' not found in colData.")
  }

  value = as.numeric(meta[[variable]])
  names(value) = as.character(meta$label)

  condition = as.character(meta$condition)
  names(condition) = as.character(meta$label)

  clusters = read.csv(file.path(cohort, clusters_path), check.names = FALSE)
  subgroup = setNames(as.character(clusters$`kmeans_k=2`), clusters$patid)

  group = ifelse(condition == "ctrl", "control", subgroup[names(condition)])
  group = factor(group, levels = c("control", "alpha", "beta"))
  names(group) = names(condition)

  if (any(is.na(group))){
    message(cohort, ": dropping ", sum(is.na(group)), " sample(s) with no control/alpha/beta label (e.g. theta).")
  }
  missing_value = is.na(value[names(group)]) & !is.na(group)
  if (any(missing_value)){
    message(cohort, ": dropping ", sum(missing_value), " sample(s) with missing ", variable, " value.")
  }

  keep = !is.na(group) & !is.na(value[names(group)])

  data.frame(cohort = cohort,
            sample = names(group)[keep],
            group = group[keep],
            value = value[names(group)][keep],
            stringsAsFactors = FALSE)}

# Boxplot of the clinical variable across control/alpha/beta, pooling all
# cohorts (jitter point shape marks which cohort each sample is from).
plot_clinical_boxplot = function(data, group_colors, variable){
  ggplot(data, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.6) +
    geom_jitter(aes(shape = cohort), width = 0.15, size = 1.5, alpha = 0.6, colour = "grey30") +
    stat_compare_means(comparisons = list(c("control", "alpha"), c("control", "beta"), c("alpha", "beta")),
                       method = "wilcox.test", label = "p.format", size = 3, tip.length = 0.01) +
    scale_fill_manual(values = group_colors) +
    scale_shape_manual(values = c(Discovery = 16, Validation = 17)) +
    labs(title = paste0(variable, " by subgroup"),
         x = NULL, y = variable, fill = NULL, shape = "Cohort") +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))}

# Single-row heatmap of the clinical variable's z-score (computed WITHIN
# each cohort before pooling), samples ordered lowest to highest, with top
# annotations for Subgroup and Cohort.
plot_clinical_heatmap = function(data, group_colors, cohort_colors, variable){
  ordered = data[order(data$z_value), ]

  sample_id = paste(ordered$cohort, ordered$sample, sep = "::")
  mat = matrix(ordered$z_value, nrow = 1, dimnames = list(paste0(variable, " (z-score)"), sample_id))

  top_anno = HeatmapAnnotation(
    Subgroup = ordered$group,
    Cohort = ordered$cohort,
    col = list(Subgroup = group_colors, Cohort = cohort_colors),
    annotation_name_side = "left",
    simple_anno_size = unit(4, "mm")
  )

  col_fun = colorRamp2(
    c(min(ordered$z_value, na.rm = TRUE), 0, max(ordered$z_value, na.rm = TRUE)),
    c("#1B5B9D", "white", "#D7191C")
  )

  Heatmap(mat,
         name = "z-score",
         col = col_fun,
         cluster_rows = FALSE,
         cluster_columns = FALSE,
         show_column_names = FALSE,
         show_row_names = TRUE,
         row_names_side = "left",
         top_annotation = top_anno,
         column_title = paste0(variable, " z-score, Discovery + Validation pooled (lowest to highest)"),
         border = TRUE)}
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--intensity", "-i"),
              type = "character", default = "02_Missing_Inspection/norm_imp_MinProb.rds",
              help = "normalised, imputed SummarizedExperiment (all patients: control + ALS), relative to each cohort folder. Only colData is used."
  ),make_option(c("--clusters", "-c"),
                type = "character", default = "08_Clustering_als/cluster_assignments_2.csv",
                help = "cluster_assignments_2.csv (patid + kmeans_k=2 alpha/beta), relative to each cohort folder."
  ),make_option(c("--variable", "-v"),
                type = "character", default = "Nfl",
                help = "colData clinical variable to plot. Default Nfl."
  ),make_option(c("--output", "-o"),
                type = "character", default = "31_NEFL",
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

if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
}else{
  seed = opt$seed
  set.seed(seed)
}

variable = opt$variable
group_colors = c(control = "grey70",
                 alpha = RColorBrewer::brewer.pal(8, "Set2")[2],
                 beta = RColorBrewer::brewer.pal(8, "Set2")[3])
cohort_colors = c(Discovery = "salmon", Validation = "#26b3ff")

################################################################################################
#PULL ALL SAMPLES TOGETHER (DISCOVERY + VALIDATION)
cohorts = c("Discovery", "Validation")
all_data = lapply(cohorts, get_clinical_data, intensity_path = opt$intensity,
                  clusters_path = opt$clusters, variable = variable)

combined = do.call(rbind, all_data)
combined$cohort = factor(combined$cohort, levels = cohorts)

# z-score within each cohort separately, so cohort-level scale differences
# don't distort the pooled heatmap ordering
combined$z_value = ave(combined$value, combined$cohort, FUN = function(x) as.numeric(scale(x)))

write.csv(combined, paste0(output_dir, "/", variable, "_all_samples.csv"), row.names = FALSE)

################################################################################################
#PLOTS
box_plot = plot_clinical_boxplot(combined, group_colors, variable)
ggsave(paste0(output_dir, "/", variable, "_boxplot.pdf"), box_plot, width = 5.5, height = 5, dpi = 300)

ht = plot_clinical_heatmap(combined, group_colors, cohort_colors, variable)
pdf(paste0(output_dir, "/", variable, "_heatmap.pdf"), width = 12, height = 3)
draw(ht)
dev.off()
