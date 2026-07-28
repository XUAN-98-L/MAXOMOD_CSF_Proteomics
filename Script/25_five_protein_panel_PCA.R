#=========================Script Description=================================
# PCA of the five ML-selected proteins (ATRN, CNTN1, PARK7, PCSK1N, and
# PTPRS) in ALS patients. Samples are coloured by subgroup assignment
# (alpha, beta, theta), with ellipses indicating the 95% confidence regions.
# The first two principal components are shown, with the proportion of
# variance explained indicated on each axis.
#
# Input:
#   <cohort>/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds
#     normalised, imputed SummarizedExperiment (ALS patients only)
#   <cohort>/08_Clustering_als/cluster_assignments_2.csv
#     patid + kmeans_k=3 (alpha/beta/theta) subgroup assignment
#     (matched to samples via colData(se)$label)
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("ggplot2"))
suppressMessages(library("RColorBrewer"))
#===========================Function Definition=============================
# PCA of a fixed protein panel, coloured by subgroup, with 95% ellipses.
# mat: proteins (rows) x samples (columns) numeric matrix (subset to panel).
# subgroup: named factor (alpha/beta/theta), names matching colnames(mat).

plot_panel_pca = function(mat, subgroup, main_title, cluster_colors){
  pca_data = as.data.frame(t(mat))
  pca_res = prcomp(pca_data, center = TRUE, scale. = TRUE)

  var_explained = summary(pca_res)$importance["Proportion of Variance", 1:2] * 100

  pc_df = data.frame(patid = rownames(pca_res$x),
                     PC1 = pca_res$x[, 1],
                     PC2 = pca_res$x[, 2],
                     subgroup = subgroup[rownames(pca_res$x)])

  lab_x = sprintf("PC1 (%.1f%% variance explained)", var_explained[1])
  lab_y = sprintf("PC2 (%.1f%% variance explained)", var_explained[2])

  plot = ggplot(pc_df, aes(x = PC1, y = PC2, colour = subgroup, fill = subgroup)) +
    stat_ellipse(type = "norm", level = 0.95, geom = "polygon", alpha = 0.15, colour = NA) +
    stat_ellipse(type = "norm", level = 0.95, linewidth = 0.6) +
    geom_point(size = 2.5, alpha = 0.8, shape = 16) +
    scale_colour_manual(values = cluster_colors) +
    scale_fill_manual(values = cluster_colors) +
    labs(title = main_title,
         subtitle = "Panel: ATRN, CNTN1, PARK7, PCSK1N, PTPRS",
         x = lab_x,
         y = lab_y,
         colour = "Subgroup",
         fill = "Subgroup") +
    theme_classic() +
    theme(axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 15, hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5),
          text = element_text(size = 14))

  return(list(plot = plot, pca_res = pca_res, scores = pc_df, var_explained = var_explained))
}


# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--cohort", "-C"),
              type = "character", default = NULL,
              help = "Discovery or Validation. Required."
  ),make_option(c("--input", "-i"),
                type = "character", default = "02_Missing_Inspection_subclusters/norm_imp_MinProb.rds",
                help = "normalised, imputed SummarizedExperiment (ALS patients only), relative to the cohort folder."
  ),make_option(c("--clusters", "-c"),
                type = "character", default = "08_Clustering_als/cluster_assignments_2.csv",
                help = "cluster_assignments_2.csv with patid + kmeans_k=3 subgroup labels, relative to the cohort folder."
  ),make_option(c("--output", "-o"),
                type = "character", default = "25_five_protein_panel_PCA",
                help = "output directory path, relative to the cohort folder."
  ),make_option(c("--seed", "-e"),
                type = "integer", default = 9,
                help = "set.seed"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
#============================================================================
if (is.null(opt$cohort) || !(opt$cohort %in% c("Discovery", "Validation"))) {
  stop("Please provide --cohort as either 'Discovery' or 'Validation'!")
}else{
  cohort = opt$cohort
}

if (is.null(opt$output)) {
  print("NO OUTPUT PATH SUPPLIED,current directory will be used!")
  output_dir <- getwd()
} else {
  output_dir <- file.path(cohort, opt$output)
  if (!file.exists(output_dir)) {
    dir.create(output_dir, recursive = T)
  }
}

if (is.null(opt$input)) {
  stop("Please provide the normalised, imputed SummarizedExperiment object file path!")
}else{
  se = readRDS(file.path(cohort, opt$input))
}

if (is.null(opt$clusters)) {
  stop("Please provide the cluster_assignments_2.csv file path!")
}else{
  cluster_assignments = read.csv(file.path(cohort, opt$clusters), check.names = FALSE)
}

if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
}else{
  seed = opt$seed
  set.seed(seed)
}

################################################################################################
#SUBSET TO THE FIVE-PROTEIN PANEL
protein_panel = c("ATRN", "CNTN1", "PARK7", "PCSK1N", "PTPRS")

present = protein_panel[protein_panel %in% rownames(se)]
missing = setdiff(protein_panel, present)
if (length(missing) > 0){
  message("Warning: the following panel proteins were not found and will be skipped: ",
         paste(missing, collapse = ", "))
}
if (length(present) < 2){
  stop("Fewer than two panel proteins were found in the input data; cannot run PCA.")
}

mat = assay(se)[present, , drop = FALSE]

################################################################################################
#SUBGROUP LABELS (alpha/beta/theta)
# SE colnames are site-indexed IDs (spinal_1, bulbar_1, ...); cluster_assignments
# uses clinical labels (CSF55, ...). Match via colData(se)$label, as elsewhere
# in this project (e.g. 03/05/09 scripts).
cluster_assignments$`kmeans_k=3` = factor(as.character(cluster_assignments$`kmeans_k=3`),
                                            levels = c("alpha", "beta", "theta"))
sample_label = as.character(colData(se)$label)
names(sample_label) = colnames(se)
subgroup = cluster_assignments$`kmeans_k=3`[match(sample_label[colnames(mat)],
                                                    cluster_assignments$patid)]
names(subgroup) = colnames(mat)
subgroup = droplevels(subgroup)

no_label = colnames(mat)[is.na(subgroup)]
if (length(no_label) > 0){
  message("Warning: ", length(no_label), " sample(s) had no subgroup label and were dropped: ",
         paste(no_label, collapse = ", "))
  mat = mat[, !is.na(subgroup), drop = FALSE]
  subgroup = subgroup[!is.na(subgroup)]
}
if (ncol(mat) < 3){
  stop("Fewer than 3 labelled samples remain after matching subgroup labels; cannot run PCA.")
}

# Drop constant proteins (cannot scale to unit variance)
protein_var = apply(mat, 1, function(x) var(x, na.rm = TRUE))
zero_var = names(protein_var)[!is.finite(protein_var) | protein_var == 0]
if (length(zero_var) > 0){
  message("Warning: dropping zero-variance protein(s): ", paste(zero_var, collapse = ", "))
  mat = mat[setdiff(rownames(mat), zero_var), , drop = FALSE]
  present = rownames(mat)
}
if (nrow(mat) < 2){
  stop("Fewer than two variable panel proteins remain; cannot run PCA.")
}

################################################################################################
#PCA + PLOT
cluster_colors = brewer.pal(n = 8, "Set2")[c(2, 3, 5)]
names(cluster_colors) = c("alpha", "beta", "theta")

res = plot_panel_pca(mat, subgroup,
                     main_title = paste0("PCA of Five-Protein ML Panel — ", cohort, " Cohort (ALS)"),
                     cluster_colors = cluster_colors)

ggsave(paste0(output_dir, "/Five_protein_panel_PCA_", cohort, ".pdf"), res$plot, width = 6, height = 5, dpi = 300)

write.csv(res$scores, paste0(output_dir, "/PCA_scores_", cohort, ".csv"), row.names = FALSE)
saveRDS(res$pca_res, paste0(output_dir, "/PCA_object_", cohort, ".rds"))

sink(paste0(output_dir, "/PCA_summary_", cohort, ".txt"))
cat("PCA of five-protein ML panel (", paste(present, collapse = ", "), ") - ", cohort, " cohort\n\n", sep = "")
if (length(missing) > 0) cat("Missing proteins (skipped): ", paste(missing, collapse = ", "), "\n\n", sep = "")
cat("Variance explained:\n")
print(round(res$var_explained, 2))
cat("\nSample counts per subgroup:\n")
print(table(res$scores$subgroup))
cat("\nSummary of PCA object:\n")
print(summary(res$pca_res))
sink()
