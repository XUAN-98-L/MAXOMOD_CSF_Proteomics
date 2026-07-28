#=========================Script Description=================================
# Mean protein expression of subtype-defining signatures vs control.
#
# Two signatures are evaluated, each combining a differential-expression
# (DE) based gene set and a WGCNA module:
#   - "alpha-up DEGs" + "turquoise module"  -> inflammation-related proteins
#   - "beta-up DEGs"  + "blue module"       -> synapto-axonal proteins
#
# Both the alpha/beta-up DEG lists and the WGCNA module memberships are
# defined ONCE from the Discovery cohort, then evaluated as FIXED protein
# sets in both cohorts independently:
#   - alpha/beta-up DEGs: Discovery/04_Vis_Differential_expression_analysis_
#     subclusters/all_important_protein_names_k2_alpha.csv / _k2_beta.csv
#     (alpha_vs_beta comparison, k=2 clustering solution, FDR <= 0.05,
#     already filtered by 04_Vis_Differential_expression_analysis_subclusters.R)
#   - modules: Discovery/10_WGCNA_subclusters/moduleLabels.rds
# Only control/alpha/beta samples are used (theta excluded; k=2 is the
# solution supported by 24_clusterboot.R).
#
# For each cohort, every protein is z-scored across that cohort's own
# samples (control + alpha + beta); the per-sample composite score for a
# set is the mean z-score across the proteins in that set. Boxplots use
# ggplot's default geom_boxplot (median/IQR + Tukey whiskers). Discovery
# results are shown on top, Validation on the bottom.
#
# Per-protein boxplots (*_protein_boxplot.pdf) show the ML-selected
# five-protein panel from selected_features.txt (default:
# 15_ML_multi_models_500_04_noAge_3sd/selected_features.txt).
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggpubr"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))
#===========================Function Definition=============================
# Define the four fixed protein sets: alpha-up/beta-up DEGs read directly
# from the k=2 subcluster DE result already filtered by
# 04_Vis_Differential_expression_analysis_subclusters.R (alpha_vs_beta_diff
# sign + FDR <= 0.05), and turquoise/blue WGCNA module membership.
get_signature_gene_sets = function(alpha_degs_path, beta_degs_path, wgcna_path){
  alpha_up = read.csv(alpha_degs_path, check.names = FALSE)$name
  beta_up  = read.csv(beta_degs_path, check.names = FALSE)$name

  module = readRDS(wgcna_path)
  turquoise = module$gene[module$color == "turquoise"]
  blue      = module$gene[module$color == "blue"]

  sets = list("alpha-up DEGs" = unique(alpha_up),
             "beta-up DEGs" = unique(beta_up),
             "turquoise module" = unique(turquoise),
             "blue module" = unique(blue))
  return(sets)}

# Load one cohort's normalised/imputed protein matrix (all patients: control
# + ALS), attach control/alpha/beta group labels (k=2 clustering; theta and
# any unlabeled samples dropped), and z-score every protein across the
# cohort's own samples.
load_cohort_data = function(cohort, intensity_path, clusters_path){
  se = readRDS(file.path(cohort, intensity_path))
  clusters = read.csv(file.path(cohort, clusters_path), check.names = FALSE)

  mat = assay(se)
  colnames(mat) = as.character(colData(se)$label)

  condition = as.character(colData(se)$condition)
  names(condition) = colnames(mat)

  subgroup = setNames(as.character(clusters$`kmeans_k=2`), clusters$patid)

  group = ifelse(condition == "ctrl", "control", subgroup[names(condition)])
  group = factor(group, levels = c("control", "alpha", "beta"))
  names(group) = names(condition)

  keep = !is.na(group)
  if (any(!keep)){
    message(cohort, ": dropping ", sum(!keep), " sample(s) with no control/alpha/beta label (e.g. theta).")
  }
  mat = mat[, keep, drop = FALSE]
  group = group[keep]

  z = t(scale(t(mat)))

  return(list(cohort = cohort, mat = mat, z = z, group = group))}

# Composite (mean z-score) per sample for each protein set, long format:
# cohort, sample, group, set, score. Also returns which set proteins were
# actually found in this cohort's data (present_genes).
build_set_scores = function(cohort_data, gene_sets){
  z = cohort_data$z
  group = cohort_data$group

  score_list = list()
  present_genes = list()
  for (set_name in names(gene_sets)){
    genes = intersect(gene_sets[[set_name]], rownames(z))
    present_genes[[set_name]] = genes
    if (length(genes) == 0){
      message(cohort_data$cohort, ": no proteins found for set '", set_name, "', skipping.")
      next
    }
    score = if (length(genes) == 1) z[genes, ] else colMeans(z[genes, , drop = FALSE], na.rm = TRUE)
    score_list[[set_name]] = data.frame(cohort = cohort_data$cohort,
                                        sample = names(score),
                                        group = group[names(score)],
                                        set = set_name,
                                        score = as.numeric(score),
                                        stringsAsFactors = FALSE)
  }
  set_long = do.call(rbind, score_list)
  rownames(set_long) = NULL
  return(list(set_long = set_long, present_genes = present_genes))}

# Protein-level long data (z-scored expression), one row per sample x
# protein x set (a protein that belongs to no set is not included).
build_protein_long = function(cohort_data, present_genes){
  z = cohort_data$z
  group = cohort_data$group

  rows = list()
  for (set_name in names(present_genes)){
    genes = present_genes[[set_name]]
    if (length(genes) == 0) next
    sub = as.data.frame(t(z[genes, , drop = FALSE]))
    sub$sample = rownames(sub)
    sub$group = group[sub$sample]
    long = tidyr::pivot_longer(sub, cols = -c(sample, group), names_to = "protein", values_to = "score")
    long$set = set_name
    long$cohort = cohort_data$cohort
    rows[[set_name]] = long
  }
  out = do.call(rbind, rows)
  rownames(out) = NULL
  return(out)}

# Protein-level long data for a flat protein list (e.g. ML selected features),
# one row per sample x protein. Uses normalised log2 intensities (not z-scored)
# so per-protein boxplots match the "Expression (log2, normalized)" scale.
build_panel_protein_long = function(cohort_data, proteins){
  mat = cohort_data$mat
  group = cohort_data$group
  genes = intersect(proteins, rownames(mat))
  missing = setdiff(proteins, genes)
  if (length(missing) > 0){
    message(cohort_data$cohort, ": ML panel protein(s) not found and skipped: ",
            paste(missing, collapse = ", "))
  }
  if (length(genes) == 0){
    return(data.frame(cohort = character(), sample = character(), group = character(),
                      protein = character(), score = numeric(), stringsAsFactors = FALSE))
  }
  # preserve the order from selected_features.txt
  genes = proteins[proteins %in% genes]
  sub = as.data.frame(t(mat[genes, , drop = FALSE]))
  sub$sample = rownames(sub)
  sub$group = group[sub$sample]
  long = tidyr::pivot_longer(sub, cols = -c(sample, group), names_to = "protein", values_to = "score")
  long$protein = factor(long$protein, levels = genes)
  long$cohort = cohort_data$cohort
  long = as.data.frame(long)
  rownames(long) = NULL
  return(long)}

# Pairwise Wilcoxon tests (control vs alpha, control vs beta, alpha vs beta)
# for every cohort x set (or cohort x set x protein) combination.
pairwise_wilcox_tests = function(long_df, by){
  comparisons = list(c("control", "alpha"), c("control", "beta"), c("alpha", "beta"))
  long_df = as.data.frame(long_df)
  keys = as.data.frame(unique(long_df[, by, drop = FALSE]), stringsAsFactors = FALSE)

  out = list()
  for (i in seq_len(nrow(keys))){
    sub = long_df
    # use [[col]][i] (scalar); keys[i, col] can return a 1-column data.frame/tibble
    for (col in by) sub = sub[as.character(sub[[col]]) == as.character(keys[[col]][i]), , drop = FALSE]
    for (comp in comparisons){
      x = sub$score[as.character(sub$group) == comp[1]]
      y = sub$score[as.character(sub$group) == comp[2]]
      pval = tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
      row = as.data.frame(keys[i, , drop = FALSE], stringsAsFactors = FALSE)
      row$comparison = paste(comp, collapse = "_vs_")
      row$n1 = length(x); row$n2 = length(y)
      row$p.value = pval
      out[[length(out) + 1]] = row
    }
  }
  res = do.call(rbind, out)
  rownames(res) = NULL
  return(res)}

# Default ggplot boxplot (median/IQR + Tukey whiskers) + jittered points +
# pairwise Wilcoxon brackets, faceted by set (and by cohort if facet_by_cohort).
plot_vs_control = function(long_df, group_colors, main_title,
                           facet_by_cohort = TRUE, set_levels){
  long_df$set = factor(long_df$set, levels = set_levels)

  p = ggplot(long_df, aes(x = group, y = score, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.9, linewidth = 0.4) +
    geom_jitter(width = 0.15, size = 0.7, alpha = 0.5, colour = "grey30") +
    stat_compare_means(comparisons = list(c("control", "alpha"), c("control", "beta"), c("alpha", "beta")),
                       method = "wilcox.test", label = "p.format", size = 2.8, tip.length = 0.01) +
    scale_fill_manual(values = group_colors) +
    labs(title = main_title, x = NULL, y = "Mean expression (z-scored)") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey95", colour = NA))

  if (facet_by_cohort){
    p = p + facet_grid(cohort ~ set)
  } else {
    p = p + facet_wrap(~ set, nrow = 1)
  }
  return(p)}
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--intensity", "-i"),
              type = "character", default = "02_Missing_Inspection/norm_imp_MinProb.rds",
              help = "normalised, imputed SummarizedExperiment (all patients: control + ALS), relative to each cohort folder."
  ),make_option(c("--clusters", "-c"),
                type = "character", default = "08_Clustering_als/cluster_assignments_2.csv",
                help = "cluster_assignments_2.csv (patid + kmeans_k=2 alpha/beta), relative to each cohort folder."
  ),make_option(c("--alpha_degs", "-a"),
                type = "character", default = "Discovery/04_Vis_Differential_expression_analysis_subclusters/all_important_protein_names_k2_alpha.csv",
                help = "Discovery alpha-up DEGs (k=2, FDR<=0.05), used to define the fixed alpha-up DEG set."
  ),make_option(c("--beta_degs", "-b"),
                type = "character", default = "Discovery/04_Vis_Differential_expression_analysis_subclusters/all_important_protein_names_k2_beta.csv",
                help = "Discovery beta-up DEGs (k=2, FDR<=0.05), used to define the fixed beta-up DEG set."
  ),make_option(c("--wgcna", "-w"),
                type = "character", default = "Discovery/10_WGCNA_subclusters/moduleLabels.rds",
                help = "Discovery WGCNA moduleLabels.rds, used to define the fixed turquoise/blue module gene sets."
  ),make_option(c("--features", "-f"),
                type = "character", default = "15_ML_multi_models_500_04_noAge_3sd/selected_features.txt",
                help = "ML-selected protein list for per-protein boxplots (one protein name per line)."
  ),make_option(c("--output", "-o"),
                type = "character", default = "26_vsCtr_proteins",
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

set_levels = c("alpha-up DEGs", "beta-up DEGs", "turquoise module", "blue module")
set_category = c("alpha-up DEGs" = "inflammation-related",
                 "beta-up DEGs" = "synapto-axonal",
                 "turquoise module" = "inflammation-related",
                 "blue module" = "synapto-axonal")
group_colors = c(control = "grey70",
                 alpha = RColorBrewer::brewer.pal(8, "Set2")[2],
                 beta = RColorBrewer::brewer.pal(8, "Set2")[3])

################################################################################################
#DEFINE FIXED PROTEIN SETS (FROM DISCOVERY) + ML PANEL FOR PER-PROTEIN PLOTS
gene_sets = get_signature_gene_sets(opt$alpha_degs, opt$beta_degs, opt$wgcna)

ml_panel = gsub('^"|"$', "", trimws(readLines(opt$features)))
ml_panel = ml_panel[nzchar(ml_panel)]
if (length(ml_panel) == 0){
  stop("No proteins found in --features file: ", opt$features)
}
message("ML panel proteins for per-protein boxplots (", length(ml_panel), "): ",
        paste(ml_panel, collapse = ", "))

sink(paste0(output_dir, "/signature_gene_sets.txt"))
cat("Fixed protein sets, defined from the Discovery cohort:\n\n", sep = "")
for (set_name in set_levels){
  cat(set_name, " (", set_category[[set_name]], "), n = ", length(gene_sets[[set_name]]), ":\n", sep = "")
  cat(paste(gene_sets[[set_name]], collapse = ", "), "\n\n")
}
cat("ML-selected panel for per-protein boxplots (from ", opt$features, "):\n", sep = "")
cat(paste(ml_panel, collapse = ", "), "\n")
sink()

################################################################################################
#PER-COHORT DATA, SET SCORES, AND PROTEIN-LEVEL SCORES
cohorts = c("Discovery", "Validation")
set_long_all = list()
ml_protein_long_all = list()

for (cohort in cohorts){
  cohort_data = load_cohort_data(cohort, opt$intensity, opt$clusters)
  built = build_set_scores(cohort_data, gene_sets)
  set_long_all[[cohort]] = built$set_long
  ml_protein_long_all[[cohort]] = build_panel_protein_long(cohort_data, ml_panel)

  # wide/long per-sample expression for the ML panel proteins (log2 normalised)
  panel_genes = intersect(ml_panel, rownames(cohort_data$mat))
  panel_genes = ml_panel[ml_panel %in% panel_genes]
  protein_wide = as.data.frame(t(cohort_data$mat[panel_genes, , drop = FALSE]))
  protein_wide$sample = rownames(protein_wide)
  protein_wide$group = as.character(cohort_data$group[protein_wide$sample])
  protein_wide = protein_wide[, c("sample", "group", panel_genes)]
  write.csv(protein_wide, paste0(output_dir, "/", cohort, "_protein_expression.csv"), row.names = FALSE)
  write.csv(ml_protein_long_all[[cohort]], paste0(output_dir, "/", cohort, "_protein_expression_long.csv"), row.names = FALSE)

  set_wide = tidyr::pivot_wider(set_long_all[[cohort]], names_from = set, values_from = score)
  write.csv(set_wide, paste0(output_dir, "/", cohort, "_set_expression.csv"), row.names = FALSE)
  write.csv(set_long_all[[cohort]], paste0(output_dir, "/", cohort, "_set_expression_long.csv"), row.names = FALSE)
}

set_long = do.call(rbind, set_long_all)
ml_protein_long = do.call(rbind, ml_protein_long_all)
set_long$cohort = factor(set_long$cohort, levels = cohorts)
ml_protein_long$cohort = factor(ml_protein_long$cohort, levels = cohorts)

################################################################################################
#STATISTICAL TESTS
group_comparison_tests = pairwise_wilcox_tests(set_long, by = c("cohort", "set"))
write.csv(group_comparison_tests, paste0(output_dir, "/group_comparison_tests.csv"), row.names = FALSE)

protein_comparison_tests = pairwise_wilcox_tests(ml_protein_long, by = c("cohort", "protein"))
write.csv(protein_comparison_tests, paste0(output_dir, "/protein_comparison_tests.csv"), row.names = FALSE)

################################################################################################
#PLOTS
# combined figure: Discovery on top, Validation on the bottom
combined_plot = plot_vs_control(set_long, group_colors,
                                main_title = "Subtype-defining proteins vs control",
                                facet_by_cohort = TRUE, set_levels = set_levels)
ggsave(paste0(output_dir, "/Combined_subtype_vs_control_boxplot.pdf"), combined_plot,
      width = 11, height = 7, dpi = 300)

# one boxplot figure per cohort (single row, 4 panels), matching each cohort separately
for (cohort in cohorts){
  cohort_long = set_long[set_long$cohort == cohort, ]
  p = plot_vs_control(cohort_long, group_colors,
                      main_title = paste0(cohort, " - subtype-defining proteins vs control"),
                      facet_by_cohort = FALSE, set_levels = set_levels)
  ggsave(paste0(output_dir, "/", cohort, "_subtype_vs_control_boxplot.pdf"), p, width = 11, height = 3.2, dpi = 300)
}

# per-protein boxplots for the ML-selected panel: one row, one panel per protein
for (cohort in cohorts){
  cohort_protein_long = ml_protein_long[ml_protein_long$cohort == cohort, ]
  if (nrow(cohort_protein_long) == 0){
    message(cohort, ": no ML panel proteins available; skipping protein boxplot.")
    next
  }
  n_proteins = nlevels(droplevels(cohort_protein_long$protein))

  p = ggplot(cohort_protein_long, aes(x = group, y = score, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.6) +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.4, colour = "grey30") +
    stat_compare_means(comparisons = list(c("control", "alpha"), c("control", "beta"), c("alpha", "beta")),
                       method = "wilcox.test", label = "p.format", size = 2.5, tip.length = 0.01, hide.ns = FALSE) +
    facet_wrap(~ protein, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = group_colors) +
    labs(title = paste0(cohort, " – selected protein expression vs control"),
         x = NULL, y = "Expression (log2, normalized)") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          strip.text = element_text(size = 10, face = "bold"),
          strip.background = element_rect(fill = "grey90", colour = NA))

  ggsave(paste0(output_dir, "/", cohort, "_protein_boxplot.pdf"), p,
         width = max(10, 2.4 * n_proteins), height = 4, dpi = 300, limitsize = FALSE)
}
