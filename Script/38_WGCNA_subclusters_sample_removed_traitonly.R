#=========================Script Description=================================
# "Trait-only" companion to Script/38_WGCNA_subclusters_sample_removed.R.
#
# This script does NOT re-run WGCNA (no pickSoftThreshold / blockwiseModules
# / GO enrichment). It reuses the module assignments (net$colors) already
# computed by 10_WGCNA_subclusters.R for a cohort (Discovery/10_WGCNA_subclusters
# or Validation/10_WGCNA_subclusters), removes the same fixed set of samples
# used elsewhere for that cohort:
#   Discovery:  CSF105, CSF68, CSF79, CSF96, CSF104, CSF109
#   Validation: CSF68, CSF79
# (see the header comment in Script/38_WGCNA_subclusters_sample_removed.R for
# where the Discovery IDs come from), and only recomputes the two
# trait-relationship-dependent outputs on the remaining samples:
#   - module eigengenes (MEs) recomputed on the sample-removed expression
#     matrix, using the ORIGINAL (full-cohort) module colors from net$colors
#   - Module-trait Pearson correlations/p-values -> WGCNA_ModuleTrait_heatmap_gradient.pdf
#     (signed -log10(p-value) diverging heatmap; same style/functions as
#     the ModuleTrait block in Script/11_WGCNA_comparison_subclusters.R)
#   - k2_boxplot.pdf (per-module eigenprotein boxplots by k2 subcluster;
#     same prepare_toplot()/calculate_pvalues()/plot_boxplots() functions as
#     Script/11_WGCNA_comparison_subclusters.R)
#
# Example:
# Rscript Script/38_WGCNA_subclusters_sample_removed_traitonly.R \
#   -i Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds \
#   -n Discovery/10_WGCNA_subclusters/WGCNA_net.rds \
#   -c Discovery/08_Clustering_als/cluster_assignments_2.csv \
#   -o Discovery/38_WGCNA_subclusters_sample_removed_traitonly \
#   --cohort Discovery -e 9
#
# Rscript Script/38_WGCNA_subclusters_sample_removed_traitonly.R \
#   -i Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds \
#   -n Validation/10_WGCNA_subclusters/WGCNA_net.rds \
#   -c Validation/08_Clustering_als/cluster_assignments_2.csv \
#   -o Validation/38_WGCNA_subclusters_sample_removed_traitonly \
#   --cohort Validation -e 9
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("tidyverse"))
suppressMessages(library("WGCNA"))
suppressMessages(library("ComplexHeatmap"))
suppressMessages(library("circlize"))
suppressMessages(library("grid"))
#===========================Function Definition=============================
# Same functions as Script/11_WGCNA_comparison_subclusters.R (relies on
# `cleanDat` and `MEs` being defined in the calling environment by the time
# these are invoked -- same convention used throughout the WGCNA scripts).
prepare_toplot <- function(MEs, cluster_file, metaData, cluster_col, ModuleTrait) {
  cluster_assignments <- read.csv(cluster_file)
  cluster_assignments <- cluster_assignments[match(colnames(cleanDat), cluster_assignments$patid), ]
  # remove kmeans_k columns from metaData
  metaData <- metaData %>% dplyr::select(-starts_with("kmeans_k"))
  factorMeta <- left_join(cluster_assignments, metaData, by = c("patid" = "label"))

  # grep the number from the cluster_col
  cluster_num <- as.numeric(gsub("\\D", "", cluster_col))
  colnames(factorMeta)[which(colnames(factorMeta) == paste0("kmeans_k.", cluster_num))] <- cluster_col

  factorMeta[[cluster_col]] <- as.factor(factorMeta[[cluster_col]])
  toplot <- t(MEs)
  toplot <- toplot[, factorMeta$patid]
  if (ModuleTrait) {
    factorMeta <- factorMeta[, c(cluster_col, "sex", "condition", "genetics")]
  } else {
    factorMeta <- factorMeta[, cluster_col, drop = FALSE]
  }

  list(toplot = toplot, factorMeta = factorMeta)
}

calculate_pvalues <- function(MEs, toplot, factorMeta, cluster_col) {
  group <- factorMeta[[cluster_col]]
  group_levels <- unique(group)
  num_groups <- length(group_levels)

  pvec <- numeric(ncol(MEs))

  for (i in seq_len(ncol(MEs))) {
    y <- MEs[, i]

    if (num_groups == 2) {
      # two group: Wilcoxon
      pval <- tryCatch({
        wilcox.test(y ~ group)$p.value
      }, error = function(e) NA)
    } else if (num_groups >= 3) {
      # multi group: Kruskal-Wallis test (non-parametric ANOVA)
      pval <- tryCatch({
        kruskal.test(y ~ group)$p.value
      }, error = function(e) NA)
    } else {
      pval <- NA
    }

    pvec[i] <- pval
  }

  names(pvec) <- colnames(MEs)
  pvec <- pvec[match(names(pvec), rownames(toplot))]
  rownames(toplot) <- paste(rownames(toplot), "\np = ", signif(pvec, 2), sep = "")
  return(toplot)
}

# Function to plot Eigenprotein value boxplots
plot_boxplots <- function(toplot, factorMeta, cluster_col, file_name) {
  pdf(file_name, onefile = FALSE)
  num_plots <- length(setdiff(seq_len(nrow(toplot)), which(colnames(MEs) == "grey")))
  layout_dim <- ceiling(sqrt(num_plots))
  par(mfrow = c(ceiling(num_plots / layout_dim), layout_dim))

  for (i in setdiff(seq_len(nrow(toplot)), which(colnames(MEs) == "grey"))) {
    boxplot(toplot[i, ] ~ factorMeta[[cluster_col]], col = colnames(MEs)[i],
            ylab = "Eigenprotein Value", main = rownames(toplot)[i], xlab = NULL)
  }
  dev.off()
  par(mfrow = c(1, 1))  # Reset plotting parameters
}
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds",
              help = "SummarizedExperiment object after normalize and imputation."
  ), make_option(c("--net", "-n"),
                 type = "character", default = "Discovery/10_WGCNA_subclusters/WGCNA_net.rds",
                 help = "WGCNA_net.rds produced by 10_WGCNA_subclusters.R for this SAME cohort. Module colors (net$colors) are reused as-is; the network itself is NOT re-run."
  ), make_option(c("--cluster_assignments", "-c"),
                 type = "character", default = "Discovery/08_Clustering_als/cluster_assignments_2.csv",
                 help = "cluster_assignments_2.csv"
  ), make_option(c("--output", "-o"),
                 type = "character", default = "38_WGCNA_subclusters_sample_removed_traitonly",
                 help = "output directory path."
  ), make_option(c("--cohort"),
                 type = "character", default = NULL,
                 help = "'Discovery', 'Validation', or 'External' -- selects the default --sample_removed list (Discovery: CSF105,CSF68,CSF79,CSF96,CSF104,CSF109; Validation: CSF68,CSF79; External: none). No default: must be set explicitly unless --sample_removed is passed instead. Ignored if --sample_removed is set explicitly."
  ), make_option(c("--sample_removed", "-r"),
                 type = "character", default = NULL,
                 help = "Comma-separated sample labels to remove. If unset, derived from --cohort. Pass \"\" (empty string) for no removal."
  ), make_option(c("--seed", "-e"),
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
    dir.create(output_dir, recursive = TRUE)
  }
}

if (is.null(opt$input)) {
  stop("Please provide the cleaned SummarizedExperiment object after normalize and imputation file path!")
} else if (!file.exists(opt$input)) {
  stop("SummarizedExperiment object after normalize and imputation does not exist!")
} else {
  cleanDat_file <- readRDS(opt$input)
}

if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
} else {
  seed <- opt$seed
  set.seed(seed)
}

if (is.null(opt$net)) {
  stop("Please provide the WGCNA_net.rds file path from 10_WGCNA_subclusters.R for this cohort!")
} else if (!file.exists(opt$net)) {
  stop("WGCNA network file does not exist!")
} else {
  net <- readRDS(opt$net)
}

if (is.null(opt$cluster_assignments)) {
  stop("Please provide the cluster_assignments file path!")
} else if (!file.exists(opt$cluster_assignments)) {
  stop("cluster_assignments file does not exist!")
} else {
  cluster_assignments_2_file <- opt$cluster_assignments
}

# Which samples to remove: explicit --sample_removed wins ("" = none);
# otherwise fall back to the fixed per-cohort lists. --cohort has NO
# default (see option_list above), so omitting both --cohort and
# --sample_removed is a hard stop, rather than silently defaulting to
# Discovery's removal list for some other cohort (e.g. External).
if (!is.null(opt$sample_removed)) {
  sample_removed <- if (!nzchar(trimws(opt$sample_removed))) {
    character(0)
  } else {
    trimws(strsplit(opt$sample_removed, ",")[[1]])
  }
} else if (identical(opt$cohort, "Discovery")) {
  sample_removed <- c("CSF105", "CSF68", "CSF79", "CSF96", "CSF104", "CSF109")
} else if (identical(opt$cohort, "Validation")) {
  sample_removed <- c("CSF68", "CSF79")
} else if (identical(opt$cohort, "External")) {
  sample_removed <- character(0)  # no fixed exclusion list defined for External
} else {
  stop("--cohort must be 'Discovery', 'Validation', or 'External' (or pass --sample_removed explicitly, e.g. --sample_removed \"\" for no removal).")
}
message("Removing samples: ", if (length(sample_removed) == 0) "(none)" else paste(sample_removed, collapse = ", "))

#################Main script##################
metaData <- as.data.frame(colData(cleanDat_file))
cleanDat <- assay(cleanDat_file)
colnames(cleanDat) <- metaData$label

# remove samples
if (length(sample_removed) > 0) {
  metaData <- metaData[!metaData$label %in% sample_removed, ]
  cleanDat <- cleanDat[, metaData$label]
}
rownames(metaData) <- metaData[, "label"]

# Module colors come from the ORIGINAL (full-cohort) net; restrict to
# proteins present in both as a safety check (same pattern as
# Script/11_WGCNA_comparison_subclusters.R's common_proteins step).
common_proteins <- intersect(names(net$colors), rownames(cleanDat))
if (length(common_proteins) < length(net$colors)) {
  message(length(net$colors) - length(common_proteins),
          " protein(s) in net$colors not found in the input matrix; they will be dropped.")
}
net_colors <- net$colors[common_proteins]
cleanDat <- cleanDat[common_proteins, ]

# Not every cohort's colData carries the same clinical trait columns (e.g.
# External only has label/ID/condition/replicate -- no age/Nfl/pNFh/etc.),
# so this is built from whichever of the usual trait columns are actually
# present, rather than assuming all of them exist. If NONE are present, the
# Module-Trait heatmap step below is skipped (k2_boxplot.pdf does not need
# numericMeta and is unaffected).
candidate_numericMeta_cols <- if (sum(!is.na(metaData$pNFh)) > 20) {
  c("age", "Nfl", "pNFh", "progression_rate", "slow_vital_capacity", "age_at_onset", "disease_duration")
} else {
  c("age", "Nfl", "progression_rate", "slow_vital_capacity", "age_at_onset", "disease_duration")
}
numericMeta_cols <- intersect(candidate_numericMeta_cols, colnames(metaData))
if (length(numericMeta_cols) < length(candidate_numericMeta_cols)) {
  message("Clinical trait column(s) not found in this cohort's metadata (skipped): ",
          paste(setdiff(candidate_numericMeta_cols, numericMeta_cols), collapse = ", "))
}
if (length(numericMeta_cols) > 0) {
  numericMeta <- metaData[, numericMeta_cols, drop = FALSE]
  numericMeta <- numericMeta[match(colnames(cleanDat), rownames(numericMeta)), , drop = FALSE]
} else {
  numericMeta <- NULL
  message("No clinical trait columns available for this cohort -- WGCNA_ModuleTrait_heatmap_gradient.pdf will be skipped (k2_boxplot.pdf is unaffected).")
}

# Recompute module eigengenes on the sample-removed matrix, using the
# ORIGINAL module colors (the network itself is not re-run).
MEList <- moduleEigengenes(t(cleanDat), colors = net_colors)
MEs <- MEList$eigengenes
MEs <- MEs[, colnames(MEs) != "MEgrey"]
colnames(MEs) <- substr(colnames(MEs), 3, 100)
rownames(MEs) <- colnames(cleanDat)

# Reorder module columns to match the module order in the ORIGINAL net
# (colnames(net$MEs), minus grey) rather than moduleEigengenes()'s own
# (alphabetical) ordering -- this is what determines row order in the
# Module-Trait heatmap and the module order in the k2 boxplot below.
net_module_order <- colnames(net$MEs)
net_module_order <- net_module_order[net_module_order != "MEgrey"]
net_module_order <- substr(net_module_order, 3, 100)
net_module_order <- intersect(net_module_order, colnames(MEs))  # drop any module absent after removal
MEs <- MEs[, net_module_order, drop = FALSE]

saveRDS(MEs, file.path(output_dir, "MEs_sample_removed.rds"))
if (!is.null(numericMeta)) {
  saveRDS(numericMeta, file.path(output_dir, "numericMeta_sample_removed.rds"))
}

####===========================Module-trait relationships (gradient heatmap)==================================
# Pearson correlation + p-value, same style as the ModuleTrait block in
# Script/11_WGCNA_comparison_subclusters.R (signed, symmetric diverging color scale).
# Skipped entirely if this cohort has no clinical trait columns at all
# (e.g. External -- see the numericMeta_cols block above).
if (!is.null(numericMeta)) {
  moduleTraitCor <- cor(MEs, numericMeta, use = "pairwise.complete.obs", method = "pearson")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(MEs))

  write.csv(moduleTraitCor, file.path(output_dir, "moduleTraitCor.csv"))
  write.csv(moduleTraitPvalue, file.path(output_dir, "moduleTraitPvalue.csv"))

  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) <- dim(moduleTraitCor)
  par(mfrow = c(1, 1))
  par(mar = c(6, 8.5, 3, 3))

  pdf(file.path(output_dir, "WGCNA_ModuleTrait_heatmap_gradient.pdf"), width = 6, height = 5.5)

  white_threshold <- -log10(0.1)
  # signed -log10(p-value)
  signed_logP <- sign(moduleTraitCor) * -log10(moduleTraitPvalue)

  # color scale fixed to -2..2 (same as Script/11_WGCNA_comparison_subclusters.R)
  max_abs_logP <- 2

  col_fun <- colorRamp2(
    c(-max_abs_logP, -white_threshold, 0, white_threshold, max_abs_logP),
    c("#4575B4", "white", "white", "white", "#D73027")
  )

  legend_breaks <- seq(-max_abs_logP, max_abs_logP, by = 1)

  ht <- Heatmap(
    matrix = signed_logP,
    col = col_fun,
    cell_fun = function(j, i, x, y, w, h, col) {
      grid.text(textMatrix[i, j], x, y, gp = gpar(fontsize = 8))
    },
    row_names_side = "left",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 12),
    column_names_gp = gpar(fontsize = 12),
    column_names_rot = 45,
    column_title = "Module-trait relationships\npearson r-value\n(p-value)",
    column_title_gp = gpar(fontsize = 12),
    border_gp = gpar(col = "black"),
    heatmap_legend_param = list(
      title = "signed -log10(p-value)",
      border = "black",
      at = legend_breaks,
      labels = as.character(legend_breaks)
    )
  )

  draw(ht)
  dev.off()
} else {
  message("Skipped WGCNA_ModuleTrait_heatmap_gradient.pdf (no clinical trait columns for this cohort).")
}

####===========================k2 eigenprotein boxplots==================================
data_k2 <- prepare_toplot(MEs, cluster_assignments_2_file, metaData, "k2", ModuleTrait = FALSE)
toplot_k2 <- calculate_pvalues(MEs, data_k2$toplot, data_k2$factorMeta, "k2")
plot_boxplots(toplot_k2, data_k2$factorMeta, "k2", file.path(output_dir, "k2_boxplot.pdf"))

message("Done. Outputs written to: ", output_dir)
