#=========================Script Description=================================
# Heatmaps of protein detection patterns across samples in the Discovery and
# Validation cohorts, before and after filtering out low-abundance proteins
# with >50% missing values across all samples (DEP::filter_proteins).
#
# Columns represent samples, annotated by Cohort (Discovery/Validation) and
# Condition (als/ctrl). Rows represent proteins, grouped by membership
# (Shared: detected in both cohorts' protein panels; Discovery only;
# Validation only). Black = detected, grey = missing (either a missing value
# within a cohort's panel, or a protein never measured in that cohort).
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("ComplexHeatmap"))
suppressMessages(library("circlize"))
#===========================Function Definition=============================
# Build a Discovery-vs-Validation protein-detection heatmap for one panel
# (e.g. before or after the completeness filter).
# se_disc / se_val: SummarizedExperiment objects (unimputed; missing values
# are NA in assay()), one per cohort.
build_missingness_heatmap = function(se_disc, se_val, panel_title){

  det_disc = !is.na(assay(se_disc))
  det_val  = !is.na(assay(se_val))

  # order samples within each cohort: ctrl first, als last (matches
  # convention used in 07_Clinical_heatmap.R)
  disc_order = colnames(se_disc)[order(colData(se_disc)$condition == "als")]
  val_order  = colnames(se_val)[order(colData(se_val)$condition == "als")]

  det_disc = det_disc[, disc_order]
  det_val  = det_val[, val_order]

  # protein membership: Shared / Discovery only / Validation only
  all_proteins = union(rownames(det_disc), rownames(det_val))
  membership = ifelse(all_proteins %in% rownames(det_disc) & all_proteins %in% rownames(det_val),
                      "Shared",
                      ifelse(all_proteins %in% rownames(det_disc), "Discovery only", "Validation only"))
  membership = factor(membership, levels = c("Shared", "Discovery only", "Validation only"))
  row_ord = order(membership)
  all_proteins = all_proteins[row_ord]
  membership = membership[row_ord]

  n_disc = length(disc_order)
  n_val = length(val_order)

  # Discovery and Validation are independent cohorts and may reuse the same
  # sample/tube IDs (e.g. both starting at "1", "2", ...). Give columns
  # cohort-prefixed, guaranteed-unique labels so name-based column lookups
  # can never accidentally collide between the two cohorts.
  disc_labels = paste0("Discovery_", disc_order)
  val_labels  = paste0("Validation_", val_order)
  all_samples = c(disc_labels, val_labels)

  # build the combined detection matrix; proteins not measured in a cohort
  # stay NA for that cohort's columns and are treated as "Missing".
  # Rows are matched by protein name (intentional: this is how Shared /
  # Discovery only / Validation only membership is defined). Columns are
  # filled by position (Discovery block first, Validation block second) to
  # avoid any risk of column-name collisions between the two cohorts.
  mat = matrix(NA_character_, nrow = length(all_proteins), ncol = length(all_samples),
              dimnames = list(all_proteins, all_samples))

  row_idx_disc = match(rownames(det_disc), all_proteins)
  mat[row_idx_disc, seq_len(n_disc)] = ifelse(det_disc, "Detected", "Missing")

  row_idx_val = match(rownames(det_val), all_proteins)
  mat[row_idx_val, n_disc + seq_len(n_val)] = ifelse(det_val, "Detected", "Missing")

  mat[is.na(mat)] = "Missing"

  # Binary detection matrix (1 = Detected, 0 = Missing) for within-group
  # row ordering. Prefer this over raw expression: cohort-only proteins
  # and incomplete measurements leave large NA blocks that make dist()/hclust
  # fail with "NA/NaN/Inf in foreign function call".
  det_mat = matrix(as.numeric(mat == "Detected"),
                   nrow = nrow(mat), ncol = ncol(mat),
                   dimnames = dimnames(mat))

  # Order proteins within each membership group (group order fixed:
  # Shared -> Discovery only -> Validation only).
  # Sort by detection rate (blacker first), not Euclidean hclust on the full
  # 0/1 pattern: hclust groups similar *which-sample* patterns, so proteins
  # with similar overall blackness but different missing samples stay scattered.
  within_group_order = unlist(lapply(levels(membership), function(g) {
    idx = which(membership == g)
    if (length(idx) <= 1) return(idx)
    sub = det_mat[idx, , drop = FALSE]
    disc_rate = rowMeans(sub[, seq_len(n_disc), drop = FALSE])
    val_rate  = rowMeans(sub[, n_disc + seq_len(n_val), drop = FALSE])
    # overall rate first so black rows gather; then cohort rates as tie-breakers
    idx[order(-(disc_rate + val_rate), -disc_rate, -val_rate)]
  }), use.names = FALSE)
  mat = mat[within_group_order, , drop = FALSE]
  membership = membership[within_group_order]
  # left annotation must follow the same row order
  # (built below from membership)

  cohort = factor(c(rep("Discovery", n_disc), rep("Validation", n_val)),
                  levels = c("Discovery", "Validation"))
  condition = factor(c(as.character(colData(se_disc)[disc_order, "condition"]),
                       as.character(colData(se_val)[val_order, "condition"])),
                     levels = c("ctrl", "als"))

  top_anno = HeatmapAnnotation(
    Cohort = cohort,
    Condition = condition,
    col = list(Cohort = c(Discovery = "salmon", Validation = "#26b3ff"),
              Condition = c(als = "#D73027", ctrl = "#4575B4")),
    annotation_name_side = "left",
    simple_anno_size = unit(4, "mm")
  )

  left_anno = rowAnnotation(
    Protein = membership,
    col = list(Protein = c("Shared" = "mediumpurple1",
                          "Discovery only" = "salmon",
                          "Validation only" = "#26b3ff")),
    annotation_name_side = "top",
    simple_anno_size = unit(4, "mm")
  )

  ht = Heatmap(mat,
              name = "Detected",
              col = c("Detected" = "black", "Missing" = "grey80"),
              column_split = cohort,
              column_title = "%s",
              column_gap = unit(2, "mm"),
              cluster_columns = FALSE,
              show_column_names = FALSE,
              row_split = membership,
              row_title = NULL,
              cluster_rows = FALSE,
              cluster_row_slices = FALSE,
              show_row_dend = FALSE,
              row_gap = unit(1, "mm"),
              show_row_names = FALSE,
              left_annotation = left_anno,
              top_annotation = top_anno,
              border = TRUE)

  return(list(heatmap = ht,
              panel_title = panel_title,
              membership = data.frame(Protein = rownames(mat),
                                      Membership = as.character(membership),
                                      stringsAsFactors = FALSE)))}
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--before", "-b"),
              type = "character", default = "01_Pre_Processing/se_abu_data.rds",
              help = "unfiltered SummarizedExperiment object, before completeness filter."
  ),make_option(c("--after", "-a"),
                type = "character", default = "01_Pre_Processing/se_abu_data_filtered.rds",
                help = "filtered SummarizedExperiment object, after completeness filter."
  ),make_option(c("--output", "-o"),
                type = "character", default = "23_Heatmap_missingness",
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

if (is.null(opt$before) | is.null(opt$after)) {
  stop("Please provide both the before- and after-completeness-filter SummarizedExperiment object file paths!")
}else{
  res_Discovery_before = readRDS(paste0("Discovery/",opt$before))
  res_Validation_before = readRDS(paste0("Validation/",opt$before))
  res_Discovery_after = readRDS(paste0("Discovery/",opt$after))
  res_Validation_after = readRDS(paste0("Validation/",opt$after))
}

if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
}else{
  seed = opt$seed
  set.seed(seed)
}

################################################################################################
#MAKE MISSINGNESS HEATMAPS
before = build_missingness_heatmap(res_Discovery_before, res_Validation_before, "Before completeness filter")
after  = build_missingness_heatmap(res_Discovery_after,  res_Validation_after,  "After completeness filter")

pdf(paste0(output_dir, "/Missingness_heatmap_before.pdf"), width = 11, height = 8)
ComplexHeatmap::draw(before$heatmap,
                     column_title = before$panel_title,
                     column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
                     heatmap_legend_side = "right",
                     annotation_legend_side = "right")
dev.off()

pdf(paste0(output_dir, "/Missingness_heatmap_after.pdf"), width = 11, height = 8)
ComplexHeatmap::draw(after$heatmap,
                     column_title = after$panel_title,
                     column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
                     heatmap_legend_side = "right",
                     annotation_legend_side = "right")
dev.off()

# Export protein membership tables (per-protein and counts)
membership_label_short = c("Shared" = "Shared",
                           "Discovery only" = "Discovery",
                           "Validation only" = "Validation")

export_membership = function(mem_df, stage, output_dir){
  protein_membership = data.frame(
    Stage = stage,
    Protein = mem_df$Protein,
    Membership = mem_df$Membership,
    stringsAsFactors = FALSE
  )
  protein_membership = protein_membership[order(protein_membership$Protein), ]
  rownames(protein_membership) = NULL

  counts = as.data.frame(table(Membership = protein_membership$Membership),
                         stringsAsFactors = FALSE)
  colnames(counts) = c("Membership", "Freq")
  counts$Membership = as.character(counts$Membership)
  counts$Membership = ifelse(counts$Membership %in% names(membership_label_short),
                             membership_label_short[counts$Membership],
                             counts$Membership)
  counts = data.frame(Stage = stage,
                      Membership = counts$Membership,
                      Freq = counts$Freq,
                      stringsAsFactors = FALSE)
  # keep Shared / Discovery / Validation order
  counts$Membership = factor(counts$Membership, levels = c("Shared", "Discovery", "Validation"))
  counts = counts[order(counts$Membership), ]
  counts$Membership = as.character(counts$Membership)
  rownames(counts) = NULL

  write.csv(protein_membership,
            file = paste0(output_dir, "/protein_membership_", stage, ".csv"),
            row.names = FALSE)
  write.csv(counts,
            file = paste0(output_dir, "/protein_membership_counts_", stage, ".csv"),
            row.names = FALSE)
}

export_membership(before$membership, "before", output_dir)
export_membership(after$membership, "after", output_dir)
