#=========================Script Description=================================
# Revised version of 13_GSEA_IC_heatmap.R.
#
# We thank the reviewer for this observation. We agree that the previous
# heatmap could be somewhat misleading. We have revised Fig. 1f and 2h
# because, in the earlier version, non-significant pathways were treated as
# missing and displayed as white. As a result, cohorts could appear strongly
# discordant when a pathway was significant in one cohort but below the FDR
# threshold in another. In the updated heatmap, we display the actual signed
# -log10(FDR) for all included pathways across all cohorts. Light color
# cells now reflect weak or non-significant effects, rather than missing
# values. This revision reduces the apparent cross-cohort differences, while
# still preserving genuine biological discordance where it exists.
#
# Significantly enriched terms (FDR < 0.05 in at least one cohort) are kept,
# grouped by GO semantic similarity, and consolidated into --cluster_number
# categories. Enrichment is shown as signed -log10(FDR): red = enriched in
# the alpha subtype, blue = enriched in the beta subtype (sign of NES).
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("msigdbr"))
suppressMessages(library("circlize"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))
suppressMessages(library("ggsci"))
suppressMessages(library("tidyverse"))
suppressMessages(library("GOSemSim"))
suppressMessages(library("dendextend"))
suppressMessages(library("ComplexHeatmap"))
#===========================Functions=================================
parse_named_paths <- function(input_string) {
  parts <- strsplit(input_string, ",")[[1]]
  split_parts <- strsplit(parts, "=")
  names_vec <- sapply(split_parts, `[[`, 1)
  paths_vec <- sapply(split_parts, `[[`, 2)
  named_paths <- setNames(paths_vec, names_vec)
  return(named_paths)
}

# read the FULL (unfiltered) GSEA result table per cohort -- 06_GSEA.R is
# run with pvalueCutoff = 1, so every tested term/NES/FDR is available here,
# not just the significant ones.
read_gsea_results <- function(path_vector, ontology = "BP") {
  ego_list <- lapply(path_vector, function(p) {
    gsea <- readRDS(p)
    gsea[[ontology]]@result
  })
  names(ego_list) <- names(path_vector)
  return(ego_list)
}

# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "Discovery=Discovery/06_GSEA/GSEA_result.rds,Validation=Validation/06_GSEA/GSEA_result.rds",
              help = "Named GSEA result RDS files, format: name1=path1,name2=path2,..."),
  make_option(c("--output", "-o"),
              type = "character", default = "27_GSEA_IC_heatmap_replot",
              help = "Output directory"
  ),make_option(c("--cluster_number", "-c"),
                type = "integer", default = 6,
                help = "Number of clusters"
  ),make_option(c("--width", "-W"),
                type = "double", default = NULL,
                help = "Plot width in inches (default: auto from number of GO terms)"
  ),make_option(c("--height", "-H"),
                type = "double", default = NULL,
                help = "Plot height in inches (default: auto from number of cohorts)"
  ),make_option(c("--seed", "-e"),
                type = "integer", default = 9,
                help = "set.seed"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
#============================================================================
if (is.null(opt$input)) {
  stop("Please provide the GSEA results!")
}else{
  gsea_paths <- parse_named_paths(opt$input)
  ego_list <- read_gsea_results(gsea_paths, ontology = "BP")
}

if (is.null(opt$output)) {
  print("NO OUTPUT PATH SUPPLIED,current directory will be used!")
  output_dir <- getwd()
} else {
  output_dir <- opt$output
  if (!file.exists(output_dir)) {
    dir.create(output_dir, recursive = T)
  }
}

if (is.null(opt$cluster_number)) {
  k_num <- 6
} else {
  k_num <- opt$cluster_number
}

if (is.null(opt$seed)) {
  set.seed(9)
} else {
  set.seed(opt$seed)
}

#============================================================================
# Set the cutoff for adjusted p-value (used only to decide WHICH terms are
# included -- not to blank out cells for cohorts where a term misses this
# threshold)
padj_cutoff <- 0.05

# keep only terms significant in at least one cohort
all_go_terms <- unique(unlist(
  lapply(ego_list, function(x) {
    sig <- x
    sig <- sig[sig$p.adjust < padj_cutoff, ]
    sig$ID
  })
))

# build GO term x Comparison matrix
logp_mat <- matrix(NA_real_, nrow = length(all_go_terms), ncol = length(ego_list))
rownames(logp_mat) <- all_go_terms

# get background GO terms
bg_map <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  distinct(gs_name, gs_exact_source)  # gs_name is GOBP_xxx, gs_exact_source is GO:xxxxxxx

# merge all_go_terms with bg_map to get GO ID
all_go_terms = as.data.frame(all_go_terms)
colnames(all_go_terms) <- "Description"
all_go_terms <- left_join(all_go_terms, bg_map, by = c("Description" = "gs_name"))

colnames(logp_mat) <- names(ego_list)

# fill the matrix: signed -log10(FDR) for EVERY included term, in EVERY
# cohort, looked up from that cohort's full GSEA result table -- regardless
# of whether the term individually cleared the FDR cutoff in that cohort.
# This is the fix for the reviewer's comment: a term that is significant in
# one cohort but not the other is no longer forced to 0/white; its actual
# (weaker) signed -log10(FDR) is shown instead.
for (i in seq_along(ego_list)) {
  ego_res <- ego_list[[i]]
  idx <- match(rownames(logp_mat), ego_res$ID)

  if (any(is.na(idx))){
    missing_terms = rownames(logp_mat)[is.na(idx)]
    message(names(ego_list)[i], ": ", length(missing_terms),
           " selected term(s) not found in this cohort's GSEA result (left as NA): ",
           paste(missing_terms, collapse = ", "))
  }

  p_adj <- ego_res$p.adjust[idx]
  nes <- ego_res$NES[idx]
  logp_mat[, i] <- -log10(p_adj) * sign(nes)
}
rownames(logp_mat) <- all_go_terms$gs_exact_source

# get GO term semantic similarity matrix
go_sim <- godata(annoDb = "org.Hs.eg.db", ont = "BP", computeIC = TRUE)
go_terms <- rownames(logp_mat)
# get IC
ic_values <- go_sim@IC[go_terms]
# calculate semantic similarity
sim_matrix <- termSim(rownames(logp_mat), rownames(logp_mat), semData = go_sim, method = "Wang")
# transform similarity matrix to distance matrix
go_dist <- as.dist(1 - sim_matrix)

# Step 1: build cluster tree + colored tree branches
col_clust <- hclust(go_dist, method = "ward.D2")
col_dend <- as.dendrogram(col_clust)

# -------------------------------
# Step 2: build GO term description labels (replace default GO ID)
term_info <- ego_list[[1]]
term_desc_map <- setNames(term_info$Description, term_info$ID)

# only keep GO terms used in the heatmap (column names)
term_desc_map <- term_desc_map[rownames(logp_mat)]

clusters <- cutree(col_clust, k = k_num)

cluster_colors <- pal_d3("category20")(k_num)
names(cluster_colors) <- as.character(1:k_num)

col_dend <- color_branches(col_dend, k = k_num, groupLabels = TRUE)

top_annot <- HeatmapAnnotation(
  Cluster = as.factor(clusters),
  col = list(Cluster = cluster_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  annotation_legend_param = list(
    title = "Cluster",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
)

# -------------------------------
# Step 4: define color bar color function (same diverging scale as before:
# red = enriched in alpha subtype, blue = enriched in beta subtype)
col_fun <- colorRamp2(
   c(floor(min(logp_mat, na.rm = TRUE)), 0, ceiling(max(logp_mat, na.rm = TRUE))),
  c("#1B5B9D", "white", "#D7191C")
)

# -------------------------------
# Step 5: draw heatmap
ht <- Heatmap(
  t(logp_mat),
  name = "Signed -log10(adj. p-value)",
  top_annotation = top_annot,
  cluster_columns = col_dend,
  cluster_rows = FALSE,
  col = col_fun,
  na_col = "grey85",
  border = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 45,
  row_names_side = "left",
  column_dend_height = unit(2, "cm"),
  row_dend_width = unit(2, "cm"),
  column_title = "GO Term Clusters by Semantic Similarity",
  heatmap_legend_param = list(
    title = bquote("Signed " * -log[10]("adj. p-value")),
    at = seq(floor(min(logp_mat, na.rm = TRUE)), ceiling(max(logp_mat, na.rm = TRUE)), by = 2),
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 9)
  )
)

ht_drawn <- draw(ht)

# get the order of GO term names after clustering
col_ordered <- column_order(ht_drawn)
clustered_GO_terms <- rownames(logp_mat)[col_ordered]

# get the final cluster color order for the dendrogram
color_order = unique(clusters[clustered_GO_terms])
col_dend <- color_branches(col_dend, k = k_num, groupLabels = TRUE, col = cluster_colors[color_order])

ht <- ComplexHeatmap::Heatmap(
  t(logp_mat),
  name = "Signed -log10(adj. p-value)",
  top_annotation = top_annot,
  cluster_columns = col_dend,
  cluster_rows = FALSE,
  col = col_fun,
  na_col = "grey85",
  border = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 45,
  row_names_side = "left",
  column_dend_height = unit(2, "cm"),
  row_dend_width = unit(2, "cm"),
  column_title = "GO Term Clusters by Semantic Similarity",
  heatmap_legend_param = list(
    title = bquote("Signed " * -log[10]("adj. p-value")),
    at = seq(floor(min(logp_mat, na.rm = TRUE)), ceiling(max(logp_mat, na.rm = TRUE)), by = 2),
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 9)
  )
)

# save eps / pdf
plot_width = if (!is.null(opt$width)) opt$width else max(ceiling(dim(logp_mat)[1] * 0.06), 6)
plot_height = if (!is.null(opt$height)) opt$height else max(ceiling(dim(logp_mat)[2] * 0.5), 4)
message("Saving heatmap: width = ", plot_width, " in, height = ", plot_height, " in")

eps_path <- file.path(output_dir, "heatmap_plot.eps")
cairo_ps(
  filename = eps_path,
  width = plot_width,
  height = plot_height,
  family = "Arial",
  onefile = FALSE,
  fallback_resolution = 600
)
draw(ht)
dev.off()

pdf_path <- file.path(output_dir, "heatmap_plot.pdf")
cairo_pdf(
  filename = pdf_path,
  width = plot_width,
  height = plot_height,
  family = "Arial",
  fallback_resolution = 600
)
draw(ht)
dev.off()

# save the heatmap result to file
df_long <- as.data.frame(logp_mat) %>%
  tibble::rownames_to_column("GO_term") %>%
  pivot_longer(-GO_term, names_to = "Comparison", values_to = "log10_pvalue")

# add cluster information
df_long$Cluster <- clusters[df_long$GO_term]

df_long_filtered <- df_long %>%
  mutate(p_value = 10^(-abs(log10_pvalue)),
         IC = ic_values[GO_term])

# based on the min IC to represent the function of each cluster
go_representatives <- df_long_filtered %>%
  group_by(Cluster) %>%
  slice_min(IC, n = 1, with_ties = FALSE) %>%
  ungroup()

all_term_info <- lapply(ego_list, function(x) x[, c("ID", "Description")]) %>%
  bind_rows() %>%
  distinct(ID, .keep_all = TRUE) %>%
  left_join(all_go_terms, by = c("Description" = "Description")) %>%
  dplyr::select(gs_exact_source, Description) %>%
  rename(ID = gs_exact_source)

term_desc_full_map <- setNames(all_term_info$Description, all_term_info$ID)

df_long_filtered <- df_long_filtered %>%
  mutate(description = term_desc_full_map[GO_term])

go_representatives <- go_representatives %>%
  mutate(description = term_desc_full_map[GO_term])

df_long_filtered = df_long_filtered %>%
  left_join(
    go_representatives %>% dplyr::select(Cluster, Cluster_summary = description),
    by = "Cluster"
  )
df_long_filtered$description = gsub("GOBP_", "", df_long_filtered$description)
df_long_filtered$Cluster_summary = gsub("GOBP_", "", df_long_filtered$Cluster_summary)

output_csv_path <- file.path(output_dir, "GSEA_results_summary.csv")
write.csv(df_long_filtered, output_csv_path, row.names = FALSE)
