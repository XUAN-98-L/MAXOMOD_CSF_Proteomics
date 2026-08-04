#=========================Script Description=================================
# This script is used for GSVA analysis
# cd /Users/xliu2942/Documents/Projects/MAXOMOD/TESTing/Pipeline
# Rscript 6GSVA.R -s GO_BP -c 5_Clustering_als/cluster_assignments_2.csv >output.log
#
# Alignment with 06_GSEA.R (so GSVA settings are comparable to the GSEA
# results): when --category is the default "C5" and no --subset is given,
# this script now runs GO Biological Process, Molecular Function, and
# Cellular Component SEPARATELY (three full analyses, one per ontology),
# the same way 06_GSEA.R loops ont = c("BP","MF","CC"). Outputs are suffixed
# _BP / _MF / _CC accordingly. Passing an explicit --subset (or a non-C5
# --category) runs once, as before, with no suffix.
# --min_sz now defaults to 10 (matching 06_GSEA.R's minGSSize = 10) and
# --max_sz defaults to 500 (matching maxGSSize = 500); both are applied via
# gsva()'s own min.sz/max.sz rather than a separate hand-rolled pre-filter.
#
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("DEP"))
# suppressMessages(library("pheatmap"))
suppressMessages(library("ComplexHeatmap"))
suppressMessages(library("circlize"))
suppressMessages(library("grid"))
suppressMessages(library("GSVA"))
suppressMessages(library("tidyverse"))
suppressMessages(library("limma"))
suppressMessages(library("msigdbr"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("GOSemSim"))
suppressMessages(library("ggsci"))
# suppressMessages(library("ggplot2"))
# suppressMessages(library("ggthemes"))
# suppressMessages(library("writexl"))
#===========================Function Definition=============================
# #all colourblind friendly options
# display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE,
#                    colorblindFriendly=TRUE)
#
# #create a list with all colours for the final paper figures
# #each list element will be a named vector
# final_colours = list(
#   volcano_plot = brewer.pal(n=10, "RdYlBu")[c(2,9)],
#   male_female = brewer.pal(n=8, "Set2")[c(1,4)],
#   heatmap_scale = rev(brewer.pal(n=11, "RdBu")),
#   clustering = brewer.pal(n=8, "Set2")[c(2,3,5)],
#   age_scale = brewer.pal(n = 9, "YlGn"),
#   disease_progression_scale = brewer.pal(n = 9, "YlOrRd"),
#   onset = c(brewer.pal(n=11, "RdYlBu")[c(4,5)], "#B3B3B3"),
#   disease_status = c(brewer.pal(n=11, "RdYlBu")[2],"#B3B3B3"),
#   genetic_testing = c("#B3B3B3", brewer.pal(n = 11, "PRGn")[c(3, 9)]),
#   center = c("purple4", "orange3"),
#   neurofilaments = brewer.pal(n = 9, "PuBu"),
#   pNFh_scale = brewer.pal(n = 9, "Purples"),
#   age_at_onset_scale = brewer.pal(n = 9, "Blues"),
#   slow_vital_capacity_scale = brewer.pal(n = 9, "Reds")
# )
#
# names(final_colours$volcano_plot) =  c("up", "down")
# names(final_colours$male_female) = c("Male", "Female")
# names(final_colours$clustering) = c("alpha", "beta", "theta")
# names(final_colours$disease_status) = c("als", "ctrl")
# names(final_colours$onset) = c("spinal", "bulbar", "ctrl")
# names(final_colours$genetic_testing) = c("not_performed", "negative", "C9orf72")
# names(final_colours$center) = c("goettingen", "munich")


# # Define a function to process each comparison dataframe (Excel topn export; unused for mean_k3 outputs)
# process_comparison <- function(comparison_df, comparison_name, GO_list, filtered_mapped_matrix, topn, output_dir) {
#
#   comparison_topn_up <- comparison_df %>% dplyr::slice_max(n = (topn+20)/2,order_by = t)  %>% as.data.frame()
#   comparison_topn_down <- comparison_df %>% dplyr::slice_min(n = (topn+20)/2,order_by = t)  %>% as.data.frame()
#   comparison_topn_df = rbind(comparison_topn_up,comparison_topn_down)
#
#   common_GO_terms <- intersect(names(GO_list), rownames(comparison_topn_df))
#   GO_list_subset <- GO_list[names(GO_list) %in% common_GO_terms]
#   assay_genes <- rownames(filtered_mapped_matrix)
#   GO_list_subset_filtered <- lapply(GO_list_subset, function(genes) {
#     intersect(genes, assay_genes)
#   })
#   GO_list_top_terms <- GO_list_subset_filtered[names(GO_list_subset_filtered) %in% common_GO_terms]
#
#   comparison_topn <- comparison_topn_df %>%
#     rownames_to_column("terms") %>%
#     mutate(genes = map(terms, ~ GO_list_top_terms[[.]])) %>%
#     select(terms, genes, everything())
#   comparison_topn$genes <- sapply(comparison_topn$genes, function(x) paste(x, collapse = ","))
#
#   write_xlsx(comparison_topn %>% as.data.frame(),
#              path = file.path(output_dir, paste0("GSVA_topn_terms_", comparison_name, "_with_matched_proteins.xlsx")))
#
#   return(comparison_topn)
# }

# Combine genes for each term that are present in the assay matrix.
get_matched_proteins = function(terms, GO_list, assay_genes){
  vapply(terms, function(term){
    genes = GO_list[[term]]
    if (is.null(genes)) return(NA_character_)
    paste(intersect(genes, assay_genes), collapse = ",")
  }, character(1))}

# Annotate cluster-means table with limma contrast stats and matched proteins.
# contrast_list: named list of data.frames, each with columns terms + limma stats.
enrich_cluster_means_table = function(cluster_means_df, contrast_list, GO_list, assay_genes){
  out = cluster_means_df
  out$matched_proteins = get_matched_proteins(out$terms, GO_list, assay_genes)
  stat_cols = c("logFC", "CI.L", "CI.R", "AveExpr", "t", "P.Value", "adj.P.Val", "B")

  for (nm in names(contrast_list)){
    ct = contrast_list[[nm]]
    if (!is.data.frame(ct)) ct = as.data.frame(ct)
    if (!("terms" %in% colnames(ct))){
      ct = ct %>% tibble::rownames_to_column("terms")
    }
    keep = c("terms", intersect(stat_cols, colnames(ct)))
    stats = ct[, keep, drop = FALSE]
    rename_cols = setdiff(colnames(stats), "terms")
    colnames(stats)[colnames(stats) != "terms"] = paste0(nm, "_", rename_cols)
    out = dplyr::left_join(out, stats, by = "terms")
  }
  return(out)}

# Format MSigDB pathway IDs for display, e.g. GOBP_CELLULAR_RESPONSE_TO_LIPID -> Cellular Response To Lipid
format_pathway_label <- function(term, max_length = 60) {
  x = gsub("^(GOBP|GOCC|GOMF|HP|REACTOME|KEGG|HALLMARK)_", "", term)
  x = gsub("_", " ", x)
  x = stringr::str_to_title(tolower(x))
  if (nchar(x) > max_length) {
    paste0(substr(x, 1, max_length - 3), "...")
  } else {
    x
  }
}

# Custom function to abbreviate long terms by adding "..." at the end
abbreviate_terms_end <- function(term, max_length = 60) {
  format_pathway_label(term, max_length = max_length)
}

# Run one full GSVA + k=3 subtype-comparison analysis for a single gene set
# collection (e.g. GO:BP), writing suffixed outputs to output_dir. Everything
# that does NOT depend on the choice of gene set collection (the expression
# matrix, sample metadata, k=3 design/contrasts) is computed once by the
# caller and passed in, so this function is only re-run for the parts that
# do depend on the gene set collection.
run_gsva_geneset_analysis = function(gs_subcat, output_suffix, label_for_messages,
                                     GO_gene_sets_full, filtered_mapped_matrix,
                                     min_sz, max_sz, topn, n_cut,
                                     plot_width, plot_height,
                                     annolabel, group_info, design, contrast_matrix,
                                     output_dir){

  message(paste0("Running GSVA for gene set collection: ", label_for_messages))

  GO_gene_sets = if (is.null(gs_subcat)) GO_gene_sets_full else GO_gene_sets_full[GO_gene_sets_full$gs_subcat == gs_subcat, ]

  # Background removal: restrict to genes present in our data, using the
  # SAME gene identifier (rownames of filtered_mapped_matrix, i.e. the
  # SummarizedExperiment's own "name" column) that is actually fed into
  # gsva() below -- and that 06_GSEA.R uses as names(geneList). This keeps
  # the gene universe underlying GSVA's gene sets identical to GSEA's, and
  # keeps GO_list consistent with gsva()'s own min.sz/max.sz size filtering.
  GO_gene_sets = GO_gene_sets[GO_gene_sets$gene_symbol %in% rownames(filtered_mapped_matrix), ]

  # msigdbr()'s raw output can contain duplicate (gs_name, gene_symbol) rows
  # (it's built from a table with additional ID columns -- e.g. multiple
  # Entrez/Ensembl IDs mapping to the same gene_symbol -- that get dropped
  # by the gs_name/gene_symbol select() above without collapsing duplicates).
  # clusterProfiler::GSEA() effectively deduplicates genes within each
  # pathway before applying minGSSize/maxGSSize; gsva()'s min.sz/max.sz does
  # not, so without distinct() here a handful of pathways can look like they
  # have >=10 genes when they actually have <10 unique genes, inflating
  # GSVA's tested-pathway count relative to GSEA's (confirmed empirically:
  # 724 pathways pass a naive, non-deduplicated overlap count on either
  # script's background, but clusterProfiler::GSEA() -- and, with this fix,
  # gsva() -- only tests 670 once genes are deduplicated per pathway).
  GO_gene_sets = dplyr::distinct(GO_gene_sets, gs_name, gene_symbol)

  GO_list <- split(GO_gene_sets$gene_symbol, GO_gene_sets$gs_name)

  message("Performing GSVA analysis")

  # https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
  #https://www.nature.com/articles/s41593-021-01006-0
  # Gene set size filtering (min/max number of proteins overlapping our
  # data) is handled entirely by gsva()'s own min.sz/max.sz here, matching
  # how 06_GSEA.R relies on clusterProfiler::GSEA()'s minGSSize/maxGSSize
  # rather than a separate hand-rolled pre-filter.
  gsva_results <- gsva(
    filtered_mapped_matrix,
    GO_list,
    method = "ssgsea",
    # Appropriate for our vst transformed data
    kcdf = "Gaussian",
    # Minimum gene set size
    min.sz = min_sz,
    # Maximum gene set size
    max.sz = max_sz,
    # Compute Gaussian-distributed scores
    mx.diff = TRUE,
    # Don't print out the progress bar
    verbose = FALSE
  )

  saveRDS(gsva_results, file.path(output_dir, paste0("gsva_results", output_suffix, ".rds")))

  # Z-score normalize each row (gene set) across samples
  zscore_normalize <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  gsva_results <- t(apply(gsva_results, 1, zscore_normalize))

  ####for kmeans k = 3
  # Fit the linear model
  fit <- lmFit(gsva_results, design)

  # Apply contrasts to the fit
  fit2 <- contrasts.fit(fit, contrast_matrix)

  # Compute empirical Bayes statistics
  fit2 <- eBayes(fit2)

  # Full limma tables (all pathways) for cluster-means annotation; significant subset for CSVs/heatmaps
  theta_vs_beta_all <- topTable(fit2, coef = "theta_vs_beta", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE)
  alpha_vs_theta_all <- topTable(fit2, coef = "alpha_vs_theta", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE)
  alpha_vs_beta_all <- topTable(fit2, coef = "alpha_vs_beta", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE)
  theta_vs_beta <- topTable(fit2, coef = "theta_vs_beta", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE, p.value = 0.05)
  alpha_vs_theta <- topTable(fit2, coef = "alpha_vs_theta", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE, p.value = 0.05)
  alpha_vs_beta <- topTable(fit2, coef = "alpha_vs_beta", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE, p.value = 0.05)

  # Keep topn subsetting: these terms define GSVA_mean_k3.csv / GSVA_mean_k3.pdf
  theta_vs_beta_up = theta_vs_beta %>% dplyr::slice_max(n = topn/2,order_by = t)  %>% as.data.frame()
  theta_vs_beta_down = theta_vs_beta %>% dplyr::slice_min(n = topn/2,order_by = t)  %>% as.data.frame()
  theta_vs_beta = rbind(theta_vs_beta_up,theta_vs_beta_down)

  alpha_vs_theta_up = alpha_vs_theta %>% dplyr::slice_max(n = topn/2, order_by = t) %>% as.data.frame()
  alpha_vs_theta_down = alpha_vs_theta %>% dplyr::slice_min(n = topn/2, order_by = t) %>% as.data.frame()
  alpha_vs_theta = rbind(alpha_vs_theta_up, alpha_vs_theta_down)

  alpha_vs_beta_up = alpha_vs_beta %>% dplyr::slice_max(n = topn/2, order_by = t) %>% as.data.frame()
  alpha_vs_beta_down = alpha_vs_beta %>% dplyr::slice_min(n = topn/2, order_by = t) %>% as.data.frame()
  alpha_vs_beta = rbind(alpha_vs_beta_up, alpha_vs_beta_down)

  ### cluster means heatmap k = 3
  cluster_means <- gsva_results %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    mutate(Cluster = annolabel$Cluster) %>%
    group_by(Cluster) %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    column_to_rownames(var = "Cluster") %>%
    t()

  write_cluster_means = cluster_means %>% as.data.frame() %>% rownames_to_column("terms")
  write_cluster_means = enrich_cluster_means_table(
    write_cluster_means,
    contrast_list = list(
      alpha_vs_beta = alpha_vs_beta_all,
      alpha_vs_theta = alpha_vs_theta_all,
      theta_vs_beta = theta_vs_beta_all
    ),
    GO_list = GO_list,
    assay_genes = rownames(filtered_mapped_matrix)
  )
  write.csv(write_cluster_means, quote = TRUE, row.names = FALSE,
            file = file.path(output_dir, paste0("GSVA_cluster_means_k3", output_suffix, ".csv")))

  sig_terms_k3 = unique(c(rownames(theta_vs_beta), rownames(alpha_vs_theta), rownames(alpha_vs_beta)))
  if (length(sig_terms_k3) == 0){
    message(label_for_messages, ": no pathway reached adj.P < 0.05 in any k=3 contrast; skipping heatmap.")
    return(invisible(NULL))
  }
  write_mean_k3 = write_cluster_means[match(sig_terms_k3, write_cluster_means$terms), , drop = FALSE]

  cluster_means = cluster_means[sig_terms_k3, , drop = FALSE]
  # Keep subtype column order for heatmap
  cluster_means = as.matrix(cluster_means[, c("alpha", "beta", "theta"), drop = FALSE])

  # Star mark: pathway in a cluster is significant vs both other clusters
  # (both pairwise limma adj.P.Val < 0.05, consistent direction)
  star_mat = matrix("", nrow = nrow(cluster_means), ncol = ncol(cluster_means),
                    dimnames = dimnames(cluster_means))
  for (term in rownames(cluster_means)) {
    ab_p = alpha_vs_beta_all[term, "adj.P.Val"]
    at_p = alpha_vs_theta_all[term, "adj.P.Val"]
    tb_p = theta_vs_beta_all[term, "adj.P.Val"]
    ab_fc = alpha_vs_beta_all[term, "logFC"]
    at_fc = alpha_vs_theta_all[term, "logFC"]
    tb_fc = theta_vs_beta_all[term, "logFC"]

    # alpha vs beta & theta: same sign of (alpha-beta) and (alpha-theta)
    if (!anyNA(c(ab_p, at_p, ab_fc, at_fc)) &&
        ab_p < 0.05 && at_p < 0.05 &&
        sign(ab_fc) == sign(at_fc) && sign(ab_fc) != 0) {
      star_mat[term, "alpha"] = "*"
    }
    # beta vs alpha & theta: same sign of (alpha-beta) and (theta-beta)
    if (!anyNA(c(ab_p, tb_p, ab_fc, tb_fc)) &&
        ab_p < 0.05 && tb_p < 0.05 &&
        sign(ab_fc) == sign(tb_fc) && sign(ab_fc) != 0) {
      star_mat[term, "beta"] = "*"
    }
    # theta vs alpha & beta: (alpha-theta) and (theta-beta) opposite signs
    if (!anyNA(c(at_p, tb_p, at_fc, tb_fc)) &&
        at_p < 0.05 && tb_p < 0.05 &&
        sign(at_fc) == -sign(tb_fc) && sign(at_fc) != 0) {
      star_mat[term, "theta"] = "*"
    }
  }

  n_cut_use = min(n_cut, nrow(cluster_means))
  if (n_cut_use < n_cut) {
    message(paste0(label_for_messages, ": requested --number=", n_cut, " exceeds nrow=", nrow(cluster_means),
                   "; using ", n_cut_use))
  }

  row_hc = hclust(dist(cluster_means), method = "ward.D2")
  row_clusters = cutree(row_hc, k = n_cut_use)

  # GO IC-based cluster labeling (same approach as 13_GSEA_IC_heatmap.R):
  # only meaningful for GO ontologies (BP/MF/CC). ont_for_IC is derived
  # from gs_subcat ("GO:BP" -> "BP", etc.); NA for non-GO collections
  # (e.g. Hallmark, Reactome) or when the whole C5 category is used
  # unsplit, in which case we skip IC and fall back to dendrogram order
  # for the cluster representative label.
  ont_for_IC = if (!is.null(gs_subcat) && grepl("^GO:", gs_subcat)) sub("^GO:", "", gs_subcat) else NA

  if (!is.na(ont_for_IC)) {
    message("Computing GO IC values for pathway clusters (", ont_for_IC, ")")
    bg_map = msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = paste0("GO:", ont_for_IC)) %>%
      dplyr::distinct(gs_name, gs_exact_source)
    go_id_map = setNames(bg_map$gs_exact_source, bg_map$gs_name)
    go_ids = unname(go_id_map[rownames(cluster_means)])
    names(go_ids) = rownames(cluster_means)

    go_sim = GOSemSim::godata(annoDb = "org.Hs.eg.db", ont = ont_for_IC, computeIC = TRUE)
    ic_values = go_sim@IC[go_ids]
    names(ic_values) = rownames(cluster_means)
  } else {
    go_ids = setNames(rep(NA_character_, nrow(cluster_means)), rownames(cluster_means))
    ic_values = setNames(rep(NA_real_, nrow(cluster_means)), rownames(cluster_means))
  }

  pathway_cluster_df = data.frame(
    terms = rownames(cluster_means),
    GO_ID = unname(go_ids),
    IC = as.numeric(ic_values),
    Cluster = as.integer(row_clusters),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Lowest-IC term within each cluster as cluster representative / name
  # (falls back to first term in dendrogram order when IC is unavailable,
  # e.g. for non-GO gene set collections)
  go_representatives = pathway_cluster_df %>%
    dplyr::group_by(Cluster) %>%
    dplyr::slice_min(order_by = IC, n = 1, with_ties = FALSE, na_rm = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Cluster_summary = terms)

  # If a cluster had all NA IC, fall back to first term in dendrogram order
  missing_rep = setdiff(unique(pathway_cluster_df$Cluster), go_representatives$Cluster)
  if (length(missing_rep) > 0) {
    fallback = pathway_cluster_df %>%
      dplyr::filter(Cluster %in% missing_rep) %>%
      dplyr::group_by(Cluster) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Cluster_summary = terms)
    go_representatives = dplyr::bind_rows(go_representatives, fallback)
  }

  pathway_cluster_df = pathway_cluster_df %>%
    dplyr::left_join(
      go_representatives %>% dplyr::select(Cluster, Cluster_summary),
      by = "Cluster"
    )

  # Attach IC / cluster info to GSVA_mean_k3.csv
  write_mean_k3 = write_mean_k3 %>%
    dplyr::left_join(pathway_cluster_df, by = "terms")
  write.csv(write_mean_k3, quote = TRUE, row.names = FALSE,
            file = file.path(output_dir, paste0("GSVA_mean_k3", output_suffix, ".csv")))

  # Row annotation / legend labels = formatted lowest-IC pathway name
  cluster_summary_map = setNames(
    vapply(go_representatives$Cluster_summary, format_pathway_label, character(1)),
    as.character(go_representatives$Cluster)
  )
  level_ids = unique(row_clusters[row_hc$order])
  row_annot = factor(
    unname(cluster_summary_map[as.character(row_clusters)]),
    levels = unname(cluster_summary_map[as.character(level_ids)])
  )
  names(row_annot) = names(row_clusters)

  cluster_colors = ggsci::pal_d3("category20")(n_cut_use)
  names(cluster_colors) = levels(row_annot)

  left_annot = rowAnnotation(
    Cluster = row_annot,
    col = list(Cluster = cluster_colors),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      title = "Cluster (min IC)",
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 7)
    ),
    width = unit(3, "mm")
  )

  col_fun = colorRamp2(
    c(min(cluster_means, na.rm = TRUE), 0, max(cluster_means, na.rm = TRUE)),
    c("navy", "white", "firebrick3")
  )

  # dendrogram + numeric row_split keeps the clustering tree; titles hidden
  pathway_heatmap_mean = Heatmap(
    cluster_means,
    name = "mean z",
    col = col_fun,
    cluster_rows = as.dendrogram(row_hc),
    row_split = n_cut_use,
    row_gap = unit(2.5, "mm"),
    row_dend_width = unit(2, "cm"),
    cluster_columns = FALSE,
    left_annotation = left_annot,
    show_row_names = TRUE,
    row_labels = sapply(rownames(cluster_means), abbreviate_terms_end),
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 10),
    row_title = NULL,
    border = TRUE,
    heatmap_legend_param = list(title = "mean GSVA z"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (star_mat[i, j] != "") {
        grid.text(star_mat[i, j], x, y, gp = gpar(fontsize = 8, col = "black"))
      }
    }
  )

  pdf_height = if (is.null(plot_height) || is.na(plot_height)) {
    max(6, nrow(cluster_means) * 0.12)
  } else {
    plot_height
  }
  pdf(file.path(output_dir, paste0("GSVA_mean_k3", output_suffix, ".pdf")), width = plot_width, height = pdf_height)
  draw(pathway_heatmap_mean,
       annotation_legend_list = list(
         Legend(title = "limma",
                labels = "* sig. vs other two\n  (adj.P < 0.05)",
                type = "points", pch = 8,
                background = "white")
       ))
  dev.off()

  invisible(NULL)
}
#===========================Command Parameters Setting=========================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "2_Missing_Inspection_als/norm_imp_MinProb.rds",
              help = "SummarizedExperiment object after normalize and imputation."
  ),make_option(c("--output", "-o"),
                type = "character", default = "6_GSVA",
                help = "output directory path."
  ),make_option(c("--seed", "-e"),
                type = "integer", default = 9,
                help = "set.seed"
  ),make_option(c("--cluster_assignments", "-c"),
                type = "character", default = "5_Clustering/cluster_assignments_2.csv",
                help = "cluster_assignments_2.csv"
  ),make_option(c("--uniprot_to_genename", "-u"),
                type = "character", default = "1_pre_processing/uniprot_to_genename.rds",
                help = "uniprot_to_genename.rds"
  ),make_option(c("--min_sz", "-m"),
                type = "integer", default = 10,
                help = "Minimum gene set size (proteins overlapping our data), passed to gsva()'s min.sz. Default 10, matching 06_GSEA.R's minGSSize."
  ),make_option(c("--max_sz", "-x"),
                type = "integer", default = 500,
                help = "Maximum gene set size (proteins overlapping our data), passed to gsva()'s max.sz. Default 500, matching 06_GSEA.R's maxGSSize."
  ),make_option(c("--category", "-g"),
                type = "character", default = "C5",
                help = "MSigDB category, default is C5 (GO terms)"
  ),make_option(c("--subset", "-s"),
                type = "character",default = NULL,
                help = "subset of MSigDB category. Default NULL: if --category is C5, runs GO Biological Process, Molecular Function, and Cellular Component SEPARATELY (matching 06_GSEA.R), with outputs suffixed _BP/_MF/_CC. Pass e.g. 'GO_BP' to run a single subset instead (output has no suffix)."),make_option(c("--topn", "-t"),
                type = "integer", default = 30,
                help = "top n pathways to be displayed in the heatmap"
  ),make_option(c("--number", "-n"),
                type = "integer", default = 3,
                help = "number of pathway clusters (cutree) for GSVA_mean_k3 heatmap row split"
  ),make_option(c("--width", "-W"),
                type = "double", default = 7,
                help = "GSVA_mean_k3.pdf width in inches (default: 7)"
  ),make_option(c("--height", "-H"),
                type = "double", default = NULL,
                help = "GSVA_mean_k3.pdf height in inches (default: auto from number of pathways)")
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
  stop("Please provide the cleaned SummarizedExperiment object after normalize and imputation file path!")
}else if (!file.exists(opt$input)) {
  stop("SummarizedExperiment object after normalize and imputation does not exist!")
}else{
  input = opt$input
  se = readRDS(input)
}

if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
}else{
  seed = opt$seed
  set.seed(seed)
}

# NOTE: --uniprot_to_genename is validated for backward CLI compatibility
# (existing calls that pass -u keep working) but the loaded object is no
# longer used: the GSVA matrix is now built from se's own "name" rownames
# directly, to match 06_GSEA.R's gene universe (see header comment).
if (is.null(opt$uniprot_to_genename)) {
  stop("Please provide the uniprot_to_genename file path!")
}else if (!file.exists(opt$uniprot_to_genename)) {
  stop("uniprot_to_genename file does not exist!")
}else{
  uniprot_to_genename = readRDS(opt$uniprot_to_genename)
}

if (is.null(opt$cluster_assignments)) {
  stop("Please provide the cluster_assignments.rds file path!")
}else if (!file.exists(opt$cluster_assignments)) {
  stop("cluster_assignments file does not exist!")
}else{
  cluster_assignments = read.csv(opt$cluster_assignments,check.names = FALSE)
}

if (is.null(opt$min_sz)) {
  stop("Please provide the min_sz for GSVA!")
}else{
  min_sz = opt$min_sz
}

if (is.null(opt$max_sz)) {
  stop("Please provide the max_sz for GSVA!")
}else{
  max_sz = opt$max_sz
}

if (is.null(opt$category)) {
  stop("Please provide the MSigDB category!")
}else{
  category = opt$category
}

if (is.null(opt$subset)) {
  subset = NULL
}else{
  subset = opt$subset
  subset = toupper(gsub("_",":",subset))
}

if (is.null(opt$topn)) {
  stop("Please provide the top n pathways to be displayed in the heatmap!")
}else{
  topn = opt$topn
}

if (is.null(opt$number)) {
  stop("Please provide --number for cutree pathway clusters in the heatmap!")
}else{
  n_cut = opt$number
  if (n_cut < 1) stop("--number must be >= 1")
}

plot_width = opt$width
if (is.null(plot_width) || is.na(plot_width) || plot_width <= 0) {
  stop("--width must be a positive number")
}
plot_height = opt$height  # NULL = auto at plot time

message(paste("msigdbr category used in GSVA is:",category))
message(paste("GSVA gene set size filter (--min_sz / --max_sz):", min_sz, "-", max_sz))
message(paste("GSVA_mean_k3 cutree clusters (--number):", n_cut))
message(paste("GSVA_mean_k3 plot width (--width):", plot_width))
if (is.null(plot_height)) {
  message("GSVA_mean_k3 plot height (--height): auto")
} else {
  message(paste("GSVA_mean_k3 plot height (--height):", plot_height))
}

GO_gene_sets_full <- msigdbr::msigdbr(
  species = "Homo sapiens", # Can change this to what species you need
  category = category
)

# Determine which gene set collection(s) to run. Mirrors 06_GSEA.R's loop
# over GO Biological Process / Molecular Function / Cellular Component when
# the default C5 category is used with no explicit --subset. An explicit
# --subset (or a non-C5 --category) runs once, unsuffixed, as before.
if (category == "C5" && is.null(subset)) {
  message("category=C5, no --subset supplied: running GO Biological Process, ",
          "Molecular Function, and Cellular Component separately, matching 06_GSEA.R")
  run_plan = list(
    list(gs_subcat = "GO:BP", output_suffix = "_BP", label = "GO Biological Process (BP)"),
    list(gs_subcat = "GO:MF", output_suffix = "_MF", label = "GO Molecular Function (MF)"),
    list(gs_subcat = "GO:CC", output_suffix = "_CC", label = "GO Cellular Component (CC)")
  )
} else if (!is.null(subset)) {
  message("Running a single gene set collection: category=", category, ", subset=", subset)
  run_plan = list(list(gs_subcat = subset, output_suffix = "", label = paste0(category, " / ", subset)))
} else {
  message("Running a single gene set collection: category=", category, " (no subset)")
  run_plan = list(list(gs_subcat = NULL, output_suffix = "", label = category))
}

################################################################################################
# One-time setup shared by every gene set collection above: expression
# matrix (name x sample) and k=3 subtype design/contrasts.
#
# Matrix rownames = rownames(se), i.e. the SummarizedExperiment's own "name"
# column -- the SAME identifier 06_GSEA.R uses as names(geneList), so both
# scripts draw gene sets from the identical gene universe (see header
# comment). This replaces an earlier version that remapped rownames to
# "gene_name" via uniprot_to_genename.rds before running gsva().
filtered_mapped_matrix = assay(se)
colnames(filtered_mapped_matrix) = se@colData$label

# Create a heatmap of the GSVA results
metadata = se@colData %>% as.data.frame() %>% left_join(cluster_assignments,by = c("label" = "patid"))
metadata = metadata %>% column_to_rownames("label")

####for kmeans k = 3
message("Performing heatmap for kmeans k=3")
annolabel <- data.frame(Cluster = metadata[,"kmeans_k=3"])
rownames(annolabel) <- rownames(metadata)

group_info <- factor(annolabel$Cluster,levels = c("alpha","beta","theta"))

# Create a design matrix for the linear model
design <- model.matrix(~ 0 + group_info)
colnames(design) <- levels(group_info)  # Name the columns according to the groups

# Create contrasts for group comparisons
contrast_matrix <- makeContrasts(
  theta_vs_beta =  theta - beta,  # Comparison between Group 1 and Group 0
  alpha_vs_theta = alpha - theta,  # Comparison between Group 2 and Group 0
  alpha_vs_beta = alpha - beta,  # Comparison between Group 2 and Group 1
  levels = design
)

################################################################################################
# Run GSVA + k=3 subtype comparison for each gene set collection in run_plan
for (plan in run_plan) {
  run_gsva_geneset_analysis(
    gs_subcat = plan$gs_subcat,
    output_suffix = plan$output_suffix,
    label_for_messages = plan$label,
    GO_gene_sets_full = GO_gene_sets_full,
    filtered_mapped_matrix = filtered_mapped_matrix,
    min_sz = min_sz,
    max_sz = max_sz,
    topn = topn,
    n_cut = n_cut,
    plot_width = plot_width,
    plot_height = plot_height,
    annolabel = annolabel,
    group_info = group_info,
    design = design,
    contrast_matrix = contrast_matrix,
    output_dir = output_dir
  )
}

#########################################################
# ####for kmeans k = 2
# message("Performing heatmap for kmeans k=2")
# annolabel <- data.frame(Cluster = metadata[,"kmeans_k=2"])
# rownames(annolabel) <- rownames(metadata)
#
# group_info <- factor(annolabel$Cluster,levels = c("alpha","beta"))
#
# design <- model.matrix(~ 0 + group_info)
# colnames(design) <- levels(group_info)
#
# fit <- lmFit(gsva_results, design)
#
# contrast_matrix <- makeContrasts(
#   alpha_vs_beta = alpha - beta,
#   levels = design
# )
#
# fit2 <- contrasts.fit(fit, contrast_matrix)
# fit2 <- eBayes(fit2)
#
# alpha_vs_beta_all <- topTable(fit2, coef = "alpha_vs_beta", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE)
# alpha_vs_beta <- topTable(fit2, coef = "alpha_vs_beta", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE, p.value = 0.05)
#
# write.csv(alpha_vs_beta %>% as.data.frame()%>%rownames_to_column("terms"),quote = F,row.names = F, file = file.path(output_dir, "GSVA_alpha_vs_beta_k2.csv"))
#
# alpha_vs_beta_topn <- process_comparison(alpha_vs_beta, "alpha_vs_beta_k2", GO_list, filtered_mapped_matrix, topn, output_dir)
#
# alpha_vs_beta_up = alpha_vs_beta %>% dplyr::slice_max(n = topn/2,order_by = t)  %>% as.data.frame()
# alpha_vs_beta_down = alpha_vs_beta %>% dplyr::slice_min(n = topn/2,order_by = t)  %>% as.data.frame()
# alpha_vs_beta = rbind(alpha_vs_beta_up,alpha_vs_beta_down)
#
# write.csv(alpha_vs_beta %>% as.data.frame()%>%rownames_to_column("terms"),quote = F,row.names = F, file = file.path(output_dir, "GSVA_topn_terms_k2.csv"))
#
# ordered_columns <- rownames(annolabel)[order(annolabel$Cluster)]
# cluster_colors <- list(Cluster = c("alpha" = "#FC8D62","beta" = "#8DA0CB"))
#
# plot_df <- gsva_results[unique(rownames(alpha_vs_beta)), ordered_columns]
# legend_breaks <- seq(from = round(min(plot_df), 1), to = round(max(plot_df), 1), by = 0.2)
# rownames(plot_df) <- sapply(rownames(plot_df), abbreviate_terms_end)
# pathway_heatmap_k2 <- pheatmap::pheatmap(plot_df,
#                                       annotation_col = annolabel,
#                                       annotation_colors = cluster_colors,
#                                       fontsize_row = 6,fontsize_col = 6, cluster_cols = FALSE, color = colorRampPalette(c("navy", "white", "firebrick3"))(100),legend_breaks = legend_breaks, border_color = NA)
#
# pdf(file.path(output_dir, "GSVA_sample_level_k2.pdf"),width = 10,height = 4)
# print(pathway_heatmap_k2)
# dev.off()
#
# ### cluster means heatmap k = 2
# cluster_means <- gsva_results %>%
#   as.data.frame() %>%
#   t() %>%
#   as.data.frame() %>%
#   mutate(Cluster = annolabel$Cluster) %>%
#   group_by(Cluster) %>%
#   summarise(across(everything(), mean, na.rm = TRUE)) %>%
#   column_to_rownames(var = "Cluster") %>%
#   t()
#
# write_cluster_means = cluster_means %>% as.data.frame() %>% rownames_to_column("terms")
# write_cluster_means = enrich_cluster_means_table(
#   write_cluster_means,
#   contrast_list = list(alpha_vs_beta = alpha_vs_beta_all),
#   GO_list = GO_list,
#   assay_genes = rownames(filtered_mapped_matrix)
# )
# write.csv(write_cluster_means, quote = TRUE, row.names = FALSE,
#           file = file.path(output_dir, "GSVA_cluster_means_k2.csv"))
#
# cluster_means = cluster_means[unique(rownames(alpha_vs_beta)),]
# rownames(cluster_means) <- sapply(rownames(cluster_means), abbreviate_terms_end)
# pathway_heatmap_mean_2 <- pheatmap::pheatmap(cluster_means,
#                                            fontsize_row = 6, cluster_cols = FALSE, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), border_color = NA)
#
# pdf(file.path(output_dir, "GSVA_mean_k2.pdf"),width = 6,height = 4)
# print(pathway_heatmap_mean_2)
# dev.off()


# #####correlation testing relation to age on set
# if (unique(is.na(colData(se)$age_at_onset))){
#   message("age_at_onset is not available in the clinical data")
# }else{
# metadata = se@colData %>% as.data.frame()
# common_samples <- intersect(colnames(gsva_results), metadata$label)
# gsva_subset <- gsva_results[, common_samples]
# metadata_subset <- metadata[metadata$label %in% common_samples, ]
#
# correlations <- apply(gsva_subset, 1, function(x) cor(x, metadata_subset$age_at_onset, method = "pearson"))
# sorted_correlations <- sort(correlations, decreasing = TRUE)
# top_10_highest <- head(sorted_correlations, 10)
# top_10_lowest <- tail(sorted_correlations, 10)
#
# top_10_df <- data.frame(
#   Term = names(top_10_highest),
#   Correlation = top_10_highest,
#   Type = rep("Positive", 10)
# )
# bottom_10_df <- data.frame(
#   Term = names(top_10_lowest),
#   Correlation = top_10_lowest,
#   Type = rep("Negative", 10)
# )
# result_df <- rbind(top_10_df, bottom_10_df)
#
# barplot = ggplot(result_df, aes(x = reorder(Term, Correlation), y = Correlation, fill = Type)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +
#   labs(
#     title = "Age at Onset",
#     x = "GSVA Terms",
#     y = "Pearson Correlation"
#   ) +
#   theme_few() +
#   scale_fill_manual(values = c("Positive" = "salmon", "Negative" = "steelblue"))
#
# pdf(file.path(output_dir, "GSVA_correlation_age_onset.pdf"),width = 10,height = 4)
# print(barplot)
# dev.off()
# }
