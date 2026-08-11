#=========================Script Description=================================
# Grouped bar plot of subcluster-vs-control GSEA results (alpha vs ctrl,
# beta vs ctrl; from 34_subcluster_vsCTR_GSEA.R), Discovery and Validation
# side by side, mirroring the "ALS vs CTR" reference figure: x = signed
# -log10(FDR), one bar pair per GO term, terms grouped into blocks by GO
# semantic similarity with a representative label per block, dotted lines
# at the FDR<0.05 threshold. Alpha vs CTR and beta vs CTR are drawn as two
# separate figures.
#
# As in 13_GSEA_IC_heatmap.R / 27_GSEA_IC_heatmap_replot.R, terms are kept
# if significant (FDR < --fdr_cutoff) in EITHER cohort, but the signed
# -log10(FDR) shown for each cohort is its actual value even where that
# term individually misses the cutoff in that cohort (not forced to 0).
# Only the top --top_n terms (by smallest FDR across cohorts) are plotted,
# and the block label per cluster is the term with the lowest information
# content (IC) in that cluster -- an automatic stand-in for a manually
# curated category name, and may be worth relabelling by hand for
# publication.
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("msigdbr"))
suppressMessages(library("dplyr"))
suppressMessages(library("GOSemSim"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))
#===========================Function Definition=============================
# Format MSigDB pathway names for display, e.g. HUMORAL_IMMUNE_RESPONSE -> Humoral Immune Response
format_pathway_label = function(term){
  x = gsub("^(GOBP|GOCC|GOMF|HP|REACTOME|KEGG|HALLMARK)_", "", term)
  x = gsub("_", " ", x)
  stringr::str_to_title(tolower(x))
}

# Build the term x cohort signed -log10(FDR) table for one comparison,
# keeping terms significant in either cohort, using each cohort's FULL
# (unfiltered) result table to look up the value regardless of whether
# that specific term is individually significant in that cohort.
build_term_matrix = function(gsea_list, ontology = "BP", fdr_cutoff = 0.05){
  ego_list = lapply(gsea_list, function(g) g[[ontology]]@result)

  all_go_terms = unique(unlist(lapply(ego_list, function(x) x$ID[x$p.adjust < fdr_cutoff])))
  if (length(all_go_terms) == 0){
    stop("No terms reached FDR < ", fdr_cutoff, " in either cohort for ontology '", ontology, "'.")
  }

  logp_mat = matrix(NA_real_, nrow = length(all_go_terms), ncol = length(ego_list),
                    dimnames = list(all_go_terms, names(ego_list)))
  for (cohort in names(ego_list)){
    idx = match(rownames(logp_mat), ego_list[[cohort]]$ID)
    p_adj = ego_list[[cohort]]$p.adjust[idx]
    nes = ego_list[[cohort]]$NES[idx]
    logp_mat[, cohort] = -log10(p_adj) * sign(nes)
  }

  min_fdr = sapply(all_go_terms, function(id){
    min(sapply(ego_list, function(x) x$p.adjust[match(id, x$ID)]), na.rm = TRUE)
  })
  description = ego_list[[1]]$Description[match(all_go_terms, ego_list[[1]]$ID)]
  # fall back to the other cohort's description if a term is absent from the first
  missing_desc = is.na(description)
  if (any(missing_desc) && length(ego_list) > 1){
    description[missing_desc] = ego_list[[2]]$Description[match(all_go_terms[missing_desc], ego_list[[2]]$ID)]
  }

  list(logp_mat = logp_mat, min_fdr = min_fdr, description = description)}

# Map msigdbr gs_name (e.g. "GOBP_XXXX") to the exact GO ID (e.g.
# "GO:0007399") needed by GOSemSim, same lookup used in
# 27_GSEA_IC_heatmap_replot.R.
map_to_go_id = function(gs_names, ontology_subcat = "GO:BP"){
  bg_map = msigdbr(species = "Homo sapiens", category = "C5", subcategory = ontology_subcat) %>%
    distinct(gs_name, gs_exact_source)
  bg_map$gs_exact_source[match(gs_names, bg_map$gs_name)]}

# Cluster the given GO IDs by semantic similarity (Wang method), and label
# each cluster with the Description of its lowest-IC (most general) member.
cluster_terms = function(go_ids, descriptions, k, ontology = "BP"){
  go_sim = godata(annoDb = "org.Hs.eg.db", ont = ontology, computeIC = TRUE)
  ic_values = go_sim@IC[go_ids]

  sim_matrix = termSim(go_ids, go_ids, semData = go_sim, method = "Wang")
  go_dist = as.dist(1 - sim_matrix)
  k = min(k, length(go_ids) - 1)
  hc = hclust(go_dist, method = "ward.D2")
  cluster = cutree(hc, k = k)

  cluster_label = sapply(sort(unique(cluster)), function(cl){
    members = go_ids[cluster == cl]
    rep_id = members[which.min(ic_values[members])]
    descriptions[go_ids == rep_id][1]
  })
  names(cluster_label) = sort(unique(cluster))

  list(cluster = cluster, cluster_label = unname(cluster_label[as.character(cluster)]))}

# Grouped horizontal bar plot: one row per term, one bar per cohort,
# rows grouped into cluster blocks (separator line + block label), dotted
# lines at the FDR<0.05 threshold.
plot_gsea_barplot = function(df, case_label, cohort_colors, fdr_cutoff = 0.05,
                             wrap_width = 35){
  thresh = -log10(fdr_cutoff)

  # Title-case + wrap pathway / block labels for display
  df$description = format_pathway_label(as.character(df$description))
  df$cluster_label = format_pathway_label(as.character(df$cluster_label))

  # order rows: by cluster (clusters ordered by mean signed value, most
  # "up" first), then within cluster by signed value
  cluster_order = df %>%
    group_by(cluster) %>%
    summarise(m = mean(signed_logp, na.rm = TRUE)) %>%
    arrange(desc(m)) %>%
    pull(cluster)
  df$cluster = factor(df$cluster, levels = cluster_order)

  term_order = df %>%
    group_by(description) %>%
    summarise(cluster = cluster[1], m = mean(signed_logp, na.rm = TRUE)) %>%
    arrange(cluster, desc(m)) %>%
    pull(description)
  df$description = factor(df$description, levels = rev(term_order))

  # y position (row index) of each cluster's rows, for separator lines and
  # a block label placed at the vertical middle of the block, near x = 0
  row_pos = data.frame(description = levels(df$description),
                       y = seq_along(levels(df$description)))
  df = left_join(df, row_pos, by = "description")

  block_info = df %>%
    distinct(description, cluster, cluster_label, y) %>%
    group_by(cluster, cluster_label) %>%
    summarise(y_min = min(y), y_max = max(y), y_mid = mean(y), .groups = "drop") %>%
    arrange(y_min)

  # Separators between blocks: draw above each block except the topmost.
  # (Do NOT use block_info[-nrow(...)]: group_by order is by cluster id,
  # not visual y-order, which previously dropped the wrong block and left
  # two adjacent groups without a dividing line.)
  block_separators = block_info %>% dplyr::filter(y_max < max(y_max))

  xrange = max(abs(df$signed_logp), na.rm = TRUE) * 1.15
  n_rows = max(df$y)
  header_y = n_rows + 1.2 # one row above the topmost term, for the Down/Up header text

  p = ggplot(df, aes(x = signed_logp, y = description, fill = cohort)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
    geom_vline(data = data.frame(x = c(-thresh, thresh)),
              aes(xintercept = x, linetype = "FDR<0.05"), colour = "grey30") +
    geom_hline(data = block_separators, aes(yintercept = y_max + 0.5),
              colour = "black", linewidth = 0.4) +
    geom_text(data = block_info, aes(x = 0, y = y_mid, label = str_wrap(cluster_label, 18)),
             inherit.aes = FALSE, size = 3, fontface = "italic") +
    annotate("text", x = -xrange * 0.5, y = header_y,
            label = paste0("Down in ", case_label), fontface = "bold", size = 4) +
    annotate("text", x = xrange * 0.5, y = header_y,
            label = paste0("Up in ", case_label), fontface = "bold", size = 4) +
    scale_fill_manual(values = cohort_colors, name = NULL) +
    scale_linetype_manual(values = c("FDR<0.05" = "dotted"), name = NULL) +
    scale_x_continuous(limits = c(-xrange, xrange)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = wrap_width),
                     expand = expansion(add = c(0.6, 2))) +
    labs(title = paste0(case_label, " vs CTR"),
         x = "signed -log10(FDR)", y = NULL) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
         axis.text.y = element_text(size = 8, lineheight = 0.9),
         legend.position = "right") +
    coord_cartesian(clip = "off")

  return(p)}
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--discovery_dir", "-D"),
              type = "character", default = "Discovery/34_subcluster_vsCTR_GSEA",
              help = "Discovery/34_subcluster_vsCTR_GSEA output folder."
  ),make_option(c("--validation_dir", "-V"),
                type = "character", default = "Validation/34_subcluster_vsCTR_GSEA",
                help = "Validation/34_subcluster_vsCTR_GSEA output folder."
  ),make_option(c("--ontology", "-t"),
                type = "character", default = "BP",
                help = "GO ontology to plot (BP/MF/CC). Default BP."
  ),make_option(c("--fdr_cutoff", "-f"),
                type = "numeric", default = 0.05,
                help = "FDR cutoff for term inclusion. Default 0.05."
  ),make_option(c("--top_n", "-n"),
                type = "integer", default = 25,
                help = "Number of top terms (by smallest FDR across cohorts) to plot. Default 25."
  ),make_option(c("--cluster_number", "-c"),
                type = "integer", default = 4,
                help = "Number of GO semantic-similarity blocks. Default 4."
  ),make_option(c("--output", "-o"),
                type = "character", default = "35_subcluster_vsCTR_GSEA_vis",
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
  set.seed(9)
} else {
  set.seed(opt$seed)
}

cohort_colors = c(Discovery = "salmon", Validation = "#26b3ff")
comparisons = c(alpha_vs_ctrl = "Alpha", beta_vs_ctrl = "Beta")
ontology_subcat = paste0("GO:", opt$ontology)

################################################################################################
for (comp in names(comparisons)){
  case_label = comparisons[[comp]]
  message("Plotting: ", comp)

  gsea_list = list(
    Discovery = readRDS(file.path(opt$discovery_dir, paste0("GSEA_result_", comp, ".rds"))),
    Validation = readRDS(file.path(opt$validation_dir, paste0("GSEA_result_", comp, ".rds")))
  )

  built = build_term_matrix(gsea_list, ontology = opt$ontology, fdr_cutoff = opt$fdr_cutoff)

  # keep the top_n terms by smallest FDR across cohorts
  keep_ids = names(sort(built$min_fdr))[seq_len(min(opt$top_n, length(built$min_fdr)))]
  logp_mat = built$logp_mat[keep_ids, , drop = FALSE]
  description = built$description[match(keep_ids, rownames(built$logp_mat))]
  description = format_pathway_label(description)

  go_ids = map_to_go_id(keep_ids, ontology_subcat = ontology_subcat)
  valid = !is.na(go_ids)
  if (sum(valid) < 3){
    warning(comp, ": fewer than 3 terms could be mapped to GO IDs for semantic similarity clustering; skipping.")
    next
  }
  logp_mat = logp_mat[valid, , drop = FALSE]
  description = description[valid]
  go_ids = go_ids[valid]

  clustered = cluster_terms(go_ids, description, k = opt$cluster_number, ontology = opt$ontology)

  df = as.data.frame(logp_mat) %>%
    mutate(description = description, cluster = clustered$cluster, cluster_label = clustered$cluster_label) %>%
    tidyr::pivot_longer(cols = c(Discovery, Validation), names_to = "cohort", values_to = "signed_logp") %>%
    mutate(cohort = factor(cohort, levels = c("Discovery", "Validation")))

  write.csv(df, paste0(output_dir, "/GSEA_barplot_data_", comp, ".csv"), row.names = FALSE)

  p = plot_gsea_barplot(df, case_label, cohort_colors, fdr_cutoff = opt$fdr_cutoff)
  # Taller figure when many terms so wrapped y-axis labels have room
  plot_height = max(8, min(opt$top_n, nrow(logp_mat)) * 0.35 + 2)
  ggsave(paste0(output_dir, "/GSEA_barplot_", comp, ".pdf"), p,
         width = 10, height = plot_height, dpi = 300)
  ggsave(paste0(output_dir, "/GSEA_barplot_", comp, ".svg"), p,
         width = 10, height = plot_height, dpi = 300)
}
