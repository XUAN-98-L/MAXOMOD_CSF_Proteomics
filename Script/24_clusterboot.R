#=========================Script Description=================================
# Cluster-wise stability assessment of the k-means subtype partition by
# non-parametric bootstrap resampling, using fpc::clusterboot
# (https://doi.org/10.1016/j.csda.2006.11.025).
#
# For each number of clusters k (default k = 2-10), k-means clustering is
# repeatedly applied to bootstrap resamples of the normalised, imputed
# protein matrix, and the stability of each original cluster is quantified
# as its mean Jaccard similarity to the most similar cluster recovered
# across resamples. Clusters with a mean Jaccard >= 0.75 are considered
# stable (>= 0.85 highly stable).
#
# Produces when --vis FALSE (default):
#   clusterboot_results.rds, mean_jaccard_by_k.csv, clusterwise_jaccard.csv,
#   clusterboot_summary.txt
# Produces when --vis TRUE (reads existing results; skips clusterboot):
#   mean_jaccard_by_k.pdf, clusterwise_stability.pdf
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("fpc"))
suppressMessages(library("ggplot2"))
#===========================Function Definition=============================
# Run fpc::clusterboot with kmeansCBI for every k in 2:kmax.
# assay_mat: samples (rows) x proteins (columns) numeric matrix, no NAs.
# Returns a named list of "clboot" objects, one per k.
run_clusterboot_range = function(assay_mat, kmax = 10, B = 500, nstart = 25,
                                 iter.max = 50, seed = 9){
  results = list()
  for (k in 2:kmax){
    message("clusterboot: k = ", k)
    set.seed(seed)
    cboot = clusterboot(assay_mat,
                        B = B,
                        bootmethod = "boot",
                        clustermethod = kmeansCBI,
                        krange = k,
                        runs = nstart,
                        iter.max = iter.max,
                        seed = seed,
                        count = FALSE)
    results[[as.character(k)]] = cboot
  }
  return(results)}

# Build tidy summary tables from a list of clusterboot results.
# Returns list(by_k = one row per k, by_cluster = one row per cluster per k)
summarise_clusterboot = function(clusterboot_results){
  by_k = data.frame(k = integer(), mean_jaccard = numeric())
  by_cluster = data.frame(k = integer(), cluster = character(),
                          size = integer(), jaccard = numeric())

  for (k_name in names(clusterboot_results)){
    k = as.integer(k_name)
    cboot = clusterboot_results[[k_name]]
    jacc = cboot$bootmean
    # original (full-data) partition sizes, ordered by cluster id 1..k
    sizes = as.integer(table(factor(cboot$partition, levels = seq_along(jacc))))

    by_k = rbind(by_k, data.frame(k = k, mean_jaccard = mean(jacc)))
    by_cluster = rbind(by_cluster,
                       data.frame(k = k,
                                 cluster = paste0("cluster_", seq_along(jacc)),
                                 size = sizes,
                                 jaccard = jacc))
  }
  by_k = by_k[order(by_k$k),]
  by_cluster = by_cluster[order(by_cluster$k),]
  return(list(by_k = by_k, by_cluster = by_cluster))}

# Fig 6: one bar per k, mean cluster-wise Jaccard
plot_mean_jaccard_by_k = function(by_k, B, cutoff){
  by_k$k = factor(by_k$k, levels = sort(unique(by_k$k)))

  ggplot(by_k, aes(x = k, y = mean_jaccard)) +
    geom_col(fill = "coral") +
    geom_hline(yintercept = cutoff, linetype = "dashed", colour = "grey30") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                       expand = expansion(mult = c(0, 0))) +
    #labs(title = "Mean cluster-wise Jaccard by k (fpc::clusterboot)",
    labs(title = NULL,
         #subtitle = paste0(B, " bootstrap runs; dashed line = cutoff (", cutoff, ")"),
         x = "Number of clusters (k)",
         y = "Mean cluster Jaccard") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 11))}

# Fig 7: per-cluster Jaccard, faceted by k; darker bars = stable (>= cutoff).
# Cluster size (original full-data partition) is annotated on each bar.
plot_clusterwise_stability = function(by_cluster, B, cutoff){
  by_cluster$k = factor(by_cluster$k,
                        levels = sort(unique(by_cluster$k)),
                        labels = paste0("k = ", sort(unique(by_cluster$k))))
  by_cluster$stable = by_cluster$jaccard >= cutoff
  if (!("size" %in% colnames(by_cluster))){
    stop("by_cluster is missing a 'size' column; rebuild summaries from clusterboot_results.rds.")
  }
  by_cluster$cluster = factor(by_cluster$cluster, levels = unique(by_cluster$cluster))

  ggplot(by_cluster, aes(x = cluster, y = jaccard, fill = stable)) +
    geom_col() +
    geom_hline(yintercept = cutoff, linetype = "dashed", colour = "grey30") +
    facet_wrap(~ k, scales = "free_x") +
    scale_fill_manual(values = c(`TRUE` = "grey40", `FALSE` = "grey75"), guide = "none") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
    geom_text(aes(label = paste0("n=", size)), vjust = -0.35, size = 2.5) +
    labs(title = NULL,
         subtitle = NULL,
         x = "Cluster",
         y = "Mean Jaccard similarity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA))}
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "02_Missing_Inspection_subclusters/norm_imp_MinProb.rds",
              help = "normalised, imputed SummarizedExperiment object used for k-means subtype clustering."
  ),make_option(c("--output", "-o"),
                type = "character", default = "24_clusterboot",
                help = "output directory path."
  ),make_option(c("--seed", "-e"),
                type = "integer", default = 9,
                help = "set.seed"
  ),make_option(c("--kmax", "-k"),
                type = "integer", default = 10,
                help = "maximum number of clusters to evaluate (range tested is 2:kmax). Default 10."
  ),make_option(c("--B", "-B"),
                type = "integer", default = 500,
                help = "number of bootstrap resampling runs for clusterboot. Default 500."
  ),make_option(c("--nstart", "-n"),
                type = "integer", default = 25,
                help = "number of random k-means initializations (kmeansCBI 'runs'). Default 25."
  ),make_option(c("--itermax", "-m"),
                type = "integer", default = 50,
                help = "maximum number of k-means iterations. Default 50."
  ),make_option(c("--cutoff", "-c"),
                type = "numeric", default = 0.75,
                help = "mean Jaccard cutoff for a cluster to be considered stable. Default 0.75."
  ),make_option(c("--vis", "-v"),
                type = "logical", default = FALSE,
                help = "if TRUE, skip clusterboot and plot from existing results in --output (clusterboot_results.rds / CSVs). Default FALSE runs clusterboot only."
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

################################################################################################
if (isTRUE(opt$vis)){
  # PLOT MODE: load saved results; do not re-run clusterboot
  rds_path = file.path(output_dir, "clusterboot_results.rds")
  by_k_path = file.path(output_dir, "mean_jaccard_by_k.csv")
  by_cluster_path = file.path(output_dir, "clusterwise_jaccard.csv")

  if (file.exists(rds_path)){
    # always rebuild by_cluster from RDS so size is available even for older CSVs
    message("Loading clusterboot results from ", rds_path)
    clusterboot_results = readRDS(rds_path)
    summary_tables = summarise_clusterboot(clusterboot_results)
    by_k = summary_tables$by_k
    by_cluster = summary_tables$by_cluster
    write.csv(by_k, by_k_path, row.names = FALSE)
    write.csv(by_cluster, by_cluster_path, row.names = FALSE)
  } else if (file.exists(by_k_path) && file.exists(by_cluster_path)){
    message("Loading summary tables from ", output_dir)
    by_k = read.csv(by_k_path, check.names = FALSE)
    by_cluster = read.csv(by_cluster_path, check.names = FALSE)
    if (!("size" %in% colnames(by_cluster))){
      stop("clusterwise_jaccard.csv has no 'size' column and clusterboot_results.rds is missing; cannot annotate cluster sizes.")
    }
  } else {
    stop("No clusterboot results found in ", output_dir,
         ". Run once with --vis FALSE to generate clusterboot_results.rds / CSVs first.")
  }

  fig6 = plot_mean_jaccard_by_k(by_k, B = opt$B, cutoff = opt$cutoff)
  ggsave(paste0(output_dir, "/mean_jaccard_by_k.pdf"), fig6, width = 7, height = 5, dpi = 300)

  fig7 = plot_clusterwise_stability(by_cluster, B = opt$B, cutoff = opt$cutoff)
  ggsave(paste0(output_dir, "/clusterwise_stability.pdf"), fig7, width = 11, height = 8, dpi = 300)

  message("Wrote plots to ", output_dir)

} else {
  # COMPUTE MODE: run clusterboot; do not plot
  if (is.null(opt$input)) {
    stop("Please provide the normalised, imputed SummarizedExperiment object file path!")
  }
  se = readRDS(opt$input)

  # samples (rows) x proteins (columns), matching 08_Clustering_subclusters.R
  assay_mat = t(as.matrix(assay(se)))

  clusterboot_results = run_clusterboot_range(assay_mat,
                                              kmax = opt$kmax,
                                              B = opt$B,
                                              nstart = opt$nstart,
                                              iter.max = opt$itermax,
                                              seed = seed)

  saveRDS(clusterboot_results, paste0(output_dir, "/clusterboot_results.rds"))

  summary_tables = summarise_clusterboot(clusterboot_results)
  write.csv(summary_tables$by_k, paste0(output_dir, "/mean_jaccard_by_k.csv"), row.names = FALSE)
  write.csv(summary_tables$by_cluster, paste0(output_dir, "/clusterwise_jaccard.csv"), row.names = FALSE)

  sink(paste0(output_dir, "/clusterboot_summary.txt"))
  cat("Mean cluster-wise Jaccard by k (", opt$B, " bootstrap runs, cutoff = ", opt$cutoff, "):\n", sep = "")
  print(summary_tables$by_k, row.names = FALSE)
  cat("\nPer-cluster Jaccard:\n")
  print(summary_tables$by_cluster, row.names = FALSE)
  sink()

  message("Skipping plots (--vis FALSE). Re-run with --vis TRUE to generate PDFs from saved results.")
}
