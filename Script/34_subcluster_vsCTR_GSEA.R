#=========================Script Description=================================
# This script is used for GSEA analysis of the subcluster vs control
# differential expression results from 32_subcluster_vsCTR.R: alpha vs ctrl
# and beta vs ctrl, each ranked and tested separately (same GSEA approach as
# 06_GSEA.R). Only the age+sex covariate, all-patients model is used
# (norm_imp_MinProb_subcluster_age_sex_cov_all_patients), matching
# 33_subcluster_vsCTR_vis.R.
# Rscript 34_subcluster_vsCTR_GSEA.R -i Discovery/32_subcluster_vsCTR/Differential_Expression_Results.rds -o Discovery/34_subcluster_vsCTR_GSEA
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("DEP"))
suppressMessages(library("clusterProfiler"))
suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library("msigdbr"))
suppressMessages(library("stringr"))
suppressMessages(library("ggplot2"))
#===========================Function Definition=============================
# Run GSEA (BP/MF/CC) for one comparison's signed -log10(p-value) ranking,
# save the raw results + CSV exports + a top-20-by-NES barplot, matching the
# outputs 06_GSEA.R produces for a single comparison.
# result: data.frame with columns "name", "<comp>_p.val", "<comp>_diff"
# comp: e.g. "alpha_vs_ctrl"; case_label: e.g. "Alpha" (used in plot legend)
run_gsea_for_comparison = function(result, comp, case_label, output_dir, seed){
  set.seed(seed)

  pval_col = paste0(comp, "_p.val")
  diff_col = paste0(comp, "_diff")
  if (!(pval_col %in% colnames(result)) || !(diff_col %in% colnames(result))){
    stop("Expected columns '", pval_col, "' and '", diff_col, "' not found for comparison '", comp, "'.")
  }

  signed_pval = -log10(result[[pval_col]]) * sign(result[[diff_col]])
  geneList = signed_pval
  names(geneList) = result$name
  geneList = geneList[!is.na(geneList)]
  geneList = sort(geneList, decreasing = TRUE)

  GSEA_result = list()
  ont = c("BP", "MF", "CC")
  for (o in ont){
    bg = msigdbr(species = "Homo sapiens", category = "C5", subcategory = o) %>%
      dplyr::select(gs_name, gene_symbol)
    bg = bg[bg$gene_symbol %in% names(geneList), ]

    GSEA_res = clusterProfiler::GSEA(geneList = geneList,
                                     nPermSimple = 100000,
                                     minGSSize = 10,
                                     maxGSSize = 500,
                                     pvalueCutoff = 1,
                                     verbose = TRUE,
                                     TERM2GENE = bg,
                                     pAdjustMethod = "BH",
                                     eps = 1e-10,
                                     seed = seed)
    GSEA_result[[o]] = GSEA_res
  }

  for (o in ont){
    sig = GSEA_result[[o]]@result[GSEA_result[[o]]@result$p.adjust < 0.05, ]
    if (nrow(sig) > 2){
      write.csv(sig, paste0(output_dir, "/", comp, "_GSEA_result_", o, "_005_cutoff.csv"))
      write.csv(GSEA_result[[o]]@result, paste0(output_dir, "/", comp, "_GSEA_result_", o, ".csv"))
    }
  }

  saveRDS(GSEA_result, paste0(output_dir, "/GSEA_result_", comp, ".rds"))

  # top 20 highest + top 20 lowest NES, BP only
  BP_top_highest = GSEA_result$BP@result %>%
    arrange(desc(NES)) %>%
    head(20)
  BP_top_lowest = GSEA_result$BP@result %>%
    arrange(NES) %>%
    head(20)
  BP_top = bind_rows(BP_top_highest, BP_top_lowest) %>%
    mutate(ID = str_remove(ID, "^GOBP_"),
          Description = str_remove(Description, "^GOBP_"))

  gsea_df = BP_top %>%
    dplyr::select(Description, NES, p.adjust) %>%
    mutate(logp = -log10(p.adjust),
          signed_logp = ifelse(NES > 0, logp, -logp)) %>%
    arrange(signed_logp) %>%
    mutate(Description = as.character(Description),
          Description_trunc = str_trunc(Description, width = 50, side = "right")) %>%
    arrange(signed_logp) %>%
    mutate(Description_trunc = factor(Description_trunc, levels = Description_trunc))

  p = ggplot(gsea_df, aes(x = signed_logp, y = Description_trunc, fill = NES > 0)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("TRUE" = "firebrick3", "FALSE" = "steelblue3"),
                      labels = c(paste0("Down in ", case_label), paste0("Up in ", case_label)),
                      name = "") +
    labs(x = expression(-log[10](adjusted~p~value)),
         y = "",
         title = paste0("GSEA: ", case_label, " vs Control")) +
    theme_classic(base_size = 12) +
    theme(axis.text.y = element_text(size = 9),
         axis.text.x = element_text(size = 10),
         axis.title.x = element_text(size = 13, face = "bold"),
         plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
         legend.position = "right")

  ggsave(paste0(output_dir, "/GSEA_BP_top20_", comp, ".pdf"), p, width = 12, height = 10)
  ggsave(paste0(output_dir, "/GSEA_BP_top20_", comp, ".svg"), p, width = 12, height = 10)
  write.csv(BP_top, paste0(output_dir, "/GSEA_BP_top20_", comp, ".csv"), row.names = FALSE)

  return(GSEA_result)}
# ===========================Command Parameters Setting=========================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "32_subcluster_vsCTR/Differential_Expression_Results.rds",
              help = "32_subcluster_vsCTR/Differential_Expression_Results.rds"
  ),make_option(c("--output", "-o"),
                type = "character", default = "34_subcluster_vsCTR_GSEA",
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
  stop("Please provide the cleaned SummarizedExperiment object after normalize and imputation file path!")
}else if (!file.exists(opt$input)) {
  stop("SummarizedExperiment object after normalize and imputation does not exist!")
}else{
  input = opt$input
  res = readRDS(input)
}

if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
}else{
  seed = opt$seed
  set.seed(seed)
}

#===========================Main Body=========================
# only the age+sex covariate, all-patients model (matches
# 32_subcluster_vsCTR.R / 33_subcluster_vsCTR_vis.R)
name_df = "norm_imp_MinProb_subcluster_age_sex_cov_all_patients"
result = res[[name_df]]

if (is.null(result)){
  stop("'", name_df, "' not found in ", input,
      ". Available comparisons: ", paste(names(res), collapse = ", "))
}

comparisons = c(alpha_vs_ctrl = "Alpha", beta_vs_ctrl = "Beta")

all_GSEA_results = list()
for (comp in names(comparisons)){
  message("Running GSEA for: ", comp)
  all_GSEA_results[[comp]] = run_gsea_for_comparison(result, comp, comparisons[[comp]], output_dir, seed)
}

saveRDS(all_GSEA_results, paste0(output_dir, "/GSEA_result_all_comparisons.rds"))
