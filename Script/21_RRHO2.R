#=========================Script Description=================================
# Rank-rank hypergeometric overlap (RRHO2) analysis of cross-cohort proteomic
# reproducibility in ALS vs. Controls.
#
# To assess whether ALS-associated protein changes were reproducible between
# independent cohorts, proteins in the discovery (y-axis) and validation
# (x-axis) datasets were each ranked by signed -log10(nominal P-value), with
# the sign indicating the direction of change in ALS relative to controls
# (positive = up-regulated in ALS, negative = down-regulated in ALS).
#
# Heatmap color intensity indicates the statistical significance of overlap
# between the ranked lists at each pair of rank thresholds (hypergeometric
# test). For display, the matrix is rotated 180 degrees so that warm colors
# in the upper-right = concordant up-up (uu) and lower-left = concordant
# down-down (dd); the other two quadrants (ud/du) are discordant. The
# returned RRHO object keeps the unflipped hypermat so gene lists stay aligned.
#
# Requires the RRHO2 package (not on CRAN):
#   remotes::install_github("RRHO2/RRHO2")
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("data.table"))
suppressMessages(library("RRHO2"))
suppressMessages(library("openxlsx"))
#===========================Function Definition=============================
# Build a signed -log10(nominal p-value) rank score for RRHO2.
# Positive = up-regulated in ALS vs Control, negative = down-regulated.
signed_log10p = function(p_val, diff){
  score = -log10(p_val)
  score[diff < 0] = -score[diff < 0]
  return(score)
}

# Run RRHO2 for a discovery/validation pair and plot the heatmap.
# list_x is plotted on the x-axis, list_y is plotted on the y-axis.
# Each is a data.frame with columns "name" and "score".
# RRHO2's native hypermat has uu at bottom-left and dd at top-right under
# image(). For plotting only, rotate 180 degrees so upper-right = uu and
# lower-left = dd. The returned object keeps the original hypermat so
# genelist_uu/dd/ud/du remain correctly indexed.
rrho2_plot = function(list_x,
                      list_y,
                      labels = c("Validation", "Discovery"),
                      main_title = "RRHO2: ALS vs Control",
                      stepsize = 1,
                      log10.ind = TRUE,
                      multipleTesting = "BH",
                      method = "hyper"){

  RRHO_obj = RRHO2_initialize(list_x, list_y,
                              labels = labels,
                              stepsize = stepsize,
                              log10.ind = log10.ind,
                              multipleTesting = multipleTesting,
                              method = method)

  # flip for display only: upper-right = uu, lower-left = dd
  RRHO_obj_plot = RRHO_obj
  RRHO_obj_plot$hypermat = RRHO_obj$hypermat[nrow(RRHO_obj$hypermat):1,
                                             ncol(RRHO_obj$hypermat):1]
  RRHO2_heatmap(RRHO_obj_plot, labels = labels, main = main_title)

  return(RRHO_obj)}
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "03_Differential_expression_analysis/Differential_Expression_Results.rds",
              help = "03_Differential_expression_analysis/Differential_Expression_Results.rds"
  ),make_option(c("--output", "-o"),
                type = "character", default = "21_RRHO2",
                help = "output directory path."
  ),make_option(c("--seed", "-e"),
                type = "integer", default = 9,
                help = "set.seed"
  ),make_option(c("--stepsize", "-s"),
                type = "integer", default = 1,
                help = "RRHO2 stepsize; smaller = finer heatmap (default: 1)"
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
}else{
  res_Discovery = readRDS(paste0("Discovery/",opt$input))
  res_Validation = readRDS(paste0("Validation/",opt$input))
}

if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
}else{
  seed = opt$seed
  set.seed(seed)
}

################################################################################################
#RUN RRHO2
inter = intersect(res_Discovery$norm_imp_MinProb_age_sex_cov_all_patients$name,
                  res_Validation$norm_imp_MinProb_age_sex_cov_all_patients$name)

res_Discovery_inter = res_Discovery$norm_imp_MinProb_age_sex_cov_all_patients[which(res_Discovery$norm_imp_MinProb_age_sex_cov_all_patients$name %in% inter),]
# sort res_Discovery by the order in inter
res_Discovery_inter = res_Discovery_inter[match(inter, res_Discovery_inter$name),]

res_Validation_inter = res_Validation$norm_imp_MinProb_age_sex_cov_all_patients[which(res_Validation$norm_imp_MinProb_age_sex_cov_all_patients$name %in% inter),]
# sort res_Validation by the order in inter
res_Validation_inter = res_Validation_inter[match(inter, res_Validation_inter$name),]

# signed -log10(nominal p-value) ranking metric, positive = up in ALS
discovery_list = data.frame(name = res_Discovery_inter$name,
                            score = signed_log10p(res_Discovery_inter$als_vs_ctrl_p.val,
                                                   res_Discovery_inter$als_vs_ctrl_diff),
                            stringsAsFactors = FALSE)

validation_list = data.frame(name = res_Validation_inter$name,
                             score = signed_log10p(res_Validation_inter$als_vs_ctrl_p.val,
                                                    res_Validation_inter$als_vs_ctrl_diff),
                             stringsAsFactors = FALSE)

# drop any proteins with missing p-value/logFC in either cohort
complete_names = intersect(discovery_list$name[complete.cases(discovery_list)],
                           validation_list$name[complete.cases(validation_list)])

discovery_list = discovery_list[discovery_list$name %in% complete_names,]
validation_list = validation_list[validation_list$name %in% complete_names,]

# x-axis = Validation, y-axis = Discovery
# stepsize=1 gives one pixel per protein rank (finest resolution)
pdf(paste0(output_dir, "/RRHO2_heatmap.pdf"), width = 7, height = 6)
RRHO_obj = rrho2_plot(validation_list, discovery_list,
                      labels = c("Validation", "Discovery"),
                      main_title = "RRHO2: ALS vs Control",
                      stepsize = opt$stepsize)
dev.off()

saveRDS(RRHO_obj, paste0(output_dir, "/RRHO2_obj.rds"))

# Summary table at the most significant UU pixel (rank thresholds from that pixel)
max_neg_log_p = max(RRHO_obj$hypermat, na.rm = TRUE)
min_hypergeometric_p = if (isTRUE(RRHO_obj$log10.ind)) 10^(-max_neg_log_p) else exp(-max_neg_log_p)

rrho2_summary = data.frame(
  comparison = "ALS_vs_Control",
  stat_type = "pvalue",
  method = RRHO_obj$method,
  max_neg_log_p = max_neg_log_p,
  min_hypergeometric_p = min_hypergeometric_p,
  # list1 = Validation (x), list2 = Discovery (y)
  rank_threshold_discovery = length(RRHO_obj$genelist_uu$gene_list2_uu),
  rank_threshold_validation = length(RRHO_obj$genelist_uu$gene_list1_uu),
  overlap_uu = length(RRHO_obj$genelist_uu$gene_list_overlap_uu),
  overlap_dd = length(RRHO_obj$genelist_dd$gene_list_overlap_dd),
  overlap_ud = length(RRHO_obj$genelist_ud$gene_list_overlap_ud),
  overlap_du = length(RRHO_obj$genelist_du$gene_list_overlap_du),
  stringsAsFactors = FALSE
)

write.table(rrho2_summary,
            file = paste0(output_dir, "/RRHO2_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(rrho2_summary,
          file = paste0(output_dir, "/RRHO2_summary.csv"),
          row.names = FALSE)

# Overlap proteins at the most significant pixel of each quadrant (uu/dd/ud/du)
overlap_proteins = list(
  uu = data.frame(protein = RRHO_obj$genelist_uu$gene_list_overlap_uu, stringsAsFactors = FALSE),
  dd = data.frame(protein = RRHO_obj$genelist_dd$gene_list_overlap_dd, stringsAsFactors = FALSE),
  ud = data.frame(protein = RRHO_obj$genelist_ud$gene_list_overlap_ud, stringsAsFactors = FALSE),
  du = data.frame(protein = RRHO_obj$genelist_du$gene_list_overlap_du, stringsAsFactors = FALSE)
)
openxlsx::write.xlsx(overlap_proteins, file = paste0(output_dir, "/RRHO2_overlap_proteins.xlsx"))
saveRDS(overlap_proteins, paste0(output_dir, "/RRHO2_overlap_proteins.rds"))
for (q in names(overlap_proteins)) {
  write.csv(overlap_proteins[[q]],
            file = paste0(output_dir, "/RRHO2_overlap_proteins_", q, ".csv"),
            row.names = FALSE)
}
