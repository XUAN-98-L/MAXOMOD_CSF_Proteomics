#=========================Script Description=================================
# Rank-rank hypergeometric overlap (RRHO2) analysis of cross-cohort proteomic
# reproducibility in ALS vs. Controls.
#
# To assess whether ALS-associated protein changes were reproducible between
# independent cohorts, proteins in the discovery (x-axis) and validation
# (y-axis) datasets were each ranked by signed -log10(nominal P-value), with
# the sign indicating the direction of change in ALS relative to controls
# (positive = up-regulated in ALS, negative = down-regulated in ALS).
#
# Discovery is plotted on the x-axis, Validation on the y-axis (Discovery is
# passed as list_x/list1, Validation as list_y/list2 below). Heatmap color
# intensity indicates the statistical significance of overlap between the
# ranked lists at each pair of rank thresholds (hypergeometric test). For
# display, the matrix is rotated 180 degrees so that warm colors in the
# upper-right = concordant up-up (uu) and lower-left = concordant down-down
# (dd). The other two quadrants are discordant:
#   top-left    = down in Discovery, up in Validation   ("du")
#   bottom-right = up in Discovery,   down in Validation ("ud")
# Because Discovery = list1 and Validation = list2 here, RRHO2's own
# internal genelist_ud/genelist_du fields (which follow [list1][list2]
# order) already match this project's [Discovery][Validation] labelling
# convention directly -- no relabelling/swapping is needed downstream (this
# was NOT the case in an earlier version of this script where Validation
# was list1/x-axis instead).
# The returned RRHO object keeps the unflipped hypermat so gene lists stay aligned.
#
# Requires the RRHO2 package (not on CRAN):
#   remotes::install_github("RRHO2/RRHO2")
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("data.table"))
suppressMessages(library("RRHO2"))
suppressMessages(library("openxlsx"))
suppressMessages(library("ggplot2"))
suppressMessages(library("grid"))
#===========================Function Definition=============================
# Quadrant-direction legend (Discovery x Validation, up/down x up/down),
# drawn to match the ACTUAL displayed orientation of the 180-degree-flipped
# heatmap above -- NOT RRHO2's native/unflipped orientation. Discovery is
# the x-axis, Validation is the y-axis. In the flipped heatmap: x-axis
# (Discovery) runs down-regulated (left) -> up-regulated (right); y-axis
# (Validation) runs down-regulated (bottom) -> up-regulated (top). So
# concordant up-up sits top-right, concordant down-down sits bottom-left,
# and the two discordant quadrants are top-left (down-Discovery/up-
# Validation, "du") and bottom-right (up-Discovery/down-Validation, "ud")
# -- matching the RRHO2_summary.csv / RRHO2_overlap_proteins_ud.csv /
# _du.csv naming below.
rrho2_quadrant_legend = function(){
  n_grad = 200

  # y-axis gradient bar (Validation): bottom = down-regulated (blue), top = up-regulated (red)
  y_bar = data.frame(
    ymin = seq(0, 8, length.out = n_grad + 1)[-(n_grad + 1)],
    ymax = seq(0, 8, length.out = n_grad + 1)[-1]
  )
  y_bar$value = seq(-1, 1, length.out = n_grad)

  # x-axis gradient bar (Discovery): left = down-regulated (blue), right = up-regulated (red)
  x_bar = data.frame(
    xmin = seq(0, 8, length.out = n_grad + 1)[-(n_grad + 1)],
    xmax = seq(0, 8, length.out = n_grad + 1)[-1]
  )
  x_bar$value = seq(-1, 1, length.out = n_grad)

  grad_colors = c("#3366CC", "white", "#CC3333")  # down (blue) - mid - up (red)

  # one row per (quadrant, cohort): arrow position/direction + quadrant label
  quadrants = data.frame(
    quadrant = rep(c("bottom-left", "top-right", "top-left", "bottom-right"), each = 2),
    cohort   = rep(c("Discovery", "Validation"), times = 4),
    cx       = rep(c(2, 6, 2, 6), each = 2),
    cy       = rep(c(2, 6, 6, 2), each = 2),
    direction = c("down", "down",  "up",   "up",
                 "down", "up",    "up",   "down"),
    concordance = rep(c("Concordant", "Concordant", "Discordant", "Discordant"), each = 2),
    stringsAsFactors = FALSE
  )
  # offset the two cohorts' arrows so they sit side by side within each quadrant
  quadrants$x = quadrants$cx + ifelse(quadrants$cohort == "Discovery", -0.7, 0.7)
  quadrants$yend = ifelse(quadrants$direction == "up", quadrants$cy + 0.9, quadrants$cy - 0.9)
  quadrants$ystart = ifelse(quadrants$direction == "up", quadrants$cy - 0.9, quadrants$cy + 0.9)
  # cohort label sits to the LEFT of its arrow, vertically centred on the
  # arrow (arrows are symmetric +-0.9 around cy, so the midpoint is just cy)
  quadrants$label_x = quadrants$x - 0.22
  quadrants$label_y = quadrants$cy

  quad_labels = data.frame(
    x = c(2, 6, 2, 6), y = c(2, 6, 6, 2) - 1.6,
    label = c("Concordant", "Concordant", "Discordant", "Discordant")
  )

  p = ggplot() +
    # axis gradient bars
    geom_rect(data = y_bar, aes(xmin = -1.3, xmax = -0.5, ymin = ymin, ymax = ymax, fill = value)) +
    geom_rect(data = x_bar, aes(ymin = -1.3, ymax = -0.5, xmin = xmin, xmax = xmax, fill = value)) +
    scale_fill_gradientn(colours = grad_colors, guide = "none") +
    # outer box + quadrant divider
    annotate("rect", xmin = 0, xmax = 8, ymin = 0, ymax = 8, fill = NA, colour = "black", linewidth = 0.6) +
    annotate("segment", x = 4, xend = 4, y = 0, yend = 8, colour = "black", linewidth = 0.4) +
    annotate("segment", x = 0, xend = 8, y = 4, yend = 4, colour = "black", linewidth = 0.4) +
    # per-cohort arrows + labels (label left-of-arrow, vertically centred)
    geom_segment(data = quadrants,
                aes(x = x, xend = x, y = ystart, yend = yend),
                arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
                linewidth = 0.4) +
    geom_text(data = quadrants, aes(x = label_x, y = label_y, label = cohort),
             size = 1.8, angle = 90) +
    geom_text(data = quad_labels, aes(x = x, y = y, label = label), size = 2.8, fontface = "italic") +
    # axis text labels (rotated to match the reference RRHO2 legend style).
    # y-axis (Validation): bottom (y=2, blue part of the bar) = Downregulated,
    # top (y=6, red part) = Upregulated. x-axis (Discovery): left (x=2,
    # blue) = Downregulated, right (x=6, red) = Upregulated.
    annotate("text", x = -2.1, y = 2, label = "Downregulated", angle = 90, colour = "#3366CC", size = 2.4, fontface = "italic") +
    annotate("text", x = -2.1, y = 6, label = "Upregulated", angle = 90, colour = "#CC3333", size = 2.4, fontface = "italic") +
    annotate("text", x = 2, y = -2.1, label = "Downregulated", colour = "#3366CC", size = 2.4, fontface = "italic") +
    annotate("text", x = 6, y = -2.1, label = "Upregulated", colour = "#CC3333", size = 2.4, fontface = "italic") +
    coord_fixed(xlim = c(-2.8, 8.5), ylim = c(-2.8, 8.5), clip = "off") +
    theme_void()

  return(p)}

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
                      labels = c("Discovery", "Validation"),
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

# x-axis = Discovery, y-axis = Validation
# stepsize=1 gives one pixel per protein rank (finest resolution)
pdf(paste0(output_dir, "/RRHO2_heatmap.pdf"), width = 7, height = 6)
RRHO_obj = rrho2_plot(discovery_list, validation_list,
                      labels = c("Discovery", "Validation"),
                      main_title = "RRHO2: ALS vs Control",
                      stepsize = opt$stepsize)
# copy the already-rendered plot to an svg device (avoids recomputing the
# rank-rank hypergeometric overlap a second time)
dev.copy(svg, filename = paste0(output_dir, "/RRHO2_heatmap.svg"), width = 7, height = 6)
dev.off()
dev.off()

saveRDS(RRHO_obj, paste0(output_dir, "/RRHO2_obj.rds"))

# Quadrant-direction legend matching the heatmap's displayed orientation
legend_plot = rrho2_quadrant_legend()
ggsave(paste0(output_dir, "/RRHO2_quadrant_legend.pdf"), legend_plot, width = 3, height = 3)
ggsave(paste0(output_dir, "/RRHO2_quadrant_legend.svg"), legend_plot, width = 3, height = 3)

# Summary table at the most significant UU pixel (rank thresholds from that pixel)
max_neg_log_p = max(RRHO_obj$hypermat, na.rm = TRUE)
min_hypergeometric_p = if (isTRUE(RRHO_obj$log10.ind)) 10^(-max_neg_log_p) else exp(-max_neg_log_p)

rrho2_summary = data.frame(
  comparison = "ALS_vs_Control",
  stat_type = "pvalue",
  method = RRHO_obj$method,
  max_neg_log_p = max_neg_log_p,
  min_hypergeometric_p = min_hypergeometric_p,
  # list1 = Discovery (x), list2 = Validation (y)
  rank_threshold_discovery = length(RRHO_obj$genelist_uu$gene_list1_uu),
  rank_threshold_validation = length(RRHO_obj$genelist_uu$gene_list2_uu),
  overlap_uu = length(RRHO_obj$genelist_uu$gene_list_overlap_uu),
  overlap_dd = length(RRHO_obj$genelist_dd$gene_list_overlap_dd),
  # RRHO2's own internal "ud"/"du" fields follow [list1=Discovery][list2=
  # Validation] order, which is now the SAME order as this project's
  # [Discovery][Validation] convention -- so, unlike an earlier version of
  # this script (when Validation was list1/x-axis), no swap is needed here:
  #   overlap_ud = up-Discovery,   down-Validation = RRHO2's genelist_ud (direct)
  #   overlap_du = down-Discovery, up-Validation   = RRHO2's genelist_du (direct)
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

# Overlap proteins at the most significant pixel of each quadrant (uu/dd/ud/du).
# ud/du are named in [Discovery][Validation] order, which now matches
# RRHO2's own internal genelist_ud/genelist_du fields directly (list1 =
# Discovery, list2 = Validation), so no swap is needed:
#   "ud" file = up-Discovery,   down-Validation = RRHO2's genelist_ud (direct)
#   "du" file = down-Discovery, up-Validation   = RRHO2's genelist_du (direct)
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
