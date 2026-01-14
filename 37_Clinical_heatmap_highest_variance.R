#=========================Script Description=================================
# This script is used for clinical heatmap visualization in the manuscript.
#Rscript 37_Clinical_heatmap_highest_variance.R --input Discovery/2_Missing_Inspection/norm_imp_MinProb.rds --output Discovery/37_Clinical_heatmap_top_variable --seed 9
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("pheatmap"))
suppressMessages(library("circlize"))
suppressMessages(library("ComplexHeatmap"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("DEP"))
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "Discovery/2_Missing_Inspection/norm_imp_MinProb.rds",
              help = "SummarizedExperiment object after normalize and imputation."
  ),make_option(c("--output", "-o"),
                type = "character", default = "37_Clinical_heatmap",
                help = "output directory path."
  ),make_option(c("--seed", "-e"),
                type = "integer", default = 9,
                help = "set.seed"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
#============================================================================
if (is.null(opt$input)) {
  stop("Please provide the cleaned SummarizedExperiment object after normalize and imputation file path!")
}else if (!file.exists(opt$input)) {
  stop("SummarizedExperiment object after normalize and imputation does not exist!")
}else{
  input = readRDS(opt$input)
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

if (is.null(opt$seed)) {
  set.seed(9)
} else {
  set.seed(opt$seed)
}

#============================================================================
expression = assay(input)
metadata = input@colData

# make sure the order of samples "ctrl first and als last"
ordered_samples <- rownames(metadata)[order(metadata$condition == "als")]

# reorder the expression matrix based on the ordered samples
mat <- expression[, ordered_samples]

# order by sample names
# anno_col <- as.data.frame(metadata[ordered_samples, c(
#   "condition",
#   "progression_group",
#   "progression_rate",
#   "slow_vital_capacity",
#   "pNFh",
#   "Nfl",
#   "sex",
#   "age",
#   "age_at_onset",
#   "genetics",
#   "onset"
# )])

# if (!identical(rownames(anno_col), colnames(mat))) {
#   stop("Error: Row names of annotation data do not match column names of the expression matrix.")
# }


# ann_colors <- list(
#   condition = c(ctrl = "white", als = "red"),
#   sex = c(Female = "lightpink", Male = "lightblue3"),
  
#   #Nfl = colorRampPalette(c("white", "skyblue", "navy"))(100),
#   Nfl = colorRampPalette(c("white", "orchid", "purple4"))(100),
#   pNFh = colorRampPalette(c("white", "orchid", "purple4"))(100),
#   progression_rate = colorRampPalette(c("white", "orangered", "darkred"))(100),
#   slow_vital_capacity = colorRampPalette(c("white", "orangered", "darkred"))(100),
#   age = colorRampPalette(c("white", "yellowgreen", "darkgreen"))(100),
#   age_at_onset = colorRampPalette(c("white", "yellowgreen", "darkgreen"))(100),
#   progression_group = c(
#     "SP" = "#1B9E77",  
#     "IP" = "#D95F02",   
#     "FP" = "#7570B3"  
#   ),
#   onset = c(spinal = "#FDAE61", bulbar = "#FEE090"),
#   genetics = c(
#     "negative" = "grey50",
#     "C9orf72" = "firebrick",
#     "not_performed" = "skyblue",
#     "VUS" = "orange"
#   )
# )

# my_breaks <- seq(-3, 3, length.out = 51)
# my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(length(my_breaks) - 1)

# #pdf(file.path(output_dir, "heatmap.pdf"), width = 12, height = 10)
# cairo_ps(file = file.path(output_dir, "clinical_heatmap.eps"), width = 12, height = 10, onefile = FALSE, bg = "transparent")
# pheatmap::pheatmap(
#   mat,
#   scale = "row",
#   color = my_palette,
#   breaks = my_breaks,
#   cluster_rows = TRUE,
#   cluster_cols = FALSE,
#   annotation_col = anno_col,
#   annotation_colors = ann_colors,
#   show_colnames = FALSE,
#   show_rownames = FALSE,
#   fontsize = 8,
#   border_color = NA,
#   treeheight_row = 0  # 不显示聚类树
# )
# dev.off()




# col annotation
anno_col <- as.data.frame(metadata[ordered_samples, c(
  "condition",
  "progression_group",
  "progression_rate",
  "slow_vital_capacity",
  "pNFh",
  "Nfl",
  "sex",
  "age",
  "age_at_onset",
  "genetics",
  "onset",
  "ECAS"
)])

if (!identical(rownames(anno_col), colnames(mat))) {
  stop("Error: Row names of annotation data do not match column names of the expression matrix.")
}

ha_colors <- list(
  condition = c(ctrl = "#4575B4", als = "#D73027"),
  sex = c(Female = "#E78AC3", Male = "#66C2A5"),
  Nfl = colorRamp2(
    quantile(anno_col$Nfl, probs = c(0, 0.5, 1), na.rm = TRUE),
    c("white", "orchid", "purple4")
  ),
  pNFh = colorRamp2(
    quantile(anno_col$pNFh, probs = c(0, 0.5, 1), na.rm = TRUE),
    c("white", "orchid", "purple4")
  ),
  progression_rate = colorRamp2(
    quantile(anno_col$progression_rate, probs = c(0, 0.5, 1), na.rm = TRUE),
    c("white", "orangered", "darkred")
  ),
  slow_vital_capacity = colorRamp2(
    quantile(anno_col$slow_vital_capacity, probs = c(0, 0.5, 1), na.rm = TRUE),
    c("white", "orangered", "darkred")
  ),
  age = colorRamp2(
    quantile(anno_col$age, probs = c(0, 0.5, 1), na.rm = TRUE),
    c("white", "yellowgreen", "darkgreen")
  ),
  age_at_onset = colorRamp2(
    quantile(anno_col$age_at_onset, probs = c(0, 0.5, 1), na.rm = TRUE),
    c("white", "yellowgreen", "darkgreen")
  ),
  progression_group = c(SP = "#1B9E77", IP = "#D95F02", FP = "#7570B3"),
  onset = c(spinal = "#FDAE61", bulbar = "#FEE090"),
  genetics = c(
    #negative = "grey50",
    negative = "white",
    C9orf72 = "firebrick",
    #not_performed = "skyblue",
    not_performed = "black",
    VUS = "orange",
    ROCK = "darkgreen",
    #SOD1 = "darkorange",
    SOD1 = "skyblue",
    SOD1_FIG4 = "darkviolet"
  ),
  ECAS = colorRamp2(
    quantile(anno_col$ECAS, probs = c(0, 0.5, 1), na.rm = TRUE),
    c("lightblue", "deepskyblue", "darkblue")
  )
)

# set the top annotation
ha <- HeatmapAnnotation(
  df = anno_col,
  col = ha_colors,
  annotation_height = unit(4, "mm"),
  #missing value color
  na_col = "black",
  annotation_name_side = "left"
)

# across samples  variance）
row_vars <- apply(mat, 1, var, na.rm = TRUE)
# Select top 100 most variable proteins
top100_genes <- names(sort(row_vars, decreasing = TRUE))[1:100]
# subset expression matrix
mat <- mat[top100_genes, ]

scaled_mat = t(scale(t(mat))) 
# limit color bar: clip to [-3, 3]
scaled_mat[scaled_mat > 3] <- 3
scaled_mat[scaled_mat < -3] <- -3

col_fun = colorRamp2(c(-3, 0, 3),c("#1B5B9D", "white", "#D7191C"))
col_fun(seq(-3, 3))

#cairo_ps(file = file.path(output_dir, "clinical_heatmap.eps"), width = 16/2, height = 14/2, onefile = FALSE, bg = "transparent")
cairo_ps(file = file.path(output_dir, "clinical_heatmap.eps"), width = 16/1.6, height = 14/2, onefile = FALSE, bg = "transparent")
ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
  scaled_mat,
  na_col = "white",
  #scale = "row",
  col = col_fun,
  #breaks = my_breaks,
  cluster_rows = TRUE,
  #cluster_cols = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_dend_width = unit(0, "mm"),
  show_row_dend = FALSE,
  top_annotation  = ha,
  heatmap_legend_param = list(
    at = c(-3, -2, -1, 0, 1, 2, 3),
    color_bar = "continuous",
    legend_height = unit(4, "cm"),
    title = NULL
  )
),heatmap_legend_side = "left",     
annotation_legend_side = "right",
padding = unit(c(5, 28, 5, 2), "mm")
)
dev.off()

#pdf(file = file.path(output_dir, "clinical_heatmap.pdf"), width = 16/2, height = 14/2, onefile = FALSE, bg = "transparent")
pdf(file = file.path(output_dir, "clinical_heatmap.pdf"), width = 16/1.6, height = 14/2, onefile = FALSE, bg = "transparent")
ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
  scaled_mat,
  na_col = "white",
  #scale = "row",
  col = col_fun,
  #breaks = my_breaks,
  cluster_rows = TRUE,
  #cluster_cols = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_dend_width = unit(0, "mm"),
  show_row_dend = FALSE,
  top_annotation  = ha,
  heatmap_legend_param = list(
    at = c(-3, -2, -1, 0, 1, 2, 3),
    color_bar = "continuous",
    legend_height = unit(4, "cm"),
    title = NULL
  )
),heatmap_legend_side = "left",     
annotation_legend_side = "right",
padding = unit(c(5, 28, 5, 2), "mm")
)
dev.off()