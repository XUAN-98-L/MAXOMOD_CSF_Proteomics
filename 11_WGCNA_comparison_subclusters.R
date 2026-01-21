#=========================Script Description=================================
# This script using WGCNA modules found from Discovery cohort to see it's expression and relationship with Validation cohort or external cohort
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("tidyverse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("clusterProfiler"))
suppressMessages(library("DEP"))
suppressMessages(library("ggpubr"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("stats"))
suppressMessages(library("WGCNA"))
suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library("openxlsx"))
suppressMessages(library("patchwork"))
suppressMessages(library("dplyr"))
suppressMessages(library("gridExtra"))
suppressMessages(library("NMF"))
suppressMessages(library("readxl"))
suppressMessages(library("GOSemSim"))
suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library("pheatmap"))
suppressMessages(library("tidyr"))
suppressMessages(library("ComplexHeatmap"))
suppressMessages(library("circlize"))
suppressMessages(library("grid"))
#===========================Function Definition=============================
# Function to prepare 'toplot' matrix
prepare_toplot <- function(MEs, cluster_file, metaData, cluster_col,ModuleTrait) {
  cluster_assignments <- read.csv(cluster_file)
  cluster_assignments <- cluster_assignments[match(colnames(cleanDat), cluster_assignments$patid), ]
  #remove kmeans_k columns from metaData
  metaData <- metaData %>% dplyr::select(-starts_with("kmeans_k"))
  factorMeta <- left_join(cluster_assignments, metaData, by = c("patid" = "label"))
  
  #grep the number from the cluster_col
  cluster_num <- as.numeric(gsub("\\D", "", cluster_col))
  colnames(factorMeta)[which(colnames(factorMeta) ==paste0("kmeans_k.",cluster_num))] <- cluster_col
  
  factorMeta[[cluster_col]] <- as.factor(factorMeta[[cluster_col]])
  #  factorMeta <- factorMeta[order(factorMeta[[cluster_col]]), ]
  toplot <- t(MEs)
  toplot <- toplot[, factorMeta$patid]
  if (ModuleTrait) {
    factorMeta <- factorMeta[, c(cluster_col, "sex", "condition", "genetics")]
  }else{
    factorMeta <- factorMeta[, cluster_col,drop = FALSE]
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
      # two group： Wilcoxon 
      pval <- tryCatch({
        wilcox.test(y ~ group)$p.value
      }, error = function(e) NA)
    } else if (num_groups >= 3) {
      # multi group：use ANOVA（linear model)
      #model <- lm(y ~ group)
      #f <- summary(model)$fstatistic
      #pval <- pf(f[1], f[2], f[3], lower.tail = FALSE)
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
              type = "character", default = "Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds",
              help = "SummarizedExperiment object after normalize and imputation."
  ),make_option(c("--output", "-o"),
                type = "character", default = "11_WGCNA_comparison",
                help = "output directory path."
  ),make_option(c("--net", "-n"),
                type = "character", default = "Discovery/23_WGCNA_updated/WGCNA_net.rds",
                help = "network object from Discovery cohort WGCNA."
  ),make_option(c("--seed", "-e"),
                type = "integer", default = 9,
                help = "set.seed"
  ),make_option(c("--cluster_assignments", "-c"), 
                type = "character", default = "Validation/08_Clustering_als/cluster_assignments_2.csv",
                help = "cluster_assignments_2.csv"
  ),make_option(c("--ModuleTrait", "-m"),
                type = "logical", default = FALSE,
                help = "Whether to plot Module-Trait relationships heatmap, default is TRUE.")
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
  cleanDat_file = readRDS(opt$input)
}

if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
}else{
  seed = opt$seed
  set.seed(seed)
}

if (is.null(opt$net)) {
  stop("Please provide the WGCNA network file path!")
}else if (!file.exists(opt$net)) {
  stop("WGCNA network file does not exist!")
}else{
  net = readRDS(opt$net)
}

if (is.null(opt$cluster_assignments)) {
  stop("Please provide the cluster_assignments file path!")
}else if (!file.exists(opt$cluster_assignments)) {
  stop("cluster_assignments file does not exist!")
}else{
  cluster_assignments_2_file = opt$cluster_assignments
}

if (is.null(opt$ModuleTrait)) {
  ModuleTrait = FALSE
}else{
  ModuleTrait = opt$ModuleTrait
}
#################Main script##################
metaData = as.data.frame(colData(cleanDat_file))
cleanDat = assay(cleanDat_file)
colnames(cleanDat) = metaData$label

## Set up a quantitative matrix for phenotype correlations
rownames(metaData) <- metaData[,"label"]

# Find common proteins between DC and VC
common_proteins = names(net$colors)[names(net$colors) %in% rownames(cleanDat)]
net_inDC = net$colors[common_proteins]
# save the common proteins and their color
common_protein_modules = as.data.frame(net_inDC)
common_protein_modules$Protein = rownames(common_protein_modules)
colnames(common_protein_modules)[1] = "module"
common_protein_modules = common_protein_modules[,c("Protein", "module")]
common_protein_modules = common_protein_modules[order(common_protein_modules$module), ]
write.csv(common_protein_modules, 
            file = paste0(output_dir, "/common_protein_modules.csv"), row.names = FALSE)


cleanDat = cleanDat[common_proteins, ]

MEList = moduleEigengenes(t(cleanDat), colors = net_inDC)
MEs = MEList$eigengenes
MEs <- MEs[, colnames(MEs) != "MEgrey"]
colnames(MEs) <- substr(colnames(MEs),3,100)

# reorder to the same order as DC
MEs = MEs[, setdiff(substr(colnames(net$MEs), 3, 100), "grey")]

# Generate plots for k2
data_k2 <- prepare_toplot(MEs, cluster_assignments_2_file, metaData, "k2", ModuleTrait = ModuleTrait)
toplot_k2 <- calculate_pvalues(MEs, data_k2$toplot, data_k2$factorMeta, "k2")
#plot_heatmap(toplot_k2, data_k2$factorMeta, paste0(output_dir, "/k2.pdf"), 
#             "Plot of Eigengene-Trait Relationships")
#plot_mean_heatmap(toplot_k2, data_k2$factorMeta, "k2", 
#                  paste0(output_dir, "/k2_mean_heatmap.pdf"))
#plot_mean_heatmap_with_star(toplot_k2, data_k2$factorMeta, "k2", 
#                            paste0(output_dir, "/k2_mean_heatmap_with_star.pdf"))
plot_boxplots(toplot_k2, data_k2$factorMeta, "k2", paste0(output_dir, "/k2_boxplot.pdf"))


# Generate plots for k3
data_k3 <- prepare_toplot(MEs, cluster_assignments_2_file, metaData, "k3", ModuleTrait = ModuleTrait)
toplot_k3 <- calculate_pvalues(MEs, data_k3$toplot, data_k3$factorMeta, "k3")
#plot_heatmap(toplot_k3, data_k3$factorMeta, paste0(output_dir, "/k3.pdf"), 
#             "Plot of Eigengene-Trait Relationships")
#plot_mean_heatmap(toplot_k3, data_k3$factorMeta, "k3", 
#                  paste0(output_dir, "/k3_mean_heatmap.pdf"))
#plot_mean_heatmap_with_star(toplot_k3, data_k3$factorMeta, "k3", 
#                            paste0(output_dir, "/k3_mean_heatmap_with_star.pdf"))
plot_boxplots(toplot_k3, data_k3$factorMeta, "k3", paste0(output_dir, "/k3_boxplot.pdf"))

############# Plot Module-Trait relationships heatmap
if (ModuleTrait){
#####plot numeric virables
if (sum(!is.na(metaData$pNFh)) >20){
  numericMeta <- metaData[,c("age","Nfl","pNFh","progression_rate","slow_vital_capacity","age_at_onset","disease_duration")]
}else{
  numericMeta <- metaData[,c("age","Nfl","progression_rate","slow_vital_capacity","age_at_onset","disease_duration")]
}

# Pearson correlation
moduleTraitCor <- cor(MEs, numericMeta, use = "pairwise.complete.obs", method = "pearson")
# calculate p-value
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(MEs))

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mfrow=c(1,1))
par(mar = c(6, 8.5, 3, 3));

# save workspace
#save.image(paste0(output_dir,"/WGCNA.RData"))
## Display the correlation values within a heatmap plot
#colvec <- rep("white",100)
#colvec[1:10] <- "red"
pdf(paste0(output_dir,"/WGCNA_ModuleTrait_heatmap_gradient.pdf"),width = 6,height = 5.5)

white_threshold <- -log10(0.1)
# Calculate -log10(p-value)
logP <- -log10(moduleTraitPvalue)

# signed -log10(p-value)
signed_logP <- sign(moduleTraitCor) * -log10(moduleTraitPvalue)

# Assume logP and textMatrix are ready
#max_logP <- floor(max(logP, na.rm = TRUE))  # rounding down, e.g., if max is 2.8, set to 2

##############################
# set the color all from -2 to 2
#max_abs_logP <- floor(max(abs(signed_logP), na.rm = TRUE))
max_abs_logP = 2

# Define color mapping function
#col_fun <- colorRamp2(
#  c(0, 1, max_logP),
#  c("white", "white", "red")
#)
#col_fun <- colorRamp2(
#  c(-max_abs_logP, 0, max_abs_logP),
#  c("blue", "white", "red")
#)
col_fun <- colorRamp2(
  c(-max_abs_logP, -white_threshold, 0, white_threshold, max_abs_logP),
  c("#4575B4", "white", "white", "white", "#D73027")
)

# Define tick positions for the legend
legend_breaks <- seq(-max_abs_logP, max_abs_logP, by = 1)  # Add tick every 1 unit

# Create the heatmap object
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

# Draw the heatmap
draw(ht)
dev.off()
}