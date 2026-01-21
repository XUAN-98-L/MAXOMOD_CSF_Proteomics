#=========================Script Description=================================
# This script is used for WGCNA analysis.
# Rscript 10_WGCNA_subclusters.R -i Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -o Discovery/10_WGCNA_subclusters -e 9 -c Discovery/08_Clustering_als/cluster_assignments_2.csv
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
prepare_toplot <- function(MEs, cluster_file, metaData, cluster_col) {
  cluster_assignments <- read.csv(cluster_file)
  cluster_assignments <- cluster_assignments[match(colnames(cleanDat), cluster_assignments$patid), ]
  factorMeta <- left_join(cluster_assignments, metaData, by = c("patid" = "label"))
  
  #grep the number from the cluster_col
  cluster_num <- as.numeric(gsub("\\D", "", cluster_col))
  colnames(factorMeta)[which(colnames(factorMeta) ==paste0("kmeans_k.",cluster_num))] <- cluster_col
  
  factorMeta[[cluster_col]] <- as.factor(factorMeta[[cluster_col]])
#  factorMeta <- factorMeta[order(factorMeta[[cluster_col]]), ]
  toplot <- t(MEs)
  toplot <- toplot[, factorMeta$patid]
  factorMeta <- factorMeta[, c(cluster_col, "sex", "condition", "genetics")]
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

# Function to plot heatmaps
plot_heatmap <- function(toplot, factorMeta, file_name, main_title, scale = "row") {
  pdf(file_name, onefile = FALSE,width = 10,height = 7)
  aheatmap(x = toplot, main = main_title, annCol = factorMeta, scale = scale,
           distfun = "correlation", hclustfun = "average", cexRow = 1, 
           col = blueWhiteRed(100), treeheight = 40, Colv = NA,
           annColors = list(
             k2 = setNames(brewer.pal(8, "Set2")[c(2, 3, 5)], c("alpha", "beta", "theta")),
             k3 = setNames(brewer.pal(8, "Set2")[c(2, 3, 5)], c("alpha", "beta", "theta")),
             sex = setNames(brewer.pal(n=8, "Set2")[c(1,4)], c("Male", "Female")),
             condition = setNames(c(brewer.pal(n=11, "RdYlBu")[c(4,5)], "#B3B3B3"), c("spinal", "bulbar", "ctrl")),
             genetics = setNames(brewer.pal(n = length(levels(factorMeta$genetics)), name = "PRGn")[seq_along(levels(factorMeta$genetics))],
  levels(factorMeta$genetics)
)))
  dev.off()
}

# Function to plot mean heatmap
plot_mean_heatmap <- function(toplot, factorMeta, cluster_col, file_name) {
  cluster_means <- toplot %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    mutate(Cluster = factorMeta[[cluster_col]]) %>%
    group_by(Cluster) %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    column_to_rownames(var = "Cluster") %>%
    t()
  
  n_rows <- nrow(cluster_means)
  pdf_height <- max(4, min(0.6 * n_rows, 20))  # Set max height to 20 inches
  
  pdf(file_name, onefile = FALSE, width = 4, height = pdf_height)
  
  aheatmap(x = cluster_means,
           main = "Mean Heatmap \n Eigengene-Trait Relationships",
           scale = "none",
           cexRow = 1,
           col = blueWhiteRed(100),
           treeheight = 40,
           Colv = NA)
  
  dev.off()
}

plot_mean_heatmap_with_star <- function(toplot, factorMeta, cluster_col, file_name) {
  # calculate means for each module
  cluster_means <- toplot %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    mutate(Cluster = factorMeta[[cluster_col]]) %>%
    group_by(Cluster) %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    column_to_rownames(var = "Cluster") %>%
    t()
  
  # extract module names and p-values
  raw_module_names <- rownames(toplot)
  p_values <- gsub(".*p = ", "", raw_module_names)
  module_names <- gsub("\\n.*", "", raw_module_names)
  rownames(cluster_means) <- module_names
  
  # build star matrix
  stars <- sapply(as.numeric(p_values), function(p) {
    if (is.na(p)) return("")
    else if (p < 0.001) return("***")
    else if (p < 0.01) return("**")
    else if (p < 0.05) return("*")
    else return("")
  })
  
  star_matrix <- matrix(rep(stars, ncol(cluster_means)), 
                        nrow = nrow(cluster_means), 
                        ncol = ncol(cluster_means),
                        byrow = FALSE)
  rownames(star_matrix) <- module_names
  colnames(star_matrix) <- colnames(cluster_means)
  
  n_rows <- nrow(cluster_means)
  pdf_height <- max(4, min(0.6 * n_rows, 20))
  
  if (dev.cur() != 1) dev.off() 
  pdf(file_name, width = 5.5, height = pdf_height)
  pheatmap::pheatmap(
    mat = cluster_means,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    display_numbers = star_matrix,
    number_color = "white",
    fontsize_number = 28,
    fontsize_row = 12,
    fontsize_col = 12,
    main = "Mean Heatmap\nEigengene–Trait Relationships",
    border_color = "black"
  )
  dev.off()
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

# function to generate top protein boxplot (sorted by kME)
plot_module_top10_boxplot <- function(module_name, dat_module, metadata, cluster_col) {
  # length of boxplot is the same as the min number in each module
  top_expr <- dat_module[1:6, colnames(cleanDat)]
  
  long_df <- as.data.frame(t(top_expr)) %>%
    rownames_to_column("Sample") %>%
    pivot_longer(-Sample, names_to = "Protein", values_to = "Expression") %>%
    left_join(metadata, by = c("Sample" = "label"))
  
  # actual group information
  group_levels <- unique(long_df[[cluster_col]])
  n_groups <- length(group_levels)
  
  # automatically generate pairwise comparison group (for more 3 groups)
  comparisons <- combn(group_levels, 2, simplify = FALSE)
  
  # based on cluster number, decide which test to use
  if (n_groups >= 2) {
    add_stats <- stat_compare_means(
      method = "wilcox.test",
      comparisons = comparisons,
      label = "p.format",        
      hide.ns = TRUE,            # don't show ns
      size = 3,
      tip.length = 0.01
    )
  } else {
    add_stats <- NULL
  }
  
  custom_colors <- setNames(brewer.pal(8, "Set2")[c(2, 3, 5)], c("alpha", "beta", "theta"))
  color_vals <- custom_colors[group_levels] 
  
  p <- ggplot(long_df, aes_string(x = cluster_col, y = "Expression", fill = cluster_col)) +
    geom_boxplot(alpha = 1, outlier.shape = NA) +  
    geom_jitter(width = 0.2, alpha = 1, size = 1)+
    facet_wrap(~Protein, ncol = 3, scales = "free_y")  +
    scale_fill_manual(values = color_vals)+
    labs(title = paste("Module", module_name, "-", cluster_col),
         x = NULL, y = "log2 Expression") +  
    #stat_compare_means(aes(label = ..p.format..), method = "wilcox.test", label = "p.format")+
    theme_classic(base_size = 11) +
    theme(strip.text = element_text(face = "bold"),
          legend.position = "none")
  
  if (!is.null(add_stats)) {
    p <- p + add_stats
  }
  
  return(p)
}


plot_top_ic_heatmap <- function(net,enrichment_list_with_IC, top_n = 5, fill_na = TRUE, output_dir = ".") {
  # get module names, remove MEgrey
  desired_order <- colnames(net$MEs)
  desired_order <- desired_order[desired_order != "MEgrey"]
  desired_order <- gsub("^ME", "", desired_order)  # remove "ME" prefix

  # extract top_n GO terms with minimum IC for each module
  go_ic_top_df <- lapply(names(enrichment_list_with_IC), function(module) {
    df <- enrichment_list_with_IC[[module]]
    if (nrow(df) == 0) return(NULL)
    df_top <- df %>%
      #arrange(IC, p.adjust) %>%
      arrange(p.adjust,IC) %>%
      slice_head(n = top_n) %>%
      mutate(Module = module)
    return(df_top)
  }) %>% bind_rows()
  
  #  build heatmap matrix (-log10 p.adjust)
  heatmap_mat <- go_ic_top_df %>%
    mutate(logpadj = -log10(p.adjust)) %>%
    dplyr::select(Module, Description, logpadj) %>%
    pivot_wider(names_from = Description, values_from = logpadj) %>%
    column_to_rownames("Module")
  
  # build the star matrix for significance
  star_matrix <- go_ic_top_df %>%
    mutate(stars = case_when(
      p.adjust < 0.001 ~ "***",
      p.adjust < 0.01 ~ "**",
      p.adjust < 0.05 ~ "*",
      TRUE ~ ""
    )) %>%
    dplyr::select(Module, Description, stars) %>%
    pivot_wider(names_from = Description, values_from = stars) %>%
    column_to_rownames("Module")
  
  # reorder the heatmap matrix and star matrix according to the desired order
  heatmap_mat <- heatmap_mat[desired_order, , drop = FALSE]
  star_matrix <- star_matrix[desired_order, , drop = FALSE]
  
  # replace NA
  if (fill_na) {
    heatmap_mat[is.na(heatmap_mat)] <- 0
    star_matrix[is.na(star_matrix)] <- ""
  }
  

  pdf_file <- file.path(output_dir, paste0("top_", top_n, "_GO_IC_heatmap.pdf"))
  
  n_rows <- nrow(heatmap_mat)
  n_cols <- ncol(heatmap_mat)
  row_height <- 1.5
  col_width <- 1
  pdf_height <- max(7, n_rows * row_height)
  pdf_width  <- max(7, n_cols * col_width)

  pdf(pdf_file, width = pdf_width, height = pdf_height)
  #pdf(pdf_file, width = 18, height = 14)
  
  pheatmap::pheatmap(
    mat = heatmap_mat,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("white", "red"))(100),
    border_color = "black",
    display_numbers = as.matrix(star_matrix),
    number_color = "white",
    fontsize_number = 28,
    fontsize = 12,
    main = paste0("Top ", top_n, " GO Terms per Module (Min IC)"),
    angle_col = 90
  )
  
  dev.off()
  
  message("Heatmap saved to: ", pdf_file)
}
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "02_Missing_Inspection_subclusters/norm_imp_MinProb.rds",
              help = "SummarizedExperiment object after normalize and imputation."
  ),make_option(c("--output", "-o"),
                type = "character", default = "10_WGCNA_subclusters",
                help = "output directory path."
  ),make_option(c("--seed", "-e"),
               type = "integer", default = 9,
               help = "set.seed"
  ),make_option(c("--cluster_assignments", "-c"), 
                type = "character", default = "08_Clustering_als/cluster_assignments_2.csv",
                help = "cluster_assignments_2.csv"
  ),make_option(c("--top_n", "-t"),
                type = "integer", default = 3,
                help = "Number of top GO terms to display in the heatmap (default: 3)."
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
  cleanDat_file = readRDS(opt$input)
}

if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
}else{
  seed = opt$seed
  set.seed(seed)
}

if (is.null(opt$cluster_assignments)) {
  stop("Please provide the cluster_assignments file path!")
}else if (!file.exists(opt$cluster_assignments)) {
  stop("cluster_assignments file does not exist!")
}else{
  cluster_assignments_2_file = opt$cluster_assignments
}

if (is.null(opt$top_n)) {
  top_n = 3
}else{
  top_n = opt$top_n
}
#################WGCNA##################
metaData = as.data.frame(colData(cleanDat_file))
cleanDat = assay(cleanDat_file)
colnames(cleanDat) = metaData$label

powers <- seq(2,15,by=1)
#sft <- pickSoftThreshold(t(cleanDat),
#                         powerVector=powers,
#                         corFnc="bicor",networkType="signed",RsquaredCut=0.8)

sft <- pickSoftThreshold(t(cleanDat),
                         powerVector=powers,
                         networkType="signed",RsquaredCut=0.85)

# Plot the results
pdf(paste0(output_dir,"/WGCNA_SoftThreshold.pdf"),height =12)
par(mfrow = c(2,1))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1.5,col="red")

# this line corresponds to using an R^2 cut-off of h
h = 0.85
abline(h=h,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=1.5,col="red")	 
dev.off()

## Run an automated network analysis
net <- blockwiseModules(t(cleanDat),power=sft$powerEstimate,deepSplit=4,minModuleSize=6,
                        mergeCutHeight=0.07, corType="pearson",networkType="signed",pamStage=TRUE,pamRespectsDendro=TRUE,reassignThreshold = 0.05,
                        verbose=3,saveTOMs=FALSE,maxBlockSize=10000)

table(net$colors)
length(table(net$colors))-1 # 7 modules

# genes belong to which module
moduleLabels = as.data.frame(net$colors)
colnames(moduleLabels) = "color"
moduleLabels$gene = rownames(moduleLabels)

## Set up a quantitative matrix for phenotype correlations
rownames(metaData) <- metaData[,"label"]

#####plot numeric virables
if (sum(!is.na(metaData$pNFh)) >20){
  numericMeta <- metaData[,c("age","Nfl","pNFh","progression_rate","slow_vital_capacity","age_at_onset","disease_duration")]
}else{
  numericMeta <- metaData[,c("age","Nfl","progression_rate","slow_vital_capacity","age_at_onset","disease_duration")]
}

numericMeta<-numericMeta[match(colnames(cleanDat),rownames(numericMeta)),]

## Plot dendrogram with module colors and trait correlations
saveRDS(net, paste0(output_dir,"/WGCNA_net.rds"))
MEList = moduleEigengenes(t(cleanDat), colors = net$colors)
MEs = MEList$eigengenes

# set pairwise.complete.obs to deal with missing values
geneSignificance <- cor(numericMeta,t(cleanDat),use="pairwise.complete.obs")
rownames(geneSignificance) <- colnames(numericMeta)
geneSigColors <- t(numbers2colors(t(geneSignificance),,signed=TRUE,lim=c(-1,1),naColor="black"))
rownames(geneSigColors) <- colnames(numericMeta)

pdf(paste0(output_dir,"/WGCNA_Dendrogram_numericVariable.pdf"))
plotDendroAndColors(dendro=net$dendrograms[[1]],
                    colors=t(rbind(net$colors,geneSigColors)),
                    cex.dendroLabels=1.2,addGuide=TRUE,
                    dendroLabels=FALSE,
                    groupLabels=c("Module Colors",colnames(numericMeta)))
dev.off()

###########Relationship between models
## Plot eigengene dendrogram/heatmap - using bicor
#MEs <- net$MEs
MEs <- net$MEs[, colnames(net$MEs) != "MEgrey"]
#plotEigengeneNetworks(MEs, "Eigengene Network", marHeatmap = c(3,4,2,2), marDendro = c(0,4,2,0),plotDendrograms = TRUE, xLabelsAngle = 90,heatmapColors=blueWhiteRed(50))

nMods <- ncol(MEs)
pdf_height <- max(6, min(nMods * 0.4 + 4, 15))
pdf(paste0(output_dir,"/WGCNA_EigengeneNetworks.pdf"), height = pdf_height)
plotEigengeneNetworks(
  MEs,
  setLabels = "Eigengene Network",
  marHeatmap = c(3, 4, 2, 2),
  marDendro = c(3, 4, 2, 2),
  plotDendrograms = TRUE,
  plotHeatmaps = TRUE,
  xLabelsAngle = 90,
  heatmapColors = blueWhiteRed(50),
  colorLabels = TRUE,coloredBarplot = TRUE,barplotMeans = TRUE
)
dev.off()

colnames(MEs) <- substr(colnames(MEs),3,100)
rownames(MEs) <- colnames(cleanDat)



## Get sigend kME values
tmpMEs <- MEs
colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
#kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc="bicor")
kMEdat <- signedKME(t(cleanDat), tmpMEs)

## Plot eigengene-trait correlations - using bicor
#MEcors <- bicorAndPvalue(MEs,numericMeta,alternative = "two.sided",use = "pairwise.complete.obs")
#moduleTraitCor <- MEcors$bicor
#moduleTraitPvalue <- MEcors$p
saveRDS(numericMeta, paste0(output_dir,"/numericMeta.rds"))

### Module-trait relationships
plot_df <- cbind(net$MEs, numericMeta)
module_names <- colnames(net$MEs)
clinical_vars <- colnames(numericMeta)
plot_list <- list()
for (mod in module_names) {
  mod_color <- gsub("ME", "", mod)  #  "MEgreen" -> "green"
  for (trait in clinical_vars) {
    df <- plot_df[, c(mod, trait)]
    df <- df[complete.cases(df), ]
    
    # Pearson and p-value calculation
    cor_test <- cor.test(df[[1]], df[[2]], method = "pearson")
    
    if (!is.na(cor_test$p.value) && cor_test$p.value < 0.1) {
      p <- ggplot(df, aes_string(x = trait, y = mod)) +
        geom_point(color = mod_color, size = 2.5, alpha = 1) +
        geom_smooth(method = "lm", se = TRUE, color = "black") +
        stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 3.5) +
        theme_classic(base_size = 12) +
        labs(
          x = trait,
          y = paste0("MEs (", mod_color, ")"),
          title = paste(mod_color, "vs", trait)
        )
      plot_list[[paste(mod, trait, sep = "_")]] <- p
    }
  }
}

n_show <- max(12, length(plot_list))
MEs_plot = do.call(gridExtra::grid.arrange, c(plot_list[1:n_show], ncol = 4))

ggsave(filename = file.path(output_dir, "MEs_vs_Clinical_Variables.pdf"), plot = MEs_plot, width = 12, height = 10, device = cairo_pdf)


# Pearson correlation
moduleTraitCor <- cor(MEs, numericMeta, use = "pairwise.complete.obs", method = "pearson")
# calculate p-value
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(MEs))

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mfrow=c(1,1))
par(mar = c(6, 8.5, 3, 3));

## Display the correlation values within a heatmap plot
colvec <- rep("white",100)
colvec[1:10] <- "red"
pdf(paste0(output_dir,"/WGCNA_ModuleTrait_heatmap.pdf"))
labeledHeatmap(Matrix = moduleTraitPvalue,
               xLabels = colnames(numericMeta),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = colvec,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(0,1),
               #main = paste("Module-trait relationships\n bicor r-value \n (p-value)"),
               main = paste("Module-trait relationships\n pearson r-value \n (p-value)"),
               cex.main=0.8)
dev.off()

pdf(paste0(output_dir,"/WGCNA_ModuleTrait_heatmap_gradient.pdf"),width = 6,height = 5.5)
# Calculate -log10(p-value)
logP <- -log10(moduleTraitPvalue)

# Assume logP and textMatrix are ready
max_logP <- floor(max(logP, na.rm = TRUE))  # rounding down, e.g., if max is 2.8, set to 2

# Define color mapping function
col_fun <- colorRamp2(
  c(0, 1, max_logP),
  c("white", "white", "red")
)

# Define tick positions for the legend
legend_breaks <- 0:max_logP  # Add tick every 1 unit

# Create the heatmap object
ht <- Heatmap(
  matrix = logP,
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
    title = expression(-log[10](p-value)),
    border = "black",
    at = legend_breaks,                     # Add tick positions
    labels = as.character(legend_breaks)    # Tick labels
  )
)

# Draw the heatmap
draw(ht)
dev.off()

######Correlation plot
#moduleTraitCor = cor(MEs,numericMeta,use = "p")
#moduleTraitPvalue = corPvalueStudent(cor(MEs,numericMeta,use = "p"),dim(cleanDat)[2])

# Function to calculate p-values for each module
#calculate_pvalues <- function(MEs, toplot, factorMeta, cluster_col) {
#  lm1 <- lm(data.matrix(MEs[colnames(toplot), ]) ~ get(cluster_col), data = factorMeta)
#  pvec <- sapply(1:ncol(MEs), function(i) {
#    f <- summary(lm1)[[i]]$fstatistic
#    pf(f[1], f[2], f[3], lower.tail = FALSE)
#  })
#  names(pvec) <- colnames(MEs)
#  pvec <- pvec[match(names(pvec), rownames(toplot))]
#  rownames(toplot) <- paste(rownames(toplot), "\np = ", signif(pvec, 2), sep = "")
#  toplot
#}

# Generate plots for k2
data_k2 <- prepare_toplot(MEs, cluster_assignments_2_file, metaData, "k2")
toplot_k2 <- calculate_pvalues(MEs, data_k2$toplot, data_k2$factorMeta, "k2")
plot_heatmap(toplot_k2, data_k2$factorMeta, paste0(output_dir, "/k2.pdf"), 
             "Plot of Eigengene-Trait Relationships")
plot_mean_heatmap(toplot_k2, data_k2$factorMeta, "k2", 
                  paste0(output_dir, "/k2_mean_heatmap.pdf"))
plot_mean_heatmap_with_star(toplot_k2, data_k2$factorMeta, "k2", 
                  paste0(output_dir, "/k2_mean_heatmap_with_star.pdf"))
plot_boxplots(toplot_k2, data_k2$factorMeta, "k2", paste0(output_dir, "/k2_boxplot.pdf"))

# Generate plots for k3
data_k3 <- prepare_toplot(MEs, cluster_assignments_2_file, metaData, "k3")
toplot_k3 <- calculate_pvalues(MEs, data_k3$toplot, data_k3$factorMeta, "k3")
plot_heatmap(toplot_k3, data_k3$factorMeta, paste0(output_dir, "/k3.pdf"), 
             "Plot of Eigengene-Trait Relationships")
plot_mean_heatmap(toplot_k3, data_k3$factorMeta, "k3", 
                  paste0(output_dir, "/k3_mean_heatmap.pdf"))
plot_mean_heatmap_with_star(toplot_k3, data_k3$factorMeta, "k3", 
                  paste0(output_dir, "/k3_mean_heatmap_with_star.pdf"))
plot_boxplots(toplot_k3, data_k3$factorMeta, "k3", paste0(output_dir, "/k3_boxplot.pdf"))

##################################
#Find hub genes
chooseTopHubInEachModule(t(cleanDat),net$colors,power = sft$powerEstimate)
#save hub genes
hub_genes = as.data.frame(chooseTopHubInEachModule(t(cleanDat),net$colors))
hub_genes$module = rownames(hub_genes)
colnames(hub_genes)[1] = "gene"
write.csv(hub_genes,paste0(output_dir,"/hub_genes.csv"),row.names = FALSE)

###ORA enrichment analysis for moduleLabels
# for each module, get the genes
moduleLabels = as.data.frame(net$colors)
colnames(moduleLabels) = "color"
moduleLabels$gene = rownames(moduleLabels)

saveRDS(moduleLabels, paste0(output_dir,"/moduleLabels.rds"))
# according to color, split moduleLabels into a list
moduleLabels_list = split(moduleLabels,moduleLabels$color)
# write into excel file with different sheets
write.xlsx(moduleLabels_list, paste0(output_dir,"/module_genelist.xlsx"))

# Perform enrichment analysis for each module and save results
enrichment_results <- list()

for (color in names(moduleLabels_list)) {
  genes <- moduleLabels[moduleLabels$color == color,]$gene
  
  
  # Perform GO enrichment analysis (Biological Process)
  go_enrich <- enrichGO(gene = genes,
                        OrgDb = org.Hs.eg.db, 
                        keyType = 'SYMBOL',
                        ont = "BP",
                        qvalueCutoff  = 0.05,
                        pAdjustMethod = "BH",
                        minGSSize = 10)
  
  # simplify the GO terms
  #go_enrich = clusterProfiler::simplify(x = go_enrich, cutoff = 0.7)
  
  if(dim(go_enrich)[1] != 0){
    go_enrich = mutate(go_enrich, logpadj = -log(p.adjust, base=10)) %>% filter(Count >= 3) 
    
    # Store results in a list
    enrichment_results[[color]] <- go_enrich@result
    
    # save pdf
    pdf(paste0(output_dir,"/",color,"_GO_enrichment.pdf"))
    p = go_enrich %>% barplot(x="-log10(p.adjust)", showCategory=20,main = color)
    print(p)
    dev.off()
  }
}

# Write results to an Excel file with each module's enrichment in separate sheets
wb <- createWorkbook()
for (sheet_name in names(enrichment_results)) {
  addWorksheet(wb, sheetName = sheet_name)
  writeData(wb, sheet = sheet_name, enrichment_results[[sheet_name]])
}
saveWorkbook(wb, paste0(output_dir, "/module_enrichment_results.xlsx"), overwrite = TRUE)


# separate results by modules, order by kME
moduleColors = net$colors
dat.res <- data.frame(cleanDat,moduleColors , kMEdat)

# export kME values
export_dat.res <- dat.res[,setdiff(colnames(dat.res), colnames(cleanDat))] %>% rownames_to_column("Protein")
write.csv(export_dat.res, paste0(output_dir,"/module_kME.csv"), row.names = FALSE)

list.cluster.dat <- lapply(levels(as.factor(moduleColors)), 
                           function(x) {dtemp = dat.res[dat.res$moduleColors == x,];
                           dtemp[order(dtemp[,paste0('kME',x)==colnames(dtemp)], decreasing=TRUE),
                                 -setdiff(grep("^kME", colnames(dtemp)), which(paste0('kME',x)==colnames(dtemp)))]} )

names(list.cluster.dat) <- levels(as.factor(moduleColors))

# for list.cluster.dat, find top 10 kME proteins and make boxplot
cluster_assignments = read.csv(cluster_assignments_2_file)
metadata = as.data.frame(colData(cleanDat_file)) %>% left_join(cluster_assignments, by = c("label" = "patid"))

# make boxplot for each module
#for (mod in names(list.cluster.dat)) {
for (mod in setdiff(names(list.cluster.dat), "grey")) {
  dat_module <- list.cluster.dat[[mod]]
  
  # kmeans_k.2
  p_k2 <- plot_module_top10_boxplot(mod, dat_module, metadata, "kmeans_k.2")
  ggsave(filename = paste0(output_dir,"/top6_module_proteins_expression-", mod, "_k2.pdf"),
         plot = p_k2, width = 8, height = 6)
  
  # kmeans_k.3
  p_k3 <- plot_module_top10_boxplot(mod, dat_module, metadata, "kmeans_k.3")
  ggsave(filename = paste0(output_dir,"/top6_module_proteins_expression-", mod, "_k3.pdf"),
         plot = p_k3, width = 8, height = 6)
}
dev.off()

####===========================IC Vis==================================
#obtain the module enrichment results from the excel file
file_path <- paste0(output_dir, "/module_enrichment_results.xlsx")
sheet_names <- excel_sheets(file_path)

#remove "grey" module from the list
sheets_to_read <- setdiff(sheet_names, "grey")

# read each sheet into a list of data frames
enrichment_list <- lapply(sheets_to_read, function(sheet) {
  as.data.frame(read_excel(file_path, sheet = sheet))
})

# Give names to the list
names(enrichment_list) <- sheets_to_read

# Check the structure
str(enrichment_list, max.level = 1)

# Filter p.adjust < 0.1
enrichment_list_filtered <- lapply(enrichment_list, function(df) {
  df %>% filter(p.adjust < 0.1)
})

# Check the number of entries after filtering for each module
sapply(enrichment_list_filtered, nrow)

# set the organism for GOSemSim
hsGO <- godata('org.Hs.eg.db', ont = "BP")  


# enrichment_list_filtered is a list of data frames, each containing GO enrichment results for a module.
# Each data frame has columns: ID, Description, p.adjust, and other relevant information.
enrichment_list_with_IC <- list()

for (module in names(enrichment_list_filtered)) {
  df <- enrichment_list_filtered[[module]]
  
  # get GO ID and calculate IC value
  go_ids <- df$ID
  ic_values <- hsGO@IC[go_ids]

  # Some GO IDs may not exist in the IC table, so handle with NA
  df$IC <- ic_values[match(df$ID, names(ic_values))]

  # Save results
  enrichment_list_with_IC[[module]] <- df
}

#save.image(file = paste0(output_dir, "/WGCNA_results.RData"))
# use the function to plot the heatmap
plot_top_ic_heatmap(net, enrichment_list_with_IC, top_n = top_n, output_dir = output_dir)
