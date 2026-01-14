#=========================Script Description=================================
# This script is used for visualize differential expression analysis between subclusters result (volcanoplot)
# cd /Users/xliu2942/Documents/Projects/MAXOMOD/TESTing/Pipeline
# Rscript 9Vis_Differential_expression_analysis_subclusters.R
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggrepel"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("ggpubr"))
suppressMessages(library("stringr"))
suppressMessages(library("pheatmap"))
suppressMessages(library("purrr"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("DEP"))
#===========================Function Definition=============================
final_colours = list(
  #red and blue for volcano plot
  volcano_plot = brewer.pal(n=10, "RdYlBu")[c(2,9)],
  male_female = brewer.pal(n=8, "Set2")[c(1,4)],
  #heatmap_scale = rev(brewer.pal(n = 11, "RdYlBu")),
  heatmap_scale = rev(brewer.pal(n=11, "RdBu")),
  clustering = brewer.pal(n=8, "Set2")[c(2,3,5)],
  age_scale = brewer.pal(n = 9, "YlGn"),
  disease_progression_scale = brewer.pal(n = 9, "YlOrRd"),
  onset = c(brewer.pal(n=11, "RdYlBu")[c(4,5)], "#B3B3B3"),
  disease_status = c(brewer.pal(n=11, "RdYlBu")[2],"#B3B3B3"),
  genetic_testing = c("#B3B3B3", brewer.pal(n = 11, "PRGn")[c(3, 9)]),
  center = c("purple4", "orange3"),
  neurofilaments = brewer.pal(n = 9, "PuBu"),
  pNFh_scale = brewer.pal(n = 9, "Purples"),
  age_at_onset_scale = brewer.pal(n = 9, "Blues"),
  slow_vital_capacity_scale = brewer.pal(n = 9, "Reds")
)

#give the vectors with discrete variables names
names(final_colours$volcano_plot) =  c("up", "down")
names(final_colours$male_female) = c("Male", "Female")
names(final_colours$clustering) = c("alpha", "beta", "theta")
names(final_colours$disease_status) = c("als", "ctrl")
names(final_colours$onset) = c("spinal", "bulbar", "ctrl")
names(final_colours$genetic_testing) = c("not_performed", "negative", "C9orf72")
names(final_colours$center) = c("goettingen", "munich")

volcano_plot <- function(df, alpha_sig, name_title,labels,case_color,ctrl_color){
  df <- df %>%
    dplyr::mutate(omic_type = case_when(x >= 0 & y >= (-log10(alpha_sig)) ~ "up",
                                        x <= (0) & y >= (-log10(alpha_sig)) ~ "down",
                                        TRUE ~ "ns")) 
  cols <- c("up" = case_color, "down" = ctrl_color, "ns" = "grey") 
  sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
  alphas <- c("up" = 0.7, "down" = 0.7, "ns" = 0.5)
  ggplot(data = df, aes(x,y)) + 
    geom_point(aes(colour = omic_type), 
               alpha = 0.5, 
               shape = 16,
               size = 3) + 
    geom_hline(yintercept = -log10(alpha_sig),
               linetype = "dashed") +
    geom_text_repel(data = filter(df, name %in% labels),
                    aes(label = name),
                    force = 1,
                    hjust = 1,
                    nudge_x = -0.3,
                    nudge_y = 1.5,
                    direction = "both",
                    max.overlaps = 20,
                    size = 4) + 
    geom_vline(xintercept = 0,linetype = "dashed")+
    scale_colour_manual(values = cols) + 
    scale_fill_manual(values = cols) + 
    labs(title = name_title,
         x = "log2(fold change)",
         y = expression(-log[10] ~ "(adjusted p-value)"),
         colour = "Differential \nExpression") +
    theme_classic() + # Select theme with a white background  
    theme(axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 15, hjust = 0.5),
          text = element_text(size = 14)) +
    annotate("text", x = 1, y = 0.5, label = paste0(sum(df$omic_type=="up"), " more abundant \n", sum(df$omic_type=="down"), " less abundant"))
}



# change size of the figure to relatively small
save_pheatmap_pdf <- function(x, filename, width=4, height=2.5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "8_Differential_expression_analysis_subclusters/res.rds",
              help = "8_Differential_expression_analysis_subclusters/res.rds."
  ),make_option(c("--output", "-o"),
                type = "character", default = "9_Vis_Differential_expression_analysis_subclusters",
                help = "output directory path."
  ),make_option(c("--seed", "-e"),
                type = "integer", default = 9,
                help = "set.seed"
  ),make_option(c("--label", "-l"),
                type = "character", default = NULL,
                help = "label proteins"
  ),make_option(c("--data", "-d"),
                type = "character", default = "2_Missing_Inspection_als/norm_imp_MinProb.rds",
                help = "2_Missing_Inspection_als/norm_imp_MinProb.rds"
  ),make_option(c("--cluster_assignments", "-c"), 
                type = "character", default = "5_Clustering_als/cluster_assignments_2.csv",
                help = "cluster_assignments_2.csv"
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

if (is.null(opt$label)) {
  message("NO label proteins provided, use default setting!")
  label = opt$label
}else{
  label = opt$label
  label <- unlist(strsplit(label, ","))
}

if (is.null(opt$data)) {
  stop("Please provide the cleaned SummarizedExperiment object after normalize and imputation file path!")
}else if (!file.exists(opt$data)) {
  stop("SummarizedExperiment object after normalize and imputation does not exist!")
}else{
  data = readRDS(opt$data)
}

if (is.null(opt$cluster_assignments)) {
  stop("Please provide the cluster_assignments file path!")
}else if (!file.exists(opt$cluster_assignments)) {
  stop("cluster_assignments file does not exist!")
}else{
  cluster_assignments = read.csv(opt$cluster_assignments,check.names = F)
}

###########################volcano plot################################
plots_FDR0.05 = plots_FDR0.1 = list()
l = 1

for(i in 1:length(res)){
  data_res = res[[i]]
  diff_index = grep("diff",colnames(data_res))
  fdr_index = grep("fdr",colnames(data_res))
  if(length(diff_index) == 1){
    title = paste0(names(res)[i], "_alpha_vs_beta")
    logFC = data_res[,diff_index]
    fdr = data_res[,fdr_index]
    df <- data.frame(x = logFC, 
                     y = -log10(fdr),
                     name = data_res$name)
    names(df) <- c("x","y","name")
    case_color = final_colours$clustering["alpha"]
    ctrl_color = final_colours$clustering["beta"]
    names(case_color) =NULL
    names(ctrl_color) =NULL
    
    top10_label_up <- df %>%
      filter(x >= 0 & y >= (-log10(0.05))) %>%  # Apply filtering conditions
      arrange(desc(x)) %>%  # Sort by 'x' in descending order
      slice_head(n = 10)
    
    top10_label_down <- df %>%
      filter(x <= 0 & y >= (-log10(0.05))) %>%  # Apply filtering conditions
      arrange(x) %>%  # Sort by 'x' in descending order
      slice_head(n = 10)
    
    plots_FDR0.05[[l]] = volcano_plot(df, 0.05 , paste0("Volcano plot clustering proteomics \nalpha = FDR 0.05\n", title),labels = c(top10_label_up$name,top10_label_down$name,label),case_color,ctrl_color)
    plots_FDR0.1[[l]] = volcano_plot(df, 0.1 , paste0("Volcano plot clustering proteomics \nalpha = FDR 0.1\n", title),labels =c(top10_label_up$name,top10_label_down$name,label),case_color,ctrl_color)
    names(plots_FDR0.05)[l] = names(plots_FDR0.1)[l] = title
    l = l+1
  }else{
    names = colnames(data_res)[diff_index]
    names = gsub("_diff", "", names)
    for(j in 1:length(names)){
      case =  str_extract(names[j], "^[^_]+")
      case_color = final_colours$clustering[case]
      names(case_color) =NULL
      ctrl = str_extract(names[j], "(?<=_vs_).*")
      ctrl_color = final_colours$clustering[ctrl]
      names(ctrl_color) = NULL
      title = paste0(names(res)[i], "_", names[j])
      logFC = data_res[,diff_index[j]]
      fdr = data_res[,fdr_index[j]]
      df <- data.frame(x = logFC, 
                       y = -log10(fdr),
                       name = data_res$name)
      names(df) <- c("x","y","name")
      
      top10_label_up <- df %>%
        filter(x >= 0 & y >= (-log10(0.05))) %>%  # Apply filtering conditions
        arrange(desc(x)) %>%  # Sort by 'x' in descending order
        slice_head(n = 10)
      
      top10_label_down <- df %>%
        filter(x <= 0 & y >= (-log10(0.05))) %>%  # Apply filtering conditions
        arrange(x) %>%  # Sort by 'x' in descending order
        slice_head(n = 10)
      
      plots_FDR0.05[[l]] = volcano_plot(df, 0.05 , paste0("Volcano plot clustering proteomics \nalpha = FDR 0.05\n", title),labels = c(top10_label_up$name,top10_label_down$name,label),case_color,ctrl_color)
      plots_FDR0.1[[l]] = volcano_plot(df, 0.1 , paste0("Volcano plot clustering proteomics \nalpha = FDR 0.1\n", title),labels = c(top10_label_up$name,top10_label_down$name,label),case_color,ctrl_color)
      names(plots_FDR0.05)[l] = names(plots_FDR0.1)[l] = title
      l = l+1
    }
  }
}


k2_plots = plots_FDR0.05[grep("k2", names(plots_FDR0.05))]
k3_plots = plots_FDR0.05[grep("k3", names(plots_FDR0.05))]

volcano_plots = c(k2_plots,k3_plots)

ggarrange(plotlist = volcano_plots, ncol = 2, nrow = 2)
title = "volcano_plots_kmeans"
ggsave(filename = paste0(output_dir,"/", title, ".pdf"), width = 12, height = 10, units = "in")


#######################Heatmap################################
important_protein_names = list()  
all_important_protein_names = list()

#k = 2 ALPHA
d = res[[1]]
d = d[d$alpha_vs_beta_diff>0,]
d = d[d$fdr<=0.05,]
write.csv(d, file = paste0(output_dir,"/","all_important_protein_names_k2_alpha.csv"),row.names = F)
all_important_protein_names$k2_alpha = d$name

#k = 2 BETA
d = res[[1]]
d = d[d$alpha_vs_beta_diff<0,]
d = d[d$fdr<=0.05,]
write.csv(d, file = paste0(output_dir,"/","all_important_protein_names_k2_beta.csv"),row.names = F)
all_important_protein_names$k2_beta = d$name


#k = 3 ALPHA
d = res[[2]]
d = d[,c("name", colnames(d)[grep("alpha", colnames(d))])]
d = d[d$alpha_vs_beta_diff>0 & d$alpha_vs_theta_diff>0,]
d = d[d$alpha_vs_beta_fdr<=0.05 & d$alpha_vs_theta_fdr<=0.05,]
write.csv(d, file = paste0(output_dir,"/","all_important_protein_names_k3_alpha.csv"),row.names = F)
all_important_protein_names$k3_alpha = d$name


#k = 3 BETA
d = res[[2]]
d = d[,c("name", colnames(d)[grep("beta", colnames(d))])]
d = d[d$alpha_vs_beta_diff<0 & d$theta_vs_beta_diff<0,]
d$alpha_vs_beta_diff = -d$alpha_vs_beta_diff
d$theta_vs_beta_diff = -d$theta_vs_beta_diff
d = d[d$alpha_vs_beta_fdr<=0.05 & d$theta_vs_beta_fdr<=0.05,]
write.csv(d, file = paste0(output_dir,"/","all_important_protein_names_k3_beta.csv"),row.names = F)
all_important_protein_names$k3_beta = d$name


#k = 3 THETA
d = res[[2]]
d = d[,c("name", colnames(d)[grep("theta", colnames(d))])]
d = d[d$alpha_vs_theta_diff<0 & d$theta_vs_beta_diff>0,]
d = d[d$alpha_vs_theta_fdr<0.05 & d$theta_vs_beta_fdr<0.05,]
write.csv(d, file = paste0(output_dir,"/","all_important_protein_names_k3_theta.csv"),row.names = F)
all_important_protein_names$k3_theta = d$name

#select only the important proteins
k2_proteins = rbind(
  cbind(protein = all_important_protein_names[["k2_alpha"]], cluster = "alpha"),
  cbind(protein = all_important_protein_names[["k2_beta"]], cluster = "beta")
)
k2_proteins = as.data.frame(as.matrix(k2_proteins))
k2_proteins$cluster = as.factor(k2_proteins$cluster)
k3_proteins = rbind(
  cbind(protein = all_important_protein_names[["k3_alpha"]], cluster ="alpha"),
  cbind(protein = all_important_protein_names[["k3_beta"]], cluster = "beta"),
  cbind(protein = all_important_protein_names[["k3_theta"]], cluster = "theta")
)
k3_proteins = as.data.frame(as.matrix(k3_proteins))
k3_proteins$cluster = as.factor(k3_proteins$cluster)

assay = t(assay(data))
rownames(assay) = colData(data)$label
assay_k2 = assay[,k2_proteins$protein]
assay_k3 = assay[,k3_proteins$protein]

cluster_assignments$`kmeans_k=2` = as.factor(cluster_assignments$`kmeans_k=2`)
cluster_assignments$`kmeans_k=3` = as.factor(cluster_assignments$`kmeans_k=3`)

ordered_k2 = cluster_assignments[order(cluster_assignments$`kmeans_k=2`),]
ordered_k3 = cluster_assignments[order(cluster_assignments$`kmeans_k=3`),]

assay_k2 = assay_k2[ordered_k2$patid,]
assay_k3 = assay_k3[ordered_k3$patid,]


assay_k2_scaled = scale(assay_k2)
assay_k3_scaled = scale(assay_k3)

assay_k2 = as.data.frame(t(assay_k2))
assay_k2_scaled = as.data.frame(t(assay_k2_scaled))
assay_k2_capped = assay_k2_scaled
assay_k2_capped[assay_k2_capped>3] = 3
assay_k2_capped[assay_k2_capped<(-3)] = -3

assay_k3 = as.data.frame(t(assay_k3))
assay_k3_scaled = as.data.frame(t(assay_k3_scaled))
assay_k3_capped = assay_k3_scaled
assay_k3_capped[assay_k3_capped>3] = 3
assay_k3_capped[assay_k3_capped<(-3)] = -3


# Organize heatmap data into a list, focusing only on the capped datasets
heatmap_data <- list(
  datasets = list(assay_k2_capped, assay_k3_capped),
  patids = list(ordered_k2$patid, ordered_k3$patid),
  patient_cluster = list(ordered_k2$`kmeans_k=2`, ordered_k3$`kmeans_k=3`),
  proteins = list(k2_proteins$protein, k3_proteins$protein),
  protein_cluster = list(k2_proteins$cluster, k3_proteins$cluster),
  titles = c("k2_capped", "k3_capped")
)

# Loop through the capped datasets to generate heatmaps
for (i in seq_along(heatmap_data$datasets)) {
  
  # Create column annotations for patient cluster assignments
  col_annotations <- data.frame(
    Cluster_assignment = as.character(heatmap_data$patient_cluster[[i]])
  )
  rownames(col_annotations) <- heatmap_data$patids[[i]]
  
  # Create row annotations for protein cluster subtypes
  row_annotations <- data.frame(
    Cluster_subtypes = heatmap_data$protein_cluster[[i]]
  )
  rownames(row_annotations) <- heatmap_data$proteins[[i]]
  
  # Define annotation colors based on final colors for clustering
  ann_colors <- list(
    Cluster_assignment = final_colours$clustering[names(final_colours$clustering)%in% unique(heatmap_data$patient_cluster[[i]])],
    Cluster_subtypes = final_colours$clustering[names(final_colours$clustering)%in% unique(heatmap_data$patient_cluster[[i]])]
  )
  
  # Set the title for the current heatmap
  title <- heatmap_data$titles[[i]]
  
  # Create the heatmap using pheatmap
  p <- pheatmap(
    heatmap_data$datasets[[i]],
    annotation_row = row_annotations,
    annotation_col = col_annotations,
    annotation_colors = ann_colors,
    annotation_names_row = FALSE,  # Disable row labels
    annotation_names_col = FALSE,  # Disable column labels
    color = final_colours$heatmap_scale,  # Set heatmap color scale
    cluster_rows = FALSE,  # Disable row clustering
    cluster_cols = FALSE,  # Disable column clustering
    fontsize = 5,  # Adjust the font size for readability
    border_color = NA,  # Remove borders
    main = title  # Set the main title of the heatmap
  )
  
  # Save the heatmap to a PDF file with different dimensions based on the dataset
  if (title == "k2_capped") {
    save_pheatmap_pdf(p, paste0(output_dir, "/heatmap_clustering_kmeans_", title, ".pdf"), width = 14, height = 14)
  } else {
    save_pheatmap_pdf(p, paste0(output_dir, "/heatmap_clustering_kmeans_", title, ".pdf"), width = 9, height = 8)
  }
}