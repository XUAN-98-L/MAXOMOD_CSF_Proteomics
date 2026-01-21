#=========================Script Description=================================
# This script is used for visualize differential expression analysis between als vs ctrl result (volcano plot)
# Rscript 05_Vis_umap.R -i Discovery/02_Missing_Inspection -o Discovery/05_Vis_umap -e 9 -d all -l TRUE
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggrepel"))
suppressMessages(library("ggpubr"))
suppressMessages(library("stringr"))
suppressMessages(library("purrr"))
suppressMessages(library("umap"))
suppressMessages(library("ggthemes"))
suppressMessages(library("cowplot"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("DEP"))
#===========================Function Definition=============================
#all colourblind friendly options
display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, 
                   colorblindFriendly=TRUE)

#create a list with all colours for the final paper figures
#each list element will be a named vector
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

UMAP_density_plot = function(data, 
                             ggtitle, 
                             legend_name, 
                             labels, 
                             file_location, 
                             colour_set = c("seagreen4", "slateblue1", "salmon"),
                             patients_label = NULL){
      # run umap function
      umap_out = umap::umap(data)
      umap_plot = as.data.frame(umap_out$layout)
      
      #add condition labels
      umap_plot$group = labels

      # plot umap
      p1 = ggplot(umap_plot) + 
        geom_point(aes(x=V1, y=V2, color = as.factor(group)), alpha = 0.75, size = 4) +
        ggtitle(ggtitle) +
          theme_few() +
          scale_colour_few() +
          scale_color_manual(name = legend_name, 
                           labels = levels(as.factor(umap_plot$group)), 
                           values = colour_set) + 
          scale_fill_manual(values=colour_set) +
        labs(x = "UMAP1", y = "UMAP2")
      
      if (is.null(patients_label)){
        print("No patient labels provided")
      }else{
      p1 = p1 + geom_text(label = patients_label, x = umap_plot$V1, y = umap_plot$V2, hjust = 0, nudge_x = 1, size = 1.5, colour = "grey")
      }
  
      xdens <- 
        axis_canvas(p1, axis = "x") + 
        geom_density(data = umap_plot, aes(x = V1, fill = group, colour = group), alpha = 0.3) +
        scale_fill_manual( values = colour_set) + 
        scale_colour_manual( values = colour_set)
      ydens <-
        axis_canvas(p1, axis = "y", coord_flip = TRUE) + 
        geom_density(data = umap_plot, aes(x = V2, fill = group, colour = group), alpha = 0.3) +
        coord_flip() +
        scale_fill_manual(values = colour_set) + 
        scale_colour_manual( values = colour_set)
      p1 = p1 %>%
        insert_xaxis_grob(xdens, grid::unit(1, "in"), position = "top") %>%
        insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
        ggdraw()
    
      # save umap with labels
      ggsave(file_location,p1, width = 11/2, height = 8/2, units = "in")
}


# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "02_Missing_Inspection",
              help = "folder name which contains SummarizedExperiment object after normalize and imputation."
  ),make_option(c("--output", "-o"),
                type = "character", default = "05_Vis_umap",
                help = "output directory path."
  ),make_option(c("--seed", "-e"),
                type = "integer", default = 9,
                help = "set.seed"
  ),make_option(c("--condition", "-d"), default = "all",
                help = "condition information, default is all (als and ctrl). Otherwise can be als or ctrl."
  ),make_option(c("--cluster_assignments", "-c"), 
                type = "character", default = NULL,
                help = "cluster_assignments folder path (cluster_assignments_2.csv)"
  ),make_option(c("--label", "-l"), 
                type = "logical",default = TRUE,
                help = "whether to add patient labels to the UMAP plot"
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
  input = opt$input
}

if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
}else{
  seed = opt$seed
  set.seed(seed)
}

if (is.null(opt$condition)) {
  stop("Please provide the condition information!")
}else{
  condition = opt$condition
  
  if (condition == "all" & is.null(opt$cluster_assignments)) {
    cluster_assignments = opt$cluster_assignments
  }else if(condition == "als" | condition == "ctrl"){
    cluster_assignments = paste0(opt$cluster_assignments,"/cluster_assignments_2.csv")
    cluster_assignments = read.csv(cluster_assignments,check.names = FALSE)
  }
}
  
  
#===========================Main Script======================================
# Load data from input directories
data <- readRDS(file.path(input, "norm_imp_MinProb.rds"))
assay_data <- t(assay(data))

if(opt$label == TRUE){
  patients_label = colData(data)$label
  label_file_name = "label"
}else{
  patients_label = NULL
  label_file_name = "No_label"
}

#=========================== Helper Function to Generate UMAP Plots =============================
generate_umap_plots <- function(assay_data, output_dir, condition, label_file_name, labels_list, colour_list, patients_label) {
  for (label_name in names(labels_list)) {
    UMAP_density_plot(assay_data, 
                      ggtitle = paste0("UMAP with ", label_name, " ", condition), 
                      legend_name = label_name, 
                      labels = labels_list[[label_name]], 
                      file_location = file.path(output_dir, paste0("UMAP_", tolower(label_name), "_", condition, "_", label_file_name, ".pdf")), 
                      colour_set = colour_list[[label_name]], 
                      patients_label = patients_label)
  }
}

#=========================== Main Script =====================================
if (condition == "als") {
  # Define labels
  labels_list <- list(
    "Condition" = data$condition,
    "Sex" = data$sex,
    "Age" = data$age_cat
  )
  
  # Make cluster_assignments the same order as data
  cluster_assignments <- cluster_assignments[match(colData(data)$label, cluster_assignments$patid),]
  labels_list$`kmeans_k=2` <- as.factor(cluster_assignments$`kmeans_k=2`)
  labels_list$`kmeans_k=3` <- as.factor(cluster_assignments$`kmeans_k=3`)
  
  # Define color schemes
  colour_list <- list(
    "Condition" = final_colours$onset,
    "Sex" = final_colours$male_female,
    "Age" = final_colours$age_scale[c(4, 7)],
    "kmeans_k=2" = final_colours$clustering[names(final_colours$clustering) %in% levels(labels_list$`kmeans_k=2`)],
    "kmeans_k=3" = final_colours$clustering[names(final_colours$clustering) %in% levels(labels_list$`kmeans_k=3`)]
  )
  
  # Generate UMAP plots
  generate_umap_plots(assay_data, output_dir, condition, label_file_name, labels_list, colour_list, patients_label)
  
} else if (condition == "ctrl") {
  # Define labels
  labels_list <- list(
    "Condition" = data$condition,
    "Age" = data$age_cat,
    "Center" = data$center
  )
  
  # Make cluster_assignments the same order as data
  cluster_assignments <- cluster_assignments[match(colData(data)$label, cluster_assignments$patid),]
  labels_list$`kmeans_k=2` <- as.factor(cluster_assignments$`kmeans_k=2`)
  labels_list$`kmeans_k=3` <- as.factor(cluster_assignments$`kmeans_k=3`)
  
  # Define color schemes
  colour_list <- list(
    "Condition" = final_colours$male_female,
    "Age" = final_colours$age_scale[c(4, 7)],
    "Center" = final_colours$center,
    "kmeans_k=2" = final_colours$clustering[names(final_colours$clustering) %in% levels(labels_list$`kmeans_k=2`)],
    "kmeans_k=3" = final_colours$clustering[names(final_colours$clustering) %in% levels(labels_list$`kmeans_k=3`)]
  )
  
  # Generate UMAP plots
  generate_umap_plots(assay_data, output_dir, condition, label_file_name, labels_list, colour_list, patients_label)
  
} else if (condition == "all") {
  # Define labels
  labels_list <- list(
    "Condition" = data$condition,
    "Sex" = data$sex,
    "Age" = data$age_cat
  )
  
  # Add "Onset" conditionally if it exists
  if (!is.null(data$onset)) {
    labels_list$Onset <- data$onset
  }
  
  if (!is.null(data$center)) {
    labels_list$Center <- data$center
  }
  
  
  # Define color schemes
  colour_list <- list(
    "Condition" = final_colours$disease_status,
    "Sex" = final_colours$male_female,
    "Age" = final_colours$age_scale[c(4, 7)]
  )
  
  if (!is.null(data$onset)) {
    colour_list$Onset <- final_colours$onset
  }
  
  if (!is.null(data$center)) {
    colour_list$Center <- final_colours$center
  }
  
  # Generate UMAP plots
  generate_umap_plots(assay_data, output_dir, condition, label_file_name, labels_list, colour_list, patients_label)
  
} else {
  stop("Please provide the correct condition information!")
}
