suppressMessages(library("RColorBrewer"))
set.seed(123)
output_dir = "External"
if (!file.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

protein_expression <- read.csv("Input/Non-Normalized Relative Abundance of Individual Proteins From Individual Samples.csv",check.names = FALSE)

# split protein Id into SwissProt_Tremble, UniProtAccession, UniProtName
suppressMessages(library("tidyverse"))

# Proteins that were defined with at least one unique peptide were included in the protein data set for further analysis
protein_expression_clean <- protein_expression %>%
  mutate(id_all_sp = str_extract_all(`Protein Group Accession`, "sp\\|[^;]+"),
         id_all_sp = sapply(id_all_sp, function(x) paste(x, collapse = ";")))
# filter protein_expression_clean which id_all_sp == "" or unique peptide < 1
protein_expression_clean <- protein_expression_clean %>%
  filter(id_all_sp != "", `Unique peptides` > 0)

protein_expression_clean <- protein_expression_clean %>%
  mutate(id_first = gsub(";.*", "", `id_all_sp`)) %>%
  tidyr::separate(id_first, into = c("SwissProt_Tremble", "UniProtAccession", "UniProtName"), sep = "\\|")

selected_columns <- grep("^(ALS.*-1|HC)", colnames(protein_expression_clean), value = TRUE)

abu_data <- protein_expression_clean %>%
  dplyr::filter(SwissProt_Tremble == "sp", UniProtName != "Biognosys_iRT") %>%
  dplyr::select(UniProtName,UniProtAccession,id_all_sp, selected_columns)

# UniProtAccession has a suffix that needs to be removed
abu_data$UniProtAccession <- gsub("-.*", "", abu_data$UniProtAccession)

#===========================Protein Name Cleaning=============================
message("Reading protein to gene mapping...")
database = "Input/HUMAN_9606_idmapping_selected.tab.gz"
prot_to = readr::read_tsv(database, col_names = FALSE)
prot_to_gene <- prot_to[,c("X1","X2", "X19")]
colnames(prot_to_gene) <- c("UniProtAccession","UniProtName","gene_id")

# there are two proteins that changed accession number in uniprot database, need to be manually fixed.
#setdiff(abu_data$UniProtAccession, prot_to_gene$UniProtAccession)
# which abu_data$UniProtAccession == "P04220" need to be change into "P01871"
# which abu_data$UniProtAccession == "P01891" need to be change into "P04439"
abu_data$UniProtAccession[abu_data$UniProtAccession == "P04220"] <- "P01871"
abu_data$UniProtName[abu_data$UniProtAccession == "P04220"] <- "IGHM_HUMAN"
abu_data$UniProtAccession[abu_data$UniProtAccession == "P01891"] <- "P04439"
abu_data$UniProtName[abu_data$UniProtAccession == "P01891"] <- "HLAA_HUMAN"

# Find intersecting UniProt accessions between 'abu_data' and 'prot_to_gene'
inter = intersect(abu_data$UniProtAccession, prot_to_gene$UniProtAccession)

# according to UniProtAccession, remove duplicated rows
abu_data <- abu_data[!duplicated(abu_data$UniProtAccession), ]

# Filter 'prot_to_gene' to include only the matching UniProt accessions found in 'abu_data'
match_ids <- prot_to_gene[which(prot_to_gene$UniProtAccession %in% inter),]
# Sort 'match_ids' by 'UniProtAccession' to ensure it matches the order in 'abu_data'
match_ids = match_ids[match(abu_data$UniProtAccession,match_ids$UniProtAccession),]
# Remove the "_HUMAN" suffix from 'UniProtName' to obtain a clean 'gene_name'
match_ids$gene_name = gsub("_HUMAN","",match_ids$UniProtName)

# define gene ID use for search
gene_ids <- match_ids$UniProtAccession
suppressMessages(library("biomaRt"))
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_info = getBM(
  attributes=c('ensembl_gene_id','uniprotswissprot','hgnc_symbol'), 
  mart = mart) 

# rename specific columns
colnames(gene_info)[colnames(gene_info) == "hgnc_symbol"] <- "SYMBOL"
colnames(gene_info)[colnames(gene_info) == "uniprotswissprot"] <- "UniProtAccession"

# keep only unique uniprot accession numbers
gene_info = gene_info %>% distinct(gene_info$UniProtAccession, .keep_all = TRUE)

df = left_join(match_ids, gene_info, by = "UniProtAccession" )

# if df$SYMBOL is "", replace it into NA
df$SYMBOL[df$SYMBOL == ""] = NA
# if df$SYMBOL is not NA, use SYMBOL instead of gene_name
df$gene_name[!is.na(df$SYMBOL)] = df$SYMBOL[!is.na(df$SYMBOL)]

#rename SYMBOL in df into better_gene_name
uniprot_to_genename <- df %>%
  dplyr::rename(better_gene_name = SYMBOL) %>%
  dplyr::select(ensembl_gene_id,gene_name,UniProtAccession,better_gene_name,UniProtName)

colnames(uniprot_to_genename)[colnames(uniprot_to_genename) == "ensembl_gene_id"] <- "gene_id"

# make rows in uniprot_to_genename the same order as in abu_data
uniprot_to_genename = uniprot_to_genename[match(abu_data$UniProtAccession,uniprot_to_genename$UniProtAccession),]

# add gene_name and uniprot_accession to abu_data
abu_data$name = uniprot_to_genename$gene_name
abu_data$name = make.unique(abu_data$name)
abu_data$ID = uniprot_to_genename$UniProtAccession
abu_data$ID = make.unique(abu_data$ID)

message("Update uniprot to genename")
saveRDS(uniprot_to_genename, file.path(output_dir, "uniprot_to_genename.rds"))

write_delim(uniprot_to_genename, file.path(output_dir, "uniprot_to_genename.txt"), delim = "\t")

saveRDS(abu_data, file.path(output_dir, "abu_data_cleaned.rds"))

#make clinical data
clin = as.data.frame(selected_columns) %>% 
  dplyr::rename(label = selected_columns) %>%
  mutate(condition = ifelse(grepl("^ALS.*-1$", label), "als", "ctrl"))

experimental.design = clin
# Add a 'replicate' column to the experimental.design data frame
experimental.design <- experimental.design %>%
  group_by(condition) %>%  # Group by 'condition'
  mutate(replicate = row_number()) %>%  # Number each occurrence within group
  ungroup() 

abu_data = abu_data[,c("UniProtName",experimental.design$label,"name","ID")]
# get abundance column numbers
abundance.columns <- which(colnames(abu_data) %in% experimental.design$label)
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("DEP"))
#construct the SE (summarized experiment)
se_abu_data <- make_se(abu_data, abundance.columns, experimental.design) 

saveRDS(se_abu_data, file.path(output_dir, "se_abu_data.rds"))

#save the missing heatmap
suppressMessages(library("visdat"))
vis_miss(as.data.frame(assay(se_abu_data)) ,show_perc = TRUE, show_perc_col = TRUE, cluster = T) 
ggsave(file.path(output_dir, "missing_vis_miss_heatmap_before.pdf"), width = 11, height = 8, units = "in") 

threshold = 0.5
se_abu_data_filtered = DEP::filter_proteins(se_abu_data,type = "fraction", min = threshold)
saveRDS(se_abu_data_filtered, file.path(output_dir, "se_abu_data_filtered.rds"))
vis_miss(as.data.frame(assay(se_abu_data_filtered)),show_perc = TRUE, show_perc_col = TRUE, cluster = T)
ggsave(file.path(output_dir, "missing_vis_miss_heatmap_after.pdf"), width = 11, height = 8, units = "in")

# Calculate the number of missing values (0 or NA) per sample (column)
missing_per_sample <- abu_data %>%
  dplyr::select(-UniProtName) %>%  # Exclude 'UniProtName' from calculations
  summarise_all(~ sum(. == 0 | is.na(.))) %>%  # Count 0s and NAs in each column
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Missing_Count") %>%  # Convert to long format
  mutate(
    Missing_Percentage = (Missing_Count / nrow(abu_data)) * 100  # Calculate the percentage of missing values
  )

# Find samples with more than 50% missing values
#missing_per_sample[which(missing_per_sample$Missing_Percentage>50),]
write.csv(missing_per_sample, file.path(output_dir, "missing_per_sample.csv"),row.names = FALSE)

# Calculate the number and percentage of missing values per protein (row)
num_numeric_cols <- abu_data %>%
  dplyr::select(where(is.numeric)) %>%
  ncol()

# Calculate the number and percentage of missing values per protein (row)
missing_summary <- abu_data %>%
  rowwise() %>%  # Row-wise operations
  mutate(
    Missing_Count = sum(c_across(where(is.numeric)) == 0 | is.na(c_across(where(is.numeric)))),  # Count 0s and NAs in each row for numeric columns
    Missing_Percentage = (Missing_Count / num_numeric_cols) * 100  # Use pre-calculated number of numeric columns
  ) %>%
  ungroup() %>%  # Ungroup to return a regular tibble
  dplyr::select(UniProtName,name,ID, Missing_Count, Missing_Percentage)

write.csv(missing_summary, file.path(output_dir, "missing_per_protein.csv"),row.names = FALSE)


# 1. Reshape `abu_data` into long format, excluding `UniProtName`, `name`, and `ID` columns
abu_data_long <- abu_data %>%
  pivot_longer(cols = -c(UniProtName, name, ID), names_to = "label", values_to = "value")

# 2. Join with `experimental.design` to add condition information
abu_data_with_conditions <- abu_data_long %>%
  left_join(experimental.design, by = "label")

# 3. Group by both `UniProtName` and `condition`, and count zeros in each group
zero_counts_all <- abu_data_with_conditions %>%
  group_by(UniProtName, condition) %>%
  summarise(zero_count = sum(value == 0), .groups = "drop")  # Count zeros in each group

# Merge zero_counts_all into missing_summary
missing_summary_combined <- missing_summary %>%
  left_join(zero_counts_all, by = "UniProtName")

write.csv(missing_summary_combined, file.path(output_dir, "missing_per_protein_with_conditions.csv"),row.names = FALSE)

# Normalization
set.seed(123)
scalar = 0.01
norm <- normalize_vsn(se_abu_data_filtered)

pdf(file.path(output_dir, "meanSdPlot.pdf"))
meanSdPlot(norm)
dev.off()

pdf(file.path(output_dir, "Normalization_Boxplot.pdf"),height=12)
plot_normalization(se_abu_data_filtered,norm)
dev.off()

# Imputation
norm_imp_MinProb <- impute(norm, fun = "MinProb", q= scalar)
saveRDS(norm, file.path(output_dir, "norm.rds"))
saveRDS(norm_imp_MinProb, file.path(output_dir, "norm_imp_MinProb.rds"))

# Differential expression
suppressMessages(library("limma"))
fomula = formula(~0 + condition)
t = test_diff(norm_imp_MinProb, type = "control", control = "ctrl", #this function performs the tests and stores the results in 't'
              test = NULL, design_formula =fomula)

res = as.data.frame(t@elementMetadata@listData)
res$fdr = p.adjust(res[,grep("p.val",colnames(res))], method="BH")

res <- left_join(res, 
                      uniprot_to_genename %>% dplyr::select(UniProtAccession, gene_name), 
                      by = c("ID" = "UniProtAccession")) #join the gene names to the results
title = "als_vs_ctrl"
write.csv(res,quote = F,row.names = F, file = paste0(output_dir,"/", title, ".csv"))

#########################Function##################
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
  slow_vital_capacity = brewer.pal(n = 9, "Reds")
)

#give the vectors with discrete variables names
names(final_colours$volcano_plot) =  c("up", "down")
names(final_colours$male_female) = c("Male", "Female")
names(final_colours$clustering) = c("alpha", "beta", "theta")
names(final_colours$disease_status) = c("als", "ctrl")
names(final_colours$onset) = c("spinal", "bulbar", "ctrl")
names(final_colours$genetic_testing) = c("not_performed", "negative", "C9orf72")
names(final_colours$center) = c("goettingen", "munich")

#function to calculate total within sum of squares
calc_SS = function(df) sum(as.matrix(dist(df)^2)) / (2 * nrow(df))
calc_TWSS = function(df, clusters){
  number_clusters = length(levels(as.factor(clusters)))
  sum_of_squares = c()
  for(i in 1:number_clusters){sum_of_squares[i] = calc_SS(df[clusters == i,])}
  return(sum(sum_of_squares))
}

#function to calculate BIC and AIC
BIC2 <- function(df, clusters){
  m = ncol(df)
  n = nrow(df)
  k = length(levels(as.factor(clusters)))
  D = calc_TWSS(df, clusters)
  return(data.frame(AIC = D + 2*m*k,
                    BIC = D + log(n)*m*k))
}

perform_clustering <- function(assay_data, seed, clin_labels) {
  set.seed(seed)
  
  silhouette_scores = TWSS_scores = AIC_scores = BIC_scores = as.data.frame(matrix()) #create empty dataframes for the different scores
  rownames(silhouette_scores) = rownames(TWSS_scores) = rownames(AIC_scores) = rownames(BIC_scores) = "hclust" #they only have one row and first row will be used for hierarchical clustering results      
  
  cluster_results <- list()
  cluster_assignments <- as.data.frame(clin_labels)
  colnames(cluster_assignments) <- "patid"
  
  cluster_methods <- c("hclust", "mclust", "kmeans", "pam")
  
  for (i in 2:10) {
    
    # Perform clustering
    #hierarchical clustering
    title = paste0("hclust_k=", i)
    dist_mat <- dist(assay_data, method = 'euclidean')
    cl <- hclust(dist_mat, method = 'ward.D')
    cluster_assignments[, title] <- cutree(cl, k = i)
    #cluster fit measures
    ss = silhouette(cluster_assignments[,title], dist(assay))
    silhouette_scores["hclust", i] = mean(ss[, 3])
    TWSS_scores["hclust", i] = calc_TWSS(assay, cluster_assignments[,title])
    AIC_scores["hclust", i] = BIC2(assay, cluster_assignments[,title])[1,1]
    BIC_scores["hclust", i] = BIC2(assay, cluster_assignments[,title])[1,2]
    
    #model based clustering
    title = paste0("mclust_k=", i)
    library(mclust)
    cl <- Mclust(assay_data, G = i)
    cluster_assignments[, title] <- cl$classification
    
    ss = silhouette(cluster_assignments[,title], dist(assay))
    silhouette_scores["mclust", i] = mean(ss[, 3])
    TWSS_scores["mclust", i] = calc_TWSS(assay, cluster_assignments[,title])
    AIC_scores["mclust", i] = BIC2(assay, cluster_assignments[,title])[1,1]
    BIC_scores["mclust", i] = BIC2(assay, cluster_assignments[,title])[1,2]   
    
    #kmeans clustering
    title = paste0("kmeans_k=", i)
    cl <- kmeans(assay_data, centers = i, nstart=25,iter.max = 50)
    cluster_assignments[, title] <- cl$cluster
    
    ss = silhouette(cluster_assignments[,title], dist(assay))
    silhouette_scores["kmeans", i] = mean(ss[, 3])
    TWSS_scores["kmeans", i] = calc_TWSS(assay, cluster_assignments[,title])
    AIC_scores["kmeans", i] = BIC2(assay, cluster_assignments[,title])[1,1]
    BIC_scores["kmeans", i] = BIC2(assay, cluster_assignments[,title])[1,2]
    
    #partitioning around medoids
    title = paste0("pam_k=", i)
    cl <- pam(assay_data, k = i)
    cluster_assignments[, title] <- cl$clustering
    
    ss = silhouette(cluster_assignments[,title], dist(assay))
    silhouette_scores["pam", i] = mean(ss[, 3])
    TWSS_scores["pam", i] = calc_TWSS(assay, cluster_assignments[,title])
    AIC_scores["pam", i] = BIC2(assay, cluster_assignments[,title])[1,1]
    BIC_scores["pam", i] = BIC2(assay, cluster_assignments[,title])[1,2]
  }
  
  cluster_assignments[, 2:ncol(cluster_assignments)] <- cluster_assignments[, 2:ncol(cluster_assignments)] - 1
  return(list(
    cluster_assignments = cluster_assignments,
    AIC_scores = AIC_scores,
    BIC_scores = BIC_scores,
    silhouette_scores = silhouette_scores,
    TWSS_scores = TWSS_scores
  ))
} 




# Define a function to process each comparison dataframe
process_comparison <- function(comparison_df, comparison_name, GO_list, filtered_mapped_matrix, topn, output_dir) {
  
  # Select top terms from the comparison dataframe
  comparison_topn_up <- comparison_df %>% dplyr::slice_max(n = (topn+20)/2,order_by = t)  %>% as.data.frame()
  comparison_topn_down <- comparison_df %>% dplyr::slice_min(n = (topn+20)/2,order_by = t)  %>% as.data.frame()
  comparison_topn_df = rbind(comparison_topn_up,comparison_topn_down)
  
  # Find common GO terms between GO_list and gsva_result
  common_GO_terms <- intersect(names(GO_list), rownames(comparison_topn_df))
  
  # Subset GO_list to keep only the terms present in common_GO_terms
  GO_list_subset <- GO_list[names(GO_list) %in% common_GO_terms]
  
  # Extract the rownames of the assay genes
  assay_genes <- rownames(filtered_mapped_matrix)
  
  # Subset genes in GO_list_subset to keep only those present in assay_genes
  GO_list_subset_filtered <- lapply(GO_list_subset, function(genes) {
    intersect(genes, assay_genes)
  })
  
  # Get top terms and genes
  GO_list_top_terms <- GO_list_subset_filtered[names(GO_list_subset_filtered) %in% common_GO_terms]
  
  # Convert GO_list_top_terms to a data frame
  go_df <- enframe(GO_list_top_terms, name = "terms", value = "genes") %>%
    mutate(genes = sapply(genes, function(gene_list) paste(gene_list, collapse = ", ")))
  
  # Create a data frame with the comparison_df terms and their corresponding genes
  comparison_topn <- comparison_topn_df %>%
    rownames_to_column("terms") %>%
    left_join(go_df, by = "terms") %>%
    dplyr::select(terms, genes, everything())
  
  # Combine the genes for each term
  comparison_topn$genes <- sapply(comparison_topn$genes, function(x) paste(x, collapse = ","))
  
  # Write the result into an Excel file
  write_xlsx(comparison_topn %>% as.data.frame(),
             path = file.path(output_dir, paste0("GSVA_topn_terms_", comparison_name, "_with_matched_proteins.xlsx")))
  
  return(comparison_topn)
}

# Custom function to abbreviate long terms by adding "..." at the end
abbreviate_terms_end <- function(term, max_length = 60) {
  if (nchar(term) > max_length) {
    paste0(substr(term, 1, max_length - 3), "...")
  } else {
    term
  }
}



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

#Helper Function to Generate UMAP Plots
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


#########################
suppressMessages(library("cluster"))
suppressMessages(library("mclust"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggthemes"))
suppressMessages(library("gridExtra"))
suppressMessages(library("ggpubr"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("ggalluvial"))
####### Clustering
seed = 123
se = norm_imp_MinProb
assay = as.data.frame(assay(se)) #make the abundancy matrix into a dataframe
#assay = assay[,clin$ID] #order patients in assay in same order as clinical table
colnames(assay) = clin$label #give them CSF labels again
#columns should be variables so we need to transpose the table
assay = as.data.frame(t(assay))

#####we only want to use ALS samples for clustering
ALS_data = se@colData %>% as.data.frame() %>% dplyr::filter(condition == "als")
assay = assay[ALS_data$label,]

#Perform different types of clustering
#use clustering function to get cluster labels
cluster_assignments_result <- perform_clustering(assay,seed = seed, clin_labels = ALS_data$label)
cluster_assignments <- cluster_assignments_result$cluster_assignments
AIC_scores = cluster_assignments_result$AIC_scores
BIC_scores = cluster_assignments_result$BIC_scores
silhouette_scores = cluster_assignments_result$silhouette_scores
TWSS_scores = cluster_assignments_result$TWSS_scores
saveRDS(cluster_assignments, file.path(output_dir, "cluster_assignments.rds"))

# fit measures
scores = list(AIC = AIC_scores, BIC = BIC_scores, silhouette = silhouette_scores, TWSS = TWSS_scores) 
plots = list() #initiate empty plot list

for(i in 1:length(scores)){ #loop that iterates through all types of scores
  colnames(scores[[i]]) = 1:10
  scores[[i]]$method = rownames(scores[[i]])
  melt = reshape2::melt(scores[[i]])
  
  plots[[i]] = ggplot(melt, aes(x=variable, y=value, group = method, colour = method)) +
    geom_line() +
    geom_point() +
    ggtitle(names(scores)[i]) +
    theme_few()
  
}
names(plots) = names(scores)

allplots <- ggarrange(plotlist=plots,
                      labels = LETTERS[1:length(plots)],
                      ncol = 2, nrow = 2)

# Save the plots
ggsave(file.path(output_dir,"fit_scores.pdf"), width = 11, height = 8, units = "in")

# Sankey plots for all methods (Sankey_plots_within_method.pdf)
cluster_assignments <- cluster_assignments_result$cluster_assignments
number_of_clusters = as.character(2:5)
plots = list()

for(i in 1:length(number_of_clusters)){
  data = cluster_assignments[,c(1, grep(number_of_clusters[i], colnames(cluster_assignments)))]
  data = reshape2::melt(data)
  colnames(data) = c("patid", "method", "cluster")
  data = as.data.frame(lapply(data, as.factor))
  
  plots[[i]] = ggplot(data,
                      aes(x = method, stratum = cluster, alluvium = patid,
                          fill = cluster, label = cluster)) +
    scale_fill_brewer(type = "qual", palette = "Set3") +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
              color = "darkgray") +
    geom_stratum() +
    theme(legend.position = "bottom") +
    ggtitle(paste0("How patients are clustered in different clustering algorithms \nk=", number_of_clusters[i])) +
    theme_few()
}

# Sankey plots within a method, different number of clusters
methods = c("kmeans", "hclust", "mclust", "pam")
plots = list()
number_of_clusters

for(i in 1:length(methods)){
  data = cluster_assignments[,c("patid", paste0(methods[i],"_k=",number_of_clusters))]
  
  data = reshape2::melt(data, id = "patid")
  colnames(data) = c("patid", "method", "cluster")
  data = as.data.frame(lapply(data, as.factor))
  
  plots[[i]] = ggplot(data,
                      aes(x = method, stratum = cluster, alluvium = patid,
                          fill = cluster, label = cluster)) +
    scale_fill_brewer(type = "qual", palette = "Set3") +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
              color = "darkgray") +
    geom_stratum() +
    theme(legend.position = "bottom") +
    ggtitle(paste0("How patients are clustered within one clustering algorithm \nmethod = ", methods[i])) +
    theme_few() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

allplots <- ggarrange(plotlist=plots,
                      labels = 1:length(plots),
                      ncol = 2, nrow = 2)

# Save the plots
ggsave(file.path(output_dir, "Sankey_plots_within_method.pdf"), width = 11*2, height = 8*2, units = "in") 

## Sankey plots for kmeans k2 and k3
cluster_assignments_2 = cluster_assignments[,c("patid", "kmeans_k=2", "kmeans_k=3")]
#set the levels for the clusters
#levels k=2
# Create a table to count the objects in each cluster
cluster_counts <- table(cluster_assignments_2$"kmeans_k=2")
# Identify the larger and smaller cluster
larger_cluster <- names(cluster_counts)[which.max(cluster_counts)]
smaller_cluster <- names(cluster_counts)[which.min(cluster_counts)]
# Dynamically assign "alpha" to the larger cluster and "beta" to the smaller cluster
cluster_assignments_2$`kmeans_k=2` = as.factor(cluster_assignments_2$`kmeans_k=2`)
levels(cluster_assignments_2$`kmeans_k=2`) <- ifelse(levels(cluster_assignments_2$`kmeans_k=2`) == larger_cluster, "alpha", "beta")
cluster_assignments_2$`kmeans_k=2` = ordered(cluster_assignments_2$`kmeans_k=2`, levels = c("alpha", "beta", "theta"))

message("kmeans = 2:")
print(table(cluster_assignments_2$`kmeans_k=2`))

#levels k=3
# If “alpha” is the majority, assign “alpha”.
# If “beta” is the majority, assign “beta”.
# If there is a mix (or tie), assign “theta”.
adjusted_clusters <- cluster_assignments_2 %>%
  group_by(`kmeans_k=3`) %>%
  summarize(
    majority_alpha = sum(`kmeans_k=2` == "alpha"),
    majority_beta = sum(`kmeans_k=2` == "beta"),
    total = n()
  ) %>%
  mutate(
    alpha_ratio = majority_alpha / total,
    beta_ratio = majority_beta / total,
    alpha_beta = alpha_ratio - beta_ratio
  ) %>%
  # Assign based on alpha_beta values
  mutate(
    kmeans_k3_adjusted = case_when(
      alpha_beta == max(alpha_beta) ~ "alpha",  # Highest alpha_beta
      alpha_beta == min(alpha_beta) ~ "beta",   # Lowest alpha_beta
      TRUE ~ "theta"                            # Remaining one
    )
  )

adjusted_clusters = adjusted_clusters[,c("kmeans_k=3","kmeans_k3_adjusted")]
cluster_assignments_2 = left_join(cluster_assignments_2,adjusted_clusters, by = "kmeans_k=3") %>% dplyr::select(-`kmeans_k=3`) 
colnames(cluster_assignments_2)[3] = "kmeans_k=3"

cluster_assignments_2$`kmeans_k=3` = as.factor(cluster_assignments_2$`kmeans_k=3`)
cluster_assignments_2$`kmeans_k=3` <- ordered(cluster_assignments_2$`kmeans_k=3`, levels = c("alpha", "beta", "theta"))

write.csv(cluster_assignments_2, file = file.path(output_dir, "cluster_assignments_2.csv"), row.names = FALSE)

message("kmeans = 2:")
print(table(cluster_assignments_2$`kmeans_k=3`))

#Sankey plot with k=2 and k=3
data = pivot_longer(cluster_assignments_2, cols = 2:dim(cluster_assignments_2)[2])
colnames(data) = c("patid", "k", "cluster")

data23 = data[data$k == "kmeans_k=2" | data$k == "kmeans_k=3",]


ggplot(data23,
       aes(x = k, stratum = cluster, alluvium = patid,
           fill = cluster, label = cluster)) +
  scale_fill_manual(values = final_colours$clustering) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom") +
  ggtitle("Patient Sankey plot within K-means, from 2 clusters to 3 clusters \nPATIENTS") +
  theme_few()


ggsave(file.path(output_dir, "Sankey_plot_kmeans_2_3_patients.pdf"), width = 11/2, height = 8/2, units = "in")

# Filter the data frame to find the patient who changed from "alpha" to "beta"
changed_cluster <- cluster_assignments_2 %>%
  filter(`kmeans_k=2` == "alpha" & `kmeans_k=3` == "beta")
print(changed_cluster)


###############################ssGSEA analyais################################
min_sz = 5
category = "C5"
subset = "GO:BP"
topn = 30
suppressMessages(library("msigdbr"))
suppressMessages(library("writexl"))
cluster_assignments = "External/cluster_assignments_2.csv"
cluster_assignments = read.csv(cluster_assignments,check.names = FALSE)
message(paste("msigdbr category used in GSVA is:",category))
GO_gene_sets <- msigdbr::msigdbr(
  species = "Homo sapiens", # Can change this to what species you need
  category = category # Only GO gene sets
) 
GO_gene_sets = GO_gene_sets[which(GO_gene_sets$gs_subcat ==subset),]

# Background removal
GO_gene_sets <- GO_gene_sets[GO_gene_sets$gene_symbol %in% rownames(se), ]

GO_list <- split(
  GO_gene_sets$gene_symbol, # The genes we want split into pathways
  GO_gene_sets$gs_name # The pathways made as the higher levels of the list
)

elementmeta = se@elementMetadata %>% as.data.frame()%>% dplyr::select(name,ID)
filtered_mapped_matrix = assay(se) %>% as.data.frame() %>% rownames_to_column("name") %>% left_join(elementmeta,by = join_by(name))

filtered_mapped_matrix = filtered_mapped_matrix %>% left_join(uniprot_to_genename %>% dplyr::select(UniProtAccession,gene_name),by = c("ID" = "UniProtAccession")) %>% dplyr::select(-ID,-name) 

filtered_mapped_matrix  = filtered_mapped_matrix %>%   column_to_rownames("gene_name") %>% as.matrix()

# Extract the row names from metadata according to the order of colnames in plot_df
colnames(filtered_mapped_matrix) = se@colData$label

message("Performing GSVA analysis")

# https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

# Extract gene names from filtered_mapped_matrix
genes_of_interest <- rownames(filtered_mapped_matrix)

# Function to calculate overlap with genes of interest
overlap_with_genes <- function(geneset) {
  length(intersect(geneset, genes_of_interest))
}

# Filter GO_list, keeping only sets with at least 10 overlapping genes
filtered_GO_list <- GO_list[sapply(GO_list, overlap_with_genes) >= 10]

### We only keep ALS samples
filtered_mapped_matrix = filtered_mapped_matrix[,ALS_data$label]

suppressMessages(library("GSVA"))
#https://rdrr.io/bioc/GSVA/man/filterGeneSets.html
gsva_results <- gsva(
  filtered_mapped_matrix,
  filtered_GO_list,
  method = "ssgsea",
  # Appropriate for our vst transformed data
  kcdf = "Gaussian",
  # Minimum gene set size
  min.sz = min_sz,
  # Maximum gene set size
  max.sz = 500,
  # Compute Gaussian-distributed scores
  mx.diff = TRUE,
  # Don't print out the progress bar
  verbose = FALSE
)

saveRDS(gsva_results, file.path(output_dir, "gsva_results.rds"))

# Define the Z-score normalization function
zscore_normalize <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Apply the Z-score normalization to each row (gene set)
gsva_results <- t(apply(gsva_results, 1, zscore_normalize))

# Create a heatmap of the GSVA results
metadata = ALS_data  %>% left_join(cluster_assignments,by = c("label" = "patid")) 
metadata = metadata %>% column_to_rownames("label")

####for kmeans k = 3
message("Performing heatmap for kmeans k=3")
annolabel <- data.frame(Cluster = metadata[,"kmeans_k=3"])
rownames(annolabel) <- rownames(metadata)

group_info <- factor(annolabel$Cluster,levels = c("alpha","beta","theta"))

# Create a design matrix for the linear model
design <- model.matrix(~ 0 + group_info)
colnames(design) <- levels(group_info)  # Name the columns according to the groups

# Fit the linear model
fit <- lmFit(gsva_results, design)

# Create contrasts for group comparisons
contrast_matrix <- makeContrasts(
  theta_vs_beta =  theta - beta,  # Comparison between Group 1 and Group 0
  alpha_vs_theta = alpha - theta,  # Comparison between Group 2 and Group 0
  alpha_vs_beta = alpha - beta,  # Comparison between Group 2 and Group 1
  levels = design
)

# Apply contrasts to the fit
fit2 <- contrasts.fit(fit, contrast_matrix)

# Compute empirical Bayes statistics
fit2 <- eBayes(fit2)

# Extract the top differentially expressed pathways for each contrast
theta_vs_beta <- topTable(fit2, coef = "theta_vs_beta", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE,p.value=0.05)
alpha_vs_theta <- topTable(fit2, coef = "alpha_vs_theta", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE,p.value=0.05)
alpha_vs_beta <- topTable(fit2, coef = "alpha_vs_beta", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE,p.value=0.05)


write.csv(theta_vs_beta %>%  as.data.frame()%>%rownames_to_column("terms") ,quote = F,row.names = F, file = file.path(output_dir, "GSVA_theta_vs_beta_k3.csv"))
write.csv(alpha_vs_theta %>% as.data.frame()%>%rownames_to_column("terms"),quote = F,row.names = F, file = file.path(output_dir, "GSVA_alpha_vs_theta_k3.csv"))
write.csv(alpha_vs_beta %>% as.data.frame()%>%rownames_to_column("terms"),quote = F,row.names = F, file = file.path(output_dir, "GSVA_alpha_vs_beta_k3.csv"))

############
# Process alpha_vs_beta, theta_vs_beta, and alpha_vs_theta
alpha_vs_beta_topn <- process_comparison(alpha_vs_beta, "alpha_vs_beta_k3", GO_list, filtered_mapped_matrix, topn, output_dir)
theta_vs_beta_topn <- process_comparison(theta_vs_beta, "theta_vs_beta_k3", GO_list, filtered_mapped_matrix, topn, output_dir)
alpha_vs_theta_topn <- process_comparison(alpha_vs_theta, "alpha_vs_theta_k3", GO_list, filtered_mapped_matrix, topn, output_dir)
############
theta_vs_beta_up = theta_vs_beta %>% dplyr::slice_max(n = topn/2,order_by = t)  %>% as.data.frame()
theta_vs_beta_down = theta_vs_beta %>% dplyr::slice_min(n = topn/2,order_by = t)  %>% as.data.frame()
theta_vs_beta = rbind(theta_vs_beta_up,theta_vs_beta_down)

# Split alpha_vs_theta into top and bottom halves based on t values
alpha_vs_theta_up = alpha_vs_theta %>% dplyr::slice_max(n = topn/2, order_by = t) %>% as.data.frame()
alpha_vs_theta_down = alpha_vs_theta %>% dplyr::slice_min(n = topn/2, order_by = t) %>% as.data.frame()
alpha_vs_theta = rbind(alpha_vs_theta_up, alpha_vs_theta_down)

# Split alpha_vs_beta into top and bottom halves based on t values
alpha_vs_beta_up = alpha_vs_beta %>% dplyr::slice_max(n = topn/2, order_by = t) %>% as.data.frame()
alpha_vs_beta_down = alpha_vs_beta %>% dplyr::slice_min(n = topn/2, order_by = t) %>% as.data.frame()
alpha_vs_beta = rbind(alpha_vs_beta_up, alpha_vs_beta_down)

#merge theta_vs_beta, alpha_vs_theta, alpha_vs_beta, and also add the contrast column
topn_terms = rbind(theta_vs_beta %>% mutate(contrast = "theta_vs_beta"),
                   alpha_vs_theta %>% mutate(contrast = "alpha_vs_theta"),
                   alpha_vs_beta %>% mutate(contrast = "alpha_vs_beta"))
topn_terms = topn_terms %>% rownames_to_column("terms")
write.csv(topn_terms,quote = F,row.names = F, file = file.path(output_dir, "GSVA_topn_terms_k3.csv"))

# Order the columns of 'plot_df' according to 'annolabel'
ordered_columns <- rownames(annolabel)[order(annolabel$Cluster)]
cluster_colors <- list(Cluster = c("alpha" = "#FC8D62","beta" = "#8DA0CB", "theta" = "#A6D854"))



plot_df <- gsva_results[unique(c(rownames(theta_vs_beta),rownames(alpha_vs_theta),rownames(alpha_vs_beta))), ordered_columns]


legend_breaks <- seq(from = round(min(plot_df), 1), to = round(max(plot_df), 1), by = 0.2)

# Extract the row names from metadata according to the order of colnames in plot_df

rownames(plot_df) <- sapply(rownames(plot_df), abbreviate_terms_end)

pathway_heatmap <- pheatmap::pheatmap(plot_df,
                                      annotation_col = annolabel,
                                      annotation_colors = cluster_colors,  
                                      fontsize_row = 6,fontsize_col = 6, cluster_cols = FALSE, color = colorRampPalette(c("navy", "white", "firebrick3"))(100),legend_breaks = legend_breaks, border_color = NA)


pdf(file.path(output_dir, "GSVA_sample_level_k3.pdf"),width = 10)
print(pathway_heatmap)
dev.off()

### cluster means heatmap k = 3
cluster_means <- gsva_results %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  mutate(Cluster = annolabel$Cluster) %>%
  group_by(Cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "Cluster") %>%
  t()

write_cluster_means = cluster_means %>% as.data.frame() %>%rownames_to_column("terms")
write.csv(write_cluster_means,quote = F,row.names = F, file = file.path(output_dir, "GSVA_cluster_means_k3.csv"))



cluster_means = cluster_means[unique(c(rownames(theta_vs_beta),rownames(alpha_vs_theta),rownames(alpha_vs_beta))),]

# Custom function to abbreviate long terms by adding "..." at the end
rownames(cluster_means) <- sapply(rownames(cluster_means), abbreviate_terms_end)
pathway_heatmap_mean <- pheatmap::pheatmap(cluster_means,
                                           fontsize_row = 6, cluster_cols = FALSE, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), border_color = NA)


pdf(file.path(output_dir, "GSVA_mean_k3.pdf"),width = 7)
print(pathway_heatmap_mean)
dev.off()


#########################################################
####for kmeans k = 2
message("Performing heatmap for kmeans k=2")
annolabel <- data.frame(Cluster = metadata[,"kmeans_k=2"])
rownames(annolabel) <- rownames(metadata)

group_info <- factor(annolabel$Cluster,levels = c("alpha","beta"))

# Create a design matrix for the linear model
design <- model.matrix(~ 0 + group_info)
colnames(design) <- levels(group_info)  # Name the columns according to the groups

# Fit the linear model
fit <- lmFit(gsva_results, design)

# Create contrasts for group comparisons
contrast_matrix <- makeContrasts(
  alpha_vs_beta = alpha - beta,  # Comparison between alpha and beta
  levels = design
)

# Apply contrasts to the fit
fit2 <- contrasts.fit(fit, contrast_matrix)

# Compute empirical Bayes statistics
fit2 <- eBayes(fit2)

# Extract the top differentially expressed pathways for each contrast
alpha_vs_beta <- topTable(fit2, coef = "alpha_vs_beta", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE,p.value=0.05)

write.csv(alpha_vs_beta %>% as.data.frame()%>%rownames_to_column("terms"),quote = F,row.names = F, file = file.path(output_dir, "GSVA_alpha_vs_beta_k2.csv"))

##########
#write into excel
alpha_vs_beta_topn <- process_comparison(alpha_vs_beta, "alpha_vs_beta_k2", GO_list, filtered_mapped_matrix, topn, output_dir)
#########

alpha_vs_beta_up = alpha_vs_beta %>% dplyr::slice_max(n = topn/2,order_by = t)  %>% as.data.frame()
alpha_vs_beta_down = alpha_vs_beta %>% dplyr::slice_min(n = topn/2,order_by = t)  %>% as.data.frame()
alpha_vs_beta = rbind(alpha_vs_beta_up,alpha_vs_beta_down)

write.csv(alpha_vs_beta %>% as.data.frame()%>%rownames_to_column("terms"),quote = F,row.names = F, file = file.path(output_dir, "GSVA_topn_terms_k2.csv"))

# Order the columns of 'test' according to 'annolabel'
ordered_columns <- rownames(annolabel)[order(annolabel$Cluster)]
cluster_colors <- list(Cluster = c("alpha" = "#FC8D62","beta" = "#8DA0CB"))



plot_df <- gsva_results[unique(rownames(alpha_vs_beta)), ordered_columns]


legend_breaks <- seq(from = round(min(plot_df), 1), to = round(max(plot_df), 1), by = 0.2)

# Extract the row names from metadata according to the order of colnames in plot_df
rownames(plot_df) <- sapply(rownames(plot_df), abbreviate_terms_end)
pathway_heatmap_k2 <- pheatmap::pheatmap(plot_df,
                                         annotation_col = annolabel,
                                         annotation_colors = cluster_colors,  
                                         fontsize_row = 6,fontsize_col = 6, cluster_cols = FALSE, color = colorRampPalette(c("navy", "white", "firebrick3"))(100),legend_breaks = legend_breaks, border_color = NA)


pdf(file.path(output_dir, "GSVA_sample_level_k2.pdf"),width = 10,height = 4)
print(pathway_heatmap_k2)
dev.off()


### cluster means heatmap k = 2
cluster_means <- gsva_results %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  mutate(Cluster = annolabel$Cluster) %>%
  group_by(Cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "Cluster") %>%
  t()

write_cluster_means = cluster_means %>% as.data.frame() %>%rownames_to_column("terms")
write.csv(write_cluster_means,quote = F,row.names = F, file = file.path(output_dir, "GSVA_cluster_means_k2.csv"))

cluster_means = cluster_means[unique(rownames(alpha_vs_beta)),]

# Custom function to abbreviate long terms by adding "..." at the end
rownames(cluster_means) <- sapply(rownames(cluster_means), abbreviate_terms_end)

pathway_heatmap_mean_2 <- pheatmap::pheatmap(cluster_means,
                                             fontsize_row = 6, cluster_cols = FALSE, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), border_color = NA)


pdf(file.path(output_dir, "GSVA_mean_k2.pdf"),width = 6,height = 4)
print(pathway_heatmap_mean_2)
dev.off()


#######Make umap plot
assay_data = t(filtered_mapped_matrix)
patients_label = NULL
label_file_name = "No_label"

# Define labels
labels_list <- list()

# Make cluster_assignments the same order as data
cluster_assignments <- cluster_assignments[match(rownames(assay_data), cluster_assignments$patid),]
labels_list$`kmeans_k=2` <- as.factor(cluster_assignments$`kmeans_k=2`)
labels_list$`kmeans_k=3` <- as.factor(cluster_assignments$`kmeans_k=3`)

# Define color schemes
colour_list <- list(
  "kmeans_k=2" = final_colours$clustering[names(final_colours$clustering) %in% levels(labels_list$`kmeans_k=2`)],
  "kmeans_k=3" = final_colours$clustering[names(final_colours$clustering) %in% levels(labels_list$`kmeans_k=3`)]
)

condition = "als"
suppressMessages(library("umap"))
suppressMessages(library("cowplot"))
set.seed(123)
generate_umap_plots(assay_data, output_dir, condition, label_file_name, labels_list, colour_list, patients_label)


#################Differential expression analysis between alpha vs beta##########
clin_ALS = clin[which(clin$condition =="als"),]

abu_data_als = abu_data[,c("UniProtName",clin_ALS$label,"name","ID")]

clin_ALS = left_join(clin_ALS,cluster_assignments,by = c("label" = "patid"))

experimental.design = clin_ALS
experimental.design = experimental.design[,-which(names(experimental.design) =="condition")]
colnames(experimental.design)[which(colnames(experimental.design) =="kmeans_k=2")] = "condition"

# Add a 'replicate' column to the experimental.design data frame
experimental.design <- experimental.design %>%
  group_by(condition) %>%  # Group by 'condition'
  mutate(replicate = row_number()) %>%  # Number each occurrence within group
  ungroup() 


# get abundance column numbers
abundance.columns <- which(colnames(abu_data_als) %in% experimental.design$label)

#make se
se_abu_data_als <- make_se(abu_data_als, abundance.columns, experimental.design) 
saveRDS(se_abu_data_als, file.path(output_dir, "se_abu_data_als.rds"))

# filter proteins
threshold = 0.5
se_abu_data_filtered_als = DEP::filter_proteins(se_abu_data_als,type = "fraction", min = threshold)
saveRDS(se_abu_data_filtered_als, file.path(output_dir, "se_abu_data_filtered_als.rds"))

# Calculate the number of missing values (0 or NA) per sample (column)
missing_per_sample <- abu_data_als %>%
  dplyr::select(-UniProtName) %>%  # Exclude 'UniProtName' from calculations
  summarise_all(~ sum(. == 0 | is.na(.))) %>%  # Count 0s and NAs in each column
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Missing_Count") %>%  # Convert to long format
  mutate(
    Missing_Percentage = (Missing_Count / nrow(abu_data_als)) * 100  # Calculate the percentage of missing values
  )

# Find samples with more than 50% missing values
#missing_per_sample[which(missing_per_sample$Missing_Percentage>50),]
write.csv(missing_per_sample, file.path(output_dir, "missing_per_sample_als.csv"),row.names = FALSE)





# Calculate the number and percentage of missing values per protein (row)
num_numeric_cols <- abu_data_als %>%
  dplyr::select(where(is.numeric)) %>%
  ncol()

# Calculate the number and percentage of missing values per protein (row)
missing_summary <- abu_data_als %>%
  rowwise() %>%  # Row-wise operations
  mutate(
    Missing_Count = sum(c_across(where(is.numeric)) == 0 | is.na(c_across(where(is.numeric)))),  # Count 0s and NAs in each row for numeric columns
    Missing_Percentage = (Missing_Count / num_numeric_cols) * 100  # Use pre-calculated number of numeric columns
  ) %>%
  ungroup() %>%  # Ungroup to return a regular tibble
  dplyr::select(UniProtName,name,ID, Missing_Count, Missing_Percentage)

write.csv(missing_summary, file.path(output_dir, "missing_per_protein_als.csv"),row.names = FALSE)



# 1. Reshape `abu_data` into long format, excluding `UniProtName`, `name`, and `ID` columns
abu_data_long <- abu_data_als %>%
  pivot_longer(cols = -c(UniProtName, name, ID), names_to = "label", values_to = "value")

# 2. Join with `experimental.design` to add condition information
abu_data_with_conditions <- abu_data_long %>%
  left_join(experimental.design, by = "label")

# 3. Group by both `UniProtName` and `condition`, and count zeros in each group
zero_counts_all <- abu_data_with_conditions %>%
  group_by(UniProtName, condition) %>%
  summarise(zero_count = sum(value == 0), .groups = "drop")  # Count zeros in each group

# Merge zero_counts_all into missing_summary
missing_summary_combined <- missing_summary %>%
  left_join(zero_counts_all, by = "UniProtName")

write.csv(missing_summary_combined, file.path(output_dir, "missing_per_protein_with_conditions_als.csv"),row.names = FALSE)



# Normalization
set.seed(123)
scalar = 0.01
norm <- normalize_vsn(se_abu_data_filtered_als)

pdf(file.path(output_dir, "meanSdPlot.pdf"))
meanSdPlot(norm)
dev.off()

pdf(file.path(output_dir, "Normalization_Boxplot.pdf"),height=12)
plot_normalization(se_abu_data_filtered_als,norm)
dev.off()

# Imputation
norm_imp_MinProb <- impute(norm, fun = "MinProb", q= scalar)
saveRDS(norm, file.path(output_dir, "norm.rds"))
saveRDS(norm_imp_MinProb, file.path(output_dir, "norm_imp_MinProb_als.rds"))




# Differential expression
suppressMessages(library("limma"))
fomula = formula(~0 + condition)
t = test_diff(norm_imp_MinProb, type = "control", control = "beta", #this function performs the tests and stores the results in 't'
              test = NULL, design_formula =fomula)

res = as.data.frame(t@elementMetadata@listData)
res$fdr = p.adjust(res[,grep("p.val",colnames(res))], method="BH")

res <- left_join(res, 
                 uniprot_to_genename %>% dplyr::select(UniProtAccession, gene_name), 
                 by = c("ID" = "UniProtAccession")) #join the gene names to the results
title = "alpha_vs_beta"
write.csv(res,quote = F,row.names = F, file = paste0(output_dir,"/", title, "_als.csv"))


####Volcano plot between alpha vs beta
suppressMessages(library("ggrepel"))
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
    geom_text_repel(data = dplyr::filter(df, name %in% labels),
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



data_res = res
diff_index = grep("diff",colnames(data_res))
fdr_index = grep("fdr",colnames(data_res))


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
  dplyr::filter(x >= 0 & y >= (-log10(0.05))) %>%  # Apply filtering conditions
  arrange(desc(x)) %>%  # Sort by 'x' in descending order
  slice_head(n = 10)

top10_label_down <- df %>%
  dplyr::filter(x <= 0 & y >= (-log10(0.05))) %>%  # Apply filtering conditions
  arrange(x) %>%  # Sort by 'x' in descending order
  slice_head(n = 10)

label = NULL
plots_FDR0.05 = volcano_plot(df, 0.05 , paste0("Volcano plot clustering proteomics \nalpha = FDR 0.05\n", title),labels = c(top10_label_up$name,top10_label_down$name,label),case_color,ctrl_color)
plots_FDR0.1 = volcano_plot(df, 0.1 , paste0("Volcano plot clustering proteomics \nalpha = FDR 0.1\n", title),labels =c(top10_label_up$name,top10_label_down$name,label),case_color,ctrl_color)

volcano_plots = list(plots_FDR0.05, plots_FDR0.1)
ggarrange(plotlist = volcano_plots, ncol = 2, nrow =1 )
title = "volcano_plots_kmeans_alpha_vs_beta"
ggsave(filename = paste0(output_dir,"/", title, ".pdf"), width =12, height = 5, units = "in")