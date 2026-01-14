
#=========================Script Description=================================
# This script is used for summarized GSEA results heatmap DC & VC.

# GSEA_result_Discovery = readRDS("/Users/xliu2942/Documents/Projects/MAXOMOD/Script/Discovery/19_GSEA_allsamples/GSEA_result.rds")
# GSEA_result_Validation = readRDS("/Users/xliu2942/Documents/Projects/MAXOMOD/Script/Validation/19_GSEA_allsamples/GSEA_result.rds")

# GSEA_result_Discovery = GSEA_result_Discovery$BP@result
# GSEA_result_Validation = GSEA_result_Validation$BP@result

# ego_list = list(GSEA_result_Discovery, GSEA_result_Validation)
# names(ego_list) = c("Discovery", "Validation")
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("msigdbr"))
suppressMessages(library("circlize"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))
suppressMessages(library("ggsci"))
suppressMessages(library("tidyverse"))
suppressMessages(library("GOSemSim"))
suppressMessages(library("dendextend"))
suppressMessages(library("ComplexHeatmap"))
#===========================Functions=================================
# 解析传入字符串为命名向量
parse_named_paths <- function(input_string) {
  # 拆分 name=path 部分
  parts <- strsplit(input_string, ",")[[1]]
  split_parts <- strsplit(parts, "=")
  names_vec <- sapply(split_parts, `[[`, 1)
  paths_vec <- sapply(split_parts, `[[`, 2)
  named_paths <- setNames(paths_vec, names_vec)
  return(named_paths)
}

read_gsea_results <- function(path_vector, ontology = "BP") {
  ego_list <- lapply(path_vector, function(p) {
    gsea <- readRDS(p)
    gsea[[ontology]]@result
  })
  names(ego_list) <- names(path_vector)
  return(ego_list)
}

# ===========================Command Parameters Setting=============================

# 设置参数
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "Discovery=/Users/xliu2942/Documents/Projects/MAXOMOD/Script/Discovery/19_GSEA_allsamples/GSEA_result.rds,Validation=/Users/xliu2942/Documents/Projects/MAXOMOD/Script/Validation/19_GSEA_allsamples/GSEA_result.rds",
              help = "Named GSEA result RDS files, format: name1=path1,name2=path2,..."),
  make_option(c("--output", "-o"),
              type = "character", default = "36_GSEA_IC_heatmap",
              help = "Output directory"
  ),make_option(c("--cluster_number", "-c"),
                type = "integer", default = 6,
                help = "Number of clusters"
  ),make_option(c("--seed", "-e"),
                type = "integer", default = 9,
                help = "set.seed"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
#============================================================================
if (is.null(opt$input)) {
  stop("Please provide the GSEA results!")
}else{
  gsea_paths <- parse_named_paths(opt$input)
  ego_list <- read_gsea_results(gsea_paths, ontology = "BP")
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

if (is.null(opt$cluster_number)) {
  k_num <- 6
} else {
  k_num <- opt$cluster_number
}

if (is.null(opt$seed)) {
  set.seed(9)
} else {
  set.seed(opt$seed)
}

#============================================================================
# Set the cutoff for adjusted p-value
padj_cutoff <- 0.05

# keep only significant GO terms
all_go_terms <- unique(unlist(
  lapply(ego_list, function(x) {
    sig <- x
    sig <- sig[sig$p.adjust < padj_cutoff, ]
    sig$ID
  })
))

# build GO term x Comparison matrix
logp_mat <- matrix(0, nrow = length(all_go_terms), ncol = length(ego_list))
rownames(logp_mat) <- all_go_terms

# 获取背景集（你原来用于 GSEA 的）
bg_map <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  distinct(gs_name, gs_exact_source)  # gs_name 是 GOBP_xxx，gs_exact_source 是 GO:xxxxxxx

# 将 merged_common 的 Description 与 bg_map 合并，获得 GO ID
all_go_terms = as.data.frame(all_go_terms)
colnames(all_go_terms) <- "Description"
all_go_terms <- left_join(all_go_terms, bg_map, by = c("Description" = "gs_name"))

#rownames(logp_mat) <- logp_mat$gs_exact_source
#logp_mat = logp_mat %>% dplyr::select(names(ego_list))
colnames(logp_mat) <- names(ego_list)

# fill the matrix: signed -log10(p.adjust) for significant GO terms
for (i in seq_along(ego_list)) {
  ego_res <- ego_list[[i]]
  # filter significant GO terms
  ego_res <- ego_res[ego_res$p.adjust < padj_cutoff, ]
  
  # calculate signed log10(p.adjust)
  signed_logp <- -log10(ego_res$p.adjust) * sign(ego_res$NES)
  
  # assign GO ID as names
  names(signed_logp) <- ego_res$ID
  
  # fill in the matrix
  logp_mat[names(signed_logp), i] <- signed_logp
}
rownames(logp_mat) <- all_go_terms$gs_exact_source

# get GO term semantic similarity matrix
go_sim <- godata(annoDb = "org.Hs.eg.db", ont = "BP", computeIC = TRUE)
go_terms <- rownames(logp_mat)
# get IC 
ic_values <- go_sim@IC[go_terms]
# calculate semantic similarity
sim_matrix <- termSim(rownames(logp_mat), rownames(logp_mat), semData = go_sim, method = "Wang")
# transform similarity matrix to distance matrix
go_dist <- as.dist(1 - sim_matrix)


# Step 1: 构建聚类树 + 彩色树枝
col_clust <- hclust(go_dist, method = "ward.D2")
col_dend <- as.dendrogram(col_clust)

# -------------------------------
# Step 2: 构建 GO term 描述标签（替代默认 GO ID）
term_info <- ego_list[[1]]
term_desc_map <- setNames(term_info$Description, term_info$ID)

# 只保留你 heatmap 中用到的 GO term（列名）
term_desc_map <- term_desc_map[rownames(logp_mat)]

clusters <- cutree(col_clust, k = k_num)



cluster_colors <- pal_d3("category20")(k_num)
names(cluster_colors) <- as.character(1:k_num)

# 可选：预览颜色条
#barplot(rep(1, k_num), col = cluster_colors, names.arg = names(cluster_colors), las = 2, main = "Cluster Colors")

col_dend <- color_branches(col_dend, k = k_num,groupLabels = TRUE)  # 可调整聚类数 k

top_annot <- HeatmapAnnotation(
  Cluster = as.factor(clusters),
  col = list(Cluster = cluster_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "left"
)

# -------------------------------
# Step 4: 定义色条颜色函数
col_fun <- colorRamp2(
   c(floor(min(logp_mat, na.rm = TRUE)), 0, ceiling(max(logp_mat, na.rm = TRUE))),
  c("#1B5B9D", "white", "#D7191C")
)

# -------------------------------
# Step 5: 绘制热图
ht <- Heatmap(
  t(logp_mat),                                 # compare group in rows
  name = "-log10(p-value)",
  top_annotation = top_annot,
  cluster_columns = col_dend,                 # GO term 
  cluster_rows = FALSE,                        
  col = col_fun,
  border = TRUE,                              
  show_row_names = TRUE,
  show_column_names = FALSE,                   
  column_names_gp = gpar(fontsize = 6),      
  column_names_rot = 45,                     
  row_names_side = "left",
  column_dend_height = unit(2, "cm"),
  row_dend_width = unit(2, "cm"),
  column_title = "GO Term Clusters by Semantic Similarity",
  heatmap_legend_param = list(title = bquote("Signed " * -log[10]("adj. p-value")), at = seq(
  floor(min(logp_mat, na.rm = TRUE)),
  ceiling(max(logp_mat, na.rm = TRUE)),
  by = 2
))
)

ht_drawn <- draw(ht)

# grep the order of GO term names
# 提取列名（即聚类后顺序）
col_ordered <- column_order(ht_drawn)

# 获取聚类后列对应的 GO term ID（即原始矩阵中的列名）
clustered_GO_terms <- rownames(logp_mat)[col_ordered]

# get the final color in the plot
color_order = unique(clusters[clustered_GO_terms])

col_dend <- color_branches(col_dend, k = k_num,groupLabels = TRUE,col = cluster_colors[color_order])

ht <- ComplexHeatmap::Heatmap(
  t(logp_mat),                                 # compare group in rows
  name = "-log10(p-value)",
  top_annotation = top_annot,
  cluster_columns = col_dend,                 # GO term 
  cluster_rows = FALSE,                        
  col = col_fun,
  border = TRUE,                              
  show_row_names = TRUE,
  show_column_names = FALSE,                  
  column_names_gp = gpar(fontsize = 6),       
  column_names_rot = 45,                      
  row_names_side = "left",
  column_dend_height = unit(2, "cm"),
  row_dend_width = unit(2, "cm"),
  column_title = "GO Term Clusters by Semantic Similarity",
  heatmap_legend_param = list(title = bquote("Signed " * -log[10]("adj. p-value")), at =seq(
  floor(min(logp_mat, na.rm = TRUE)),
  ceiling(max(logp_mat, na.rm = TRUE)),
  by = 2
)
)
)

# save eps file
eps_path <- file.path(output_dir, "heatmap_plot.eps")
cairo_ps(
  filename = eps_path,
  width = max(ceiling(dim(logp_mat)[1] * 0.06), 6),
  height = max(ceiling(dim(logp_mat)[2] * 0.5), 4),
  family = "Arial",
  onefile = FALSE,
  fallback_resolution = 600
)
draw(ht)
dev.off()  

pdf_path <- file.path(output_dir, "heatmap_plot.pdf")
cairo_pdf(
  filename = pdf_path,
  width = max(ceiling(dim(logp_mat)[1] * 0.06), 6),
  height = max(ceiling(dim(logp_mat)[2] * 0.5), 4),
  family = "Arial",
  fallback_resolution = 600
)
draw(ht)
dev.off()

# save the heatmap result to file
df_long <- as.data.frame(logp_mat) %>%
  tibble::rownames_to_column("GO_term") %>%
  pivot_longer(-GO_term, names_to = "Comparison", values_to = "log10_pvalue")

# add cluster information
df_long$Cluster <- clusters[df_long$GO_term]

df_long_filtered <- df_long %>%
  #filter(log10_pvalue > 0) %>% 
  mutate(p_value = 10^(-log10_pvalue),
         IC = ic_values[GO_term])

#write.csv(df_long_filtered, file = "GO_term_enrichment_results.csv", row.names = FALSE)

# based on the min IC to represent the function of each cluster
go_representatives <- df_long_filtered %>%
  group_by(Cluster) %>%
  slice_min(IC, n = 1, with_ties = FALSE) %>%
  ungroup()

all_term_info <- lapply(ego_list, function(x) x[, c("ID", "Description")]) %>%
  bind_rows() %>%
  distinct(ID, .keep_all = TRUE) %>%
  left_join(all_go_terms, by = c("Description" = "Description")) %>%
  dplyr::select(gs_exact_source, Description) %>%
  rename(ID = gs_exact_source)

term_desc_full_map <- setNames(all_term_info$Description, all_term_info$ID)


df_long_filtered <- df_long_filtered %>%
  mutate(description = term_desc_full_map[GO_term])

go_representatives <- go_representatives %>%
  mutate(description = term_desc_full_map[GO_term])

# add Cluster_summury to df_long_filtered, using go_representatives$description
df_long_filtered =  df_long_filtered %>%
  left_join(
    go_representatives %>% dplyr::select(Cluster, Cluster_summary = description),
    by = "Cluster"
  )
df_long_filtered$description = gsub("GOBP_", "", df_long_filtered$description)
df_long_filtered$Cluster_summary = gsub("GOBP_", "", df_long_filtered$Cluster_summary)
# save the final data frame to a CSV file
output_csv_path <- file.path(output_dir, "GSEA_results_summary.csv")
write.csv(df_long_filtered, output_csv_path, row.names = FALSE)