#=========================Script Description=================================
# This script is used for GSVA analysis
# cd /Users/xliu2942/Documents/Projects/MAXOMOD/TESTing/Pipeline
# Rscript 19GSEA_allsamples.R >output.log
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("DEP"))
suppressMessages(library("clusterProfiler"))
suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library("msigdbr"))
suppressMessages(library("stringr"))
#===========================Function Definition=============================
#===========================Command Parameters Setting=========================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "8_Differential_expression_analysis_subclusters/res.rds"
  ),make_option(c("--output", "-o"),
                type = "character", default = "19_GSEA_subclusters_k2",
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
result = res$k2
# use the signed p-value as the gene list (Use origincal pvalue to avoid warning: 1: In preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam,  :There are ties in the preranked stats (33.27% of the list).)
# result$signed_pval <- -log10(result$fdr) * sign(result$als_vs_ctrl_diff)
result$signed_pval <- -log10(result$alpha_vs_beta_p.val) * sign(result$alpha_vs_beta_diff)
#geneList <- result$als_vs_ctrl_diff
geneList <- result$signed_pval
names(geneList) <- result$name
geneList <- sort(geneList, decreasing = TRUE)

GSEA_result = list()

ont = c("BP", "MF", "CC") 
for (i in ont) {
  
  #create background according to the ontology used for the analysis    
  bg <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = i) %>% 
    dplyr::select(gs_name, gene_symbol)
  bg <- bg[bg$gene_symbol %in% names(geneList), ]
  
  #the gsea analysis with cut-off   
  GSEA_res = clusterProfiler::GSEA(geneList=geneList, #performing the gsea
              nPermSimple = 100000, 
              minGSSize = 10, #minimum gene set size
              maxGSSize = 500, #maximum gene set size
              #pvalueCutoff = 0.05, 
              pvalueCutoff = 1, 
              verbose = TRUE, 
              TERM2GENE = bg, #background
              pAdjustMethod = "BH", #benjamini hochberg correction
              eps = 1e-10, 
              seed = 9)
  # simplify the GO terms
  #GSEA_res = clusterProfiler::simplify(x = GSEA_res, cutoff = 0.7)
  
  GSEA_result[[i]] = GSEA_res
}




if (nrow(GSEA_result$BP@result[GSEA_result$BP@result$p.adjust < 0.05,] ) >2){
GSEA_result_BP_005_cutoff = GSEA_result$BP@result[GSEA_result$BP@result$p.adjust < 0.05,] 
write.csv(GSEA_result_BP_005_cutoff, paste0(output_dir,"/GSEA_result_BP_005_cutoff.csv"))
write.csv(GSEA_result$BP@result, paste0(output_dir,"/GSEA_result_BP.csv"))
}
if (nrow(GSEA_result$MF@result[GSEA_result$MF@result$p.adjust < 0.05,] ) >2){
write.csv(GSEA_result$MF@result[GSEA_result$MF@result$p.adjust < 0.05,],paste0(output_dir,"/GSEA_result_MF_005_cutoff.csv"))
write.csv(GSEA_result$MF@result,paste0(output_dir,"/GSEA_result_MF.csv"))
}
if (nrow(GSEA_result$CC@result[GSEA_result$CC@result$p.adjust < 0.05,] ) >2){
write.csv(GSEA_result$CC@result[GSEA_result$CC@result$p.adjust < 0.05,],paste0(output_dir,"/GSEA_result_CC_005_cutoff.csv"))
write.csv(GSEA_result$CC@result,paste0(output_dir,"/GSEA_result_CC.csv"))
}

saveRDS(GSEA_result, paste0(output_dir,"/GSEA_result.rds"))


# plot
# select top 10 highest NES and lowest NES
if (dim(GSEA_result$BP@result)[1] > 40){
print("Top 20")
BP_top_highest <- GSEA_result$BP@result %>% 
  arrange(desc(NES)) %>%  # Arrange in descending order to get highest NES
  head(20)

BP_top_lowest <- GSEA_result$BP@result %>% 
  arrange(NES) %>%  # Arrange in ascending order to get lowest NES
  head(20)

top = 20
} else if (dim(GSEA_result$BP@result)[1] > 20){
print("Top 10")
BP_top_highest <- GSEA_result$BP@result %>% 
  arrange(desc(NES)) %>%  # Arrange in descending order to get highest NES
  head(10)

BP_top_lowest <- GSEA_result$BP@result %>%
  arrange(NES) %>%  # Arrange in ascending order to get lowest NES
  head(10)

top = 10
} else {
  print("Enriched terms less than 10")
  print("Skip the plot!")
}
  
if (top %in% c(10,20)){
# Combine the two data frames
BP_top <- bind_rows(BP_top_highest, BP_top_lowest)

BP_top <- BP_top %>%
  mutate(
    ID = str_remove(ID, "^GOBP_"),
    Description = str_remove(Description, "^GOBP_")
  )

library(ggplot2)
gsea_df <- BP_top %>%
  dplyr::select(Description, NES, "p.adjust") %>%
  dplyr::mutate(
    padj.group = cut(.$p.adjust, breaks = c(-Inf, 0.05, Inf), labels = c(1, 0)),
    nes.group = cut(.$NES, breaks = c(-Inf, 0, Inf), labels = c(0, 1)),
    comb.group = paste0(padj.group, nes.group)
  ) %>%
  dplyr::mutate(group = dplyr::case_when(
    comb.group == "10" ~ "CTRL", # CTRL is nes<0 & p<0.05
    comb.group == "11" ~ "ALS", # ALS is nes>0 & p<0.05
    TRUE ~ "C"
  )) %>% # C is grey
  dplyr::arrange(NES) %>%
  dplyr::mutate(index = 1:dplyr::n())

colour <- c("\\#4575B4", "\\#D73027", "\\#A9A9A9")
colour <- stringr::str_remove_all(colour, ".*#") %>% paste0("#", .)

# p <- ggplot(gsea_df, aes(x = index, y = NES, fill = group)) +
#   geom_bar(stat = "identity", width = 0.8) +
#   scale_fill_manual(values = c("CTRL" = colour[1], "ALS" = colour[2], "C" = colour[3])) +
#   scale_x_discrete(expand = expansion(add = .5)) +
#   scale_y_continuous(breaks = seq(
#     floor(min(gsea_df$NES)), ceiling(max(gsea_df$NES)),
#     ceiling((ceiling(max(gsea_df$NES)) - floor(min(gsea_df$NES))) / 6)
#   )) +
#   coord_flip()+
#   theme_classic() + 
#   theme(axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         plot.title = element_text(size = 15, hjust = 0.5)) +
#   geom_text(
#     data = subset(gsea_df, NES > 0),
#     aes(x = index, y = 0, label = paste0("  ", Description)),
#     size = 3, hjust = 1
#   ) +
#   geom_text(
#     data = subset(gsea_df, NES < 0),
#     aes(x = index, y = 0, label = paste0("  ", Description)),
#     size = 3, hjust = 0
#   )+labs(x = "", y = "Normalized enrichment score") 

# add column signed log10(padj)
gsea_df <- gsea_df %>%
  mutate(
    logp = -log10(p.adjust),
    signed_logp = ifelse(NES > 0, logp, -logp)
  )


# order based on signed_logp 
gsea_df <- gsea_df %>%
  arrange(signed_logp) %>%
  dplyr::mutate(Description = factor(Description, levels = Description))

# truncate the description, avoid too long description
max_char <- 50
gsea_df <- gsea_df %>%
  mutate(
    Description = as.character(Description),  # change into character, avoid invalid factor level, NA generated 
    Description_trunc = str_trunc(Description, width = max_char, side = "right")
  ) %>%
  arrange(signed_logp) %>%
  mutate(Description_trunc = factor(Description_trunc, levels = Description_trunc))

p <- ggplot(gsea_df, aes(x = signed_logp, y = Description_trunc, fill = NES > 0)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("TRUE" = "firebrick3", "FALSE" = "steelblue3"),
                    labels = c("Down in alpha", "Up in alpha"),
                    name = "") +
  labs(
    x = expression(-log[10](adjusted~p~value)),
    y = "",
    title = "GSEA"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right"
  )

ggsave(paste0(output_dir,"/GSEA_BP_top",top,".pdf"), p, width = 12, height = 10)
write.csv(BP_top, paste0(output_dir,"/GSEA_BP_top",top,".csv"),row.names = FALSE)
}