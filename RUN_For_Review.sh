# new msigdbr_26.1.0, old msigdbr 7.5.1
# packageVersion("msigdbr")

Rscript Script/12_Scatterplot_FDR.R -i 03_Differential_expression_analysis/Differential_Expression_Results.rds -o 12_Scatterplot_FDR -e 9

Rscript Script/12_Scatterplot_logFC.R -i 03_Differential_expression_analysis/Differential_Expression_Results.rds -o 12_Scatterplot_logFC -e 9

Rscript Script/21_RRHO2.R -i 03_Differential_expression_analysis/Differential_Expression_Results.rds -o 21_RRHO2 -e 9

Rscript Script/22_Scatter_intensity.R -i 02_Missing_Inspection/norm_imp_MinProb.rds -o 22_Scatter_intensity -e 9

Rscript Script/23_Heatmap_missingness.R -b 01_Pre_Processing/se_abu_data.rds -a 01_Pre_Processing/se_abu_data_filtered.rds -o 23_Heatmap_missingness -e 9

Rscript Script/24_clusterboot.R -i Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -o Discovery/24_clusterboot -e 9

Rscript Script/24_clusterboot.R -i Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -o Discovery/24_clusterboot -e 9 -v TRUE

Rscript Script/25_five_protein_panel_PCA.R --cohort Discovery
Rscript Script/25_five_protein_panel_PCA.R --cohort Validation

Rscript Script/26_vsCtr_proteins.R -o 26_vsCtr_proteins -e 9

Rscript Script/05_Vis_umap.R -i Discovery/02_Missing_Inspection \
  -o Discovery/05_Vis_umap -d all -e 9 -l FALSE -v batchid

Rscript Script/05_Vis_umap.R -i Discovery/02_Missing_Inspection \
  -o Discovery/05_Vis_umap -d all -e 9 -l TRUE -v batchid

Rscript Script/06_GSEA.R -i Discovery/03_Differential_expression_analysis/Differential_Expression_Results.rds -o Discovery/06_GSEA

Rscript Script/06_GSEA.R -i Validation/03_Differential_expression_analysis/Differential_Expression_Results.rds -o Validation/06_GSEA

Rscript Script/27_GSEA_IC_heatmap_replot.R -i Discovery=Discovery/06_GSEA/GSEA_result.rds,Validation=Validation/06_GSEA/GSEA_result.rds -o 27_GSEA_IC_heatmap_replot -c 6

Rscript Script/06_GSEA_subclusters.R -i Discovery/03_Differential_expression_analysis_subclusters/res.rds -o Discovery/06_GSEA_subclusters_k2

Rscript Script/06_GSEA_subclusters.R -i Validation/03_Differential_expression_analysis_subclusters/res.rds -o Validation/06_GSEA_subclusters_k2

Rscript Script/27_GSEA_IC_heatmap_replot.R -i Discovery=Discovery/06_GSEA_subclusters_k2/GSEA_result.rds,Validation=Validation/06_GSEA_subclusters_k2/GSEA_result.rds,External=External/26_External_GSEA_subclusters/GSEA_result.rds -o 27_GSEA_IC_heatmap_subclusters_all_cohorts_replot -c 6

Rscript Script/28_Protein_numbers.R -o 28_Protein_numbers -e 9

Rscript Script/29_missing_per_protein.R -o 29_missing_per_protein -e 9

Rscript Script/30_GSVA.R -i Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -c Discovery/08_Clustering_als/cluster_assignments_2.csv -e 9 -o Discovery/30_GSVA -u Discovery/01_Pre_Processing_als/uniprot_to_genename.rds -s GO_BP --number 5 --width 10

Rscript Script/30_GSVA.R -i Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -c Validation/08_Clustering_als/cluster_assignments_2.csv -e 9 -o Validation/30_GSVA -u Validation/01_Pre_Processing_als/uniprot_to_genename.rds -s GO_BP --number 3 --width 7

Rscript Script/31_NEFL.R -o 31_NEFL -e 9

Rscript Script/32_subcluster_vsCTR.R \
  -i Discovery/02_Missing_Inspection/norm_imp_MinProb.rds \
  -c Discovery/08_Clustering_als/cluster_assignments_2.csv \
  -u Discovery/01_Pre_Processing/uniprot_to_genename.rds \
  -o Discovery/32_subcluster_vsCTR -e 9

Rscript Script/32_subcluster_vsCTR.R \
  -i Validation/02_Missing_Inspection/norm_imp_MinProb.rds \
  -c Validation/08_Clustering_als/cluster_assignments_2.csv \
  -u Validation/01_Pre_Processing/uniprot_to_genename.rds \
  -o Validation/32_subcluster_vsCTR -e 9


Rscript Script/33_subcluster_vsCTR_vis.R -i Discovery/32_subcluster_vsCTR/Differential_Expression_Results.rds -o Discovery/33_subcluster_vsCTR_vis -e 9

Rscript Script/33_subcluster_vsCTR_vis.R -i Validation/32_subcluster_vsCTR/Differential_Expression_Results.rds -o Validation/33_subcluster_vsCTR_vis -e 9

Rscript Script/34_subcluster_vsCTR_GSEA.R -i Discovery/32_subcluster_vsCTR/Differential_Expression_Results.rds -o Discovery/34_subcluster_vsCTR_GSEA -e 9

Rscript Script/34_subcluster_vsCTR_GSEA.R -i Validation/32_subcluster_vsCTR/Differential_Expression_Results.rds -o Validation/34_subcluster_vsCTR_GSEA -e 9

Rscript Script/35_subcluster_vsCTR_GSEA_vis.R -D Discovery/34_subcluster_vsCTR_GSEA -V Validation/34_subcluster_vsCTR_GSEA -o 35_subcluster_vsCTR_GSEA_vis -e 9

Rscript Script/36_MA_plot.R -f 0.05 -o 36_MA_plot -e 9
Rscript Script/36_MA_plot.R -s p.val -f 0.05 -o 36_MA_plot_pval -e 9

Rscript Script/37_Scatterplot_subtype_Ctrl.R -i 32_subcluster_vsCTR/Differential_Expression_Results.rds -o 37_Scatterplot_subtype_Ctrl

# # remove sample from WGCNA
# Rscript Script/38_WGCNA_subclusters_sample_removed.R -i Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -o Discovery/38_WGCNA_subclusters_sample_removed -e 9 -c Discovery/08_Clustering_als/cluster_assignments_2.csv

# opt$input = "Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds"
# opt$output = "Discovery/38_WGCNA_subclusters_sample_removed"
# opt$seed = 9
# opt$cluster_assignments = "Discovery/08_Clustering_als/cluster_assignments_2.csv"

# Rscript Script/39_WGCNA_comparison_subclusters_sample_removed.R -i Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -o Discovery/39_WGCNA_comparison_subclusters_sample_removed --net Discovery/38_WGCNA_subclusters_sample_removed/WGCNA_net.rds -e 9 -c Discovery/08_Clustering_als/cluster_assignments_2.csv --ModuleTrait TRUE --cohort Discovery

# Rscript Script/39_WGCNA_comparison_subclusters_sample_removed.R -i Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -o Validation/39_WGCNA_comparison_subclusters_sample_removed --net Discovery/38_WGCNA_subclusters_sample_removed/WGCNA_net.rds -e 9 -c Validation/08_Clustering_als/cluster_assignments_2.csv --ModuleTrait TRUE --cohort Validation

# # opt$input = "Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds"
# # opt$output = "Validation/39_WGCNA_comparison_subclusters_sample_removed"
# # opt$net = "Discovery/38_WGCNA_subclusters_sample_removed/WGCNA_net.rds"
# # opt$seed = 9
# # opt$cluster_assignments = "Validation/08_Clustering_als/cluster_assignments_2.csv"
# # opt$ModuleTrait = TRUE
# # opt$cohort = "Validation"

# Rscript Script/11_WGCNA_comparison_subclusters.R -i External/norm_imp_MinProb_als.rds -o External/39_WGCNA_comparison_subclusters_sample_removed --net Discovery/38_WGCNA_subclusters_sample_removed/WGCNA_net.rds -e 9 -c External/cluster_assignments_2.csv

# Rscript Script/40_WGCNA_Enrich.R -w Discovery/38_WGCNA_subclusters_sample_removed -i Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -c Discovery/08_Clustering_als/cluster_assignments_2.csv -o Discovery/40_WGCNA_Enrich -t 3 --cell_width 0.6 --cell_height 0.8 --textbox_k 5 --textbox_width 4 --top_k_proteins 5 --textbox_k 3

# # For original 
# Rscript Script/40_WGCNA_Enrich.R -w Discovery/10_WGCNA_subclusters -i Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -c Discovery/08_Clustering_als/cluster_assignments_2.csv -o Discovery/40_WGCNA_Enrich_original -t 3 --cell_width 0.6 --cell_height 0.8 --textbox_k 5 --textbox_width 4 --top_k_proteins 5 --textbox_k 3 --row_label_wrap 32 --sample_removed "" 

# # opt$wgcna_dir = "Discovery/38_WGCNA_subclusters_sample_removed"
# # opt$input = "Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds"
# # opt$cluster_assignments = "Discovery/08_Clustering_als/cluster_assignments_2.csv"
# # opt$output = "Discovery/40_WGCNA_Enrich"
# # opt$top_n = 3

# Rscript Script/41_WGCNA_Enrich_ALL.R -w Discovery/38_WGCNA_subclusters_sample_removed -i Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -c Discovery/08_Clustering_als/cluster_assignments_2.csv -o Discovery/41_WGCNA_Enrich_ALL -t 3

# opt$wgcna_dir = "Discovery/38_WGCNA_subclusters_sample_removed"
# opt$input = "Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds"
# opt$cluster_assignments = "Discovery/08_Clustering_als/cluster_assignments_2.csv"
# opt$output = "Discovery/41_WGCNA_Enrich_ALL"
# opt$top_n = 3

# Trait-only re-run of Module-Trait relationships + k2 boxplot on top of the
# ORIGINAL (non-sample-removed) 10_WGCNA_subclusters network, after removing
# the same samples used in 38_WGCNA_subclusters_sample_removed.R. Does NOT
# re-run WGCNA itself.
Rscript Script/38_WGCNA_subclusters_sample_removed_traitonly.R -i Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -n Discovery/10_WGCNA_subclusters/WGCNA_net.rds -c Discovery/08_Clustering_als/cluster_assignments_2.csv -o Discovery/38_WGCNA_subclusters_sample_removed_traitonly --cohort Discovery -e 9

Rscript Script/38_WGCNA_subclusters_sample_removed_traitonly.R -i Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -n Discovery/10_WGCNA_subclusters/WGCNA_net.rds -c Validation/08_Clustering_als/cluster_assignments_2.csv -o Validation/38_WGCNA_subclusters_sample_removed_traitonly --cohort Validation -e 9

Rscript Script/38_WGCNA_subclusters_sample_removed_traitonly.R -i External/norm_imp_MinProb_als.rds -o External/39_WGCNA_comparison_subclusters_sample_removed --net Discovery/10_WGCNA_subclusters/WGCNA_net.rds -e 9 -c External/cluster_assignments_2.csv --cohort External

Rscript Script/42_Convert_Protein.R

Rscript Script/43_Predict_ELA_LS_Subtype.R

Rscript Script/44_Pre_Processing_ELA_LS.R -i 42_Convert_Protein/ELA-LS_Report_with_gene_names.csv -o 44_Pre_Processing_ELA_LS -e 9 -t 0.5 -q 0.01