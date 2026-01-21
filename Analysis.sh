# # nohup bash Analysis.sh > Analysis.log 2>&1 &
# ########################################################### Data preparation (Discovery cohort) ###########################################################
# # Discovery cohort data preparation
# echo "Discovery cohort data preparation"
# Rscript 00_Discovery_Data_Prepare.R

# ########################################################### Data analysis (Discovery cohort, All samples) ###########################################################
# # All samples
# Rscript 01_Pre_Processing.R -i Discovery/data/proteinGroups.txt -c Discovery/data/clinical_info.csv -o Discovery/01_Pre_Processing -e 9 -d Input/HUMAN_9606_idmapping_selected.tab.gz

# Rscript 02_Missing_Inspection.R -i Discovery/01_Pre_Processing/se_abu_data_filtered.rds -o Discovery/02_Missing_Inspection -e 9 -s 0.01

# Rscript 03_Differential_expression_analysis.R -i Discovery/02_Missing_Inspection/norm_imp_MinProb.rds -o Discovery/03_Differential_expression_analysis -e 9 -u Discovery/01_Pre_Processing/uniprot_to_genename.rds

# Rscript 04_Vis_Differential_expression_analysis.R -i Discovery/03_Differential_expression_analysis/Differential_Expression_Results.rds -o Discovery/04_Vis_Differential_expression_analysis -e 9
# Rscript 04_Vis_Differential_expression_analysis_FDR_0102.R -i Discovery/03_Differential_expression_analysis/Differential_Expression_Results.rds -o Discovery/04_Vis_Differential_expression_analysis_FDR_0102 -e 9

# echo "Rscript 05_Vis_umap.R"
# Rscript 05_Vis_umap.R -i Discovery/02_Missing_Inspection -o Discovery/05_Vis_umap -e 9 -d all -l TRUE
# Rscript 05_Vis_umap.R -i Discovery/02_Missing_Inspection -o Discovery/05_Vis_umap -e 9 -d all -l FALSE

# echo "Rscript Rscript 06_GSEA.R"
# Rscript 06_GSEA.R -i Discovery/03_Differential_expression_analysis/Differential_Expression_Results.rds -o Discovery/06_GSEA

# Rscript 07_Clinical_heatmap.R --input Discovery/02_Missing_Inspection/norm_imp_MinProb.rds --output Discovery/07_Clinical_heatmap --seed 9

# ########################################################### Data analysis (Discovery cohort, ALS only) ###########################################################
# # als only
# Rscript 01_Pre_Processing.R -i Discovery/data/proteinGroups.txt -c Discovery/data/clinical_info.csv -o Discovery/01_Pre_Processing_als -e 9 --subset als -d Input/HUMAN_9606_idmapping_selected.tab.gz

# Rscript 02_Missing_Inspection.R  -i Discovery/01_Pre_Processing_als/se_abu_data_filtered.rds -o Discovery/02_Missing_Inspection_subclusters -e 9 -s 0.01

# Rscript 08_Clustering_subclusters.R -i Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -o Discovery/08_Clustering_als -e 9

# Rscript 03_Differential_expression_analysis_subclusters.R -i Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -o Discovery/03_Differential_expression_analysis_subclusters -e 9 -c Discovery/08_Clustering_als/cluster_assignments_2.csv

# Rscript 04_Vis_Differential_expression_analysis_subclusters.R -i Discovery/03_Differential_expression_analysis_subclusters/res.rds -o Discovery/04_Vis_Differential_expression_analysis_subclusters -d Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -c Discovery/08_Clustering_als/cluster_assignments_2.csv -e 9

# echo "Rscript 05_Vis_umap.R"
# Rscript 05_Vis_umap.R -i Discovery/02_Missing_Inspection_subclusters -o Discovery/05_Vis_umap -d als -e 9 -l TRUE -c Discovery/08_Clustering_als

# Rscript 05_Vis_umap.R -i Discovery/02_Missing_Inspection_subclusters -o Discovery/05_Vis_umap -d als -e 9 -l FALSE -c Discovery/08_Clustering_als

# Rscript 06_GSEA_subclusters.R -i Discovery/03_Differential_expression_analysis_subclusters/res.rds -o Discovery/06_GSEA_subclusters_k2

# # sex
# Rscript 09_Clinical_features_subclusters.R --input Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds --cluster Discovery/08_Clustering_als/cluster_assignments_2.csv --output Discovery/09_Clinical_features_subclusters --var sex

# # progression group
# Rscript 09_Clinical_features_subclusters.R --input Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds --cluster Discovery/08_Clustering_als/cluster_assignments_2.csv --output Discovery/09_Clinical_features_subclusters --var progression_group

# #onset
# Rscript 09_Clinical_features_subclusters.R --input Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds --cluster Discovery/08_Clustering_als/cluster_assignments_2.csv --output Discovery/09_Clinical_features_subclusters --var condition

# #age
# Rscript 09_Clinical_features_subclusters.R --input Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds --cluster Discovery/08_Clustering_als/cluster_assignments_2.csv --output Discovery/09_Clinical_features_subclusters --var age

# #age at onset
# Rscript 09_Clinical_features_subclusters.R --input Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds --cluster Discovery/08_Clustering_als/cluster_assignments_2.csv --output Discovery/09_Clinical_features_subclusters --var age_at_onset

# #Nfl
# Rscript 09_Clinical_features_subclusters.R --input Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds --cluster Discovery/08_Clustering_als/cluster_assignments_2.csv --output Discovery/09_Clinical_features_subclusters --var Nfl

# #pNFh
# Rscript 09_Clinical_features_subclusters.R --input Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds --cluster Discovery/08_Clustering_als/cluster_assignments_2.csv --output Discovery/09_Clinical_features_subclusters --var pNFh

# echo "Rscript 10_WGCNA_subclusters.R"
# Rscript 10_WGCNA_subclusters.R -i Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -o Discovery/10_WGCNA_subclusters -e 9 -c Discovery/08_Clustering_als/cluster_assignments_2.csv

# ########################################################### Data preparation (Validation cohort) ###########################################################
# # Validation cohort data preparation
# echo "Validation cohort data preparation"
# Rscript 00_Validation_Data_Prepare.R

# ########################################################### Data analysis (Validation cohort, All samples) ###########################################################
# # All samples
# Rscript 01_Pre_Processing.R -i Validation/data/proteinGroups.txt -c Validation/data/clinical_info.csv -o Validation/01_Pre_Processing -e 9 -d Input/HUMAN_9606_idmapping_selected.tab.gz
# echo "Rscrtipt 02_Missing_Inspection.R "

# Rscript 02_Missing_Inspection.R  -i Validation/01_Pre_Processing/se_abu_data_filtered.rds -o Validation/02_Missing_Inspection -e 9 -s 0.01

# Rscript 03_Differential_expression_analysis.R -i Validation/02_Missing_Inspection/norm_imp_MinProb.rds -o Validation/03_Differential_expression_analysis -e 9 -u Validation/01_Pre_Processing/uniprot_to_genename.rds

# Rscript 04_Vis_Differential_expression_analysis.R -i Validation/03_Differential_expression_analysis/Differential_Expression_Results.rds -o Validation/04_Vis_Differential_expression_analysis -e 9
# Rscript 04_Vis_Differential_expression_analysis_FDR_0102.R -i Validation/03_Differential_expression_analysis/Differential_Expression_Results.rds -o Validation/04_Vis_Differential_expression_analysis_FDR_0102 -e 9

# Rscript 05_Vis_umap.R -i Validation/02_Missing_Inspection -o Validation/05_Vis_umap -e 9 -d all -l TRUE
# Rscript 05_Vis_umap.R -i Validation/02_Missing_Inspection -o Validation/05_Vis_umap -e 9 -d all -l FALSE

# Rscript 06_GSEA.R -i Validation/03_Differential_expression_analysis/Differential_Expression_Results.rds -o Validation/06_GSEA

# Rscript 07_Clinical_heatmap.R --input Validation/02_Missing_Inspection/norm_imp_MinProb.rds --output Validation/07_Clinical_heatmap --seed 9

# ########################################################### Data analysis (Validation cohort, ALS only) ###########################################################
# # als only
# Rscript 01_Pre_Processing.R -i Validation/data/proteinGroups.txt -c Validation/data/clinical_info.csv -o Validation/01_Pre_Processing_als -e 9 --subset als -d Input/HUMAN_9606_idmapping_selected.tab.gz

# Rscript 02_Missing_Inspection.R  -i Validation/01_Pre_Processing_als/se_abu_data_filtered.rds -o Validation/02_Missing_Inspection_subclusters -e 9 -s 0.01

# Rscript 08_Clustering_subclusters.R -i Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -o Validation/08_Clustering_als -e 9 --reverse TRUE

# Rscript 03_Differential_expression_analysis_subclusters.R -i Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -o Validation/03_Differential_expression_analysis_subclusters -e 9 -c Validation/08_Clustering_als/cluster_assignments_2.csv

# Rscript 04_Vis_Differential_expression_analysis_subclusters.R -i Validation/03_Differential_expression_analysis_subclusters/res.rds -o Validation/04_Vis_Differential_expression_analysis_subclusters -d Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -c Validation/08_Clustering_als/cluster_assignments_2.csv -e 9

# Rscript 05_Vis_umap.R -i Validation/02_Missing_Inspection_subclusters -o Validation/05_Vis_umap -d als -e 9 -l TRUE -c Validation/08_Clustering_als
# Rscript 05_Vis_umap.R -i Validation/02_Missing_Inspection_subclusters -o Validation/05_Vis_umap -d als -e 9 -l FALSE -c Validation/08_Clustering_als

# Rscript 06_GSEA_subclusters.R -i Validation/03_Differential_expression_analysis_subclusters/res.rds -o Validation/06_GSEA_subclusters_k2

# # sex
# Rscript 09_Clinical_features_subclusters.R --input Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds --cluster Validation/08_Clustering_als/cluster_assignments_2.csv --output Validation/09_Clinical_features_subclusters --var sex

# # progression group
# Rscript 09_Clinical_features_subclusters.R --input Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds --cluster Validation/08_Clustering_als/cluster_assignments_2.csv --output Validation/09_Clinical_features_subclusters --var progression_group

# #onset
# Rscript 09_Clinical_features_subclusters.R --input Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds --cluster Validation/08_Clustering_als/cluster_assignments_2.csv --output Validation/09_Clinical_features_subclusters --var condition

# #age
# Rscript 09_Clinical_features_subclusters.R --input Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds --cluster Validation/08_Clustering_als/cluster_assignments_2.csv --output Validation/09_Clinical_features_subclusters --var age

# #age at onset
# Rscript 09_Clinical_features_subclusters.R --input Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds --cluster Validation/08_Clustering_als/cluster_assignments_2.csv --output Validation/09_Clinical_features_subclusters --var age_at_onset

# #Nfl
# Rscript 09_Clinical_features_subclusters.R --input Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds --cluster Validation/08_Clustering_als/cluster_assignments_2.csv --output Validation/09_Clinical_features_subclusters --var Nfl

# # WGCNA using module detected in Discovery, to see it's performance in Discovery
# Rscript 11_WGCNA_comparison_subclusters.R -i Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -o Discovery/11_WGCNA_comparison --net Discovery/10_WGCNA_subclusters/WGCNA_net.rds -e 9 -c Discovery/08_Clustering_als/cluster_assignments_2.csv --ModuleTrait TRUE    

# # WGCNA using module detected in Discovery, to see it's performance in Validation
# Rscript 11_WGCNA_comparison_subclusters.R -i Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -o Validation/11_WGCNA_comparison --net Discovery/10_WGCNA_subclusters/WGCNA_net.rds -e 9 -c Validation/08_Clustering_als/cluster_assignments_2.csv --ModuleTrait TRUE

# ########################################################### Between Validation and Discovery cohort ###########################################################
# #scatterplot
# Rscript 12_Scatterplot_FDR.R -i 03_Differential_expression_analysis/Differential_Expression_Results.rds -o 12_Scatterplot_FDR -e 9

# Rscript 12_Scatterplot_FDR_subclusters.R -i 03_Differential_expression_analysis_subclusters/res.rds -o 12_Scatterplot_FDR_subclusters

# Rscript 12_Scatterplot_FDR_0102.R -i 03_Differential_expression_analysis/Differential_Expression_Results.rds -o 12_Scatterplot_FDR_01 -e 9 --FDR_Cutoff 0.1

# Rscript 12_Scatterplot_FDR_0102.R -i 03_Differential_expression_analysis/Differential_Expression_Results.rds -o 12_Scatterplot_FDR_02 -e 9 --FDR_Cutoff 0.2

# # GSEA heatmap
# Rscript 13_GSEA_IC_heatmap.R -i Discovery=Discovery/06_GSEA/GSEA_result.rds,Validation=Validation/06_GSEA/GSEA_result.rds -o 13_GSEA_IC_heatmap -c 6

# Rscript 13_GSEA_IC_heatmap.R -i Discovery=Discovery/06_GSEA_subclusters_k2/GSEA_result.rds,Validation=Validation/06_GSEA_subclusters_k2/GSEA_result.rds -o 13_GSEA_IC_heatmap_subclusters -c 6

# Rscript 14_Pathway_between_cohort_GSEA.R
# Rscript 14_Pathway_between_cohort_GSEA_subclusters.R

# #input file for Discovery should be "Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds"
# #input file for Validation should be "Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds"
# Rscript 15_ML_multi_models.R --n_boot 500 --feature_freq_cutoff 0.1 --sd_threshold 3 --output 15_ML_multi_models_500_01_noAge_3sd
# Rscript 15_ML_multi_models.R --n_boot 500 --feature_freq_cutoff 0.2 --sd_threshold 3 --output 15_ML_multi_models_500_02_noAge_3sd
# Rscript 15_ML_multi_models.R --n_boot 500 --feature_freq_cutoff 0.3 --sd_threshold 3 --output 15_ML_multi_models_500_03_noAge_3sd
# Rscript 15_ML_multi_models.R --n_boot 500 --feature_freq_cutoff 0.4 --sd_threshold 3 --output 15_ML_multi_models_500_04_noAge_3sd
# Rscript 15_ML_multi_models.R --n_boot 500 --feature_freq_cutoff 0.5 --sd_threshold 3 --output 15_ML_multi_models_500_05_noAge_3sd
# Rscript 15_ML_multi_models.R --n_boot 500 --feature_freq_cutoff 0.6 --sd_threshold 3 --output 15_ML_multi_models_500_06_noAge_3sd

# ML Step by Step
Rscript 16_ML_multi_models_Step_by_step.R --n_boot 500 --topn 10 --sd_threshold 3 --output 15_ML_multi_models_Step_by_step

# AUC line plot
Rscript 17_AUC_comparision.R --output_dir 17_AUC_comparision

# # For Discovery
# Rscript 40Correlation_LASSO_protein_vs_Clinical.R --input Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -s 15_ML_multi_models/selected_features.txt -c Discovery/08_Clustering_als/cluster_assignments_2.csv -g Discovery/03_Differential_expression_analysis_subclusters/Differential_expression_analysis_for_k2.csv -o 40_Correlation_LASSO_protein_vs_Clinical_Discovery 
# # For Validation
# Rscript 40Correlation_LASSO_protein_vs_Clinical.R --input Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -s 15_ML_multi_models/selected_features.txt -c Validation/08_Clustering_als/cluster_assignments_2.csv -g Validation/03_Differential_expression_analysis_subclusters/Differential_expression_analysis_for_k2.csv -o 40_Correlation_LASSO_protein_vs_Clinical_Validation 

########################################################### External Cohort ###########################################################
Rscript 18_External.R

Rscript 18_GSEA_External.R -i External/als_vs_ctrl.csv -o External/18_GSEA_External
Rscript 18_GSEA_External_subclusters.R -i External/alpha_vs_beta_als.csv -o External/18_GSEA_External_subclusters_k2

Rscript 18_External_protein_expression.R -i External/norm_imp_MinProb_als.rds --features 15_ML_multi_models_500_04_noAge_3sd/selected_features.txt -o External/18_External_protein_expression/

# WGCNA using module detected in Discovery, to see it's performance in External Cohort
Rscript 11_WGCNA_comparison_subclusters.R -i External/norm_imp_MinProb_als.rds -o External/11_WGCNA_comparison --net Discovery/10_WGCNA_subclusters/WGCNA_net.rds -e 9 -c External/cluster_assignments_2.csv

Rscript 13_GSEA_IC_heatmap.R -i Discovery=Discovery/06_GSEA_subclusters_k2/GSEA_result.rds,Validation=Validation/06_GSEA_subclusters_k2/GSEA_result.rds,External=External/18_GSEA_External_subclusters_k2/GSEA_result.rds -o 13_GSEA_IC_heatmap_subclusters_all_cohorts -c 6

########################################################### Brain data ###########################################################
Rscript 19_Brain.R

# All samples
Rscript 03_Differential_expression_analysis.R -i Brain/02_Missing_Inspection/norm_imp_MinProb.rds -o Brain/03_Differential_expression_analysis -e 9 -u Brain/01_Pre_Processing/uniprot_to_genename.rds

Rscript 04_Vis_Differential_expression_analysis.R -i Brain/03_Differential_expression_analysis/Differential_Expression_Results.rds -o Brain/04_Vis_Differential_expression_analysis -e 9

Rscript 05_Vis_umap.R -i Brain/02_Missing_Inspection -o Brain/05_Vis_umap -e 9 -d all -l TRUE

Rscript 05_Vis_umap.R -i Brain/02_Missing_Inspection -o Brain/05_Vis_umap -e 9 -d all -l FALSE

Rscript 06_GSEA.R -i Brain/03_Differential_expression_analysis/Differential_Expression_Results.rds -o Brain/06_GSEA

# als only
# this script is used to cluster the als samples, without NFL information
Rscript 08_Clustering_subclusters_brain.R -i Brain/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -o Brain/08_Clustering_als -e 9

Rscript 03_Differential_expression_analysis_subclusters.R -i Brain/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -o Brain/03_Differential_expression_analysis_subclusters -e 9 -c Brain/08_Clustering_als/cluster_assignments_2.csv

# when k =3, there is no alpha or theta subcluster specific proteins, therefore, cannot draw the expression heatmap. This will return an error.
Rscript 04_Vis_Differential_expression_analysis_subclusters.R -i Brain/03_Differential_expression_analysis_subclusters/res.rds -o Brain/04_Vis_Differential_expression_analysis_subclusters -d Brain/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -c Brain/08_Clustering_als/cluster_assignments_2.csv -e 9

Rscript 05_Vis_umap.R -i Brain/02_Missing_Inspection_subclusters -o Brain/05_Vis_umap -d als -e 9 -l TRUE -c Brain/08_Clustering_als
Rscript 05_Vis_umap.R -i Brain/02_Missing_Inspection_subclusters -o Brain/05_Vis_umap -d als -e 9 -l FALSE -c Brain/08_Clustering_als

Rscript 06_GSEA_subclusters.R -i Brain/03_Differential_expression_analysis_subclusters/res.rds -o Brain/06_GSEA_subclusters_k2

Rscript 13_GSEA_IC_heatmap.R -i Brain=Brain/06_GSEA_subclusters_k2/GSEA_result.rds -o Brain/13_GSEA_IC_heatmap_subclusters -c 6