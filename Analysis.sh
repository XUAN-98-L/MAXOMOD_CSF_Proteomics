########################################################### Data preparation (Discovery cohort) ###########################################################
# Discovery cohort data preparation
echo "Discovery cohort data preparation"
Rscript Discovery_Data_prepare_updated_20250731.R

########################################################### Data analysis (Discovery cohort, All samples) ###########################################################
# All samples
Rscript 1Pre_processing_updated_20250731.R -i Discovery/data/proteinGroups.txt -c Discovery/data/clinical_info.csv -o Discovery/1_pre_processing -e 9 -d Input/HUMAN_9606_idmapping_selected.tab.gz

Rscript 2Missing_Inspection.R -i Discovery/1_pre_processing/se_abu_data_filtered.rds -o Discovery/2_Missing_Inspection -e 9 -s 0.01

Rscript 3Differential_expression_analysis.R -i Discovery/2_Missing_Inspection/norm_imp_MinProb.rds -o Discovery/3_Differential_Expression_Analysis -e 9 -u Discovery/1_pre_processing/uniprot_to_genename.rds

Rscript 12Vis_Differential_expression_analysis.R -i Discovery/3_Differential_Expression_Analysis/Differential_Expression_Results.rds -o Discovery/12_Vis_Differential_expression_analysis -e 9
Rscript 12Vis_Differential_expression_analysis_FDR_0102.R -i Discovery/3_Differential_Expression_Analysis/Differential_Expression_Results.rds -o Discovery/12_Vis_Differential_expression_analysis_FDR_0102 -e 9

echo "Rscript 13Vis_umap.R"
Rscript 13Vis_umap.R -i Discovery/2_Missing_Inspection -o Discovery/13_Vis_umap -e 9 -d all -l TRUE
Rscript 13Vis_umap.R -i Discovery/2_Missing_Inspection -o Discovery/13_Vis_umap -e 9 -d all -l FALSE

echo "Rscript Rscript 19GSEA_allsamples.R"
Rscript 19GSEA_allsamples.R -i Discovery/3_Differential_Expression_Analysis/Differential_Expression_Results.rds -o Discovery/19_GSEA_allsamples

Rscript 37_Clinical_heatmap_highest_variance.R --input Discovery/2_Missing_Inspection/norm_imp_MinProb.rds --output Discovery/37_Clinical_heatmap_top_variable --seed 9

########################################################### Data analysis (Discovery cohort, ALS only) ###########################################################
# als only
Rscript 1Pre_processing_updated_20250731.R -i Discovery/data/proteinGroups.txt -c Discovery/data/clinical_info.csv -o Discovery/1_pre_processing_als -e 9 --subset als -d Input/HUMAN_9606_idmapping_selected.tab.gz

Rscript 2Missing_Inspection.R -i Discovery/1_pre_processing_als/se_abu_data_filtered.rds -o Discovery/2_Missing_Inspection_als -e 9 -s 0.01

Rscript 5Clustering_updated_20250506.R -i Discovery/2_Missing_Inspection_als/norm_imp_MinProb.rds -o Discovery/5_Clustering_als -e 9

Rscript 8Differential_expression_analysis_subclusters.R -i Discovery/2_Missing_Inspection_als/norm_imp_MinProb.rds -o Discovery/8_Differential_expression_analysis_subclusters -e 9 -c Discovery/5_Clustering_als/cluster_assignments_2.csv

Rscript 9Vis_Differential_expression_analysis_subclusters.R -i Discovery/8_Differential_expression_analysis_subclusters/res.rds -o Discovery/9_Vis_Differential_expression_analysis_subclusters -d Discovery/2_Missing_Inspection_als/norm_imp_MinProb.rds -c Discovery/5_Clustering_als/cluster_assignments_2.csv -e 9

echo "Rscript 13Vis_umap.R"
Rscript 13Vis_umap.R -i Discovery/2_Missing_Inspection_als -o Discovery/13_Vis_umap -d als -e 9 -l TRUE -c Discovery/5_Clustering_als

Rscript 13Vis_umap.R -i Discovery/2_Missing_Inspection_als -o Discovery/13_Vis_umap -d als -e 9 -l FALSE -c Discovery/5_Clustering_als

echo "Rscript 19GSEA_subclusters.R"
Rscript 19GSEA_subclusters.R -i Discovery/8_Differential_expression_analysis_subclusters/res.rds -o Discovery/19_GSEA_subclusters_k2

echo "Rscript 23WGCNA.R"
Rscript 23WGCNA.R -i Discovery/2_Missing_Inspection_als/norm_imp_MinProb.rds -o Discovery/23_WGCNA -e 9 -c Discovery/5_Clustering_als/cluster_assignments_2.csv

# sex
Rscript 38Clinical_features.R --input Discovery/2_Missing_Inspection_als/norm_imp_MinProb.rds --cluster Discovery/5_Clustering_als/cluster_assignments_2.csv --output Discovery/38_Clinical_features --var sex

# progression group
Rscript 38Clinical_features.R --input Discovery/2_Missing_Inspection_als/norm_imp_MinProb.rds --cluster Discovery/5_Clustering_als/cluster_assignments_2.csv --output Discovery/38_Clinical_features --var progression_group

#onset
Rscript 38Clinical_features.R --input Discovery/2_Missing_Inspection_als/norm_imp_MinProb.rds --cluster Discovery/5_Clustering_als/cluster_assignments_2.csv --output Discovery/38_Clinical_features --var condition

#age
Rscript 38Clinical_features.R --input Discovery/2_Missing_Inspection_als/norm_imp_MinProb.rds --cluster Discovery/5_Clustering_als/cluster_assignments_2.csv --output Discovery/38_Clinical_features --var age

#age at onset
Rscript 38Clinical_features.R --input Discovery/2_Missing_Inspection_als/norm_imp_MinProb.rds --cluster Discovery/5_Clustering_als/cluster_assignments_2.csv --output Discovery/38_Clinical_features --var age_at_onset

#Nfl
Rscript 38Clinical_features.R --input Discovery/2_Missing_Inspection_als/norm_imp_MinProb.rds --cluster Discovery/5_Clustering_als/cluster_assignments_2.csv --output Discovery/38_Clinical_features --var Nfl

#pNFh
Rscript 38Clinical_features.R --input Discovery/2_Missing_Inspection_als/norm_imp_MinProb.rds --cluster Discovery/5_Clustering_als/cluster_assignments_2.csv --output Discovery/38_Clinical_features --var pNFh

########################################################### Data preparation (Validation cohort) ###########################################################
# Validation cohort data preparation
echo "VAlidation cohort data preparation"
Rscript Validation_Data_prepare_updated_20250731.R

########################################################### Data analysis (Validation cohort, All samples) ###########################################################
# All samples
Rscript 1Pre_processing_updated_20250731.R -i Validation/data/proteinGroups.txt -c Validation/data/clinical_info.csv -o Validation/1_pre_processing -e 9 -d Input/HUMAN_9606_idmapping_selected.tab.gz
echo "Rscrtipt 2Missing_Inspection.R"

Rscript 2Missing_Inspection.R -i Validation/1_pre_processing/se_abu_data_filtered.rds -o Validation/2_Missing_Inspection -e 9 -s 0.01

Rscript 3Differential_expression_analysis.R -i Validation/2_Missing_Inspection/norm_imp_MinProb.rds -o Validation/3_Differential_Expression_Analysis -e 9 -u Validation/1_pre_processing/uniprot_to_genename.rds

Rscript 12Vis_Differential_expression_analysis.R -i Validation/3_Differential_Expression_Analysis/Differential_Expression_Results.rds -o Validation/12_Vis_Differential_expression_analysis -e 9
Rscript 12Vis_Differential_expression_analysis_FDR_0102.R -i Validation/3_Differential_Expression_Analysis/Differential_Expression_Results.rds -o Validation/12_Vis_Differential_expression_analysis_FDR_0102 -e 9

Rscript 13Vis_umap.R -i Validation/2_Missing_Inspection -o Validation/13_Vis_umap -e 9 -d all -l TRUE
Rscript 13Vis_umap.R -i Validation/2_Missing_Inspection -o Validation/13_Vis_umap -e 9 -d all -l FALSE

Rscript 19GSEA_allsamples.R -i Validation/3_Differential_Expression_Analysis/Differential_Expression_Results.rds -o Validation/19_GSEA_allsamples

Rscript 37_Clinical_heatmap_highest_variance.R --input Validation/2_Missing_Inspection/norm_imp_MinProb.rds --output Validation/37_Clinical_heatmap_top_variable --seed 9

########################################################### Data analysis (Validation cohort, ALS only) ###########################################################
# als only
Rscript 1Pre_processing_updated_20250731.R -i Validation/data/proteinGroups.txt -c Validation/data/clinical_info.csv -o Validation/1_pre_processing_als -e 9 --subset als -d Input/HUMAN_9606_idmapping_selected.tab.gz

Rscript 2Missing_Inspection.R -i Validation/1_pre_processing_als/se_abu_data_filtered.rds -o Validation/2_Missing_Inspection_als -e 9 -s 0.01

Rscript 5Clustering_Validation_updated_20250331.R -i Validation/2_Missing_Inspection_als/norm_imp_MinProb.rds -o Validation/5_Clustering_als -e 9

Rscript 8Differential_expression_analysis_subclusters.R -i Validation/2_Missing_Inspection_als/norm_imp_MinProb.rds -o Validation/8_Differential_expression_analysis_subclusters -e 9 -c Validation/5_Clustering_als/cluster_assignments_2.csv

Rscript 9Vis_Differential_expression_analysis_subclusters.R -i Validation/8_Differential_expression_analysis_subclusters/res.rds -o Validation/9_Vis_Differential_expression_analysis_subclusters -d Validation/2_Missing_Inspection_als/norm_imp_MinProb.rds -c Validation/5_Clustering_als/cluster_assignments_2.csv -e 9

Rscript 13Vis_umap.R -i Validation/2_Missing_Inspection_als -o Validation/13_Vis_umap -d als -e 9 -l TRUE -c Validation/5_Clustering_als
Rscript 13Vis_umap.R -i Validation/2_Missing_Inspection_als -o Validation/13_Vis_umap -d als -e 9 -l FALSE -c Validation/5_Clustering_als

Rscript 19GSEA_subclusters.R -i Validation/8_Differential_expression_analysis_subclusters/res.rds -o Validation/19_GSEA_subclusters_k2

# WGCNA using module detected in Discovery, to see it's performance in Discovery
Rscript 23_WGCNA_comparison.R -i Discovery/2_Missing_Inspection_als/norm_imp_MinProb.rds -o Discovery/23_WGCNA_comparison --net Discovery/23_WGCNA/WGCNA_net.rds -e 9 -c Discovery/5_Clustering_als/cluster_assignments_2.csv --ModuleTrait TRUE    

# WGCNA using module detected in Discovery, to see it's performance in Validation
Rscript 23_WGCNA_comparison.R -i Validation/2_Missing_Inspection_als/norm_imp_MinProb.rds -o Validation/23_WGCNA_comparison --net Discovery/23_WGCNA/WGCNA_net.rds -e 9 -c Validation/5_Clustering_als/cluster_assignments_2.csv --ModuleTrait TRUE

# sex
Rscript 38Clinical_features.R --input Validation/2_Missing_Inspection_als/norm_imp_MinProb.rds --cluster Validation/5_Clustering_als/cluster_assignments_2.csv --output Validation/38_Clinical_features --var sex

# progression group
Rscript 38Clinical_features.R --input Validation/2_Missing_Inspection_als/norm_imp_MinProb.rds --cluster Validation/5_Clustering_als/cluster_assignments_2.csv --output Validation/38_Clinical_features --var progression_group

#onset
Rscript 38Clinical_features.R --input Validation/2_Missing_Inspection_als/norm_imp_MinProb.rds --cluster Validation/5_Clustering_als/cluster_assignments_2.csv --output Validation/38_Clinical_features --var condition

#age
Rscript 38Clinical_features.R --input Validation/2_Missing_Inspection_als/norm_imp_MinProb.rds --cluster Validation/5_Clustering_als/cluster_assignments_2.csv --output Validation/38_Clinical_features --var age

#age at onset
Rscript 38Clinical_features.R --input Validation/2_Missing_Inspection_als/norm_imp_MinProb.rds --cluster Validation/5_Clustering_als/cluster_assignments_2.csv --output Validation/38_Clinical_features --var age_at_onset

#Nfl
Rscript 38Clinical_features.R --input Validation/2_Missing_Inspection_als/norm_imp_MinProb.rds --cluster Validation/5_Clustering_als/cluster_assignments_2.csv --output Validation/38_Clinical_features --var Nfl

########################################################### Between Validation and Discovery cohort ###########################################################
#UMAP between duplicated samples
#Rscript /Users/xliu2942/Documents/Projects/MAXOMOD/Script/24UMAP_duplicated_samples.R

#scatterplot
echo "Rscript 20scatterplot_FDR_signedp_withCI.R"
Rscript 20scatterplot_FDR_signedp_withCI.R -i 3_Differential_Expression_Analysis/Differential_Expression_Results.rds -o 20_scatterplot_FDR_signedp_withCI -e 9

Rscript 20scatterplot_FDR_signedp_withCI_alpha_vs_beta.R -i 8_Differential_expression_analysis_subclusters/res.rds

Rscript 20scatterplot_FDR_signedp_withCI_FDR_0102.R -i 3_Differential_Expression_Analysis/Differential_Expression_Results.rds -o 20_scatterplot_FDR_signedp_withCI_FDR_01 -e 9 --FDR_Cutoff 0.1

Rscript 20scatterplot_FDR_signedp_withCI_FDR_0102.R -i 3_Differential_Expression_Analysis/Differential_Expression_Results.rds -o 20_scatterplot_FDR_signedp_withCI_FDR_02 -e 9 --FDR_Cutoff 0.2

# GSEA heatmap
Rscript 36GSEA_IC_heatmap.R -i Discovery=Discovery/19_GSEA_allsamples/GSEA_result.rds,Validation=Validation/19_GSEA_allsamples/GSEA_result.rds -o 36_GSEA_IC_heatmap_ALS_VS_CTRL -c 6

Rscript 36GSEA_IC_heatmap.R -i Discovery=Discovery/19_GSEA_subclusters_k2/GSEA_result.rds,Validation=Validation/19_GSEA_subclusters_k2/GSEA_result.rds -o 36_GSEA_IC_heatmap_alpha_VS_beta -c 6

Rscript 25Pathway_between_cohort_GSEA_als_vs_ctrl.R
Rscript 25Pathway_between_cohort_GSEA.R

#input file for Discovery should be "Discovery/2_Missing_Inspection_als/norm_imp_MinProb.rds"
#input file for Validation should be "Validation/2_Missing_Inspection_als/norm_imp_MinProb.rds"
# ML Step by Step
Rscript 31ML_multi_models_Step_by_step_test_copy.R --n_boot 500 --topn 10 --sd_threshold 3 --output 31_ML_multi_models_step_by_step_test_copy_500

Rscript 31ML_multi_models_Script_20250811_noAge_errorline_updated_sd.R --n_boot 500 --feature_freq_cutoff 0.4 --sd_threshold 3 --output 31_ML_multi_models_20250915_500_04_noAge_3sd

# AUC line plot
Rscript 42AUC_comparision_noAge_sd.R

# For Discovery
Rscript 40Correlation_LASSO_protein_vs_Clinical.R --input Discovery/2_Missing_Inspection_als/norm_imp_MinProb.rds -s 31_ML_multi_models_20250915_500_04_noAge_3sd/selected_features.txt -c Discovery/5_Clustering_als/cluster_assignments_2.csv -g Discovery/8_Differential_expression_analysis_subclusters/Differential_expression_analysis_for_k2.csv -o 40_Correlation_LASSO_protein_vs_Clinical_Discovery 
# For Validation
Rscript 40Correlation_LASSO_protein_vs_Clinical.R --input Validation/2_Missing_Inspection_als/norm_imp_MinProb.rds -s 31_ML_multi_models_20250915_500_04_noAge_3sd/selected_features.txt -c Validation/5_Clustering_als/cluster_assignments_2.csv -g Validation/8_Differential_expression_analysis_subclusters/Differential_expression_analysis_for_k2.csv -o 40_Correlation_LASSO_protein_vs_Clinical_Validation 

########################################################### External Cohort ###########################################################
Rscript 26External_data.R

Rscript 26External_GSEA_allsamples.R -i 26External_data/als_vs_ctrl.csv -o 26External_data/26_External_GSEA_allsamples
Rscript 26External_GSEA_subclusters.R -i 26External_data/alpha_vs_beta_als.csv -o 26External_data/26_External_GSEA_subclusters_k2

Rscript 26External_protein_expression.R -i 26External_data/norm_imp_MinProb_als.rds --features 31_ML_multi_models_20250915_500_04_noAge_3sd/selected_features.txt -o 26External_data/26External_protein_expression/

# WGCNA using module detected in Discovery, to see it's performance in External Cohort
Rscript 23_WGCNA_comparison.R -i 26External_data/norm_imp_MinProb_als.rds -o 26External_data/23_WGCNA_comparison --net Discovery/23_WGCNA/WGCNA_net.rds -e 9 -c 26External_data/cluster_assignments_2.csv

Rscript 36GSEA_IC_heatmap.R -i Discovery=Discovery/19_GSEA_subclusters_k2/GSEA_result.rds,Validation=Validation/19_GSEA_subclusters_k2/GSEA_result.rds,External_data=26External_data/26_External_GSEA_subclusters_k2/GSEA_result.rds -o 36_GSEA_IC_heatmap_alpha_VS_beta_all_cohorts -c 6

########################################################### Brain data ###########################################################
Rscript 28_Brain_data_preprocessing.R

# All samples
Rscript 3Differential_expression_analysis.R -i 28_Brain_data/2_Missing_Inspection/norm_imp_MinProb.rds -o 28_Brain_data/3_Differential_Expression_Analysis -e 9 -u 28_Brain_data/1_pre_processing/uniprot_to_genename.rds

Rscript 12Vis_Differential_expression_analysis.R -i 28_Brain_data/3_Differential_Expression_Analysis/Differential_Expression_Results.rds -o 28_Brain_data/12_Vis_Differential_expression_analysis -e 9

Rscript 13Vis_umap.R -i 28_Brain_data/2_Missing_Inspection -o 28_Brain_data/13_Vis_umap -e 9 -d all -l TRUE

Rscript 13Vis_umap.R -i 28_Brain_data/2_Missing_Inspection -o 28_Brain_data/13_Vis_umap -e 9 -d all -l FALSE

Rscript 19GSEA_allsamples.R -i 28_Brain_data/3_Differential_Expression_Analysis/Differential_Expression_Results.rds -o 28_Brain_data/19_GSEA_allsamples

# als only
# this script is used to cluster the als samples, without NFL information
Rscript 5Clustering_updated_20250512_brain_clustering.R -i 28_Brain_data/2_Missing_Inspection_als/norm_imp_MinProb.rds -o 28_Brain_data/5_Clustering_als -e 9

Rscript 8Differential_expression_analysis_subclusters.R -i 28_Brain_data/2_Missing_Inspection_als/norm_imp_MinProb.rds -o 28_Brain_data/8_Differential_expression_analysis_subclusters -e 9 -c 28_Brain_data/5_Clustering_als/cluster_assignments_2.csv

# when k =3, there is no alpha or theta subcluster specific proteins, therefore, cannot draw the expression heatmap. This will return an error.
Rscript 9Vis_Differential_expression_analysis_subclusters.R -i 28_Brain_data/8_Differential_expression_analysis_subclusters/res.rds -o 28_Brain_data/9_Vis_Differential_expression_analysis_subclusters -d 28_Brain_data/2_Missing_Inspection_als/norm_imp_MinProb.rds -c 28_Brain_data/5_Clustering_als/cluster_assignments_2.csv -e 9

Rscript 13Vis_umap.R -i 28_Brain_data/2_Missing_Inspection_als -o 28_Brain_data/13_Vis_umap -d als -e 9 -l TRUE -c 28_Brain_data/5_Clustering_als
Rscript 13Vis_umap.R -i 28_Brain_data/2_Missing_Inspection_als -o 28_Brain_data/13_Vis_umap -d als -e 9 -l FALSE -c 28_Brain_data/5_Clustering_als

Rscript 19GSEA_subclusters.R -i 28_Brain_data/8_Differential_expression_analysis_subclusters/res.rds -o 28_Brain_data/19_GSEA_subclusters_k2

Rscript 36GSEA_IC_heatmap.R -i Brain=28_Brain_data/19_GSEA_subclusters_k2/GSEA_result.rds -o 28_Brain_data/36_GSEA_IC_heatmap_alpha_VS_beta -c 6