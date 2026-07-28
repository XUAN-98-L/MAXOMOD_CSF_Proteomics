#=========================Script Description=================================
# Per-protein missingness summary, with per-condition breakdown, for the
# supplementary tables:
#   Supp Table 14 - missing_per_protein_with_conditions_DC_ALSvsCTR.csv
#   Supp Table 15 - missing_per_protein_with_conditions_DC_alphavsbeta.csv
#   Supp Table 16 - missing_per_protein_with_conditions_VC_ALSvsCTR.csv
#   Supp Table 17 - missing_per_protein_with_conditions_VC_alphavsbeta.csv
#
# Reuses the missing-value counting logic from 01_Pre_Processing.R (a
# protein/sample entry is "missing" if its raw intensity is 0 or NA), on the
# same abu_data_cleaned.rds tables that script produces:
#   <cohort>/01_Pre_Processing/abu_data_cleaned.rds      (all patients, ALSvsCTR)
#   <cohort>/01_Pre_Processing_als/abu_data_cleaned.rds  (ALS patients only)
#
# ALSvsCTR condition (als/ctrl) comes from clinical_cleaned.csv. For the ALS
# subset, 01_Pre_Processing.R itself would label condition as onset
# (spinal/bulbar); per request, we instead label samples by their k=2
# subtype (alpha/beta) from 08_Clustering_als/cluster_assignments_2.csv.
#
# Output columns: UniProtName, name, ID, Missing_Count, Missing_Percentage,
# condition, Missing_Count_Per_Condition, Missing_Percentage_Per_Condition.
# Missing_Count/Missing_Percentage are computed across ALL samples in the
# table; Missing_Count_Per_Condition/Missing_Percentage_Per_Condition are
# computed within just the samples of that row's condition.
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
#===========================Function Definition=============================
# per-protein missingness (overall + per condition) for one abu_data table.
# abu_data: data.frame with columns UniProtName, name, ID, and one column
# per sample (raw intensity, 0 = missing), as saved by 01_Pre_Processing.R.
# condition_map: named character vector, names = sample column names
# (patid), values = condition label for that sample.
compute_missing_per_protein_with_conditions = function(abu_data, condition_map){
  sample_cols = intersect(names(condition_map), colnames(abu_data))
  missing_ids = setdiff(names(condition_map), colnames(abu_data))
  if (length(missing_ids) > 0){
    message(length(missing_ids), " sample(s) in the condition map were not found in abu_data and are skipped: ",
           paste(missing_ids, collapse = ", "))
  }

  num_samples = length(sample_cols)

  mat = as.matrix(abu_data[, sample_cols, drop = FALSE])
  storage.mode(mat) = "numeric"
  missing_mat = (mat == 0) | is.na(mat)

  missing_summary = data.frame(
    UniProtName = abu_data$UniProtName,
    name = abu_data$name,
    ID = abu_data$ID,
    Missing_Count = rowSums(missing_mat),
    Missing_Percentage = rowSums(missing_mat) / num_samples * 100,
    stringsAsFactors = FALSE
  )

  conditions = unique(condition_map[sample_cols])
  per_condition_list = list()
  for (cond in conditions){
    cond_cols = sample_cols[condition_map[sample_cols] == cond]
    n_cond = length(cond_cols)
    cond_missing = missing_mat[, cond_cols, drop = FALSE]
    per_condition_list[[cond]] = data.frame(
      UniProtName = abu_data$UniProtName,
      condition = cond,
      Missing_Count_Per_Condition = rowSums(cond_missing),
      Missing_Percentage_Per_Condition = rowSums(cond_missing) / n_cond * 100,
      stringsAsFactors = FALSE
    )
  }
  per_condition = do.call(rbind, per_condition_list)

  out = merge(missing_summary, per_condition, by = "UniProtName", all.x = TRUE)
  out = out[, c("UniProtName", "name", "ID", "Missing_Count", "Missing_Percentage",
               "condition", "Missing_Count_Per_Condition", "Missing_Percentage_Per_Condition")]
  out = out[order(out$UniProtName, out$condition), ]
  rownames(out) = NULL
  return(out)}
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--full_dir", "-f"),
              type = "character", default = "01_Pre_Processing",
              help = "01_Pre_Processing output dir (all patients), relative to each cohort folder."
  ),make_option(c("--als_dir", "-a"),
                type = "character", default = "01_Pre_Processing_als",
                help = "01_Pre_Processing output dir (ALS patients only), relative to each cohort folder."
  ),make_option(c("--clusters", "-c"),
                type = "character", default = "08_Clustering_als/cluster_assignments_2.csv",
                help = "cluster_assignments_2.csv (patid + kmeans_k=2 alpha/beta), relative to each cohort folder."
  ),make_option(c("--output", "-o"),
                type = "character", default = "29_missing_per_protein",
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

if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
}else{
  seed = opt$seed
  set.seed(seed)
}

cohorts = c(Discovery = "DC", Validation = "VC")

manifest = data.frame(Supp_Table = character(), file = character(), description = character(),
                      stringsAsFactors = FALSE)

################################################################################################
for (cohort in names(cohorts)){
  prefix = cohorts[[cohort]]

  #-------------------------------------------------------------------------
  # ALS vs CTR (all patients)
  abu_data_full = readRDS(file.path(cohort, opt$full_dir, "abu_data_cleaned.rds"))
  clinical_full = read.csv(file.path(cohort, opt$full_dir, "clinical_cleaned.csv"), check.names = FALSE)
  condition_map_full = setNames(as.character(clinical_full$disease), as.character(clinical_full$patid))

  table_ALSvsCTR = compute_missing_per_protein_with_conditions(abu_data_full, condition_map_full)
  file_ALSvsCTR = paste0("missing_per_protein_with_conditions_", prefix, "_ALSvsCTR.csv")
  write.csv(table_ALSvsCTR, file.path(output_dir, file_ALSvsCTR), row.names = FALSE)

  #-------------------------------------------------------------------------
  # alpha vs beta (ALS patients only, k=2 subtype instead of onset)
  abu_data_als = readRDS(file.path(cohort, opt$als_dir, "abu_data_cleaned.rds"))
  clusters = read.csv(file.path(cohort, opt$clusters), check.names = FALSE)
  condition_map_ab = setNames(as.character(clusters$`kmeans_k=2`), as.character(clusters$patid))

  table_alphavsbeta = compute_missing_per_protein_with_conditions(abu_data_als, condition_map_ab)
  file_alphavsbeta = paste0("missing_per_protein_with_conditions_", prefix, "_alphavsbeta.csv")
  write.csv(table_alphavsbeta, file.path(output_dir, file_alphavsbeta), row.names = FALSE)

  manifest = rbind(manifest,
                   data.frame(Supp_Table = NA, file = file_ALSvsCTR,
                             description = paste0(cohort, " cohort, ALS vs Control")),
                   data.frame(Supp_Table = NA, file = file_alphavsbeta,
                             description = paste0(cohort, " cohort, alpha vs beta subtype (k=2)")))
}

manifest = manifest[order(match(manifest$file, c(
  "missing_per_protein_with_conditions_DC_ALSvsCTR.csv",
  "missing_per_protein_with_conditions_DC_alphavsbeta.csv",
  "missing_per_protein_with_conditions_VC_ALSvsCTR.csv",
  "missing_per_protein_with_conditions_VC_alphavsbeta.csv"
))), ]
manifest$Supp_Table = paste("Supp Table", 13 + seq_len(nrow(manifest)))
write.csv(manifest, file.path(output_dir, "table_manifest.csv"), row.names = FALSE)
