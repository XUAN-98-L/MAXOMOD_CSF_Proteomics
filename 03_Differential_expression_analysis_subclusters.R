#=========================Script Description=================================
# This script is used for differential expression analysis between subclusters
# Rscript 03_Differential_expression_analysis_subclusters.R -i Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds -o Discovery/03_Differential_expression_analysis_subclusters -e 9 -c Discovery/08_Clustering_als/cluster_assignments_2.csv
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("purrr"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("DEP"))
suppressMessages(library("data.table"))
#===========================Function Definition=============================
test_diff = function(se, type = c("control", "all", "manual"),
                     control = NULL, test = NULL,
                     design_formula = formula(~ 0 + condition)){
  assertthat::assert_that(inherits(se, "SummarizedExperiment"), 
                          is.character(type), class(design_formula) == "formula")
  type <- match.arg(type)
  col_data <- colData(se)
  raw <- assay(se)
  if (any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '", 
         deparse(substitute(se)), "'\nRun make_unique() and make_se() to obtain the required columns", 
         call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '", 
         deparse(substitute(se)), "'\nRun make_se() or make_se_parse() to obtain the required columns", 
         call. = FALSE)
  }
  if (any(is.na(raw))) {
    warning("Missing values in '", deparse(substitute(se)), 
            "'")
  }
  if (!is.null(control)) {
    assertthat::assert_that(is.character(control), length(control) == 
                              1)
    if (!control %in% unique(col_data$condition)) {
      stop("run test_diff() with a valid control.\nValid controls are: '", 
           paste0(unique(col_data$condition), collapse = "', '"), 
           "'", call. = FALSE)
    }
  }
  variables <- terms.formula(design_formula) %>% attr(., "variables") %>% 
    as.character() %>% .[-1]
  if (any(!variables %in% colnames(col_data))) {
    stop("run make_diff() with an appropriate 'design_formula'")
  }
  if (variables[1] != "condition") {
    stop("first factor of 'design_formula' should be 'condition'")
  }
  for (var in variables) {
    temp <- factor(col_data[[var]])
    assign(var, temp)
  }
  #design <- model.matrix(design_formula, data = environment())
  design <- model.matrix(design_formula, data = col_data)
  colnames(design) <- gsub("condition", "", colnames(design))
  conditions <- as.character(unique(condition))
  if (type == "all") {
    cntrst <- apply(utils::combn(conditions, 2), 2, paste, 
                    collapse = " - ")
    if (!is.null(control)) {
      flip <- grep(paste("^", control, sep = ""), cntrst)
      if (length(flip) >= 1) {
        cntrst[flip] <- cntrst[flip] %>% gsub(paste(control, 
                                                    "- ", sep = " "), "", .) %>% paste(" - ", control, 
                                                                                       sep = "")
      }
    }
  }
  if (type == "control") {
    if (is.null(control)) 
      stop("run test_diff(type = 'control') with a 'control' argument")
    cntrst <- paste(conditions[!conditions %in% control], 
                    control, sep = " - ")
  }
  if (type == "manual") {
    if (is.null(test)) {
      stop("run test_diff(type = 'manual') with a 'test' argument")
    }
    assertthat::assert_that(is.character(test))
    if (any(!unlist(strsplit(test, "_vs_")) %in% conditions)) {
      stop("run test_diff() with valid contrasts in 'test'", 
           ".\nValid contrasts should contain combinations of: '", 
           paste0(conditions, collapse = "', '"), "', for example '", 
           paste0(conditions[1], "_vs_", conditions[2]), 
           "'.", call. = FALSE)
    }
    cntrst <- gsub("_vs_", " - ", test)
  }
  message("Tested contrasts: ", paste(gsub(" - ", "_vs_", cntrst), 
                                      collapse = ", "))
  fit <- limma::lmFit(raw, design = design)
  made_contrasts <- limma::makeContrasts(contrasts = cntrst, levels = design)
  contrast_fit <- limma::contrasts.fit(fit, made_contrasts)
  if (any(is.na(raw))) {
    for (i in cntrst) {
      covariates <- strsplit(i, " - ") %>% unlist
      single_contrast <- limma::makeContrasts(contrasts = i, levels = design[, 
                                                                             covariates])
      single_contrast_fit <- limma::contrasts.fit(fit[, covariates], 
                                                  single_contrast)
      contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[, 
                                                                         1]
      contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[, 
                                                                             1]
    }
  }
  eB_fit <- limma::eBayes(contrast_fit)
  retrieve_fun <- function(comp, fit = eB_fit) {
    res <- limma::topTable(fit, sort.by = "t", coef = comp, number = Inf, 
                           confint = TRUE)
    res <- res[!is.na(res$t), ]
    fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
    res$qval <- fdr_res$qval
    res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    res <- tibble::rownames_to_column(res)
    return(res)
  }
  limma_res <- map_df(cntrst, retrieve_fun)
  table <- limma_res %>% dplyr::select(rowname, logFC, CI.L, CI.R, 
                                       P.Value, qval, comparison) %>% dplyr::mutate(comparison = gsub(" - ", 
                                                                                                      "_vs_", comparison)) %>% tidyr::gather(variable, value, -c(rowname, 
                                                                                                                                                                 comparison)) %>% dplyr::mutate(variable = recode(variable, logFC = "diff", 
                                                                                                                                                                                                                  P.Value = "p.val", qval = "p.adj")) %>% tidyr::unite(temp, comparison, 
                                                                                                                                                                                                                                                                       variable) %>% tidyr::spread(temp, value)
  rowData(se) <- merge(rowData(se, use.names = FALSE), table, 
                                   by.x = "name", by.y = "rowname", all.x = TRUE, sort = FALSE)
  return(se)
}
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "02_Missing_Inspection_subclusters/norm_imp_MinProb.rds",
              help = "SummarizedExperiment object after normalize and imputation."
  ),make_option(c("--output", "-o"),
                type = "character", default = "03_Differential_expression_analysis_subclusters",
                help = "output directory path."
  ),make_option(c("--seed", "-e"),
                type = "integer", default = 9,
                help = "set.seed"
  ),make_option(c("--cluster_assignments", "-c"), 
                type = "character", default = "08_Clustering_als/cluster_assignments_2.csv",
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
  data = readRDS(input)
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
  cluster_assignments = read.csv(opt$cluster_assignments,check.names = F)
}
#===========================Differential expression analysis=============================
covariates_f =  ~0 + condition + age + sex
res = list() #empty list to save the results

#control = "beta"

# make cluster_assignments  the same oder of colData(data)$label
cluster_assignments = cluster_assignments[match(colData(data)$label, cluster_assignments$patid),]
#===========================Differential expression analysis=============================

############################################################################################################
#for k = 2
message("Differential expression analysis for k=2")
title = "Differential_expression_analysis_for_k2"
l = "k2"
se_abu_data_ALS = data
# add cluster information
colData(se_abu_data_ALS)$cluster = as.factor(cluster_assignments$`kmeans_k=2`)
colData(se_abu_data_ALS)$condition = as.character(colData(se_abu_data_ALS)$cluster)
experimental.design <- colData(se_abu_data_ALS) %>% as.data.frame() %>% 
  group_by(condition) %>%  # Group by 'condition'
  mutate(replicate = ave(condition, condition, FUN = seq_along)) 
colData(se_abu_data_ALS)$replicate = experimental.design$replicate
colData(se_abu_data_ALS)$ID = paste0(colData(se_abu_data_ALS)$condition, "_", colData(se_abu_data_ALS)$replicate)
# rename colData in se_abu_data_ALS using cluster information as replicates
colnames(se_abu_data_ALS) = colData(se_abu_data_ALS)$ID

#perform DEx
d = se_abu_data_ALS
t = test_diff(d, type = "all", control = "beta", 
              test = NULL, design_formula = formula(covariates_f))
res[[l]] = as.data.frame(t@elementMetadata@listData)
pval_columns = grep("p.val", colnames(res[[l]]))
pval_names = colnames(res[[l]])[pval_columns]
fdr_names = gsub("p.val", "fdr", pval_names)
if(length(pval_columns)==1){
  res[[l]]$fdr = p.adjust(res[[l]][,pval_columns], method="BH") #calculate fdr with BH adjustment
}else{
  #if we have 3 clusters we have three different comparisons
  res[[l]]$fdr.1 = p.adjust(res[[l]][,pval_columns[1]], method="BH") 
  res[[l]]$fdr.2 = p.adjust(res[[l]][,pval_columns[2]], method="BH")
  res[[l]]$fdr.3 = p.adjust(res[[l]][,pval_columns[3]], method="BH")
  colnames(res[[l]])[(ncol(res[[l]])-2) : ncol(res[[l]])] = fdr_names
}
write.csv(res[[l]],quote = F,row.names = F, file = paste0(output_dir,"/", title, ".csv"))


############################################################################################################
#for k = 3
message("Differential expression analysis for k=3")
title = "Differential_expression_analysis_for_k3"
l = "k3"
se_abu_data_ALS = data
# add cluster information
colData(se_abu_data_ALS)$cluster = as.factor(cluster_assignments$`kmeans_k=3`)
colData(se_abu_data_ALS)$condition = as.character(colData(se_abu_data_ALS)$cluster)
experimental.design <- colData(se_abu_data_ALS) %>% as.data.frame() %>% 
  group_by(condition) %>%  # Group by 'condition'
  mutate(replicate = ave(condition, condition, FUN = seq_along)) 
colData(se_abu_data_ALS)$replicate = experimental.design$replicate
colData(se_abu_data_ALS)$ID = paste0(colData(se_abu_data_ALS)$condition, "_", colData(se_abu_data_ALS)$replicate)
# rename colData in se_abu_data_ALS using cluster information as replicates
colnames(se_abu_data_ALS) = colData(se_abu_data_ALS)$ID

#perform DEx
d = se_abu_data_ALS
t = test_diff(d, type = "all", control = "beta",
              test = NULL, design_formula = formula(covariates_f))
res[[l]] = as.data.frame(t@elementMetadata@listData)
pval_columns = grep("p.val", colnames(res[[l]]))
pval_names = colnames(res[[l]])[pval_columns]
fdr_names = gsub("p.val", "fdr", pval_names)
if(length(pval_columns)==1){
  res[[l]]$fdr = p.adjust(res[[l]][,pval_columns], method="BH") #calculate fdr with BH adjustment
}else{
  #if we have 3 clusters we have three different comparisons
  res[[l]]$fdr.1 = p.adjust(res[[l]][,pval_columns[1]], method="BH") 
  res[[l]]$fdr.2 = p.adjust(res[[l]][,pval_columns[2]], method="BH")
  res[[l]]$fdr.3 = p.adjust(res[[l]][,pval_columns[3]], method="BH")
  colnames(res[[l]])[(ncol(res[[l]])-2) : ncol(res[[l]])] = fdr_names
}
write.csv(res[[l]],quote = F,row.names = F, file = paste0(output_dir,"/", title, ".csv"))

saveRDS(res, file = paste0(output_dir,"/","res.rds")) #save the results as an RDS file
