#=========================Script Description=================================
# Differential expression analysis of ALS subtypes vs Control: alpha vs
# ctrl and beta vs ctrl. Structured the same way as
# 03_Differential_expression_analysis.R (same test_diff() function, same
# covariate/patient-stratification grid), but the two-level condition
# (als/ctrl) is replaced by a three-level condition (ctrl/alpha/beta):
# control samples keep condition = ctrl, ALS samples are relabelled to
# their k=2 subtype (alpha/beta) from cluster_assignments_2.csv. Theta and
# unlabeled ALS samples are dropped, since only alpha/beta subtypes exist
# at k=2.
#
# test_diff(type = "control", control = "ctrl") on this three-level factor
# produces exactly two contrasts: "alpha - ctrl" and "beta - ctrl" (not
# alpha vs beta), matching what 03_Differential_expression_analysis.R does
# for als vs ctrl.
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
              type = "character", default = "02_Missing_Inspection/norm_imp_MinProb.rds",
              help = "normalised, imputed SummarizedExperiment (all patients: control + ALS)."
  ),make_option(c("--cluster_assignments", "-c"),
                type = "character", default = "08_Clustering_als/cluster_assignments_2.csv",
                help = "cluster_assignments_2.csv (patid + kmeans_k=2 alpha/beta subtype)."
  ),make_option(c("--output", "-o"),
                type = "character", default = "32_subcluster_vsCTR",
                help = "output directory path."
  ),make_option(c("--seed", "-e"),
                type = "integer", default = 9,
                help = "set.seed"
  ),make_option(c("--uniprot_to_genename", "-u"),
                type = "character", default = "01_Pre_Processing/uniprot_to_genename.rds",
                help = "uniprot_to_genename.rds"
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

if (is.null(opt$uniprot_to_genename)) {
  stop("Please provide the uniprot_to_genename file path!")
}else if (!file.exists(opt$uniprot_to_genename)) {
  stop("uniprot_to_genename file does not exist!")
}else{
  uniprot_to_genename = readRDS(opt$uniprot_to_genename)
}

if (is.null(opt$cluster_assignments)) {
  stop("Please provide the cluster_assignments_2.csv file path!")
}else if (!file.exists(opt$cluster_assignments)) {
  stop("cluster_assignments_2.csv does not exist!")
}else{
  cluster_assignments = read.csv(opt$cluster_assignments, check.names = FALSE)
}
#===========================Relabel condition: ctrl / alpha / beta=========================
subgroup = setNames(as.character(cluster_assignments$`kmeans_k=2`), cluster_assignments$patid)

new_condition = ifelse(colData(data)$condition == "ctrl", "ctrl", subgroup[colData(data)$label])
keep = !is.na(new_condition)
if (any(!keep)){
  message("Dropping ", sum(!keep), " sample(s) with no ctrl/alpha/beta label (e.g. theta): ",
         paste(colData(data)$label[!keep], collapse = ", "))
}

data = data[, keep]
colData(data)$condition = factor(new_condition[keep], levels = c("ctrl", "alpha", "beta"))

message("Sample counts by condition:")
print(table(colData(data)$condition))
#===========================Differential expression analysis=============================
covariates = c("no_cov", "age_cov", "sex_cov", "age_sex_cov") #a vector for the names of the covariates
covariates_f = c(~0 + condition,  #the actual formulas to test with different combinations of covariates
                 ~0 + condition + age, #because this function only takes categorical covariates i'm using the categorical age variable
                 ~0 + condition + sex,
                 ~0 + condition + age + sex)
patients = c("all_patients", "only_female", "only_male") #names for the analysis without stratification, with only females, and with only males
patients_f = c(NA, "Female", "Male") #the actual values that need to be used in the function to do the stratification
res = list() #empty list to save the results

control = "ctrl"

#loop to perform the different types of differential expression analysis
l = 1
for(i in 1:length(covariates)){ #a loop for the different covariates
  for(j in 1:length(patients)){ #a loop for sex stratification and no sex stratification

    title = paste0("norm_imp_MinProb_subcluster_",covariates[i], "_", patients[j]) #construct a title where all the varaibles are stored
    d = data #select the k'th dataset
    if(j >1) { d = d[,d$sex == patients_f[j]] } #if j>1 we will apply sex stratification, and this line of code will select either only females or only males
    if(i < 3 | j == 1){ #this if() command is to prevent the differential expression analysis from running the test on stratified subsets with sex stratification because the sample sizes are too small for this
      print(title)
      print(dim(d))
      t = test_diff(d, type = "control", control = control, #this function performs the tests (alpha_vs_ctrl, beta_vs_ctrl) and stores the results in 't'
                    test = NULL, design_formula = formula(covariates_f[[i]]))
      res[[l]] = as.data.frame(t@elementMetadata@listData) #the statistics are in @elementMetadata@listData and we will save this part of the results to our 'res' results list
      pval_columns = grep("p.val", colnames(res[[l]]))
      pval_names = colnames(res[[l]])[pval_columns]
      fdr_names = gsub("p.val", "fdr", pval_names)
      if(length(pval_columns)==1){
        res[[l]]$fdr = p.adjust(res[[l]][,pval_columns], method="BH") #we add an FDR column by using Benjamini Hochberg correction on the "p.val" column
      }else{
        #alpha_vs_ctrl and beta_vs_ctrl each get their own FDR column
        for(p in seq_along(pval_columns)){
          res[[l]][[fdr_names[p]]] = p.adjust(res[[l]][,pval_columns[p]], method="BH")
        }
      }
      print(dim(res[[l]]))
      names(res)[l] = title
      res[[l]] <- left_join(res[[l]],
                            uniprot_to_genename %>% dplyr::select(UniProtAccession, gene_name),
                            by = c("ID" = "UniProtAccession")) #join the gene names to the results
      write.csv(res[[l]],quote = F,row.names = F, file = paste0(output_dir,"/", title, ".csv")) #write the results in a CSV file to the results folder. Each loop will create a CSV file
      l = l+1
    }}}
saveRDS(res, file = paste0(output_dir,"/","Differential_Expression_Results.rds")) #save the results as an RDS file
