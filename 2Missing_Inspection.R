#=========================Script Description=================================
# This script is used for make pre-processing of the data, including normalization and imputation
# cd /Users/xliu2942/Documents/Projects/MAXOMOD/TESTing/Pipeline
# Rscript 2Missing_Inspection.R -i 1_pre_processing/se_abu_data_filtered.rds -o 2_Missing_Inspection -e 9 -s 0.01 >output.log
#===========================Loading Packages==================================
suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("DEP"))
#===========================Function Definition================================

#===========================Command Parameters Setting=========================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "1_pre_processing/se_abu_data_filtered.rds",
              help = "filtered SummarizedExperiment object file path."
  ),make_option(c("--output", "-o"),
                type = "character", default = "2_Missing_Inspection",
                help = "output directory path."
  ),make_option(c("--seed", "-e"),
                type = "integer", default = 9,
                help = "set.seed"
  ),make_option(c("--scalar", "-s"),
                type = "numeric", default = 0.01,
                help = "A scalar used to determine a low expression value to be used for missing data imputation. 0 < q < 1, in this case q should be set to a low value. The default value is q = 0.01."
  ),make_option(c("--normalization", "-n"),
                type = "logical", default = TRUE,
                help = "Whether to perform vsn normalization. Default is TRUE."
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
  stop("Please provide the cleaned proteinGroups file path!")
}else if (!file.exists(opt$input)) {
  stop("The cleaned proteinGroups file does not exist!")
}else{
  input = opt$input
  se_abu_data_filtered = readRDS(input)
}

if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
}else{
  seed = opt$seed
  set.seed(seed)
}

#===========================Data Pre-processing==============================

# Normalization
if (opt$normalization){
  norm <- normalize_vsn(se_abu_data_filtered)
}else{
  norm <- se_abu_data_filtered
}

#norm <- normalize_vsn(se_abu_data_filtered)
pdf(file.path(output_dir, "meanSdPlot.pdf"))
meanSdPlot(norm)
dev.off()

pdf(file.path(output_dir, "Normalization_Boxplot.pdf"),height=12)
plot_normalization(se_abu_data_filtered,norm)
dev.off()

# Imputation
norm_imp_MinProb <- impute(norm, fun = "MinProb", q= opt$scalar)

saveRDS(norm, file.path(output_dir, "norm.rds"))
saveRDS(norm_imp_MinProb, file.path(output_dir, "norm_imp_MinProb.rds"))
