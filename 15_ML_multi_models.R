#=========================Script Description=================================
# This script is used for ML model building and evaluation using the Discovery and validation cohort
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("glmnet")) # LASSO regression
suppressMessages(library("pROC")) # AUC calculation
suppressMessages(library("caret")) # Confusion matrix, performance metrics
suppressMessages(library("boot")) # Bootstrap sampling
suppressMessages(library("DEP"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("randomForest"))
suppressMessages(library("e1071"))# Support Vector Machine (SVM)
suppressMessages(library("xgboost"))
suppressMessages(library("nnet")) # Logistic Regression support
suppressMessages(library("ggplot2"))
suppressMessages(library("ggpubr"))
suppressMessages(library("gridExtra"))
suppressMessages(library("tidyverse"))
suppressMessages(library("dplyr")) # For data manipulation
suppressMessages(library("RColorBrewer")) # For data manipulation
suppressMessages(library("ggrepel")) # For better text label placement in plots
suppressMessages(library("nortest"))
#===========================Function Definition=============================
# normalise expression values 
scale_manual <- function(df) {
  as.data.frame(apply(df, 2, function(x) (x - min(x)) / diff(range(x))))
}

# Calculate ROC and return the set of points (with endpoints added) + AUC ）
get_roc_points <- function(y_true, y_score, run_id) {
  roc_i <- pROC::roc(response = y_true, predictor = y_score, quiet = TRUE)
  Sp <- roc_i$specificities
  Sn <- roc_i$sensitivities
  eps <- 1e-12
  add_start <- !(any(abs(Sp - 1) < eps & abs(Sn - 0) < eps))
  add_end   <- !(any(abs(Sp - 0) < eps & abs(Sn - 1) < eps))
  if (add_start) { Sp <- c(1, Sp); Sn <- c(0, Sn) }
  if (add_end)   { Sp <- c(Sp, 0); Sn <- c(Sn, 1) }
  df <- data.frame(Sp = Sp, Sn = Sn, n = seq_along(Sn), run = run_id)
  df$auc <- as.numeric(pROC::auc(roc_i))
  df
}

get_probs <- function(fit, newx, y_true) {
  probs <- predict(fit, newx, type = "prob")
  pos <- if (is.factor(y_true)) levels(y_true)[2] else colnames(probs)[2]
  if (!pos %in% colnames(probs)) pos <- colnames(probs)[2]
  probs[[pos]]
}
get_class <- function(fit, newx) {
  as.character(predict(fit, newx, type = "raw"))
}

runML_with_plot <- function(discovery_df, validation_df,
                            scale_func = scale_manual,
                            BS_number = 100, cv = 10,
                            plot_mean_roc = TRUE, output_dir = NULL) {
  
  #models <- c("logit","rf","svm_radial","svm_linear","xgb")
  models <- c("rf","svm_radial","svm_linear","xgb","ridge","lasso")
  
  auc_list <- setNames(lapply(models, function(m) numeric(BS_number)), models)
  rocc_store <- setNames(vector("list", length(models)), models)  # ROC points storage
  names(rocc_store) <- models
  for (m in models) rocc_store[[m]] <- data.frame()  # Sp, Sn, n, run
  
  all_models <- vector("list", BS_number)       # Model each iteration
  all_predictions <- vector("list", BS_number)  # Predicted probabilities and class labels on the validation set for each iteration
  
  train_control <- trainControl(method = "cv",
                                number = cv,
                                classProbs = TRUE,
                                summaryFunction = twoClassSummary,
                                savePredictions = TRUE)
  
  # bootstrap iterations
  for (i in seq_len(BS_number)) {
    message("Running bootstrap iteration ", i)
    
    # Stratified bootstrap
    class_levels <- levels(discovery_df$class)
    max_n <- max(table(discovery_df$class))
    sampled_indices <- unlist(lapply(class_levels, function(cl) {
      sample(which(discovery_df$class == cl), size = max_n, replace = TRUE)
    }))
    
    train_bs <- discovery_df[sampled_indices, ]
    train_scaled <- scale_func(train_bs[, -ncol(train_bs)])
    train_scaled$class <- train_bs$class
    
    # Validation
    X_valid_scaled <- scale_func(validation_df[, -ncol(validation_df)])
    y_valid <- validation_df$class
    
    # Models
    #logit_model <- train(class ~ ., data = train_scaled, method = "glm", family = "binomial",
    #                     metric = "ROC", trControl = train_control, tuneLength = 5)

    p <- ncol(train_scaled) 
    rf_grid <- expand.grid(mtry = 1:p)
    
    rf_model <- train(class ~ ., data = train_scaled, method = "rf",
                      metric = "ROC", trControl = train_control, tuneLength = 5, tuneGrid = rf_grid, ntree = 500)
    svm_radial_model <- train(class ~ ., data = train_scaled, method = "svmRadial",
                              metric = "ROC", trControl = train_control, tuneLength = 5)
    svm_linear_model <- train(class ~ ., data = train_scaled, method = "svmLinear",
                              metric = "ROC", trControl = train_control, tuneLength = 5)
    xgb_model <- train(class ~ ., data = train_scaled, method = "xgbTree",
                       metric = "ROC", trControl = train_control, tuneLength = 5)
    
    # ridge regression
    ridgeGrid <- expand.grid(alpha = 0,
                         lambda = 10^seq(-3, 1, length = 100)) 
    ridge_model <- train(class ~ ., data = train_scaled, method = "glmnet",
                         metric = "ROC", trControl = train_control, tuneLength = 5, tuneGrid = ridgeGrid)

    # lasso regression
    lassoGrid <- expand.grid(alpha = 1,
                         lambda = 10^seq(-3, 1, length = 100)) 
    lasso_model <- train(class ~ ., data = train_scaled, method = "glmnet",
                         metric = "ROC", trControl = train_control, tuneLength = 5, tuneGrid = lassoGrid)


    # save the model
    all_models[[i]] <- list(
      #logit = logit_model,
      rf = rf_model,
      svm_radial = svm_radial_model,
      svm_linear = svm_linear_model,
      xgb = xgb_model,
      ridge = ridge_model,
      lasso = lasso_model
    )
    
    # Prediction: probability + class.
    #p_logit <- get_probs(logit_model, X_valid_scaled, y_valid)
    p_ridge <- get_probs(ridge_model, X_valid_scaled, y_valid)
    p_lasso <- get_probs(lasso_model, X_valid_scaled, y_valid)
    p_rf    <- get_probs(rf_model,    X_valid_scaled, y_valid)
    p_svr   <- get_probs(svm_radial_model, X_valid_scaled, y_valid)
    p_svl   <- get_probs(svm_linear_model, X_valid_scaled, y_valid)
    p_xgb   <- get_probs(xgb_model,   X_valid_scaled, y_valid)
    
    #c_logit <- get_class(logit_model, X_valid_scaled)
    c_ridge <- get_class(ridge_model, X_valid_scaled)
    c_lasso <- get_class(lasso_model, X_valid_scaled)
    c_rf    <- get_class(rf_model,    X_valid_scaled)
    c_svr   <- get_class(svm_radial_model, X_valid_scaled)
    c_svl   <- get_class(svm_linear_model, X_valid_scaled)
    c_xgb   <- get_class(xgb_model,   X_valid_scaled)
    
    # Store the prediction results of the current iteration.
    all_predictions[[i]] <- data.frame(
      run = i,
      truth = as.character(y_valid),
      #logit_prob = p_logit,  logit_pred = c_logit,
      rf_prob    = p_rf,     rf_pred    = c_rf,
      svm_radial_prob = p_svr, svm_radial_pred = c_svr,
      svm_linear_prob = p_svl, svm_linear_pred = c_svl,
      xgb_prob   = p_xgb,    xgb_pred   = c_xgb,
      ridge_prob = p_ridge, ridge_pred = c_ridge,
      lasso_prob = p_lasso, lasso_pred = c_lasso,
      stringsAsFactors = FALSE
    )
    
    # ROC points + AUC
    #df_logit <- get_roc_points(y_valid, p_logit, i)
    df_rf    <- get_roc_points(y_valid, p_rf,    i)
    df_svr   <- get_roc_points(y_valid, p_svr,   i)
    df_svl   <- get_roc_points(y_valid, p_svl,   i)
    df_xgb   <- get_roc_points(y_valid, p_xgb,   i)
    df_ridge <- get_roc_points(y_valid, p_ridge, i)
    df_lasso <- get_roc_points(y_valid, p_lasso, i)
    
    #rocc_store$logit      <- rbind(rocc_store$logit,      df_logit[, c("Sp","Sn","n","run")])
    rocc_store$rf         <- rbind(rocc_store$rf,         df_rf[,    c("Sp","Sn","n","run")])
    rocc_store$svm_radial <- rbind(rocc_store$svm_radial, df_svr[,   c("Sp","Sn","n","run")])
    rocc_store$svm_linear <- rbind(rocc_store$svm_linear, df_svl[,   c("Sp","Sn","n","run")])
    rocc_store$xgb        <- rbind(rocc_store$xgb,        df_xgb[,   c("Sp","Sn","n","run")])
    rocc_store$ridge      <- rbind(rocc_store$ridge,      df_ridge[, c("Sp","Sn","n","run")])
    rocc_store$lasso      <- rbind(rocc_store$lasso,      df_lasso[, c("Sp","Sn","n","run")])

    #auc_list$logit[i]      <- df_logit$auc[1]
    auc_list$rf[i]         <- df_rf$auc[1]
    auc_list$svm_radial[i] <- df_svr$auc[1]
    auc_list$svm_linear[i] <- df_svl$auc[1]
    auc_list$xgb[i]        <- df_xgb$auc[1]
    auc_list$ridge[i]      <- df_ridge$auc[1]
    auc_list$lasso[i]      <- df_lasso$auc[1]
  }
  
  # Aggregate: mean ROC + SD, with fallback to include (0,0) and (1,1) points.
  summary_df <- bind_rows(lapply(names(rocc_store), function(m){
    rocc <- rocc_store[[m]]
    Sp_mean <- aggregate(Sp ~ n, rocc, mean)$Sp
    Sn_mean <- aggregate(Sn ~ n, rocc, mean)$Sn
    Sp_sd   <- aggregate(Sp ~ n, rocc, sd)$Sp
    Sn_sd   <- aggregate(Sn ~ n, rocc, sd)$Sn
    
    out <- data.frame(
      method   = m,
      FPR      = 1 - Sp_mean,
      TPR      = Sn_mean,
      FPR_sd   = Sp_sd,
      TPR_sd   = Sn_sd,
      mean_auc = mean(auc_list[[m]], na.rm = TRUE),
      sd_auc   = sd(auc_list[[m]],   na.rm = TRUE)
    )
    
    eps <- 1e-12
    if (!any(out$FPR <= eps & out$TPR <= eps)) {
      out <- rbind(
        data.frame(method=m,FPR=0,TPR=0,FPR_sd=0,TPR_sd=0,
                   mean_auc=out$mean_auc[1], sd_auc=out$sd_auc[1]),
        out
      )
    }
    if (!any(abs(out$FPR-1)<=eps & abs(out$TPR-1)<=eps)) {
      out <- rbind(
        out,
        data.frame(method=m,FPR=1,TPR=1,FPR_sd=0,TPR_sd=0,
                   mean_auc=out$mean_auc[1], sd_auc=out$sd_auc[1])
      )
    }
    out[order(out$FPR, out$TPR), ]
  }))
  
  # Draw the mean ROC curve with shaded bands (mean ± 0.95 × SD)
  if (plot_mean_roc) {
    auc_text <- summary_df |>
      group_by(method) |>
      summarise(txt = sprintf("%s: %.3f±%.3f",
                              unique(method), unique(mean_auc), unique(sd_auc)),
                .groups = "drop") |>
      pull(txt) |>
      paste(collapse = "\n")
    
    p <- ggplot(summary_df, aes(x = FPR, y = TPR, color = method, fill = method)) +
      geom_ribbon(aes(ymin = pmax(0, TPR - 0.95*TPR_sd),
                      ymax = pmin(1, TPR + 0.95*TPR_sd),
                      xmin = pmax(0, FPR - 0.95*FPR_sd),
                      xmax = pmin(1, FPR + 0.95*FPR_sd)),
                  alpha = 0.25, color = NA) +
      geom_line(size = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.3) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      labs(x = "False Positive Rate (1 - Specificity)",
           y = "True Positive Rate (Sensitivity)",
           title = "Mean ± SD and the 95% confidence interval shaded",
           subtitle = paste0("AUC (mean±SD):\n", auc_text)) +
      scale_color_brewer(palette = "Set2") +  # line
      scale_fill_brewer(palette = "Set2") +   # shading
      theme_classic(base_size = 12) +
      theme(panel.grid.minor = element_blank(),
            legend.title = element_blank())
    print(p)
    
    if (!is.null(output_dir)) {
      ggsave(file.path(output_dir, "mean_ROC_plot.png"), plot = p, width = 8, height = 6, dpi = 300)
      ggsave(file.path(output_dir, "mean_ROC_plot.pdf"), plot = p, width = 8, height = 6)
      message("Mean ROC plot saved to: ", normalizePath(output_dir))
    }
  }
  
  invisible(list(
    auc_list        = auc_list,
    rocc_store      = rocc_store,
    roc_summary     = summary_df,
    all_models      = all_models,
    all_predictions = all_predictions
  ))
}

# calculate varImp and draw barplot（mean ± SD）
plot_model_importance <- function(ml_results, model_key, scale_importance = FALSE, output_dir = NULL) {
  
  # from ML_results, extract variable importance for the specified model_key
  imp_list <- lapply(ml_results$all_models, function(models) {
    fit <- models[[model_key]]
    imp <- caret::varImp(fit, scale = scale_importance)$importance
    #imp_df <- data.frame(Protein = rownames(imp), Importance = imp$Overall, stringsAsFactors = FALSE)
    imp_df <- data.frame(
      Protein = rownames(imp),
      Importance = rowMeans(imp, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    return(imp_df)
  })
  
  # into one df
  imp_df <- bind_rows(imp_list, .id = "Run")
  
  # mean and sd of importance for each protein
  imp_summary <- imp_df %>%
    group_by(Protein) %>%
    summarise(mean_importance = mean(Importance, na.rm = TRUE),
              sd_importance = sd(Importance, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(desc(mean_importance))
  
  p <- ggplot(imp_summary, aes(x = reorder(Protein, mean_importance), 
                               y = mean_importance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_errorbar(aes(ymin = mean_importance - sd_importance,
                      ymax = mean_importance + sd_importance),
                  width = 0.4, color = "black") +
    coord_flip() +
    labs(title = paste("Variable Importance (", model_key, ")", sep = ""),
         x = "Protein", 
         y = "Mean Importance ± SD") +
    theme_classic(base_size = 12)
  
  # output
  #print(p)
  if (!is.null(output_dir)) {
    ggsave(file.path(output_dir, paste0("VarImp_", model_key, ".pdf")),
           plot = p, width = max(7, nrow(imp_summary)/4), height = max(6,length(ml_results$all_models[[1]][[model_key]]$coefnames)/2))
  }
  
  # return the summary data frame
  return(imp_summary)
}
# ===========================Command Parameters Setting=============================
option_list <- list(
make_option(c("--output", "-o"),
                type = "character", default = "15_ML_multi_models",
                help = "output directory path."
  ),make_option(c("--seed", "-e"),
                type = "integer", default = 9,
                help = "set.seed"
  ),make_option(c("--n_boot", "-b"),
                type = "integer", default = 10,
                help = "Number of bootstrapping."
  ),make_option(c("--feature_freq_cutoff", "-f"),
                type = "numeric", default = 0.4,
                help = "Feature frequency cutoff in LASSO (bootstrapping)."
  ),make_option(c("--lambda", "-l"),
                type = "character", default = "min",
                help = "lambda value for LASSO regression, options are 'min' or '1se'. Default is 'min'."
  ),make_option(c("--sd_threshold", "-s"),
                type = "integer", default = 3,
                help = "Standard deviation threshold for highly age correlated proteins. Default is 3."
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

if (is.null(opt$n_boot)) {
  stop("Please provide the number of bootstrapping!")
}else{
  n_boot = opt$n_boot
}

if (is.null(opt$sd_threshold)) {
  stop("Please provide the standard deviation threshold for highly age correlated proteins!")
}else{
  sd_threshold = opt$sd_threshold
}

if (is.null(opt$feature_freq_cutoff)) {
  stop("Please provide the feature frequency cutoff in LASSO!")
}else{
  feature_freq_cutoff = opt$feature_freq_cutoff
}

if (is.null(opt$lambda)) {
  print("NO LAMBDA VALUE SUPPLIED, default lambda value will be used!")
  lambda_value <- "min"
} else {
  lambda_value <- opt$lambda
}

#============================================================================
Discovery = readRDS("Discovery/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds")
Validation = readRDS("Validation/02_Missing_Inspection_subclusters/norm_imp_MinProb.rds")

inter = intersect(Discovery@elementMetadata$name,Validation@elementMetadata$name)

# based on DEP list in discovery cohort, find differential expressed proteins
df = read.csv("Discovery/03_Differential_expression_analysis_subclusters/Differential_expression_analysis_for_k2.csv")
DEP = df$name[df$fdr < 0.05]
DEP = intersect(DEP,inter)
inter = DEP

Discovery_als_ctrl = readRDS("Discovery/02_Missing_Inspection/norm_imp_MinProb.rds")
Validation_als_ctrl = readRDS("Validation/02_Missing_Inspection/norm_imp_MinProb.rds")
assay_All = assay(Discovery_als_ctrl)
colnames(assay_All) = Discovery_als_ctrl$label
assay_All = assay_All[inter,]
# which protein are highly correlated with age (pearson correlation > 0.4 or < -0.4)
meta = as.data.frame(colData(Discovery_als_ctrl))
age_cor = c()
for (i in 1:nrow(assay_All)){
  cor_i = cor.test(as.numeric(assay_All[i, ]),meta$age,method = "pearson")
  age_cor = c(age_cor,cor_i$estimate)
}
names(age_cor) = rownames(assay_All)

#first do a normal distribution test
shapiro.test(age_cor)
# library(nortest)
ad.test(age_cor)

# this is normally distributed
# highly age correlated proteins are more than mean + 3*sd or less than mean - 3*sd
#sd_threshold <- 3
sd_threshold <- sd_threshold
m <- mean(age_cor, na.rm = TRUE)
s <- sd(age_cor, na.rm = TRUE)
thr_low  <- m - sd_threshold*s
thr_high <- m + sd_threshold*s

highly_age_cor = names(age_cor[age_cor > thr_high | age_cor < thr_low])
message("Number of highly age correlated proteins: ", length(highly_age_cor))
message("Highly age correlated proteins: ", paste(highly_age_cor, collapse = ", "))

age_cor = as.data.frame(age_cor)
age_cor$Protein = rownames(age_cor)
out <- subset(
  age_cor,
  age_cor > thr_high | age_cor < thr_low
)

plot = ggplot(age_cor, aes(x = age_cor)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.02, 
                 fill = "steelblue", color = "black") +
  geom_density(color = "grey", size = 1) +
  # for two linesm add label
  geom_text(x = thr_high, y = 0, label = paste0("mean + ", sd_threshold, "*sd: ", round(thr_high,3)), vjust = -0.5, hjust = 1, color = "darkred") +
  geom_text(x = thr_low, y = 0, label = paste0("mean - ", sd_threshold, "*sd: ", round(thr_low,3)), vjust = -0.5, hjust = 0, color = "darkred") +
  geom_vline(xintercept = c(thr_high, thr_low), linetype = "dashed", color = "darkred") +
  geom_rug(alpha = 0.3) +
  geom_text_repel(
    data = subset(age_cor, age_cor > thr_high | age_cor < thr_low),
    aes(label = Protein, y = 0),
    nudge_y = 0.2,
    angle = 90,
    size = 3,
    max.overlaps = Inf
  ) +
  scale_x_continuous(breaks = seq(-1, 1, by = 0.05)) +   # every 0.05 on x-axis
  theme_classic() +
  labs(x = "Correlation with Age", y = "Density",
       title = "Distribution of Age Correlations")
ggsave(filename = file.path(output_dir, "age_correlation_distribution.pdf"), plot = plot, width = 8, height = 6)
ggsave(filename = file.path(output_dir, "age_correlation_distribution.png"), plot = plot, width = 8, height = 6)

# remove highly correlated proteins with age
inter = setdiff(inter,highly_age_cor)

# the expression data is log2 transformed, so we need to convert it back
data = as.data.frame(2^assay(Discovery))
colnames(data) = Discovery$label

data = data[inter,]

cluster_result = read.csv("Discovery/08_Clustering_als/cluster_assignments_2.csv")
colnames(cluster_result)[2] = "k2"
# if cluster_result$k2 is alpha, change it to 1
cluster_result$k2 = as.numeric(cluster_result$k2 == "alpha")

# Define response variable (binary outcome)
data = as.data.frame(t(data))
data$outcome = as.factor(cluster_result$k2)

# Define response variable (binary outcome)
y <- as.factor(data$outcome) 

# Define predictor variables (protein features)
X <- as.matrix(data[, -which(names(data) == "outcome")])

# normalization
sample_names = rownames(X)
X = as.matrix(scale_manual(as.data.frame(X)))
rownames(X) = sample_names

# Bootstrap LASSO Regression
# Number of bootstrap iterations
# n_boot <- 1000
# Store selected features from each bootstrap
selected_features <- matrix(0, ncol = ncol(X), nrow = n_boot)
colnames(selected_features) <- colnames(X)

cv_fit_result = list()
best_lambda_result = list()
lambda_result = list()
lambda_all = list()
feat_names <- colnames(X)
coef_names <- c("(Intercept)", feat_names)

# Model coefficients at the selected λ (min/1se) for each iteration
coef_mat <- matrix(NA_real_, nrow = n_boot, ncol = length(coef_names),
                   dimnames = list(paste0("run_", seq_len(n_boot)), coef_names))


# Perform bootstrap sampling
for (i in 1:n_boot) {
  # Bootstrap sample
  index <- sample(1:nrow(X), replace = TRUE)
  X_boot <- X[index, ]
  y_boot <- y[index]
  
  # LASSO with 10-fold cross-validation
  cv_fit <- cv.glmnet(X_boot, y_boot, alpha = 1, family = "binomial", nfolds = 10)
  cv_fit_result[[i]] = cv_fit
  
  # grep the lambda value
  lambda.min <- cv_fit$lambda.min
  lambda.1se <- cv_fit$lambda.1se
  
  #Rebuild the model using a specified λ value (to select genes based on the λ).
  model_lasso_min <- glmnet(X_boot, y_boot, alpha = 1, lambda = lambda.min,family = "binomial")
  model_lasso_1se <- glmnet(X_boot, y_boot, alpha = 1, lambda = lambda.1se,family = "binomial")
  
  if (lambda_value == "min") {
    lasso_fit <- model_lasso_min
    lambda_result[[i]] = lasso_fit
    best_lambda_result[[i]] = lambda.min
  } else if (lambda_value == "1se") {
    lasso_fit <- model_lasso_1se
    lambda_result[[i]] = lasso_fit
    best_lambda_result[[i]] = lambda.1se
  } else {
    stop("Invalid lambda value specified. Use 'min' or '1se'.")
  }
  
  lasso_fit_all <- glmnet(X_boot, y_boot, alpha = 1, family = "binomial")
  lambda_all[[i]] = lasso_fit_all

  # Store selected features (nonzero coefficients)
  selected_features[i, ] <- as.numeric(coef(lasso_fit)[-1] != 0)
  
  # ===== New: Save the coefficients at the selected lambda =====
  # coef(lasso_fit) returns a sparse matrix ((p+1) x 1); convert it to a vector and align by name
  beta_vec <- as.numeric(coef(lasso_fit))
  # the row names are already "(Intercept)" and feat_names; 
  # here we explicitly set the names to ensure proper alignment
  names(beta_vec) <- coef_names
  coef_mat[i, ] <- beta_vec
}

saveRDS(cv_fit_result, file = file.path(output_dir, "cv_fit_result.rds"))
saveRDS(best_lambda_result, file = file.path(output_dir, "best_lambda_result.rds"))
saveRDS(lambda_result, file = file.path(output_dir, "lambda_result.rds"))
saveRDS(lambda_all, file = file.path(output_dir, "lambda_all.rds"))
saveRDS(coef_mat, file = file.path(output_dir, "coef_mat.rds"))

# turn back to the original data (protein expression without log2 or min-max)
X <- as.matrix(data[, -which(names(data) == "outcome")])

# Compute selection frequency of each feature
feature_freq <- colMeans(selected_features)
selected_final <- names(feature_freq[feature_freq > feature_freq_cutoff])  # Features selected in >50% (feature_freq_cutoff) of bootstraps

# Print the most important proteins
print(selected_final)
write.table(selected_final, file = file.path(output_dir, "selected_features.txt"), row.names = FALSE, col.names = FALSE)


#validation_data <- as.data.frame(assay(Validation))
validation_data = as.data.frame(2^assay(Validation))
colnames(validation_data) <- Validation$label

validation_data = validation_data[inter,]
validation_data = as.data.frame(t(validation_data))

validation_cluster_result = read.csv("Validation/08_Clustering_als/cluster_assignments_2.csv")
colnames(validation_cluster_result)[2] = "k2"
# if cluster_result$k2 is alpha, change it to 1
validation_cluster_result$k2 = as.numeric(validation_cluster_result$k2 == "alpha")

# Define response variable (binary outcome)
validation_data$outcome <- as.factor(validation_cluster_result$k2)

# Subset validation data with selected features
X_valid <- validation_data[, selected_final]

X_selected <- X[, selected_final]
#The glm algorithm can be used in caret as follows:
discovery_selected = as.data.frame(X_selected)
discovery_selected$class = y
discovery_selected$class = ifelse(discovery_selected$class == 1, "alpha", "beta")

validation_selected = as.data.frame(X_valid)
validation_selected$class = validation_data$outcome
validation_selected$class = ifelse(validation_selected$class == 1, "alpha", "beta")

# Data Preprocessing: Ensure 'class' variable is a factor
discovery_selected$class <- as.factor(discovery_selected$class)
validation_selected$class <- as.factor(validation_selected$class)

# Ensure 'alpha' is set as the reference level (Caret requires the "positive" class to be first)
discovery_selected$class <- relevel(discovery_selected$class, ref = "alpha")

# split data
X_train <- discovery_selected[, -ncol(discovery_selected)]  # select all protein expression data
y_train <- discovery_selected$class  # classification
X_valid <- validation_selected[, -ncol(validation_selected)]
y_valid <- validation_selected$class

######Plot
# Generate Boxplots for Feature Expression in Discovery Cohort
plot_discovery <- log2(X_train)  # Log2 transform for better visualization
plot_discovery <- as.data.frame(plot_discovery)
plot_discovery$outcome <- y_train  # Assign class labels

# Convert to long format for ggplot
plot_long_discovery <- data.table::melt(plot_discovery, id.vars = "outcome", variable.name = "Protein", value.name = "Expression")

# Create the boxplot for Discovery Cohort
boxplot_discovery <- ggplot(plot_long_discovery, aes(x = as.factor(outcome), y = Expression, fill = as.factor(outcome))) +
  geom_boxplot(alpha = 1, outlier.shape = NA) +  
  geom_jitter(width = 0.2, alpha = 1, size = 1) +  
  facet_wrap(~Protein, scales = "free_y") +  
  stat_compare_means(aes(label = ..p.format..), method = "wilcox.test", label = "p.format") +  
  theme_classic() +
  labs(title = "Protein Expression in Discovery Cohort", 
       x = "ALS Subtype", y = "log2 Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = brewer.pal(8, "Set2")[c(2, 3)]) 

# Save Discovery Cohort boxplot
ggsave(filename = file.path(output_dir, "boxplot_discovery.pdf"), plot = boxplot_discovery, width = 12, height = 8)

# Generate Boxplots for Feature Expression in Validation Cohort
plot_validation <- log2(X_valid)  
plot_validation <- as.data.frame(plot_validation)
plot_validation$outcome <- y_valid  # Assign class labels

# Convert to long format for ggplot
plot_long_validation <- data.table::melt(plot_validation, id.vars = "outcome", variable.name = "Protein", value.name = "Expression")

# Create the boxplot for Validation Cohort
boxplot_validation <- ggplot(plot_long_validation, aes(x = as.factor(outcome), y = Expression, fill = as.factor(outcome))) +
  geom_boxplot(alpha = 1, outlier.shape = NA) +  
  geom_jitter(width = 0.2, alpha = 1, size = 1) +  
  facet_wrap(~Protein, scales = "free_y") +  
  stat_compare_means(aes(label = ..p.format..), method = "wilcox.test", label = "p.format") +  
  theme_classic() +
  labs(title = "Protein Expression in Validation Cohort", 
       x = "ALS Subtype", y = "log2 Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = brewer.pal(8, "Set2")[c(2, 3)]) 

# Save Validation Cohort boxplot
ggsave(filename = file.path(output_dir, "boxplot_validation.pdf"), plot = boxplot_validation, width = 12, height = 8)



# Plot Feature Selection Stability
feature_freq_df <- data.frame(Protein = names(feature_freq), Frequency = feature_freq) %>% filter(Frequency >= 0.1)

# Calculate standard error for error bars
feature_se <- apply(selected_features, 2, function(x) sd(x) / sqrt(length(x)))
feature_se_df <- data.frame(Protein = names(feature_se), SE = feature_se)
feature_freq_df <- dplyr::left_join(feature_freq_df, feature_se_df, by = "Protein")

write.csv(feature_freq_df, file = file.path(output_dir, "feature_selection_stability.csv"), row.names = FALSE)

feature_plot <- ggplot(feature_freq_df, aes(x = reorder(Protein, -Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = Frequency - SE, ymax = Frequency + SE), 
                width = 0.2, color = "black", size = 0.5) +
  geom_hline(yintercept = feature_freq_cutoff, linetype = "dashed", color = "red", size = 1) +  # Add dashed red line at feature_freq_cutoff
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Feature Selection Stability", x = "Protein", y = "Selection Frequency") +
  theme_classic()

# Save Feature Selection Stability plot
# base on the dim of the feature_freq_df, we can adjust the width and height
#ggsave(filename = file.path(output_dir, "feature_selection_stability.pdf"), plot = feature_plot, width = 10, height = 6)
ggsave(filename = file.path(output_dir, "feature_selection_stability.pdf"), plot = feature_plot, width = 10, height =  max(6, nrow(feature_freq_df) / 8))


#############Mean beta coefficients for selected variables in Lasso regression#############
coef_no_intercept <- coef_mat[, setdiff(colnames(coef_mat), "(Intercept)"), drop = FALSE]
# only keep proteins in selected_final
sel_in_mat <- intersect(selected_final, colnames(coef_no_intercept))
if (length(sel_in_mat) == 0) stop("None of selected_final found in coef_mat columns.")

# calculate mean for beta coefficients
beta_mean <- colMeans(coef_no_intercept[, sel_in_mat, drop = FALSE], na.rm = TRUE)

beta_df <- data.frame(
  Protein = names(beta_mean),
  beta    = as.numeric(beta_mean),
  stringsAsFactors = FALSE
)

DEP_results <- df[df$name %in% sel_in_mat, c("name", "alpha_vs_beta_diff", "fdr")]
DEP_results <- DEP_results[match(beta_df$Protein, DEP_results$name), ]

beta_df <- dplyr::left_join(beta_df, DEP_results, by = c("Protein" = "name")) %>%
  dplyr::select(Protein, beta, alpha_vs_beta_diff, fdr)

beta_df$log10_fdr <- -log10(beta_df$fdr)
beta_df$color <- ifelse(beta_df$beta > 0, "up", "down")

# sort beta_df by beta in descending order
beta_df$Protein <- factor(beta_df$Protein, levels = beta_df$Protein[order(beta_df$beta, decreasing = TRUE)])

# based on alpha_vs_beta_diff, set label colors
beta_df$label_color <- ifelse(beta_df$alpha_vs_beta_diff > 0, "#D73027", "#4575B4")
label_colors <- setNames(beta_df$label_color, as.character(beta_df$Protein))

p <- ggplot(beta_df, aes(x = Protein, y = beta,
                         size = log10_fdr,
                         fill = alpha_vs_beta_diff)) + geom_point(shape = 21, color = "black", stroke = 0.2) +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0,
                       name = expression(Log[2]~fold~change)) +
  scale_size_continuous(name = expression(-log[10]~adjusted~p~value)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "right"
  ) +
  ylab(expression(Mean~beta~coefficient)) +
  xlab("") +
  guides(fill = guide_colourbar(order = 1),
         size = guide_legend(order = 2))

# save plot to output directory
ggsave( file.path(output_dir, "mean_beta_coefficients.pdf"), plot = p, width = 8, height = 6, dpi = 300)
ggsave( file.path(output_dir, "mean_beta_coefficients.png"), plot = p, width = 8, height = 6, dpi = 300)
################mean ROC plot
results <- runML_with_plot(discovery_selected, validation_selected,
                           BS_number = n_boot, 
                           cv = 10, output_dir = output_dir)

# Save the results
saveRDS(results, file = file.path(output_dir, "ML_results.rds"))

# Variable importance plots
# svm_linear
svm_linear_imp <- plot_model_importance(results, model_key = "svm_linear", scale_importance = TRUE, output_dir = output_dir)
# svm_radial
svm_radial_imp <- plot_model_importance(results, model_key = "svm_radial", scale_importance = TRUE, output_dir = output_dir)
# ridge
ridge_imp <- plot_model_importance(results, model_key = "ridge", scale_importance = TRUE, output_dir = output_dir)
# lasso
lasso_imp <- plot_model_importance(results, model_key = "lasso", scale_importance = TRUE, output_dir = output_dir)
# logit
#logit_imp <- plot_model_importance(results, model_key = "logit", scale_importance = TRUE, output_dir = output_dir)
# rf
rf_imp <- plot_model_importance(results, model_key = "rf", scale_importance = TRUE, output_dir = output_dir)
# xgb
xgb_imp <- plot_model_importance(results, model_key = "xgb", scale_importance = TRUE, output_dir = output_dir)