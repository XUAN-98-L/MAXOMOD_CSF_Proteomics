#=========================Script Description=================================
# This script is used for Clustering
# cd /Users/xliu2942/Documents/Projects/MAXOMOD/TESTing/Pipeline
# Rscript 5Clustering.R  >output.log
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("DEP"))
suppressMessages(library("cluster"))
suppressMessages(library("mclust"))
suppressMessages(library("reshape2"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggthemes"))
suppressMessages(library("ggalluvial"))
suppressMessages(library("tidyverse"))
suppressMessages(library("nnet"))
suppressMessages(library("stats"))
suppressMessages(library("dunn.test"))
suppressMessages(library("gridExtra"))
suppressMessages(library("ggpubr"))
suppressMessages(library("RColorBrewer"))
#===========================Function Definition=============================
#all colourblind friendly options
display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, 
                   colorblindFriendly=TRUE)

#create a list with all colours for the final paper figures
#each list element will be a named vector
final_colours = list(
  #red and blue for volcano plot
  volcano_plot = brewer.pal(n=10, "RdYlBu")[c(2,9)],
  male_female = brewer.pal(n=8, "Set2")[c(1,4)],
  #heatmap_scale = rev(brewer.pal(n = 11, "RdYlBu")),
  heatmap_scale = rev(brewer.pal(n=11, "RdBu")),
  clustering = brewer.pal(n=8, "Set2")[c(2,3,5)],
  age_scale = brewer.pal(n = 9, "YlGn"),
  disease_progression_scale = brewer.pal(n = 9, "YlOrRd"),
  onset = c(brewer.pal(n=11, "RdYlBu")[c(4,5)], "#B3B3B3"),
  disease_status = c(brewer.pal(n=11, "RdYlBu")[2],"#B3B3B3"),
  genetic_testing = c("#B3B3B3", brewer.pal(n = 11, "PRGn")[c(3, 9)]),
  center = c("purple4", "orange3"),
  Nfl = brewer.pal(n = 9, "PuBu"),
  pNFh_scale = brewer.pal(n = 9, "Purples"),
  age_at_onset_scale = brewer.pal(n = 9, "Blues"),
  slow_vital_capacity = brewer.pal(n = 9, "Reds")
)

#give the vectors with discrete variables names
names(final_colours$volcano_plot) =  c("up", "down")
names(final_colours$male_female) = c("Male", "Female")
names(final_colours$clustering) = c("alpha", "beta", "theta")
names(final_colours$disease_status) = c("als", "ctrl")
names(final_colours$onset) = c("spinal", "bulbar", "ctrl")
names(final_colours$genetic_testing) = c("not_performed", "negative", "C9orf72")
names(final_colours$center) = c("goettingen", "munich")


#function to calculate total within sum of squares
calc_SS = function(df) sum(as.matrix(dist(df)^2)) / (2 * nrow(df))
calc_TWSS = function(df, clusters){
  number_clusters = length(levels(as.factor(clusters)))
  sum_of_squares = c()
  for(i in 1:number_clusters){sum_of_squares[i] = calc_SS(df[clusters == i,])}
  return(sum(sum_of_squares))
}

#function to calculate BIC and AIC
BIC2 <- function(df, clusters){
  m = ncol(df)
  n = nrow(df)
  k = length(levels(as.factor(clusters)))
  D = calc_TWSS(df, clusters)
  return(data.frame(AIC = D + 2*m*k,
                    BIC = D + log(n)*m*k))
}

perform_clustering <- function(assay_data, seed, clin_labels) {
  set.seed(seed)
  
  silhouette_scores = TWSS_scores = AIC_scores = BIC_scores = as.data.frame(matrix()) #create empty dataframes for the different scores
  rownames(silhouette_scores) = rownames(TWSS_scores) = rownames(AIC_scores) = rownames(BIC_scores) = "hclust" #they only have one row and first row will be used for hierarchical clustering results      
  
  cluster_results <- list()
  cluster_assignments <- as.data.frame(clin_labels)
  colnames(cluster_assignments) <- "patid"
  
  cluster_methods <- c("hclust", "mclust", "kmeans", "pam")
  
  for (i in 2:10) {
    
    # Perform clustering
    #hierarchical clustering
    title = paste0("hclust_k=", i)
    dist_mat <- dist(assay_data, method = 'euclidean')
    cl <- hclust(dist_mat, method = 'ward.D')
    cluster_assignments[, title] <- cutree(cl, k = i)
    #cluster fit measures
    ss = silhouette(cluster_assignments[,title], dist(assay))
    silhouette_scores["hclust", i] = mean(ss[, 3])
    TWSS_scores["hclust", i] = calc_TWSS(assay, cluster_assignments[,title])
    AIC_scores["hclust", i] = BIC2(assay, cluster_assignments[,title])[1,1]
    BIC_scores["hclust", i] = BIC2(assay, cluster_assignments[,title])[1,2]
    
    #model based clustering
    title = paste0("mclust_k=", i)
    library(mclust)
    cl <- Mclust(assay_data, G = i)
    cluster_assignments[, title] <- cl$classification
    
    ss = silhouette(cluster_assignments[,title], dist(assay))
    silhouette_scores["mclust", i] = mean(ss[, 3])
    TWSS_scores["mclust", i] = calc_TWSS(assay, cluster_assignments[,title])
    AIC_scores["mclust", i] = BIC2(assay, cluster_assignments[,title])[1,1]
    BIC_scores["mclust", i] = BIC2(assay, cluster_assignments[,title])[1,2]   
    
    #kmeans clustering
    title = paste0("kmeans_k=", i)
    cl <- kmeans(assay_data, centers = i, nstart=25,iter.max = 50)
    cluster_assignments[, title] <- cl$cluster
    
    ss = silhouette(cluster_assignments[,title], dist(assay))
    silhouette_scores["kmeans", i] = mean(ss[, 3])
    TWSS_scores["kmeans", i] = calc_TWSS(assay, cluster_assignments[,title])
    AIC_scores["kmeans", i] = BIC2(assay, cluster_assignments[,title])[1,1]
    BIC_scores["kmeans", i] = BIC2(assay, cluster_assignments[,title])[1,2]
    
    #partitioning around medoids
    title = paste0("pam_k=", i)
    cl <- pam(assay_data, k = i)
    cluster_assignments[, title] <- cl$clustering
    
    ss = silhouette(cluster_assignments[,title], dist(assay))
    silhouette_scores["pam", i] = mean(ss[, 3])
    TWSS_scores["pam", i] = calc_TWSS(assay, cluster_assignments[,title])
    AIC_scores["pam", i] = BIC2(assay, cluster_assignments[,title])[1,1]
    BIC_scores["pam", i] = BIC2(assay, cluster_assignments[,title])[1,2]
  }
  
  cluster_assignments[, 2:ncol(cluster_assignments)] <- cluster_assignments[, 2:ncol(cluster_assignments)] - 1
  return(list(
    cluster_assignments = cluster_assignments,
    AIC_scores = AIC_scores,
    BIC_scores = BIC_scores,
    silhouette_scores = silhouette_scores,
    TWSS_scores = TWSS_scores
  ))
} 
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "2_Missing_Inspection_als/norm_imp_MinProb.rds",
              help = "SummarizedExperiment object after normalize and imputation."
  ),make_option(c("--output", "-o"),
                type = "character", default = "5_Clustering",
                help = "output directory path."
  ),make_option(c("--seed", "-e"),
                type = "integer", default = 9,
                help = "set.seed"
  ),make_option(c("--disease", "-d"), 
                type = "logical", default = TRUE,
                help = "Clustering base on als patients or ctrl group, default is TRUE (als patients)."
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
  se = readRDS(input)
}

if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
}else{
  seed = opt$seed
  set.seed(seed)
}

if (is.null(opt$disease)) {
  stop("Please provide the disease status!")
}else{
  disease = opt$disease
}

# Clustering data loading
clin = as.data.frame(se@colData@listData) #take the clinical data from the summarized experiment
#we can find the onset variable in the patient ID


if(disease){
  clin$onset = clin$condition
  clin$onset = as.factor(clin$onset)
  clin$sex = as.factor(clin$sex)
}else{
  clin$sex = clin$condition
  clin$sex = as.factor(clin$sex)
}


assay = as.data.frame(assay(se)) #make the abundancy matrix into a dataframe
assay = assay[,clin$ID] #order patients in assay in same order as clinical table
colnames(assay) = clin$label #give them CSF labels again
#columns should be variables so we need to transpose the table
assay = as.data.frame(t(assay))

#Perform different types of clustering
#use clustering function to get cluster labels
cluster_assignments_result <- perform_clustering(assay,seed = seed, clin_labels = clin$label)
cluster_assignments <- cluster_assignments_result$cluster_assignments
AIC_scores = cluster_assignments_result$AIC_scores
BIC_scores = cluster_assignments_result$BIC_scores
silhouette_scores = cluster_assignments_result$silhouette_scores
TWSS_scores = cluster_assignments_result$TWSS_scores
saveRDS(cluster_assignments, file.path(output_dir, "cluster_assignments.rds"))

# fit measures
scores = list(AIC = AIC_scores, BIC = BIC_scores, silhouette = silhouette_scores, TWSS = TWSS_scores) 
plots = list() #initiate empty plot list

for(i in 1:length(scores)){ #loop that iterates through all types of scores
  colnames(scores[[i]]) = 1:10
  scores[[i]]$method = rownames(scores[[i]])
  melt = reshape2::melt(scores[[i]])
  
  plots[[i]] = ggplot(melt, aes(x=variable, y=value, group = method, colour = method)) +
    geom_line() +
    geom_point() +
    ggtitle(names(scores)[i]) +
    theme_few()
  
}
names(plots) = names(scores)

allplots <- ggarrange(plotlist=plots,
                      labels = LETTERS[1:length(plots)],
                      ncol = 2, nrow = 2)

# Save the plots
ggsave(file.path(output_dir,"fit_scores.pdf"), width = 11, height = 8, units = "in")


# Sankey plots for all methods (Sankey_plots_within_method.pdf)
cluster_assignments <- cluster_assignments_result$cluster_assignments
number_of_clusters = as.character(2:5)
plots = list()

for(i in 1:length(number_of_clusters)){
  data = cluster_assignments[,c(1, grep(number_of_clusters[i], colnames(cluster_assignments)))]
  data = reshape2::melt(data)
  colnames(data) = c("patid", "method", "cluster")
  data = as.data.frame(lapply(data, as.factor))
  
  plots[[i]] = ggplot(data,
                      aes(x = method, stratum = cluster, alluvium = patid,
                          fill = cluster, label = cluster)) +
    scale_fill_brewer(type = "qual", palette = "Set3") +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
              color = "darkgray") +
    geom_stratum() +
    theme(legend.position = "bottom") +
    ggtitle(paste0("How patients are clustered in different clustering algorithms \nk=", number_of_clusters[i])) +
    theme_few()
}

# Sankey plots within a method, different number of clusters
methods = c("kmeans", "hclust", "mclust", "pam")
plots = list()
number_of_clusters

for(i in 1:length(methods)){
  data = cluster_assignments[,c("patid", paste0(methods[i],"_k=",number_of_clusters))]
  
  data = reshape2::melt(data, id = "patid")
  colnames(data) = c("patid", "method", "cluster")
  data = as.data.frame(lapply(data, as.factor))
  
  plots[[i]] = ggplot(data,
                      aes(x = method, stratum = cluster, alluvium = patid,
                          fill = cluster, label = cluster)) +
    scale_fill_brewer(type = "qual", palette = "Set3") +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
              color = "darkgray") +
    geom_stratum() +
    theme(legend.position = "bottom") +
    ggtitle(paste0("How patients are clustered within one clustering algorithm \nmethod = ", methods[i])) +
    theme_few() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

allplots <- ggarrange(plotlist=plots,
                      labels = 1:length(plots),
                      ncol = 2, nrow = 2)

# Save the plots
ggsave(file.path(output_dir, "Sankey_plots_within_method.pdf"), width = 11*2, height = 8*2, units = "in") 

## Sankey plots for kmeans k2 and k3
cluster_assignments_2 = cluster_assignments[,c("patid", "kmeans_k=2", "kmeans_k=3")]
#set the levels for the clusters
#levels k=2
# Create a table to count the objects in each cluster
cluster_counts <- table(cluster_assignments_2$"kmeans_k=2")
# Identify the larger and smaller cluster
larger_cluster <- names(cluster_counts)[which.max(cluster_counts)]
smaller_cluster <- names(cluster_counts)[which.min(cluster_counts)]
# Dynamically assign "alpha" to the larger cluster and "beta" to the smaller cluster
cluster_assignments_2$`kmeans_k=2` = as.factor(cluster_assignments_2$`kmeans_k=2`)
levels(cluster_assignments_2$`kmeans_k=2`) <- ifelse(levels(cluster_assignments_2$`kmeans_k=2`) == larger_cluster, "alpha", "beta")
cluster_assignments_2$`kmeans_k=2` = ordered(cluster_assignments_2$`kmeans_k=2`, levels = c("alpha", "beta", "theta"))

message("kmeans = 2:")
print(table(cluster_assignments_2$`kmeans_k=2`))

#levels k=3
# If “alpha” is the majority, assign “alpha”.
# If “beta” is the majority, assign “beta”.
# If there is a mix (or tie), assign “theta”.
adjusted_clusters <- cluster_assignments_2 %>%
  group_by(`kmeans_k=3`) %>%
  summarize(
    majority_alpha = sum(`kmeans_k=2` == "alpha"),
    majority_beta = sum(`kmeans_k=2` == "beta"),
    total = n()
  ) %>%
  mutate(
    alpha_ratio = majority_alpha / total,
    beta_ratio = majority_beta / total,
    alpha_beta = alpha_ratio - beta_ratio
  ) %>%
  # Assign based on alpha_beta values
  mutate(
    kmeans_k3_adjusted = case_when(
      alpha_beta == max(alpha_beta) ~ "alpha",  # Highest alpha_beta
      alpha_beta == min(alpha_beta) ~ "beta",   # Lowest alpha_beta
      TRUE ~ "theta"                            # Remaining one
    )
  )

adjusted_clusters = adjusted_clusters[,c("kmeans_k=3","kmeans_k3_adjusted")]
cluster_assignments_2 = left_join(cluster_assignments_2,adjusted_clusters, by = "kmeans_k=3") %>% dplyr::select(-`kmeans_k=3`) 
colnames(cluster_assignments_2)[3] = "kmeans_k=3"

cluster_assignments_2$`kmeans_k=3` = as.factor(cluster_assignments_2$`kmeans_k=3`)
cluster_assignments_2$`kmeans_k=3` <- ordered(cluster_assignments_2$`kmeans_k=3`, levels = c("alpha", "beta", "theta"))

write.csv(cluster_assignments_2, file = file.path(output_dir, "cluster_assignments_2.csv"), row.names = FALSE)

message("kmeans = 2:")
print(table(cluster_assignments_2$`kmeans_k=3`))

#Sankey plot with k=2 and k=3
data = pivot_longer(cluster_assignments_2, cols = 2:dim(cluster_assignments_2)[2])
colnames(data) = c("patid", "k", "cluster")

data23 = data[data$k == "kmeans_k=2" | data$k == "kmeans_k=3",]


ggplot(data23,
       aes(x = k, stratum = cluster, alluvium = patid,
           fill = cluster, label = cluster)) +
  scale_fill_manual(values = final_colours$clustering) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom") +
  ggtitle("Patient Sankey plot within K-means, from 2 clusters to 3 clusters \nPATIENTS") +
  theme_few()


ggsave(file.path(output_dir, "Sankey_plot_kmeans_2_3_patients.pdf"), width = 11/2, height = 8/2, units = "in")

# Filter the data frame to find the patient who changed from "alpha" to "beta"
changed_cluster <- cluster_assignments_2 %>%
  filter(`kmeans_k=2` == "alpha" & `kmeans_k=3` == "beta")
print(changed_cluster)

#CHECK ASSUMPTIONS: are the clinical variables normally distributed?
message("Checking assumptions: age")
pdf(file.path(output_dir, "age_histogram.pdf"))
hist(clin$age)
dev.off()
shapiro.test(clin$age) #YES

message("Checking assumptions: Nfl")
pdf(file.path(output_dir, "Nfl_histogram.pdf"))
hist(clin$Nfl)
dev.off()
shapiro.test(clin$Nfl) #NO

if(disease){
  message("Checking assumptions: age_at_onset")
  pdf(file.path(output_dir, "age_at_onset_histogram.pdf"))
  hist(clin$age_at_onset)
  dev.off()
  shapiro.test(clin$age_at_onset) #YES'
  
  message("Checking assumptions: progression_rate")
  pdf(file.path(output_dir, "progression_rate_histogram.pdf"))
  hist(clin$progression_rate)
  dev.off()
  shapiro.test(clin$progression_rate) #NO
  
  message("Checking assumptions: slow_vital_capacity")
  pdf(file.path(output_dir, "Slow_vital_capacity_histogram.pdf"))
  hist(clin$slow_vital_capacity)
  dev.off()
  shapiro.test(clin$slow_vital_capacity) #NO
  
  message("Checking assumptions: disease_duration")
  pdf(file.path(output_dir, "Disease_duration_histogram.pdf"))
  hist(clin$disease_duration)
  dev.off()
  shapiro.test(clin$disease_duration) #NO
  
  message("Checking assumptions: pNFh")
  pdf(file.path(output_dir, "pNFh_histogram.pdf"))
  hist(clin$pNFh)
  dev.off()
  shapiro.test(clin$pNFh) 
}else{
  message("This is not a disease dataset, so no need to check the assumptions for age_at_onset")
}


#Test difference clinical characteristics between clusters

#for each clustering algorithm, for each number of clusters, test if there are differences in the clinical variables
if(disease){
  clin_vars = c( "sex", "Nfl", "progression_rate",  "age_at_onset","age","onset","slow_vital_capacity","disease_duration","limb","pNFh","progression_group") # "age", "sex", "Nfl",, "progression_rate", "onset" 
}else{
  clin_vars = c( "sex", "Nfl",  "age") # "age", "sex", "Nfl",, "progression_rate", "onset" 
}

pval_clin = c()
models_clin = list()
l = 1

for(i in 2:ncol(cluster_assignments_2)){
  for(j in 1:length(clin_vars)){
    var = clin_vars[j]
    title = paste0(colnames(cluster_assignments_2)[i],"_",var)
    
    if(grepl("k=2", title, fixed = TRUE)){
      if(var %in%c("Nfl","progression_rate","disease_duration","pNFh")){
        #for k = 2
        model <- kruskal.test( clin[,var] ~ cluster_assignments_2[,i])
        d = as.data.frame(rbind(model[[2]], round(model[[3]],2)))
        rownames(d) = c("degrees of freedom", "p-value")
        colnames(d) = "Kruskal-wallis test"
        models_clin[[l]] = d
        pval_clin[l] = model[["p.value"]]
      }else{
        #for k = 2
        model <- glm(cluster_assignments_2[,i] ~ clin[,var], family = "binomial")
        models_clin[[l]] = summary(model)
        pval_clin[l] = models_clin[[l]]$coefficients[2,"Pr(>|z|)"]
      }
      
    }else{
      if(var %in%c("Nfl","progression_rate","disease_duration","pNFh")){
        #for k > 2
        model <- kruskal.test( clin[,var] ~ cluster_assignments_2[,i])
        d = as.data.frame(rbind(model[[2]], round(model[[3]],2)))
        rownames(d) = c("degrees of freedom", "p-value")
        colnames(d) = "Kruskal-wallis test"
        models_clin[[l]] = d
        pval_clin[l] = model[["p.value"]]
      }else{
        # for k > 2
        model <- multinom(cluster_assignments_2[,i] ~ clin[,var])
        models_clin[[l]] = summary(model)
        z <- summary(model)$coefficients/summary(model)$standard.errors
        p <- (1 - pnorm(abs(z), 0, 1)) * 2
        pval_clin[l] = min(p[,2])
      }
    }
    names(pval_clin)[l] = title
    names(models_clin)[l] = title
    l = l+1
  }
}

#for each clustering algorithm, for each number of clusters, test if there are differences in the clinical variables
colours = list(
  sex = c(Female = "lightpink", Male ="lightblue3"),
  Nfl = "royalblue",
  pNFh = "purple",
  progression_rate = "orange3",
  progression_group = c(
    "SP" = "#1B9E77",  
    "IP" = "#D95F02",   
    "FP" = "#7570B3"  
  ),
  age_at_onset = c("theta" = "#A6D854", "beta" = "#8DA0CB", "alpha" = "#FC8D62"),
  onset = c("spinal" = "#FDAE61", "bulbar" = "#FEE090", "ctrl" = "#B3B3B3"),
  slow_vital_capacity = "#A6D854",
  disease_duration = "#8DA0CB",
  limb = c("upper" = "#6A5ACD", "lower" = "#FF7F50","both" = "#A6D854")
)
plots = list()



l = 0
for(i in 2:ncol(cluster_assignments_2)){
  for(j in 1:length(clin_vars)){
    l = l+1
    var = clin_vars[j]
    title = paste0(colnames(cluster_assignments_2)[i],"_",var)
    
    #how many clusters?
    if(grepl("k=2", colnames(cluster_assignments_2)[i])){
      if(var %in%c("Nfl","progression_rate","disease_duration","pNFh")){
        coeff = models_clin[[title]]
        coeff = coeff[2,,drop = FALSE]}else{
          #with k=2
          coeff = round(models_clin[[title]]$coefficients,4)
          coeff = coeff[2,"Pr(>|z|)",drop = FALSE]
        }
    }else{
      if(var %in%c("Nfl","progression_rate","disease_duration","pNFh")){
        coeff = models_clin[[title]]
      }else{
        #with k>2
        #https://stats.stackexchange.com/questions/63222/getting-p-values-for-multinom-in-r-nnet-package
        z <- models_clin[[title]]$coefficients/models_clin[[title]]$standard.errors
        # 2-tailed Wald z tests to test significance of coefficients
        p <- (1 - pnorm(abs(z), 0, 1)) * 2
        coeff = as.data.frame(models_clin[[title]]$coefficients)
        coeff$pvalue = p
        coeff = round(coeff, 4)
        coeff = as.data.frame(as.matrix(coeff))
        colnames(coeff) = c("Intercept", var, "p-value intercept", paste0("p-value ", var))
        coeff = coeff[,paste0("p-value ", var),drop = FALSE]
      }
      
    }
    linesize = ifelse(pval_clin[title]<=0.05, 2, 0)
    
    #is it a factor?
    if(is.factor(clin[,var])){
      
      #factor variables
      df = data.frame(clin[,var], droplevels(cluster_assignments_2[,i]))
      df = na.omit(df)
      colnames(df) = c("variable", "cluster")
      df = as.data.frame(proportions(table(df), margin = 2))
      
      p = ggplot(df, aes(fill=variable, y=Freq, x=cluster)) + 
        geom_bar(position="fill", stat="identity") + 
        annotation_custom(tableGrob(coeff, theme= ttheme_minimal(base_size = 8)), 
                          #xmin = 1.8, xmax = 2, 
                          ymin = 0, ymax = 0.2) +
        scale_fill_manual(values=colours[[var]]) +
        labs(title = title,
             x = "cluster",
             y = paste0("proportions of ", var),
             fill = var) +
        theme_few() +
        theme(panel.background = element_rect( size=linesize))
      plots[[l]] = p
      names(plots)[l] = title
      
    }else{
      #numeric variables
      df = data.frame(clin[,var], as.factor(cluster_assignments_2[,i]))
      colnames(df) = c("variable", "cluster")
      df = na.omit(df)
      
      p = ggplot(df, aes(y=variable, x=cluster, fill=cluster)) + 
        # geom_boxplot(fill = colours[[var]]) +
        geom_violin(trim = FALSE, fill = "white")+
        labs(title = title,
             x = "cluster",
             y = var) +
        # guides(fill="none") +
        theme_few() +
        theme(panel.background = element_rect( size=linesize)) +
        annotation_custom(tableGrob(coeff, theme = ttheme_minimal(base_size = 8)), 
                          ymin = min(df$variable), 
                          ymax = min(df$variable)+10) +
        geom_dotplot(binaxis='y', stackdir='center', dotsize=1, color = "black", aes(fill = cluster))+scale_fill_manual( values = colours$age_at_onset)
      
      
      plots[[l]] = p
      names(plots)[l] = title
    }
  }}

# mean age at onset
if(disease){
  mean_age_at_onset = as.data.frame(tapply(clin$age_at_onset, cluster_assignments_2$`kmeans_k=3`, mean, na.rm = TRUE))
  colnames(mean_age_at_onset)[1] = "mean_age_at_onset"
  write.table(mean_age_at_onset, file.path(output_dir, "mean_age_at_onset.txt"), row.names = FALSE, col.names = FALSE,quote = FALSE)
}else{
  message("This is not a disease dataset, so no need to calculate the mean age at onset")
}

# Loop through each plot in the list and save it
for (plot_name in names(plots)) {
  # Create the file path using the plot name
  file_path <- file.path(paste0(output_dir,"/",plot_name, ".pdf"))
  
  # Save the plot using ggsave
  ggsave(filename = file_path, plot = plots[[plot_name]], width = 4, height = 4)
}


