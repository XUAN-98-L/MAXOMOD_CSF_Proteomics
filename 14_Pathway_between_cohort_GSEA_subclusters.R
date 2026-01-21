#=========================Script Description=================================
# This script is used for pathway scatterplot between discovery and validation
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggrepel"))
suppressMessages(library("tibble"))
#============================================================================
Discovery = readRDS("Discovery/06_GSEA_subclusters_k2/GSEA_result.rds")
Discovery = Discovery$BP@result

Validation = readRDS("Validation/06_GSEA_subclusters_k2/GSEA_result.rds")
Validation = Validation$BP@result

#Making signed p value plot
#MAKE PLOT WITH FDR
inter = intersect(Discovery$Description,
                  Validation$Description)
# sort by the order in inter
Discovery_inter = Discovery[match(inter, Discovery$Description),]
Validation_inter = Validation[match(inter, Validation$Description),]

library(data.table)


output_dir <- "14_Pathway_between_cohort_GSEA_subclusters"
if (!file.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}
#######################Make plot when p value cutoff 0.05#######
data = as.data.table(cbind(Discovery_inter$Description, 
                           Discovery_inter$p.adjust, 
                           Validation_inter$p.adjust))
# x is Discovery
# y is Validation
colnames(data) = c("name", "x", "y")
data$x = as.numeric(data$x)
data$y = as.numeric(data$y)
data$x = -log10(data$x)
data$y = -log10(data$y)

data$x[Discovery_inter$NES <0] = -data$x[Discovery_inter$NES <0]
data$y[Validation_inter$NES <0] = -data$y[Validation_inter$NES <0]

cut_off = -log10(0.05)
main_title = "alpha vs beta (GSEA p.adjust)"
max.overlaps = 10
lab_x = "signed -log10(p.adjust) for Discovery"
lab_y = "signed -log10(p.adjust) for Validation"
labels_T_F = T
annotate_YN = T
text_y = "significant in Validation"
text_x = "significant in Discovery"

data$omic_type = rep("ns", nrow(data))
data$omic_type[abs(data$y) >= cut_off] = text_y
data$omic_type[abs(data$x) >= cut_off] = text_x
data$omic_type[(abs(data$x) >= cut_off) & (abs(data$y) >= cut_off)] = "significant in both"
cols <- c("x" = "salmon", "y" = "#26b3ff", "ns" = "grey", "significant in both" = "mediumpurple1") 
attributes(cols)$names[1] = text_x
attributes(cols)$names[2] = text_y

p = ggplot(data, aes(x, y, text = name)) +  # Include 'text' in aes for interactive tooltip
  geom_point(aes(colour = omic_type),
             alpha = 0.5,
             shape = 16,
             size = 2) +
  geom_point(data = filter(data, abs(y) >= cut_off | abs(x) >= cut_off),
             aes(colour = omic_type), 
             alpha = 0.5, 
             shape = 16,
             size = 3) + 
  geom_hline(yintercept = cut_off, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -cut_off, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = cut_off, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = -cut_off, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey80") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey80") +
  scale_colour_manual(values = cols) + 
  scale_fill_manual(values = cols) +
  labs(title = main_title,
       x = lab_x,
       y = lab_y,
       colour = "Differential \nExpression") +
  theme_classic() + # Select theme with a white background  
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 15, hjust = 0.5),
        text = element_text(size = 14))+ 
  annotate("text", x = 3, y = -5, label = 
             paste0(sum(data$omic_type == text_y), " ", text_y, "\n", 
                    sum(data$omic_type == text_x), " ", text_x, "\n", 
                    sum(data$omic_type == "significant in both"), " significant in both"), size = 8 / .pt)+
  # annotate significant points and prevent label overlap
  geom_text_repel(data = filter(data, abs(y) >= cut_off | abs(x) >= cut_off), 
                  aes(label = name), 
                  size = 1.5, 
                  min.segment.length = 0, # ensure all labels have connection lines
                  direction = "both",    # allow labels to be adjusted around
                  segment.color = "grey50",  # connection line color
                  segment.size = 0.3,    # connection line width
                  max.overlaps = Inf,    # show all labels
                  box.padding = 0.3,     # text-point distance
                  point.padding = 0.2) 


ggsave(paste0(output_dir,"/scatterplot_FDR_005.pdf"), p, width = 10, height = 8, units = "in", dpi = 300)

sink(paste0(output_dir, "/FDR_cor.txt"))
cor.test(data$x, data$y, method = "pearson")
sink()