###Final R Script to make figures


rm(list = ls())
gc()


####FIGURE 1------------
library(Seurat)
library(ggplot2)
library(ggrepel)



# Load the Seurat object
seurat_obj <- readRDS("C:/Users/ostrobel/Indiana University/O365-Strobel_Sequencing_Team - Documents/General/Analysis/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/2024_09_20_FINAL_SUERAT_OBJECT.rds")


seurat_obj_pseudo <- readRDS("C:/Users/ostrobel/Indiana University/O365-Strobel_Sequencing_Team - Documents/General/Analysis/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/2024_10_11_FINAL_ALL_pseudo.rds")



# Convert seurat_clusters to numeric if it is stored as a factor
numeric_clusters <- as.numeric(as.character(seurat_obj$seurat_clusters))

# Define the clusters for each general cell type
cluster_assignments <- list(
  "Immune" = c(2, 6, 7, 10, 11, 17, 19, 21, 24),
  "Stromal" = c(0, 1, 5, 12, 16, 20, 22, 23),
  "Seminal Vesicle" = c(13),
  "Epithelial" = c(3, 4, 8, 9, 14, 15, 18)
)

# Create a vector to store the cell type for each cluster
cell_type_general <- rep(NA, length(unique(numeric_clusters)))

# Assign the general cell type based on clusters and print assignments
for (cell_type in names(cluster_assignments)) {
  clusters <- cluster_assignments[[cell_type]]
  cat("Assigning", cell_type, "to clusters:", clusters, "\n")
  cell_type_general[clusters + 1] <- cell_type  # Adjust for 0-based indexing
}

# Assign this to a new metadata column in your Seurat object
seurat_obj$Cell_Type_General <- factor(cell_type_general[numeric_clusters + 1])
View(seurat_obj@meta.data)


saveRDS(seurat_obj, file = "C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/2024_09_20_FINAL_SUERAT_OBJECT.rds")


#####DimPlot overall assignments
# Load RColorBrewer
library(RColorBrewer)

# Adjust the DimPlot to use the Dark2 palette
fig1_dimplot <- DimPlot(seurat_obj, group.by = "Cell_Type_General", reduction = "umap.harmonyintegration", label = TRUE, pt.size = 0.01) + 
  ggtitle("Cell Type General") +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Dark2")  # Use Dark2 color palette



fig_1_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure1"


# Define the file path to save the plot
file_name <- file.path(fig_1_path, "DimPlot_Cell_Type_General_Dark2.png")

# Save the plot to the specified folder
ggsave(filename = file_name, plot = fig1_dimplot, width = 4, height = 4, dpi = 300)

# Output confirmation
cat("DimPlot saved to:", file_name, "\n")


FeaturePlot(seurat_obj, features = "Epcam", reduction = "umap.harmonyintegration")



#####DimPlot unsupervized
fig1_dimplot <- DimPlot(seurat_obj, group.by = "seurat_clusters", reduction = "umap.harmonyintegration", label = TRUE, pt.size = 0.01) + 
  ggtitle("Cell Type General") +
  theme(legend.position = "none")


# Define the file path to save the plot
file_name <- file.path(fig_1_path, "DimPlot_unsupervized.png")

# Save the plot to the specified folder
ggsave(filename = file_name, plot = fig1_dimplot, width = 4, height = 4, dpi = 300)

# Output confirmation
cat("DimPlot saved to:", file_name, "\n")


###Condition Plot
# Assuming the condition information is stored in the "condition" metadata column
cell_counts <- table(seurat_obj$condition)

# Print the number of cells for each condition
cat("Number of cells in each condition:\n")
print(cell_counts)
#control 41174
#experimental 41525





##Downsampling to try the condition plot 
library(RColorBrewer)
# Define the group of interest (e.g., 'condition' column should contain 'Control' and 'Experimental')
group_column <- "condition"

# Specify the number of cells you want to sample from each group
n_control <- 7000  # Adjust this to your desired sample size
n_experimental <- 7054  # Adjust this to your desired sample size

# Identify cells from each group
control_cells <- WhichCells(seurat_obj, expression = condition == "Control")
experimental_cells <- WhichCells(seurat_obj, expression = condition == "Experimental")

# Randomly sample cells from each group based on the specified sample size
set.seed(123)  # For reproducibility
control_sample <- sample(control_cells, min(n_control, length(control_cells)))
experimental_sample <- sample(experimental_cells, min(n_experimental, length(experimental_cells)))

# Combine the downsampled cells
balanced_cells <- c(control_sample, experimental_sample)

# Subset the Seurat object to include only the downsampled cells
seurat_obj_downsampled <- subset(seurat_obj, cells = balanced_cells)

# Define custom colors for control and experimental groups
custom_colors <- c("Control" = "navyblue", "Experimental" = "orange")

# Plot the UMAP with the downsampled data, split by condition and custom colors
balanced_umap <- DimPlot(seurat_obj_downsampled, group.by = "condition", reduction = "umap.harmonyintegration", cols = custom_colors, pt.size = 0.1) +
  ggtitle("Manually Adjusted UMAP - Control vs Experimental")

# Display the plot
print(balanced_umap)

# Save the plot
ggsave(filename = file.path(fig_1_path, "Manually_Adjusted_UMAP_Condition.png"), plot = balanced_umap, width = 5, height = 4, dpi = 300, device = "png")


####Dotplot
library(ggplot2)

# Define your genes in the desired order for the x-axis
genes_to_plot <- c("Vim", "Dcn", "Col6a1", "Acta2",
                   "Svs2", "Svs5", "Svs4", "Svs1",   
                   "Ptprc", "Cd3d", "Cd3g", "Cd4",   
                   "Epcam", "Krt8", "Krt18", "Krt19") 

# Create the DotPlot with the specified gene order
dotplot <- DotPlot(seurat_obj, features = genes_to_plot, group.by = "Cell_Type_General", dot.scale = 8) +
  scale_colour_gradient(low = "grey", high = "darkblue") +  # Apply color gradient
  ggtitle("Gene Expression by Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),   # Rotate x-axis labels
        axis.text.y = element_text(size = 10))               # Ensure y-axis labels show as expected


# Display the plot
print(dotplot)

# Save the plot
ggsave(filename = file.path(fig_1_path, "CellType_GeneExpression_DotPlot35.png"), plot = dotplot, width = 10, height = 3.5)


###4x3 FeaturePlot
library(Seurat)
library(ggplot2)
library(patchwork)

# Define the genes to plot
genes_to_plot <- c("Svs2", "Svs5", "Svs4", "Svs1",   
                   "Ptprc", "Cd3d", "Cd3g", "Cd4",   
                   "Vim", "Dcn", "Acta2", "Col6a1",
                   "Epcam", "Krt8", "Krt18", "Krt19")


# Create the individual FeaturePlots without x and y axis labels
plots <- lapply(genes_to_plot, function(gene) {
  FeaturePlot(seurat_obj, features = gene) + 
    ggtitle(gene) +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank())
})

# Combine the plots into a 4x3 grid to accommodate all 12 genes
combined_plot <- wrap_plots(plots, ncol = 3, nrow = 4)

# Display the combined plot
print(combined_plot)


# Save the plot
ggsave(filename = file.path(fig_1_path, "4_3FeaturePlot.png"), plot = combined_plot, width = 12, height = 16)




####bar plot
library(dplyr)
library(ggplot2)


View(seurat_obj@meta.data)

# Assuming you have a metadata column "Cell_Type_General" for cell types and "orig.ident" for sample identities

# Step 1: Calculate percentage of each cell type per sample
cell_type_percentages <- seurat_obj@meta.data %>%
  group_by(orig.ident, Cell_Type_General) %>%
  summarise(count = n()) %>%
  mutate(percent = count / sum(count) * 100)

# Step 2: Plot the 100% stacked bar plot
ggplot(cell_type_percentages, aes(x = orig.ident, y = percent, fill = Cell_Type_General)) +
  geom_bar(stat = "identity", position = "fill") +  # position = "fill" makes it a 100% stacked bar plot
  labs(x = "Sample", y = "Percentage", title = "Cell Type Composition per Sample") +
  scale_y_continuous(labels = scales::percent_format()) +  # Show y axis as percentages
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x axis labels if needed

ggsave(filename = file.path(fig_1_path, "CellType_Composition_BarPlot.png"), width = 4, height = 4, dpi = 300)

#try2
library(dplyr)
library(ggplot2)

# Define the order of cell types
cell_type_order <- c("Seminal Vesicle", "Immune","Epithelial","Stromal" )

# Step 1: Set Cell_Type_General as a factor with the desired order
seurat_obj@meta.data$Cell_Type_General <- factor(seurat_obj@meta.data$Cell_Type_General, levels = cell_type_order)

# Step 2: Calculate percentage of each cell type per sample
cell_type_percentages <- seurat_obj@meta.data %>%
  group_by(orig.ident, Cell_Type_General) %>%
  summarise(count = n()) %>%
  mutate(percent = count / sum(count) * 100)

# Step 3: Plot the 100% stacked bar plot with the desired order
ggplot(cell_type_percentages, aes(x = orig.ident, y = percent, fill = Cell_Type_General)) +
  geom_bar(stat = "identity", position = "fill") +  # position = "fill" makes it a 100% stacked bar plot
  labs(x = "Sample", y = "Percentage", title = "Cell Type Composition per Sample") +
  scale_y_continuous(labels = scales::percent_format()) +  # Show y axis as percentages
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x axis labels if needed
  scale_fill_manual(values = c("Stromal" = "#e7298a",  # Light blue
                               "Immune" = "#d95f02",      # Light green
                               "Epithelial" = "#1b9e77",  # Light orange
                               "Seminal Vesicle" = "#7570b3"))  # Light pink



# Step 4: Save the plot
ggsave(filename = file.path(fig_1_path, "CellType_Composition_BarPlot_try345.png"), width = 4, height = 4, dpi = 300)

View(seurat_obj@meta.data)

###making the plot with raw numbers so that i can combine the two later
library(tidyverse)

# Create the data frame with raw counts
cell_counts <- tribble(
  ~orig.ident, ~`Seminal Vesicle`, ~Immune, ~Epithelial, ~Stromal,
  "C3", 403, 3059, 3383, 6393,
  "C4", 526, 2386, 5420, 5238,
  "C5", 301, 2639, 5600, 5826,
  "E3", 129, 4119, 2278, 6401,
  "E4", 117, 4997, 2239, 4400,
  "E5", 120, 5660, 3321, 7744
)

cell_long$Cell_Type_General <- factor(cell_long$Cell_Type_General,
                                      levels = c("Seminal Vesicle", "Immune", "Epithelial", "Stromal"))


# Reorder the factor levels to control stack order (bottom to top)
cell_long$Cell_Type_General <- factor(cell_long$Cell_Type_General,
                                      levels = c("Seminal Vesicle", "Immune", "Epithelial", "Stromal"))

# Plot
p<-ggplot(cell_long, aes(x = orig.ident, y = percent, fill = Cell_Type_General)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Sample", y = "Percentage", title = "Cell Type Composition per Sample") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(
    "Stromal" = "#e7298a",
    "Immune" = "#d95f02",
    "Epithelial" = "#1b9e77",
    "Seminal Vesicle" = "#7570b3"
  ))

ggsave(
  filename = "Figure1_cell_type_composition.svg",
  plot = p,  # or omit this if it's the last plot
  path = "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure1",
  width = 4, height = 3, dpi = 300
)


# Save as SVG
ggsave(
  filename = file.path(fig_1_path, "CellType_Composition_Plot_minus_Immune.svg"),
  plot = p,
  width = 4, height = 4, dpi = 300
)




library(tibble)
library(readr)

# ---- Create tibble ----
cell_counts <- tribble(
  ~orig.ident, ~`Seminal Vesicle`, ~Immune, ~Epithelial, ~Stromal,
  "C3", 403, 3059, 3383, 6393,
  "C4", 526, 2386, 5420, 5238,
  "C5", 301, 2639, 5600, 5826,
  "E3", 129, 4119, 2278, 6401,
  "E4", 117, 4997, 2239, 4400,
  "E5", 120, 5660, 3321, 7744
)

# ---- Save as CSV ----
write_csv(
  cell_counts,
  "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure2/components/cell_counts.csv"
)



# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# Define the percentages for each cell type in each sample
data <- data.frame(
  Cell_Type_General = c("Stromal", "Epithelial", "Seminal Vesicle"),
  C3 = c(0.628057766, 0.332350919, 0.039591315),
  C4 = c(0.468347639, 0.484620887, 0.047031474),
  C5 = c(0.496802251, 0.477530485, 0.025667264),
  E3 = c(0.726725704, 0.25862852, 0.014645777),
  E4 = c(0.651272943, 0.331409118, 0.01731794),
  E5 = c(0.692355834, 0.296915512, 0.010728654)
)

# Reshape the data into a long format for ggplot2
data_long <- tidyr::gather(data, key = "Sample", value = "Percentage", -Cell_Type_General)

# Reorder Cell_Type_General so that Seminal Vesicle is first in the order (appearing at the top of the bars)
data_long$Cell_Type_General <- factor(data_long$Cell_Type_General, levels = c("Seminal Vesicle", "Epithelial", "Stromal"))

# Filter out only the Stromal cells for label display
stromal_data <- data_long %>% filter(Cell_Type_General == "Stromal")

# Plot the percentages as a 100% stacked bar plot
plotctcps1 <- ggplot(data_long, aes(x = Sample, y = Percentage, fill = Cell_Type_General)) +
  geom_bar(stat = "identity", position = "fill") +  # position = "fill" makes it a 100% stacked bar plot
  labs(x = "Sample", y = "Percentage", title = "Cell Type Composition per Sample") +
  scale_y_continuous(labels = scales::percent_format()) +  # Show y axis as percentages
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x axis labels if needed
  scale_fill_manual(values = c("Stromal" = "#e7298a",  # Dark pink
                               "Epithelial" = "#1b9e77",  # Dark green
                               "Seminal Vesicle" = "#7570b3")) +  # Light purple
  # Add labels for the Stromal cell percentages
  geom_text(data = stromal_data, aes(label = scales::percent(Percentage, accuracy = 0.1)), 
            position = position_fill(vjust = 0.3), angle=65,color = "white", size = 3.5)

# Save as SVG to hardcoded path
ggsave(
  filename = "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure1/CellType_Composition_Plot_minus_Immune.svg",
  plot = plotctcps1,
  width = 4, height = 3, dpi = 300
)


library(Seurat)
library(dplyr)

seurat_obj <- PrepSCTFindMarkers(seurat_obj, verbose = TRUE)

all_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top10 <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

htmp DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()


library(Seurat)
library(dplyr)

# make sure cluster identities are set
Idents(seurat_obj) <- "seurat_clusters"

# 1. Find top markers
all_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top10 <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# 2. Collapse to cluster-average expression
avg_exp_10 <- AverageExpression(seurat_obj, features = top10$gene, return.seurat = TRUE)

# 3. Heatmap of cluster means
htmp_10 <- DoHeatmap(
            avg_exp_10,
            features = top10$gene,
          )

htmp_10


htmp1 <- DoHeatmap(seurat_obj, 
                   features = top5$gene,
                   cells = sample(colnames(seurat_obj), 10000)
                   ) + NoLegend()
htmp1

top25

# 5. output folder
out_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure1"

# save as PNG
png(file.path(out_dir, "heatmap1.png"), width = 2400, height = 2100, res = 300)
print(htmp1)
dev.off()



####FIGURE 2------------
fib_subset_path <- "C:\\Users\\ostrobel\\Indiana University\\Jerde, Travis J - JERDE LAB (J Drive)\\Strobel\\Single Cell Analysis Stuff\\R_Analysis_Strobel\\Combined_Sample_Analysis\\Combined_analysis\\combined_obj_try4\\Harmony_SCT\\Fib_subset\\2024_10_07_fib_subset_FINAL.rds"


fib_subset <- readRDS(fib_subset_path)



####overall clusters
#View(fib_subset@meta.data)
fig_2_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure2"


# Create a DimPlot for the fib_subset Seurat object
fib_dimplot <- DimPlot(fib_subset, group.by = "seurat_clusters", reduction = "umap", label = TRUE, pt.size = 0.01)

# Display the plot
print(fib_dimplot)

# Save the DimPlot as a PNG file
ggsave(filename = file.path(fig_2_path, "Fib_Subset_DimPlot.png"), plot = fib_dimplot, width = 4, height = 4, dpi = 300)


###Condition Plot
# Assuming the condition information is stored in the "condition" metadata column
cell_counts <- table(seurat_obj$condition)

# Print the number of cells for each condition
cat("Number of cells in each condition:\n")
print(cell_counts)
#control 41174
#experimental 41525





##Downsampling to try the condition plot 
library(RColorBrewer)
# Define the group of interest (e.g., 'condition' column should contain 'Control' and 'Experimental')
group_column <- "condition"

# Specify the number of cells you want to sample from each group
n_control <- 7000  # Adjust this to your desired sample size
n_experimental <- 7054  # Adjust this to your desired sample size

# Identify cells from each group
control_cells <- WhichCells(seurat_obj, expression = condition == "Control")
experimental_cells <- WhichCells(seurat_obj, expression = condition == "Experimental")

# Randomly sample cells from each group based on the specified sample size
set.seed(123)  # For reproducibility
control_sample <- sample(control_cells, min(n_control, length(control_cells)))
experimental_sample <- sample(experimental_cells, min(n_experimental, length(experimental_cells)))

# Combine the downsampled cells
balanced_cells <- c(control_sample, experimental_sample)

# Subset the Seurat object to include only the downsampled cells
seurat_obj_downsampled <- subset(seurat_obj, cells = balanced_cells)

# Define custom colors for control and experimental groups
custom_colors <- c("Control" = "navyblue", "Experimental" = "orange")

# Plot the UMAP with the downsampled data, split by condition and custom colors
balanced_umap <- DimPlot(seurat_obj_downsampled, group.by = "condition", reduction = "umap.harmonyintegration", cols = custom_colors, pt.size = 0.1) +
  ggtitle("Manually Adjusted UMAP - Control vs Experimental")

# Display the plot
print(balanced_umap)

# Save the plot
ggsave(filename = file.path(fig_1_path, "Manually_Adjusted_UMAP_Condition.png"), plot = balanced_umap, width = 5, height = 4, dpi = 300, device = "png")






###percent composition change. Pulled values from excel 
# Load necessary libraries
library(ggplot2)

# Create a dataframe with your data
sample_data <- data.frame(
  Fib_cluster = c(3, 4, 0, 5, 6, 2, 1, 7),
  pct_control_avg = c(
    mean(c(0.015915119, 0.010282546, 0.00750405)),
    mean(c(0.014539739, 0.010550787, 0.015946107)),
    mean(c(0.451223106, 0.316434192, 0.355248572)),
    mean(c(0.004027901, 0.003755365, 0.007845144)),
    mean(c(0.001866588, 0.000983548, 0.001620193)),
    mean(c(0.023086747, 0.018776824, 0.017907393)),
    mean(c(0.115924944, 0.106759657, 0.088684233)),
    mean(c(0.001473622, 0.000804721, 0.002046559))
  ),
  pct_experimental_avg = c(
    mean(c(0.058015441, 0.01731794, 0.020742065)),
    mean(c(0.015440509, 0.02294257, 0.026642825)),
    mean(c(0.546889192, 0.461219657, 0.486276263)),
    mean(c(0.005449591, 0.004144464, 0.009566384)),
    mean(c(0.000681199, 0.002220249, 0.002324542)),
    mean(c(0.022593097, 0.024718769, 0.020384443)),
    mean(c(0.07708901, 0.117821196, 0.124988824)),
    mean(c(0.000567666, 0.000888099, 0.001430487))
  )
)

# Calculate fold change as log2 ratio (log2(Experimental/Control))
sample_data$log2_fold_change <- log2(sample_data$pct_experimental_avg / sample_data$pct_control_avg)
# 
# # Create a dataframe with individual sample data (for the dots)
# sample_data <- data.frame(
#   Fib_cluster = rep(df$Fib_cluster, each = 3),
#   pct_control = c(0.015915119, 0.010282546, 0.00750405,
#                   0.014539739, 0.010550787, 0.015946107,
#                   0.451223106, 0.316434192, 0.355248572,
#                   0.004027901, 0.003755365, 0.007845144,
#                   0.001866588, 0.000983548, 0.001620193,
#                   0.023086747, 0.018776824, 0.017907393,
#                   0.115924944, 0.106759657, 0.088684233,
#                   0.001473622, 0.000804721, 0.002046559),
#   pct_experimental = c(0.058015441, 0.01731794, 0.020742065,
#                        0.015440509, 0.02294257, 0.026642825,
#                        0.546889192, 0.461219657, 0.486276263,
#                        0.005449591, 0.004144464, 0.009566384,
#                        0.000681199, 0.002220249, 0.002324542,
#                        0.022593097, 0.024718769, 0.020384443,
#                        0.07708901, 0.117821196, 0.124988824,
#                        0.000567666, 0.000888099, 0.001430487)
# )
# 
# # Calculate the log2 fold change for individual points
# sample_data$log2_fold_change <- log2(sample_data$pct_experimental / sample_data$pct_control)

# Create the bar plot with points
fold_change_plot <- ggplot(df, aes(x = factor(Fib_cluster), y = log2_fold_change)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  
  labs(title = "Log2 Fold Change per Fibroblast Cluster", 
       x = "Fibroblast Cluster", y = "Log2 Fold Change (Experimental / Control)") +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")  # Line at 0 to indicate no change

# Display the plot
print(fold_change_plot)

# Save the plot to the Figure 2 folder
ggsave(file.path(fig_2_path, "Log2_Fold_Change_Per_Fibroblast_Cluster.png"), plot = fold_change_plot, width = 4, height = 4)


# Create the bar plot with color gradient
fold_change_plot <- ggplot(df, aes(x = factor(Fib_cluster), y = log2_fold_change, fill = log2_fold_change)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_gradient2(low = "red", mid = "skyblue", high = "darkblue", midpoint = 0) +  # Color scale for gradient
  labs(title = "Log2 Fold Change per Fibroblast Cluster", 
       x = "Fibroblast Cluster", y = "Log2 Fold Change (Experimental / Control)") +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")  # Line at 0 to indicate no change

# Display the plot
print(fold_change_plot)

  
# Save the plot to the Figure 2 folder
ggsave(file.path(fig_2_path, "Log2_Fold_Change_Per_Fibroblast_Clusterv3.png"), plot = fold_change_plot, width = 4, height = 4)

  
  
  
  
  
  
  
  
  
  

###trying another version of making this plot
# Assuming your data is structured in a dataframe like this:
# Control and experimental percentage values along with cluster IDs.

# Create the dataframe from your image data
data <- data.frame(
  Fib_cluster = c(3, 4, 0, 5, 6, 2, 1, 7),
  pct.control3 = c(0.015915119, 0.014539739, 0.451223106, 0.004027901, 0.001866588, 0.023086747, 0.115924944, 0.001473622),
  pct.control4 = c(0.010282546, 0.010550787, 0.316434192, 0.003755365, 0.000983548, 0.018776824, 0.106759657, 0.000804721),
  pct.control5 = c(0.00750405, 0.015946107, 0.355248572, 0.007845144, 0.001620193, 0.017907393, 0.088684233, 0.002046559),
  pct.experimental3 = c(0.058015441, 0.015440509, 0.546889192, 0.005449591, 0.000681199, 0.022593097, 0.07708901, 0.000567666),
  pct.experimental4 = c(0.01731794, 0.02294257, 0.461219657, 0.004144464, 0.002220249, 0.024718769, 0.117821196, 0.000888099),
  pct.experimental5 = c(0.020742065, 0.026642825, 0.486276263, 0.009566384, 0.002324542, 0.020384443, 0.124988824, 0.001430487)
)

# Reshape the data to long format using tidyr
library(tidyr)
library(ggplot2)

long_data <- data %>%
  pivot_longer(cols = starts_with("pct"),
               names_to = c("condition", "animal"),
               names_pattern = "pct.(.*)([345])",
               values_to = "percent") %>%
  mutate(condition = ifelse(condition == "control", "Control", "Experimental"),
         Fib_cluster = as.factor(Fib_cluster),
         animal = factor(animal, labels = c("Animal 3", "Animal 4", "Animal 5")))

# Plot the data, adding different shapes for animals and colors for conditions
fold_change_plot <- ggplot(long_data, aes(x = Fib_cluster, y = percent, color = condition, shape = animal)) +
  geom_point(size = 3) + # Adjust point size
  scale_color_manual(values = c("Control" = "orange", "Experimental" = "navyblue")) + # Set colors
  theme_minimal() +
  ggtitle("Fold Change per Fibroblast Cluster by Animal") +
  ylab("Percent") +
  xlab("Fibroblast Cluster") +
  theme(legend.position = "right") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") # Optional dashed line at 0

# Display the plot
print(fold_change_plot)


# Save the plot to your Figure 2 folder
ggsave(file.path(fig_2_path, "Fold_Change_Per_Cluster_Per_Animal.png"), plot = fold_change_plot, width = 6, height = 4)


###Final way to make this stupid plot. 
# Data provided
data <- data.frame(
  Fib_cluster = as.factor(c(3, 4, 0, 5, 6, 2, 1, 7)),
  pct_diff = c(2.914706324, 1.608445349, 1.345089231, 1.293606704, 1.186563526, 1.121848887, 1.041025399, 0.69584377)
)

# Assign colors based on increase or decrease
data$change <- ifelse(data$pct_diff > 1, "increase", "decrease")

# Create the bar plot with log2 scale for y-axis
library(ggplot2)

fold_change_plot <- ggplot(data, aes(x = Fib_cluster, y = log2(pct_diff), fill = change)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("increase" = "blue", "decrease" = "red")) +
  ylab("Log2 Fold Change (Experimental/Control)") +
  xlab("Fibroblast Cluster") +
  ggtitle("Log2 Fold Change per Fibroblast Cluster") +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")  # Line at y = 0 (log2(1) = 0) for no change

# Save the plot to the Figure2 folder
ggsave(file.path(fig_2_path, "Log2_Fold_Change_Per_Cluster.png"), plot = fold_change_plot, width = 4, height = 4)

# Display the plot
print(fold_change_plot)


####i need to actually make this so that the colors match. 
# Farben definieren (8 Cluster, 8 Farben)
set1_colors <- c("#8DA0CB", "#4682B4", "#66C2A5", "#FFD700", 
                 "#FF69B4", "#8A2BE2", "#3c3c3c", "#E41A1C")

# Data
data <- data.frame(
  Fib_cluster = as.factor(c(3, 4, 0, 5, 6, 2, 1, 7)),
  pct_diff = c(2.914706324, 1.608445349, 1.345089231, 
               1.293606704, 1.186563526, 1.121848887, 
               1.041025399, 0.69584377)
)

library(ggplot2)

fold_change_plot <- ggplot(data, aes(x = Fib_cluster, 
                                     y = log2(pct_diff), 
                                     fill = Fib_cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = set1_colors) +   # nur set1_colors
  ylab("Log2 Fold Change (Experimental/Control)") +
  xlab("Fibroblast Cluster") +
  ggtitle("Log2 Fold Change per Fibroblast Cluster") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

# Speichern
ggsave(file.path(fig_2_path, "Log2_Fold_Change_Per_Clusterv19.svg"), 
       plot = fold_change_plot, width = 4, height = 4)

print(fold_change_plot)

print(fig_2_path)


###assigning cluster identities. 
# Define the cluster identities
cluster_identities <- c(
  "Prostate Fibroblast",    # Cluster 0
  "Urethral/Ductal Fibroblast",  # Cluster 1
  "Cxcl13+ Fibroblast",     # Cluster 2
  "Immune Reactive",                 # Cluster 3
  "Myofibroblast",          # Cluster 4
  "Endothelial",            # Cluster 5
  "Schwann Cells",          # Cluster 6
  "Smoc2+ Fibroblast"       # Cluster 7
)

# Convert seurat_clusters to numeric and subtract 1 for correct indexing
numeric_clusters <- as.numeric(as.character(fib_subset$seurat_clusters))

# Assign the identities to a new metadata column based on the numeric clusters
fib_subset$Fib_Cluster_Identities <- factor(cluster_identities[numeric_clusters + 1])

# View the first few rows to check the assignment
head(fib_subset@meta.data)

###trying to label the legend
# Define the cluster names in the format you want
cluster_names3 <- c(
  "0 - Prostate Fibroblast",    # Cluster 0
  "1 - Urethral/Ductal Fibroblast",  # Cluster 1
  "2 - Cxcl13+ Fibroblast",     # Cluster 2
  "3 - Immune Reactive",                 # Cluster 3
  "4 - Myofibroblast",          # Cluster 4
  "5 - Endothelial",            # Cluster 5
  "6 - Schwann Cells",          # Cluster 6
  "7 - Smoc2+ Fibroblast"       # Cluster 7
)

# Assign the cluster names to a new metadata column
fib_subset$cluster_labels3 <- factor(fib_subset$seurat_clusters, labels = cluster_names3)

View(fib_subset@meta.data)
# 
# set1_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#F781BF", "#A65628", "#FFFF33")
# 
# set1_colors <- viridis::viridis(8)
# 
# set1_colors <- c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC")
# set1_colors <- c("#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", ", "#3288BD")
# set1_colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3") ###descent
# set1_colors <- cividis::cividis(8)
# set1_colors <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5")
set1_colors <- c("#8DA0CB", "#4682B4", "#66C2A5", "#FFD700", "#FF69B4", "#8A2BE2","#3c3c3c", "#E41A1C") ###using this going forward

#4D4D4D (Dark Grey)
#5C5C5C (Slightly Lighter Dark Grey)
#3C3C3C (Darker Grey)
#2E2E2E (Very Dark Grey)
#6E6E6E (Medium Dark Grey)
#707070 (Another Medium Dark Grey)
#323232 (Almost Black Grey)
#282828 (Very Deep Grey)

# Create a DimPlot for the fib_subset Seurat object using the new cluster labels
fib_dimplot <- DimPlot(fib_subset, group.by = "cluster_labels3", reduction = "umap", label = FALSE) +
  # scale_color_manual(values = set1_colors) +  # Apply Set1 hex colors, excluding the first
  guides(color = guide_legend(title = "Clusters", override.aes = list(size = 4))) +  # Increase legend dot size
  theme(legend.position = "right")  # Position legend on the right side

# Display the plot
print(fib_dimplot)

# Save the DimPlot as a PNG file
ggsave(filename = file.path(fig_2_path, "Fib_Subset_DimPlot_Set1_Colorsv56433345.png"), plot = fib_dimplot, width = 6, height = 4, dpi = 300)



# Create a DimPlot for the fib_subset Seurat object using the new cluster labels
fib_dimplot <- DimPlot(fib_subset, group.by = "SCT_snn_res.0.1", reduction = "umap", label = TRUE) +
  theme(legend.position = "none") 

# Display the plot
print(fib_dimplot)

# Save the DimPlot as a PNG file
ggsave(filename = file.path(fig_2_path, "Fib_Subset_DimPlot_Set1_Colorsv56433344_nolegend_origcolor.png"), plot = fib_dimplot, width = 4, height = 4, dpi = 300)


# Create a DimPlot for the fib_subset Seurat object using SCT_snn_res.0.1 as the label
fib_dimplot <- DimPlot(fib_subset, group.by = "SCT_snn_res.0.1", reduction = "umap", label = TRUE, 
                       label.size = 4, repel = FALSE) +
  scale_color_manual(values = set1_colors) +  # Use custom colors
  theme(legend.position = "none")  # Remove legend
# Save the DimPlot as a PNG file
ggsave(filename = file.path(fig_2_path, "Fib_Subset_DimPlot_Set1_Colorsv564333441_nolegend_.png"), plot = fib_dimplot, width = 4, height = 4, dpi = 300)


# Display the plot
print(fib_dimplot)



###Dot plot for Fib Genes


genes_to_plot <- c(
  "Vim", "Dcn", "Fbn1", # Overall
  "Cxcl13", "Ifi207", "Sparcl1", # Cxcl13+ Fibroblast (Cluster 2)
  "Pecam1", "Ccl21a", "Cldn5", # Endothelial (Cluster 5)
  "Ccl5", "Cd3g", "Ptprc", # Immune (Cluster 3)
  "C3", "Ebf1", "Igf1", # Prostate Fibroblast (Cluster 0)
  "Plp1", "Cadm2", "Mbp", # Schwann Cells (Cluster 6)
  "Smoc2", "Cpe", "Rasgrp2",  # Smoc2+ Fibroblast (Cluster 7)
  "Acta2", "Tagln", "Plcb1", # Smooth Muscle (Cluster 4)
  "Lgr5", "Mfap4", "Ctla2a" # Urethral/Ductal Fibroblast (Cluster 1)
)





# Define genes for each cluster
genes_to_plot <- c(
  "Vim", "Dcn", "Fbn1", # Overall
  "C3", "Ebf1", "Igf1", # Prostate Fibroblast (Cluster 0)
  "Lgr5", "Mfap4", "Ctla2a", # Urethral/Ductal Fibroblast (Cluster 1)
  "Cxcl13", "Ifi207", "Sparcl1", # Cxcl13+ Fibroblast (Cluster 2)
  "Ccl5", "Cd52", "Hcst", # Immune (Cluster 3)
  "Acta2", "Tagln", "Plcb1", # Smooth Muscle (Cluster 4)
  "Pecam1", "Ccl21a", "Cldn5", # Endothelial (Cluster 5)
  "Plp1", "Cadm2", "Mbp", # Schwann Cells (Cluster 6)
  "Smoc2", "Cpe", "Rasgrp2"  # Smoc2+ Fibroblast (Cluster 7)
)

# Create DotPlot with customized dark blue color
dot_plot <- DotPlot(fib_subset, features = genes_to_plot, group.by = "cluster_labels3", dot.scale = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotate x-axis labels
  ggtitle("Dot Plot of Selected Genes Across Fibroblast Cluster Identities") +
  theme(legend.position = "right") +  # Adjust legend position if necessary
  scale_color_gradient(low = "grey", high = "darkblue")  # Customize color gradient

# Display the plot
print(dot_plot)

# Reverse the factor levels for cluster_labels3 to have Cluster 0 at the top
fib_subset$cluster_labels3 <- factor(fib_subset$cluster_labels3, levels = rev(levels(fib_subset$cluster_labels3)))



# Define genes for each cluster
genes_to_plot <- c(
  "Vim", "Dcn", "Fbn1", # Overall
  "Fn1", "C3", "Igf1", # Prostate Fibroblast (Cluster 0)
  "Ctla2a", "Lgr5", "Mfap4",  # Urethral/Ductal Fibroblast (Cluster 1)
  "Cxcl13", "Ifi207", "Sparcl1", # Cxcl13+ Fibroblast (Cluster 2)
  "Ccl5", "Cd52", "Hcst", # Immune (Cluster 3)
  "Acta2", "Tagln", "Plcb1", # Smooth Muscle (Cluster 4)
  "Pecam1", "Ccl21a", "Cldn5", # Endothelial (Cluster 5)
  "Plp1", "Cadm2", "Mbp", # Schwann Cells (Cluster 6)
  "Smoc2", "Cpe", "Rasgrp2"  # Smoc2+ Fibroblast (Cluster 7)
)

# Create DotPlot with customized dark blue color
dot_plot <- DotPlot(fib_subset, features = genes_to_plot, group.by = "cluster_labels3", dot.scale = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotate x-axis labels
  ggtitle("Dot Plot of Selected Genes Across Fibroblast Cluster Identities") +
  theme(legend.position = "right") +  # Adjust legend position if necessary
  scale_color_gradient(low = "grey", high = "darkblue")  # Customize color gradient

# Display the plot
print(dot_plot)


# Save the DotPlot
ggsave(filename = file.path(fig_2_path, "DotPlot_Cluster_Genesv2.png"), plot = dot_plot, width = 12, height = 4, dpi = 300)


library(ggplot2)  # Make sure ggplot2 is loaded
# 
# # Create DotPlot, using Fib_Cluster_Identities for cluster names
# dot_plot <- DotPlot(fib_subset, features = genes_to_plot, group.by = "Fib_Cluster_Identities", dot.scale = 8) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotate x-axis labels
#   ggtitle("Dot Plot of Selected Genes Across Fibroblast Cluster Identities") +
#   
#   # Customize the color scale: blue for downregulated, orange for upregulated
#   scale_color_gradientn(colors = c("blue", "darkgrey", "orange")) +
#   
#   theme(legend.position = "right")  # Adjust legend position if necessary
# 
# # Display the plot
# print(dot_plot)
# 
# # Save the DotPlot
# ggsave(filename = file.path(fig_2_path, "DotPlot_Cluster_Genesv3.png"), plot = dot_plot, width = 12, height = 4, dpi = 300)








###Condition Plot
# Assuming the condition information is stored in the "condition" metadata column
cell_counts <- table(fib_subset$condition)

# Print the number of cells for each condition
cat("Number of cells in each condition:\n")
print(cell_counts)
#control 17457        
#experimental 18545 





# Define the group of interest (e.g., 'condition' column should contain 'Control' and 'Experimental')
group_column <- "condition"

# Specify the number of cells you want to sample from each group
n_control <- 7000  # Adjust this to your desired sample size
n_experimental <- 7440  # Adjust this to your desired sample size

# Identify cells from each group
control_cells <- WhichCells(fib_subset, expression = condition == "Control")
experimental_cells <- WhichCells(fib_subset, expression = condition == "Experimental")

# Randomly sample cells from each group based on the specified sample size
set.seed(123)  # For reproducibility
control_sample <- sample(control_cells, min(n_control, length(control_cells)))
experimental_sample <- sample(experimental_cells, min(n_experimental, length(experimental_cells)))

# Combine the downsampled cells
balanced_cells <- c(control_sample, experimental_sample)

# Subset the Seurat object to include only the downsampled cells
fib_subset_downsampled <- subset(fib_subset, cells = balanced_cells)

# Define custom colors for control and experimental groups
custom_colors <- c("Control" = "navyblue", "Experimental" = "orange")

# Plot the UMAP with the downsampled data, split by condition and custom colors
balanced_umap <- DimPlot(fib_subset_downsampled, group.by = "condition", reduction = "umap", cols = custom_colors, pt.size = 0.1) +
  ggtitle("Cell Distribution - Condition") +  theme(legend.position = "right") 

# Display the plot
print(balanced_umap)

# Save the plot
#ggsave(filename = file.path(fig_2_path, "Manually_Adjusted_UMAP_Condition_legend.png"), plot = balanced_umap, width = 4, height = 4, dpi = 300, device = "png")
ggsave(filename = file.path(fig_2_path, "Manually_Adjusted_UMAP_Condition_legend.svg"), plot = balanced_umap, width = 4, height = 2.6, dpi = 300, device = "svg")







# Plot the UMAP with the downsampled data, split by condition and custom colors
balanced_umap <- DimPlot(fib_subset_downsampled, group.by = "condition", reduction = "umap", cols = custom_colors, pt.size = 0.01) +
  ggtitle("Cell Distribution - Condition") +
  theme(legend.position = "none")  # Remove legend

# Display the plot
print(balanced_umap)


ggsave(filename = file.path(fig_2_path, "Manually_Adjusted_UMAP_Condition_nolegend.svg"), plot = balanced_umap, width = 2.6, height = 2.6, dpi = 300, device = "svg")



###making plot to show composition based on cell type. This is normalized only to fibroblasts. 


library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# Create the data frame
data <- data.frame(
  Fib_cluster = c(3, 4, 0, 5, 6, 2, 1, 7),
  pct.C1 = c(0.015915119, 0.014539739, 0.451223106, 0.004027901, 0.001866588, 0.023086747, 0.115924944, 0.001473622),
  pct.C2 = c(0.010282546, 0.010550787, 0.316434192, 0.003755365, 0.000983548, 0.018776824, 0.106759657, 0.000804721),
  pct.C3 = c(0.00750405, 0.015946107, 0.355248572, 0.007845144, 0.001620193, 0.017907393, 0.088684233, 0.002046559),
  pct.E1 = c(0.058015441, 0.015440509, 0.546889192, 0.005449591, 0.000681199, 0.022593097, 0.07708901, 0.000567666),
  pct.E2 = c(0.01731794, 0.02294257, 0.461219657, 0.004144464, 0.002220249, 0.024718769, 0.117821196, 0.000888099),
  pct.E3 = c(0.020742065, 0.026642825, 0.486276263, 0.009566384, 0.002324542, 0.020384443, 0.124988824, 0.001430487)
)

# Transform data to long format for ggplot
data_long <- data %>%
  pivot_longer(cols = starts_with("pct"), 
               names_to = "Sample", 
               values_to = "Percentage") %>%
  mutate(Sample = gsub("pct.", "", Sample),  # Clean up Sample names
         Fib_cluster = factor(Fib_cluster, levels = rev(sort(unique(Fib_cluster)))))  # Order clusters 0 to 7 bottom to top

# Define custom colors and labels for each Fib_cluster
custom_colors <- c(
  "0" = "#8da0cb",  # Prostate Fibroblast
  "1" = "#4682b4",  # Urethral/Ductal Fibroblast
  "2" = "#66c2a5",  # Cxcl13+ Fibroblast
  "3" = "#ffd700",  # Immune Reactive
  "4" = "#ff69b4",  # Myofibroblast
  "5" = "#8a2be2",  # Endothelial
  "6" = "#3c3c3c",  # Schwann Cells
  "7" = "#e41a1c"   # Smoc2+ Fibroblast
)

# Custom labels for Fib_cluster
cluster_labels <- c(
  "0" = "0 - Prostate Fibroblast",
  "1" = "1 - Urethral/Ductal Fibroblast",
  "2" = "2 - Cxcl13+ Fibroblast",
  "3" = "3 - Immune Reactive",
  "4" = "4 - Myofibroblast",
  "5" = "5 - Endothelial",
  "6" = "6 - Schwann Cells",
  "7" = "7 - Smoc2+ Fibroblast"
)

# Plot
plot1 <- ggplot(data_long, aes(x = Sample, y = Percentage, fill = Fib_cluster)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked 100% bar plot
  labs(x = "Sample", y = "Percentage", title = "Stromal Cell Type Composition per Sample", fill = "Fib Cluster") +
  scale_y_continuous(labels = percent_format()) +  # Show y axis as percentages
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x axis labels if needed
  scale_fill_manual(values = custom_colors, labels = cluster_labels)  # Use custom colors and labels

# Display the plot
print(plot1)

# Save the plot
output_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure2/CellType_Stromal_Composition_BarPlot.png"
ggsave(filename = output_path, plot = plot1, width = 4.5, height = 2.6, dpi = 300)



# Convert Fib_cluster to a factor with Cluster 0 at the bottom
data_long$Fib_cluster <- factor(data_long$Fib_cluster, levels = c("0", "1", "2", "3", "4", "5", "6", "7"))

# Plot
plot1 <- ggplot(data_long, aes(x = Sample, y = Percentage, fill = Fib_cluster)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked 100% bar plot
  labs(x = "Sample", y = "Percentage", title = "Stromal Cell Type Composition per Sample", fill = "Fib Cluster") +
  scale_y_continuous(labels = scales::percent_format()) +  # Show y axis as percentages
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),  # Rotate and set x-axis labels to black
    axis.text.y = element_text(color = "black"),  # Set y-axis labels to black
    axis.title.x = element_text(color = "black"),  # Set x-axis title to black
    axis.title.y = element_text(color = "black")   # Set y-axis title to black
  ) +
  scale_fill_manual(values = custom_colors, labels = cluster_labels)  # Use custom colors and labels

# Display the plot
print(plot1)


output_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure2/CellType_Stromal_Composition_BarPlotv3.png"
ggsave(filename = output_path, plot = plot1, width = 4.5, height = 2.6, dpi = 300)



# Reverse the factor levels for Fib_cluster to have 0 at the top of the legend
data_long <- data_long %>%
  mutate(Fib_cluster = factor(Fib_cluster, levels = rev(sort(unique(Fib_cluster)))))  # Reverse order of clusters for legend

# Plot
plot1 <- ggplot(data_long, aes(x = Sample, y = Percentage, fill = Fib_cluster)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked 100% bar plot
  labs(x = "Sample", y = "Percentage", title = "Stromal Cell Type Composition per Sample", fill = "Fib Cluster") +
  scale_y_continuous(labels = scales::percent_format()) +  # Show y axis as percentages
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),  # Rotate and set x-axis labels to black
    axis.text.y = element_text(color = "black"),  # Set y-axis labels to black
    axis.title.x = element_text(color = "black"),  # Set x-axis title to black
    axis.title.y = element_text(color = "black")   # Set y-axis title to black
  ) +
  scale_fill_manual(values = custom_colors, labels = cluster_labels)  # Use custom colors and labels

# Display the plot
print(plot1)

# Save the plot
output_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure2/CellType_Composition_BarPlot_v3.svg"
ggsave(filename = output_path, plot = plot1, width = 4.5, height = 2.6, dpi = 300, device = "svg")




####Supps for figure 2--------
library(Seurat)
library(ggplot2)
library(patchwork)


fib_subset_path <- "C:\\Users\\ostrobel\\Indiana University\\Jerde, Travis J - JERDE LAB (J Drive)\\Strobel\\Single Cell Analysis Stuff\\R_Analysis_Strobel\\Combined_Sample_Analysis\\Combined_analysis\\combined_obj_try4\\Harmony_SCT\\Fib_subset\\2024_10_07_fib_subset_FINAL.rds"


fib_subset <- readRDS(fib_subset_path)

#loading in my object 
seurat_obj <- readRDS("C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/2024_09_20_FINAL_SUERAT_OBJECT.rds")


genes_of_interest <- c("Vim")

save_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Supplementals/Supps for Figure 2"



# Loop through each gene in the list
for (gene_of_interest in genes_of_interest) {
  
  # 1) FeaturePlot for the full Seurat object (top)
  plot1 <- FeaturePlot(
    object    = seurat_obj, 
    features  = gene_of_interest, 
    reduction = "umap.harmonyintegration"
  ) + 
    ggtitle(paste0(gene_of_interest, " - All Cells"))
  
  # 2) FeaturePlot for the fibroblast subset (bottom)
  plot2 <- FeaturePlot(
    object    = fib_subset, 
    features  = gene_of_interest, 
    reduction = "umap"
  ) + 
    ggtitle(paste0(gene_of_interest, " - Stromal Subset"))
  
  # 3) Stack the two plots vertically (ncol = 1)
  combined_plot <- wrap_plots(plot1, plot2, ncol = 1)
  
  # 4) Save the combined plot
  file_name <- file.path(save_path, paste0(gene_of_interest, "_plot.png"))
  ggsave(filename = file_name, plot = combined_plot, width = 4.5, height = 9)
  
  # 5) Confirmation message
  cat("Saved plot for", gene_of_interest, "to", file_name, "\n")
}
View(seurat_obj@meta.data)


library(Seurat)
library(ggplot2)
library(patchwork)

# Define your genes
row1_genes <- c("Vim", "Dcn", "Fbn1")
row2_genes <- c("C3", "Ebf1", "Igf1")

# A helper function to create the top/bottom (vertical) layout for one gene
make_combined_plot <- function(gene_of_interest) {
  plot_top <- FeaturePlot(
    object    = seurat_obj, 
    features  = gene_of_interest, 
    reduction = "umap.harmonyintegration",
    cols      = c("grey", "darkblue"),
  ) + ggtitle(paste0(gene_of_interest, " - All Cells"))
  
  plot_bottom <- FeaturePlot(
    object    = fib_subset,
    features  = gene_of_interest,
    reduction = "umap",
    cols      = c("grey", "darkblue"),
  ) + ggtitle(paste0(gene_of_interest, " - Stromal Subset"))
  
  # Stack top and bottom vertically
  wrap_plots(plot_top, plot_bottom, ncol = 1)
}

# Create a combined plot (top+bottom) for each gene
p1 <- make_combined_plot("Vim")
p2 <- make_combined_plot("Dcn")
p3 <- make_combined_plot("Fbn1")
p4 <- make_combined_plot("C3")
p5 <- make_combined_plot("Ebf1")
p6 <- make_combined_plot("Igf1")

# Arrange these 6 vertical combos (p1-p6) in 2 rows x 3 columns
final_plot <- wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 3)
#final_plot
# Save to file (adjust dimensions as needed)
ggsave(
  filename = file.path(save_path, "Combined_FeaturePlots_2rows_3columns.png"),
  plot     = final_plot,
  width    = 13.5, 
  height   = 18
)




"Vim", "Dcn", "Fbn1", # Overall
"C3", "Ebf1", "Igf1", # Prostate Fibroblast (Cluster 0)
"Lgr5", "Mfap4", "Ctla2a", # Urethral/Ductal Fibroblast (Cluster 1)
"Cxcl13", "Ifi207", "Sparcl1", # Cxcl13+ Fibroblast (Cluster 2)
"Ccl5", "Cd52", "Hcst" # Immune (Cluster 3)
"Acta2", "Tagln", "Plcb1", # Smooth Muscle (Cluster 4)
"Pecam1", "Ccl21a", "Cldn5", # Endothelial (Cluster 5)
"Plp1", "Cadm2", "Mbp", # Schwann Cells (Cluster 6)
"Smoc2", "Cpe", "Rasgrp2"  # Smoc2+ Fibroblast (Cluster 7)



library(Seurat)
library(ggplot2)
library(patchwork)

# 1. Define a helper function for a single gene
#    Top: full Seurat object, Bottom: fib_subset
make_combined_plot <- function(gene_of_interest) {
  top_plot <- FeaturePlot(
    object    = seurat_obj, 
    features  = gene_of_interest, 
    reduction = "umap.harmonyintegration",
    cols      = c("grey", "darkblue")
  ) + ggtitle(paste0(gene_of_interest, " - All Cells"))
  
  bottom_plot <- FeaturePlot(
    object    = fib_subset,
    features  = gene_of_interest,
    reduction = "umap",
    cols      = c("grey", "darkblue")
  ) + ggtitle(paste0(gene_of_interest, " - Stromal Subset"))
  
  # Stack top and bottom vertically
  wrap_plots(top_plot, bottom_plot, ncol = 1)
}

# 2. Define the 9 sets of genes (3 genes per set)
#    and name them based on the comment after the # in your list
gene_sets <- list(
  "Overall"       = c("Vim", "Dcn", "Fbn1"),
  "Cluster0"      = c("C3", "Ebf1", "Igf1"),
  "Cluster1"      = c("Lgr5", "Mfap4", "Ctla2a"),
  "Cluster2"      = c("Cxcl13", "Ifi207", "Sparcl1"),
  "Cluster3"      = c("Ccl5", "Cd52", "Hcst"),
  "Cluster4"      = c("Acta2", "Tagln", "Plcb1"),
  "Cluster5"      = c("Pecam1", "Ccl21a", "Cldn5"),
  "Cluster6"      = c("Plp1", "Cadm2", "Mbp"),
  "Cluster7"      = c("Smoc2", "Cpe", "Rasgrp2")
)

# 3. Loop over each gene set to create a single-row panel (3 genes wide)
for(set_name in names(gene_sets)) {
  
  # Extract the 3 genes in this set
  genes_current <- gene_sets[[set_name]]
  
  # Build a vertical stack for each gene
  plot_list <- lapply(genes_current, make_combined_plot)
  
  # Combine the three vertical stacks in one row
  # => 1 row x 3 columns
  final_plot <- wrap_plots(plot_list, ncol = 3)
  
  # 4. Save each combined figure with a distinct filename
  #    Height = 9 (as requested); adjust width as desired
  ggsave(
    filename = file.path(save_path, paste0(set_name, ".svg")),
    plot     = final_plot,
    device   = "svg",
    width    = 13.5,   # or adjust to your preference
    height   = 9
  )
  
  cat("Saved", set_name, "plot as", paste0(set_name, ".svg"), "\n")
}



library(Seurat)
library(ggplot2)
library(patchwork)

# ---- Save path ----
save_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/S6"

# 1. Define a helper function for a single gene
make_combined_plot <- function(gene_of_interest) {
  top_plot <- FeaturePlot(
    object    = seurat_obj, 
    features  = gene_of_interest, 
    reduction = "umap.harmonyintegration",
    cols      = c("grey", "darkblue")
  ) + ggtitle(paste0(gene_of_interest, " - All Cells"))
  
  bottom_plot <- FeaturePlot(
    object    = fib_subset,
    features  = gene_of_interest,
    reduction = "umap",
    cols      = c("grey", "darkblue")
  ) + ggtitle(paste0(gene_of_interest, " - Stromal Subset"))
  
  wrap_plots(top_plot, bottom_plot, ncol = 1)
}

# 2. Define only Cluster 3 genes
genes_cluster3 <- c("Ccl5", "Cd52", "Hcst")

# 3. Build plots for each Cluster 3 gene
plot_list <- lapply(genes_cluster3, make_combined_plot)

# 4. Combine into one row (3 genes wide)
final_plot <- wrap_plots(plot_list, ncol = 3)

# 5. Save output as SVG
ggsave(
  filename = file.path(save_path, "Cluster3.svg"),
  plot     = final_plot,
  device   = "svg",
  width    = 13.5,
  height   = 9
)

# 6. Save as PNG
ggsave(
  filename = file.path(save_path, "Cluster3.png"),
  plot     = final_plot,
  device   = "png",
  width    = 13.5,
  height   = 9,
  dpi      = 300
)

cat("Saved Cluster 3 plot as Cluster3.svg in", save_path, "\n")






gahhh

genes_of_interest <- c("Cxcl13")


# Define the folder where you want to save the plots
save_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Cxcl13_plots_by_sample"

# Get unique orig.ident values
ident_values <- unique(seurat_obj$orig.ident)

# Loop through each gene
for (gene_of_interest in genes_of_interest) {
  
  # Loop through each orig.ident
  for (ident_value in ident_values) {
    
    # Subset the main Seurat object by this orig.ident
    seurat_obj_subset <- subset(
      x = seurat_obj,
      subset = orig.ident == ident_value
    )
    
    # Subset the fibroblast object by the same orig.ident
    fib_subset_subset <- subset(
      x = fib_subset,
      subset = orig.ident == ident_value
    )
    
    # 1) FeaturePlot for the (subsetted) full object
    plot1 <- FeaturePlot(
      object    = seurat_obj_subset, 
      features  = gene_of_interest, 
      reduction = "umap.harmonyintegration"
    ) + 
      ggtitle(
        paste0(gene_of_interest, " - All Cells - ", ident_value)
      )
    
    # 2) FeaturePlot for the (subsetted) fib subset
    plot2 <- FeaturePlot(
      object    = fib_subset_subset, 
      features  = gene_of_interest, 
      reduction = "umap"
    ) + 
      ggtitle(
        paste0(gene_of_interest, " - Stromal Subset - ", ident_value)
      )
    
    # 3) Combine the two plots vertically
    combined_plot <- wrap_plots(plot1, plot2, ncol = 1)
    
    # 4) Construct a filename with gene + orig.ident
    file_name <- file.path(
      save_path,
      paste0(gene_of_interest, "_", ident_value, "_plot.png")
    )
    
    # Save the combined plot
    ggsave(
      filename = file_name,
      plot     = combined_plot,
      width    = 4.5,
      height   = 9
    )
    
    # 5) Message for confirmation
    cat("Saved plot for gene:", gene_of_interest,
        "and orig.ident:", ident_value, 
        "->", file_name, "\n")
  }
}





####Figure3 DEG analysis-----------
# Load necessary libraries
library(dplyr)
library(DESeq2)
library(fgsea)
library(ggplot2)
library(ggrepel)
library(openxlsx)

fib_subset_pseudo_path <- "C:\\Users\\ostrobel\\Indiana University\\Jerde, Travis J - JERDE LAB (J Drive)\\Strobel\\Single Cell Analysis Stuff\\R_Analysis_Strobel\\Combined_Sample_Analysis\\Combined_analysis\\combined_obj_try4\\Harmony_SCT\\Fib_subset\\2024_10_07_fib_sub_psudo_FINAL.rds"

fib_subset.pseudo <- readRDS(fib_subset_pseudo_path)

# Path to the GMT file for GSEA
gmt_file <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/DE_excel_docs_for_analysis/Mouse_GOBP_AllPathways_noPFOCR_no_GO_iea_October_01_2024_symbol.gmt"

# Load the pathways using the GMT file
pathways <- gmtPathways(gmt_file)

# Define the cluster to run the analysis (Cluster 0)
clusters_to_run <- 0

# Specify the output directory
output_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure4"

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through only the specified cluster(s)
for (cluster in clusters_to_run) {
  cat("Processing cluster:", cluster, "\n")
  
  # Create directory for the current cluster
  cluster_dir <- file.path(output_dir, paste0("Cluster", cluster))
  if (!dir.exists(cluster_dir)) {
    dir.create(cluster_dir, recursive = TRUE)
  }
  
  # Subset the pseudobulk object for the current cluster
  cluster_data <- subset(fib_subset.pseudo, idents = cluster)
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = GetAssayData(cluster_data, layer = "counts"),
    colData = cluster_data@meta.data,
    design = ~ condition  # Compare conditions within this cluster
  )
  
  # Filter genes with fewer than 10 total counts in at least 3 samples
  smallestGroupSize <- 3
  keep_genes <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds <- dds[keep_genes, ]
  
  # Run DESeq2 analysis
  dds <- DESeq(dds)
  
  # Extract results and convert to a data frame
  de_results <- results(dds, contrast = c("condition", "Experimental", "Control"))
  de_results_df <- as.data.frame(de_results)
  
  # Add gene_id as a column
  de_results_df$gene_id <- rownames(de_results_df)
  
  # Print the summary for this cluster's results
  cat("Summary for cluster", cluster, ":\n")
  print(summary(de_results))
  
  # Prepare the ranking metric for GSEA, using the stat score
  ranks <- de_results$stat
  names(ranks) <- rownames(de_results)
  
  # Filter out non-finite values from ranks
  ranks <- ranks[is.finite(ranks)]
  
  # Run GSEA using fgsea
  fgsea_results <- fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500)
  
  # Filter pathways by adjusted p-value (significant pathways)
  significant_pathways <- fgsea_results %>%
    dplyr::filter(padj < 0.05)
  
  # Sort pathways by NES (largest to smallest)
  sorted_pathways <- significant_pathways %>%
    dplyr::arrange(desc(NES))
  
  # Convert the `leadingEdge` list column to a character string for CSV export
  sorted_pathways$leadingEdge <- sapply(sorted_pathways$leadingEdge, function(x) paste(x, collapse = ";"))
  
  # Save GSEA results as CSV
  gsea_file_path <- file.path(cluster_dir, paste0("V2fgsea_significant_pathways_Cluster_", cluster, ".csv"))
  write.csv(sorted_pathways, file = gsea_file_path, row.names = FALSE)
  cat("GSEA results saved for cluster", cluster, "to:", gsea_file_path, "\n")
  
  # Select the top 25 most significant genes for labeling
  top_genes <- de_results_df %>%
    dplyr::filter(padj < 0.05) %>%
    arrange(padj) %>%
    head(15)
  
  # Save DESeq2 results as Excel, ensuring gene_id is included as a column
  excel_file_path <- file.path(cluster_dir, paste0("V2DESeq2_Results_Cluster_", cluster, ".xlsx"))
  write.xlsx(de_results_df, file = excel_file_path, rowNames = FALSE)
  cat("DESeq2 results saved for cluster", cluster, "to:", excel_file_path, "\n")
  
  # Generate the volcano plot for the current cluster
  volcano_plot <- ggplot(de_results_df, aes(x = log2FoldChange, y = -log10(padj), 
                                            color = ifelse(padj < 0.05, 
                                                           ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"), 
                                                           "Not Significant"))) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
    geom_text_repel(data = top_genes, aes(label = gene_id), size = 3, box.padding = 0.2, point.padding = 0.2) +
    labs(title = paste("Volcano Plot: DESeq2 Analysis - Cluster", cluster), 
         x = "Log2 Fold Change", y = "-log10 Adjusted P-value", color = "Legend") + 
    theme_minimal()
  
  # Save the volcano plot for the current cluster
  volcano_plot_path <- file.path(cluster_dir, paste0("V2Volcano_Plot_Cluster_", cluster, ".png"))
  ggsave(filename = volcano_plot_path, plot = volcano_plot, device = "png", width = 5, height = 5)
  cat("Volcano plot saved for cluster", cluster, "to:", volcano_plot_path, "\n")
}

cat("Cluster 0 processed successfully.\n")





# Assuming you have 'fgsea_results' from the previous analysis
# First, filter for significant pathways and sort by NES
library(dplyr)
library(ggplot2)

# Number of top pathways to plot
n_top_pathways <- 20

# Select the top n pathways based on NES
top_pathways <- fgsea_results %>%
  dplyr::filter(padj < 0.05) %>%       # Only include significant pathways
  dplyr::arrange(desc(NES)) %>%        # Sort by NES in descending order
  dplyr::slice(1:n_top_pathways)       # Select the top pathways

# Create the bar plot
pathway_plot <- ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES)) +  # Reorder by NES
  geom_bar(stat = "identity", aes(fill = NES > 0), show.legend = FALSE) +        # Positive/Negative NES coloring
  scale_fill_manual(values = c("blue", "orange")) +                             # Custom colors for positive/negative NES
  coord_flip() +                                                                # Flip coordinates for horizontal bars
  labs(title = "Top Pathway Enrichment by NES", x = "Pathway", y = "Normalized Enrichment Score (NES)") +
  theme_minimal()                                                               # Use minimal theme

# Show the plot
print(pathway_plot)

# Save the plot
ggsave(filename = file.path("C:/path/to/save", "Top_Pathway_Enrichment.png"), plot = pathway_plot, width = 8, height = 6, dpi = 300)


Cxcl9_plot <- FeaturePlot(fib_subset, features = "Cxcl9", reduction = "umap")
Cxcl9_plot

Cxcl10_plot <- FeaturePlot(fib_subset, features = "Cxcl10", reduction = "umap")
Cxcl10_plot


###Cxcl9 and 10 volplot
# Load required libraries
library(ggplot2)
library(readxl)
library(ggrepel)
library(dplyr)

# Define the file path to read in the DESeq2 results
file_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/DESeq2_outputs_allfib/Cluster0/DESeq2_Results_Cluster_0.xlsx"

# Read the DESeq2 results from the Excel file
de_results_df <- read_excel(file_path)

# Ensure that the gene names are included in the dataframe
if (!"gene_id" %in% colnames(de_results_df)) {
  stop("The 'gene_id' column is missing from the DESeq2 results.")
}

# Identify the top 15 upregulated and top 5 downregulated significant genes
top_genes_up <- de_results_df %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  arrange(padj) %>%
  head(15)

top_genes_down <- de_results_df %>%
  filter(padj < 0.05 & log2FoldChange < 0) %>%
  arrange(padj) %>%
  head(5)

# Combine the genes to label (top up and downregulated)
genes_to_label <- unique(c(top_genes_up$gene_id, top_genes_down$gene_id, "Cxcl9", "Cxcl10"))

# Generate the volcano plot
volcano_plot <- ggplot(de_results_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"), "Not Significant")), alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  
  # Add gene labels with lines pointing to the dots
  geom_text_repel(data = subset(de_results_df, gene_id %in% genes_to_label), 
                  aes(label = gene_id), 
                  size = 3, box.padding = 0.4, point.padding = 0.4, 
                  nudge_x = 0.15, nudge_y = 1,  # Nudge the labels slightly to create space for the lines
                  segment.color = 'black', segment.size = 0.7) +  # Add line with segment color and size
  
  labs(title = "Volcano Plot: Prostate Fibroblasts", 
       x = "Log2 Fold Change", y = "-log10 Adjusted P-value", color = "Legend") +
  theme_minimal()

# Define the file path to save the plot
output_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/DESeq2_outputs_allfib/Cluster0/Volcano_Plot_TopGenes_Cluster_0.png"

# Display the plot
print(volcano_plot)

# Save the volcano plot as a PNG file
ggsave(filename = output_path, plot = volcano_plot, width = 8, height = 6, dpi = 300)









####trying to make some volcano plots and running deseq2
# Load necessary libraries
library(Seurat)
library(DESeq2)
library(openxlsx)

# Load the Seurat object
fib_subset_pseudo_path <- "C:\\Users\\ostrobel\\Indiana University\\Jerde, Travis J - JERDE LAB (J Drive)\\Strobel\\Single Cell Analysis Stuff\\R_Analysis_Strobel\\Combined_Sample_Analysis\\Combined_analysis\\combined_obj_try4\\Harmony_SCT\\Fib_subset\\2024_10_07_fib_sub_psudo_FINAL.rds"
fib_subset.pseudo <- readRDS(fib_subset_pseudo_path)

# Define the cluster to run the analysis (Cluster 0)
clusters_to_run <- 0

# Specify the output directory
output_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure4"

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through the specified cluster(s)
for (cluster in clusters_to_run) {
  cat("Processing cluster:", cluster, "\n")
  
  # Create directory for the current cluster
  cluster_dir <- file.path(output_dir, paste0("Cluster", cluster))
  if (!dir.exists(cluster_dir)) {
    dir.create(cluster_dir, recursive = TRUE)
  }
  
  # Subset the pseudobulk object for the current cluster
  cluster_data <- subset(fib_subset.pseudo, idents = cluster)
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = GetAssayData(cluster_data, layer = "counts"),
    colData = cluster_data@meta.data,
    design = ~ condition
  )
  
  # Filter genes with fewer than 10 total counts in at least 3 samples
  smallestGroupSize <- 3
  keep_genes <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds <- dds[keep_genes, ]
  
  # Run DESeq2 analysis
  dds <- DESeq(dds)
  
  # Extract results and convert to a data frame
  de_results <- results(dds, contrast = c("condition", "Experimental", "Control"))
  de_results_df <- as.data.frame(de_results)
  
  # Add gene_id and label_gene columns
  de_results_df$gene_id <- rownames(de_results_df)
  de_results_df$label_gene <- FALSE  # Default column for gene labeling in the volcano plot
  
  # Save DESeq2 results as Excel
  excel_file_path <- file.path(cluster_dir, paste0("DESeq2_Results_Cluster_", cluster, ".xlsx"))
  write.xlsx(de_results_df, file = excel_file_path, rowNames = FALSE)
  cat("DESeq2 results saved for cluster", cluster, "to:", excel_file_path, "\n")
}

cat("Cluster processing completed and results saved to Excel.\n")













# Load required libraries
library(ggplot2)
library(readxl)
library(dplyr)

# Define the file path to read in the DESeq2 results Excel file
file_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure4/Cluster0/DESeq2_Results_Cluster_0.xlsx"

# Read the DESeq2 results from the Excel file
de_results_df <- read_excel(file_path)

# Ensure no missing values in log2FoldChange and padj columns
de_results_df <- de_results_df %>%
  filter(!is.na(log2FoldChange) & !is.na(padj))

# Generate a basic volcano plot
volcano_plot <- ggplot(de_results_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, 
                                ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"), 
                                "Not Significant")), 
             alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  labs(title = "Volcano Plot", 
       x = "Log2 Fold Change", y = "-log10 Adjusted P-value", color = "Legend") +
  theme_minimal()

# Display the plot
print(volcano_plot)





# Load required libraries
library(ggplot2)
library(readxl)
library(ggrepel)
library(dplyr)

# Define the file path to read in the DESeq2 results Excel file
file_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure4/Cluster0/DESeq2_Results_Cluster_0.xlsx"

# Read the DESeq2 results from the Excel file
de_results_df <- read_excel(file_path)

# Ensure no missing values in log2FoldChange, padj, and label_gene columns
de_results_df <- de_results_df %>%
  filter(!is.na(log2FoldChange) & !is.na(padj))

# Select genes to label based on the label_gene column being TRUE
genes_to_label <- de_results_df %>%
  filter(label_gene == TRUE)

# Generate the volcano plot with labeled genes
volcano_plot <- ggplot(de_results_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, 
                                ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"), 
                                "Not Significant")), 
             alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_text_repel(data = genes_to_label, aes(label = gene_id), 
                  size = 3, box.padding = 0.3, point.padding = 0.2, 
                  segment.color = "black") +
  labs(title = "Volcano Plot with Selected Gene Labels", 
       x = "Log2 Fold Change", y = "-log10 Adjusted P-value", color = "Legend") +
  theme_minimal()

# Display the plot
print(volcano_plot)








# Load required libraries
library(ggplot2)
library(readxl)
library(ggrepel)
library(dplyr)

# Define the file path to read in the DESeq2 results
file_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/DESeq2_outputs_allfib/Cluster0/DESeq2_Results_Cluster_0.xlsx"

# Read the DESeq2 results from the Excel file
de_results_df <- read_excel(file_path)

# Ensure that the gene names are included in the dataframe
if (!"gene_id" %in% colnames(de_results_df)) {
  stop("The 'gene_id' column is missing from the DESeq2 results.")
}

# Filter significant genes with padj < 0.05
significant_genes <- de_results_df %>%
  filter(padj < 0.05)

# Sort by adjusted p-value and log2 fold change
sorted_genes <- significant_genes %>%
  arrange(padj, desc(log2FoldChange))

# Select the top 15 upregulated (positive log2FoldChange) significant genes
top_15_upregulated <- sorted_genes %>%
  filter(log2FoldChange > 1) %>%
  head(15)

# Select the top 10 downregulated (negative log2FoldChange) significant genes
top_10_downregulated <- sorted_genes %>%
  filter(log2FoldChange < -1) %>%
  head(4)

# Combine the gene lists to label
genes_to_label <- c(top_15_upregulated$gene_id, top_10_downregulated$gene_id)

# Generate the volcano plot with log2FC threshold of -1 to 1 for "Not Significant"
volcano_plot <- ggplot(de_results_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, 
                                ifelse(log2FoldChange > 1, "Upregulated", "Downregulated"), 
                                "Not Significant")), alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  
  # Add gene labels with lines pointing to the dots for both upregulated and downregulated genes
  geom_text_repel(data = subset(de_results_df, gene_id %in% genes_to_label), 
                  aes(label = gene_id), 
                  size = 3, box.padding = 0.4, point.padding = 0.4, 
                  nudge_x = 0.15, nudge_y = 1,  # Nudge the labels slightly to create space for the lines
                  segment.color = 'black', segment.size = 0.7, max.overlaps = 30) +  # Add line with segment color and size
  
  labs(title = "Volcano Plot: Prostate Fibroblasts", 
       x = "Log2 Fold Change", y = "-log10 Adjusted P-value", color = "Legend") +
  theme_minimal()

# Define the file path to save the plot
output_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/DESeq2_outputs_allfib/Cluster0/Volcano_Plot_Top15_Up_Top10_Down.png"

# Display the plot
print(volcano_plot)




# Save the volcano plot as a PNG file
ggsave(filename = output_path, plot = volcano_plot, width = 6, height = 4, dpi = 300)








####committee Meeting image code 2024_10_21
# Load necessary libraries
library(Seurat)
library(ggplot2)

# Define file paths
seurat_object_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Grad School/Year 4/Sequencing Documents November 2023/Chrm_394_Arrizabalaga_Jerde_scRNAseq2plus8_Apr2023_Seurat_integrative_analysis_10102023.rds"
save_plot_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Grad School/2024/Presentations and such/2024_10_21_Com_meeting/FeaturePlot_Ptprc.png"

# Step 2: Read in the Seurat object
Strobel <- readRDS(seurat_object_path)

DefaultAssay(Strobel) <- "RNA"


# Step 3: Create a feature plot for "Ptprc"
feature_plot <- FeaturePlot(Strobel, features = "Ptprc", reduction = "umap") +
  ggtitle("Feature Plot of Ptprc") + 
  theme_minimal()

# Step 4: Save the feature plot as a PNG file
ggsave(filename = save_plot_path, plot = feature_plot, width = 4, height = 4, dpi = 300)

# Step 1: Subset cells in cluster 6
cluster_6_cells <- subset(Strobel, seurat_clusters == 6)

# Step 2: Check if cells express Ptprc and add it as a metadata column
cluster_6_cells$Ptprc_Expressed <- FetchData(cluster_6_cells, vars = "Ptprc") > 0

# Step 3: Calculate the percentage of cells expressing Ptprc, split by condition
percentage_expression <- cluster_6_cells@meta.data %>%
  group_by(condition) %>%
  summarise(
    Total_Cells = n(),
    Ptprc_Positive_Cells = sum(Ptprc_Expressed),
    Percentage_Expressing = (Ptprc_Positive_Cells / Total_Cells) * 100
  )

# Display the results
print(percentage_expression)




FeaturePlot(seurat_obj, features="Stat6", reduction = "umap.harmonyintegration", split.by = "condition")





###2024_10_28 Figure



##jerde figures------


# Minimal code to create and display a FeaturePlot for FN1
FeaturePlot(
  object    = seurat_obj,
  features  = "Fn1",
  reduction = "umap.harmonyintegration"
) + 
  ggtitle("FN1 Expression")


FeaturePlot(
  object    = fib_subset,
  features  = "Fn1",
  reduction = "umap.integration"
) + 
  ggtitle("FN1 Expression")


#To do:
#show where ifnar1, ifnar2, ifngr1, ifngr2 are located in all cells and fibroblasts.
#Show number of Fn1 positive cells, control vs experimental within cluster 0.
#maybe look at subtypes within cluster0?


#Only break up cluster 0 if it's a few very digestiable genes that make it different. 



library(Seurat)
library(patchwork)  # For wrap_plots
library(ggplot2)

# Define the gene list
genes_of_interest <- c(
  "Cxcl9", "H2-Aa", "H2-K1", "Cxcl10", "Cd74", "H2-Ab1", "Psmb8", 
  "B2m", "Ifit1", "H2-D1", "Ifi47", "Iigp1", "Gbp2", "Psmb9", 
  "Zbp1", "Cd274", "AW112010", "H2-Eb1", "H2-T23", "Isg15", 
  "Gbp3", "Bst2", "Igtp", "Gbp5", "Fn1", "Ifnar1", "Ifnar2", "Ifngr1", "Ifngr2"
)

genes_of_interest <- c(
  "Fn1", "Ifnar1", "Ifnar2", "Ifngr1", "Ifngr2"
)


# Create a list of FeaturePlots, one for each gene
plot_list <- lapply(genes_of_interest, function(gene) {
  FeaturePlot(
    object    = fib_subset,
    features  = gene,
    reduction = "umap"  # change if your reduction is named differently
  ) + 
    ggtitle(gene)
})

# Combine all FeaturePlots into one big figure
combined_plots <- wrap_plots(plot_list, ncol = 4)

# Save to the specified folder on your system
ggsave(
  filename = file.path(
    "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Jerde_images",
    "fib_subset_all_genes.png"
  ),
  plot  = combined_plots,
  width = 16,   # Adjust as needed
  height = 20   # Adjust as needed
)

cat("Saved FeaturePlots to 'C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Jerde_images'.\n")






library(Seurat)
library(ggplot2)
library(patchwork)

# Define your list of genes
genes_of_interest <- c(
  "Cxcl9", "H2-Aa", "H2-K1", "Cxcl10", "Cd74", "H2-Ab1", "Psmb8",
  "B2m", "Ifit1", "H2-D1", "Ifi47", "Iigp1", "Gbp2", "Psmb9",
  "Zbp1", "Cd274", "AW112010", "H2-Eb1", "H2-T23", "Isg15",
  "Gbp3", "Bst2", "Igtp", "Gbp5", "Fn1", "Ifnar1", "Ifnar2",
  "Ifngr1", "Ifngr2"
)



library(Seurat)
library(ggplot2)
library(patchwork)

# List of genes
genes_of_interest <- c("Cxcl9", "Fn1")  # Add more if needed

# Output directory
save_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/2025_02_12_Jerde_data"

# Loop over genes
for (gene in genes_of_interest) {
  
  # TOP PLOT: fib_subset, split by condition, force 2 columns
  plot_top <- FeaturePlot(
    object   = fib_subset,
    features = gene,
    reduction = "umap",   # adjust if your fib_subset uses a different reduction
    split.by  = "condition",
    ncol      = 2
  ) + ggtitle(paste(gene, "- fib_subset (Top)"))
  
  # BOTTOM PLOT: seurat_obj, split by condition, force 2 columns
  plot_bottom <- FeaturePlot(
    object   = seurat_obj,
    features = gene,
    reduction = "umap.harmonyintegration",  # or "umap" if that's your main embedding
    split.by  = "condition",
    ncol      = 2
  ) + ggtitle(paste(gene, "- seurat_obj (Bottom)"))
  
  # STACK THEM (top over bottom) using patchwork
  stacked_plot <- plot_top + plot_bottom + plot_layout(ncol = 1)
  
  # File name for each gene
  out_file <- file.path(save_dir, paste0(gene, "_fibTop_allBottom_byCondition.png"))
  
  # Save the stacked plot
  ggsave(
    filename = out_file,
    plot     = stacked_plot,
    width    = 10,
    height   = 8
  )
  
  cat("Saved stacked-by-condition plot for gene:", gene, "->", out_file, "\n")
}











###this makes a stacked plot with whole object top, stromal bottom, split by condition-----




# Load necessary packages
library(Seurat)
library(patchwork)

# Define colors (dark blue for high expression)
dark_blue <- c("lightgray", "darkblue") 

# 1) fib_subset FeaturePlot (top)
plot_top <- FeaturePlot(
  object   = fib_subset,
  features = gene_of_interest,
  reduction = "umap",     
  split.by  = "condition",
  ncol      = 2,
  cols      = dark_blue  # Set dark blue color gradient
) + ggtitle(paste(gene_of_interest, "- fib_subset (top)"))

# 2) seurat_obj FeaturePlot (bottom)
plot_bottom <- FeaturePlot(
  object   = seurat_obj,
  features = gene_of_interest,
  reduction = "umap.harmonyintegration",  
  split.by  = "condition",
  ncol      = 2,
  cols      = dark_blue  # Set dark blue color gradient
) + ggtitle(paste(gene_of_interest, "- all cells (bottom)"))

# Stack them
stacked_plot <- plot_top / plot_bottom  # Stack vertically

# Define file path
file_path <- file.path(folder, "Fn1_STACKED_byCondition.png")

# Save the plot
ggsave(
  filename = file_path,
  plot     = stacked_plot,
  width    = 8,
  height   = 8
)



# Load necessary packages
library(Seurat)
library(patchwork)

# List of genes
genes_of_interest <- c("Ccl5")


# Define colors (dark blue for high expression)
dark_blue <- c("lightgray", "darkblue") 

# Loop through each gene in genes_of_interest
for (gene_of_interest in genes_of_interest) {
  
  # 1) fib_subset FeaturePlot (top)
  plot_top <- FeaturePlot(
    object   = fib_subset,
    features = gene_of_interest,
    reduction = "umap",     
    split.by  = "condition",
    ncol      = 2,
    cols      = dark_blue  # Set dark blue color gradient
  ) + ggtitle(paste(gene_of_interest, "- fib_subset (top)"))
  
  # 2) seurat_obj FeaturePlot (bottom)
  plot_bottom <- FeaturePlot(
    object   = seurat_obj,
    features = gene_of_interest,
    reduction = "umap.harmonyintegration",  
    split.by  = "condition",
    ncol      = 2,
    cols      = dark_blue  # Set dark blue color gradient
  ) + ggtitle(paste(gene_of_interest, "- all cells (bottom)"))
  
  # Stack them
  stacked_plot <- plot_top / plot_bottom  # Stack vertically
  
  # Define dynamic file path with gene name
  file_path <- file.path(save_dir, paste0(gene_of_interest, "_STACKED_byCondition.png"))
  
  # Save the plot
  ggsave(
    filename = file_path,
    plot     = stacked_plot,
    width    = 8,
    height   = 8
  )
}






#time to make a few more for fun 
###immune sig analysis!-----
# ====== Pakete ======
library(readxl)
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(Matrix)
library(reshape2)

# ====== Deine Inputs ======
obj <- fib_subset2
Idents(obj) <- "seurat_clusters"
cl_of_interest <- "3"

de_file <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure4/Cluster3/V2DESeq2_Results_Cluster_3.xlsx"
out_dir <- dirname(de_file)

# ====== DESeq2 laden & Ranking ======
res_cl3 <- read_xlsx(de_file)
stopifnot(all(c("gene_id","log2FoldChange","padj") %in% colnames(res_cl3)))
res_cl3 <- res_cl3 %>%
  mutate(padj = ifelse(is.na(padj), 1, padj),
         rank = sign(log2FoldChange) * -log10(pmax(padj, 1e-300)))

# ====== Module ======
isg  <- c("Isg15","Ifi6","Ifit1","Ifit2","Ifit3","Ifi44","Oasl1","Oasl2","Rsad2",
          "Mx1","Bst2","Xaf1","Stat1","Stat2","Irf7","Irf1","Gbp2","Gbp4","Gbp5")
chem <- c("Cxcl9","Cxcl10","Cxcl11","Ccl2","Ccl7","Ccl5","Cxcl13","Ccl19","Ccl21")
mhc2 <- c("Cd74","H2-Aa","H2-Ab1","H2-Eb1","Ciita")
prr  <- c("Ddx58","Ifih1","Tlr3","Tlr2","Nfkbia","Relb","Ccl20","Icam1")

# ====== Shortlist aus DE + Modulen ======
res_up <- res_cl3 %>% filter(log2FoldChange > 0, padj < 0.05)

pick_top <- function(df_up, mod_genes, k=5) {
  cand <- df_up %>%
    filter(toupper(gene_id) %in% toupper(mod_genes)) %>%
    arrange(desc(rank))
  head(cand$gene_id, k)
}

shortlist <- unique(c(
  pick_top(res_up, isg, 5),
  pick_top(res_up, chem, 5),
  pick_top(res_up, mhc2, 3),
  pick_top(res_up, prr, 3)
))

# Falls zu kurz: mit Top-Up auffüllen
if (length(shortlist) < 10) {
  extra <- setdiff(res_up %>% arrange(desc(rank)) %>% pull(gene_id), shortlist)
  shortlist <- unique(c(shortlist, head(extra, 10 - length(shortlist))))
}

# Auf Gene beschränken, die im Seurat-Objekt vorhanden sind
shortlist <- shortlist[toupper(shortlist) %in% toupper(rownames(obj))]
# Optional: Auf max. 20 Gene begrenzen, damit die Figur „clean“ bleibt
shortlist <- head(shortlist, 20)

message("Shortlist Größe: ", length(shortlist))
print(shortlist)

# ====== RAM-sicheres DotPlot-Äquivalent ======
DefaultAssay(obj) <- "SCT"
# Nur die benötigten Gene als sparse Submatrix holen
mat <- GetAssayData(obj, assay = "SCT", layer = "data")[shortlist, , drop = FALSE]  # dgCMatrix, klein
clusters <- Idents(obj)
cl_levels <- levels(clusters)

# Mittelwert und %>0 je Cluster (sparse-freundlich)
avg_by_cluster <- sapply(cl_levels, function(cl) {
  idx <- which(clusters == cl)
  Matrix::rowMeans(mat[, idx, drop = FALSE])
})
pct_by_cluster <- sapply(cl_levels, function(cl) {
  idx <- which(clusters == cl)
  Matrix::rowMeans(mat[, idx, drop = FALSE] > 0) * 100
})

avg_df <- reshape2::melt(as.matrix(avg_by_cluster), varnames = c("gene","cluster"), value.name = "avg")
pct_df <- reshape2::melt(as.matrix(pct_by_cluster), varnames = c("gene","cluster"), value.name = "pct")
dot_df <- dplyr::left_join(avg_df, pct_df, by = c("gene","cluster"))

# Plot (wie Seurat::DotPlot: Farbe=avg, Größe=pct)
p_dot_safe <- ggplot(dot_df, aes(x = cluster, y = gene)) +
  geom_point(aes(size = pct, color = avg)) +
  scale_color_viridis_c(name = "Avg exp") +
  scale_size(range = c(1, 8), name = "% expr") +
  coord_flip() +
  theme_classic(base_size = 11) +
  labs(title = "Immune-like shortlist across clusters (RAM-safe)",
       x = "Cluster", y = "Gene")

p_dot_safe

ggsave(file.path(out_dir, "dotplot_cluster3_shortlist_RAMsafe.png"),
       p_dot_safe, width = 9.5, height = 6, dpi = 300)


View(cl3_obj@meta.data)

# ====== Heatmap: Cluster 3 pro Sample (pseudobulk-like) ======
cl3_obj <- subset(obj, idents = cl_of_interest)
avg_cl3 <- AverageExpression(cl3_obj, assays = "SCT", layer = "data", group.by = "orig.ident")$SCT
samples <- colnames(avg_cl3)
meta <- data.frame(condition = cl3_obj$condition[match(samples, cl3_obj$orig.ident)],
                   row.names = samples)

# --- Auswahl, Matrix & Row-Z-Scaling ---
genes_heat <- intersect(shortlist, rownames(avg_cl3))
mat2 <- as.matrix(avg_cl3[genes_heat, , drop = FALSE])

# Row-Z (Zeilenweise skalieren); NAs auf 0
mat2_z <- t(scale(t(mat2)))
mat2_z[is.na(mat2_z)] <- 0

# --- Farbskala: 0 = weiß, symmetrisch um 0 ---
n_col <- 101  # ungerade Anzahl -> mittlere Farbe entspricht genau 0
lim   <- max(abs(range(mat2_z)), na.rm = TRUE)
breaks <- seq(-lim, lim, length.out = n_col + 1)
cols   <- colorRampPalette(c("navy", "white", "goldenrod"))(n_col)

# --- PNG export wie gehabt ---
png(file.path(out_dir, "heatmap_cluster3_shortlist_by_sample.png"),
    width = 8.5, height = 6, units = "in", res = 300)
pheatmap(mat2_z,
         cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = meta,
         color = cols, breaks = breaks,
         main = paste0("Cluster ", cl_of_interest, " — immune-like shortlist (row-z)")
)
dev.off()

library(svglite)

svglite(file.path(out_dir, "heatmap_cluster3_shortlist_by_sample.svg"),
        width = 8.5, height = 6)

pheatmap(
  mat2_z,
  cluster_rows = TRUE, cluster_cols = TRUE,
  annotation_col = meta,
  color  = cols,
  breaks = breaks,
  main   = paste0("Cluster ", cl_of_interest, " — immune-like shortlist (row-z)")
)

dev.off()


# --- Zusätzlich: SVG export (0 bleibt exakt weiß) ---
svg(file.path(out_dir, "heatmap_cluster3_shortlist_by_sample.svg"),
    width = 8.5, height = 6)
pheatmap(mat2_z,
         cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = meta,
         color = cols, breaks = breaks,
         main = paste0("Cluster ", cl_of_interest, " — immune-like shortlist (row-z)")
)
dev.off()




# ====== Optional: Cluster-Average Heatmap nur Immun-Gene (klein halten) ======
immune_union <- unique(c(isg, chem, mhc2, prr))
imm_in_obj  <- immune_union[toupper(immune_union) %in% toupper(rownames(obj))]
imm_in_obj  <- head(imm_in_obj, 60)  # cap zur Sicherheit

avg_all <- AverageExpression(obj, assays = "SCT", layer = "data",
                             group.by = "seurat_clusters", features = imm_in_obj)$SCT
mat3 <- as.matrix(avg_all)
mat3_z <- t(scale(t(mat3))); mat3_z[is.na(mat3_z)] <- 0

# Min/Max holen
lim <- max(abs(mat3_z), na.rm = TRUE)

# Breaks symmetrisch um 0
bk <- seq(-lim, lim, length.out = 101)

pheatmap(
  mat3_z,
  filename = file.path(out_dir, "heatmap_immune_signatures_by_cluster.png"),
  width = 9, height = 6, dpi = 300,
  cluster_rows = TRUE, cluster_cols = TRUE,
  show_rownames = TRUE, show_colnames = TRUE,
  color = colorRampPalette(c("navy","white","goldenrod"))(100),
  breaks = bk,
  main = "Immune signatures across clusters (SCT, row-z)"
)


# ====== Shortlist export ======
shortlist_df <- res_cl3 %>%
  filter(toupper(gene_id) %in% toupper(shortlist)) %>%
  mutate(module = case_when(
    toupper(gene_id) %in% toupper(isg)  ~ "ISG",
    toupper(gene_id) %in% toupper(chem) ~ "Chemokines",
    toupper(gene_id) %in% toupper(mhc2) ~ "MHCII",
    toupper(gene_id) %in% toupper(prr)  ~ "PRR_NFkB",
    TRUE ~ "Other"
  )) %>%
  arrange(desc(rank)) %>%
  select(gene_id, module, log2FoldChange, padj, rank)

write.csv(shortlist_df,
          file.path(out_dir, "immune_like_shortlist_cluster3.csv"),
          row.names = FALSE)



library(Seurat)
library(pheatmap)

# 1) Immun-Genliste (wie zuvor)
isg  <- c("Isg15","Ifi6","Ifit1","Ifit2","Ifit3","Ifi44","Oasl1","Oasl2","Rsad2",
          "Mx1","Bst2","Xaf1","Stat1","Stat2","Irf7","Irf1","Gbp2","Gbp4","Gbp5")
chem <- c("Cxcl9","Cxcl10","Cxcl11","Ccl2","Ccl7","Ccl5","Cxcl13","Ccl19","Ccl21")
mhc2 <- c("Cd74","H2-Aa","H2-Ab1","H2-Eb1","Ciita")
prr  <- c("Ddx58","Ifih1","Tlr3","Tlr2","Nfkbia","Relb","Ccl20","Icam1")
gset <- unique(c(isg,chem,mhc2,prr))
gset <- gset[toupper(gset) %in% toupper(rownames(obj))]

DefaultAssay(obj) <- "SCT"

# 2) Cluster-Mittelwerte separat für Control/Experimental
obj_ctrl <- subset(obj, subset = condition == "Control")
obj_exp  <- subset(obj, subset = condition == "Experimental")

avg_ctrl <- AverageExpression(obj_ctrl, assays="SCT", layer="data",
                              group.by="seurat_clusters", features=gset)$SCT
avg_exp  <- AverageExpression(obj_exp,  assays="SCT", layer="data",
                              group.by="seurat_clusters", features=gset)$SCT

# 3) gemeinsame Cluster auswählen & Δ bilden
common_cl <- intersect(colnames(avg_ctrl), colnames(avg_exp))
mat_delta <- as.matrix(avg_exp[gset, common_cl, drop=FALSE] -
                         avg_ctrl[gset, common_cl, drop=FALSE])

# 4) Farbskala um 0 zentrieren (0 = weiß)
lim <- max(abs(mat_delta), na.rm=TRUE)
bk  <- seq(-lim, lim, length.out=101)

png(file.path(out_dir, "heatmap_cluster3_delta.png"),
    width = 9, height = 6, units = "in", res = 300)

pheatmap(
  mat_delta,
  cluster_rows = TRUE, cluster_cols = TRUE,
  color  = colorRampPalette(c("navy","white","goldenrod"))(100),
  breaks = bk,
  main   = "Immune Gene Expression Change During Infection (per cluster)"
)

dev.off()

library(svglite)

svglite(file.path(out_dir, "heatmap_cluster3_delta.svg"),
        width = 9, height = 6)

pheatmap(
  mat_delta,
  cluster_rows = TRUE, cluster_cols = TRUE,
  color  = colorRampPalette(c("navy","white","goldenrod"))(100),
  breaks = bk,
  main   = "Immune Gene Expression Change During Infection (per cluster)"
)

dev.off()




# === Packages ===
library(readxl)
library(pheatmap)
library(tools)

# === Pfade ===
base_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure4"
out_dir  <- file.path(base_dir, "DESeq2Heatmaps")
if (!dir.exists(out_dir)) dir.create(out_dir)

# === Immun-Genliste (Mouse) ===
isg  <- c("Isg15","Ifi6","Ifit1","Ifit2","Ifit3","Ifi44","Oasl1","Oasl2","Rsad2",
          "Mx1","Bst2","Xaf1","Stat1","Stat2","Irf7","Irf1","Gbp2","Gbp4","Gbp5")
chem <- c("Cxcl9","Cxcl10","Cxcl11","Ccl2","Ccl7","Ccl5","Cxcl13","Ccl19","Ccl21")
mhc2 <- c("Cd74","H2-Aa","H2-Ab1","H2-Eb1","Ciita")
prr  <- c("Ddx58","Ifih1","Tlr3","Tlr2","Nfkbia","Relb","Ccl20","Icam1")
gset <- unique(c(isg, chem, mhc2, prr))

# === Helper: erste passende Spalte (case-insensitive) finden ===
first_match_col <- function(df, candidates) {
  cn <- tolower(colnames(df))
  cand <- tolower(candidates)
  hit <- which(cn %in% cand)
  if (length(hit) == 0) return(NA_character_)
  colnames(df)[hit[1]]
}

# === Cluster-Ordner erkennen (z. B. "Cluster3" oder "Cluster 3") ===
all_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
cl_dirs  <- all_dirs[ grepl("^Cluster\\s*\\d+$", basename(all_dirs), ignore.case = TRUE) ]

if (length(cl_dirs) == 0) {
  message("Keine Cluster-Ordner gefunden. Verfügbar unter Figure4: ",
          paste(basename(all_dirs), collapse=", "))
}

# === Farbschema ===
cols <- colorRampPalette(c("navy","white","goldenrod"))(100)

for (d in cl_dirs) {
  cl_base <- basename(d)                        # z.B. "Cluster3" oder "Cluster 3"
  cl_num  <- sub(".*?(\\d+)$", "\\1", cl_base) # extrahiere Ziffern am Ende
  # Excel-Dateien, die DESeq2 & Cluster-Nummer enthalten (case-insensitive)
  xls <- list.files(d, pattern = paste0("(?i)DESeq2.*Cluster[_ ]?", cl_num, ".*\\.(xlsx|xls)$"),
                    full.names = TRUE)
  if (length(xls) == 0) {
    message("Überspringe ", cl_base, " (keine passende DESeq2-Datei gefunden)")
    next
  }
  # Neueste Datei nehmen (nach Änderungszeit)
  f <- xls[ which.max(file.info(xls)$mtime) ]
  message("Lese: ", f)
  
  # Excel einlesen (erstes Sheet)
  res <- read_xlsx(f)
  
  # Spalten zuordnen (robust gg. unterschiedliche Namen)
  gene_col <- first_match_col(res, c("gene_id","gene","symbol","genes","gene_symbol"))
  lfc_col  <- first_match_col(res, c("log2foldchange","log2_fc","lfc","logfc"))
  padj_col <- first_match_col(res, c("padj","qvalue","q_value","padjusted","fdr"))
  
  if (any(is.na(c(gene_col, lfc_col, padj_col)))) {
    message("Überspringe ", cl_base, " (Spalten nicht gefunden: gene/log2FC/padj). Spalten sind: ",
            paste(colnames(res), collapse=", "))
    next
  }
  
  # Filter auf Immun-Gene (case-insensitive)
  res$..gene <- as.character(res[[gene_col]])
  keep <- toupper(res$..gene) %in% toupper(gset)
  df <- res[keep, c(gene_col, lfc_col, padj_col)]
  colnames(df) <- c("gene", "log2FoldChange", "padj")
  
  if (nrow(df) == 0) {
    message("Keine Immun-Gene in ", cl_base, " gefunden.")
    next
  }
  
  rownames(df) <- df$gene
  # Heatmap-Matrix: eine Spalte (dieser Cluster), Zeilen = Gene
  mat <- as.matrix(df["log2FoldChange"])
  colnames(mat) <- cl_base
  
  # Skala um 0 zentrieren
  lim <- max(abs(mat), na.rm = TRUE); lim <- ifelse(is.finite(lim), lim, 1)
  bk  <- seq(-lim, lim, length.out = 101)
  
  # Sterne für Signifikanz
  stars <- ifelse(df$padj < 0.05 & !is.na(df$padj), "*", "")
  
  # Ausgabe-Datei
  svg_file <- file.path(out_dir, paste0("Heatmap_", cl_base, ".svg"))
  
  # SVG speichern
  svg(svg_file, width = 4.2, height = max(4, 0.28 * nrow(mat)))
  pheatmap(mat,
           cluster_rows = TRUE, cluster_cols = FALSE,
           color = cols, breaks = bk,
           display_numbers = matrix(stars, nrow = nrow(mat), ncol = 1,
                                    dimnames = list(rownames(mat), cl_base)),
           number_color = "black",
           main = paste0("DESeq2 log2FC — ", cl_base))
  dev.off()
  
  message("Gespeichert: ", svg_file)
}



####making plot comparing IPA data---------
# =========================
# Panel B – IPA grouped dotplot
# =========================

# Pakete
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)

theme_set(theme_minimal(base_size = 12))

# ---- Pfade / Cluster ----
base_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/IPA_analysis"
clusters <- paste0("Fib", 0:4)  # ggf. anpassen (z.B. 0:7)

outfile_png <- file.path(base_dir, "Fig4B_IPA_grouped_dotplot_Fib0-4.png")
outfile_csv <- file.path(base_dir, "Fig4B_IPA_grouped_values_Fib0-4.csv")
outfile_map <- file.path(base_dir, "Fig4B_IPA_grouped_mapping_Fib0-4.csv")

# ---- Einlesen einer IPA-Datei (Header ist immer in Zeile 2) ----
read_ipa_file <- function(path_xls, cl){
  stopifnot(file.exists(path_xls))
  df <- read_excel(path_xls, skip = 1)      # Zeile 1 = Copyright, also skip=1
  names(df) <- str_trim(names(df))          # Header trimmen (entfernt führende Leerzeichen)
  
  # Header müssen exakt vorhanden sein (nach Trim):
  must <- c("Ingenuity Canonical Pathways", "-log(p-value)", "Ratio", "z-score", "Molecules")
  if (!all(must %in% names(df))) {
    stop("Unerwartete Spaltennamen in: ", path_xls,
         "\nGefunden: ", paste(names(df), collapse=" | "),
         "\nErwartet: ", paste(must, collapse=" | "))
  }
  
  # Datenspalten extrahieren und p, padj berechnen
  out <- df %>%
    transmute(
      cluster   = gsub("^Fib", "", cl),
      pathway   = `Ingenuity Canonical Pathways`,
      logp      = as.numeric(`-log(p-value)`),
      ratio     = as.numeric(Ratio),
      z         = suppressWarnings(as.numeric(`z-score`)),  # "#NUM!" -> NA
      molecules = Molecules
    ) %>%
    mutate(
      p    = 10^(-logp),
      padj = p.adjust(p, "BH")
    )
  out
}

# ---- Mapping: IPA-Pfade -> Biologie-Gruppen (priorisiert; erstes Match gewinnt) ----
groups_tbl <- tibble::tribble(
  ~group,                    ~pattern,
  "Interferon (JAK/STAT)",  "INTERFERON|\\bIFN\\b|\\bSTAT\\b|\\bJAK\\b",
  "Antigen presentation",   "ANTIGEN|\\bMHC\\b|\\bHLA\\b|\\bTAP\\b",
  "NFkB/TNF/IL1",           "NF.?KB|\\bTNF\\b|IL-?1",
  "Chemokine/IL6",          "CHEMOKINE|CXCL|CCL|IL-?6|CYTOKINE|CYTOKINE STORM",
  "PRR/IRF (RIG-I/MDA5)",   "PATTERN|\\bPRR\\b|RIG|MDA5|IFIH1|\\bIRF\\b",
  "T/NK/Adaptive/Checkpoint","\\bTH1\\b|\\bTH2\\b|\\bTCR\\b|PD-?1|PD-?L1|\\bNK\\b|CHECKPOINT",
  "Ubiquitination/ISGylation","UBIQUITIN|ISG(Y|)LATION|\\bISG\\b",
  "Complement/Acute phase", "COMPLEMENT|ACUTE.*PHASE",
  "TGFb/Fibrosis/SMAD",     "\\bTGF\\b|\\bSMAD\\b|FIBROS",
  "VEGF/Angiogenesis",      "\\bVEGF\\b|ANGIO",
  "IGF/IGFBP",              "\\bIGF\\b|IGFBP",
  "Stress/HSP/GR",          "STRESS|HSPA|HSPB|GLUCOCORTICOID",
  "MAPK/WNT/EGFR/ERK",      "MAPK|\\bWNT\\b|\\bEGFR\\b|\\bERK\\b"
)

assign_group <- function(pathway_name){
  U <- toupper(as.character(pathway_name))
  idx <- which(stringr::str_detect(U, groups_tbl$pattern))
  if (length(idx) == 0) "Other" else groups_tbl$group[idx[1]]  # Priorität: erstes Match
}

# ---- Alle Cluster laden ----
ipa_all <- dplyr::bind_rows(lapply(clusters, function(cl){
  read_ipa_file(file.path(base_dir, cl, "all pathways.xls"), cl)
}))

# Vollständiges Mapping (für Supplement)
ipa_map <- ipa_all %>%
  mutate(group = vapply(pathway, assign_group, FUN.VALUE = character(1)))

# ---- Repräsentativster Pfad je Cluster×Gruppe wählen ----
# 1) wenn z vorhanden: nimm max |z|
ipa_sum_z <- ipa_map %>%
  filter(is.finite(z)) %>%
  group_by(cluster, group) %>%
  slice_max(order_by = abs(z), n = 1, with_ties = FALSE) %>%
  ungroup()

# 2) falls z fehlt: bestes padj als Fallback
ipa_sum_p <- ipa_map %>%
  anti_join(ipa_sum_z, by = c("cluster","group")) %>%
  group_by(cluster, group) %>%
  slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
  mutate(z = NA_real_) %>%
  ungroup()

ipa_sum <- bind_rows(ipa_sum_z, ipa_sum_p)

# C3-Spezifität berechnen (einmal)
spec_tbl <- ipa_sum %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(spec_C3 = {
    z_c3  <- z[cluster=="3"]; z_oth <- z[cluster!="3"]
    if (length(z_c3)==0 || all(is.na(z_c3))) {
      (-log10(padj[cluster=="3"])) - median(-log10(padj[cluster!="3"]), na.rm=TRUE)
    } else {
      z_c3 - median(z_oth, na.rm=TRUE)
    }
  }, .groups="drop")
spec_tbl
# Level-Reihenfolge aus spec_tbl
group_levels <- spec_tbl %>%
  dplyr::arrange(spec_C3) %>%           # oder desc(spec_C3) wenn du C3-stark oben willst
  dplyr::pull(group)

# Einmal joinen, dann Levels setzen (KEIN zweites left_join mehr!)
ipa_sum <- ipa_sum %>%
  dplyr::left_join(spec_tbl, by = "group") %>%
  dplyr::mutate(
    group   = factor(group, levels = group_levels),
    cluster = factor(cluster, levels = sort(unique(cluster)))
  )
print(n=70, ipa_sum)




#plot
# |z| in 4 Kategorien mappen: ≤1, (1,2], (2,3], >3
dd <- ipa_sum %>%
  mutate(
    z_cat = cut(abs(z),
                breaks = c(-Inf, 1, 2, 3, Inf),
                labels = c("1","2","3","4"),
                right = TRUE)  # 1=≤1, 2=(1,2], 3=(2,3], 4=>3
  )

# ---- Plot mit manuellen Punktgrößen ----
p <- ggplot(dd, aes(x = cluster, y = group)) +
  geom_point(
    aes(size  = z_cat,                 # diskrete Größe
        color = -log10(padj)),
    alpha = 0.95
  ) +
  geom_text(
    aes(label = ifelse(is.na(z), "", ifelse(z >= 0, "+", "−"))),
    size = 3, nudge_x = 0.18, alpha = .7
  ) +
  scale_color_viridis_c(name = "-log10(adj p)") +
  # HIER legst du die Größen fest (in mm). Passe die Zahlen nach Geschmack an:
  scale_size_manual(
    name   = "|z|-Kategorie",
    values = c("1" = 1,   # ≤1
               "2" = 3,   # 1–2
               "3" = 7,  # 2–3
               "4" = 12), # >3
    breaks = c("1","2","3","4"),
    labels = c("≤1","1–2","2–3",">3"),
    drop   = FALSE
  ) +
  labs(x = "Fibroblast-Cluster", y = NULL,
       title = "IPA grouped pathways per cluster (repräsentativer Pfad je Gruppe)") +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"),
        plot.title  = element_text(face = "bold"))

p


# ---- Speichern ----
ggsave(outfile_png, p, width = 9, height = 6.2, dpi = 300)
write.csv(ipa_sum %>% arrange(group, as.numeric(as.character(cluster))), outfile_csv, row.names = FALSE)
write.csv(ipa_map %>% arrange(as.numeric(as.character(cluster)), group, padj), outfile_map, row.names = FALSE)

message("Gespeichert:\n  ", outfile_png, "\n  ", outfile_csv, "\n  ", outfile_map)


# globales Minimum & Maximum der z-Spalte
min_z <- min(ipa_sum$z, na.rm = TRUE)
max_z <- max(ipa_sum$z, na.rm = TRUE)
min_z; max_z


print(n=70, ipa_sum)


##try2
# =========================
# Fig. 4B – IPA grouped dotplot (manuelle Größen-Skala 1/3/7/12)
# =========================

# Pakete
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)

theme_set(theme_minimal(base_size = 12))

# ---- Pfade / Cluster ----
base_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/IPA_analysis"
clusters <- paste0("Fib", 0:4)  # ggf. anpassen

outfile_png <- file.path(base_dir, "Fig4B_IPA_grouped_dotplot_Fib0-4_manualsize.png")
outfile_csv <- file.path(base_dir, "Fig4B_IPA_grouped_values_Fib0-4.csv")
outfile_map <- file.path(base_dir, "Fig4B_IPA_grouped_mapping_Fib0-4.csv")

# ---- Einlesen einer IPA-Datei (Zeile 1 Copyright, Zeile 2 Header) ----
read_ipa_file <- function(path_xls, cl){
  stopifnot(file.exists(path_xls))
  df <- read_excel(path_xls, skip = 1)
  names(df) <- str_trim(names(df))  # Header trimmen
  
  must <- c("Ingenuity Canonical Pathways", "-log(p-value)", "Ratio", "z-score", "Molecules")
  if (!all(must %in% names(df))) {
    stop("Unerwartete Spaltennamen in: ", path_xls,
         "\nGefunden: ", paste(names(df), collapse=" | "),
         "\nErwartet: ", paste(must, collapse=" | "))
  }
  
  df %>%
    transmute(
      cluster      = gsub("^Fib", "", cl),
      pathway      = `Ingenuity Canonical Pathways`,
      logp         = as.numeric(`-log(p-value)`),
      ratio        = as.numeric(Ratio),
      z            = suppressWarnings(as.numeric(`z-score`)),  # "#NUM!" -> NA
      molecules    = Molecules,
      p            = 10^(-as.numeric(`-log(p-value)`)),
      padj         = p.adjust(10^(-as.numeric(`-log(p-value)`)), "BH")
    )
}

# ---- Mapping: IPA-Pfade -> Biologie-Gruppen (priorisiert; erstes Match gewinnt) ----
groups_tbl <- tibble::tribble(
  ~group,                    ~pattern,
  "Interferon (JAK/STAT)",  "INTERFERON|\\bIFN\\b|\\bSTAT\\b|\\bJAK\\b",
  "Antigen presentation",   "ANTIGEN|\\bMHC\\b|\\bHLA\\b|\\bTAP\\b",
  "NFkB/TNF/IL1",           "NF.?KB|\\bTNF\\b|IL-?1",
  "Chemokine/IL6",          "CHEMOKINE|CXCL|CCL|IL-?6|CYTOKINE|CYTOKINE STORM",
  "PRR/IRF (RIG-I/MDA5)",   "PATTERN|\\bPRR\\b|RIG|MDA5|IFIH1|\\bIRF\\b",
  "T/NK/Adaptive/Checkpoint","\\bTH1\\b|\\bTH2\\b|\\bTCR\\b|PD-?1|PD-?L1|\\bNK\\b|CHECKPOINT",
  "Ubiquitination/ISGylation","UBIQUITIN|ISG(Y|)LATION|\\bISG\\b",
  "Complement/Acute phase", "COMPLEMENT|ACUTE.*PHASE",
  "TGFb/Fibrosis/SMAD",     "\\bTGF\\b|\\bSMAD\\b|FIBROS",
  "VEGF/Angiogenesis",      "\\bVEGF\\b|ANGIO",
  "IGF/IGFBP",              "\\bIGF\\b|IGFBP",
  "Stress/HSP/GR",          "STRESS|HSPA|HSPB|GLUCOCORTICOID",
  "MAPK/WNT/EGFR/ERK",      "MAPK|\\bWNT\\b|\\bEGFR\\b|\\bERK\\b"
)
assign_group <- function(pathway_name){
  U <- toupper(as.character(pathway_name))
  idx <- which(stringr::str_detect(U, groups_tbl$pattern))
  if (length(idx) == 0) "Other" else groups_tbl$group[idx[1]]
}

# ---- Daten laden & gruppieren ----
ipa_all <- dplyr::bind_rows(lapply(clusters, function(cl){
  read_ipa_file(file.path(base_dir, cl, "all pathways.xls"), cl)
}))

# Vollständiges Mapping (für Supplement/Transparenz)
ipa_map <- ipa_all %>%
  mutate(group = vapply(pathway, assign_group, FUN.VALUE = character(1)))

# Repräsentativster Pfad je Cluster×Gruppe:
# 1) mit z: max |z|
ipa_sum_z <- ipa_map %>%
  filter(is.finite(z)) %>%
  group_by(cluster, group) %>%
  slice_max(order_by = abs(z), n = 1, with_ties = FALSE) %>%
  ungroup()
# 2) ohne z: bestes padj (Fallback)
ipa_sum_p <- ipa_map %>%
  anti_join(ipa_sum_z, by = c("cluster","group")) %>%
  group_by(cluster, group) %>%
  slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
  mutate(z = NA_real_) %>%
  ungroup()

ipa_sum <- bind_rows(ipa_sum_z, ipa_sum_p)

# ---- C3-Spezifität (nur für Sortierung der y-Achse) ----
spec_tbl <- ipa_sum %>%
  group_by(group) %>%
  summarise(spec_C3 = {
    z_c3  <- z[cluster=="3"]; z_oth <- z[cluster!="3"]
    if (length(z_c3)==0 || all(is.na(z_c3))) {
      (-log10(padj[cluster=="3"])) - median(-log10(padj[cluster!="3"]), na.rm=TRUE)
    } else {
      z_c3 - median(z_oth, na.rm=TRUE)
    }
  }, .groups="drop")

group_levels <- spec_tbl %>% arrange(desc(spec_C3)) %>% pull(group)

# ---- Zwei Namensspalten: pathway_name (original), plot_y (für die Grafik) ----
# Mapping alt -> neu (Group-Label -> hübsches y-Achsenlabel)
label_map <- c(
  "Interferon (JAK/STAT)"    = "Interferon / JAK–STAT",
  "Antigen presentation"     = "Antigenpräsentation",
  "NFkB/TNF/IL1"             = "NFκB / TNF / IL-1",
  "Chemokine/IL6"            = "Chemokine / IL-6",
  "PRR/IRF (RIG-I/MDA5)"     = "PRR / IRF (RIG-I/MDA5)",
  "T/NK/Adaptive/Checkpoint" = "T/NK/Adaptive/Checkpoint",
  "Ubiquitination/ISGylation"= "Ubiquitination / ISGyl.",
  "Complement/Acute phase"   = "Komplement / Akute Phase",
  "TGFb/Fibrosis/SMAD"       = "TGFβ / Fibrose / SMAD",
  "VEGF/Angiogenesis"        = "VEGF / Angiogenese",
  "IGF/IGFBP"                = "IGF / IGFBP",
  "Stress/HSP/GR"            = "Stress / HSP / GR",
  "MAPK/WNT/EGFR/ERK"        = "MAPK / WNT / EGFR / ERK"
)
plot_y_levels <- dplyr::recode(group_levels, !!!label_map, .default = group_levels)

ipa_sum <- ipa_sum %>%
  dplyr::left_join(spec_tbl, by = "group") %>%
  dplyr::mutate(
    pathway_name = pathway,                              # Original-IPA-Name
    group        = factor(group, levels = group_levels), # sortiert nach C3
    plot_y       = factor(                               # hübsche y-Labels (mit gleicher Sortierung)
      dplyr::recode(as.character(group), !!!label_map,
                    .default = as.character(group)),
      levels = plot_y_levels
    ),
    cluster      = factor(stringr::str_trim(as.character(cluster)),
                          levels = sort(unique(stringr::str_trim(as.character(cluster)))))
  )

# ---- Plot-Daten: diskrete Größenklassen 1..4; NA -> "1" (sichtbar, klein) ----
dd <- ipa_sum %>%
  mutate(
    padj_plot = pmax(padj, .Machine$double.xmin),         # schützt vor Inf
    z_cat_raw = cut(abs(z), breaks = c(-Inf, 1, 2, 3, Inf),
                    labels = c("1","2","3","4"), right = TRUE),
    z_cat     = factor(ifelse(is.na(z_cat_raw), "1", as.character(z_cat_raw)),
                       levels = c("1","2","3","4"))
  )

# ---- Plot ----
p <- ggplot(dd, aes(x = cluster, y = plot_y)) +
  geom_point(
    aes(size = z_cat, color = -log10(padj_plot)),
    alpha = 0.95
  ) +
  geom_text(
    aes(label = ifelse(is.na(z), "", ifelse(z >= 0, "+", "−"))),
    size = 3, nudge_x = 0.18, alpha = .7
  ) +
  scale_color_viridis_c(name = "-log10(adj p)") +
  scale_size_manual(
    name   = "|z|-Kategorie",
    values = c("1" = 1,   # ≤1
               "2" = 3,   # 1–2
               "3" = 7,   # 2–3
               "4" = 12), # >3
    breaks = c("1","2","3","4"),
    labels = c("≤1","1–2","2–3",">3"),
    drop   = FALSE
  ) +
  labs(
    x = "Fibroblast-Cluster", y = NULL,
    title = "IPA grouped pathways per cluster (repräsentativer Pfad je Gruppe)"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold"),
    plot.title  = element_text(face = "bold")
  )

p

# ---- Speichern ----
ggsave(outfile_png, p, width = 9, height = 6.2, dpi = 300)
write.csv(ipa_sum %>% arrange(group, as.numeric(as.character(cluster))), outfile_csv, row.names = FALSE)
write.csv(ipa_map %>% arrange(as.numeric(as.character(cluster)), group, padj), outfile_map, row.names = FALSE)

message("Gespeichert:\n  ", outfile_png, "\n  ", outfile_csv, "\n  ", outfile_map)



DefaultAssay(fib_subset) <- "SCT"
library(Seurat)
library(ggplot2)
# you can plot raw counts as well
vlnplt <- VlnPlot(fib_subset, 
                  features = c("Ccl5", "Cd52", "Cxcl13"), 
                  log = FALSE,
                  group.by = "condition",
                  layer = "scale.data")



out_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures"
ggsave(filename = file.path(out_dir, "violin_plot.png"),
       plot = vlnplt,
       width = 10, height = 8, dpi = 300)


View(fib_subset@meta.data)
View(fib_subset)

Assays(fib_subset)
VlnPlot(fib_subset, 
        features = c("Cd52", "Ccl5", "Cxcl13"),
        log = FALSE,
        assay = "RNA",
        group.by = "condition",
        layer = "data",
        pt.size = 0.1)
VlnPlot(
  fib_subset,
  features = c("Cxcl13"),
  assay   = "RNA",
  layer   = "data",
  group.by = "condition",
  pt.size  = 0.01
)

library(patchwork)

library(patchwork)
library(Seurat)
library(ggplot2)

p1 <- VlnPlot(
  fib_subset, features = "Cd52", assay = "RNA", layer = "data",
  group.by = "condition", pt.size = 0.01
) + ylim(0, 7) + NoLegend()

p2 <- VlnPlot(
  fib_subset, features = "Ccl5", assay = "RNA", layer = "data",
  group.by = "condition", pt.size = 0.01
) + ylim(0, 7) + NoLegend()

p3 <- VlnPlot(
  fib_subset, features = "Hcst", assay = "RNA", layer = "data",
  group.by = "condition", pt.size = 0.01
) + ylim(0, 7) + NoLegend()

# combine
combined <- p2 | p1

# save paths
out_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure3"
out_png <- file.path(out_dir, "violin_Cd52_Ccl5V4.png")
out_svg <- file.path(out_dir, "violin_Cd52_Ccl5V4.svg")

# save PNG
ggsave(out_png, plot = combined, width = 4.5, height = 3, units = "in", dpi = 300)

# save SVG
ggsave(out_svg, plot = combined, width = 4.5, height = 3, units = "in")









View(fib_subset@meta.data)
library(Seurat)
library(dplyr)

# pull expression for genes of interest + cluster IDs
dat <- FetchData(
  fib_subset,
  vars = c("Ccl5", "Cd52", "seurat_clusters")
)

# define whether a gene is expressed (>0 counts in RNA/data layer)
dat <- dat %>%
  mutate(
    Ccl5_pos  = Ccl5  > 0,
    Cd52_pos  = Cd52  > 0,
    Both_pos  = Ccl5_pos & Cd52_pos
  )

# count totals per cluster
counts <- dat %>%
  group_by(seurat_clusters) %>%
  summarise(
    total_cells    = n(),
    Ccl5_positive  = sum(Ccl5_pos),
    Cd52_positive  = sum(Cd52_pos),
    Both_positive  = sum(Both_pos),
    .groups = "drop"
  )

counts
669/1225
775+849
# total number of cells in the object
total_cells <- ncol(fib_subset)
total_cells



































# === Packages ===
library(readxl)
library(pheatmap)

# --- helper: find first matching column name (case-insensitive) ---
first_match_col <- function(df, candidates) {
  cn <- tolower(colnames(df))
  hit <- which(cn %in% tolower(candidates))
  if (length(hit) == 0) return(NA_character_)
  colnames(df)[hit[1]]
}

# --- main function ---
preview_deseq2_heatmap <- function(
    genes,
    base_dir,
    fontsize_row = 10,
    cluster_rows = FALSE   # FALSE = keep input order, TRUE = cluster
) {
  # find cluster directories like "Cluster3"
  all_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  cl_dirs  <- all_dirs[grepl("^Cluster\\s*\\d+$", basename(all_dirs), ignore.case = TRUE)]
  
  cols <- colorRampPalette(c("navy","white","goldenrod"))(100)
  mats <- list()
  
  for (d in cl_dirs) {
    cl_base <- basename(d)
    cl_num  <- sub(".*?(\\d+)$", "\\1", cl_base)
    xls <- list.files(
      d,
      pattern = paste0("(?i)DESeq2.*Cluster[_ ]?", cl_num, ".*\\.(xlsx|xls)$"),
      full.names = TRUE
    )
    if (length(xls) == 0) next
    f <- xls[which.max(file.info(xls)$mtime)]
    
    res <- as.data.frame(read_xlsx(f))  # tibble → data.frame
    gene_col <- first_match_col(res, c("gene","gene_id","symbol","genes"))
    lfc_col  <- first_match_col(res, c("log2foldchange","logfc","lfc"))
    if (any(is.na(c(gene_col, lfc_col)))) next
    
    res$..gene <- as.character(res[[gene_col]])
    
    # build df with ALL requested genes
    df <- data.frame(
      gene   = genes,
      log2FC = NA_real_,
      stringsAsFactors = FALSE
    )
    
    for (g in genes) {
      idx <- which(toupper(res$..gene) == toupper(g))[1]  # first match if duplicate
      if (!is.na(idx)) {
        df[df$gene == g, "log2FC"] <- res[[lfc_col]][idx]
      }
    }
    
    # safe rownames (handles duplicates automatically)
    rownames(df) <- make.unique(df$gene)
    
    mat <- as.matrix(df["log2FC"])
    colnames(mat) <- gsub("Cluster", "", cl_base, ignore.case = TRUE)
    
    mats[[cl_base]] <- mat
  }
  
  if (length(mats) == 0) {
    message("No data found.")
    return(NULL)
  }
  
  # combine all clusters into one matrix
  M <- do.call(cbind, mats)
  lim <- max(abs(M), na.rm = TRUE); if (!is.finite(lim)) lim <- 1
  bk  <- seq(-lim, lim, length.out = 101)
  
  print("Final rownames in heatmap matrix:")
  print(rownames(M))
  
  pheatmap(
    M,
    color = cols, breaks = bk,
    cluster_rows = cluster_rows, cluster_cols = FALSE,
    labels_row = rownames(M),   # always display gene names
    fontsize_row = fontsize_row,
    na_col = "grey90",
    main = "DESeq2 log2FC (preview)"
  )
}

# ---------------- Example usage ----------------
genes_to_view <- c("Isg15","Ifit1","Ifit3","Cxcl9","Cxcl10","Ccl5","Cxcl13")

preview_deseq2_heatmap(
  genes   = genes_to_view,
  base_dir = "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure4",
  fontsize_row = 10,
  cluster_rows = FALSE   # FALSE = keep your order, TRUE = cluster
)








# === Packages ===
library(readxl)
library(pheatmap)
library(svglite)

# --- helper: find first matching column name (case-insensitive) ---
first_match_col <- function(df, candidates) {
  cn <- tolower(colnames(df))
  hit <- which(cn %in% tolower(candidates))
  if (length(hit) == 0) return(NA_character_)
  colnames(df)[hit[1]]
}

# --- function to build matrix for a given gene list ---
make_matrix_for_genes <- function(genes, base_dir) {
  all_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  cl_dirs  <- all_dirs[grepl("^Cluster\\s*\\d+$", basename(all_dirs), ignore.case = TRUE)]
  
  mats <- list()
  
  for (d in cl_dirs) {
    cl_base <- basename(d)
    cl_num  <- sub(".*?(\\d+)$", "\\1", cl_base)
    xls <- list.files(
      d,
      pattern = paste0("(?i)DESeq2.*Cluster[_ ]?", cl_num, ".*\\.(xlsx|xls)$"),
      full.names = TRUE
    )
    if (length(xls) == 0) next
    f <- xls[which.max(file.info(xls)$mtime)]
    
    res <- as.data.frame(read_xlsx(f))  # tibble → data.frame
    gene_col <- first_match_col(res, c("gene","gene_id","symbol","genes"))
    lfc_col  <- first_match_col(res, c("log2foldchange","logfc","lfc"))
    if (any(is.na(c(gene_col, lfc_col)))) next
    
    res$..gene <- as.character(res[[gene_col]])
    
    # build df with ALL requested genes
    df <- data.frame(
      gene   = genes,
      log2FC = NA_real_,
      stringsAsFactors = FALSE
    )
    for (g in genes) {
      idx <- which(toupper(res$..gene) == toupper(g))[1]
      if (!is.na(idx)) {
        df[df$gene == g, "log2FC"] <- res[[lfc_col]][idx]
      }
    }
    
    rownames(df) <- make.unique(df$gene)
    mat <- as.matrix(df["log2FC"])
    colnames(mat) <- gsub("Cluster", "", cl_base, ignore.case = TRUE)
    
    mats[[cl_base]] <- mat
  }
  
  if (length(mats) == 0) return(NULL)
  
  M <- do.call(cbind, mats)
  
  # drop rows or columns that are all NA
  M <- M[rowSums(is.na(M)) < ncol(M), , drop = FALSE]
  M <- M[, colSums(is.na(M)) < nrow(M), drop = FALSE]
  
  if (nrow(M) == 0 || ncol(M) == 0) return(NULL)
  
  return(M)
}

# --- function to plot and save heatmap ---
save_heatmap <- function(M, title, out_dir, fontsize_row = 9) {
  lim <- max(abs(M), na.rm = TRUE); if (!is.finite(lim)) lim <- 1
  bk  <- seq(-lim, lim, length.out = 101)
  cols <- colorRampPalette(c("navy","white","goldenrod"))(100)
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # PNG
  png(file.path(out_dir, paste0(title, ".png")),
      width = 9, height = 6, units = "in", res = 300)
  pheatmap(M,
           color = cols, breaks = bk,
           cluster_rows = TRUE, cluster_cols = TRUE,
           labels_row = rownames(M),
           fontsize_row = fontsize_row,
           na_col = "grey90",
           main = title)
  dev.off()
  
  # SVG
  svglite(file.path(out_dir, paste0(title, ".svg")),
          width = 9, height = 6)
  pheatmap(M,
           color = cols, breaks = bk,
           cluster_rows = TRUE, cluster_cols = TRUE,
           labels_row = rownames(M),
           fontsize_row = fontsize_row,
           na_col = "grey90",
           main = title)
  dev.off()
}

# === Define gene sets ===
gene_sets <- list(
  "1_Chemokines_Cytokines" = c("Ccl5","Ccl4","Cxcl9","Cxcl10","Cxcl11","Ccl2","Ccl7","Il6","Lif"),
  "2_IFN_JAK_STAT_ISGs"    = c("Isg15","Ifit1","Ifit2","Ifit3","Gbp2","Gbp3","Gbp5","Ly6e","Bst2",
                               "Rsad2","Mx1","Oas1a","Oas1b","Oasl1","Irf1"),
  "3_Antigen_Processing"   = c("H2-K1","H2-D1","B2m","Tap1","Tap2","Psmb8","Psmb9","Psmb10",
                               "Cd74","Ciita","H2-Aa","H2-Ab1","H2-Eb1","H2-T23"),
  "4_PRR_Inflammation"     = c("Ddx58","Ifih1","Zbp1","Tlr3","Tlr4","Nfkbia","Tnfaip3","Tnfaip2",
                               "Nfkbiz","Icam1","Vcam1"),
  "5_Immune_Interface"     = c("Cd52","Cd55","Cd47","Cd200","Lbp","H2-T23","Pdcd1lg2","Cd274","Fas"),
  "6_Fibroblast_Identity"  = c("Col1a1","Dcn","Pdgfra","Pi16","Fbn1","Col14a1","Ltbp2","Adamts5","Cd248")
)

# === Run all and save ===
base_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure4"
out_dir  <- file.path(base_dir, "Exploratory_Heatmaps")

for (set_name in names(gene_sets)) {
  genes <- gene_sets[[set_name]]
  M <- make_matrix_for_genes(genes, base_dir)
  if (!is.null(M)) {
    save_heatmap(M, paste0("Exploratory_", set_name), out_dir)
  } else {
    message("No usable data found for ", set_name)
  }
}











# === Packages ===
library(readxl)
library(pheatmap)
library(svglite)

# --- helper: find first matching column name (case-insensitive) ---
first_match_col <- function(df, candidates) {
  cn <- tolower(colnames(df))
  hit <- which(cn %in% tolower(candidates))
  if (length(hit) == 0) return(NA_character_)
  colnames(df)[hit[1]]
}

# --- function to build matrix for a given gene list ---
make_matrix_for_genes <- function(genes, base_dir) {
  all_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  cl_dirs  <- all_dirs[grepl("^Cluster\\s*\\d+$", basename(all_dirs), ignore.case = TRUE)]
  
  mats <- list()
  
  for (d in cl_dirs) {
    cl_base <- basename(d)
    cl_num  <- sub(".*?(\\d+)$", "\\1", cl_base)
    xls <- list.files(
      d,
      pattern = paste0("(?i)DESeq2.*Cluster[_ ]?", cl_num, ".*\\.(xlsx|xls)$"),
      full.names = TRUE
    )
    if (length(xls) == 0) next
    f <- xls[which.max(file.info(xls)$mtime)]
    
    res <- as.data.frame(read_xlsx(f))
    gene_col <- first_match_col(res, c("gene","gene_id","symbol","genes"))
    lfc_col  <- first_match_col(res, c("log2foldchange","logfc","lfc"))
    if (any(is.na(c(gene_col, lfc_col)))) next
    
    res$..gene <- as.character(res[[gene_col]])
    
    df <- data.frame(
      gene   = genes,
      log2FC = NA_real_,
      stringsAsFactors = FALSE
    )
    for (g in genes) {
      idx <- which(toupper(res$..gene) == toupper(g))[1]
      if (!is.na(idx)) {
        df[df$gene == g, "log2FC"] <- res[[lfc_col]][idx]
      }
    }
    
    rownames(df) <- make.unique(df$gene)
    mat <- as.matrix(df["log2FC"])
    colnames(mat) <- gsub("Cluster", "", cl_base, ignore.case = TRUE)
    mats[[cl_base]] <- mat
  }
  
  if (length(mats) == 0) return(NULL)
  
  M <- do.call(cbind, mats)
  M <- M[rowSums(is.na(M)) < ncol(M), , drop = FALSE]
  M <- M[, colSums(is.na(M)) < nrow(M), drop = FALSE]
  
  if (nrow(M) == 0 || ncol(M) == 0) return(NULL)
  
  return(M)
}

# --- define your categories ---
gene_groups <- list(
  "Chemokines"          = c("Cxcl9","Ccl5","Ccl4"),
  "ISGs"                = c("Rsad2","Ifit2","Ifit3","Isg14"),
  "Antigen_Presentation"= c("H2-Ab1","H2-Eb1","H2-Aa","Cd74","Ciita"),
  "PRRs"                = c("Zbp1","Tlr3"),
  "Immune_Checkpoints"  = c("Cd274","Fas","Cd200","Cd52"),
  "ECM_Identity"        = c("Fbn1","Pi16","Adamts5")
)

# flatten into one gene list, but keep group info
all_genes <- unlist(gene_groups, use.names = FALSE)
gene_to_group <- rep(names(gene_groups), lengths(gene_groups))
names(gene_to_group) <- all_genes

# === Run ===
base_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure4"
out_file <- file.path(base_dir, "Exploratory_CombinedHeatmap.png")

M <- make_matrix_for_genes(all_genes, base_dir)

if (!is.null(M)) {
  # order rows by group, keep your input order
  row_order <- all_genes[all_genes %in% rownames(M)]
  M <- M[row_order, , drop = FALSE]
  
  # color scale
  lim <- max(abs(M), na.rm = TRUE); if (!is.finite(lim)) lim <- 1
  bk  <- seq(-lim, lim, length.out = 101)
  cols <- colorRampPalette(c("navy","white","goldenrod"))(100)
  
  # add group annotation for rows
  row_anno <- data.frame(Category = gene_to_group[rownames(M)])
  rownames(row_anno) <- rownames(M)
  
  # save PNG
  png(out_file, width = 9, height = 7, units = "in", res = 300)
  pheatmap(M,
           color = cols, breaks = bk,
           cluster_rows = FALSE, cluster_cols = TRUE,
           labels_row = rownames(M),
           fontsize_row = 10,
           na_col = "grey90",
           annotation_row = row_anno,
           main = "Exploratory Combined Heatmap")
  dev.off()
  
  message("Saved heatmap: ", out_file)
} else {
  message("No usable data for these genes")
}






list.files("C:/Users/ostrobel/Indiana University/O365-Strobel_Sequencing_Team - Documents/General/Analysis/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/Fib_subset/Pathways_info",
           pattern = "reactome", full.names = TRUE)




####figure 5 attempts ---------------


# =========================
# Global fibroblast reprogramming (IPA -> Families) Fig 5A
# =========================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tibble)
  library(purrr)
})

theme_set(theme_minimal(base_size = 12))

# ---- Paths ----
infile  <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/IPA_analysis/all_fibs/All_pathway_export.xls"
out_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/IPA_analysis/Figure5A"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

outfile_png <- file.path(out_dir, "Fig5A_IPA_global_barplot.png")
outfile_svg <- file.path(out_dir, "Fig5A_IPA_global_barplot.svg")
outfile_csv <- file.path(out_dir, "Fig5A_IPA_global_summary.csv")
outfile_map <- file.path(out_dir, "Fig5A_IPA_global_mapping.csv")

# ---- Options ----
agg_method <- "median"  # "median" or "stouffer"

# ---- Read IPA file ----
read_ipa_file <- function(path_xls){
  stopifnot(file.exists(path_xls))
  df <- read_excel(path_xls, skip = 1)
  names(df) <- str_trim(names(df))
  df %>%
    transmute(
      pathway   = `Ingenuity Canonical Pathways`,
      logp      = suppressWarnings(as.numeric(`-log(p-value)`)),
      ratio     = suppressWarnings(as.numeric(Ratio)),
      z         = suppressWarnings(as.numeric(`z-score`)),
      molecules = Molecules,
      p         = 10^(-suppressWarnings(as.numeric(`-log(p-value)`))),
      padj      = p.adjust(10^(-suppressWarnings(as.numeric(`-log(p-value)`))), "BH")
    )
}

ipa_all <- read_ipa_file(infile)

# ---- Group definitions ----
groups_tbl <- tribble(
  ~group,                           ~pattern,
  
  # IFN axis
  "Interferon / JAK–STAT",          "(?i)interferon|\\bIFN\\b|STAT1|\\bSTAT\\b|\\bJAK\\b|ISG15|ISGylation|\\bOAS\\b|\\bPKR\\b|cGAS|STING",
  
  # Antigen/APC
  "Antigen presentation",           "(?i)Antigen Presentation|\\bMHC\\b|\\bHLA\\b|H[-–]?2|TAP\\d?|CIITA|NLRC5|class\\s*I|class\\s*II|Dendritic Cell Maturation|Tumor Microenvironment|Costimulation|B Cell Receptor|T Cell Receptor",
  
  # Pro-inflammatory
  "NF-κB / TNF / IL-1",             "(?i)NF.?KB|REL[AB]?|\\bTNF\\b|IL-?1\\b|Interleukin-?\\s*1\\b|\\bNIK\\b|Macrophage Classical|Macrophage Alternative|MSP-RON|TAK1",
  
  # Cytokines
  "Cytokines / Chemokines (IL-6)",  "(?i)chemokine|CXCL\\d+|CCL\\d+|IL-?6\\b|Interleukin-?\\s*6\\b|IL-?8\\b|Interleukin-?\\s*8\\b|IL-?10\\b|Interleukin-?\\s*10\\b|IL-?3\\b|Interleukin-?\\s*3\\b|IL-?4\\b|Interleukin-?\\s*4\\b|IL-?13\\b|Interleukin-?\\s*13\\b|IL-?17\\b|Interleukin-?\\s*17\\b|IL-?27\\b|Interleukin-?\\s*27\\b|IL-?12\\b|Interleukin-?\\s*12\\b|IL-?15\\b|Interleukin-?\\s*15\\b|Cytokine Storm|Hypercytokinemia|Cachexia",
  
  # PRRs
  "PRR / TLR sensors",              "(?i)pattern recognition|\\bPRR\\b|RIG|MDA5|IRF\\d*|TLR\\d*|Toll-?like\\s+Receptor|Cytosolic sensors|Coronavirus Pathogenesis|C-?type\\s+lectin|\\bCLR\\b|\\bNOD1\\b|\\bNOD2\\b|HMGB1|TREM1",
  
  # T/NK/Checkpoint
  "T/NK / Checkpoint",              "(?i)\\bTCR\\b|TH1|TH2|TH17|\\bNK\\b|Natural\\s+Killer|checkpoint|PD-?1|CTLA-?4|ICOS|T cell exhaustion|PKC.?θ|FCER[I1]|Fc\\s*epsilon\\s*receptor|IL-?2\\s*Expression|Interleukin-?\\s*2\\s*Expression|CD40|OX40|Immunoregulatory interactions|Innate and Adaptive|\\bNFAT\\b",
  
  # Complement/acute phase
  "Complement / Acute phase",       "(?i)complement|acute phase|LXR/RXR|Neutrophil degranulation|Hematoma Resolution",
  
  # Fibrosis/dev
  "TGFb/Fibrosis/SMAD",             "(?i)\\bTGF\\b|fibros|SMAD\\d*|\\bECM\\b|collagen|Wound Healing|RUNX1|RUNX2|RUNX3|\\bID1\\b|NOTCH4",
  
  # Angiogenesis, IGF
  "VEGF/Angiogenesis",              "(?i)angiogenesis|\\bVEGF\\b|endothelial|\\bCSF3\\b|SCF-?KIT|Thrombopoietin|Growth Hormone",
  "IGF/IGFBP",                      "(?i)\\bIGF\\b|IGFBP\\d*|Prolactin",
  
  # Stress/HSP
  "Stress/HSP/GR",                  "(?i)stress|HSP\\d*|heat shock|hypoxia|\\bHIF\\b|KEAP1|NFE2L2|oxidative|\\bROS\\b|metabolism|polyamines|folate|\\bNAD\\b|DHCR24|Protein Ubiquitination|Deubiquitination|Post-?translational\\s+protein\\s+phosphorylation|Neddylation|\\bFAT10\\b|BAG2",
  
  # Growth/dev cascades
  "MAPK/WNT/EGFR/ERK",              "(?i)MAPK\\d*|\\bERK\\b|\\bJNK\\b|\\bp38\\b|\\bWNT\\b|\\bEGFR\\b|\\bRAF\\b|beta-?catenin|G\\s*alpha|\\bROBO\\b|\\bPTEN\\b|\\bCDC42\\b|Hedgehog|RAR Activation|platelet\\s+cytosolic\\s*Ca2",
  
  # Others
  "Cell death",                     "(?i)apoptosis|pyroptosis|necroptosis|ferroptosis|Death Receptor|Immunogenic Cell Death|Autophagy|Microautophagy|Chaperone Mediated Autophagy",
  "Inflammation",                   "(?i)inflammasome|nitric oxide|\\biNOS\\b|Eicosanoid",
  "Autoimmune signatures",          "(?i)lupus|rheumatoid|psoriasis|autoimmune|Type I Diabetes|Multiple Sclerosis|Atherosclerosis",
  "Neurodegeneration / CNS",        "(?i)huntington|parkinson|alzheimer|Neuroinflammation|Amyloid",
  "Cell cycle / DNA repair",        "(?i)cell cycle|DNA damage|checkpoint|mitotic|\\bS Phase\\b|\\bG1\\b|\\bG2\\b|Replication|Synthesis of DNA|Pre-?Initiation|Metaphase|Anaphase",
  "Endocytosis / trafficking",      "(?i)endocytosis|clathrin|virus entry|ABC-?family|Tubby|O-?linked\\s+glycosylation|Heparan\\s+Sulfate|HOTAIR"
)

# ---- Family map ----
family_map <- tribble(
  ~group,                          ~family,
  "Interferon / JAK–STAT",         "IFN/ISG & JAK–STAT",
  "NF-κB / TNF / IL-1",            "Proinflammatory cytokines & chemokines",
  "Cytokines / Chemokines (IL-6)", "Proinflammatory cytokines & chemokines",
  "Complement / Acute phase",      "Proinflammatory cytokines & chemokines",
  "PRR / TLR sensors",             "PRR/TLR sensors",
  "Antigen presentation",          "Antigen presentation & T/NK crosstalk",
  "T/NK / Checkpoint",             "Antigen presentation & T/NK crosstalk",
  "TGFb/Fibrosis/SMAD",            "ECM & tissue remodeling",
  "VEGF/Angiogenesis",             "ECM & tissue remodeling",
  "MAPK/WNT/EGFR/ERK",             "Growth factor & developmental signaling",
  "IGF/IGFBP",                     "Growth factor & developmental signaling",
  "Stress/HSP/GR",                 "Stress & metabolism",
  "Endocytosis / trafficking",     "Stress & metabolism",
  "Cell cycle / DNA repair",       "Cell cycle & DNA damage",
  "Cell death",                    "Cell death",
  "Inflammation",                  "Inflammation",
  "Autoimmune signatures",         "Disease signatures",
  "Neurodegeneration / CNS",       "Disease signatures"
)

family_priority <- c(
  "IFN/ISG & JAK–STAT",
  "PRR/TLR sensors",
  "Proinflammatory cytokines & chemokines",
  "Antigen presentation & T/NK crosstalk",
  "ECM & tissue remodeling",
  "Growth factor & developmental signaling",
  "Stress & metabolism",
  "Cell cycle & DNA damage",
  "Cell death",
  "Inflammation",
  "Disease signatures"
)

group_priority <- c(
  "Interferon / JAK–STAT",
  "PRR / TLR sensors",
  "NF-κB / TNF / IL-1",
  "Cytokines / Chemokines (IL-6)",
  "Complement / Acute phase",
  "Antigen presentation",
  "T/NK / Checkpoint",
  "TGFb/Fibrosis/SMAD",
  "VEGF/Angiogenesis",
  "IGF/IGFBP",
  "MAPK/WNT/EGFR/ERK",
  "Stress/HSP/GR",
  "Cell cycle / DNA repair",
  "Cell death",
  "Inflammation",
  "Autoimmune signatures",
  "Neurodegeneration / CNS",
  "Endocytosis / trafficking"
)


# ---- Helpers ----
assign_groups_all <- function(pathway_name, groups_tbl) {
  if (is.na(pathway_name) || pathway_name == "") return(character(0))
  groups_tbl %>%
    filter(str_detect(pathway_name, pattern)) %>%
    pull(group) %>%
    unique()
}

assign_group_priority <- function(pathway_name, groups_tbl, priority) {
  hits <- assign_groups_all(pathway_name, groups_tbl)
  if (length(hits) == 0) return("Unassigned")
  hits[order(match(hits, priority))][1]
}

estimate_pathway_size <- function(molecules_str, ratio) {
  if (is.na(ratio) || ratio <= 0) return(NA_real_)
  overlap <- length(unlist(strsplit(molecules_str, "\\s*[,;]\\s*")))
  if (overlap == 0) return(NA_real_)
  overlap / ratio
}

stouffer_z <- function(z_vec, w_vec = NULL) {
  z_vec <- z_vec[is.finite(z_vec)]
  if (length(z_vec) == 0) return(NA_real_)
  if (is.null(w_vec)) w_vec <- rep(1, length(z_vec))
  sum(w_vec * z_vec) / sqrt(sum(w_vec^2))
}

# ---- Map & aggregate ----
ipa_map <- ipa_all %>%
  mutate(
    group_all = purrr::map(pathway, assign_groups_all, groups_tbl = groups_tbl),
    group     = purrr::map_chr(pathway, assign_group_priority, groups_tbl = groups_tbl, priority = group_priority),
    groups_all = purrr::map_chr(group_all, ~ ifelse(length(.x) == 0, "", paste(.x, collapse = "; ")))
  ) %>%
  select(-group_all) %>%
  left_join(family_map, by = "group") %>%
  mutate(
    family = ifelse(is.na(family), "Unassigned", family),
    family = factor(family, levels = c(family_priority, "Unassigned"))
  )


write.csv(ipa_map, outfile_map, row.names = FALSE)

ipa_sig <- ipa_map %>%
  filter(padj <= 0.05, !is.na(z), is.finite(z)) %>%
  rowwise() %>%
  mutate(
    pathway_size_est = estimate_pathway_size(molecules, ratio),
    w_stouffer       = ifelse(is.finite(pathway_size_est), sqrt(pathway_size_est), NA_real_)
  ) %>%
  ungroup()

if (agg_method == "median") {
  ipa_sum_family <- ipa_sig %>%
    group_by(family, .drop = FALSE) %>%
    summarise(z_agg = median(z, na.rm = TRUE), padj_min = min(padj, na.rm = TRUE), n_paths = n(), .groups = "drop")
} else {
  ipa_sum_family <- ipa_sig %>%
    group_by(family, .drop = FALSE) %>%
    summarise(z_agg = stouffer_z(z, w_stouffer), padj_min = min(padj, na.rm = TRUE), n_paths = n(), .groups = "drop")
}

ipa_sum_family <- ipa_sum_family %>%
  mutate(padj_plot = pmax(padj_min, .Machine$double.xmin),
         neglog10FDR = -log10(padj_plot)) %>%
  filter(!is.na(family) & family != "Unassigned") %>%
  mutate(family = factor(as.character(family), levels = family_priority)) %>%
  arrange(match(as.character(family), family_priority))

write.csv(ipa_sum_family, outfile_csv, row.names = FALSE)

# ---- Final bar plot (gray) ----
dd <- ipa_sum_family %>%
  arrange(desc(z_agg)) %>%
  mutate(family = factor(family, levels = rev(family)))

p <- ggplot(dd, aes(x = family, y = z_agg)) +
  geom_col(fill = "grey70", color = "black") +
  coord_flip() +
  labs(
    x = NULL, y = "Aggregated Z-score",
    title = "Global fibroblast reprogramming with infection",
    subtitle = "Families aggregated from IPA terms (padj ≤ 0.05)"
  ) +
  theme(
    axis.text.y = element_text(face = "bold"),
    plot.title  = element_text(face = "bold")
  )

ggsave(outfile_png, p, width = 6, height = 5, dpi = 300)
ggsave(outfile_svg, p, width = 6, height = 5, dpi = 300)

message("Saved outputs:\n  ", outfile_png,
        "\n  ", outfile_csv,
        "\n  ", outfile_map)







###now run this by cluster-------
# ============================
# Figure 5B: Fibroblast clusters (Fib0–Fib4) IPA -> Zoom-in categories

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tibble)
  library(purrr)
})

theme_set(theme_minimal(base_size = 12))

# ---- Paths ----
base_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/IPA_analysis"
clusters <- paste0("Fib", 0:4)   # Fib0–Fib4
out_dir  <- file.path(base_dir, "Figure5B")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

outfile_png   <- file.path(out_dir, "Fig5B_IPA_heatmap.png")
outfile_csv   <- file.path(out_dir, "Fig5B_clusters_summary.csv")
outfile_full  <- file.path(out_dir, "Fig5B_full_mapping.csv")   # legacy
outfile_perpw <- file.path(out_dir, "Fig5B_full_mapping_per_pathway.csv")  # NEW

# ---- Options ----
agg_method <- "median"  # "median" (robust) or "stouffer"

# ---- Read IPA file ----
read_ipa_file <- function(path_xls){
  stopifnot(file.exists(path_xls))
  df <- read_excel(path_xls, skip = 1)
  names(df) <- str_trim(names(df))
  must <- c("Ingenuity Canonical Pathways", "-log(p-value)", "Ratio", "z-score", "Molecules")
  if (!all(must %in% names(df))) {
    stop("Unexpected columns in: ", path_xls,
         "\nFound: ", paste(names(df), collapse=" | "),
         "\nExpected: ", paste(must, collapse=" | "))
  }
  df %>%
    transmute(
      pathway   = `Ingenuity Canonical Pathways`,
      logp      = suppressWarnings(as.numeric(`-log(p-value)`)),
      ratio     = suppressWarnings(as.numeric(Ratio)),
      z         = suppressWarnings(as.numeric(`z-score`)),
      molecules = Molecules,
      p         = 10^(-suppressWarnings(as.numeric(`-log(p-value)`))),
      padj      = p.adjust(10^(-suppressWarnings(as.numeric(`-log(p-value)`))), "BH")
    )
}

# ---- Group definitions (Figure 5B) ----
groups_tbl <- tribble(
  ~group,                           ~pattern,
  
  # IFN axis
  "Interferon / JAK–STAT",          "(?i)interferon|\\bIFN\\b|STAT1|\\bSTAT\\b|\\bJAK\\b|ISG15|ISGylation|\\bOAS\\b|\\bPKR\\b|cGAS|STING",
  
  # Antigen/APC
  "Antigen presentation",           "(?i)Antigen Presentation|\\bMHC\\b|\\bHLA\\b|H[-–]?2|TAP\\d?|CIITA|NLRC5|class\\s*I|class\\s*II|Dendritic Cell Maturation|Tumor Microenvironment|Costimulation|B Cell Receptor|T Cell Receptor",
  
  # Pro-inflammatory
  "NF-κB / TNF / IL-1",             "(?i)NF.?KB|REL[AB]?|\\bTNF\\b|IL-?1\\b|Interleukin-?\\s*1\\b|\\bNIK\\b|Macrophage Classical|Macrophage Alternative|MSP-RON|TAK1",
  
  # Cytokines
  "Cytokines / Chemokines (IL-6)",  "(?i)chemokine|CXCL\\d+|CCL\\d+|IL-?6\\b|Interleukin-?\\s*6\\b|IL-?8\\b|Interleukin-?\\s*8\\b|IL-?10\\b|Interleukin-?\\s*10\\b|IL-?3\\b|Interleukin-?\\s*3\\b|IL-?4\\b|Interleukin-?\\s*4\\b|IL-?13\\b|Interleukin-?\\s*13\\b|IL-?17\\b|Interleukin-?\\s*17\\b|IL-?27\\b|Interleukin-?\\s*27\\b|IL-?12\\b|Interleukin-?\\s*12\\b|IL-?15\\b|Interleukin-?\\s*15\\b|Cytokine Storm|Hypercytokinemia|Cachexia",
  
  # PRRs
  "PRR / TLR sensors",              "(?i)pattern recognition|\\bPRR\\b|RIG|MDA5|IRF\\d*|TLR\\d*|Toll-?like\\s+Receptor|Cytosolic sensors|Coronavirus Pathogenesis|C-?type\\s+lectin|\\bCLR\\b|\\bNOD1\\b|\\bNOD2\\b|HMGB1|TREM1",
  
  # T/NK/Checkpoint
  "T/NK / Checkpoint",              "(?i)\\bTCR\\b|TH1|TH2|TH17|\\bNK\\b|Natural\\s+Killer|checkpoint|PD-?1|CTLA-?4|ICOS|T cell exhaustion|PKC.?θ|FCER[I1]|Fc\\s*epsilon\\s*receptor|IL-?2\\s*Expression|Interleukin-?\\s*2\\s*Expression|CD40|OX40|Immunoregulatory interactions|Innate and Adaptive|\\bNFAT\\b",
  
  # Complement/acute phase
  "Complement / Acute phase",       "(?i)complement|acute phase|LXR/RXR|Neutrophil degranulation|Hematoma Resolution",
  
  # Fibrosis/dev
  "TGFb/Fibrosis/SMAD",             "(?i)\\bTGF\\b|fibros|SMAD\\d*|\\bECM\\b|collagen|Wound Healing|RUNX1|RUNX2|RUNX3|\\bID1\\b|NOTCH4",
  
  # Angiogenesis, IGF
  "VEGF/Angiogenesis",              "(?i)angiogenesis|\\bVEGF\\b|endothelial|\\bCSF3\\b|SCF-?KIT|Thrombopoietin|Growth Hormone",
  "IGF/IGFBP",                      "(?i)\\bIGF\\b|IGFBP\\d*|Prolactin",
  
  # Stress/HSP
  "Stress/HSP/GR",                  "(?i)stress|HSP\\d*|heat shock|hypoxia|\\bHIF\\b|KEAP1|NFE2L2|oxidative|\\bROS\\b|metabolism|polyamines|folate|\\bNAD\\b|DHCR24|Protein Ubiquitination|Deubiquitination|Post-?translational\\s+protein\\s+phosphorylation|Neddylation|\\bFAT10\\b|BAG2",
  
  # Growth/dev cascades
  "MAPK/WNT/EGFR/ERK",              "(?i)MAPK\\d*|\\bERK\\b|\\bJNK\\b|\\bp38\\b|\\bWNT\\b|\\bEGFR\\b|\\bRAF\\b|beta-?catenin|G\\s*alpha|\\bROBO\\b|\\bPTEN\\b|\\bCDC42\\b|Hedgehog|RAR Activation|platelet\\s+cytosolic\\s*Ca2",
  
  # Others
  "Cell death",                     "(?i)apoptosis|pyroptosis|necroptosis|ferroptosis|Death Receptor|Immunogenic Cell Death|Autophagy|Microautophagy|Chaperone Mediated Autophagy",
  "Inflammation",                   "(?i)inflammasome|nitric oxide|\\biNOS\\b|Eicosanoid",
  "Autoimmune signatures",          "(?i)lupus|rheumatoid|psoriasis|autoimmune|Type I Diabetes|Multiple Sclerosis|Atherosclerosis",
  "Neurodegeneration / CNS",        "(?i)huntington|parkinson|alzheimer|Neuroinflammation|Amyloid",
  "Cell cycle / DNA repair",        "(?i)cell cycle|DNA damage|checkpoint|mitotic|\\bS Phase\\b|\\bG1\\b|\\bG2\\b|Replication|Synthesis of DNA|Pre-?Initiation|Metaphase|Anaphase",
  "Endocytosis / trafficking",      "(?i)endocytosis|clathrin|virus entry|ABC-?family|Tubby|O-?linked\\s+glycosylation|Heparan\\s+Sulfate|HOTAIR"
)

# ---- Broad families (reuse from Fig 5A) ----
family_map <- tribble(
  ~group,                          ~family,
  "Interferon / JAK–STAT",         "IFN/ISG & JAK–STAT",
  "NF-κB / TNF / IL-1",            "Proinflammatory cytokines & chemokines",
  "Cytokines / Chemokines (IL-6)", "Proinflammatory cytokines & chemokines",
  "Complement / Acute phase",      "Proinflammatory cytokines & chemokines",
  "PRR / TLR sensors",             "PRR/TLR sensors",
  "Antigen presentation",          "Antigen presentation & T/NK crosstalk",
  "T/NK / Checkpoint",             "Antigen presentation & T/NK crosstalk",
  "TGFb/Fibrosis/SMAD",            "ECM & tissue remodeling",
  "VEGF/Angiogenesis",             "ECM & tissue remodeling",
  "MAPK/WNT/EGFR/ERK",             "Growth factor & developmental signaling",
  "IGF/IGFBP",                     "Growth factor & developmental signaling",
  "Stress/HSP/GR",                 "Stress & metabolism",
  "Endocytosis / trafficking",     "Stress & metabolism",
  "Cell cycle / DNA repair",       "Cell cycle & DNA damage",
  "Cell death",                    "Cell death",
  "Inflammation",                  "Inflammation",
  "Autoimmune signatures",         "Disease signatures",
  "Neurodegeneration / CNS",       "Disease signatures"
)

# ---- Zoom-in categories (Figure 5B) ----
zoom_map <- tribble(
  ~group,                          ~zoom_in_category,
  "Interferon / JAK–STAT",         "Interferon / JAK–STAT",
  "NF-κB / TNF / IL-1",            "NF-κB / TNF / IL-1",
  "Cytokines / Chemokines (IL-6)", "Cytokines / Chemokines (IL-6)",
  "Complement / Acute phase",      "Complement / Acute phase",
  "Antigen presentation",          "Antigen presentation",
  "T/NK / Checkpoint",             "T/NK / Checkpoint",
  "PRR / TLR sensors",             "PRR / TLR sensors"
)

zoom_priority <- c(
  "Interferon / JAK–STAT",
  "PRR / TLR sensors",
  "NF-κB / TNF / IL-1",
  "Cytokines / Chemokines (IL-6)",
  "Complement / Acute phase",
  "Antigen presentation",
  "T/NK / Checkpoint"
)


# ---- Aggregation helpers ----
estimate_pathway_size <- function(molecules_str, ratio){
  if (is.na(ratio) || ratio <= 0) return(NA_real_)
  overlap <- length(unlist(strsplit(molecules_str, "\\s*[,;]\\s*")))
  if (overlap == 0) return(NA_real_)
  overlap / ratio
}
stouffer_z <- function(z_vec, w_vec = NULL){
  z_vec <- z_vec[is.finite(z_vec)]
  if (length(z_vec) == 0) return(NA_real_)
  if (is.null(w_vec)) w_vec <- rep(1, length(z_vec))
  sum(w_vec * z_vec) / sqrt(sum(w_vec^2))
}

# ---- Loop over clusters ----
per_pathway <- list()
per_cluster <- list()

for (cl in clusters) {
  infile <- file.path(base_dir, cl, "all pathways.xls")
  ipa <- read_ipa_file(infile) %>%
    mutate(
      group = map_chr(pathway, function(pw){
        hit <- groups_tbl %>% filter(str_detect(pw, pattern)) %>% pull(group)
        if (length(hit) == 0) "Unassigned" else hit[1]
      })
    ) %>%
    left_join(family_map, by = "group") %>%
    left_join(zoom_map,   by = "group") %>%
    filter(padj <= 0.05, !is.na(z), is.finite(z)) %>%
    mutate(cluster = cl)
  
  # save unmatched for debugging
  ipa %>%
    filter(group == "Unassigned") %>%
    distinct(pathway) %>%
    write.csv(file.path(out_dir, paste0(cl, "_unmatched.csv")), row.names = FALSE)
  
  # save pathway-level detail
  per_pathway[[cl]] <- ipa
  
  # weights
  ipa <- ipa %>%
    rowwise() %>%
    mutate(pathway_size_est = estimate_pathway_size(molecules, ratio),
           w_stouffer = ifelse(is.finite(pathway_size_est), sqrt(pathway_size_est), NA_real_)) %>%
    ungroup()
  
  # aggregate by category
  if (agg_method=="median") {
    fam <- ipa %>% group_by(zoom_in_category) %>%
      summarise(z_agg = median(z, na.rm=TRUE), .groups="drop")
  } else {
    fam <- ipa %>% group_by(zoom_in_category) %>%
      summarise(z_agg = stouffer_z(z, w_stouffer), .groups="drop")
  }
  fam$cluster <- cl
  per_cluster[[cl]] <- fam
}

# combine
all_pathways <- bind_rows(per_pathway)   # pathway-level detail
all_results  <- bind_rows(per_cluster)   # aggregated summary

# save pathway-level detail
write.csv(all_pathways, outfile_perpw, row.names = FALSE)

# save aggregated full mapping (legacy name)
write.csv(all_results, outfile_full, row.names = FALSE)

# ---- Reorder & restrict to Fib3-driven categories ----
fib3_order <- all_results %>%
  filter(cluster == "Fib3", !is.na(zoom_in_category)) %>%
  arrange(desc(z_agg)) %>%
  pull(zoom_in_category) %>%
  unique()

fib3_order <- rev(fib3_order)

dd_fib3 <- all_results %>%
  filter(!is.na(zoom_in_category)) %>%
  filter(zoom_in_category %in% fib3_order) %>%
  mutate(
    zoom_in_category = factor(zoom_in_category, levels = fib3_order),
    cluster = factor(cluster, levels = clusters)
  )

write.csv(dd_fib3, outfile_csv, row.names = FALSE)

# ---- Plot heatmap ----
p <- ggplot(dd_fib3, aes(x = cluster, y = zoom_in_category, fill = z_agg)) +
  geom_tile(color = "grey70") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, name = "Aggregated z"
  ) +
  labs(
    x = NULL, y = NULL,
    title = "Fibroblast cluster reprogramming with infection",
    subtitle = "Zoom-in categories ordered by Fib3 aggregated z-score"
  ) +
  theme(
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    plot.title  = element_text(face = "bold")
  )



outfile_png <- file.path(out_dir, "Fig5B_IPA_zoom_dotplot.png")
outfile_svg <- file.path(out_dir, "Fig5B_IPA_zoom_dotplot.svg")
ggsave(outfile_png, p, width = 7, height = 6, dpi = 300)
ggsave(outfile_svg, p, width = 5, height = 4, device = "svg")


p


###Figure3C----
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# ---- Define modules ----
chemokine_genes <- c("Cxcl9", "Ccl5", "Ccl4", "Cxcl10")
mhcII_genes    <- c("H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74", "Ciita")

# ---- Add module scores ----
fib_subset <- AddModuleScore(fib_subset, features = list(chemokine_genes), name = "Chemokine")
fib_subset <- AddModuleScore(fib_subset, features = list(mhcII_genes),    name = "MHCII")

# Module score column names
chemokine_col <- "Chemokine1"
mhcII_col     <- "MHCII1"

# ---- Output path ----
fig_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/IPA_analysis/Figure5C"

# ---- Function for violin plots (no points) ----
make_violin <- function(obj, feature, title, filename){
  p <- VlnPlot(
    obj,
    features = feature,
    group.by = "condition",
    pt.size = 0
  ) +
    theme_minimal(base_size = 14) +
    labs(title = title, y = "Module Score") +
    theme(legend.position = "none")
  
  ggsave(file.path(fig_path, filename), p, width = 4, height = 4, dpi = 300)
  return(p)
}

# ---- All fibroblasts ----
p_all_chemokine <- make_violin(fib_subset, chemokine_col, "Chemokine module (all fibroblasts)", "allfib_chemokine.png")
p_all_mhcII     <- make_violin(fib_subset, mhcII_col, "MHC-II module (all fibroblasts)", "allfib_mhcII.png")

# ---- Cluster 3 only ----
fib_clust3 <- subset(fib_subset, seurat_clusters == "3")

p_cl3_chemokine <- make_violin(fib_clust3, chemokine_col, "Chemokine module (cluster 3)", "clust3_chemokine.png")
p_cl3_mhcII     <- make_violin(fib_clust3, mhcII_col, "MHC-II module (cluster 3)", "clust3_mhcII.png")

# ---- Optional preview ----
(p_all_chemokine | p_all_mhcII) / (p_cl3_chemokine | p_cl3_mhcII)



### ---- Quantitative summaries and per-mouse plots ----

library(dplyr)
library(ggplot2)

# --- column names created by AddModuleScore
ifng_col <- "IFNG_STAT11"
mhcII_col <- "MHCII1"

# --- identify mouse/animal column
mouse_col <- case_when(
  "mouse_id" %in% colnames(fib_subset@meta.data) ~ "mouse_id",
  "orig.ident" %in% colnames(fib_subset@meta.data) ~ "orig.ident",
  TRUE ~ "sample"
)

# --- make sure condition and cluster are present
meta <- fib_subset@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  select(cell, !!mouse_col, condition, seurat_clusters, !!ifng_col, !!mhcII_col)

# --- per-mouse median scores by cluster
scores_mouse <- meta %>%
  group_by(!!sym(mouse_col), condition, seurat_clusters) %>%
  summarise(
    IFNG_STAT1_median = median(.data[[ifng_col]], na.rm = TRUE),
    MHCII_median      = median(.data[[mhcII_col]], na.rm = TRUE),
    .groups = "drop"
  )

# --- per-cluster grand medians across mice (for ranking)
scores_cluster <- scores_mouse %>%
  group_by(seurat_clusters, condition) %>%
  summarise(
    IFNG_STAT1_overall = median(IFNG_STAT1_median),
    MHCII_overall      = median(MHCII_median),
    .groups = "drop"
  ) %>%
  mutate(Execution_Index = MHCII_overall - IFNG_STAT1_overall)

# --- export summaries
out_dir <- file.path(fig_path, "summary_scores")
dir.create(out_dir, showWarnings = FALSE)

write.csv(scores_mouse,   file.path(out_dir, "IFNG_MHCII_scores_per_mouse.csv"),   row.names = FALSE)
write.csv(scores_cluster, file.path(out_dir, "IFNG_MHCII_scores_per_cluster.csv"), row.names = FALSE)

# --- print concise summary to console
print(scores_cluster %>% arrange(desc(Execution_Index)))

# --- overlay per-mouse points on violins (All fibroblasts example)
make_violin_with_points <- function(obj, feature, title, filename, score_df, module){
  p <- VlnPlot(obj, features = feature, group.by = "condition", pt.size = 0) +
    geom_jitter(
      data = score_df %>%
        group_by(condition) %>%
        summarise(med = median(if (module == "ifng") IFNG_STAT1_median else MHCII_median)),
      aes(x = condition, y = med, color = condition),
      width = 0.1, size = 2
    ) +
    theme_minimal(base_size = 14) +
    labs(title = title, y = "Module Score") +
    theme(legend.position = "none")
  ggsave(file.path(fig_path, filename), p, width = 4, height = 4, dpi = 300)
  return(p)
}

# --- recreate All fibroblast violins with per-mouse points
p_all_ifng_mouse  <- make_violin_with_points(
  fib_subset, ifng_col, "IFN/STAT1 module (All fibroblasts, per mouse)", "allfib_ifng_stat1_permouse.png",
  scores_mouse %>% filter(seurat_clusters %in% unique(fib_subset$seurat_clusters)), "ifng"
)
p_all_mhcII_mouse <- make_violin_with_points(
  fib_subset, mhcII_col, "MHC-II module (All fibroblasts, per mouse)", "allfib_mhcII_permouse.png",
  scores_mouse %>% filter(seurat_clusters %in% unique(fib_subset$seurat_clusters)), "mhc"
)
(p_all_ifng_mouse | p_all_mhcII_mouse)

# --- per-cluster paired plots with per-mouse points (optional loop)
for (cl in 0:4) {
  obj_cl <- subset(fib_subset, seurat_clusters == as.character(cl))
  df_cl  <- filter(scores_mouse, seurat_clusters == as.character(cl))
  
  p1 <- VlnPlot(obj_cl, features = ifng_col, group.by = "condition", pt.size = 0) +
    geom_point(data = df_cl, aes(x = condition, y = IFNG_STAT1_median, color = condition),
               position = position_jitter(width = 0.15), size = 2) +
    theme_minimal(base_size = 14) +
    labs(title = paste("IFN/STAT1 module - Cluster", cl), y = "Per-mouse median score") +
    theme(legend.position = "none")
  
  p2 <- VlnPlot(obj_cl, features = mhcII_col, group.by = "condition", pt.size = 0) +
    geom_point(data = df_cl, aes(x = condition, y = MHCII_median, color = condition),
               position = position_jitter(width = 0.15), size = 2) +
    theme_minimal(base_size = 14) +
    labs(title = paste("MHC-II module - Cluster", cl), y = "Per-mouse median score") +
    theme(legend.position = "none")
  
  combo <- p1 | p2
  ggsave(file.path(fig_path, paste0("cluster", cl, "_ifng_mhcII_permouse.png")), combo,
         width = 5, height = 3, units = "in", dpi = 300)
  ggsave(file.path(fig_path, paste0("cluster", cl, "_ifng_mhcII_permouse.svg")), combo,
         width = 5, height = 3, units = "in")
}

# Δscore = infected - control median per cluster
delta_tbl <- scores_mouse %>%
  tidyr::pivot_wider(
    id_cols = c(seurat_clusters),
    names_from = condition,
    values_from = c(IFNG_STAT1_median, MHCII_median)
  ) %>%
  mutate(
    dIFNG  = IFNG_STAT1_median_Experimental - IFNG_STAT1_median_Control,
    dMHCII = MHCII_median_Experimental - MHCII_median_Control
  )

print(delta_tbl %>% arrange(desc(dIFNG)))














### Figure 3C — IFN/STAT1 vs MHC-II module analysis ----
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# ---- 1. Define gene modules ----
ifng_stat1_core <- c(
  "Stat1","Irf1","Ifit1","Ifit2","Ifit3",
  "Gbp2","Gbp3","Gbp5","Cxcl9","Cxcl10",
  "Socs1","Bst2","Oas1a","Oas2"
)

mhcII_genes <- c(
  "H2-Aa", "H2-Ab1", "Cd74",
  "H2-DMa", "H2-DMb1",
  "Ciita",
  "Ctss", "Ctsl", "Lgmn"
)

# ---- 2. Add module scores (SCT-normalized) ----
fib_subset <- AddModuleScore(fib_subset, features = list(ifng_stat1_core), name = "IFNG_STAT1")
fib_subset <- AddModuleScore(fib_subset, features = list(mhcII_genes),    name = "MHCII")

ifng_col <- "IFNG_STAT11"
mhcII_col <- "MHCII1"

# ---- 3. Set output path ----
fig_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/IPA_analysis/Figure5C"
out_dir  <- file.path(fig_path, "summary_scores")
dir.create(out_dir, showWarnings = FALSE)

# ---- 4. Identify mouse/animal column ----
mouse_col <- case_when(
  "mouse_id" %in% colnames(fib_subset@meta.data) ~ "mouse_id",
  "orig.ident" %in% colnames(fib_subset@meta.data) ~ "orig.ident",
  TRUE ~ "sample"
)

# ---- 5. Build meta table ----
meta <- fib_subset@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  select(cell, !!mouse_col, condition, seurat_clusters, !!ifng_col, !!mhcII_col)

# ---- 6. Per-mouse medians by cluster ----
scores_mouse <- meta %>%
  group_by(!!sym(mouse_col), condition, seurat_clusters) %>%
  summarise(
    IFNG_STAT1_median = median(.data[[ifng_col]], na.rm = TRUE),
    MHCII_median      = median(.data[[mhcII_col]], na.rm = TRUE),
    .groups = "drop"
  )

# ---- 7. Cluster-level medians & execution index ----
scores_cluster <- scores_mouse %>%
  group_by(seurat_clusters, condition) %>%
  summarise(
    IFNG_STAT1_overall = median(IFNG_STAT1_median),
    MHCII_overall      = median(MHCII_median),
    .groups = "drop"
  ) %>%
  mutate(Execution_Index = MHCII_overall - IFNG_STAT1_overall)

# ---- 8. Export numeric summaries ----
write.csv(scores_mouse,   file.path(out_dir, "IFNG_MHCII_scores_per_mouse.csv"),   row.names = FALSE)
write.csv(scores_cluster, file.path(out_dir, "IFNG_MHCII_scores_per_cluster.csv"), row.names = FALSE)
print(scores_cluster %>% arrange(desc(Execution_Index)))

# ---- 9. Δscore = Experimental − Control per cluster (unpaired comparison) ----
# Compute condition-level means and medians per cluster
delta_summary <- scores_mouse %>%
  group_by(seurat_clusters, condition) %>%
  summarise(
    mean_IFNG  = mean(IFNG_STAT1_median, na.rm = TRUE),
    mean_MHCII = mean(MHCII_median, na.rm = TRUE),
    median_IFNG  = median(IFNG_STAT1_median, na.rm = TRUE),
    median_MHCII = median(MHCII_median, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # pivot wider so we have Control vs Experimental columns
  pivot_wider(
    id_cols = seurat_clusters,
    names_from = condition,
    values_from = c(mean_IFNG, mean_MHCII, median_IFNG, median_MHCII)
  ) %>%
  # calculate delta (Experimental − Control)
  mutate(
    dIFNG_mean    = mean_IFNG_Experimental  - mean_IFNG_Control,
    dMHCII_mean   = mean_MHCII_Experimental - mean_MHCII_Control,
    dIFNG_median  = median_IFNG_Experimental  - median_IFNG_Control,
    dMHCII_median = median_MHCII_Experimental - median_MHCII_Control
  )

# ---- Export numeric results ----
write.csv(delta_summary, file.path(out_dir, "Delta_scores_per_cluster.csv"), row.names = FALSE)

# ---- Console summary ----
cat("\nΔ IFN/STAT1 (Experimental − Control) per cluster:\n")
print(delta_summary %>% arrange(desc(dIFNG_mean)))

# ---- Optional visualization: bar plots of Δscores ----
# IFN/STAT1
p_delta_ifng <- ggplot(delta_summary, aes(x = seurat_clusters, y = dIFNG_mean)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = round(dIFNG_mean, 3)), vjust = -0.5, size = 3) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Δ IFN/STAT1 module (Experimental − Control)",
    x = "Cluster", y = "Mean Δ module score"
  )

# MHC-II
p_delta_mhcII <- ggplot(delta_summary, aes(x = seurat_clusters, y = dMHCII_mean)) +
  geom_col(fill = "darkorange") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = round(dMHCII_mean, 3)), vjust = -0.5, size = 3) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Δ MHC-II module (Experimental − Control)",
    x = "Cluster", y = "Mean Δ module score"
  )

# Save plots
ggsave(file.path(fig_path, "delta_IFNG_STAT1_percluster.png"), p_delta_ifng, width = 5, height = 4, dpi = 300)
ggsave(file.path(fig_path, "delta_MHCII_percluster.png"), p_delta_mhcII, width = 5, height = 4, dpi = 300)

# Optional side-by-side preview
p_delta_ifng | p_delta_mhcII

# ---- 10. Functions for visualization ----

## (a) Plain violin plot
make_violin <- function(obj, feature, title, filename) {
  p <- VlnPlot(obj, features = feature, group.by = "condition", pt.size = 0) +
    theme_minimal(base_size = 14) +
    labs(title = title, y = "Module Score") +
    theme(legend.position = "none")
  ggsave(file.path(fig_path, filename), p, width = 4, height = 4, dpi = 300)
  p
}

## (b) Violin + per-mouse points
make_violin_with_points <- function(obj, feature, title, filename, score_df, module){
  # choose correct y variable
  yvar <- if (module == "ifng") "IFNG_STAT1_median" else "MHCII_median"
  
  # base violin from Seurat
  p <- VlnPlot(obj, features = feature, group.by = "condition", pt.size = 0) +
    # add external points — explicitly set data and mapping
    ggplot2::geom_point(
      data = score_df,
      mapping = ggplot2::aes(x = condition, y = .data[[yvar]], color = condition),
      position = ggplot2::position_jitter(width = 0.15),
      size = 2,
      inherit.aes = FALSE           # <-- prevents ggplot from inheriting "fill=ident" etc.
    ) +
    theme_minimal(base_size = 14) +
    labs(title = title, y = "Per-mouse median module score") +
    theme(legend.position = "none")
  
  # save to file
  ggsave(file.path(fig_path, filename), p, width = 4, height = 4, dpi = 300)
  return(p)
}


# ---- 11. All fibroblasts ----
p_all_ifng  <- make_violin(fib_subset, ifng_col, "IFN/STAT1 module (All fibroblasts)", "allfib_ifng_stat1.png")
p_all_mhcII <- make_violin(fib_subset, mhcII_col, "MHC-II module (All fibroblasts)",   "allfib_mhcII.png")
(p_all_ifng | p_all_mhcII)

# ---- 12. Overlay per-mouse medians on violins ----
p_all_ifng_mouse  <- make_violin_with_points(fib_subset, ifng_col,
                                             "IFN/STAT1 module (All fibroblasts, per mouse)", "allfib_ifng_stat1_permouse.png",
                                             scores_mouse, "ifng"
)
p_all_mhcII_mouse <- make_violin_with_points(fib_subset, mhcII_col,
                                             "MHC-II module (All fibroblasts, per mouse)", "allfib_mhcII_permouse.png",
                                             scores_mouse, "mhc"
)
(p_all_ifng_mouse | p_all_mhcII_mouse)

# ---- 13. Per-cluster violins with per-mouse points (clusters 0–4) ----
for (cl in 0:4) {
  obj_cl <- subset(fib_subset, seurat_clusters == as.character(cl))
  df_cl  <- filter(scores_mouse, seurat_clusters == as.character(cl))
  
  # IFN/STAT1 module
  p1 <- VlnPlot(obj_cl, features = ifng_col, group.by = "condition", pt.size = 0) +
    ggplot2::geom_point(
      data = df_cl,
      mapping = ggplot2::aes(x = condition, y = IFNG_STAT1_median, color = condition),
      position = ggplot2::position_jitter(width = 0.15),
      size = 2,
      inherit.aes = FALSE    # <-- prevents Seurat-internal aesthetics from carrying over
    ) +
    theme_minimal(base_size = 14) +
    labs(title = paste("IFN/STAT1 module - Cluster", cl),
         y = "Per-mouse median score") +
    theme(legend.position = "none")
  
  # MHC-II module
  p2 <- VlnPlot(obj_cl, features = mhcII_col, group.by = "condition", pt.size = 0) +
    ggplot2::geom_point(
      data = df_cl,
      mapping = ggplot2::aes(x = condition, y = MHCII_median, color = condition),
      position = ggplot2::position_jitter(width = 0.15),
      size = 2,
      inherit.aes = FALSE
    ) +
    theme_minimal(base_size = 14) +
    labs(title = paste("MHC-II module - Cluster", cl),
         y = "Per-mouse median score") +
    theme(legend.position = "none")
  
  combo <- p1 | p2
  
  # Save plots
  ggsave(file.path(fig_path, paste0("cluster", cl, "_ifng_mhcII_permouse.png")),
         combo, width = 5, height = 3, units = "in", dpi = 300)
  ggsave(file.path(fig_path, paste0("cluster", cl, "_ifng_mhcII_permouse.svg")),
         combo, width = 5, height = 3, units = "in")
}


# ---- 14. ΔIFN/STAT1 summary bar plots (use unpaired delta_summary) ----

# IFN/STAT1 Δ plot
p_delta_ifng <- ggplot(delta_summary, aes(x = seurat_clusters, y = dIFNG_mean)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = round(dIFNG_mean, 3)), vjust = -0.5, size = 3) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Δ IFN/STAT1 module (Experimental − Control)",
    x = "Cluster",
    y = "Mean Δ IFN/STAT1 module score"
  )

# MHC-II Δ plot
p_delta_mhcII <- ggplot(delta_summary, aes(x = seurat_clusters, y = dMHCII_mean)) +
  geom_col(fill = "darkorange") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = round(dMHCII_mean, 3)), vjust = -0.5, size = 3) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Δ MHC-II module (Experimental − Control)",
    x = "Cluster",
    y = "Mean Δ MHC-II module score"
  )

# save both
ggsave(file.path(fig_path, "delta_IFNG_STAT1_percluster.png"), p_delta_ifng, width = 5, height = 4, dpi = 300)
ggsave(file.path(fig_path, "delta_MHCII_percluster.png"), p_delta_mhcII, width = 5, height = 4, dpi = 300)

# optional side-by-side preview
p_delta_ifng | p_delta_mhcII

ggplot(delta_summary, aes(x=seurat_clusters, y=dIFNG_mean)) +
  geom_col(fill="steelblue") +
  geom_hline(yintercept=0, lty="dashed") +
  theme_minimal(base_size=13) +
  labs(title="Δ IFN/STAT1 (Experimental − Control)", x="Cluster", y="Mean Δ module score")





# ---- Figure 5C: responder fractions (per animal), then Δ across conditions ----
library(dplyr)
library(tidyr)
library(ggplot2)

# 1) define "high" using control-derived thresholds (robust to size)
# threshold = Control median + 1*mad per CLUSTER (computed within each cluster)
ctrl_stats <- fib_subset@meta.data %>%
  mutate(IFNG = IFNG_STAT11, MHCII = MHCII1) %>%
  group_by(seurat_clusters) %>%
  summarise(
    IFNG_thr  = median(IFNG[condition=="Control"], na.rm=TRUE) + mad(IFNG[condition=="Control"], na.rm=TRUE),
    MHCII_thr = median(MHCII[condition=="Control"], na.rm=TRUE) + mad(MHCII[condition=="Control"], na.rm=TRUE),
    .groups="drop"
  )

md <- fib_subset@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  mutate(IFNG = IFNG_STAT11, MHCII = MHCII1) %>%
  left_join(ctrl_stats, by="seurat_clusters") %>%
  mutate(
    IFNG_high  = IFNG  > IFNG_thr,
    MHCII_high = MHCII > MHCII_thr,
    DP_high    = IFNG_high & MHCII_high  # double-positive executional responders
  )

# 2) per-animal fractions within each cluster
mouse_col <- if ("mouse_id" %in% colnames(md)) "mouse_id" else if ("orig.ident" %in% colnames(md)) "orig.ident" else "sample"

frac_mouse <- md %>%
  group_by(!!sym(mouse_col), condition, seurat_clusters) %>%
  summarise(
    frac_IFNG_high = mean(IFNG_high,  na.rm=TRUE),
    frac_DP_high   = mean(DP_high,    na.rm=TRUE),
    n_cells        = dplyr::n(),
    .groups="drop"
  )

# 3) compute Δ fractions (Experimental − Control) per cluster
delta_frac <- frac_mouse %>%
  pivot_wider(
    id_cols = c(!!sym(mouse_col), seurat_clusters),
    names_from = condition,
    values_from = c(frac_IFNG_high, frac_DP_high)
  ) %>%
  mutate(
    d_frac_IFNG = frac_IFNG_high_Experimental - frac_IFNG_high_Control,
    d_frac_DP   = frac_DP_high_Experimental   - frac_DP_high_Control
  )

# 4) summarise across animals for plotting (mean ± sem if desired)
sum_frac <- delta_frac %>%
  group_by(seurat_clusters) %>%
  summarise(
    mean_d_frac_IFNG = mean(d_frac_IFNG, na.rm=TRUE),
    mean_d_frac_DP   = mean(d_frac_DP,   na.rm=TRUE),
    .groups="drop"
  )

# 5) plots for Figure 5C
p5c_left <- ggplot(sum_frac, aes(x=seurat_clusters, y=mean_d_frac_IFNG)) +
  geom_col(fill="steelblue") +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(x="Cluster", y="Δ fraction IFN-high (Exp − Ctrl)", title="IFN/STAT1 responders") +
  theme_minimal(base_size=13)

p5c_right <- ggplot(sum_frac, aes(x=seurat_clusters, y=mean_d_frac_DP)) +
  geom_col(fill="firebrick3") +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(x="Cluster", y="Δ fraction IFN+MHC-II double-positive", title="IFN + MHC-II executors") +
  theme_minimal(base_size=13)

ggsave(file.path(fig_path, "Fig5C_IFNG_responderFraction.png"), p5c_left,  width=5, height=4, dpi=300)
ggsave(file.path(fig_path, "Fig5C_DoublePositive_responderFraction.png"), p5c_right, width=5, height=4, dpi=300)

# (optional) export the per-mouse table for supplement
out_dir <- file.path(fig_path, "summary_scores")
dir.create(out_dir, showWarnings=FALSE)
write.csv(delta_frac, file.path(out_dir, "PerMouse_DeltaResponderFractions.csv"), row.names=FALSE)






library(Seurat)
library(ggplot2)

# ---- Updated map ----
celltype_map <- list(
  Immune     = c(10, 6, 7, 2, 19, 21, 24, 11, 17),
  Stromal    = c(12, 20, 23, 1, 5, 6, 22, 16, 0),
  Epithelial = c(8, 18, 3, 14, 4, 15, 9),
  Seminal    = c(13)
)

# assign
seurat_obj$celltype <- NA_character_
for (ct in names(celltype_map)) {
  seurat_obj$celltype[seurat_obj$seurat_clusters %in% celltype_map[[ct]]] <- ct
}

# convert to factor
seurat_obj$celltype <- factor(
  seurat_obj$celltype,
  levels = c("Immune", "Stromal", "Epithelial", "Seminal")
)

# sanity check
table(seurat_obj$seurat_clusters, seurat_obj$celltype)

# ---- DotPlot ----
DotPlot(seurat_obj,
        features = c("Ccr5", "Cxcr3"),
        group.by = "celltype",
        assay = "RNA") +
  RotatedAxis() +
  ggtitle("CCR5 / CXCR3 expression across cell types")





DotPlot(seurat_obj,
        features = c("Ccr5", "Cxcr3"),
        group.by = "celltype",
        split.by = "condition",
        assay = "RNA") +
  RotatedAxis() +
  theme(legend.position = "right")




library(tidyr)
library(ggplot2)

# Generate DotPlot and extract data
dp <- DotPlot(seurat_obj,
              features = c("Ccr5", "Cxcr3"),
              group.by = "celltype",
              split.by = "condition",
              assay = "RNA")

plot_data <- dp$data %>%
  separate(id, into = c("celltype", "condition"), sep = "_")

plot_data$celltype <- factor(plot_data$celltype,
                             levels = c("Immune", "Stromal", "Epithelial", "Seminal"))

# Rebuild plot with gradient and legend
plotccrcxcr <- ggplot(plot_data,
               aes(x = features.plot,
               y = celltype,
               size = pct.exp,
              color = avg.exp.scaled)) +
  geom_point() +
  facet_wrap(~condition) +
  scale_size(range = c(0, 8), name = "Percent Expressed") +
  scale_color_gradient(low = "lightgrey", high = "darkblue", name = "Scaled Expression") +
  RotatedAxis() +
  ggtitle("CCR5 / CXCR3 expression across cell types") +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")



# Define save path
save_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure5"

# PNG
ggsave(filename = "CCR5_CXCR3_dotplot.png",
       plot = plotccrcxcr,
       path = save_path,
       width = 6, height = 4, dpi = 300)

# SVG
ggsave(filename = "CCR5_CXCR3_dotplot.svg",
       plot = plotccrcxcr,
       path = save_path,
       width = 6, height = 4)





##Figure 2C----
# === Libraries ===
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

# === Zentrales Theme für alle Fonts ===
my_theme <- theme_minimal(base_size = 10, base_family = "Arial") +
  theme(
    # Fonts
    plot.title      = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title.x    = element_text(size = 10, face = "plain"),
    axis.title.y    = element_text(size = 10, face = "plain"),
    axis.text.x     = element_text(size = 10, face = "plain"),
    axis.text.y     = element_text(size = 10, face = "plain"),
    legend.title    = element_text(size = 10, face = "plain"),
    legend.text     = element_text(size = 10, face = "plain"),
    strip.text      = element_text(size = 10, face = "bold"),
    
    # Keine Linien
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_blank(),
    panel.border     = element_blank()
  )



# === Pfade ===
input_file <- "C:/Users/ostrobel/Indiana University/O365-IN-Phtx-Jerde Lab Main Storage - Strobel/Imaging/AnalysisV2/Full_analysis/Cd52_images_analyzed2/CD52_summary_per_image.xlsx"
output_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure2"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# === Daten laden ===
df <- read_excel(input_file, sheet = "Summed_values") %>%
  rename_with(~trimws(.x)) %>%
  mutate(
    Group = trimws(as.character(Group)),
    DPI   = as.integer(DPI)
  )

value_col <- "fib/total"

# === globale y-Achse berechnen ===
ymin <- 0
ymax <- max(df[[value_col]], na.rm = TRUE) * 1.25

# === Plotfunktion ===
plot_dpi <- function(df_in, dpi_value, ylim) {
  sub <- df_in %>% filter(DPI == dpi_value)
  
  control_vals  <- sub %>% filter(Group == "Control") %>% pull(value_col)
  infected_vals <- sub %>% filter(Group == "Infected") %>% pull(value_col)
  
  means <- c(mean(control_vals, na.rm = TRUE),
             mean(infected_vals, na.rm = TRUE))
  sems  <- c(sd(control_vals, na.rm = TRUE)/sqrt(length(control_vals)),
             sd(infected_vals, na.rm = TRUE)/sqrt(length(infected_vals)))
  
  # Welch t-test
  t_res <- t.test(control_vals, infected_vals, var.equal = FALSE)
  p_val <- t_res$p.value
  
  # Daten für Balken
  bar_df <- data.frame(
    Group = c("Control", "Infected"),
    Mean  = means,
    SEM   = sems
  )
  
  # Punkte mit Jitter
  point_df <- sub %>%
    mutate(
      Group = factor(Group, levels = c("Control", "Infected")),
      xjitter = as.numeric(Group) + rnorm(n(), sd = 0.05)
    )
  
  
  g <- ggplot(bar_df, aes(x = Group, y = Mean, fill = Group)) +
    geom_col(color = "black", alpha = 0.8) +
    geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2) +
    geom_point(data = point_df, aes(x = xjitter, y = !!sym(value_col)),
               inherit.aes = FALSE, size = 2, alpha = 0.7) +
    scale_fill_manual(values = c("skyblue", "salmon")) +
    labs(title = paste0(dpi_value, " DPI Fibroblasts"),
         y = "Fibroblasts / All Cells",
         x = "Group") +
    ylim(ymin, ylim) +
    annotate("text", x = 1.5, y = ylim * 0.95,
             label = paste0("p = ", signif(p_val, 3)), size = 4) +
    my_theme
  
  
  return(g)
}

# === Gesamtplot (2 nebeneinander) ===
p1 <- plot_dpi(df, 28, ymax)
p2 <- plot_dpi(df, 60, ymax)

combined <- plot_grid(p1, p2, ncol = 2, align = "hv")

combined

# === Speichern ===
out_svg <- file.path(output_dir, "fibroblast_ratio_28vs60DPI_bar.svg")
out_png <- file.path(output_dir, "fibroblast_ratio_28vs60DPI_bar.png")

ggsave(out_svg, combined, width = 6, height = 2.55)
ggsave(out_png, combined, width = 6, height = 2.55, dpi = 300)

message("✅ saved: ", out_svg)
message("✅ saved: ", out_png)





##Figure 4C----
# === Libraries ===
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

# === Zentrales Theme für alle Fonts ===
my_theme <- theme_minimal(base_size = 10, base_family = "Arial") +
  theme(
    # Fonts
    plot.title      = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title.x    = element_text(size = 10, face = "plain"),
    axis.title.y    = element_text(size = 10, face = "plain"),
    axis.text.x     = element_text(size = 10, face = "plain"),
    axis.text.y     = element_text(size = 10, face = "plain"),
    legend.title    = element_text(size = 10, face = "plain"),
    legend.text     = element_text(size = 10, face = "plain"),
    strip.text      = element_text(size = 10, face = "bold"),
    
    # Keine Linien
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_blank(),
    panel.border     = element_blank()
  )



# === Pfade ===
input_file <- "C:/Users/ostrobel/Indiana University/O365-IN-Phtx-Jerde Lab Main Storage - Strobel/Imaging/AnalysisV2/Full_analysis/Cd52_images_analyzed2/CD52_summary_per_image.xlsx"
output_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure4"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# === Daten laden ===
df <- read_excel(input_file, sheet = "Summed_values_renamed") %>%
  rename_with(~trimws(.x)) %>%
  mutate(
    Group = trimws(as.character(Group)),
    DPI   = as.integer(DPI)
  )

value_col <- "Cd52_fib"

# === globale y-Achse berechnen ===
ymin <- 0
ymax <- max(df[[value_col]], na.rm = TRUE) * 1.25

# === Plotfunktion ===
plot_dpi <- function(df_in, dpi_value, ylim) {
  sub <- df_in %>% filter(DPI == dpi_value)
  
  control_vals  <- sub %>% filter(Group == "Control") %>% pull(value_col)
  infected_vals <- sub %>% filter(Group == "Infected") %>% pull(value_col)
  
  means <- c(mean(control_vals, na.rm = TRUE),
             mean(infected_vals, na.rm = TRUE))
  sems  <- c(sd(control_vals, na.rm = TRUE)/sqrt(length(control_vals)),
             sd(infected_vals, na.rm = TRUE)/sqrt(length(infected_vals)))
  
  # Welch t-test
  t_res <- t.test(control_vals, infected_vals, var.equal = FALSE)
  p_val <- t_res$p.value
  
  # Daten für Balken
  bar_df <- data.frame(
    Group = c("Control", "Infected"),
    Mean  = means,
    SEM   = sems
  )
  
  # Punkte mit Jitter
  point_df <- sub %>%
    mutate(
      Group = factor(Group, levels = c("Control", "Infected")),
      xjitter = as.numeric(Group) + rnorm(n(), sd = 0.05)
    )
  
  
  g <- ggplot(bar_df, aes(x = Group, y = Mean, fill = Group)) +
    geom_col(color = "black", alpha = 0.8) +
    geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2) +
    geom_point(data = point_df, aes(x = xjitter, y = !!sym(value_col)),
               inherit.aes = FALSE, size = 2, alpha = 0.7) +
    scale_fill_manual(values = c("skyblue", "salmon")) +
    labs(title = paste0(dpi_value, " DPI Fibroblasts"),
         y = "Fibroblasts / All Cells",
         x = "Group") +
    ylim(ymin, ylim) +
    annotate("text", x = 1.5, y = ylim * 0.95,
             label = paste0("p = ", signif(p_val, 3)), size = 4) +
    my_theme
  
  
  return(g)
}

# === Gesamtplot (2 nebeneinander) ===
p1 <- plot_dpi(df, 28, ymax)
p2 <- plot_dpi(df, 60, ymax)

combined <- plot_grid(p1, p2, ncol = 2, align = "hv")

combined

# === Speichern ===
out_svg <- file.path(output_dir, "CD52_ratio_28vs60DPI_bar.svg")
out_png <- file.path(output_dir, "CD52_ratio_28vs60DPI_bar.png")

ggsave(out_svg, combined, width = 6, height = 2.55)
ggsave(out_png, combined, width = 6, height = 2.55, dpi = 300)

message("✅ saved: ", out_svg)
message("✅ saved: ", out_png)





# === Pfade ===
input_file <- "C:/Users/ostrobel/Indiana University/O365-IN-Phtx-Jerde Lab Main Storage - Strobel/Imaging/AnalysisV2/Full_analysis/Cxcl13_images_analyzed/Cxcl13_summary_per_image.xlsx"
output_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure4"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# === Daten laden ===
df <- read_excel(input_file, sheet = "Summed_values_renamed") %>%
  rename_with(~trimws(.x)) %>%
  mutate(
    Group = trimws(as.character(Group)),
    DPI   = as.integer(DPI)
  )

value_col <- "Cxcl13_fib"

# === globale y-Achse berechnen ===
ymin <- 0
ymax <- max(df[[value_col]], na.rm = TRUE) * 1.25

# === Plotfunktion ===
plot_dpi <- function(df_in, dpi_value, ylim) {
  sub <- df_in %>% filter(DPI == dpi_value)
  
  control_vals  <- sub %>% filter(Group == "Control") %>% pull(value_col)
  infected_vals <- sub %>% filter(Group == "Infected") %>% pull(value_col)
  
  means <- c(mean(control_vals, na.rm = TRUE),
             mean(infected_vals, na.rm = TRUE))
  sems  <- c(sd(control_vals, na.rm = TRUE)/sqrt(length(control_vals)),
             sd(infected_vals, na.rm = TRUE)/sqrt(length(infected_vals)))
  
  # Welch t-test
  t_res <- t.test(control_vals, infected_vals, var.equal = FALSE)
  p_val <- t_res$p.value
  
  # Daten für Balken
  bar_df <- data.frame(
    Group = c("Control", "Infected"),
    Mean  = means,
    SEM   = sems
  )
  
  # Punkte mit Jitter
  point_df <- sub %>%
    mutate(
      Group = factor(Group, levels = c("Control", "Infected")),
      xjitter = as.numeric(Group) + rnorm(n(), sd = 0.05)
    )
  
  
  g <- ggplot(bar_df, aes(x = Group, y = Mean, fill = Group)) +
    geom_col(color = "black", alpha = 0.8) +
    geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2) +
    geom_point(data = point_df, aes(x = xjitter, y = !!sym(value_col)),
               inherit.aes = FALSE, size = 2, alpha = 0.7) +
    scale_fill_manual(values = c("skyblue", "salmon")) +
    labs(title = paste0(dpi_value, " DPI Fibroblasts"),
         y = "Fibroblasts / All Cells",
         x = "Group") +
    ylim(ymin, ylim) +
    annotate("text", x = 1.5, y = ylim * 0.95,
             label = paste0("p = ", signif(p_val, 3)), size = 4) +
    my_theme
  
  
  return(g)
}

# === Gesamtplot (2 nebeneinander) ===
p1 <- plot_dpi(df, 28, ymax)
p2 <- plot_dpi(df, 60, ymax)

combined <- plot_grid(p1, p2, ncol = 2, align = "hv")

combined

# === Speichern ===
out_svg <- file.path(output_dir, "Cxcl13_ratio_28vs60DPI_bar.svg")
out_png <- file.path(output_dir, "Cxcl13_ratio_28vs60DPI_bar.png")

ggsave(out_svg, combined, width = 6, height = 2.55)
ggsave(out_png, combined, width = 6, height = 2.55, dpi = 300)

message("✅ saved: ", out_svg)
message("✅ saved: ", out_png)






library(msigdbr)

library(dplyr)

# All Hallmark gene sets for mouse
m_df <- msigdbr(species = "Mus musculus", category = "H")

# Browse available pathway names
unique(m_df$gs_name)


# pull all Hallmark pathways for mouse
hallmark <- msigdbr(species = "Mus musculus", category = "H")

# IFN α and γ response sets
ifn_alpha <- hallmark %>% filter(gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE")
ifn_gamma <- hallmark %>% filter(gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE")

# combine and deduplicate
ifn_genes <- unique(c(ifn_alpha$gene_symbol, ifn_gamma$gene_symbol))
ifn_genes <- intersect(ifn_genes, rownames(fib_subset))


ifn_extra <- hallmark %>% filter(gs_name == "HALLMARK_IL6_JAK_STAT3_SIGNALING")
ifn_genes <- unique(c(ifn_genes, ifn_extra$gene_symbol))





# Re-run the GO pull first
go_bp <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")

# Check what MHC-related sets exist
unique(go_bp$gs_name[grep("ANTIGEN|MHC", go_bp$gs_name, ignore.case = TRUE)])

# Select the one you want
mhcII_go <- go_bp %>%
  filter(gs_name == "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II")

# Fix capitalization mismatch
mhc_genes_raw <- unique(mhcII_go$gene_symbol)

# Make everything consistent
mhc_genes <- intersect(toupper(rownames(fib_subset)), toupper(mhc_genes_raw))

# Optional: print what overlapped
mhc_genes








library(readxl)
library(dplyr)
library(stringr)
library(writexl)
library(tidyr)

# ---- File paths ----
input_file  <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/DE_excel_docs_for_analysis/2024_10_07_contvsinfected_gene_exp.xlsx"
output_file <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/2025_10_08_Figure5/IFNg_fibroblast_genes_DE_hits.xlsx"

# ---- Target genes ----
target_genes <- c("Irf1","Irf5","Irf7","Stat1","Gbp2","Gbp5","Irgb6","Ciita","H2-Ab1","H2-DMa","Psmb8")

# ---- Identify cluster sheets ----
sheets <- excel_sheets(input_file)
cluster_sheets <- sheets[str_detect(sheets, regex("^Cluster\\d+$", ignore_case = TRUE))]

# ---- Read, clean, and filter ----
results <- lapply(cluster_sheets, function(sheet) {
  df <- read_excel(input_file, sheet = sheet, col_names = TRUE)
  df <- df[, 1:6]
  colnames(df) <- c("gene_id","log2FC.single","p_val_adj.single",
                    "log2FC.bulk","p_val_adj.bulk","cluster")
  df %>%
    filter(gene_id %in% target_genes) %>%
    mutate(sheet = sheet)
})

combined_df <- bind_rows(results)

# ---- Create wide format ----
wide_df <- combined_df %>%
  select(gene_id, cluster, log2FC.bulk, p_val_adj.bulk) %>%
  mutate(cluster = paste0("cluster", cluster)) %>%
  pivot_wider(
    names_from = cluster,
    values_from = c(log2FC.bulk, p_val_adj.bulk),
    names_glue = "{.value}_{cluster}"
  ) %>%
  arrange(gene_id)

# ---- Verify content ----
cat("\nRows in Combined:", nrow(combined_df),
    "\nRows in Wide_View:", nrow(wide_df), "\n")

# ---- Write both tabs ----
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write_xlsx(
  x = list(
    "Combined"  = combined_df,
    "Wide_View" = wide_df
  ),
  path = output_file
)

message("✅  Saved two-sheet Excel file to: ", output_file)





library(readxl)
library(dplyr)
library(stringr)
library(writexl)
library(tidyr)

# ---- File paths ----
input_file  <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/DE_excel_docs_for_analysis/2024_10_07_contvsinfected_gene_exp.xlsx"
output_file <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure5/2025_10_08_IFNg_fibroblast_genes_DE_hits_SC.xlsx"

# ---- Target genes ----
target_genes <- c("Irf1","Irf5","Irf7","Stat1","Gbp2","Gbp5","Irgb6","Ciita","H2-Ab1","H2-DMa","Psmb8")

# ---- Identify cluster sheets ----
sheets <- excel_sheets(input_file)
cluster_sheets <- sheets[str_detect(sheets, regex("^Cluster\\d+$", ignore_case = TRUE))]

# ---- Read, clean, and filter ----
results <- lapply(cluster_sheets, function(sheet) {
  df <- read_excel(input_file, sheet = sheet, col_names = TRUE)
  df <- df[, 1:6]
  colnames(df) <- c("gene_id","log2FC.single","p_val_adj.single",
                    "log2FC.bulk","p_val_adj.bulk","cluster")
  df %>%
    filter(gene_id %in% target_genes) %>%
    mutate(sheet = sheet)
})

combined_df <- bind_rows(results)

# ---- Create wide format using SINGLE-CELL columns ----
wide_df <- combined_df %>%
  select(gene_id, cluster, log2FC.single, p_val_adj.single) %>%
  mutate(cluster = paste0("cluster", cluster)) %>%
  pivot_wider(
    names_from = cluster,
    values_from = c(log2FC.single, p_val_adj.single),
    names_glue = "{.value}_{cluster}"
  ) %>%
  arrange(gene_id)

# ---- Verify content ----
cat("\nRows in Combined:", nrow(combined_df),
    "\nRows in Wide_View:", nrow(wide_df), "\n")

# ---- Write both tabs ----
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write_xlsx(
  x = list(
    "Combined"  = combined_df,
    "Wide_View" = wide_df
  ),
  path = output_file
)

message("✅  Saved single-cell IFN-γ gene hits to: ", output_file)












library(readxl)
library(dplyr)
library(stringr)
library(writexl)
library(tidyr)

# ---- File paths ----
input_file  <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/DE_excel_docs_for_analysis/2024_10_07_contvsinfected_gene_exp.xlsx"
output_file <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/2025_10_08_Figure5/Cluster3_signature_bulk_DE_hits.xlsx"

# ---- Target genes (cluster 3 upregulated list) ----
target_genes <- c(
  "Cxcl9","Cxcl10","Isg15","Ifit1","Gbp2","Ifi47","Iigp1","Ly6e",
  "H2-Ab1","Zbp1","H2-K1","Cd3d","Ifi27l2a","Hcst","H2-D1","Plac8"
)

# ---- Identify cluster sheets ----
sheets <- excel_sheets(input_file)
cluster_sheets <- sheets[str_detect(sheets, regex("^Cluster\\d+$", ignore_case = TRUE))]

# ---- Read, clean, and filter ----
results <- lapply(cluster_sheets, function(sheet) {
  df <- read_excel(input_file, sheet = sheet, col_names = TRUE)
  df <- df[, 1:6]
  colnames(df) <- c("gene_id","log2FC.single","p_val_adj.single",
                    "log2FC.bulk","p_val_adj.bulk","cluster")
  df %>%
    filter(gene_id %in% target_genes) %>%
    mutate(sheet = sheet)
})

combined_df <- bind_rows(results)

# ---- Create wide format (bulk columns) ----
wide_df <- combined_df %>%
  select(gene_id, cluster, log2FC.bulk, p_val_adj.bulk) %>%
  mutate(cluster = paste0("cluster", cluster)) %>%
  pivot_wider(
    names_from = cluster,
    values_from = c(log2FC.bulk, p_val_adj.bulk),
    names_glue = "{.value}_{cluster}"
  ) %>%
  arrange(gene_id)

# ---- Verify content ----
cat("\nRows in Combined:", nrow(combined_df),
    "\nRows in Wide_View:", nrow(wide_df), "\n")

# ---- Write both tabs ----
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write_xlsx(
  x = list(
    "Combined"  = combined_df,
    "Wide_View" = wide_df
  ),
  path = output_file
)

message("✅  Saved cluster 3 gene hits (bulk DE) to: ", output_file)




View(fib_subset@meta.data)





library(Seurat)
library(ggplot2)
library(patchwork)

# ---- Paths ----
output_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure5/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Gene list ----
genes_to_plot <- c("Stat1", "Gbp2", "Cxcl9", "Cxcl10", "Cxcl13", "Ccl5")

# ---- Helper function: make violin plots for one cluster ----
make_violin_cluster <- function(obj, cluster_id, genes, output_dir) {
  # Subset to cluster of interest
  subset_obj <- subset(obj, subset = seurat_clusters == cluster_id)
  
  # Extract max expression value across all genes (from full object, SCT assay, data layer)
  expr_values <- FetchData(obj, vars = genes, layer = "data", assay = "SCT")
  ymax <- max(expr_values)
  
  # Generate violin plots for each gene
  plots <- lapply(genes, function(gene) {
    VlnPlot(
      subset_obj,
      features = gene,
      group.by = "condition",
      pt.size = 0,
      assay = "SCT",
      layer = "data",
      cols = c("#9fd8ef", "#fb998e")  # blue = Control, red = Experimental
    ) +
      ggtitle(paste0(gene, " (Cluster ", cluster_id, ")")) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11)
      ) +
      ylim(0, ymax)
  })
  
  # Combine plots (2 rows x 3 columns)
  combined_plot <- wrap_plots(plots, ncol = 3)
  
  # Save figure
  png_filename <- file.path(output_dir, paste0("Cluster", cluster_id, "_IFN_genes_violin_SCT.png"))
  png(png_filename, width = 3200, height = 1800, res = 300)
  print(combined_plot)
  dev.off()
  
  message("✅ Saved: ", png_filename)
}

# ---- Run for clusters 2 and 3 ----
make_violin_cluster(fib_subset, cluster_id = 2, genes = genes_to_plot, output_dir = output_dir)
make_violin_cluster(fib_subset, cluster_id = 3, genes = genes_to_plot, output_dir = output_dir)









####actual figure code------
library(Seurat)
library(ggplot2)
library(patchwork)

# ---- Output-Ordner ----
output_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure5/Violin_plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Gene-Liste ----
genes_to_plot <- c("Stat1", "Gbp2", "Cxcl9", "Cxcl10", "Cxcl13", "Ccl5")

# ---- Hilfsfunktion: Violinplots pro Cluster ----
make_violin_cluster <- function(obj, cluster_id, genes, output_dir) {
  # Objekt auf Cluster beschränken
  subset_obj <- subset(obj, subset = seurat_clusters == cluster_id)
  
  # Plots generieren (automatische y-Skala)
  plots <- lapply(genes, function(gene) {
    VlnPlot(
      subset_obj,
      features = gene,
      group.by = "condition",
      pt.size = 0,
      assay = "SCT",
      layer = "data",
      cols = c("#9fd8ef", "#fb998e")  # blau = Control, rot = Experimental
    ) +
      ggtitle(paste0(gene, " (Cluster ", cluster_id, ")")) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9)
      )
  })
  
  # 2 Reihen x 3 Spalten kombinieren
  combined_plot <- wrap_plots(plots, ncol = 3)
  
  # Als SVG speichern (keine Unterordner)
  svg_filename <- file.path(output_dir, paste0("Cluster", cluster_id, "_IFN_genes_violin_SCT.svg"))
  svg(svg_filename, width = 10, height = 6)
  print(combined_plot)
  dev.off()
  
  message("✅ Saved: ", svg_filename)
}

# ---- Für Cluster 0–4 laufen lassen ----
for (cl in 0:4) {
  make_violin_cluster(fib_subset, cluster_id = cl, genes = genes_to_plot, output_dir = output_dir)
}




library(Seurat)
library(ggplot2)

# ---- Output-Ordner ----
output_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure5/Violin_plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Gene-Liste ----
genes_to_plot <- c("Stat1", "Gbp2", "Cxcl9", "Cxcl10", "Cxcl13", "Ccl5", "Gbp1")

# ---- Plot-Größe (in Zoll) ----
plot_width  <- 3
plot_height <- 3

# ---- Globales Maximum für Y-Achse (über alle Cluster und Gene) ----
expr_values_global <- FetchData(fib_subset, vars = genes_to_plot, layer = "data", assay = "SCT")
ymax_global <- max(expr_values_global, na.rm = TRUE)
message("🔹 Globales Y-Maximum: ", round(ymax_global, 2))

# ---- Hilfsfunktion ----
make_violin_cluster <- function(obj, cluster_id, genes, output_dir, width, height, ymax_global) {
  subset_obj <- subset(obj, subset = seurat_clusters == cluster_id)
  
  for (gene in genes) {
    p <- VlnPlot(
      subset_obj,
      features = gene,
      group.by = "condition",
      pt.size = 0,
      assay = "SCT",
      layer = "data",
      cols = c("#9fd8ef", "#fb998e")  # blau = Control, rot = Experimental
    ) +
      ggtitle(paste0(gene, " (Cluster ", cluster_id, ")")) +
      theme_classic(base_size = 8) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)
      ) +
      ylim(0, ymax_global)
    
    # Dateiname
    svg_filename <- file.path(output_dir, paste0("Cluster", cluster_id, "_", gene, "_violin_SCT.svg"))
    
    # SVG speichern
    svg(svg_filename, width = width, height = height)
    print(p)
    dev.off()
    
    message("✅ Saved: ", svg_filename)
  }
}

# ---- Für Cluster 0–4 laufen lassen ----
for (cl in 2:2) {
  make_violin_cluster(
    fib_subset,
    cluster_id = cl,
    genes = genes_to_plot,
    output_dir = output_dir,
    width = plot_width,
    height = plot_height,
    ymax_global = ymax_global
  )
}










library(Seurat)
library(ggplot2)
library(patchwork)

# ---- Output-Ordner ----
output_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure5/Violin_plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Gene-Liste ----
genes_to_plot <- c("Stat1", "Gbp2", "Cxcl9", "Cxcl10", "Cxcl13", "Ccl5")

# ---- 1. Globales Maximum über alle Cluster bestimmen ----
expr_values_global <- FetchData(fib_subset, vars = genes_to_plot, layer = "data", assay = "SCT")
ymax_global <- max(expr_values_global, na.rm = TRUE)
message("🔹 Globales Y-Maximum: ", round(ymax_global, 2))

# ---- Hilfsfunktion: Violinplots pro Cluster ----
make_violin_cluster <- function(obj, cluster_id, genes, output_dir, ymax_global) {
  subset_obj <- subset(obj, subset = seurat_clusters == cluster_id)
  
  plots <- lapply(genes, function(gene) {
    VlnPlot(
      subset_obj,
      features = gene,
      group.by = "condition",
      pt.size = 0,
      assay = "SCT",
      layer = "data",
      cols = c("#9fd8ef", "#fb998e")  # blau = Control, rot = Experimental
    ) +
      ggtitle(paste0(gene, " (Cluster ", cluster_id, ")")) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)
      ) +
      ylim(0, ymax_global)
  })
  
  combined_plot <- wrap_plots(plots, ncol = 3)
  
  svg_filename <- file.path(output_dir, paste0("Cluster", cluster_id, "_IFN_genes_violin_SCT.svg"))
  svg(svg_filename, width = 10, height = 6)
  print(combined_plot)
  dev.off()
  
  message("✅ Saved: ", svg_filename)
}

# ---- Für Cluster 0–4 laufen lassen ----
for (cl in 0:4) {
  make_violin_cluster(fib_subset, cluster_id = cl, genes = genes_to_plot, output_dir = output_dir, ymax_global = ymax_global)
}







library(readxl)
library(writexl)
library(dplyr)
library(stringr)
library(purrr)

# ---- Paths ----
input_path  <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Table S5/Table S5.xlsx"
output_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Table S5/Table_S5_selected_genes.xlsx"

# ---- Genes to keep ----
genes_to_keep <- c("Stat1", "Gbp2", "Cxcl9", "Cxcl10", "Cxcl13", "Ccl5")

# Helper: pick first column whose name matches any pattern
pick_col <- function(nms, patterns) {
  idx <- which(Reduce(`|`, lapply(patterns, function(p) str_detect(nms, regex(p, ignore_case = TRUE)))))
  if (length(idx) == 0) return(NA_integer_)
  idx[[1]]
}

# Standardize one sheet's columns and force cluster from sheet name
read_and_standardize <- function(sheet) {
  df_raw <- read_excel(input_path, sheet = sheet, col_names = TRUE)
  nms <- names(df_raw)
  
  # Find columns by flexible patterns (handles duplicated headers like gene_id...1, etc.)
  i_gene   <- pick_col(nms, c("^gene[_ ]?id$", "^gene$", "gene_id\\.\\.\\."))       # gene_id, gene, gene_id...1
  i_l2fc_s <- pick_col(nms, c("^log2fc\\.single$", "log2fc[_ ]?single", "l2fc.*single"))
  i_padj_s <- pick_col(nms, c("^p[_ ]?val[_ ]?adj\\.single$", "adj.*single"))
  i_l2fc_b <- pick_col(nms, c("^log2fc\\.bulk$", "log2fc[_ ]?bulk", "l2fc.*bulk"))
  i_padj_b <- pick_col(nms, c("^p[_ ]?val[_ ]?adj\\.bulk$", "adj.*bulk"))
  i_cluster<- pick_col(nms, c("^cluster$", "cluster\\.\\.\\."))
  
  # Build clean tibble using available columns
  df <- tibble(
    gene_id          = if (!is.na(i_gene))   df_raw[[i_gene]]   else NA,
    log2FC.single    = if (!is.na(i_l2fc_s)) df_raw[[i_l2fc_s]] else NA,
    p_val_adj.single = if (!is.na(i_padj_s)) df_raw[[i_padj_s]] else NA,
    log2FC.bulk      = if (!is.na(i_l2fc_b)) df_raw[[i_l2fc_b]] else NA,
    p_val_adj.bulk   = if (!is.na(i_padj_b)) df_raw[[i_padj_b]] else NA,
    cluster          = if (!is.na(i_cluster)) df_raw[[i_cluster]] else NA
  )
  
  # Force cluster from sheet name (e.g., "Cluster0" -> 0); keeps existing if numbers fail
  cluster_from_sheet <- suppressWarnings(as.integer(str_extract(sheet, "\\d+")))
  if (!is.na(cluster_from_sheet)) {
    df$cluster <- cluster_from_sheet
  }
  
  # Clean gene_id (trim)
  df$gene_id <- as.character(df$gene_id)
  df$gene_id <- str_trim(df$gene_id)
  
  # Drop rows with no gene_id
  df <- df %>% filter(!is.na(gene_id) & gene_id != "")
  
  df
}

# ---- Read all sheets, standardize, combine ----
sheet_names <- excel_sheets(input_path)
all_data <- map_dfr(sheet_names, read_and_standardize)

# ---- Filter to requested genes ----
filtered_data <- all_data %>%
  filter(gene_id %in% genes_to_keep) %>%
  mutate(cluster = suppressWarnings(as.integer(cluster))) %>%
  arrange(cluster, factor(gene_id, levels = genes_to_keep))

# ---- Write single Excel ----
write_xlsx(filtered_data, output_path)

message("✅ Wrote: ", output_path)

.
















library(Seurat)

genes_ifn_jakstat <- c(
  "Stat1","Stat2","Irf1","Irf7","Irf9",
  "Gbp2","Gbp5","Ifit1","Ifit2","Ifit3","Isg15",
  "Rsad2","Oas1a","Oas2",
  "Cxcl9","Cxcl10","Ccl5",
  "Ciita","Cd74"
)
DefaultAssay(fib_subset) <- "SCT"

fib_subset <- AddModuleScore(
  fib_subset,
  features = list(genes_ifn_jakstat),
  name = "IFN_JAKSTAT",
  nbin = 24
)



DefaultAssay(fib_subset) <- "SCT"

fib_subset <- AddModuleScore(
  fib_subset,
  features = list(genes_ifn_jakstat),
  name = "IFN_JAKSTAT",
  nbin = 24
)



library(ggplot2)

# consistent y-scale across clusters
ylims <- range(fib_subset@meta.data$IFN_JAKSTAT1, na.rm = TRUE)

VlnPlot(
  fib_subset,
  features = "IFN_JAKSTAT1",
  group.by = "seurat_clusters",
  pt.size = 0
) +
  scale_y_continuous(limits = ylims) +
  ylab("IFN/JAK-STAT module score") +
  ggtitle("Interferon & JAK-STAT pathway activity across fibroblast clusters")


FeaturePlot(fib_subset, features = "IFN_JAKSTAT1", min.cutoff = "q05", max.cutoff = "q95")


library(ggplot2)
library(Seurat)

# --- subset to fibroblast clusters 0–4 ---
fib_0to4 <- subset(fib_subset, idents = c(0, 1, 2, 3, 4))

# --- consistent y-limits across all ---
ylims <- range(fib_0to4@meta.data$IFN_JAKSTAT1, na.rm = TRUE)

# --- violin plot ---
p1 <- VlnPlot(
  fib_0to4,
  features = "IFN_JAKSTAT1",
  group.by = "seurat_clusters",
  pt.size = 0
) +
  scale_y_continuous(limits = ylims) +
  ylab("IFN/JAK-STAT module score") +
  ggtitle("Interferon & JAK-STAT pathway activity (Clusters 0–4)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# --- UMAP overlay for those same clusters ---
p2 <- FeaturePlot(
  fib_0to4,
  features = "IFN_JAKSTAT1",
  min.cutoff = "q05",
  max.cutoff = "q95"
) +
  ggtitle("IFN/JAK-STAT module score (UMAP, Clusters 0–4)")

# --- print both ---
p1
p2



library(Seurat)
library(ggplot2)

# --- subset to fibroblast clusters 0–4 ---
fib_0to4 <- subset(fib_subset, idents = c(0, 1, 2, 3, 4))

# --- consistent y-axis limits across all ---
ylims <- range(fib_0to4@meta.data$IFN_JAKSTAT1, na.rm = TRUE)

# --- violin plot: clusters on x, split by condition ---
p_ifn <- VlnPlot(
  fib_0to4,
  features = "IFN_JAKSTAT1",
  group.by = "seurat_clusters",
  split.by = "condition",      # <- assumes metadata column named "condition"
  pt.size = 0
) +
  scale_y_continuous(limits = ylims) +
  scale_fill_manual(values = c("Control" = "#4DBBD5", "Experimental" = "#E64B35")) +
  ylab("IFN/JAK-STAT module score") +
  ggtitle("IFN/JAK-STAT activity in fibroblast clusters\nControl vs Infected") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  )

p_ifn












pC <- VlnPlot(
  subset(fib_subset, idents = 0:4),
  features = "IFN_JAKSTAT1",
  group.by = "seurat_clusters",
  split.by = "condition",
  pt.size = 0
) +
  scale_fill_manual(values=c("Control"="#4DBBD5","Experimental"="#E64B35")) +
  ylab("IFN/JAK-STAT module score") +
  ggtitle("Cluster-specific IFN activation during infection")

pC


pD <- FeaturePlot(
  subset(fib_subset, idents = 0:4),
  features = "IFN_JAKSTAT1",
  min.cutoff="q05", max.cutoff="q95"
)

pD



genes <- c("Stat1","Irf1","Gbp2","Cxcl9","Cxcl10","Ccl5")
pE <- DotPlot(
  subset(fib_subset, idents=0:4),
  features = genes,
  group.by = "seurat_clusters",
  split.by = "condition"
) + RotatedAxis() + ggtitle("Key IFN response genes across fibroblast clusters")

library(dplyr)
fib_df <- fib_subset@meta.data
delta <- fib_df %>%
  filter(seurat_clusters %in% 0:4) %>%
  group_by(seurat_clusters, condition) %>%
  summarise(mean_IFN = mean(IFN_JAKSTAT1), .groups="drop") %>%
  tidyr::pivot_wider(names_from = condition, values_from = mean_IFN) %>%
  mutate(delta = Experimental - Control)

pF <- ggplot(delta, aes(x=seurat_clusters, y=delta, fill=seurat_clusters)) +
  geom_bar(stat="identity") +
  ylab("Δ IFN module score (Infected - Control)") +
  ggtitle("Magnitude of infection-induced IFN response per cluster")




pE
pF







library(Seurat)
library(dplyr)
library(ggplot2)
DefaultAssay(fib_subset) <- "SCT"

# ---- subset to fibroblast clusters 0–4 ----
fib_0to4 <- subset(fib_subset, idents = 0:4)

# =========================
# G) Expression level views
# =========================
# Violin split by condition (same y across clusters)
ylims_ccl5 <- range(FetchData(fib_0to4, "Ccl5", slot = "data"), na.rm = TRUE)

pG1 <- VlnPlot(fib_0to4, features = "Ccl5",
               group.by = "seurat_clusters", split.by = "condition", pt.size = 0) +
  scale_y_continuous(limits = ylims_ccl5) +
  labs(y = "Ccl5 (SCT data)", title = "Ccl5 expression by cluster and condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

# UMAP overlay (clusters 0–4 only)
pG2 <- FeaturePlot(fib_0to4, features = "Ccl5", min.cutoff = "q05", max.cutoff = "q95") +
  ggtitle("Ccl5 on UMAP (fibroblasts 0–4)")

# ======================================
# H) % Ccl5+ cells per cluster/condition
# ======================================
# define positivity threshold (robust; adjust if you prefer another cutoff)
ccl5_vals <- FetchData(fib_0to4, vars = "Ccl5", slot = "data")[,1]
thr <- quantile(ccl5_vals[ccl5_vals > 0], 0.25, na.rm = TRUE)  # lower-quartile of nonzeros

fib_0to4$Ccl5_pos <- FetchData(fib_0to4, "Ccl5", slot = "data")[,1] > thr

# % positive by cluster x condition
pct_df <- fib_0to4@meta.data |>
  transmute(seurat_clusters, condition, Ccl5_pos = as.logical(Ccl5_pos)) |>
  group_by(seurat_clusters, condition) |>
  summarise(pct_pos = mean(Ccl5_pos) * 100, n = n(), .groups = "drop")

pH1 <- ggplot(pct_df, aes(x = seurat_clusters, y = pct_pos, fill = condition)) +
  geom_col(position = position_dodge(width = 0.8)) +
  labs(y = "% Ccl5⁺ cells", x = "Cluster",
       title = "%Ccl5⁺ by cluster (Control vs Infected)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

# Per-sample view (helps show robustness mouse-by-mouse)
stopifnot("orig.ident" %in% colnames(fib_0to4@meta.data))
pct_sample <- fib_0to4@meta.data |>
  transmute(seurat_clusters, condition, orig.ident, Ccl5_pos = as.logical(Ccl5_pos)) |>
  group_by(seurat_clusters, condition, orig.ident) |>
  summarise(pct_pos = mean(Ccl5_pos) * 100, .groups = "drop")

pH2 <- ggplot(pct_sample,
              aes(x = seurat_clusters, y = pct_pos, color = condition, group = orig.ident)) +
  geom_point(position = position_jitter(width = 0.15, height = 0)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black") +
  labs(y = "% Ccl5⁺ cells (per sample)", title = "Per-sample %Ccl5⁺ by cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

# =================================
# I) Coordination with IFN program
# =================================
# DotPlot (control vs infected) for Ccl5 + IFN genes
genes_coor <- c("Ccl5","Cxcl9","Cxcl10","Stat1","Irf1","Gbp2")
pI <- DotPlot(subset(fib_0to4, idents = 0:4),
              features = genes_coor,
              group.by = "seurat_clusters",
              split.by = "condition") +
  RotatedAxis() +
  ggtitle("Ccl5 and IFN-responsive genes across clusters")

# Optional: within-cluster correlation (AvgExpression)
avg <- AverageExpression(fib_0to4, assays = "SCT",
                         features = genes_coor, group.by = "seurat_clusters")$SCT
cor_mat <- cor(t(avg))  # rows = genes
# (Plot with your preferred heatmap tool if desired)

# ============================
# J) Δ %Ccl5⁺ (Inf – Ctrl)
# ============================
delta <- pct_df |>
  tidyr::pivot_wider(names_from = condition, values_from = pct_pos) |>
  mutate(delta_pct = Experimental - Control)

pJ <- ggplot(delta, aes(x = seurat_clusters, y = delta_pct)) +
  geom_col() +
  labs(y = "Δ %Ccl5⁺ (Infected – Control)",
       title = "Infection-induced gain of Ccl5⁺ fibroblasts per cluster")








library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

Idents(fib_subset) <- "seurat_clusters"
cl3 <- subset(fib_subset, idents = 3)
DefaultAssay(cl3) <- "SCT"

# --- define "STAT1-positive" cells ---
stat1_vals <- FetchData(cl3, "Stat1", slot = "data")[,1]
thr <- quantile(stat1_vals[stat1_vals > 0], 0.25, na.rm = TRUE)   # adjust if you prefer 0.1 or 0.5 quantile
cl3$STAT1_pos <- stat1_vals > thr

# --- % STAT1⁺ per mouse and condition ---
df_stat1 <- cl3@meta.data %>%
  group_by(orig.ident, condition) %>%
  summarise(pct_pos = mean(STAT1_pos) * 100, .groups = "drop")

# paired plot
p5d <- ggplot(df_stat1, aes(x = condition, y = pct_pos, group = orig.ident)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  labs(y = "% STAT1⁺ fibroblasts (cluster 3)",
       x = NULL,
       title = "Fraction of STAT1⁺ cluster-3 fibroblasts per mouse") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

# paired Wilcoxon
df_wide <- pivot_wider(df_stat1, names_from = condition, values_from = pct_pos)
p_val_5d <- wilcox.test(df_wide$Experimental, df_wide$Control, paired = TRUE)$p.value
p_val_5d


df_wide



library(ggplot2)
library(dplyr)

# df_stat1 already contains: orig.ident, condition, pct_pos

p5d <- ggplot(df_stat1,
              aes(x = condition, y = pct_pos, color = condition)) +
  geom_jitter(width = 0.15, size = 3) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2, width = 0.4, color = "black") +
  labs(
    x = NULL,
    y = "% STAT1⁺ fibroblasts (cluster 3)",
    title = "Cluster 3 STAT1⁺ fibroblast fraction per mouse"
  ) +
  scale_color_manual(values = c("Control" = "#4DBBD5", "Experimental" = "#E64B35")) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_text(size = 12)
  )

p5d















library(Seurat)
library(ggplot2)

# --- restrict to fibroblast clusters 0–4 ---
fib_0to4 <- subset(fib_subset, idents = c(0,1,2,3,4))

# --- set global y-limits once ---
ylims <- range(fib_0to4@meta.data$IFN_JAKSTAT1, na.rm = TRUE)

# --- violin plot grouped by cluster, split by sample ID ---
VlnPlot(
  fib_0to4,
  features = "IFN_JAKSTAT1",
  group.by = "seurat_clusters",
  split.by = "orig.ident",      # or whatever column stores your sample IDs
  pt.size = 0
) +
  scale_y_continuous(limits = ylims) +
  ylab("IFN/JAK-STAT module score") +
  ggtitle("IFN/JAK-STAT activity across fibroblast clusters (per sample)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

.








library(dplyr)
library(ggplot2)

Idents(fib_subset) <- "seurat_clusters"
cl3 <- subset(fib_subset, idents = 3)

# calculate mean IFN module per mouse
df_ifn <- cl3@meta.data %>%
  group_by(orig.ident, condition) %>%
  summarise(mean_IFN = mean(IFN_JAKSTAT1, na.rm = TRUE), .groups = "drop")

p5d_ifn <- ggplot(df_ifn, aes(x = condition, y = mean_IFN, color = condition)) +
  geom_jitter(width = 0.15, size = 3) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2, width = 0.4, color = "black") +
  labs(
    x = NULL,
    y = "Mean IFN/JAK-STAT module (cluster 3)",
    title = "Per-mouse IFN module in cluster 3 fibroblasts"
  ) +
  scale_color_manual(values = c("Control" = "#4DBBD5", "Experimental" = "#E64B35")) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_text(size = 12)
  )

p5d_ifn






library(Seurat)
library(ggplot2)
library(dplyr)

# subset to cluster 3 fibroblasts
Idents(fib_subset) <- "seurat_clusters"
cl3 <- subset(fib_subset, idents = 3)
DefaultAssay(cl3) <- "SCT"

# ---- (5D) IFN/JAK-STAT module violin ----
p5d_ifn <- VlnPlot(
  cl3,
  features = "IFN_JAKSTAT1",
  group.by = "condition",
  pt.size = 0
)  +
  scale_fill_manual(values = c("Control" = "#4DBBD5", "Experimental" = "#E64B35")) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_text(size = 12)
  )

# optional: per-mouse summary (mean IFN score)
df_ifn <- cl3@meta.data %>%
  group_by(orig.ident, condition) %>%
  summarise(mean_IFN = mean(IFN_JAKSTAT1, na.rm = TRUE), .groups = "drop")

p5d_ifn_permouse <- ggplot(df_ifn, aes(x = condition, y = mean_IFN, color = condition)) +
  geom_jitter(width = 0.15, size = 3) +
  
  scale_color_manual(values = c("Control" = "#4DBBD5", "Experimental" = "#E64B35")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

p5d_ifn
p5d_ifn_permouse




library(Seurat)
library(ggplot2)

# --- subset to fibroblast clusters 0–4 ---
Idents(fib_subset) <- "seurat_clusters"
fib_0to4 <- subset(fib_subset, idents = 0:4)

# --- violin plot for IFN/JAK–STAT module across clusters, split by condition ---
p_ifn_clusters <- VlnPlot(
  fib_0to4,
  features = "IFN_JAKSTAT1",
  group.by = "seurat_clusters",
  split.by = "condition",
  pt.size = 0
) +
  scale_fill_manual(values = c("Control" = "#4DBBD5", "Experimental" = "#E64B35")) +
  ylab("IFN/JAK–STAT module score") +
  xlab("Fibroblast clusters (0–4)") +
  ggtitle("Interferon/JAK–STAT module across fibroblast clusters") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11),
    legend.position = "top",
    legend.title = element_blank()
  )

# --- define save path and dimensions ---
fig5_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure5"
if (!dir.exists(fig5_dir)) dir.create(fig5_dir, recursive = TRUE)

save_path <- file.path(fig5_dir, "Fig5_IFNmodule_clusters0to4_violin.svg")

# --- save as SVG (vector format for Illustrator editing) ---
ggsave(
  filename = save_path,
  plot = p_ifn_clusters,
  width = 6.5,   # inches
  height = 4.5,  # inches
  dpi = 300
)

message("✅ Figure saved as SVG at: ", save_path)






library(dplyr)
library(ggplot2)

Idents(fib_subset) <- "seurat_clusters"
fib_0to4 <- subset(fib_subset, idents = 0:4)

df <- fib_0to4@meta.data %>%
  group_by(seurat_clusters, condition) %>%
  summarise(mean_ifn = mean(IFN_JAKSTAT1, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = condition, values_from = mean_ifn) %>%
  mutate(delta = Experimental - Control)

p_delta <- ggplot(df, aes(x = seurat_clusters, y = delta, fill = delta)) +
  geom_col(width = 0.6) +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35") +
  labs(
    x = "Cluster",
    y = expression(Delta*" IFN/JAK–STAT module (Exp - Ctrl)"),
    title = "Change in IFN/JAK–STAT activity during infection"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

ggsave("C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure5/Fig5_IFNmodule_delta_clusters0to4.svg",
       p_delta, width = 5.5, height = 4, dpi = 300)



genes_ck <- c("Ccl5","Cxcl9","Cxcl10","Ccl2","Ccl7","Ccl8")
DotPlot(fib_subset, features = genes_ck, group.by = "seurat_clusters") +
  scale_color_gradient(low="grey90", high="#E64B35") +
  theme_classic() +
  labs(title="Fibroblast chemokine expression across clusters",
       x="Cluster", y="Gene")





library(Seurat)
library(ggplot2)

# --- subset to fibroblast clusters 0–4 ---
Idents(fib_subset) <- "seurat_clusters"
fib_0to4 <- subset(fib_subset, idents = 0:4)
DefaultAssay(fib_0to4) <- "SCT"

# --- Violin plot for Ccl2 and Ccl5 across clusters, split by condition ---
p_ccl <- VlnPlot(
  fib_0to4,
  features = c("Ccl2", "Ccl5"),
  group.by = "seurat_clusters",
  split.by = "condition",
  pt.size = 0
) +
  scale_fill_manual(values = c("Control" = "#4DBBD5", "Experimental" = "#E64B35")) +
  ylab("Expression (SCT)") +
  xlab("Fibroblast clusters (0–4)") +
  ggtitle("Ccl2 and Ccl5 expression across fibroblast clusters") +
  theme_classic() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 14),
    axis.text.x  = element_text(size = 11, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 11),
    legend.position = "top",
    legend.title    = element_blank()
  )

# --- Save as SVG for Illustrator editing ---
fig5_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure5"
if (!dir.exists(fig5_dir)) dir.create(fig5_dir, recursive = TRUE)

save_path <- file.path(fig5_dir, "Fig5F_Ccl2_Ccl5_violin_clusters0to4.svg")

ggsave(
  filename = save_path,
  plot = p_ccl,
  width = 6.5,   # inches
  height = 5,    # inches
  dpi = 300
)

message("✅ Saved violin plot at: ", save_path)







library(dplyr); library(ggplot2); Idents(fib_subset) <- "seurat_clusters"
fib_0to4 <- subset(fib_subset, idents = 0:4)
df_delta <- fib_0to4@meta.data %>%
  group_by(seurat_clusters, condition) %>%
  summarise(mean_ifn = mean(IFN_JAKSTAT1, na.rm = TRUE),
            sem_ifn = sd(IFN_JAKSTAT1, na.rm = TRUE)/sqrt(dplyr::n()), .groups="drop") %>%
  tidyr::pivot_wider(names_from = condition, values_from = c(mean_ifn, sem_ifn)) %>%
  mutate(delta = mean_ifn_Experimental - mean_ifn_Control)

p_delta <- ggplot(df_delta, aes(x = seurat_clusters, y = delta, fill = delta)) +
  geom_col(width = 0.6) +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35") +
  labs(x = "Cluster", y = expression(Delta*" IFN/JAK–STAT module (Exp–Ctrl)"),
       title = "Infection-induced change in IFN/JAK–STAT activity") +
  theme_classic() + theme(plot.title = element_text(hjust=0.5), legend.position="none")

ggsave("C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure5/Fig5D_Delta_IFNmodule.svg",
       p_delta, width = 5.5, height = 4, dpi = 300)




genes_mhcII <- c(
  "Ciita","Cd74","H2-Aa","H2-Ab1","H2-Eb1",
  "H2-DMa","H2-DMb1","Ctsl","Lgmn"
)


library(Seurat)
library(ggplot2)

Idents(fib_subset) <- "seurat_clusters"
fib_0to4 <- subset(fib_subset, idents = 0:4)
DefaultAssay(fib_0to4) <- "SCT"

# Violin plot across clusters 0–4, split by condition
p_mhcII <- VlnPlot(
  fib_0to4,
  features = "MHCII_MS1",   # or whatever AddModuleScore column name you used
  group.by = "seurat_clusters",
  split.by = "condition",
  pt.size = 0
) +
  scale_fill_manual(values = c("Control"="#4DBBD5", "Experimental"="#E64B35")) +
  ylab("MHC II module score") +
  xlab("Fibroblast clusters (0–4)") +
  ggtitle("Antigen-presentation (MHC II) module across fibroblast subsets") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust=0.5, size=14),
    axis.text.x  = element_text(size=11, angle=45, hjust=1),
    axis.text.y  = element_text(size=11),
    legend.position="top",
    legend.title=element_blank()
  )

# Save as SVG
fig5_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure5"
save_path <- file.path(fig5_dir, "Fig5D_MHCII_module_clusters0to4.svg")

ggsave(
  filename = save_path,
  plot = p_mhcII,
  width = 6.5,
  height = 4.5,
  dpi = 300
)

message("✅ Saved MHCII violin at: ", save_path)


View(fib_subset@meta.data)




library(Seurat)
library(dplyr)
library(tidyr)

# Pick your cluster
cl <- 3
genes <- c("Ccl5","Cxcl13","Cxcl10","Stat1","Gbp2")

obj <- subset(fib_subset, subset = seurat_clusters == cl)

# Pull SCT data layer
mat <- GetAssayData(obj, assay = "SCT", layer = "data")

meta <- obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::select(cell, orig.ident, condition)

# Long table: one row per cell per gene
df_long <- as.data.frame(t(as.matrix(mat[genes, , drop = FALSE]))) %>%
  tibble::rownames_to_column("cell") %>%
  left_join(meta, by = "cell") %>%
  pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expr")

# Summaries per replicate (orig.ident)
rep_summ <- df_long %>%
  group_by(condition, orig.ident, gene) %>%
  summarise(
    n_cells = n(),
    mean_expr = mean(expr),
    pct_expr = mean(expr > 0) * 100,
    mean_expr_nonzero = ifelse(sum(expr > 0) > 0, mean(expr[expr > 0]), NA_real_),
    .groups = "drop"
  )

rep_summ


# ---- Output path ----
out_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure5/Troubleshooting"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(out_dir, "Cluster3_replicate_summary_expression.csv")

# ---- Summarize + export ----
rep_summ %>%
  group_by(condition, gene) %>%
  summarise(
    mean_of_means      = mean(mean_expr),
    mean_pct_expr      = mean(pct_expr),
    mean_expr_nonzero  = mean(mean_expr_nonzero, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(gene, condition) %>%
  write.csv(out_file, row.names = FALSE)

message("✅ CSV saved to: ", out_file)


























library(Seurat)
library(dplyr)

#genes <- c("Ccl5","Cxcl13","Cxcl10","Stat1","Gbp2")
genes <- c("Stat1")

DefaultAssay(fib_subset) <- "SCT"  # change to what you THINK you're plotting
Idents(fib_subset) <- "seurat_clusters"

genes_to_check <- genes

obj2 <- subset(fib_subset, idents = 3)

# Confirm counts per condition inside cluster 2
print(table(obj2$condition))

# Function to summarize per gene
summ_gene <- function(obj, gene, assay_use = DefaultAssay(obj), layer_use = "data") {
  x <- FetchData(obj, vars = gene, assay = assay_use, layer = layer_use)[,1]
  tibble(
    mean_all = mean(x),
    pct_expr = mean(x > 0) * 100,
    mean_nonzero = ifelse(sum(x > 0) > 0, mean(x[x > 0]), NA_real_)
  )
}

# Summaries by condition
res <- lapply(genes_to_check, function(g) {
  bind_rows(
    lapply(split(Cells(obj2), obj2$condition), function(cell_ids) {
      tmp <- subset(obj2, cells = cell_ids)
      summ_gene(tmp, g) %>% mutate(condition = unique(tmp$condition))
    }) %>% bind_rows()
  ) %>% mutate(gene = g)
}) %>% bind_rows()

print(res)

# Now compare SCT vs RNA quickly for one gene
g <- genes_to_check[1]
p1 <- VlnPlot(obj2, features = g, group.by = "condition", assay = "SCT", layer = "data")
p2 <- VlnPlot(obj2, features = g, group.by = "condition", assay = "RNA", layer = "data")
p1 | p2




library(Seurat)
library(dplyr)

Idents(fib_subset) <- "seurat_clusters"
obj2 <- subset(fib_subset, idents = 3)

genes <- genes_to_plot

summ_one_gene <- function(obj, gene, assay_use = "RNA", layer_use = "data") {
  df <- FetchData(obj, vars = c("condition", gene), assay = assay_use, layer = layer_use)
  colnames(df)[2] <- "expr"
  df %>%
    group_by(condition) %>%
    summarise(
      n_cells = n(),
      mean_all = mean(expr),
      pct_expr = mean(expr > 0) * 100,
      mean_nonzero = ifelse(sum(expr > 0) > 0, mean(expr[expr > 0]), NA_real_),
      .groups = "drop"
    ) %>%
    mutate(gene = gene)
}

res <- bind_rows(lapply(genes, \(g) summ_one_gene(obj2, g, assay_use="RNA", layer_use="data")))
print(res)



library(Seurat)
library(dplyr)

Idents(fib_subset) <- "seurat_clusters"
obj2 <- subset(fib_subset, idents = 2)

summ_gene <- function(obj, gene, assay_use = "RNA", layer_use = "data") {
  df <- FetchData(obj, vars = c("condition", gene), assay = assay_use, layer = layer_use)
  colnames(df)[2] <- "expr"
  df %>%
    group_by(condition) %>%
    summarise(
      n_cells = n(),
      mean_all = mean(expr),
      pct_expr = mean(expr > 0) * 100,
      mean_nonzero = ifelse(sum(expr > 0) > 0, mean(expr[expr > 0]), NA_real_),
      .groups = "drop"
    )
}

summ_gene(obj2, "Stat1", assay_use = "RNA", layer_use = "data")



library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

file_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure6/2026_01_14_Hu_slides_quant.xlsx"

df <- read_excel(file_path)
colnames(df)

df <- df %>%
  select(
    Slide,
    inflamed,
    `CD52+ Fib Score`
  )

colnames(df)

df <- df %>%
  rename(
    CD52_score = `CD52+ Fib Score`
  )

df <- df %>%
  mutate(
    Slide = as.factor(Slide),
    inflamed = factor(inflamed, levels = c("N", "Y")),
    CD52_score = as.numeric(CD52_score)
  )
summary(df)
table(df$inflamed, useNA = "ifany")
table(df$CD52_score, useNA = "ifany")
head(df)
tail(df)


df %>%
  group_by(Slide, inflamed) %>%
  summarise(mean_score = mean(CD52_score, na.rm = TRUE))

df_patient <- df %>%
  group_by(Slide, inflamed) %>%
  summarise(
    mean_score = mean(CD52_score, na.rm = TRUE),
    .groups = "drop"
  
df_wide


wilcox_res <- wilcox.test(
  df_wide$Inflamed_mean,
  df_wide$NonInflamed_mean,
  paired = TRUE,
  alternative = "greater",
  exact = FALSE
)

wilcox_res

median_diff <- median(
  df_wide$Inflamed_mean - df_wide$NonInflamed_mean
)

median_diff



df_long <- df_wide %>%
  tidyr::pivot_longer(
    cols = c(NonInflamed_mean, Inflamed_mean),
    names_to = "Region",
    values_to = "MeanScore"
  )

plot <- ggplot(df_long, aes(x = Region, y = MeanScore, group = Slide)) +
  geom_point(size = 3) +
  geom_line(alpha = 0.6) +
  theme_classic() +
  labs(
    x = "",
    y = "Mean CD52+ score per region"
  )

out_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Paper_Figures/Figure6"

ggsave(
  filename = file.path(out_dir, "CD52_inflamed_vs_noninflamed_paired.png"),
  plot = plot,
  width = 4,
  height = 3,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "CD52_inflamed_vs_noninflamed_paired.svg"),
  plot = plot,
  width = 4,
  height = 2
)




# Folder containing the FASTQ files
dir_path <- "D:/Sequencing Stuff Strobel/Full Experiment Data/More fastq files"

# Get file names only
files <- list.files(dir_path)

# Extract the sample prefix (C2–C5 or E2–E5)
sample_id <- sub("^([CE][2-5]).*", "\\1", files)

# Split files by sample
split_files <- split(files, sample_id)

# Convert list to a dataframe where each row is a sample
max_len <- max(sapply(split_files, length))

df <- do.call(rbind, lapply(split_files, function(x) {
  length(x) <- max_len
  x
}))

df <- as.data.frame(df, stringsAsFactors = FALSE)
df$sample <- rownames(df)
df <- df[, c("sample", setdiff(names(df), "sample"))]
rownames(df) <- NULL

# Output CSV in the same folder
output_path <- file.path(dir_path, "fastq_files_by_sample.csv")
write.csv(df, output_path, row.names = FALSE)
