#trying to modify the code a bit of my full analysis.
library(Seurat)
rm(list=ls())
gc()
set.seed(123456)


###QC graphs function code---------

library(Seurat)
library(ggplot2)

# Define the QC graphs function with a specified Seurat object and its name
# this function does the following:
# creates a folder with seurat object name (i change the name as i run things on the object)
# makes a sub folder within the seurat object folder titled QC_plots
# saves the QC plots in that folder
QC_graphs_function <- function(seurat_obj, seurat_object_name) {
  
  # Create a parent folder with the specified Seurat object name
  if (!dir.exists(seurat_object_name)) {
    dir.create(seurat_object_name)
  }
  
  # Create a QC_plots folder inside the parent folder
  qc_plots_folder <- file.path(seurat_object_name, "QC_plots")
  if (!dir.exists(qc_plots_folder)) {
    dir.create(qc_plots_folder)
  }
  
  # Plot 1: Violin plot with points
  #plot1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hemo"), pt.size = 0.1) + NoLegend()
 # ggsave(filename = file.path(qc_plots_folder, "Violin_plot_withpoints.svg"), plot = plot1, device = "svg", width = 10, height = 8)
  
  # Plot 2: Violin plot without points
  plot2 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hemo"), pt.size = 0) + NoLegend()
  ggsave(filename = file.path(qc_plots_folder, "Violin_plots_without_points.svg"), plot = plot2, device = "svg", width = 10, height = 8)
  
  # Plot 3: Scatter plot of nCount_RNA vs. percent.mt
  plot3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  ggsave(filename = file.path(qc_plots_folder, "nCount_RNA_vs_percent_mt.svg"), plot = plot3, device = "svg", width = 5, height = 4)
  
  # Plot 4: Scatter plot of nCount_RNA vs. percent.ribo
  plot4 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.ribo")
  ggsave(filename = file.path(qc_plots_folder, "nCount_RNA_vs_percent_ribo.svg"), plot = plot4, device = "svg", width = 5, height = 4)
  
  # Plot 5: Scatter plot of nCount_RNA vs. percent.hemo
  plot5 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.hemo")
  ggsave(filename = file.path(qc_plots_folder, "nCount_RNA_vs_percent_hemo.svg"), plot = plot5, device = "svg", width = 5, height = 4)
  
  plot6 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  ggsave(filename = file.path(qc_plots_folder, "nCount_RNA_vs_nFeature_RNA.svg"), plot = plot6, device = "svg", width = 5, height = 4)
  
  # Plot 7: Density plot of nUMI per sample
  metadata <- seurat_obj@meta.data
  plot7 <- metadata %>% 
    ggplot(aes(color = orig.ident, x = nCount_RNA, fill = orig.ident)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
  ggsave(filename = file.path(qc_plots_folder, "nUMI_density_per_sample.svg"), plot = plot7, device = "svg", width = 10, height = 8)
  
  # Plot 8: Distribution of genes detected per cell via histogram
  plot8 <- metadata %>% 
    ggplot(aes(color = orig.ident, x = nFeature_RNA, fill = orig.ident)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 300)
  ggsave(filename = file.path(qc_plots_folder, "nGene_density_per_sample.svg"), plot = plot8, device = "svg", width = 10, height = 8)
  
  # Plot 9: Distribution of genes detected per cell via boxplot
  plot9 <- metadata %>% 
    ggplot(aes(x = orig.ident, y = log10(nFeature_RNA), fill = orig.ident)) + 
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle("NCells vs NGenes")
  ggsave(filename = file.path(qc_plots_folder, "nGene_boxplot_per_sample.svg"), plot = plot9, device = "svg", width = 10, height = 8)
  
  # Print message
  cat("QC plots saved in the folder:", qc_plots_folder, "\n")
}

# Example usage:
# QC_graphs_function(seurat_obj = combined_obj_try2_filt_mt, seurat_object_name = "combined_obj_try2_filt_mt")


###starting on the analysis--------
#time go through and analyze each object individually quickly. Kind of like what the core did. 
library(Seurat)
library(dplyr)
library(patchwork)

# Set the working directory
setwd("C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis")

# Define the path to the data folder
data_folder <- "C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/data"


#C3_h5 <- Read10X_h5(file.path(data_folder, 'C3_filtered_feature_bc_matrix.h5'))
#C3_unfiltered <- CreateSeuratObject(counts = C3_h5, project='C3')



# List of sample identifiers
samples <- c("C3", "C4", "C5","E3", "E4", "E5")

# Initialize an empty list to store the Seurat objects
seurat_objects <- list()

# Loop through each sample to read the data and create Seurat objects
for (sample in samples) {
  # Read the .h5 file
  h5_data <- Read10X_h5(file.path(data_folder, paste0(sample, "_filtered_feature_bc_matrix.h5")))
  
  # Create a Seurat object
  seurat_obj <- CreateSeuratObject(counts = h5_data, project = sample)
  
  # Store the Seurat object in the list
  seurat_objects[[sample]] <- seurat_obj
}

# Access each Seurat object using the sample name
C3_unfiltered <- seurat_objects[["C3"]]
C4_unfiltered <- seurat_objects[["C4"]]
C5_unfiltered <- seurat_objects[["C5"]]
E3_unfiltered <- seurat_objects[["E3"]]
E4_unfiltered <- seurat_objects[["E4"]]
E5_unfiltered <- seurat_objects[["E5"]]




# Set the base working directory
setwd("C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/")

# Create the "combined_obj" folder if it doesn't already exist
if (!dir.exists("combined_obj_try4")) {
  dir.create("combined_obj_try4")
}

# Set the working directory to the new "C3" folder
setwd(file.path(getwd(), "combined_obj_try4"))

# Print the current working directory to confirm
print(getwd())


combined_obj_try4 <- merge(C3_unfiltered, y = c(C4_unfiltered, C5_unfiltered, E3_unfiltered, E4_unfiltered, E5_unfiltered), add.cell.ids = c("C3", "C4", "C5", "E3", "E4", "E5" ))


# Create a new column 'condition' based on 'orig.ident'
combined_obj_try4@meta.data$condition <- ifelse(grepl("^C", combined_obj_try4@meta.data$orig.ident), "Control", "Experimental")


combined_obj_try4_unf <- combined_obj_try4

View(combined_obj_try4_unf@meta.data)
table(combined_obj_try4_unf$orig.ident)



##change object name 3 times
# Run the standard Seurat pipeline
combined_obj_try4_unf[["percent.mt"]] <- PercentageFeatureSet(combined_obj_try4_unf, pattern = "^mt-")
combined_obj_try4_unf[["percent.ribo"]] <- PercentageFeatureSet(combined_obj_try4_unf, pattern = '^Rpl|^Rps')
combined_obj_try4_unf[["percent.hemo"]] <- PercentageFeatureSet(combined_obj_try4_unf, pattern = '^Hba|^Hbb')

QC_graphs_function(seurat_obj = combined_obj_try4_unf, seurat_object_name = "combined_obj_try4_unf")



View(combined_obj_try4_unf)

combined_obj_try4_unf
#Note that since the data is split into layers, normalization and variable feature identification is performed for each batch independently (a consensus set of variable features is automatically identified). https://satijalab.org/seurat/articles/seurat5_integration#:~:text=Once%20integrative%20analysis%20is%20complete,performing%20any%20differential%20expression%20analysis.




##here i am keeping track of the number of cells removed at each step to keep an idea
# Define the file name for the CSV file within the current working directory
filtering_csv_path <- "cell_filtering_info_per_sample.csv"

# Get the initial cell count per sample from your Seurat object
initial_cell_count_per_sample <- table(combined_obj_try4_unf$orig.ident)

# Create a data frame with the step as the first column and cell counts as the rest
filtering_info <- data.frame(Step = "Read_in", t(as.matrix(initial_cell_count_per_sample)))

# Save the initial information to the CSV file
write.csv(filtering_info, file = filtering_csv_path, row.names = FALSE)

# Print confirmation message
cat("Filtering information saved to:", filtering_csv_path, "\n")


###sidequest to count number of times we have a cell with either Krt5 or Krt14
num_cells_with_Krt5_Krt14_Trp63 <- sum(FetchData(combined_obj_try4_unf, vars = c("Krt5", "Krt14", 'Trp63')) > 0)
num_cells_with_Krt5 <- sum(FetchData(combined_obj_try4_unf, vars = c("Krt5")) > 0)
num_cells_with_Krt14 <- sum(FetchData(combined_obj_try4_unf, vars = c("Krt14")) > 0)
num_cells_with_Trp63 <- sum(FetchData(combined_obj_try4_unf, vars = c("Trp63")) > 0)



num_cells_with_Krt5_Krt14_Trp63 #256
num_cells_with_Krt5 #55
num_cells_with_Krt14 #36
num_cells_with_Trp63 #174



# 
# ###sidequest to quickly go through and normalize data and make a featureplot of Krt5, Krt14, Trp63. 
# #spoilers... i still dont think that i have basel cells. 
# # Normalize the data
# combined_obj_try4_unf <- NormalizeData(combined_obj_try4_unf)
# combined_obj_try4_unf <- FindVariableFeatures(combined_obj_try4_unf, selection.method = "vst", nfeatures = 3000)
# combined_obj_try4_unf <- ScaleData(combined_obj_try4_unf)
# combined_obj_try4_unf <- RunPCA(combined_obj_try4_unf, features = VariableFeatures(object = combined_obj_try4_unf))
# combined_obj_try4_unf <- FindNeighbors(combined_obj_try4_unf, dims = 1:30)
# combined_obj_try4_unf <- FindClusters(combined_obj_try4_unf, resolution = 0.5)
# combined_obj_try4_unf <- RunUMAP(combined_obj_try4_unf, reduction = "pca", dims = 1:30)
# FeaturePlot(combined_obj_try4_unf, features = c("Krt5", "Krt14", 'Trp63'), raster = FALSE)
# ggsave("combined_obj_try4_unf_featureplot_Krt5_Krt14_Trp63_2.svg", width = 10, height = 10, units = "in", dpi = 300)



###filtering my seurat object very basic min cells 3 and min features 200--------------------------------------------------

# Set the working directory
#setwd("C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis")

# Define the path to the data folder
data_folder <- "C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/data"

# List of sample identifiers
samples <- c("C3", "C4", "C5", "E3", "E4", "E5")

# Initialize an empty list to store the Seurat objects
seurat_objects <- list()

# Loop through each sample to read the data and create Seurat objects
for (sample in samples) {
  # Read the .h5 file
  h5_data <- Read10X_h5(file.path(data_folder, paste0(sample, "_filtered_feature_bc_matrix.h5")))
  
  # Create a Seurat object
  seurat_obj <- CreateSeuratObject(counts = h5_data, project = sample, min.cells = 3, min.features = 200)
  
  # Store the Seurat object in the list
  seurat_objects[[sample]] <- seurat_obj
}


# Access each Seurat object using the sample name
C3_filtered <- seurat_objects[["C3"]]
C4_filtered <- seurat_objects[["C4"]]
C5_filtered <- seurat_objects[["C5"]]
E3_filtered <- seurat_objects[["E3"]]
E4_filtered <- seurat_objects[["E4"]]
E5_filtered <- seurat_objects[["E5"]]



combined_obj_try4_filt <- merge(C3_filtered, y = c(C4_filtered, C5_filtered,E3_filtered, E4_filtered, E5_filtered), add.cell.ids = c("C3", "C4", "C5","E3", "E4", "E5" ))


# Create a new column 'condition' based on 'orig.ident'
combined_obj_try4_filt@meta.data$condition <- ifelse(grepl("^C", combined_obj_try4_filt@meta.data$orig.ident), "Control", "Experimental")

head(combined_obj_try4_filt@meta.data)

#adding cell count to my csv
# Get the cell count per sample after filtering
cell_count_after_filtering <- table(combined_obj_try4_filt$orig.ident)

# Create a data frame for this filtering step
filtering_info_after_filtering <- data.frame(Step = "After_Filtering", t(as.matrix(cell_count_after_filtering)))

# Append this information to the original CSV file
write.table(filtering_info_after_filtering, file = filtering_csv_path, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)

# Print confirmation message
cat("Filtering information after filtering step saved to:", filtering_csv_path, "\n")



##change object name 6 times
# Run the standard Seurat pipeline
combined_obj_try4_filt[["percent.mt"]] <- PercentageFeatureSet(combined_obj_try4_filt, pattern = "^mt-")
combined_obj_try4_filt[["percent.ribo"]] <- PercentageFeatureSet(combined_obj_try4_filt, pattern = '^Rpl|^Rps')
combined_obj_try4_filt[["percent.hemo"]] <- PercentageFeatureSet(combined_obj_try4_filt, pattern = '^Hba|^Hbb')


#Change object name twice
#redoing QC metrics
#Plot the QC metrics
#QC_graphs_function(seurat_obj = combined_obj_try4_filt, seurat_object_name = "combined_obj_try4_filt")








###adding and filtering by 3MAD mitochondrial counts---------------------------------------------------------------------

#some good info on filtering https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html


combined_obj_try4_filt

# Sample: C3 Threshold for mitochondrial content: 6.699135 
# Sample: C4 Threshold for mitochondrial content: 7.610519 
# Sample: C5 Threshold for mitochondrial content: 7.902271 
# Sample: E3 Threshold for mitochondrial content: 5.496056 
# Sample: E4 Threshold for mitochondrial content: 5.725418 
# Sample: E5 Threshold for mitochondrial content: 6.566773 


# Get the metadata
metadata <- combined_obj_try4_filt@meta.data
head(metadata)
# Initialize an empty list to store thresholds for each sample
thresholds <- list()

# Loop through each unique sample in the metadata to calculate thresholds
for (sample in unique(metadata$orig.ident)) {
  
  # Subset metadata for the current sample
  sample_metadata <- metadata[metadata$orig.ident == sample, ]
  
  # Correctly access 'percent.mt' from the metadata of Seurat object
  percent_mt_sample <- sample_metadata$percent.mt
  
  # Calculate median and MAD directly using the correct data
  median_percent_mt <- median(percent_mt_sample, na.rm = TRUE)
  mad_percent_mt <- mad(percent_mt_sample, na.rm = TRUE)
  
  # Calculate the threshold with 3 MAD
  threshold_mito <- median_percent_mt + 3 * mad_percent_mt
  
  # Store the threshold in the list
  thresholds[[sample]] <- threshold_mito
  
  # Print the threshold for the current sample
  cat("Sample:", sample, "Threshold for mitochondrial content:", threshold_mito, "\n")
}


# Initialize a list to store the filtered Seurat objects
filtered_objects <- list()

# Loop through each unique sample to filter by respective threshold
for (sample in unique(metadata$orig.ident)) {
  
  # Subset the cells in the Seurat object for the current sample
  sample_seurat_obj <- subset(combined_obj_try4_filt, subset = orig.ident == sample)
  
  # Apply the filter based on the respective threshold
  threshold_mito <- thresholds[[sample]]
  sample_filtered <- subset(sample_seurat_obj, subset = percent.mt < threshold_mito)
  
  # Store the filtered Seurat object in the list
  filtered_objects[[sample]] <- sample_filtered
  
  # Print the number of cells remaining after filtering for the current sample
  cat("Sample:", sample, "Number of cells after filtering:", ncol(sample_filtered), "\n")
}

# Combine all the filtered objects back into one Seurat object
combined_obj_try4_filt <- merge(filtered_objects[[1]], y = filtered_objects[-1])

# Print confirmation message
cat("Combined filtered Seurat object created with", ncol(combined_obj_try4_filt), "cells.\n")




# Get the cell count per sample after filtering
cell_count_after_filtering <- table(combined_obj_try4_filt$orig.ident)

# Create a data frame for this filtering step
filtering_info_after_filtering <- data.frame(Step = "After_3MAD_mt_filtering", t(as.matrix(cell_count_after_filtering)))

# Append this information to the original CSV file
write.table(filtering_info_after_filtering, file = filtering_csv_path, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)

# Print confirmation message
cat("Filtering information after filtering step saved to:", filtering_csv_path, "\n")

head(combined_obj_try4_filt@meta.data)

#Change object name twice
#redoing QC metrics
#Plot the QC metrics
QC_graphs_function(seurat_obj = combined_obj_try4_filt, seurat_object_name = "combined_obj_try2_filt_mt")



###remove 10% hemo cells-----------------------------------------------------------------------------------------------
ncol(combined_obj_try4_filt)

# Define the threshold for hemoglobin genes (10%)
hemo_threshold <- 10

# Filter the Seurat object to remove cells with more than 10% hemoglobin genes
combined_obj_try4_filt <- subset(combined_obj_try4_filt, subset = percent.hemo <= hemo_threshold)

# Print a message with the number of remaining cells
cat("Cells remaining after filtering for hemoglobin genes:", ncol(combined_obj_try4_filt), "\n")


QC_graphs_function(seurat_obj = combined_obj_try4_filt, seurat_object_name = "combined_obj_try4_filt_hemo")


# Get the cell count per sample after filtering
cell_count_after_filtering <- table(combined_obj_try4_filt$orig.ident)

# Create a data frame for this filtering step
filtering_info_after_filtering <- data.frame(Step = "After_10pct_hemo", t(as.matrix(cell_count_after_filtering)))

# Append this information to the original CSV file
write.table(filtering_info_after_filtering, file = filtering_csv_path, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)

# Print confirmation message
cat("Filtering information after filtering step saved to:", filtering_csv_path, "\n")



#saveRDS(combined_obj_try4_filt, "combined_obj_try4_filt.rds")


###normalize data with log transformation (for cell cycle score) ---------------------------------------------

##change object name once, and make new object
# Normalize the data

 combined_obj_try4_filt <- NormalizeData(combined_obj_try4_filt, normalization.method = "LogNormalize", scale.factor = 10000)






###cell cycle scoring my object so that it's all done for later-----------------------------------------------------



# Split the merged Seurat object by sample
samples <- SplitObject(combined_obj_try4_filt, split.by = "orig.ident")



#this is my gene list for cell cycle. It comes with Seurat. There is reference in the function. 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


#time to format our gene lists properly.....
# Function to capitalize the first letter and make the rest lowercase
capitalize_first_letter <- function(x) {
  paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))
}

# Apply the function to the original vector
s.genes1 <- sapply(s.genes, capitalize_first_letter)

# Remove the names from the new vector
s.genes1 <- unname(s.genes1)



# Function to capitalize the first letter and make the rest lowercase
capitalize_first_letter <- function(x) {
  paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))
}
# Apply the function to the original vector
g2m.genes1 <- sapply(g2m.genes, capitalize_first_letter)

# Remove the names from the new vector
g2m.genes1 <- unname(g2m.genes1)






# Split the merged Seurat object by sample
samples <- SplitObject(combined_obj_try4_filt, split.by = "orig.ident")

# Initialize an empty list to store the Seurat objects after cell cycle scoring
scored_samples <- list()

# Loop through each sample, perform cell cycle scoring, and store the results
for (i in seq_along(samples)) {
  # Run cell cycle scoring on each sample
  samples[[i]] <- CellCycleScoring(samples[[i]], s.features = s.genes1, g2m.features = g2m.genes1, set.ident = FALSE)
  
  # Store the scored Seurat object in the list
  scored_samples[[i]] <- samples[[i]]
}

# Merge the scored Seurat objects back together
combined_obj_try4_filt <- merge(scored_samples[[1]], y = scored_samples[-1], add.cell.ids = names(samples), project = "Merged_Seurat")

# View the metadata to confirm that the cell cycle scores have been added
head(combined_obj_try4_filt@meta.data)



#saveRDS(combined_obj_try4_filt, "combined_obj_try4_filt_final.rds")

###3MAD Filtering-----------------------------------------------------------------------------------------------

# Enhanced MAD filtering function with detailed output
#designed to filter ncount and nfeature by 3MAD. Then it runs qc metrics and prints them out for me.
filter_3MAD <- function(seurat_obj, metric_name) {
  # Extract the metric data from the Seurat object's metadata
  metric_data <- seurat_obj@meta.data[[metric_name]]
  
  # Calculate median, MAD, and thresholds
  median_metric <- median(metric_data, na.rm = TRUE)
  mad_metric <- mad(metric_data, na.rm = TRUE)
  threshold_lower <- median_metric - 3 * mad_metric
  threshold_upper <- median_metric + 3 * mad_metric
  
  # Determine min and max values in the data
  min_metric_before <- min(metric_data, na.rm = TRUE)
  max_metric_before <- max(metric_data, na.rm = TRUE)
  
  # Filter the Seurat object based on the thresholds
  initial_cell_count <- ncol(seurat_obj)
  filtered_obj <- subset(seurat_obj, cells = which(metric_data >= threshold_lower & metric_data <= threshold_upper))
  cells_removed <- initial_cell_count - ncol(filtered_obj)
  
  # Determine min and max values after filtering
  filtered_metric_data <- filtered_obj@meta.data[[metric_name]]
  min_metric_after <- min(filtered_metric_data, na.rm = TRUE)
  max_metric_after <- max(filtered_metric_data, na.rm = TRUE)
  
  # Print detailed information
  cat("Processing metric:", metric_name, "\n")
  cat("Minimum value before filtering:", min_metric_before, "\n")
  cat("Maximum value before filtering:", max_metric_before, "\n")
  cat("Minimum value after filtering:", min_metric_after, "\n")
  cat("Maximum value after filtering:", max_metric_after, "\n")
  cat("Median:", median_metric, "\n")
  cat("MAD:", mad_metric, "\n")
  cat("Lower threshold (3MAD):", threshold_lower, "\n")
  cat("Upper threshold (3MAD):", threshold_upper, "\n")
  cat("Initial number of cells:", initial_cell_count, "\n")
  cat("Number of cells removed:", cells_removed, "\n")
  cat("Number of cells after filtering:", ncol(filtered_obj), "\n\n")
  
  # Return the filtered object and the number of cells removed
  list(filtered_obj = filtered_obj, cells_removed = cells_removed, initial_cell_count = initial_cell_count)
}



# Function to save filtering information per sample to the CSV
save_filtering_info <- function(seurat_obj, step_name, csv_path) {
  # Get the cell count per sample after filtering
  cell_count_after_filtering <- table(seurat_obj$orig.ident)
  
  # Create a data frame for this filtering step
  filtering_info_after_filtering <- data.frame(Step = step_name, t(as.matrix(cell_count_after_filtering)))
  
  # Append this information to the original CSV file
  write.table(filtering_info_after_filtering, file = csv_path, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
  
  # Print confirmation message
  cat("Filtering information for", step_name, "saved to:", csv_path, "\n")
}

# Apply 3MAD filtering for nFeature_RNA
result_nfeature <- filter_3MAD(combined_obj_try4_filt, "nFeature_RNA")
combined_obj_try4_filt <- result_nfeature$filtered_obj
save_filtering_info(combined_obj_try4_filt, "Filter_by_3MAD_nFeature_RNA", filtering_csv_path)

# Apply 3MAD filtering for nCount_RNA
result_ncount <- filter_3MAD(combined_obj_try4_filt, "nCount_RNA")
combined_obj_try4_filt <- result_ncount$filtered_obj
save_filtering_info(combined_obj_try4_filt, "Filter_by_3MAD_nCount_RNA", filtering_csv_path)



#running the QC function on my final object now
QC_graphs_function(seurat_obj = combined_obj_try4_filt, seurat_object_name = "combined_obj_try4_filt_3MAD")





###running doublet finder----------------------------------------------------------------------------------------

library(scDblFinder)
library(Seurat)
library(SingleCellExperiment)
set.seed(123456)

# Split the merged Seurat object by sample
samples <- SplitObject(combined_obj_try4_filt, split.by = "orig.ident")

# Initialize a list to store filtered Seurat objects after doublet removal
filtered_samples <- list()

# Loop through each sample, run scDblFinder, and remove doublets
for (i in seq_along(samples)) {
  
  # Print which sample is being processed
  cat("Processing sample:", names(samples)[i], "\n")
  
  # Convert Seurat object to SingleCellExperiment
  sce <- as.SingleCellExperiment(samples[[i]])
  
  # Detect doublets
  sce <- scDblFinder(sce)
  
  # Convert the doublet information to Seurat-compatible format
  seurat_obj <- samples[[i]]
  seurat_obj$doublet <- colData(sce)$scDblFinder.class
  
  # Save doublet information
  doublet_table <- table(seurat_obj$doublet)
  doublet_df <- as.data.frame(doublet_table)
  doublet_csv_path <- paste0("doublet_detection_info_", names(samples)[i], ".csv")
  write.csv(doublet_df, file = doublet_csv_path, row.names = FALSE)
  cat("Doublet detection table saved to:", doublet_csv_path, "\n")
  
  # Subset the Seurat object to retain only singlets
  seurat_obj <- subset(seurat_obj, subset = doublet == "singlet")
  
  # Store the filtered object
  filtered_samples[[i]] <- seurat_obj
  
  # Print confirmation message for the sample
  cat("Finished processing sample:", names(samples)[i], "\n\n")
}


###change name here once
# Merge the filtered Seurat objects back together
combined_obj_try4_filt <- merge(filtered_samples[[1]], y = filtered_samples[-1], add.cell.ids = names(samples), project = "Merged_Seurat")
combined_obj_try4_filt
head(combined_obj_try4_filt@meta.data)


##change name here once
# Get the cell count per sample after filtering
cell_count_after_filtering <- table(combined_obj_try4_filt$orig.ident)

# Create a data frame for this filtering step
filtering_info_after_filtering <- data.frame(Step = "After_doublet_removal", t(as.matrix(cell_count_after_filtering)))

# Append this information to the original CSV file
write.table(filtering_info_after_filtering, file = filtering_csv_path, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)

# Print confirmation message
cat("Filtering information after filtering step saved to:", filtering_csv_path, "\n")




QC_graphs_function(seurat_obj = combined_obj_try4_filt, seurat_object_name = "combined_obj_try4_filt_db")


saveRDS(combined_obj_try4_filt, "combined_obj_try4_filt_unmerged.rds")

beepr::beep(sound = 8)












###SCTransform with no variables; determining if i have to integrate; CONCLUSION: I have to integrate-------

names(combined_obj_try4_filt@assays$RNA@layers)


head(combined_obj_try4_filt@meta.data)

combined_obj_try4_filt <- readRDS("combined_obj_try4_filt_unmerged.rds")

SCT_novars <- readRDS('C:\\Users\\ostrobel\\Indiana University\\Jerde, Travis J - JERDE LAB (J Drive)\\Strobel\\Single Cell Analysis Stuff\\R_Analysis_Strobel\\Combined_Sample_Analysis\\Combined_analysis\\combined_obj_try4\\SCTransformTry1\\SCT_no_vars\\try3_SCT_novar_nojoin.rds')

###regressing no variables
# Required Libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(openxlsx)  # For writing Excel files

setwd('C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4')

# Define your base directory
getwd()

# Set your base directory to the current working directory
base_dir <- getwd()
base_dir

#dir.create(file.path(base_dir, "SCTransformTry1"), showWarnings = TRUE)
# Navigate to the "SCTransform" folder within the current working directory
sctransform_dir <- file.path(base_dir, "SCTransformTry1")


# Create a subfolder named "SCT_no_vars" within the "SCTransform" folder
sct_no_vars_dir <- file.path(sctransform_dir, "SCT_no_vars")
#dir.create(sct_no_vars_dir, showWarnings = TRUE)

# Print the base directory to confirm
cat("SCT_no_vars directory set to:", sct_no_vars_dir, "\n")

SCTnovars_dir <- file.path(sctransform_dir, "SCT_no_vars")


##start here!!!!!!
seurat_obj <- combined_obj_try4_filt

# Apply SCTransform without regression
seurat_obj <- SCTransform(seurat_obj, verbose = TRUE)
?SCTransform
head(seurat_obj@meta.data)


# Standard workflow: PCA, UMAP, Clustering
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)


# Generate UMAP and save
umap_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = c("condition", "seurat_clusters"))
ggsave(filename = file.path(SCTnovars_dir, "umap_by_cond_clust.svg"), plot = umap_plot, device = "svg", width = 16, height = 8)

# Create DimPlot with PCA and save
pca_plot <- DimPlot(seurat_obj, reduction = "pca")
ggsave(filename = file.path(SCTnovars_dir, "pca.svg"), plot = pca_plot, device = "svg", width = 8, height = 6)



cell_cycle_plot <- DimPlot(seurat_obj,
                           reduction = "pca",
                           group.by= "Phase",
                           split.by = "Phase")
ggsave(filename = file.path(SCTnovars_dir, "cell_cycle_pca.svg"), plot = cell_cycle_plot, device = "svg", width = 18, height = 6)
cell_cycle_plot
# Create ElbowPlot and save
elbow_plot <- ElbowPlot(seurat_obj, ndims = 50)
ggsave(filename = file.path(SCTnovars_dir, "elbow.svg"), plot = elbow_plot, device = "svg", width = 8, height = 6)

# Create VariableFeaturePlot and save
plot1 <- VariableFeaturePlot(object = seurat_obj)
top10 <- head(VariableFeatures(seurat_obj), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave(filename = file.path(SCTnovars_dir, "var_feat.svg"), plot = plot1 + plot2, device = "svg", width = 8, height = 6)

# Generate feature plots for each key marker gene
features <- c("Ptprc", "Vim", "Epcam", "Acta2", "Cdk1", "Mki67", "Svs2")
for (feature in features) {
  if (feature %in% rownames(seurat_obj)) {
    plot <- FeaturePlot(seurat_obj, features = feature, reduction = "umap")
    ggsave(filename = file.path(SCTnovars_dir, paste0(feature, ".svg")), 
           plot = plot, device = "svg", width = 8, height = 6)
  } else {
    cat("Feature", feature, "not found in the Seurat object.\n")
  }
}

# Generate Excel table with cluster and sample counts
clusters <- Idents(seurat_obj)
cluster_table <- table(clusters, seurat_obj$orig.ident)
cluster_df <- as.data.frame.matrix(cluster_table)

cluster_df$Total_Control <- rowSums(cluster_df[, grep("^C", colnames(cluster_df))])
cluster_df$Total_Experimental <- rowSums(cluster_df[, grep("^E", colnames(cluster_df))])

cluster_df <- cbind(seurat_clusters = rownames(cluster_df), cluster_df)

excel_file_path <- file.path(base_path, "2024_10_07clust_counts.xlsx")
write.xlsx(cluster_df, file = excel_file_path, rowNames = FALSE)

cat("Excel table saved to:", excel_file_path, "\n")



# Save the Seurat object in the respective folder
saveRDS(seurat_obj, file = "try3_SCT_novar_nojoin.rds")

beepr::beep(5)

##it looks ok. Some looks a little weird with how the Epcam and Vim cells are together in areas but split by batches. Definitely need to integrate. I will also try the conventional method next. 


###time to determine which feature selection method I want to use and see if i can get away without integrating.; CONCLUSION: I have to integrate ----------------------------------------------------------------------



##this is my code to test which feature selection method works best. There are 3: dispersion, mean.var, vst. I anticipate VST, but one should try all. 
#this code takes my current base directory, reads in the folder of my current seurat object name, and then creates 3 folders (one for each selection method). WIthin that folder will be the key marker genes plots, dimplot, feature plot, pca, umap, and then a top_genes_plots folder. That folder contains featuremaps and violin plots of the top 5 genes for each cluster. 

#overall this code works really well and is very useful. I will use it for the big object as well.


#I NOW RUN JOIN LAYERS BECAUSE EVERYTHING IS NORMALIZED.
#in the seurat tutorial https://satijalab.org/seurat/articles/seurat5_integration this is done after integration because integration needs separate layers. I am hoping to avoid integration. Let's see. That is something else to explore.


try4_filt <- combined_obj_try4_filt
# Required Libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Define your base directory
base_dir <- getwd()
base_dir

# Function to generate UMAP, other plots, and an Excel table with cluster and sample counts
generate_umap_and_plots <- function(seurat_obj, selection_method) {
  
  # Create a folder based on the Seurat object name
  seurat_object_name <- deparse(substitute(seurat_obj))
  object_folder <- file.path(base_dir, seurat_object_name)
  
  # Create the main object folder if it doesn't exist
  dir.create(object_folder, showWarnings = FALSE)
  
  # Create subfolders based on the selection method within the Seurat object folder
  method_folder <- file.path(object_folder, selection_method)
  dir.create(method_folder, showWarnings = FALSE)
  
  # Standard workflow: PCA, UMAP, Clustering
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = selection_method)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  
  # Generate UMAP split by condition and colored by Seurat clusters
  umap_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = c("condition", "seurat_clusters"))
  ggsave(filename = file.path(method_folder, "umap_plot_by_condition_clusters.svg"), plot = umap_plot, device = "svg", width = 16, height = 8)
  
  # Create DimPlot with PCA and save
  pca_plot <- DimPlot(seurat_obj, reduction = "pca")
  ggsave(filename = file.path(method_folder, "pca_plot.svg"), plot = pca_plot, device = "svg", width = 8, height = 6)
  
  # Create ElbowPlot and save
  elbow_plot <- ElbowPlot(seurat_obj, ndims = 50)
  ggsave(filename = file.path(method_folder, "elbow_plot.svg"), plot = elbow_plot, device = "svg", width = 8, height = 6)
  
  # Create VariableFeaturePlot and save
  plot1 <- VariableFeaturePlot(seurat_obj)
  top10 <- head(VariableFeatures(seurat_obj), 10)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  ggsave(filename = file.path(method_folder, "variable_feature_plot.svg"), plot = plot1 + plot2, device = "svg", width = 8, height = 6)
  
  # Generate feature plots for each key marker gene
  features <- c("Ptprc", "Vim", "Epcam", "Acta2", "Cdk1", "Mki67", "Spink1")
  for (feature in features) {
    if (feature %in% rownames(seurat_obj)) {
      plot <- FeaturePlot(seurat_obj, features = feature, reduction = "umap")
      ggsave(filename = file.path(method_folder, paste0(feature, ".svg")), 
             plot = plot, device = "svg", width = 8, height = 6)
    } else {
      cat("Feature", feature, "not found in the Seurat object.\n")
    }
  }
  
  # Generate Excel table with cluster and sample counts
  clusters <- Idents(seurat_obj)
  cluster_table <- table(clusters, seurat_obj$orig.ident)
  cluster_df <- as.data.frame.matrix(cluster_table)
  
  cluster_df$Total_Control <- rowSums(cluster_df[, grep("^C", colnames(cluster_df))])
  cluster_df$Total_Experimental <- rowSums(cluster_df[, grep("^E", colnames(cluster_df))])
  
  cluster_df <- cbind(seurat_clusters = rownames(cluster_df), cluster_df)
  
  excel_file_path <- file.path(method_folder, "cluster_counts.xlsx")
  write.xlsx(cluster_df, file = excel_file_path, rowNames = FALSE)
  
  cat("Excel table saved to:", excel_file_path, "\n")
  
  # Save the Seurat object in the respective folder
  saveRDS(seurat_obj, file = file.path(method_folder, paste0(seurat_object_name, "_", selection_method, ".rds")))
}



# Generate and save plots for each selection method
generate_umap_and_plots(try4_filt, selection_method = "vst")
generate_umap_and_plots(try4_filt, selection_method = "mean.var.plot")
generate_umap_and_plots(try4_filt, selection_method = "dispersion")

##it might be worth trying a resolution of 0.4. Haven't messed with that yet thoguh. 

beepr::beep(5)












###trying integration with Harmony; THIS WAS MY METHOD FORWARD---------

#seurat_obj <- readRDS("C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/SCTransformTry1/SCT_no_vars/try3_SCT_novar_nojoin.rds")



#loading my finished object for future things
seurat_obj <- readRDS("C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/try4_HarmonyIntegration.rds")


DimPlot(seurat_obj, reduction = "umap.harmonyintegration", group.by = "seurat_clusters", label = TRUE)


# Print a summary of the object to confirm successful loading
print(seurat_obj)



library(harmony)

getwd()
seurat_obj <- IntegrateLayers(object=seurat_obj, method=HarmonyIntegration, normalization.method="SCT", orig.reduction="pca", new.reduction="integrated.harmony")
seurat_obj

#Determine number of PCs to use
ElbowPlot(seurat_obj, ndims=50)

Find neighbors with 30 PCs
seurat_obj <- FindNeighbors(seurat_obj, reduction="integrated.harmony", dims=1:30)
# 
seurat_obj_0.4 <- FindClusters(seurat_obj, resolution=0.4, cluster.name="harmony_clusters")
seurat_obj_0.4 <- RunUMAP(seurat_obj_0.4, reduction="integrated.harmony", dims=1:30, reduction.name="umap.harmony")
DimPlot(seurat_obj_0.4, reduction="umap.harmony", group.by=c("harmony_clusters"), label=TRUE)
DimPlot(seurat_obj_0.4, reduction="umap.harmony", group.by=c("condition"), label=TRUE)

SCTHarm_dir <- "C:\\Users\\ostrobel\\Indiana University\\Jerde, Travis J - JERDE LAB (J Drive)\\Strobel\\Single Cell Analysis Stuff\\R_Analysis_Strobel\\Combined_Sample_Analysis\\Combined_analysis\\combined_obj_try4\\Harmony_SCT"


# Generate UMAP and save
#umap_plot <- DimPlot(seurat_obj_0.4, reduction = "umap.harmony", group.by = c("condition", "harmony_clusters"))
#ggsave(filename = file.path(SCTHarm_dir, "png"), plot = umap_plot, device = "png", width = 16, height = 8)


# Generate UMAP colored by Seurat clusters and save
umap_seurat_clusters_plot <- DimPlot(seurat_obj, reduction = 'umap.harmonyintegration', group.by = "seurat_clusters", label = TRUE)
ggsave(filename = file.path(SCTHarm_dir, "umap_by_seurat_clusters.png"), plot = umap_seurat_clusters_plot, device = "png", width = 8, height = 8)

# Generate UMAP colored by condition and the method-specific clusters and save
umap_cond_clusters_plot <- DimPlot(seurat_obj, reduction = 'umap.harmonyintegration', group.by = c("condition", "harmonyintegration_clusters"), label = TRUE)
ggsave(filename = file.path(SCTHarm_dir, "umap_by_cond_clust.png"), plot = umap_cond_clusters_plot, device = "png", width = 16, height = 8)


# Generate UMAP colored by sampleID s and save
umap_sample_clusters_plot <- DimPlot(seurat_obj, reduction = 'umap.harmonyintegration', group.by = 'orig.ident', label = TRUE)
ggsave(filename = file.path(SCTHarm_dir, "umap_by_sample_clust.png"), plot = umap_sample_clusters_plot, device = "png", width = 8, height = 8)
View(seurat_obj@meta.data)
cell_cycle_plot <- DimPlot(seurat_obj,
                           reduction = "umap.harmonyintegration",
                           group.by= "Phase",
                           split.by = "Phase")
ggsave(filename = file.path(SCTHarm_dir, "cell_cycle_pca.svg"), plot = cell_cycle_plot, device = "png", width = 18, height = 6)
cell_cycle_plot
# Create ElbowPlot and save
elbow_plot <- ElbowPlot(seurat_obj_0.4, ndims = 50)
ggsave(filename = file.path(SCTHarm_dir, "elbow.png"), plot = elbow_plot, device = "png", width = 8, height = 6)


# Generate feature plots for each key marker gene
features <- c("Ptprc", "Vim", "Epcam", "Acta2", "Cdk1", "Mki67", "Svs2")
for (feature in features) {
  if (feature %in% rownames(seurat_obj)) {
    plot <- FeaturePlot(seurat_obj, features = feature, reduction = "umap.harmonyintegration")
    ggsave(filename = file.path(SCTHarm_dir, paste0(feature, ".png")), 
           plot = plot, device = "png", width = 5, height = 5)
  } else {
    cat("Feature", feature, "not found in the Seurat object.\n")
  }
}





# Generate Excel table with cluster and sample counts
clusters <- Idents(seurat_obj)
cluster_table <- table(clusters, seurat_obj$orig.ident)
cluster_df <- as.data.frame.matrix(cluster_table)

cluster_df$Total_Control <- rowSums(cluster_df[, grep("^C", colnames(cluster_df))])
cluster_df$Total_Experimental <- rowSums(cluster_df[, grep("^E", colnames(cluster_df))])

cluster_df <- cbind(seurat_clusters = rownames(cluster_df), cluster_df)

excel_file_path <- file.path(SCTHarm_dir, "clust_countsv2.xlsx")
write.xlsx(cluster_df, file = excel_file_path, rowNames = FALSE)

cat("Excel table saved to:", excel_file_path, "\n")

# SCTHarm_dir
# # Save the Seurat object in the SCTHarm_dir folder
# saveRDS(seurat_obj, file = file.path(SCTHarm_dir, "try3_SCT_novar_nojoin.rds"))

beepr::beep(5)


#SCTHarm_obj <- readRDS("C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/try3_SCT_novar_nojoin.rds")

# getwd()
# saveRDS(SCTHarm_obj, "2024_09_20_FINAL_SUERAT_OBJECT.rds")
# 
# seurat_obj <- readRDS("C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/2024_09_20_FINAL_SUERAT_OBJECT.rds")

seurat_clusters_plot <- DimPlot(seurat_obj, reduction = "umap.harmonyintegration", group.by = "seurat_clusters", label = TRUE)
ggsave(filename = file.path(SCTHarm_dir, "umap_by_seurat_clusters.png"), plot = seurat_clusters_plot, device = "png", width = 8, height = 8)


# List all reductions in the Seurat object
names(seurat_obj@reductions)

# Keep only the 'umap' reduction and remove all other reductions
seurat_obj@reductions <- seurat_obj@reductions["umap.harmonyintegration"]

# Verify that only the desired reduction is kept
names(seurat_obj@reductions)


#saving a loupe file
# Assuming seurat_obj is the object you loaded earlier using readRDS()
library(loupeR)

# Modify the create_loupe_from_seurat function to use seurat_obj
create_loupe_from_seurat(
  seurat_obj,
  output_dir = 'C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT',
  output_name = 'Strobel_Final_Loupev2.cloupe'
)
DimPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE)




###continuing with my analysis. Looking at Fibroblasts/Stromal--------

# Set the base directory
base_path <- "C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT"

# Set the working directory to the base path
setwd(base_path)

# Create subfolders
subfolders <- c("Fib_subset", "Epi_subset", "Immune_subset")
for (subfolder in subfolders) {
  dir.create(file.path(base_path, subfolder), showWarnings = FALSE)
}

# Generate paths for each subfolder
fib_subset_path <- file.path(base_path, "Fib_subset")
epi_subset_path <- file.path(base_path, "Epi_subset")
immune_subset_path <- file.path(base_path, "Immune_subset")

# Print paths to confirm
cat("Fib subset path:", fib_subset_path, "\n")
cat("Epi subset path:", epi_subset_path, "\n")
cat("Immune subset path:", immune_subset_path, "\n")


#loading in my object 
seurat_obj <- readRDS("C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/2024_09_20_FINAL_SUERAT_OBJECT.rds")


library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)


##subsetting and analyzing my fibroblasts
# 1. Create a subset of the original Seurat object for the specified clusters
clusters_of_interest <- c(0,1,5,12,16,20,22,23)
fib_subset <- subset(seurat_obj, idents = clusters_of_interest)

# 2. Store the original cluster assignments in a new metadata column
fib_subset$original_cluster <- Idents(seurat_obj)
DefaultAssay(fib_subset) <- "RNA"

# 3. Perform clustering on the subsetted data
fib_subset <- SCTransform(fib_subset, verbose = TRUE)  # Normalize and scale data
fib_subset <- RunPCA(fib_subset, verbose = TRUE)       # Run PCA
elbow_plot <- ElbowPlot(fib_subset, ndims = 50)

# Define the file path for the Elbow Plot
elbow_plot_path <- file.path(fib_subset_path, "2024_10_07_elbow_plot.png")

# Save the Elbow Plot as a PNG file
ggsave(filename = elbow_plot_path, plot = elbow_plot, device = "png", width = 10, height = 8)


fib_subset <- RunUMAP(fib_subset, dims = 1:30, verbose = TRUE)  # Run UMAP
fib_subset <- FindNeighbors(fib_subset, dims = 1:30)    # Find neighbors
fib_subset <- FindClusters(fib_subset, resolution = 0.1) # Recluster with a new resolution

# Visualization of the new clusters with reference to original clusters
plot1 <- DimPlot(fib_subset, group.by = "original_cluster", label = TRUE)  # Original clusters
plot2 <- DimPlot(fib_subset, label = TRUE)  # New clusters
combined_plot <- plot1 + plot2

# Save combined plot
ggsave(file.path(fib_subset_path, "2024_10_07_combined_plot.png"), plot = combined_plot, device = "png", width = 16, height = 8)

#new clusters
plot3 <- DimPlot(fib_subset, group.by = "orig.ident", label = TRUE) 
ggsave(file.path(fib_subset_path, "2024_10_07_sample_orig_plot.png"), plot = plot3, device = "png", width = 6, height = 6)


# Generate Excel table with cluster and sample counts
clusters <- Idents(fib_subset)
cluster_table <- table(clusters, fib_subset$orig.ident)
cluster_df <- as.data.frame.matrix(cluster_table)

cluster_df$Total_Control <- rowSums(cluster_df[, grep("^C", colnames(cluster_df))])
cluster_df$Total_Experimental <- rowSums(cluster_df[, grep("^E", colnames(cluster_df))])

cluster_df <- cbind(seurat_clusters = rownames(cluster_df), cluster_df)

# Define the correct file path
excel_file_path <- file.path(fib_subset_path, "2024_10_07_fib_clust_counts.xlsx")

# Save the data frame as an Excel file
write.xlsx(cluster_df, file = excel_file_path, rowNames = FALSE)

cat("Excel table saved to:", excel_file_path, "\n")



# Define the path for saving the Seurat object
seurat_save_path <- file.path(fib_subset_path, "2024_10_07_fib_subset.rds")

# Save the Seurat object
saveRDS(fib_subset, file = seurat_save_path)

cat("Seurat object saved to:", seurat_save_path, "\n")

##Couple more plots
# Define the features you want to plot
features <- c("Ptprc", "Vim", "Epcam", "Acta2", "Cdk1", "Mki67")

# Loop through each feature to create and save individual FeaturePlots
for (feature in features) {
  feature_plot <- FeaturePlot(fib_subset, features = feature, pt.size = 0.5)
  
  # Define the path for saving each FeaturePlot
  feature_plot_path <- file.path(fib_subset_path, paste0("2024_10_07_feature_plot_", feature, ".png"))
  
  # Save each FeaturePlot as PNG
  ggsave(filename = feature_plot_path, plot = feature_plot, device = "png", width = 6, height = 6)
  
  cat("FeaturePlot for", feature, "saved to:", feature_plot_path, "\n")
}


# Create UMAP plot split by condition
umap_plot <- DimPlot(fib_subset, reduction = "umap", group.by = "condition", pt.size = 0.5)

# Define the file path for saving the plot
umap_plot_path <- file.path(fib_subset_path, "2024_10_07_umap_split_by_condition.png")

# Save the UMAP plot as a PNG file
ggsave(filename = umap_plot_path, plot = umap_plot, device = "png", width = 6, height = 6)

cat("UMAP plot split by condition saved to:", umap_plot_path, "\n")



fib_subset2 <- PrepSCTFindMarkers(fib_subset)

# Find and save markers
markers <- FindAllMarkers(fib_subset2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Group markers by cluster and find the top 10 markers based on p-value for each cluster
top_markers_per_cluster <- top_n(markers[order(markers$p_val), ], n = 10, wt = markers$p_val)

# Arrange by cluster and p_val
top_markers_per_cluster <- top_markers_per_cluster[order(top_markers_per_cluster$cluster, top_markers_per_cluster$p_val), ]


write.csv(top_markers_per_cluster, "C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/Fib_subset/2024_10_07_Fibs_Top_Markers_Per_Cluster0.1.csv", row.names = FALSE)





library(DESeq2)
library(Seurat)
library(ggplot2)
library(presto)
library(dplyr)
library(openxlsx)




# View the metadata to confirm the structure
View(fib_subset2@meta.data)

# Check if the two columns are exactly the same
identical(fib_subset2@meta.data$seurat_clusters, fib_subset2@meta.data$SCT_snn_res.0.1)

# Set the default assay to RNA to focus on the original data
DefaultAssay(fib_subset2) <- "RNA"

# Aggregate expression data by original identity, condition, and Seurat clusters
fib_subset.pseudo <- AggregateExpression(fib_subset, return.seurat=TRUE, group.by=c("orig.ident", "condition", "seurat_clusters"))
View(fib_subset.pseudo@meta.data)

# Set identifiers for the new pseudo object
Idents(fib_subset.pseudo) <- "seurat_clusters"


# Define the path for saving the Seurat object
seurat_save_path <- file.path(fib_subset_path, "2024_10_07_fib_subset_pseudo.rds")

# Save the Seurat object
saveRDS(fib_subset.pseudo, file = seurat_save_path)





# 
# # Differential expression analysis using DESeq2 for cluster 0
# DE.bulk.0 <- FindMarkers(fib_subset.pseudo, ident.1="0", min.cells.group=0, test.use="DESeq2")
# 
# # View results and details about cells in the pseudo object
# View(DE.bulk.0)
# tail(Cells(fib_subset.pseudo))



# Define the original and pseudo objects for marker finding
obj_original <- fib_subset2
obj.pseudo <- fib_subset.pseudo



# Set identities based on Seurat clusters
Idents(obj_original) <- "seurat_clusters"
Idents(obj.pseudo) <- "seurat_clusters"

# Set the RNA assay as the default assay
DefaultAssay(object = obj_original) <- "RNA"
DefaultAssay(object = obj.pseudo) <- "RNA"

obj_original <- JoinLayers(obj_original, assay = "RNA")
obj.pseudo <- JoinLayers(obj.pseudo, assay = "RNA")


# Loop through clusters to find markers, comparing single cell and pseudo-bulk data
n <- 7 #this is the number of your LAST cluster. Pointer makes it so that we start at 0
keep <- data.frame()
wb_cluster_DE <- createWorkbook()
ptr <- 0

for(i in 0:n){
  print(paste("Processing cluster", i))
  
  DE_temp <- FindMarkers(obj_original, ident.1=paste(i), test.use="wilcox", slot="data")
  DE_temp <- DE_temp[c(2:5)]
  colnames(DE_temp) <- c("log2FC.single", "pct.current_cluster", "pct.other_clusters", "padj.single")
  DE_temp$gene_id <- rownames(DE_temp)
  
  bulk_temp <- FindMarkers(object=obj.pseudo, ident.1=paste(i), min.cells.group=0, test.use="DESeq2")
  bulk_temp <- bulk_temp[c(2,5)]
  colnames(bulk_temp) <- c("log2FC.bulk","padj.bulk")
  bulk_temp$gene_id <- rownames(bulk_temp)
  
  temp <- merge(DE_temp,bulk_temp,by="gene_id", all = TRUE)
  
  temp2 <- temp
  temp2$cluster <- paste(i)
  temp.sig.either <- subset(temp2, padj.single < 0.05 | padj.bulk < 0.05)
  temp.sig.both <- subset(temp2, padj.single < 0.05 & padj.bulk < 0.05)
  
  keep <- rbind(keep, temp.sig.both)
  
  addWorksheet(wb_cluster_DE, paste("cluster", i, "all", sep="."))
  addWorksheet(wb_cluster_DE, paste("cluster", i, "sig", "either", sep="."))
  addWorksheet(wb_cluster_DE, paste("cluster", i, "sig", "both", sep="."))
  
  writeData(wb_cluster_DE, sheet=(ptr+1), temp, rowNames=FALSE)
  writeData(wb_cluster_DE, sheet=(ptr+2), temp.sig.either, rowNames=FALSE)
  writeData(wb_cluster_DE, sheet=(ptr+3), temp.sig.both, rowNames=FALSE)
  
  ptr <- ptr+3
}

# Save the workbook summarizing differential expression results
saveWorkbook(wb_cluster_DE, file.path(fib_subset_path, "2024_10_07_DE_summary_fib_subset.xlsx"), overwrite = TRUE)


# Export a summary of markers for manual annotation
keep.pos <- subset(keep, log2FC.single > 0 & log2FC.bulk > 0)
id <- data.frame("gene_id"=keep.pos[duplicated(keep.pos[,1]),1])
keep.pos.unique <- dplyr::anti_join(keep.pos, id, by="gene_id")
keep.pos.overlapping <- merge(keep.pos, id, by="gene_id")

wb_cluster_summary <- createWorkbook()
addWorksheet(wb_cluster_summary, "all_markers")
addWorksheet(wb_cluster_summary, "all_positive_markers")
addWorksheet(wb_cluster_summary, "all_positive_markers_overlap")
addWorksheet(wb_cluster_summary, "all_positive_markers_unique")
writeData(wb_cluster_summary, sheet=1, keep, rowNames=FALSE)
writeData(wb_cluster_summary, sheet=2, keep.pos, rowNames=FALSE)
writeData(wb_cluster_summary, sheet=3, keep.pos.overlapping, rowNames=FALSE)
writeData(wb_cluster_summary, sheet=4, keep.pos.unique, rowNames=FALSE)

# Save the workbook summarizing differential expression results
saveWorkbook(wb_cluster_DE, file.path(fib_subset_path, "2024_10_07_Fib_susbet_markers_for_man_annot.xlsx"), overwrite = TRUE)



###Differential gene experession bulk and single cell control vs experimental FIBS-----------

# Set the base path for saving the volcano plots
volcano_plot_dir <- file.path(fib_subset_path, "Volcano_Plots_2024_10_07")

# Create the directory if it doesn't exist
if (!dir.exists(volcano_plot_dir)) {
  dir.create(volcano_plot_dir)
}


library(ggplot2)

# Ensure the correct identities are set based on the Seurat clusters
Idents(obj_original) <- "seurat_clusters"
Idents(obj.pseudo) <- "seurat_clusters"





head(obj_original@meta.data)
head(obj.pseudo@meta.data)
# Create a new column 'condition' based on 'orig.ident'
obj_original@meta.data$condition <- ifelse(grepl("^C", obj_original@meta.data$orig.ident), "Control", "Experimental")

head(obj_original@meta.data)
head(obj.pseudo@meta.data)




library(ggplot2)
library(ggrepel)


# Initialize an empty list to store results for each cluster
results_list <- list()

# Loop through each cluster
clusters <- unique(Idents(obj_original))
for (cluster in clusters) {
  print(paste("Processing cluster", cluster))
  
  # Subset the Seurat objects for the current cluster
  single_cell_cluster <- subset(obj_original, idents = cluster)
  pseudo_bulk_cluster <- subset(obj.pseudo, idents = cluster)
  
  # Ensure the correct identities are set based on the condition
  Idents(single_cell_cluster) <- single_cell_cluster$condition
  Idents(pseudo_bulk_cluster) <- pseudo_bulk_cluster$condition
  
  # Perform differential expression analysis for single-cell data
  DE_single <- FindMarkers(single_cell_cluster, ident.1 = "Experimental", ident.2 = "Control", test.use = "wilcox", slot = "data")
  DE_single <- DE_single %>% select(log2FC.single = avg_log2FC, p_val_adj.single = p_val_adj)
  DE_single$gene_id <- rownames(DE_single)
  
  # Perform differential expression analysis for pseudo-bulk data
  DE_bulk <- FindMarkers(pseudo_bulk_cluster, ident.1 = "Experimental", ident.2 = "Control", min.cells.group = 0, test.use = "DESeq2")
  DE_bulk <- DE_bulk %>% select(log2FC.bulk = avg_log2FC, p_val_adj.bulk = p_val_adj)
  DE_bulk$gene_id <- rownames(DE_bulk)
  
  # Merge single-cell and pseudo-bulk results
  merged_results <- merge(DE_single, DE_bulk, by = "gene_id", all = TRUE)
  merged_results$cluster <- cluster
  
  # Store the results in the list
  results_list[[as.character(cluster)]] <- merged_results
  
  # Select top 25 upregulated and downregulated genes for labeling
  top_upregulated <- merged_results %>% filter(p_val_adj.bulk < 0.05) %>% arrange(desc(log2FC.bulk)) %>% head(25)
  top_downregulated <- merged_results %>% filter(p_val_adj.bulk < 0.05) %>% arrange(log2FC.bulk) %>% head(25)
  top_genes <- bind_rows(top_upregulated, top_downregulated)
  
  # Create a volcano plot for the cluster
  volcano_plot <- ggplot(merged_results, aes(x = log2FC.bulk, y = -log10(p_val_adj.bulk))) +
    geom_point(aes(color = ifelse(p_val_adj.bulk < 0.05 & log2FC.bulk > 0, "Upregulated",
                                  ifelse(p_val_adj.bulk < 0.05 & log2FC.bulk < 0, "Downregulated", "Not Significant")))) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
    labs(title = paste("Volcano Plot - Cluster", cluster),
         x = "Log2 Fold Change (Bulk)",
         y = "-Log10 Adjusted P-value",
         color = "Gene Regulation") +
    geom_text_repel(data = top_genes, aes(label = gene_id), size = 3, max.overlaps = 20) +
    theme_minimal()
  
  # Save the volcano plot
  ggsave(filename = file.path(volcano_plot_dir, paste0("Volcano_Cluster_", cluster, ".png")),
         plot = volcano_plot, device = "png", width = 8, height = 6)
  
  print(paste("Volcano plot saved for cluster", cluster))
}

# Combine all clusters into a single data frame
final_results <- bind_rows(results_list)

# Save the results to an Excel file
output_file <- file.path(fib_subset_path, "2024_10_07_Fib_DE_Comp_Cont_vs_Exp.xlsx")
write.xlsx(final_results, file = output_file, rowNames = FALSE)

cat("Results saved to:", output_file, "\n")


#### ptprc positive number in each cluster by condition-------------
# Define a threshold for PTPRC positivity (you may adjust this based on your data)
ptprc_threshold <- 0  # Assuming any expression greater than 0 is considered positive

# Identify PTPRC-positive cells
ptprc_positive_cells <- WhichCells(fib_subset, expression = Ptprc > ptprc_threshold)

# Subset the Seurat object to only include PTPRC-positive cells
fib_subset_ptprc_positive <- subset(fib_subset, cells = ptprc_positive_cells)

# Get the number of PTPRC-positive cells in each cluster, split by condition
ptprc_positive_counts_by_condition <- table(Idents(fib_subset_ptprc_positive), fib_subset_ptprc_positive$condition)

# Convert the result to a data frame for easy viewing
ptprc_positive_counts_by_condition_df <- as.data.frame(ptprc_positive_counts_by_condition)
colnames(ptprc_positive_counts_by_condition_df) <- c("Cluster", "Condition", "PTPRC_Positive_Cells")

# Print the result
print(ptprc_positive_counts_by_condition_df)


# Define the file path for saving
excel_file_path <- file.path(fib_subset_path, "Ptprc_Positive_Cells.xlsx")

# Save the result to an Excel file
write.xlsx(ptprc_positive_counts_by_condition_df, file = excel_file_path, rowNames = FALSE)

cat("Excel file saved to:", excel_file_path, "\n")





###determining the number of cells in my object, and minus immune--------
# Generate Excel table with cluster and sample counts
clusters <- Idents(fib_subset)
cluster_table <- table(clusters, fib_subset$orig.ident)
cluster_df <- as.data.frame.matrix(cluster_table)

cluster_df$Total_Control <- rowSums(cluster_df[, grep("^C", colnames(cluster_df))])
cluster_df$Total_Experimental <- rowSums(cluster_df[, grep("^E", colnames(cluster_df))])

cluster_df <- cbind(seurat_clusters = rownames(cluster_df), cluster_df)

excel_file_path <- file.path(fib_subset_path, "2024_10_07Fib_clust_counts.xlsx")
write.xlsx(cluster_df, file = excel_file_path, rowNames = FALSE)

cat("Excel table saved to:", excel_file_path, "\n")

#loading in my original object again to make sure it's good
# seurat_obj <- readRDS("C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/2024_09_20_FINAL_SUERAT_OBJECT.rds")


# Get the total number of cells in Control condition
total_control_cells <- sum(seurat_obj$condition == "Control")

# Get the total number of cells in Experimental condition
total_experimental_cells <- sum(seurat_obj$condition == "Experimental")

# Print the results
print(paste("Total number of cells in Control condition:", total_control_cells))
print(paste("Total number of cells in Experimental condition:", total_experimental_cells))


# Specify the clusters to exclude
clusters_to_exclude <- c(2,6,7,10,11,17,19,21,24)

# Subset the Seurat object to exclude the specified clusters
filtered_imm_subset <- subset(seurat_obj, idents = setdiff(unique(Idents(seurat_obj)), clusters_to_exclude))

# Get the total number of cells in Control condition after excluding clusters
total_control_cells_filtered <- sum(filtered_imm_subset$condition == "Control")

# Get the total number of cells in Experimental condition after excluding clusters
total_experimental_cells_filtered <- sum(filtered_imm_subset$condition == "Experimental")

# Print the results
print(paste("Total number of cells in Control condition after excluding clusters:", total_control_cells_filtered))

#33121"
print(paste("Total number of cells in Experimental condition after excluding clusters:", total_experimental_cells_filtered))
#26590"











####running DESEq2 and GSEA on ALL OF MY FIBROBLASTS COMBINED--------------

library(DESeq2)
library(fgsea)
library(org.Mm.eg.db)
library(openxlsx)
library(dplyr)
library(tidyr)
library(Seurat)


fib_subset_path <- "C:\\Users\\ostrobel\\Indiana University\\Jerde, Travis J - JERDE LAB (J Drive)\\Strobel\\Single Cell Analysis Stuff\\R_Analysis_Strobel\\Combined_Sample_Analysis\\Combined_analysis\\combined_obj_try4\\Harmony_SCT\\Fib_subset\\2024_10_07_fib_subset_FINAL.rds"

fib_subset_pseudo_path <- "C:\\Users\\ostrobel\\Indiana University\\Jerde, Travis J - JERDE LAB (J Drive)\\Strobel\\Single Cell Analysis Stuff\\R_Analysis_Strobel\\Combined_Sample_Analysis\\Combined_analysis\\combined_obj_try4\\Harmony_SCT\\Fib_subset\\2024_10_07_fib_sub_psudo_FINAL.rds"

fib_subset <- readRDS(fib_subset_path)

fib_subset.pseudo <- readRDS(fib_subset_pseudo_path)

#Aggregate expression data by original identity, condition, and Seurat clusters
fib_subset.pseudo_all <- AggregateExpression(fib_subset, return.seurat=TRUE, group.by=c("orig.ident", "condition"))
fib_subset.pseudo_all@meta.data

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
my_pseudo_all <- fib_subset.pseudo_all






# Create a DESeqDataSet object from the Seurat object
dds <- DESeqDataSetFromMatrix(
  countData = GetAssayData(my_pseudo_all, layer = "counts"),
  colData = my_pseudo_all@meta.data,
  design = ~ condition  
)

dds

rowSums(counts(dds))
# Run the DESeq2 analysis
dds <- DESeq(dds)

# Extract results including the stat column
de_results <- results(dds, contrast = c("condition", "Experimental", "Control"))

# Summarize the results
summary(de_results)

head(de_results)


# Convert DESeq2 results object to a data frame
de_results_df <- as.data.frame(de_results)
head(de_results_df)


de_results_df$gene <- rownames(de_results_df)  # Add gene names as a column
head(de_results_df)

# Remove rows with NA values 
###should i remove NA Values? IDFK. https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-are-some-p-values-set-to-na
de_results_df <- de_results_df %>%
  drop_na()



# Assuming de_results_df is your data frame with the 'rank_metric' column already calculated
# Count the total number of rows
total_count <- nrow(de_results_df)

# Count how many rows have the same 'rank_metric' score
duplicate_counts <- de_results_df %>%
  group_by(stat) %>%               # Group by the rank_metric column
  filter(n() > 1) %>%                     # Filter for groups with more than 1 row (duplicates)
  summarise(count = n())                  # Summarise to get counts of each duplicate score

# Calculate the total number of rows with duplicate scores
total_duplicates <- sum(duplicate_counts$count)

# Calculate the percentage overlap
percent_overlap <- (total_duplicates / total_count) * 100

# Display results
cat("Total rows:", total_count, "\n")
cat("Rows with duplicate scores:", total_duplicates, "\n")
cat("Percentage overlap:", percent_overlap, "%\n")




# Run GSEA
pathways_file <- "C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Core_Object_Analysis/GSEA/MSigDB/M2/m2.all.v2023.2.Mm.symbols.gmt"

# Load pathways
pathways <- gmtPathways(pathways_file)



# Prepare the ranking metric for GSEA, using the stat score
ranks <- de_results_df$stat
names(ranks) <- rownames(de_results_df)

length(ranks)
#13363
head(ranks)

# # Filter out non-finite values from ranks
# ranks <- ranks[is.finite(ranks)]
# length(ranks)
# #19228

# Run GSEA using fgsea
fgsea_results <- fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500)
#.94% were ties


# Filter pathways by adjusted p-value (significant pathways)
significant_pathways <- fgsea_results %>%
  dplyr::filter(padj < 0.05)

# Sort pathways by NES (largest to smallest)
sorted_pathways <- significant_pathways %>%
  dplyr::arrange(desc(NES))

# Convert the `leadingEdge` list column to a character string for CSV export
sorted_pathways$leadingEdge <- sapply(sorted_pathways$leadingEdge, function(x) paste(x, collapse = ";"))

# Get today's date in the format YYYY_MM_DD
today_date <- format(Sys.Date(), "%Y_%m_%d")


deseq2_dir <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/DESeq2_outputs_allfib"

# Define the file path to save the results, including today's date
save_path <- file.path(deseq2_dir, paste0("ALLFIB_fgsea_sign_paths.csv"))

# Save the results as a CSV file
write.csv(sorted_pathways, file = save_path, row.names = FALSE)

# Print a message to confirm saving
cat("Filtered GSEA results saved to:", save_path, "\n")

# Create volcano plot directory if it doesn't exist
volcano_plot_dir <- file.path(deseq2_dir, "Volcano_Plots")
if (!dir.exists(volcano_plot_dir)) {
  dir.create(volcano_plot_dir, recursive = TRUE)
}
# Generate volcano plot
library(ggplot2)
library(ggrepel)

# Determine the top 50 genes to label
top_genes <- de_results_df %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(50)

# Generate the volcano plot
volcano_plot <- ggplot(de_results_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05 & log2FoldChange > 0, "Upregulated",
                                ifelse(padj < 0.05 & log2FoldChange < 0, "Downregulated", "Not Significant")))) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "Gene Regulation") +
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes)), size = 3, max.overlaps = 20) +
  theme_minimal()

# Save the volcano plot
volcano_plot_path <- file.path(volcano_plot_dir, "Volcano_Plot_50.png")
ggsave(filename = volcano_plot_path, plot = volcano_plot, device = "png", width = 8, height = 6)

# Print a message to confirm saving
cat("Volcano plot saved to:", volcano_plot_path, "\n")










# Define your output folder path
output_folder <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/DESeq2_outputs_allfib/ALLFIB"

# Prepare the ranking metric for GSEA, using the stat score
ranks <- de_results$stat
names(ranks) <- rownames(de_results)

# Filter out non-finite values from ranks
ranks <- ranks[is.finite(ranks)]

# Enrichment Plot for the Top Pathway
top_pathway <- sorted_pathways$pathway[1]  # Select the top enriched pathway
enrichment_plot <- plotEnrichment(pathways[[top_pathway]], ranks) +
  labs(title = top_pathway, x = "Rank", y = "Enrichment Score")

# Save the enrichment plot
enrichment_plot_path <- file.path(output_folder, "GSEA_Enrichment_Plot.png")
ggsave(filename = enrichment_plot_path, plot = enrichment_plot, width = 8, height = 6)
cat("Enrichment plot saved to:", enrichment_plot_path, "\n")

# Bar Plot of Top Pathways by NES
top_n_pathways <- 25  # Adjust the number of pathways you want to visualize
bar_plot <- ggplot(sorted_pathways[1:top_n_pathways, ], aes(x = reorder(pathway, NES), y = NES, fill = NES)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 10 Enriched Pathways", x = "Pathway", y = "Normalized Enrichment Score") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "grey", midpoint = 0) +
  theme_minimal()

# Save the bar plot
bar_plot_path <- file.path(output_folder, "Top_10_Enriched_Pathways.png")
ggsave(filename = bar_plot_path, plot = bar_plot, width = 8, height = 6)
cat("Bar plot saved to:", bar_plot_path, "\n")

# Dot Plot for Pathway Significance
dot_plot <- ggplot(sorted_pathways, aes(x = NES, y = reorder(pathway, NES), size = -log10(padj), color = NES)) +
  geom_point() +
  scale_color_gradient2(low = "blue", high = "red", mid = "grey", midpoint = 0) +
  labs(title = "Significant Pathways", x = "Normalized Enrichment Score", y = "Pathway") +
  theme_minimal()

# Save the dot plot
dot_plot_path <- file.path(output_folder, "Dot_Plot_Significant_Pathways.png")
ggsave(filename = dot_plot_path, plot = dot_plot, width = 8, height = 6)
cat("Dot plot saved to:", dot_plot_path, "\n")




####GSEA and DESEq2 on individual clusters----------

summary(de_results)

###trying this a different way
# Loop through each cluster
clusters <- unique(fib_subset.pseudo$seurat_clusters)

for (cluster in clusters) {
  cat("Processing cluster:", cluster, "\n")
  
  # Create directory for the current cluster
  cluster_dir <- file.path(deseq2_dir, paste0("Cluster", cluster))
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
  
  #https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
  # Filter genes with fewer than 10 total counts
  #keep_genes <- rowSums(counts(dds)) >= 10
  smallestGroupSize <- 3
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  
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
  
  # Run GSEA
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
    head(25)
  
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
         x = "Log2 Fold Change", y = "-log10 Adjusted P-value", color = "Legend") + # Change the legend title to "Legend"
    theme_minimal()
  
  
  # Save the volcano plot for the current cluster
  volcano_plot_path <- file.path(cluster_dir, paste0("V2Volcano_Plot_Cluster_", cluster, ".png"))
  ggsave(filename = volcano_plot_path, plot = volcano_plot, device = "png", width = 8, height = 6)
  cat("Volcano plot saved for cluster", cluster, "to:", volcano_plot_path, "\n")
}

cat("All clusters processed successfully.\n")






##combining the sheets into one
# Load necessary libraries
library(openxlsx)
library(dplyr)

# Define the parent folder for saving the combined Excel file
combined_file_path <- file.path(deseq2_dir, "V2DESeq2_Combined_Results_All_Clusters.xlsx")

# Create a new workbook to store all the clusters' data
wb_combined <- createWorkbook()

# Loop through each cluster, read the Excel files, and add them to the combined workbook
for (cluster in clusters) {
  cat("Adding results for cluster:", cluster, "to the combined workbook.\n")
  
  # Define the path to the individual cluster Excel file
  excel_file_path <- file.path(deseq2_dir, paste0("Cluster", cluster, "/V2DESeq2_Results_Cluster_", cluster, ".xlsx"))
  
  # Read the individual cluster's Excel file
  cluster_data <- read.xlsx(excel_file_path)
  
  # Add the data to a new sheet in the combined workbook
  sheet_name <- paste0("Cluster_", cluster)
  
  # Make sure the sheet name does not exceed 31 characters (Excel's limit)
  sheet_name <- substr(sheet_name, 1, 31)
  
  # Add a new worksheet with the cluster data
  addWorksheet(wb_combined, sheetName = sheet_name)
  writeData(wb_combined, sheet = sheet_name, cluster_data)
}

# Save the combined workbook in the parent directory
saveWorkbook(wb_combined, combined_file_path, overwrite = TRUE)
cat("Combined DESeq2 results saved to:", combined_file_path, "\n")





###verifying the number of counts in my object
# Load necessary libraries
library(Seurat)
library(dplyr)
library(tibble)

# Define the function
get_gene_counts_by_cluster <- function(seurat_obj, gene_of_interest, cluster_of_interest) {
  
  # Extract the counts matrix
  counts_matrix <- GetAssayData(seurat_obj, layer = "counts")
  
  # Check if the gene exists in the counts matrix
  if(gene_of_interest %in% rownames(counts_matrix)) {
    
    # Extract counts for the gene of interest
    gene_counts <- counts_matrix[gene_of_interest, ]
    
    # Extract metadata and convert rownames to a column
    meta_data <- seurat_obj@meta.data %>%
      rownames_to_column(var = "cell")
    
    # Combine counts and metadata into a single dataframe
    result <- data.frame(cell = colnames(counts_matrix), counts = gene_counts) %>%
      left_join(meta_data, by = "cell")
    
    # Filter by the defined cluster
    result_filtered <- result %>%
      filter(seurat_clusters == cluster_of_interest)
    
    # Summarize counts by orig.ident and sort
    counts_summary <- result_filtered %>%
      group_by(orig.ident) %>%
      summarise(total_counts = sum(counts)) %>%
      arrange(desc(total_counts))  # Sort in descending order of total counts
    
    # Return the result
    return(counts_summary)
    
  } else {
    cat("The gene '", gene_of_interest, "' is not found in the counts matrix.\n", sep = "")
  }
}

# Example call: Replace "fib_subset.pseudo" with your Seurat object
get_gene_counts_by_cluster(fib_subset.pseudo, "H2-Aa", 0)






###Differential gene experession bulk and single cell control vs experimental FIBS-----------

# Set the base path for saving the volcano plots
volcano_plot_dir <- file.path(epi_subset_path, "Volcano_Plots")

# Create the directory if it doesn't exist
if (!dir.exists(volcano_plot_dir)) {
  dir.create(volcano_plot_dir)
}


library(ggplot2)

# Ensure the correct identities are set based on the Seurat clusters
Idents(obj_original) <- "seurat_clusters"
Idents(obj.pseudo) <- "seurat_clusters"





head(obj_original@meta.data)
head(obj.pseudo@meta.data)
# Create a new column 'condition' based on 'orig.ident'
obj_original@meta.data$condition <- ifelse(grepl("^C", obj_original@meta.data$orig.ident), "Control", "Experimental")

head(obj_original@meta.data)
head(obj.pseudo@meta.data)




library(ggplot2)
library(ggrepel)
library(dplyr)
library(Seurat)
library(openxlsx)

# Initialize an empty list to store results for each cluster
results_list <- list()

# Loop through each cluster
clusters <- unique(Idents(obj_original))
for (cluster in clusters) {
  print(paste("Processing cluster", cluster))
  
  # Subset the Seurat objects for the current cluster
  single_cell_cluster <- subset(obj_original, idents = cluster)
  pseudo_bulk_cluster <- subset(obj.pseudo, idents = cluster)
  
  # Ensure the correct identities are set based on the condition
  Idents(single_cell_cluster) <- single_cell_cluster$condition
  Idents(pseudo_bulk_cluster) <- pseudo_bulk_cluster$condition
  
  # Perform differential expression analysis for single-cell data
  DE_single <- FindMarkers(single_cell_cluster, ident.1 = "Experimental", ident.2 = "Control", test.use = "wilcox", slot = "data")
  DE_single <- DE_single %>% dplyr::select(log2FC.single = avg_log2FC, p_val_adj.single = p_val_adj)
  DE_single$gene_id <- rownames(DE_single)
  
  # Perform differential expression analysis for pseudo-bulk data
  DE_bulk <- FindMarkers(pseudo_bulk_cluster, ident.1 = "Experimental", ident.2 = "Control", min.cells.group = 0, test.use = "DESeq2")
  DE_bulk <- DE_bulk %>% dplyr::select(log2FC.bulk = avg_log2FC, p_val_adj.bulk = p_val_adj)
  DE_bulk$gene_id <- rownames(DE_bulk)
  
  # Merge single-cell and pseudo-bulk results
  merged_results <- merge(DE_single, DE_bulk, by = "gene_id", all = TRUE)
  merged_results$cluster <- cluster
  
  # Store the results in the list
  results_list[[as.character(cluster)]] <- merged_results
  
  # Select top 25 upregulated and downregulated genes for labeling
  top_upregulated <- merged_results %>% filter(p_val_adj.bulk < 0.05) %>% arrange(desc(log2FC.bulk)) %>% head(25)
  top_downregulated <- merged_results %>% filter(p_val_adj.bulk < 0.05) %>% arrange(log2FC.bulk) %>% head(25)
  top_genes <- bind_rows(top_upregulated, top_downregulated)
  
  # Create a volcano plot for the cluster
  volcano_plot <- ggplot(merged_results, aes(x = log2FC.bulk, y = -log10(p_val_adj.bulk))) +
    geom_point(aes(color = ifelse(p_val_adj.bulk < 0.05 & log2FC.bulk > 0, "Upregulated",
                                  ifelse(p_val_adj.bulk < 0.05 & log2FC.bulk < 0, "Downregulated", "Not Significant")))) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
    labs(title = paste("Volcano Plot - Cluster", cluster),
         x = "Log2 Fold Change (Bulk)",
         y = "-Log10 Adjusted P-value",
         color = "Gene Regulation") +
    geom_text_repel(data = top_genes, aes(label = gene_id), size = 3, max.overlaps = 20) +
    theme_minimal()
  
  # Save the volcano plot
  ggsave(filename = file.path(volcano_plot_dir, paste0("Volcano_Cluster_", cluster, ".png")),
         plot = volcano_plot, device = "png", width = 8, height = 6)
  
  print(paste("Volcano plot saved for cluster", cluster))
}

# Combine all clusters into a single data frame
final_results <- bind_rows(results_list)

# Save the results to an Excel file
output_file <- file.path(epi_subset_path, "DE_Comparison_Control_vs_Experimental.xlsx")
write.xlsx(final_results, file = output_file, rowNames = FALSE)

cat("Results saved to:", output_file, "\n")








#TRY2 to include percents
# Load necessary libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(Seurat)
library(openxlsx)

# Create a directory for volcano plots if it doesn't exist
volcano_plot_dir <- file.path(epi_subset_path, "DESeq2/Volcano_Plots")
if (!dir.exists(volcano_plot_dir)) {
  dir.create(volcano_plot_dir, recursive = TRUE)
}

# Initialize an empty list to store results for each cluster
results_list <- list()

# Loop through each cluster
clusters <- unique(Idents(obj_original))
for (cluster in clusters) {
  print(paste("Processing cluster", cluster))
  
  # Subset the Seurat objects for the current cluster
  single_cell_cluster <- subset(obj_original, idents = cluster)
  pseudo_bulk_cluster <- subset(obj.pseudo, idents = cluster)
  
  # Ensure the correct identities are set based on the condition
  Idents(single_cell_cluster) <- single_cell_cluster$condition
  Idents(pseudo_bulk_cluster) <- pseudo_bulk_cluster$condition
  
  # Calculate the percentage of cells expressing each gene in Control and Experimental for single-cell data
  pct_express_single <- single_cell_cluster@assays$RNA@data %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene_id") %>%
    gather(key = "cell", value = "expression", -gene_id) %>%
    inner_join(as.data.frame(single_cell_cluster@meta.data) %>%
                 rownames_to_column(var = "cell") %>%
                 select(cell, condition), by = "cell") %>%
    group_by(gene_id, condition) %>%
    summarise(pct_express = sum(expression > 0) / n() * 100) %>%
    tidyr::pivot_wider(names_from = condition, values_from = pct_express) %>%
    rename(pct_express_control_single = Control, pct_express_experimental_single = Experimental)
  
  # Calculate the percentage of cells expressing each gene in Control and Experimental for pseudo-bulk data
  pct_express_bulk <- pseudo_bulk_cluster@assays$RNA@data %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene_id") %>%
    gather(key = "cell", value = "expression", -gene_id) %>%
    inner_join(as.data.frame(pseudo_bulk_cluster@meta.data) %>%
                 rownames_to_column(var = "cell") %>%
                 select(cell, condition), by = "cell") %>%
    group_by(gene_id, condition) %>%
    summarise(pct_express = sum(expression > 0) / n() * 100) %>%
    tidyr::pivot_wider(names_from = condition, values_from = pct_express) %>%
    rename(pct_express_control_bulk = Control, pct_express_experimental_bulk = Experimental)
  
  # Perform differential expression analysis for single-cell data
  DE_single <- FindMarkers(single_cell_cluster, ident.1 = "Experimental", ident.2 = "Control", test.use = "wilcox", slot = "data")
  DE_single <- DE_single %>% select(log2FC.single = avg_log2FC, p_val_adj.single = p_val_adj)
  DE_single$gene_id <- rownames(DE_single)
  
  # Perform differential expression analysis for pseudo-bulk data
  DE_bulk <- FindMarkers(pseudo_bulk_cluster, ident.1 = "Experimental", ident.2 = "Control", min.cells.group = 0, test.use = "DESeq2")
  DE_bulk <- DE_bulk %>% select(log2FC.bulk = avg_log2FC, p_val_adj.bulk = p_val_adj)
  DE_bulk$gene_id <- rownames(DE_bulk)
  
  # Merge single-cell and pseudo-bulk results with percentage data
  merged_results <- DE_single %>%
    full_join(DE_bulk, by = "gene_id") %>%
    full_join(pct_express_single, by = "gene_id") %>%
    full_join(pct_express_bulk, by = "gene_id")
  merged_results$cluster <- cluster
  
  # Store the results in the list
  results_list[[as.character(cluster)]] <- merged_results
  
  # Select top 25 upregulated and downregulated genes for labeling
  top_upregulated <- merged_results %>% filter(p_val_adj.bulk < 0.05) %>% arrange(desc(log2FC.bulk)) %>% head(25)
  top_downregulated <- merged_results %>% filter(p_val_adj.bulk < 0.05) %>% arrange(log2FC.bulk) %>% head(25)
  top_genes <- bind_rows(top_upregulated, top_downregulated)
  
  # Create a volcano plot for the cluster
  volcano_plot <- ggplot(merged_results, aes(x = log2FC.bulk, y = -log10(p_val_adj.bulk))) +
    geom_point(aes(color = ifelse(p_val_adj.bulk < 0.05 & log2FC.bulk > 0, "Upregulated",
                                  ifelse(p_val_adj.bulk < 0.05 & log2FC.bulk < 0, "Downregulated", "Not Significant")))) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
    labs(title = paste("Volcano Plot - Cluster", cluster),
         x = "Log2 Fold Change (Bulk)",
         y = "-Log10 Adjusted P-value",
         color = "Gene Regulation") +
    geom_text_repel(data = top_genes, aes(label = gene_id), size = 3, max.overlaps = 20) +
    theme_minimal()
  
  # Save the volcano plot
  ggsave(filename = file.path(volcano_plot_dir, paste0("Volcano_Cluster_", cluster, ".png")),
         plot = volcano_plot, device = "png", width = 8, height = 6)
  
  print(paste("Volcano plot saved for cluster", cluster))
}

# Combine all clusters into a single data frame
final_results <- bind_rows(results_list)

# Save the results to an Excel file
output_file <- file.path(fib_subset_path, "DE_Comparison_Control_vs_Experimental_with_Percent_Expression.xlsx")
write.xlsx(final_results, file = output_file, rowNames = FALSE)

cat("Results saved to:", output_file, "\n")


#### ptprc positive number in each cluster by condition-------------
# Define a threshold for PTPRC positivity (you may adjust this based on your data)
ptprc_threshold <- 0  # Assuming any expression greater than 0 is considered positive

# Identify PTPRC-positive cells
ptprc_positive_cells <- WhichCells(fib_subset, expression = Ptprc > ptprc_threshold)

# Subset the Seurat object to only include PTPRC-positive cells
fib_subset_ptprc_positive <- subset(fib_subset, cells = ptprc_positive_cells)

# Get the number of PTPRC-positive cells in each cluster, split by condition
ptprc_positive_counts_by_condition <- table(Idents(fib_subset_ptprc_positive), fib_subset_ptprc_positive$condition)

# Convert the result to a data frame for easy viewing
ptprc_positive_counts_by_condition_df <- as.data.frame(ptprc_positive_counts_by_condition)
colnames(ptprc_positive_counts_by_condition_df) <- c("Cluster", "Condition", "PTPRC_Positive_Cells")

# Print the result
print(ptprc_positive_counts_by_condition_df)


# Define the file path for saving
excel_file_path <- file.path(fib_subset_path, "Ptprc_Positive_Cells.xlsx")

# Save the result to an Excel file
write.xlsx(ptprc_positive_counts_by_condition_df, file = excel_file_path, rowNames = FALSE)

cat("Excel file saved to:", excel_file_path, "\n")






###statistically determining if cell number changes in my experiment-----------
# Load necessary libraries
library(readxl)
library(dplyr)

# Define file path
clust_counts_path <- "C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/clust_counts.xlsx"

# Read the Excel file
data <- read_excel(clust_counts_path)

# Initialize a list to store Fisher's Exact Test results
fishers_results <- list()

# Loop through each cluster and perform Fisher's Exact Test
for (i in 1:nrow(data)) {
  # Get control and experimental cell counts for the current cluster
  control_counts <- as.numeric(data[i, c("C3", "C4", "C5")])
  experimental_counts <- as.numeric(data[i, c("E3", "E4", "E5")])
  
  # Create a contingency table
  contingency_table <- matrix(c(sum(control_counts), sum(experimental_counts), 
                                sum(data$Total_Control[i]) - sum(control_counts), 
                                sum(data$Total_Experimental[i]) - sum(experimental_counts)), 
                              nrow = 2, byrow = TRUE)
  
  # Perform Fisher's Exact Test
  fisher_test <- fisher.test(contingency_table)
  
  # Store the p-value and cluster information
  fishers_results[[i]] <- list(cluster = data$seurat_clusters[i], p_value = fisher_test$p.value)
}

# Convert the results to a data frame
fishers_results_df <- do.call(rbind, lapply(fishers_results, as.data.frame))
colnames(fishers_results_df) <- c("Cluster", "P_Value")

# View the results
print(fishers_results_df)

# Save the results to a CSV file
write.csv(fishers_results_df, "C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/fishers_exact_results.csv", row.names = FALSE)

###Trying to figure out what marker genes i could use for Fib cell types-----------------
# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

fib_subset_path <- "C:\\Users\\ostrobel\\Indiana University\\Jerde, Travis J - JERDE LAB (J Drive)\\Strobel\\Single Cell Analysis Stuff\\R_Analysis_Strobel\\Combined_Sample_Analysis\\Combined_analysis\\combined_obj_try4\\Harmony_SCT\\Fib_subset\\2024_10_07_fib_subset_FINAL.rds"

# fib_subset_pseudo_path <- "C:\\Users\\ostrobel\\Indiana University\\Jerde, Travis J - JERDE LAB (J Drive)\\Strobel\\Single Cell Analysis Stuff\\R_Analysis_Strobel\\Combined_Sample_Analysis\\Combined_analysis\\combined_obj_try4\\Harmony_SCT\\Fib_subset\\2024_10_07_fib_sub_psudo_FINAL.rds"

fib_subset <- readRDS(fib_subset_path)

#loading in my object 
seurat_obj <- readRDS("C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/2024_09_20_FINAL_SUERAT_OBJECT.rds")



# Define the list of marker genes
#cluster0 markers
genes_of_interest <- c("Fn1", "Metrnl", "Ebf2", "Scara5", "Ifi202b", "Sult1e1", "Gfpt2", "Ccl11", "Ebf1", "Gpc3", 
                       "Dpep1", "Adamts5", "Tnfaip2", "C3", "Slfn5", "Pcolce2")


#cluster3 markers
genes_of_interest <- c("Ccl5", "Cd52", "Hcst", "Nkg7", "Ptpn18", "Cd3g", "Fcer1g", "Trbc2", "Cd3d", "Ms4a4b", "Ptprc", "Rac2", "Arhgdib", "AW112010", "Tyrobp", "Ltb")


genes_of_interest <- c("Dcn")

#cluster0 DE
genes_of_interest <- c("Cxcl10", "Isg15", "Gbp2", "Iigp1", "Ifit1", "Bst2", "Ly6e", "Gbp3", 
                       "H2-K1", "Plac8", "H2-D1", "Psmb8", "Tnfaip2", "B2m", "Ifi27l2a", 
                       "Socs1", "H2-T23", "Adamts5", "Fbn1", "Ctla2a", "Irgm1", "Cd248", 
                       "Ebpl", "Pi16", "Cd55", "Lbp", "Col14a1")


genes_of_interest <- c("H2-Aa", "H2-Ab1", "Cxcl9", "Gm20513", "Cd74", "Gbp8", "H2-DMb1", "Gbp10", "Gbp4", "Il18bp", "H2-Eb1", "Gbp6", "Ifi47", "Cd274", "Zbp1", "H2-DMa", "Ifi44", "Olfr56", "Rsad2", "Cxcl10", "Igtp")

genes_of_interest <- c("Psmb10", "Psmb8", "H2-Aa") 



# Define the save path
save_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Fib_dimplots_potential_markers"

# Loop through each gene in the list
for (gene_of_interest in genes_of_interest) {
  
  # Create the FeaturePlot for the full Seurat object, split by condition
  plot1 <- FeaturePlot(seurat_obj, features = gene_of_interest, split.by = "condition", reduction = "umap.harmonyintegration") + 
    ggtitle(paste0(gene_of_interest, " - Full Seurat Experimental"))
  
  # Create the FeaturePlot for the fibroblast subset, split by condition
  plot2 <- FeaturePlot(fib_subset, features = gene_of_interest, split.by = "condition", reduction = "umap") + 
    ggtitle(paste0(gene_of_interest, " - Fibroblast Subset Experimental"))
  
  # Combine the two plots into a 2x2 grid
  combined_plot <- wrap_plots(plot1, plot2, ncol = 1)
  
  # Save the combined plot as a PNG file
  file_name <- file.path(save_path, paste0(gene_of_interest, "_plot.png"))  # Ensure .png extension is included
  ggsave(filename = file_name, plot = combined_plot, width = 10, height = 10)  # 2 plots, each 5x5 inches
  
  # Output to confirm file saved
  cat("Saved plot for", gene_of_interest, "to", file_name, "\n")
}





FeaturePlot(seurat_obj, features = "Svs1", reduction = "umap.harmonyintegration", split.by = "condition")
FeaturePlot(fib_subset, features = "Svs1", reduction = "umap", split.by = "condition")

FeaturePlot(seurat_obj, features = "Lsp1", reduction = "umap.harmonyintegration", split.by = "condition")
FeaturePlot(fib_subset, features = "Lsp1", reduction = "umap", split.by = "condition")

FeaturePlot(seurat_obj, features = "Gm2a", reduction = "umap.harmonyintegration", split.by = "condition")
FeaturePlot(fib_subset, features = "Gm2a", reduction = "umap", split.by = "condition")

FeaturePlot(seurat_obj, features = "Cd74", reduction = "umap.harmonyintegration", split.by = "condition")
FeaturePlot(fib_subset, features = "Cd74", reduction = "umap", split.by = "condition")

FeaturePlot(seurat_obj, features = "H2-Aa", reduction = "umap.harmonyintegration", split.by = "condition")
FeaturePlot(fib_subset, features = "H2-Aa", reduction = "umap", split.by = "condition")





Reductions(fib_subset)





# Define the base save path
base_save_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Fib_dimplots_potential_markers/Cluster0vs1"

# Create the folders for overlap, cluster0 specific, and cluster3 specific
overlap_dir <- file.path(base_save_path, "overlap")
cluster0_dir <- file.path(base_save_path, "cluster0_specific")
cluster3_dir <- file.path(base_save_path, "cluster3_specific")

dirs <- c(overlap_dir, cluster0_dir, cluster3_dir)
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# Define the three gene lists
overlap_genes <- c("Cxcl10", "H2-K1", "Psmb8", "Ifit1", "Gbp2", "Iigp1", "B2m", "Gbp3", "H2-D1", "Isg15", "H2-T23", "Bst2", "Ly6e", "Plac8", "Socs1")
cluster0_specific <- c("Psmb9", "Tap1", "C2", "Gbp5", "AW112010", "Srgn", "Tap2", "Irgm2", "Psmb10", "Ifih1", "Gm4951", "Tmem140", "H2-Q4")
cluster3_specific <- c("Tnfaip2", "Ifi27l2a", "Adamts5", "Fbn1", "Ctla2a", "Cd248", "Ebpl", "Pi16", "Cd55", "Lbp", "Col14a1", "Tgfbr2")

# Function to generate and save FeaturePlots for a list of genes
generate_and_save_plots <- function(genes, save_dir) {
  for (gene_of_interest in genes) {
    
    # Create the FeaturePlot for the full Seurat object, split by condition
    plot1 <- FeaturePlot(seurat_obj, features = gene_of_interest, split.by = "condition", reduction = "umap.harmonyintegration") + 
      ggtitle(paste0(gene_of_interest, " - Full Seurat Experimental"))
    
    # Create the FeaturePlot for the fibroblast subset, split by condition
    plot2 <- FeaturePlot(fib_subset, features = gene_of_interest, split.by = "condition", reduction = "umap") + 
      ggtitle(paste0(gene_of_interest, " - Fibroblast Subset Experimental"))
    
    # Combine the two plots into a 2x1 grid
    combined_plot <- wrap_plots(plot1, plot2, ncol = 1)
    
    # Save the combined plot as a PNG file
    file_name <- file.path(save_dir, paste0(gene_of_interest, "_plot.png"))  # Ensure .png extension is included
    ggsave(filename = file_name, plot = combined_plot, width = 10, height = 10)  # 2 plots, each 5x5 inches
    
    # Output to confirm file saved
    cat("Saved plot for", gene_of_interest, "to", file_name, "\n")
  }
}

# Run the function for each list
generate_and_save_plots(overlap_genes, overlap_dir)
generate_and_save_plots(cluster0_specific, cluster0_dir)
generate_and_save_plots(cluster3_specific, cluster3_dir)

cat("All plots generated and saved successfully.\n")


FeaturePlot(seurat_obj, features = "H2-Aa", reduction = "umap.harmonyintegration", split.by = "condition")
FeaturePlot(fib_subset, features = "H2-Aa", reduction = "umap", split.by = "condition")



fib_subset2 <- PrepSCTFindMarkers(fib_subset)

fib_markers <- FindAllMarkers(fib_subset2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

seurat_obj2 <- PrepSCTFindMarkers(seurat_obj)
all_markers <- FindAllMarkers(seurat_obj2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)



# View the different SCT models applied
seurat_obj[["SCT"]]@SCTModel.list


head(all_markers)



###DE gene expression ALL CELLS findmarkers and cont vs inf-----------------
library(DESeq2)
library(Seurat)
library(ggplot2)
library(presto)
library(dplyr)
library(openxlsx)

seurat_obj <- readRDS("C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/2024_09_20_FINAL_SUERAT_OBJECT.rds")


# Set the default assay to RNA to focus on the original data
DefaultAssay(seurat_obj) <- "RNA"

# Aggregate expression data by original identity, condition, and Seurat clusters
seurat_obj.pseudo <- AggregateExpression(seurat_obj, return.seurat=TRUE, group.by=c("orig.ident", "condition", "seurat_clusters"))
View(seurat_obj.pseudo@meta.data)

# Set identifiers for the new pseudo object
Idents(seurat_obj.pseudo) <- "seurat_clusters"


# Define the save path
save_path <- "C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/2024_10_11_FINAL_ALL_pseudo.rds"


output_path <- "C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/"



# Save the Seurat object
saveRDS(seurat_obj.pseudo, file = save_path)


# Define the original and pseudo objects for marker finding
obj_original <- seurat_obj
obj.pseudo <- seurat_obj.pseudo



# Set identities based on Seurat clusters
Idents(obj_original) <- "seurat_clusters"
Idents(obj.pseudo) <- "seurat_clusters"

# Set the RNA assay as the default assay
DefaultAssay(object = obj_original) <- "RNA"
DefaultAssay(object = obj.pseudo) <- "RNA"

obj_original <- JoinLayers(obj_original, assay = "RNA")
obj.pseudo <- JoinLayers(obj.pseudo, assay = "RNA")

unique_clusters <- unique(Idents(seurat_obj))
print(unique_clusters)

# Loop through clusters to find markers, comparing single cell and pseudo-bulk data
n <- 24 #this is the number of your LAST cluster. Pointer makes it so that we start at 0
keep <- data.frame()
wb_cluster_DE <- createWorkbook()
ptr <- 0

for(i in 0:n){
  print(paste("Processing cluster", i))
  
  DE_temp <- FindMarkers(obj_original, ident.1=paste(i), test.use="wilcox", slot="data")
  DE_temp <- DE_temp[c(2:5)]
  colnames(DE_temp) <- c("log2FC.single", "pct.current_cluster", "pct.other_clusters", "padj.single")
  DE_temp$gene_id <- rownames(DE_temp)
  
  bulk_temp <- FindMarkers(object=obj.pseudo, ident.1=paste(i), min.cells.group=0, test.use="DESeq2")
  bulk_temp <- bulk_temp[c(2,5)]
  colnames(bulk_temp) <- c("log2FC.bulk","padj.bulk")
  bulk_temp$gene_id <- rownames(bulk_temp)
  
  temp <- merge(DE_temp,bulk_temp,by="gene_id", all = TRUE)
  
  temp2 <- temp
  temp2$cluster <- paste(i)
  temp.sig.either <- subset(temp2, padj.single < 0.05 | padj.bulk < 0.05)
  temp.sig.both <- subset(temp2, padj.single < 0.05 & padj.bulk < 0.05)
  
  keep <- rbind(keep, temp.sig.both)
  
  addWorksheet(wb_cluster_DE, paste("cluster", i, "all", sep="."))
  addWorksheet(wb_cluster_DE, paste("cluster", i, "sig", "either", sep="."))
  addWorksheet(wb_cluster_DE, paste("cluster", i, "sig", "both", sep="."))
  
  writeData(wb_cluster_DE, sheet=(ptr+1), temp, rowNames=FALSE)
  writeData(wb_cluster_DE, sheet=(ptr+2), temp.sig.either, rowNames=FALSE)
  writeData(wb_cluster_DE, sheet=(ptr+3), temp.sig.both, rowNames=FALSE)
  
  ptr <- ptr+3
}

# Save the workbook summarizing differential expression results
saveWorkbook(wb_cluster_DE, file.path(output_path, "2024_10_10_ALL_Markers.xlsx"), overwrite = TRUE)


# Export a summary of markers for manual annotation
keep.pos <- subset(keep, log2FC.single > 0 & log2FC.bulk > 0)
id <- data.frame("gene_id"=keep.pos[duplicated(keep.pos[,1]),1])
keep.pos.unique <- dplyr::anti_join(keep.pos, id, by="gene_id")
keep.pos.overlapping <- merge(keep.pos, id, by="gene_id")

wb_cluster_summary <- createWorkbook()
addWorksheet(wb_cluster_summary, "all_markers")
addWorksheet(wb_cluster_summary, "all_positive_markers")
addWorksheet(wb_cluster_summary, "all_positive_markers_overlap")
addWorksheet(wb_cluster_summary, "all_positive_markers_unique")
writeData(wb_cluster_summary, sheet=1, keep, rowNames=FALSE)
writeData(wb_cluster_summary, sheet=2, keep.pos, rowNames=FALSE)
writeData(wb_cluster_summary, sheet=3, keep.pos.overlapping, rowNames=FALSE)
writeData(wb_cluster_summary, sheet=4, keep.pos.unique, rowNames=FALSE)

# Save the workbook summarizing differential expression results
saveWorkbook(wb_cluster_DE, file.path(output_path, "2024_10_10_ALL_markers_for_man_annot.xlsx"), overwrite = TRUE)







# Set the base path for saving the volcano plots
volcano_plot_dir <- "C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/VolPlots"

library(ggplot2)

# Ensure the correct identities are set based on the Seurat clusters
Idents(obj_original) <- "seurat_clusters"
Idents(obj.pseudo) <- "seurat_clusters"





head(obj_original@meta.data)
head(obj.pseudo@meta.data)
# Create a new column 'condition' based on 'orig.ident'
obj_original@meta.data$condition <- ifelse(grepl("^C", obj_original@meta.data$orig.ident), "Control", "Experimental")

head(obj_original@meta.data)
head(obj.pseudo@meta.data)




library(ggplot2)
library(ggrepel)


# Initialize an empty list to store results for each cluster
results_list <- list()

# Loop through each cluster
clusters <- unique(Idents(obj_original))
for (cluster in clusters) {
  print(paste("Processing cluster", cluster))
  
  # Subset the Seurat objects for the current cluster
  single_cell_cluster <- subset(obj_original, idents = cluster)
  pseudo_bulk_cluster <- subset(obj.pseudo, idents = cluster)
  
  # Ensure the correct identities are set based on the condition
  Idents(single_cell_cluster) <- single_cell_cluster$condition
  Idents(pseudo_bulk_cluster) <- pseudo_bulk_cluster$condition
  
  # Perform differential expression analysis for single-cell data
  DE_single <- FindMarkers(single_cell_cluster, ident.1 = "Experimental", ident.2 = "Control", test.use = "wilcox", slot = "data")
  DE_single <- DE_single %>% select(log2FC.single = avg_log2FC, p_val_adj.single = p_val_adj)
  DE_single$gene_id <- rownames(DE_single)
  
  # Perform differential expression analysis for pseudo-bulk data
  DE_bulk <- FindMarkers(pseudo_bulk_cluster, ident.1 = "Experimental", ident.2 = "Control", min.cells.group = 0, test.use = "DESeq2")
  DE_bulk <- DE_bulk %>% select(log2FC.bulk = avg_log2FC, p_val_adj.bulk = p_val_adj)
  DE_bulk$gene_id <- rownames(DE_bulk)
  
  # Merge single-cell and pseudo-bulk results
  merged_results <- merge(DE_single, DE_bulk, by = "gene_id", all = TRUE)
  merged_results$cluster <- cluster
  
  # Store the results in the list
  results_list[[as.character(cluster)]] <- merged_results
  
  # Select top 25 upregulated and downregulated genes for labeling
  top_upregulated <- merged_results %>% filter(p_val_adj.bulk < 0.05) %>% arrange(desc(log2FC.bulk)) %>% head(25)
  top_downregulated <- merged_results %>% filter(p_val_adj.bulk < 0.05) %>% arrange(log2FC.bulk) %>% head(25)
  top_genes <- bind_rows(top_upregulated, top_downregulated)
  
  # Create a volcano plot for the cluster
  volcano_plot <- ggplot(merged_results, aes(x = log2FC.bulk, y = -log10(p_val_adj.bulk))) +
    geom_point(aes(color = ifelse(p_val_adj.bulk < 0.05 & log2FC.bulk > 0, "Upregulated",
                                  ifelse(p_val_adj.bulk < 0.05 & log2FC.bulk < 0, "Downregulated", "Not Significant")))) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
    labs(title = paste("Volcano Plot - Cluster", cluster),
         x = "Log2 Fold Change (Bulk)",
         y = "-Log10 Adjusted P-value",
         color = "Gene Regulation") +
    geom_text_repel(data = top_genes, aes(label = gene_id), size = 3, max.overlaps = 20) +
    theme_minimal()
  
  # Save the volcano plot
  ggsave(filename = file.path(volcano_plot_dir, paste0("Volcano_Cluster_", cluster, ".png")),
         plot = volcano_plot, device = "png", width = 8, height = 6)
  
  print(paste("Volcano plot saved for cluster", cluster))
}

# Combine all clusters into a single data frame
final_results <- bind_rows(results_list)

# Save the results to an Excel file
output_file <- file.path(output_path, "2024_10_11_ALL_Cont_vs_Exp.xlsx")
write.xlsx(final_results, file = output_file, rowNames = FALSE)

cat("Results saved to:", output_file, "\n")


FeaturePlot(fib_subset, features = "Cd3g", reduction = "umap", split.by = "condition")



###comparing strand BPH DEG to Mouse Fib DEG;;;;;   making both stromal and all cells feature plots!!!!!!!!!!----------------
#i manually coverted teh genes. https://www.informatics.jax.org/marker/MGI:99450

#This is my list of human overlap: CXCL9,CD74,GBP6,IL18BP,HLA-DMA,RSAD2,CXCL10,CD274,ZBP1,PSMB8,IFIT1
#mouse
# Cxcl9
# Cd74
# Gbp6
# Il18bp
# H2-DMa
# Rsad2
# Cxcl10
# Cd274
# Zbp1
# Psmb8
# Ifit1


fib_subset_path <- "C:\\Users\\ostrobel\\Indiana University\\Jerde, Travis J - JERDE LAB (J Drive)\\Strobel\\Single Cell Analysis Stuff\\R_Analysis_Strobel\\Combined_Sample_Analysis\\Combined_analysis\\combined_obj_try4\\Harmony_SCT\\Fib_subset\\2024_10_07_fib_subset_FINAL.rds"


fib_subset <- readRDS(fib_subset_path)

#loading in my object 
seurat_obj <- readRDS("C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/2024_09_20_FINAL_SUERAT_OBJECT.rds")


#Now time to make featureplots
#genes_of_interest <- c("Cxcl9", "Cd74", "Gbp6", "Il18bp", "H2-DMa", "Rsad2", "Cxcl10", "Cd274", "Zbp1", "Psmb8", "Ifit1")

genes_of_interest <- c("Vim")

genes_of_interest <- c("Fn1", "Metrnl", "Ebf2", "Scara5", "Ifi202b", "Sult1e1", "Gfpt2", "Ccl11", "Ebf1", "Gpc3")



# Define the save path
#save_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Fib_dimplots_potential_markers/Overlap_BPH_Mouse"
save_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Fib_dimplots_potential_markers"



# Check if the directory exists, if not create it
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

# Now you can use save_path to save files in that directory

# Loop through each gene in the list
for (gene_of_interest in genes_of_interest) {
  
  # Create the FeaturePlot for the full Seurat object, split by condition
  plot1 <- FeaturePlot(seurat_obj, features = gene_of_interest, split.by = "condition", reduction = "umap.harmonyintegration") + 
    ggtitle(paste0(gene_of_interest, " - Full Seurat Experimental"))
  
  # Create the FeaturePlot for the fibroblast subset, split by condition
  plot2 <- FeaturePlot(fib_subset, features = gene_of_interest, split.by = "condition", reduction = "umap") + 
    ggtitle(paste0(gene_of_interest, " - Fibroblast Subset Experimental"))
  
  # Combine the two plots into a 2x2 grid
  combined_plot <- wrap_plots(plot1, plot2, ncol = 1)
  
  # Save the combined plot as a PNG file
  file_name <- file.path(save_path, paste0(gene_of_interest, "_plot.png"))  # Ensure .png extension is included
  ggsave(filename = file_name, plot = combined_plot, width = 10, height = 10)  # 2 plots, each 5x5 inches
  
  # Output to confirm file saved
  cat("Saved plot for", gene_of_interest, "to", file_name, "\n")
}


#making an excel where i compile the information of the two files
# Load necessary libraries
library(readxl)
library(dplyr)
library(openxlsx)

# Define the file paths to your Excel sheets
file1_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Strand_Human/20240911_DESeq2_Results_Strand_Fib.xlsx"
file2_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/DE_excel_docs_for_analysis/2024_10_10_DESeq2_results_fibs.xlsx"

# Define the output file path (same as the previous one)
output_path <- "C:/Users/ostrobel/OneDrive - Indiana University/Desktop/Strand_Human/Filtered_DESeq2_Results_Strand_Fib.xlsx"

# Define the list of genes to find
genes_of_interest <- c("Cxcl9", "Cd74", "Gbp6", "Il18bp", "H2-DMa", "Rsad2", "Cxcl10", "Cd274", "Zbp1", "Psmb8", "Ifit1","CXCL9", "CD74", "GBP6", "IL18BP", "HLA-DMA", "RSAD2", "CXCL10", "CD274", "ZBP1", "PSMB8", "IFIT1")

# Read the Excel files
df1 <- read_excel(file1_path)
df2 <- read_excel(file2_path)

# Filter both DataFrames to only include rows with those genes
filtered_df1 <- df1 %>% filter(gene_id %in% genes_of_interest)
filtered_df2 <- df2 %>% filter(gene %in% genes_of_interest)

# Combine both filtered DataFrames
combined_filtered_df <- bind_rows(filtered_df1, filtered_df2)


# Check if output file exists
if (file.exists(output_path)) {
  # Load existing workbook and get sheet names
  wb <- loadWorkbook(output_path)
  sheet_names <- getSheetNames(output_path)
  
  # If the "Filtered_DESeq2_Results" sheet does not exist, create it
  if (!"Filtered_DESeq2_Results" %in% sheet_names) {
    addWorksheet(wb, "Filtered_DESeq2_Results")
  }
} else {
  # Create a new workbook and add the sheet
  wb <- createWorkbook()
  addWorksheet(wb, "Filtered_DESeq2_Results")
}

# Write the combined filtered data to the existing or new sheet
writeData(wb, "Filtered_DESeq2_Results", combined_filtered_df, startRow = 1, colNames = TRUE)

# Save the updated or new Excel file
saveWorkbook(wb, output_path, overwrite = TRUE)

# Output to confirm file saved
cat("Combined filtered data saved to:", output_path, "\n")









# Example character strings
string1 <- c("IGKC", "CXCL13", "C1QTNF3", "IGHG1", "PNMT", "NRK", "KCNQ2", "PSRC1", "IL1B", "KRT222", "PCSK1N", "ENTPD3", "FNDC5", "CCL4", "CRYGD", "PGM5-AS1", "HAPLN1", "IFNG", "RANBP3L", "CCL4L2", "CALCA", "CHGB", "TAF7L", "ENHO", "ODAD2", "IGHM", "GZMK", "GRIK1", "RGS13", "CD48", "RASGRF1", "SLC7A3", "PLEK", "NELL2", "RRAD", "NKG7", "RNF112", "RGS1", "FAM171A2", "HSPA1L", "PAGE4", "F8", "PRKCB", "FAM43B", "CCR7", "CCL3L3", "KAZALD1", "EXOC3L1", "LINC01239", "CPA3", "DHRS2", "FAM13C", "TSBP1-AS1", "ATP1A2", "LRRC17", "MS4A6A", "POSTN", "LINC00882", "CLGN", "TAGAP", "MT1M", "B3GALT2", "RAC2", "ZFYVE28", "HLA-DRA", "SCT", "C1QL1", "CD2", "LCP1", "TGFBR3L", "CCL3", "SSX2IP", "NELL1", "HCST", "CDH10", "HMGCLL1", "GSTCD", "C1QC", "EMID1", "HRC", "CD53", "KCNH2", "CD69", "SRRM4", "EAF2", "ANGPT1", "GPR65", "PENK", "HSPA6", "TRHDE-AS1", "SMOC1", "TNFSF4", "HSD11B2", "MPP1", "RGS2", "PCLO", "ADCY5", "TPSB2", "ICAM4", "THBS1-IT1", "ANKRD36", "HLA-DPA1", "CD37", "FYB1", "MVK", "CMTM5", "EGR2", "LINC02814", "FCGR2A", "LRRC34", "JAKMIP2", "EPCAM-DT", "NAALAD2", "NPW", "CADPS", "RUNX3", "HLA-DQA1", "PLA2G7", "SPEG", "PPP1R1B", "RSPO2", "CILK1", "NDC80", "ERMARD", "PCDHB16", "FRZB", "LINC00941", "FILIP1L", "PGM5", "RGS7BP", "LINC02147", "BTG2", "PARPBP", "MAGI2", "MMP25-AS1", "CYB5D1", "MT1A", "TRUB1", "TSPAN18", "BMP5", "BCL2A1", "APOL4", "EPHX2", "MDK", "GASAL1", "KIAA0408", "ENPP2", "SORT1", "IKZF4", "TRBC1", "IGF1", "MAN1C1", "CST7", "RASL11A", "ADH1C", "CRISPLD1", "IDNK", "ADRA2A", "ROGDI", "PTGIS", "CPPED1", "HAS1", "ABCG2", "ITGA4", "HIVEP2-DT", "CA8", "CHRDL1", "ALOX5AP", "ZNF618", "LY6G6D", "LINC00877", "PGAP2", "NFKBID", "ZNF571", "FILIP1", "MYLK", "HELB", "CNN1", "CEP78", "ZNF460", "SRD5A2", "GSTM2", "QPRT", "IER5L", "RCAN2", "MAL", "GATA6", "CLIC3", "CHAF1A", "PDE4C", "SEMA3E", "GATA6-AS1", "FRRS1L", "PRLR", "PREX2", "UBA7", "ETV1", "HLF", "TMEM88", "NUDT10", "PEBP4", "SAMD12", "GPR183", "ZNF891", "TSPAN33", "CCND2-AS1", "OLFML1", "NAF1", "LRRC4", "JUND", "IFIT1", "RND1", "CCDC8", "H19", "TMEM25", "PARP10", "BTG3-AS1", "CNTN1", "COL1A1", "MNS1", "LINC00205", "APPAT", "PAK3", "KCNG1", "TICAM2-AS1", "P2RX2", "COX4I2", "IFI44L", "IL10RA", "NCAM1", "HLA-DRB1", "DES", "EZH2", "ISL1-DT", "F10", "FAM184A", "CROCC", "NME5", "IFI44", "GXYLT2", "RHOJ", "COL23A1", "RBFOX1", "LAPTM5", "ANKRD37", "SAMD11", "IRF1", "TRIM69", "IRAG1", "PPM1L", "PDE5A", "NUSAP1", "COL4A6", "MXRA5", "ZNF331", "CADPS2", "SCRG1", "EDIL3", "PLCL1", "CNTN4", "RASL12", "KDELR3", "GSTM5", "RIGI", "ROBO2", "ANKRD29", "PRECSIT", "FAM222A", "SAMSN1", "ITPR1", "ITM2A", "LTC4S", "DUSP2", "KIAA0586", "PCDH10", "TCEAL7", "ANK3", "SAP30-DT", "KIAA0040", "JUNB", "PRRT2", "RAMP2", "SEMA6A", "ZNF10", "CEP152", "ZNF141", "EIF4A3", "COL13A1", "LUM", "MYH10", "ITIH5", "MMP23B", "ATP1B2", "MAP9", "FOXC1", "RBP1", "GADD45B", "BATF3", "IFI6", "SGMS1-AS1", "SHTN1", "COPZ2", "LINC01082", "COL14A1", "ALDH7A1", "CPM", "SPOCK3", "PTGER2", "NFIL3", "CCDC34", "HIF3A", "TBL1X", "CMBL", "MAGED2", "TNFRSF10D", "IL11RA", "SLC66A2", "C4orf33", "MR1", "ZNF711", "MIR3936HG", "ADAM22", "NR4A3", "CNTFR", "SGCD", "RCAN1", "MEG9", "GAMT", "HIC1", "GFOD1", "GKAP1", "HSPB6", "FZD7"
)
string2 <- c("HLA-DQA1", "HLA-DQB1", "CXCL9", "CD74", "GBP6", "HLA-DMB", "GBP6", "GBP6", "IL18BP", "HLA-DRB1", "GBP6", "CD274", "ZBP1", "HLA-DMA", "IFI44", "OR2V1", "RSAD2", "CXCL10", "IRGM", "HLA-A", "GBP6", "PSMB8", "SLFN12", "IFIT1"
)

# Find the overlap (common elements)
common_elements <- intersect(string1, string2)

# Find elements unique to string1
unique_to_string1 <- setdiff(string1, string2)

# Find elements unique to string2
unique_to_string2 <- setdiff(string2, string1)

# Print results
cat("Common elements:\n", common_elements, "\n\n")
cat("Unique to string1:\n", unique_to_string1, "\n\n")
cat("Unique to string2:\n", unique_to_string2, "\n\n")









###looking at STRAND Fib data

FeaturePlot(fib_subset, features = "Vim", reduction = "umap")
FeaturePlot(fib_subset, features = "Dcn", reduction = "umap")
FeaturePlot(fib_subset, features = "Acta2", reduction = "umap")
FeaturePlot(fib_subset, features = "Pecam", reduction = "umap")
FeaturePlot(fib_subset, features = "C3", reduction = "umap")
FeaturePlot(fib_subset, features = "Lgr5", reduction = "umap")
FeaturePlot(fib_subset, features = "Wnt2", reduction = "umap")
FeaturePlot(fib_subset, features = "Uchl1", reduction = "umap")
FeaturePlot(fib_subset, features = "Sca-1", reduction = "umap")
FeaturePlot(fib_subset, features = "Cxcl13", reduction = "umap")
FeaturePlot(fib_subset, features = "Cxcl9", reduction = "umap")
FeaturePlot(fib_subset, features = "Cxcl10", reduction = "umap")


# Define gene groups
prostate_fibs_genes <- c("C3", "Ebf1", "Gpx3", "Sult1e1", "Igf1")
urethral_fibs_genes <- c("Lgr5", "Apoe", "Osr1", "Sfrp2", "Mfap4")
ductal_fibs_genes <- c("Wnt2", "Rorb", "Wif1", "Ifitm1", "Srd5a2")

# Function to generate FeaturePlots
generate_featureplots <- function(gene_list, group_name) {
  plots <- lapply(gene_list, function(gene) {
    FeaturePlot(fib_subset, features = gene, reduction = "umap") + 
      ggtitle(paste(group_name, "-", gene))  # Add group and gene name in title
  })
  
  combined_plot <- patchwork::wrap_plots(plots, ncol = 3)  # Combine the plots in a 3x2 grid
  return(combined_plot)
}

# Generate FeaturePlots for each group
prostate_fibs_plot <- generate_featureplots(prostate_fibs_genes, "Prostate Fibroblasts")
urethral_fibs_plot <- generate_featureplots(urethral_fibs_genes, "Urethral Fibroblasts")
ductal_fibs_plot <- generate_featureplots(ductal_fibs_genes, "Ductal Fibroblasts")

# Display the plots
print(prostate_fibs_plot)
print(urethral_fibs_plot)
print(ductal_fibs_plot)

# Save the plots to the Figure2 folder
ggsave(file.path(fig_2_path, "Prostate_Fibroblasts_FeaturePlots.png"), plot = prostate_fibs_plot, width = 12, height = 8)
ggsave(file.path(fig_2_path, "Urethral_Fibroblasts_FeaturePlots.png"), plot = urethral_fibs_plot, width = 12, height = 8)
ggsave(file.path(fig_2_path, "Ductal_Fibroblasts_FeaturePlots.png"), plot = ductal_fibs_plot, width = 12, height = 8)





#peri-epi fibs
#APOD, PTGDS, PTGS2 and MMP2 

#ifib
#C7, CCK, PCOLCE2 and GSN






# Load required libraries
library(readxl)
library(openxlsx)
library(dplyr)

# Define the file path to read in the Excel file
file_path <- "C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/2024_10_10_ALL_markers_for_man_annot.xlsx"

# Load the names of all sheets
sheet_names <- excel_sheets(file_path)

# Split sheets into cluster-related and non-cluster-related
cluster_sheets <- grep("^cluster\\.[0-9]+\\.sig\\.either$", sheet_names, value = TRUE)
non_cluster_sheets <- setdiff(sheet_names, cluster_sheets)

# Sort cluster sheets in ascending order based on the number
cluster_sheets <- cluster_sheets[order(as.numeric(gsub("[^0-9]", "", cluster_sheets)))]

# Combine sorted cluster sheets and non-cluster sheets
sorted_sheets <- c(cluster_sheets, sort(non_cluster_sheets))

# Create a new workbook to save the modified sheets
wb <- createWorkbook()

# Define a basic table style to format the table in each sheet
table_style <- createStyle(
  fontName = "Arial", fontSize = 10, border = "TopBottomLeftRight",
  fgFill = "#D3D3D3", halign = "center", textDecoration = "Bold"
)

# Loop through each sorted sheet
for (sheet in sorted_sheets) {
  # Read the sheet into a data frame
  df <- read_excel(file_path, sheet = sheet)
  
  # Check if "pct.diff" column exists, if not, create it
  if (!"pct.diff" %in% colnames(df)) {
    if ("pct.current_cluster" %in% colnames(df) & "pct.other_clusters" %in% colnames(df)) {
      df$pct.diff <- df$pct.current_cluster - df$pct.other_clusters
    } else {
      warning(paste("The sheet", sheet, "does not have the necessary columns 'pct.current_cluster' and 'pct.other_clusters'."))
      next
    }
  }
  
  # Sort by pct.diff column, largest to smallest
  df <- df %>%
    arrange(desc(pct.diff))
  
  # Add the sheet back to the workbook
  addWorksheet(wb, sheet)
  
  # Write the sorted table into the sheet
  writeData(wb, sheet = sheet, df)
  
  # Insert a formatted table
  addStyle(wb, sheet = sheet, style = table_style, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
  
  # Add auto filter for the table
  addFilter(wb, sheet = sheet, row = 1, cols = 1:ncol(df))
}

# Define the output path
output_path <- "C:/Users/ostrobel/Indiana University/Jerde, Travis J - JERDE LAB (J Drive)/Strobel/Single Cell Analysis Stuff/R_Analysis_Strobel/Combined_Sample_Analysis/Combined_analysis/combined_obj_try4/Harmony_SCT/2024_10_10_ALL_markers_for_man_annot_modified.xlsx"

# Save the workbook with the sorted sheets and formatted tables
saveWorkbook(wb, output_path, overwrite = TRUE)

cat("Processing complete. Modified file saved at:", output_path)



###making loupe for fibroblast-----------

