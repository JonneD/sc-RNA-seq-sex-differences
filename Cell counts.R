### 2021-05-06
# Project SCS 46 patients - Cell counts of dataset
# Jonne Damen
# R version 4.0.2 (2020-06-22) -- "Taking Off Again"
# R Studio version 1.3.1056 "Water Lily" (5a4dee98, 2020-07-07) for macOS

# Comparing cell counts between sexes per cell cluster
#--------------------------------------------------------------------------
### Settings
## Set working directory
setwd("/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Comparing counts/")

## Load packages
library(Seurat)   #3.2.2
library(RColorBrewer)   #1.1-2
library(ggplot2)    #3.3.3

## Load data
# Seurat object 
seurdata <- readRDS("/Volumes/DATASHUR_P2/Project SCS 46 patients/Data/20210409_scRNAseq_46_patients.RDS")
#--------------------------------------------------------------------------
### Give clusters the correct names

## What do we want to name the clusters? (in order from cluster 0 to 22)
celltypes <- c("CD8A+ T Cells I", "CD4+ T Cells I", "CD16+ NK Cells", "Resident-like Macrophages", "Mixed Cells I",
               "CD4+ T Cells II", "Inflammatory Macrophages", "CD8A+ T Cells II", "Unknown Cells I", "Smooth Muscle Cells",
               "Foam Cells", "Endothelial Cells I", "CD56+ NK Cells", "Endothelial Cells II", "B Cells",
               "Mixed Cells II", "Unknown Cells II", "CD1c+ Dendritic Cells", "Mast Cells", "Plasmacytoid Dendritic Cells", 
               "Plasma B Cells", "Dividing T Cells", "Neutrophils")

## Make list of clusters for iterating over
# Make list of clusters
list.of.clusters <- as.list(celltypes)
# Give corresponding cluster numbers to list
names(list.of.clusters) <- c(0:22)

## Add new.ident column to metadata
seurdata@meta.data$new.ident <- as.character(seurdata@meta.data$seurat_clusters)

## Add new cluster names to new.ident
# For iterating
numbers <- c(1:23)
# Add cluster names to corresponding cluster
for(number in numbers){
  seurdata@meta.data$new.ident[seurdata@meta.data$new.ident == names(list.of.clusters)[number]] <- list.of.clusters[number]
}

## Set factor levels of new.ident column
seurdata@meta.data$new.ident <- factor(seurdata@meta.data$new.ident, levels = c("Smooth Muscle Cells", "Endothelial Cells I", "Endothelial Cells II",
                                       "CD8A+ T Cells I", "CD8A+ T Cells II", "CD4+ T Cells I", "CD4+ T Cells II",
                                       "Dividing T Cells", "Resident-like Macrophages", "Inflammatory Macrophages", 
                                       "Foam Cells", "CD1c+ Dendritic Cells", "Plasmacytoid Dendritic Cells",
                                       "Neutrophils", "CD16+ NK Cells", "CD56+ NK Cells", "B Cells", "Plasma B Cells",
                                       "Mast Cells", "Mixed Cells I", "Mixed Cells II", "Unknown Cells I", "Unknown Cells II"))

## UMAP plot of entire dataset
# Set new.idents column as active ident
Idents(seurdata) <- "new.ident"
# Make UMAP plot
UMAPPlot(seurdata)
#--------------------------------------------------------------------------
## Put metadata in global environment
meta.data <- seurdata@meta.data
#--------------------------------------------------------------------------
### Counts per cell cluster

## List of cells per cell cluster
list.of.cells <- lapply(celltypes, function(cluster){
  cl.of.interest <- rownames(meta.data[meta.data$new.ident == cluster,])
  return(cl.of.interest)
})

## Make dataframe of cell counts per cell cluster for plotting
df <- data.frame(
  cellcounts = lengths(list.of.cells), # Get cell counts per cluster
  cluster = celltypes, # Cluster name corresponding to cell counts
  whole.dataset = rep("Whole dataset", times = 23)) # For filling in x aesthetic in ggplot

## Set cluster levels
df$cluster <- factor(df$cluster, levels = c("Smooth Muscle Cells", "Endothelial Cells I", "Endothelial Cells II",
                                            "CD8A+ T Cells I", "CD8A+ T Cells II", "CD4+ T Cells I", "CD4+ T Cells II",
                                            "Dividing T Cells", "Resident-like Macrophages", "Inflammatory Macrophages", 
                                            "Foam Cells", "CD1c+ Dendritic Cells", "Plasmacytoid Dendritic Cells",
                                            "Neutrophils", "CD16+ NK Cells", "CD56+ NK Cells", "B Cells", "Plasma B Cells",
                                            "Mast Cells", "Mixed Cells I", "Mixed Cells II", "Unknown Cells I", "Unknown Cells II"))

## Pick colors for plot
# How many colors are necessary
nb.cols <- 23
# Create color palette with sufficient amount of colors
color <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)

## Plot of cell counts
ggplot(df, aes(x = whole.dataset, y = cellcounts)) +
  geom_col(aes(fill = cluster), color = "black") + # Make a bar chart of cell counts per cluster
  ggtitle("Cell counts") + # Add title
  labs(y = "Cell count") + # Add label on y-axis
  theme(axis.title.x=element_blank()) + # Remove label on x-axis
  scale_fill_manual(values = color) # Add predetermined colors to color each cluster
#--------------------------------------------------------------------------
### Counts per cell cluster per sex

## Metadata of male cells
male.md <- subset(seurdata@meta.data, seurdata@meta.data$Sex == "male")

## List of male cells per cell cluster
list.of.male.cells <- lapply(celltypes, function(cluster){
  cl.of.interest <- rownames(male.md[male.md$new.ident == cluster,])
  return(cl.of.interest)
})

## Metadata of female cells
female.md <- subset(seurdata@meta.data, seurdata@meta.data$Sex == "female")

## List of female cells per cell cluster
list.of.female.cells <- lapply(celltypes, function(cluster){
  cl.of.interest <- rownames(female.md[female.md$new.ident == cluster,])
  return(cl.of.interest)
})

## Make dataframe
df <- data.frame(
  cellcounts = lengths(c(list.of.male.cells, list.of.female.cells)), # Get cell counts per cluster
  cluster = rep(celltypes, times = 2), # Cluster name corresponding to cell counts
  sex = rep(c("male", "female"), each = 23)) # For filling in x aesthetic in ggplot

## Set cluster levels
df$cluster <- factor(df$cluster, levels = c("Smooth Muscle Cells", "Endothelial Cells I", "Endothelial Cells II",
                                            "CD8A+ T Cells I", "CD8A+ T Cells II", "CD4+ T Cells I", "CD4+ T Cells II",
                                            "Dividing T Cells", "Resident-like Macrophages", "Inflammatory Macrophages", 
                                            "Foam Cells", "CD1c+ Dendritic Cells", "Plasmacytoid Dendritic Cells",
                                            "Neutrophils", "CD16+ NK Cells", "CD56+ NK Cells", "B Cells", "Plasma B Cells",
                                            "Mast Cells", "Mixed Cells I", "Mixed Cells II", "Unknown Cells I", "Unknown Cells II"))

## Plot of cellcounts by sex
ggplot(df, aes(x = sex, y = cellcounts)) +
  geom_col(aes(fill = cluster), color = "black") + # Make a bar chart of cell counts per cluster
  ggtitle("Cell counts per sex") + # Add title
  labs(x = "Sex", y = "Cell count") + # Add label on x- and y-axis
  scale_fill_manual(values = color) # Add predetermined colors to color each cluster

# Plot of relative contribution of cluster to sex total
ggplot(df, aes(x= sex, y = cellcounts)) +
  geom_col(aes(fill = cluster), position = "fill", colour = "black") + # Make a bar chart of relative cell counts per cluster
  ggtitle("Relative contribution of cluster to total") + # Add title
  labs(x = "Sex", y = "Fraction of total") + # Add label on x- and y-axis
  scale_fill_manual(values = color) # Add predetermined colors to color each cluster

