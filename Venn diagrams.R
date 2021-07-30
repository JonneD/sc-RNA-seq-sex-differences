### 2021-04-16
# Project 46 patients - Venn diagrams
# Jonne Damen
# R version 4.0.2 (2020-06-22) -- "Taking Off Again"
# R Studio version 1.3.1056 "Water Lily" (5a4dee98, 2020-07-07) for macOS

# Venn diagrams of DEGs identified in single cell and pseudobulk analysis
#---------------------------------------------------------------------------
### Settings
## Set working directory
setwd("/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Visualize DEG/")

## load packages
library(Seurat)   #3.2.3
library(DESeq2)   #1.28.1
library(ggplot2)  #3.3.3
library(ggvenn)   #0.1.9

## Load data
# Seurat object 
seurdata <- readRDS("/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Correct counts for sc analysis/20210420_seurdata_corrected_cellcounts.RDS")

# sc DEG analysis results
sc.DEG <- readRDS(file = "/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Calculate DEG/20210419.sc.DEG.RDS")

# Pseudobulk DEG analysis results
pb.DEG <- readRDS(file = "/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Construct pseudobulk/20210712 DEG pseudobulk.RDS")
#--------------------------------------------------------------------------
### Put metadata in global environment
meta.data <- seurdata@meta.data
#--------------------------------------------------------------------------
### Parameter of interest
parameter <- "Sex"
#--------------------------------------------------------------------------
## Set parameter of interest as active ident in seurdata
# Make new column with parameter and new.ident per cell
seurdata[["comparison.idents"]] <- paste(seurdata@meta.data[[parameter]], seurdata@meta.data$new.ident)

# Set active ident to "comparison.idents"
Idents(seurdata) <- "comparison.idents"
#--------------------------------------------------------------------------
### Compare sc results with pseudobulk

## Clusters that are tested in the pseudobulk 
celltypes.pb <- c("CD8A+ T Cells I", "CD4+ T Cells I", "CD16+ NK Cells", "Resident-like Macrophages", "Mixed Cells I",
                  "CD4+ T Cells II", "Inflammatory Macrophages", "CD8A+ T Cells II", "Smooth Muscle Cells",
                  "Foam Cells", "Endothelial Cells I", "CD56+ NK Cells", "Endothelial Cells II", "B Cells",
                  "Mixed Cells II",  "CD1c+ Dendritic Cells", "Mast Cells")

## List of parameters: parameter variables and celltype combined
list.of.parameters <- lapply(celltypes.pb, function(cluster){
  paste(unique(seurdata@meta.data[[parameter]]), cluster, sep = " ")
})

names(list.of.parameters) <- celltypes.pb

## Possible combinations of comparisons per cell cluster
# How many groups per cell cluster can be compared?
choices <- choose(length(unique(seurdata@meta.data[[parameter]])), 2)
# How many unique combinations of groups can be made?
combinations <- combn(length(unique(seurdata@meta.data[[parameter]])), 2)


## Testing names: what groups are compared to each other?
testing.names <- sapply(celltypes.pb, function(cluster){
  
  # What comparisons are possible
  comparisons <- lapply(1:choices, function(pick){
    choice1 <- combinations[, pick][1]
    choice2 <- combinations[, pick][2]
    return(c(choice1, choice2))
    
  })
  
  # List of every possible combination
  lapply(1:length(comparisons), function(pair){
    cluster1 <- list.of.parameters[[cluster]][comparisons[[pair]][1]]
    cluster2 <- list.of.parameters[[cluster]][comparisons[[pair]][2]]
    paste(cluster1, cluster2, sep = " vs. ")
  })
})


## Separate strings of groups that are compared to each other
testing.clusters <- sapply(celltypes.pb, function(cluster){
  strsplit(testing.names[[cluster]], split = " vs. ")
})

names(testing.clusters) <- testing.names


## Get list of sc and pb DEG names
sc.pb.genes <- lapply(testing.names, function(cluster){
  
  # Select sc DEG names
  sc <- rownames(sc.DEG[[cluster]])
  
  # Get significant genes in pb
  pb.sign <- as.data.frame(pb.DEG[[cluster]][pb.DEG[[cluster]]$pvalue < 0.05,])
  
  # Select pb DEG names
  pb <- rownames(pb.sign)
  
  # Put into list
  res <- list("sc padj < 0.05" = sc,
              "pb p-value < 0.05" = pb)
  return(res)
  
})
names(sc.pb.genes) <- testing.names


## Venn diagram of SMC comparison
ggvenn(sc.pb.genes[["male Smooth Muscle Cells vs. female Smooth Muscle Cells"]], show_percentage = F, text_size = 7) + 
  ggtitle("SMC sc and pb gene overlap") +
  theme(plot.title = element_text(size = 20)) +
  theme(plot.title = element_text(hjust = 0.5))

## Venn diagram of ECI comparison
ggvenn(sc.pb.genes[["male Endothelial Cells I vs. female Endothelial Cells I"]], show_percentage = F, text_size = 7) + 
  ggtitle("ECI sc and pb gene overlap") +
  theme(plot.title = element_text(size = 20)) +
  theme(plot.title = element_text(hjust = 0.5))

## Venn diagram of ECII comparison
ggvenn(sc.pb.genes[["male Endothelial Cells II vs. female Endothelial Cells II"]], show_percentage = F, text_size = 7) + 
  ggtitle("ECII sc and pb gene overlap") +
  theme(plot.title = element_text(size = 20)) +
  theme(plot.title = element_text(hjust = 0.5))



