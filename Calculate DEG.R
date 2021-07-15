### 2021-04-12
# Project SCS 46 patients - Calculate DEG
# Jonne Damen
# R version 4.0.2 (2020-06-22) -- "Taking Off Again"
# R Studio version 1.3.1056 "Water Lily" (5a4dee98, 2020-07-07) for macOS

# Calculating differentially expressed genes for every cell cluster between selected parameters
#---------------------------------------------------------------------------
### Settings
## Set working directory
setwd("/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Calculate DEG")

## Load packages
library(Seurat)   #3.2.3

## Load data
# Seurat object
seurdata <- readRDS("/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Correct counts for sc analysis/20210420_seurdata_corrected_cellcounts.RDS")
#--------------------------------------------------------------------------
### Put metadata in global environment
meta.data <- seurdata@meta.data
#--------------------------------------------------------------------------
## Parameter of interest
parameter <- "Sex"
#--------------------------------------------------------------------------
### Set parameter of interest as active ident in seurdata
## Make new column with parameter and new.ident per cell
seurdata[["comparison.idents"]] <- paste(seurdata@meta.data[[parameter]], seurdata@meta.data$new.ident)

## Set active ident to "comparison.idents"
Idents(seurdata) <- "comparison.idents"
#--------------------------------------------------------------------------
### All the celltypes in the Seurat object
celltypes <- unique(seurdata@meta.data$new.ident)
#--------------------------------------------------------------------------
### Calculate DEG between clinical parameter variables per cell cluster

## List of parameters: parameter variables and celltype combined
list.of.parameters <- lapply(celltypes, function(cluster){
  paste(unique(seurdata@meta.data[[parameter]]), cluster)
})

names(list.of.parameters) <- celltypes

## Possible combinations of comparisons per cell cluster
# How many sets per cell cluster can be compared?
choices <- choose(length(unique(seurdata@meta.data[[parameter]])), 2)
# How many unique combinations of groups can be made?
combinations <- combn(length(unique(seurdata@meta.data[[parameter]])), 2)


## Settings for DE analysis function
# Minimum percent of cells in which a gene is detected
minimum.gene.expression <- 10

# Minimum X-fold difference (log-scale) between groups of cells
logfoldchange <- 0.25

# Filter DE on significant adjusted p-value
padj <- 0.05

# Slot to pull data from
slotchoice <- "counts"

# Which statistical test to use
statistical.test <- "poisson"


## Run sc DE analysis for every cluster comparison
DEG <- sapply(celltypes, function(cluster){
  
  # What comparisons per cell cluster are possible
  comparisons <- lapply(1:choices, function(pick){
    choice1 <- combinations[, pick][1]
    choice2 <- combinations[, pick][2]
    return(c(choice1, choice2))
    
  })
  
  # Find cluster names for the comparison per celltype
  lapply(1:length(comparisons), function(pair){
    cluster1 <- list.of.parameters[[cluster]][comparisons[[pair]][1]]
    cluster2 <- list.of.parameters[[cluster]][comparisons[[pair]][2]]
    print(c(cluster1, cluster2))
    
    # Only run function if there are 3 or more cells per cluster
    if(length(which(seurdata@meta.data[["comparison.idents"]] == cluster1)) > 2 &
       length(which(seurdata@meta.data[["comparison.idents"]] == cluster2)) > 2) {
      
      
      # Function for calculating DEG between selected cell clusters
      pairwise.comparison <- FindMarkers(object = seurdata,
                                         ident.1 = cluster1 ,
                                         ident.2 = cluster2,
                                         slot = slotchoice,
                                         logfc.threshold = logfoldchange, #default 0.25
                                         min.pct = minimum.gene.expression/100,
                                         test.use = statistical.test) # default wilcox
      
      # Select only values with adjusted p-value of certain value
      pairwise.comparison <- pairwise.comparison[(pairwise.comparison$p_val_adj <= padj), ]
      
      ## Add column of gene names and column in which cluster the DEG is more highly expressed
      # Positive logfoldchange corresponds to cluster1
      cluster1.subset <- subset(pairwise.comparison, avg_logFC > 0) 
      
      # If there are avg_logFC > 0 the genes correspond to cluster1
      if (dim(cluster1.subset)[1] != 0) {
        cluster1.subset$cluster <- cluster1
        cluster1.subset$gene <- rownames(cluster1.subset)
      }
      
      
      # Negative logfoldchange corresponds to cluster2
      cluster2.subset <- subset(pairwise.comparison, avg_logFC < 0) 
      
      # If there are avg_logFC < 0 the genes correspond to cluster2
      if (dim(cluster2.subset)[1] != 0) {
        cluster2.subset$cluster <- cluster2
        cluster2.subset$gene <- rownames(cluster2.subset)
      }
      
      
      # Return result
      result <- rbind(cluster1.subset, cluster2.subset) 
      return(result)
      
    }
  })
})


## Set names for DEG: which cluster comparisons are made?
testing.names <- sapply(celltypes, function(cluster){
  
  # What comparisons per cell cluster are possible
  comparisons <- lapply(1:choices, function(pick){
    choice1 <- combinations[, pick][1]
    choice2 <- combinations[, pick][2]
    return(c(choice1, choice2))
    
  })
  
  # List of every possible comparison
  lapply(1:length(comparisons), function(pair){
    cluster1 <- list.of.parameters[[cluster]][comparisons[[pair]][1]]
    cluster2 <- list.of.parameters[[cluster]][comparisons[[pair]][2]]
    paste(cluster1, cluster2, sep = " vs. ")
  })
})

## Set names for DEG
names(DEG) <- testing.names

## If applicable: remove all comparisons that are NULL values in DEG list
DEG <- DEG[lengths(DEG) != 0]

## Save DEG file
#saveRDS(DEG, file = "20210419.sc.DEG.RDS")

