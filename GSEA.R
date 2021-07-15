### 2021-04-16
# Project SCS 46 patients - GSEA
# Jonne Damen
# R version 4.0.2 (2020-06-22) -- "Taking Off Again"
# R Studio version 1.3.1056 "Water Lily" (5a4dee98, 2020-07-07) for macOS

# Gene Set Enrichment Analysis (Functional Class Scoring) for every cell cluster comparison
#--------------------------------------------------------------------------
### Settings
## Set working directory
setwd("/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Enrichment analyses/")

## Load packages
library(Seurat)   #3.2.3
library(DESeq2)   #1.28.1
library(clusterProfiler) #3.16.1
#library(DOSE)   #3.14.0
library(ggplot2)    #3.3.3
library(org.Hs.eg.db)   #3.11.4

## Load data
# Seurat object
seurdata <- readRDS("/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Correct counts for sc analysis/20210420_seurdata_corrected_cellcounts.RDS")
#--------------------------------------------------------------------------
### Put metadata in global environment
meta.data <- seurdata@meta.data
#--------------------------------------------------------------------------
### Parameter of interest
parameter <- "Sex"
#--------------------------------------------------------------------------
### Set parameter of interest as active ident in seurdata
## Make new column with parameter and new.ident per cell
seurdata[["comparison.idents"]] <- paste(seurdata@meta.data[[parameter]], seurdata@meta.data$new.ident)

## Set active ident to "comparison.idents"
Idents(seurdata) <- "comparison.idents"
#--------------------------------------------------------------------------
### All the celltypes in the pseudobulk
celltypes.pb <- c("CD8A+ T Cells I", "CD4+ T Cells I", "CD16+ NK Cells", "Resident-like Macrophages", "Mixed Cells I",
                  "CD4+ T Cells II", "Inflammatory Macrophages", "CD8A+ T Cells II", "Smooth Muscle Cells",
                  "Foam Cells", "Endothelial Cells I", "CD56+ NK Cells", "Endothelial Cells II", "B Cells",
                  "Mixed Cells II",  "CD1c+ Dendritic Cells", "Mast Cells")
#--------------------------------------------------------------------------
### Prepare GSEA input - use logfoldchanges of pseudobulk

## Load pseudobulk DEG results
pb.DEG <- readRDS(file = "/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Construct pseudobulk/20210712 DEG pseudobulk.RDS")

## List of parameters: parameter variables and celltype combined for pseudobulk celltypes
list.of.parameters <- lapply(celltypes.pb, function(cluster){
  paste(unique(seurdata@meta.data[[parameter]]), cluster, sep = " ")
})

names(list.of.parameters) <- celltypes.pb

## Possible combinations of comparisons per cell cluster
# How many groups per cell cluster can be compared?
choices <- choose(length(unique(seurdata@meta.data[[parameter]])), 2)
# How many unique combinations of groups can be made?
combinations <- combn(length(unique(seurdata@meta.data[[parameter]])), 2)

## Get comparisons in pseudobulk
testing.names <- sapply(celltypes.pb, function(cluster){
  
  # What comparisons are possible
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

## Extract pseudo-bulk results per celltype
pb.comparisons <- lapply(testing.names, function(comparison){
  data.frame(pb.DEG[[comparison]])
})

## Get foldchanges from pseudo-bulk results per celltype
foldchanges <- lapply(celltypes.pb, function(celltype){
  pb.comparisons[[celltype]]$log2FoldChange
})
names(foldchanges) <- celltypes.pb

## Find entrez ID names
entrez.names <- lapply(celltypes.pb, function(celltype){
  mapIds(org.Hs.eg.db, rownames(pb.comparisons[[celltype]]), "ENTREZID", "SYMBOL")  
})
names(entrez.names) <- celltypes.pb

## Add corresponding entrez ID names and clean up logfolchanges
for(celltype in celltypes.pb){
  names(foldchanges[[celltype]]) <- entrez.names[[celltype]] # Name logfolchange vectors with entrez ID names
  foldchanges[[celltype]] <- na.omit(foldchanges[[celltype]]) # Omit any NA values
  foldchanges[[celltype]] <- foldchanges[[celltype]][foldchanges[[celltype]]!=0] # Remove the 0 values
  foldchanges[[celltype]] <- sort(foldchanges[[celltype]], decreasing = T) # Sort the list in decreasing order
  foldchanges[[celltype]] <- foldchanges[[celltype]][!is.na(names(foldchanges[[celltype]]))] # Remove NA gene names
}
#--------------------------------------------------------------------------
### Run GSEA

## Set seed
set.seed(123)

## GSEA using gene sets from KEGG pathways
gseaKEGG <- lapply(celltypes.pb, function(celltype){
  gseKEGG(geneList = foldchanges[[celltype]],  # Ordered list of genes
          organism = "hsa",  # Organism = Homo sapiens
          eps = 0.0,  # Boundary for calculating the p-value
          pvalueCutoff = 1,  # Cutoff of p-value
          verbose = F)  # Print message or not
})
names(gseaKEGG) <- celltypes.pb

#--------------------------------------------------------------------------
### Dotplots of enrichment results of interest 

## SMC enrichment results
# Extract SMC enrichment results
GSEAresult <- gseaKEGG[["Smooth Muscle Cells"]]
# Return metric of 1 or -1 when the normalized enrichment score is positive or negative, respectively
GSEAresult@result$metric <- sign(GSEAresult@result$NES)
# When metric is positive the pathway is enriched towards males
GSEAresult@result$metric[GSEAresult@result$metric == "1"] <- "male"
# When metric is negative the pathway is enriched towards females
GSEAresult@result$metric[GSEAresult@result$metric == "-1"] <- "female"
# Make dotplot of results and split by sex
dotplot(GSEAresult, split = "metric", showCategory = 10, title = "SMC enriched pathways") + facet_grid(~metric)


## ECI enrichment results
# Extract ECI enrichment results
GSEAresult <- gseaKEGG[["Endothelial Cells I"]]
# Return metric of 1 or -1 when the normalized enrichment score is positive or negative, respectively
GSEAresult@result$metric <- sign(GSEAresult@result$NES)
# When metric is positive the pathway is enriched towards males
GSEAresult@result$metric[GSEAresult@result$metric == "1"] <- "male"
# When metric is negative the pathway is enriched towards females
GSEAresult@result$metric[GSEAresult@result$metric == "-1"] <- "female"
# Make dotplot of results and split by sex
dotplot(GSEAresult, split = "metric", showCategory = 10, title = "ECI enriched pathways") + facet_grid(~metric)


## ECII enrichment results
# Extract ECII enrichment results
GSEAresult <- gseaKEGG[["Endothelial Cells II"]]
# Return metric of 1 or -1 when the normalized enrichment score is positive or negative, respectively
GSEAresult@result$metric <- sign(GSEAresult@result$NES)
# When metric is positive the pathway is enriched towards males
GSEAresult@result$metric[GSEAresult@result$metric == "1"] <- "male"
# When metric is negative the pathway is enriched towards females
GSEAresult@result$metric[GSEAresult@result$metric == "-1"] <- "female"
# Make dotplot of results and split by sex
dotplot(GSEAresult, split = "metric", showCategory = 10, title = "ECII enriched pathways") + facet_grid(~metric)

