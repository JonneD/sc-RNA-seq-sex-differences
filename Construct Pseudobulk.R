### 2021-07-12
# Project SCS 46 patients - Construct pseudobulk
# Jonne Damen
# R version 4.0.2 (2020-06-22) -- "Taking Off Again"
# R Studio version 1.3.1056 "Water Lily" (5a4dee98, 2020-07-07) for macOS

# Constructing pseudobulk and calculating DEG between selected parameter variables
#---------------------------------------------------------------------------
### Settings
## Set working directory
setwd("/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Construct pseudobulk/")

## Load packages
library(Seurat)   #3.2.2
library(DESeq2)   #1.28.1

## Load data
# Seurat object 46 patients
seurdata <- readRDS("/Volumes/DATASHUR_P2/Project SCS 46 patients/Data/20210409_scRNAseq_46_patients.RDS")
#--------------------------------------------------------------------------
### Give clusters in seurdata the correct names

## What do we want to name the clusters? (in correct order)
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
seurdata@meta.data$new.ident <- factor(seurdata@meta.data$new.ident, levels = celltypes)
#--------------------------------------------------------------------------
### Put metadata in global environment
meta.data <- seurdata@meta.data
#--------------------------------------------------------------------------
### Parameter of interest
parameter <- "Sex"
#--------------------------------------------------------------------------
### All the celltypes in the Seurat object
celltypes <- unique(seurdata@meta.data$new.ident)
#--------------------------------------------------------------------------
### Construct pseudobulk per patient per celltype

## Make dataframe 'genes' of RNA expression of every cell
genes <- as.data.frame(seurdata@assays$SCT@counts)

## Get all patient names
patients <- unique(seurdata@meta.data$Patient)

## Minimum amount of total cells of a patient
minimum.total.cells <- 49

## Minimum amount of cells per cell cluster of a patient
cell.minimum <- 4

## Create pseudobulk of single cell RNA expression data (per patient per celltype)
patient.pseudobulk <- lapply(patients, function(patient){
  print(patient)
  
  # Get cells corresponding to patient
  cells <- seurdata@meta.data[seurdata@meta.data$Patient == patient,]
  
  # Only include patients that have more than x amount of cells
  if(length(rownames(cells)) > minimum.total.cells){
    
    # Get cell names for every celltype separately (per patient)
    per.celltype <- lapply(celltypes, function(cluster){
      print(cluster)
      
      cells1 <- rownames(cells[cells$new.ident == cluster,])
      print(cells1)
      
      # Get sums of RNA counts per celltype for every patient (only if x or more cells per celltype per patient)
      if (length(cells1) > cell.minimum) {
        counts <- data.frame(rowSums(genes[,cells1]))
        
        return(counts)
      }
    })
    
    # Name list elements
    names(per.celltype) <- celltypes
    return(per.celltype)
  }
})

## Make matrix of pseudobulk counts
# How many patients are in patient.pseudobulk list
numbers <- c(1:length(patients))

# Remove all cell clusters that are NULL values in patient.pseudobulk list
patient.pseudobulk <- lapply(numbers, function(n){
  patient.pseudobulk[[n]][lapply(patient.pseudobulk[[n]], length)!= 0]
})

# Name every patient.pseudobulk list element with the corresponding patient name
names(patient.pseudobulk) <- patients

# Remove patients with 0 pseudobulk cell clusters from patient.pseudobulk list
patient.pseudobulk <- patient.pseudobulk[lengths(patient.pseudobulk) != 0]

# Rounded matrix to put into DESeq
patient.pseudobulk <- round(as.matrix(as.data.frame(patient.pseudobulk)))
#--------------------------------------------------------------------------
### Make coldata corresponding to pseudobulk matrix

## Get annotation data for every column in patient.pseudobulk (use same function as used for constructing the pseudobulk)
coldata <- lapply(patients, function(patient){
  print(patient)
  
  # Get cells corresponding to patient
  cells <- seurdata@meta.data[seurdata@meta.data$Patient == patient,]
  
  # Only include patients that have more than x amount of cells
  if(length(rownames(cells)) > minimum.total.cells){
    
    # Get cell names for every celltype separately (per patient)
    lapply(celltypes, function(cluster){
      print(cluster)
      
      cells1 <- rownames(cells[cells$new.ident == cluster,])
      print(cells1)
      
      # Get patient number, parameter of interest, cell name, and batch corresponding to the patient.pseudobulk columns
      if (length(cells1) > cell.minimum) {
        patient.number <- meta.data[cells1, which(colnames(seurdata@meta.data) == "Patient")][1]
        parameter.oi <- meta.data[cells1, which(colnames(seurdata@meta.data) == parameter)][1]
        cell.name <- as.character(meta.data[cells1, which(colnames(seurdata@meta.data) == "new.ident")][1])
        
        return(c(patient.number, parameter.oi, cell.name))
      }
    })
  }
})

## Make dataframe of column data
# Remove all cell clusters that are NULL values in coldata list
coldata <- lapply(numbers, function(n){
  coldata[[n]][lapply(coldata[[n]], length)!= 0]
})

# Remove patients with 0 pseudobulk cell clusters from coldata list
coldata <- coldata[lengths(coldata) != 0]

# Get coldata in correct dataframe structure
coldata <- as.data.frame(t(as.data.frame(coldata)))

# Rename columns of coldata
colnames(coldata) <- c("patient.number", "parameter.oi", "celltype")

## Paste celltype and parameter.oi in new column in coldata to put as design in DeSeqDataSetFromMatrix function
coldata$celltype.and.parameter <- paste(coldata$celltype, coldata$parameter.oi)
#--------------------------------------------------------------------------
### Name columns of pseudobulk matrix object with patient number and celltype
## Paste patient number and celltype
patient.and.celltype <- paste(coldata$patient.number, coldata$celltype, sep = " ")

## Name columns of pseudobulk matrix
colnames(patient.pseudobulk) <- patient.and.celltype
#--------------------------------------------------------------------------
### Make DESeqDataSet from pseudobulk matrix and coldata
dds <- DESeqDataSetFromMatrix(countData = patient.pseudobulk, colData = coldata, design = ~ celltype.and.parameter)
#-------------------------------------------------------------------------
### Run DESeq function: DEG analysis of pseudobulk
dds <- DESeq(dds)

## Save file
#saveRDS(dds, file = "20210712 pseudobulk.RDS")
#--------------------------------------------------------------------------
## If applicable: Load the DeSeqDataSet (dds file)
#dds <- readRDS(file = "/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Construct pseudobulk/20210712 pseudobulk.RDS")
#--------------------------------------------------------------------------
### Pseudobulk DE analysis

## Check all cell clusters present in the dds object
resultsNames(dds)

## Which cell clusters have enough cell counts for pseudobulk and can be compared to each other (check in resultsNames(dds))
celltypes.pb <- c("CD8A+ T Cells I", "CD4+ T Cells I", "CD16+ NK Cells", "Resident-like Macrophages", "Mixed Cells I",
                  "CD4+ T Cells II", "Inflammatory Macrophages", "CD8A+ T Cells II", "Smooth Muscle Cells",
                  "Foam Cells", "Endothelial Cells I", "CD56+ NK Cells", "Endothelial Cells II", "B Cells",
                  "Mixed Cells II",  "CD1c+ Dendritic Cells", "Mast Cells")

## List of parameters: parameter variables and celltype combined for pseudobulk celltypes
list.of.parameters <- lapply(celltypes.pb, function(cluster){
  paste(cluster, unique(seurdata@meta.data[[parameter]]), sep = " ")
})

names(list.of.parameters) <- celltypes.pb


## Possible combinations of comparisons per cell cluster
# How many groups per cell cluster can be compared?
choices <- choose(length(unique(seurdata@meta.data[[parameter]])), 2)
# How many unique combinations of groups can be made?
combinations <- combn(length(unique(seurdata@meta.data[[parameter]])), 2)


## Calculate pseudobulk DE results between clusters
DE.results <- sapply(celltypes.pb, function(cluster){
  
  # What comparisons per cell cluster are possible
  comparisons <- lapply(1:choices, function(pick){
    choice1 <- combinations[, pick][1]
    choice2 <- combinations[, pick][2]
    return(c(choice1, choice2))
    
  })
  
  # Print every possible comparison
  lapply(1:length(comparisons), function(pair){
    cluster1 <- list.of.parameters[[cluster]][comparisons[[pair]][1]]
    cluster2 <- list.of.parameters[[cluster]][comparisons[[pair]][2]]
    print(c(cluster1, cluster2))
    
    # Run DE results formula for comparison
    res <- results(dds, contrast = c("celltype.and.parameter", cluster1, cluster2))
    return(res)
  })
})

## List of testing clusters ordered differently: parameter first, then cluster
list.of.parameters <- lapply(celltypes.pb, function(cluster){
  paste(unique(seurdata@meta.data[[parameter]]), cluster, sep = " ")
})

names(list.of.parameters) <- celltypes.pb

## Set names for DE.results: which comparisons are made?
testing.names <- sapply(celltypes.pb, function(cluster){
  
  # What comparisons are possible
  comparisons <- lapply(1:choices, function(pick){
    choice1 <- combinations[, pick][1]
    choice2 <- combinations[, pick][2]
    return(c(choice1, choice2))
    
  })
  
  # Print every possible comparison
  lapply(1:length(comparisons), function(pair){
    cluster1 <- list.of.parameters[[cluster]][comparisons[[pair]][1]]
    cluster2 <- list.of.parameters[[cluster]][comparisons[[pair]][2]]
    paste(cluster1, cluster2, sep = " vs. ")
  })
})

## Set names for DEG
names(DE.results) <- testing.names

## Save DE.results object
#saveRDS(DE.results, file = "20210712 DEG pseudobulk.RDS")





