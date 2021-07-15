### 2021-05-22
# Project SCS 46 patients - Macrophage and T cell versus structural cell ratio
# Jonne Damen
# R version 4.0.2 (2020-06-22) -- "Taking Off Again"
# R Studio version 1.3.1056 "Water Lily" (5a4dee98, 2020-07-07) for macOS

# Calculate the macrophage and T cell versus structural cell ratios between parameters of interest
#--------------------------------------------------------------------------
### Settings
## Set working directory
setwd("/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Calculate ratios/")

## Load packages
library(Seurat)   #3.2.2
library(RColorBrewer)   #1.1-2
library(ggplot2)    #3.3.3
library(ggpubr)   #0.4.0

## Load data
# Seurat object
seurdata <- readRDS("/Volumes/DATASHUR_P2/Project SCS 46 patients/Data/20210409_scRNAseq_46_patients.RDS")
#--------------------------------------------------------------------------
### Make new column in metadata

## Add new.ident column from seurat_clusters
seurdata@meta.data$new.ident <- as.character(seurdata@meta.data$seurat_clusters)
#--------------------------------------------------------------------------
### Take together cell types
## Cluster 3, 6, and 10 are macrophages
macrophages <- c("3", "6", "10")
# Rename clusters
for(cluster in macrophages){
  seurdata@meta.data$new.ident[seurdata@meta.data$new.ident == cluster] <- "Macrophage"
}

## Cluster 0, 1, 5, and 7 are T cells
tcells <- c("0", "1", "5", "7")
# Rename clusters
for(cluster in tcells){
  seurdata@meta.data$new.ident[seurdata@meta.data$new.ident == cluster] <- "T cell"
}

## Cluster 9, 11, 13 are structural cells
structural <- c("9", "11", "13")
# Rename clusters
for(cluster in structural){
  seurdata@meta.data$new.ident[seurdata@meta.data$new.ident == cluster] <- "structural"
}
#--------------------------------------------------------------------------
## Put metadata in global environment
meta.data <- seurdata@meta.data
#--------------------------------------------------------------------------
### Parameter of interest
parameter <- "Sex"
#--------------------------------------------------------------------------
### Calculate ratios
## Of which celltypes do we want to calculate ratios to structural cells?
cell.options <- c("Macrophage", "T cell")

## Function for calculating ratios
all.ratios <- lapply(cell.options, function(cl){
  
  ## Which celltypes to calculate ratio for
  celltypes <- c(cl, "structural") 
  
  ## Prepare dataframe
  # Patient numbers
  patients <- unique(meta.data$Patient)
  
  # Patient total counts
  total.counts <- lapply(patients, function(p){
    length(rownames(meta.data[meta.data$Patient == p,]))
  })
  names(total.counts) <- patients
  
  # Minimum amount of total cells of a patient
  minimum.total.cells <- 49
  
  # Patients with x or more cells
  patients <- unlist(lapply(patients, function(p){
    if(total.counts[[p]] > minimum.total.cells){
      return(p)
    }
  }))
  
  # List of metadata of every patient
  list.of.patients.md <- lapply(patients, function(patient){
    patient <- subset(meta.data, meta.data$Patient == patient)
    return(patient)
  })
  names(list.of.patients.md) <- patients
  
  # Get rownames of cells of every cluster for every patient
  cells.per.patient <- lapply(patients, function(patient){
    lapply(celltypes, function(cluster){
      cl.of.interest <- rownames(list.of.patients.md[[patient]][list.of.patients.md[[patient]]$new.ident == cluster,])
      return(cl.of.interest)
    })
  })
  names(cells.per.patient) <- patients
  
  # Count rownames of all cells to get counts per patient
  counts.per.patient <- lapply(patients, function(patient){
    counts <- lengths(cells.per.patient[[patient]])
    return(counts)
  })
  names(counts.per.patient) <- patients
  
  
  # Get the relative counts per patient
  # How many cell clusters are there
  numbers <- 1:length(celltypes)
  
  # Calculate relative counts
  relative.counts.per.patient <- lapply(patients, function(patient){
    sapply(numbers, function(num){
      counts.per.patient[[patient]][num] / sum(counts.per.patient[[patient]]) * 100
    })
  })
  names(relative.counts.per.patient) <- patients
  
  # Parameter of interest
  parameter.oi <- sapply(patients, function(patient){
    param <- list.of.patients.md[[patient]][1, which(colnames(seurdata@meta.data) == parameter)]
    return(param)
  })
  
  
  # Dataframe for plotting
  df <- data.frame(
    counts <- unlist(relative.counts.per.patient),
    cluster <- rep(celltypes, times = length(patients)),
    param <- rep(parameter.oi, each = length(celltypes))
  )
  
  # Names of dataframe
  names(df) <- c("cellcounts", "cluster", "sex")
  
  # Make factors of clusters and set order of factors
  df$cluster <- factor(df$cluster, levels = celltypes)
  
  ## Make boxplot
  relative.count.boxplot <- ggplot(df, aes(x = cluster, y = cellcounts, fill = sex)) +
    geom_boxplot(outlier.size = 0.5) +   # Make boxplot and adjust outlier size
    geom_point(position = position_dodge(width = 0.75), aes(group = sex), size = 0.5) +    # Every datapoint as a dot
    stat_compare_means(method = "wilcox.test") +    # Statistical test
    labs(title = "Relative cell counts of patients per cell type separated by sex", y = "Relative cell counts")    # Title and y-label
  
  ## Return boxplot
  #return(relative.count.boxplot)
  
  ## Add patient number to df
  df$patient <- rep(patients, each = length(celltypes))
  
  ## Calculate macrophage/structural ratio per patient
  MS.ratio <- lapply(patients, function(patient){
    subset(df$cellcounts, df$patient == patient)[1] / subset(df$cellcounts, df$patient == patient)[2]
  })
  
  ## Make dataframe
  MS.ratio <- as.data.frame(do.call(rbind, MS.ratio))
  names(MS.ratio) <- "ratio"
  
  ## Add patient column to dataframe
  MS.ratio$patient <- patients
  
  ## Add sex column to dataframe
  MS.ratio$sex <- parameter.oi
  
  ## Non normally distributed data so boxplot with Wilcoxon rank sum test
  ratio <- ggplot(MS.ratio, aes(x = sex, y = ratio)) + 
    geom_boxplot() + # Make boxplot
    geom_point(position = position_dodge(width = 0.75)) + # Every datapoint as a dot
    labs(title = paste(cl, "/structural cell ratio per sex", sep = "")) + # Title
    stat_compare_means(method = "wilcox.test") # Statistical test
  
  return(ratio)
})

