### 2021-04-20
# Project SCS 46 patients - Corrected counts for sc analysis
# Jonne Damen
# R version 4.0.2 (2020-06-22) -- "Taking Off Again"
# R Studio version 1.3.1056 "Water Lily" (5a4dee98, 2020-07-07) for macOS

### Correcting patient contributions to cell clusters:
# Determining patient contributions to cell clusters
# 1 patient should maximally contribute 20% to a cell cluster
# Random sample of cells for a patient to make cell contribution of 20%
#--------------------------------------------------------------------------
### Settings
## Set working directory
setwd("/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Correct counts for sc analysis/")

## Load packages
library(Seurat)   #3.2.2
library(ggplot2)  #3.3.3

## Load data
# Seurat object 
seurdata <- readRDS("/Volumes/DATASHUR_P2/Project SCS 46 patients/Data/20210409_scRNAseq_46_patients.RDS")
#--------------------------------------------------------------------------
### Give clusters the correct names

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

## Set active ident of seurdata
Idents(seurdata) <- seurdata@meta.data$new.ident
#--------------------------------------------------------------------------

### ROUND 1 of sampling excess cells

#--------------------------------------------------------------------------
## Put metadata in global environment
meta.data <- seurdata@meta.data
#--------------------------------------------------------------------------
### Order cell clusters from biggest to smallest

## Patient numbers in the Seurat object
patients <- unique(meta.data$Patient)

## List of metadata of every patient
list.of.patients.md <- lapply(patients, function(patient){
  patient <- subset(meta.data, meta.data$Patient == patient)
  return(patient)
})
names(list.of.patients.md) <- patients

## Determine how big the cell clusters are
# Get cell names for every cluster for every patient
cells.per.cluster <- lapply(celltypes, function(cluster){
  lapply(patients, function(patient){
    cl.of.interest <- rownames(list.of.patients.md[[patient]][list.of.patients.md[[patient]]$new.ident == cluster,])
    return(cl.of.interest)
  })
})
names(cells.per.cluster) <- celltypes 

# Vector of sum of all the cells in every cluster
counts.per.cluster <- unlist(lapply(celltypes, function(cluster){
  sum(lengths(cells.per.cluster[[cluster]]))
}))

# Order from biggest to smallest count
correct.order <- order(counts.per.cluster, decreasing = T)

# Order celltypes by correct.order
celltypes <- sapply(correct.order, function(place){
  celltypes[place]
})
#--------------------------------------------------------------------------
### Barplot of relative cell counts per patient for every cell cluster

## Count the amount of cells per patient for every cell cluster
counts.per.cluster <- lapply(celltypes, function(cluster){
  counts <- lengths(cells.per.cluster[[cluster]])
  return(counts)
})
names(counts.per.cluster) <- celltypes

## Give counts corresponding patient names
for(cluster in celltypes){
  names(counts.per.cluster[[cluster]]) <- patients
}

## Calculate the relative counts per patient
# How many patients are there
numbers <- 1:length(patients)

# Calculate relative counts per cell cluster per patient
relative.counts.per.cluster <- lapply(celltypes, function(cluster){
  sapply(numbers, function(num){
    counts.per.cluster[[cluster]][num] / sum(counts.per.cluster[[cluster]]) * 100
  })
})
names(relative.counts.per.cluster) <- celltypes

## Give relative counts corresponding patient names
for(cluster in celltypes){
  names(relative.counts.per.cluster[[cluster]]) <- patients
}

## Barplot of relative counts
# Unlist and make dataframe of relative counts per cell cluster per patient
relative.counts <- data.frame(relative.count = unlist(relative.counts.per.cluster),
                              patient = rep(patients, times = length(celltypes)),
                              cluster = rep(celltypes, each = length(patients)))

# Make factors of clusters and set correct order of factor
relative.counts$cluster <- factor(relative.counts$cluster, levels = celltypes)

## Plot relative counts of patients per cell cluster in barplot
ggplot(relative.counts, aes(x = cluster, y = relative.count)) +
  geom_col(aes(fill = patient), colour = "black") + # Make a bar chart of relative cell counts per cluster per patient
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Put text on x-axis in 45 degree angle
  ggtitle("Initial patient contribution to celltype") # Add title
#--------------------------------------------------------------------------
### Select which patients to sample cells from and how many cells to sample

## Which patients contribute more than 20 percent to which cell cluster (determine position in vector)
patient.position <- lapply(celltypes, function(cluster){
  which(relative.counts.per.cluster[[cluster]] > 20)
})
names(patient.position) <- celltypes

## Find out how many cells to sample
# Remove the 0 values in patient.position
patient.position <- patient.position[lengths(patient.position) > 0L]

# How many elements are in patient.position
num <- 1:sum(lengths(patient.position))

# Vector of which clusters have a high contributing patient
select.clusters <- sapply(num, function(number){
  over.20.cluster.and.patient <- strsplit(names(unlist(patient.position)), split = "[.]") # Split unlisted names of patient.position into cell clusters and patient
  over.20.cluster <- over.20.cluster.and.patient[[number]][1] # Get cluster name of high contributing patient
  return(over.20.cluster)
})         

# Make vector of patient.position (get positions of high contributing patients in counts.per.cluster vector)
patient.position <- unlist(patient.position)

# Get cell count of the clusters of high contributing patients
cell.count <- sapply(num, function(number){
  counts.per.cluster[[select.clusters[number]]][patient.position[number]]
})
names(cell.count) <- select.clusters

# Get counts.per.cluster corresponding to selected clusters
selected.counts <- sapply(select.clusters, function(cluster){
  list(counts.per.cluster[[cluster]])
})

# How many cells to remove from dataset
to.sample <- lapply(num, function(number){
  a <- sum(selected.counts[[number]]) * 0.2 # What is 20% of the total amount of cells in the cell cluster of interest
  b <- cell.count[[number]] - a # How many cells to remove from the cell cluster of the high contributing patient
  c <- sum(selected.counts[[number]]) - b # How many cells are left in the total cell cluster when the patient's cells are removed
  d <- c * 0.2 # What is 20% of the new amount of cells in the cell cluster of interest
  for.sampling <- round(cell.count[[number]] - d) # How many cells should be removed from the cell cluster of the high contributing patient
})

#--------------------------------------------------------------------------
### Sample cells to remove and construct new Seurat object

## Extract patient names of patients with high contributing clusters
select.patients <- sapply(num, function(number){
  over.20.cluster.and.patient <- strsplit(names(unlist(patient.position)), split = "[.]") # Split unlisted names of patient.position into cell clusters and patient
  over.20.patient <- over.20.cluster.and.patient[[number]][2] # Get high contributing patient name
  return(over.20.patient)
})   

## List of patient names, clusters, and amount of elements
list.of.info <- list(select.patients, select.clusters, num)

## Random sampling of cells to exclude
set.seed(123)
to.remove <- unlist(sapply(num, function(n){
    sample(rownames(subset(meta.data, meta.data$Patient == list.of.info[[1]][n] & meta.data$new.ident == list.of.info[[2]][n])), # Determine patient to sample and cluster to sample from list.of.info
      size = to.sample[[list.of.info[[3]][n]]], # Determine amount of cells to sample from this patient cluster
      replace = F) # Once a cell has been sampled, take out of the pool for sampling
}))

# All the cell names in the Seurat object
all.cells <- rownames(meta.data)

# Select cells to keep in Seurat object by removing the to.remove cells
keep.cells <- all.cells[!all.cells %in% to.remove]

## Subset Seurat object with desired cells
corrected.counts.seurdata <- subset(seurdata, cells = keep.cells)
seurdata <- corrected.counts.seurdata
#--------------------------------------------------------------------------

### ROUND 2 of sampling excess cells

#--------------------------------------------------------------------------
## Put metadata in global environment
meta.data <- seurdata@meta.data
#--------------------------------------------------------------------------
### Order cell clusters from biggest to smallest

## Patient numbers in the Seurat object
patients <- unique(meta.data$Patient)

## List of metadata of every patient
list.of.patients.md <- lapply(patients, function(patient){
  patient <- subset(meta.data, meta.data$Patient == patient)
  return(patient)
})
names(list.of.patients.md) <- patients

## Determine how big the cell clusters are
# Get cell names for every cluster for every patient
cells.per.cluster <- lapply(celltypes, function(cluster){
  lapply(patients, function(patient){
    cl.of.interest <- rownames(list.of.patients.md[[patient]][list.of.patients.md[[patient]]$new.ident == cluster,])
    return(cl.of.interest)
  })
})
names(cells.per.cluster) <- celltypes 

# Vector of sum of all the cells in every cluster
counts.per.cluster <- unlist(lapply(celltypes, function(cluster){
  sum(lengths(cells.per.cluster[[cluster]]))
}))

# Order from biggest to smallest count
correct.order <- order(counts.per.cluster, decreasing = T)

# Order celltypes by correct.order
celltypes <- sapply(correct.order, function(place){
  celltypes[place]
})
#--------------------------------------------------------------------------
### Barplot of relative cell counts per patient for every cell cluster

## Count the amount of cells per patient for every cell cluster
counts.per.cluster <- lapply(celltypes, function(cluster){
  counts <- lengths(cells.per.cluster[[cluster]])
  return(counts)
})
names(counts.per.cluster) <- celltypes

## Give counts corresponding patient names
for(cluster in celltypes){
  names(counts.per.cluster[[cluster]]) <- patients
}

## Calculate the relative counts per patient
# How many patients are there
numbers <- 1:length(patients)

# Calculate relative counts per cell cluster per patient
relative.counts.per.cluster <- lapply(celltypes, function(cluster){
  sapply(numbers, function(num){
    counts.per.cluster[[cluster]][num] / sum(counts.per.cluster[[cluster]]) * 100
  })
})
names(relative.counts.per.cluster) <- celltypes

## Give relative counts corresponding patient names
for(cluster in celltypes){
  names(relative.counts.per.cluster[[cluster]]) <- patients
}

## Barplot of relative counts
# Unlist and make dataframe of relative counts per cell cluster per patient
relative.counts <- data.frame(relative.count = unlist(relative.counts.per.cluster),
                              patient = rep(patients, times = length(celltypes)),
                              cluster = rep(celltypes, each = length(patients)))

# Make factors of clusters and set correct order of factor
relative.counts$cluster <- factor(relative.counts$cluster, levels = celltypes)

## Plot relative counts of patients per cell cluster in barplot
ggplot(relative.counts, aes(x = cluster, y = relative.count)) +
  geom_col(aes(fill = patient), colour = "black") + # Make a bar chart of relative cell counts per cluster per patient
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Put text on x-axis in 45 degree angle
  ggtitle("Patient contribution to celltype after 1 round of sampling") # Add title
#--------------------------------------------------------------------------
### Select which patients to sample cells from and how many cells to sample

## Which patients contribute more than 20 percent to which cell cluster (determine position in vector)
patient.position <- lapply(celltypes, function(cluster){
  which(relative.counts.per.cluster[[cluster]] > 20)
})
names(patient.position) <- celltypes

## Find out how many cells to sample
# Remove the 0 values in patient.position
patient.position <- patient.position[lengths(patient.position) > 0L]

# How many elements are in select
num <- 1:sum(lengths(patient.position))

# Vector of which clusters have a high contributing patient
select.clusters <- sapply(num, function(number){
  over.20.cluster.and.patient <- strsplit(names(unlist(patient.position)), split = "[.]") # Split unlisted names of patient.position into cell clusters and patient
  over.20.cluster <- over.20.cluster.and.patient[[number]][1] # Get cluster name of high contributing patient
  return(over.20.cluster)
})         

# Make vector of patient.position (get positions of high contributing patients in counts.per.cluster vector)
patient.position <- unlist(patient.position)

# Get cell count of the clusters of high contributing patients
cell.count <- sapply(num, function(number){
  counts.per.cluster[[select.clusters[number]]][patient.position[number]]
})
names(cell.count) <- select.clusters

# Get counts.per.cluster corresponding to selected clusters
selected.counts <- sapply(select.clusters, function(cluster){
  list(counts.per.cluster[[cluster]])
})

# How many cells to remove from dataset
to.sample <- lapply(num, function(number){
  a <- sum(selected.counts[[number]]) * 0.2 # What is 20% of the total amount of cells in the cell cluster of interest
  b <- cell.count[[number]] - a # How many cells to remove from the cell cluster of the high contributing patient
  c <- sum(selected.counts[[number]]) - b # How many cells are left in the total cell cluster when the patient's cells are removed
  d <- c * 0.2 # What is 20% of the new amount of cells in the cell cluster of interest
  for.sampling <- round(cell.count[[number]] - d) # How many cells should be removed from the cell cluster of the high contributing patient
})

#--------------------------------------------------------------------------
### Sample cells to remove and construct new Seurat object

## Extract patient names of patients with high contributing clusters
select.patients <- sapply(num, function(number){
  over.20.cluster.and.patient <- strsplit(names(unlist(patient.position)), split = "[.]") # Split unlisted names of patient.position into cell clusters and patient
  over.20.patient <- over.20.cluster.and.patient[[number]][2] # Get high contributing patient name
  return(over.20.patient)
})   

## List of patient names, clusters, and amount of elements
list.of.info <- list(select.patients, select.clusters, num)

## Random sampling of cells to exclude
set.seed(123)
to.remove <- unlist(sapply(num, function(n){
  sample(rownames(subset(meta.data, meta.data$Patient == list.of.info[[1]][n] & meta.data$new.ident == list.of.info[[2]][n])), # Determine patient to sample and cluster to sample from list.of.info
         size = to.sample[[list.of.info[[3]][n]]], # Determine amount of cells to sample from this patient cluster
         replace = F) # Once a cell has been sampled, take out of the pool for sampling
}))

# All the cell names in the Seurat object
all.cells <- rownames(meta.data)

# Selecting cells to keep in Seurat object by removing the to.remove cells
keep.cells <- all.cells[!all.cells %in% to.remove]

## Subset Seurat object with desired cells
corrected.counts.seurdata <- subset(seurdata, cells = keep.cells)
seurdata <- corrected.counts.seurdata
#--------------------------------------------------------------------------

### CHECK NEWLY CREATED SEURAT OBJECT CELL COUNT RATIO'S

#--------------------------------------------------------------------------
## Put metadata in global environment
meta.data <- seurdata@meta.data
#--------------------------------------------------------------------------
### Order cell clusters from biggest to smallest

## Patient numbers in the Seurat object
patients <- unique(meta.data$Patient)

## List of metadata of every patient
list.of.patients.md <- lapply(patients, function(patient){
  patient <- subset(meta.data, meta.data$Patient == patient)
  return(patient)
})
names(list.of.patients.md) <- patients

## Determine how big the cell clusters are
# Get cell names for every cluster for every patient
cells.per.cluster <- lapply(celltypes, function(cluster){
  lapply(patients, function(patient){
    cl.of.interest <- rownames(list.of.patients.md[[patient]][list.of.patients.md[[patient]]$new.ident == cluster,])
    return(cl.of.interest)
  })
})
names(cells.per.cluster) <- celltypes 

# Vector of sum of all the cells in every cluster
counts.per.cluster <- unlist(lapply(celltypes, function(cluster){
  sum(lengths(cells.per.cluster[[cluster]]))
}))

# Order from biggest to smallest count
correct.order <- order(counts.per.cluster, decreasing = T)

# Order celltypes by correct.order
celltypes <- sapply(correct.order, function(place){
  celltypes[place]
})
#--------------------------------------------------------------------------
### Barplot of relative cell counts per patient for every cell cluster

## Count the amount of cells per patient for every cell cluster
counts.per.cluster <- lapply(celltypes, function(cluster){
  counts <- lengths(cells.per.cluster[[cluster]])
  return(counts)
})
names(counts.per.cluster) <- celltypes

## Give counts corresponding patient names
for(cluster in celltypes){
  names(counts.per.cluster[[cluster]]) <- patients
}

## Calculate the relative counts per patient
# How many patients are there
numbers <- 1:length(patients)

# Calculate relative counts per cell cluster per patient
relative.counts.per.cluster <- lapply(celltypes, function(cluster){
  sapply(numbers, function(num){
    counts.per.cluster[[cluster]][num] / sum(counts.per.cluster[[cluster]]) * 100
  })
})
names(relative.counts.per.cluster) <- celltypes

## Give relative counts corresponding patient names
for(cluster in celltypes){
  names(relative.counts.per.cluster[[cluster]]) <- patients
}

## Barplot of relative counts
# Unlist and make dataframe of relative counts per cell cluster per patient
relative.counts <- data.frame(relative.count = unlist(relative.counts.per.cluster),
                              patient = rep(patients, times = length(celltypes)),
                              cluster = rep(celltypes, each = length(patients)))

# Make factors of clusters and set correct order of factor
relative.counts$cluster <- factor(relative.counts$cluster, levels = celltypes)

## Plot relative counts of patients per cell cluster in barplot
ggplot(relative.counts, aes(x = cluster, y = relative.count)) +
  geom_col(aes(fill = patient), colour = "black") + # Make a bar chart of relative cell counts per cluster per patient
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Put text on x-axis in 45 degree angle
  ggtitle("Final patient contribution to celltype") # Add title
#--------------------------------------------------------------------------
## Save new sampled Seurat object
#saveRDS(seurdata, "20210420_seurdata_corrected_cellcounts.RDS")
