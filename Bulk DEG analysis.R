### 2021-05-17
# Project SCS 46 patients - Calculate DEG in bulkdata
# Jonne Damen
# R version 4.0.2 (2020-06-22) -- "Taking Off Again"
# R Studio version 1.3.1056 "Water Lily" (5a4dee98, 2020-07-07) for macOS

# Calculating differentially expressed genes between sexes in bulk data and projecting results onto single cell data
#--------------------------------------------------------------------------
### Settings
## Set working directory
setwd("/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Bulk to single cell/")

## Load package
library(Seurat) # 3.2.3
library(DESeq2) # 1.28.1

## Load data
# Bulk RNA seq data
bulkdata <- readRDS("/Volumes/DATASHUR_P2/Project SCS 46 patients/Data/bulk_data_for_Jonne.RDS")
# Update Seurat object to be compatible with current Seurat version
bulkdata <- UpdateSeuratObject(bulkdata)
#--------------------------------------------------------------------------
### Put metadata in global environment
meta.data <- bulkdata@meta.data
#--------------------------------------------------------------------------
### Calculate DEG between sexes

## Get bulk gene expression
countdata <- as.matrix(bulkdata@assays[["RNA"]]@counts)

## Determine descriptive data of bulk gene expression columns
coldata <- data.frame(sex = bulkdata@meta.data$sex, 
                      cluster = bulkdata@meta.data$seurat_clusters)
# Paste cluster and sex together
coldata$cluster.and.sex <- paste(coldata$cluster, coldata$sex, sep = "")

## Construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ cluster.and.sex)

## Run DESeq
dds <- DESeq(dds)

## Save dds object
#saveRDS(dds, file = "Bulk-cluster.and.sex.RDS")
## Load dds object
#dds <- readRDS("/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Bulk to single cell/Bulk-cluster.and.sex.RDS")
#--------------------------------------------------------------------------
### Get DEG results for every bulk cluster

## Determine all clusters in the object
clusters <- unique(bulkdata@meta.data$seurat_clusters)

## Determine all cluster and sex combinations
list.of.options <- lapply(clusters, function(cluster){
  paste(cluster, unique(bulkdata@meta.data$sex), sep = "")
})
names(list.of.options) <- clusters


## Get bulk DEG results for all options
list.of.results <- lapply(clusters, function(cluster){
  cluster1 <- list.of.options[[cluster]][1] # Determine cluster 1 for contrast
  cluster2 <- list.of.options[[cluster]][2] # Determine cluster 2 for contrast
  res <- results(dds, contrast = c("cluster.and.sex", cluster1, cluster2)) # Extract results between clusters
  res <- na.omit(res) # Omit NA results
  res <- as.data.frame(subset(res, padj < 0.1)) # Subset by adjusted p-value < 0.1
  res$gene <- rownames(res) # Add column with gene name
  return(res)
})
names(list.of.results) <- clusters
#--------------------------------------------------------------------------
### Look up DEG in scRNAseq data
## Load scRNAseq data
seurdata <- readRDS("/Volumes/DATASHUR_P2/Project SCS 46 patients/Data/20210409_scRNAseq_46_patients.RDS")

### Give clusters of scRNAseq data the correct names
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

## Set new.idents column as active ident
Idents(seurdata) <- seurdata@meta.data$new.ident
#--------------------------------------------------------------------------
### Look up genes from bulk data in single cell data

## Get all DEGs of bulk clusters
genes <- unlist(sapply(clusters, function(cluster){
  list.of.results[[cluster]][["gene"]]
}))

## Remove gene that is not in scRNAseq dataset
genes <- genes[! genes %in% "PSMC4.1"]

## Set directory for results
#setwd("/Volumes/JONNE/Documenten/BoD/Stage/Thesis/5. supplementary/")

## Create result folder
#result.folder <- paste(Sys.Date(), "result", sep = " ")
#dir.create(result.folder, showWarnings = FALSE)

## Make violin plots per DEG - projected onto single cell data
#pdf(paste(result.folder, "/", Sys.Date(), "_plots_bulk_to_sc.pdf", sep = ""), width = 15, height = 10)
lapply(genes, function(gene){
  VlnPlot(seurdata, features = gene, split.by = "Sex")
})
#dev.off()

## Make violin plots per DEG - bulk expression
#pdf(paste(result.folder, "/", Sys.Date(), "_plots_bulk_in_bulk.pdf", sep = ""), width = 15, height = 10)
lapply(genes, function(gene){
  VlnPlot(bulkdata, features = gene, split.by = "sex")
})
#dev.off()


