### 2021-04-16
# Side project 46 patients - Visualize DEG
# Jonne Damen
# R version 4.0.2 (2020-06-22) -- "Taking Off Again"
# R Studio version 1.3.1056 "Water Lily" (5a4dee98, 2020-07-07) for macOS

# Dotplots of average expression of DEGs significant in both sc and pseudobulk analysis for every variable
#---------------------------------------------------------------------------
### Settings
## Set working directory
setwd("/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Visualize DEG/")

## load packages
library(Seurat)   #3.2.3
library(DESeq2)   #1.28.1
library(ggplot2)  #3.3.3

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


## Look up sc DEG in pseudobulk
comparison <- lapply(testing.names, function(cluster){
  
  # Select sc DEG names
  sc <- rownames(sc.DEG[[cluster]])
  
  # Look up sc DEG names in pseudobulk  
  pb <- pb.DEG[[cluster]][sc,]
  return(pb)
})

names(comparison) <- testing.names


## Order comparison by adjusted p-value
comparison <- lapply(testing.names, function(cluster){
  pb.ordered <- as.data.frame(comparison[[cluster]][order(comparison[[cluster]]$padj),])
  return(pb.ordered)
})

names(comparison) <- testing.names


## Replace NA values with 1
for (cluster in testing.names){
  comparison[[cluster]][is.na(comparison[[cluster]])] <- 1
}


## Select comparison DEG by significant p-value and make dataframe
comparison.sign <- lapply(testing.names, function(cluster){
  pb.sign <- as.data.frame(comparison[[cluster]][comparison[[cluster]]$pvalue < 0.05,])
  
  ## Add column of gene names and column in which cluster the DEG is more highly expressed
  # Positive logfoldchange corresponds to cluster1
  cluster1.subset <- subset(pb.sign, log2FoldChange > 0) 
  
  # If there are avg_logFC > 0 the genes correspond to cluster1
  if (dim(cluster1.subset)[1] != 0) {
    cluster1.subset$cluster <- testing.clusters[[cluster]][1]
    cluster1.subset$gene <- rownames(cluster1.subset)
  }
  
  
  # Negative logfoldchange corresponds to cluster2
  cluster2.subset <- subset(pb.sign, log2FoldChange < 0) 
  
  # If there are avg_logFC < 0 the genes correspond to cluster2
  if (dim(cluster2.subset)[1] != 0) {
    cluster2.subset$cluster <- testing.clusters[[cluster]][2]
    cluster2.subset$gene <- rownames(cluster2.subset)
  }
  
  result <- rbind(cluster1.subset, cluster2.subset) 
  
  return(result)
})

names(comparison.sign) <- testing.names

#--------------------------------------------------------------------------
### Plots of significant DEG for every cell cluster per variable

## Preparing data
# Significant gene names for every cluster comparison
gene.names <- lapply(testing.names, function(comparison){
  rownames(comparison.sign[[comparison]])
})

names(gene.names) <- testing.names

# Average expression of significant genes per cluster comparison
avg.exp <- lapply(testing.names, function(comparison){
  if(length(gene.names[[comparison]]) != 0){
    AverageExpression(seurdata, assay = "SCT", features = gene.names[[comparison]], slot = "data", return.seurat = T)
  }
})

names(avg.exp) <- testing.names

# Make dataframes of average expression of SCT scale data slot per cluster comparison
avg.exp <- lapply(testing.names, function(cluster){
  if(length(avg.exp[[cluster]][["SCT"]]) != 0){
    exp <- avg.exp[[cluster]][["SCT"]]@scale.data[, c(strsplit(cluster, split = " vs. ")[[1]][1], strsplit(cluster, split = " vs. ")[[1]][2])]
  }
})

names(avg.exp) <- testing.names

# Combine both variable columns of average expression into one large column per cluster comparison
avg.exp <- lapply(testing.names, function(cluster){
  if(length(avg.exp[[cluster]]) > 2){
  data.frame(c(avg.exp[[cluster]][, strsplit(cluster, split = " vs. ")[[1]][1]], 
               avg.exp[[cluster]][, strsplit(cluster, split = " vs. ")[[1]][2]]))
  }
})

names(avg.exp) <- testing.names


# Name avg.exp list elements
for (i in c(1:length(avg.exp))){
  if(length(avg.exp[[i]] != 0)){
    names(avg.exp[[i]]) <- "avg.exp"
  }
}  

# Percentage expression of significant genes per cluster comparison
pct <- lapply(testing.names, function(comparison){
  data.frame(c(sc.DEG[[comparison]][gene.names[[comparison]], "pct.1"], sc.DEG[[comparison]][gene.names[[comparison]], "pct.2"]))
})

names(pct) <- testing.names

# Name pct list elements
for (i in c(1:length(avg.exp))){
  if(length(avg.exp[[i]] != 0)){
    names(pct[[i]]) <- "percent"
  }
}  

# Remove avg.exp elements that are NULL
avg.exp <- avg.exp[lengths(avg.exp) != 0]

## Dataframe of significant genes per cluster comparison
df <- lapply(names(avg.exp), function(comparison){
  data.frame(
    rep((rownames(sc.DEG[[comparison]][gene.names[[comparison]],])), times = length(testing.clusters[[comparison]])), # Column of gene names
    pct[[comparison]],  # Column of percentages of cells in which genes are expressed
    avg.exp[[comparison]],  # Column of average expression of genes
    rep(c(strsplit(comparison, split = " vs. ")[[1]][1],
          strsplit(comparison, split = " vs. ")[[1]][2]), 
        each = length(pct[[comparison]][["percent"]])/2))  # In which parameter is the gene more highly expressed?
})

names(df) <- names(avg.exp)

# Name df list elements
for (i in c(1:length(df))){
  if(length(df[[i]] == 4)){
    names(df[[i]]) <- c("gene", "pct", "avg.exp", "parameter")
  }
} 

## Dotplot/bubble plot of significant genes per cluster comparison
dotplots <- lapply(names(avg.exp), function(cluster){
  if(length(df[[cluster]] == 4)){
    ggplot(data = df[[cluster]], aes(x = parameter, y = gene)) +
      geom_point(aes(colour = avg.exp, size = pct)) +   # Make into dotplot
      labs(colour = "Average Expression") +   # Label scale bar
      scale_color_gradient2(low = "#fbfdff", mid = "#a7a7de", high = "#08085d") +   # Colors of scale bar
      scale_size_area("Expressed in fraction of cells", breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) +   # Size of the dots for fraction of cells in which they are expressed
      theme(panel.background = element_blank()) +    # Remove background
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Text on x-axis in 45 degree angle
      ggtitle(label = paste(cluster))   # Add title
  }
})

names(dotplots) <- names(avg.exp)

#--------------------------------------------------------------------------
### Customize plots

## Plot of SMCs with genes colored by sex
# Get dataframe of genes and sex in which they are more highly expressed
genes.color <- comparison.sign[["male Smooth Muscle Cells vs. female Smooth Muscle Cells"]][, c("gene", "cluster")]

# Order genes
genes.color <- genes.color[order(genes.color$cluster),]

## Make sure ggplot order matches gene order
df[["male Smooth Muscle Cells vs. female Smooth Muscle Cells"]]$gene <- factor(df[["male Smooth Muscle Cells vs. female Smooth Muscle Cells"]]$gene, levels = rownames(genes.color))

# Plot
ggplot(data = df[["male Smooth Muscle Cells vs. female Smooth Muscle Cells"]], aes(x = parameter, y = gene)) +
  geom_point(aes(colour = avg.exp, size = pct)) +   # Make into dotplot
  labs(colour = "Avg Exp") +   # Label scale bar
  scale_color_gradient2(low = "#fbfdff", mid = "#a7a7de", high = "#08085d") +   # Colors of scale bar
  scale_size_area("Fraction", breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) +   # Size of the dots for fraction of cells in which they are expressed
  theme(panel.background = element_blank()) +    # Remove background
  theme(axis.text.y = element_text(colour = ifelse(genes.color$cluster == "male Smooth Muscle Cells", "#4ec900", "#FF0000"), size = 12)) +  # Give color to genes corresponding to sex
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) + # Text on x-axis in 45 degree angle
  scale_x_discrete(labels = c("female SMCs", "male SMCs")) + # Custom x group labels
  ggtitle(label = "female vs. male Smooth Muscle Cells")   # Add title



## Plot of ECIs with genes colored by sex
# Get dataframe of genes and sex in which they are more highly expressed
genes.color <- comparison.sign[["male Endothelial Cells I vs. female Endothelial Cells I"]][, c("gene", "cluster")]

# Order genes
genes.color <- genes.color[order(genes.color$cluster),]

## Make sure ggplot order matches gene order
df[["male Endothelial Cells I vs. female Endothelial Cells I"]]$gene <- factor(df[["male Endothelial Cells I vs. female Endothelial Cells I"]]$gene, levels = rownames(genes.color))

# Plot
ggplot(data = df[["male Endothelial Cells I vs. female Endothelial Cells I"]], aes(x = parameter, y = gene)) +
  geom_point(aes(colour = avg.exp, size = pct)) +   # Make into dotplot
  labs(colour = "Avg Exp") +   # Label scale bar
  scale_color_gradient2(low = "#fbfdff", mid = "#a7a7de", high = "#08085d") +   # Colors of scale bar
  scale_size_area("Fraction", breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) +   # Size of the dots for fraction of cells in which they are expressed
  theme(panel.background = element_blank()) +    # Remove background
  theme(axis.text.y = element_text(colour = ifelse(genes.color$cluster == "male Endothelial Cells I", "#4ec900", "#FF0000"), size = 12)) +  # Give color to genes corresponding to sex
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) + # Text on x-axis in 45 degree angle
  scale_x_discrete(labels = c("female ECIs", "male ECIs")) + # Custom x group labels
  ggtitle(label = "female vs. male Endothelial Cells I")   # Add title




## Plot of ECIIs with genes colored by sex
# Get dataframe of genes and sex in which they are more highly expressed
genes.color <- comparison.sign[["male Endothelial Cells II vs. female Endothelial Cells II"]][, c("gene", "cluster")]

# Order genes
genes.color <- genes.color[order(genes.color$cluster),]

## Make sure ggplot order matches gene order
df[["male Endothelial Cells II vs. female Endothelial Cells II"]]$gene <- factor(df[["male Endothelial Cells II vs. female Endothelial Cells II"]]$gene, levels = rownames(genes.color))

# Plot
ggplot(data = df[["male Endothelial Cells II vs. female Endothelial Cells II"]], aes(x = parameter, y = gene)) +
  geom_point(aes(colour = avg.exp, size = pct)) +   # Make into dotplot
  labs(colour = "Avg Exp") +   # Label scale bar
  scale_color_gradient2(low = "#fbfdff", mid = "#a7a7de", high = "#08085d") +   # Colors of scale bar
  scale_size_area("Fraction", breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) +   # Size of the dots for fraction of cells in which they are expressed
  theme(panel.background = element_blank()) +    # Remove background
  theme(axis.text.y = element_text(colour = ifelse(genes.color$cluster == "male Endothelial Cells II", "#4ec900", "#FF0000"), size = 12)) +  # Give color to genes corresponding to sex
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) + # Text on x-axis in 45 degree angle
  scale_x_discrete(labels = c("female ECIIs", "male ECIIs")) + # Custom x group labels
  ggtitle(label = "female vs. male Endothelial Cells II")   # Add title
