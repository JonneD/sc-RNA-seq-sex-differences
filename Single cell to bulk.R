### 2021-06-29
# Project SCS 46 patients - Single cell to bulk
# Jonne Damen
# R version 4.0.2 (2020-06-22) -- "Taking Off Again"
# R Studio version 1.3.1056 "Water Lily" (5a4dee98, 2020-07-07) for macOS

### Looking up DEGs of the single cell analysis in the bulk data
#---------------------------------------------------------------------------
### Settings
## Set working directory
setwd("/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Single cell to bulk/")

## Load package
library(Seurat) # 3.2.3
library(ggplot2) # 3.3.3
library(ggpubr) # 0.4.0

## Load data
# Bulk RNA seq data
bulkdata <- readRDS("/Volumes/DATASHUR_P2/Project SCS 46 patients/Data/bulk_data_for_Jonne.RDS")
#--------------------------------------------------------------------------
## Put metadata in global environment
meta.data <- bulkdata@meta.data

## Create result folder
#result.folder <- paste(Sys.Date(), "result", sep = " ")
#dir.create(result.folder, showWarnings = FALSE)
#--------------------------------------------------------------------------
### Project sc genes onto bulk split by sex

## Overview of male biased gene bulk expression
# What genes are male biased in the sc analysis in the SMC and EC clusters and are present in the bulkdata?
genes <- c("ZNF480", "ZNF471", "ZFP36L1", "ZFP36", "XAF1", "TMEM212", "RAB12",
           "ORAI2", "ODF2L", "JUNB", "HMGB1", "FOSB", "FOS", "F5", "CFLAR",
           "SOCS3", "PTGDS", "NR2F2", "KLF15", "ENG", "EGR1", "RPS18", "PTMA",
           "KLF6", "HLA-DRB5", "COL1A1")

# Get bulk expression data of male biased genes
exp.data <- FetchData(bulkdata, vars = genes)
# Make dataframe of expression data
exp.data <- data.frame(unlist(exp.data))
# Name column of dataframe
names(exp.data) <- "expression"
# Add column of sex of bulk patient
exp.data$sex <- rep(meta.data$sex, times = length(genes))
# Add column of gene name
exp.data$gene <- rep(genes, each = length(meta.data$sex))

# Make boxplots of bulk RNA expression of selected male biased genes
#pdf(paste(result.folder, "/", Sys.Date(), "_plot_male_biased.pdf", sep = ""), width = 15, height = 15)
ggplot(exp.data, aes(x = gene, y = expression, fill = sex)) + 
  geom_boxplot(outlier.size = 0.5) + # Make boxplot
  geom_point(position = position_dodge(width = 0.75), size = 0.5) + # Add individual observations to boxplot
  facet_wrap(~gene, scale = "free") + # Make a separate plot for every gene
  stat_compare_means(method = "wilcox.test") + # Compare means of observations
  ggtitle("bulk RNA expression of male biased genes from scRNAseq analysis") # Add title
#dev.off()



## Overview of female biased gene bulk expression

# What genes are female biased in the sc analysis in the SMC and EC clusters and are present in the bulkdata?
genes <- c("TYROBP", "TLN1", "IGLL5", "CD74", "PLS3", "HLA-B", "HLA-A",
           "GNAQ", "DKK3", "B2M", "ZBTB20", "TAGLN", "IL2RG", "IFITM2",
           "IFITM1", "HDGF")

# Get bulk expression data of female biased genes
exp.data <- FetchData(bulkdata, vars = genes)
# Make dataframe of expression data
exp.data <- data.frame(unlist(exp.data))
# Name column of dataframe
names(exp.data) <- "expression"
# Add column of sex of bulk patient
exp.data$sex <- rep(meta.data$sex, times = length(genes))
# Add column of gene name
exp.data$gene <- rep(genes, each = length(meta.data$sex))

# Make boxplots of bulk RNA expression of selected female biased genes
#pdf(paste(result.folder, "/", Sys.Date(), "_plot_female_biased.pdf", sep = ""), width = 15, height = 15)
ggplot(exp.data, aes(x = gene, y = expression, fill = sex)) + 
  geom_boxplot(outlier.size = 0.5) + # Make boxplot
  geom_point(position = position_dodge(width = 0.75), size = 0.5) + # Add individual observations to boxplot
  facet_wrap(~gene, scale = "free") + # Make a separate plot for every gene
  stat_compare_means(method = "wilcox.test") + # Compare means of observations
  ggtitle("bulk RNA expression of female biased genes from scRNAseq analysis") # Add title
#dev.off()

#--------------------------------------------------------------------------
### Project sc genes onto bulk clusters split by sex

## Overview of male biased gene bulk expression
# What genes are male biased in the sc analysis in the SMC and EC clusters and are present in the bulkdata?
genes <- c("ZNF480", "ZNF471", "ZFP36L1", "ZFP36", "XAF1", "TMEM212", "RAB12",
           "ORAI2", "ODF2L", "JUNB", "HMGB1", "FOSB", "FOS", "F5", "CFLAR",
           "SOCS3", "PTGDS", "NR2F2", "KLF15", "ENG", "EGR1", "RPS18", "PTMA",
           "KLF6", "HLA-DRB5", "COL1A1")

# Make plot of gene expression per bulk cluster for every gene separately
#pdf(paste(result.folder, "/", Sys.Date(), "_plots_male_per_bulk_cluster.pdf", sep = ""))
lapply(genes, function(gene){
  # Get bulk expression data of male biased genes
  exp.data <- FetchData(bulkdata, vars = gene)
  # Make dataframe of expression data
  exp.data <- data.frame(unlist(exp.data))
  # Name column of dataframe
  names(exp.data) <- "expression"
  # Add column of sex of bulk patient
  exp.data$sex <- rep(meta.data$sex, times = length(gene))
  # Add column of gene name
  exp.data$gene <- rep(gene, each = length(meta.data$sex))
  # Add column of bulk cluster
  exp.data$cluster <- rep(meta.data$seurat_clusters, times = length(gene))  
  
  # Make boxplots of bulk RNA expression of every gene per bulk cluster
  ggplot(exp.data, aes(x = gene, y = expression, fill = sex)) + 
    geom_boxplot(outlier.size = 0.5) + # Make boxplot
    geom_point(position = position_dodge(width = 0.75), size = 0.5) + # Add individual observations to boxplot
    facet_wrap(~cluster, scale = "free") + # Make a separate plot for every cluster
    stat_compare_means(method = "wilcox.test") + # Compare means of observations
    ggtitle(paste("bulk RNA expression of", gene, "per bulk cluster", sep = " ")) # Add title
})
#dev.off()

## Overview of female biased gene bulk expression
# What genes are female biased in the sc analysis in the SMC and EC clusters and are present in the bulkdata?
genes <- c("TYROBP", "TLN1", "IGLL5", "CD74", "PLS3", "HLA-B", "HLA-A",
           "GNAQ", "DKK3", "B2M", "ZBTB20", "TAGLN", "IL2RG", "IFITM2",
           "IFITM1", "HDGF")

#pdf(paste(result.folder, "/", Sys.Date(), "_plots_female_per_bulk_cluster.pdf", sep = ""))
lapply(genes, function(gene){
  # Get bulk expression data of male biased genes
  exp.data <- FetchData(bulkdata, vars = gene)
  # Make dataframe of expression data
  exp.data <- data.frame(unlist(exp.data))
  # Name column of dataframe
  names(exp.data) <- "expression"
  # Add column of sex of bulk patient
  exp.data$sex <- rep(meta.data$sex, times = length(gene))
  # Add column of gene name
  exp.data$gene <- rep(gene, each = length(meta.data$sex))
  # Add column of bulk cluster
  exp.data$cluster <- rep(meta.data$seurat_clusters, times = length(gene)) 
  
  # Make boxplots of bulk RNA expression of every gene per bulk cluster
  ggplot(exp.data, aes(x = gene, y = expression, fill = sex)) + 
    geom_boxplot(outlier.size = 0.5) + # Make boxplot
    geom_point(position = position_dodge(width = 0.75), size = 0.5) + # Add individual observations to boxplot
    facet_wrap(~cluster, scale = "free") + # Make a separate plot for every cluster
    stat_compare_means(method = "wilcox.test") + # Compare means of observations
    ggtitle(paste("bulk RNA expression of", gene, "per bulk cluster", sep = " ")) # Add title
})
#dev.off()
