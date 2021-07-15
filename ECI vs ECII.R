### 2021-05-12
# Project SCS 46 patients - ECI vs ECII
# Jonne Damen
# R version 4.0.2 (2020-06-22) -- "Taking Off Again"
# R Studio version 1.3.1056 "Water Lily" (5a4dee98, 2020-07-07) for macOS

# Finding out what separates ECI and ECII
#--------------------------------------------------------------------------
### Settings
## Set working directory
setwd("/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/ECI vs ECII")

## Load packages
library(Seurat)   #3.2.3

## Load data
seurdata <- readRDS("/Volumes/DATASHUR_P2/Project SCS 46 patients/Final scripts/Correct counts for sc analysis/20210420_seurdata_corrected_cellcounts.RDS")
#--------------------------------------------------------------------------
### Calculate differential gene expression

## Set correct idents
Idents(seurdata) <- seurdata@meta.data$new.ident

## Calculate DEG between ECI and ECII
difference <- FindMarkers(seurdata, ident.1 = "Endothelial Cells I", ident.2 = "Endothelial Cells II")
#--------------------------------------------------------------------------
### DEG of ECI
## Average logfoldchange of > 0 corresponds to ECI
ECI <- difference[difference$avg_logFC > 0,]

## Filter on adjusted p-value
ECI <- ECI[ECI$p_val_adj < 0.05,]

## Get gene names
ECI.names <- rownames(ECI)

## Save gene names
#write.table(ECI.names, file = "ECI vs ECII.txt", sep = "\t", row.names = F, quote = F)
#--------------------------------------------------------------------------
### DEG of ECII
## Average logfoldchange of < 0 corresponds to ECII
ECII <- difference[difference$avg_logFC < 0,]

## Filter on adjusted p-value
ECII <- ECII[ECII$p_val_adj < 0.05,]

## Get gene names
ECII.names <- rownames(ECII)

## Save gene names
#write.table(ECII.names, file = "ECII vs ECI.txt", sep = "\t", row.names = F, quote = F)
#--------------------------------------------------------------------------
### Pathway (overrepresentation) analysis

## Load packages
library(ReactomePA) #1.32.0
library(biomaRt)    #2.44.4
library(org.Hs.eg.db) #3.12.0

## Select dataset and BioMart database to pull information from
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#--------------------------------------------------------------------------
### ECI pathway analysis
## Get gene names in Entrez ID
DEGs_entrez <- getBM( # Get information from the BioMart database
  filters = "hgnc_symbol", # What format are our gene names in?
  attributes = c("hgnc_symbol", 'ensembl_gene_id', 'entrezgene_id'), # What formats do we want to obtain for our gene names
  values= ECI.names, # For which genes do we want to convert gene names?
  mart = mart # From which part of the BioMart database do we want to pull information?
)

## Set random seed
set.seed(456)

## Pathway Enrichment Analysis of ECI
PA.ECI <- enrichPathway(gene= DEGs_entrez$entrezgene_id ,pvalueCutoff=0.05, readable=T)

## Dotplot of enriched pathways
dotplot(PA.ECI,
        showCategory= 10, 
        title = "Pathway Enrichment of ECI")
#--------------------------------------------------------------------------
### ECII pathway analysis
## Get gene names in Entrez ID
DEGs_entrez <- getBM( # Get information from the BioMart database
  filters = "hgnc_symbol", # What format are our gene names in
  attributes = c("hgnc_symbol", 'ensembl_gene_id', 'entrezgene_id'), # What formats do we want to obtain for our gene names
  values= ECII.names, # For which genes do we want to convert gene names
  mart = mart # From which part of the BioMart database do we want to pull information?
)

## Pathway Enrichment Analysis of ECII
PA.ECII <- enrichPathway(gene= DEGs_entrez$entrezgene_id ,pvalueCutoff=0.05, readable=T)

## Dotplot of enriched pathways
dotplot(PA.ECII,
        showCategory= 10, 
        title = "Pathway Enrichment of ECII")

