library(scCustomize)
library(ComplexHeatmap)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(tidyverse)
library(limma) 
library(edgeR) 
library("sva")
library(dplyr)
library("writexl")
library(gplots)
library('xlsx')
library('readxl')
library(ggplot2)

##### Read in the KTA data.
scRNAseqFile <- "C:/Users/MrBes/Documents/R/KidneyTissueAtlas/521c5b34-3dd0-4871-8064-61d3e3f1775a_PREMIERE_Alldatasets_08132021.h5Seurat"

# Define a list of genes of interest
genes <- c("ADGRL1", "PROM1", "PLS3", "EPHB2", "NECTIN1", "PDCD1LG2",
           "SCGB1A1", "FAM151A", "LCN1", "CLIP2", "PLXDC2", "PTPRJ")

# Convert gene names to uppercase
genes <- toupper(genes)

# Single Cell
# Load Seurat object from the H5Seurat file
scData <- LoadH5Seurat(scRNAseqFile, assays = "data")
scMetadata <- scData@meta.data

# Define the order of cell types for plotting
scOrder <- c('EC', 'POD', 'PEC', 'PT', 'DTL', 'ATL/TAL', 'TAL', 'DCT', 'CNT', 'PC', 'IC', 'Immune', 'Interstitial')
scData@meta.data$subclass.l1 <- factor(scData@meta.data$subclass.l1, scOrder)

# Remove genes not found in scRNA-seq
genes <- intersect(genes, row.names(scData))

# Create a dataframe to store mean intensity values for each gene within each cell type
df <- data.frame(matrix(nrow = length(genes), ncol = length(scOrder)))
names(df) <- scOrder
rownames(df) <- genes

# Iterate through cell types and genes to establish the mean Intensity for each gene within each cell type
for (j in 1:length(scOrder)) {
  ndxCellType <- scData@meta.data$subclass.l1 == scOrder[j] & scData@meta.data$sampletype == "LD"
  dfSmall <- scData[genes, ndxCellType]
  meanVals <- rowMeans(dfSmall)
  df[, j] <- meanVals
}

geneSums <- rowSums(df)

# Normalize the mean intensity values
dfUniq <- data.frame(matrix(nrow = length(genes), ncol = length(scOrder)))
names(dfUniq) <- scOrder
rownames(dfUniq) <- genes
for (j in 1:length(scOrder)) {
  for (jj in 1:length(genes)) {
    dfUniq[jj, j] <- df[jj, j] / geneSums[jj]
  }
}

# Prepare data for heatmap
dataScaled <-  as.matrix(dfUniq)
toKeep <- !rowSums(is.na(dataScaled))
genes <- rownames(dataScaled)[toKeep]
dataScaled <- dataScaled[toKeep,]
dataScaled <- t(scale(t(dataScaled)))

# Supervised clustering
dfMaxUniq <- data.frame(matrix(nrow = length(genes), ncol = 3))
colnames(dfMaxUniq) <- c("Gene", "CellType", "MaxUniq")

for (j in 1:length(genes)) {
  thisGene <- genes[j]
  thisNdx <- rownames(dataScaled) == thisGene
  thisRow <- dataScaled[thisNdx,]
  thisMaxNdx <- which.max(thisRow)
  
  thisMax <- thisRow[thisMaxNdx]
  thisCellType <- as.character(colnames(dataScaled)[thisMaxNdx])
  dfMaxUniq[j, 1] <- thisGene
  dfMaxUniq[j, 2] <- thisCellType
  dfMaxUniq[j, 3] <- thisMax
}

# Sort genes by the maximum intensity within each cell type
for (j in 1:length(scOrder)) {
  ndxs <- which(dfMaxUniq$CellType == scOrder[j])
  dfSub <- dfMaxUniq[ndxs,]
  dfSub <- dfSub[order(dfSub$MaxUniq, decreasing = TRUE),]
  if (j == 1) {
    dfAll <- dfSub
  } else {
    dfAll <- rbind(dfAll,dfSub)
  }
}
geneOrder <- rev(dfAll$Gene)
dataScaled <- dataScaled[geneOrder,]

# Define a color palette
my_palette <- colorRampPalette(c("#5A5AA0", 'white', "#FFB600"))(n = 100)

# Plot heatmap with heatmap.2
par(cex.main=1) # Shrink title fonts on plot

# Plot heatmap
setwd('C:\\Users\\MrBes\\Documents\\R\\KidneyTissueAtlas\\covidAKI\\Uniqueness')
pdf(file = 'scUnique_Healthy_scaledRows_supervised_nephronorder_CovidFeatures.pdf',height = 7, width = 7)
gplots::heatmap.2(dataScaled,                     # Tidy, normalized data
                  Colv = FALSE,                   # Experiments clusters in cols
                  Rowv = FALSE,                   # Protein clusters in rows
                  revC = TRUE,                    # Flip plot to match pheatmap
                  density.info = "histogram",    # Plot histogram of data and color key
                  trace = "none",                # Turn off trace lines from heat map
                  col = my_palette,              # Use custom color scheme
                  cexRow = .7, cexCol = 1,
                  keysize = 1,
                  margins = c(10, 7))
dev.off()