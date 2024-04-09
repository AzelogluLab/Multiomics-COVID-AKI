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
# Define the path to the H5Seurat file containing single-cell RNA-seq data
scRNAseqFile <- "C:/Users/MrBes/Documents/R/KidneyTissueAtlas/521c5b34-3dd0-4871-8064-61d3e3f1775a_PREMIERE_Alldatasets_08132021.h5Seurat"

# Single Cell
# Load Seurat object from the H5Seurat file
scData <- LoadH5Seurat(scRNAseqFile, assays = "data")
scMetadata <- scData@meta.data

# Define the order of cell types for plotting
scOrder <- c('EC','POD','PEC','PT','DTL','ATL/TAL','TAL','DCT','CNT','PC','IC','Immune','Interstitial')
scData@meta.data$subclass.l1 <- factor(scData@meta.data$subclass.l1, scOrder)


## Generate the genes
# Function to retrieve genes associated with a given pathway from a specified Excel file
getGeneList <- function(pathwayName) {
  dfFO <- read_excel("C:/Users/MrBes/Documents/R/KidneyTissueAtlas/covidAKI/Functional overlap_KPMP.xlsx", col_names = TRUE)
  thisPathwayNdx <- which(dfFO$...1 == pathwayName)
  theseGenesUnf <- dfFO[thisPathwayNdx, 2:16] # Exclude organoid-related genes
  allGenes <- c()
  for (j in 1:15) {
    parts <- unlist(strsplit(theseGenesUnf[[j]], "\\|\\|"))
    theseGenesFor <- unlist(strsplit(parts[2], ", "))
    allGenes <- c(allGenes, theseGenesFor)
  }
  # Reduce to genes found in 2 or more sources (appear more than once in the list)
  geneCounts <- table(allGenes)
  genes <- names(geneCounts[geneCounts >=2])
  return(genes)
}

# Retrieve genes associated with three pathways
genes1 <- getGeneList("GO:0043062~extracellular structure organization")
genes2 <- getGeneList("GO:0030198~extracellular matrix organization")
genes3 <- getGeneList("GO:0022617~extracellular matrix disassembly")

# Combine genes from all pathways and remove duplicates
genesToConsider <- c(genes1, genes2, genes3)
genesToConsider <- unique(genesToConsider)

genesToRemove <- c()

if (length(genesToRemove) > 0) {
  # Remove specified genes
  for (j in 1:length(genesToRemove)) {
    thisNdx <- which(genesToConsider == genesToRemove[j])
    genesToConsider <- genesToConsider[-thisNdx]
  }
}

# Plot DotPlot to visualize gene expression across cell types
A <- DotPlot(subset(scData, subset = diseasetype == "LivingDonor"), 
             features = genesToConsider, group.by = 'subclass.l1', scale.by = 'size') + 
  RotatedAxis() + coord_flip() + scale_colour_gradient2(low = "#5A5AA0", mid = "white", high = "#FFB600")

# Retrieve data from the DotPlot
dfDot <- A$data

# Prepare dataframe to store maximum scaled expression and percentage expressed for each gene
dfMaxScaledExp <- data.frame(matrix(nrow = length(genesToConsider), ncol = 4))
colnames(dfMaxScaledExp) <- c("Gene", "CellType", "MaxScaledExp", "PctExpressed")

# Iterate through genes to determine maximum scaled expression and percentage expressed
for (j in 1:length(genesToConsider)) {
  thisGene <- genesToConsider[j]
  dfTemp <- dfDot[dfDot$features.plot == thisGene,]
  thisMaxNdx <- which.max(dfTemp$avg.exp.scaled)
  thisMax <- dfTemp$avg.exp.scaled[thisMaxNdx]
  thisCellType <- as.character(dfTemp$id[thisMaxNdx])
  thisPctExp <- dfTemp$pct.exp[thisMaxNdx]
  dfMaxScaledExp[j,1] <- thisGene
  dfMaxScaledExp[j,2] <- thisCellType
  dfMaxScaledExp[j,3] <- thisMax
  dfMaxScaledExp[j,4] <- thisPctExp
}

# Sort genes by maximum scaled expression within each cell type
dfMaxScaledExp$CellType <- factor(dfMaxScaledExp$CellType, scOrder)
dfMaxScaledExp <- dfMaxScaledExp[order(dfMaxScaledExp$CellType),]
for (j in 1:length(scOrder)) {
  thisCellType <- scOrder[j]
  ndxs <- which(dfMaxScaledExp$CellType == thisCellType)
  dfSub <- dfMaxScaledExp[ndxs,]
  uniqExp <- sort(unique(dfSub$MaxScaledExp), decreasing = TRUE)
  for (jExp in 1:length(uniqExp)) {
    jExpNdx <- dfSub$MaxScaledExp == uniqExp[jExp]
    dfSubSub <- dfSub[jExpNdx,]
    dfSubSub <- dfSubSub[order(dfSubSub$PctExpressed, decreasing = TRUE),]
    if (jExp == 1 & j == 1) {
      dfAll <- dfSubSub
    }
    else {
      dfAll <- rbind(dfAll, dfSubSub)
    }
  }
}
geneListSupervised <- rev(dfAll$Gene)

# Prepare data for heatmap
genes <- genesToConsider
df <- data.frame(matrix(nrow = length(genes), ncol = length(scOrder)))
names(df) <- scOrder
rownames(df) <- genes

# Iterate through cell types and genes to establish the mean intensity for each gene within each cell type
for (j in 1:length(scOrder)) {
  ndxCellType <- scData@meta.data$subclass.l1 == scOrder[j] & scData@meta.data$sampletype == "LD"
  dfSmall <- scData[genes, ndxCellType]
  meanVals <- rowMeans(dfSmall)
  df[,j] <- meanVals
}
geneSums <- rowSums(df)

# Normalize the mean intensity values
dfUniq <- data.frame(matrix(nrow = length(genes), ncol = length(scOrder)))
names(dfUniq) <- scOrder
rownames(dfUniq) <- genes
for (j in 1:length(scOrder)) {
  for(jj in 1:length(genes)) {
    dfUniq[jj,j] <- df[jj,j] / geneSums[jj]
  }
}

dataScaled <- as.matrix(dfUniq)
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
  dfMaxUniq[j,1] <- thisGene
  dfMaxUniq[j,2] <- thisCellType
  dfMaxUniq[j,3] <- thisMax
}

# Sort genes by maximum scaled expression within each cell type
for (j in 1:length(scOrder)) {
  ndxs <- which(dfMaxUniq$CellType == scOrder[j])
  dfSub <- dfMaxUniq[ndxs,]
  dfSub <- dfSub[order(dfSub$MaxUniq, decreasing = TRUE),]
  if (j == 1) {
    dfAll <- dfSub
  } else {
    dfAll <- rbind(dfAll, dfSub)
  }
}
geneOrder <- rev(dfAll$Gene)
dataScaled <- dataScaled[geneOrder,]

my_palette <- colorRampPalette(c("#5A5AA0", 'white', "#FFB600"))(n=100)

# Plot heatmap with heatmap.2
par(cex.main=1) # Shrink title fonts on plot

# Plot heatmap
setwd('C:\\Users\\MrBes\\Documents\\R\\KidneyTissueAtlas\\covidAKI\\Uniqueness')
pdf(file = 'scUnique_Healthy_ECMRelated_supervised_nephronorder_withOrganoids.pdf', height = 7, width = 7)
gplots::heatmap.2(dataScaled,                    # Tidy, normalized data
                  Colv = FALSE,                 # Experiments clusters in columns
                  Rowv = FALSE,                 # Protein clusters in rows
                  revC = TRUE,                  # Flip plot to match pheatmap
                  density.info = "histogram",   # Plot histogram of data and color key
                  trace = "none",               # Turn off trace lines from heat map
                  col = my_palette,             # Use custom color scheme
                  cexRow = .7, cexCol = 1,
                  keysize = .75,
                  margins = c(10, 7))
dev.off()