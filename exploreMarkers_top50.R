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

# Read dataframes containing genes
df1 <- read_excel("C:\\Users\\MrBes\\Documents\\R\\KidneyTissueAtlas\\covidAKI\\Uniqueness\\DEG_AKI_pro.xlsx")
df2 <- read_excel("C:\\Users\\MrBes\\Documents\\R\\KidneyTissueAtlas\\covidAKI\\Uniqueness\\covid_proteomics_Composite_Outcome_Y_vs_N_unadj_limma_all.xlsx")

# Order dataframes by p-values
df1 <- df1[order(df1$P.Value, decreasing = FALSE), ]
df2 <- df2[order(df2$p, decreasing = FALSE), ]

# Combine top genes from both dataframes and remove specific genes
covidGenes <- c(df1$Gene[1:50], df2$Symbol[1:50])
toRemove <- which(covidGenes %in% c("PRSS2", "PI3", "CPLX2", "REG3A"))
covidGenesOrdered <- covidGenes[-toRemove]

# Read ordered gene list
covidGenesOrder <- read_excel("C:/Users/MrBes/Documents/R/KidneyTissueAtlas/covidAKI/top50clustergeneorder.xlsx")
covidGenesOrder <- toupper(covidGenesOrder$Gene)

# Read additional gene lists
thisListOrdered <- read_excel("C:/Users/MrBes/Documents/R/KidneyTissueAtlas/covidAKI/MigrationGenes.xlsx")
thisListOrdered <- toupper(thisListOrdered$Gene)
covidGenesOrdered <- thisListOrdered

thisListOrderedOxygen <- read_excel("C:/Users/MrBes/Documents/R/KidneyTissueAtlas/covidAKI/OxygenGenesOrderSC.xlsx")
thisListOrderedOxygen <- toupper(thisListOrderedOxygen$Gene)
covidGenesOrdered <- thisListOrderedOxygen

thisList <- read_excel("C:\\Users\\MrBes\\Documents\\R\\KidneyTissueAtlas\\covidAKI\\CellAdhesionGenes.xlsx")
thisList <- toupper(thisList$Gene)
covidGenesOrdered <- thisList

# Plot clustered dot plot
p <- Clustered_DotPlot(subset(scData, subset = diseasetype == "LivingDonor"),
                       features = covidGenesHealthy, group.by = 'subclass.l1', flip = FALSE,
                       cluster_ident = FALSE)
pdf(file = 'scRNA_COVID_Healthy_clustered_2.pdf', 7, 15)
p[[2]]
dev.off()

# Plot DotPlot with ordered genes
p <- Clustered_DotPlot(scData,
                       features = covidGenesHealthy, group.by = 'subclass.l1', flip = FALSE,
                       cluster_ident = FALSE)
A <- p[[2]]
covidGenesOrdered <- covidGenes[A@ht_list$matrix_4@row_order]
dev.off()

# Plot DotPlot with top features
pdf(file = 'scRNA_COVID_Healthy_top_features_nephronorder_4.pdf', 7, 4)
A <- DotPlot(subset(scData, subset = diseasetype == "LivingDonor"), features = covidGenesOrdered, group.by = 'subclass.l1', scale.by = 'size') + RotatedAxis() + coord_flip() + scale_colour_gradient2(low = "#5A5AA0", mid = "white", high = "#FFB600")
dev.off()

# Use info from dotplot to perform supervised clustering
dfDot <- A$data
dfMaxScaledExp <- data.frame(matrix(nrow = length(covidGenesOrdered), ncol = 4))
colnames(dfMaxScaledExp) <- c("Gene","CellType","MaxScaledExp","PctExpressed")
for (j in 1:length(covidGenesOrdered)) {
  thisGene <- covidGenesOrdered[j]
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
dfMaxScaledExp$CellType <- factor(dfMaxScaledExp$CellType, scOrder)
dfMaxScaledExp <- dfMaxScaledExp[order(dfMaxScaledExp$CellType),]
for (j in 1:length(scOrder)) {
  thisCellType <- scOrder[j]
  ndxs <- which(dfMaxScaledExp$CellType == thisCellType)
  dfSub <- dfMaxScaledExp[ndxs,]
  # Add another loop to break ties in expression with % expressed. 
  uniqExp <- sort(unique(dfSub$MaxScaledExp), decreasing = TRUE)
  for (jExp in 1:length(uniqExp)) {
    jExpNdx <- dfSub$MaxScaledExp == uniqExp[jExp]
    dfSubSub <- dfSub[jExpNdx,]
    dfSubSub <- dfSubSub[order(dfSubSub$PctExpressed, decreasing = TRUE),]
    if (jExp == 1 & j == 1) {
      dfAll <- dfSubSub
    } else {
      dfAll <- rbind(dfAll, dfSubSub)
    }
  }
}
geneListSupervised <- rev(dfAll$Gene)

# Plot DotPlot with supervised genes
pdf(file = 'scRNA_COVID_Healthy_top50supervised_nephronorder_4.pdf', 7, 19)
DotPlot(subset(scData, subset = diseasetype == "LivingDonor"), features = geneListSupervised, group.by = 'subclass.l1', scale.by = 'size') + RotatedAxis() + coord_flip() + scale_colour_gradient2(low = "#5A5AA0", mid = "white", high = "#FFB600")
dev.off()

# Plot DotPlots with specific gene lists
pdf(file = 'scRNA_COVID_Healthy_Migration_supervised_nephronorder_4.pdf', 7, 18)
DotPlot(subset(scData, subset = diseasetype == "LivingDonor"), features = geneListSupervised, group.by = 'subclass.l1', scale.by = 'size') + RotatedAxis() + coord_flip() + scale_colour_gradient2(low = "#5A5AA0", mid = "white", high = "#FFB600")
dev.off()

pdf(file = 'scRNA_COVID_Healthy_Oxygen_supervised_nephronorder_4.pdf', 7, 7)
DotPlot(subset(scData, subset = diseasetype == "LivingDonor"), features = geneListSupervised, group.by = 'subclass.l1', scale.by = 'size') + RotatedAxis() + coord_flip() + scale_colour_gradient2(low = "#5A5AA0", mid = "white", high = "#FFB600")
dev.off()

pdf(file = 'scRNA_COVID_Healthy_CellAdhesion_supervised_nephronorder_4.pdf', 7, 7.5)
DotPlot(subset(scData, subset = diseasetype == "LivingDonor"), features = geneListSupervised, group.by = 'subclass.l1', scale.by = 'size') + RotatedAxis() + coord_flip() + scale_colour_gradient2(low = "#5A5AA0", mid = "white", high = "#FFB600")
dev.off()
