
#! R/4.0.3

Args <- commandArgs(trailingOnly = T)

library(CellChat)

ExprF <- Args[1]
AnnF <- Args[2]
DBS <- Args[3]

#Expr <- read.delim(file = ExprF, row.names = 1, as.is = T, check.names = F)
Expr <- readRDS(file = ExprF)
Ann <- read.delim(file = AnnF, row.names = 1, as.is = T)

cellchat <- createCellChat(object = as.matrix(Expr), meta = Ann, group.by = "cell_type")

levels(cellchat@idents)

groupSize <- as.numeric(table(cellchat@idents))


#CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

if(DBS == "mouse"){

	CellChatDB <- CellChatDB.mouse

}else if(DBS == "homo"){

	CellChatDB <- CellChatDB.human

}else{

	CellChatDB <- readRDS(file = DBS)

}

dplyr::glimpse(CellChatDB$interaction)

#pdf("CellChatDB.pdf")
#
#	showDatabaseCategory(CellChatDB)
#
#dev.off()


CellChatDB.use <- CellChatDB

cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

saveRDS(cellchat, file = "cellchat.rds")

future::plan("multiprocess", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)


saveRDS(cellchat, file = "cellchat.rds")


pdf("Overview.pdf")

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

dev.off()


pdf("OverviewCount.pdf")

mat <- cellchat@net$count

par(mfrow = c(3,4), xpd=TRUE)

for (i in 1:nrow(mat)) {

  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))

  mat2[i, ] <- mat[i, ]

  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], vertex.label.cex = 0.3)

}

dev.off()


pdf("OverviewWeight.pdf")

mat <- cellchat@net$weight

par(mfrow = c(3,4), xpd=TRUE)

for (i in 1:nrow(mat)) {

  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))

  mat2[i, ] <- mat[i, ]

  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], vertex.label.cex = 0.3)

}

dev.off()

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
#netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

pdf(paste0("All", ".pca.pdf"))
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
#gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
#gg1 + gg2
gg1

dev.off()


pdf(paste0("All", ".pathway.pdf"), width = 14)

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 4)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 4)
ht1 + ht2

write.table(ht1@matrix, file = "outgoing.xls", sep = "\t", quote = F, col.names = NA)

write.table(ht2@matrix, file = "incoming.xls", sep = "\t", quote = F, col.names = NA)

dev.off()


saveRDS(cellchat, file = "cellchat.rds")

Resu <- subsetCommunication(cellchat)

write.table(Resu, file = "TalkTalk.xls", sep = "\t", quote = F, row.names = F)

## Access all the signaling pathways showing significant communications
#pathways.show.all <- cellchat@netP$pathways
#
## check the order of cell identity to set suitable vertex.receiver
#levels(cellchat@idents)
#
#vertex.receiver = seq(1,4)
#
#for (i in 1:length(pathways.show.all)) {
#
#  # Visualize communication network associated with both signaling pathway and individual L-R pairs
#
#  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
#
#  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
#
#  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
#
#  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
#}
