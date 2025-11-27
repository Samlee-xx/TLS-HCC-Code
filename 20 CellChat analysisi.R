
###### CellChat analysis ######
library(Seurat)
library(dplyr)
library(CellChat)
library(patchwork)

chol = readRDS('...\\Bile Duct.rds')
neutro = readRDS('...\\01 Neu seurat.rds')

DefaultAssay(chol)   <- "RNA"
DefaultAssay(neutro) <- "RNA"
chol_c0    <- subset(chol, idents = 0)   
neutro_all <- neutro                  
chol_c0$cell_group    <- "Chol_C0"
neutro_all$cell_group <- "Neutro"

obj <- merge(chol_c0, neutro_all)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 2000)
obj <- ScaleData(obj, features = VariableFeatures(obj))

Idents(obj) <- "cell_group"
expr <- GetAssayData(obj, assay = "RNA", slot = "data")  
meta <- data.frame(group = Idents(obj))
rownames(meta) <- colnames(obj)
cellchat <- createCellChat(object = expr, meta = meta, group.by = "group")

CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat)                 
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)         
cellchat <- filterCommunication(cellchat, min.cells = 10) 
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)               

setwd('...\\12-Cellchat neu chol')
saveRDS(cellchat, "CellChat_cholC0_vs_neutro.rds")

cellchat = readRDS("CellChat_cholC0_vs_neutro.rds")
unique(cellchat@netP$pathway.name)
netVisual_bubble(cellchat, sources.use = "Chol_C0", targets.use = "Neutro", remove.isolate = TRUE)
head(cellchat@DB$interaction$pathway_name)

netVisual_bubble(
  cellchat,
  sources.use = "Chol_C0",
  targets.use = "Neutro",
  signaling = "CXCL",
  remove.isolate = TRUE
)

pdf("bubble_CholC0_to_Neutro.pdf", width = 3, height = 2)
netVisual_bubble(
  cellchat,
  sources.use = "Chol_C0",
  targets.use = "Neutro",
  signaling = "CXCL",
  remove.isolate = TRUE
)
dev.off()
