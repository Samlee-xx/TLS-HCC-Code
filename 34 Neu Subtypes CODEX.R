
setwd('...\\03 Neu')

total.meta <- readRDS('...\\01 Total Cell Metadata.rds')
data <- readRDS('...\\Batch Removed data.rds')

neu.meta <- subset(total.meta,celltype == 'Neutrophil')
neu.data <- data[,neu.meta$cellid]
neu.cluster <- CreateSeuratObject(counts = neu.data)
neutro.marker <- c('CD14','CD66b','MPO','S100A8-9','S100A12','HLA.DR','CD74',
                  'MIF','ISG15','PD.L1','IDO1','iNOS','MMP9','APOC3','COX2',
                  'CXCL8','NE','IFNG','IL1B','CD38','Ki67')

VariableFeatures(neu.cluster) <- neutro.marker
neu.cluster <- ScaleData(neu.cluster)
neu.cluster <- RunPCA(object = neu.cluster, npcs =10, verbose = TRUE)
neu.cluster <- RunUMAP(object = neu.cluster, dims = 1:10, verbose = TRUE)
neu.cluster <- FindNeighbors(object = neu.cluster, dims = 1:10, verbose = TRUE)
neu.cluster <- FindClusters(object = neu.cluster, verbose = TRUE, resolution =0.8)
DimPlot(neu.cluster)

pheatmap(cor(t(neu.cluster@assays$RNA@scale.data)))
a = DoHeatmap(neu.cluster, features = neutro.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '01 DoHeatmap Neutro.pdf',height = 5,width = 10)
print(a)
dev.off()

new.cluster.ids <- c('Neu_C0_S100A8',
                     'Neu_C0_S100A8',
                     'Neu_C4_S100A12',
                     'Neu_C2_MIF',
                     'Neu_C2_MIF',
                     'Neu_Others',
                     'Neu_C5_APC',
                     'Neu_Others',
                     'Neu_C1_APOC3',
                     'Neu_C5_APC',
                     'Neu_Others',
                     'Neu_C4_S100A12', 
                     'Neu_C3_PD.L1',
                     'Neu_Others' ,
                     'Neu_Others',
                     'Neu_C0_S100A8',
                     'Neu_C2_MIF', 
                     'Neu_Others',
                     'Neu_Others'
)


names(new.cluster.ids) <- levels(neu.cluster)
neu.cluster<- RenameIdents(neu.cluster, new.cluster.ids)
neu.cluster$celltype <- neu.cluster@active.ident
Idents(neu.cluster) <- neu.cluster$seurat_clusters
a <- DimPlot(neu.cluster,split.by   = 'celltype')

a <- DoHeatmap(subset(neu.cluster, downsample=2000), features = neutro.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '02 DoHeatmap Neu subtypes marker.pdf',height = 10,width = 20)
print(a)
dev.off()

saveRDS(neu.cluster,'01 Neu.RDS')
neu.cluster <- readRDS('01 Neu.RDS')


