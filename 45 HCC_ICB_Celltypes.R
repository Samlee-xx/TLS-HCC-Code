
library(dplyr)
library(reshape2)
library(stringr)
library(sva)
library(Seurat)
setwd('')
data <- read.csv('...\\cell_table_arcsinh_transformed.csv')

data.sub <- filter(data, cell_size <=200)
data.sub <- filter(data, cell_size >=10)

data.sub$sample = str_extract(data.sub$fov,'^ICI_\\d{2}')
data.sub$cellid = paste0(data.sub$fov,"_",data.sub$label)
data.sub$batch = str_extract(data.sub$fov, "^ICI_\\d{2}_ROI\\d+")
table(data.sub$batch)

colnames(data.sub)
count.data = data.sub[,c(2:24,27:30,32,36:46)]
count.data = t(count.data)
combat_data <- ComBat(count.data, batch = data.sub$batch)
colnames(combat_data) =data.sub$cellid

data.obj = CreateSeuratObject(counts =  combat_data)
colnames(data.obj)
rownames(data.obj)
celltype.marker = c('Ecadherin','PanKeratin', #tumor
                    'LYVE1','CD34',# Endothelial
                    'SMA','Collagen' ,# Fibroblast
                    'CD45','CD45RO','CD20','CD3','CD4','CD8a','CD56','CD15','CD68' # Immune marker
)


VariableFeatures(data.obj) <- rownames(data.obj)
data.obj = ScaleData(data.obj)
data.obj@assays$RNA@scale.data[1:15,1:5]
data.obj <- RunPCA(object = data.obj, npcs = 10, verbose = TRUE)
data.obj <- RunUMAP(object = data.obj, dims = 1:10, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
data.obj <- FindNeighbors(object = data.obj, dims = 1:10, verbose = TRUE)
data.obj <- FindClusters(object = data.obj, verbose = TRUE, resolution = 0.5)
DimPlot(data.obj)

data.obj.sub = subset(data.obj,seurat_clusters %in% c(0:25)) #remove cluster less 1000
dim(data.obj.sub)
data.obj.sub = ScaleData(data.obj.sub)
data.obj.sub <- RunPCA(object = data.obj.sub, npcs = 10, verbose = TRUE)
data.obj.sub <- RunUMAP(object = data.obj.sub, dims = 1:10, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
data.obj.sub <- FindNeighbors(object = data.obj.sub, dims = 1:10, verbose = TRUE)
data.obj.sub <- FindClusters(object = data.obj.sub, verbose = TRUE, resolution = 0.5)
DimPlot(data.obj.sub)
table(data.obj.sub$seurat_clusters)

a = DoHeatmap(subset(data.obj.sub, downsample=2000), features = celltype.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '01 DoHeatmap marker.pdf',height = 7,width = 20)
print(a)
dev.off()

immune.cluster = subset(data.obj.sub,seurat_clusters %in% c(0,2,4,6,15,18))
stroma.cluster = subset(data.obj.sub,seurat_clusters %in% c(1,10,11,17,19,23,25))
epith.cluster = subset(data.obj.sub,seurat_clusters %in% c(3,7,9,12,13,14,16,20,21,24,26,27,28))
others.cluster = subset(data.obj.sub,seurat_clusters %in% c(5,8,22,29,30))
dim(others.cluster)
table(others.cluster$seurat_clusters)
## others cluster re-defined

VariableFeatures(others.cluster) <- rownames(others.cluster)
others.cluster = ScaleData(others.cluster)
#data.obj@assays$RNA@scale.data[1:15,1:5]
others.cluster <- RunPCA(object = others.cluster, npcs = 10, verbose = TRUE)
others.cluster <- RunUMAP(object = others.cluster, dims = 1:10, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
others.cluster <- FindNeighbors(object = others.cluster, dims = 1:10, verbose = TRUE)
others.cluster <- FindClusters(object = others.cluster, verbose = TRUE, resolution = 0.5)
DimPlot(others.cluster)

a = DoHeatmap(subset(others.cluster, downsample=2000), features = celltype.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '02 DoHeatmap others marker.pdf',height = 7,width = 20)
print(a)
dev.off()

other.immune = subset(others.cluster,seurat_clusters %in% c(0,3,11,12))
other.stroma = subset(others.cluster,seurat_clusters %in% c(7,9))
other.epith = subset(others.cluster,seurat_clusters %in% c(5,6,10))
other.other = subset(others.cluster,seurat_clusters %in% c(1,2,4,8,13,14))

## immune total defined celltype 
immune.total = data.obj.sub[,c(colnames(immune.cluster),colnames(other.immune))]
immune.total = ScaleData(immune.total)
immune.total <- RunPCA(object = immune.total, npcs = 10, verbose = TRUE)
immune.total <- RunUMAP(object = immune.total, dims = 1:10, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
immune.total <- FindNeighbors(object = immune.total, dims = 1:10, verbose = TRUE)
immune.total <- FindClusters(object = immune.total, verbose = TRUE, resolution = 0.5)
DimPlot(immune.total)
rownames(immune.total)

immune.marker = c('CD3','CD4','CD8a','CD20','CD56','CD68','CD15','CD163','HLA.DR')
a = DoHeatmap(subset(immune.total, downsample=2000), features = immune.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '03 DoHeatmap IMMUNE marker.pdf',height = 7,width = 20)
print(a)
dev.off()

new.cluster.ids <- c('Macrophage','Macrophage','CD4T','Neutrophil','CD4T','CD8T','Macrophage',
                     'Neutrophil','immune_others','B','CD8T','Macrophage','Macrophage','CD4T','Neutrophil','immune_others'
                     
)

names(new.cluster.ids) <- levels(immune.total)
immune.total<- RenameIdents(immune.total, new.cluster.ids)
immune.total$celltype = immune.total@active.ident
DimPlot(immune.total,group.by = 'celltype')
saveRDS(immune.total,'01 Immune.rds')

### Stroma cell
stroma.total = data.obj.sub[,c(colnames(stroma.cluster),colnames(other.stroma))]
stroma.total = ScaleData(stroma.total)
stroma.total <- RunPCA(object = stroma.total, npcs = 10, verbose = TRUE)
stroma.total <- RunUMAP(object = stroma.total, dims = 1:10, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
stroma.total <- FindNeighbors(object = stroma.total, dims = 1:10, verbose = TRUE)
stroma.total <- FindClusters(object = stroma.total, verbose = TRUE, resolution = 0.5)
DimPlot(stroma.total)
table(stroma.total$seurat_clusters)

rownames(stroma.total)
stroma.marker = c('CD34','Collagen','LYVE1','SMA')

a = DoHeatmap(subset(stroma.total, downsample=2000), features = stroma.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '04 DoHeatmap Stroma marker.pdf',height = 7,width = 20)
print(a)
dev.off()

stroma.total$celltype = case_when(stroma.total$seurat_clusters %in% c(0,1,2,6,7,12,13,15,18)~'Endothelial',
                                  stroma.total$seurat_clusters %in% c(3,4,5,8,9,11,14,16,17,19,20,21,22,
                                                                      23,24,25,26,27,28,29,30,31,32,33,34,
                                                                      35,36,37,38,39)~'Fibroblast',
                                  stroma.total$seurat_clusters %in% c(10)~'Lymphatic')
DimPlot(stroma.total,group.by = 'celltype')
saveRDS(stroma.total,'02 Stroma.rds')
stroma.total = readRDS('02 Stroma.rds')

### Epthelial 
epith.total = data.obj.sub[,c(colnames(epith.cluster),colnames(other.epith))]
epith.total = ScaleData(epith.total)
epith.total <- RunPCA(object = epith.total, npcs = 10, verbose = TRUE)
epith.total <- RunUMAP(object = epith.total, dims = 1:10, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
epith.total <- FindNeighbors(object = epith.total, dims = 1:10, verbose = TRUE)
epith.total <- FindClusters(object = epith.total, verbose = TRUE, resolution = 0.5)
DimPlot(epith.total)
table(epith.total$seurat_clusters)
rownames(epith.total)
epith.marker = c('PanKeratin','Ecadherin','CD15')

a = DoHeatmap(subset(epith.total, downsample=2000), features = epith.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '05 DoHeatmap Epith marker.pdf',height = 7,width = 20)
print(a)
dev.off()

DimPlot(epith.total,label = TRUE)

epith.total$celltype = ifelse(epith.total$seurat_clusters == 12,'Bile Duct','Heptocyte')
DimPlot(epith.total,group.by = 'celltype')
saveRDS(epith.total,'03 Epithelial.rds')
epith.total = readRDS('03 Epithelial.rds')

table(epith.total$celltype)
table(immune.total$celltype)
table(stroma.total$celltype)
phenotype.df = data.frame(cellid = c(colnames(epith.total),
                                     colnames(stroma.total),
                                     colnames(immune.total)),
                          celltype = c(epith.total$celltype,
                                       stroma.total$celltype,
                                       immune.total$celltype))
table(phenotype.df$celltype)

data.sub = data.sub %>% left_join(phenotype.df)
colnames(data.sub)
metadata.total = data.sub[,c(47:67)]
table(metadata.total$celltype)
metadata.total$celltype = case_when(metadata.total$celltype == 1 ~ 'Macrophage',
                                    metadata.total$celltype == 2 ~ 'CD4T',
                                    metadata.total$celltype == 3 ~ 'Neutrophil',
                                    metadata.total$celltype == 4 ~ 'CD8T',
                                    metadata.total$celltype == 5 ~ 'Immune_others',
                                    metadata.total$celltype == 6 ~ 'B',
                                    metadata.total$celltype == 'Bile Duct' ~ 'Bile Duct',
                                    metadata.total$celltype == 'Endothelial' ~ 'Endothelial',
                                    metadata.total$celltype == 'Fibroblast' ~ 'Fibroblast',
                                    metadata.total$celltype == 'Heptocyte' ~ 'Heptocyte',
                                    metadata.total$celltype == 'Lymphatic' ~ 'Lymphatic')
metadata.total[is.na(metadata.total$celltype),'celltype'] = 'Unassign'

write.csv(metadata.total,'...\\Celltype total metadata.csv')
