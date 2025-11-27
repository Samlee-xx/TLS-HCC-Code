

setwd('...\\08 Data Merge')
base.path = '...\\03 Segmentation\\'
slide.name = list.files('...\\03 Segmentation')[1:9]

PATH1='01 124989\\cell_table\\cell_table_arcsinh_transformed.csv'
codex.pre1 = read.csv(paste0(base.path,PATH1))
cellid1 = paste0(codex.pre1$fov,'_',codex.pre1$label)
rownames(codex.pre1)=cellid1
codex.pre1$slide = 'Sample_124989'
table(codex.pre1$fov)

PATH2='02 86130\\cell_table\\cell_table_arcsinh_transformed.csv'
codex.pre2 = read.csv(paste0(base.path,PATH2))
cellid2 = paste0(codex.pre2$fov,'_',codex.pre2$label)
rownames(codex.pre2)=cellid2
codex.pre2$slide = 'Sample_86130'
table(codex.pre2$fov)

PATH3='03 136759\\cell_table\\cell_table_arcsinh_transformed.csv'
codex.pre3 = read.csv(paste0(base.path,PATH3))
cellid3 = paste0(codex.pre3$fov,'_',codex.pre3$label)
rownames(codex.pre3)=cellid3
codex.pre3$slide = 'sample_136759'
table(codex.pre3$fov)

PATH4='04 134205\\cell_table\\cell_table_arcsinh_transformed.csv'
codex.pre4 = read.csv(paste0(base.path,PATH4))
cellid4 = paste0(codex.pre4$fov,'_',codex.pre4$label)
rownames(codex.pre4)=cellid4
codex.pre4$slide = 'Sample_134205'
table(codex.pre4$fov)

PATH5='05 133698\\cell_table\\cell_table_arcsinh_transformed.csv'
codex.pre5 = read.csv(paste0(base.path,PATH5))
cellid5 = paste0(codex.pre5$fov,'_',codex.pre5$label)
rownames(codex.pre5)=cellid5
codex.pre5$slide = 'Sample_133698'
table(codex.pre5$fov)

PATH6='06 123942\\cell_table\\cell_table_arcsinh_transformed.csv'
codex.pre6 = read.csv(paste0(base.path,PATH6))
cellid6 = paste0(codex.pre6$fov,'_',codex.pre6$label)
rownames(codex.pre6)=cellid6
codex.pre6$slide = 'Sample_123942'
table(codex.pre6$fov)

PATH7='07 136254\\cell_table\\cell_table_arcsinh_transformed.csv'
codex.pre7 = read.csv(paste0(base.path,PATH7))
cellid7 = paste0(codex.pre7$fov,'_',codex.pre7$label)
rownames(codex.pre7)=cellid7
codex.pre7$slide = 'Sample_136254'
table(codex.pre7$fov)

PATH8='08 134024\\cell_table\\cell_table_arcsinh_transformed.csv'
codex.pre8 = read.csv(paste0(base.path,PATH8))
cellid8 = paste0(codex.pre8$fov,'_',codex.pre8$label)
rownames(codex.pre8)=cellid8
codex.pre8$slide = 'Sample_134024'
table(codex.pre8$fov)

PATH9='09 176552\\cell_table\\cell_table_arcsinh_transformed.csv'
codex.pre9 = read.csv(paste0(base.path,PATH9))
cellid9 = paste0(codex.pre9$fov,'_',codex.pre9$label)
rownames(codex.pre9)=cellid9
codex.pre9$slide = 'Sample_176552'
table(codex.pre9$fov)

codex.total = rbind(codex.pre1,codex.pre2,codex.pre3,codex.pre4,
                    codex.pre5,codex.pre6,codex.pre7,codex.pre8,codex.pre9)

all.feature = colnames(codex.total)
select.feature = all.feature[c(1:52,60,61,68,69)] 
feature.filtered = codex.total[,select.feature]

write.csv(feature.filtered,'Merge_total.csv')

feature.filtered = read.csv('...\\Merge_total.csv')
colnames(feature.filtered)

size.prevent = filter(feature.filtered, cell_size < 800 , cell_size > 30) 
colnames(size.filtered)[1] = 'cellid'
metadata = size.filtered[,c(1,2,53:57)] 
saveRDS(metadata,'Filtered size Metadata.rds')


size.prevent = size.prevent[,c(3:52)]
colnames(size.prevent)
rownames(size.prevent) = metadata$cellid
size.prevent = t(size.prevent)
combat_data <- ComBat(size.prevent, batch = metadata$slide)
saveRDS(combat_data,'Batch Removed data.rds')
combat_data = readRDS('Batch Removed data.rds')


rownames(combat_data)
combat_data =combat_data[-25,]
tls = CreateSeuratObject(counts =  combat_data)
colnames(tls)
rownames(tls)
celltype.marker = c('APOC3','COX2','Pan.Cytokeratin', #tumor
                    'MECA79','CD31', 'Podoplanin',# Endothelial
                    'aSMA', # Fibroblast
                    'CD45','CD45RO', # Immune marker
                    'CD14','CD66b','CD68','CD11c','MPO', # Myeloid
                    'CD20', 'CD21','CD79a', # B cell
                    'CD3e','CD56','CD4','CD8'  # NKT 
)

VariableFeatures(tls) <- celltype.marker
tls = ScaleData(tls)
tls@assays$RNA@scale.data[1:21,1:5]
tls <- RunPCA(object = tls, npcs = 20, verbose = TRUE)
tls <- RunUMAP(object = tls, dims = 1:20, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
tls <- FindNeighbors(object = tls, dims = 1:20, verbose = TRUE)
tls <- FindClusters(object = tls, verbose = TRUE, resolution = 0.5)
DimPlot(tls)
table(tls$seurat_clusters)
saveRDS(tls,'Seurat TLS.rds')  


