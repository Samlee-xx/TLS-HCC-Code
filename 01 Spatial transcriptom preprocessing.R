
###### code for spatial transcriptom analysis ######
###01 merge tls data ----
setwd('...\\Spata_segmentation')
PATH='...\\Spata_segmentation'
slice.name=list.files()
for (i in 1:length(slice.name)){
  sec.path=paste0(PATH,'\\',slice.name[i])
  slice.file.name=list.files(path=sec.path)
  TLS.name=grep('TLS_',slice.file.name,value=TRUE)
  TLS.head.name=gsub('.rds','',TLS.name)
  for (g in TLS.name){
    TLS.head.name=gsub('.rds','',g)
    assign(TLS.head.name,readRDS(paste0(sec.path,'\\',g)))
    Seurat_obj=get(TLS.head.name)
    assign(paste0(TLS.head.name,'_matrix'),Seurat_obj@assays[["RNA"]]@counts)
    matrix=get(paste0(TLS.head.name,'_matrix'))
    assign(paste0(TLS.head.name,'_rowsum'),as.data.frame(rowSums(matrix)))
    tem=get(paste0(TLS.head.name,'_rowsum'))
    colnames(tem)=paste0(slice.name[i],'_',TLS.head.name)
    assign(paste0(TLS.head.name,'_final'),tem)
  }
  combine_file=c()
  for (m in 1:length(TLS.name)) {
    combine_file[m]=paste0('TLS_',m,'_final')
  }
  if (length(combine_file)>=3){
    combined_final=cbind(TLS_1_final,TLS_2_final)
    for (f in combine_file[3:length(combine_file)]){
      combined_tem=get(f)
      combined_final=cbind(combined_final,combined_tem)
    }
    assign(paste0(slice.name[i],'_final'),combined_final)
  } else if (length(combine_file)==2){
    combined_final=cbind(TLS_1_final,TLS_2_final)
    assign(paste0(slice.name[i],'_final'),combined_final)
  } else if (length(combine_file)==1){
    assign(paste0(slice.name[i],'_final'),TLS_1_final)
  }
}

total.slice.file.name=c()
for (i in 1:length(slice.name)){
  total.slice.file.name[i]=paste0(slice.name[i],'_final')
}

library(dplyr)
`084717_A_final`$group=rownames(`084717_A_final`)
`084717_B_final`$group=rownames(`084717_B_final`)
total.count=`084717_A_final` %>% full_join(`084717_B_final`,by='group')


for (i in 1:length(slice.name)){
  tem=total.slice.file.name[i]
  tem=get(tem)
  tem$group=rownames(tem)
  assign(total.slice.file.name[i],tem)
}

for (i in total.slice.file.name[3:length(total.slice.file.name)]){
  tem=get(i)
  total.count=total.count %>% full_join(tem,by='group')
}

rownames(total.count)=total.count$group
total.count=total.count[,-16]
total.count[is.na(total.count)]=0
saveRDS(total.count,'...\\Total_tls_count.rds')

### 02 Remove batch ----
library(sva)
library(stringr)
total.count=readRDS('B:\\ST\\Analysis\\01_Total_TLS\\Total_tls_count_20220814.rds')
batch.df=as.data.frame(colnames(total.count))
for (i in 1:dim(batch.df)[1]){
  a=batch.df$`colnames(total.count)`[i]
  batch.df$patient[i]=str_sub(a,1,6)
}
table(batch.df$patient)

total.count=total.count[,1:166]
batch.level=c(rep(1,64),rep(2,15),rep(3,9),rep(4,19),rep(5,32),rep(6,16),rep(7,11))
adjusted=ComBat_seq(as.matrix(total.count),
                    batch = batch.level)
saveRDS(adjusted,'...\\Batch_remove_tls_count_sample166.rds')

### 03 TPM calculate and clustering ----
library(dplyr)
library(IOBR)
library(Seurat)
gene_length=read.csv('...\\gene_length.csv')
TLS.count=readRDS('...\\Batch_remove_tls_count_sample166.rds')

gene_length=gene_length %>% distinct(gene_name,.keep_all = T)
rownames(gene_length)=gene_length$gene_name

TLS.rowname=rownames(TLS.count)
total.genename=gene_length$gene_name
commom_name=intersect(TLS.rowname,total.genename)
TLS.count.fliter=TLS.count[commom_name,]
gene_length.fliter=gene_length[commom_name,]

TPM=count2tpm(as.matrix(TLS.count.fliter),effLength = gene_length.fliter,
              id='gene_name',
              length = 'gene_length',
              gene_symbol = 'gene_name')

TPM=TPM[,-139]
colnames(TPM)[117]='226766_C_L_TLS_1'

tls<- CreateSeuratObject(counts = TPM, project = "TLS", min.cells = 0, min.features = 0)
tls@assays[["RNA"]]@scale.data=TPM
all.genes <- rownames(tls)
tls <- RunPCA(tls, features = all.genes)
tls<- FindNeighbors(tls, dims = 1:5)
tls <- FindClusters(tls, resolution = 0.5)
tls <- RunUMAP(tls, dims = 1:15)
DimPlot(tls, reduction = "umap")

tls <- NormalizeData(tls, normalization.method = "LogNormalize", scale.factor = 10000)
tls <- FindVariableFeatures(tls, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tls)
tls <- ScaleData(tls, features = all.genes)
marker=FindAllMarkers(tls,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

