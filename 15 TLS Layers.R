###### TLS Layers ######
library(dplyr)
library(reshape2)

spatil.layer = read.csv('...\\TLS Spatial Layer labels.csv')
total.meta = readRDS('...\\TA_TLS_Metadata.rds')

spatil.layer <- spatil.layer[!duplicated(spatil.layer$CellID), ]
rownames(total.meta) = total.meta$cellid
total.meta = total.meta[spatil.layer$CellID,]

spatil.layer$Level = as.numeric(spatil.layer$Level)
spatil.layer$tls.id = total.meta$tls_id
tls.id = unique(spatil.layer$tls.id)

spatil.layer$Boundary_Class = spatil.layer$Level

for (i in tls.id) {
  sub.delanuay = subset(spatil.layer, tls.id == i)
  if (max(sub.delanuay$Level) <= 3) {
    sub.delanuay$Boundary_Class = case_when(
      sub.delanuay$Level == 1 ~ 'Edge',
      sub.delanuay$Level == 2 ~ 'Transition',
      sub.delanuay$Level == 3 ~ 'Core'
    )
  } else {
    each_slice = max(sub.delanuay$Level) %/% 3
    sub.delanuay$Boundary_Class = case_when(
      sub.delanuay$Level %in% c(1:each_slice) ~ 'Edge',
      sub.delanuay$Level %in% c((1 + each_slice):(2 * each_slice)) ~ 'Transition',
      sub.delanuay$Level %in% c((1 + 2 * each_slice):(3 * each_slice+3)) ~ 'Core'
    )
  }
  spatil.layer[spatil.layer$tls.id == i, ] = sub.delanuay
}


spatil.layer$celltype = total.meta$celltype
spatil.layer$subtype  = total.meta$subtype

celltype.table.meta = readRDS('...\\TA TLS Clustering.rds')
t.enriched.tls = subset(celltype.table.meta,group == 'T-enriched')
t.enriched.tls = t.enriched.tls$tls.id

b.enriched.tls = subset(celltype.table.meta,group == 'B-enriched')
b.enriched.tls = b.enriched.tls$tls.id


spatil.layer.sub = subset(spatil.layer,tls.id %in% t.enriched.tls)
class.boundary = dcast(spatil.layer.sub,celltype~Boundary_Class)
rownames(class.boundary) = class.boundary$celltype
class.boundary = class.boundary[,-1]
class.boundary = class.boundary[,c('Edge','Transition','Core')]
rownames(class.boundary)


a = class.boundary
p.value = class.boundary
roe = class.boundary
for (i in 1:dim(a)[1]){
  for (j in 1:dim(a)[2]){
    onecelltype_onetissue=a[i,j]
    onetissue_othercelltype=sum(a[i,])-onecelltype_onetissue
    onecelltype_othertissue=sum(a[,j])-onecelltype_onetissue
    othercelltype_ohertissue=sum(a)-onecelltype_onetissue - onetissue_othercelltype - onecelltype_othertissue
    x=matrix(c(onecelltype_onetissue,onecelltype_othertissue,onetissue_othercelltype,othercelltype_ohertissue),
             nrow = 2,ncol = 2)
    b=chisq.test(x)$expected[1,1]
    m=onecelltype_onetissue/b
    roe[i,j]=m
    p.value[i,j] = chisq.test(x)$p.value
  }
} 
colnames(roe)

bk <- c(seq(0.5,1,by=0.01),seq(1.001,1.6,by=0.01))
library(pheatmap)
rownames(roe)
roe =roe[c('B','T-B','CD4T','CD8T','DC','Neutrophil','Macrophage',
           'Endothelial','HEV','Lymphatic','Fibroblast','Bile_Duct','Unassign'),]
a = pheatmap(roe,cluster_rows = FALSE,cluster_cols = FALSE,
             #gaps_col = c(5),
             #gaps_row = c(9,17),
             cellwidth = 25, cellheight = 20,
             color = c(colorRampPalette(colors = c("white","#F7E6D7"))(length(bk)/2),
                       colorRampPalette(colors = c("#F7E6D7","#FA151D"))(length(bk)/2)),
             #legend_breaks=seq(-1,0,1),
             breaks=bk)

pdf(file = 'b enriched tls celltype.pdf',height = 7,width = 4)
print(a)
dev.off()
