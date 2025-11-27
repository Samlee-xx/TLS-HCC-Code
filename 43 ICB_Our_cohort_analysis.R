
###### Task 1 ######

setwd('...\\task1\\04-TLS Dbscan')
dbscan.out = read.csv('DBSCAN2.csv')
dbscan.out = subset(dbscan.out,!(cluster_label == '-1'))
dbscan.out$tls.label = paste0(dbscan.out$fov,'_',dbscan.out$cluster_label)
unique(dbscan.out$tls.label)
dbscan.out = subset(dbscan.out,!(tls.label %in% c('TLS-B3_0')))
unique(dbscan.out$tls.label)

celltype_count <- dbscan.out %>%
  group_by(tls.label, celltype) %>%
  summarise(count = n(), .groups = 'drop')
tls.cell.count = celltype_count %>% group_by(tls.label) %>% summarise(sum(count))
colnames(tls.cell.count)[2] ='CellNum'
celltype_count.filter <- celltype_count %>%
  group_by(tls.label) %>%
  mutate(percent = count / sum(count)*100)

wide_percent <- dcast(celltype_count.filter, tls.label ~ celltype, value.var = "percent", fill = 0)
wide_percent = wide_percent %>% as.data.frame() 
rownames(wide_percent) = wide_percent$tls.label
wide_percent = wide_percent[,-1]
colnames(wide_percent)

wide_percent$t.cell = wide_percent$CD4T+wide_percent$CD8T
colnames(wide_percent)
wide_percent = wide_percent[,-c(3,4)]
wide_percent$Endothelial = wide_percent$Endothelial+wide_percent$lymphatic
wide_percent = wide_percent[,-c(5)]
pheatmap::pheatmap(wide_percent)

library(pheatmap)
a <- pheatmap(as.matrix(wide_percent))
library(pheatmap)
row_hc <- a$tree_row
plot(row_hc, main = "Row Clustering Dendrogram from pheatmap")
groups <- cutree(row_hc, k = 2)
wide_percent$group <- groups[rownames(wide_percent)]
table(wide_percent$group)
wide_percent$group = ifelse(wide_percent$group %in% c(2),'B-TLS','T-TLS')
table(wide_percent$group)

wide_percent$tb = log((wide_percent$B+1)/(wide_percent$t.cell+1))
library(ggpubr)
colnames(wide_percent)
a = ggplot(wide_percent,aes(x=group,y=bile.duct,fill=group))+geom_boxplot(outlier.shape = NA,width=0.5)+
  geom_jitter(width =0.3,size=1,color = "grey50", alpha = 0.6)+ 
  theme(panel.background = element_blank(),axis.line = element_line())+ 
  stat_compare_means(label="p.format")+    
  stat_boxplot(geom = "errorbar",width=0.5)+   
  theme(legend.position = "bottom") +
  ylab('Percentage')+
  scale_fill_manual(values=c('#c095e4','#72b043'))

###### Task 2 ######
library(RedeR)
library(RColorBrewer)
library(igraph)
library(TreeAndLeaf)
library(pheatmap)
load('...\\Task02.rds')
task2.ct = subset(wide_percent.filter,B_cell+CD4T+CD8T+T_cell>=30)
task2.ct$tcell = task2.ct$CD4T+task2.ct$CD8T+task2.ct$T_cell
colnames(task2.ct)

task2.ct = task2.ct[,-c(3,4,8)]
task2.ct = task2.ct[-match(c('M6_cluster4','O6_cluster1','O6_cluster7',
                             'Q6_cluster1','Q6_cluster3','Q6_cluster4'), rownames(task2.ct), nomatch = 0), , drop=FALSE]
task2.ct$fov =  sub("_.*", "", rownames(task2.ct))
task2.ct = subset(task2.ct,fov %in% c('C4','D4','E7','L5','M5','N5','O5',
                                      'G1','G8','I8','K8','E3','J2','K2',
                                      'F5','G5','K7','Q8','H4','J4','K4',
                                      'C2','L4','M4','O4','I1','J1','K1','M1'))
task2.ct = task2.ct[,-8]

tls.cell.count = tls.cell.count %>% as.data.frame()
rownames(tls.cell.count) = tls.cell.count$tls.label
tls.cell.count = tls.cell.count[rownames(task2.ct),]

a <- pheatmap(as.matrix(task2.ct), silent = TRUE)
library(pheatmap)
row_hc <- a$tree_row
plot(row_hc, main = "Row Clustering Dendrogram from pheatmap")
groups <- cutree(row_hc, k = 5)
task2.ct$group <- groups[rownames(task2.ct)]
table(task2.ct$group)
task2.ct$group = ifelse(task2.ct$group %in% c(4,5),'B-TLS','T-TLS')
table(task2.ct$group)

colnames(task2.ct)
task2.ct$tb = log((task2.ct$B_cell+1)/(task2.ct$tcell+1))
tls.cell.count$tb = task2.ct$tb
tls.cell.count$b = task2.ct$B_cell

colnames(task2.ct)
a = ggplot(task2.ct,aes(x=group,y=tcell,fill=group))+geom_boxplot(outlier.shape = NA,width=0.5)+
  geom_jitter(width =0.3,size=1,color = "grey50", alpha = 0.6)+ 
  theme(panel.background = element_blank(),axis.line = element_line())+ 
  stat_compare_means(label="p.format")+    
  stat_boxplot(geom = "errorbar",width=0.5)+   
  theme(legend.position = "bottom") +
  ylab('Percentage')+
  scale_fill_manual(values=c('#c095e4','#72b043'))
