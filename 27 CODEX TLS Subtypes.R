
#### TLS Subtypes ####

load('...\\HD TLS CODEX Celltype.RDATA')
rownames(tls.S4) = paste0('S4_',c("N_1_2","N_2_1","N_3","N_4","N_5","N_6","P_1","P_2","P_3","T_1"))
rownames(tls.S2) = paste0('S2_',c("P_1" ,"P_2", "P_3", "P_4", "P_5" ,"P_6" ,"T_1"))
rownames(tls.S1) = paste0('S1_',c( "P_3", "P_4","T_2",'T_3'))
rownames(tls.S5) = paste0('S5_',c("N_1","P_1","P_2","P_3","P_4","P_5","T_1","T_2","T_3","T_4","T_5"))
rownames(tls.S6) = paste0('S6_',rownames(tls.S6))

tls.S4$size = apply(tls.S4,1,sum)
tls.S2$size = apply(tls.S2,1,sum)
tls.S1$size = apply(tls.S1,1,sum)
tls.S5$size = apply(tls.S5,1,sum)
tls.S6$size = apply(tls.S6,1,sum)

tls.S4$t_cluster = tls.S4$CD4T+tls.S4$CD8T
tls.S2$t_cluster = tls.S2$CD4T+tls.S2$CD8T
tls.S1$t_cluster = tls.S1$CD4T+tls.S1$CD8T
tls.S5$t_cluster = tls.S5$CD4T+tls.S5$CD8T
tls.S6$t_cluster = tls.S6$CD4T+tls.S6$CD8T

tls.S4 = tls.S4[,c('B','t_cluster','size')]
tls.S2 = tls.S2[,c('B','t_cluster','size')]
tls.S1 = tls.S1[,c('B','t_cluster','size')]
tls.S5 = tls.S5[,c('B','t_cluster','size')]
tls.S6 = tls.S6[,c('B','t_cluster','size')]

tls.size.total = rbind(tls.S4,
                       tls.S2,
                       tls.S1,
                       tls.S5,
                       tls.S6)

hd.tls.meta = readRDS('...\\HD Subtype 0711.rds')
ta.tls.hd = hd.tls.meta[c(5:18,20:29,38:47),]
ta.tls.hd = ta.tls.hd[,-9]

ta.tls.hd$t.total = ta.tls.hd$CD4T+ta.tls.hd$CD8T
ta.tls.hd = ta.tls.hd[,-c(3,4)]

rownames(ta.tls.hd) %in% rownames(tls.size.total)
tls.size.total = tls.size.total[rownames(ta.tls.hd),]
tls.size.total$b.per = ta.tls.hd$B
tls.size.total$t.per = ta.tls.hd$t.total
tls.size.total$tb = log((tls.size.total$b.per+1)/(tls.size.total$t.per+1))

a <- pheatmap::pheatmap(as.matrix(ta.tls.hd), silent = TRUE)
library(pheatmap)
row_hc <- a$tree_row
plot(row_hc, main = "Row Clustering Dendrogram from pheatmap")
groups <- cutree(row_hc, k = 2)
ta.tls.hd$group <- groups[rownames(ta.tls.hd)]
ta.tls.hd$group = ifelse(ta.tls.hd$group == 1,'T-TLS','B-TLS')
table(ta.tls.hd$group)

library(RedeR)
library(RColorBrewer)
library(igraph)
library(TreeAndLeaf)
tal <- treeAndLeaf(row_hc)
tal <- att.mapv(tal, tls.size.total, refcol = 0)

pal <- c(rev(brewer.pal(5, "Greens")),brewer.pal(5, "Reds"))
tal <- att.setv(g = tal, from = "tb", to = "nodeColor",
                cols = pal, nquant = 40)
tal <- att.setv(g = tal, from = "size", to = "nodeSize")
tal <- att.adde(tal, "edgeWidth", value =15)
tal <- att.addv(tal, "nodeColor", value = "black", index=!V(tal)$isLeaf) 

rdp <- RedPort()
calld(rdp)
resetd(rdp)
addGraph(obj = rdp, g = tal, gzoom=20)
relax(rdp, p1=25, p2=200, p3=10, p4=100, p5=10, ps=TRUE)
addLegend.color(obj = rdp, tal, title = "B cell enrichment",
                position = "bottomright")
addLegend.size(obj = rdp, tal, title = "TLS Number")

## Each Branch plot
groups <- cutree(row_hc, k = 2)
ta.tls.hd$group <- groups[rownames(ta.tls.hd)]
ta.tls.hd$group = ifelse(ta.tls.hd$group == 1,'T-TLS','B-TLS')
tls.size.total$group =ta.tls.hd$group

red_line <- median(tls.size.total$tb)
a = ggplot(tls.size.total, aes_string(x = 'group', y = 'tb')) +
  geom_boxplot(outlier.shape = 21, width = 0.5, aes(fill = ifelse(..middle.. > red_line, "#C25E72", "#D9D9D9"))) +
  theme_minimal() +
  geom_jitter(width = 0.2, alpha = 0.5, size =1, color = "#919191") +
  theme(
    axis.title = element_text(size = 13, color = "black", face = "bold"),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.line.x = element_line(color = "black", size = 0.5),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  labs(x = '', y = 'Percentage') +
  stat_compare_means(aes(label = "p.format"), method = "wilcox.test", label = "p.format") +
  geom_hline(yintercept = red_line, color = 'red', linetype = 'dashed') +
  scale_fill_identity()

pdf(file = '...\\tb Ratio.pdf',height = 2.5,width = 2)
print(a)
dev.off()

setwd('...\\07 HD TLS Subtpye')
ct = readRDS('CODEX-identified ta-TLS group.rds') 

a = ggplot(ct,aes(x=group,y=Neutrophil,fill=group))+geom_boxplot(outlier.shape = NA,width=0.5)+
  geom_jitter(width =0.3,size=1,color = "grey50", alpha = 0.6)+ 
  #scale_shape_manual(values = c(15, 16))+
  theme(panel.background = element_blank(),axis.line = element_line())+ 
  stat_compare_means(label="p.format")+    
  stat_boxplot(geom = "errorbar",width=0.5)+   
  theme(legend.position = "bottom") +
  ylab('Percentage')+
  scale_fill_manual(values=c('#c095e4','#72b043'))
pdf(file = 'Neutrophil.pdf',height = 4,width = 3)
print(a)
dev.off()
