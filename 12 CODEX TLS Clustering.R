
###### TLS Clustering ######
### 01 Prepare Data ----
total.meta = readRDS('...\\TA_TLS_Metadata.rds')
celltype.table = total.meta[,c('tls_id','celltype')]
celltype.table = dcast(celltype.table,tls_id~celltype)
rownames(celltype.table) = celltype.table$tls_id
colnames(celltype.table)
celltype.table = celltype.table[,-c(1)]

tls.size = apply(celltype.table, 1, sum)
celltype.table = celltype.table[,-13]

normalize = function(x){
  x/sum(x)
}

celltype.table = apply(celltype.table, 1, normalize)
celltype.table = celltype.table %>% t() %>% as.data.frame()
celltype.table = celltype.table*100

back.info = str_split(rownames(celltype.table),'_',simplify = TRUE)
back.info = back.info %>% as.data.frame()
celltype.table$location = back.info$V2

back.info$stage = case_when(back.info$V1 =='123942'~'Early',
                            back.info$V1 =='124989'~'Early',
                            back.info$V1 =='133698'~'Late',
                            back.info$V1 =='134024'~'Middle',
                            back.info$V1 =='134205'~'Early',
                            back.info$V1 =='136254'~'Middle',
                            back.info$V1 =='136759'~'Middle',
                            back.info$V1 =='176552'~'Late',
                            back.info$V1 =='86130'~'Early')
celltype.table$stage = back.info$stage

celltype.table.meta = data.frame(B = celltype.table$B,size = tls.size,`T` = celltype.table$CD4T+celltype.table$CD8T,
                                 stage = celltype.table$stage)

rownames(celltype.table.meta) = rownames(celltype.table)
celltype.table.meta$size = scale(celltype.table.meta$size)
str(celltype.table.meta)
colnames(celltype.table.meta)


colnames(celltype.table)
celltype.table = celltype.table[,-c(13,14)]

celltype.table.meta$tls.id = rownames(celltype.table.meta)
rownames(celltype.table.meta) = c(1:127)
rownames(celltype.table) = c(1:127)

celltype.table.meta$tb = log((celltype.table.meta$B+1)/(celltype.table.meta$`T`+1))
celltype.table$group = celltype.table.meta$group
write.csv(celltype.table,'...\\TA-TLS Metadata.csv')

### 02 Tree Plot ----

hc <- hclust(dist(celltype.table), "complete")
celltype.table.meta$group = cutree(hc, h = 57)
unique(celltype.table.meta$group)
celltype.table.meta$fixed_size  = 40
tal <- treeAndLeaf(hc)
tal <- att.mapv(tal, celltype.table.meta, refcol = 0)

pal <- c(rev(brewer.pal(5, "Greens")),brewer.pal(5, "Reds"))
pal = brewer.pal(6,'Set1')
colnames(celltype.table)
tal <- att.setv(g = tal, from = "tb", to = "nodeColor",
                cols = pal, nquant = 40)

tal <- att.setv(g = tal, from = "group", to = "nodeColor",
                cols = pal, nquant = 20)

tal <- att.setv(g = tal, from = "fixed_size", to = "nodeSize")
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

### 03 Clustering Plot each Branch ----

table(celltype.table.meta$group)
celltype.table.meta$group = case_when(celltype.table.meta$group == '1'~ 'Branch1',
                                      celltype.table.meta$group == '2'~ 'Branch2',
                                      celltype.table.meta$group == '3'~ 'Branch3',
                                      celltype.table.meta$group == '4'~ 'Branch4',
                                      celltype.table.meta$group == '5'~ 'Branch5')

table(celltype.table.meta$group)

colnames(celltype.table.meta)
red_line <- median(celltype.table.meta$tb)

a = ggplot(celltype.table.meta, aes_string(x = 'group', y = 'tb')) +
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
  stat_compare_means(aes(label = "p.format"), method = "anova", label = "p.format") +
  geom_hline(yintercept = red_line, color = 'red', linetype = 'dashed') +
  scale_fill_identity()

celltype.table.meta$group = case_when(celltype.table.meta$group %in% c('Branch1','Branch4')~ 'B-enriched',
                                      celltype.table.meta$group %in% c('Branch2','Branch3','Branch5')~ 'T-enriched')
saveRDS(celltype.table.meta,'TA TLS Clustering.rds')

