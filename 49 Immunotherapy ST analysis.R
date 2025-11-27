
#### Spatial Transcriptom gse238264 analysis ----

###gene signature
b.tls.up <- readRDS('E:\\01-TLS\\01-ST\\B_TLS_Sig.rds')
t.tls.up <- readRDS('E:\\01-TLS\\01-ST\\T_TLS_Sig.rds')
b_genes <- b.tls.up$gene
t_genes <- t.tls.up$gene

### spatial data
base_dir  <- "...\\GSE238264\\GSE238264_RAW"

samples <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
samples <- samples[grepl("^HCC\\d+(R|NR)$", samples, ignore.case = TRUE)]

res <- data.frame()
for (s in samples) {
  obj <- Load10X_Spatial(file.path(base_dir, s), assay = "Spatial")
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- AddModuleScore(obj, features = list(intersect(b_genes, rownames(obj))),
                        name = "b.enriched", nbin = 10, verbose = FALSE)
  obj <- AddModuleScore(obj, features = list(intersect(t_genes, rownames(obj))),
                        name = "t.enriched", nbin = 10, verbose = FALSE)
  
  res <- rbind(res, data.frame(
    sample = s,
    group  = ifelse(grepl("NR$", s, ignore.case = TRUE), "NR", "R"),
    b.enriched = median(obj$`b.enriched1`, na.rm = TRUE),
    t.enriched = median(obj$`t.enriched1`, na.rm = TRUE)
  ))
}
print(res)
colnames(res)
a <- ggplot(res,aes(x=group,y=t.enriched,fill=group))+geom_boxplot(outlier.shape = NA,width=0.5)+
  geom_jitter(width =0.3,size=1,color = "grey50", alpha = 0.6)+ 
  #scale_shape_manual(values = c(15, 16))+
  theme(panel.background = element_blank(),axis.line = element_line())+ 
  stat_compare_means(label="p.format")+    
  stat_boxplot(geom = "errorbar",width=0.5)+   
  #labs(title = "B_C0_IGHA1")+   
  theme(legend.position = "bottom") +
  ylab('Percentage')+
  scale_fill_manual(values=c('#c095e4','#72b043'))

pdf(file = 'total t tls.pdf',
    height =3.5,width =2.5)
print(a)
dev.off()

hcc2r <- Load10X_Spatial('...\\GSE238264\\GSE238264_RAW\\HCC2R')
hcc2r <- AddModuleScore(hcc2r,
                     features = list(b_genes),
                     name = 'b.enriched',nbin = 10)
hcc2r <- AddModuleScore(hcc2r,
                     features = list(t_genes),
                     name = 't.enriched',nbin = 10)
a <- SpatialFeaturePlot(object = hcc2r, features = c('t.enriched1'),
                       alpha = c(0.1, 1), ncol = 1,min.cutoff =0)

a <- SpatialFeaturePlot(object = hcc2r, features = c('b.enriched1'),
                       alpha = c(0.1, 1), ncol = 1,max.cutoff =75)

pdf(file = 'Response hcc2r t enriched.pdf')
print(a)
dev.off()

