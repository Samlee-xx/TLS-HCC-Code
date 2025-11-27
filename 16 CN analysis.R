
## cell neighborhood analysis ----
setwd('...\\25-CN Analysis')
total.meta = readRDS('...\\TA_TLS_Metadata.rds')

library(dplyr)
results_list <- list()
fovs <- unique(total.meta$tls_id)
for (f in fovs) {
  message("Processing FOV: ", f)
  df <- total.meta %>% filter(tls_id == f)
  
  fov_results <- lapply(1:nrow(df), function(i) {
    this_cell <- df[i, ]
    
    distances <- sqrt((df$centroid.0 - this_cell$centroid.0)^2 + 
                        (df$centroid.1 - this_cell$centroid.1)^2)
    
    neighbors <- df[distances <= 100 & distances > 0, ]
    
    if (nrow(neighbors) == 0) return(NULL)
    
    celltype_counts <- table(neighbors$celltype)
    out <- as.data.frame(t(celltype_counts))
    out$cellid <- this_cell$cellid
    return(out)
  })
  
  fov_df <- bind_rows(fov_results)
  results_list[[f]] <- fov_df
}

final_result <- bind_rows(results_list)
colnames(final_result)
final_result = final_result[,-1]
final_result_long = dcast(final_result,cellid~Var2,value.var = 'Freq')

rownames(final_result_long) = final_result_long$cellid
final_result_long = final_result_long[,-1]
final_result_long[is.na(final_result_long)] = 0

final_prop <- final_result_long / rowSums(final_result_long)

wss <- numeric()
max_k <- 15  

for (k in 1:max_k) {
  kmeans_result <- kmeans(final_prop, centers = k, nstart = 10)
  wss[k] <- kmeans_result$tot.withinss
}

plot(1:max_k, wss, type = "b", pch = 19,
     xlab = "Number of clusters (K)",
     ylab = "Total within-cluster sum of squares",
     main = "Elbow Method for Determining K")

K =8
set.seed(123) 
kmeans_result <- kmeans(final_prop, centers = K, nstart = 20)
table(kmeans_result$cluster)

final_prop$group = kmeans_result$cluster
final.plot = final_prop%>%group_by(group) %>% summarise_each(funs =mean)
final.plot = final.plot[,-1]
rownames(final.plot) = paste0('CN_',rownames(final.plot))
pheatmap(scale(final.plot),cluster_rows = FALSE,cluster_cols = FALSE,cellwidth = 15,
         cellheight = 15,
         colorRampPalette(c('blue',"yellow", "red"))(100))
pheatmap(scale(final.plot),cluster_rows = FALSE,cluster_cols = FALSE,cellwidth = 15,
         cellheight = 15)

breaks <- c(seq(-2, 0, length.out = 50), seq(0.01, 3, length.out = 50))
colors <- colorRampPalette(c("#51a16a", "white", "red"))(99)

a = pheatmap(
  scale(final.plot),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 15,
  cellheight = 15,
  color = colors,
  breaks = breaks
)
a

pdf(file = 'Distance 100 k 8 CN.pdf',width = 8,height = 8)
print(a)
dev.off()

sum(!rownames(final_prop) == total.meta$cellid)
total.meta = total.meta[rownames(final_prop),]
total.meta$group = kmeans_result$cluster
table(total.meta$group)

total.meta$group = case_when(total.meta$group == 1~'CN_CD8T_DC',
                             total.meta$group == 2~'CN_Bile_duct_Fibroblast',
                             total.meta$group == 3~'CN_CD4T_B',
                             total.meta$group == 4~'CN_Lymphatic_Unassign',
                             total.meta$group == 5~'CN_CD4T',
                             total.meta$group == 6~'CN_B',
                             total.meta$group == 7~'CN_Macrophage',
                             total.meta$group == 8~'CN_Neutrophil_Bile_duct')

total.meta$group = paste0('CN_',total.meta$group)

table(total.meta$group)
ta.tls = readRDS('...\\TA TLS Clustering.rds')
ta.tls = ta.tls[,c('tls.id','group')]
colnames(ta.tls)[1] = 'tls_id'
colnames(ta.tls)[2] = 'tls_group'
total.meta = total.meta%>% left_join(ta.tls)

group_count <- total.meta %>%
  dplyr::group_by(tls_group, group) %>%
  dplyr::summarise(count = n(), .groups = "drop")

roe = dcast(group_count,tls_group~group) %>% as.data.frame()
rownames(roe) = roe$tls_group
roe = roe[,-1]

a=roe
p.value = roe
for (i in 1:dim(a)[1]){
  for (j in 1:dim(a)[2]){
    onecelltype_onetissue=a[i,j]
    onetissue_othercelltype=sum(a[i,])-onecelltype_onetissue
    onecelltype_othertissue=sum(a[,j])-onecelltype_onetissue
    othercelltype_ohertissue=sum(a) - onecelltype_onetissue - onetissue_othercelltype - onecelltype_othertissue
    x=matrix(c(onecelltype_onetissue,onecelltype_othertissue,onetissue_othercelltype,othercelltype_ohertissue),
             nrow = 2,ncol = 2)
    b=chisq.test(x)$expected[1,1]
    m=onecelltype_onetissue/b
    roe[i,j]=m
    p.value[i,j] = chisq.test(x)$p.value
  }
} 

a = pheatmap(
  as.matrix(roe),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 15,
  cellheight = 15,
  color = colorRampPalette(c("gray90", "white", "red"))(100)
)
pdf(file = 'ROE cn tb tls.pdf',width = 8,height = 8)
print(a)
dev.off()
