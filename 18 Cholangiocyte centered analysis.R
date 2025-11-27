
###### The co-occurrence of Bile duct and neutrophil ######

total.cell <- readRDS('.../01 Total Cell Metadata-0318.rds')
ta.tls.meta <- readRDS('.../TA_TLS_Metadata.rds')

fovs <- unique(ta.tls.meta$fov)
total.cell <- subset(total.cell, fov %in% fovs)

ta.bile.duct <- ta.tls.meta[ta.tls.meta$celltype == 'Bile_Duct', 'cellid']
bile.cells <- total.cell[total.cell$cellid %in% ta.bile.duct, ]
total.cell$centroid.0 <- as.numeric(total.cell$centroid.0)
total.cell$centroid.1 <- as.numeric(total.cell$centroid.1)

radii <- seq(50, 300, by = 10)
fov_list <- split(total.cell, total.cell$fov)
score_list <- lapply(1:nrow(bile.cells), function(i) {
  cell <- bile.cells[i, ]
  fov_data <- fov_list[[as.character(cell$fov)]]
  scores <- sapply(radii, function(r) {
    dists <- sqrt((fov_data$centroid.0 - cell$centroid.0)^2 + (fov_data$centroid.1 - cell$centroid.1)^2)
    in_radius <- fov_data[dists <= r, ]
    if (nrow(in_radius) == 0) return(NA)
    p_local <- sum(in_radius$new.celltype == "Neutrophil") / nrow(in_radius)
    p_global <- sum(fov_data$new.celltype == "Neutrophil") / nrow(fov_data)
    if (p_global == 0) return(NA) else return(p_local / p_global)
  })
  names(scores) <- paste0("r", radii)
  return(scores)
})
score_df <- do.call(rbind, score_list)
rownames(score_df) <- bile.cells$cellid
write.csv(score_df, ".../bile_duct_cooccurrence_score.csv")

score_df <- score_df[complete.cases(score_df), ]           
score_df <- score_df[apply(score_df, 1, function(x) all(is.finite(x))), ]

row_sds <- apply(score_df, 1, sd)
sum(row_sds == 0)

score_df <- score_df[apply(score_df, 1, sd) != 0, ]
score_df2 <- t(apply(score_df, 1, scale))  
pheatmap(score_df2,cluster_cols = FALSE)

row_dist <- dist(score_df2, method = "euclidean")  
row_clust <- hclust(row_dist, method = "ward.D2")  

cluster_k <- 4
row_clusters <- cutree(row_clust, k = cluster_k)

score_df2 = as.data.frame(score_df2)
score_df2$cluster <- as.factor(row_clusters)
score_df2$id <- rownames(score_df2)

score_df2 <- score_df2[order(score_df2$cluster), ]
heatmap_data <- score_df2[, c(1:26)]
annotation_row <- data.frame(cluster = as.factor(score_df2$cluster))
rownames(annotation_row) <- rownames(score_df2)
a = pheatmap(heatmap_data,
             annotation_row = annotation_row,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             gaps_row = cumsum(table(score_df2$cluster)),  
             show_rownames = FALSE,
             show_colnames = TRUE,
             scale = "none",
             color = colorRampPalette(c("#b0ead7", "white", "#b6042a"))(100)) 

pdf(file = '...\\pheatmap cluster.pdf',width = 4,height = 5)
print(a)
dev.off()

df_long <- melt(score_df2, id.vars = c("id", "cluster"))
a = ggplot(df_long, aes(x = variable, y = value, group = id, color = cluster)) +
  geom_line(alpha = 0.4) +
  facet_wrap(~ cluster, scales = "free_y") +
  theme_classic() +
  labs(title = "Sample trends by cluster", x = "Feature", y = "Scaled value")

pdf(file = '...\\Trend cluster.pdf',width = 6,height = 5)
print(a)
dev.off()

score_df2.result = score_df2[,c('cluster','id')]
colnames(score_df2.result)[2] = 'cellid'
bile.cells = bile.cells %>% left_join(score_df2.result)
ta.tls = readRDS('...\\TA TLS Clustering.rds')
colnames(ta.tls)[5] = 'tls_id'
ta.tls = ta.tls[,c('group','tls_id')]

bile.cells = bile.cells %>% left_join(ta.tls)
tls_ids_to_keep <- bile.cells %>%
  count(tls_id) %>%
  filter(n > 10) %>%
  pull(tls_id)
filtered_data <- bile.cells %>% filter(tls_id %in% tls_ids_to_keep)
filtered_data = filtered_data[!is.na(filtered_data$cluster),] 
filtered_data$cluster = case_when(filtered_data$cluster == 1~'Cluster1',
                                  filtered_data$cluster == 2~'Cluster2',
                                  filtered_data$cluster == 3~'Cluster3',
                                  filtered_data$cluster == 4~'Cluster4')
unique(filtered_data$slide)
write.csv(filtered_data,'...\\Bile duct type.csv')

filtered_data2 = filtered_data
filtered_data2$cluster = ifelse(filtered_data$cluster %in% c("Cluster1","Cluster2","Cluster3"),'co-localization','Bystander')

proportion_df <- filtered_data2 %>%
  group_by(tls_id, cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tls_id) %>%
  mutate(proportion = count / sum(count))%>%
  filter(cluster %in% c("co-localization") ) 

group_info <- filtered_data2 %>%
  select(tls_id, group) %>%
  distinct()

final_df <- left_join(proportion_df, group_info, by = "tls_id")
colnames(final_df)
a = ggboxplot(final_df, x = "group", y = "proportion", color = "group", palette = "jco",add = "jitter") + 
  stat_compare_means(method = 'wilcox.test',label = "p.format")+
  theme(legend.position = "top") +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank())

pdf(file = '...\\co localzation bystander.pdf',width =3,height = 4)
print(a)
dev.off()
save.image('...\\Bile duct neutro relationship.rdata')
