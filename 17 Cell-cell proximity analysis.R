
###### Cell paired analysis ####

### 01 Prepare data ----

library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(tibble)

set.seed(123)
setwd('...\\21-Distance Calculate')
ta.meta <- readRDS('...\\TA_TLS_Metadata.rds')
df <- ta.meta

ta.tls = ct.tls.meanta.tls = readRDS('...\\TA TLS Clustering.rds')
colnames(ta.tls)[5] = 'tls_id'
ta.tls = ta.tls[,c('group','tls_id')]

df = df %>% left_join(ta.tls)
df.b = subset(df,group == 'B-enriched')

radius <- 100
celltypes <- unique(df.b$celltype)
result_all <- data.frame()

fov_groups <- split(df.b, df.b$tls_id)
for (target_celltype in celltypes) {
  message("Processing center celltype: ", target_celltype)
  
  result <- data.frame(fov = character(),
                       center = character(),
                       neighbor = character(),
                       z_score = numeric(),
                       stringsAsFactors = FALSE)
  
  for (fov_id in names(fov_groups)) {
    fov_data <- fov_groups[[fov_id]]
    sample_size <- round(nrow(fov_data) * 0.2)
    background_cells <- fov_data[sample(1:nrow(fov_data), sample_size), ]
    
    background_neighbor_matrix <- sapply(1:nrow(background_cells), function(i) {
      bg_x <- background_cells$centroid.0[i]
      bg_y <- background_cells$centroid.1[i]
      dist <- sqrt((fov_data$centroid.0 - bg_x)^2 + (fov_data$centroid.1 - bg_y)^2)
      sapply(celltypes, function(ct) {
        sum(dist <= radius & fov_data$celltype == ct)
      })
    })
    rownames(background_neighbor_matrix) <- celltypes
    background_mean <- rowMeans(background_neighbor_matrix)
    background_sd <- apply(background_neighbor_matrix, 1, sd)
    centers <- fov_data %>% filter(celltype == target_celltype)
    if (nrow(centers) == 0) next
    
    neighbor_compositions <- sapply(1:nrow(centers), function(i) {
      cx <- centers$centroid.0[i]
      cy <- centers$centroid.1[i]
      dist <- sqrt((fov_data$centroid.0 - cx)^2 + (fov_data$centroid.1 - cy)^2)
      sapply(celltypes, function(ct) {
        sum(dist <= radius & fov_data$celltype == ct)
      })
    })
    rownames(neighbor_compositions) <- celltypes
    
    for (ct in celltypes) {
      idx <- which(celltypes == ct)
      neigh_vals <- neighbor_compositions[idx, ]
      if (is.null(dim(neigh_vals))) {
        neigh_vals <- matrix(neigh_vals, nrow = 1)
      }
      z_score <- (mean(neigh_vals) - background_mean[ct]) / background_sd[ct]
      result <- rbind(result,
                      data.frame(fov = fov_id,
                                 center = target_celltype,
                                 neighbor = ct,
                                 z_score = z_score))
    }
  }
  
  result_all <- rbind(result_all, result)
}

#write.csv(result_all,'...\\Paired Cell Distance B-TLS.csv')
result_all = read.csv('...\\Paired Cell Distance B-TLS.csv')







### 02 Calculate and Plot ----

result_all = read.csv('E:\\01-TLS\\14-Tables\\Figure2\\Paired Cell Distance T-TLS.csv')
result_all = read.csv('E:\\01-TLS\\14-Tables\\Figure2\\Paired Cell Distance B-TLS.csv')
# center cell:neutrophil neighbor cell: bile duct cell
zscore_matrix <- result_all %>%
  group_by(center, neighbor) %>%
  summarise(avg_z = mean(z_score, na.rm = TRUE)) %>%
  pivot_wider(names_from = neighbor, values_from = avg_z) %>%
  column_to_rownames("center")

library(igraph)
library(ggraph)
library(tidyverse)

target_celltype <- "Bile_Duct"
colnames(zscore_matrix)

df_net <- data.frame(
  from = target_celltype,
  to = colnames(zscore_matrix),
  weight = as.numeric(zscore_matrix[target_celltype, ])
) %>% filter(from != to)


df_net <- data.frame(
  from = target_celltype,
  to = colnames(zscore_matrix),
  weight = as.numeric(zscore_matrix[,target_celltype])
) %>% filter(from != to)

celltypes <- unique(c(df_net$from, df_net$to))

n <- nrow(df_net)
angles <- seq(-pi/2, pi/2, length.out = n)  
radius <- 1 

node_df <- data.frame(
  name = c(target_celltype, df_net$to),
  x = c(0, radius * cos(angles)),
  y = c(0, radius * sin(angles))
)

edge_df <- df_net %>%
  left_join(node_df, by = c("from" = "name")) %>%
  rename(x_from = x, y_from = y) %>%
  left_join(node_df, by = c("to" = "name")) %>%
  rename(x_to = x, y_to = y)

names(celltype_colors) <- c(df_net$to, target_celltype)

node_w <- edge_df %>%
  dplyr::select(to, weight) %>% dplyr::distinct() %>%
  dplyr::rename(name = to)
node_plot <- node_df %>%
  dplyr::left_join(node_w, by = "name") %>%
  dplyr::mutate(weight = ifelse(is.na(weight), 0, weight)) 

a = ggplot() +
  geom_curve(data = edge_df,
             aes(x = x_from, y = y_from, xend = x_to, yend = y_to,
                 colour = weight, size = abs(weight)),
             curvature = 0.3, lineend = "round", show.legend = TRUE) +
  geom_point(data = node_plot,
             aes(x = x, y = y, fill = weight),
             size = 6, shape = 21, color = "black", stroke = 0.8, show.legend = FALSE) +
  geom_text(data = node_df, aes(x = x, y = y, label = name),
            hjust = ifelse(node_df$name == target_celltype, 1.2, -0.2),
            size = 4, color = "black") +
  scale_colour_gradient2(low = "#43b0f1", mid = "white", high = "#ca0020",
                         midpoint = 0, limits = c(-0.15, 0.25), name = "Z-score") +
  scale_fill_gradient2(low = "#43b0f1", mid = "white", high = "#ca0020",
                       midpoint = 0, limits = c(-0.15, 0.25)) +
  scale_size(range = c(0,1), name = "|Z|") +
  theme_void() +
  ggtitle(paste("Cell-cell spatial z-score (center:", target_celltype, ")"))


