library(ggplot2)
library(dplyr)
library(pheatmap)
library(viridis)
library(colorspace)
# Load data and functions
sapply((paste0("../Core/", list.files("../Core/"))), source)
load("ests.rda")

# Look at rank
sapply(cvslis, "[[", "r")

### Intragroup statistics
Ls <- lapply(slis, "[[", 2) # Get Ls
intra_df <- c()
for(i in seq_along(Ls)){
  r <- slis[[i]]$rank
  L <- Ls[[i]]
  
  # Eigendecomp
  eigL <- eigen(L)
  df <- as.data.frame(eigL$vectors[,1:2])
  df <- cbind(df, group = names(Ls)[i])
  
  # Norms
  frob_norm <- norm(L, "F")
  spec_norm <- norm(L, "2")
  
  # Clustering
  clusts <- table(kmeans(eigL$vectors[,1:r], r, 1000)$cluster)
  min_clust <- min(clusts)
  max_clust <- max(clusts)
  
  # Append
  row <- c(frob_norm, spec_norm, min_clust, max_clust)
  intra_df <- rbind(intra_df, row)
}
intra_df <- as.data.frame(intra_df)
intra_df <- cbind(names(slis), intra_df)
colnames(intra_df) <- c("group", "frob_norm", "spec_norm", 
                        "min_clust", "max_clust")
rownames(intra_df) <- NULL
print(intra_df)

# Labels and palette
AAL <- read.csv("AAL.txt")
vox_ctrs <- read.csv("voxel_centers.txt", header = FALSE)

# For each row in vox_ctrs, find the closest row in AAL
vox_ctrs$AAL_Label <- apply(vox_ctrs, 1, function(vox_row) {
  # Calculate distances to all points in AAL
  distances <- mapply(eucl_dist,
                      x1 = vox_row[1], y1 = vox_row[2], z1 = vox_row[3],
                      x2 = AAL$X, y2 = AAL$Y, z2 = AAL$Z)
  
  # Find the index of the minimum distance
  nearest_index <- which.min(distances)
  
  # Return the corresponding AAL label
  return(AAL$Group[nearest_index])
})

labels <- data.frame(Group = vox_ctrs$AAL_Label)
palette <- make_palette(length(unique(labels$Group))) # Make colours for all
group_colors <- list(Group = (setNames(palette, unique(labels$Group))))

### Intergroup statistics
inter_df <- matrix(c(1, 2, 1, 2, 3, 3), 3, 2) # Groups to compare
inter_df <- as.data.frame(inter_df)
sin_theta <- frob_norm <- spec_norm <- adj_rand_indx <- c()
for(i in 1:nrow(inter_df)){
  set.seed(123)
  r <- slis[[i]]$rank
  g1 <- inter_df[i,1]; g2 <- inter_df[i,2]
  L1 <- Ls[[g1]]; L2 <- Ls[[g2]]
  eigL1 <- eigen(L1); eigL2 <- eigen(L2)
  
  sin_theta <- c(sin_theta,
                 sintheta(eigL1$vectors[,1], eigL2$vectors[,1]))
  
  frob_norm <- c(frob_norm, norm(L1 - L2, type = "F"))
  spec_norm <- c(spec_norm, norm(L1 - L2, type = "2"))
  
  clust1 <- kmeans(eigL1$vectors[,1:r], r, 10000)$cluster
  clust2 <- kmeans(eigL2$vectors[,1:r], r, 10000)$cluster
  adj_rand_indx <- c(adj_rand_indx, ari(clust1, clust2))
  
  png(paste0(g1, g2, "_L.png"), width = 1600, height = 1200, res = 200)
  diffL <- (L1 - L2)^2
  rownames(diffL) <- colnames(diffL) <- rownames(labels)
  pheatmap(diffL,
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           color = viridis(256),
           breaks = seq(0.01, 0.02, length.out = 257),
           legend = FALSE,
           show_rownames = FALSE, 
           show_colnames = FALSE,
           border_color = NA,
           annotation_col = labels,
           annotation_colors = group_colors,
           cellheight=3, cellwidth=3)
  dev.off()
}
inter_df <- cbind(inter_df, sin_theta, frob_norm, spec_norm, adj_rand_indx)
print(inter_df)