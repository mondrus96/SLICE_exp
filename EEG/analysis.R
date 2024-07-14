library(ggplot2)
library(dplyr)

# Load data and functions
sapply((paste0("../Core/", list.files("../Core/"))), source)
load("ests.rda")

# Look at rank
sapply(cvslis, "[[", "r")

# Intragroup statistics
Ls <- lapply(slis, "[[", 2) # Get Ls
eigvecs_df <- list()
intra_df <- c()
for(i in seq_along(Ls)){
  r <- slis[[i]]$rank
  L <- Ls[[i]]
  
  # Eigendecomp
  eigL <- eigen(L)
  df <- as.data.frame(eigL$vectors[,1:2])
  df <- cbind(df, group = names(Ls)[i])
  eigvecs_df <- rbind(eigvecs_df, df)
  
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

# Scatter plot of eigenvectors
ggplot(eigvecs_df, aes(x = V1, y = V2, color = group)) +
  geom_point(size = 3) +
  labs(title = "Scatter Plot of Eigenvectors", x = "X1", y = "X2") +
  theme_minimal()

# Intergroup statistics
inter_df <- matrix(c(1, 2, 1, 2, 3, 3), 3, 2) # Groups to compare
inter_df <- as.data.frame(inter_df)
sin_theta <- frob_norm <- spec_norm <- adj_rand_indx <- c()
for(i in 1:nrow(inter_df)){
  g1 <- inter_df[i,1]; g2 <- inter_df[i,2]
  L1 <- Ls[[g1]]; L2 <- Ls[[g2]]
  eigL1 <- eigen(L1); eigL2 <- eigen(L2)
  
  sin_theta <- c(sin_theta,
                 sintheta(eigL1$vectors[,1], eigL2$vectors[,1]))
  
  frob_norm <- c(frob_norm, norm(L1 - L2, type = "F"))
  spec_norm <- c(spec_norm, norm(L1 - L2, type = "2"))
  
  clust1 <- kmeans(eigL1$vectors[,1:r], r, 1000)$cluster
  clust2 <- kmeans(eigL2$vectors[,1:r], r, 1000)$cluster
  adj_rand_indx <- c(adj_rand_indx, ari(clust1, clust2))
}
inter_df <- cbind(inter_df, sin_theta, frob_norm, spec_norm, adj_rand_indx)
print(inter_df)
save(inter_df, intra_df, eigvecs_df, file = paste0(rho, "out.rda"))
