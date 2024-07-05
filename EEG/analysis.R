# Load functions
sapply((paste0("../Core/", list.files("../Core/"))), source)

# Define output
Xs = vector("list", 3)
conds = c("Famous", "Unfamiliar", "Scrambled")
names(Xs) = conds
# Loop over subjects
for(i in 1:16){
  # Loop over conditions
  for(j in conds){
    X = read.csv(paste0("./sub-", sprintf("%02d", i),"/meeg/", j, "_hyp.txt"))
    colnames(X) = 1:ncol(X)
    Xs[[j]] = rbind(Xs[[j]], X)
  }
}

sli1 = slice(cor(Xs$Famous), 0.1, 5)
heatmap(sli1$L, Rowv = NA, Colv = NA)
heatmap(1*(sli1$S != 0), Rowv = NA, Colv = NA)

sli2 = slice(cor(Xs$Unfamiliar), 0.1, 5)
heatmap(sli2$L, Rowv = NA, Colv = NA)
heatmap(1*(sli2$S != 0), Rowv = NA, Colv = NA)

sli3 = slice(cor(Xs$Scrambled), 0.1, 5)
heatmap(sli3$L, Rowv = NA, Colv = NA)
heatmap(1*(sli3$S != 0), Rowv = NA, Colv = NA)

norm(sli1$L - sli2$L, "F")
norm(sli2$L - sli3$L, "F")
norm(sli1$L - sli3$L, "F")

norm(sli1$L, "2")
norm(sli2$L, "2")
norm(sli3$L, "2")

eigL1 = eigen(sli1$L)
eigL2 = eigen(sli2$L)
eigL3 = eigen(sli3$L)

plot(eigL1$vectors[,1], eigL1$vectors[,2])
plot(eigL2$vectors[,1], eigL2$vectors[,2])
plot(eigL3$vectors[,1], eigL3$vectors[,2])

# Create a data frame for each group
df1 <- data.frame(x = eigL1$vectors[,1], y = eigL1$vectors[,2], group = "eigL1")
df2 <- data.frame(x = eigL2$vectors[,1], y = eigL2$vectors[,2], group = "eigL2")
df3 <- data.frame(x = eigL3$vectors[,1], y = eigL3$vectors[,2], group = "eigL3")

# Combine the data frames into one
combined_df <- rbind(df1, df2, df3)
#combined_df <- rbind(df1, df3)

# Plot the data using ggplot2
ggplot(combined_df, aes(x = x, y = y, color = group)) +
  geom_point() +
  labs(title = "Scatter Plot of Eigenvectors", x = "First Eigenvector", y = "Second Eigenvector") +
  theme_minimal()


cor.test(eigL1$vectors[,1], eigL2$vectors[,1])

cor.test(eigL1$vectors[,1], eigL3$vectors[,1])

cor.test(eigL3$vectors[,1], eigL2$vectors[,1])

# Example data
set.seed(123)
# Assuming you have a matrix 'eigenvectors' where each row is a subject's eigenvector
eigenvectors <- matrix(rnorm(160), nrow = 40, ncol = 4)
# Group assignments for each subject
groups <- factor(rep(conds, each = 16))

# Perform t-SNE
eigenvectors = c()
for(i in 1:16){
  print(i)
  indx = ((100*(i-1))+1):(i*100)
  for(j in 1:3){
    sli = slice(cor(Xs[[j]][indx,]), 0.01, 2)
    L = sli$L
    eigL = eigen(L)
    eigenvectors = rbind(eigenvectors, c(j, eigL$vectors[,1]))
  }
}
set.seed(123)
tsne_result <- Rtsne(eigenvectors, dims = 2, perplexity = 15, verbose = TRUE, max_iter = 10000)
# Extract the t-SNE coordinates
tsne_data <- data.frame(tsne_result$Y)
tsne_data$Group <- groups

# Plot t-SNE results
ggplot(tsne_data, aes(x = X1, y = X2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "t-SNE of Eigenvector Data", x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
  theme_minimal()
