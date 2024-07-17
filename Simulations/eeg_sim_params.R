# Load functions
library(pheatmap)
library(viridis)
library(Matrix)
sapply((paste0("../Core/", list.files("../Core/"))), source)

# Define output
Xs <- c()
# Loop over subjects
for(i in 1:16){
  # Loop over conditions
  X <- read.csv(paste0("../EEG/sub-", sprintf("%02d", i),"/meeg/Famous.txt"))
  colnames(X) <- 1:ncol(X)
  Xs <- rbind(Xs, X)
}

# Get correlation matrix and its inverse
cormat <- cor(Xs)
invcormat <- solve(makePD(cormat))

# Function to find the largest k off-diagonal elements
sparsify <- function(matrix, k) {
  upp_tri <- which(upper.tri(matrix), arr.ind = TRUE)
  off_diag <- matrix[upp_tri]
  
  inds <- order(abs(off_diag), decreasing = TRUE)[1:k]
  S <- matrix(0, nrow(matrix), ncol(matrix))

  max_vals <- off_diag[inds]
  row_inds <- upp_tri[inds, 1]
  col_inds <- upp_tri[inds, 2]
  
  S[cbind(row_inds, col_inds)] <- max_vals
  S[cbind(col_inds, row_inds)] <- max_vals  # Because the matrix is symmetric
  return(S)
}

# Get S
S_star <- sparsify(invcormat, k=400)
S_star <- makePD(S_star)
png("eeg_S_star.png")
pheatmap(1*(S_star != 0),
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "black"))(256),
         show_rownames = FALSE, 
         show_colnames = FALSE, 
         legend = FALSE,
         border_color = NA)
dev.off()
upp_tri <- S_star[upper.tri(S_star)]
hist(upp_tri[upp_tri != 0])
diag(S_star)[1]

# Get L
load("../EEG/ests.rda")
eigL <- eigen(cormat)
png("elbow.png", res = 150, width = 800, height = 1000)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(eigL$values[1:40], xlab = "rank", ylab = "eigenvalue", 
     type = "l", cex.lab = 1.5, cex.axis = 1.5)
r <- cvslis$Famous$r
abline(v = r, col = "red")
dev.off()
dev.off()

L_star <- eigL$vectors[,1:r] %*% 
  diag(eigL$values[1:r]) %*% 
  t(eigL$vectors[,1:r])
png("eeg_L_star.png")
pheatmap(L_star,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         color = viridis(256),
         show_rownames = FALSE, 
         show_colnames = FALSE, 
         legend = FALSE,
         border_color = NA)
dev.off()

# Define output
Lout <- list()
set.seed(123)
Lout$L <- L_star
Lout$z <- kmeans(eigL$vectors[,1:r], r, nstart = 1000)$cluster

# Save parameter matrices
save(S_star, Lout, file = "eeg_sim_params.rda")