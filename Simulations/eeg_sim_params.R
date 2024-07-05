# Load functions
sapply((paste0("../Core/", list.files("../Core/"))), source)

# Define output
Xs <- c()
# Loop over subjects
for(i in 1:16){
  # Loop over conditions
  X <- read.csv(paste0("../EEG/sub-", sprintf("%02d", i),"/meeg/Famous_hyp.txt"))
  colnames(X) <- 1:ncol(X)
  Xs <- rbind(Xs, X)
}

# Get correlation matrix and its inverse
cormat <- cor(Xs)
diag(cormat) <- diag(cormat) + 0.1
invcormat <- solve(cormat)

# Get sparse component
heatmap(invcormat, Rowv = NA, Colv = NA)

# Load required package
library(Matrix)

# Function to find the largest k off-diagonal elements
sparsify <- function(matrix, k=100) {
  upp_tri <- which(upper.tri(matrix), arr.ind = TRUE)
  off_diag <- matrix[upp_tri]
  
  inds <- order(-off_diag)[1:k]
  S <- matrix(0, nrow(matrix), ncol(matrix))

  max_vals <- off_diag[inds]
  row_inds <- upp_tri[inds, 1]
  col_inds <- upp_tri[inds, 2]
  
  S[cbind(row_inds, col_inds)] <- max_vals
  S[cbind(col_inds, row_inds)] <- max_vals  # Because the matrix is symmetric
  S <- makePD(S)
  return(S)
}

# Get S
S_star <- sparsify(invcormat, k=100)
heatmap(1*(S_star != 0), Rowv = NA, Colv = NA)

# Get L
eigL <- eigen(cormat)
plot(eigL$values[1:40], type = "l")
abline(v = 6, col = "red")
L_star <- eigL$vectors[,1:6] %*% diag(1/eigL$values[1:6]) %*% t(eigL$vectors[,1:6])
heatmap(L_star, Rowv = NA, Colv = NA)

# Define output
Lout <- list()
set.seed(123)
Lout$L <- L_star
Lout$z <- kmeans(eigL$vectors[,1:6], 6, nstart = 1000)$cluster

# Save parameter matrices
save(S_star, Lout, file = "eeg_sim_params.rda")