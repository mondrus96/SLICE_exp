library(MASS)
library(ggplot2)
library(reshape2)
sapply((paste0("../Core/", list.files("../Core/"))), source)

### Generate Ss and Ls ###
sapply((paste0("../Core/", list.files("../Core/"))), source)
set.seed(123)

pobs <- 100 # Number of observed variables for S
plat <- 4 # Number of latent variables for L
n <- 150 # Number of observations
simtype <- "rand"

S_star <- Smat(pobs, 2, 1.5)
S_star[S_star < 0.01] <- 0 # True sparse component

Lout <- Lrand(pobs, plat, 1.5)
L_star <- Lout$L; z_star <- Lout$z # True latent component; True cluster labels

Sigma_star <- solve(S_star + L_star) # True Sigma

X <- mvrnorm(n, rep(0, pobs), Sigma_star) # Run multivariate normal
Sigma = cov(X)

### Estimate w/ diff Sests ###
outs <- vector("list", 4)
outs[[1]] <- list(S = S_star, L = L_star)
models <- c("glasso", "clime", "gscad")
for(i in seq_along(models)){
  out <- slice(Sigma, 0.12, plat, models[i])
  out$S[abs(out$S) < 1e-5] <- 0
  outs[[i + 1]] <- list(S = out$S, L = out$L)
}
names(outs) <- c("true", models)

### Generate plots ###
for(i in seq_along(outs)){
  png(paste0(names(outs)[i], "_S.png"), width = 600, height = 600)
  heatmap(x = 1*(outs[[i]]$S != 0),
          Rowv = NA,  # Disable row clustering
          Colv = NA,  # Disable column clustering
          symm = TRUE, # Indicates the matrix is symmetric
          scale = "none", # Don't scale the data
          labRow = "", # Remove x-axis label
          labCol = "", # Remove y-axis label
          col = colorRampPalette(c("white", "black"))(256)) # Black & White color scheme
  dev.off()
  png(paste0(names(outs)[i], "_L.png"), width = 600, height = 600)
  heatmap(x = outs[[i]]$L,
          Rowv = NA,  # Disable row clustering
          Colv = NA,  # Disable column clustering
          symm = TRUE, # Indicates the matrix is symmetric
          scale = "none", # Don't scale the data
          labRow = "", # Remove x-axis label
          labCol = "", # Remove y-axis label
          col = colorRampPalette(c("white", "black"))(256)) # Black & White color scheme
  dev.off()
}
