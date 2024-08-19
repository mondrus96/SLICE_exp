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

# Generate sample covariance
X <- mvrnorm(1000, rep(0, pobs), Sigma_star)
Sigma <- cov(X)

### Estimate w/ diff Sests ###
outs <- vector("list", 4)
outs[[1]] <- list(S = S_star, L = L_star)
models <- c("glasso", "clime", "gscad")
rhos <- c(0.05, 0.07, 0.05)
for(i in seq_along(models)){
  out <- slice(Sigma, rhos[i], plat, Sest = models[i])
  out$S[abs(out$S) < 1e-5] <- 0
  outs[[i + 1]] <- list(S = out$S, L = out$L)
}
names(outs) <- c("true", models)

### Generate plots ###
for(i in seq_along(outs)){
  png(paste0(names(outs)[i], "_S.png"), width = 600, height = 600)
  pheatmap(1*(outs[[i]]$S != 0), 
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           color = colorRampPalette(c("white", "black"))(256),
           show_rownames = FALSE, 
           show_colnames = FALSE, 
           legend = FALSE,
           border_color = NA)
  dev.off()
  
  png(paste0(names(outs)[i], "_L.png"), width = 600, height = 600)
  pheatmap(outs[[i]]$L, 
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           color = colorRampPalette(c("white", "black"))(256),
           show_rownames = FALSE, 
           show_colnames = FALSE, 
           legend = FALSE,
           border_color = NA)
  dev.off()
}
