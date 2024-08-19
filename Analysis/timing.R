# Timing comparison
library(R.utils)
sapply((paste0("../Core/", list.files("../Core/"))), source)

# Settings
pobs <- c(100, 250, 500, 1000, 2500)
plat <- 7
rho <- 0.01
gamma <- 0.001

# Run SLICE function
S_star <- Smat(500, 2, 1.5); S_star[S_star < 0.01] <- 0 # Sparse
Lout <- Lrand(500, plat, 1.5); L_star <- Lout$L # Latent
Sigma <- Matrix::chol2inv(Matrix::chol(S_star + L_star)) # True Sigma
rm(L_star, Lout, S_star)

# Check profile of function
Rprof("slice_profile.out")
out <- slice(Sigma, rho, plat)
Rprof(NULL)
summary <- summaryRprof("slice_profile.out")
print(summary)

# Compare timing between estimators
models <- c("SLICE_glasso", "SLICE_huge_glasso", "rcLVGLASSO", "nnLVGLASSO")
out <- matrix(NA, nrow = length(pobs), ncol = length(models))
colnames(out) <- models; rownames(out) <- pobs

runtiming <- function(Sigma, estimator, ...){
  withTimeout({
    startT <- Sys.time()
    estimator(Sigma, ...)
    endT <- Sys.time()
    return(difftime(endT, startT, units = "secs"))
  }, timeout = 600, onTimeout = "warning")  # Timeout set 10 minutes
}

set.seed(123)
for(i in seq_along(pobs)){
  # Define simulation
  S_star <- Smat(pobs[i], 2, 1.5); S_star[S_star < 0.01] <- 0 # Sparse
  Lout <- Lrand(pobs[i], plat, 1.5); L_star <- Lout$L # Latent
  Sigma <- Matrix::chol2inv(Matrix::chol(S_star + L_star)) # True Sigma
  rm(L_star, Lout, S_star)
  
  # Loop through models
  for(j in seq_along(models)){
    if(models[j] == "SLICE_glasso"){
      T <- runtiming(Sigma, slice, rho, plat, Sest = "glasso")
    } else if(models[j] == "SLICE_huge_glasso"){
      T <- runtiming(Sigma, slice, rho, plat, Sest = "huge_glasso")
    } else if(models[j] == "rcLVGLASSO"){
      T <- runtiming(Sigma, rclvg, rho, plat)
    } else if(models[j] == "nnLVGLASSO"){
      T <- runtiming(Sigma, nnlvg, rho, gamma)
    }
    T <- ifelse(is.null(T), NA, as.numeric(T))
    out[i, j] <- T
    print(out)
  }
}