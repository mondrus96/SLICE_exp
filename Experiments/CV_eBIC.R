library(MASS)
sapply((paste0("../Models/", list.files("../Models/"))), source)
sapply((paste0("../Simulations/", list.files("../Simulations/"))), source)

set.seed(123)
pobs <- 150 # Number of observed variables for S
plat <- 4 # Number of latent variables for L
n <- 10000

S_star <- Smat(pobs, 2, 1.5)
Lout <- Lexp(pobs, plat, 1.5)
L_star <- Lout$L; z_star <- Lout$z
Sigma_star <- solve(S_star + L_star)
X <- mvrnorm(n, rep(0, pobs), Sigma_star)
Sigma = cov(X)

lambdas = c(1e-4, 1e-3, 1e-2, 0.1, 0.2)
rs = 2:8

CVout <-cv.slice(X, 5, lambdas, rs) # Cross validation

ebicmat <- matrix(NA, length(rs), length(lambdas))
rownames(ebicmat) <- rs; colnames(ebicmat) <- lambdas
for(i in 1:length(rs)){
  print(paste0("rank: ", rs[i]))
  for(j in 1:length(lambdas)){
    print(paste0("lambda: ", lambdas[j]))
    
    out <- slice(Sigma, lambdas[j], rs[i]) # Run method
    likl <- logL(Sigma, out$S, out$L)
    
    k_S <- sum(out$S[upper.tri(out$S)] != 0); k_L <- rs[i] + (pobs*rs[i]) - (rs[i]*(rs[i] - 1)/2)
    ebicmat[i, j] <- ebic(likl, pobs, n, k_S + k_L)
  }
}
best <- which(ebicmat == min(ebicmat), arr.ind = TRUE)
