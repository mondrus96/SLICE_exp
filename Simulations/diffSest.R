sapply((paste0("../Core/", list.files("../Core/"))), source)

pobs <- 100 # Number of observed variables for S
plat <- 4 # Number of latent variables for L
n <- 500 # Number of observations
simtype <- "rand"
iters <- 1:100

S_star <- Smat(pobs, 2, 1.5)
S_star[S_star < 0.01] <- 0 # True sparse component

Lout <- Lrand(pobs, plat, 1.5)
L_star <- Lout$L; z_star <- Lout$z # True latent component; True cluster labels

Sigma_star <- solve(S_star + L_star) # True Sigma

sest <- c("glasso", "clime", "gscad") # Different sparse estimators
mainlist <- list()
mainlist$glasso <- mainlist$clime <- mainlist$gscad <- vector("list", 100)
for(i in iters){
  set.seed(i*123)
  X <- mvrnorm(n, rep(0, pobs), Sigma = Sigma_star) # Finite sample data
  Sigma <- cov(X) # Sample Sigma
  
  # Loop through different estimators
  for(j in sest){
    print(j)
    cvout <- cv.slice(X, lambdas = logseq(1e-5, 0.05, 5), Sest = j)
    mainlist[[j]][[i]] <- slice(Sigma, cvout$lambda, cvout$r)
  }
  save(mainlist, file = "diffSest.rda")
}