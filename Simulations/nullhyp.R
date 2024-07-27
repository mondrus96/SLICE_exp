library(MASS)
sapply((paste0("../Core/", list.files("../Core/"))), source)

# Loop through the two settings
testout <- c()
for(i in 1:2){
  for(j in 1:100){
    # Settings
    set.seed(j*123)
    pobs <- 150
    S_star <- Smat(pobs, 2, 1.5); S_star[S_star < 0.01] <- 0 # Sparse
    if(i == 1){
      L_star <- matrix(0, pobs, pobs) 
    } else if(i == 2){
      Lout <- Lrand(pobs, 4, 1.5) # 4 latent variables
      L_star <- Lout$L
    }
    Sigma <- Matrix::chol2inv(Matrix::chol(S_star + L_star)) # True Sigma
    n <- 10000
    X <- mvrnorm(n, rep(0, pobs), Sigma)
    
    # Run cross validation
    cvout <- cv.slice(X)
    sliout <- slice(cov(X), cvout$rho, cvout$r)
    likl_A <- logL(Sigma, sliout$S + sliout$L)
    likl_0 <- logL(Sigma, sliout$S)
    df <- sum(sliout$S[upper.tri(sliout$S)] != 0) + 
      (cvout$r*pobs-(cvout$r*(cvout$r-1)/2))
    testout <- rbind(testout, c(i, unlist(L.ratio.test(likl_0, likl_A, df))))
    print(testout)
  }
}
colnames(testout)[1] <- "setting"