library(pracma)

cv.nnlvg = function(X, folds = 3, rhos = logseq(1e-5, 0.1, 5), gammas = logseq(1e-5, 0.1, 5)){
  # X = input data matrix
  # k = number of folds to perform for CV
  # rhos = list of rhos values to try
  # gammas = list of gamma values to try
  
  n <- nrow(X) # number of samples
  cvmat <- matrix(NA, length(gammas), length(rhos)); rownames(cvmat) <- gammas; colnames(cvmat) <- rhos
  
  # Go over grid of rhos and gammas
  for(i in 1:length(gammas)){
    print(paste0("gamma: ", gammas[i]))
    for(j in 1:length(rhos)){
      print(paste0("rho: ", rhos[j]))
      
      ind <- sample(1:folds, n, replace = TRUE) # Define indices
      mulogL <- c()
      for(k in 1:folds){
        train <- X[ind != k,]; test <- X[ind == k,] # Train and test sets
        
        out <- nnlvg(cov(train), rhos[j], gammas[i]) # Run method
        S <- out$S; L <- out$L
        
        likl <- logL(cov(test), S + L) # Append to mulogL
        mulogL <- c(mulogL, likl)
      }
      cvmat[i, j] <- mean(mulogL)
    }
  }
  best <- which(cvmat == max(cvmat, na.rm=TRUE), arr.ind = TRUE)
  if(nrow(best) > 2){
    best <- best[1,]
  }
  
  return(list(cvmat = cvmat, maxlogL = max(cvmat), 
              rho = rhos[best[2]], gamma = gammas[best[1]]))
}