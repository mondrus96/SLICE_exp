# For cross-validation
cv.lvg = function(X, folds = 3, lambdas = c(1e-4, 1e-3, 1e-2, 0.1), gammas = c(1e-4, 1e-3, 1e-2, 0.1)){
  # X = input data matrix
  # k = number of folds to perform for CV
  # lambdas = list of lambdas values to try
  # gammas = list of gamma values to try
  
  n <- nrow(X) # number of samples
  cvmat <- matrix(NA, length(gammas), length(lambdas)); rownames(cvmat) <- gammas; colnames(cvmat) <- lambdas
  
  # Go over grid of lambdas and gammas
  for(i in 1:length(gammas)){
    print(paste0("gamma: ", gammas[i]))
    for(j in 1:length(lambdas)){
      print(paste0("lambda: ", lambdas[j]))
      
      ind <- sample(1:folds, n, replace = TRUE) # Define indices
      mulogL <- c()
      for(k in 1:folds){
        train <- X[ind != k,]; test <- X[ind == k,] # Train and test sets
        
        out <- lvglasso(cov(train), lambdas[j], gammas[i]) # Run method
        S <- out$S; L <- out$L
        
        logL <- logL(cov(test), S, L) # Append to mulogL
        mulogL <- c(mulogL, logL)
      }
      cvmat[i, j] <- mean(mulogL)
    }
  }
  best <- which(cvmat == max(cvmat), arr.ind = TRUE)
  
  return(list(cvmat = cvmat, maxlogL = max(cvmat), lambda = lambdas[best[2]], gamma = gammas[best[1]]))
}
