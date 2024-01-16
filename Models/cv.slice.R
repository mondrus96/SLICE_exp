# For cross-validation
cv.slice = function(X, folds = 3, lambdas = c(1e-4, 1e-3, 1e-2, 0.1, 0.2, 0.4), rs = c(2:8), sest = "glasso"){
  # X = input data matrix
  # k = number of folds to perform for CV
  # lambdas = list of lambdas values to try
  # rs = list of rs to try
  # sest = sparse estimator - either glasso or clime
  
  n <- nrow(X) # number of samples
  cvmat <- matrix(NA, length(rs), length(lambdas)); rownames(cvmat) <- rs; colnames(cvmat) <- lambdas
  
  # Go over grid of lambdas and rs
  for(i in 1:length(rs)){
    print(paste0("rank: ", rs[i]))
    for(j in 1:length(lambdas)){
      print(paste0("lambda: ", lambdas[j]))
      
      ind <- sample(1:folds, n, replace = TRUE) # Define indices
      mulogL <- c()
      for(k in 1:folds){
        train <- X[ind != k,]; test <- X[ind == k,] # Train and test sets
        
        out <- slice(cov(train), lambdas[j], rs[i], sest) # Run method
        S <- out$S; L <- out$L
        
        logL <- logL(cov(test), S, L) # Append to mulogL
        mulogL <- c(mulogL, logL)
      }
      cvmat[i, j] <- mean(mulogL)
    }
  }
  best <- which(cvmat == max(cvmat), arr.ind = TRUE)
  
  return(list(cvmat = cvmat, maxlogL = max(cvmat), lambda = lambdas[best[1]], r = rs[best[2]]))
}
