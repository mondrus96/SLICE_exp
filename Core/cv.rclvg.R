library(pracma)

cv.rclvg = function(X, folds = 3, rhos = logseq(1e-5, 0.1, 5), rs = c(2:6)){
  # X = input data matrix, or Sigma
  # k = number of folds to perform for CV
  # rhos = list of rhos values to try
  # rs = list of rs to try
  # sest = sparse estimator - either glasso or clime
  
  p <- ncol(X); n <- nrow(X) # dimensions of X
  cvmat <- matrix(NA, length(rs), length(rhos)); rownames(cvmat) <- rs; colnames(cvmat) <- rhos
  
  # Go over grid of rhos and rs
  for(i in 1:length(rs)){
    print(paste0("rank: ", rs[i]))
    for(j in 1:length(rhos)){
      print(paste0("rho: ", rhos[j]))
      
      ind <- sample(1:folds, n, replace = TRUE) # Define indices
      mulogL <- c()
      for(k in 1:folds){
        train <- X[ind != k,]; test <- X[ind == k,] # Train and test sets
        
        out <- rclvg(cov(train), rhos[j], rs[i]) # Run method
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
              rho = rhos[best[2]], r = rs[best[1]]))
}
