library(glasso)

ebic.tg = function(X, lambdas = logseq(1e-5, 1e-3, 10), taus = logseq(1e-5, 0.9, 10)){
  # X = input data matrix
  # lambda = vector of lambdas to try
  # tau = vector of thresholds to try
  
  p <- ncol(X); n <- nrow(X) # dimensions of X
  Sigma <- cov(X)
  ebicmat <- matrix(NA, length(taus), length(lambdas)); rownames(ebicmat) <- taus; colnames(ebicmat) <- lambdas
  
  # Go over grid of taus and lambdas
  for(i in 1:length(taus)){
    print(paste0("tau: ", taus[i]))
    for(j in 1:length(lambdas)){
      print(paste0("lambda: ", lambdas[j]))
      
      # Fit glasso
      gout <- glasso(Sigma, lambdas[j])
      
      # Fit threshold
      S <- gout$wi
      S[abs(S) <= taus[i]] <- 0

      # Calculate EBIC
      likl <- logL(Sigma, S, 0)
      k <- sum(S[upper.tri(S)] != 0)
      ebicmat[i, j] <- ebic(likl, p, n, k)
    } 
  }
  best <- which(ebicmat == min(ebicmat, na.rm=TRUE), arr.ind = TRUE)
  if(nrow(best) > 2){
    best <- best[1,]
  }
  
  # Refit
  Theta <- glasso(Sigma, lambdas[best[2]])$wi
  S[abs(S) < taus[best[1]]] <- 0
  
  return(list(ebicvec = ebicmat, minebic = min(ebicmat), S = S,
              lambda = lambdas[best[2]], tau = taus[best[1]]))
}
