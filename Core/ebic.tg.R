library(glasso)

ebic.tg = function(X, lambda = logseq(1e-5, 0.1, 5), tau = logseq(1e-5, 0.1, 5)){
  # X = input data matrix
  # lambda = vector of lambdas to try
  # tau = vector of thresholds to try
  
  p <- ncol(X); n <- nrow(X) # dimensions of X
  Sigma <- cov(X)
  ebicmat <- matrix(NA, length(tau), length(lambdas)); rownames(ebicmat) <- tau; colnames(ebicmat) <- lambdas
  
  # Go over grid of lambdas
  for(i in 1:length(lambdas)){
    print(paste0("lambda: ", lambdas[i]))
    for(j in 1:length(tau)){
      print(paste0("tau: ", tau[j]))
      
      # Fit glasso
      Theta <- glasso(Sigma, lambdas[i])$wi
      
      # Fit threshold
      S <- Theta
      S[abs(S) < tau[j]] <- 0

      # Calculate EBIC
      likl <- logL(Sigma, S, 0)
      k <- sum(S[upper.tri(S)] != 0)
      ebicmat[i, j] <- ebic(likl, p, n, k)
    } 
  }
  best <- which(ebicmat == min(ebicmat), arr.ind = TRUE)
  
  # Refit
  Theta <- glasso(Sigma, lambdas[best[2]])$wi
  S[abs(S) < tau[best[1]]] <- 0
  
  return(list(ebicvec = ebicmat, minebic = min(ebicmat), S = S,
              lambda = lambdas[best[2]], tau = tau[best[1]]))
}
