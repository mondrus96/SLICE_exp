library(glasso)

ebic.tg = function(X, lambdas = NULL, taus = logseq(1e-5, 0.2, 20)){
  # X = input data matrix
  # lambda = vector of lambdas to try
  # tau = vector of thresholds to try
  
  p <- ncol(X); n <- nrow(X) # dimensions of X
  Sigma <- cov(X)
  
  o <- sqrt(log(p)/n)
  lambdas <- seq(o-(0.05*o), o+(0.05*o), length.out = 9)
  
  # Get best lambda
  ebicvec <- c()
  for(j in 1:length(lambdas)){
    print(paste0("lambda: ", lambdas[j]))
    
    # Fit glasso
    Theta <- glasso(Sigma, lambdas[j])$wi
    
    # Calculate EBIC
    likl <- logL(Sigma, Theta)
    k <- sum(Theta[upper.tri(Theta)] != 0)
    ebicvec <- c(ebicvec, ebic(likl, p, n, k))
  }
  lambda <- lambdas[which(ebicvec == min(ebicvec, na.rm=TRUE), arr.ind = TRUE)]
  if(length(lambda) > 2){
    lambda <- lambda[1]
  }
  
  # Get best tau
  ebicvec <- c()
  for(i in 1:length(taus)){
    print(paste0("tau: ", taus[i]))
    
    # Fit threshold
    Theta <- glasso(Sigma, lambda)$wi
    Theta[abs(Theta) <= taus[i]] <- 0
    
    # Calculate EBIC
    likl <- logL(Sigma, Theta)
    k <- sum(Theta[upper.tri(Theta)] != 0)
    ebicvec <- c(ebicvec, ebic(likl, p, n, k))
  }
  tau <- taus[which(ebicvec == min(ebicvec, na.rm=TRUE), arr.ind = TRUE)]
  if(length(tau) > 2){
    tau <- tau[1]
  }
  
  # Refit
  Theta <- glasso(Sigma, lambda)$wi
  Theta[abs(Theta) < tau] <- 0
  
  return(list(S = Theta, lambda = lambda, tau = tau))
}
