library(glasso)

ebic.tg = function(X, rhos = logseq(1e-5, 0.01, 5), taus = logseq(1e-5, 0.2, 5)){
  # X = input data matrix
  # rho = vector of rhos to try
  # tau = vector of thresholds to try
  
  p <- ncol(X); n <- nrow(X) # dimensions of X
  Sigma <- cov(X)
  
  #o <- sqrt(log(p)/n)
  #rhos <- seq(o-(0.05*o), o+(0.05*o), length.out = 9)
  
  # Get best rho
  ebicvec <- c()
  for(j in 1:length(rhos)){
    print(paste0("rho: ", rhos[j]))
    
    # Fit glasso
    Theta <- glasso(Sigma, rhos[j])$wi
    
    # Calculate EBIC
    likl <- logL(Sigma, Theta)
    k <- sum(Theta[upper.tri(Theta)] != 0)
    ebicvec <- c(ebicvec, ebic(likl, p, n, k))
  }
  rho <- rhos[which(ebicvec == min(ebicvec, na.rm=TRUE), arr.ind = TRUE)]
  if(length(rho) > 2){
    rho <- rho[1]
  }
  
  # Get best tau
  ebicvec <- c()
  for(i in 1:length(taus)){
    print(paste0("tau: ", taus[i]))
    
    # Fit threshold
    Theta <- glasso(Sigma, rho)$wi
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
  Theta <- glasso(Sigma, rho)$wi
  Theta[abs(Theta) < tau] <- 0
  
  return(list(S = Theta, rho = rho, tau = tau))
}
