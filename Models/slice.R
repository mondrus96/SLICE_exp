library(clime)
library(glasso)

# Main slice estimator
slice <- function(Sigma, lambda, rank, sest = "glasso", tol = 1e-4, randinit = FALSE, maxiter = 100){
  # Sigma = the input covariance matrix
  # lambda = regularization parameter for clime/graphical lasso
  # rank = rank
  # sest = sparse estimator - glasso or clime
  
  p <- ncol(Sigma) # Make Sigma PD
  Sigma <- makePD(Sigma)
  invSigma <- solve(Sigma)
  
  # Initial estimate of L
  if(randinit == TRUE){
    L <- matrix(rnorm(p * p, 0, 0.01), p, p)
  } else{
    L <- matrix(0, p, p)
  }
  E <- invSigma - L
  S <- matrix(0, p, p)
  
  for(i in 1:maxiter){
    # Make expectation PD
    if(!isPD(E)){
      E <- makePD(E)
    }
    
    Sold <- S # Sparse step
    if(sest == "glasso"){
      S <- glasso(solve(E), lambda)$wi 
    } else if(sest == "clime"){
      S <- clime(solve(E), lambda, Sigma = TRUE)$Omegalist[[1]]
    }
    S <- (S + t(S))/2
    
    Lold <- L # Latent step
    L <- invSigma - S 
    eigL <- eigen(L)
    L <- eigL$vectors[,1:rank] %*% diag(eigL$values[1:rank]) %*% t(eigL$vectors[,1:rank])
    L <- (L + t(L))/2
    
    E <- invSigma - L # Define new expectation
    E <- (E + t(E))/2
    
    # Convergence check
    paramdiff = (norm(S - Sold) < tol) && (norm(L - Lold) < tol)
    objdiff = norm(invSigma - (S + L)) < tol
    if(paramdiff || objdiff){
      break
    }
  }
  return(list(S = S, L = L))
}