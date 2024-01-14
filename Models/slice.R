# Load libraries
library(clime)
library(glasso)

# Main slice estimator code
slice = function(sigma, lambda, rank, sest = "glasso", tol = 1e-4, randinit = FALSE, maxiter = 1000){
  # sigma = the input covariance matrix
  # lambda = regularization parameter for clime/graphical lasso
  # rank = rank
  
  p <- ncol(sigma) # Make Sigma PD
  sigma <- makePD(sigma)
  invsigma <- solve(sigma)
  
  # Initial estimate of L
  if(randinit == TRUE){
    L <- matrix(rnorm(p * p, 0, 0.01), p, p)
  } else{
    L <- matrix(0, p, p)
  }
  E <- invsigma - L
  S <- matrix(0, p, p)
  
  for(i in 1:maxiter){
    print(i)
    # Make expectation PD
    if(!isPD(E)){
      E <- makePD(E)
    }
    
    Sold <- S # Sparse step
    if(sest == "glasso"){
      S <- glasso(solve(E), lambda)$wi 
    }
    S <- (S + t(S))/2
    
    Lold <- L # Latent step
    L <- invsigma - S 
    eigL <- eigen(L)
    L <- eigL$vectors[,1:rank] %*% diag(eigL$values[1:rank]) %*% t(eigL$vectors[,1:rank])
    L <- (L + t(L))/2
    
    E <- invsigma - L # Define new expectation
    E <- (E + t(E))/2
    
    # Convergence check
    if(norm(S - Sold) < tol && (norm(L - Lold) < tol)){
      break
    }
  }
  return(list(S = S, L = L))
}