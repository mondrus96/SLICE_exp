library(clime)
library(glasso)
library(RSpectra)

# Main slice estimator
slice <- function(Sigma, lambda, rank, sest = "glasso",
                  tol = 1e-3, maxiter = 100, randinit = FALSE){
  # Sigma = the input covariance matrix
  # lambda = regularization parameter for clime/graphical lasso
  # rank = rank
  # sest = sparse estimator - glasso or clime
  
  if(!sest %in% c("glasso", "clime")){
    stop(paste(sest, "is not a valid sparse model"))
  }
  
  p <- ncol(Sigma) # Make Sigma PD
  Sigma <- makePD(Sigma)
  invSigma <- solve(Sigma)
  
  if(randinit == TRUE){
    L <- matrix(rnorm(p * p, 0, 0.01), p, p) # random initialization L
  } else if(randinit == FALSE){
    L <- matrix(0, p, p) # zero initialization L
  }
  E <- invSigma - L # Expectation
  S <- matrix(0, p, p) # Empty S
  
  for(i in 1:maxiter){
    if(!isPD(E)){
      E <- makePD(E) # Make expectation PD
    }
    
    Sold <- S # Sparse step
    if(sest == "glasso"){
      S <- glasso(solve(E), lambda, thr = tol, maxit = maxiter)$wi  
    } else if(sest == "clime"){
      S <- clime(solve(E), lambda, sigma = TRUE)$Omegalist[[1]]
    } else if(sest){
      
    }
    S <- (S + t(S))/2
    
    Lold <- L # Latent step
    L <- invSigma - S 
    tsvdL <- svds(L, rank)
    L <- tsvdL$u %*% diag(tsvdL$d) %*% t(tsvdL$v)
    L <- (L + t(L))/2
    
    E <- invSigma - L # Define new expectation
    E <- (E + t(E))/2
    
    # Convergence check
    deltaS = mean(abs(S - Sold)); deltaL = mean(abs(L - Lold))
    deltalogL = suppressWarnings(abs(logL(Sigma, S, L) - logL(Sigma, Sold, Lold)))
    if((deltaS < tol) && (deltaL < tol) || 
       ifelse(is.na(deltalogL < tol), FALSE, deltalogL < tol)){
      break
    }
  }
  return(list(S = S, L = L, lambda = lambda, rank = rank, 
              misc = list(converged = ifelse(i != maxiter, TRUE, FALSE), iters = i,
              deltaS = deltaS, deltaL = deltaL, deltalogL = deltalogL)))
}