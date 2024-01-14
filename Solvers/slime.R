slime = function(Sigma, lambda, rank, TOL = 1e-4, rand.init = FALSE){
  # Sigma = the input covariance matrix
  # lambda = regularization parameter for clime/graphical lasso
  # rank = rank
  
  # Make Sigma PD
  p = ncol(Sigma)
  Sigma = makePD(Sigma)
  invSigma = solve(Sigma)
  
  # Initial estimate of L
  if(rand.init == TRUE){
    L = matrix(rnorm(p * p, 0, 0.001), p, p)
  } else{
    L = matrix(0, p, p)
  }
  E = invSigma - L
  S = matrix(0, p, p)
  
  for(i in 1:100){
    print(i)
    # Make expectation PD
    if(!isPD(E)){
      E = makePD(E)
    }
    
    # Maximize
    Sold = S
    S = glasso(solve(E), lambda)$wi
    S = (S + t(S))/2
    
    # Define new L
    Lold = L
    L = invSigma - S
    eigL = eigen(L)
    L = eigL$vectors[,1:rank] %*% diag(eigL$values[1:rank]) %*% t(eigL$vectors[,1:rank])
    L = (L + t(L))/2
    
    # Define new expected value
    E = invSigma - L
    E = (E + t(E))/2
    
    # Convergence check
    if(norm(S - Sold) < TOL && (norm(L - Lold) < TOL)){
      break
    }
  }
  return(list(S = S, L = L))
}