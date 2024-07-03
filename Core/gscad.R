library(glasso)

# Function for calculating graphical scad
gscad <- function(Sigma, lambda, a = 3.7){
  # Sigma = input covariance matrix
  # lambda = regularization parameter
  # a = concavity parameter, 3.7 by default
  
  # Initial L1 estimate
  L1est <- glasso(Sigma,lambda)$wi
  
  # Generate rhomat
  p <- nrow(Sigma)
  lambdamat <- scadrightderv(abs(L1est), a, lambda)
  SCADest <- glasso(Sigma, lambdamat)$wi
  return(SCADest)
}

# Helper function
scadrightderv <- function(lamhat, a, lambda){
  return(pmax(lambda*((lamhat<=lambda)+pmax(a*lambda-lamhat,0)*(lamhat>lambda)/(a-1)/lambda),1e-4))
}