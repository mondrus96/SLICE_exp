library(glasso)

# Function for calculating graphical scad
gscad <- function(Sigma, rho, a = 3.7){
  # Sigma = input covariance matrix
  # rho = regularization parameter
  # a = concavity parameter, 3.7 by default
  
  # Initial L1 estimate
  L1est <- glasso(Sigma,rho)$wi
  
  # Generate rhomat
  p <- nrow(Sigma)
  rhomat <- scadrightderv(abs(L1est), a, rho)
  SCADest <- glasso(Sigma, rhomat)$wi
  return(SCADest)
}

# Helper function
scadrightderv <- function(lamhat, a, rho){
  return(pmax(rho*((lamhat<=rho)+pmax(a*rho-lamhat,0)*(lamhat>rho)/(a-1)/rho),1e-4))
}