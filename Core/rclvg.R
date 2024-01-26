library(corpcor)
library(glasso)

### Main rank constrained lvglasso function
rclvg <- function(Sigma, lambda, nLatents, tol = 1e-3, maxiter = 100){
  # To make life easier
  pobs <- ncol(Sigma)
  ptot <- pobs + nLatents
  O <- c(rep(TRUE, pobs), rep(FALSE, nLatents)); H <- !O
  
  # Initialize
  K <- matrix(0, ptot, ptot)
  K[O,O] <- solve(makePD(Sigma)); K[O,H] <- K[H,H] <- K[H,O] <- 1
  
  iter <- 1
  Kold <- K
  Sold <- K[O,O]; Lold <- K[O,H] %*% K[H,H] %*% K[H,O] # Define S and L
  
  for(iter in 1:maxiter){
    expS <- as.matrix(forceSymmetric(Estep(Sigma, K, O, H))) # E step
    K <- as.matrix(forceSymmetric(Mstep(expS, O, lambda))) # M step
    K <- forcePositive(K) 
    
    S <- K[O,O]; L <- K[O,H] %*% K[H,H] %*% K[H,O] # Define S and L
    
    # Check for convergence
    deltaK <- norm(K - Kold)
    deltalogL <- suppressWarnings(abs(logL(Sigma, S + L) - logL(Sigma, Sold + Lold)))
    if(deltaK < tol || ifelse(is.na(deltalogL < tol), FALSE, deltalogL < tol)){
      break
    } else {
      Kold <- K; Sold <- S; Lold <- L
    }
  }
  return(list(S = S, L = L))
}

# E-step in the optimization algorithm:
Estep <- function(Sigma, Kcur, O, H){
  # Current estimate of S
  Scur <- corpcor::pseudoinverse(Kcur)
  
  # Expected Sigma_OH
  Sigma_OH <- Sigma %*% corpcor::pseudoinverse(Scur[O,O]) %*% Scur[O, H]
  
  # Expected Sigma_H
  Sigma_H <- Scur[H, H] - Scur[H,O] %*% 
    corpcor::pseudoinverse(Scur[O,O]) %*% Scur[O,H] + Scur[H,O] %*% 
    corpcor::pseudoinverse(Scur[O,O]) %*% Sigma %*% 
    corpcor::pseudoinverse(Scur[O,O]) %*% Scur[O,H]
  
  # Construct expected Sigma
  expSigma <- rbind(cbind(Sigma,Sigma_OH),cbind(t(Sigma_OH), Sigma_H))  
  
  return(expSigma)
}

# M-step in the optimization algorithm:
Mstep <- function(expS, O, lambda){
  if (!isPD(expS)){
    expS <- as.matrix(forcePositive(expS))
  }
  # Rho matrix
  n <- nrow(expS)
  RhoMat <- matrix(lambda, n, n)
  RhoMat[!O,] <- 0
  RhoMat[,!O] <- 0
  
  K <- glasso(expS, RhoMat, penalize.diagonal=FALSE)$wi
  return(K)  
}

forcePositive <- function(x){
  if(any(eigen(x)$values<0)){
    cov2cor(x - (min(eigen(x)$values)-.1) * diag(nrow(x)))
  } else {
    x
  }
}

forceSymmetric <- function(x){
  return((x + t(x))/2)
}