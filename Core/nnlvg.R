library(RSpectra)

nnlvg <- function(Sigma, rho, gamma, lambda = 1, maxiter = 100, tol = 1e-3){
  p <- ncol(Sigma)
  
  R <- matrix(0, p, p) # for logdet + trace step
  S <- matrix(0, p, p) # sparse component
  L <- matrix(0, p, p) # latent component
  U <- matrix(0, p, p) # dual variable
  
  history <- list(objval = c(), r_norm = c(), s_norm = c(), eps_pri = c(), eps_dual = c())
  
  for (i in 1:maxiter){
    # R-update
    Rold <- R
    eigout <- eigen(lambda * (S - U) - Sigma)
    es <- eigout$values
    xi <- (es + sqrt(es^2 + 4*lambda)) / (2*lambda)
    R <- eigout$vectors %*% diag(xi) %*% t(eigout$vectors)
    R <- (R + t(R))/2
    
    # S-update
    Sold <- S
    S <- L1_shr(R + U, rho / lambda)
    S <- (S + t(S))/2
    
    # L-update
    Lold <- L
    L <- nucl_shr(R - S, gamma / lambda, tol)
    L <- (L + t(L))/2
    
    # Update dual variable
    U <- U + (R - S - L)
    
    # Diagnostics, reporting, termination checks
    history$objval <- c(history$objval, objective(Sigma, R, S, L, rho, gamma))
    history$r_norm <- c(history$r_norm, norm(R - S + L, type = "F"))
    history$s_norm <- c(history$s_norm, lambda*(norm((R - Rold), type = "F")))
    history$eps_pri <- c(history$eps_pri, sqrt(p*p) * tol + tol * max(norm(R, type = "F"), norm(S - L, type = "F")))
    history$eps_dual <- c(history$eps_dual, sqrt(p*p) * tol + tol * norm(lambda * U, type = "F"))
    
    if (history$r_norm[i] < history$eps_pri[i] && history$s_norm[i] < history$eps_dual[i]) {
      break
    }
  }
  
  list(S = S, L = L, history = history)
}

objective <- function(Sigma, R, S, L, rho, gamma){
  return(sum(diag(Sigma %*% R)) - log(det(R)) + rho * sum(abs(S)) + gamma * nucl_norm(L))
}

L1_shr <- function(a, kappa){
  return(sign(a) * pmax(0, abs(a) - kappa))
}

nucl_shr <- function(a, kappa, tol){
  a = eigen(a, symmetric = TRUE)
  return(a$vectors %*% diag(pmax(0, a$values - kappa)) %*% t(a$vectors))
}

nucl_norm <- function(mat){
  s <- eigen(mat)$values  # Get the singular values
  return(sum(abs(s)))  # Return the sum of the singular values
}
