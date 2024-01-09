admmlatentglasso <- function(S, lambda, gamma, rho, MAX_ITER = 1000, ABSTOL = 1e-4, RELTOL = 1e-2) {
  n <- ncol(S)
  
  X <- matrix(0, n, n)
  Z <- matrix(0, n, n)
  U <- matrix(0, n, n)
  L <- matrix(0, n, n) ###
  
  history <- list(objval = rep(NA, MAX_ITER), r_norm = rep(NA, MAX_ITER), 
                  s_norm = rep(NA, MAX_ITER), eps_pri = rep(NA, MAX_ITER), eps_dual = rep(NA, MAX_ITER))
  
  for (k in 1:MAX_ITER) {
    
    # X-update
    Xold <- X
    eigout <- eigen(rho * (Z - U) - S)
    es <- eigout$values
    xi <- (es + sqrt(es^2 + 4*rho)) / (2*rho)
    X <- eigout$vectors %*% diag(xi) %*% t(eigout$vectors)
    X = (X + t(X))/2
    
    # Z-update
    Zold <- Z
    Z <- L1_shr(X + U, lambda / rho)
    Z = (Z + t(Z))/2
    
    # L-update ###
    Lold <- L 
    L <- nucl_shr(X - Z, gamma / rho)
    L = (L + t(L))/2
    
    # Update dual variable
    U <- U + (X - Z - L) ###
    
    # Diagnostics, reporting, termination checks
    history$objval[k] <- objective(S, X, Z, L, lambda, gamma) ###
    history$r_norm[k] <- norm(X - Z + L, type = "F") ###
    history$s_norm[k] <- rho*(norm((X - Xold), type = "F")) ###
    history$eps_pri[k] <- sqrt(p*p) * ABSTOL + RELTOL * max(norm(X, type = "F"), norm(Z - L, type = "F")) ###
    history$eps_dual[k] <- sqrt(p*p) * ABSTOL + RELTOL * norm(rho * U, type = "F")
    
    #cat(sprintf("%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n", k, 
    #            history$r_norm[k], history$eps_pri[k], history$s_norm[k], 
    #            history$eps_dual[k], history$objval[k]))
    
    if (history$r_norm[k] < history$eps_pri[k] && history$s_norm[k] < history$eps_dual[k]) {
      break
    }
  }
  
  list(Z = Z, L = L, history = history)
}

objective <- function(S, X, Z, L, lambda, gamma) {
  return(sum(diag(S %*% X)) - log(det(X)) + lambda * sum(abs(Z)) + gamma * nucl_norm(L))
}

L1_shr <- function(a, kappa) {
  return(sign(a) * pmax(0, abs(a) - kappa))
}

nucl_shr <- function(a, kappa) {
  b = eigen(a)
  return(b$vectors %*% diag(pmax(0, b$values - kappa)) %*% t(b$vectors))
}

nucl_norm <- function(mat) {
  s <- eigen(mat)$values  # Get the singular values
  return(sum(abs(s)))  # Return the sum of the singular values
}
