admmlatentglasso <- function(Sigma, lambda, gamma, rho, MAX_ITER = 1000, ABSTOL = 1e-4, RELTOL = 1e-2) {
  p <- ncol(Sigma)
  
  R <- matrix(0, p, p)
  S <- matrix(0, p, p)
  L <- matrix(0, p, p)
  U <- matrix(0, p, p)
  
  history <- list(objval = rep(NA, MAX_ITER), r_norm = rep(NA, MAX_ITER), 
                  s_norm = rep(NA, MAX_ITER), eps_pri = rep(NA, MAX_ITER), eps_dual = rep(NA, MAX_ITER))
  
  for (k in 1:MAX_ITER) {
    
    # R-update
    Rold <- R
    eigout <- eigen(rho * (S - U) - Sigma)
    es <- eigout$values
    xi <- (es + sqrt(es^2 + 4*rho)) / (2*rho)
    R <- eigout$vectors %*% diag(xi) %*% t(eigout$vectors)
    R = (R + t(R))/2
    
    # S-update
    Sold <- S
    S <- L1_shr(R + U, lambda / rho)
    S = (S + t(S))/2
    
    # L-update ###
    Lold <- L 
    L <- nucl_shr(R - S, gamma / rho)
    L = (L + t(L))/2
    
    # Update dual variable
    U <- U + (R - S - L) ###
    
    # Diagnostics, reporting, termination checks
    history$objval[k] <- objective(Sigma, R, S, L, lambda, gamma) ###
    history$r_norm[k] <- norm(R - S + L, type = "F") ###
    history$s_norm[k] <- rho*(norm((R - Rold), type = "F")) ###
    history$eps_pri[k] <- sqrt(p*p) * ABSTOL + RELTOL * max(norm(R, type = "F"), norm(S - L, type = "F")) ###
    history$eps_dual[k] <- sqrt(p*p) * ABSTOL + RELTOL * norm(rho * U, type = "F")
    
    if (history$r_norm[k] < history$eps_pri[k] && history$s_norm[k] < history$eps_dual[k]) {
      break
    }
  }
  
  list(S = S, L = L, history = history)
}

objective <- function(Sigma, R, S, L, lambda, gamma) {
  return(sum(diag(Sigma %*% R)) - log(det(R)) + lambda * sum(abs(S)) + gamma * nucl_norm(L))
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
