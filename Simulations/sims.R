# This function is a simulation type where there is an exponential decay from the main diagonal
Smat <- function(p, decay, init_val){
  # Fill the matrix with exponential decay values
  S <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j) {
        # Main diagonal element
        S[i, j] <- init_val
      } else {
        # Off-diagonal elements decay exponentially
        distance <- abs(i - j)
        S[i, j] <- init_val*exp(-decay*distance)
      }
    }
  }
  perm <- sample(p) # Permute the matrix
  S <- S[perm, perm]
  
  return(S)
}

# This function is for simulating uniform community structure
Lunif <- function(p, r, init_val, sd = 0){
  Z <- matrix(0, p, r) # Define Z
  for(i in 1:p){
    Z[i,sample(1:r, 1)] <- 1
  }
  z <- apply(Z, 1, function(row) which(row == 1)) # true clusters
  
  L <- Z %*% t(Z) 
  L <- L * init_val

  return(list(L = L, z = z))
}

# This function is for simulating exponential community structure
Lexp <- function(p, r, init_val){
  probs <- exp(-1*(0:(r - 1))) # Get probabilities
  probs <- probs/sum(probs)
  
  Z <- matrix(0, p, r)
  for(i in 1:p){
    Z[i,sample(1:r, 1, prob = probs)] <- 1
  }
  z <- apply(Z, 1, function(row) which(row == 1))
  
  L <- Z %*% t(Z)
  L <- L * init_val
  
  return(list(L = L, z = z))
}

# For creating crescent shape
Lcres <- function(p, sd = 0.1){
  theta_out <- seq(0, pi, length.out = p %/% 2) # Outer circle
  outer_circ_x <- cos(theta_out)
  outer_circ_y <- sin(theta_out)
  
  theta_in <- seq(0, pi, length.out = p - (p %/% 2)) # Inner circle
  inner_circ_x <- 1 - cos(theta_in)
  inner_circ_y <- 1 - sin(theta_in) - 0.5
  
  X <- rbind(cbind(outer_circ_x, outer_circ_y), cbind(inner_circ_x, inner_circ_y)) # Combine coordinates
  colnames(X) = c("x", "y")
  z <- c(rep(0, p/2), rep(1, p/2)) # Cluster membership
  
  X <- X + matrix(rnorm(prod(dim(X)), mean = 0, sd = sd), nrow = nrow(X), ncol = ncol(X))
  
  indices <- sample(1:nrow(X)) # Permute indices
  X <- X[indices,]
  z <- z[indices]
  
  L = X %*% t(X) # Generate L
  
  return(list(X = X, L = L, z = z))
}

# For creating circles, one within the other
Lcirc <- function(p, sd = 0.02){
  linspace_out <- seq(0, 2 * pi, length.out = p %/% 2)
  outer_circ_x <- cos(linspace_out)
  outer_circ_y <- sin(linspace_out)

  linspace_in <- seq(0, 2 * pi, length.out = p - (p %/% 2))
  inner_circ_x <- cos(linspace_in) * 0.8
  inner_circ_y <- sin(linspace_in) * 0.8
  
  X <- rbind(cbind(outer_circ_x, outer_circ_y), cbind(inner_circ_x, inner_circ_y))
  colnames(X) = c("x", "y")
  z <- c(rep(0, p %/% 2), rep(1, p - (p %/% 2)))
  
  X <- X + matrix(rnorm(prod(dim(X)), mean = 0, sd = sd), nrow = nrow(X), ncol = ncol(X))
  
  indices <- sample(1:nrow(X)) # Permute indices
  X <- X[indices,]
  z <- z[indices]
  
  L = X %*% t(X) # Generate L
  
  return(list(X = X, L = L, z = z))
}