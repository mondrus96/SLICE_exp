# This function is a simulation type where there is an exponential decay from the main diagonal
expdecay = function(p, decay, init_val){
  # Create an empty matrix
  Theta = matrix(0, p, p)
  
  # Fill the matrix with exponential decay values
  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j) {
        # Main diagonal element
        Theta[i, j] = init_val
      } else {
        # Off-diagonal elements decay exponentially
        distance = abs(i - j)
        Theta[i, j] = init_val*exp(-decay*distance)
      }
    }
  }
  
  # Permute the matrix
  perm = sample(p)
  Theta = Theta[perm, perm]
  
  # Return the output
  return(Theta)
}

# This function is for simulating community data
Lclust = function(p, r, init_val, sd){
  # Define Z
  Z = matrix(0, p, r)
  for(i in 1:p){
    Z[i,round(runif(1, 1, r))] = 1
  }
  
  # Test w/ L
  L = Z %*% t(Z)
  L = L * init_val
  
  # Return output
  return(list(L = L, Z = Z))
}