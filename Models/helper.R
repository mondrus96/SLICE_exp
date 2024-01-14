# Extended Bayes Information Criterion
ebic = function(likl, p, n, k, gamma = 0.5){
  return(-2*likl + k*log(n) + 2*gamma*log(p/k))
}

# Bayes Information Criterion
bic = function(likl, n, k){
  return(-2*likl + k*log(n))
}

# Log likelihood function
logL = function(Sigma, S, L){
  return(log(det(S + L)) - sum(diag(Sigma %*% (S + L))))
}

# Make a matrix positive definite by adding a small value to diagonal
makePD = function(mat){
  p = ncol(mat)
  eigvals = eigen(mat, only.values=T)$values
  perturb = max(max(eigvals) - p*min(eigvals), 0)/(p-1)
  mat = mat+diag(p)*perturb
  return(mat)
}

# Check if a matrix is PD through Cholesky decomposition
isPD = function(mat){
  tryCatch({
    chol(mat)
    return(TRUE)
  }, error = function(e){
    return(FALSE)
  })
}