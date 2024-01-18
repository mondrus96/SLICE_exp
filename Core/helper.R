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
  eigvals = suppressWarnings(eigs(mat, ncol(mat), opts = list(retvec = FALSE))$values)
  perturb = max(max(eigvals) - p*min(eigvals), 0)/(p-1)
  mat = mat+diag(p)*perturb
  return(mat)
}

# Check if a matrix is PD through Cholesky decomposition (faster than full eigendecomp)
isPD = function(mat){
  tryCatch({
    chol(mat)
    return(TRUE)
  }, error = function(e){
    return(FALSE)
  })
}

# Get sin angle between two vectors
sintheta <- function(v, v_hat){
  return(sqrt(1 - sum(v_hat * v)^2))
}

# Calculates the mutual information score
nmi <- function(labels_true, labels_pred) {
  # Adapted from sklearn \url{https://scikit-learn.org/stable/modules/generated/sklearn.metrics.normalized_mutual_info_score.html}
  
  contingency <- table(labels_true, labels_pred) # Create contingency table
  contingency_sum <- sum(contingency) # Total sum of the contingency table
  row_sums <- rowSums(contingency) # Row and column sums
  col_sums <- colSums(contingency)
  
  mi <- 0 # Initialize mutual information
  for (i in seq_len(nrow(contingency))) {
    for (j in seq_len(ncol(contingency))) {
      if (contingency[i, j] > 0) {
        log_contingency_nm <- log(contingency[i, j])
        contingency_nm <- contingency[i, j] / contingency_sum
        outer <- row_sums[i] * col_sums[j]
        log_outer <- -log(outer) + log(sum(row_sums)) + log(sum(col_sums))
        mi_contrib <- contingency_nm * (log_contingency_nm - log(contingency_sum)) + contingency_nm * log_outer
        mi <- mi + mi_contrib
      }
    }
  }
  
  if (mi == 0) { # No need to proceed if mi == 0
    return(0)
  }
  
  h_true <- entropy(labels_true); h_pred <- entropy(labels_pred) # Calculate entropies
  norm <- mean(h_true, h_pred) # Normalization
  return(as.numeric(mi/norm))
}

# Calculates the adjusted rand index
ari <- function(labels_true, labels_pred){
  # Adapted from sklearn \url{https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_rand_score.html}
  
  contingency <- pair_confusion_matrix(labels_true, labels_pred)
  num <- sum(diag(contingency))
  den <- sum(contingency)
  
  if (num == den || den == 0){ 
    return(1)
  } else {
    return(num/den)
  }
}

# Calculate cluster entropy
entropy <- function(labels) {
  # Adapted from sklearn \url{https://github.com/scikit-learn/scikit-learn/blob/main/sklearn/metrics/cluster/_supervised.py#L1261}
  
  if (length(labels) == 0) {
    return(1)
  }
  
  label_counts <- table(labels) # Count the occurrence of each label
  if (length(label_counts) == 1) { # If there's only one label, entropy is 0
    return(0)
  }
  
  total_count <- sum(label_counts) # Calculate probabilities
  probs <- label_counts / total_count
  entropy <- -sum(probs * log(probs)) # Calculate entropy
  return(entropy)
}

# Matching labels through contingency matrix
pair_confusion_matrix <- function(labels_true, labels_pred) {
  # Adapted from sklearn \url{https://github.com/scikit-learn/scikit-learn/blob/main/sklearn/metrics/cluster/_supervised.py#L190}
  
  labels_true <- factor(labels_true) # Ensure labels are factors and have the same levels
  labels_pred <- factor(labels_pred, levels = levels(labels_true))
  
  contingency <- table(labels_true, labels_pred) # Compute contingency matrix
  
  n_samples <- sum(contingency)
  n_c <- rowSums(contingency)
  n_k <- colSums(contingency)
  contingency_values <- as.vector(contingency) # Convert contingency matrix to vector of its values
  sum_squares <- sum(contingency_values^2) # Compute sum of squares of contingency matrix
  
  cmat <- matrix(numeric(4), nrow = 2) # Initialize and compute the pair confusion matrix
  cmat[2, 2] <- sum_squares - n_samples
  cmat[1, 2] <- sum(contingency %*% n_k) - sum_squares
  cmat[2, 1] <- sum(t(contingency) %*% n_c) - sum_squares
  cmat[1, 1] <- n_samples^2 - cmat[1, 2] - cmat[2, 1] - sum_squares
  
  return(cmat) # Return contingency table
}