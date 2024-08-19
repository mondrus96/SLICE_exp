library(ggplot2)
sapply((paste0("../Core/", list.files("../Core/"))), source)

# Define simulation
set.seed(123)
p <- 100
L_star <- Lspir(p)$L
S_star <- Smat(p, 2, 1.5)
S_star[S_star < 0.01] <- 0 # True sparse component
Sigma_star <- solve(S_star + L_star)

# Get first two eigenvalues
eigL_star <- eigen(L_star)

# Define the function to calculate the distance from the true point
f <- function(x, y, Sigma_star, S_star, eigL_star){
  vecsL_star <- eigL_star$vectors[,1:2]
  likl <- logL(Sigma_star, 
               S_star + 
                 (vecsL_star %*% diag(c(x, y)) %*% t(vecsL_star)))
  return(likl)
}

# Function for creating arrows df
arrowdf <- function(arrow_data){
  df <- data.frame(
    x = arrow_data[,1][-length(arrow_data[,1])],  # Starting x points, exclude the last point
    y = arrow_data[,2][-length(arrow_data[,2])],  # Starting y points, exclude the last point
    xend = arrow_data[,1][-1],  # Ending x points, exclude the first point
    yend = arrow_data[,2][-1])  # Ending y points, exclude the first point
  return(df)
}

# Generate grid
x <- seq(0, 100, length.out = 100)
y <- seq(0, 15, length.out = 100)
grid <- expand.grid(x = x, y = y)
grid$z <- mapply(f, grid$x, grid$y, MoreArgs = list(Sigma_star = Sigma_star, 
                                        S_star = S_star, 
                                        eigL_star = eigL_star))

# Get data for arrows
sli_eigs <- sliceLs(Sigma_star, 0.1, 2)
nnlvg_eigs <- nnlvgLs(Sigma_star, 0.01, 0.01)
rclvg_eigs <- rclvgLs(Sigma_star, 0.1, 2)
arrows_df1 <- arrowdf(sli_eigs)
arrows_df2 <- arrowdf(nnlvg_eigs)
arrows_df3 <- arrowdf(rclvg_eigs)

# Start plotting
plt <- ggplot(grid, aes(x = x, y = y)) +  # Define only x and y in the global aes
  geom_contour(aes(z = z), color = "black", bins = 20) +  # Contour lines
  labs(x = expression(lambda[1]), y = expression(lambda[2])) +
  coord_fixed(5) +  # Keep aspect ratio of x and y the same
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())  # Remove minor gridlines

# Add arrows to the plot
plt <- plt + 
  geom_segment(data = arrows_df1, 
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(type = "closed", length = unit(2, "mm"), angle = 15),  # Adjust angle and length
               color = "#0C6291", alpha = 0.5) + 
  geom_segment(data = arrows_df2, 
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(type = "closed", length = unit(2, "mm"), angle = 15),  # Adjust angle and length
               color = "#A63446", alpha = 0.5) +
  geom_segment(data = arrows_df3, 
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(type = "closed", length = unit(2, "mm"), angle = 15),  # Adjust angle and length
               color = "#34A646", alpha = 0.5)
ggsave("grad.png", plot = plt, width = 5, height = 4, dpi = 600)  

# Print the plot
print(plt)
print(nrow(arrows_df1))
print(nrow(arrows_df2))
print(nrow(arrows_df3))

### HELPER FUNCTIONS ###
# Slice Ls
sliceLs <- function(Sigma, lambda, rank, tol = 1e-3, maxiter = 100){
  
  p <- ncol(Sigma) # Make Sigma PD
  Sigma <- makePD(Sigma)
  invSigma <- Matrix::chol2inv(Matrix::chol(Sigma))
  
  L <- matrix(rnorm(p^2), p, p)
  L <- L %*% t(L)
  E <- invSigma - L # Expectation
  S <- 0 # Empty S
  alleigs <- c()
  
  deltaS <- deltaL <- deltalogL <- c() 
  for(i in 1:maxiter){
    if(!isPD(E)){
      E <- makePD(E) # Make expectation PD
    }
    
    # Sparse step
    Sold <- S
    S <- huge(Matrix::chol2inv(Matrix::chol(E)), lambda, method = "glasso", verbose = FALSE)$icov[[1]]
    
    # Latent step
    Lold <- L
    tsvdL <- svds(invSigma - S, rank)
    L <- tsvdL$v %*% diag(tsvdL$d) %*% t(tsvdL$v)
    L <- (L + t(L))/2
    
    # New expectation
    E <- invSigma - L
    E <- (E + t(E))/2
    
    deltaS <- c(deltaS, sqrt(sum((S - Sold)^2))); deltaL <- c(deltaL, sqrt(sum((L - Lold)^2))) # Convergence check
    if(i == 1){
      deltalogL <- c(deltalogL, suppressWarnings(abs(logL(Sigma, S + L) - logL(Sigma, Diagonal(p, x = tol)))))
    } else{
      deltalogL <- c(deltalogL, suppressWarnings(abs(logL(Sigma, S + L) - logL(Sigma, Sold + Lold)))) 
    }
    if((deltaS[i] < tol) && (deltaL[i] < tol) || 
       ifelse(is.na(deltalogL[i] < tol), FALSE, deltalogL[i] < tol)){
      break
    }
    # Get eigendecomposition
    eigL <- eigen(L)
    alleigs <- rbind(alleigs, eigL$values[1:2])
  }
  return(alleigs)
}

# nnLVGLASSO Ls
nnlvgLs <- function(Sigma, lambda, gamma, rho = 1, maxiter = 100, tol = 1e-3){
  p <- ncol(Sigma)
  
  R <- matrix(0, p, p) # for logdet + trace step
  S <- matrix(0, p, p) # sparse component
  L <- matrix(0, p, p) # latent component
  U <- matrix(0, p, p) # dual variable
  alleigs <- c()
  
  history <- list(objval = c(), r_norm = c(), s_norm = c(), eps_pri = c(), eps_dual = c())
  
  for (i in 1:maxiter){
    # R-update
    Rold <- R
    eigout <- eigen(rho * (S - U) - Sigma)
    es <- eigout$values
    xi <- (es + sqrt(es^2 + 4*rho)) / (2*rho)
    R <- eigout$vectors %*% diag(xi) %*% t(eigout$vectors)
    R <- (R + t(R))/2
    
    # S-update
    Sold <- S
    S <- L1_shr(R + U, lambda / rho)
    S <- (S + t(S))/2
    
    # L-update
    Lold <- L
    L <- nucl_shr(R - S, gamma / rho, tol)
    L <- (L + t(L))/2
    
    # Update dual variable
    U <- U + (R - S - L)
    
    # Diagnostics, reporting, termination checks
    history$objval <- c(history$objval, objective(Sigma, R, S, L, lambda, gamma))
    history$r_norm <- c(history$r_norm, norm(R - S + L, type = "F"))
    history$s_norm <- c(history$s_norm, rho*(norm((R - Rold), type = "F")))
    history$eps_pri <- c(history$eps_pri, sqrt(p*p) * tol + tol * max(norm(R, type = "F"), norm(S - L, type = "F")))
    history$eps_dual <- c(history$eps_dual, sqrt(p*p) * tol + tol * norm(rho * U, type = "F"))
    
    if (history$r_norm[i] < history$eps_pri[i] && history$s_norm[i] < history$eps_dual[i]) {
      break
    }
    # Get eigendecomposition
    eigL <- eigen(L)
    alleigs <- rbind(alleigs, eigL$values[1:2])
  }
  return(alleigs)
}

### Main rank constrained lvglasso function
rclvgLs <- function(Sigma, lambda, nLatents, tol = 1e-3, maxiter = 100){
  # To make life easier
  pobs <- ncol(Sigma)
  ptot <- pobs + nLatents
  O <- c(rep(TRUE, pobs), rep(FALSE, nLatents)); H <- !O
  
  # Initialize
  K <- matrix(rnorm(ptot*ptot), ptot, ptot)
  K[O,O] <- solve(makePD(Sigma)); K[O,H] <- K[H,H] <- K[H,O] <- 1
  
  iter <- 1
  Kold <- K
  Sold <- K[O,O]; Lold <- K[O,H] %*% K[H,H] %*% K[H,O] # Define S and L
  alleigs <- c()
  
  for(iter in 1:maxiter){
    if(!isPD(K)){
      K <- makePD(K)
    }
    expS <- as.matrix(forceSymmetric(Estep(Sigma, K, O, H))) # E step
    K <- as.matrix(forceSymmetric(Mstep(expS, O, lambda))) # M step
    
    S <- K[O,O]; L <- K[O,H] %*% K[H,H] %*% K[H,O] # Define S and L
    
    # Check for convergence
    deltaK <- norm(K - Kold)
    deltalogL <- suppressWarnings(abs(logL(Sigma, S + L) - logL(Sigma, Sold + Lold)))
    if(deltaK < tol || ifelse(is.na(deltalogL < tol), FALSE, deltalogL < tol)){
      break
    } else {
      Kold <- K; Sold <- S; Lold <- L
    }
    eigL <- eigen(L)
    alleigs <- rbind(alleigs, eigL$values[1:2])
  }
  return(alleigs)
}
