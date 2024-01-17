library(MASS)
sapply((paste0("../Models/", list.files("../Models/"))), source)
sapply((paste0("../Simulations/", list.files("../Simulations/"))), source)

set.seed(123)
pobs <- 200 # Number of observed variables for S
plat <- 4 # Number of latent variables for L
n <- 10000

S_star <- Smat(pobs, 2, 1.5)
S_star[S_star < 0.01] = 0
Lout <- Lrand(pobs, plat, 1.5)
L_star <- Lout$L; z_star <- Lout$z
Sigma_star <- solve(S_star + L_star)
X <- mvrnorm(n, rep(0, pobs), Sigma_star)
Sigma = cov(X)

#lambdas = c(1e-4, 1e-3, 1e-2, 0.1, 0.2)
#rs = 2:8

CVout <- cv.slice(X) # Cross validation
out <- slice(cov(X), CVout$lambda, CVout$r) # Run method

heatmap(1*(out$S != 0), Rowv = NA, Colv = NA)
heatmap(1*(S_star != 0), Rowv = NA, Colv = NA)

heatmap(out$L, Rowv = NA, Colv = NA)
heatmap(L_star, Rowv = NA, Colv = NA)

ebicmat <- matrix(NA, length(rs), length(lambdas))
rownames(ebicmat) <- rs; colnames(ebicmat) <- lambdas
for(i in 1:length(rs)){
  print(paste0("rank: ", rs[i]))
  for(j in 1:length(lambdas)){
    print(paste0("lambda: ", lambdas[j]))
    
    out <- slice(Sigma_star, lambdas[j], rs[i]) # Run method
    likl <- logL(Sigma_star, out$S, out$L)
    
    k_S <- sum(out$S[upper.tri(out$S)] != 0); k_L <- rs[i] + (pobs*rs[i]) - (rs[i]*(rs[i] - 1)/2)
    ebicmat[i, j] <- ebic(likl, pobs, n, k_S + k_L)
  }
}
best <- which(ebicmat == min(ebicmat), arr.ind = TRUE)

plot_spiral <- function(turns, point_density = 100) {
  # The total number of points is turns multiplied by point density
  total_points <- turns * point_density
  
  # Create a sequence of angles for the spiral
  theta <- seq(0, 2 * pi * turns, length.out = total_points)
  
  # The radius grows as the angle grows
  radius <- seq(0, turns, length.out = total_points)
  
  # Convert polar coordinates (radius, theta) to Cartesian coordinates (x, y)
  X <- 0.3*cbind(radius * cos(theta), radius * sin(theta))
  
  # Plot the spiral
  plot(X, type = 'p', main = "2D Spiral", xlab = "X", ylab = "Y", asp = 1)
}

# Example usage: plot a spiral with 5 turns
plot_spiral(3)

