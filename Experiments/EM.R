# EM ALGORITHM?
source("../Simulations/sims.R")
source("../Models/slice.R")
source("../Models/cv.slice.R")
source("../Models/helper.R")
library(clime)
library(glasso)
library(Matrix)

# Set seed
set.seed(123)

pobs = 200 # Number of observed variables for S
plat = 4 # Number of latent variables for L
n = 10000 # Number of observations

init_S = 1.5 # Initial value for S
init_L = 1.5 # Initial value for L

# Get sparse S_star components
S_star = Smat(pobs, 1.5, init_S)
S_star[S_star < 0.01] = 0

# Define L
clustout = Lmat(pobs, plat, init_L)
L_star = clustout$L; Z_star = clustout$Z
heatmap(L_star)

# Add noise
noise = matrix(rnorm(pobs * pobs, 0, 0.5), pobs, pobs)
noise = noise %*% t(noise)

# Define true Sigma, Sigma_star
Sigma_star = solve(S_star + L_star)

# Define true Sigma, Sigma_star
library(MASS)
X = mvrnorm(n, rep(0, pobs), Sigma = Sigma_star)
sigma = cov(X)

# Run slice
cvout = cv.slice(X, lambdas = c(0.1, 0.2), rs = 4:5)
out = slice(sigma, 0.1, 4)
S = out$S; L = out$L

heatmap(1*(S != 0), Rowv = NA, Colv = NA)
heatmap(1*(S_star != 0), Rowv = NA, Colv = NA)
heatmap(L, Rowv = NA, Colv = NA)
heatmap(L_star, Rowv = NA, Colv = NA)

eigL = eigen(L)$vectors[,1:plat]
eigL_star = eigen(L_star)$vectors[,1:plat]

predclust = kmeans(eigL, plat, 10000)$cluster
trueclust = kmeans(eigL_star, plat, 10000)$cluster
NMI(predclust, trueclust)

v_hat = eigL[,1]
v = eigL_star[,1]

dot_product <- sum(v_hat * v)  # Calculate dot product
sin_theta <- sqrt(1 - dot_product^2)  # Calculate sine of the angle
sin_theta

norm(L - L_star)
norm(Omega - S_star)


library(ggplot2)

# Function to generate 'circles' data
generate_circles <- function(n_samples, noise) {
  set.seed(42)
  t <- runif(n_samples, 0, 2 * pi)
  r <- sqrt(runif(n_samples, 0, 1))
  x <- c(r * cos(t), (r + 0.5) * cos(t))
  y <- c(r * sin(t), (r + 0.5) * sin(t))
  x <- x + rnorm(length(x), sd=noise)
  y <- y + rnorm(length(y), sd=noise)
  data.frame(x = x, y = y)
}

# Function to generate 'moons' data
generate_moons <- function(n_samples, noise) {
  set.seed(42)
  t <- runif(n_samples / 2, 0, pi)
  x <- c(cos(t), 1 - cos(t))
  y <- c(sin(t), 1 - sin(t))
  x <- x + rnorm(length(x), sd=noise)
  y <- y + rnorm(length(y), sd=noise)
  data.frame(x = x, y = y)
}

# Generate the datasets
n_samples <- 300
circles_data <- generate_circles(n_samples, noise = 0.05)
moons_data <- generate_moons(n_samples, noise = 0.05)

# Plot the datasets using ggplot2
ggplot() +
  geom_point(data = circles_data, aes(x = x, y = y), color = "blue") +
  ggtitle("Circles") +
  theme_minimal() +
  xlim(c(-1.5, 2.5)) + ylim(c(-1.5, 1.5))

ggplot() +
  geom_point(data = moons_data, aes(x = x, y = y), color = "red") +
  ggtitle("Moons") +
  theme_minimal() +
  xlim(c(-0.5, 2.5)) + ylim(c(-0.5, 1.5))
