# EM ALGORITHM?
source("../Simulations/gen_mat.R")
source("../Solvers/ilrglasso.R")
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
S_star = expdecay(pobs, 1.5, init_S)
S_star[S_star < 0.01] = 0

# Define L
clustout = Lclust(pobs, plat, init_L, latsd)
L_star = clustout$L; Z_star = clustout$Z

# Add noise
noise = matrix(rnorm(pobs * pobs, 0, 0.5), pobs, pobs)
noise = noise %*% t(noise)

# Define true Sigma, Sigma_star
Sigma_star = solve(S_star + L_star)

# Define true Sigma, Sigma_star
library(MASS)
data = mvrnorm(n, rep(0, pobs), Sigma = Sigma_star)
Sigma = cov(data)

TOL = 1e-8
nfolds = 3
lambdas = c(0.0001, 0.001, 0.01, 0.1, 0.2)
rs = c(2, 3, 4, 5, 6)
modselec = matrix(NA, length(rs), length(lambdas)); rownames(modselec) = rs; colnames(modselec) = lambdas

ebics = bics = CVs = modselec

# Loop through different rank and lambdas
for(k in 1:length(rs)){
  print(paste0("rank: ", rs[k]))
  for(j in 1:length(lambdas)){
    print(paste0("lambda: ", lambdas[j]))
    # Run estimation
    indices = sample(1:nfolds, n, replace = TRUE)
    meanlikl = c()
    for(l in 1:nfolds){
      # Define train and test
      train = data[indices != l,]
      test = data[indices == l,]
      
      # Run ilrglasso
      out = ilrglasso(cov(train), lambdas[j], rs[k])
      S = out$S; L = out$L
      
      # Append to meanlikl
      loglikl = loglikelihood(cov(test), S, L)
      meanlikl = c(meanlikl, loglikl)
    }
    CVs[k, j] = mean(meanlikl)
    
    # Over whole dataset for ebic and bic  
    out = ilrglasso(Sigma, lambdas[j], rs[k])
    S = out$S; L = out$L
    
    # Get likelihood
    loglikl = loglikelihood(Sigma, S, L)
    S_params = sum(S[upper.tri(S)] != 0)
    L_params = rs[k] + (rs[k]*(pobs-1))/2
    bics[k, j] = bic(loglikl, n, (L_params + S_params))
    ebics[k, j] = ebic(loglikl, pobs, n, (L_params + S_params))
  }
}
loglikl = loglikelihood(Sigma, S, L)
S_params = sum(S[upper.tri(S)] != 0)
L_params = r + (r*(pobs-1))/2
bic(loglikl, n, (L_params + S_params))


out = ilrglasso(Sigma, 0.01, 4)
S = out$S; L = out$L

heatmap(1*(S != 0), Rowv = NA, Colv = NA)
heatmap(1*(S_star != 0), Rowv = NA, Colv = NA)
heatmap(L, Rowv = NA, Colv = NA)
heatmap(L_star, Rowv = NA, Colv = NA)

eigL = eigen(L)$vectors[,1:plat]
eigL_star = eigen(L_star)$vectors[,1:plat]

library(randnet)

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