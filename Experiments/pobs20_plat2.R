# Load glasso functions and matrix generation
source("gen_mat.R")
source("../Solvers/ADMMlatentglasso.R")

# Load glasso library
library(glasso)
library(MASS)
library(randnet)
library(ggplot2)
library(dplyr)
library(tidyr)

# Set seed
set.seed(123)

pobs = 100 # Intralayer size
plat = 3 # Number of latent variables for L

### 1 Simulate covariance ###
reps = 100 # Number of simulation repetitions
ns = c(seq(25, 300, 25), NA) # Number of samples

latmu = 1.5
latsd = 0

decay_rate = 0.1 # Decay rate for sparse matrix Theta
init_val = 1.5 # Initial value for Theta

df = c()
for(i in 1:reps){
  print(paste0("Sim #", i))
  
  # Get sparse theta components
  Thetalist = multidecay(pobs, 5, 1.5, 0.01)
  Theta = Thetalist[[1]]
  
  # Define L
  clustout = Lclust(sum(p), plat, latmu, latsd)
  L = clustout$L; Z = clustout$Z
  
  # Define true Sigma, Sigma_star
  Sigma_star = solve(Theta + L)
  
  for(j in 1:length(ns)){
    # Create simulated data
    if (!is.na(ns[j])){
      X = mvrnorm(ns[j], rep(0, p), Sigma_star)
      Sigma = cov(X) 
    } else {
      Sigma = Sigma_star
    }
    
    # Get latent glasso estimate
    lg = admmlatentglasso(Sigma, 0.05, 0.01, 0.5)
    
    # Get glasso estimate
    g = glasso(Sigma, 0.05)$wi
    
    # Threshold glasso estimate
    tg = g
    tg[tg < 0.2] = 0
    
    # Make Sigma PD
    invSigma = solve(Sigma + diag(rep(0.1, dim(Sigma)[1])))
    
    # Truncated eigendecomp
    iltg = eigen(invSigma - tg)
    iltg = iltg$vectors[,1:plat] %*% diag(iltg$values[1:plat]) %*% t(iltg$vectors[,1:plat])
    
    # Community detection
    lg_clust = kmeans(eigen(lg$L)$vectors[,1:plat], plat, 1000)$cluster
    iltg_clust = kmeans(eigen(iltg)$vectors[,1:plat], plat, 1000)$cluster
    true_clust = kmeans(eigen(L)$vectors[,1:plat], plat, 1000)$cluster
    
    # Compare the two
    iltg_nmi = NMI(true_clust, iltg_clust)
    lg_nmi = NMI(true_clust, lg_clust)
    
    # Add to results
    df = rbind(df, c(ns[j], iltg_nmi, lg_nmi))
  }
  group_means <- as.data.frame(df) %>%
    group_by(V1) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  print(group_means)
}

# Calculating mean for each group
df = as.data.frame(df)
colnames(df) = c("n", "lrtgl", "lvgl")
group_means <- as.data.frame(df) %>%
  group_by(n) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
group_means[,2:3] = log(group_means[,2:3])

# Convert the dataframe to long format
group_means_long <- group_means %>%
  pivot_longer(cols = c("lrtgl", "lvgl"), names_to = "Variable", values_to = "Value")

# Plot the data
ggplot(group_means_long, aes(x = n, y = Value, color = Variable)) +
  geom_line() +
  theme_minimal() +
  labs(x = "n", y = "log(NMI)", color = "Method")