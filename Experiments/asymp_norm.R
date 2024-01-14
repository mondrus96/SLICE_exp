# Load glasso functions and matrix generation
source("../Simulations/gen_mat.R")
source("../Solvers/ADMMlatentglasso.R")

# Load glasso library
library(glasso)
library(MASS)
library(randnet)
library(ggplot2)
library(dplyr)
library(plot3D)
library(tidyr)
library(RColorBrewer)

# Set seed
set.seed(123)

pobs = 100 # Number of observed variables for S
plat = 5 # Number of latent variables for L

reps = 100 # Number of simulation repetitions
ns = NA # Number of samples

decay_rate = 0.1 # Decay rate for sparse matrix Theta
init_S = seq(0.5, 1.5, 0.25) # Initial value for S
init_L = seq(0.5, 1.5, 0.25) # Initial value for L

df = data.frame()
for(i in 1:reps){
  print(paste0("Sim #", i))
  for(k in 1:length(init_S)){
    print(paste0("init_S ", init_S[k]))
    for(j in 1:length(init_L)){
      print(paste0("init_L ", init_L[j]))
      
      # Get sparse S_star components
      S_star = expdecay(pobs, 2, init_S[k])
      S_star[S_star < 0.01] = 0
      
      # Define L
      clustout = Lclust(pobs, plat, init_L[j], latsd)
      L_star = clustout$L; Z_star = clustout$Z
      
      # Define true Sigma, Sigma_star
      Sigma = solve(S_star + L_star)
        
      # Get latent glasso estimate
      lvg = admmlatentglasso(Sigma, 0.1, 0.01, 0.5)
      
      # Get glasso estimate
      g = glasso(Sigma, 0.1)$wi
      
      # Threshold glasso estimate
      tg = g
      tg[tg < 0.2] = 0
      
      # Make Sigma PD
      if(is.na(ns[j])){
        SigmaPD = Sigma
      } else{
        SigmaPD = Sigma + diag(rep(0.2, ncol(Sigma))) 
      }
      invSigmaPD = solve(SigmaPD) # Invert
      
      # Truncated eigendecomp
      ilr = eigen(invSigmaPD - tg)
      ilr = ilr$vectors[,1:plat] %*% diag(ilr$values[1:plat]) %*% t(ilr$vectors[,1:plat])
      
      # Evaluate Frobenius norm
      ilr_fnorm = norm(ilr - L_star, type = "F")
      lvg_fnorm = norm(lvg$L - L_star, type = "F")
      
      # Evaluate trace
      ilr_trace = norm(diag(diag(ilr) - diag(L_star)), type = "F")
      lvg_trace = norm(diag(diag(lvg$L) - diag(L_star)), type = "F")
      
      # Community detection
      lvg_clust = kmeans(eigen(lvg$L)$vectors[,1:plat], plat, 1000)$cluster
      ilr_clust = kmeans(eigen(ilr)$vectors[,1:plat], plat, 1000)$cluster
      true_clust = kmeans(eigen(L_star)$vectors[,1:plat], plat, 1000)$cluster
      
      # Compare the clusters
      ilr_nmi = NMI(true_clust, ilr_clust)
      lvg_nmi = NMI(true_clust, lvg_clust)
      
      # Add to results
      df = rbind(df, c(init_S[k], init_L[j], ilr_nmi, lvg_nmi, ilr_trace, lvg_trace, ilr_fnorm, lvg_fnorm))
    }
  }
}
# Rename columns
colnames(df) = c("init_S", "init_L", "ilr_nmi", "lvg_nmi", "ilr_trace", "lvg_trace", "ilr_fnorm", "lvg_fnorm")
save(df, file = "asymp.rda")

# Calculating mean for each group by 'n' and 'plat'
group_means <- df %>%
  group_by(init_S, init_L) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

# Calcualte the trace and fnorm for each one
ilr_trace <- with(group_means, tapply(ilr_trace, list(init_L, init_S), mean))
lvg_trace <- with(group_means, tapply(lvg_trace, list(init_L, init_S), mean))
ilr_fnorm <- with(group_means, tapply(ilr_fnorm, list(init_L, init_S), mean))
lvg_fnorm <- with(group_means, tapply(lvg_fnorm, list(init_L, init_S), mean))
