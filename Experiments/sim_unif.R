# Uniform cluster sizes
sapply((paste0("../Models/", list.files("../Models/"))), source)
sapply((paste0("../Simulations/", list.files("../Simulations/"))), source)

# Set seed
set.seed(123)

pobs <- 200 # Number of observed variables for S
plats <- 2:8 # Number of latent variables for L
ns <- seq(50, 300, 50) # Number of observations

init_S <- 1.5 # Initial value for S
init_L <- 1.5 # Initial value for L

df = c() # Results
# Loop through different number of latent variables
runsim = function(pobs, plats, ns, init_S, init_L){
  df = c()
  for(i in 1:length(plats)){
    print(paste0("plat: ", plats[i]))
    # Loop through different values of ns
    for(j in 1:length(ns)){
      print(paste0("n: ", ns[j]))
      # Loop through 100 iterations of simulation
      subdf = foreach(k = 1:100, .combine = rbind, .packages = c("MASS")) %dopar% {
        set.seed(k*123)
        
        S_star <- Smat(pobs, 1.5, init_S)
        S_star[S_star < 0.01] <- 0 # True sparse component
        
        Lout <- Luni(pobs, plats[i], init_L) 
        L_star <- Lout$L; z_star <- Lout$z # True latent component; True cluster labels
        
        Sigma_star <- solve(S_star + L_star) # True Sigma
        X <- mvrnorm(ns[j], rep(0, pobs), Sigma = Sigma_star) # Finite sample data
        Sigma <- cov(X) # Sample Sigma
        
        cvlvg <- cv.lvg(X) # Run lvglasso
        lvg <- lvglasso(Sigma, cvlvg$lambda, cvlvg$gamma)
        
        cvsli <- cv.slice(X) # Run slice
        sli <- slice(Sigma, cvsli$lambda, cvsli$r)
        
        eigLlvg <- eigen(lvg$L) # Eigendecomp
        eigLsli <- eigen(sli$L)
        eigL_star <- eigen(L_star)
        
        z_lvg <- kmeans(eigLlvg$vectors[,1:plats[i]], plats[i], 1000)$cluster # Cluster predictions
        z_sli <- kmeans(eigLsli$vectors[,1:plats[i]], plats[i], 1000)$cluster
        
        lvg_nmi <- nmi(z_star, z_lvg) # NMI
        sli_nmi <- nmi(z_star, z_sli) 
        
        lvg_ari <- ari(z_star, z_lvg) # ARI
        sli_ari <- ari(z_star, z_sli)
        
        lvg_sin <- sintheta(eigL_star$vectors[,1], eigLlvg$vectors[,1]) # Sin Theta
        sli_sin <- sintheta(eigL_star$vectors[,1], eigLsli$vectors[,1])
        
        return(c(plats[i], ns[j], lvg_nmi, sli_nmi, 
                 lvg_ari, sli_ari, lvg_sin, sli_sin))
      }
      df = rbind(df, subdf)
      colnames(df) = c("plat", "n", "lvg_nmi", "sli_nmi",
                       "lvg_ari", "sli_ari", "lvg_sin", "sli_sin")
      rownames(df) = NULL
      write.table(df, file = paste0("sim", simtype), row.names = FALSE)
    }
  }
}