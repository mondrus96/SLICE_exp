library(MASS)

# Function for running simulations
runsim = function(simtype, pobs, plat = NA, ns, init_S, init_L = NULL){
  df = c()
  # Loop through different values of ns
  for(i in 1:length(ns)){
    print(paste0("n: ", ns[i]))
    # Loop through 100 iterations of simulation
    for(j in 1:100){
      S_star <- Smat(pobs, 2, init_S)
      S_star[S_star < 0.01] <- 0 # True sparse component
      
      if (simtype == "unif"){
        Lout <- Lunif(pobs, plat, init_L)
      } else if (simtype == "exp"){
        Lout <- Lexp(pobs, plat, init_L)
      } else if (simtype == "cres"){
        Lout <- Lcres(pobs)
      } else if (simtype == "circ"){
        Lout <- Lcirc(pobs)
      }
      
      L_star <- Lout$L; z_star <- Lout$z # True latent component; True cluster labels
      
      Sigma_star <- solve(S_star + L_star) # True Sigma
      X <- mvrnorm(ns[i], rep(0, pobs), Sigma = Sigma_star) # Finite sample data
      Sigma <- cov(X) # Sample Sigma
      
      cvlvg <- cv.lvg(X) # Run lvglasso
      lvg <- lvglasso(Sigma, cvlvg$lambda, cvlvg$gamma)
      
      cvsli <- cv.slice(X) # Run slice
      sli <- slice(Sigma, cvsli$lambda, cvsli$r)
      
      eigLlvg <- eigen(lvg$L) # Eigendecomp
      eigLsli <- eigen(sli$L)
      eigL_star <- eigen(L_star)
      
      z_lvg <- kmeans(eigLlvg$vectors[,1:plat], plat, 1000)$cluster # Cluster predictions
      z_sli <- kmeans(eigLsli$vectors[,1:plat], plat, 1000)$cluster
      
      lvg_nmi <- nmi(z_star, z_lvg) # NMI
      sli_nmi <- nmi(z_star, z_sli) 
      
      lvg_ari <- ari(z_star, z_lvg) # ARI
      sli_ari <- ari(z_star, z_sli)
      
      lvg_sin <- sintheta(eigL_star$vectors[,1], eigLlvg$vectors[,1]) # Sin Theta
      sli_sin <- sintheta(eigL_star$vectors[,1], eigLsli$vectors[,1])
      
      df = rbind(df, c(plat, ns[i], lvg_nmi, sli_nmi, 
               lvg_ari, sli_ari, lvg_sin, sli_sin))
    }
    colnames(df) = c("plat", "n", "lvg_nmi", "sli_nmi",
                     "lvg_ari", "sli_ari", "lvg_sin", "sli_sin")
    rownames(df) = NULL
    write.table(df, file = paste0("sim_", simtype, ".txt"), row.names = FALSE)
  }
}
