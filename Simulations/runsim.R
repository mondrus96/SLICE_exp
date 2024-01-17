library(MASS)

# Function for running simulations
runsim = function(simtype, pobs, plat = NULL, n, iters){
  df = c()
  # Loop through 100 iterations of simulation
  for(i in iters){
    set.seed(123*i)
    
    S_star <- Smat(pobs, 2, 1.5)
    S_star[S_star < 0.01] <- 0 # True sparse component
    
    if (simtype == "exp"){
      Lout <- Lexp(pobs, plat, 1.5)
    } else if (simtype == "cres"){
      Lout <- Lcres(pobs)
      plat <- 2
    } else if (simtype == "circ"){
      Lout <- Lcirc(pobs)
      plat <- 2
    }
    
    L_star <- Lout$L; z_star <- Lout$z # True latent component; True cluster labels
    
    Sigma_star <- solve(S_star + L_star) # True Sigma
    if(n == "inf"){
      Sigma <- Sigma_star
    } else if(is.numeric(n)){
      X <- mvrnorm(n, rep(0, pobs), Sigma = Sigma_star) # Finite sample data
      Sigma <- cov(X) # Sample Sigma 
    }
    
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
    colnames(df) = c("plat", "n", "lvg_nmi", "sli_nmi",
                     "lvg_ari", "sli_ari", "lvg_sin", "sli_sin")
    rownames(df) = NULL
    write.table(df, file = paste0("sim", simtype, "_plat", 
                                  plat, "_n", n, ".txt"), row.names = FALSE)
  }
}