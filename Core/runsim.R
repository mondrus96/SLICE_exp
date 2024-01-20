library(MASS)

runsim = function(simtype, pobs, plat = NULL, n, iters){
  df = c()
  # Loop through 100 iterations of simulation
  for(i in iters){
    print(paste0("SIM ITER ", i))
    set.seed(123*i)
    
    S_star <- Smat(pobs, 2, 1.5)
    S_star[S_star < 0.01] <- 0 # True sparse component
    
    if (simtype == "exp"){
      Lout <- Lexp(pobs, plat, 1.5)
    } else if (simtype == "rand"){
      Lout <- Lrand(pobs, plat, 1.5)
    } else if (simtype == "cres"){
      Lout <- Lcres(pobs, 0.1)
      plat <- 2
    } else if (simtype == "circ"){
      Lout <- Lcirc(pobs, 0.05)
      plat <- 2
    } else if (simtype == "spir"){
      Lout <- Lspir(pobs, 0.05)
      plat <- 2
    }
    
    L_star <- Lout$L; z_star <- Lout$z # True latent component; True cluster labels
    
    Sigma_star <- solve(S_star + L_star) # True Sigma
    X <- mvrnorm(n, rep(0, pobs), Sigma = Sigma_star) # Finite sample data
    Sigma <- cov(X) # Sample Sigma 
    
    cvlvg <- cv.lvg(X) # Cross validation
    cvsli <- cv.slice(X)
    
    lvg <- lvglasso(Sigma, cvlvg$lambda, cvlvg$gamma) # Best fit
    sli <- slice(Sigma, cvsli$lambda, cvsli$r)
    
    eigLlvg <- eigen(lvg$L) # Eigendecomp
    eigLsli <- eigen(sli$L)
    eigL_star <- eigen(L_star)
    
    lvg_sin <- sintheta(eigL_star$vectors[,1], eigLlvg$vectors[,1]) # Sin angle
    sli_sin <- sintheta(eigL_star$vectors[,1], eigLsli$vectors[,1])
    
    lvg_fnorm <- norm(L_star - lvg$L) # Frobenius norm
    sli_fnorm <- norm(L_star - sli$L)
    
    if(simtype %in% c("exp", "rand")){
      z_lvg <- kmeans(eigLlvg$vectors[,1:plat], plat, 1000)$cluster # Cluster predictions
      z_sli <- kmeans(eigLsli$vectors[,1:plat], plat, 1000)$cluster
      
      lvg_nmi <- nmi(z_star, z_lvg) # NMI
      sli_nmi <- nmi(z_star, z_sli) 
      
      lvg_ari <- ari(z_star, z_lvg) # ARI
      sli_ari <- ari(z_star, z_sli)
      
      df = rbind(df, c(plat, n, lvg_nmi, sli_nmi, 
                       lvg_ari, sli_ari, lvg_sin, sli_sin,
                       lvg_fnorm, sli_fnorm))
      colnames(df) = c("plat", "n", "lvg_nmi", "sli_nmi",
                       "lvg_ari", "sli_ari", "lvg_sin", "sli_sin",
                       "lvg_fnorm", "sli_fnorm")
    } else{
      df = rbind(df, c(plat, n, lvg_sin, sli_sin,
                       lvg_fnorm, sli_fnorm))
      colnames(df) = c("plat", "n", "lvg_sin", "sli_sin",
                       "lvg_fnorm", "sli_fnorm")
    }
    
    rownames(df) = NULL
    write.table(df, file = paste0("sim", simtype, "_plat", 
                                  plat, "_n", n, "_iters", min(iters), "to", max(iters), ".txt"), row.names = FALSE)
  }
}