library(MASS)

runsim <- function(simtype, method, pobs, plat = NULL, n, iters){
  S_hats <- L_hats <- S_stars <- L_stars <- z_stars <- vector("list", length(iters))
  # Loop through 100 iterations of simulation
  for(i in 1:length(iters)){
    print(paste0("SIM ITER ", iters[i]))
    set.seed(123*iters[i])
    
    dir_path <- paste0("./", method) # Make directory if it doesn't exist
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    S_star <- Smat(pobs, 2, 1.5)
    S_star[S_star < 0.01] <- 0 # True sparse component
    
    if (simtype == "rand"){
      Lout <- Lrand(pobs, plat, 1.5)
    } else if (simtype == "cres"){
      Lout <- Lcres(pobs, 0.1)
      plat <- 2
    } else if (simtype == "spir"){
      Lout <- Lspir(pobs, 0.05)
      plat <- 2
    } else if (simtype == "eeg"){
      load("eeg_sim_params.rda")
    }
    
    L_star <- Lout$L; z_star <- Lout$z # Latent component; cluster labels
    
    Sigma_star <- solve(S_star + L_star) # Sigma
    X <- mvrnorm(n, rep(0, pobs), Sigma = Sigma_star) # Finite sample data
    Sigma <- cov(X) # Sample Sigma
    
    # Model selection and parameter estimation
    if(method == "SLICE"){
      if(simtype == "eeg"){
        cvsli <- cv.slice(X, rs = 5:10)
      } else{
        cvsli <- cv.slice(X) 
      }
      sli <- slice(Sigma, cvsli$rho, cvsli$r)
      S <- sli$S; L <- sli$L
    } else if(method == "SLICE_GSCAD"){
      if(simtype == "eeg"){
        cvsli <- cv.slice(X, rs = 5:10, Sest = "gscad")
      } else{
        cvsli <- cv.slice(X, Sest = "gscad")
      }
      sli <- slice(Sigma, cvsli$rho, cvsli$r, Sest = "gscad")
      S <- sli$S; L <- sli$L
    } else if(method == "SLICE_CLIME"){
      if(simtype == "eeg"){
        cvsli <- cv.slice(X, rs = 5:10, Sest = "clime")
      } else{
        cvsli <- cv.slice(X, Sest = "clime")
      }
      sli <- slice(Sigma, cvsli$rho, cvsli$r, Sest = "clime")
      S <- sli$S; L <- sli$L
    } else if(method == "nnLVGLASSO"){
      cvnnlvg <- cv.nnlvg(X)
      nnlvg <- nnlvg(Sigma, cvnnlvg$rho, cvnnlvg$gamma)
      S <- nnlvg$S; L <- nnlvg$L
    } else if(method == "rcLVGLASSO"){
      if(simtype == "eeg"){
        cvrclvg <- cv.rclvg(X, rs = 5:10)
      } else{
        cvrclvg <- cv.rclvg(X)
      }
      rclvg <- rclvg(Sigma, cvrclvg$rho, cvrclvg$r)
      S <- rclvg$S; L <- rclvg$L
    } else if(method == "tGLASSO"){
      S <- ebic.tg(X)$S
      L <- NULL
    }
    
    # Append to all outputs
    S_hats[[i]] <- S; L_hats[[i]] <- L
    S_stars[[i]] <- S_star; L_stars[[i]] <- L_star; z_stars[[i]] <- z_star
    
    # Save as rda
    save(S_hats, L_hats, S_stars, L_stars, z_stars, 
         file = paste0(method, "/", method, "_sim", simtype, "_plat", 
                       plat, "_n", n, "_iters", min(iters), 
                       "to", max(iters), ".rda"))
  }
}