library(ggplot2)
library(dplyr)
library(Matrix)

### Getting Summary Values ###
sapply((paste0("../Core/", list.files("../Core/"))), source)
# Get file names
models <- c("SLICE", "SLICE_CLIME", "SLICE_GSCAD",
            "nnLVGLASSO", "tGLASSO", "rcLVGLASSO")
allsims <- data.frame()
for(model in models){
  print(model)
  files <- list.files(paste0("../Simulations/", model, "/"), "n10000")
  for(j in 1:length(files)){
    load(paste0("../Simulations/", model, "/", files[j]))
    simtype <- strsplit(strsplit(files[j], paste0(model, "_"))[[1]][2], "_plat")[[1]][1]
    plat <- as.numeric(strsplit(strsplit(files[j], "plat")[[1]][2], "_")[[1]][1])
    
    # Loop through each one
    indx <- which(!sapply(S_stars, is.null))
    iters <- strsplit(strsplit(strsplit(files[j], "iters")[[1]][2], ".rda")[[1]][1], "to")[[1]]
    iters <- as.numeric(iters[1]):as.numeric(iters[2])
    iters <- iters[indx]
    F1 <- TP <- TN <- sin_theta <- frob_norm <- spec_norm <- ARI <- rep(NA, length(iters))
    for(k in 1:length(iters)){
      # Sparse metrics
      F1[k] <- F1score(S_stars[[k]], S_hats[[k]])
      TP[k] <- TPrate(S_stars[[k]], S_hats[[k]])
      TN[k] <- TNrate(S_stars[[k]], S_hats[[k]])
      
      # Low rank metrics
      if(model != "tGLASSO"){
        eigL_hat <- eigen(L_hats[[k]])
        eigL_star <- eigen(L_stars[[k]])
        
        # Sin theta
        sin_theta[k] <- sintheta(eigL_star$vectors[,1], eigL_hat$vectors[,1])
        
        # Frobenius norm
        frob_norm[k] <- norm(L_stars[[k]] - L_hats[[k]])
        
        # The spectral norm
        eigDiff <- eigen(L_stars[[k]] - L_hats[[k]])
        spec_norm[k] <- eigDiff$values[1]
        
        # Clustering
        if(simtype == "simrand"){
          # Figure out rank of matrix
          r <- rankMatrix(L_hats[[k]])
          
          # Do kmeans clustering
          z_hat <- kmeans(eigL_hat$vectors[,1:plat], max(z_stars[[k]]), 
                          iter.max = 100, nstart = 1000)$cluster
          ARI[k] <- ari(z_stars[[k]], z_hat)
        }
      }
    }
    allsims <- rbind(allsims, 
                     cbind(model, simtype, iters, TP, TN, F1, 
                           sin_theta, frob_norm, spec_norm, ARI)) # Add to df
  }
}

# Get summary values
rm(list = setdiff(ls(), "allsims"))
allsims <- allsims %>%
  mutate(TP = as.numeric(TP),
         TN = as.numeric(TN),
         F1 = as.numeric(F1),
         sin_theta = as.numeric(sin_theta),
         frob_norm = as.numeric(frob_norm),
         spec_norm = as.numeric(spec_norm),
         ARI = as.numeric(ARI))

df_mean <- allsims %>%
  group_by(model, simtype) %>%
  summarize(across(c(TP, TN, F1, 
                     sin_theta, frob_norm, spec_norm, ARI), 
                   mean, na.rm = TRUE), .groups = 'drop')
write.table(df_mean, "highn_mean.txt")
df_sd <- allsims %>%
  group_by(model, simtype) %>%
  summarize(across(c(TP, TN, F1, sin_theta, frob_norm, spec_norm, ARI), 
                   sd, na.rm = TRUE), .groups = 'drop')
write.table(df_sd, "highn_sd.txt")