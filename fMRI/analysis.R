library(fabisearch)
library(ggplot2)
library(reshape2)
library(viridis)
library(rgl)
sapply((paste0("../Core/", list.files("../Core/"))), source)

# Get data for all models
models <- c("SLICE", "nnLVGLASSO", "tGLASSO", "rcLVGLASSO")
visits <- c("visit2", "visit3")
df <- c()
for(model in models){
  for(visit in visits){
    df <- rbind(df, cbind(model, visit, read.table(paste0(model, "_", visit, ".txt")))) 
  }
}
mean_df <- aggregate(. ~ model + visit, data = df, FUN = mean)

# Load data
Xall <- read.table("visit2.txt")
AAL <- read.table("AAL.txt")

# Get first subjets data
X <- Xall[1:190,]

# Set seed for estimation
set.seed(123)
models <- c("SLICE", "nnLVGLASSO", "rcLVGLASSO")
for(model in models){
  if(model == "SLICE"){
    cvout <- cv.slice(X, rhos = 0.5) # SLICE
    out <- slice(cov(X), cvout$rho, cvout$r)
  } else if(model == "nnLVGLASSO"){
    cvout <- cv.nnlvg(X, rhos = 0.3) # nnLVGLASSO
    out <- nnlvg(cov(X), cvout$rho, cvout$gamma)
  } else if(model == "rcLVGLASSO"){
    cvout <- cv.rclvg(X, rhos = 0.3) # rcLVGLASSO
    out <- rclvg(cov(X), cvout$rho, cvout$r)
  } else if(model == "tGLASSO"){
    out <- ebic.tg(X, taus = logseq(1e-5, 0.2, 20)) # tGLASSO
  }
  # Latent
  rownames(out$L) <- colnames(out$L) <- 1:nrow(out$L)
  rownames(out$L)[!(1:length(rownames(out$L)) %% 10 == 0)] <- ""
  colnames(out$L)[!(1:length(colnames(out$L)) %% 10 == 0)] <- ""
  png(paste0(model, "_L.png"), width = 600, height = 600)
  pheatmap(out$L, 
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           color = viridis(256),
           legend = FALSE,
           border_color = NA,
           fontsize_row = 20,
           fontsize_col = 20)
  dev.off()
  
  # Sparse
  net.3dplot(1*(out$S != 0), coordROIs = cbind("c",AAL[,3:5]))
  rglwidget()
}