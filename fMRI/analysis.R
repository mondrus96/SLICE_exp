library(fabisearch)
library(ggplot2)
library(reshape2)
library(viridis)
library(rgl)
library(pheatmap)
library(colorspace)
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
labels <- sub(".*\\s([^_]+)_.*", "\\1", AAL$Labels)
labels <- data.frame(Group = labels)

# Get first subjets data
X <- Xall[1:190,]

# Set seed for estimation
set.seed(123)
palette <- make_palette(28) # Make colours for all
group_colors <- list(Group = (setNames(palette, unique(labels$Group))))
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
  rownames(out$L) <- colnames(out$L) <- rep("", nrow(out$L))
  png(paste0(model, "_L.png"), width = 750, height = 600)
  pheatmap(out$L,
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           color = viridis(256),
           legend = FALSE,
           border_color = NA,
           annotation_col = labels,
           annotation_colors = group_colors,
           cellheight=4, cellwidth=4)
  dev.off()
  
  # Sparse
  net.3dplot(1*(out$S != 0), coordROIs = cbind("c",AAL[,3:5]))
  rglwidget()
}
