library(fabisearch)
library(ggplot2)
library(reshape2)
library(viridis)
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

# Set seed for t-SNE
set.seed(123)
models <- c("SLICE", "nnLVGLASSO", "rcLVGLASSO")
for(model in models){
  if(model == "SLICE"){
    cvout <- cv.slice(X, lambdas = 0.5) # SLICE
    out <- slice(cov(X), cvout$lambda, cvout$r)
  } else if(model == "nnLVGLASSO"){
    cvout <- cv.nnlvg(X, lambdas = 0.3) # nnLVGLASSO
    out <- nnlvg(cov(X), cvout$lambda, cvout$gamma)
  } else if(model == "rcLVGLASSO"){
    cvout <- cv.rclvg(X, lambdas = 0.3) # rcLVGLASSO
    out <- rclvg(cov(X), cvout$lambda, cvout$r)
  } else if(model == "tGLASSO"){
    out <- ebic.tg(X, taus = logseq(1e-5, 0.2, 20)) # tGLASSO
  }
  Llong <- melt(out$L)
  p <- ggplot(Llong, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c() +  # Optional: Use a Viridis color scale
    labs(fill = "Value") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  ggsave(paste0(model, "_L.png"), p)
  net.3dplot(1*(out$S != 0), coordROIs = cbind("c",AAL[,3:5]))
  rglwidget()
}
