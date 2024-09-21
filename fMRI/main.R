# Load functions
sapply((paste0("../Core/", list.files("../Core/"))), source)

# Load data
Xname <- "visit3"
Xall <- read.table(paste0(Xname, ".txt"))
AAL <- read.table("AAL.txt")

# Define the coords variable - the coordintes relative to the variables
coords <- AAL[,3:5]
distmat <- as.matrix(dist(coords))
distmat[lower.tri(distmat)] <- 0

# Define subjects and length of experiment
subj <- 25; tlen <- 190

# Accessory function
wdist <- function(S, distmat){
  return(sum((abs(S)*distmat))/sum(S[upper.tri(S)] != 0))
}

# Set seed
set.seed(123)

# Loop through subjects
df <- Fnorm <- wd <- sparsity <- c()
for(i in 1:subj){
  print(i) # Current subject
  
  inds <- (((i-1)*tlen)+1):(i*tlen) # Get the length of indices
  X1 <- as.matrix(Xall[inds,])
  Sig1 <- cov(X1)
  
  if(model == "SLICE"){
    cvout <- cv.slice(X1, rhos = 0.5) # SLICE
    out <- slice(Sig1, cvout$rho, cvout$r)
  } else if(model == "nnLVGLASSO"){
    cvout <- cv.nnlvg(X1, rhos = 0.3) # nnLVGLASSO
    out <- nnlvg(Sig1, cvout$rho, cvout$gamma)
  } else if(model == "rcLVGLASSO"){
    cvout <- cv.rclvg(X1, rhos = 0.3) # rcLVGLASSO
    out <- rclvg(Sig1, cvout$rho, cvout$r)
  } else if(model == "tGLASSO"){
    out <- ebic.tg(X1, taus = logseq(1e-5, 0.2, 20)) # tGLASSO
    out$L <- NULL; out2 <- ebic.tg(X2)
  }

  wd <- c(wd, wdist(out$S, distmat))
  sparsity <- c(sparsity, sum(out$S[upper.tri(out$S)] == 0)/sum(upper.tri(out$S)))
  
  df <- cbind(Fnorm, wd, sparsity)
  print(colMeans(df))
  write.table(df, file = paste0(model, "_", Xname, ".txt"))
}
