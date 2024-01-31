# Load functions
sapply((paste0("../Core/", list.files("../Core/"))), source)

# Load data
Xall <- read.table("visit2.txt")
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
  midp <- length(inds)/2 # Find the midpoint
    
  X <- as.matrix(Xall[inds,])
  X1 <- X[1:midp,]; X2 <- X[(midp+1):length(inds),]
  Sig1 <- cov(X1); Sig2 <- cov(X2);
  
  if(model == "SLICE"){
    cvout <- cv.slice(X1, lambdas = seq(0.1, 0.5, 0.1)) # SLICE
    out <- slice(Sig1, cvout$lambda, cvout$r); out2 <- slice(Sig2, cvout$lambda, cvout$r)
  } else if(model == "nnLVGLASSO"){
    cvout <- cv.nnlvg(X1, lambdas = seq(0.1, 0.5, 0.1)) # nnLVGLASSO
    out <- nnlvg(Sig1, cvout$lambda, cvout$gamma); out2 <- nnlvg(Sig2, cvout$lambda, cvout$lambda)
  } else if(model == "rcLVGLASSO"){
    cvout <- cv.rclvg(X1, lambdas = seq(0.1, 0.5, 0.1)) # rcLVGLASSO
    out <- rclvg(Sig1, cvout$lambda, cvout$r); out2 <- rclvg(Sig2, cvout$lambda, cvout$r)
  } else if(model == "tGLASSO"){
    out <- ebic.tg(X1) # tGLASSO
    out$L <- NULL; out2 <- ebic.tg(X2)
  }
  
  if(model == "tGLASSO"){
    Fnorm = 0
  } else{
    Fnorm <- c(Fnorm, norm(out$L - out2$L, type = "F")) 
  }
  wd <- c(wd, wdist(out$S, distmat))
  sparsity <- c(sparsity, sum(out$S[upper.tri(out$S)] == 0)/sum(upper.tri(out$S)))
  
  df <- cbind(Fnorm, wd, sparsity)
  print(colMeans(df))
  write.table(df, file = paste0(model, ".txt"))
}