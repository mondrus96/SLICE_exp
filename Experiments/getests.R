# Load functions
sapply((paste0("../Core/", list.files("../Core/"))), source)

# Load data
visit1 <- read.table("gordonvisit2.txt"); visit2 <- read.table("gordonvisit3.txt")
subj <- 25; tlen <- 197

# Set seed
set.seed(123)

# Loop through subjects
allS <- allL <- vector("list", subj)
for(i in 1:subj){
  print(i) # Current subject
  
  inds <- (1+((i-1)*tlen)):(i*tlen) # Define the first and second visits for subject
  X1 <- as.matrix(visit1[inds, 2:ncol(visit1)]); X2 <- as.matrix(visit2[inds, 2:ncol(visit2)])
  Sig1 <- cov(X1); Sig2 <- cov(X2); Xall <- rbind(X1, X2)
  
  cvsli <- cv.slice(Xall) # SLICE
  sli1 <- slice(Sig1, cvsli$lambda, cvsli$r); sli2 <- slice(Sig2, cvsli$lambda, cvsli$r)
  
  cvnnlvg <- cv.nnlvg(Xall) # nnLVGLASSO
  nnlvg1 <- nnlvg(Sig1, cvnnlvg$lambda, cvnnlvg$gamma); nnlvg2 <- nnlvg(Sig2, cvnnlvg$lambda, cvnnlvg$gamma)
  
  cvrclvg <- cv.rclvg(Xall) # rcLVGLASSO
  rclvg1 <- rclvg(Sig1, cvrclvg$lambda, cvrclvg$r); rclvg2 <- rclvg(Sig2, cvrclvg$lambda, cvrclvg$r)
  
  V1 <- list(sli1$S, nnlvg1$S, rclvg1$S); V2 <- list(sli2$S, nnlvg2$S, rclvg2$S)
  names(V1) <- names(V2) <- c("SLICE", "nnLVGLASSO", "rcLVGLASSO")
  allS[[i]] <- list(V1, V2)
  
  V1 <- list(sli1$L, nnlvg1$L, rclvg1$L); V2 <- list(sli2$L, nnlvg2$L, rclvg2$L)
  names(V1) <- names(V2) <- c("SLICE", "nnLVGLASSO", "rcLVGLASSO")
  allL[[i]] <- list(V1, V2)
  
  save(allS, allL, file = "ests.rda") # Save results
}