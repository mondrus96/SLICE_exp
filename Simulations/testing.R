sapply((paste0("../Core/", list.files("../Core/"))), source)

pobs <- 100 # Number of observed variables for S
plat <- 4 # Number of latent variables for L
n <- 100000 # Number of observations

set.seed(12345)

S_star <- Smat(pobs, 2, 1)
S_star[S_star < 0.01] <- 0

Lout <- Lrand(pobs, plat, 1.5)
L_star <- Lout$L; z_star <- Lout$z

Sigma_star <- solve(S_star + L_star)
X <- mvrnorm(n, rep(0, pobs), Sigma_star)

cvsli <- cv.slice(X)
cvlvg <- cv.lvg(X)

sliout <- slice(cov(X), cvsli$lambda, cvsli$r)
lvgout <- lvglasso(cov(X), cvlvg$lambda, cvlvg$gamma)

# Function for calculating F1 score
F1score <- function(true, pred) {
  TP <- sum(true == 1 & pred == 1)
  FP <- sum(true == 0 & pred == 1)
  FN <- sum(true == 1 & pred == 0)
  
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  return(2 * precision * recall / (precision + recall))
}

true <- 1*(S_star[upper.tri(S_star)] != 0)
sli <- 1*(sliout$S[upper.tri(sliout$S)] != 0)
lvg <- 1*(lvgout$S[upper.tri(lvgout$S)] != 0)

F1score(true, sli)
F1score(true, lvg)
