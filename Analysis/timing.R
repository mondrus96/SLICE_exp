# Timing comparison
sapply((paste0("../Core/", list.files("../Core/"))), source)
pobs <- 150
plat <- 7
n <- 10000

set.seed(123)
S_star <- Smat(150, 2, 1.5); S_star[S_star < 0.01] <- 0 # Sparse
Lout <- Lrand(pobs, plat, 1.5); L_star <- Lout$L # Latent
Sigma_star <- solve(S_star + L_star) # True Sigma
X <- mvrnorm(n, rep(0, pobs), Sigma = Sigma_star) # Finite sample data
Sigma <- cov(X) # Sample Sigma

cvsli <- cv.slice(X, rs = 7)
cvrclvg <- cv.rclvg(X, rs = 7)
cvnnlvg <- cv.nnlvg(X)

startT = Sys.time() # SLICE
slice(Sigma, cvsli$lambda, 7)
endT = Sys.time()
print(paste("SLICE:", endT - startT))

startT = Sys.time() # rcLVGLASSO
rclvg(Sigma, cvrclvg$lambda, 7)
endT = Sys.time()
print(paste("rcLVGLASSO:", endT - startT))

startT = Sys.time() # nnLVGLASSO
nnlvg(Sigma, cvnnlvg$lambda, cvnnlvg$gamma)
endT = Sys.time()
print(paste("nnLVGLASSO:", endT - startT))