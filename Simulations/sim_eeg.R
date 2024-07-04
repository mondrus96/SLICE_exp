sapply((paste0("../Core/", list.files("../Core/"))), source)

pobs <- 5000 # Number of observed variables for S
#n <- 10000 # Number of observations
simtype <- "cres"

S_star <- Smat(pobs, 2, 1.5)
S_star[S_star < 0.01] <- 0 # True sparse component

if (simtype == "exp"){
  Lout <- Lexp(pobs, plat, 1.5)
} else if (simtype == "rand"){
  Lout <- Lrand(pobs, plat, 1.5)
} else if (simtype == "cres"){
  Lout <- Lcres(pobs, 0.1)
  plat <- 2
} else if (simtype == "circ"){
  Lout <- Lcirc(pobs, 0.05)
  plat <- 2
} else if (simtype == "spir"){
  Lout <- Lspir(pobs, 0.05)
  plat <- 2
}

L_star <- Lout$L; z_star <- Lout$z # True latent component; True cluster labels
Sigma_star <- Matrix::chol2inv(Matrix::chol(S_star + L_star)) # True Sigma

start <- Sys.time()
sli <- slice(Sigma_star, 0.05, 2, Sest = "huge_glasso")
end <- Sys.time()
end - start
sum(sli$S[upper.tri(sli$S)] != 0)