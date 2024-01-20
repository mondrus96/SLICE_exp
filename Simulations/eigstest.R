sapply((paste0("../Core/", list.files("../Core/"))), source)

pobs <- 150 # Number of observed variables for S
plat <- 5 # Number of latent variables for L
n <- 150 # Number of observations
simtype <- "rand"
iters <- 24

runsim(simtype, pobs, plat, n, iters)