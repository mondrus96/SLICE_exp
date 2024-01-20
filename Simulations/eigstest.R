sapply((paste0("../Core/", list.files("../Core/"))), source)

pobs <- 150 # Number of observed variables for S
plat <- 3 # Number of latent variables for L
n <- 75 # Number of observations
simtype <- "rand"
iters <- 76:100

runsim(simtype, pobs, plat, n, iters)
