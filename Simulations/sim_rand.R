sapply((paste0("../Core/", list.files("../Core/"))), source)

pobs <- 150 # Number of observed variables for S
plat <- 4 # Number of latent variables for L
n <- 10000 # Number of observations
simtype <- "rand"
iters <- 1:100

runsim(simtype, pobs, plat, n, iters)