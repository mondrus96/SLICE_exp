# Uniform cluster sizes
sapply((paste0("../Models/", list.files("../Models/"))), source)
sapply((paste0("../Simulations/", list.files("../Simulations/"))), source)

set.seed(123)
pobs <- 150 # Number of observed variables for S
plat <- 4 # Number of latent variables for L
n <- 10000 # Number of observations
simtype <- "rand"
iters <- 1:100

runsim(simtype, pobs, plat, n, iters)
