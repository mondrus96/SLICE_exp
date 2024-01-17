# Uniform cluster sizes
sapply((paste0("../Models/", list.files("../Models/"))), source)
sapply((paste0("../Simulations/", list.files("../Simulations/"))), source)

set.seed(123)
pobs <- 150 # Number of observed variables for S
plat <- 4 # Number of latent variables for L
n <- "inf" # Number of observations
simtype <- "exp"

runsim(simtype, pobs, plat, n, 1:50)