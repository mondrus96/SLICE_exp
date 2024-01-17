# Uniform cluster sizes
sapply((paste0("../Models/", list.files("../Models/"))), source)
sapply((paste0("../Simulations/", list.files("../Simulations/"))), source)

pobs <- 150 # Number of observed variables for S
n <- 10000 # Number of observations
simtype <- "circ"
iters <- 1:100

runsim(simtype, pobs, plat, n, iters)