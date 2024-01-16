# Uniform cluster sizes
sapply((paste0("../Models/", list.files("../Models/"))), source)
sapply((paste0("../Simulations/", list.files("../Simulations/"))), source)

set.seed(123)
pobs <- 200 # Number of observed variables for S
plat <- 2 # Number of latent variables for L
ns <- seq(50, 300, 50) # Number of observations

init_S <- 1.5 # Initial value for S

simtype <- "cres"
runsim(simtype, pobs, plat, ns, init_S, init_L)