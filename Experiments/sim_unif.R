# Uniform cluster sizes
sapply((paste0("../Models/", list.files("../Models/"))), source)
sapply((paste0("../Simulations/", list.files("../Simulations/"))), source)

set.seed(123)
pobs <- 200 # Number of observed variables for S
plats <- 2:8 # Number of latent variables for L
ns <- seq(50, 300, 50) # Number of observations

init_S <- 1.5 # Initial value for S
init_L <- 1.5 # Initial value for L

simtype = "unif"
runsim(simtype, pobs, plats, ns, init_S, init_L)