sapply((paste0("../Core/", list.files("../Core/"))), source)

pobs <- 150 # Number of observed variables for S
n <- 10000 # Number of observations
simtype <- "cres"
iters <- 1:100

runsim(simtype, pobs, plat, n, iters)