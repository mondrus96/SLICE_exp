args <- commandArgs(trailingOnly = TRUE)
model <- as.character(args[1])
simtype <- as.character(args[2])

sapply((paste0("../Core/", list.files("../Core/"))), source)
load("eeg_sim_params.rda")

if(simtype == "eeg"){
  pobs <- ncol(S_star) # Number of observed variables for S
  plat <- max(Lout$z)
} else if(simtype == "rand"){
  pobs <- 150
  plat <- 4
} else{
  pobs <- 150
  plat <- 2
}
n <- 10000 # Number of observations
iters <- 1:100

runsim(simtype, model, pobs, plat, n, iters)