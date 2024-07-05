sapply((paste0("../Core/", list.files("../Core/"))), source)
load("eeg_sim_params.rda")

pobs <- ncol(S_star) # Number of observed variables for S
plat <- max(Lout$z)
n <- 10000 # Number of observations
simtype <- "eeg"
iters <- 1:100

models <- c("SLICE", "SLICE_GSCAD", "SLICE_CLIME", "nnLVGLASSO", "tGLASSO", "rcLVGLASSO")
for(i in models){
  runsim(simtype, i, pobs, plat, n, iters)
}