sapply((paste0("../Core/", list.files("../Core/"))), source)

pobs <- 150 # Number of observed variables for S
n <- 10000 # Number of observations
simtype <- "spir"
iters <- 1:100

models <- c("SLICE", "nnLVGLASSO", "tGLASSO", "rcLVGLASSO")
for(i in models){
  runsim(simtype, i, pobs, plat, n, iters)
}