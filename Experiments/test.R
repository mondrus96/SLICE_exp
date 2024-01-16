# Uniform cluster sizes
sapply((paste0("../Models/", list.files("../Models/"))), source)
sapply((paste0("../Simulations/", list.files("../Simulations/"))), source)

pobs <- 50 # Number of observed variables for S
plats <- 2 # Number of latent variables for L
ns <- seq(50, 300, 50) # Number of observations

init_S <- 1.5 # Initial value for S

simtype <- "circ"

S_star <- Smat(pobs, 2, init_S)
Lout <- Lcirc(pobs)
L_star <- Lout$L

Sigma = solve(S_star + L_star)
test1 = slice(Sigma, 0.01, 4)
test2 = slice(Sigma, 0.01, 4, start = "zero")
