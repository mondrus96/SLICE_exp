#!/bin/bash
#SBATCH --account=def-ijc               # account associated with job
#SBATCH --nodes=1                       # number of node MUST be 1
#SBATCH --cpus-per-task=1               # number of processes
#SBATCH --mem=0                         # memory usage
#SBATCH --time=0-00:30                  # time (DD-HH:MM)
#SBATCH --mail-user=mondrus@ualberta.ca # Send email updates to you or someone else
#SBATCH --mail-type=ALL                 # Send an email in all cases (job started, job ended, job aborted)

module load gcc/9.3.0 r/4.3.1
export R_LIBS=~/local/R_libs/

export n=1000
export plat=4
export start=1
export end=25

Rscript -e "
# Uniform cluster sizes
sapply((paste0("../Models/", list.files("../Models/"))), source)
sapply((paste0("../Simulations/", list.files("../Simulations/"))), source)

pobs <- 150 # Number of observed variables for S
plat <- as.integer(Sys.getenv('plat')) # Number of latent variables for L
n <- as.integer(Sys.getenv('n')) # Number of observations
simtype <- "rand"
start <- as.integer(Sys.getenv('start'))
end <- as.integer(Sys.getenv('end'))
iters <- start:end

runsim(simtype, pobs, plat, n, iters)
"