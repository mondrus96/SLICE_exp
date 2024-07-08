#!/bin/bash
#SBATCH --account=def-ijc               # account associated with job
#SBATCH --nodes=1                       # number of node MUST be 1
#SBATCH --cpus-per-task=1               # number of processes
#SBATCH --mem=0                         # memory usage
#SBATCH --mail-user=mondrus@ualberta.ca # Send email updates to you or someone else
#SBATCH --mail-type=ALL                 # Send an email in all cases (job started, job ended, job aborted)

module load gcc/9.3.0 r/4.3.1
export R_LIBS=~/local/R_libs/

echo "Received arguments: $1 $2"

Rscript -e "
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
" $1 $2