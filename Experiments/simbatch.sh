#!/bin/bash
#SBATCH --account=def-ijc               # account associated with job
#SBATCH --nodes=1                       # number of node MUST be 1
#SBATCH --cpus-per-task=1               # number of processes
#SBATCH --mem=0                         # memory usage
#SBATCH --time=0-02:59                  # time (DD-HH:MM)
#SBATCH --mail-user=mondrus@ualberta.ca # Send email updates to you or someone else
#SBATCH --mail-type=ALL                 # Send an email in all cases (job started, job ended, job aborted)

module load gcc/9.3.0 r/4.3.1
export R_LIBS=~/local/R_libs/

Rscript -e "
# Uniform cluster sizes
sapply((paste0('../Models/', list.files('../Models/'))), source)
sapply((paste0('../Simulations/', list.files('../Simulations/'))), source)

simtype <- 'rand'
pobs <- 150 # Number of observed variables for S
args <- commandArgs(trailingOnly = TRUE)
plat <- as.integer(args[1])
n <- as.integer(args[2])
start <- as.integer(args[3])
end <- as.integer(args[4])
iters <- start:end

runsim(simtype, pobs, plat, n, iters)
" $plat $n $start $end