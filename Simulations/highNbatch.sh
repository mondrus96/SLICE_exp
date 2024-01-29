#!/bin/bash
#SBATCH --account=def-ijc               # account associated with job
#SBATCH --nodes=1                       # number of node MUST be 1
#SBATCH --cpus-per-task=1               # number of processes
#SBATCH --mem=0                         # memory usage
#SBATCH --mail-user=mondrus@ualberta.ca # Send email updates to you or someone else
#SBATCH --mail-type=ALL                 # Send an email in all cases (job started, job ended, job aborted)

module load gcc/9.3.0 r/4.3.1
export R_LIBS=~/local/R_libs/

echo "Received arguments: $1 $2 $3 $4"

Rscript -e "
sapply((paste0('../Core/', list.files('../Core/'))), source)

args <- commandArgs(trailingOnly = TRUE)

pobs <- 150 # Number of observed variables for S
n <- 10000 # Number of observations
simtype <- as.character(args[1])
plat <- as.integer(args[2])
start <- as.integer(args[3])
end <- as.integer(args[4])

runsim(simtype, 'rcLVGLASSO', pobs, plat, n, start:end) 
" $1 $2 $3 $4