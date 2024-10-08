#!/bin/bash
#SBATCH --account=def-ijc               # account associated with job
#SBATCH --nodes=1                       # number of node MUST be 1
#SBATCH --cpus-per-task=1               # number of processes
#SBATCH --mem=0                         # memory usage
#SBATCH --mail-user=mondrus@ualberta.ca # Send email updates to you or someone else
#SBATCH --mail-type=ALL                 # Send an email in all cases (job started, job ended, job aborted)

module load gcc/9.3.0 r/4.3.1
export R_LIBS=~/local/R_libs/

echo "Received arguments: $1 $2 $3 $4 $5"

Rscript -e "
sapply((paste0('../Core/', list.files('../Core/'))), source)

simtype <- 'rand'
pobs <- 150 # Number of observed variables for S
args <- commandArgs(trailingOnly = TRUE)
plat <- as.integer(args[1])
n <- as.integer(args[2])
start <- as.integer(args[3])
end <- as.integer(args[4])
model <- as.character(args[5])

print(paste('plat:', plat))
print(paste('n:', n))
print(paste('start:', start))
print(paste('end:', end))

iters <- start:end

runsim(simtype, model, pobs, plat, n, iters)
" $1 $2 $3 $4 $5