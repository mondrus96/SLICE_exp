#!/bin/bash
#SBATCH --account=def-ijc               # account associated with job
#SBATCH --nodes=1                       # number of node MUST be 1
#SBATCH --cpus-per-task=1               # number of processes
#SBATCH --mem=0                         # memory usage
#SBATCH --time=0-23:59                  # time (DD-HH:MM)
#SBATCH --mail-user=mondrus@ualberta.ca # Send email updates to you or someone else
#SBATCH --mail-type=ALL                 # Send an email in all cases (job started, job ended, job aborted)

module load StdEnv/2023 gcc/13.3 r/4.4.0
export R_LIBS=~/local/R_libs/
R CMD BATCH --no-save --no-restore sim_eeg.R