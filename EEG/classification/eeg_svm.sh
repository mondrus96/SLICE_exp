#!/bin/bash
#SBATCH --account=def-ijc               # account associated with job
#SBATCH --nodes=1                       # number of node MUST be 1
#SBATCH --cpus-per-task=32              # number of processes
#SBATCH --mem=0                         # memory usage
#SBATCH --mail-user=mondrus@ualberta.ca # Send email updates to you or someone else
#SBATCH --mail-type=ALL                 # Send an email in all cases (job started, job ended, job aborted)
#SBATCH --time=0-11:59                  # Total run time 

module load StdEnv/2023
module load gcc/12.3 r/4.3.1
export R_LIBS=~/local/R_libs/

Rscript eeg_svm.R