#!/bin/bash

sims=("cres" "spir" "rand")
plats=(2 2 4)

# Loop over n
for i in {0..2}; do
    sim=${sims[i]}
    plat=${plats[i]}
    # Loop over iterations
    for j in {1..10}; do
        # Define start and end
        start=$((1 + 10 * (j - 1)))
        end=$((10 * j))

        # Submit the job to SLURM with Rscript arguments
        sbatch --job-name="sim${sim}_batch${j}" --time=0-23:59 highNbatch.sh $sim $plat $start $end
        sleep 1
    done
done