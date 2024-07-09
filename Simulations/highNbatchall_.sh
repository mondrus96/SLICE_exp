#!/bin/bash
models=("SLICE_CLIME" "rcLVGLASSO")
simtypes=("cres" "spir" "rand" "eeg")
iters=({1..10})  # This creates a sequence from 1 to 10

# Loop over models
for model in "${models[@]}"; do
    # Loop over simtypes
    for simtype in "${simtypes[@]}"; do
        # Loop over iterations
        for iter in "${iters[@]}"; do
            # Submit the job to SLURM with Rscript arguments
            sbatch --job-name="sim${simtype}_model${model}_iter${iter}" --time=0-11:59 highNbatch.sh $model $simtype $iter
            sleep 1
        done
    done
done