#!/bin/bash
models=("SLICE" "SLICE_GSCAD" "SLICE_CLIME" "nnLVGLASSO" "rcLVGLASSO")
simtypes=("cres" "spir" "rand" "eeg")

# Loop over models
for model in models; do
    # Loop over simtypes
    for simtype in simtypes; do
        # Submit the job to SLURM with Rscript arguments
        sbatch --job-name="sim${simtype}_model${model}" --time=0-11:59 highNbatch.sh $model $simtype
        sleep 1
    done
done