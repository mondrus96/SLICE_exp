#!/bin/bash

models=("SLICE" "tGLASSO" "nnLVGLASSO" "rcLVGLASSO")
#rank=(3 4 5 6)
#n=(75 150 225 300 375)
rank=(4)
n=(375)

# Loop over model
for model in "${models[@]}"; do
    # Loop over rank
    for plat in "${rank[@]}"; do
        # Loop over n
        for num in "${n[@]}"; do
            # Loop over iterations
            for i in {1..4}; do
                # Define start and end
                start=$((1 + 25 * (i - 1)))
                end=$((25 * i))

                # Submit the job to SLURM with Rscript arguments
                sbatch --job-name="${model}_plat${plat}_n${num}_batch${i}" simbatch.sh $plat $num $start $end $model
                sleep 1
            done
        done
    done
done