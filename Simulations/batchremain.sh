#!/bin/bash

# Read the file line by line
while IFS=" " read -r model plat num i
do
    # Process each column as needed
    echo "Submitting: $model, $plat, $num, $i"

    # Define start and end
    start=$((1 + 25 * (i - 1)))
    end=$((25 * i))

    # Submit the job to SLURM with Rscript arguments
    sbatch --job-name="${model}_plat${plat}_n${num}_batch${i}" --time=0-23:00 simbatch.sh $plat $num $start $end $model
    sleep 1
done < "remain.txt"