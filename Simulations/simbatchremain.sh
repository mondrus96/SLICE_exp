#!/bin/bash

# Define all the remaining ones
reamining=(
    "3,300,1"
    "3,225,4"
    "3,375,2"
    "4,225,3"
    "6,225,1"
    "5,375,3"
    "4,300,2"
    "6,375,1"
    "4,225,4"
    "3,75,4"
    "6,150,2"
    "6,300,3"
    "4,300,3"
    "6,300,4"
    "5,225,1"
    "5,150,1"
)

# Loop over rank
for params in "${reamining[@]}"; do
    # Define params
    IFS=',' read -r plat num i <<< "$params"
    echo "Processing: $plat, $num, $i"

    # Define start and end
    start=$((1 + 25 * (i - 1)))
    end=$((25 * i))

    # Submit the job to SLURM with Rscript arguments
    sbatch --job-name="plat${plat}_n${num}_batch${i}" simbatch.sh $plat $num $start $end
    sleep 1
done