#!/bin/bash
#SBATCH --job-name=shallow_water
#SBATCH --account=edu25.dd2356
#SBATCH -p main
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --output=shallow_%j.out

gcc -O3 -fopenmp main.c -o shallow

export OMP_NUM_THREADS=16

srun ./shallow
