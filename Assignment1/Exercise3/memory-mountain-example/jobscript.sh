#!/bin/bash
#SBATCH -A edu25.dd2356
#SBATCH -J memoryMountainTest
#SBATCH -t 0:10:00
#SBATCH -p main
#SBATCH -N 1
srun -n 1 ./mountain.out > results.txt
