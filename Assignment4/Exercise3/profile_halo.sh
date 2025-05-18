#!/bin/bash -l
#SBATCH --job-name=profile_halo
#SBATCH --account=edu25.dd2356
#SBATCH --partition=main
#SBATCH --nodes=2             
#SBATCH --ntasks-per-node=16
#SBATCH --time=00:05:00
#SBATCH --output=prof_%j.txt

module purge
module load PDC/23.12 
module load score-p/8.4-cpeGNU   

scorep mpicc -O3 ex1_halo.c -o ex1_halo.out

export SCOREP_ENABLE_PROFILING=true
export SCOREP_ENABLE_TRACING=false

srun ./ex1_halo.out
