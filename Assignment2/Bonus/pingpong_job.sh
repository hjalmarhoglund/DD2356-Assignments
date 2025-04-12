#!/bin/bash -l
#SBATCH -A edu25.dd2356
#SBATCH -p main
#SBATCH -t 00:05:00   
#SBATCH -N 2         
#SBATCH -n 2        
#SBATCH -J pingpong
#SBATCH --ntasks-per-node=1

srun -N2 --ntasks=2 ./mpi_ping_pong
