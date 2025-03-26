#!/bin/bash -l
# The -l above is required to get the full environment with modules

#SBATCH -J mpi_hello_world_batch_job 
#SBATCH -t 00:02:00
#SBATCH -A edu25.DD2356
# Number of nodes
#SBATCH -p main 
#SBATCH --nodes=1
#SBATCH -e error_file.e

# Run the executable file 
srun -n 128 ./mpi_hello_world.out > mpi_hello_world_output
