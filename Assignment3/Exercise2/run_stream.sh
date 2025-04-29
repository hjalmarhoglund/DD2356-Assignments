#!/bin/bash
#SBATCH --account=edu25.dd2356
#SBATCH --job-name=stream_bm
#SBATCH --output=stream_bm_%j.out
#SBATCH -p main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=00:20:00
#SBATCH --mem=100G

# compile
gcc -O3 -fopenmp -DSTREAM_ARRAY_SIZE=20000000 stream.c -o stream

threads=(1 32 64 128)

for t in "${threads[@]}"; do
  export OMP_NUM_THREADS=$t
  > results_${t}.txt
  for i in {1..5}; do
    srun --cpus-per-task=$t ./stream | grep Copy: | awk '{print $2}' >> results_${t}.txt
  done
done
