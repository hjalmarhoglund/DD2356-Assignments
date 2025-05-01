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

numactl --show
echo "Cores available: $(nproc)"

# compile
gcc -O3 -fopenmp -DSTREAM_ARRAY_SIZE=20000000 stream.c -o stream

export OMP_NUM_THREADS=128
schedules=(static dynamic guided)

for schedule in "${schedules[@]}"; do
  export OMP_SCHEDULE=$schedule
  > results_128_${schedule}.txt

  for i in {1..5}; do
    srun --cpus-per-task=128 ./stream | awk '/Copy:/ {print $2}' >> results_128_${schedule}.txt
  done
done
