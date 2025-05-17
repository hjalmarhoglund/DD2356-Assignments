#!/bin/bash -l
#SBATCH --job-name=scaling_test
#SBATCH --account=edu25.dd2356
#SBATCH --partition=main
#SBATCH --nodes=1             # one node can host 32 ranks
#SBATCH --time=00:05:00
#SBATCH --output=scaling_run_%j.txt

# build files
cc -O3 original.c -o original.out
cc -O3 halo.c     -o halo.out

# helper function
run_and_time () {
    ranks=$1
    exe=$2
    echo ">>> $exe on $ranks rank(s)"
    srun -n "$ranks" "./$exe"            \
        | grep "Time ="                  \
        | awk -v e="$exe" -v p="$ranks" '{print e, p, $3}' \
        >> times.txt
}

echo "#exe  ranks  seconds" > times.txt

# serial baseline
run_and_time 1 original.out

# parallel runs
for ranks in 2 4 8 16 32; do
    run_and_time "$ranks" halo.out
done

# results
echo "Finished.  Collected timings:"
cat times.txt
