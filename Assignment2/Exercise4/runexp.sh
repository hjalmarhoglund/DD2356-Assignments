#!/bin/bash

perf stat -e cycles,instructions,duration_time,L1-dcache-loads,L1-dcache-load-misses,L1-dcache-stores ./matmul64N.out > res64N.txt
perf stat -e cycles,instructions,duration_time,L1-dcache-loads,L1-dcache-load-misses,L1-dcache-stores ./matmul64O.out > res64O.txt
perf stat -e cycles,instructions,duration_time,L1-dcache-loads,L1-dcache-load-misses,L1-dcache-stores ./matmul1000N.out > res1000N.txt
perf stat -e cycles,instructions,duration_time,L1-dcache-loads,L1-dcache-load-misses,L1-dcache-stores ./matmul1000O.out > res1000O.txt

