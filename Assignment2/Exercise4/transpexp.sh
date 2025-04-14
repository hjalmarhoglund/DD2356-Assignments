#!/bin/bash

perf stat -e cycles,instructions,duration_time,L1-dcache-loads,L1-dcache-load-misses,L1-dcache-stores ./transp64.out > restransp64.txt
perf stat -e cycles,instructions,duration_time,L1-dcache-loads,L1-dcache-load-misses,L1-dcache-stores ./transp128.out > restransp128.txt
perf stat -e cycles,instructions,duration_time,L1-dcache-loads,L1-dcache-load-misses,L1-dcache-stores ./transp2048.out > restransp2048.txt

