import numpy as np
import matplotlib.pyplot as plt

NUM_WARMUP = 3
NUM_REAL = 10

versions = ["DFT_serial", "DFT_omp", "DFT_omp_schedule", "DFT_omp_reduction", "DFT_swap", "DFT_manual"]
NUM_VERSIONS = len(versions)

threads = [1,32,64,128]
NUM_THREADS = len(threads)

data = np.zeros((NUM_THREADS, NUM_VERSIONS, NUM_REAL))

for tr in range(NUM_THREADS):
    # Read first line
    input()

    for ver in range(NUM_VERSIONS):
        for _ in range(NUM_WARMUP):
            input()
            input() # Read test output
        for i in range(NUM_REAL):
            t = float(input().split(' ')[4])
            data[tr,ver,i] = t
            input() # Read test output

plt.figure()
base = []

for ver in range(NUM_VERSIONS):
    base.append(np.mean(data[0,ver,:]))

base = np.array(base)

for ver in range(NUM_VERSIONS):
    vals = np.mean(data[:,ver,:], axis=1)
    speedup = base[ver] / vals
    plt.plot(threads,speedup,'o-',label=versions[ver])

plt.legend()
plt.title("Average speedup for N=10000")
plt.xlabel("Number of threads")
plt.ylabel("Speedup factor")
plt.xticks(threads)
plt.savefig("DFTW-speedup.png")
#plt.grid(True, which='both', linestyle='--')
plt.show()
