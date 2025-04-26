import numpy as np

NUM_WARMUP = 3
NUM_REAL = 10

versions = ['serial', 'omp_sum', 'omp_critical', 'omp_local_sum', 'omp_reduction']
NUM_VERSIONS = len(versions)

threads = [1,2,4,8,16,20,24,28,32,64,128]
NUM_THREAD_CONFS = len(threads)

# Set up matrix to store the output
data = np.zeros((NUM_VERSIONS, NUM_THREAD_CONFS, NUM_REAL))

# For each number of threads
for t in range(NUM_THREAD_CONFS):
    input() # Read "SETTING OMP=..."
    for ver in range(NUM_VERSIONS):
        input() # Read "=== START ... ==="
        for _ in range(NUM_WARMUP):
            input()
        for r in range(NUM_REAL):
            v = float(input().split(' ')[3])
            data[ver,t,r] = v
        input() # Read "=== END ... ==="

# Print for N in 1 ... 32
for ver in range(3):
    print("\t\\begin{tabular}{|c|c|c|}")
    print("\t\t\\hline\n\t\t\\# thrds & $\\mu$ (s) & $\\sigma$ \\\\\n\t\t\\hline")
    for n in range(9):
        s = f"\t\t${threads[n]}$ & ${np.mean(data[ver,n,:]):.2e}$ & ${np.std(data[ver,n,:]):.2e}$ \\\\"
        for i in range(1,8):
            subr = "{-" + str(i) + "}"
            s = s.replace(f"e-0{i}", f" \\cdot 10^{subr}")
        print(s)
    print("\t\t\\hline")
    print("\t\\end{tabular}\n\t\\quad")

