import numpy as np

NUM_WARMUP = 3
NUM_REAL = 10

versions = ['serial', '\\texttt{omp\\_sum}', '\\texttt{omp\\_critical}']
NUM_VERSIONS = len(versions)

threads = [1,2,4,8,16,20,24,28,32]
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
handle_first = NUM_VERSIONS
print("\t\\adjustbox{max width=\\textwidth}{")
print("\t\\begin{tabular}{ccccccc}")
print("\t\t\\toprule")
print("\t\t\\multicolumn{1}{c}{}", end='')
for ver in range(handle_first):
    s = " & \\multicolumn{2}{c}{" + versions[ver] + '}'
    print(s, end='')
print(" \\\\")
print("\t\t", end='')
for i in range(1,handle_first+1):
    s = " \\cmidrule(rl){" + str(2*i) + "-" + str(2*i+1) + "}"
    print(s, end='')
print()
print("\t\t\\# threads", end='')
for _ in range(handle_first):
    print(" & $\\mu$ (s) & $\\sigma$", end='')
print(" \\\\")
print("\t\t", end='')
for i in range(1,handle_first*2+2):
    s = "\\cmidrule(rl){" + str(i) + '-' + str(i) + "}"
    print(s, end=' ')
print()
for n in range(9):
    print(f"\t\t${threads[n]}$", end='')
    for ver in range(handle_first):
        s = f" & ${np.mean(data[ver,n,:]):.2e}$ & ${np.std(data[ver,n,:]):.2e}$"
        for i in range(1,8):
            subr = "{-" + str(i) + "}"
            s = s.replace(f"e-0{i}", f" \\cdot 10^{subr}")
        print(s, end='')
    print(" \\\\")
print("\t\t\\bottomrule")
print("\t\\end{tabular}\n\t}")

