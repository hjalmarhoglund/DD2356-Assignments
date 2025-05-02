import numpy as np

NUM_WARMUP = 3
NUM_REAL = 10

versions = ["\\texttt{DFT\\_serial}", "\\texttt{DFT\\_omp}", "\\texttt{DFT\\_omp\\_schedule}", "\\texttt{DFT\\_omp\\_reduction}", "\\texttt{DFT\\_swap}", "\\texttt{DFT\\_manual}"]
NUM_VERSIONS = len(versions)

data = np.zeros((NUM_VERSIONS, NUM_REAL))

# Read first line
input()

for ver in range(NUM_VERSIONS):
    for _ in range(NUM_WARMUP):
        input()
        input() # Read test output
    for i in range(NUM_REAL):
        t = float(input().split(' ')[4])
        data[ver,i] = t
        input() # Read test output

s = ''
for ver in range(NUM_VERSIONS):
    s += f'\t\t{versions[ver]} & ${np.mean(data[ver,:]):.2e}$ & ${np.std(data[ver,:]):.2e}$ \\\\\n'

for i in range(10):
    subr = ' \\cdot 10^{' + str(-i) + '}'
    s = s.replace(f'e-0{i}', subr)
for i in range(10):
    subr = ' \\cdot 10^{' + str(i) + '\\phantom{-}}'
    s = s.replace(f'e+0{i}', subr)

print(s, end='')

