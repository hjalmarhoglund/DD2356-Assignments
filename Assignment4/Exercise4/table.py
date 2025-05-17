import numpy as np

NUM_RUNS = 10

PROC_VALS = [4,8,16,32,64]
NUM_PROCS_CONFS = len(PROC_VALS)

N_VALS = [1000, 2000, 4000, 8000, 16000]
NUM_N_CONFS = len(N_VALS)

def getprocconf():
    d = np.zeros((NUM_N_CONFS, NUM_RUNS))
    for i in range(NUM_N_CONFS):
        #input() # Read "N = ..."
        for j in range(NUM_RUNS):
            d[i,j] = float(input().split(' ')[3])
    return d


data = np.zeros((NUM_PROCS_CONFS, NUM_N_CONFS, NUM_RUNS))

for n in range(NUM_PROCS_CONFS):
    data[n,:,:] = getprocconf()

avg = np.zeros((NUM_PROCS_CONFS, NUM_N_CONFS))
std = np.zeros((NUM_PROCS_CONFS, NUM_N_CONFS))

for i in range(NUM_PROCS_CONFS):
    for j in range(NUM_N_CONFS):
        avg[i,j] = np.average(data[i,j,:])
        std[i,j] = np.std(data[i,j,:])

print("\\begin{table}[H]")
print("\t\\centering")
print("\t\\adjustbox{max width=\\textwidth}{")
print("\t\\begin{tabular}{ccccccccccc}")
print("\t\t\\toprule")
print("\t\t\\multicolumn{1}{c}{} & \\multicolumn{2}{c}{$N = 1000$} & \\multicolumn{2}{c}{$N = 2000$} & \\multicolumn{2}{c}{$N = 4000$} & \\multicolumn{2}{c}{$N = 8000$} & \\multicolumn{2}{c}{$N = 16000$} \\\\")
print("\t\t\\cmidrule(rl){2-3} \\cmidrule(rl){4-5} \\cmidrule(rl){6-7} \\cmidrule(rl){8-9} \\cmidrule(rl){10-11} ")
print("\t\t\\# procs & $\\mu$ (s) & $\\sigma$ & $\\mu$ (s) & $\\sigma$ & $\\mu$ (s) & $\\sigma$ & $\\mu$ (s) & $\\sigma$ & $\\mu$ (s) & $\\sigma$ \\\\ ")
print("\t\t\\cmidrule(rl){1-1} \\cmidrule(rl){2-2} \\cmidrule(rl){3-3} \\cmidrule(rl){4-4} \\cmidrule(rl){5-5} \\cmidrule(rl){6-6} \\cmidrule(rl){7-7} \\cmidrule(rl){8-8} \\cmidrule(rl){9-9} \\cmidrule(rl){10-10} \\cmidrule(rl){11-11}")

for ni in range(NUM_PROCS_CONFS):
    print(f"\t\t${PROC_VALS[ni]}$", end='')
    for n in range(NUM_N_CONFS):
        s = f" & ${avg[ni,n]:.2e}$ & ${std[ni,n]:.2e}$"
        for i in range(1,8):
            subr = "{-" + str(i) + "}"
            s = s.replace(f"e-0{i}", f" \\cdot 10^{subr}")
        print(s, end='')
    print(" \\\\")

print("\t\t\\bottomrule")
print("\t\\end{tabular}")
print("\t}")
print("\t\\caption{Running time and standard deviation for different number of processes and matrix sizes. Average of $10$ runs.}")
print("\t\\label{tab:matvectime}")
print("\\end{table}")
