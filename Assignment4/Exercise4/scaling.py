import numpy as np
import matplotlib.pyplot as plt

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

# Strong scaling
plt.figure()
plt.plot(np.array(PROC_VALS), avg[:,-1])
plt.xlabel("Number of processes")
plt.xticks(np.array(PROC_VALS))
plt.ylabel("Running time (seconds)")
plt.title("Strong scaling of Matrix-Vector multiplication")
plt.savefig("matvecstrong.png")
