import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("row_sums_output.txt")
plt.figure()
plt.plot(data, marker='o', label='Row Sums')
plt.xlabel("Row Index")
plt.ylabel("Sum Value")
plt.title("Parallel Row Sum Computation with MPI")
plt.legend()
plt.savefig('par.png')

data = np.loadtxt("ser_row_sums_output.txt")
plt.figure()
plt.plot(data, marker='o', label='Row Sums')
plt.xlabel("Row Index")
plt.ylabel("Sum Value")
plt.title("Serial Row Sum Computation")
plt.legend()
plt.savefig('serial.png')
