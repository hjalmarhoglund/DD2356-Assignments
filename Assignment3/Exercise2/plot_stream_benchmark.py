import numpy as np
import matplotlib.pyplot as plt

threads = [1, 32, 64, 128]
means = []
stds = []

for t in threads:
    data = np.loadtxt(f"results_{t}.txt")
    means.append(data.mean())
    stds.append(data.std())

plt.errorbar(threads, means, yerr=stds, fmt='o-')
plt.xlabel('Number of Threads')
plt.ylabel('Copy Bandwidth (MB/s)')
plt.title('STREAM Copy Bandwidth vs Threads')
plt.grid(True)
plt.savefig('stream_benchmark_plot.png')
plt.show()
