import numpy as np
import matplotlib.pyplot as plt

schedules = ['static', 'dynamic', 'guided']
means = []
stds = []

for sched in schedules:
    data = np.loadtxt(f"results_128_{sched}.txt")
    means.append(data.mean())
    stds.append(data.std(ddof=1))

x = np.arange(len(schedules))

plt.errorbar(x, means, yerr=stds, fmt='o-', capsize=5)
plt.xticks(x, schedules)
plt.xlabel('OpenMP Schedule')
plt.ylabel('Copy Bandwidth (MB/s)')
plt.title('STREAM Copy Bandwidth - Different Schedules (128 threads)')
plt.grid(True, linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.savefig('stream_benchmark_schedules_plot.png')
plt.show()
