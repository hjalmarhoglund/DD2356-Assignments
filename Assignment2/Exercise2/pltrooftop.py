import numpy as np
import matplotlib.pyplot as plt

operational_intensity = np.array([0.042, 0.042, 0.042, 0.042])  # FLOPs per byte
performance = np.array([0.5, 0.68, 0.73, 0.40]) * 10 ** 9 # FLOP/s
#operational_intensity = np.array([0.1, 1, 10, 50])  # FLOPs per byte
#performance = np.array([10, 50, 100, 150]) * 10 ** 9 # FLOP/s

memory_bandwidth = 979.2149 * 10 ** 9 # bytes/s
compute_peak = 1440 * 10 ** 9 # FLOP/s

x_for_slope = np.pow(np.linspace(0.1,10,1000),10)

memory_line = memory_bandwidth * x_for_slope

plt.figure(figsize=(8,6))
plt.loglog(operational_intensity, performance, marker='o', label='Measured Data')
plt.axhline(y=compute_peak, color='r', linestyle='--', label='Peak FLOP/s')
plt.loglog(x_for_slope, memory_line, linestyle='dotted', color='gray', label='Memory Bandwidth Limit')
plt.xlabel('Operational Intensity (FLOPs/Byte)')
plt.ylabel('Performance (FLOP/s)')
plt.title('Roofline Model')
plt.legend()
plt.grid(True, which='both', linestyle='--')
plt.show()

