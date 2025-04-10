import numpy as np
import matplotlib.pyplot as plt

# Example data (replace with actual measurements)
operational_intensity = np.array([0.1, 1, 10, 50])  # FLOPs per byte
performance = np.array([10, 50, 100, 150])  # GFLOP/s

plt.figure(figsize=(8,6))
plt.loglog(operational_intensity, performance, marker='o', label='Measured Data')
plt.axhline(y=200, color='r', linestyle='--', label='Peak FLOP/s')
# Use workaround from here
# https://stackoverflow.com/questions/77964915/how-to-use-axline-with-log-scale-axis
x_for_slope = np.linspace(0.1,10,1000)
plt.plot(x_for_slope, x_for_slope, linestyle='dotted', color='gray', label='Memory Bandwidth Limit')
plt.xlabel('Operational Intensity (FLOPs/Byte)')
plt.ylabel('Performance (GFLOP/s)')
plt.title('Roofline Model')
plt.legend()
plt.grid(True, which='both', linestyle='--')
plt.show()

