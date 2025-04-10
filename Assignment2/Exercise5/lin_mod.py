import numpy as np
import matplotlib.pyplot as plt
# replace with actual data
message_sizes = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024])  # in Bytes
rtt_times = np.array([10, 12, 16, 22, 30, 42, 55, 70, 90, 120, 160])  # in microseconds

coeffs = np.polyfit(message_sizes, rtt_times, 1)
latency_extrapolated = coeffs[1]  # y-intercept

plt.plot(message_sizes, rtt_times, 'o', label='Measured RTT')
plt.plot(message_sizes, coeffs[0]*message_sizes + coeffs[1], '--', label='Linear Fit')
plt.xlabel('Message Size (Bytes)')
plt.ylabel('Round-Trip Time (microseconds)')
plt.title(f'Extrapolated Latency at Zero Message Size: {latency_extrapolated:.2f} µs')
plt.legend()
plt.show()
