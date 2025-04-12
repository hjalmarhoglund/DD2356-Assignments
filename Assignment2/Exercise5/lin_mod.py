import numpy as np
import matplotlib.pyplot as plt
# replace with actual data
message_sizes = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024])  # in Bytes
rtt_times = np.array([5.4302, 4.8921, 4.9133, 4.9492, 5.0105, 6.2017, 6.3169, 6.369, 6.6887, 7.4462, 8.0062]) # in microseconds

coeffs = np.polyfit(message_sizes, rtt_times, 1)
latency_extrapolated = coeffs[1]  # y-intercept

plt.plot(message_sizes, rtt_times, 'o', label='Measured RTT')
plt.plot(message_sizes, coeffs[0]*message_sizes + coeffs[1], '--', label='Linear Fit')
plt.xlabel('Message Size (Bytes)')
plt.ylabel('Round-Trip Time (microseconds)')
plt.title(f'Extrapolated Latency at Zero Message Size: {latency_extrapolated:.2f} µs')
plt.legend()
plt.show()
