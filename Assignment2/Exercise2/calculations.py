import numpy as np

times = np.array([
[0.002025, 0.002029, 0.002029, 0.002015, 0.002005, 0.002033, 0.001983, 0.002017, 0.002001, 0.002012],
[0.014617, 0.014838, 0.014754, 0.014563, 0.014514, 0.014579, 0.014546, 0.014763, 0.014653, 0.014472],
[0.136685, 0.137124, 0.136628, 0.137246, 0.136997, 0.136434, 0.136941, 0.136398, 0.136785, 0.136634],
[2.858346, 2.552086, 2.517934, 2.460361, 2.442166, 2.432661, 2.457380, 2.405495, 2.418833, 2.446097],
])

avgtimes = np.mean(times, 1)

print("AVERAGE (seconds)")
print(avgtimes)
print("STANDARD DEVIATION")
print(np.std(times, 1))



N = np.array([1e6, 1e7, 1e8, 1e9])
print("FLOP/s")
print(N / avgtimes)
print("OI")
print(N / (24 * N))
