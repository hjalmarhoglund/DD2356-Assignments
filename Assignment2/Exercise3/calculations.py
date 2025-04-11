nrows = [1e2, 1e4, 1e6, 1e8]
nnzs = [460, 49600, 4996000, 499960000]
exectimes = [0.000001, 0.000058, 0.007771, 0.806055]

clock_freq = 2_250_000_000

print("Expected vs real running time")
for i in range(4):
    print(f"{2 * nnzs[i] / clock_freq:.2e} VS {exectimes[i]:.2e}")

print("Read bandwidth (MB/s)")
for i in range(4):
    print(f"{8*(2.5 * nnzs[i] + 0.5 * nrows[i]) / (1e6 * exectimes[i]):.2e}")

