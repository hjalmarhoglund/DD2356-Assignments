import numpy as np
import matplotlib.pyplot as plt
import glob

files = sorted(glob.glob("mpi_output_*.txt"))
i = 0
for file in files:
    data = np.loadtxt(file)
    plt.imshow(data, cmap='binary', origin='lower')
    plt.title(f"Game of Life - {file}")
    plt.pause(0.08)
    i += 1

plt.show()
