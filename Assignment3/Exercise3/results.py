import numpy as np

serial = np.array([0.009178, 0.009189, 0.009167, 0.009191, 0.009126, 0.009146, 0.009300, 0.009170, 0.009137, 0.009191])

print("SERIAL")
print(f'Average: {np.mean(serial):.2e}')
print(f'Std dev: {np.std(serial):.2e}')

omp_sum = np.array([0.008245, 0.008216, 0.008146, 0.008192, 0.008183, 0.008196, 0.008261, 0.008196, 0.008153, 0.008147])

print("OMP SUM")
print(f'Average: {np.mean(omp_sum):.2e}')
print(f'Std dev: {np.std(omp_sum):.2e}')

omp_critical = np.array([0.648281, 0.645557, 0.645193, 0.645970, 0.647752, 0.645404, 0.647858, 0.645948, 0.669996, 0.665490])

print("OMP CRITICAL")
print(f'Average: {np.mean(omp_critical):.2e}')
print(f'Std dev: {np.std(omp_critical):.2e}')
