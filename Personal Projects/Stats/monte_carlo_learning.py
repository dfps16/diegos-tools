import numpy as np
import matplotlib.pyplot as plt

# Initiate analysis variables
sim_count = 10000000

A = np.random.uniform(1, 5, sim_count)
B = np.random.uniform(2, 6, sim_count)

duration = A + B

plt.figure(figsize=(3.5, 2.5), dpi=300)
plt.hist(duration, density=True)
plt.axvline(9, color='r')
plt.show()
print((duration > 9).sum()/sim_count)
