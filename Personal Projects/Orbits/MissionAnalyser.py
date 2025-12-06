import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


data = pd.read_csv("LunaAlt.txt", sep="   ")

time, altitude = data.columns

epoch = np.array(data[time])
alt = np.array(data[altitude])

plt.figure(1)
plt.plot(epoch, alt)
plt.show()
