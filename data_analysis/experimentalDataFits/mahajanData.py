import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pathlib
import matplotlib


xFN, yFN = np.loadtxt("mahajanBadData.dat", unpack=True, delimiter=",")

voltageData = 1./xFN
currentData = np.exp(yFN) * voltageData**2

print ("Field:" , sorted(voltageData))

print ("Current Density:", sorted(currentData))

plt.plot(voltageData, currentData)

plt.grid()
plt.savefig("mahajan.png")
