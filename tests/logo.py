#! /usr/bin/python3

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

x1,y1 = np.loadtxt("SNbar.dat", unpack=True)
# np.savetxt("SNbar.dat", np.c_[x1, y1])
# plt.plot(x1,np.sqrt(y1),"*", color = '#b95000', linewidth=2)



# print (x1, np.sqrt(y1))
x2, y2 = np.loadtxt("KXbar.dat", unpack=True)

plt.fill(x2,np.sqrt(y2),'#004f20',linewidth=2)
plt.fill(x1,np.sqrt(y1),'#b95000',linewidth=2)
plt.savefig("log.")
plt.show()

