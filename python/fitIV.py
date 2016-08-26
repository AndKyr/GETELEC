import numpy as np
import sys
import os
import matplotlib.pyplot as plt

mainpath,filename = os.path.split(os.path.realpath(__file__))
emissionpath,mainfolder = os.path.split(mainpath)
pythonpath = emissionpath + '/python'
sys.path.append(pythonpath)

import getelec_mod as gt

xdata, ydata =  np.loadtxt('data/IV_He_330.csv',delimiter=',',unpack=True)

F0 = [0.005, 0.05, 0.15]
W0 = [1.5, 1.7, 1.9]
R0 = [999.9, 1000., 1000.1]
gamma0 = [1., 1.01, 1.02]
Temp0 = [630., 630.01, 630.02]


fit= gt.fitML(xdata,ydata, F0, W0, R0, gamma0, Temp0)

popt = fit.x
print fit
print 'Floc = ', popt[0] / xdata

yopt = gt.MLplot(xdata, popt[0], popt[1], popt[2], popt[3], popt[4])
yshift = max(yopt) - max(ydata)
print 'sigma Aeff = ', np.exp(-yshift)

plt.plot(xdata, ydata,'r*')

plt.plot(xdata,yopt - yshift,'b-')
plt.show()


