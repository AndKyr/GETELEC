#! /usr/bin/python

#fitting parameters: minimum, initial guess and maximum values
F0 = [0.1, 0.51, 2.]
R0 = [100., 101., 102.]
gamma0 = [1., 1.2, 1.3]
Temp0 = [299., 300., 300.]
W0 = [5.2, 5.21, 5.22]


#ploting parameters
markers = 'o'
colors = 'b'
font = 45
lw = 5
mw = 15

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib
import sys

mainpath,filename = os.path.split(os.path.realpath(__file__))
emissionpath,mainfolder = os.path.split(mainpath)
pythonpath = emissionpath + '/python'
sys.path.append(pythonpath)
import getelec_mod as gt

matplotlib.rcParams["font.family"] = 'Times New Roman'
matplotlib.rcParams["font.size"] = font
matplotlib.rcParams["axes.labelsize"] = font
matplotlib.rcParams["xtick.labelsize"] = font
matplotlib.rcParams["ytick.labelsize"] = font
matplotlib.rcParams["legend.fontsize"] = font
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

filename  = sys.argv[1]
fig = plt.figure(figsize=(20,15))
ax = fig.gca()
ax.grid()
ax.set_xlabel(r"$1/F\ [ \mathrm{nm/V} ]$")
ax.set_ylabel(r"$I \ [ \mathrm{nA} ]$")


Vdata, Idata =  np.loadtxt(filename,unpack=True)
xdata = 1./Vdata
ydata = np.log(Idata)

fit= gt.fitML(xdata,ydata, F0, W0, R0, gamma0, Temp0)
popt = fit.x
yopt = gt.MLplot(xdata, popt[0], popt[1], popt[2], popt[3], popt[4])
yshift = max(yopt) - max(ydata)
    
print 'beta = %10.3e, W = %10.3f, R = %10.3f, gamma = %10.3f, Temp = %10.3f, sigmaAeff = %10.3e' \
        % (popt[0], popt[1],  popt[2], popt[3], popt[4], 1e-9*np.exp(-yshift))
                
xth = np.linspace(min(xdata),max(xdata),100)
yth = np.exp(gt.MLplot(xth, popt[0], popt[1], popt[2], popt[3], popt[4]) - yshift)           
ax.semilogy(xdata / popt[0],Idata,markers, \
            label = r'Data', markersize = mw, \
            mec = colors, mfc = 'none', mew = 2)
ax.semilogy(xth / popt[0],yth,c=colors, linewidth = lw, label = 'fitting')

ax.legend(loc="best")
fig.tight_layout()
plt.show()

