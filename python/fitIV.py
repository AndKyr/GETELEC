import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib



mainpath,filename = os.path.split(os.path.realpath(__file__))
emissionpath,mainfolder = os.path.split(mainpath)
pythonpath = emissionpath + '/python'
sys.path.append(pythonpath)

import getelec_mod as gt


F0 = [1., 5., 14.]
R0 = [1., 5., 50.]
gamma0 = [1., 10., 100.]
Temp0 = [299.99, 300., 300.01]

filenames  = ['IV_Guerrera.csv', 'IV_scaling.csv', \
            'IV_Spindt.csv', 'IV_CERN_2004_S22.csv']
Wi =        [4.05,                4.5,               4.35,                4.5]
markers = ['s', '+', 'd', 'o']

colors = ['b', 'g', 'r', 'k']
                
#ploting parameters
font = 45
lw = 5
mw = 15
matplotlib.rcParams["font.family"] = 'Times New Roman'
matplotlib.rcParams["font.size"] = font
matplotlib.rcParams["axes.labelsize"] = font
matplotlib.rcParams["xtick.labelsize"] = font
matplotlib.rcParams["ytick.labelsize"] = font
matplotlib.rcParams["legend.fontsize"] = font
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']


fig = plt.figure(figsize=(20,15))
ax = fig.gca()
ax.grid()
ax.set_xlabel(r"$1/F\ [ \mathrm{nm/V} ]$")
ax.set_ylabel(r"$I \ [ \mathrm{nA} ]$")


i=0                
for fname in filenames:
    Vdata, Idata =  np.loadtxt('data/' + fname,delimiter=',',unpack=True)
    xdata = 1./Vdata
    ydata = np.log(Idata)
    W0 = list(np.array([1.-1e-4, 1., 1.+1e-4]) * Wi[i])

    fit= gt.fitML(xdata,ydata, F0, W0, R0, gamma0, Temp0)
    popt = fit.x
    yopt = gt.MLplot(xdata, popt[0], popt[1], popt[2], popt[3], popt[4])
    yshift = max(yopt) - max(ydata)
    
    print 'beta = %10.3e, W = %10.3f, R = %10.3f, gamma = %10.3f, Temp = %10.3f, sigmaAeff = %10.3e' \
            % (popt[0], popt[1],  popt[2], popt[3], popt[4], 1e-9*np.exp(-yshift))
                
    xth = np.linspace(min(xdata),max(xdata),100)
    yth = np.exp(gt.MLplot(xth, popt[0], popt[1], popt[2], popt[3], popt[4]) - yshift)           
    ax.semilogy(xdata / popt[0],Idata,markers[i], \
                label = r'Data set (' + str(i+1) + r')', markersize = mw, \
                mec = colors[i], mfc = 'none', mew = 2)
    ax.semilogy(xth / popt[0],yth,c=colors[i], linewidth = lw, label = None)
    i = i + 1

ax.set_ylim([1.e-3, 2.e8])
ax.legend(loc="best")
fig.tight_layout()
#fig.savefig('/home/kyritsak/Documents/LaTeX/papers/getelec_paper/eps/figure_fittings.eps')
#fig.savefig('/home/kyritsak/Documents/LaTeX/papers/getelec_paper/png/figure_fittings.png')
plt.show()

