#! /usr/bin/python3

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

data  = {'Guerrera': (np.array([2.413e+02, 2.511e+02, 2.622e+02, 2.706e+02, 2.803e+02, 2.915e+02, 2.999e+02, 3.096e+02, 3.208e+02, 3.305e+02, 3.403e+02, 3.515e+02, 3.612e+02, 3.710e+02, 3.808e+02, 3.891e+02, 4.003e+02, 4.100e+02, 4.184e+02, 4.338e+02, 4.491e+02, 4.644e+02, 4.826e+02, 4.993e+02]), \
                    np.array([8.719e-01, 1.582e+00, 2.967e+00, 5.038e+00, 8.555e+00, 1.406e+01, 2.309e+01, 3.670e+01, 5.643e+01, 8.678e+01, 1.249e+02, 1.797e+02, 2.502e+02, 3.371e+02, 4.540e+02, 6.116e+02, 7.459e+02, 9.720e+02, 1.267e+03, 1.764e+03, 2.376e+03, 3.096e+03, 4.310e+03, 5.617e+03])) \
        }
data['ETH_scaling'] = (np.array([6.443e+00, 6.327e+00, 6.188e+00, 6.029e+00, 5.895e+00, 5.774e+00, 5.650e+00, 5.496e+00, 5.384e+00, 5.270e+00, 5.154e+00, 5.043e+00, 4.903e+00, 4.775e+00, 4.685e+00, 4.604e+00, 4.505e+00, 4.421e+00, 4.361e+00, 4.269e+00, 4.180e+00, 4.127e+00, 4.064e+00, 4.006e+00, 3.924e+00, 3.852e+00]), \
                      np.array([2.998e+02, 2.359e+02, 1.697e+02, 1.164e+02, 8.714e+01, 6.288e+01, 4.623e+01, 3.095e+01, 2.189e+01, 1.515e+01, 1.068e+01, 7.385e+00, 4.553e+00, 2.876e+00, 1.999e+00, 1.394e+00, 9.257e-01, 6.180e-01, 4.712e-01, 3.005e-01, 1.958e-01, 1.407e-01, 1.047e-01, 7.498e-02, 4.599e-02, 3.137e-02])) 

data['Spindt'] = (np.array([56.923, 60.055, 64.011, 67.967, 73.077, 79.011, 86.923, 95.989, 107.030, 119.070, 122.030, 126.980, 130.930, 134.070]), \
                       np.array([1.000e+03, 2.970e+03, 1.000e+04, 2.970e+04, 1.000e+05, 2.970e+05, 1.000e+06, 2.970e+06, 1.000e+07, 2.970e+07, 3.981e+07, 6.051e+07, 8.111e+07, 1.000e+08]))
data['CERN_2004'] =  (np.array([1.238e-01, 1.192e-01, 1.117e-01, 1.058e-01, 1.011e-01, 9.618e-02, 9.282e-02, 8.918e-02, 8.490e-02, 7.980e-02, 7.862e-02, 7.563e-02, 7.320e-02, 7.092e-02, 6.878e-02, 6.620e-02, 6.459e-02, 6.356e-02]), \
                       np.array([1.814e+03, 9.700e+02, 5.192e+02, 2.278e+02, 1.018e+02, 5.039e+01, 2.564e+01, 1.293e+01, 5.429e+00, 1.688e+00, 9.456e-01, 5.051e-01, 1.964e-01, 9.531e-02, 3.333e-02, 1.431e-02, 6.666e-03, 3.526e-03]))
Wi = {'Guerrera': 4.05, 'ETH_scaling': 4.5, 'Spindt': 4.35, 'CERN_2004': 4.5}
markers = {'Guerrera': 's', 'ETH_scaling': '+', 'Spindt': 'd', 'CERN_2004': 'o'}
#markers = ['s', '+', 'd', 'o']
colors = {'Guerrera': 'b', 'ETH_scaling': 'g', 'Spindt': 'r', 'CERN_2004': 'k'}
                
#ploting parameters
font = 45
lw = 5
mw = 18
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

                
for key in data.keys():
    Vdata = data[key][0]
    Idata =  data[key][1]
    xdata = 1./Vdata
    ydata = np.log(Idata)
    #print Vdata
    #print Idata
    W0 = list(np.array([1.-1e-4, 1., 1.+1e-4]) * Wi[key])

    fit= gt.fitML(xdata,ydata, F0, W0, R0, gamma0, Temp0)
    popt = fit.x
    yopt = gt.MLplot(xdata, popt[0], popt[1], popt[2], popt[3], popt[4])
    yshift = max(yopt) - max(ydata)
    
    print('beta = %10.3e, W = %10.3f, R = %10.3f, gamma = %10.3f, Temp = %10.3f, sigmaAeff = %10.3e' \
            % (popt[0], popt[1],  popt[2], popt[3], popt[4], 1e-9*np.exp(-yshift)))
                
    xth = np.linspace(min(xdata),max(xdata),100)
    yth = np.exp(gt.MLplot(xth, popt[0], popt[1], popt[2], popt[3], popt[4]) - yshift)           
    ax.semilogy(xdata / popt[0],Idata,markers[key], \
                label = key, markersize = mw, \
                mec = colors[key], mfc = 'none', mew = 2)
    ax.semilogy(xth / popt[0],yth,c=colors[key], linewidth = lw, label = None)

ax.set_ylim([1.e-3, 2.e8])
ax.set_yticks([1.e8, 1.e6, 1.e4, 1.e2, 1.e0, 1.e-2])
ax.legend(loc="best")
fig.tight_layout()
#fig.savefig('/home/kyritsak/Documents/LaTeX/papers/getelec_paper/eps/figure_fittings.eps')
#fig.savefig('/home/kyritsak/Documents/LaTeX/papers/getelec_paper/png/figure_fittings.png')
plt.show()

