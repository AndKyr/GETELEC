#!/usr/bin/python

"""This plots the spectra as outputed in spectra.csv from getelec.""" 
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mb
import os

font = 30
mb.rcParams["font.family"] = "Serif"
mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 2.




fig1 = plt.figure()
ax1 = fig1.gca()
ax1.set_xlabel(r"$E [eV]$")
ax1.set_ylabel(r"$ j(E) | [A nm{-2} (eV)^{-1}]$")

fig2 = plt.figure()
ax2 = fig2.gca()
ax2.set_xlabel(r"$E [eV]$")
ax2.set_ylabel(r"$ G(E) $")


fields = [6., 8., 10.]

colors = ['b', 'r', 'g', 'k', 'm']



for i in range(len(fields)):
    field = fields[i]
    appstr = os.popen('./bin/current.exe %f 4.5 300. | grep "Gamow\|@Ef" '%field).read().split()
    Gam = float(appstr[2])
    dGam = float(appstr[-2])
    E,JE,J,fj,G = np.loadtxt("spectra.csv", delimiter = ',' ,unpack=True)
    maxJ = max(fj)
    Eplot = E[fj > maxJ * 0.01]
    Jplot = fj[fj > maxJ * 0.01]
    Gplot = G[fj > maxJ * 0.01]
    Gapprox = -Eplot * dGam + Gam
    Japprox = np.log(1. + np.exp(-Eplot/0.0258)) / (1. + np.exp(Gapprox))

    ax1.plot(Eplot, Jplot, colors[i]+'-', label = "F=%3.2f, full"%field)
    ax1.plot(Eplot, Japprox,  colors[i]+'--', label = "F=%3.2f, approx"%field)
    ax2.plot(Eplot, Gplot, colors[i]+'-', label = "F=%3.2f, full"%field)
    ax2.plot(Eplot, Gapprox, colors[i]+'--', label = "F=%3.2f, approx"%field)
    print "Field = %f, Japprox = %f, Jfull = %f"%(field, np.trapz(Japprox,Eplot), np.trapz(Jplot,Eplot))
    
ax1.set_yscale('linear')    
ax1.legend(loc = "best")

ax2.set_yscale('linear')    
ax2.legend(loc = "best")
plt.show()
