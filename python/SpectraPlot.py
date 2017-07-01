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


files = os.listdir("./")
sfiles = [myfile for myfile in files if "spectra" in myfile]




for sfile in sfiles:
    E,JE,J, G = np.loadtxt(sfile, delimiter = ',' ,unpack=True)
    ax1.plot(E,J,label = sfile)
    ax2.plot(E,G,'b-', label = sfile)
    
ax1.set_yscale('log')    
ax1.legend(loc = "best")

ax2.set_yscale('linear')    
ax2.legend(loc = "best")
plt.show()



