#!/usr/bin/python

"""This tests plots the behaviour of the emission current density for high fields
    according to three different approximations. The FN approximation (Miller-Good)
    Version of the FN equation, the GTF approximation (Jensen) and the full
    calculation by GETELEC"""
import numpy as np
import getelec_mod as gt

import matplotlib.pyplot as plt
import matplotlib as mb

font = 25
# mb.rcParams["font.family"] = "Serif"
mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 2.

Npoints = 256

F = np.linspace(1,12)
Jfn = np.copy(F)
Jgtf = np.copy(F)
Jfull = np.copy(F)
Jrld= np.copy(F)

fig = plt.figure()
ax = fig.gca()
ax.grid()

T = [600,3200]
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

for j in range(len(T)):
    this = gt.emission_create(W = 4.5, R = 5000., Temp = T[j])

    for i in range(len(F)):
        this.F = F[i]
        this.approx = 2
        this.cur_dens()
        Jfull[i] = this.Jem
        
        this.approx = -1
        this.cur_dens()
        Jfn[i] = this.Jem
        
        this.approx = 0
        this.cur_dens()
        Jgtf[i] = this.Jem
        
        this.approx = -2
        this.cur_dens()
        Jrld[i] = this.Jem




    ax.set_xlabel(r"$F[GV/m]$")
    ax.set_ylabel(r"$J[A / nm^2]$")
    ax.set_ylim([1.e-10, 1.e-4])
    # ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid()

    ax.plot(F,Jfn, label = r"FN, T = %dK"%T[j], c=colors[j], linestyle = "-." )
    ax.plot(F,Jgtf,label = r"GTF, T = %dK "%T[j], c=colors[j], linestyle = "--")
    ax.plot(F,Jrld,label = r"RLD, T = %dK"%T[j], c=colors[j], linestyle = ":")
    ax.plot(F,Jfull,label = r"Num, T = %dK"%T[j], c=colors[j])
    ax.legend(loc = 'best')
    

plt.show()



