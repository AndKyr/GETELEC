#!/usr/bin/python

"""This tests plots the behaviour of the emission current density for high fields
    according to three different approximations. The FN approximation (Miller-Good)
    Version of the FN equation, the GTF approximation (Jensen) and the full
    calculation by GETELEC"""
import numpy as np
import getelec_mod as gt

import matplotlib.pyplot as plt
import matplotlib as mb

font = 30
mb.rcParams["font.family"] = "Serif"
mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 2.

Npoints = 128

F = np.logspace(np.log10(5.), np.log10(15.), Npoints)
Jfn = np.copy(F)
Jgtf = np.copy(F)
Jfull = np.copy(F)

fig1 = plt.figure()
ax1 = fig1.gca()


for work in [3.5,4.,4.5]:
    this = gt.emission_create(W = work, R = 500., Temp = 1000.)

    for i in range(len(F)):
        this.F = F[i]
        this.approx = 1
        this.cur_dens()
        Jfull[i] = this.Jem
        
        this.approx = -1
        this.cur_dens()
        Jfn[i] = this.Jem
        
        this.approx = 0
        this.cur_dens()
        Jgtf[i] = this.Jem




    ax1.set_xlabel(r"$F[GV/m]$")
    ax1.set_ylabel(r"$J[K]$")
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    ax1.plot(F,Jfn,label = r"FN")
    ax1.plot(F,Jgtf,label = r"GTF")
    ax1.plot(F,Jfull,label = r"full")
    ax1.legend()

plt.show()



