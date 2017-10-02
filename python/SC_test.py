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

Npoints = 32

F = np.logspace(np.log10(5.), np.log10(25.), Npoints)
J_SC = np.copy(F)
F_p = np.copy(F)

fig1 = plt.figure()
ax1 = fig1.gca()
fig2 = plt.figure()
ax2 = fig2.gca()



for i in range(len(F)):
    J_SC[i], theta = gt.emit_SC(F[i], 4.5, approx = 2)
    F_p[i] = theta * F[i]




ax1.set_xlabel(r"$F[GV/m]$")
ax2.set_xlabel(r"$F[GV/m]$")
ax1.set_ylabel(r"$J[A/nm^2]$")
ax1.set_yscale('log')
ax1.plot(F,J_SC)

ax2.set_ylabel(r"$F_p[GV/m]$")

ax2.plot(F,F_p)

plt.show()



