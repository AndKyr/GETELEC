#! /usr/bin/python
import numpy as np
import getelec_mod as gt

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mb

font = 30
mb.rcParams["font.family"] = "Serif"
mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 2.5
fsize = (18,10)

Npoints = 2048
F = np.logspace(np.log10(1.), np.log10(15.), Npoints)

Jem = np.copy(F)

this = gt.emission_create(W = 4.5, R = 5000., approx = 2)

for j in range(len(F)):
    this.F = F[j]
    this.cur_dens()
    Jem[j] = this.Jem
    

np.savetxt("J-F.dat", np.transpose(np.array([F,Jem])), delimiter = " ")

fig1 = plt.figure(figsize=fsize)
ax1 = fig1.gca()

ax1.set_xlabel(r"$F[GV/m]$")
ax1.set_ylabel(r"$J[A/nm^2]$")

ax1.semilogy(F,Jem)
ax1.grid()

plt.show()
