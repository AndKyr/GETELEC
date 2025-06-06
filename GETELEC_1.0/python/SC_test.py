#!/usr/bin/python

"""This tests plots the behaviour of the emission current density for high fields
    according to three different calculations. The FN equation (Murphy-Good)
    the Child Langmuir-Limit, and the full self-consistent calculation GETELEC"""
import numpy as np
import getelec_mod as getelec_old

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mb

font = 30
# mb.rcParams["font.family"] = "Serif"
mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 2.5
fsize = (18,10)


Npoints = 128
# Fticks = np.array([.7, 1.,2.,3.,5.,8.])

F = np.logspace(np.log10(5.), np.log10(20.), Npoints)

Jfn = np.copy(F)
Jsc = np.copy(F)

Vappl = 10000

this = getelec_old.emission_create(R = 5000., voltage = Vappl, approx = -1)

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

fig = plt.figure(figsize=fsize)
ax = fig.gca()

for j in range(len(F)):
    this.F = F[j]
    
    this.cur_dens()
    Jfn[j] = this.Jem
    

    this.cur_dens_SC()
    Jsc[j] = this.Jem
    

    
ax.semilogy(1/F, Jfn, '--', label = 'FN')#label = 'T = %.0f K, Full Numerical'%T[i])
ax.semilogy(1/F, Jsc, label = 'Self consistent')#label = 'T = %.0f K, Full Numerical'%T[i])
ax.semilogy(1/F,getelec_old.J_ChildL(Vappl, F / Vappl), '--', label = 'CL') #linestyle = '--', label = 'T = %.0f K, F-N'%T[i])

ax.grid()
ax.set_xlabel(r"$1/F[nm/V]$")
ax.set_ylabel(r"$J[A/nm^2]$")
ax.legend()

plt.savefig('/media/sf_VM-shared/SC_comparison.png')

plt.show()



