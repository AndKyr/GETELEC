#! /usr/bin/python
import numpy as np
import getelec_mod as gt

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

cMaphot = mb.colors.ListedColormap(['white', 'red'])
cMapcold = mb.colors.ListedColormap(['white', 'blue'])
cMap = plt.cm.plasma
# cMap = plt.cm.copper
Npoints = 128
Jmin = 1.e-26
Fticks = np.array([1.,2.,3.,5.,8.])

T = 2000.

R = [4., 8., 1000.]

F = np.logspace(np.log10(1.), np.log10(9.), 64)

Jem = np.copy(F)
# Jrld = np.copy(Jem)
# Jfn = np.copy(Jem)
Jgtf = np.copy(Jem)

this = gt.emission_create(gamma = 30)

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
this.Temp = T
fig = plt.figure(figsize=fsize)
ax = fig.gca()

for i in range(len(R)):
    this.R = R[i]
    if (R[i] > 500):
        lbl = r'SN barrier (R = $\infty$)'
    else:
        lbl = r'$R = $%g nm'%R[i]
    for j in range(len(F)):

        this.F = F[j]
        
        this.approx = 2
        this.mode = 0
        this.cur_dens()
        Jem[j] = this.Jem
        
        this.approx = 0
        this.mode = -1
        this.cur_dens()
        Jgtf[j] = this.Jem
        
    ax.loglog(F, Jem, c = colors[i], label = lbl)#label = 'T = %.0f K, Full Numerical'%T[i])
    ax.loglog(F, Jgtf, '--', c = colors[i])#, label = 'General GTF, R=%g nm'%R[i])#label = 'T = %.0f K, Full Numerical'%T[i])


ax.set_xlabel(r"$F[GV/m]$")
ax.set_ylabel(r"$J [A/nm^2]$")
ax.xaxis.set_major_formatter(mb.ticker.FormatStrFormatter('%g'))
# ax1.yaxis.set_major_formatter(mb.ticker.FormatStrFormatter('%0.0f'))
ax.xaxis.set_ticks(Fticks)
# ax.set_ylim([1.e-13, 1.e-5])
# ax.set_xlim([1., 9.])
# ax1.yaxis.set_ticks(Tticks)
ax.grid()
ax.legend()


fig.savefig('J-F_regimes_comparison.svg')
fig.savefig('J-F_regimes_comparison.png')

plt.show()
