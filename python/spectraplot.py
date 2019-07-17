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
fsize = (16,10)

cMaphot = mb.colors.ListedColormap(['white', 'red'])
cMapcold = mb.colors.ListedColormap(['white', 'blue'])
cMap = plt.cm.plasma
# cMap = plt.cm.copper

plotG = False

T = 300.

F = np.logspace(np.log10(2.), np.log10(8.),2.)

this = gt.emission_create(R = 5000.)
sfile = "output/spectra.csv"

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']



if (plotG):
    ax2 = ax.twinx()
    ax2.set_yscale("log")


for i in range(len(F)):
    this.Temp = T
    this.F = F[i]    
    this.approx = 2
    this.cur_dens()
    
    E,JE,J, Nj, G = np.loadtxt(sfile, delimiter = ',' ,unpack=True)
    
     
    fig = plt.figure(figsize=fsize)
    fig.tight_layout()
    ax = fig.gca()
    ax.plot(E, 100* J / max(J), label = 'T = %.0f K, F = %.0f V/nm'%(T, F[i]))
    
    if (plotG):
        ax2.plot(E, 1./(1. + np.exp(G)), c = colors[i], label = 'T = %.0f K, F = %.0f V/nm'%(T[i], F[i]))

    ax.set_xlabel(r"$E-E_F$ [eV]")
    ax.set_ylabel(r"j [% of max]")
    ax.set_xlim([-2.,0.5])
    ax.grid()

    fig.savefig("spectra_%.3d.png"%i)
    plt.close()

# plt.show()
