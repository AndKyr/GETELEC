#! /usr/bin/python3
import numpy as np
import getelec_mod as gt

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mb

font = 20
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

this = gt.emission_create(R = 5000., gamma = 30)
sfile = "output/spectra.csv"

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


if (plotG):
    ax2 = ax.twinx()
    ax2.set_yscale("log")


fig = plt.figure(figsize=fsize)
fig.tight_layout()

fig2 = plt.figure(figsize=fsize)
fig2.tight_layout()

ax = fig.gca()
ax2 = fig2.gca()


F = [5., 5.98, 8.92, 10.58]

R = [5000., 4., 1., 0.7]

style = ['-', '--', '.', ':']

for i in range(len(F)):
    this.Temp = T
    this.F = F[i] 
    this.R = R[i]   
    this.approx = 2
    this.cur_dens()
    
    E,JE,J, Nj, G = np.loadtxt(sfile,unpack=True)
    
    if this.R > 100:
        Rstr = r'$\infty$'
    else:
        Rstr = '%g nm'%R[i]
        
    ax.plot(E, 100 * J / np.max(J) , style[i], \
            label = r'R = %s, F = %g V/nm, J = %.2g A/nm$^2$'%(Rstr, F[i], this.Jem))
        
    ax2.plot(E, J , style[i], \
            label = r'R = %s, F = %g V/nm, J = %.2g A/nm$^2$'%(Rstr, F[i], this.Jem))
    
    if (plotG):
        ax = plt.gca()
        ax2 = ax.twinx()
        ax2.plot(E, 1./(1. + np.exp(G)), c = colors[i], label = 'R = %.0f nm, F = %g V/nm'%(R[i], F[i]))



ax.set_xlabel(r"$E-E_F$ [eV]")
ax.set_ylabel(r"j [% of max]")
ax.set_xlim([-1.2,0.4])
ax.legend()
ax.grid()

ax2.set_xlabel(r"$E-E_F$ [eV]")
ax2.set_ylabel(r"j [A nm$^{-2}$ eV$^{-1}$]")
ax2.set_xlim([-1.25,0.5])
ax2.legend()
ax2.grid()

# fig.savefig("spectra_%.3d.png"%i)
# plt.close()

plt.show()
