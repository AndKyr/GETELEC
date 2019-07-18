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

Npoints = 256

# Temps = [1.e-2, 300, 1000, 2000]
Temps = [1.e-2]

Xfn = np.linspace(0.12, 0.35, 256)
F = 1./Xfn

Jem = np.copy(F)

this = gt.emission_create(W = 4.5, R = 5000., approx = 2)

fig1 = plt.figure(figsize=fsize)
ax1 = fig1.gca()

ax1.set_xlabel(r"$1/F$ [m GV$^{-1}$]")
ax1.set_ylabel(r"$J$ [A nm$^{-2}$]")

ax2 = ax1.twinx()
ax2.set_ylabel(r'$\log_{10}(J/F^2)$', color = 'red')
ax2.tick_params(axis='y', labelcolor='red')
# ax2.set_yscale('linear')

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


for i in range(len(Temps)):
    this.Temp = Temps[i]
    if (this.Temp < 10.):
        this.approx = -1
    else:
        this.approx = 2

    for j in range(len(F)):
        this.F = F[j]
        this.cur_dens()
        Jem[j] = this.Jem
    
    ax1.semilogy(Xfn,Jem, 'k-', label = r'M-L plot')
    ax1.semilogy([Xfn[0], Xfn[-1]], [Jem[0], Jem[-1]], 'k--', label = 'Straight line ML')
    ax2.plot(Xfn, np.log10(Jem / F**2),'r-',  label = r'F-N plot')
    ax2.plot([Xfn[0], Xfn[-1]], [ np.log10(Jem[0] / F[0]**2), np.log10(Jem[-1] / F[-1]**2)],'r--', label = 'Straight line FN')
    

# np.savetxt("J-F.dat", np.transpose(np.array([F,Jem])), delimiter = " ")

ax2.set_ylim([-15, -7.5])
ax1.grid()
ax1.legend(loc = 3)
ax2.legend(loc = 1)

plt.savefig("ML-FN_plot_compare.png")

plt.show()
