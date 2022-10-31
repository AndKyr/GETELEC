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

Temps = [1.e-2, 300, 800, 1500]

Xfn = np.linspace(0.12, 0.35, Npoints)
F = 1./Xfn

Jem = np.copy(F)

this = gt.emission_create(W = 4.5, R = 5000., approx = 2)

fig1 = plt.figure(figsize=fsize)
ax1 = fig1.gca()

ax1.set_xlabel(r"$1/F$ [m GV$^{-1}$]")
ax1.set_ylabel(r"$J$ [A nm$^{-2}$]")

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
    
    ax1.semilogy(Xfn,Jem, label = r'T = %d K'%this.Temp)
    
# for i in range(len(Temps)):
    # this.Temp = Temps[i]
    # if (this.Temp < 10.):
        # this.approx = -1
    # else:
        # this.approx = -1

    # for j in range(len(F)):
        # this.F = F[j]
        # this.cur_dens()
        # Jem[j] = this.Jem
    
    # ax1.semilogy(Xfn,Jem, '--', color = colors[i], label = r'T = %d K'%this.Temp)
    

# np.savetxt("J-F.dat", np.transpose(np.array([F,Jem])), delimiter = " ")

ax1.grid()
ax1.legend()

plt.savefig("JFplot_Tparam.svg")

plt.savefig("JFplot_Tparam.png")

plt.show()
