#! /usr/bin/python
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
fsize = (18,10)

Npoints = 256

Works = [3.5, 4., 4.5, 5.]
T = 300.



Fmin = [2., 2.5, 3., 3.5]
Fmax = [5., 6., 7., 8.]

xdata = np.linspace(0.15, 0.35, 128)
Jem = np.copy(xdata)

this = gt.emission_create(W = 4.5, R = 5000., approx = 2)

fig1 = plt.figure(figsize=fsize)
ax1 = fig1.gca()

fig2 = plt.figure(figsize=fsize)
ax2 = fig2.gca()

fig3 = plt.figure(figsize=fsize)
ax3 = fig3.gca()

ax1.set_xlabel(r"$1/F$ [m GV$^{-1}$]")
ax1.set_ylabel(r"$J$ [A nm$^{-2}$]")


colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

for i in range(len(Works)):
    this.W = Works[i]
    xdata = np.linspace(1./Fmax[i], 1./Fmin[i], 128)
    F = 1./xdata

    for j in range(len(F)):
        this.F = F[j]
        this.cur_dens()
        Jem[j] = this.Jem

    ydata = np.log(Jem)
    ax1.semilogy(xdata,Jem, label = r'W = %g eV'%this.W)

    fit = gt.fitML(1./F, np.log(Jem), W0 = [4.5-1.e-8, 4.5, 4.5+1.e-8], R0 = [5000, 5001, 5002], gamma0=[10-1.e-8, 10., 10+1.e-8], Temp0=[T-1.e-8, T, T+1.e-8])
    popt = fit.x

  
    yopt = gt.MLplot(xdata, popt[0], popt[1], popt[2], popt[3], popt[4], approx=2)
    yshift = max(yopt) - max(ydata)
    # print(popt, np.exp(yshift))



    xth = xdata/ popt[0]
    Jth = np.exp(gt.MLplot(xdata, popt[0], popt[1], popt[2], popt[3], popt[4]) - yshift)   
    rmse = np.mean(((Jth - Jem) / Jem)**2)**0.5
    ax2.semilogy(xth, Jem, label = r'W = %g eV, $\beta$=%.4g, S = %.4g, REMSE= %.3f %s'%(this.W, popt[0]**-1, np.exp(yshift), 100 * rmse, "%"))

    der1 = np.gradient(np.log(Jem)) / np.gradient(xth)
    der2 = np.gradient(np.log(Jth)) / np.gradient(xth)
    rmse = np.mean(((der1 - der2) / der2)**2)**0.5
    ax3.plot(xth[1:-1], der1[1:-1], label = "W = %.2g, RMSE = %.1f %s"%(this.W, 100 * rmse, '%') )

    # if (i == 2):
    #     ax3.plot(xth[1:-1], der2[1:-1])

ax3.set_xlabel(r"$(\beta \cdot F)^{-1}$ [m GV$^{-1}$]")
ax3.set_ylabel(r"$\partial \logJ / \partial(F^-1)$ [GV/m]")
ax3.grid()


ax1.legend()
ax2.legend()
ax3.legend()

ax2.set_xlabel(r"$(\beta \cdot F)^{-1}$ [m GV$^{-1}$]")
ax2.set_ylabel(r"$S \cdot J$ [A nm$^{-2}$]")
ax1.grid()
ax2.grid()

fig1.savefig("Works_plot_compare.png")
fig2.savefig("Works_scaled.png")
fig3.savefig("derivative_scaled.png")

plt.show()
