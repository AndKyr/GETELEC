#! /usr/bin/python
import numpy as np
import getelec_mod as gt

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mb

mb.rcParams["font.family"] = "Serif"
mb.rcParams["font.size"] = 38
mb.rcParams["axes.labelsize"] = 38
mb.rcParams["xtick.labelsize"] = 38
mb.rcParams["ytick.labelsize"] = 38
mb.rcParams["legend.fontsize"] = 38

t = np.logspace(np.log10(500.), np.log10(5000.), 50)
f = np.logspace(np.log10(1.), np.log10(12.), 50)

#t = np.linspace(300, 5000., 50)
#f = np.linspace(0.5, 12., 50)

T,F = np.meshgrid(t,f)
Jem = np.copy(T)

this = gt.emission_create()

for i in range(len(t)):
    for j in range(len(f)):
        this.approx = 0
        this.Temp = T[i,j]
        this.F = F[i,j]
        this.cur_dens()
        Jem[i,j] = this.Jem
        #print this.F, this.Temp, this.Jem

fig = plt.figure()
ax = fig.gca()

ax.set_xlabel(r"$F[GV/m]$")
ax.set_ylabel(r"$T[K]$")
ax.set_title(r"G error")



# Plot the surface.
cs = ax.contourf(F,T,Jem,50, locator=mb.ticker.LogLocator(numticks = 50))

ax.set_xscale('log')
ax.set_yscale('log')


## Add a color bar which maps values to colors.
fig.colorbar(cs, )

plt.show()
