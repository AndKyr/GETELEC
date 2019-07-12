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
Ncontours = 20
Jmin = 1.e-26
Fticks = np.array([.2, .3, .5, .7, 1.,2.,3.,5.,8.])
Tticks = np.array([300.,500.,1000.,2000.,4000,6000, 10000])
good_lim = 1.2
#Jlevels = np.array([1.e-6, 1.e-8, 1.e-10, 1.e-12, 1.e-14, 1.e-16])

T = [2500]

F = np.logspace(np.log10(0.2), np.log10(9.), 64)

Jem = np.copy(F)
Jrld = np.copy(Jem)
Jfn = np.copy(Jem)

this = gt.emission_create(R = 5000.)

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

fig = plt.figure(figsize=fsize)
ax = fig.gca()

for i in range(len(T)):
    for j in range(len(F)):
        this.Temp = T[i]
        this.F = F[j]
        
        this.approx = 2
        this.cur_dens()
        Jem[j] = this.Jem
        
        this.approx = -1
        this.cur_dens()
        Jfn[j] = this.Jem
        
        this.approx = -2
        this.cur_dens()
        Jrld[j] = this.Jem
    
    ax.loglog(F, Jem, c = colors[i], label = 'T = %.0f K, Full Numerical'%T[i])
    ax.loglog(F,Jfn, c = colors[i], linestyle = '--', label = 'T = %.0f K, F-N'%T[i])
    ax.loglog(F,Jrld,c = colors[i], linestyle = ':', label = 'T = %.0f K, R-L-D'%T[i])


ax.set_xlabel(r"$F[GV/m]$")
ax.set_ylabel(r"$J [A/nm^2]$")
ax.xaxis.set_major_formatter(mb.ticker.FormatStrFormatter('%0.0g'))
# ax1.yaxis.set_major_formatter(mb.ticker.FormatStrFormatter('%0.0f'))
ax.xaxis.set_ticks(Fticks)
ax.set_ylim([1.e-13, 1.e-5])
ax.set_xlim([.2, 9.])
# ax1.yaxis.set_ticks(Tticks)
ax.grid()
ax.legend()


fig.savefig('J-F_regimes_comparison.svg')
fig.savefig('J-F_regimes_comparison.png')
# ax1.contour(F,T,good_fn,levels = np.array([-100., good_lim]), cmap = cMapcold, linestyle = "--")
# ax1.contour(F,T,good_rld,levels = np.array([-100., good_lim]), cmap = cMaphot, linestyle = "--")
# ax1.annotate('Cold FE regime', xy = (4,600))
# ax1.annotate('Thermionic regime', xy = (1.01,8000))
# ax1.annotate('Intermediate regime', xy = (1.5,1500))
# fig1.savefig('numerical_regions.svg')
# fig1.savefig('numerical_regions.png')

# fig2 = plt.figure(figsize=fsize)
# ax2 = fig2.gca()

# ax2.set_xlabel(r"$F[GV/m]$")
# ax2.set_ylabel(r"$T[K]$")
# ax2.set_title(r"Fowler-Nordheim")

# # Plot the surface.
# cs = ax2.contourf(F,T,Jfn,50, vmin = 1.e-18, vmax=1.e-6, locator=mb.ticker.LogLocator(numticks = Ncontours),  cmap = cMap)

# ax2.set_xscale('log')
# ax2.set_yscale('log')
# ax2.xaxis.set_major_formatter(mb.ticker.FormatStrFormatter('%0.0f'))
# ax2.yaxis.set_major_formatter(mb.ticker.FormatStrFormatter('%0.0f'))
# ax2.xaxis.set_ticks(Fticks)
# ax2.yaxis.set_ticks(Tticks)
# ax2.grid()
# cbar2 = fig2.colorbar(cs)
# cbar2.set_label(r"$J [A/nm^2]$")
# fig2.savefig('F-N.svg')
# fig2.savefig('F-N.png')

# fig3 = plt.figure(figsize = fsize)
# ax3 = fig3.gca()

# ax3.set_xlabel(r"$F[GV/m]$")
# ax3.set_ylabel(r"$T[K]$")
# ax3.set_title(r"RLD")

# # Plot the surface.
# cs = ax3.contourf(F,T,Jrld,50, locator=mb.ticker.LogLocator(numticks = Ncontours),  cmap = cMap)

# ax3.set_xscale('log')
# ax3.set_yscale('log')
# ax3.xaxis.set_major_formatter(mb.ticker.FormatStrFormatter('%0.0f'))
# ax3.yaxis.set_major_formatter(mb.ticker.FormatStrFormatter('%0.0f'))
# ax3.xaxis.set_ticks(Fticks)
# ax3.yaxis.set_ticks(Tticks)
# ax3.grid()
# cbar3 = fig3.colorbar(cs)
# cbar3.set_label(r"$J [A/nm^2]$")
# fig3.savefig('R-L-D.svg')
# fig3.savefig('R-L-D.png')

plt.show()
