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

cMaphot = mb.colors.ListedColormap(['white', 'red'])
cMapcold = mb.colors.ListedColormap(['white', 'blue'])
#cMap = plt.cm.gist_heat
cMap = plt.cm.copper
Npoints = 200
Ncontours = 20
Jmin = 1.e-26
Fticks = np.array([1.,2.,3.,4.,5.,6.,8.,10.])
Tticks = np.array([300.,500.,1000.,2000.,4000,6000, 10000])
good_lim = 1.2
#Jlevels = np.array([1.e-6, 1.e-8, 1.e-10, 1.e-12, 1.e-14, 1.e-16])

t = np.logspace(np.log10(300.), np.log10(10000.), Npoints)
f = np.logspace(np.log10(1.), np.log10(11.), Npoints)

#t = np.linspace(300, 5000., 50)
#f = np.linspace(0.5, 12., 50)

T,F = np.meshgrid(t,f)
Jem = np.copy(T)
Jrld = np.copy(Jem)
Jfn = np.copy(Jem)
good_rld = np.copy(Jrld)
good_fn = np.copy(Jrld)

this = gt.emission_create(R = 500.)

for i in range(len(t)):
    for j in range(len(f)):
        this.Temp = T[i,j]
        this.F = F[i,j]
        
        this.approx = 1
        this.cur_dens()
        Jem[i,j] = this.Jem
        
        this.approx = -1
        this.cur_dens()
        Jfn[i,j] = this.Jem
        
        this.approx = -2
        this.cur_dens()
        Jrld[i,j] = this.Jem
        
        #if (Jem[i,j] / Jrld[i,j] < good_lim): good_rld[i,j] = 1.
        #if (Jem[i,j] / Jfn[i,j] < good_lim): good_fn[i,j] = 1.
        
    
good_rld = Jem / Jrld
good_fn = Jem / Jfn    
    
Jem[Jem < Jmin] = Jmin
Jrld[Jrld < Jmin] = Jmin
Jfn[Jfn < Jmin] = Jmin   

     

fig1 = plt.figure(figsize=fsize)
ax1 = fig1.gca()

ax1.set_xlabel(r"$F[GV/m]$")
ax1.set_ylabel(r"$T[K]$")
ax1.set_title(r"Full numerical")

# Plot the surface.

cs = ax1.contour(F,T,Jem,50, locator=mb.ticker.LogLocator(numticks = Ncontours), cmap = cMap)

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.xaxis.set_major_formatter(mb.ticker.FormatStrFormatter('%0.0f'))
ax1.yaxis.set_major_formatter(mb.ticker.FormatStrFormatter('%0.0f'))
ax1.xaxis.set_ticks(Fticks)
ax1.yaxis.set_ticks(Tticks)
ax1.grid()
cbar1 = fig1.colorbar(cs)
cbar1.set_label(r"$J [A/nm^2]$")
fig1.savefig('numerical.png')
ax1.contourf(F,T,good_fn,levels = np.array([-100., good_lim]), cmap = cMapcold)
ax1.contourf(F,T,good_rld,levels = np.array([-100., good_lim]), cmap = cMaphot)
ax1.annotate('Cold FE regime', xy = (4,600))
ax1.annotate('Thermionic RLD regime', xy = (1.01,8000))
ax1.annotate('Intermediate regime', xy = (1.5,1500))
fig1.savefig('numerical_regions.png')

fig2 = plt.figure(figsize=fsize)
ax2 = fig2.gca()

ax2.set_xlabel(r"$F[GV/m]$")
ax2.set_ylabel(r"$T[K]$")
ax2.set_title(r"Fowler-Nordheim")

# Plot the surface.
cs = ax2.contour(F,T,Jfn,50, vmin = 1.e-18, vmax=1.e-6, locator=mb.ticker.LogLocator(numticks = Ncontours),  cmap = cMap)

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.xaxis.set_major_formatter(mb.ticker.FormatStrFormatter('%0.0f'))
ax2.yaxis.set_major_formatter(mb.ticker.FormatStrFormatter('%0.0f'))
ax2.xaxis.set_ticks(Fticks)
ax2.yaxis.set_ticks(Tticks)
ax2.grid()
cbar2 = fig2.colorbar(cs)
cbar2.set_label(r"$J [A/nm^2]$")
fig2.savefig('F-N.png')

fig3 = plt.figure(figsize = fsize)
ax3 = fig3.gca()

ax3.set_xlabel(r"$F[GV/m]$")
ax3.set_ylabel(r"$T[K]$")
ax3.set_title(r"RLD")

# Plot the surface.
cs = ax3.contour(F,T,Jrld,50, locator=mb.ticker.LogLocator(numticks = Ncontours),  cmap = cMap)

ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.xaxis.set_major_formatter(mb.ticker.FormatStrFormatter('%0.0f'))
ax3.yaxis.set_major_formatter(mb.ticker.FormatStrFormatter('%0.0f'))
ax3.xaxis.set_ticks(Fticks)
ax3.yaxis.set_ticks(Tticks)
ax3.grid()
cbar3 = fig3.colorbar(cs)
cbar3.set_label(r"$J [A/nm^2]$")
fig3.savefig('R-L-D.png')

#plt.show()
