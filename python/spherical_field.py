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
fsize = (15,10)

cMap = plt.cm.plasma
Npoints = 1024
Ncontours = 20

Radius = 1.

""" plotting hte contours"""
# x = np.linspace(-3, 3, Npoints)
# X, Z = np.meshgrid(x,z)
# R = np.sqrt(X**2 + Z**2)

# phi =  Z * (1. - (Radius / R)**3)

# phi[R < Radius] = -1.e-7

# fig1 = plt.figure(figsize = fsize)
# fig1.tight_layout()
# ax1 = fig1.gca()
# ax1.set_xlabel(r"$x$ [$\mu$m]")
# ax1.set_ylabel(r"$z$ [$\mu$m]")
# ax1.set_title(r'$F_{a} = 1 GV/m$')
# cs = ax1.contour(X,Z,phi,20, cmap = cMap)
# ax1.set_xlim([-3.1, 3.1])
# ax1.set_ylim([-0.1, 6.2])
# ax1.axis('equal')
# cbar1 = fig1.colorbar(cs)
# cbar1.set_label(r"$\Phi$ [kV]")
# ax1.grid()
# plt.savefig('hemisphere_contours.svg')

""" Comparison between top and side"""
# z = np.linspace(0, 4., Npoints)
# xtop = 0.
# xside = 10.
# rtop =  np.sqrt(xtop**2 + z**2)
# rside =  np.sqrt(xside**2 + z**2)
# phitop = z * (1. - (Radius / rtop)**3)
# phiside = z * (1. - (Radius / rside)**3)
# phitop[rtop < Radius] = 0.


# fig1 = plt.figure(figsize = fsize)
# fig1.tight_layout()
# ax1 = fig1.gca()
# ax1.set_xlabel(r"$z$ [$\mu$m]")
# ax1.set_ylabel(r"$\Phi$ [kV]")

# ax1.plot(z,phitop, label = 'apex')
# ax1.plot(z,phiside, label = 'side')
# ax1.set_xlim([-0.1, 4.])
# ax1.set_ylim([-0.1, 4.])

# ax1.legend()

# ax1.grid()
# plt.savefig('hemisphere_lines.svg')


# ax1.set_xlim([0.9, 1.5])
# ax1.set_ylim([0, 1.5])

# plt.savefig('hemisphere_lines_zoom1.svg')

# ax1.set_xlim([0.9, 1.1])
# ax1.set_ylim([0, 1.1])

# plt.savefig('hemisphere_lines_zoom2.svg')


""" comparison of enhancement with different radius"""


""" Comparison between top and side"""


z = np.linspace(Radius, 3 * Radius, Npoints)
xtop = 0.
rtop =  np.sqrt(xtop**2 + z**2)
phitop = z * (1. - (Radius / rtop)**3)
phitop[rtop < Radius] = 0.


fig1 = plt.figure(figsize = fsize)
fig1.tight_layout()
ax1 = fig1.gca()
ax1.set_xlabel(r"$z$ [nm]")
ax1.set_ylabel(r"$\Phi$ [V]")

for Radius in [5, 10, 50, 100, 1000]:
    z = np.linspace(Radius, 3 * Radius, Npoints)
    xtop = 0.
    rtop =  np.sqrt(xtop**2 + z**2)
    phitop = z * (1. - (Radius / rtop)**3)
    phitop[rtop < Radius] = 0.
    phitop = z * (1. - (Radius / rtop)**3)
    phitop[rtop < Radius] = 0.
    ax1.plot(z - Radius,phitop, label = r'R = %g nm'%Radius)

ax1.set_xlim([-0.1, 4.])
ax1.set_ylim([-0.1, 10.])

ax1.legend()

ax1.grid()
# plt.savefig('hemisphere_lines.svg')


# ax1.set_xlim([0.9, 1.5])
# ax1.set_ylim([0, 1.5])

# plt.savefig('hemisphere_lines_zoom1.svg')

# ax1.set_xlim([0.9, 1.1])
# ax1.set_ylim([0, 1.1])

plt.savefig('hemisphere_lines_radii.svg')




plt.show()





