import gamowCalculator as gc
import matplotlib.pyplot as plt
import numpy as np
import os

import getelec_wrap as gt
font = 25
import matplotlib as mb
mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 1.5
mb.rcParams["text.usetex"] = True
figureSize = [10,10]


emitter = gt.GetelecInterface()
emitter.set_kT(0.5)

emitter.run(["CurrentDensity", "TotalEnergyDistribution"], writeFlag=True)


files = [s for s in os.listdir(".") if "odeSolution_i_0_E_" in s]
energies = [float(s[18:-4]) for s in files]

energies = np.sort(np.array(energies))

cmap = mb.colormaps["jet"]

norm = mb.colors.Normalize(vmin=min(energies), vmax=max(energies))

fig, axs = plt.subplots(3, 1, sharex=True, figsize=figureSize)
ax1, ax2, ax3 = axs
for i in range(1, len(energies), 3):
    filename = "odeSolution_i_0_E_%f.dat"%energies[i]
    y1 = np.loadtxt(filename)
    color = cmap(norm(energies[i]))


    ax1.plot(y1[:,0], y1[:,1], color = color)
    ax2.plot(y1[:,0], y1[:,2], color = color)
    ax3.plot(y1[:,0], y1[:,3], color = color)

# ax1.grid()
# ax2.grid()
# ax3.grid()

fig.colorbar(plt.cm.ScalarMappable(cmap = "jet", norm=norm), ax = axs, label = r"$E \textrm{ [eV]}$")


ax3.set_xlabel(r"$z \textrm{ [nm]}$")
ax1.set_ylabel(r"$\Re(s')$")
ax2.set_ylabel(r"$\Im(s')$")
ax3.set_ylabel(r"$\Re(s)$")

plt.savefig("solutionFamily.pdf")


filename = ["splineSolution_i_0.dat", "splineNodes_i_0.dat"]

y1 = np.loadtxt(filename[0])
nodes = np.loadtxt(filename[1])

fig, ax1 = plt.subplots(1, 1, figsize=figureSize)

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


ax1.plot(y1[:,0], y1[:,1], label=r"$\Re(s')$", color = colors[0])
ax1.plot(nodes[:,0], nodes[:,1], 'o',  color = colors[0])

ax1.plot(y1[:,0], y1[:,2], label=r"$\Im(s')$", color = colors[1])
ax1.plot(nodes[:,0], nodes[:,3], 'o', color = colors[1])

ax1.plot(y1[:,0], y1[:,3], label=r"$\Re(s)$", color = colors[2])
ax1.plot(nodes[:,0], nodes[:,5], 'o', label = r"$\textrm{Spline nodes}$", color = colors[2])


ax1.legend()
ax1.set_xlabel(r"$E-E_F \textrm{ [eV]}$")
ax1.set_ylabel(r"$\mathbf{s}(z_0)$")
ax1.grid()
plt.savefig("splineSolution.pdf")
plt.show()

