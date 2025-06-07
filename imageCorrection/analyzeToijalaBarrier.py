import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pathlib
import matplotlib

font = 15
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams["font.family"] = "a"
matplotlib.rcParams["font.size"] = font
matplotlib.rcParams["axes.labelsize"] = font
matplotlib.rcParams["xtick.labelsize"] = font
matplotlib.rcParams["ytick.labelsize"] = font
matplotlib.rcParams["legend.fontsize"] = font
# matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['figure.figsize'] = 10, 6
matplotlib.rcParams['lines.markersize'] = 10
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

x, potential = np.loadtxt("potentialToijala.dat", unpack=True, delimiter=",")

xPBE, potentialPBE = np.loadtxt("potentialPBE.dat", unpack=True, delimiter=",")
sortInds = np.argsort(xPBE)

x0 = 16.25
Field = .2
Q = 3.6
WorkFunction = 4.76
SNbarrier = -Field *(x - x0) + WorkFunction - Q/(x-x0)

def barrierIntegrate (x, V):
    return np.trapz(np.sqrt(V[V>0]), x[V>0])

# goodinds = np.where(np.logical_and(x > 17.5, x < 21.8))
goodinds = np.where(x > x0 + .4)


# plt.semilogy(x[goodinds], correction[goodinds] )

plt.plot(x - x0, potential, label="mix")
plt.plot(x[goodinds] - x0, SNbarrier[goodinds], label="SN barrier")
plt.plot(xPBE[sortInds] - x0, potentialPBE[sortInds], label="PBE")
plt.plot(x[x > x0] - x0, WorkFunction - Field*(x[x > x0] - x0), label="Triangular" )
# plt.plot(x[potential > 0], np.sqrt(potential[potential > 0]), label="mix")
# plt.plot(x[np.logical_and(SNbarrier > 0, x > x0)], np.sqrt(SNbarrier[np.logical_and(SNbarrier > 0, x > x0)]), label="SN barrier")
# plt.plot(xPBE[sortInds] - x0, potentialPBE[sortInds], label="PBE corrected")


print("Dmix = ", np.exp(-barrierIntegrate(x, potential)), "D_SN = ", np.exp(-barrierIntegrate(x[x > x0],SNbarrier[x > x0])), "DPBE = ", np.exp(-barrierIntegrate(xPBE[sortInds], potentialPBE[sortInds])))

plt.plot
plt.grid()
plt.legend()
plt.xlabel(r"$x \textrm{[} \AA \textrm{]}$")
plt.ylabel(r"$\textrm{Barrier [eV]}$")
plt.savefig("toijalaPotential.png")
plt.show()