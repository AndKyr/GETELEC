import matplotlib.pyplot as plt
import numpy as np

import getelec_wrap as gt
font = 45
import matplotlib as mb
mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 2.0
mb.rcParams["text.usetex"] = True
figureSize = [15,10]
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

from gamowCalculator import GamowCalculator




gamowCalculator = GamowCalculator(solverType="WKB", XCdataFile="", minimumPotential=10.)
gamowCalculator.findBarrierMax()
gamowWKB = gamowCalculator.calculateGamowWKB(4.5)
print(gamowWKB)



waveVectors = np.linspace(1.e-5, 40., 256)
energies = np.linspace(-4.5, -4.5 + 1.e-5, 256)
kineticEnergies = waveVectors**2 * gt.Globals.hbarSqrOver2m

emitter = gt.GetelecInterface()
coefficients = emitter.calculateTransmissionCoefficients(energies, waveVectors, writeFiles=False)
fig, ax1 = plt.subplots(figsize = figureSize, tight_layout=True)
ax2 = ax1.twinx()

ax1.plot(waveVectors, np.abs(coefficients)**2 / waveVectors, label=r"$D$", color = colors[0])
ax1.plot([min(waveVectors), max(waveVectors)], [np.exp(-gamowWKB), np.exp(-gamowWKB)], label = r"$D - \textrm{JWKB}$", color = colors[3])
ax3 = ax2.twiny()

ax3.plot( kineticEnergies, np.angle(coefficients) / np.pi, label=r"$\arg(d)$", color = colors[1])



ax1.set_xlabel(r"$k \textrm{[nm}^{-1}\textrm{]}$")
ax1.set_ylabel(r"$D$")
ax2.set_ylabel(r"$\arg(d) / \pi$")
ax3.set_xlabel(r"$\frac{\hbar^2 k_z^2}{2m} \textrm{ [eV]}$")
# ax1.set_ylim([0., 0.05])
ax3.set_ylim([-0.5, 0])

ax1.yaxis.label.set_color(colors[0])
ax2.yaxis.label.set_color(colors[1])
ax1.tick_params(axis='y', colors=colors[0])
ax3.tick_params(axis='y', colors=colors[1])
ax2.tick_params(axis='y', colors=colors[1])

ax1.spines["left"].set_color(colors[0])
ax3.spines["right"].set_color(colors[1])
ax2.spines["right"].set_color(colors[1])



# ax1.grid()
# ax2.grid()
ax1.legend(loc="lower right")
ax3.legend(loc="lower left")

plt.savefig("wavevectorDependence.svg")
plt.show()

