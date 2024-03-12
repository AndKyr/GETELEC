import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pathlib
import matplotlib
import scipy.optimize as opt


filePath = pathlib.Path(os.path.realpath(__file__))
srcPath = filePath.parent.parent.parent.joinpath("src")
sys.path.insert(0, str(srcPath))

import getelec as gt
# import IVDataFitter



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

Voltages = np.array([1500, 2000, 2300, 2400, 2500])

fig, ax = plt.subplots(len(Voltages), 1, sharex=True, figsize=[10,12])

energyShift = 0.7915

workFunction = 2.83
beta = 0.00064

# workFunction = 3.
# beta = 0.0009



gtModel = gt.GETELECModel(workFunction=workFunction, temperature=1800., field=1.5)
currentDensitiesExp = np.array([19., 69., 148., 186., 240.])
currentDensitiesCalc = np.copy(currentDensitiesExp)

for i in range(len(Voltages)):
    energy, spectra = np.loadtxt("kimSpectra1800_V%d.dat"%Voltages[i], unpack=True, delimiter=",")
    spectra -= i * 7
    spectra /= np.max(spectra)

    gtModel.setParameters(field=Voltages[i] * beta)
    gtModel.run(calculateCurrent=True, calculateSpectrum=True)
    currentDensitiesCalc[i] = gtModel.getCurrentDensity()[0]
    print(gtModel.field, currentDensitiesCalc[i])
    spectrum = gtModel.getElectronSpectrum()

    topOfBarrier = workFunction - 3.79 * np.sqrt(0.1 * beta * Voltages[i])
    
    energy -= energyShift
    ax[i].plot(energy, spectra, ".", markersize=5, color=colors[i],\
               label=r"Exp. $V=%.2g \textrm{ kV}, I =%d \textrm{ }\mu \textrm{A/sr}$"%(1.e-3*Voltages[i], currentDensitiesExp[i]))

    # ax[i].plot([0,0], [0,1], "k")
    ax[i].plot([workFunction, workFunction], [0,1], "k:")
    ax[i].plot([topOfBarrier, topOfBarrier], [0,1], "k--")

    ax[i].plot(spectrum["energy"][0], spectrum["electronCount"][0] / np.max(spectrum["electronCount"][0]), "-", color=colors[i], \
                label=r"Theory $J =%.2f \textrm{ nA/nm}^2$"%(1.e9*currentDensitiesCalc[i]))
    ax[i].legend(loc="upper left")
    ax[i].set_xlim([-1.,3.])
    ax[i].grid()

ax[-1].set_xlabel(r"$E-E_F$ [eV]")
ax[-1].text(workFunction+0.05, 0.3, r"$\phi$")
ax[-1].text(topOfBarrier+0.05, 0.3, r"$U_m$")
# ax[-1].text(0.05, 0.3, r"$E_F$")

plt.savefig("spectraZrO.pdf")

plt.figure()
plt.semilogy(beta * Voltages,  currentDensitiesExp / currentDensitiesExp[-1], label="Experimental")
plt.semilogy(beta * Voltages, currentDensitiesCalc / currentDensitiesCalc[-1], label="calculation")
plt.xlabel(r"$F= \beta \cdot V$ [GV/m]")
plt.ylabel(r"Current density [normalized]")
plt.legend()
plt.grid()
plt.savefig("iVplot.pdf")

plt.show()


