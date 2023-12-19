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

fig, ax = plt.subplots(len(Voltages), 1, sharex=True)

energyShift = 0.7915

workFunction = 2.78
beta = 0.00064



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
    ax[i].plot(energy, spectra, ".", color=colors[i])

    ax[i].plot([0,0], [0,1], "k")
    ax[i].plot([workFunction, workFunction], [0,1], "b")
    ax[i].plot([topOfBarrier, topOfBarrier], [0,1], "r")

    ax[i].plot(spectrum["energy"][0], spectrum["electronCount"][0] / np.max(spectrum["electronCount"][0]), "-", color=colors[i])

plt.savefig("spectraZrO.png")

plt.figure()
plt.semilogy(beta * Voltages,  currentDensitiesExp / currentDensitiesExp[0], label="Experimental")
plt.semilogy(beta * Voltages, currentDensitiesCalc / currentDensitiesCalc[0], label="calculation")
plt.legend()
plt.savefig("iVplot.png")

plt.show()


