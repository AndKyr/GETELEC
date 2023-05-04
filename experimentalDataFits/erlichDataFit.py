import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pathlib
import matplotlib

filePath = pathlib.Path(os.path.realpath(__file__))
srcPath = filePath.parent.parent.joinpath("src")
print(str(srcPath))
sys.path.insert(0, str(srcPath))
import IVDataFitter



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

ivFitter = IVDataFitter.IVDataFitter()

xFN, yFN = np.loadtxt("erlichData.dat", unpack=True, delimiter=",")

voltageData = 1.e4/xFN
currentData = np.exp(yFN) * voltageData**2 * 1.e2



ivFitter.setIVdata(voltageData=voltageData, currentData=currentData)
ivFitter.setParameterRange(workFunction=5.25)
ivFitter.fitIVCurve()
fittedCurrent = ivFitter.getOptCurrentCurve()

ivFitter.printFittingData()
plt.semilogy(ivFitter.fittingParameters["fieldConversionFactor"] * voltageData, fittedCurrent, label = r"$\beta=%.2g \textrm{ nm}^{-1}, \sigma=%.2g, \langle \delta J \rangle = %.2g$"%(ivFitter.fittingParameters["fieldConversionFactor"], ivFitter.preFactor, ivFitter.getFittingError()))

plt.semilogy(ivFitter.fittingParameters["fieldConversionFactor"] * voltageData, fittedCurrent / ivFitter.preFactor, label = r"$\beta=%.2g \textrm{ nm}^{-1}, \sigma=%.2g, \langle \delta J \rangle = %.2g$"%(ivFitter.fittingParameters["fieldConversionFactor"], 1, ivFitter.getFittingError()))


plt.errorbar(ivFitter.fittingParameters["fieldConversionFactor"] * voltageData, currentData, yerr=currentData * .3, fmt="." , label="Experimental Erlich 1978")
plt.xlabel(r"$\textrm{local field [GV/m]} $")
plt.ylabel(r"$\textrm{current density [A/nm}^2\textrm{]} $")
plt.legend()
plt.grid()
plt.savefig("fittedIVCurveErlich.png")
