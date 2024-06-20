import numpy as np 
import matplotlib.pyplot as plt
from itertools import cycle
import getelec as gt

font = 15
import matplotlib as mb
mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 1.5
mb.rcParams["text.usetex"] = True
mb.rcParams["lines.markersize"] = 10
figureSize = [16,10]
colors = cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])

gt._setTabulationPath("./tabulated/1D_1024")


Nfields = 32
Ntemperatures = 16

vectorOfFields = 1./np.linspace(1./16., 1./2., Nfields)
vectorOfTemperatures = np.logspace(2., 4.3, Ntemperatures)

emitter = gt.ConductionBandEmitter(workFunction=4.5)

mappable = plt.cm.ScalarMappable(norm ="log", cmap="jet")
mappable.set_clim(min(vectorOfTemperatures), max(vectorOfTemperatures))
vectorOfCurrentDensities = np.copy(vectorOfFields)


fig = plt.figure()
ax = fig.gca()

for i in range(len(vectorOfTemperatures)):
    # fig1 = plt.figure()

    for j in range(len(vectorOfFields)):

        emitter.barrier.setParameters(field=vectorOfFields[j])
        emitter.setParameters(kT = gt.Globals.BoltzmannConstant * vectorOfTemperatures[i])
        
        # E, TED = emitter.totalEnergySpectrumArrays(numberOfPoints=512)
        # vectorOfCurrentDensities[j] = np.trapz(TED, E)

        try:
            vectorOfCurrentDensities[j] = emitter.currentDensity()
        except:
            emitter.calculateTotalEnergySpectrum()
            vectorOfCurrentDensities[j]= emitter.currentDensityFromTED()

        # fig1.gca().semilogy(E, TED, label = "F=%g"%vectorOfFields[j])
        # fig1.gca().set_title("T = %g"%vectorOfTemperatures[i])
    # plt.legend()
    # plt.show()
        # if (j == 1):
        #     plt.semilogy(E, TED, label = "T = %g"%vectorOfTemperatures[i])
        
    # plt.show()

    # print(vectorOfCurrentDensities)
        
        
    if (np.any(vectorOfCurrentDensities < 1.e-30)):
        continue
        
    poly, cov = np.polyfit(1./vectorOfFields, np.log(vectorOfCurrentDensities), 4, cov=True)

    print(poly)
    # print(np.diag(cov)**0.5)
    
    ax.semilogy(1./vectorOfFields, vectorOfCurrentDensities , ".", c = mappable.to_rgba(vectorOfTemperatures[i]))
    ax.semilogy(1./vectorOfFields, np.exp(np.polyval(poly, 1./vectorOfFields)), c = mappable.to_rgba(vectorOfTemperatures[i]))


cbar = plt.colorbar(mappable=mappable, ax = ax)
cbar.ax.set_ylabel("Temperature[K]")
ax.set_xlabel("1/F [nm/V]")
ax.set_ylabel(r"J [A/nm$^2$]")
plt.show()