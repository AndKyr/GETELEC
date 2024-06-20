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


Nfields = 16
Ntemperatures = 4

vectorOfFields = 1./np.linspace(1./16., 1./2., Nfields)
vectorOfTemperatures = np.logspace(2., 4.3, Ntemperatures)

emitter = gt.ConductionBandEmitter(workFunction=4.5)

mappable = plt.cm.ScalarMappable(norm ="log", cmap="jet")
mappable.set_clim(min(vectorOfTemperatures), max(vectorOfTemperatures))
vectorOfCurrentDensities = np.copy(vectorOfFields)

for i in range(len(vectorOfTemperatures)):
    emitter.setParameters(kT = gt.Globals.BoltzmannConstant * vectorOfTemperatures[i])
    
    for j in range(len(vectorOfFields)):
        emitter.barrier.setParameters(field=vectorOfFields[j])
        emitter.calculateTotalEnergySpectrum()
        
        E, TED = emitter.totalEnergySpectrumArrays()
        vectorOfCurrentDensities[j] = np.trapz(TED, E)
        if (j == 1):
            plt.plot(E, TED)
        
    # plt.show()

    print(vectorOfCurrentDensities)
        
        
    if (np.any(vectorOfCurrentDensities < 1.e-30)):
        continue
    
    poly = np.polyfit(1./vectorOfFields, np.log(vectorOfCurrentDensities), 6)
    
    # plt.semilogy(1./vectorOfFields, vectorOfCurrentDensities , ".", c = mappable.to_rgba(vectorOfTemperatures[i]))
    # plt.semilogy(1./vectorOfFields, np.exp(np.polyval(poly, 1./vectorOfFields)), c = mappable.to_rgba(vectorOfTemperatures[i]))


# cbar = plt.colorbar(mappable=mappable, ax = plt.gca())
# cbar.ax.set_ylabel("Temperature[K]")
# plt.xlabel("1/F [nm/V]")
# plt.ylabel(r"J [A/nm$^2$]")
plt.show()