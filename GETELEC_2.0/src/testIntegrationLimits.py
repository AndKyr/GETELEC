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


Nfields = 8
Ntemperatures = 8

vectorOfFields = 1./np.linspace(1./16., 1./2., Nfields)
vectorOfTemperatures = np.logspace(2.5, 4.3, Ntemperatures)



vectorOfCurrentDensities = np.copy(vectorOfFields)
emitter = gt.ConductionBandEmitter(workFunction=4.5)


# itests = [2, 2, 3, 3, 4, 5]
# jtests = [6, 7, 2, 3, 1, 0]

# for i in range(len(itests)):
#     emitter.barrier.setParameters(field=vectorOfFields[jtests[i]])
#     emitter.setParameters(kT = gt.Globals.BoltzmannConstant * vectorOfTemperatures[itests[i]])
#     emitter.calculateTotalEnergySpectrum()
#     E, TED = emitter.totalEnergySpectrumArrays(numberOfPoints=512)

#     plt.semilogy(E, TED, label = "F=%g, T=%g, %s"%(emitter.barrier._field, emitter.kT / gt.Globals.BoltzmannConstant, emitter.regime))
#     plt.legend()
#     # plt.savefig("integrationLimitTest.png")
#     plt.show()




for i in range(len(vectorOfTemperatures)):

    for j in range(len(vectorOfFields)):

        emitter.barrier.setParameters(field=vectorOfFields[j])
        emitter.setParameters(kT = gt.Globals.BoltzmannConstant * vectorOfTemperatures[i])
        emitter.calculateTotalEnergySpectrum()
        
        E, TED = emitter.totalEnergySpectrumArrays(numberOfPoints=512)
        vectorOfCurrentDensities[j] = np.trapz(TED, E)

        plt.semilogy(E, TED, label = "n=%g, %s"%(emitter.regimeDeterminingParameter,emitter.regime))
        plt.title("i = %g"%i)
    plt.legend()
    plt.show()
        