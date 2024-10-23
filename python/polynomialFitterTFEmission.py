import numpy as np 
import matplotlib.pyplot as plt
from itertools import cycle
import getelec as gt
import scipy

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


def polynomialCoefficientArray(x:float, y:float, Nx:int, Ny:int) -> np.ndarray:
    powersOfX = x**np.arange(Nx)
    powersOfY = y**np.arange(Ny)
    return np.outer(powersOfX, powersOfY).flatten()


#defining getelec parameters
gt._setTabulationPath("./tabulated/1D_1024")
emitter = gt.ConductionBandEmitter(workFunction=4.5)

#getting grid of temperature-field variables
Nfields = 32
Ntemperatures = 32
minField = 3.
maxField = 16.

minTemperature = 1000.
maxTemperature = 3000.

vectorOfFields = 1./np.linspace(1./maxField, 1./minField, Nfields)
# vectorOfTemperatures = 1/np.linspace(1./maxTemperature, 1/minTemperature, Ntemperatures)
vectorOfTemperatures = np.logspace(np.log10(minTemperature), np.log10(maxTemperature), Ntemperatures)


#calculate data for all grid

NpolyTemperature = 6
NpolyField = 5
Amatrix = []
bVector = []
currentDensities = np.ones((Ntemperatures, Nfields))
for i in range(Ntemperatures):
    for j in range(Nfields):
        emitter.barrier.setParameters(field=vectorOfFields[j])
        emitter.setParameters(kT = gt.Globals.BoltzmannConstant * vectorOfTemperatures[i])
        currentDensities[i,j] = emitter.currentDensity()
        #add row to the Amatrix and element to the b vector (Ax=b)
        Amatrix.append(polynomialCoefficientArray(np.log(vectorOfTemperatures[i]), 1/vectorOfFields[j], NpolyTemperature, NpolyField))
        bVector.append(np.log(currentDensities[i,j]))

np.save("TFTabulation/currentDensitiesForTfrom_%d_to_%d"%(minTemperature, maxTemperature), currentDensities)
np.save("TFTabulation/fieldsForTfrom_%d_to_%d"%(minTemperature, maxTemperature), vectorOfFields)
np.save("TFTabulation/tempearaturesForTfrom_%d_to_%d"%(minTemperature, maxTemperature), vectorOfTemperatures)


#solve the least squares problem and get polynomial coefficients
leastSquareResults =  scipy.linalg.lstsq(np.array(Amatrix), np.array(bVector))
polynomialCoefficientMatrix = np.reshape(leastSquareResults[0], (NpolyTemperature, NpolyField))

np.save("TFTabulation/polynomialCoefficientsForTfrom_%d_to_%d"%(minTemperature, maxTemperature), polynomialCoefficientMatrix)
print("least square residual = ", leastSquareResults[1])

#test validity by plotting polynomial vs calculated values
currentDensitiesPolynomial = np.copy(currentDensities)

for i in range(Ntemperatures):
    for j in range(Nfields):
        currentDensitiesPolynomial[i,j] = np.exp(np.polynomial.polynomial.polyval2d(np.log(vectorOfTemperatures[i]), 1/vectorOfFields[j], polynomialCoefficientMatrix))

plt.figure()
plt.loglog(currentDensities.flatten(), currentDensitiesPolynomial.flatten(), ".")
plt.loglog(np.array([np.min(currentDensities), np.max(currentDensities)]), np.array([np.min(currentDensities), np.max(currentDensities)]))
print("\n", NpolyTemperature, NpolyField)
print("mean of relative error", 1-np.mean(currentDensitiesPolynomial / currentDensities))

print("std of ratio: ", np.std(currentDensitiesPolynomial / currentDensities))
plt.savefig("poly2Dvalidation.png")
plt.show()
plt.close()
mappable = plt.cm.ScalarMappable(norm ="log", cmap="jet")
mappable.set_clim(min(vectorOfTemperatures), max(vectorOfTemperatures))


plt.figure()


for i in range(len(vectorOfTemperatures)):   
    plt.semilogy(1./vectorOfFields, currentDensities[i] , ".", c = mappable.to_rgba(vectorOfTemperatures[i]))
    plt.semilogy(1./vectorOfFields, currentDensitiesPolynomial[i], c = mappable.to_rgba(vectorOfTemperatures[i]))


cbar = plt.colorbar(mappable=mappable, ax = plt.gca())
cbar.ax.set_ylabel("Temperature[K]")
plt.xlabel("1/F [nm/V]")
plt.ylabel(r"J [A/nm$^2$]")
# plt.savefig("TFJ.png")
plt.show()
plt.close()
