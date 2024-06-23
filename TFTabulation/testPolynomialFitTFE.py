import numpy as np 
import matplotlib.pyplot as plt
import scipy


def polynomialCoefficientArray(x:float, y:float, Nx:int, Ny:int) -> np.ndarray:
    powersOfX = x**np.arange(Nx)
    powersOfY = y**np.arange(Ny)
    return np.outer(powersOfX, powersOfY).flatten()

NpolyTemperature = 6
NpolyField = 5
minTemperatures = [400, 1000, 3000]
maxTemperatures = [1000, 3000, 20000]

for iSegment in range(3):
    minTemperature = minTemperatures[iSegment]
    maxTemperature = maxTemperatures[iSegment]
    
    vectorOfFields = np.load("fieldsForTfrom_%d_to_%d.npy"%(minTemperature, maxTemperature))
    vectorOfTemperatures = np.load("tempearaturesForTfrom_%d_to_%d.npy"%(minTemperature, maxTemperature))

    #calculate data for all grid

    currentDensities = np.load("currentDensitiesForTfrom_%d_to_%d.npy"%(minTemperature, maxTemperature))

    Amatrix = []
    bVector = []
    for i in range(len(vectorOfTemperatures)):
        for j in range(len(vectorOfFields)):
            #add row to the Amatrix and element to the b vector (Ax=b)
            Amatrix.append(polynomialCoefficientArray(np.log(vectorOfTemperatures[i]), 1/vectorOfFields[j], NpolyTemperature, NpolyField))
            bVector.append(np.log(currentDensities[i,j]))


    #solve the least squares problem and get polynomial coefficients
    leastSquareResults =  scipy.linalg.lstsq(np.array(Amatrix), np.array(bVector))
    polynomialCoefficientMatrix = np.reshape(leastSquareResults[0], (NpolyTemperature, NpolyField))
    np.savetxt("polynomialCoefficientsForTfrom_%d_to_%d.dat"%(minTemperature, maxTemperature), polynomialCoefficientMatrix)

    #test validity by plotting polynomial vs calculated values
    currentDensitiesPolynomial = np.copy(currentDensities)
    for i in range(len(vectorOfTemperatures)):
        for j in range(len(vectorOfFields)):
            currentDensitiesPolynomial[i,j] = np.exp(np.polynomial.polynomial.polyval2d(np.log(vectorOfTemperatures[i]), 1/vectorOfFields[j], polynomialCoefficientMatrix))

    print("mean of relative error", 1-np.mean(currentDensitiesPolynomial / currentDensities))
    print("std of relative error: ", np.std(currentDensitiesPolynomial / currentDensities))

    #plot the results
    plt.figure()
    mappable = plt.cm.ScalarMappable(norm ="log", cmap="jet")
    mappable.set_clim(min(vectorOfTemperatures), max(vectorOfTemperatures))

    for i in range(len(vectorOfTemperatures)):   
        plt.semilogy(1./vectorOfFields, currentDensities[i] , ".", c = mappable.to_rgba(vectorOfTemperatures[i]))
        plt.semilogy(1./vectorOfFields, currentDensitiesPolynomial[i], c = mappable.to_rgba(vectorOfTemperatures[i]))

    cbar = plt.colorbar(mappable=mappable, ax = plt.gca())
    cbar.ax.set_ylabel("Temperature[K]")
    plt.xlabel("1/F [nm/V]")
    plt.ylabel(r"J [A/nm$^2$]")
    plt.savefig("TFJPlotForTfrom_%d_to_%d.png"%(minTemperature, maxTemperature))
    plt.show()
    plt.close()
