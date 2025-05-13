"""
@file getelec_wrap.py
@brief Python wrapper for the Getelec C++ library using ctypes.
This module provides a Python interface to the Getelec C++ library, allowing for the calculation of various electronic properties.
Classes:
    GetelecInterface: A class to interface with the Getelec C++ library.
Usage example:
    getelec.setRadius([5., 6., 7., 8., 9.])
    getelec.run()
    densities = getelec.getCurrentDensity()
"""

import ctypes
import numpy as np
import os
import scipy.interpolate as spi

"""
@class GetelecInterface
@brief A class to interface with the Getelec C++ library.
This class provides methods to set parameters, run simulations, and retrieve results from the Getelec C++ library.
@param libPath Path to the shared library file.
@param configPath Path to the configuration file.
@param barrierType Type of barrier to use in calculations.
Methods:
    __init__(self, libPath=None, configPath:str="", barrierType="modifiedSN"): Initializes the interface and loads the shared library.
    calculateTransmissionForEnergy(self, energy, params_index): Calculates the transmission coefficient for a given energy.
    calculateTransmissionForEnergies(self, energies, params_index): Calculates the transmission coefficients for multiple energies.
    calculateTransmissionForManyEnergies(self, energies, params_index): Calculates the transmission coefficients for many energies.
    _toCtypesArray(self, numpy_array): Converts a numpy array to a ctypes array.
    setField(self, fields): Sets the field values.
    setRadius(self, radii): Sets the radius values.
    setGamma(self, gammas): Sets the gamma values.
    set_kT(self, kTs): Sets the kT values.
    setWorkFunction(self, workFunctions): Sets the work function values.
    setBandDepth(self, bandDepths): Sets the band depth values.
    setEffectiveMass(self, effectiveMasses): Sets the effective mass values.
    run(self, calculate_spectra=False): Runs the simulation.
    getCurrentDensity(self): Retrieves the current density values.
    getNottinghamHeat(self): Retrieves the Nottingham heat values.
    getSpectra(self): Retrieves the spectra values.
    calculateBarrierValues(self, energies, params_index=0): Calculates the barrier values for given energies.
    calculateCurrentDensity(self): Runs the simulation and retrieves the current density values.
    calculateNottinghamHeat(self): Runs the simulation and retrieves the Nottingham heat values.
    getBarrierIntegrationLimits(self, params_index=0): Retrieves the barrier integration limits.
    getBarrierPlotData(self, Npoints=512, params_index=0): Retrieves the barrier plot data.
    setRandomParams(self, numberOfInputs): Sets random parameters for the simulation.
    setParameters(self, **params): Sets multiple parameters for the simulation.
    __del__(self): Cleans up the object and deletes the Getelec instance.
"""

class GetelecInterface:


    flagMap = {
        "CurrentDensity": 1,
        "NottinghamHeat": 2,
        "TotalEnergyDistribution": 4,
        "NormalEnergyDistribution": 8,
        "ParallelEnergyDistribution": 16,
        "TotalEnergyDistributionDerivatives": 32,
    }

    def __init__(self, libPath=None, configPath:str="", barrierType="modifiedSN", fields = 5., radii = 1.e5, gammas = 10., kTs = 0.025, workFunctions = 4.5, bandDepths = 10., effectiveMasses = 1., numberOfThreads = 0, seed = 1987):
        if libPath is None:
            libPath = os.path.join(os.path.dirname(__file__), "../build/libgetelec.so")
        
        # Load the shared library
        self.lib = ctypes.CDLL(libPath)

        # Initialize the Getelec instance
        self.obj = self.lib.Getelec_new_with_config(ctypes.c_char_p(configPath.encode('utf-8')), ctypes.c_char_p(barrierType.encode('utf-8')), ctypes.c_int(numberOfThreads), ctypes.c_void_p(), ctypes.c_int(seed))

        # Define argument and return types for methods
        self.lib.Getelec_setField.restype = None
        self.lib.Getelec_setField.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]

        self.lib.Getelec_getField.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_getField.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t)]

        self.lib.Getelec_setRadius.restype = None
        self.lib.Getelec_setRadius.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]

        self.lib.Getelec_getRadius.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_getRadius.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t)]

        self.lib.Getelec_setGamma.restype = None
        self.lib.Getelec_setGamma.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double)]

        self.lib.Getelec_getGamma.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_getGamma.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t)]

        self.lib.Getelec_setkT.restype = None
        self.lib.Getelec_setkT.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]

        self.lib.Getelec_getkT.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_getkT.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t)]

        self.lib.Getelec_setWorkFunction.restype = None
        self.lib.Getelec_setWorkFunction.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]

        self.lib.Getelec_getWorkFunction.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_getWorkFunction.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t)]

        self.lib.Getelec_setBandDepth.restype = None
        self.lib.Getelec_setBandDepth.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]

        self.lib.Getelec_getBandDepth.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_getBandDepth.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t)]

        self.lib.Getelec_setEffectiveMass.restype = None
        self.lib.Getelec_setEffectiveMass.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]

        self.lib.Getelec_getEffectiveMass.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_getEffectiveMass.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t)]

        self.lib.Getelec_setRandomParameters.restype = None
        self.lib.Getelec_setRandomParameters.argtypes = [ctypes.c_void_p, ctypes.c_size_t]

        self.lib.Getelec_run.restype = ctypes.c_size_t
        self.lib.Getelec_run.argtypes = [ctypes.c_void_p, ctypes.c_uint32]

        self.lib.Getelec_getCurrentDensities.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_getCurrentDensities.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t)]

        self.lib.Getelec_getNottinghamHeats.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_getNottinghamHeats.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t)]

        self.lib.Getelec_getSpectraEnergies.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_getSpectraEnergies.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t), ctypes.c_char]

        self.lib.Getelec_getSpectraValues.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_getSpectraValues.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t), ctypes.c_char]

        self.lib.Getelec_getSpectraDerivatives.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_getSpectraDerivatives.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t)]

        self.lib.Getelec_calculateTransmissionProbability.restype = ctypes.c_double
        self.lib.Getelec_calculateTransmissionProbability.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_size_t]

        self.lib.Getelec_calculateTransmissionProbabilities.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_calculateTransmissionProbabilities.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_size_t]

        self.lib.Getelec_interpolateTransmissionProbabilities.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_interpolateTransmissionProbabilities.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_size_t]

        self.lib.Getelec_getCalculationStatus.restype = ctypes.c_uint
        self.lib.Getelec_getCalculationStatus.argtypes = [ctypes.c_void_p]
    
        self.lib.Getelec_delete.restype = None
        self.lib.Getelec_delete.argtypes = [ctypes.c_void_p]

        self.lib.Getelec_getBarrierIntegrationLimits.restype = None
        self.lib.Getelec_getBarrierIntegrationLimits.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]

        self.lib.Getelec_getBarrierValues.restype = None
        self.lib.Getelec_getBarrierValues.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_size_t]

        self.setField(fields)
        self.setRadius(radii)
        self.setGamma(gammas)
        self.set_kT(kTs)
        self.setWorkFunction(workFunctions)
        self.setBandDepth(bandDepths)
        self.setEffectiveMass(effectiveMasses)
        self.resultDictionary = {}

    def calculateTransmissionProbability(self, energy, waveVector = -1., params_index = 0):
        return self.lib.Getelec_calculateTransmissionProbability(self.obj, energy, waveVector, params_index)

    def calculateTransmissionProbabilities(self, energies, waveVectors, params_index):
        energies_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(energies, dtype=np.float64)))
        waveVectors_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(waveVectors, dtype=np.float64)))
        result_ptr = self.lib.Getelec_calculateTransmissionProbabilities(self.obj, energies_ctypes, waveVectors_ctypes, size, params_index)
        return np.ctypeslib.as_array(result_ptr, shape=(size,))

    def calculateTransmissionForManyEnergies(self, energies, waveVectors, params_index=0):
        energies_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(energies, dtype=np.float64)))
        waveVectors_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(waveVectors, dtype=np.float64)))
        result_ptr = self.lib.Getelec_interpolateTransmissionProbabilities(self.obj, energies_ctypes, waveVectors_ctypes, size, params_index)
        return np.ctypeslib.as_array(result_ptr, shape=(size,))

    def _toCtypesArray(self, numpy_array):
        return numpy_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(numpy_array)

    def setField(self, fields):
        self.fields = fields
        fields_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(fields, dtype=np.float64)))
        self.lib.Getelec_setField(self.obj, fields_ctypes, size)

    def setRadius(self, radii):
        self.radii = radii
        radii_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(radii, dtype=np.float64)))
        self.lib.Getelec_setRadius(self.obj, radii_ctypes, size)
    
    def setGamma(self, gammas):
        self.gammas = gammas
        gammas_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(gammas, dtype=np.float64)))
        self.lib.Getelec_setGamma(self.obj, gammas_ctypes, size)
    
    def set_kT(self, kTs):
        self.kTs = kTs
        kTs_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(kTs, dtype=np.float64)))
        self.lib.Getelec_setkT(self.obj, kTs_ctypes, size)
    
    def setWorkFunction(self, workFunctions):
        self.workFunctions = workFunctions
        workFunctions_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(workFunctions, dtype=np.float64)))
        self.lib.Getelec_setWorkFunction(self.obj, workFunctions_ctypes, size)
    
    def setBandDepth(self, bandDepths):
        self.bandDepths = bandDepths
        bandDepths_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(bandDepths, dtype=np.float64)))
        self.lib.Getelec_setBandDepth(self.obj, bandDepths_ctypes, size)
    
    def setEffectiveMass(self, effectiveMasses):
        self.effectiveMasses = effectiveMasses
        effectiveMasses_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(effectiveMasses, dtype=np.float64)))
        self.lib.Getelec_setEffectiveMass(self.obj, effectiveMasses_ctypes, size)
    
    def getInputsDict(self):
        return {"fields": self.fields, "radii": self.radii, "gammas": self.gammas, "kTs": self.kTs, "workFunctions": self.workFunctions, "bandDepths": self.bandDepths, "effectiveMasses": self.effectiveMasses}

    def _syncInputParameters(self):
        size = ctypes.c_size_t()
        fields_ptr = self.lib.Getelec_getField(self.obj, ctypes.byref(size))
        self.fields = np.ctypeslib.as_array(fields_ptr, shape=(size.value,))

        radii_ptr = self.lib.Getelec_getRadius(self.obj, ctypes.byref(size))
        self.radii = np.ctypeslib.as_array(radii_ptr, shape=(size.value,))

        gammas_ptr = self.lib.Getelec_getGamma(self.obj, ctypes.byref(size))
        self.gammas = np.ctypeslib.as_array(gammas_ptr, shape=(size.value,))

        kTs_ptr = self.lib.Getelec_getkT(self.obj, ctypes.byref(size))
        self.kTs = np.ctypeslib.as_array(kTs_ptr, shape=(size.value,))

        workFunctions_ptr = self.lib.Getelec_getWorkFunction(self.obj, ctypes.byref(size))
        self.workFunctions = np.ctypeslib.as_array(workFunctions_ptr, shape=(size.value,))

        bandDepths_ptr = self.lib.Getelec_getBandDepth(self.obj, ctypes.byref(size))
        self.bandDepths = np.ctypeslib.as_array(bandDepths_ptr, shape=(size.value,))

        effectiveMass_ptr = self.lib.Getelec_getEffectiveMass(self.obj, ctypes.byref(size))
        self.effectiveMasses = np.ctypeslib.as_array(effectiveMass_ptr, shape=(size.value,))

    def getCalculationStatus(self):
        combinedFlagsInt = self.lib.Getelec_getCalculationStatus(self.obj)

        flags = []
        for flagName, flagValue in self.flagMap.items():
            if combinedFlagsInt & flagValue:
                flags.append(flagName)

        return flags

    def run(self, calculationFlags):
        combinedFlagsInt:int = 0
        for flagStr in calculationFlags:
            if flagStr in self.flagMap:
                combinedFlagsInt |= self.flagMap[flagStr]
            else:
                raise ValueError(f"Invalid flag: {flagStr}")

        self.numberOfIterations = self.lib.Getelec_run(self.obj, combinedFlagsInt)

    def getCurrentDensity(self):
        if ("CurrentDensity" not in self.getCalculationStatus()):
            return 0
        size = ctypes.c_size_t()
        densities_ptr = self.lib.Getelec_getCurrentDensities(self.obj, ctypes.byref(size))
        result = np.ctypeslib.as_array(densities_ptr, shape=(size.value,))
        self.resultDictionary["CurrentDensity"] = result
        return result
        
    def getNottinghamHeat(self):
        if ("NottinghamHeat" not in self.getCalculationStatus()):
            return 0
        size = ctypes.c_size_t()
        densities_ptr = self.lib.Getelec_getNottinghamHeats(self.obj, ctypes.byref(size))
        result = np.ctypeslib.as_array(densities_ptr, shape=(size.value,)) 
        self.resultDictionary["NottinghamHeat"] = result
        return result
    
    def extractAllSpectra(self):
        self.totalEnergyDistributions = []
        self.normalEnergyDistributions = []
        self.parallelEnergyDistributions = []
        self.totalEnergyDistributionSplines = []
        calculationStatus = self.getCalculationStatus()

        if ("TotalEnergyDistribution" in calculationStatus):
            for i in range(self.numberOfIterations):
                (energies, values) = self._extractSpectra(i, 'T')
                self.totalEnergyDistributions.append((energies, values))
                if ("TotalEnergyDistributionDerivatives" in calculationStatus):
                    energies2, derivatives = self._extractSpectra(i, 'D')
                    assert np.allclose(energies, energies2), "abssicae for TED values and derivatives don't match"
                    self.totalEnergyDistributionSplines.append(spi.CubicHermiteSpline(energies, values, derivatives))
            self.resultDictionary["TotalEnergyDistributions"] = self.totalEnergyDistributions
            if ("TotalEnergyDistributionDerivatives" in calculationStatus):
                self.resultDictionary["TotalEnergyDistributionSplines"] = self.totalEnergyDistributionSplines


        if ("NormalEnergyDistribution" in calculationStatus):
            for i in range(self.numberOfIterations):
                self.normalEnergyDistributions.append(self._extractSpectra(i, "N"))
            self.resultDictionary["NormalEnergyDistributions"] = self.normalEnergyDistributions

        if ("ParallelEnergyDistribution" in calculationStatus):
            for i in range(self.numberOfIterations):
                self.parallelEnergyDistributions.append(self._extractSpectra(i, "P"))
            self.resultDictionary["ParallelEnergyDistributions"] = self.parallelEnergyDistributions

    def getResultDictionary(self):
        self.getCurrentDensity()
        self.getNottinghamHeat()
        self.extractAllSpectra()
        return self.resultDictionary
    
    def getAllDataDictionary(self):
        inputDict = self.getInputsDict()
        resultDict = self.getResultDictionary()
        inputDict.update(resultDict)
        return inputDict

        
    def _extractSpectra(self, i:int, spectraType:str) -> tuple[np.ndarray, np.ndarray]:
        size = ctypes.c_size_t()
        energies_ptr = self.lib.Getelec_getSpectraEnergies(self.obj, ctypes.c_size_t(i), ctypes.byref(size), ctypes.c_char(spectraType.encode("utf-8")))
        values_ptr = self.lib.Getelec_getSpectraValues(self.obj, ctypes.c_size_t(i) , ctypes.byref(size), ctypes.c_char(spectraType.encode("utf-8")))
        energies = np.ctypeslib.as_array(energies_ptr, shape=(size.value,))
        values = np.ctypeslib.as_array(values_ptr, shape=(size.value,))
        sortInds = np.argsort(energies)
        energies = energies[sortInds]
        values = values[sortInds]
        return (energies, values)
                    
    def calculateBarrierValues(self, energies, params_index = 0):
        energyArray = np.atleast_1d(np.asarray(energies, dtype=np.float64))
        potentialArray = np.copy(energyArray)
        energies_ctypes, size = self._toCtypesArray(energyArray)
        potential_ctypes, size = self._toCtypesArray(potentialArray)

        self.lib.Getelec_getBarrierValues(self.obj, energies_ctypes, potential_ctypes, size, params_index)
        return np.ctypeslib.as_array(potential_ctypes, shape=(size,))
    
    def calculateCurrentDensity(self):
        self.run(["CurrentDensity"])
        return self.getCurrentDensity()
    
    def calculateNottinghamHeat(self):
        self.run(["CurrentDensity", "NottinghamHeat"])
        return self.getNottinghamHeat()
    
    def getBarrierIntegrationLimits(self, params_index = 0):
        lowerLimit = ctypes.c_double()
        upperLimit = ctypes.c_double()
        self.lib.Getelec_getBarrierIntegrationLimits(self.obj, ctypes.byref(lowerLimit), ctypes.byref(upperLimit), ctypes.c_size_t(params_index))
        return lowerLimit.value, upperLimit.value
    
    def getBarrierPlotData(self, Npoints = 512, params_index = 0):
        lowerLimit, upperLimit = self.getBarrierIntegrationLimits(params_index=params_index)
        energies = np.linspace(lowerLimit, upperLimit, Npoints)
        values = self.calculateBarrierValues(energies, params_index=params_index)
        return energies, values
        
    
    def setRandomParams(self, numberOfInputs):
        self.lib.Getelec_setRandomParameters(self.obj, ctypes.c_size_t(numberOfInputs))
        self._syncInputParameters()

    def setParameters(self, **params) -> None:
        for key, value in params.items():
            if key == "field":
                self.setField(value)
            elif key == "radius":
                self.setRadius(value)
            elif key == "gamma":
                self.setGamma(value)
            elif key == "kT":
                self.set_kT(value)
            elif key == "temperature":
                self.set_kT(value * 8.617333262145e-5)
            elif key == "workFunction":
                self.setWorkFunction(value)
            elif key == "bandDepth":
                self.setBandDepth(value)
            elif key == "effectiveMass":
                self.setEffectiveMass(value)

    def __del__(self):
        # Cleanup the object
        if hasattr(self, "lib") and self.lib:
            self.lib.Getelec_delete(self.obj)


if (__name__ == "__main__"):
    # Example usage:
    # Assuming the shared library is compiled as "getelec.so"
    getelec = GetelecInterface(configPath="getelec.cfg")
    getelec.setRadius([5., 6., 7. , 8., 9.])

    workFunction = 4.5
    bandDepth = 10.
    getelec.setBandDepth([bandDepth])
    getelec.setWorkFunction([workFunction])
    # getelec.setRandomInputs(5)
    # getelec.set_radius(5.)

    getelec.run(calculate_spectra=True)
    densities = getelec.getCurrentDensity()

    spectra = getelec.getSpectra()

    import matplotlib.pyplot as plt
    fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True)
    i = 0
    for spline in spectra:
        energiesPlot = np.linspace(min(spline.x), max(spline.x), 512)
        spectraPlot = spline(energiesPlot)

        barrierX, potentialEnergy = getelec.getBarrierPlotData(params_index=i)
        potentialEnergy += workFunction

        ax1.plot(barrierX, potentialEnergy)

        ax2.plot(spectraPlot, energiesPlot)        
        i += 1
    
    ax1.plot([-1, barrierX[-1]], [workFunction, workFunction], "k")
    ax1.text(barrierX[-1], workFunction, "Vacuum Level", ha="right")

    ax1.plot([-1., barrierX[-1]], [-bandDepth, -bandDepth], "k")
    ax1.text(barrierX[-1], -bandDepth, "Band bottom", ha="right")

    ax1.plot([-1., barrierX[-1]], [0., 0.], "k")
    ax1.text(barrierX[-1], 0., "Fermi level", ha="right")
    ax2.plot([min(spectraPlot), max(spectraPlot)],[0., 0.], "k")
    ax1.plot([0., 0.], [-bandDepth, workFunction], "k")
    
    ax2.set_xlabel("TED")
    ax2.set_ylabel("Energy (eV)")


    ax1.set_xlabel("Barrier position (nm)")
    ax1.set_ylabel("Potential energy (eV)")
    ax1.set_ylim(-bandDepth - .5, workFunction + .5)



    getelec.setRadius(1.e5)
    xFN = np.linspace(.1, 1., 512)
    getelec.setField(1./xFN)
    getelec.run()
    currentDensities = getelec.getCurrentDensity()
    plt.figure()
    plt.semilogy(xFN, currentDensities)
    plt.show()



