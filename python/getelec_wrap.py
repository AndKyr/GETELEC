import ctypes
import numpy as np
import os
import scipy.interpolate as spi

class GetelecInterface:
    def __init__(self, libPath=None, configPath:str="", barrierType="modifiedSN"):
        if libPath is None:
            libPath = os.path.join(os.path.dirname(__file__), "../build/libgetelec.so")
        
        # Load the shared library
        self.lib = ctypes.CDLL(libPath)

        # Initialize the Getelec instance
        self.obj = self.lib.Getelec_new_with_config(ctypes.c_char_p(configPath.encode('utf-8')), ctypes.c_char_p(barrierType.encode('utf-8')))

        # Define argument and return types for methods
        self.lib.Getelec_setField.restype = None
        self.lib.Getelec_setField.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]

        self.lib.Getelec_setRadius.restype = None
        self.lib.Getelec_setRadius.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]

        self.lib.Getelec_setGamma.restype = None
        self.lib.Getelec_setGamma.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double)]

        self.lib.Getelec_setGamma.restype = None
        self.lib.Getelec_setGamma.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]

        self.lib.Getelec_setkT.restype = None
        self.lib.Getelec_setkT.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]

        self.lib.Getelec_setWorkFunction.restype = None
        self.lib.Getelec_setWorkFunction.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]

        self.lib.Getelec_setBandDepth.restype = None
        self.lib.Getelec_setBandDepth.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]

        self.lib.Getelec_setEffectiveMass.restype = None
        self.lib.Getelec_setEffectiveMass.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]

        self.lib.Getelec_run.restype = ctypes.c_size_t
        self.lib.Getelec_run.argtypes = [ctypes.c_void_p, ctypes.c_bool]

        self.lib.Getelec_getCurrentDensities.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_getCurrentDensities.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t)]

        self.lib.Getelec_getNottinghamHeats.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_getNottinghamHeats.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t)]

        self.lib.Getelec_getSpectraEnergies.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_getSpectraEnergies.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t)]

        self.lib.Getelec_getSpectraValues.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_getSpectraValues.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t)]

        self.lib.Getelec_getSpectraDerivatives.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_getSpectraDerivatives.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t)]

        self.lib.Getelec_calculateTransmissionCoefficientForEnergy.restype = ctypes.c_double
        self.lib.Getelec_calculateTransmissionCoefficientForEnergy.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_size_t]

        self.lib.Getelec_calculateTransmissionCoefficientForEnergies.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_calculateTransmissionCoefficientForEnergies.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_size_t]

        self.lib.Getelec_calculateTransmissionCoefficientForManyEnergies.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_calculateTransmissionCoefficientForManyEnergies.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_size_t]

        self.lib.Getelec_delete.restype = None
        self.lib.Getelec_delete.argtypes = [ctypes.c_void_p]

        self.lib.Getelec_getBarrierIntegrationLimits.restype = None
        self.lib.Getelec_getBarrierIntegrationLimits.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]

        self.lib.Getelec_getBarrierValues.restype = None
        self.lib.Getelec_getBarrierValues.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_size_t]

    def calculateTransmissionForEnergy(self, energy, params_index):
        return self.lib.Getelec_calculateTransmissionCoefficientForEnergy(self.obj, energy, params_index)

    def calculateTransmissionForEnergies(self, energies, params_index):
        energies_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(energies, dtype=np.float64)))
        result_ptr = self.lib.Getelec_calculateTransmissionCoefficientForEnergies(self.obj, energies_ctypes, size, params_index)
        return np.ctypeslib.as_array(result_ptr, shape=(size,))

    def calculateTransmissionForManyEnergies(self, energies, params_index):
        energies_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(energies, dtype=np.float64)))
        result_ptr = self.lib.Getelec_calculateTransmissionCoefficientForManyEnergies(self.obj, energies_ctypes, size, params_index)
        return np.ctypeslib.as_array(result_ptr, shape=(size,))

    def _toCtypesArray(self, numpy_array):
        return numpy_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(numpy_array)

    def setField(self, fields):
        fields_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(fields, dtype=np.float64)))
        self.lib.Getelec_setField(self.obj, fields_ctypes, size)

    def setRadius(self, radii):
        radii_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(radii, dtype=np.float64)))
        self.lib.Getelec_setRadius(self.obj, radii_ctypes, size)
    
    def setGamma(self, gammas):
        gammas_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(gammas, dtype=np.float64)))
        self.lib.Getelec_setGamma(self.obj, gammas_ctypes, size)
    
    def set_kT(self, kTs):
        kTs_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(kTs, dtype=np.float64)))
        self.lib.Getelec_setkT(self.obj, kTs_ctypes, size)
    
    def setWorkFunction(self, workFunctions):
        workFunctions_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(workFunctions, dtype=np.float64)))
        self.lib.Getelec_setWorkFunction(self.obj, workFunctions_ctypes, size)
    
    def setBandDepth(self, bandDepths):
        bandDepths_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(bandDepths, dtype=np.float64)))
        self.lib.Getelec_setBandDepth(self.obj, bandDepths_ctypes, size)
    
    def setEffectiveMass(self, effectiveMasses):
        effectiveMasses_ctypes, size = self._toCtypesArray(np.atleast_1d(np.asarray(effectiveMasses, dtype=np.float64)))
        self.lib.Getelec_setEffectiveMass(self.obj, effectiveMasses_ctypes, size)
    

    def run(self, calculate_spectra=False):
        self.numberOfIterations = self.lib.Getelec_run(self.obj, calculate_spectra)

    def getCurrentDensity(self):
        size = ctypes.c_size_t()
        densities_ptr = self.lib.Getelec_getCurrentDensities(self.obj, ctypes.byref(size))
        return np.ctypeslib.as_array(densities_ptr, shape=(size.value,))
    
    def getNottinghamHeat(self):
        size = ctypes.c_size_t()
        densities_ptr = self.lib.Getelec_getNottinghamHeats(self.obj, ctypes.byref(size))
        return np.ctypeslib.as_array(densities_ptr, shape=(size.value,)) 
    
    def getSpectra(self):
        self.spectra = []
        for i in range(self.numberOfIterations):
            size = ctypes.c_size_t()

            energies_ptr = self.lib.Getelec_getSpectraEnergies(self.obj, ctypes.c_size_t(i), ctypes.byref(size))
            energies = np.ctypeslib.as_array(energies_ptr, shape=(size.value,))

            values_ptr = self.lib.Getelec_getSpectraValues(self.obj, ctypes.c_size_t(i) , ctypes.byref(size))
            values = np.ctypeslib.as_array(values_ptr, shape=(size.value,))

            derivatives_ptr = self.lib.Getelec_getSpectraDerivatives(self.obj, ctypes.c_size_t(i) , ctypes.byref(size))
            derivatives = np.ctypeslib.as_array(derivatives_ptr, shape=(size.value,))

            validIndices = np.where(np.isfinite(energies) & np.isfinite(values) & np.isfinite(derivatives) & (values > max(values) * 1e-5))

            self.spectra.append(spi.CubicHermiteSpline(energies[validIndices], values[validIndices], derivatives[validIndices]))
        return self.spectra
    
    def calculateBarrierValues(self, energies, params_index = 0):
        energyArray = np.atleast_1d(np.asarray(energies, dtype=np.float64))
        potentialArray = np.copy(energyArray)
        energies_ctypes, size = self._toCtypesArray(energyArray)
        potential_ctypes, size = self._toCtypesArray(potentialArray)

        self.lib.Getelec_getBarrierValues(self.obj, energies_ctypes, potential_ctypes, size, params_index)
        return np.ctypeslib.as_array(potential_ctypes, shape=(size,))
    
    def calculateCurrentDensity(self):
        self.run()
        return self.getCurrentDensity()
    
    def calculateNottinghamHeat(self):
        self.run()
        return self.getNottinghamHeat()
    
    
    def getBarrierIntegrationLimits(self, params_index = 0):
        lowerLimit = ctypes.c_double()
        upperLimit = ctypes.c_double()
        self.lib.Getelec_getBarrierIntegrationLimits(self.obj, ctypes.byref(lowerLimit), ctypes.byref(upperLimit))
        return lowerLimit.value, upperLimit.value
    
    def getBarrierPlotData(self, Npoints = 512):
        lowerLimit, upperLimit = self.getBarrierIntegrationLimits()
        energies = np.linspace(lowerLimit, upperLimit, Npoints)
        values = self.calculateBarrierValues(energies)
        return energies, values
        
    
    def setRandomParams(self, numberOfInputs):
        fields = np.random.uniform(1., 12., numberOfInputs)
        radii = 1./np.random.uniform(0., 1., numberOfInputs)
        gammas = np.random.uniform(1., 100., numberOfInputs)
        kTs = np.random.uniform(0., 1., numberOfInputs)
        workFunctions = np.random.uniform(1., 6., numberOfInputs)
        bandDepths = np.random.uniform(1., 15., numberOfInputs)
        effectiveMasses = np.random.uniform(0.01, 2., numberOfInputs)
        self.setField(fields)
        self.setRadius(radii)
        self.setGamma(gammas)
        self.set_kT(kTs)
        self.setWorkFunction(workFunctions)
        self.setBandDepth(bandDepths)
        self.setEffectiveMass(effectiveMasses)

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
    getelec = GetelecInterface()
    getelec.setRadius([5., 6., 7. , 8., 9.])
    # getelec.setRandomInputs(5)
    # getelec.set_radius(5.)

    getelec.run(calculate_spectra=True)
    print(getelec.getCurrentDensity())

    spectra = getelec.getSpectra()

    import matplotlib.pyplot as plt
    for spline in spectra:
        xPlot = np.linspace(min(spline.x), max(spline.x), 512)
        yPlot = spline(xPlot)
        plt.plot(xPlot, yPlot)
        print(np.trapz(yPlot, xPlot))

    energies, values = getelec.getBarrierPlotData()

    plt.figure()
    plt.plot(energies, values)

    plt.show()
