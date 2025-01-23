import ctypes
import numpy as np

class GetelecInterface:
    def __init__(self, library_path):
        # Load the shared library
        self.lib = ctypes.CDLL(library_path)

        # Initialize the Getelec instance
        self.obj = self.lib.Getelec_new()

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

        self.lib.Getelec_calculateTransmissionCoefficientForEnergy.restype = ctypes.c_double
        self.lib.Getelec_calculateTransmissionCoefficientForEnergy.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_size_t]

        self.lib.Getelec_calculateTransmissionCoefficientForEnergies.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_calculateTransmissionCoefficientForEnergies.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_size_t]

        self.lib.Getelec_calculateTransmissionCoefficientForManyEnergies.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_calculateTransmissionCoefficientForManyEnergies.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_size_t]

    def calculate_transmission_coefficient_for_energy(self, energy, params_index):
        return self.lib.Getelec_calculateTransmissionCoefficientForEnergy(self.obj, energy, params_index)

    def calculate_transmission_coefficient_for_energies(self, energies, params_index):
        energies_ctypes, size = self._to_ctypes_array(np.atleast_1d(np.asarray(energies, dtype=np.float64)))
        result_ptr = self.lib.Getelec_calculateTransmissionCoefficientForEnergies(self.obj, energies_ctypes, size, params_index)
        return np.ctypeslib.as_array(result_ptr, shape=(size,))

    def calculate_transmission_coefficient_for_many_energies(self, energies, params_index):
        energies_ctypes, size = self._to_ctypes_array(np.atleast_1d(np.asarray(energies, dtype=np.float64)))
        result_ptr = self.lib.Getelec_calculateTransmissionCoefficientForManyEnergies(self.obj, energies_ctypes, size, params_index)
        return np.ctypeslib.as_array(result_ptr, shape=(size,))

    def _to_ctypes_array(self, numpy_array):
        return numpy_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(numpy_array)

    def set_field(self, fields):
        fields_ctypes, size = self._to_ctypes_array(np.atleast_1d(np.asarray(fields, dtype=np.float64)))
        self.lib.Getelec_setField(self.obj, fields_ctypes, size)

    def set_radius(self, radii):
        radii_ctypes, size = self._to_ctypes_array(np.atleast_1d(np.asarray(radii, dtype=np.float64)))
        self.lib.Getelec_setRadius(self.obj, radii_ctypes, size)
    
    def set_gamma(self, gammas):
        gammas_ctypes, size = self._to_ctypes_array(np.atleast_1d(np.asarray(gammas, dtype=np.float64)))
        self.lib.Getelec_setGamma(self.obj, gammas_ctypes, size)
    
    def set_kT(self, kTs):
        kTs_ctypes, size = self._to_ctypes_array(np.atleast_1d(np.asarray(kTs, dtype=np.float64)))
        self.lib.Getelec_setkT(self.obj, kTs_ctypes, size)
    
    def set_workFunction(self, workFunctions):
        workFunctions_ctypes, size = self._to_ctypes_array(np.atleast_1d(np.asarray(workFunctions, dtype=np.float64)))
        self.lib.Getelec_setWorkFunction(self.obj, workFunctions_ctypes, size)
    
    def set_bandDepth(self, bandDepths):
        bandDepths_ctypes, size = self._to_ctypes_array(np.atleast_1d(np.asarray(bandDepths, dtype=np.float64)))
        self.lib.Getelec_setBandDepth(self.obj, bandDepths_ctypes, size)
    
    def set_effectiveMass(self, effectiveMasses):
        effectiveMasses_ctypes, size = self._to_ctypes_array(np.atleast_1d(np.asarray(effectiveMasses, dtype=np.float64)))
        self.lib.Getelec_setEffectiveMass(self.obj, effectiveMasses_ctypes, size)
    

    def run(self, calculate_spectra=False):
        self.numberOfIterations = self.lib.Getelec_run(self.obj, calculate_spectra)

    def get_current_densities(self):
        size = ctypes.c_size_t()
        densities_ptr = self.lib.Getelec_getCurrentDensities(self.obj, ctypes.byref(size))
        return np.ctypeslib.as_array(densities_ptr, shape=(size.value,))
    
    def get_nottingham_heats(self):
        size = ctypes.c_size_t()
        densities_ptr = self.lib.Getelec_getNottinghamHeats(self.obj, ctypes.byref(size))
        return np.ctypeslib.as_array(densities_ptr, shape=(size.value,)) 
    
    def get_spectra(self):
        spectra = []
        for i in range(self.numberOfIterations):
            size = ctypes.c_size_t()

            energies_ptr = self.lib.Getelec_getSpectraEnergies(self.obj, ctypes.c_size_t(i), ctypes.byref(size))
            energies = np.ctypeslib.as_array(energies_ptr, shape=(size.value,))

            values_ptr = self.lib.Getelec_getSpectraValues(self.obj, ctypes.c_size_t(i) , ctypes.byref(size))
            values = np.ctypeslib.as_array(values_ptr, shape=(size.value,))
            spectra.append((energies, values))
        return spectra
    
    def setRandomInputs(self, numberOfInputs):
        fields = np.random.uniform(1., 12., numberOfInputs)
        radii = 1./np.random.uniform(0., 1., numberOfInputs)
        gammas = np.random.uniform(1., 100., numberOfInputs)
        kTs = np.random.uniform(0., 1., numberOfInputs)
        workFunctions = np.random.uniform(1., 6., numberOfInputs)
        bandDepths = np.random.uniform(1., 15., numberOfInputs)
        effectiveMasses = np.random.uniform(0.01, 2., numberOfInputs)
        self.set_field(fields)
        self.set_radius(radii)
        self.set_gamma(gammas)
        self.set_kT(kTs)
        self.set_workFunction(workFunctions)
        self.set_bandDepth(bandDepths)
        self.set_effectiveMass(effectiveMasses)


    def __del__(self):
        # Cleanup the object
        if hasattr(self, "lib") and self.lib:
            self.lib.Getelec_delete(self.obj)

# Example usage:
# Assuming the shared library is compiled as "getelec.so"
getelec = GetelecInterface("build/libgetelec.so")
getelec.setRandomInputs(5)
# getelec.set_radius(5.)

for i in range(5):
    print(getelec.calculate_transmission_coefficient_for_energy(0.1, i))

getelec.run(calculate_spectra=True)
print(getelec.get_current_densities())

spectra = getelec.get_spectra()

import matplotlib.pyplot as plt
for energies, values in spectra:
    plt.plot(energies, values)

plt.show()
