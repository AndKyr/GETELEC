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

    def _to_ctypes_array(self, numpy_array):
        return numpy_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(numpy_array)

    def set_field(self, fields):
        fields_ctypes, size = self._to_ctypes_array(np.asarray(fields, dtype=np.float64))
        self.lib.Getelec_setField(self.obj, fields_ctypes, size)

    def set_radius(self, radii):
        radii_ctypes, size = self._to_ctypes_array(np.asarray(radii, dtype=np.float64))
        self.lib.Getelec_setRadius(self.obj, radii_ctypes, size)

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

    def __del__(self):
        # Cleanup the object
        if hasattr(self, "lib") and self.lib:
            self.lib.Getelec_delete(self.obj)

# Example usage:
# Assuming the shared library is compiled as "getelec.so"
getelec = GetelecInterface("build/libgetelec.so")
getelec.set_field(np.linspace(5., 6., 8))
# getelec.set_radius([1000.0, 2000.0])
getelec.run(calculate_spectra=True)
print(getelec.get_current_densities())

spectra = getelec.get_spectra()

import matplotlib.pyplot as plt
for energies, values in spectra:
    plt.plot(energies, values)

plt.show()
