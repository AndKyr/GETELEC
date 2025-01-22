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

        self.lib.Getelec_run.restype = None
        self.lib.Getelec_run.argtypes = [ctypes.c_void_p, ctypes.c_bool]

        self.lib.Getelec_getCurrentDensities.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.Getelec_getCurrentDensities.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t)]

    def _to_ctypes_array(self, numpy_array):
        return numpy_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(numpy_array)

    def set_field(self, fields):
        fields_ctypes, size = self._to_ctypes_array(np.asarray(fields, dtype=np.float64))
        self.lib.Getelec_setField(self.obj, fields_ctypes, size)

    def set_radius(self, radii):
        radii_ctypes, size = self._to_ctypes_array(np.asarray(radii, dtype=np.float64))
        self.lib.Getelec_setRadius(self.obj, radii_ctypes, size)

    def run(self, calculate_spectra=False):
        self.lib.Getelec_run(self.obj, calculate_spectra)

    def get_current_densities(self):
        size = ctypes.c_size_t()
        densities_ptr = self.lib.Getelec_getCurrentDensities(self.obj, ctypes.byref(size))
        return np.ctypeslib.as_array(densities_ptr, shape=(size.value,))

    def __del__(self):
        # Cleanup the object
        if hasattr(self, "lib") and self.lib:
            self.lib.Getelec_delete(self.obj)

# Example usage:
# Assuming the shared library is compiled as "getelec.so"
getelec = GetelecInterface("build/libgetelec.so")
getelec.set_field([5.0, 6.0, 7.0])
# getelec.set_radius([1000.0, 2000.0])
getelec.run()
print(getelec.get_current_densities())
