import ctypes as ct
import numpy as np

ct.cdll.LoadLibrary("lib/libgetelec.so")
getelec = ct.CDLL("lib/libgetelec.so")

class Emission(ct.Structure):
    _fields_ = [("F", ct.c_double),
                ("W", ct.c_double),
                ("R", ct.c_double),
                ("gamma", ct.c_double),
                ("Temp", ct.c_double),
                ("Jem", ct.c_double),
                ("heat", ct.c_double),
                ("xr", ct.POINTER(ct.c_double)),
                ("Vr", ct.POINTER(ct.c_double)),
                ("regime", ct.c_char),
                ("sharp", ct.c_char),
                ("Nr", ct.c_int),
                ("full", ct.c_int),
                ("mode", ct.c_int),
                ("ierr", ct.c_int)
                ]
                
    def cur_dens(self):
        return getelec.cur_dens_c(ct.byref(self))
    
    def print_data(self):
        return getelec.print_data_c(ct.byref(self))
        
    def plot_data(self):
        return getelec.plot_data(ct.byref(self))

 
this = Emission(5.,4.7,5.,10.,500.,0.,0.,None,None,'t','t', 0, 0, 0, 0)



this.cur_dens()
print this.Jem, this.heat
this.print_data()
