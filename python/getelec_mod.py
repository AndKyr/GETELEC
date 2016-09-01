import ctypes as ct
import numpy as np
import scipy.optimize as opt
import os

pythonpath,filename = os.path.split(os.path.realpath(__file__))
emissionpath,pythonfolder = os.path.split(pythonpath)
libpath = emissionpath + '/lib/libgetelec.so'
#libslatecpath = emissionpath + '/lib/libslatec.so'
#ct.cdll.LoadLibrary(libslatecpath)
print libpath
ct.cdll.LoadLibrary(libpath)
getelec = ct.CDLL(libpath)

class Emission(ct.Structure):
    """This contains all the interchange information for emission calculation
    with getelec. To create it use the standard constructor. The xr and Vr attributes
    are passed via a numpy array as np.ctypes.as_ctypes(xr)."""
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
    
    def print_data(self,full = False):#print the data 
        if (full):
            return getelec.print_data_c(ct.byref(self),ct.c_int(1))
        else:
            return getelec.print_data_c(ct.byref(self),ct.c_int(0))
    
    def plot_data(self):
        return getelec.plot_data_c(ct.byref(self))

def emit (F = 5., W = 4.5, R = 5., gamma = 10., Temp = 300.):
    this =   Emission(F,W,R,gamma,Temp,0.,0.,None,None,'r','r',0,1,0,0)
    this.cur_dens()
    if (this.ierr == 0):
        return (this.Jem, this.heat)
    else:
        print 'Error:', this.ierr
        this.print_data()
        return (this.Jem, this.heat)
        
def MLplot (xfn, beta = 1., W = 4.5, R = 5., gamma = 10., Temp = 300.):
    yfn = np.empty([len(xfn)])
    for i in range(len(xfn)):
        F = beta / xfn[i]
        out = emit(F,W,R,gamma,Temp)
        yfn[i] = np.log(out[0])    
    return yfn

def MLerror(p, xML, yML):
    ycalc = MLplot(xML, p[0], p[1], p[2], p[3], p[4])
    yshift = max(ycalc) - max(yML)
    return ycalc - yML - yshift
    
def fitML (xML, yML, F0 = [0.05, 5., 18.], W0 = [2.5, 4., 5.5], \
            R0 = [1., 5., 30.], gamma0 = [2., 10., 100.], \
            Temp0 = [100., 300., 1080. ]):
    """ tries to fit xML, yML experimental data to the GETELEC model """
    meanFmac = np.mean(1./xML)
    minFmac = min(1./xML)
    maxFmac = max(1./xML)
    beta0 = [F0[0]/minFmac, F0[1] / meanFmac, F0[2] / maxFmac]
    
    popt = opt.least_squares (fun = MLerror, args = (xML, yML), \
                x0 =     [beta0[1], W0[1], R0[1], gamma0[1], Temp0[1]], \
                bounds =([beta0[0], W0[0], R0[0], gamma0[0], Temp0[0]],\
                        [beta0[2], W0[2], R0[2], gamma0[2], Temp0[2]]), \
                method = 'trf', jac = '3-point',xtol = 1.e-15 )
    return popt



    
