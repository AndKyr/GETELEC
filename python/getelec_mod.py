"""This module provides with all the interfcace classes and functions to call
GETELEC software from python and perform electron emissino calculations"""


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
        """Calculate current density"""
        return getelec.cur_dens_c(ct.byref(self))
    
    def print_data(self,full = False):#print the data
        """Calculate current density and print data""" 
        if (full):
            return getelec.print_data_c(ct.byref(self),ct.c_int(1))
        else:
            return getelec.print_data_c(ct.byref(self),ct.c_int(0))
    
    def plot_data(self):
        """Calculate current density and plot the barrier"""
        return getelec.plot_data_c(ct.byref(self))
        
def emission_create(F = 5., W = 4.5, R = 5., gamma = 10., Temp = 300., \
                Jem = 0., heat = 0., xr = np.array([]), Vr = np.array([]), \
                regime = 'r', sharp = 'r', full = 1, mode = 0, ierr = 0):
                    
    """Creates an Emission class object and calculates everyting. 
    Input xr ,Vr are given as numpy arrays."""
    x_c = xr.ctypes.data_as(ct.POINTER(ct.c_double))
    V_c = Vr.ctypes.data_as(ct.POINTER(ct.c_double))
    Nr = xr.size
    assert (Vr.size == Nr),"Check sizes of xr and Vr"
    this = Emission(F,W,R,gamma,Temp,Jem,heat,x_c,V_c,regime,sharp,Nr,full,mode,ierr)
    this.cur_dens()
    return this

def emit (F = 5., W = 4.5, R = 5., gamma = 10., Temp = 300.):
    """Calculate the current density and the Nottingham heating for specific set
    of input parameters"""
    this =   emission_create(F,W,R,gamma,Temp)
    if (this.ierr == 0):
        return (this.Jem, this.heat)
    else:
        print 'Error:', this.ierr
        this.print_data()
        return (this.Jem, this.heat)
        
def MLplot (xfn, beta = 1., W = 4.5, R = 5., gamma = 10., Temp = 300.):
    """Take input data xfn=1/V (x axis of FN plot) and calculate the current density
    for a specific set of parameters and F = beta * V. Output in logarithmic scale """
    yfn = np.empty([len(xfn)])
    for i in range(len(xfn)):
        F = beta / xfn[i]
        out = emit(F,W,R,gamma,Temp)
        yfn[i] = np.log(out[0])    
    return yfn

def MLerror(p, xML, yML):
    """calculate the y-shift between theoretical and experimental"""
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
    

def theta_SC(J,V,F):
    
    k = 1.904e5
    z = k * V**0.5 * J / F**2
    out = 1. - 1.33333 * z + 3. * z**2 - 8 * z**3
    return min(max(out,0.01), 1.)

def emit_SC(F = 10., W = 4.5, R = 5., gamma = 10., Temp = 300., \
                V_appl = 5.e3, err_fact = 0.2):
    """Calculate the current density and the Nottingham heating for specific set
    of input parameters, taking into account the space_charge effect"""
    F_p = F
    this = emission_create(F_p,W,R,gamma,Temp)
    theta_old = F_p / F
    
    for j in range(20):
        for i in range(20):
            this.F = F_p
            this.cur_dens()
            J_p = this.Jem
            if (i == 0 and j==0): J0 = J_p
            theta_new = theta_SC(J_p, V_appl, F_p)
            error =  (theta_new - theta_old)
            theta_new = theta_old + error * err_fact
            if (abs(error) < 1.e-5): break 
            theta_old = theta_new 
            #print "F = %f, J = %e, theta = %e" %(F_p, J_p, theta_new)
            F_p = F * theta_new
        err_fact = err_fact * 0.5
    
    #print 'J0 = %e, Jfinal = %e, reduction = %f' %(J0, J_p, J_p/J0)
    return J_p,theta_new



    
