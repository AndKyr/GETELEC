"""This module provides with all the interfcace classes and functions to call
GETELEC software from python and perform electron emissino calculations"""


import ctypes as ct
import numpy as np
import scipy.optimize as opt
import scipy.integrate as ig
import os

pythonpath,filename = os.path.split(os.path.realpath(__file__))
emissionpath,pythonfolder = os.path.split(pythonpath)
libpath = emissionpath + '/lib/dynamic/libgetelec.so'
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
                ("regime", ct.c_int),
                ("sharp", ct.c_int),
                ("Nr", ct.c_int),
                ("approx", ct.c_int),
                ("mode", ct.c_int),
                ("ierr", ct.c_int),
                ("voltage", ct.c_double),
                ("theta", ct.c_double)
                ]
    
    def cur_dens(self):
        """Calculate current density"""
        return getelec.cur_dens_c(ct.byref(self))
        
    def cur_dens_SC(self):
        return getelec.cur_dens_SC(ct.byref(self))
        
    
    def print_data(self,full = False):#print the data
        """Calculate current density and print data""" 
        if (full):
            return getelec.print_data_c(ct.byref(self),ct.c_int(1))
        else:
            return getelec.print_data_c(ct.byref(self),ct.c_int(0))
    
    def print_C_data(self):
        return getelec.print_C_data(ct.byref(self))
    
    def plot_data(self):
        """Calculate current density and plot the barrier"""
        return getelec.plot_data_c(ct.byref(self))
    
    def print_data(self,full = False):#print the data
        """Calculate current density and print data""" 
        if (full):
            return getelec.print_data_c(ct.byref(self),ct.c_int(1))
        else:
            return getelec.print_data_c(ct.byref(self),ct.c_int(0))
    
    def print_C_data(self):
        return getelec.print_C_data(ct.byref(self))
        
         
    
def emission_create(F = 5., W = 4.5, R = 5000., gamma = 1., Temp = 300., \
                Jem = 0., heat = 0., xr = np.array([]), Vr = np.array([]), \
                regime = 0, sharp = 1, approx = 1, mode = 0, ierr = 0, voltage = 500):
                    
    """Creates an Emission class object and calculates everyting. 
    Input xr ,Vr are given as numpy arrays."""
    x_c = xr.ctypes.data_as(ct.POINTER(ct.c_double))
    V_c = Vr.ctypes.data_as(ct.POINTER(ct.c_double))
    Nr = len(xr)
    assert (Vr.size == Nr),"Check sizes of xr and Vr"
    this = Emission(F,W,R,gamma,Temp,Jem,heat,x_c,V_c,regime,sharp,Nr, \
            approx,mode,ierr, voltage)
    this.cur_dens()
    return this

    # ode system for planar geometry    
def planar(y, x, kJ):
    dy = y[1]
    ddy = kJ / np.sqrt(y[0])
    return np.array([dy, ddy])

#ode system for cylindrical geometry    
def cylindrical(y, x, kJ):
    dy = y[1]
    ddy = (kJ / np.sqrt(y[0]) - y[1]) / x
    return np.array([dy, ddy])
    
#ode system for spherical geometry
def spherical(y, x, kJ):
    dy = y[1]
    ddy = (kJ / np.sqrt(y[0]) - 2 * y[1] * x) / x**2
    return np.array([dy, ddy])
    
""" This class accommodates all 1D space charge calculations"""        
class SpaceCharge():
    kappa = 1.904e5
        
    #initializes the class
    def __init__(self, Voltage = 100, Distance = 10, geo = "plane", 
                Npoints = 256, emitter = None, Radius = 1.):
        if (emitter == None):
            print "SpaceCharge class created with default emitter"
            self.emitter = emission_create()
        else:
            self.emitter =  emitter
        
        self.voltage = Voltage
        self.distance = Distance
        self.geometry = geo
        self.radius = Radius
                                                                                  
        if (self.geometry == "plane"):
            self.system = planar
            self.Radius = 1.
            self.beta = 1. / self.distance
        elif(self.geometry == "cylinder"):
            self.system = cylindrical
            self.beta = 1. / (self.radius * np.log(1 + self.distance / self.radius))
        elif (self.geometry == "sphere"):
            self.system = spherical
            self.beta = (self.radius + self.distance) / (self.radius * self.distance)
        else:
            print "Error: Wrong geometry. Should be either, plane, cylinder or sphere"
        
        self.Np = Npoints
        
        
        
    # solves the ode initial value problem with cathode field F and finds potential at d
    def Vanode(self, F):
        ## calculate current density for given cathode field F
        self.emitter.F = F
        self.emitter.cur_dens()
        self.Jc = self.emitter.Jem
        
        if (self.geometry == "plane"):
            factor = self.Jc * self.kappa
        elif (self.geometry == "sphere"):
            factor = self.kappa * self.Jc * self.radius**2
        elif (self.geometry == "cylinder"):
            factor = self.kappa * self.Jc * self.radius
        
        ## solve ode and calculate V(d)
        xs = np.linspace(self.radius, self.radius + self.distance, self.Np)
        y = ig.odeint(self.system, np.array([1.e-10,F]), xs, (factor,))
        return y[-1, 0]
        
    # uses the bisection method to calculate F(V, d)    
    def calcJF(self, V):
        rtol = 1.e-4
        Fmin = 0.
        
        Fmax = self.beta * V
        for i in range(50):
            Fi = .5 * (Fmax + Fmin)
            Vi = self.Vanode(Fi)
            if (Vi < (1 - rtol) * V):
                Fmin = Fi
            elif((Vi > (1 + rtol) * V)):
                Fmax = Fi
            else:
                break
            
            print "Fi = %.3f, Vi = %.3f, Ji = %e"%(Fi, Vi, self.Jc)
        
        return Fi

        
            
            
    

def emit (F = 5., W = 4.5, R = 5., gamma = 10., Temp = 300., verbose = False):
    """Calculate the current density and the Nottingham heating for specific set
    of input parameters"""
    this =   emission_create(F,W,R,gamma,Temp)
    if (this.ierr == 0 or isnan(this.Jem) or isnan(this.heat) ):
        if (verbose): this.print_data()
        return (this.Jem, this.heat)
    else:
        print 'Error:', this.ierr
        this.print_data()
        return (this.Jem, this.heat)
        
def MLplot (xfn, beta = 1., W = 4.5, R = 5., gamma = 10., Temp = 300.):
    """Take input data xfn=1/V (x axis of FN plot) and calculate the current density
    for a specific set of parameters and F = beta * V. Output in logarithmic scale """
    yfn = np.empty([len(xfn)])
    this = emission_create(5,W,R,gamma,Temp)
    
    for i in range(len(xfn)):
        this.F = beta / xfn[i]
        this.cur_dens()
        yfn[i] = np.log(this.Jem)    
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
    if (F<=0.01): return 1.
    z = k * V**0.5 * J / F**2
    poly = np.array([9.*z**2, -3., -4.*z+3.])
    #print "poly = ", poly
    rts = np.roots(poly)
    rdist = abs(rts - (2./3))
    theta = rts[np.argmin(rdist)]
#    print "zeta = %e, theta = %e"%(z, theta)
    return theta

def emit_SC(F = 10., W = 4.5, R = 500., gamma = 10., Temp = 300., \
                V_appl = 5.e3, err_fact = 0.5, approx = 1):
    """Calculate the current density and the Nottingham heating for specific set
    of input parameters, taking into account the space_charge effect"""
    F_p = min(F, 15.)
    this = emission_create(F_p,W,R,gamma,Temp, approx = approx)
    theta_old = F_p / F
    
    for j in range(50):
        this.F = F_p
        this.cur_dens()
        J_p = this.Jem
        if (j==0): 
            J0 = J_p
        theta_new = theta_SC(J_p, V_appl, F_p)
        print "F = %f, J = %e, theta = %e" %(F_p, J_p, theta_new)
        error =  (theta_new - theta_old)
        theta_new = theta_old + error * err_fact
        if (abs(error) < 1.e-6): break 
        theta_old = theta_new 

        F_p = F * theta_new
    
    #print 'J0 = %e, Jfinal = %e, reduction = %f' %(J0, J_p, J_p/J0)
    return J_p,theta_new
    
def theta_approx(J,V,F):
    
    k = 1.904e5
    if (F<=0.01): return 1.
    z = k * V**0.5 * J / F**2

    out = 1. - 1.33333 * z + 3. * z**2 - 8 * z**3
    if (out > 0 and out <= 1.):
        return out
    else:
        return 1.
        
def thetaroots(z):
    poly = np.array([9.*z**2, -3., -4.*z+3.])
    #print "poly = ", poly
    out = np.zeros(2)
    rts = np.roots(poly)
    if len(rts == 2):
        return rts
    else:
        out[:] = rts[0]
        return out
        
def J_ChildL(V,beta):
    kappa = 1.904e5
    return (4./9) * beta**2 * V**1.5 / kappa
    
    
        




    
