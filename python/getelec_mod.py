"""This module provides with all the interfcace classes and functions to call
GETELEC software from python and perform electron emissino calculations"""


import ctypes as ct
import numpy as np
import scipy.optimize as opt
import scipy.integrate as ig
import os
import matplotlib.pyplot as plt

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
        
    
    #calculates the total current, assuming a spherical tip, where the field falls as F(theta = 0)(cos(kappa * theta))
    def total_current(self, kappa):
        theta = np.linspace(0,np.pi/2, 45)
        Fi = self.F * np.cos(kappa * theta)
        Ji = np.copy(Fi)
        
        for i in range(len(Fi)):
            self.F = Fi[i]
            self.cur_dens()
            Ji[i] = self.Jem
        
        self.F = Fi[0]
        
        # plt.plot(theta, Ji)
        # plt.plot(theta, 2 * np.pi * self.R**2 * np.sin(theta) * Ji)
        # plt.show()
        
        I = np.trapz(2 * np.pi * self.R**2 * np.sin(theta) * Ji, theta)
        Area = I / Ji[0]
        
        return I, Area
        
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




class MultiEmitter():
    
    def __init__(self, em = emission_create(), N = 1, mu_h = 100., std_h = 10., \
                    mu_r = 10., std_r = 1.):
        self.emitter = em
        self.N_emitters = N
        self.mu_height = mu_h
        self.sigma_height = std_h
        self.mu_radii = mu_r
        self.sigma_radii = std_r
        
        
    def get_emitters(self):
        mu = np.log(self.mu_height**2/np.sqrt(self.sigma_height**2 + self.mu_height**2))
        sigma = np.sqrt(np.log(1+self.sigma_height**2/self.mu_height**2))
        
        self.heights = np.random.lognormal(mu, sigma, self.N_emitters)
        
        mu = np.log(self.mu_radii**2/np.sqrt(self.sigma_radii**2 + self.mu_radii**2))
        sigma = np.sqrt(np.log(1+self.sigma_radii**2/self.mu_radii**2))
        self.radii = np.random.lognormal(mu, sigma, self.N_emitters)
        self.betas = self.heights / self.radii
        self.currents = np.copy(self.betas) * 0.
        
           
    def emit(self, Ffar):
        for i in range(len(self.betas)):
            self.emitter.F = self.betas[i] * Ffar
            self.emitter.R = self.radii[i]
            self.currents[i] =  self.emitter.total_current(0.5)[0]
            
        return sum(self.currents)    
        
        
    def plot_histograms(self):
        
        plt.figure()
        plt.hist(self.radii, 100)
        plt.xlabel("Radii")
        plt.ylabel("count")
        
        plt.figure()
        plt.hist(self.heights, 100)
        plt.xlabel("Heights")
        plt.ylabel("count")
        
        plt.figure()
        plt.hist(self.betas, 100)
        plt.xlabel(r"$\beta$")
        plt.ylabel("count")
        
        
        plt.figure()
        counts,Iedges = np.histogram(np.log(self.currents), 20, range = (np.log(1.e-9), np.log(max(self.currents))))
        
        print counts
        print Iedges
        Ibins = .5*(np.exp(Iedges[1:]) + np.exp(Iedges[:-1]))
        plt.semilogx(Ibins, counts * Ibins)
        plt.xlabel(r"I [Amps])")
        plt.ylabel("current contribution")
        plt.show()
        
    
    def print_data_out(self):
        
        # for i in range(len(betas)):
            # print heights[i], radii[i], betas[i]
            
        print "heights = %.3g +- %.3g"%(np.mean(self.heights), np.std(self.heights))
        print "radii = %.3g +- %.3g", np.mean(self.radii), "std<radii> =  ", np.std(self.radii)
        print "betas = %.3g +- %.3g"%(np.mean(self.betas),np.std(self.betas))
        print "maxbeta = %g"%(np.max(betas))
        

        
  
 
        
        
        
        
        


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
    
def IFplot(xfn, beta = 1., W = 4.5, R = 5., gamma = 10., Temp = 300., kappa = 0.5):
    """Take input data xfn=1/V (x axis of FN plot) and calculate the current 
        for a specific set of parameters and F = beta * V. Output in logarithmic scale """
    yfn = np.empty([len(xfn)])
    this = emission_create(5,W,R,gamma,Temp)
    
    for i in range(len(xfn)):
        this.F = beta / xfn[i]
        yfn[i] = np.log(this.total_current(kappa)[0])    
    return yfn  
    
    
def IFerror(p, xML, yML):
    """calculate the y-shift between theoretical and experimental"""
    ycalc = IFplot(xML, p[0], p[1], p[2], p[3], p[4], p[5])
    yshift = max(ycalc) - max(yML)
    return ycalc - yML - yshift
    
def fit_IFplot (xML, yML, F0 = [0.05, 5., 18.], W0 = [2.5, 4., 5.5], \
            R0 = [1., 5., 30.], gamma0 = [2., 10., 100.], \
            Temp0 = [100., 300., 1080. ], kappa = [0.4, 0.5, 0.6]):
    """ tries to fit xML, yML experimental data to the GETELEC model """
    meanFmac = np.mean(1./xML)
    minFmac = min(1./xML)
    maxFmac = max(1./xML)
    beta0 = [F0[0]/minFmac, F0[1] / meanFmac, F0[2] / maxFmac]
    
    popt = opt.least_squares (fun = IFerror, args = (xML, yML), \
                x0 =     [beta0[1], W0[1], R0[1], gamma0[1], Temp0[1], kappa[1]], \
                bounds =([beta0[0], W0[0], R0[0], gamma0[0], Temp0[0], kappa[0]],\
                        [beta0[2], W0[2], R0[2], gamma0[2], Temp0[2], kappa[2]]), \
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
    
    
        




    
