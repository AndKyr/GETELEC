"""This module provides with all the interfcace classes and functions to call
GETELEC software from python and perform electron emissino calculations"""


import ctypes as ct
import numpy as np
from numpy.lib.polynomial import RankWarning
import scipy.optimize as opt
import scipy.integrate as ig
import os
import matplotlib.pyplot as plt
import scipy.interpolate as intrp
import json

import warnings

from io import StringIO 
import sys

pythonpath,filename = os.path.split(os.path.realpath(__file__))
emissionpath,pythonfolder = os.path.split(pythonpath)
libpath = emissionpath + '/lib/dynamic/libgetelec.so'
#libslatecpath = emissionpath + '/lib/libslatec.so'
#ct.cdll.LoadLibrary(libslatecpath)
print(libpath)
ct.cdll.LoadLibrary(libpath)
getelec = ct.CDLL(libpath)

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout


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
        with Capturing() as output:
            if (full):
                getelec.print_data_c(ct.byref(self),ct.c_int(1))
            else:
                getelec.print_data_c(ct.byref(self),ct.c_int(0))
        print ("output = ", output)
        return output
    
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
        theta = np.linspace(0,np.pi/2, 64)
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

def get_gamow_line(F, R, gamma, Npoints):
    Wmin = np.array([1.])
    Wmax = np.copy(Wmin)
    Gamow = np.zeros(Npoints)
    getelec.export_gamow(ct.c_double(F), ct.c_double(R), ct.c_double(gamma), ct.c_int(Npoints), \
        ct.c_void_p(Wmin.ctypes.data), ct.c_void_p(Wmax.ctypes.data), ct.c_void_p(Gamow.ctypes.data))
    W = np.linspace(Wmin[0], Wmax[0], Npoints)
    return W, Gamow
    
def make_gamow_array(Frange, Rmin, Nf, Nr, Ngam, Npoly):
    warnings.filterwarnings("error")
    Finv = np.linspace(1/Frange[1], 1./ Frange[0], Nf)
    Rinv = np.linspace(1.e-5, 1/Rmin, Nr)
    gaminv = np.linspace(1.e-5, .999, Ngam)
    outarray = np.ones([Ngam, Nr, Nf, Npoly])

    for i in range(Ngam):
        for j in range(Nr):
            for k in range(Nf):
                F = 1/Finv[k]
                R = 1/Rinv[j]
                gamma = 1/gaminv[i]
                W,G = get_gamow_line(F, R , gamma , 128)
                try:
                    poly = np.polyfit(W, G, Npoly-1)
                except(RankWarning):
                    print("Rank Warning for F = %g, R = %g, gamma = %g"%(F,R,gamma))
                    plt.plot(W,G)
                    plt.show()
                outarray[i,j,k,:] = poly
    return outarray

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



def calc_json(json_str):
    """Creates an Emission class object and calculates everyting. 
    from json imput object."""

    try:
        F = float(json_str["Field"])
    except(KeyError):
        F = 5.

    try:
        W = float(json_str["Work_function"])
    except(KeyError):
        W = 4.5

    try:
        R = float(json_str["Radius"])
    except(KeyError):
        R = 2000.

    try:
        gamma = float(json_str["gamma"])
    except(KeyError):
        gamma = 20.
    try:
        Temp = float(json_str["Temperature"])
    except(KeyError):
        Temp = 300.

    try:
        voltage = float(json_str["voltage"])
    except(KeyError):
        voltage=0.
    
    try:
        approx = int(json_str["approximation"])
    except(KeyError):
        approx = 0


    this = emission_create(F,W,R,gamma,Temp, approx=approx, voltage=voltage)

    if (voltage > 0):
        this.cur_dens_SC()
    else:
        this.cur_dens()
                    

    outdata = {'current_density': this.Jem, 'heat': this.heat, \
                'regime': this.regime, 'sharpness': this.sharp, \
                'Ierror': this.ierr}
    return json.dumps(outdata)
    

    
class Tabulator():
    
    def __init__(self, emitter, Nf = 256, Nr = 256):
        self.emitter = emitter
        if (self.check_tabulated()):
            print ("loading data from files")
            self.load()
        else:
            print ("tabulating J, F, R")
            self.tabulate_JFR(Nf, Nr)

    def tabulate_JFR(self, Nf = 256, Nr = 256):
        self.Finv = np.linspace(1./20., 1., Nf)
        self.Rinv = np.linspace(1/100., 2., Nr)
        
        self.Jmesh = np.zeros((Nr, Nf))
        
        for i in range(Nr):
            self.emitter.R = 1./self.Rinv[i]
            for j in range(Nf):
                self.emitter.F = 1./self.Finv[j]
                self.emitter.Voltage = self.emitter.F * self.emitter.R * 0.2
                self.emitter.cur_dens_SC()
                self.Jmesh[i,j] = self.emitter.Jem
                if (self.emitter.Jem <= 0.):
                    self.emitter.print_data()
        try:
            os.mkdir("tabulated")
        except:
            pass
            
        np.save("tabulated/current_density", self.Jmesh)
        np.save("tabulated/Finv", self.Finv)
        np.save("tabulated/Rinv", self.Rinv)
        np.save("tabulated/checkcase",np.array([1./self.Finv[0], 1./self.Rinv[0], self.Jmesh[0,0]]))
        self.finterp = intrp.interp2d(self.Finv, self.Rinv, np.log(self.Jmesh))
    
    def check_tabulated(self):
        print ("checking tabulation")
        try:
            checks = np.load("tabulated/checkcase.npy")
        except(IOError):
            print ("tabulation check file not found")
            return False
            
        self.emitter.F = checks[0]
        self.emitter.R = checks[1]
        
        self.emitter.Voltage = self.emitter.F * self.emitter.R * 0.2
        self.emitter.cur_dens_SC()
        
        print ("read F, R, J = ", checks)
        print ("recalculated J = ", self.emitter.Jem)
        
        return self.emitter.Jem == checks[2]
        
    def load(self):
        self.Jmesh = np.load("tabulated/current_density.npy")
        self.Finv = np.load("tabulated/Finv.npy")
        self.Rinv = np.load("tabulated/Rinv.npy")
        self.finterp = intrp.interp2d(self.Finv, self.Rinv, np.log(self.Jmesh))
        
    def get_cur_dens(self, F, R):
        # assert (np.shape(F) == np.shape(R)), "F and R do not have the same shapes"
        
        Rcopy = np.copy(R)
        Fcopy = np.copy(F)
        
        Rcopy[Rcopy > 1./self.Rinv[0]] = 1./self.Rinv[0] 
        Rcopy[Rcopy < 1./self.Rinv[-1]] = 1./self.Rinv[-1]

        Fcopy[Fcopy > 1./self.Finv[0]] = 1./self.Finv[0]
        Fcopy[Fcopy < 1./self.Finv[-1]] = 1./self.Finv[-1]
        
        return np.exp(self.finterp(1./Fcopy, 1./Rcopy))
        
    #calculates the total current, assuming a spherical tip, where the field falls as F(theta = 0)(cos(kappa * theta))
    def total_current(self, kappa):
        theta = np.linspace(0,np.pi/2, 64)
        Fi = self.emitter.F * np.cos(kappa * theta)
        
        Ji = self.get_cur_dens(Fi, self.emitter.R)
               
        I = 2 * np.pi * self.emitter.R**2 * np.trapz(np.sin(theta) * Ji, theta)
        Area = I / Ji[0]
        
        return I, Area

class MultiEmitter():
    
    def __init__(self, em = emission_create(), N = 1, mu_h = 100., std_h = 10., \
                    mu_r = 10., std_r = 1.):
        self.emitter = em
        self.N_emitters = N
        self.mu_height = mu_h
        self.sigma_height = std_h
        self.mu_radii = mu_r
        self.sigma_radii = std_r
        self.tb = Tabulator(self.emitter, 256, 128)
        
        
    def set_emitters(self, radii, heights = np.array([]), betas = np.array([])):
        self.radii = radii
        
        if(not heights.any()):
            self.betas = betas
            self.heights = betas * radii
        else:
            self.betas = heights / radii
            self.heights = heights
        
        self.N_emitters = len(self.betas)
        
        self.currents = np.copy(self.betas) * 0.
        self.areas = np.copy(self.currents)
        
        
    def load_emitters(self, filename):
        data =  np.loadtxt(filename, unpack=True)
        self.heights = data[:,0]
        self.radii = data[:,1]
        self.betas = self.heights / self.radii
        
        self.currents = np.copy(self.betas) * 0.
        self.areas = np.copy(self.currents)
        
    def get_emitters(self, max_beta = 1.e50, min_beta = 1., h_dist = 'normal', r_dist = 'lognormal', addmax = True):
        
        if (h_dist == 'lognormal'):
            mu = np.log(self.mu_height**2/np.sqrt(self.sigma_height**2 + self.mu_height**2))
            sigma = np.sqrt(np.log(1+self.sigma_height**2/self.mu_height**2))
            heights = np.random.lognormal(mu, sigma, self.N_emitters)
        elif(h_dist == 'normal'):
            mu = self.mu_height
            sigma = self.sigma_height
            heights = np.random.normal(mu, sigma, self.N_emitters)
        elif(h_dist == 'uniform'):
            low = self.mu_height - np.sqrt(12) * self.sigma_height * 0.5
            high = self.mu_height + np.sqrt(12) * self.sigma_height * 0.5
            heights = np.random.uniform(low, high, self.N_emitters)
        elif(h_dist == 'exponential'):
            heights = np.random.exponential(self.mu_height, self.N_emitters)
        else:
            print ("wrong distribution for heights")
        
        if (r_dist == "lognormal"):
            mu = np.log(self.mu_radii**2/np.sqrt(self.sigma_radii**2 + self.mu_radii**2))
            sigma = np.sqrt(np.log(1+self.sigma_radii**2/self.mu_radii**2))
            radii = np.random.lognormal(mu, sigma, self.N_emitters)
        elif(r_dist == 'normal'):
            mu = self.mu_radii
            sigma = self.sigma_radii
            radii = np.random.normal(mu, sigma, self.N_emitters)
        elif(r_dist == 'uniform'):
            low = self.mu_radii - np.sqrt(12) * self.sigma_radii * 0.5
            high = self.mu_radii + np.sqrt(12) * self.sigma_radii * 0.5
            radii = np.random.uniform(low, high, self.N_emitters)
        elif(r_dist == 'exponential'):
            curvatures = np.random.exponential(1./self.mu_radii, self.N_emitters)
            radii = 1./curvatures
        else:
            print ("wrong distribution for radii"    )
            
        betas = heights / radii
        
        isgood = np.where(np.logical_and(betas < max_beta, betas > min_beta))
        
        
        self.radii = radii[isgood]
        self.heights = heights[isgood]
        self.betas = betas[isgood]
        self.N_emitters = len(self.betas)
        
        if (addmax):
            self.betas = np.append(self.betas, max_beta)
            self.radii = np.append(self.radii, min(self.radii))
            self.heights = np.append(self.heights, self.betas[-1] * self.radii[-1])
            self.N_emitters += 1
            
        
        self.currents = np.copy(self.betas) * 0.
        self.areas = np.copy(self.currents)
        
        
    def get_emitters_beta(self, min_beta = 50, max_beta = 150, beta_dist = 'uniform', r_dist = 'lognormal'):

        if(beta_dist == 'normal'):
            mu = 0.5 * (min_beta + max_beta)
            sigma = (max_beta - min_beta) * 0.25
            self.betas = np.random.normal(mu, sigma, self.N_emitters)
        elif(beta_dist == 'uniform'):
            low = self.mu_height - np.sqrt(12) * self.sigma_height * 0.5
            high = self.mu_height + np.sqrt(12) * self.sigma_height * 0.5
            self.betas = np.random.uniform(min_beta, max_beta, self.N_emitters)
        else:
            print ("wrong distribution for betas")
            
        if (r_dist == "lognormal"):
            mu = np.log(self.mu_radii**2/np.sqrt(self.sigma_radii**2 + self.mu_radii**2))
            sigma = np.sqrt(np.log(1+self.sigma_radii**2/self.mu_radii**2))
            self.radii = np.random.lognormal(mu, sigma, self.N_emitters)
        elif(r_dist == 'normal'):
            mu = self.mu_radii
            sigma = self.sigma_radii
            self.radii = np.random.normal(mu, sigma, self.N_emitters)
        elif(r_dist == 'uniform'):
            low = self.mu_radii - np.sqrt(12) * self.sigma_radii * 0.5
            high = self.mu_radii + np.sqrt(12) * self.sigma_radii * 0.5
            self.radii = np.random.uniform(low, high, self.N_emitters)
        else:
            print ("wrong distribution for radii")
            
        self.heights = self.betas * self.radii
        self.currents = np.copy(self.betas) * 0.
        self.areas = np.copy(self.currents)
        
           
    def emit(self, Ffar):
        for i in range(len(self.betas)):
            self.emitter.F = self.betas[i] * Ffar
            self.emitter.R = self.radii[i]
            self.currents[i], self.areas[i] =  self.emitter.total_current(0.5)
            
        self.current_densities = self.currents / self.areas
        return sum(self.currents)   
        
    def emit_tab(self, Ffar):        
        for i in range(len(self.betas)):
            self.tb.emitter.F = self.betas[i] * Ffar
            self.tb.emitter.R = self.radii[i]
            self.currents[i], self.areas[i] =  self.tb.total_current(0.5)
        self.current_densities = self.currents / self.areas
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
        # counts,Iedges = np.histogram(np.log(self.currents), 20, range = (np.log(1.e-9), np.log(max(self.currents))))
        
        # print counts
        # print Iedges
        # Ibins = .5*(np.exp(Iedges[1:]) + np.exp(Iedges[:-1]))
        # plt.semilogx(Ibins, counts * Ibins)
        # plt.xlabel(r"I [Amps])")
        # plt.ylabel("current contribution")
        # plt.show()
        
    
    def print_data_out(self):
        
        # for i in range(len(betas)):
            # print heights[i], radii[i], betas[i]
            
        print ("heights = %.3g +- %.3g"%(np.mean(self.heights), np.std(self.heights)))
        print ("radii = %.3g +- %.3g", np.mean(self.radii), "std<radii> =  ", np.std(self.radii))
        print ("betas = %.3g +- %.3g"%(np.mean(self.betas),np.std(self.betas)))
        print ("maxbeta = %g"%(np.max(self.betas)))
        
        


def emit (F = 5., W = 4.5, R = 5., gamma = 10., Temp = 300., verbose = False):
    """Calculate the current density and the Nottingham heating for specific set
    of input parameters"""
    this =   emission_create(F,W,R,gamma,Temp)
    if (this.ierr == 0 or np.isnan(this.Jem) or np.isnan(this.heat) ):
        if (verbose): this.print_data()
        return (this.Jem, this.heat)
    else:
        print ('Error:', this.ierr)
        this.print_data()
        return (this.Jem, this.heat)
        
def MLplot (xfn, beta = 1., W = 4.5, R = 5., gamma = 10., Temp = 300., approx = 2):
    """Take input data xfn=1/V (x axis of FN plot) and calculate the current density
    for a specific set of parameters and F = beta * V. Output in logarithmic scale """
    yfn = np.empty([len(xfn)])
    this = emission_create(5,W,R,gamma,Temp, approx=approx)
    
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
                V_appl = 5.e3, err_fact = 0.5, approx = 1, prefactor = 1.):
    """Calculate the current density and the Nottingham heating for specific set
    of input parameters, taking into account the space_charge effect"""
    F_p = min(F, 15.)
    this = emission_create(F_p,W,R,gamma,Temp, approx = approx)
    theta_old = F_p / F
    
    for j in range(50):
        this.F = F_p
        this.cur_dens()
        J_p = this.Jem * prefactor
        theta_new = theta_SC(J_p, V_appl, F_p)
        # print ("F = %f, J = %e, theta = %e" %(F_p, J_p, theta_new))
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
    
    
        




    
