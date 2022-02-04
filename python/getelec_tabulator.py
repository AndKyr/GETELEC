
import ctypes as ct
import numpy as np
from numpy.lib.polynomial import RankWarning
from scipy.integrate.quadpack import IntegrationWarning
import scipy.optimize as opt
import scipy.integrate as ig
import os
import matplotlib.pyplot as plt
import scipy.interpolate as intrp
import json
import datetime
#import warnings

from io import StringIO 
import sys
import getelec_mod as gt

pythonpath,filename = os.path.split(os.path.realpath(__file__))
emissionpath,pythonfolder = os.path.split(pythonpath)
libpath = emissionpath + '/lib/dynamic/libgetelec.so'
#libslatecpath = emissionpath + '/lib/libslatec.so'
#ct.cdll.LoadLibrary(libslatecpath)
print(libpath)
ct.cdll.LoadLibrary(libpath)
getelec = ct.CDLL(libpath)

integrator = ct.CDLL(pythonpath + '/libintegrator.so') #use absolute path
integrator.intfun.restype = ct.c_double
integrator.intfun.argtypes = (ct.c_int, ct.c_double)
integrator.intfun_Pn.restype = ct.c_double
integrator.intfun_Pn.argtypes = (ct.c_int, ct.c_double)
# integrator.intfun_Pn.argtypes = (ct.c_int, ct.c_void_p)
integrator.Gfun.argtypes = (ct.c_int, ct.c_void_p)
integrator.Gfun.restype = ct.c_double
integrator.intfun_dbg.argtypes = (ct.c_int, ct.c_void_p)
integrator.intfun_dbg.restype = ct.c_double

Npoly = 5
NGi = 128
zs = 1.6183e-4
kBoltz = 8.6173324e-5 

m = 9.1093837015e-31 # free electron mass, same for holes
#me = 1.08*m # effective electron mass @300K, https://doi.org/10.1142/S0219749909005833
mp = 0.54*m # effective hole mass @300k, https://doi.org/10.1142/S0219749909005833
#mp = m
me = m



"""
Tabulator is class designed to reduce the complexity of GETELEC, from a 4D calculation [J(E,F,R,gammow)] to a 3D calculation [J(E,F,R)], where gammow is
either a) loaded from an existing file or b) approximated to a polynomial, which is then saved for later use. Thus reducing the computational time to calculate J and Pn.

There are 6 functions within the class Tabulator:
    1) __init__()
        This initialising function calls other functions within Tabulator to either load the polynomials describing Gammow or to calculate them 
        As inputs this function uses:
            Nf, Nr and Ngamma number of data points, given field, tip radius and gamma ranges, for which we want to evaluate Gammow
            Frange which defines the range for the field we want to apply to our emitter (in V/nm)
            Rmin is the minimum value for the tip radius (in nm)
            
    2) load()
        This function looks for files where gammow has been pre-calculated and stored for future use. Thus saving the computational time of calculating it.
        Gtab, Finv, Rinv and gaminv are loaded and ready to use
        Then "RegularGridInterpolator" is use to asign values of Gtab to the corresponding Finv, Rinv and gaminv combination.
            These relations are stored in interpolator for later use
        
    3) tabulate()
        This function is called if there are nor pre-calculated values for Finv, Rinv, gaminv and Gtab.
        It sets Finv, Rinv, gaminv for arbitrary values we give, from which the values for the field (F), tip radius (R) and gamma will be obtained
        For each value of F, R and gamma, a C function calculates Gammow (G) and the energy range (w) for which G is evaluated
            wmin is defined as the energy at the top of the potential barrier
            wmax is defined as the energy where G gets a values high enough, here G=50, so the tunnelling probability is ~0%
        Then it tries to fit G(W) to a polynomial. Such polynomial and the energy limits are stored in Gtab
                    G(W)|                                                     .
                        |                                                   .
                        |                                                .
                        |                                            .
                        |                                       .
                        |                               . 
                        |                    .               
                        | .                                  
                        |_________________|_________________________|________________W (energy)
                         w=0(Evac)        wmin (Top barrier)        wmax (G=50)
    
    4) Save()
        It saves the calculations we have made, if any, in order to be reused in future simulations
    
    5) get_Gtab()  
        Calls the interpolator (see 2)) to access the tabulated polynomials describing Gammow
"""

class Tabulator():
    def __init__(self, Nf = 256, Nr = 128, Ngamma = 32):
        self.Nf = Nf
        self.Nr = Nr
        self.Ngamma = Ngamma
        if(self.load()):
            pass
        else:
            self.Frange = [0.5, 20.]
            self.Rmin = 0.5
            print ("tabulating G, F, R, gamma")
            self.tabulate()
            self.save()

    def load(self):
        try:
            self.Gtab = np.load("tabulated/Gtable.npy")
            self.Finv = np.load("tabulated/Finv.npy")
            self.Rinv = np.load("tabulated/Rinv.npy")
            self.gaminv = np.load("tabulated/gammainv.npy")
            if(np.shape(self.Gtab) == tuple([self.Ngamma, self.Nr, self.Nf, Npoly+2])):
                self.interpolator = intrp.RegularGridInterpolator((self.gaminv, self.Rinv, self.Finv), self.Gtab)
                return True
            else:
                print("tabulation filed cannot be interpolated")
                return False
        except(IOError):
            print ("tabulation files not found")
            return False
    
    def tabulate(self):
        #warnings.filterwarnings("error")
        self.Finv = np.linspace(1/self.Frange[1], 1./ self.Frange[0], self.Nf)
        self.Rinv = np.linspace(1.e-3, 1/self.Rmin, self.Nr)
        self.gaminv = np.linspace(1.e-3, .99, self.Ngamma)
        self.Gtab = np.ones([self.Ngamma, self.Nr, self.Nf, Npoly+2])

        for i in range(self.Ngamma):
            for j in range(self.Nr):
                for k in range(self.Nf):
                    
                    F = 1/self.Finv[k]
                    R = 1/self.Rinv[j]
                    gamma = 1/self.gaminv[i]
                    
                    Wmin = np.array([1.])
                    Wmax = np.copy(Wmin)
                    G = np.zeros(NGi)
                    
                    getelec.export_gamow(ct.c_double(F), ct.c_double(R), ct.c_double(gamma), ct.c_int(NGi), \
                        ct.c_void_p(Wmin.ctypes.data), ct.c_void_p(Wmax.ctypes.data), ct.c_void_p(G.ctypes.data))
                    
                    W = np.linspace(Wmin[0], Wmax[0], NGi)
                    try:
                        poly = np.polyfit(W, G, Npoly-1)
                    except(RankWarning):
                        print("Rank Warning for F = %g, R = %g, gamma = %g"%(F,R,gamma))
                        plt.plot(W,G)
                        plt.show()
                    self.Gtab[i,j,k,:] = np.append(poly, [W[0], W[-1]])
    
    def save(self):
        try:
            os.mkdir("tabulated")
        except:
            pass
            
        np.save("tabulated/Gtable", self.Gtab)
        np.save("tabulated/Finv", self.Finv)
        np.save("tabulated/Rinv", self.Rinv)
        np.save("tabulated/gammainv", self.gaminv)
    
    def load_Gtab(self, F, R, gamma):
        outintr = self.interpolator.__call__([1/gamma, 1./R, 1./F])[0]
        return outintr
    
"""
Emitter is a class that takes as input the field, the tip radius, gamma, the electron energy and the temperate, 
along with tabulated values of Gammow from Tabulator, to calculate the current density of a metallic field emitter

There are 11 functions within the class Emitter. We could bundle them as
    A) Global functions, which are common to the two ways we have implemented to calculate J.
        These are the functions are:
        __init__()
        set()
        interpolate()
        cur_metal_dens()
            get_lims()
    B) Those proper of J calculation mode one:
            integrate_quad()
    C) Those proper of J calculation mode two:
            integrate_lin()
                lFD()
                transmission()
                    Gamow()
         
    1) __init__()
        This function initializes the the class by assigining the values from Tabulator() to itself
    
    2) set()
        It takes as inputs the values for the field, the tip radius, and the and the gamma coefficient for which we want to calculate J
    
    3) interpolate()
        With the F, R and gamma values from set(), this function calls the class Tabulator fucntion load_Gtab
        It loads the tabulated coefficients for Gammow and its energy limits, as calculated before (see class Tabulator 3))
        Then it takes the derivative of Gammow and evalutes it at the energy limits
                    G(W)|                                                    .       
                        |                                                 .      
                        |                                              .    
                        |                                           .
                        |                                       .
                        |                               . 
                        | ................  .                                           
                        |_________________|_________________________|________________W (energy)
                         w=0(Evac)        wmin (Top barrier)        wmax (G=50)

    
    4) cur_dens_metal()
        This function takes as inputs the work and the temperate in (eV)
        It calls two functions to calculate J
            get lims() (see below 5))
            integrate_lin() (see below 6)) or integrate_quad()  (see below 10))
        Then current density is the returned value
        
    5) get_lims()
        This function calculates the energy limits of the integral to calculate J
        NOTE: This is to be completed by Andreas
    
    6) integrate_lin()
        This function calculates J by inself
            NOTE: it might be redundant on future versions when a definite J calculation mode is selected
        It calls two functions to calculate the integrands
            lFD() (see below 7))
            transmission() (see below 8))
            
    7) lFD()
        This function calculates the electron supply for different regimes as a function of the electron energy E and the emitter temperature kT
        If the energy E is large compared with the thermal energy, we will have field emission 
        If the energy E is small compared with the thermal energy, we will have thermal emission
        If E and kt are similar, we will have an intermediate regime as described by the GFT theory
        
    8) transmission()
        It provides an expression for the transmission coefficient D by evaluating the Gammow coefficeint
        To evaluate Gammow, this fucntion calls the fucntion Gamow() (see below 9))
    
    9) Gamow()
        This function provides specific values for Gammow given a set of energies
        It takes the polynomial coefficients from Tabulator() and evaluate them at the different energies
        Then it returs such values to trasnmission()
        
    10) integrate_quad()
        It integrates using a function called from GETELEC to calculate J
    
    11) plot_quad()
        WARNING: I don't quite understand what is this function about
            
"""
 
class Emitter():
    def __init__(self, tabula):
        self.tabula = tabula
   
    def set(self, F, R, gamma):
        self.F = F
        self.R = R
        self.gamma = gamma
    
    def interpolate(self):
        data = self.tabula.load_Gtab(self.F, self.R, self.gamma)
        self.Wmin = data[-2]
        self.Wmax = data[-1]
        self.Gpoly = data[:Npoly]
        self.dG = np.polyder(self.Gpoly)
        self.dGmin = np.polyval(self.dG, self.Wmin)
        self.dGmax = np.polyval(self.dG, self.Wmax)
        
    def cur_dens_semi(self, Ec, Ef, Ev, kT):
        self.get_lims_semi(Ec, Ef, Ev, kT)

        Jc = self.integrate_lin_semi_conduction_band(Ec, Ef, kT)
        Jv = self.integrate_lin_semi_valence_band(Ef, Ev, kT)

        return Jc + Jv

    def get_lims_semi(self, Ec, Ef, Ev, kT):
        self.Eclow = (Ec-Ef)
        self.Echigh = max(self.Eclow, 0) + 10 * kT
        self.Evhigh = (Ev-Ef)
        self.Evlow = min(self.Evhigh, 0) - 10 / self.dGmax

    def integrate_lin_semi_conduction_band(self, Ec, Ef, kT): #calculates the integers for J from the conduction band - Eq (7) 10.1103/PhysRev.125.67
        ec = np.linspace(self.Eclow, self.Echigh, 128)
        #print(ec)
        ecm = (me/m) * ec
        ef = Ef-Ec

        Jcinteg = self.lFD(ec, kT) * (self.transmission(-Ef-ec) - ((1-(me/m)) * self.transmission(-Ef-ec+ecm)))
        #print(self.transmission(-Ef-ec))
        #print((1-(me/m)) * self.transmission(-Ef-ec+ecm)))
        return zs * kT * np.sum(Jcinteg) * (ec[1]-ec[0]) 
        
    def integrate_lin_semi_valence_band(self, Ef, Ev, kT): #calculates the integers for J from the valence band - Eq. (A1) 10.1103/PhysRev.135.A794
        ev = np.linspace(self.Evlow, self.Evhigh, 128)
        evm = (mp/m) * ev
        efv = Ef-Ev

        Evhigh_eff = Ev * (1 + mp/m)

        ev_eff = np.linspace(Ev, Evhigh_eff, 128)

        ValenceEffect = np.log(1+np.exp((Ev-efv)/2)) * np.sum(self.transmission(-Ef-ev)) * (ev_eff[1] - ev_eff[0])
        
        Jvinteg = self.lFD(ev, kT) * (self.transmission(-Ef-ev) - ((1+(mp/m)) * self.transmission(-Ef-ev-evm)))

        return zs * kT * ((np.sum(Jvinteg) * (ev[1]-ev[0])) + ValenceEffect)

    def get_Pn_semi(self, Ec, Ef, Ev, kT):
        self.get_lims_semi(Ec, Ef, Ev, kT)

        Pnc = self.Pn_semi_conduction_band(Ec, Ef, kT)
        Pnv = self.Pn_semi_valence_band(Ef, Ev, kT)
        return -Pnc + Pnv

    def Pn_semi_conduction_band(self, Ec, Ef, kT):
        ec = np.linspace(self.Eclow, self.Echigh, 128)
        
        G_heat_C = np.copy(ec)
        dif_c = ec[1]-ec[0]
        G_heat_C[0] = 0
        ecm = (me/m) * ec 
        ef = Ef-Ec

        Dc = (self.transmission(-Ef-ec) - ((1-(me/m)) * self.transmission(-Ef-ec+ecm)))

        for i in range(len(ec)-1):
            G_heat_C[i+1] = G_heat_C[i]+(dif_c*(Dc[i+1]+Dc[i])/2)
        
        heat_integ_c = (ec*G_heat_C)/(1+np.exp((ec)/kT))

        return zs * np.sum(heat_integ_c) * dif_c
    
    def Pn_semi_valence_band(self, Ef, Ev, kT):
        ev = np.linspace(self.Evlow, self.Evhigh, 128)

        G_heat_V = np.copy(ev)
        dif_v = ev[1]-ev[0]
        G_heat_V[0] = 0
        evm = (mp/m) * ev
        efv = Ef-Ev

        Evhigh_eff = Ev * (1 + mp/m)

        ev_eff = np.linspace(Ev, Evhigh_eff, 128)

        D_ValenceEffect = self.transmission(-Ef-ev) * (ev_eff[1] - ev_eff[0])

        heat_integ_ValenceEffect = (dif_v/(1+np.exp((Ev-efv)/kT)))*np.sum(D_ValenceEffect)

        Dv = (self.transmission(-Ef-ev) - ((1+(mp/m)) * self.transmission(-Ef-ev-evm)))

        for i in range(len(ev)-1):
            G_heat_V[i+1] = G_heat_V[i]+(dif_v*(Dv[i+1]+Dv[i])/2)

        heat_integ_v = (ev*G_heat_V)/(1+np.exp((ev)/kT))

        return zs * ((np.sum(heat_integ_v) * dif_v) + heat_integ_ValenceEffect)

    def cur_dens_metal(self, Work, kT, plot = False):
        self.get_lims(Work, kT)
        return self.integrate_quad(Work, kT)
        #return self.integrate_lin(Work, kT, plot)

    def get_lims(self, Work, kT):
        """finds the limits of integration"""
        self.maxbeta = np.polyval(self.dG, min(Work, self.Wmax))
        if (self.maxbeta * kT < 1.05): #field regime
            Wcenter = 0.
            self.Ehigh = 10 /(1/kT - .85*self.maxbeta)
            self.Elow = -10 / self.maxbeta
        elif (self.dGmin * kT > .95):
            Wcenter = Work - self.Wmin
            self.Ehigh = Wcenter + 10 * kT
            self.Elow = Wcenter - 10 / self.dGmin
        else:
            rootpoly = np.copy(self.dG)
            rootpoly[-1] -= 1./kT
            rts = np.roots(rootpoly)
            realroots = np.real(rts[np.nonzero(np.imag(rts) == 0)])
            Wcenter = realroots[np.nonzero(np.logical_and(realroots > self.Wmin, realroots < Work))][0]
            #print (Wcenter)
            self.Ehigh = Work - Wcenter + 10 * kT
            self.Elow = Work - Wcenter - 25 * kT
    
    def get_Pn(self, Work, kT):
        E = np.linspace(self.Elow, self.Ehigh, 128)
        D = self.transmission(Work-E)
        G_heat = np.copy(E)
        dif = (E[1]-E[0])
        G_heat[0] = 0
    
        for i in range(len(E)-1):
            G_heat[i+1] = G_heat[i]+(dif*(D[i+1]+D[i])/2)

        heat_intg = (E*G_heat)/(1+np.exp(E/kT))

        return -zs * np.sum(heat_intg) * dif
    
    def integrate_lin(self, Work, kT, plot = False):
        E = np.linspace(self.Elow, self.Ehigh, 128)
        integ = self.lFD(E, kT) * self.transmission(Work - E)

        if(plot):
            print("Whigh, Wlow = ", Work - self.Elow, Work - self.Ehigh)
            print ("Wmin, Wmax, Work = ", self.Wmin, self.Wmax, Work)
            print("maxbeta, minbeta, dGmax, beta_T: ", self.maxbeta, self.dGmin, self.dGmax, 1/ kT)
            plt.plot(E,integ)
            ax = plt.gca()
            ax2 = ax.twinx()
            ax2.plot(E,self.Gamow(Work - E), 'r-')
            ax.grid()
            plt.savefig("Jcur.png")
            # plt.show()
        return zs * kT * np.sum(integ) * (E[1] - E[0])


    def integrate_lin_Pn(self, Work, kT, plot = False):
        E = np.linspace(self.Elow, self.Ehigh, 128)
        integ = np.copy(E)
        for i in range(len(integ)):
            args = np.array([E[i], Work, kT, self.Wmin, self.Wmax, self.dGmin, self.dGmax] + list(self.Gpoly))
            N = ct.c_int(len(args))
            cargs = ct.c_void_p(args.ctypes.data)
            integ[i] = integrator.intfun_Pn(N, cargs)

        if(plot):
            print("Whigh, Wlow = ", Work - self.Elow, Work - self.Ehigh)
            print ("Wmin, Wmax, Work = ", self.Wmin, self.Wmax, Work)
            print("maxbeta, minbeta, dGmax, beta_T: ", self.maxbeta, self.dGmin, self.dGmax, 1/ kT)
            plt.plot(E,integ)
            ax = plt.gca()
            ax2 = ax.twinx()
            ax2.plot(E,self.Gamow(Work - E), 'r-')
            ax.grid()
            plt.savefig("Pnspectra.png")
            # plt.show()
        return zs * kT * np.sum(integ) * (E[1] - E[0])

    def lFD(self, E, kT):
    #assert E.size() == kT.size()
        if (isinstance(E, np.ndarray)):
            L = np.copy(E)

            highE = E > 10. * kT
            lowE = E < -10. * kT
            midE = np.logical_not(highE) * np.logical_not(lowE)

            L[highE] = np.exp(-E[highE] / kT)
            L[lowE] = - E[lowE] / kT
            L[midE] = np.log(1. + np.exp(-E[midE] / kT))
            return L
        else:
            if E > 10 * kT:
                return np.exp(-E / kT)
            elif(E < -10 * kT):
                return -E / kT
            else:
                return np.log(1. + np.exp(-E / kT))

    def transmission(self, W):
        G = self.Gamow(W)
        if (isinstance(W, np.ndarray)):
            D = np.copy(G)
            D[G < 15.] = 1 / (1 + np.exp(G[G < 15.]))
            D[G > 15] = np.exp(-G[G > 15.])
            return D
        else:
            if (G < 15.):
                return 1 / (1 + np.exp(G))
            else:
                return np.exp(-G)

    def Gamow(self, W):
    #used integrate_lin
        if (isinstance(W, np.ndarray)):
            G = np.copy(W)
            highW = W > self.Wmax
            lowW = W < self.Wmin
            G = np.polyval(self.Gpoly, W)
            G[highW] = np.polyval(self.Gpoly, self.Wmax) + self.dGmax * (W[highW] - self.Wmax)
            G[lowW] = np.polyval(self.Gpoly, self.Wmin) + self.dGmin * (W[lowW] - self.Wmin)
            return G
        else:
            if (W > self.Wmin and W < self.Wmax):
                return np.polyval(self.Gpoly, W)
            elif(W > self.Wmax):
                return np.polyval(self.Gpoly, self.Wmax) + self.dGmax * (W - self.Wmax)
            else:
                return np.polyval(self.Gpoly, self.Wmin) + self.dGmin * (W - self.Wmin)

    def integrate_quad(self, Work, kT):
        args = tuple([Work] + [kT] + [self.Wmin] + [self.Wmax] + [self.dGmin] + [self.dGmax] + list(self.Gpoly) )
        try:
            integ, abserr, info = ig.quad(integrator.intfun, self.Elow, self.Ehigh, args, full_output = 1)
            self.integ_points = info["alist"]
        except(IntegrationWarning):
            integ = 0.
        return zs * kT * integ

    def integrate_quad_Nottingham(self, Work, kT):
        args = tuple([Work] + [kT] + [self.Wmin] + [self.Wmax] + [self.dGmin] + [self.dGmax] + list(self.Gpoly) )
        try:
            integ, abserr = ig.quad(integrator.intfun_Pn, self.Elow, self.Ehigh, args, full_output = 0)
        except(IntegrationWarning):
            integ = 0.
        return -zs * kT * integ

    def plot_quad(self, Work, kT):
        E = self.integ_points
        integ = np.copy(E)
        gamow = np.copy(E)
        for i in range(len(E)):
            args = np.array([E[i], Work, kT, self.Wmin, self.Wmax, self.dGmin, self.dGmax] + list(self.Gpoly))
            N = ct.c_int(len(args))
            cargs = ct.c_void_p(args.ctypes.data)
            gamow[i] = integrator.Gfun(N, cargs)
            integ[i] = integrator.intfun_dbg(N, cargs)
        
        plt.plot(E, integ)
        plt.grid()
        ax = plt.gca()
        ax2 = ax.twinx()
        ax2.plot(E,gamow, 'r-')
        plt.show()

"""
Main programmen

Tha class Tabulator is innitiated to:
    a) Load the precalculated polynomials for Gammmow
    b) Calculate such polynomials

F, R, gamma, W and T are randomly assigned to evaluate the performance of Tabulator Vs. GETELEC
This randomization is done to prevent any biased values from the programmers when evaluating the code performance

The class Emitter is innitiated with tab as the argument to establish a connection between both classes
J is calculated from tabulator --> Ji
J is calculated from GETELEC --> JgetElem --> Jget

Then the error is calculated and displayed for comparison
"""
      
tab = Tabulator()
tab.save() #this might be redundant becuase it is already included within the __init__

Fmax = 1/tab.Finv[0]
Fmin = 1/tab.Finv[-1]
Rmax = 1/tab.Rinv[0]
Rmin = 1/tab.Rinv[-1]
gammax = 1/tab.gaminv[0]
gammin = 1/tab.gaminv[-1]

Np = 8096

Fi = np.random.rand(Np) * (Fmax - Fmin) + Fmin
Ri = np.random.rand(Np) * (Rmax - Rmin) + Rmin
gami = np.random.rand(Np) * (gammax - gammin) + gammin
Ji = np.copy(Fi)
Wi = np.random.rand(Np) * (7.5 - 2.5) + 2.5
Ti = np.random.rand(Np) * (3000 - 100) + 100
kT = Ti * kBoltz
Ji = np.copy(Fi)
Pi = np.copy(Fi)
Jget = np.copy(Ji)
Pget = np.copy(Ji)
Ji_semi = np.copy(Ji)
Eci = -np.random.rand(Np) * (4.5 - 3.5) + 3.5
Efi = -np.random.rand(Np) * (4.0 - 3.0) + 3.0
Evi = Eci + Efi 

emit = Emitter(tab)

"""
print("calculating from tabulator")
#tab_start = datetime.datetime.now()
for i in range(len(Fi)):
    emit.set(Fi[i], Ri[i], gami[i])
    emit.interpolate()
    Ji[i] = emit.cur_dens_metal(Wi[i], kT[i])
    Pi[i] = emit.get_Pn(Wi[i], kT[i])
    #Pi[i] = emit.integrate_quad_Nottingham(Wi[i], kT[i])
   #Ji_semi[i] = emit.cur_dens_semi(Eci[i], Efi[i], Evi[i], kT[i]) 
    
#tab_end = datetime.datetime.now()

print("calculating from getelec")
#get_start = datetime.datetime.now()
em = gt.emission_create(approx=2)
for i in range(len(Fi)):   
    em.F = Fi[i]
    em.W = Wi[i]
    em.Temp = Ti[i]
    em.gamma = gami[i]
    em.R = Ri[i]
    em.cur_dens()
    Jget[i] = em.Jem
    Pget[i] = em.heat
#get_end = datetime.datetime.now()


abserr = abs(Ji - Jget)
relerr = abserr / Jget
bad = np.where(np.logical_and(relerr > 0.5, abserr > 1.e-25))[0]

Pn_abserr = abs(Pi - Pget)
Pn_relerr = Pn_abserr / Pget
Pbad = np.where(np.logical_and(Pn_relerr > 0.5, Pn_abserr > 1.e-25))[0]



print("bad = ", bad)
print("rms error in J= ", np.sqrt(np.mean(relerr[abserr > 1.e-25]**2)))
print("rms error in Pn= ", np.sqrt(np.mean(Pn_relerr[Pn_abserr > 1.e-25]**2)))

#print("Tab running time", tab_end-tab_start)
#print("get running time", get_end-get_start)

# for i in bad:
#     print("Jget, Ji : ", Jget[i], Ji[i])
#     emit.set(Fi[i], Ri[i], gami[i])
#     emit.interpolate()
#     emit.get_lims(Wi[i], kT[i])
#     emit.integrate_quad(Wi[i], kT[i])
#     emit.integrate_quad_Nottingham(W[i], kT[i])

fig = plt.figure(figsize=(16,6))
plt.loglog(Ji, Jget, '.')
plt.loglog(abs(Pi), abs(Pget), '.')
plt.loglog([1.e-50, 1.], [1.e-50, 1.])
plt.grid()
plt.title("J comparison")
plt.savefig("J comparison.png")
#plt.show()

fig2 = plt.figure(figsize=(16,6))
plt.loglog(Pi, Pget, '.')
plt.loglog(abs(Pi), abs(Pget), '.')
plt.loglog([1.e-50, 1.], [1.e-50, 1.])
plt.grid()
plt.title("Pn comparison")
plt.savefig("Pn comparison.png")
#plt.show()

fig3 = plt.figure(figsize=(16,6))
plt.loglog(Ji_semi, Ji_semi, '.')
plt.loglog(abs(Ji_semi), Ji_semi, '.')
plt.loglog([1.e-50, 1.], [1.e-50, 1.])
plt.grid()
plt.title("J from semiconductors")
plt.savefig("J from semiconductors.png")"""

"""test single value of new Nottingham integration"""

W = 4.5
Temp = 300.
kT = kBoltz * Temp
emit.set(5., 10., 10.)
emit.interpolate()
emit.get_lims(W, kT)


#Pn = emit.integrate_quad_Nottingham(W, kT)
Pn_metal = emit.get_Pn(W, kT)
J_metal = emit.cur_dens_metal(W, kT)

em = gt.emission_create(approx=2) 
em.F = 5.
em.W = 4.5
em.Temp = 300.
em.gamma = 10.
em.R = 10.
em.cur_dens()
#print("Pget = ", em.heat, "Ptab = ", Pn)

Ec = -3.5
Ef = -4.5
Ev = Ec-1.1

J_semi = emit.cur_dens_semi(Ec, Ef, Ev, kT)
Pn_semi = emit.get_Pn_semi(Ec, Ef, Ev, kT)
print(J_metal, J_semi)
print(Pn_metal, Pn_semi)