
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
import getelec_mod as gt

pythonpath,filename = os.path.split(os.path.realpath(__file__))
emissionpath,pythonfolder = os.path.split(pythonpath)
libpath = emissionpath + '/lib/dynamic/libgetelec.so'
#libslatecpath = emissionpath + '/lib/libslatec.so'
#ct.cdll.LoadLibrary(libslatecpath)
print(libpath)
ct.cdll.LoadLibrary(libpath)
getelec = ct.CDLL(libpath)

Npoly = 5
NGi = 128
zs = 1.6183e-4

def get_gamow_line(F, R, gamma, Npoints = NGi):
    Wmin = np.array([1.])
    Wmax = np.copy(Wmin)
    Gamow = np.zeros(Npoints)
    getelec.export_gamow(ct.c_double(F), ct.c_double(R), ct.c_double(gamma), ct.c_int(Npoints), \
        ct.c_void_p(Wmin.ctypes.data), ct.c_void_p(Wmax.ctypes.data), ct.c_void_p(Gamow.ctypes.data))
    W = np.linspace(Wmin[0], Wmax[0], Npoints)
    return W, Gamow

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

    def tabulate(self):
        warnings.filterwarnings("error")
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
                    W,G = get_gamow_line(F, R , gamma)
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
        
    def load(self):
        try:
            self.Gtab = np.load("tabulated/Gtable.npy")
            self.Finv = np.load("tabulated/Finv.npy")
            self.Rinv = np.load("tabulated/Rinv.npy")
            self.gaminv = np.load("tabulated/gammainv.npy")
            if(np.shape(self.Gtab) == tuple([self.Ngamma, self.Nr, self.Nf, Npoly+2])):
                return True
            else:
                return False
        except(IOError):
            print ("tabulation files not found")
            return False
    
    def get_Gpoly(self, F, R, gamma):
        return intrp.interpn((self.gaminv, self.Rinv, self.Finv), self.Gtab, (1/gamma, 1./R, 1./F))[0]

class Emitter():
    def __init__(self, tabula):
        self.tabula = tabula

    def set(self, F, R, gamma):
        self.F = F
        self.R = R
        self.gamma = gamma
    
    def interpolate(self):
        data = self.tabula.get_Gpoly(self.F, self.R, self.gamma)
        self.Wmin = data[-2]
        self.Wmax = data[-1]
        self.Gpoly = data[:Npoly]
        self.dG = np.polyder(self.Gpoly)

    def lFD(self, E, kT):
    #assert E.size() == kT.size()
        if (isinstance(E, np.ndarray)):
            L = np.copy(E)

            highE = E > 10. * kT
            lowE = E < -10. * kT
            midE = np.logical_not(highE) * np.logical_not(lowE)

            L[highE] = np.exp(-E[highE] / kT)
            L[lowE] = - E[lowE] / kT
            L[midE]=np.log(1. + np.exp(-E[midE] / kT))
            return L
        else:
            if E > 10 * kT:
                return np.exp(-E / kT)
            elif(E < -10 * kT):
                return -E / kT
            else:
                return np.log(1. + np.exp(-E / kT))

    def Gamow(self, W):

        if (isinstance(W, np.ndarray)):
            G = np.copy(W)
            highW = W > self.Wmax
            lowW = W < self.Wmin
            G = np.polyval(self.Gpoly, W)
            if (len(highW + lowW) > 0):
                G[highW] = np.polyval(self.Gpoly, self.Wmax) + np.polyval(self.dG, self.Wmax) * (W[highW] - self.Wmax)
                G[lowW] = np.polyval(self.Gpoly, self.Wmin) + np.polyval(self.dG, self.Wmin) * (W[lowW] - self.Wmin)
            return G
        else:
            if (W > self.Wmin and W < self.Wmax):
                return np.polyval(self.Gpoly, W)
            elif(W > self.Wmax):
                return np.polyval(self.Gpoly, self.Wmax) + np.polyval(self.dG, self.Wmax) * (W - self.Wmax)
            else:
                return np.polyval(self.Gpoly, self.Wmin) + np.polyval(self.dG, self.Wmin) * (W - self.Wmin)

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

    def cur_dens_metal(self, Work, kT, plot = False):
        self.get_lims(Work, kT)
        return self.integrate_lin(Work, kT, plot)

    def get_lims(self, Work, kT):

        self.maxbeta = np.polyval(self.dG, min(Work, self.Wmax))
        self.dG_Wmin = np.polyval(self.dG, self.Wmin)
        if (self.maxbeta * kT < 1.05): #field regime
            Wcenter = 0.
            self.Ehigh = 10 /(1/kT - .85*self.maxbeta)
            self.Elow = -10 / self.maxbeta
        elif (self.dG_Wmin * kT > .95):
            Wcenter = Work - self.Wmin
            self.Ehigh = Wcenter + 10 * kT
            self.Elow = Wcenter - 10 / self.dG_Wmin
        else:
            rootpoly = np.copy(self.dG)
            rootpoly[-1] -= 1./kT
            Wcenter = np.roots(rootpoly)
            Wcenter = np.real(Wcenter[np.nonzero(np.logical_and(Wcenter > self.Wmin, Wcenter < Work, np.imag(Wcenter) == 0.))])[0]
            #print (Wcenter)
            self.Ehigh = Work - Wcenter + 20 * kT
            self.Elow = Work - Wcenter - 10 / kT

    def integrate_lin(self, Work, kT, plot = False):
        E = np.linspace(self.Elow, self.Ehigh, 128)
        integ = self.lFD(E, kT) * self.transmission(Work - E)

        if(plot):
            print("poly = ", self.Gpoly, "Elow, Ehigh = ", self.Elow, self.Ehigh)
            print("maxbeta, minbeta, beta_T: ", self.maxbeta, self.dG_Wmin,1/ kT)
            plt.plot(E,integ)
            ax = plt.gca()
            ax2 = ax.twinx()
            ax2.plot(E,self.calculate_Gamow(Work - E), 'r-')
            ax.grid()
            plt.show()
        return zs * kT * np.sum(integ) * (E[1] - E[0])

    def integrate_quad(self, Work, kT):
        args = tuple([Work] + [kT] + + [self.Wmin] + [self.Wmax] + list(self.Gpoly) + list(self.dG) )
        return zs * kT * ig.quad(intfun, self.Elow, self.Ehigh, args)


kBoltz = 8.6173324e-5       
        

tab = Tabulator()
tab.save()

Fmax = 1/tab.Finv[0]
Fmin = 1/tab.Finv[-1]
Rmax = 1/tab.Rinv[0]
Rmin = 1/tab.Rinv[-1]
gammax = 1/tab.gaminv[0]
gammin = 1/tab.gaminv[-1]

Np = 4096

Fi = np.random.rand(Np) * (Fmax - Fmin) + Fmin
Ri = np.random.rand(Np) * (Rmax - Rmin) + Rmin
gami = np.random.rand(Np) * (gammax - gammin) + gammin
Ji = np.copy(Fi)
Wi = np.random.rand(Np) * (7.5 - 2.5) + 2.5
Ti = np.random.rand(Np) * (3000 - 100) + 200
kT = Ti * kBoltz
Jget = np.copy(Ji)

em = Emitter(tab)

print("calculating from tabulator")
for i in range(len(Fi)):
    em.set(Fi[i], Ri[i], gami[i])
    em.interpolate()
    Ji[i] = em.cur_dens_metal(Wi[i], kT[i])

print("calculating from getelec")

em = gt.emission_create(approx=2)
for i in range(len(Fi)):   
    em.F = Fi[i]
    em.W = Wi[i]
    em.Temp = Ti[i]
    em.gamma = gami[i]
    em.R = Ri[i]
    em.cur_dens()
    Jget[i] = em.Jem

relerr = abs(Jget/Ji - 1.) 
abserr = abs(Ji - Jget)

bad = np.where(np.logical_and(relerr > 0.3, abserr > 1.e-25))[0]

print("bad = ", bad)
print("rms error = ", np.sqrt(np.mean(relerr[abserr > 1.e-25])))
for i in bad:
    print("Jget, Ji : ", Jget[i], Ji[i])
#     tab.cur_dens_metal(Fi[i], Ri[i], gami[i], Wi[i], kT[i], plot = True)

# plt.loglog(Ji, Jget, '.')
# plt.loglog([1.e-50, 1.], [1.e-50, 1.])
# plt.grid()
# plt.show()