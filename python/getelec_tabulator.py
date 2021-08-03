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

def lFD(E, kT):
    #assert E.size() == kT.size()
    if (isinstance(E, np.ndarray)):
        L = np.copy(E)

        highE = E > 10. * kT
        lowE = E < -10. * kT
        midE = np.logical_not(highE) * np.logical_not(lowE)

        L[highE] = np.exp(-E[highE] / kT)
        L[lowE] = - E[lowE] / kT
        L[midE]=np.log(1. + np.exp(-E[midE] / kT))
    else:
        if E > 10 * kT:
            return np.exp(-E / kT)
        elif(E < -10 * kT):
            return -E / kT
        else:
            return np.log(1. + np.exp(-E / kT))

class Tabulator():

    
    def __init__(self, Nf = 256, Nr = 128, Ngamma = 32):
        if(self.load()):
            pass
        else:
            self.Frange = [0.5, 15.]
            self.Rmin = 0.5
            self.Nf = Nf
            self.Nr = Nr
            self.Ngamma = Ngamma
            print ("tabulating G, F, R, gamma")
            self.tabulate()
            self.save()
    


    def tabulate(self):
        warnings.filterwarnings("error")
        self.Finv = np.linspace(1/self.Frange[1], 1./ self.Frange[0], self.Nf)
        self.Rinv = np.linspace(1.e-3, 1/self.Rmin, self.Nr)
        self.gaminv = np.linspace(1.e-3, .99, self.Ngamma)
        self.Gtab = np.ones([self.Ngamma, self.Nr, self.Nf, Npoly])

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
                    self.Gtab[i,j,k,:] = poly
    

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
            return True
        except(IOError):
            print ("tabulation files not found")
            return False
    
    def get_Gpoly(self, F, R, gamma):
        return intrp.interpn((self.gaminv, self.Rinv, self.Finv), self.Gtab, (1/gamma, 1./R, 1./F))[0]

    def cur_dens_metal(self, F, R, gamma, Work, kT):
        poly = self.get_Gpoly(F, R, gamma)
        integrand = lambda E: lFD(E, kT) / (1. + np.exp(np.polyval(poly, Work - E)))
        return ig.quad(integrand, -10., Work)




        
        


tab = Tabulator(64,32,8)

poly = tab.get_Gpoly(5., 5., 50.)
print(poly)
print(tab.cur_dens_metal(5.,5.,50,4.5,0.025))
