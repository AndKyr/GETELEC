
from array import array
import ctypes as ct
from importlib.metadata import distribution
from pickle import TRUE
from tokenize import Number
import numpy as np
from numpy.lib.polynomial import RankWarning
from pkg_resources import Distribution, WorkingSet
from scipy.integrate import IntegrationWarning
import scipy.optimize as opt
import scipy.integrate as ig
import os
import matplotlib.pyplot as plt
import scipy.interpolate as intrp
import scipy.constants as const
import json
import datetime
#import warnings

from io import StringIO 
import sys
import getelec_mod as gt
#from interp3d import interp_3d

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

class Tabulator:
    # region field declarations
    _number_of_gamma_values: int
    _number_of_radius_values: int
    _number_of_field_values: int
    # endregion
    
    # region initialization
    def __init__(self, Nf=256, Nr=128, Ngamma=32):
        self._number_of_field_values = Nf
        self._number_of_radius_values = Nr
        self._number_of_gamma_values = Ngamma
        was_table_loaded_from_files = self._Load_Table_from_Files()
        if not was_table_loaded_from_files:
            self._range_of_field_values = [0.5, 20.]
            self._minimum_radius = 0.5
            print("tabulating G, F, R, gamma")
            self._Calculate_Table()
            self._Save_Table_to_Files()

    def _Load_Table_from_Files(self) -> bool:
        try:
            self.Gtab = np.load("tabulated/Gtable.npy")
            self.Finv = np.load("tabulated/Finv.npy")
            self.Rinv = np.load("tabulated/Rinv.npy")
            self.gaminv = np.load("tabulated/gammainv.npy")
            if (np.shape(self.Gtab) == tuple([self._number_of_gamma_values, self._number_of_radius_values, self._number_of_field_values, Npoly + 2])):
                self.interpolator = intrp.RegularGridInterpolator((self.gaminv, self.Rinv, self.Finv), self.Gtab)
                #self.interpolator = interp_3d.Interp3D(self.Gtab, self.gaminv, self.Rinv, self.Finv)
                return True
            else:
                print("tabulation filed cannot be interpolated")
                return False
        except(IOError):
            print("tabulation files not found")
            return False
    
    def _Calculate_Table(self):
        # warnings.filterwarnings("error")
        self.Finv = np.linspace(1 / self._range_of_field_values[1], 1. / self._range_of_field_values[0], self._number_of_field_values)
        self.Rinv = np.linspace(1.e-3, 1 / self._minimum_radius, self._number_of_radius_values)
        self.gaminv = np.linspace(1.e-3, .99, self._number_of_gamma_values)
        self.Gtab = np.ones([self._number_of_gamma_values, self._number_of_radius_values, self._number_of_field_values, Npoly + 2])
        
        for i in range(self._number_of_gamma_values):
            for j in range(self._number_of_radius_values):
                for k in range(self._number_of_field_values):
                    
                    F = 1 / self.Finv[k]
                    R = 1 / self.Rinv[j]
                    gamma = 1 / self.gaminv[i]
                    
                    Wmin = np.array([1.])
                    Wmax = np.copy(Wmin)
                    G = np.zeros(NGi)
                    
                    getelec.export_gamow(
                        ct.c_double(F), ct.c_double(R), ct.c_double(gamma), ct.c_int(NGi), \
                        ct.c_void_p(Wmin.ctypes.data), ct.c_void_p(Wmax.ctypes.data), ct.c_void_p(G.ctypes.data)
                    )
                    
                    W = np.linspace(Wmin[0], Wmax[0], NGi)
                    try:
                        poly = np.polyfit(W, G, Npoly - 1)
                    except(RankWarning):
                        print("Rank Warning for F = %g, R = %g, gamma = %g" % (F, R, gamma))
                        plt.plot(W, G)
                        plt.show()
                    self.Gtab[i, j, k, :] = np.append(poly, [W[0], W[-1]])
    
    def _Save_Table_to_Files(self):
        try:
            os.mkdir("tabulated")
        except:
            pass
        
        np.save("tabulated/Gtable", self.Gtab)
        np.save("tabulated/Finv", self.Finv)
        np.save("tabulated/Rinv", self.Rinv)
        np.save("tabulated/gammainv", self.gaminv)
    # endregion
    
    # region user methods
    def Gammow_Exponent_for_Parameters(self, F, R, gamma):
        # Using __call_ because of scipy reasons
        outintr = self.interpolator.__call__([1 / gamma, 1. / R, 1. / F])[0]
        return outintr
    # endregion

class Emitter:
    # region field declarations
    _energy: float
    _kT: float
    _field: float
    _radius: float
    _gamma: float
    _energy_top_barrier: float
    _energy_bottom_barrier:float
    _gammow_coefficients: array
    _gammow_derivative_coefficients: array
    _gammow_derivative_top_barrier: float
    _gammow_derivative_bottom_barrier: float
    # endregion
   
    # region initialization
    def __init__(self, tabulator: Tabulator):
        self.tabulator = tabulator
    # endregion
    
    # region user methods
    def Define_Barrier_Parameters(self, field: float, radius: float, gamma: float):
        self._field = field
        self._radius = radius
        self._gamma = gamma
        
    def Interpolate_Gammow(self):
        data = self.tabulator.Gammow_Exponent_for_Parameters(self._field, self._radius, self._gamma)
        self._energy_top_barrier = data[-2] #Wmin
        self._energy_bottom_barrier = data[-1] #Wmax
        self._gammow_coefficients = data[:Npoly] #Gpoly
        self._gammow_derivative_coefficients = np.polyder(self._gammow_coefficients) #dG
        self._gammow_derivative_top_barrier = np.polyval(self._gammow_derivative_coefficients, self._energy_top_barrier) #dGmin
        self._gammow_derivative_bottom_barrier = np.polyval(self._gammow_derivative_coefficients, self._energy_bottom_barrier) #dGmax
    
    def Transmission_Coefficient(self, energy):
        """Calculates the transmission coefficient"""
        gamow_coefficient = self.Gamow(energy)
        if (isinstance(energy, np.ndarray)):
            transmission_coefficient = np.copy(gamow_coefficient)
            transmission_coefficient[gamow_coefficient < 15.] = 1 / (1 + np.exp(gamow_coefficient[gamow_coefficient < 15.]))
            transmission_coefficient[gamow_coefficient > 15] = np.exp(-gamow_coefficient[gamow_coefficient > 15.])
            return transmission_coefficient
        else:
            if (gamow_coefficient < 15.):
                return 1 / (1 + np.exp(gamow_coefficient))
            else:
                return np.exp(-gamow_coefficient)
    
    def Gamow(self, energy):
        if (isinstance(energy, np.ndarray)):
            gammow = np.copy(energy)
            high_energy = energy > self._energy_bottom_barrier
            low_energy = energy < self._energy_top_barrier
            gammow = np.polyval(self._gammow_coefficients, energy)
            gammow[high_energy] = np.polyval(self._gammow_coefficients, self._energy_bottom_barrier) + self._gammow_derivative_bottom_barrier * (energy[high_energy] - self._energy_bottom_barrier)
            gammow[low_energy] = np.polyval(self._gammow_coefficients, self._energy_top_barrier) + self._gammow_derivative_top_barrier * (energy[low_energy] - self._energy_top_barrier)
            return gammow
        else:
            if (energy > self._energy_top_barrier and energy < self._energy_bottom_barrier):
                return np.polyval(self._gammow_coefficients, energy)
            elif(energy > self._energy_bottom_barrier):
                return np.polyval(self._gammow_coefficients, self._energy_bottom_barrier) + self._gammow_derivative_bottom_barrier * (energy - self._energy_bottom_barrier)
            else:
                return np.polyval(self._gammow_coefficients, self._energy_top_barrier) + self._gammow_derivative_top_barrier * (energy - self._energy_top_barrier) 

    def Log_Fermi_Dirac_Distribution(self, energy, kT):
    #assert E.size() == kT.size()
        if (isinstance(energy, np.ndarray)):
            distribution = np.copy(energy)

            high_energy = energy > 10. * kT
            low_energy = energy < -10. * kT
            mid_energy = np.logical_not(high_energy) * np.logical_not(low_energy)

            distribution[high_energy] = np.exp(-energy[high_energy] / kT)
            distribution[low_energy] = - energy[low_energy] / kT
            distribution[mid_energy] = np.log(1. + np.exp(-energy[mid_energy] / kT))
            return distribution
        else:
            if energy > 10 * kT:
                return np.exp(-energy / kT)
            elif(energy < -10 * kT):
                return -energy / kT
            else:
                return np.log(1. + np.exp(-energy / kT))
    
    def Fermi_Dirac_Distribution(self, energy, kT):

        distribution = 1/(1+np.exp(energy/kT))
        
        return distribution
    # endregion
   
class Metal_Emitter:
    # region field declarations
    emitter: Emitter
    kT: float
    workfunction: float
    energy: array
    _maxbeta: float
    # endregion
    
    # region initialization
    def __init__(self, tabulator: Tabulator):
        self.emitter = Emitter(tabulator)
    # endregion   
    
    # region user methods 
    def Define_Emitter_Parameters(self, workfunction: float, kT: float):
        """Defines main emitter characteristics"""
        self.workfunction = workfunction
        self.kT = kT
        self.energy = self.Integration_Limits()
        
    def Integration_Limits(self):
        """Finds the limits of integration"""
        resolution = NGi #128
        self._maxbeta = np.polyval(self.emitter._gammow_derivative_coefficients, min(self.workfunction, self.emitter._energy_bottom_barrier))
        
        if (self._maxbeta * self.kT < 1.05): #field regime
            workfunction_center = 0.
            energy_top = 10 /(1/self.kT - .85*self._maxbeta)
            energy_bottom = -10 / self._maxbeta
            return  np.linspace(energy_bottom, energy_top, resolution)
        
        elif (self.emitter._gammow_derivative_top_barrier * self.kT > .95):
            workfunction_center = self.workfunction - self.emitter._energy_top_barrier
            energy_top = workfunction_center + 10 * self.kT
            energy_bottom = workfunction_center - 10 / self.emitter._gammow_derivative_top_barrier
            return np.linspace(energy_bottom, energy_top, resolution)
        
        else:
            rootpoly = np.copy(self.emitter._gammow_derivative_coefficients)
            rootpoly[-1] -= 1./self.kT
            realroots = np.roots(rootpoly)
            workfunction_center = realroots[np.nonzero(np.logical_and(realroots > self.emitter._energy_top_barrier, realroots < self.workfunction))][0]
            #print (Wcenter)
            energy_top = self.workfunction - workfunction_center + 10 * self.kT
            energy_bottom = self.workfunction - workfunction_center - 25 * self.kT
            return np.linspace(energy_bottom, energy_top, resolution)
     
    def Current_Density(self):
        """Calculates the field emitted current density from metal surfaces"""
        args = tuple([self.workfunction] + [self.kT] + [self.emitter._energy_top_barrier] + [self.emitter._energy_bottom_barrier] + [self.emitter._gammow_derivative_top_barrier] + [self.emitter._gammow_derivative_bottom_barrier] + list(self.emitter._gammow_coefficients))
        try:
            integ, abserr, info = ig.quad(integrator.intfun, self.energy[0], self.energy[-1], args, full_output = 1)
            self.integ_points = info["alist"]
        except(IntegrationWarning):
            integ = 0.
        return zs * self.kT * integ
    
    def Energy_Distribution(self):
        """Calculates the energy distribution of field emitted electrons from metals"""
        transmission_in_energy_z = self.emitter.Transmission_Coefficient(self.workfunction - self.energy)

        cumulative_transmission = np.insert(np.cumsum((transmission_in_energy_z[1:]+transmission_in_energy_z[:-1])*(self.energy[1]-self.energy[0])/2), 0, 0)
        
        electron_number = self.emitter.Fermi_Dirac_Distribution(self.energy, self.kT) * cumulative_transmission

        return  self.energy, electron_number
  
    def Nottingham_Heat(self):
        """Calculates the Nottingham heat resulted from field emitted electrons from metals"""
        args = tuple([self.workfunction] + [self.kT] + [self.emitter._energy_top_barrier] + [self.emitter._energy_bottom_barrier] + [self.emitter._gammow_derivative_top_barrier] + [self.emitter._gammow_derivative_bottom_barrier] + list(self.emitter._gammow_coefficients))
        try:
            integ, abserr = ig.quad(integrator.intfun_Pn, self.energy[0], self.energy[-1], args, full_output = 0)
        except(IntegrationWarning):
            integ = 0.
        return -zs * self.kT * integ
    
    def Current_Density_educational(self, plot = False):
        """Calculates the field emitted current density from metals on a clear equation to code manner"""
        integ = self.emitter.Log_Fermi_Dirac_Distribution(self.energy, self.kT) * self.emitter.Transmission_Coefficient(self.workfunction - self.energy)
        
        if(plot):
            print("Whigh, Wlow = ", self.workfunction - self.energy[0], self.workfunction - self.energy[-1])
            print ("Wmin, Wmax, Work = ", self.emitter._energy_bottom_barrier, self.emitter._energy_top_barrier, self.workfunction)
            print("maxbeta, minbeta, dGmax, beta_T: ", self._maxbeta, self.emitter._gammow_derivative_bottom_barrier, self.emitter._gammow_derivative_top_barrier, 1/ self.kT)
            plt.plot(self.energy,integ)
            ax = plt.gca()
            ax2 = ax.twinx()
            ax2.plot(self.energy,self.emitter.Gamow(self.workfunction - self.energy), 'r-')
            ax.grid()
            plt.savefig("Jcur.png")
            # plt.show()
        return zs * self.kT * np.sum(integ) * (self.energy[1] - self.energy[0])
    
    def Energy_Distribution_educational(self):
        """Calculates the energy distribution of field emitted electrons from metals on a clear equation to code manner"""
        transmission_in_energy_z = self.emitter.Transmission_Coefficient(self.workfunction - self.energy)

        cumulative_transmission = np.copy(self.energy)
        cumulative_transmission[0] = 0

        for i in range(len(self.energy)-1):
            cumulative_transmission[i+1] = cumulative_transmission[i]+((self.energy[1]-self.energy[0])*(transmission_in_energy_z[i+1]+transmission_in_energy_z[i])/2)
        
        energy_distribution = self.emitter.Fermi_Dirac_Distribution(self.energy, self.kT) * cumulative_transmission

        return  self.energy, energy_distribution
    
    def Nottingham_Heat_educational(self):
        """Calculates the Nottingham heat resulted from field emitted electrons from metals on a clear equation to code manner"""
        transmission_in_energy_z = self.emitter.Transmission_Coefficient(self.workfunction - self.energy)

        cumulative_transmission = np.copy(self.energy)
        cumulative_transmission[0] = 0

        for i in range(len(self.energy)-1):
            cumulative_transmission[i+1] = cumulative_transmission[i]+((self.energy[1]-self.energy[0])*(transmission_in_energy_z[i+1]+transmission_in_energy_z[i])/2)
        
        heat_components = self.emitter.Fermi_Dirac_Distribution(self.energy, self.kT) * cumulative_transmission * self.energy
        
        return -zs * np.sum(heat_components) * (self.energy[1] - self.energy[0])
    #endregion
    
class Semiconductor_Emitter:
    # region field declarations
    emitter: Emitter
    _max_energy_conduction_band: array
    _max_energy_valence_band: array
    _energy_conduction_band: array
    _energy_valence_band: array
    _Ec: float
    _Ef: float
    _Eg: float
    _Ev: float
    _kT: float
    _Eclow: float
    _Echigh: float
    _Evhigh: float
    _Evlow: float
    _m:float #free mass of the electrons
    _me: float #electron effective mass
    _mp: float #hole effective mass
    # endregion
    
    # region initialization
    def __init__(self, tabulator: Tabulator):
        self.emitter = Emitter(tabulator)
    # endregion

    # region user methods
    def Define_Semiconductor_Emitter_Parameters(self, Ec, Ef, Eg, kT, m, me, mp):
        """Defines main emitter characteristics (in order): botton of the conduction band, fermi level, band gap, temperature, free electrons mass, effective electron mass, effective hole mass"""
        self._Ec = -Ec
        self._Ef = -Ef
        self._Eg = -Eg
        self._Ev = self._Ec+self._Eg
        self._kT = kT
        self._m  = m
        self._me = me
        self._mp = mp
        self.Integration_Limits_for_Semiconductors()
    
    def Integration_Limits_for_Semiconductors(self):
        """Finds the limits of integration"""
        resolution = NGi #128
        self._Eclow = (self._Ec-self._Ef)
        self._Echigh = max(self._Eclow, 0) + 20 * self._kT
        self._Evhigh = (self._Ev-self._Ef)
        self._Evlow = min(self._Evhigh, 0) - 20 / self.emitter._gammow_derivative_bottom_barrier
        self._energy_conduction_band = np.linspace(self._Eclow, self._Echigh, resolution)
        self._energy_valence_band = np.linspace(self._Evlow, self._Evhigh, resolution)
      
    def Current_Density_from_Semiconductors(self):
        """Calculates the field emitted current density from semiconductor surfaces"""
        current_conduction = self.Current_from_Conduction_Band()
        current_valence = self.Current_from_Valence_Band()
        current_total = current_conduction+current_valence
        
        return current_conduction, current_valence, current_total

    def Current_from_Conduction_Band(self): 
        """Calculates the contribution of the conduction band to the field emitted current density from semiconductor surfaces - Andreas' equation"""
        a = self._me/self._m
        b = 1-a 
        
        jc_integ = self.emitter.Log_Fermi_Dirac_Distribution(self._energy_conduction_band, self._kT) * (self.emitter.Transmission_Coefficient(-self._Ef-self._energy_conduction_band) - (b*self.emitter.Transmission_Coefficient(-self._Ef-a*self._energy_conduction_band)))
        
        return zs * self._kT * np.sum(jc_integ) * (self._energy_conduction_band[1]-self._energy_conduction_band[0]) 
    
    def Current_from_Valence_Band(self):
        """Calculates the contribution of the valence band to the field emitted current density from semiconductor surfaces - Andreas' equation"""
        a = self._mp/self._m
        b = 1+a
        
        jv_integ = self.emitter.Log_Fermi_Dirac_Distribution(self._energy_valence_band, self._kT) * (self.emitter.Transmission_Coefficient(-self._Ef-self._energy_valence_band) - (b*self.emitter.Transmission_Coefficient(-self._Ef-b*self._energy_valence_band+a*self._Evhigh)))
        
        return zs * self._kT * ((np.sum(jv_integ) * (self._energy_valence_band[1]-self._energy_valence_band[0])))
    
    def Distributed_Current_Density_from_Semiconductors(self):
        #conduction band
        a = self._me/self._m
        b = 1-a 
        
        c_transmission_in_z = self.emitter.Transmission_Coefficient(-self._Ef-self._energy_conduction_band) - (b*self.emitter.Transmission_Coefficient(-self._Ef-a*self._energy_conduction_band))

        c_cumulative_transmission = np.insert(np.cumsum((c_transmission_in_z[1:]+c_transmission_in_z[:-1])*(self._energy_conduction_band[1]-self._energy_conduction_band[0])/2), 0, 0)
        
        c_number_of_electrons = np.sum(self.emitter.Fermi_Dirac_Distribution(self._energy_conduction_band, self._kT) * c_cumulative_transmission )
        
        total_c = np.sum(zs * c_number_of_electrons * (self._energy_conduction_band[1]-self._energy_conduction_band[0]))
        
        #valence band
        c = self._mp/self._m
        d = 1+c
        
        v_transmission_in_z = self.emitter.Transmission_Coefficient(-self._Ef-self._energy_valence_band) - (d*self.emitter.Transmission_Coefficient(-self._Ef-d*self._energy_valence_band+c*self._Evhigh))
        
        v_cumulative_transmission = np.insert(np.cumsum((v_transmission_in_z[1:]+v_transmission_in_z[:-1])*(self._energy_valence_band[1]-self._energy_valence_band[0])/2), 0, 0)
            
        v_number_of_electrons = self.emitter.Fermi_Dirac_Distribution(self._energy_valence_band, self._kT) * v_cumulative_transmission 
        
        total_v = np.sum(zs * v_number_of_electrons * (self._energy_valence_band[1]-self._energy_valence_band[0]))
        
        total = total_c + total_v
        
        return total_c, total_v, total
    
    def Energy_Distribution_from_Semiconductors(self):
        """Calculates the energy distribution of field emitted electrons from semiconductors"""
        energy_space_c, electrons_from_c = self.Energy_Distribution_from_Conduction_Band()
        energy_space_v, electrons_from_v = self.Energy_Distribution_from_Valence_Band()

        return  energy_space_c, electrons_from_c, energy_space_v, electrons_from_v

    def Energy_Distribution_from_Conduction_Band(self):
        """Calculates the contribution of the conduction band to the energy distribution of field emitted electrons from semiconductors """
        a = self._me/self._m
        b = 1-a 
        
        transmission_in_z = self.emitter.Transmission_Coefficient(-self._Ef-self._energy_conduction_band) - (b*self.emitter.Transmission_Coefficient(-self._Ef-a*self._energy_conduction_band))

        cumulative_transmission = np.insert(np.cumsum((transmission_in_z[1:]+transmission_in_z[:-1])*(self._energy_conduction_band[1]-self._energy_conduction_band[0])/2), 0, 0)
        
        number_of_electrons = self.emitter.Fermi_Dirac_Distribution(self._energy_conduction_band, self._kT) * cumulative_transmission 

        return self._energy_conduction_band , number_of_electrons
    
    def Energy_Distribution_from_Valence_Band(self): #NEED REVISION
        """Calculates the contribution of the valence band to the energy distribution of field emitted electrons from semiconductors"""
        a = self._mp/self._m
        b = 1+a
        
        transmission_in_z = self.emitter.Transmission_Coefficient(-self._Ef-self._energy_valence_band) - (b*self.emitter.Transmission_Coefficient(-self._Ef-b*self._energy_valence_band+a*self._Evhigh))
        
        cumulative_transmission = np.insert(np.cumsum((transmission_in_z[1:]+transmission_in_z[:-1])*(self._energy_valence_band[1]-self._energy_valence_band[0])/2), 0, 0)
            
        number_of_electrons = self.emitter.Fermi_Dirac_Distribution(self._energy_valence_band, self._kT) * cumulative_transmission 
      
        return self._energy_valence_band, number_of_electrons
    
    def Nottingham_Heat_from_Semiconductors(self):
        """Calculates the Nottingham heat resulted from field emitted electrons from semicoductors"""
        nottingham_conduction = self.Nottingham_Heat_from_Conduction_Band()
        nottigham_valence = self.Nottingham_Heat_from_Valence_Band()
        nottigham_replacement = self.Nottingham_Heat_fron_Replacement_Electrons()
        nottingham_total = nottingham_conduction + nottigham_valence + nottigham_replacement   
        
        return nottingham_conduction, nottigham_valence, nottingham_total
    
    def Nottingham_Heat_from_Conduction_Band(self): #NEED MODIFICATIONS
        """Calculates the contribution of the conduction band to the Nottingham heat resulted from field emitted electrons from semicoductors"""
        """the fuction will be replaced by the full code from the function, if we agree the function works well"""
        energy_distribution = self.Energy_Distribution_from_Conduction_Band()*self._energy_conduction_band
        
        return -zs* self._kT * np.sum(energy_distribution)
    
    def Nottingham_Heat_from_Valence_Band(self): #NEED MOFICATIONS
        """Calculates the contribution of the valence band to the Nottingham heat resulted from field emitted electrons from semicoductors"""
        """the fuction will be replaced by the full code from the function, if we agree the function works well"""
        energy_distribution = self.Energy_Distribution_from_Valence_Band()*self._energy_valence_band
        
        return -zs * self._kT * np.sum(energy_distribution)
    
    def Nottingham_Heat_fron_Replacement_Electrons(self):
        """To be calculated from COMSOL"""
        return 0
    
    def Energy_Distribution_from_Conduction_Band_educational(self):
        """Calculates the contribution of the conduction band to the energy distribution of field emitted electrons from semiconductors on a clear equation to code manner"""
        a = self._me/self._m
        b = 1-a 
        
        transmission_in_z = self.emitter.Transmission_Coefficient(-self._Ef-self._energy_conduction_band) - (b*self.emitter.Transmission_Coefficient(-self._Ef-a*self._energy_conduction_band))

        cumulative_transmission = np.copy(self._energy_conduction_band)
        cumulative_transmission[0] = 0

        for i in range(len(self._energy_conduction_band)-1):
            cumulative_transmission[i+1] = cumulative_transmission[i] + ((self._energy_conduction_band[1]-self._energy_conduction_band[0])*(transmission_in_z[i+1]+transmission_in_z[i])/2)
        
        energy_distribution = self.emitter.Fermi_Dirac_Distribution(self._energy_conduction_band, self._kT) * cumulative_transmission 

        return self._energy_conduction_band , energy_distribution
    # endregion

# region One data point calculation routine
    # This routine calculates the current density from semiconductors (two methods) and metals, as well as the plotting of the energy distributions
"""     
Npoly = 5
NGi = 512
zs = 1.6183e-4
kBoltz = 8.6173324e-5 

tab = Tabulator()

Workfunction = 4.5
Ec = 3
Ef = 4.5
Eg = 0.7
m = 9.1093837015e-31 
me = 1.08*m # effective electron mass @300K, https://doi.org/10.1142/S0219749909005833
mp = 0.54*m # effective hole mass @300k, https://doi.org/10.1142/S0219749909005833
Temp = 300.
kT = kBoltz * Temp

metal_emitter = Metal_Emitter(tab)

metal_emitter.emitter.Define_Barrier_Parameters(5., 10., 10.)
metal_emitter.emitter.Interpolate_Gammow()

metal_emitter.Define_Emitter_Parameters(Workfunction, kT)

j_metal = metal_emitter.Current_Density()
pn_metal = metal_emitter.Nottingham_Heat()
energy_space_metal, distribution_metal = metal_emitter.Energy_Distribution()


semiconductor_emitter = Semiconductor_Emitter(tab)

semiconductor_emitter.emitter.Define_Barrier_Parameters(5., 10., 10.)
semiconductor_emitter.emitter.Interpolate_Gammow()

semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(Ec, Ef, Eg, kT, m, me, mp)

j_c, j_v, j_total = semiconductor_emitter.Current_Density_from_Semiconductors()

pn_c, pn_v, pn_total = semiconductor_emitter.Nottingham_Heat_from_Semiconductors()

energy_c, distribution_c, energy_v, distribution_v = semiconductor_emitter.Energy_Distribution_from_Semiconductors()

j_c2, j_v2, j_total2 = semiconductor_emitter.Distributed_Current_Density_from_Semiconductors()

c_abs_error = abs(j_c2-j_c)
c_rel_error = c_abs_error/j_c

v_abs_error = abs(j_v2-j_v)
v_rel_error = v_abs_error/j_v

total_abs_error = abs(j_total2-j_total)
total_rel_error = total_abs_error/j_total

print("current_c", j_c, j_c2, "relative error:", c_rel_error)
print("current_v", j_v, j_v2, "relative error:", v_rel_error)
print("current_semi", j_total, j_total2, "relative error:", total_rel_error)
print("current_metal", j_metal)


fig = plt.figure(figsize=(16,6))
#plt.plot(energy_space_metal, distribution_metal/max(distribution_metal))
plt.plot(energy_c, distribution_c)#/max(distribution_c))
plt.plot(energy_v, distribution_v)#/max(distribution_v))
plt.grid("True")
plt.title("Energy distributions")
plt.savefig("Energy distributions.png")
"""
# endregion

# region Multiple data point calculation routine - metal
"""
Npoly = 5
NGi = 128
zs = 1.6183e-4
kBoltz = 8.6173324e-5 

tab = Tabulator()

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
Wi = np.random.rand(Np) * (7.5 - 2.5) + 2.5
Ti = np.random.rand(Np) * (3000 - 100) + 100
kT = Ti * kBoltz
Ji = np.copy(Fi)
Pi = np.copy(Fi)
Jget = np.copy(Ji)
Pget = np.copy(Ji)

metal_emitter = Metal_Emitter(tab)

print("calculating from tabulator")
#tab_start = datetime.datetime.now()
for i in range(len(Fi)):
    metal_emitter.emitter.Define_Barrier_Parameters(Fi[i], Ri[i], gami[i])
    metal_emitter.emitter.Interpolate_Gammow()
    metal_emitter.Define_Emitter_Parameters(Wi[i], kT[i])
    Ji[i] = metal_emitter.Current_Density()
    Pi[i] = metal_emitter.Nottingham_Heat()
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


#for i in bad:
#    print("Jget, Ji : ", Jget[i], Ji[i])
#    emit.set(Fi[i], Ri[i], gami[i])
#    emit.interpolate()
#    emit.get_lims(Wi[i], kT[i])
#    emit.integrate_quad(Wi[i], kT[i])
#    emit.integrate_quad_Nottingham(W[i], kT[i])

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
"""
#endregion

# region Multiple data point calculation routine - semiconductor

Npoly = 5
NGi = 128
zs = 1.6183e-4
kBoltz = 8.6173324e-5 

tab = Tabulator()

Ef = 4.5
Eg = 1.1
Temp = 300.
m = 9.1093837015e-31 
me = 1.08*m # effective electron mass @300K, https://doi.org/10.1142/S0219749909005833
mp = 0.54*m # effective hole mass @300k, https://doi.org/10.1142/S0219749909005833
kT = kBoltz * Temp

semiconductor_emitter = Semiconductor_Emitter(tab)

semiconductor_emitter.emitter.Define_Barrier_Parameters(5., 10., 10.)
semiconductor_emitter.emitter.Interpolate_Gammow()

metal_emitter = Metal_Emitter(tab)

metal_emitter.emitter.Define_Barrier_Parameters(5., 10., 10.)
metal_emitter.emitter.Interpolate_Gammow()

resolution = 100

j_c = np.zeros(resolution) 
j_v = np.zeros(resolution) 
j_total = np.zeros(resolution) 
j_c2 = np.zeros(resolution) 
j_v2 = np.zeros(resolution) 
j_total2 = np.zeros(resolution) 
pn_c = np.zeros(resolution) 
pn_v = np.zeros(resolution) 
pn_total = np.zeros(resolution) 
Bottom_Ec = np.zeros(resolution) 
j_metal = np.zeros(resolution) 
pn_metal = np.zeros(resolution)
c_rel_error = np.zeros(resolution)
v_rel_error = np.zeros(resolution)
total_rel_error = np.zeros(resolution)

for i in range(resolution):
    Ec = i*0.05+2

    semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(Ec, Ef, Eg, kT, m, me, mp)

    metal_emitter.Define_Emitter_Parameters(Ec, kT)


    j_c[i], j_v[i], j_total[i] = semiconductor_emitter.Current_Density_from_Semiconductors()
    
    j_c2[i], j_v2[i], j_total2[i] = semiconductor_emitter.Distributed_Current_Density_from_Semiconductors()
    
    pn_c[i], pn_v[i], pn_total[i] = semiconductor_emitter.Nottingham_Heat_from_Semiconductors()
    
    energy_c, distribution_c, energy_v, distribution_v = semiconductor_emitter.Energy_Distribution_from_Semiconductors()
    
    
    j_metal[i] = metal_emitter.Current_Density()
    
    pn_metal[i] = metal_emitter.Nottingham_Heat()
    
    energy_space_metal, distribution_metal = metal_emitter.Energy_Distribution()
    
    
    Bottom_Ec[i] = Ec
    
    
    c_abs_error = abs(j_c2[i]-j_c[i])
    c_rel_error[i] = c_abs_error/j_c[i]

    v_abs_error = abs(j_v2[i]-j_v[i])
    v_rel_error[i] = v_abs_error/j_v[i]

    total_abs_error = abs(j_total2[i]-j_total[i])
    total_rel_error[i] = total_abs_error/j_total[i]
    
    
    title = "Normalised energy distrution, Ec = " + str(Ec)
    
    save_png = "Normalised energy distrution, Ec = " + str(Ec) + ".png"
    
    #fig = plt.figure(figsize=(16,6))
    #plt.plot(energy_space_metal, distribution_metal/max(distribution_metal))
    #plt.plot(energy_c, distribution_c/max(distribution_c))
    #plt.plot(energy_v, distribution_v/max(distribution_v))
    #plt.grid("True")
    #plt.title(title)
    #plt.savefig(save_png)

c_mean = np.mean(c_rel_error)
c_std = np.std(c_rel_error)
c_max = max(c_rel_error)

v_mean = np.mean(v_rel_error)
v_std = np.std(v_rel_error)
v_max = max(v_rel_error)


total_mean = np.mean(total_rel_error)
total_std = np.std(total_rel_error)
total_max = max(total_rel_error)


fig3 = plt.figure(figsize=(16,6))
plt.semilogy(Bottom_Ec, j_c)
plt.semilogy(Bottom_Ec, j_v)
plt.semilogy(Bottom_Ec, j_total)
plt.semilogy(Bottom_Ec, j_metal)
plt.grid()
plt.title("Current from semiconductor")
plt.savefig("Current from semiconductor.png")

fig4 = plt.figure(figsize=(16,6))
plt.plot(Bottom_Ec, pn_c)
plt.plot(Bottom_Ec, pn_v)
plt.plot(Bottom_Ec, pn_total)
plt.plot(Bottom_Ec, pn_metal)
plt.yscale("symlog", linscale=-1)
plt.grid('True')
plt.title("Nottinghah Heat from semiconductor")
plt.savefig("Nottinghah Heat from semiconductor.png")

file = open("Field_Emission_Data.txt", "w+")

file.write("Simulation temperature = %f\r\n\r\n" % Temp)

for i in range(len(Bottom_Ec)):
    file.write("\r\nSemiconductor Ec = %d" % (Bottom_Ec[i]))
    file.write(" Metal Ef = %d\r\n" % (Bottom_Ec[i]))
    file.write("    Current from semiconductor %e" % (j_total[i]))
    file.write("    Current from Metal = %e\r\n" % (j_metal[i]))
    file.write("        Semiconductor current relative error = %e\r\n" % (total_rel_error[i]))
    file.write("    Heat from semiconductor = %e" % (pn_total[i]))
    file.write("    Heat from metal = %e\r\n" % (pn_metal[i]))
    
file.write("\r\nMax relative error = %e\r\n" % (total_max))
  
file.close()

# endregion