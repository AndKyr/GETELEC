
from array import array
import ctypes as ct
from importlib.metadata import distribution
from tokenize import Number
import numpy as np
from numpy.lib.polynomial import RankWarning
from pkg_resources import Distribution, WorkingSet
from scipy.integrate.quadpack import IntegrationWarning
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
from interp3d import interp_3d

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
me = 1.08*m # effective electron mass @300K, https://doi.org/10.1142/S0219749909005833
mp = 0.54*m # effective hole mass @300k, https://doi.org/10.1142/S0219749909005833
#mp = m
#me = m

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
            rts = np.roots(rootpoly)
            realroots = np.real(rts[np.nonzero(np.imag(rts) == 0)])
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
        
        energy_distribution = self.emitter.Fermi_Dirac_Distribution(self.energy, self.kT) * cumulative_transmission

        return  self.energy, energy_distribution
  
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
    # endregion
    
    # region initialization
    def __init__(self, tabulator: Tabulator):
        self.emitter = Emitter(tabulator)
    # endregion

    # region user methods
    def Define_Semiconductor_Emitter_Parameters(self, Ec, Ef, Eg, kT):
        """Defines main emitter characteristics"""
        self._Ec = -Ec
        self._Ef = -Ef
        self._Eg = -Eg
        self._Ev = self._Ec+self._Eg
        self._kT = kT
        self.Integration_Limits_for_Semiconductors()
    
    def Integration_Limits_for_Semiconductors(self):
        """Finds the limits of integration"""
        resolution = NGi #128
        self._Eclow = (self._Ec-self._Ef)
        self._Echigh = max(self._Eclow, 0) + 20 * self._kT
        self._Evhigh = (self._Ev-self._Ef)
        self._Evlow = min(self._Evhigh, 0) - 20 / self.emitter._gammow_derivative_bottom_barrier
        self._energy_conduction_band = np.linspace(self._Eclow, self._Echigh, resolution)
        self._max_energy_conduction_band = (me/m) * self._energy_conduction_band
        self._energy_valence_band = np.linspace(self._Evlow, self._Evhigh, resolution)
        self._max_energy_valence_band = (mp/m) * self._energy_valence_band
           
    def Current_Density_from_Semiconductors(self):
        """Calculates the field emitted current density from semiconductor surfaces"""
        current_conduction = self.Current_from_Conduction_Band()
        current_valence = self.Current_from_Valence_Band()
        current_total = current_conduction+current_valence
        
        return current_conduction, current_valence, current_total

    def Current_from_Conduction_Band(self): 
        """Calculates the contribution of the conduction band to the field emitted current density from semiconductor surfaces - Eq (7) 10.1103/PhysRev.125.67"""
        jc_integ = self.emitter.Log_Fermi_Dirac_Distribution(self._energy_conduction_band, self._kT) * (self.emitter.Transmission_Coefficient(-self._Ef-self._energy_conduction_band) - ((1-(me/m)) * self.emitter.Transmission_Coefficient(-self._Ef-self._energy_conduction_band-self._max_energy_conduction_band)))
        
        return zs * self._kT * np.sum(jc_integ) * (self._energy_conduction_band[1]-self._energy_conduction_band[0]) 
    
    def Current_from_Valence_Band(self):
        """Calculates the contribution of the valence band to the field emitted current density from semiconductor surfaces - Eq. (A1) 10.1103/PhysRev.135.A794"""
        max_energy_v_effect = self._energy_valence_band[-1]*(1+(mp/m))
        
        energy_v_effect = np.linspace(max_energy_v_effect, self._energy_valence_band[-1], 128)

        valence_effect = self.emitter.Log_Fermi_Dirac_Distribution(self._energy_valence_band[-1],self._kT) * np.sum(self.emitter.Transmission_Coefficient(-self._Ef-energy_v_effect)) * (energy_v_effect[1] - energy_v_effect[0])
        
        jv_integ = self.emitter.Log_Fermi_Dirac_Distribution(self._energy_valence_band, self._kT) * (self.emitter.Transmission_Coefficient(-self._Ef-self._energy_valence_band) + ((1+(mp/m)) * self.emitter.Transmission_Coefficient(-self._Ef-self._energy_valence_band-self._max_energy_valence_band)))
        
        return zs * self._kT * ((np.sum(jv_integ) * (self._energy_valence_band[1]-self._energy_valence_band[0])) + np.sum(valence_effect))
  
    def Energy_Distribution_from_Semiconductors(self):
        """Calculates the energy distribution of field emitted electrons from semiconductors"""
        energy_space_c, energy_distribution_c = self.Energy_Distribution_from_Conduction_Band()
        energy_space_v, energy_distribution_v = self.Energy_Distribution_from_Valence_Band()
        return  energy_space_c, energy_distribution_c, energy_space_v, energy_distribution_v

    def Energy_Distribution_from_Conduction_Band(self):
        """Calculates the contribution of the conduction band to the energy distribution of field emitted electrons from semiconductors """
        transmission_in_z = self.emitter.Transmission_Coefficient(-self._Ef-self._energy_conduction_band) - ((1-(me/m)) * self.emitter.Transmission_Coefficient(-self._Ef-self._energy_conduction_band-self._max_energy_conduction_band))

        cumulative_transmission = np.insert(np.cumsum((transmission_in_z[1:]+transmission_in_z[:-1])*(self._energy_conduction_band[1]-self._energy_conduction_band[0])/2), 0, 0)
        
        energy_distribution = self.emitter.Fermi_Dirac_Distribution(self._energy_conduction_band, self._kT) * cumulative_transmission 

        return self._energy_conduction_band , energy_distribution
    
    def Energy_Distribution_from_Valence_Band(self): #NEED MODIFICATIONS #FOR LOOP
        """Calculates the contribution of the valence band to the energy distribution of field emitted electrons from semiconductors"""
        a = mp/m
        b = 1+a
        
        transmission_in_z = self.emitter.Transmission_Coefficient(-self._Ef-self._energy_valence_band) - (1+(mp/m)) * self.emitter.Transmission_Coefficient(-self._Ef-b*self._energy_valence_band+a*self._Evhigh)  
        #Transmission_in_Ex = np.abs(Transmission_in_Ex)
        
        cumulative_transmission = np.insert(np.cumsum((transmission_in_z[1:]+transmission_in_z[:-1])*(self._energy_valence_band[1]-self._energy_valence_band[0])/2), 0, 0)
            
        energy_distribution = self.emitter.Fermi_Dirac_Distribution(self._energy_valence_band, self._kT) * cumulative_transmission 
      
        return self._energy_valence_band, energy_distribution
    
    def Nottingham_Heat_from_Semiconductors(self):
        """Calculates the Nottingham heat resulted from field emitted electrons from semicoductors"""
        nottingham_conduction = self.Nottingham_Heat_from_Conduction_Band()
        nottigham_valence = self.Nottingham_Heat_from_Valence_Band()
        nottingham_total = nottingham_conduction + nottigham_valence   
        
        return nottingham_conduction, nottigham_valence, nottingham_total
    
    def Nottingham_Heat_from_Conduction_Band(self): #NEED MODIFICATIONS
        """Calculates the contribution of the conduction band to the Nottingham heat resulted from field emitted electrons from semicoductors"""
        energy_distribution = self.Energy_Distribution_from_Conduction_Band()
        
        return -zs* self._kT * np.sum(energy_distribution)
    
    def Nottingham_Heat_from_Valence_Band(self): #NEED MOFICATIONS
        """Calculates the contribution of the valence band to the Nottingham heat resulted from field emitted electrons from semicoductors"""
        energy_distribution = self.Energy_Distribution_from_Valence_Band()
        
        return -zs * self._kT * np.sum(energy_distribution)
    
    def Energy_Distribution_from_Conduction_Band_educational(self):
        """Calculates the contribution of the conduction band to the energy distribution of field emitted electrons from semiconductors on a clear equation to code manner"""
        transmission_in_ez = self.emitter.Transmission_Coefficient(-self._Ef-self._energy_conduction_band) - ((1-(me/m)) * self.emitter.Transmission_Coefficient(-self._Ef-self._energy_conduction_band-self._max_energy_conduction_band))

        cumulative_transmission = np.copy(self._energy_conduction_band)
        cumulative_transmission[0] = 0

        for i in range(len(self._energy_conduction_band)-1):
            cumulative_transmission[i+1] = cumulative_transmission[i] + ((self._energy_conduction_band[1]-self._energy_conduction_band[0])*(transmission_in_ez[i+1]+transmission_in_ez[i])/2)
        
        energy_distribution = self.emitter.Fermi_Dirac_Distribution(self._energy_conduction_band, self._kT) * cumulative_transmission 

        return self._energy_conduction_band , energy_distribution
    # endregion

# region One data point calculation routine
tab = Tabulator()

work = 4.5
Ec = 6
Ef = 4.5
Eg = 1.1
Temp = 300.
kT = kBoltz * Temp

metal_emitter = Metal_Emitter(tab)

metal_emitter.emitter.Define_Barrier_Parameters(5., 10., 10.)
metal_emitter.emitter.Interpolate_Gammow()

metal_emitter.Define_Emitter_Parameters(work, kT)

j_metal = metal_emitter.Current_Density()
pn_metal = metal_emitter.Nottingham_Heat()
energy_space_metal, distribution_metal = metal_emitter.Energy_Distribution()


semiconductor_emitter = Semiconductor_Emitter(tab)

semiconductor_emitter.emitter.Define_Barrier_Parameters(5., 10., 10.)
semiconductor_emitter.emitter.Interpolate_Gammow()

semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(Ec, Ef, Eg, kT)

j_c, j_v, j_total = semiconductor_emitter.Current_Density_from_Semiconductors()

pn_c, pn_v, pn_total = semiconductor_emitter.Nottingham_Heat_from_Semiconductors()

energy_c, distribution_c, energy_v, distribution_v = semiconductor_emitter.Energy_Distribution_from_Semiconductors()

print("current_semi", j_total)
print("current_metal", j_metal)


fig = plt.figure(figsize=(16,6))
plt.plot(energy_space_metal, distribution_metal/max(distribution_metal))
plt.plot(energy_c, distribution_c/max(distribution_c))
plt.plot(energy_v, distribution_v/max(distribution_v))
plt.grid("True")
plt.title("Energy distributions")
plt.savefig("Energy distributions.png")
# endregion


# region Multiple data point calculation routine
"""tab = Tabulator()

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
Energy_C = np.copy(Ji)
Electrons_C = np.copy(Ji)
Energy_Ca = np.copy(Ji)
Electrons_Ca = np.copy(Ji)
Energy_V = np.copy(Ji)
Electrons_V = np.copy(Ji)
Energy_M = np.copy(Ji)
Electrons_M = np.copy(Ji)
Eci = -np.random.rand(Np) * (4.5 - 3.5) + 3.5
Efi = -np.random.rand(Np) * (4.0 - 3.0) + 3.0
Evi = Eci + Efi

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
#plt.show()"""

"""test single value of new Nottingham integration"""
"""
W = 4.5
Temp = 300.
kT = kBoltz * Temp
emit.set(5., 10., 10.)
emit.interpolate()
emit.get_lims(W, kT)


#Pn = emit.integrate_quad_Nottingham(W, kT)
Pn_metal = emit.integrate_quad_Nottingham(W, kT)
J_metal = emit.cur_dens_metal(W, kT)

em = gt.emission_create(approx=2) 
em.F = 5.
em.W = 4.5
em.Temp = 300.
em.gamma = 10.
em.R = 10.
em.cur_dens()
#print("Pget = ", em.heat, "Ptab = ", Pn)

Ec = -6
Ef = -4.5
Ev = Ec-0.5

J_c, J_v, J_semi = emit.Current_from_Semiconductors(Ec, Ef, Ev, kT)
Pn_c, Pn_v, Pn_semi = emit.Nottinghah_Heat_from_Semiconductors(Ec, Ef, Ev, kT)
Energy_C, Electrons_C, Energy_V, Electrons_V = emit.Energy_Distribution_from_Semiconductors(Ef, kT)
Energy_M, Electrons_M = emit.Energy_Distribution_from_Metals(W, kT)
#Energy_Ca, Electrons_Ca = emit.Energy_Distribution_from_Conduction_Band2(Ef, kT, 4, 4.5E7)
#print(J_metal, J_semi)
#print(Pn_metal, Pn_semi)
#print(Energy_C, Electrons_C)

fig = plt.figure(figsize=(16,6))
plt.plot(Energy_C, Electrons_C/max(Electrons_C))
#plt.plot(Energy_Ca, Electrons_Ca/max(Electrons_Ca))
plt.plot(Energy_V, Electrons_V/max(Electrons_V))
#plt.plot(Energy_M, Electrons_M/max(Electrons_M))
plt.grid("True")
plt.title("Energy distributions")
plt.savefig("Energy distributions.png")

n = 400

Current_Conduction_Band = np.zeros(n)
Current_Valence_Band = np.zeros(n)
Total_Semiconductor_Current = np.zeros(n)

Total_Metal_Current = np.zeros(n)

NottinghamHeat_Conduction_Band = np.zeros(n)
NottinghamHeat_Valence_Band = np.zeros(n)
Total_Semiconductor_NottinghamHeat = np.zeros(n)

Energy_Conduction_Band =  np.ones(n)
Electrons_Conduction_Band = np.ones(n)

ConductionBand = np.zeros(n)

for i in range(n):
    Ec = -i*0.01-2
    Ef = -4.5
    Ev = Ec-1.1

    Current_Conduction_Band[i], Current_Valence_Band[i],Total_Semiconductor_Current[i] = emit.Current_from_Semiconductors(Ec, Ef, Ev, kT)
    
    NottinghamHeat_Conduction_Band[i], NottinghamHeat_Valence_Band[i], Total_Semiconductor_NottinghamHeat[i] = emit.Nottinghah_Heat_from_Semiconductors(Ec, Ef, Ev, kT)
    
    Total_Metal_Current[i] = emit.cur_dens_metal(-Ec, kT)
    
    ConductionBand[i] = Ec

fig3 = plt.figure(figsize=(16,6))
plt.semilogy(ConductionBand, Current_Conduction_Band)
plt.semilogy(ConductionBand, Current_Valence_Band)
plt.semilogy(ConductionBand, Total_Semiconductor_Current)
plt.semilogy(ConductionBand, Total_Metal_Current)
plt.grid()
plt.title("Current from semiconductor")
plt.savefig("Current from semiconductor.png")

fig4 = plt.figure(figsize=(16,6))
plt.plot(ConductionBand, NottinghamHeat_Conduction_Band)
plt.plot(ConductionBand-1.1, NottinghamHeat_Valence_Band)
plt.plot(ConductionBand, Total_Semiconductor_NottinghamHeat)
plt.yscale("symlog", linscale=-1)
plt.grid('True')
#plt.title("Nottinghah Heat from semiconductor")
#plt.savefig("Nottinghah Heat from semiconductor.png")"""
# endregion