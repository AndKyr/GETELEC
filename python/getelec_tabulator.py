
from array import array
import ctypes as ct
#from importlib.metadata import distribution
#from pickle import TRUE
#from tkinter.ttk import LabeledScale
#from tokenize import Number
#from turtle import color
import numpy as np
from numpy.lib.polynomial import RankWarning
#from pkg_resources import Distribution, WorkingSet
from scipy.integrate import IntegrationWarning
#import scipy.optimize as opt
import scipy.integrate as ig
import os
import matplotlib.pyplot as plt
import scipy.interpolate as intrp
#import scipy.constants as const
#import matplotlib as mb
#import json
import time
#import warnings

#from io import StringIO 
#import sys
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
#integrator.intfun_Pn.argtypes = (ct.c_int, ct.c_void_p)
integrator.Gfun.argtypes = (ct.c_int, ct.c_void_p)
integrator.Gfun.restype = ct.c_double
integrator.intfun_dbg.argtypes = (ct.c_int, ct.c_void_p)
integrator.intfun_dbg.restype = ct.c_double

Npoly = 5
NGi = 512
zs = 1.6183e-4

class Tabulator:
    """
    Tabulator is class designed to reduce the computational effort of GETELEC, thus enabling to allocate computationa resources into solving other problems
    like the calculation of a semiconductor's band structure. This functionality is based on the fact that the potential barrier of an emitter is independed 
    from the emitter's material. So this barrier can be calcualted once and to be later called when needed.
    
    Tabulator first looks through the barrier maps (function: _Load_Table_from_File) to see if there is a pre-calculated barrier. If there is one, it passes 
    it to the rest of the model for further computation. If there is not, it uses the maps from whose it would infer one using tabulation 
    (function: Gammow_Exponent_for_Parameters). In the case that the maps do not exist, Tabulator calculate (function: _Calculate_Table) and save 
    (function: _Save_Table_to_Files) these maps - a process than can take several hours. 

    There are 5 functions within the class Tabulator. Four of which are private to the class (see the _ before the function name) and one which it public (no
    _ before the name) and will be inheritated by class Emitter.
    
    1) __init__()
        Initialises Tabulator by 1) loading the "resolution", 2) loading the maps and 3) if maps not found, calling the functions to calculate and sabe the maps
            
    2) _Load_Table_from_Files()
        Looks for the files where the precaculate barriers are stored. Then it uses interpolation methods to make the most accurate barrier for the given 
        input (electric field, tip radius and gamma exponent). Gtab is stores the polinomial that gives its shape to the barrier.
        
    3) _Calculate_Table()
        This function is called if the barrier maps are not found. For each field, tip radius and gamma combitation, Gtab is calculated by polynomial fitting 
        to the data points extracted from getelec.export_gamow
        For each value of the field, the radius and gamma, the Gammow coeffient (G) and the energy range (W), for such G, are evaluated
            Wmin is defined as the energy at the top of the potential barrier. Point from which G goes to 0 and tranmission cofficient equals 1
            Wmax is defined as the energy where G gets a values high enough, here G=50, so the tunnelling probability is effectively 0
        Then it tries to fit G(W) to a polynomial. Such polynomial and the energy limits are stored in Gtab
        
        Here a depiction of how the gammow coefficient changes with energy, as one goes from the vacuum level (E = 0 eV) to the bottom of the 
        energy spectrum, and associated tramission coefficient.
                    G(W)|                                                     .
                        |                                                    .
                        |                                                 .
                        |                                             .
                        |                                        .
                        |                                .
                        |                     .
                        |.
                        |_________________|_________________________|________________W (energy)
                         W=0(Vacuum)      Wmin (Top barrier)        Wmax (G=50)
    
    
               D(W) = 1 |.................                                                      
                        |                 .                                  
                        |                  .                                
                        |                    .                            
                        |                        .                       
                        |                             .               
                        |                                    .                
                        |                                           .                 
                      0 |_________________|_________________________|._______________W (energy)
                         W=0(Vacuum)      Wmin (Top barrier)        Wmax (G=50)
    
    
    
    
    4) _Save_Table_to_Files()
        It saves the barrier maps, if any, in order to be reused in future simulations
    
    5) Gammow_Exponent_for_Parameters()
        It is the public function that takes the three describing barrier parameters (field, radius, gamma) and calls the interpolator 
        (_Load_Table_from_Files()) to calculate the polinomial that describes the barrier, as well as the energy limis (Wmin and Wmax)
        of it
"""

    # region field declarations
    _number_of_gamma_values: int
    _number_of_radius_values: int
    _number_of_field_values: int
    # endregion
    
    # region initialization
    def __init__(self, Nf=256, Nr=128, Ngamma=32):
        """Initialises Tabulator by 1) loading the "resolution", 2) loading the maps and 3) if maps not found, calling the functions to calculate and sabe the maps

        Args:
            Nf (int, optional): Resolution for the electric field. Defaults to 256.
            Nr (int, optional): Resolution for the tip radius. Defaults to 128.
            Ngamma (int, optional): Resolution for gamma. Defaults to 32.
        """
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
        """Looks for the files where the precaculate barriers are stored. Then it uses interpolation methods to make the most accurate barrier for the given 
        input (electric field, tip radius and gamma exponent). Gtab is stores the polinomial that gives its shape to the barrier.
        """
        try:
            self.Gtab = np.load("tabulated/Gtable.npy")
            #self.Gtab1 = np.load("tabulated/Gtable1.npy")
            #self.Gtab2 = np.load("tabulated/Gtable2.npy")
            #self.Gtab3 = np.load("tabulated/Gtable3.npy")
            #self.Gtab4 = np.load("tabulated/Gtable4.npy")
            #self.Gtab5 = np.load("tabulated/Gtable5.npy")
            #self.Gtab6 = np.load("tabulated/Gtable6.npy")
            #self.Gtab7 = np.load("tabulated/Gtable7.npy")
            self.Finv = np.load("tabulated/Finv.npy")
            self.Rinv = np.load("tabulated/Rinv.npy")
            self.gaminv = np.load("tabulated/gammainv.npy")
            #if (np.shape(self.Gtab1) == tuple([self._number_of_gamma_values, self._number_of_radius_values, self._number_of_field_values])):
            if (np.shape(self.Gtab) == tuple([self._number_of_gamma_values, self._number_of_radius_values, self._number_of_field_values, Npoly + 2])):
                self.interpolator = intrp.RegularGridInterpolator((self.gaminv, self.Rinv, self.Finv), self.Gtab)
                #self.interpolator1 = interp_3d.Interp3D(self.Gtab1, self.gaminv, self.Rinv, self.Finv)
                #self.interpolator2 = interp_3d.Interp3D(self.Gtab2, self.gaminv, self.Rinv, self.Finv)
                #self.interpolator3 = interp_3d.Interp3D(self.Gtab3, self.gaminv, self.Rinv, self.Finv)
                #self.interpolator4 = interp_3d.Interp3D(self.Gtab4, self.gaminv, self.Rinv, self.Finv)
                #self.interpolator5 = interp_3d.Interp3D(self.Gtab5, self.gaminv, self.Rinv, self.Finv)
                #self.interpolator6 = interp_3d.Interp3D(self.Gtab6, self.gaminv, self.Rinv, self.Finv)
                #self.interpolator7 = interp_3d.Interp3D(self.Gtab7, self.gaminv, self.Rinv, self.Finv)
                return True
            else:
                print("tabulation filed cannot be interpolated")
                return False
        except(IOError):
            print("tabulation files not found")
            return False
        
    def _Calculate_Table(self):
        """Looks for the files where the precaculate barriers are stored. Then it uses interpolation methods to make the most accurate barrier for the given 
        input (electric field, tip radius and gamma exponent). Gtab is stores the polinomial that gives its shape to the barrier.
        """
        # warnings.filterwarnings("error")
        self.Finv = np.linspace(1 / self._range_of_field_values[1], 1. / self._range_of_field_values[0], self._number_of_field_values)
        self.Rinv = np.linspace(1.e-3, 1 / self._minimum_radius, self._number_of_radius_values)
        self.gaminv = np.linspace(1.e-3, .99, self._number_of_gamma_values)
        self.Gtab = np.ones([self._number_of_gamma_values, self._number_of_radius_values, self._number_of_field_values, Npoly + 2])
        #self.Gtab1 = np.ones([self._number_of_gamma_values, self._number_of_radius_values, self._number_of_field_values])
        #self.Gtab2 = np.ones([self._number_of_gamma_values, self._number_of_radius_values, self._number_of_field_values])
        #self.Gtab3 = np.ones([self._number_of_gamma_values, self._number_of_radius_values, self._number_of_field_values])
        #self.Gtab4 = np.ones([self._number_of_gamma_values, self._number_of_radius_values, self._number_of_field_values])
        #self.Gtab5 = np.ones([self._number_of_gamma_values, self._number_of_radius_values, self._number_of_field_values])
        #self.Gtab6 = np.ones([self._number_of_gamma_values, self._number_of_radius_values, self._number_of_field_values])
        #self.Gtab7 = np.ones([self._number_of_gamma_values, self._number_of_radius_values, self._number_of_field_values])

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
                    #self.Gtab1[i, j, k] = (poly[0])
                    #self.Gtab2[i, j, k] = (poly[1])
                    #self.Gtab3[i, j, k] = (poly[2])
                    #self.Gtab4[i, j, k] = (poly[3])
                    #self.Gtab5[i, j, k] = (poly[4])
                    #self.Gtab6[i, j, k] = (W[0])
                    #self.Gtab7[i, j, k] = (W[-1])
    
    def _Save_Table_to_Files(self):
        """It saves the barrier maps, if any, in order to be reused in future simulations
        """
        try:
            os.mkdir("tabulated")
        except:
            pass
        
        np.save("tabulated/Gtable", self.Gtab)
        #np.save("tabulated/Gtable1", self.Gtab1)
        #np.save("tabulated/Gtable2", self.Gtab2)
        #np.save("tabulated/Gtable3", self.Gtab3)
        #np.save("tabulated/Gtable4", self.Gtab4)
        #np.save("tabulated/Gtable5", self.Gtab5)
        #np.save("tabulated/Gtable6", self.Gtab6)
        #np.save("tabulated/Gtable7", self.Gtab7)
        np.save("tabulated/Finv", self.Finv)
        np.save("tabulated/Rinv", self.Rinv)
        np.save("tabulated/gammainv", self.gaminv)
    # endregion
    
    # region user methods
    def Gammow_Exponent_for_Parameters(self, F: float, R: float, gamma: float):
        """It is the public function that takes the three describing barrier parameters (field, radius, gamma) and calls the interpolator 
        (_Load_Table_from_Files()) to calculate the polinomial that describes the barrier, as well as the energy limis (Wmin and Wmax)
        of it

        Args:
            F (float): Value of the electric field
            R (float): Value of the emitter tip radius
            gamma (float): Value of gamma

        Returns:
            array: Contains the polynomial which describes the barrier, as well as the top (transmission = 1) and bottom (transmission is efectively 0) 
            limits of this barrier
        """
        # Using __call_ because of scipy reasons
        try:
            outintr = self.interpolator.__call__([1 / gamma, 1. / R, 1. / F], method="linear")[0]
            #outintr1 = self.interpolator1((1./gamma, 1./R, 1./F))
            #outintr2 = self.interpolator2((1./gamma, 1./R, 1./F))
            #outintr3 = self.interpolator3((1./gamma, 1./R, 1./F))
            #outintr4 = self.interpolator4((1./gamma, 1./R, 1./F))
            #outintr5 = self.interpolator5((1./gamma, 1./R, 1./F))
            #outintr6 = self.interpolator6((1./gamma, 1./R, 1./F))
            #outintr7 = self.interpolator7((1./gamma, 1./R, 1./F))
            #outintr = np.array([outintr1,outintr2,outintr3,outintr4,outintr5,outintr6,outintr7])

        except:
            # This except is here in case the field, from comosl, is very low to prevent errors. This outintr results in a current equal to 0 A/cm2
            print("WARNING: Field is too high/low for GETELEC. Review model/results")
            outintr = np.array([1, 1., 500, 500, 0, 2., 2.0])
            
        return outintr
    # endregion

class Emitter:
    """
    This class gathers those attributes that are common to all emitters regardless of the material they are made of, and inheritage from tabulator
    the potential barrier. Thus making an emitter. Later Emitter will be composed with the classes Metal_Emitter and Semiconductor_Emitter to complete
    the model for the emitter
    
    This class 7 classes, of which one is private and the rest public (called later by Metal_Emitter or Semiconductor_Emitter):
    1) __init__()
        Initialises the class by inheritating Tabulator and setting it up as a private field to be used by Emitter
        
    2) Define_Barrier_Parameters()
        This function will be called by the user to define the main parameters of the barrier they want to work with. These values will be fed to Tabulator
        from which the total barrier shape and magnitude will be tabulated
        
    3) Interpolate_Gammow()
        Here a function from Tabulator is called in order to interpolate the potential barrier. Then the top and bottom limits, as well as the polynomial
        describing the barrier's shape, are extracted. The polinomial coefficients are derivated and evaluated at the top and bottom of the barrier to analyse
        the barrier slope at such points - which will be required to later calculate the integration limits for the currrent density and heat intgrals
    
    4) Tranmission_Coefficient()
        Calculates the probability of an electron tunelling (being transmitted) through the potential barrier by means of the Kemble formula whitin the JWKB approximation
        
    5) Gamow()
        Evaluates the gammow coefficients in order to provide the exponential for the Kemble formula
    
    6) Log_Fermi_Dirac_Distribuition()
        Returs the natural logarith of the Fermi Diract distribution
        
    7) Fermi_Dirac_Distribution()
        Returns the Fermi Dirac distribution
    """
    
    
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
        """Initialises the class by inheritating Tabulator and setting it up as a private field to be used by Emitter

        Args:
            tabulator (Tabulator): Object that contains the attributes to interpolate the potential barrier given the field, radius and gamma as inputs
        """
        self._tabulator = tabulator
    # endregion
    
    # region user methods
    def Define_Barrier_Parameters(self, field: float, radius: float, gamma: float):
        """This function will be called by the user to define the main parameters of the barrier they want to work with. These values will be fed to Tabulator
        from which the total barrier shape and magnitude will be tabulated

        Args:
            field (float): Electric field (V/A)
            radius (float): Emitter tip raidus (nm)
            gamma (float): gamma coefficient (number)
        """
        
        self._field = field
        self._radius = radius
        self._gamma = gamma
        
    def Interpolate_Gammow(self):
        """Here a function from Tabulator is called in order to interpolate the potential barrier. Then the top and bottom limits, as well as the polynomial
        describing the barrier's shape, are extracted. The polinomial coefficients are derivated and evaluated at the top and bottom of the barrier to analyse
        the barrier slope at such points - which will be required to later calculate the integration limits for the currrent density and heat intgrals
        """
        
        data = self._tabulator.Gammow_Exponent_for_Parameters(self._field, self._radius, self._gamma)
        self._energy_top_barrier = data[-2] #Wmin
        self._energy_bottom_barrier = data[-1] #Wmax
        self._gammow_coefficients = data[:Npoly] #Gpoly
        self._gammow_derivative_coefficients = np.polyder(self._gammow_coefficients) #dG
        self._gammow_derivative_top_barrier = np.polyval(self._gammow_derivative_coefficients, self._energy_top_barrier) #dGmin
        self._gammow_derivative_bottom_barrier = np.polyval(self._gammow_derivative_coefficients, self._energy_bottom_barrier) #dGmax
    
    def Transmission_Coefficient(self, energy):
        """Calculates the probability of an electron tunelling (being transmitted) through the potential barrier by means of the Kemble formula whitin the JWKB approximation

        Args:
            energy (array): Energy, perpendicular to the emitting surface, values (eV) for which the transmission cofficients are calculated (eV)

        Returns:
            array: Probability of an electron with a certain kinetic energy perpendicular to the emitting surface to tunnerl through the potential barrier (number)
        """
        
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
        """Evaluates the gammow coefficients in order to provide the exponential for the Kemble formula

        Args:
            energy (array): Energy, perpendicular to the emitting surface, values (eV) for which the transmission cofficients are calculated (eV)

        Returns:
            array: Kemble's formula exponent value for a given energy (number)
        """
        
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
        """Returs the natural logarith of the Fermi Diract distribution.
        The energy is divided in 3 terms to prevent the computer from overfloading with large exponents 

        Args:
            energy (array): Energy, perpendicular to the emitting surface, values (eV) for which the transmission cofficients are calculated (eV)
            kT (float): Energy from temperature (eV)

        Returns:
            array: Natural logarith of the Fermi Diract distribution (number)
        """
        
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
        """Returns the Fermi Dirac distribution
        
        Args:
            energy (array): Energy, perpendicular to the emitting surface, values (eV) for which the transmission cofficients are calculated (eV)
            kT (float): Energy from temperature (eV)

        Returns:
            array: Fermi Diract distribution (number)
        """
        
        distribution = 1/(1+np.exp(energy/kT))
        
        return distribution
    # endregion
   
class Metal_Emitter:
    """This class will be composed along with Emitter(Tabulator) in order to create a metal field emitter. Then different functions will be executed to 
    calculate the density, nottingham heat and energy density of the electrons being emitted.
    
    1) __init__()
        Initialises the function by adding the atributes of a "tabulated" barrier emitter to metals
        
    2) Define_Emitter_Parameters()
        Takes the two relevant material properties (workfucntion and temperature) and makes them available for the rest of the class
        
    3) _Integration_Limits()
        Finds the limits of integration
    
    4) Current_Density()
        Calculates the field emitted current density from metal surfaces
     
    5) Energy_Distribution()
        Calculates the energy distribution of field emitted electrons from metals
       
    6) Nottingham_Heat()
        Calculates the Nottingham heat resulted from field emitted electrons from metals
        
    7) Current_Density_educational()
        Calculates the field emitted current density from metals on a clear equation to code manner
        
    8) Energy_Distribution_educational()
        Calculates the energy distribution of field emitted electrons from metals on a clear equation to code manner
        
    9) Nottingham_Heat_educational()
        Calculates the Nottingham heat resulted from field emitted electrons from metals on a clear equation to code manner
    """
    
    # region field declarations
    emitter: Emitter
    kT: float
    workfunction: float
    energy: array
    _maxbeta: float
    # endregion
    
    # region initialization
    def __init__(self, tabulator: Tabulator):
        """Initialises the function by adding the atributes of a "tabulated" barrier emitter to metals

        Args:
            tabulator (Tabulator): Object that contains the attributes of a x material emitter (potential barrier, transmission coefficiens and electron supply functions)
        """
        
        self.emitter = Emitter(tabulator)
    # endregion   
    
    # region user methods 
    def Define_Emitter_Parameters(self, workfunction: float, kT: float):
        """Defines main emitter characteristics

        Args:
            workfunction (float): Average minimum energy required to put an electron, whitin a material, out in vacuum (some few nm away from the material surface) (eV)
            kT (float): Energy from temperature (eV)
        """
        
        self.workfunction = workfunction
        self.kT = kT
        self.energy = self._Integration_Limits()
        
    def _Integration_Limits(self):
        """Finds the limits of integration. ANDREAS PLEASE MAKE A COMMENT WHY THE 3 DIFFERENT INSTANCES AND WHY THE VALUES
        
        Returns:
            array: Energy space for which the electron emission and Notigham heat are going to be evaluated in (eV)
        """
        
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
        """Calculates the field emitted current density from metal surfaces

        Returns:
            float: Emitted current density being emitted from an infinitesimal flat surface area (electrons/area*time)
        """
        
        args = tuple([self.workfunction] + [self.kT] + [self.emitter._energy_top_barrier] + [self.emitter._energy_bottom_barrier] + [self.emitter._gammow_derivative_top_barrier] + [self.emitter._gammow_derivative_bottom_barrier] + list(self.emitter._gammow_coefficients))
        try:
            integ, abserr, info = ig.quad(integrator.intfun, self.energy[0], self.energy[-1], args, full_output = 1)
            self.integ_points = info["alist"]
        except(IntegrationWarning):
            integ = 0.
        return zs * self.kT * integ
    
    def Energy_Distribution(self):
        """Calculates the energy distribution of field emitted electrons from metals

        Returns:
            energy (array): Energy space for which the electron emission has been evalated (eV)
            electron_number (array): Number of electrons being emitter with a certain energy (number)
        """
        
        transmission_in_energy_z = self.emitter.Transmission_Coefficient(self.workfunction - self.energy)

        cumulative_transmission = np.insert(np.cumsum((transmission_in_energy_z[1:]+transmission_in_energy_z[:-1])*(self.energy[1]-self.energy[0])/2), 0, 0)
        
        electron_number = self.emitter.Fermi_Dirac_Distribution(self.energy, self.kT) * cumulative_transmission

        return  self.energy, electron_number
  
    def Nottingham_Heat(self):
        """Calculates the Nottingham heat resulted from field emitted electrons from metals

        Returns:
            float: Resulted Nottigham from the emission of electrons from an infinitesinal area (J)
        """
        
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
    """This class will be composed along with Emitter(Tabulator) in order to create a semiconductor field emitter. Then different functions will be executed to 
    calculate the density, nottingham heat and energy density of the electrons being emitted
    
    1) __init__()
        Initialises the function by adding the atributes of a "tabulated" barrier emitter to metals
    
    2) Define_Semiconductor_Emitter_Parameters()
        Takes the material properties (workfucntion and temperature) and makes them available for the rest of the class
        
    3) _Integration_Limits_for_Semiconductors()
        Finds the limits of integration
    
    4) Current_Density_from_Semicondductors()
        Returns the field emitted current density from semiconducting surfaces, by calling two other functions
    
    4) Current_from_Conduction_Band()
        Calculates the conduction band component of the emitted current
        
    5) Current_from_Valence_Band()
        Calculates the valence band component of the emitter current
     
    6) Distributed_Current_Density_from_Semiconductors()
        Calculates the current density being emitter from a semiconducting flat surface by calculating the emitted electron
        energy distrubution and integrating over that distribution
    
    7) Energy_Distribution_from_Semiconductors()
        Returns the energy distribution of the electrons being emitted from the conduction and valence bands
        
    8) Energy_Distribution_from_Conduction_Band()
        Calculates teh energy distribution of those electrons emitted from the conduction band
        
    9) Energy_Distribution_from_Valence_Band()
        Calculates teh energy distribution of those electrons emitted from the valence band
       
    10) Nottingham_Heat_from_Semiconductors()
        Returns the Nottigham heat resulted from electrons being emitter fromt the conduction and valence band
        
    11) Nottingham_Heat_from_Conduction_Band()
        Calculates the Nottingham heat resulted from field emitted electrons from the conduction band
    
    12) Nottingham_Heat_from_Valence_Band()
        Calculates the Nottingham heat resulted from field emitted electrons from the valence band
    
    13) Nottingham_Heat_from_Replacement_Electrons()
        UNDER CONSTRUCTION
        It will return the energy of those electrons that replace the emitted ones
    
    14) Energy_Distribution_from_Conduction_Band_educational()
        Calculates the contribution of the conduction band to the energy distribution of field emitted electrons from semiconductors on a clear equation to code manner
    """
    
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
        """Initialises the function by adding the atributes of a "tabulated" barrier emitter to metals

        Args:
            tabulator (Tabulator): Object that contains the attributes of a x material emitter (potential barrier, transmission coefficiens and electron supply functions)
        """
        
        self.emitter = Emitter(tabulator)
    # endregion

    # region user methods
    def Define_Semiconductor_Emitter_Parameters(self, Ec: float, Ef: float, Eg: float, kT: float, m: float, me:float, mp: float):
        """Defines the main emitter parameters:,, band gap, temperature, free electrons mass, effective electron mass, effective hole mass

        Args:
            Ec (float): botton of the conduction band (eV)
            Ef (float): fermi level (eV)
            Eg (float): band gap (eV)
            kT (float): temperature (eV)
            m (float): free electron mass (kg)
            me (float): effective electron mass (number)
            mp (float): effective hole mass (number)
        """
        
        self._Ec = -Ec
        self._Ef = -Ef
        self._Eg = -Eg
        self._Ev = self._Ec+self._Eg
        self._kT = kT
        self._m  = m
        self._me = me
        self._mp = mp
        self._Integration_Limits_for_Semiconductors()
    
    def _Integration_Limits_for_Semiconductors(self):
        """Finds the limits of integration"""
        resolution = NGi #128
        self._Eclow = (self._Ec-self._Ef)
        self._Echigh = max(self._Eclow, 0) + 20 * self._kT
        self._Evhigh = (self._Ev-self._Ef)
        self._Evlow = min(self._Evhigh, 0) - 20 / self.emitter._gammow_derivative_bottom_barrier
        self._energy_conduction_band = np.linspace(self._Eclow, self._Echigh, resolution)
        self._energy_valence_band = np.linspace(self._Evlow, self._Evhigh, resolution)
      
    def Current_Density_from_Semiconductors(self):
        """Returns the field emitted current density from semiconducting surfaces, by calling two other functions

        Returns:
            current_conduction (float): current density emitted from the conduction band (electrons/area*time)
            current_valence (float): current density emitted from the valence band (electrons/area*time)
            current_total (float): total current emitted from a semiconducting infinitisimal and flat surface (sum of the conduction and valence contributions) (electrons/area*time)
        """
        
        current_conduction = self.Current_from_Conduction_Band()
        current_valence = self.Current_from_Valence_Band()
        current_total = current_conduction+current_valence
        
        return current_conduction, current_valence, current_total

    def Current_from_Conduction_Band(self): 
        """ Calculates the conduction band component of the emitted current - Stratton's 1962 equation

        Returns:
            float: Emitted current density from the conduction band (electrons/area*time)
        """
        
        a = self._me/self._m
        b = 1-a 
        
        jc_integ = self.emitter.Log_Fermi_Dirac_Distribution(self._energy_conduction_band, self._kT) * (self.emitter.Transmission_Coefficient(-self._Ef-self._energy_conduction_band) - (b*self.emitter.Transmission_Coefficient(-self._Ef-a*self._energy_conduction_band)))
        
        return zs * self._kT * np.sum(jc_integ) * (self._energy_conduction_band[1]-self._energy_conduction_band[0]) 
    
    def Current_from_Valence_Band(self):
        """Calculates the conduction band component of the emitted current - Andreas' equation

        Returns:
            float: Emitted current density from the valence band (electrons/area*time)
        """
        
        a = self._mp/self._m
        b = 1+a
        
        jv_integ = self.emitter.Log_Fermi_Dirac_Distribution(self._energy_valence_band, self._kT) * (self.emitter.Transmission_Coefficient(-self._Ef-self._energy_valence_band) - (b*self.emitter.Transmission_Coefficient(-self._Ef-b*self._energy_valence_band+a*self._Evhigh)))
        
        return zs * self._kT * ((np.sum(jv_integ) * (self._energy_valence_band[1]-self._energy_valence_band[0])))
    
    def Distributed_Current_Density_from_Semiconductors(self):
        """Calculates the current density being emitter from a semiconducting flat surface by calculating the emitted electron
        energy distrubution and integrating over that distribution
        
        This function was implemented for validation purposes and it is not optimised!

        Returns:
            float: Emitted current density from semiconductor emitter with infinitesimal area (electrons/area*time)
        """
        
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
        """Returns the energy distribution of the electrons being emitted from the conduction and valence bands

        Returns:
            energy_space_c (array): Energy values in the conduction band from which electrons are emitted (eV)
            electrons_from_c (array): Number of emitted electrons from the conduction band for each energy value (number)
            energy_space_v (array): Energy values in the valence band from which electrons are emitted (eV)
            electrons_from_v (array): Number of emitted electrons from the valence band for each energy value (number)
        """
        
        energy_space_c, electrons_from_c = self.Energy_Distribution_from_Conduction_Band()
        energy_space_v, electrons_from_v = self.Energy_Distribution_from_Valence_Band()

        return  energy_space_c, electrons_from_c, energy_space_v, electrons_from_v

    def Energy_Distribution_from_Conduction_Band(self):
        """Calculates teh energy distribution of those electrons emitted from the conduction band

        Returns:
            _energy_conduction_band (array): Energy values in the conduction band from which electrons are emitted (eV)
            number_of_electrons (array): Number of emitted electrons from the conduction band for each energy value (number)
        """
        
        a = self._me/self._m
        b = 1-a 
        
        transmission_in_z = self.emitter.Transmission_Coefficient(-self._Ef-self._energy_conduction_band) - (b*self.emitter.Transmission_Coefficient(-self._Ef-a*self._energy_conduction_band))

        cumulative_transmission = np.insert(np.cumsum((transmission_in_z[1:]+transmission_in_z[:-1])*(self._energy_conduction_band[1]-self._energy_conduction_band[0])/2), 0, 0)
        
        number_of_electrons = self.emitter.Fermi_Dirac_Distribution(self._energy_conduction_band, self._kT) * cumulative_transmission 

        return self._energy_conduction_band , number_of_electrons
    
    def Energy_Distribution_from_Valence_Band(self): #NEED REVISION
        """Calculates teh energy distribution of those electrons emitted from the valence band

        Returns:
            _energy_conduction_band (array): Energy values in the valence band from which electrons are emitted (eV)
            number_of_electrons (array): Number of emitted electrons from the valence band for each energy value (number)
        """
        
        a = self._mp/self._m
        b = 1+a
        
        transmission_in_z = self.emitter.Transmission_Coefficient(-self._Ef-self._energy_valence_band) - (b*self.emitter.Transmission_Coefficient(-self._Ef-b*self._energy_valence_band+a*self._Evhigh))
        
        cumulative_transmission = np.insert(np.cumsum((transmission_in_z[1:]+transmission_in_z[:-1])*(self._energy_valence_band[1]-self._energy_valence_band[0])/2), 0, 0)
            
        number_of_electrons = self.emitter.Fermi_Dirac_Distribution(self._energy_valence_band, self._kT) * cumulative_transmission 
      
        return self._energy_valence_band, number_of_electrons
    
    def Nottingham_Heat_from_Semiconductors(self):
        """Returns the Nottigham heat resulted from electrons being emitter fromt the conduction and valence band

        Returns:
            nottigham_conduction (float): Conduction band contribution to Nottingham heat (J)
            nottigham_valence (float): Valence band contribution to Nottingham heat (J)
            nottigham_total (float): Total Nottingham heat (J)
        """
        
        nottingham_conduction = self.Nottingham_Heat_from_Conduction_Band()
        nottigham_valence = self.Nottingham_Heat_from_Valence_Band()
        nottigham_replacement = self.Nottingham_Heat_from_Replacement_Electrons()
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
    
    def Nottingham_Heat_from_Replacement_Electrons(self):
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

def current_metal_emitter(Field, Radius, Gamma, Workfunction, Temperature):
 
    tab = Tabulator()

    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature
    
    metal_emitter = Metal_Emitter(tab)

    j_metal = np.copy(Field)
    
    for i in range(len(Field)):

        metal_emitter.emitter.Define_Barrier_Parameters(Field[i], Radius[i], Gamma[i])
        metal_emitter.emitter.Interpolate_Gammow()
    
        metal_emitter.Define_Emitter_Parameters(Workfunction[i], kT[i])
    
        j_metal[i] = metal_emitter.Current_Density()
        
    return j_metal

def heat_metal_emitter(Field, Radius, Gamma, Workfunction, Temperature):
    
    tab = Tabulator()
    
    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature
    
    metal_emitter = Metal_Emitter(tab)
    
    nh_metal = np.copy(Field)

    for i in range(len(Field)):
        metal_emitter.emitter.Define_Barrier_Parameters(Field[i], Radius[i], Gamma[i])
        metal_emitter.emitter.Interpolate_Gammow()
    
        metal_emitter.Define_Emitter_Parameters(Workfunction[i], kT[i])
    
        nh_metal[i] = metal_emitter.Nottingham_Heat()

    return nh_metal

def current_semiconductor_emitter(Field:float ,Ec:float, Ef:float, Eg:float):
    debug = 1
    tab = Tabulator()
    
    kBoltz = 8.6173324e-5 
    Temperature = 300
    kt = kBoltz * Temperature
    
    mass_e = 9.1093837015e-31 
    effec_e = 0.98*mass_e
    effec_p = 0.59*mass_e
    radius = 20
    gamma = 10
    
    semiconductor_emitter = Semiconductor_Emitter(tab)

    j_total = np.copy(Field)
    j_c = np.copy(Field)
    j_v = np.copy(Field)


    Radius = np.ones(len(Field))*radius
    Gamma = np.ones(len(Field))*gamma
    kT = np.ones(len(Field))*kt
    m = np.ones(len(Field))*mass_e
    me = np.ones(len(Field))*effec_e
    mp = np.ones(len(Field))*effec_p

    for i in range(len(Field)):

        semiconductor_emitter.emitter.Define_Barrier_Parameters(Field[i], Radius[i], Gamma[i])
        semiconductor_emitter.emitter.Interpolate_Gammow()

        semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(Ec[i], Ef[i], Eg[i], kT[i], m[i], me[i], mp[i])
        
        j_c[i], j_v[i], j_total[i] = semiconductor_emitter.Current_Density_from_Semiconductors()
        #pn_c, pn_v, pn_total = semiconductor_emitter.Nottingham_Heat_from_Semiconductors()
    #energy_c, distribution_c, energy_v, distribution_v = semiconductor_emitter.Energy_Distribution_from_Semiconductors()
    
    #if debug == 1:
    #y = time.time()

    file = open("Field_Emission_COMSOL_Debugging_Data.txt", "a")
    file.write("Array length = %d" % len(Field))
    file.write("\r\nArray length = %d" % len(Ec))
    file.write("\r\nCurrent density        Field    Radius    Gamma     Ec    Ef    Eg     kT            m                    me                   mp\r\n")
        

    for i in range(len(Field)):
            #file.write("\r\n%e %e %e %e %e %e %e %e %e %e " % (j_total[i]) (Field[i]) (Radius[i]) (Gamma[i]) (Ec[i]) (Ef[i]) (Eg[i]) (m[i]) (me[i]) (mp[i]))
        L = [str(j_total[i]),"   ",str(Field[i]), "     ",str(Radius[i]), "     ",str(Gamma[i]),"     ", str(Ec[i]),"   ", str(Ef[i]),"   ", str(Eg[i]), "   ", str(kT[i]/kBoltz), "   ",str(m[i]), "   ",str(me[i]),"   ", str(mp[i])] 
        file.writelines(L) 
        file.write("\r\n")
    file.close()

    return j_total

"""Add Pn for semi and e- energy distributions"""

data = current_metal_emitter(np.ones(10)*10, np.ones(10)*20,np.ones(10)*10,np.ones(10)*4.5, np.ones(10)*300)
end = time.time()
print(data)

