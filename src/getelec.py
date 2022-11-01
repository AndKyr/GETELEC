
from array import array
import ctypes as ct
import numpy as np
from scipy.integrate import IntegrationWarning
import scipy.integrate as ig
import os
import scipy.ndimage
import matplotlib.pyplot as plt

pythonpath,filename = os.path.split(os.path.realpath(__file__))

integrator = ct.CDLL(pythonpath + '/libintegrator.so') #use absolute path
integrator.intfun.restype = ct.c_double
integrator.intfun.argtypes = (ct.c_int, ct.c_double)
integrator.intfun_Pn.restype = ct.c_double
integrator.intfun_Pn.argtypes = (ct.c_int, ct.c_double)
integrator.Gfun.argtypes = (ct.c_int, ct.c_void_p)
integrator.Gfun.restype = ct.c_double
integrator.intfun_dbg.argtypes = (ct.c_int, ct.c_void_p)
integrator.intfun_dbg.restype = ct.c_double

zs = 1.6183e-4

class Interpolator:
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
            
    2) _loadTables()
        Looks for the files where the precaculate barriers are stored. Then it uses interpolation methods to make the most accurate barrier for the given 
        input (electric field, tip radius and gamma exponent). Gtab is stores the polinomial that gives its shape to the barrier.
        
    3) interpolateForValues()
        This function is called if the barrier maps are not found. For each field, tip radius and gamma combitation, Gtab is calculated by polynomial fitting 
        to the data points extracted from getelec.export_gamow
        For each value of the field, the radius and gamma, the Gammow coeffient (G) and the energy range (W), for such G, are evaluated
            Wmin is defined as the energy at the top of the potential barrier. Point from which G goes to 0 and tranmission cofficient equals 1
            Wmax is defined as the energy where G gets a values high enough, here G=50, so the tunnelling probability is effectively 0
        Then it tries to fit G(W) to a polynomial. Such polynomial and the energy limits are stored in Gtab
    """
    #define variable types
    _gamowTable : np.ndarray
    _dimension : int
    _inverseFieldLimits : np.ndarray(shape=(2,), dtype=np.float64)
    _inverseRadiusLimits : np.ndarray(shape=(2,), dtype=np.float64)
    _inverseGammaLimits : np.ndarray(shape=(2,), dtype=np.float64)
    _Nfields : int
    _Ngamma : int
    _Npolynomial : int 
    _Nradius : int  

    def __init__(self, NField = 256, NRadius = 128, Ngamma = 1, Npolynomial = 4, NGamow = 128):
        """Initialises Tabulator by loading the tables from files. If maps not found, calls tabulation scripts
        from old getelec

        Args:
            tabulationFolder (int, optional): folder where to fine tabulated values
            NField, NRadius, Ngamma (int, optional): Number of points to tabulate for electric field, radius and gamma.
            They are relevant only in case the tabulated files are not found and the tabulation script from oldGETELEC
            is invokes. 
        """

        if self._loadTables():
            pass
        else:
            tabulationScript = pythonpath + "/../oldGETELEC/python/tabulateGamow.py"
            command = "python3 %s -nf %d -nr %d -ng %d -np %d -nG %d"%(tabulationScript, NField, NRadius, Ngamma, Npolynomial, NGamow)
            print("tabulation not found. producing tabulator by running:")
            print(command)
            os.system(command)
            self._loadTables()

    def _loadTables(self) -> bool:
        """Loads the table of the Gamow factor from the files where it has been stored.
        """
        try:
            self._gamowTable = np.load("tabulated/GamowTable.npy")
            limits = np.load("tabulated/tabLimits.npy")
        except(IOError):
            print("tabulation files not found")
            return False

        #Get interpolation limits
        self._inverseFieldLimits = limits[:2]
        self._inverseRadiusLimits = limits[2:4]
        self._inverseGammaLimits = limits[4:]
        self._dimension = len(np.shape(self._gamowTable)) - 1
        
        #make sure that dimensions and limits are correct
        assert self._dimension <= 3 and self._dimension >= 1, \
            'The tabulated Gamow table has wrong shape: {}'.format(np.shape(self._gamowTable))
        assert self._inverseFieldLimits[0] < self._inverseFieldLimits[1], \
            'Limits of tabulated Field: {} are wrong'.format(self._inverseFieldLimits)

        #set sizes and assert that everything is consistent
        self._Npolynomial = np.shape(self._gamowTable)[-1] - 2       
        self._Nfields = np.shape(self._gamowTable)[-2]

        if(self._dimension >= 2):
            assert self._inverseRadiusLimits[0] < self._inverseRadiusLimits[1], \
                'Limits of tabulated Radius: {} are wrong'.format(self._inverseRadiusLimits)
            self._Nradius = np.shape(self._gamowTable)[-3]
        else:
            assert self._inverseRadiusLimits[0] == self._inverseRadiusLimits[1], \
                'Limits of tabulated Radius: {} are wrong'.format(self._inverseRadiusLimits)
            self._Nradius = 1

        if (self._dimension == 3):
            self._Ngamma = np.shape(self._gamowTable[0])
            assert self._inverseGammaLimits[0] < self._inverseGammaLimits[1], \
                'Limits of tabulated Gamma: {} are wrong'.format(self._inverseGammaLimits)
        else:
            assert self._inverseGammaLimits[0] == self._inverseGammaLimits[1], \
                'Limits of tabulated Gamma: {} are wrong'.format(self._inverseGammaLimits)
            self._Ngamma = 1
        
        return True

    
    # region user methods
    def interpolateForValues(self, field : float, radius : float = 1.e4, gamma : float = 10., interpolationOrder = 2.) -> np.ndarray:
        """Interpolate for given values of field, radius and gamma. Returns an array of numbers. 
        The first 4 numbers are the polynomial coefficient"""
        
        paramCoordinates = np.arange(self._Npolynomial + 2)
        fieldCoorinates = np.ones(len(paramCoordinates)) * (self._Nfields - 1) * (1./field - self._inverseFieldLimits[0]) / (self._inverseFieldLimits[1] - self._inverseFieldLimits[0])

        if (self._dimension >= 2):
            radiusCoordinates = np.ones(len(paramCoordinates)) * (self._Nradius - 1) * (1./radius - self._inverseRadiusLimits[0]) / (self._inverseRadiusLimits[1] - self._inverseFieldLimits[0])
        else:
            radiusCoordinates = np.zeros(len(paramCoordinates))
            if (radius < 100.):
                print("WARNING: Using tabulation without radius (planar approximation) with input radius: {} nm < 100 nm.".format(radius))
        
        if (self._dimension == 3):
            gammaCoordinates = np.ones(len(paramCoordinates)) * (self._Ngamma - 1) * (1./gamma - self._inverseGammaLimits[0]) / (self._inverseGammaLimits[1] - self._inverseGammaLimits[0])
        else:
            gammaCoordinates = np.zeros(len(paramCoordinates))
            if (gamma != 10.):
                print("WARNING: Using tabulation without gamma with gamma: {} nm != 10. Tabulation is done with gamma = 10".format(gamma))


        if (self._dimension == 1):
            interpolationCoordinates = np.array([fieldCoorinates, paramCoordinates])
        elif(self._dimension == 2):
            interpolationCoordinates = np.array([radiusCoordinates, fieldCoorinates, paramCoordinates])
            print(radiusCoordinates, fieldCoorinates, paramCoordinates)
        elif (self._dimension == 3):
            interpolationCoordinates = np.array([gammaCoordinates, radiusCoordinates, fieldCoorinates, paramCoordinates])

        return scipy.ndimage.map_coordinates(self._gamowTable, interpolationCoordinates, order = interpolationOrder)

class Barrier(Interpolator):

    """
    This class gathers those attributes that are common to all emitters regardless of the material they are made of.
    Inherits from Interpolator, as it uses its attributes to calculate the Gamow(EnergyDepth) polynomial coefficients and the interval of validity.
    
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

    """
    # region field declarations
    _field: float
    _radius: float
    _gamma: float
    _data: array
    _energy: array
    energy_top_barrier: float
    _maxEnergyDepth:float
    _gamowPolynomialCoefficients: array
    _gamowDerivativeAtMinBarrierDepth: float
    _gamowDerivativeAtMaxEnergyDepth: float
    _gamowDerivativePolynomialCoefficients: array
    # endregion
    
    # region initialization
    def __init__(self, field : float, radius : float, gamma : float):
        super().__init__()
        self.setBarrierParameters(field, radius, gamma)
   
    def _calculateParameters(self):
        """Here a function from Tabulator is called in order to interpolate the potential barrier. Then the top and bottom limits, as well as the polynomial
        describing the barrier's shape, are extracted. The polinomial coefficients are derivated and evaluated at the top and bottom of the barrier to analyse
        the barrier slope at such points - which will be required to later calculate the integration limits for the currrent density and heat intgrals
        """
        
        data = self.interpolateForValues(self._field, self._radius, self._gamma)
        self._minEnergyDepth = data[-2] #Wmin: This is the top of the barrier
        self._maxEnergyDepth = data[-1] #Wmax: The maximum barrier depth for which the polynomial is valid
        self._gamowPolynomialCoefficients = data[:self._Npolynomial] #Gpoly
        self._gamowDerivativePolynomialCoefficients = np.polyder(self._gamowPolynomialCoefficients) #dG/dW polynomial coefficients
        self._gamowDerivativeAtMinBarrierDepth = np.polyval(self._gamowDerivativePolynomialCoefficients, self._minEnergyDepth) #dG/dW @ Wmin
        self._gamowDerivativeAtMaxEnergyDepth = np.polyval(self._gamowDerivativePolynomialCoefficients, self._maxEnergyDepth) #dG/dW @ Wmax

    def calculateGamowOutOfBounds(self, energyDepth:array):

        if (1./self._radius < self._inverseRadiusLimits[0]):
            pass


    def calculateGamowForEnergy(self, energyDepth:array):
        """Evaluates the gammow coefficients in order to provide the exponential for the Kemble formula

        Parameters:
            energy (array): Incoming energy, perpendicular to the emitting surface, values (eV) for which the transmission cofficients are calculated (eV)

        Returns:
            array: Kemble's formula Gamow exponent value for a given energy (number)
        """
        self._calculateParameters()
        
        if (isinstance(energyDepth, np.ndarray)):
            Gamow = np.copy(energyDepth)
            highEnergyDepth = energyDepth > self._maxEnergyDepth
            lowEnergyDepth = energyDepth < self._minEnergyDepth
            Gamow = np.polyval(self._gamowPolynomialCoefficients, energyDepth)
            Gamow[highEnergyDepth] = np.polyval(self._gamowPolynomialCoefficients, self._maxEnergyDepth) + \
                self._gamowDerivativeAtMaxEnergyDepth * (energyDepth[highEnergyDepth] - self._maxEnergyDepth)
            Gamow[lowEnergyDepth] = np.polyval(self._gamowPolynomialCoefficients, self._minEnergyDepth) + \
                self._gamowDerivativeAtMinBarrierDepth * (energyDepth[lowEnergyDepth] - self._minEnergyDepth)
            return Gamow
        else:
            if (energyDepth > self._minEnergyDepth and energyDepth < self._maxEnergyDepth):
                return np.polyval(self._gamowPolynomialCoefficients, energyDepth)
            elif(energyDepth > self._maxEnergyDepth):
                return np.polyval(self._gamowPolynomialCoefficients, self._maxEnergyDepth) + \
                    self._gamowDerivativeAtMaxEnergyDepth * (energyDepth - self._maxEnergyDepth)
            else:
                return np.polyval(self._gamowPolynomialCoefficients, self._minEnergyDepth) + \
                    self._gamowDerivativeAtMinBarrierDepth * (energyDepth - self._minEnergyDepth) 

    def setBarrierParameters(self, field:float, radius:float = 1.e4, gamma:float = 10.):
        """Sets the main parameters of the barrier.

        Parameters:
            field (float): Electric field (V/nm)
            radius (float): Emitter tip raidus (nm). Default value 1.e4 in case it is omitted for planar approximation
            gamma (float): gamma coefficient (number). Default value 10 in case it is omitted
        """
        
        self._field = field
        self._radius = radius
        self._gamma = gamma
    
    def transmissionCoefficient(self, energyDepth:array):
        """Calculates the probability of an electron tunelling (being transmitted) through the potential barrier by 
        means of the Kemble formula within the JWKB approximation

        Parameters:
            energyDepth (array): Energy depth of the barrier, i.e. incoming electron energy perpendicular to the emitting surface, 
            measured downwards from vacuum level. The Energy values (eV) for which the transmission cofficients are calculated

        Returns:
            array: Probability of an electron with a certain kinetic energy perpendicular to the emitting surface to tunnerl through the potential barrier (number)
        """
        
        gamow = self.calculateGamowForEnergy(energyDepth)
        if (isinstance(energyDepth, np.ndarray)):
            transmissionCoefficient = np.copy(gamow)
            transmissionCoefficient[gamow < 15.] = 1 / (1 + np.exp(gamow[gamow < 15.]))
            transmissionCoefficient[gamow > 15] = np.exp(-gamow[gamow > 15.])
            return transmissionCoefficient
        else:
            if (gamow < 15.):
                return 1 / (1 + np.exp(gamow))
            else:
                return np.exp(-gamow)
    
    
    def plotGamow(self, minBarrierDepth : float = 0., maxBarrierDepth : float = 0., show = False, saveFile = ""):
        """Plots the Gamow exponent for energy depth between minBarrierDepth and maxBarrierDepth. IF saveFile is given,
        it will save the file to the given output. If show = True, it will show the output plot."""
        
        barrierDepths = np.linspace(minBarrierDepth, maxBarrierDepth, 256)
        gamowValues = self.calculateGamowForEnergy(barrierDepths)
        plt.plot(barrierDepths, gamowValues)
        plt.xlabel("Barrier Depth [eV]")
        plt.ylabel("Gamow Factor")
        if (show):
            plt.show()
        
        if (saveFile):
            plt.savefig(saveFile)

    def plotTransmissionCoefficient(self, minBarrierDepth : float = 0., maxBarrierDepth : float = 0., show = False, saveFile = ""):
        "Similar to plotGamow, but plots the transmission coefficient"
        barrierDepths = np.linspace(minBarrierDepth, maxBarrierDepth, 256)
        TValues = self.transmissionCoefficient(barrierDepths)
        plt.plot(barrierDepths, TValues)
        plt.xlabel("Barrier Depth [eV]")
        plt.ylabel("Transmission Coefficient")
        if (show):
            plt.show()
        
        if (saveFile):
            plt.savefig(saveFile)
    

class Supply:
    # region field declaration
    kT: array
    distribution:array
    energy:array
    _high_energy:array
    _low_energy:array
    _mid_energy:array
    # endregion

    # region inisialisation
    def __init__(self):
        pass
    # endregion 
    
    # region user methods
    def Log_Fermi_Dirac_Distribution(self, energy:array, kT:array):
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

            _high_energy = energy > 10. * kT
            _low_energy = energy < -10. * kT
            _mid_energy = np.logical_not(_high_energy) * np.logical_not(_low_energy)

            distribution[_high_energy] = np.exp(-energy[_high_energy] / kT)
            distribution[_low_energy] = - energy[_low_energy] / kT
            distribution[_mid_energy] = np.log(1. + np.exp(-energy[_mid_energy] / kT))
            return distribution
        else:
            if energy > 10 * kT:
                return np.exp(-energy / kT)
            elif(energy < -10 * kT):
                return -energy / kT
            else:
                return np.log(1. + np.exp(-energy / kT))
    
    def Fermi_Dirac_Distribution(self, energy:array, kT:array):
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

class Emitter:
    
    # region initialization
    def __init__(self, barrier: Barrier, supply:Supply):
        """Initialises the class by inheritating Tabulator and setting it up as a private field to be used by Emitter

        Args:
            tabulator (Tabulator): Object that contains the attributes to interpolate the potential barrier given the field, radius and gamma as inputs
        """
        self.barrier = barrier
        self.supply = supply
    # endregion

   
class Metal_Emitter(Emitter):
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
    def __init__(self, barrier:Barrier, supply:Supply):
        """Initialises the function by adding the atributes of a "tabulated" barrier emitter to metals

        Args:
            tabulator (Tabulator): Object that contains the attributes of a x material emitter (potential barrier, transmission coefficiens and electron supply functions)
        """
        
        super().__init__(barrier, supply)
        
    def _Integration_Limits(self):
        """Finds the limits of integration. ANDREAS PLEASE MAKE A COMMENT WHY THE 3 DIFFERENT INSTANCES AND WHY THE VALUES
        
        Returns:
            array: Energy space for which the electron emission and Notigham heat are going to be evaluated in (eV)
        """
        
        resolution = NGi #128
        self._maxbeta = np.polyval(self.barrier._gamowDerivativePolynomialCoefficients, min(self.workfunction, self.barrier._maxEnergyDepth))
        
        if (self._maxbeta * self.kT < 1.05): #field regime
            workfunction_center = 0.
            energy_top = 10 /(1/self.kT - .85*self._maxbeta)
            energy_bottom = -10 / self._maxbeta
            return  np.linspace(energy_bottom, energy_top, resolution)
        
        elif (self.barrier._gamowDerivativeAtMinBarrierDepth * self.kT > .95):
            workfunction_center = self.workfunction - self.barrier._minEnergyDepth
            energy_top = workfunction_center + 10 * self.kT
            energy_bottom = workfunction_center - 10 / self.barrier._gamowDerivativeAtMinBarrierDepth
            return np.linspace(energy_bottom, energy_top, resolution)
        
        else:
            rootpoly = np.copy(self.barrier._gamowDerivativePolynomialCoefficients)
            rootpoly[-1] -= 1./self.kT
            realroots = np.roots(rootpoly)
            workfunction_center = realroots[np.nonzero(np.logical_and(realroots > self.barrier._minEnergyDepth, realroots < self.workfunction))][0]
            #print (Wcenter)
            energy_top = self.workfunction - workfunction_center + 10 * self.kT
            energy_bottom = self.workfunction - workfunction_center - 25 * self.kT
            return np.linspace(energy_bottom, energy_top, resolution)
     
    # endregion   
    
    # region user methods 
    def Define_Metal_Emitter_Parameters(self, workfunction:float, kT:float):
        """Defines main emitter characteristics

        Args:
            workfunction (float): Average minimum energy required to put an electron, whitin a material, out in vacuum (some few nm away from the material surface) (eV)
            kT (float): Energy from temperature (eV)
        """
        
        self.workfunction = workfunction
        self.kT = kT
        self.energy = self._Integration_Limits()
    
    def Current_Density_from_Metals(self):
        """Calculates the field emitted current density from metal surfaces

        Returns:
            float: Emitted current density being emitted from an infinitesimal flat surface area (electrons/area*time)
        """
        
        args = tuple([self.workfunction] + [self.kT] + [self.barrier._minEnergyDepth] + [self.barrier._maxEnergyDepth] + [self.barrier._gamowDerivativeAtMinBarrierDepth] + [self.barrier._gamowDerivativeAtMaxEnergyDepth] + list(self.barrier._gamowPolynomialCoefficients))
        try:
            integ, abserr, info = ig.quad(integrator.intfun, self.energy[0], self.energy[-1], args, full_output = 1)
            self.integ_points = info["alist"]
        except(IntegrationWarning):
            integ = 0.
        return zs * self.kT * integ
    
    def Energy_Distribution_for_Metals(self):
        """Calculates the energy distribution of field emitted electrons from metals

        Returns:
            energy (array): Energy space for which the electron emission has been evalated (eV)
            electron_number (array): Number of electrons being emitter with a certain energy (number)
        """
        
        transmission_in_energy_z = self.barrier.transmissionCoefficient(self.workfunction - self.energy)

        cumulative_transmission = np.insert(np.cumsum((transmission_in_energy_z[1:]+transmission_in_energy_z[:-1])*(self.energy[1]-self.energy[0])/2), 0, 0)
        
        electron_number = self.supply.Fermi_Dirac_Distribution(self.energy, self.kT) * cumulative_transmission

        return  self.energy, electron_number
  
    def Nottingham_Heat_from_Metals(self):
        """Calculates the Nottingham heat resulted from field emitted electrons from metals

        Returns:
            float: Resulted Nottigham from the emission of electrons from an infinitesinal area (J)
        """

        args = tuple([self.workfunction] + [self.kT] + [self.barrier._minEnergyDepth] + [self.barrier._maxEnergyDepth] + [self.barrier._gamowDerivativeAtMinBarrierDepth] + [self.barrier._gamowDerivativeAtMaxEnergyDepth] + list(self.barrier._gamowPolynomialCoefficients))
        try:
            integ, abserr = ig.quad(integrator.intfun_Pn, self.energy[0], self.energy[-1], args, full_output = 0)
        except:
            return 0.
        
        return -zs * self.kT * integ    
    
    def Current_Density_from_Metals_educational(self):
        """Calculates the field emitted current density from metals on a clear equation to code manner"""
        integ = self.supply.Log_Fermi_Dirac_Distribution(self.energy, self.kT) * self.Transmission_Coefficient(self.workfunction - self.energy)
        
        """if(plot):
            print("Whigh, Wlow = ", self.workfunction - self.energy[0], self.workfunction - self.energy[-1])
            print ("Wmin, Wmax, Work = ", self.barrier.energy_bottom_barrier, self.barrier._energy_top_barrier, self.workfunction)
            print("maxbeta, minbeta, dGmax, beta_T: ", self._maxbeta, self.barrier.gammow_derivative_bottom_barrier, self.barrier.gammow_derivative_top_barrier, 1/ self.kT)
            plt.plot(self.energy,integ)
            ax = plt.gca()
            ax2 = ax.twinx()
            ax2.plot(self.energy,self.Gamow(self.workfunction - self.energy), 'r-')
            ax.grid()
            plt.savefig("Jcur.png")
            # plt.show()"""
        return zs * self.kT * np.sum(integ) * (self.energy[1] - self.energy[0])
    
    def Energy_Distribution_of_Metals_educational(self):
        """Calculates the energy distribution of field emitted electrons from metals on a clear equation to code manner"""
        transmission_in_energy_z = self.Transmission_Coefficient(self.workfunction - self.energy)

        cumulative_transmission = np.copy(self.energy)
        cumulative_transmission[0] = 0

        for i in range(len(self.energy)-1):
            cumulative_transmission[i+1] = cumulative_transmission[i]+((self.energy[1]-self.energy[0])*(transmission_in_energy_z[i+1]+transmission_in_energy_z[i])/2)
        
        energy_distribution = self.supply.Fermi_Dirac_Distribution(self.energy, self.kT) * cumulative_transmission

        return  self.energy, energy_distribution
    
    def Nottingham_Heat_from_Metals_educational(self):
        """Calculates the Nottingham heat resulted from field emitted electrons from metals on a clear equation to code manner"""
        transmission_in_energy_z = self.Transmission_Coefficient(self.workfunction - self.energy)

        cumulative_transmission = np.copy(self.energy)
        cumulative_transmission[0] = 0

        for i in range(len(self.energy)-1):
            cumulative_transmission[i+1] = cumulative_transmission[i]+((self.energy[1]-self.energy[0])*(transmission_in_energy_z[i+1]+transmission_in_energy_z[i])/2)
        
        heat_components = self.supply.Fermi_Dirac_Distribution(self.energy, self.kT) * cumulative_transmission * self.energy
        
        return -zs * np.sum(heat_components) * (self.energy[1] - self.energy[0])
    #endregion
    
class Semiconductor_Emitter(Emitter):
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
    def __init__(self, barrier: Barrier, supply:Supply):
        """Initialises the function by adding the atributes of a "tabulated" barrier emitter to metals

        Args:
            tabulator (Tabulator): Object that contains the attributes of a x material emitter (potential barrier, transmission coefficiens and electron supply functions)
        """
        
        super().__init__(barrier, supply)
    # endregion

    # region user methods
    def Define_Semiconductor_Emitter_Parameters(self, Ec:float, Ef:float, Eg:float, kT:float, m:float, me:float, mp:float):
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
        self._Evlow = min(self._Evhigh, 0) - 20 / self.barrier._gamowDerivativeAtMaxEnergyDepth
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
        
        jc_integ = self.supply.Log_Fermi_Dirac_Distribution(self._energy_conduction_band, self._kT) * (self.barrier.transmissionCoefficient(-self._Ef-self._energy_conduction_band) - (b*self.barrier.transmissionCoefficient(-self._Ef-a*self._energy_conduction_band)))
        
        return zs * self._kT * np.sum(jc_integ) * (self._energy_conduction_band[1]-self._energy_conduction_band[0]) 
    
    def Current_from_Valence_Band(self):
        """Calculates the conduction band component of the emitted current - Andreas' equation

        Returns:
            float: Emitted current density from the valence band (electrons/area*time)
        """
        
        a = self._mp/self._m
        b = 1+a
        
        jv_integ = self.supply.Log_Fermi_Dirac_Distribution(self._energy_valence_band, self._kT) * (self.barrier.transmissionCoefficient(-self._Ef-self._energy_valence_band) - (b*self.barrier.transmissionCoefficient(-self._Ef-b*self._energy_valence_band+a*self._Evhigh)))
        
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
        
        c_transmission_in_z = self.Transmission_Coefficient(-self._Ef-self._energy_conduction_band) - (b*self.Transmission_Coefficient(-self._Ef-a*self._energy_conduction_band))

        c_cumulative_transmission = np.insert(np.cumsum((c_transmission_in_z[1:]+c_transmission_in_z[:-1])*(self._energy_conduction_band[1]-self._energy_conduction_band[0])/2), 0, 0)
        
        c_number_of_electrons = np.sum(self.supply.Fermi_Dirac_Distribution(self._energy_conduction_band, self._kT) * c_cumulative_transmission )
        
        total_c = np.sum(zs * c_number_of_electrons * (self._energy_conduction_band[1]-self._energy_conduction_band[0]))
        
        #valence band
        c = self._mp/self._m
        d = 1+c
        
        v_transmission_in_z = self.Transmission_Coefficient(-self._Ef-self._energy_valence_band) - (d*self.Transmission_Coefficient(-self._Ef-d*self._energy_valence_band+c*self._Evhigh))
        
        v_cumulative_transmission = np.insert(np.cumsum((v_transmission_in_z[1:]+v_transmission_in_z[:-1])*(self._energy_valence_band[1]-self._energy_valence_band[0])/2), 0, 0)
            
        v_number_of_electrons = self.supply.Fermi_Dirac_Distribution(self._energy_valence_band, self._kT) * v_cumulative_transmission 
        
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
        
        transmission_in_z = self.barrier.transmissionCoefficient(-self._Ef-self._energy_conduction_band) - (b*self.barrier.transmissionCoefficient(-self._Ef-a*self._energy_conduction_band))

        cumulative_transmission = np.insert(np.cumsum((transmission_in_z[1:]+transmission_in_z[:-1])*(self._energy_conduction_band[1]-self._energy_conduction_band[0])/2), 0, 0)
        
        number_of_electrons = self.supply.Fermi_Dirac_Distribution(self._energy_conduction_band, self._kT) * cumulative_transmission 

        return self._energy_conduction_band , number_of_electrons
    
    def Energy_Distribution_from_Valence_Band(self): #NEED REVISION
        """Calculates teh energy distribution of those electrons emitted from the valence band

        Returns:
            _energy_conduction_band (array): Energy values in the valence band from which electrons are emitted (eV)
            number_of_electrons (array): Number of emitted electrons from the valence band for each energy value (number)
        """
        
        a = self._mp/self._m
        b = 1+a
        
        transmission_in_z = self.barrier.transmissionCoefficient(-self._Ef-self._energy_valence_band) - (b*self.barrier.transmissionCoefficient(-self._Ef-b*self._energy_valence_band+a*self._Evhigh))
        
        cumulative_transmission = np.insert(np.cumsum((transmission_in_z[1:]+transmission_in_z[:-1])*(self._energy_valence_band[1]-self._energy_valence_band[0])/2), 0, 0)
            
        number_of_electrons = self.supply.Fermi_Dirac_Distribution(self._energy_valence_band, self._kT) * cumulative_transmission 
      
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
        energy, distribution = self.Energy_Distribution_from_Conduction_Band()
        
        heat = energy
        
        return -zs* self._kT * np.sum(heat*distribution)
    
    def Nottingham_Heat_from_Valence_Band(self): #NEED MOFICATIONS
        """Calculates the contribution of the valence band to the Nottingham heat resulted from field emitted electrons from semicoductors"""
        """the fuction will be replaced by the full code from the function, if we agree the function works well"""
        energy, distribution = self.Energy_Distribution_from_Valence_Band()#*self._energy_valence_band

        heat = energy

        return -zs * self._kT * np.sum(heat*distribution)

    def Nottingham_Heat_from_Replacement_Electrons(self):
        """To be calculated from COMSOL"""
        return 0
    
    def Energy_Distribution_from_Conduction_Band_educational(self):
        """Calculates the contribution of the conduction band to the energy distribution of field emitted electrons from semiconductors on a clear equation to code manner"""
        a = self._me/self._m
        b = 1-a 
        
        transmission_in_z = self.Transmission_Coefficient(-self._Ef-self._energy_conduction_band) - (b*self.Transmission_Coefficient(-self._Ef-a*self._energy_conduction_band))

        cumulative_transmission = np.copy(self._energy_conduction_band)
        cumulative_transmission[0] = 0

        for i in range(len(self._energy_conduction_band)-1):
            cumulative_transmission[i+1] = cumulative_transmission[i] + ((self._energy_conduction_band[1]-self._energy_conduction_band[0])*(transmission_in_z[i+1]+transmission_in_z[i])/2)
        
        energy_distribution = self.supply.Fermi_Dirac_Distribution(self._energy_conduction_band, self._kT) * cumulative_transmission 

        return self._energy_conduction_band , energy_distribution
    # endregiondata


def current_metal_emitter(Field:array, Radius:array, Gamma:array, Workfunction:array, Temperature:array):

    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature
    
    metal_emitter = Metal_Emitter(Barrier(),Supply())

    j_metal = np.copy(Field)
    
    for i in range(len(Field)):

        metal_emitter.barrier.setBarrierParameters(Field[i], Radius[i], Gamma[i])
        metal_emitter.barrier._calculateParameters()
       
        metal_emitter.Define_Metal_Emitter_Parameters(Workfunction[i], kT[i])
    
        j_metal[i] = metal_emitter.Current_Density_from_Metals()
        
    return j_metal

def heat_metal_emitter(Field:array, Radius:array, Gamma:array, Workfunction:array, Temperature:array):
    
    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature
    
    metal_emitter = Metal_Emitter(Barrier(),Supply())
    
    nh_metal = np.copy(Field)

    for i in range(len(Field)):
        metal_emitter.barrier.setBarrierParameters(Field[i], Radius[i], Gamma[i])
        metal_emitter.barrier._calculateParameters()
    
        metal_emitter.Define_Metal_Emitter_Parameters(Workfunction[i], kT[i])
    
        nh_metal[i] = metal_emitter.Nottingham_Heat_from_Metals()

    return nh_metal

def heat_metal_emitter2(Field:array, Radius:array, Gamma:array, Workfunction:array, Temperature:array):

    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature
    
    metal_emitter = Metal_Emitter(Barrier(),Supply())
    
    nh_metal = np.copy(Field)

    for i in range(len(Field)):
        metal_emitter.barrier.setBarrierParameters(Field[i], Radius[i], Gamma[i])
        metal_emitter.barrier._calculateParameters()
    
        metal_emitter.Define_Metal_Emitter_Parameters(Workfunction[i], kT[i])
    
        nh_metal[i] = metal_emitter.Nottingham_Heat_from_Metals_educational()

    return nh_metal

def current_semiconductor_emitter(Field:array, Radius:array, Gamma:array ,Ec:array, Ef:array, Eg:array, Temperature:array):

    kBoltz = 8.6173324e-5 

    kT = kBoltz * Temperature
    
    mass_e = 9.1093837015e-31 
    effec_e = 0.98*mass_e
    effec_p = 0.59*mass_e
    
    semiconductor_emitter = Semiconductor_Emitter(Barrier(),Supply())

    j_total = np.copy(Field)
    j_c = np.copy(Field)
    j_v = np.copy(Field)


    m = np.ones(len(Field))*mass_e
    me = np.ones(len(Field))*effec_e
    mp = np.ones(len(Field))*effec_p

    for i in range(len(Field)):

        semiconductor_emitter.barrier.setBarrierParameters(Field[i], Radius[i], Gamma[i])
        semiconductor_emitter.barrier._calculateParameters()

        semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(Ec[i], Ef[i], Eg[i], kT[i], m[i], me[i], mp[i])
        
        j_c[i], j_v[i], j_total[i] = semiconductor_emitter.Current_Density_from_Semiconductors()
 
    """
    file = open("Field_Emission_COMSOL_Debugging_Data.txt", "a")
    file.write("Array length = %d" % len(Field))
    file.write("\r\nArray length = %d" % len(Ec))
    file.write("\r\nCurrent density        Field    Radius    Gamma     Ec    Ef    Eg     kT            m                    me                   mp\r\n")
        

    for i in range(len(Field)):
            #file.write("\r\n%e %e %e %e %e %e %e %e %e %e " % (j_total[i]) (Field[i]) (Radius[i]) (Gamma[i]) (Ec[i]) (Ef[i]) (Eg[i]) (m[i]) (me[i]) (mp[i]))
        L = [str(j_total[i]),"   ",str(Field[i]), "     ",str(Radius[i]), "     ",str(Gamma[i]),"     ", str(Ec[i]),"   ", str(Ef[i]),"   ", str(Eg[i]), "   ", str(kT[i]/kBoltz), "   ",str(m[i]), "   ",str(me[i]),"   ", str(mp[i])] 
        file.writelines(L) 
        file.write("\r\n")
    file.close()"""
    
    return j_total

def heat_semiconductor_emitter(Field:array, Radius:array, Gamma:array ,Ec:array, Ef:array, Eg:array, Temperature:array):
    
    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature
    
    mass_e = 9.1093837015e-31 
    effec_e = 0.98*mass_e
    effec_p = 0.59*mass_e
    
    semiconductor_emitter = Semiconductor_Emitter(Barrier(),Supply())

    pn_total = np.copy(Field)
    pn_c = np.copy(Field)
    pn_v = np.copy(Field)


    m = np.ones(len(Field))*mass_e
    me = np.ones(len(Field))*effec_e
    mp = np.ones(len(Field))*effec_p

    for i in range(len(Field)):

        semiconductor_emitter.barrier.setBarrierParameters(Field[i], Radius[i], Gamma[i])
        semiconductor_emitter.barrier._calculateParameters()

        semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(Ec[i], Ef[i], Eg[i], kT[i], m[i], me[i], mp[i])
        
        pn_c[i], pn_v[i], pn_total[i] = semiconductor_emitter.Nottingham_Heat_from_Semiconductors()

    """
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
    """
    return pn_total



if (__name__ == "__main__"):
    bar = Barrier(4., 10, 10.)
    print(bar.calculateGamowForEnergy(4.5))
    # bar.plotGamow(0., 10., show = True)

    # bar.plotGamow(0., 10., show = True)


