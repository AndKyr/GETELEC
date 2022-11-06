from array import array
import ctypes as ct
import numpy as np
from scipy.integrate import IntegrationWarning
import scipy.integrate as ig
import os
import scipy.ndimage
import matplotlib.pyplot as plt

pythonpath,filename = os.path.split(os.path.realpath(__file__))

BoltzmannConstant = 8.617333262e-5

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

    def __init__(self, preloadedGamowTable : np.ndarray = None, preloadedLimits : np.ndarray = None, dataFolder = "tabulated"):
        """Initialises Tabulator by loading the tables from files. If maps not found, calls tabulation scripts
        from old getelec

        Parameters:
            tabulationFolder (int, optional): folder where to fine tabulated values
            NField, NRadius, Ngamma (int, optional): Number of points to tabulate for electric field, radius and gamma.
            They are relevant only in case the tabulated files are not found and the tabulation script from oldGETELEC
            is invoked. 
        """
        if (preloadedGamowTable is None or preloadedLimits is None):
            self._isTableLoaded = self._loadTablesFromFileIfPossible(dataFolder)
            if not self._isTableLoaded: #all the script that creates tables and try to load again
                self.calculateAndSaveTable()
                self._isTableLoaded = self._loadTablesFromFileIfPossible(dataFolder)
        else:
            self._gamowTable = preloadedGamowTable
            self._limits = preloadedLimits
            self._isTableLoaded = True
        
        self._initializeTabulationVariables()

    def calculateAndSaveTable(self, NField = 256, NRadius = 128, Ngamma = 1, Npolynomial = 4, NGamow = 128):
        """Runs old Getelec script that produces tabulation and saves it at tabulated/ folder."""

        tabulationScript = pythonpath + "/../oldGETELEC/python/tabulateGamow.py"
        command = "python3 %s -nf %d -nr %d -ng %d -np %d -nG %d"%(tabulationScript, NField, NRadius, Ngamma, Npolynomial, NGamow)
        print("Producing tabulator by running:")
        print(command)
        os.system(command)

    def _loadTablesFromFileIfPossible(self, dataFolder = "tabulated") -> bool:
        """Loads the table of the Gamow factor from the files where it has been stored.

        Returns:
            True if successful loading. False if table files not found
        """
        try:
            self._gamowTable = np.load(dataFolder + "/GamowTable.npy")
            self._limits = np.load(dataFolder + "/tabLimits.npy")
            return True
        except(IOError):
            print("tabulation files not found")
            return False

    def _initializeTabulationVariables(self) -> None:
        """Initializes class variables that facilitate tabulation operations
        
        Raises:
            Assertion error if tabularized data have inconsistencies.
        """

        assert self._isTableLoaded, "Attempting to initialize tabulation variables without having loaded the tables"

        #Get interpolation limits
        self._inverseFieldLimits = self._limits[:2]
        self._inverseRadiusLimits = self._limits[2:4]
        self._inverseGammaLimits = self._limits[4:]
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

    
    # region user methods
    def interpolateForValues(self, field:float, radius:float = 1.e4, gamma:float = 10., interpolationOrder:int = 1) -> np.ndarray:
        """
        Interpolate for given values of field, radius and gamma. Returns an array of numbers. 
        The first 4 numbers are the polynomial coefficient
        
        Parameters:
            field: field value for interpolation
            radius: radius value for interpolation. Default 1.e4, in case it is called for planar approximation
            gamma: gamma value for interpolation. Default 10. , in case it is called for gamma-is-irrelevant approximation
            interpolatorOrder: order of the spline used to interpolate values. Default is 2

        Returns:
            Array of Npolynomial + 2 values interpolated. The first Npolynomial are the polynomial coefficients and the last 2 are 
            the energyDepth limits of validity of the polynomial fitting.
        """
        
        paramCoordinates = np.arange(self._Npolynomial + 2)
        fieldCoorinates = np.ones(len(paramCoordinates)) * (self._Nfields - 1) * (1./field - self._inverseFieldLimits[0]) / \
            (self._inverseFieldLimits[1] - self._inverseFieldLimits[0])

        if (self._dimension >= 2):
            radiusCoordinates = np.ones(len(paramCoordinates)) * (self._Nradius - 1) * (1./radius - self._inverseRadiusLimits[0]) / \
                (self._inverseRadiusLimits[1] - self._inverseFieldLimits[0])
        else:
            radiusCoordinates = np.zeros(len(paramCoordinates))
            if (radius < 100.):
                print("WARNING: Using tabulation without radius (planar approximation) with input radius: {} nm < 100 nm.".format(radius))
        
        if (self._dimension == 3):
            gammaCoordinates = np.ones(len(paramCoordinates)) * (self._Ngamma - 1) * (1./gamma - self._inverseGammaLimits[0]) / \
                (self._inverseGammaLimits[1] - self._inverseGammaLimits[0])
        else:
            gammaCoordinates = np.zeros(len(paramCoordinates))
            if (gamma != 10.):
                print("WARNING: Using tabulation without gamma with gamma: {} nm != 10. Tabulation is done with gamma = 10".format(gamma))


        if (self._dimension == 1):
            interpolationCoordinates = np.array([fieldCoorinates, paramCoordinates])
        elif(self._dimension == 2):
            interpolationCoordinates = np.array([radiusCoordinates, fieldCoorinates, paramCoordinates])
        elif (self._dimension == 3):
            interpolationCoordinates = np.array([gammaCoordinates, radiusCoordinates, fieldCoorinates, paramCoordinates])

        return scipy.ndimage.map_coordinates(self._gamowTable, interpolationCoordinates, order = interpolationOrder, mode='nearest')

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
    maxEnergyDepth:float
    _gamowPolynomialCoefficients: array
    gamowDerivativeAtMinBarrierDepth: float
    _gamowDerivativeAtMaxEnergyDepth: float
    gamowDerivativePolynomial: array
    # endregion
    
    # region initialization
    def __init__(self, field:float = 5., radius:float = 1.e4, gamma:float = 10., \
        preloadedGamowTable:np.ndarray = None, preloadedLimits:np.ndarray = None, tabulationFolder = 'tabulated') -> None:

        super().__init__(preloadedGamowTable, preloadedLimits, tabulationFolder)
        self.setBarrierParameters(field, radius, gamma)
   
    def _calculateParameters(self) -> None:
        """Gets the parameters necessary to evaluate the Gamow factor. """
        
        data = self.interpolateForValues(self._field, self._radius, self._gamma)
        self.minEnergyDepth = data[-2] #Wmin: This is the top of the barrier
        self.maxEnergyDepth = data[-1] #Wmax: The maximum barrier depth for which the polynomial is valid
        self._gamowPolynomialCoefficients = data[:self._Npolynomial] #Gpoly
        self.gamowDerivativePolynomial = np.polyder(self._gamowPolynomialCoefficients) #dG/dW polynomial coefficients
        self.gamowDerivativeAtMinBarrierDepth = np.polyval(self.gamowDerivativePolynomial, self.minEnergyDepth) #dG/dW @ Wmin
        self._gamowDerivativeAtMaxEnergyDepth = np.polyval(self.gamowDerivativePolynomial, self.maxEnergyDepth) #dG/dW @ Wmax
    
    def _extendFieldBounds(self, energyDepth : array) -> array:
        """
        In case self._field is out of bounds of the Table, linearly extrapolates the value of Gamow for all given energies
        
        Parameters:
        energyDepth (array): input energy values for which extension for Gamow is to be provided

        Returns:
        array of Delta Gamow to be added to the values calculated at the edge
        """

        deltaInverseField = np.diff(self._inverseFieldLimits) / (self._Nfields - 1)

        if (1./self._field < self._inverseFieldLimits[0]):
            adgacentInverseFieldValue = self._inverseFieldLimits[0] + deltaInverseField
            newBarrier = Barrier(1./adgacentInverseFieldValue, self._radius, self._gamma, preloadedGamowTable=self._gamowTable, preloadedLimits=self._limits)
            gamowDerivative = (self.Gamow - newBarrier.calculateGamow(energyDepth)) / deltaInverseField
            return gamowDerivative * (1 / self._field - self._inverseFieldLimits[0])
        elif (1./self._field > self._inverseFieldLimits[1]):
            adgacentInverseFieldValue = self._inverseFieldLimits[1] - deltaInverseField
            newBarrier = Barrier(1./adgacentInverseFieldValue, self._radius, self._gamma, preloadedGamowTable=self._gamowTable, preloadedLimits=self._limits)
            gamowDerivative = (self.Gamow - newBarrier.calculateGamow(energyDepth)) / deltaInverseField
            return gamowDerivative * (1 / self._field - self._inverseFieldLimits[1])
        else:
            return 0. 

    def calculateGamow(self, energyDepth:array) -> array:
        """Evaluates the gammow coefficients in order to provide the exponential for the Kemble formula

        Parameters:
            energy (array): Incoming energy, perpendicular to the emitting surface, values (eV) for which the transmission cofficients are calculated (eV)

        Returns:
            array: Kemble's formula Gamow exponent value for a given energy (number)
        """
        self._calculateParameters()
        
        if (isinstance(energyDepth, np.ndarray)):
            self.Gamow = np.copy(energyDepth)
            highEnergyDepth = energyDepth > self.maxEnergyDepth
            lowEnergyDepth = energyDepth < self.minEnergyDepth
            self.Gamow = np.polyval(self._gamowPolynomialCoefficients, energyDepth)
            self.Gamow[highEnergyDepth] = np.polyval(self._gamowPolynomialCoefficients, self.maxEnergyDepth) + \
                self._gamowDerivativeAtMaxEnergyDepth * (energyDepth[highEnergyDepth] - self.maxEnergyDepth)
            self.Gamow[lowEnergyDepth] = np.polyval(self._gamowPolynomialCoefficients, self.minEnergyDepth) + \
                self.gamowDerivativeAtMinBarrierDepth * (energyDepth[lowEnergyDepth] - self.minEnergyDepth)
        else:
            if (energyDepth > self.minEnergyDepth and energyDepth < self.maxEnergyDepth):
                self.Gamow = np.polyval(self._gamowPolynomialCoefficients, energyDepth)
            elif(energyDepth > self.maxEnergyDepth):
                self.Gamow = np.polyval(self._gamowPolynomialCoefficients, self.maxEnergyDepth) + \
                    self._gamowDerivativeAtMaxEnergyDepth * (energyDepth - self.maxEnergyDepth)
            else:
                self.Gamow = np.polyval(self._gamowPolynomialCoefficients, self.minEnergyDepth) + \
                    self.gamowDerivativeAtMinBarrierDepth * (energyDepth - self.minEnergyDepth) 

        if (1./self._field < self._inverseFieldLimits[0] or 1./self._field > self._inverseFieldLimits[1]):
            self.Gamow += self._extendFieldBounds(energyDepth)
        
        return self.Gamow

    def setBarrierParameters(self, field:float, radius:float = 1.e4, gamma:float = 10.) -> None:
        """Sets the main parameters of the barrier.

        Parameters:
            field (float): Electric field (V/nm)
            radius (float): Emitter tip raidus (nm). Default value 1.e4 in case it is omitted for planar approximation
            gamma (float): gamma coefficient (number). Default value 10 in case it is omitted
        """
        
        self._field = field
        self._radius = radius
        self._gamma = gamma
    
    def transmissionCoefficient(self, energyDepth : array) -> array:
        """Calculates the probability of an electron tunelling (being transmitted) through the potential barrier by 
        means of the Kemble formula within the JWKB approximation

        Parameters:
            energyDepth (array): Energy depth of the barrier, i.e. incoming electron energy perpendicular to the emitting surface, 
            measured downwards from vacuum level. The Energy values (eV) for which the transmission cofficients are calculated

        Returns:
            array: Probability of an electron with a certain kinetic energy perpendicular to the emitting surface to tunnerl through the potential barrier (number)
        """
        
        gamow = self.calculateGamow(energyDepth)
        if (isinstance(energyDepth, np.ndarray)):
            transmissionCoefficient = np.copy(gamow)
            highGamowArea = gamow > 18.
            lowGamowArea = gamow <= 18.
            transmissionCoefficient[lowGamowArea] = 1 / (1 + np.exp(gamow[lowGamowArea]))
            transmissionCoefficient[highGamowArea] = np.exp(-gamow[highGamowArea])
            return transmissionCoefficient
        else:
            if (gamow < 18.):
                return 1 / (1 + np.exp(gamow))
            else:
                return np.exp(-gamow)
      
    def plotGamow(self, minBarrierDepth : float = 0., maxBarrierDepth : float = 0., show = False, saveFile = "") -> None:
        """Plots the Gamow exponent for energy depth between minBarrierDepth and maxBarrierDepth. IF saveFile is given,
        it will save the file to the given output. If show = True, it will show the output plot."""
        
        barrierDepths = np.linspace(minBarrierDepth, maxBarrierDepth, 256)
        gamowValues = self.calculateGamow(barrierDepths)
        plt.plot(barrierDepths, gamowValues)
        plt.xlabel("Barrier Depth [eV]")
        plt.ylabel("Gamow Factor")
        if (show):
            plt.show()
        
        if (saveFile):
            plt.savefig(saveFile)

    def plotTransmissionCoefficient(self, minBarrierDepth : float = 0., maxBarrierDepth : float = 0., show = False, saveFile = "") -> None:
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
   
    #the limit for which different approximations of the Fermi Dirac functions apply
    exponentLimit =  - 0.5 * np.log(np.finfo(float).eps)
    
    
    
    # endregion

    # region inisialisation
    def __init__(self):
        pass
    # endregion 
    
    # region user methods
    def logFermiDiracFunction(self, energy:array, kT:array):
        """Returs the natural logarith of the Fermi Diract distribution.
        The energy is divided in 3 terms to prevent the computer from overfloading with large exponents 

        Args:
            energy (array): Energy, perpendicular to the emitting surface, values (eV) for which the transmission cofficients are calculated (eV)
            kT (float): Energy from temperature (eV)

        Returns:
            array: Natural logarith of the Fermi Diract distribution (number)
        """
     
        if (isinstance(energy, np.ndarray)):
            functionValue = np.copy(energy)

            highEnergyArea = energy > self.exponentLimit * kT
            lowEnergyArea = energy < - self.exponentLimit * kT
            midEnergyArea = np.logical_not(highEnergyArea) * np.logical_not(lowEnergyArea)

            functionValue[highEnergyArea] = np.exp(-energy[highEnergyArea] / kT)
            functionValue[lowEnergyArea] = - energy[lowEnergyArea] / kT + np.exp(energy[lowEnergyArea] / kT)
            functionValue[midEnergyArea] = np.log(1. + np.exp(-energy[midEnergyArea] / kT))
            return functionValue
        else:
            if energy > self.exponentLimit * kT:
                return np.exp(-energy / kT)
            elif(energy < - self.exponentLimit * kT):
                return -energy / kT + np.exp(energy / kT)
            else:
                return np.log(1. + np.exp(-energy / kT))
    
    def FermiDiracFunction(self, energy:array, kT:array):
        """Returns the Fermi Dirac distribution
        
        Args:
            energy (array): Total Electron Energy (eV)
            kT (float): Energy from temperature (eV)

        Returns:
            array: Fermi Diract distribution (number)
        """
        
        if (isinstance(energy, np.ndarray)):
            functionValue = np.copy(energy)

            highEnergyArea = energy > self.exponentLimit * kT
            lowEnergyArea = energy < - self.exponentLimit * kT
            midEnergyArea = np.logical_not(highEnergyArea) * np.logical_not(lowEnergyArea)

            # use approximation 1/(1+exp(x)) = exp(-x) + O(exp(-2*x)) for x->+inf, i.e. exp(x) >> 1
            functionValue[highEnergyArea] = np.exp(-energy[highEnergyArea] / kT)
            
            #use approximation 1/(1+exp(x)) = 1 - exp(x) + O(exp(2*x)) for x->-inf, i.e. exp(x) << 1
            functionValue[lowEnergyArea] = 1. - np.exp(energy[lowEnergyArea] / kT)
            
            functionValue[midEnergyArea] = 1. / (1. + np.exp(energy[midEnergyArea] / kT))
            return functionValue
        else:
            if energy > self.exponentLimit * kT:
                return np.exp(-energy / kT)
            elif(energy < - self.exponentLimit * kT):
                return 1. - np.exp(energy / kT)
            else:
                return 1. / (1. + np.exp(energy / kT))
        
    def plotFermiDiracFunctions(self):
        """Plots the Fermi Dirac functions. Used for testing"""
        # x = np.linspace(-1.2 * self.exponentLimit, 1.2 * self.exponentLimit, 256)
        # x = np.linspace(-self.exponentLimit- .01, - self.exponentLimit + .01, 256)
        x = np.linspace(-5, 5, 256)
        fermiDirac = self.FermiDiracFunction(x, 1.)
        logFermiDirac = self.logFermiDiracFunction(x, 1.)
        
        figure, axes = plt.subplots(2,1, sharex=True)
        
        
        axes[0].plot(x, fermiDirac, label="Fermi-Dirac")
        axes[0].plot(x, logFermiDirac, label = "Log Fermi-Dirac")
        axes[0].set_ylabel("f(E/kT)")
        axes[0].grid()
        
        axes[1].semilogy(x, fermiDirac, label="Fermi-Dirac")
        axes[1].semilogy(x, logFermiDirac, label = "Log Fermi-Dirac")
        axes[1].set_xlabel("E/kT")
        axes[1].set_ylabel("f(E/kT)")
        axes[1].grid()
        plt.savefig("fermiDiracFunctions.png")


class Emitter:
    fastIntegrator = ct.CDLL(pythonpath + '/libintegrator.so') #use absolute path
    fastIntegrator.currentDensityPerNormalEnergy.restype = ct.c_double
    fastIntegrator.currentDensityPerNormalEnergy.argtypes = (ct.c_int, ct.c_double)
    fastIntegrator.nottinghamHeatInegrand.restype = ct.c_double
    fastIntegrator.nottinghamHeatInegrand.argtypes = (ct.c_int, ct.c_double)
    
    SommerfeldConstant:float = 1.618311e-4 
    # region initialization
    def __init__(self, barrier: Barrier, supply: Supply):
        """Initialises the class by inheritating Tabulator and setting it up as a private field to be used by Emitter

        Args:
            tabulator (Tabulator): Object that contains the attributes to interpolate the potential barrier given the field, radius and gamma as inputs
        """
        self.barrier = barrier
        self.supply = supply
    # endregion

   
class MetalEmitter(Emitter):
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
    kT: float
    workFunction: float
    energy: array
    _centerOfSpectrumEnergy: float
    _GamowDerivativeAtFermi: float
    
    # endregion
    
    # region initialization
    def __init__(self, barrier:Barrier, supply:Supply):
        """Initialises the function by adding the atributes of a "tabulated" barrier emitter to metals

        Args:
            tabulator (Tabulator): Object that contains the attributes of a x material emitter (potential barrier, transmission coefficiens and electron supply functions)
        """
        
        super().__init__(barrier, supply)
        
    def _calculateIntegrationLimits(self, decayCutoff = 10.):
        """Finds the limits of integration, based on the regimes described by Jensen's GTF theory (see http://dx.doi.org/10.1063/1.4940721 for details)
            The limits are saved on self._highEnergyLimit and self._lowEnergyLimit
        """
        self.barrier._calculateParameters()
        #get the maximum (for energies above Fermi level) derivative of Gamow. It is capped at Wmax
        self._GamowDerivativeAtFermi = np.polyval(self.barrier.gamowDerivativePolynomial, min(self.workFunction, self.barrier.maxEnergyDepth))
        
        if (self._GamowDerivativeAtFermi * self.kT < 0.002): # Very low temperature and quad() misbehaves. set some help
            self._highEnergyLimit = 0.
            
            #get the polynomial P for which the equation P(E)=0 finds the place where E * dG/dE = 1
            polynomialForRoots = np.append(self.barrier.gamowDerivativePolynomial, -1.) #get W*dG/dW
            polynomialForRoots[1:] -= self.workFunction * self.barrier.gamowDerivativePolynomial # add -phi * dG/dw
            realroots = np.roots(polynomialForRoots)
            self._centerOfSpectrumEnergy = self.workFunction - realroots[np.nonzero(np.logical_and(realroots > self.workFunction, realroots < self.barrier.maxEnergyDepth))][0]
            self._lowEnergyLimit = self._centerOfSpectrumEnergy - decayCutoff / self._GamowDerivativeAtFermi 
        elif (self._GamowDerivativeAtFermi * self.kT < 0.95): #field regime. The spectrum is centered around Fermi Level
            self._highEnergyLimit = decayCutoff / (1. / self.kT - self._GamowDerivativeAtFermi) # decayCutoff divided by decay rate
            self._lowEnergyLimit = - decayCutoff / self._GamowDerivativeAtFermi 
            return        
        elif (self.barrier.gamowDerivativeAtMinBarrierDepth * self.kT > 1.05): #thermal regime. The spctrum is centered around the top of the barrier
            self._centerOfSpectrumEnergy = self.workFunction - self.barrier.minEnergyDepth #top of the barrier (Um)
            self._highEnergyLimit = self._centerOfSpectrumEnergy + decayCutoff * self.kT 
            self._lowEnergyLimit = self._centerOfSpectrumEnergy - 10 / self.barrier.gamowDerivativeAtMinBarrierDepth  
            return     
        else: #intermediate regime. The spectrum center is approximately where dG/dE = 1/kT
            #get the polynomial P for which the equation P(w)=0 finds the place where dG/dw = 1/kT
            polynomialForRoots = np.copy(self.barrier.gamowDerivativePolynomial)
            polynomialForRoots[-1] -= 1./self.kT
            
            #find meaningful root of the polynomial
            realroots = np.roots(polynomialForRoots)
            # value of w where dG/dw = 1/kT
            self._centerOfSpectrumEnergy = self.workFunction - realroots[np.nonzero(np.logical_and(realroots > self.barrier.minEnergyDepth, realroots < self.workFunction))][0]
            
            # find energy limits
            self._highEnergyLimit =  self._centerOfSpectrumEnergy + decayCutoff * self.kT
            self._lowEnergyLimit = self._centerOfSpectrumEnergy - 2.5 * decayCutoff * self.kT
            return

    def setParameters(self, workfunction:float, kT:float):
        """Defines main emitter characteristics

        Args:
            workfunction (float): Average minimum energy required to put an electron, whitin a material, out in vacuum (some few nm away from the material surface) (eV)
            kT (float): Energy from temperature (eV)
        """
        
        self.workFunction = workfunction
        self.kT = kT
        self._calculateIntegrationLimits()
        
    def currentDensityPerNormalEnergy(self, energy, saveSpectrum = False):
        integrand = self.barrier.transmissionCoefficient(self.workFunction - energy) * self.supply.logFermiDiracFunction(energy, self.kT)
        if saveSpectrum:
            self._currentDensityPerEnergy.append(integrand)
            self._spectralEnergyPoints.append(energy)
        return integrand
    
    def currentDensityFast(self):
        """Calculates the field emitted current density from metal surfaces

        Returns:
            float: Emitted current density being emitted from an infinitesimal flat surface area (electrons/area*time)
        """
        
        args = tuple([self.workFunction, self.kT, self.barrier.minEnergyDepth, self.barrier.maxEnergyDepth, self.barrier.gamowDerivativeAtMinBarrierDepth, \
            self.barrier._gamowDerivativeAtMaxEnergyDepth] + list(self.barrier._gamowPolynomialCoefficients))
        try:
            integral, abserr = ig.quad(self.fastIntegrator.currentDensityPerNormalEnergy, self._lowEnergyLimit, self._highEnergyLimit, args)
        except(IntegrationWarning):
            integral = 0.
            
        return self.SommerfeldConstant * self.kT * integral
    
    def currentDensitySlow(self, saveSpectrum = False):
        """Calculates the field emitted current density from metal surfaces

        Returns:
            float: Emitted current density being emitted from an infinitesimal flat surface area (electrons/area*time)
        """
        
        try:
            if (self._highEnergyLimit == 0.):
                breakPoints =  np.array([0., 0.5 * self._centerOfSpectrumEnergy, self._centerOfSpectrumEnergy, 0.75 * self._centerOfSpectrumEnergy + 0.25 * self._lowEnergyLimit ])
                integral, abserr, infodict= ig.quad(self.currentDensityPerNormalEnergy, self._lowEnergyLimit, self._highEnergyLimit, args=(saveSpectrum), full_output=True, points=breakPoints)
                remainder = self.kT * self.barrier.transmissionCoefficient(self.workFunction) * (0.8224670334241132 + 0.9015426773696957 * self._GamowDerivativeAtFermi)
                integral += remainder 
                
            else:
                integral, abserr= ig.quad(self.currentDensityPerNormalEnergy, self._lowEnergyLimit, self._highEnergyLimit, args=(saveSpectrum))
                
        except(IntegrationWarning):
            integral = 0.
            
        return self.SommerfeldConstant * self.kT * integral
    
    def Energy_Distribution_for_Metals(self):
        """Calculates the energy distribution of field emitted electrons from metals

        Returns:
            energy (array): Energy space for which the electron emission has been evalated (eV)
            electron_number (array): Number of electrons being emitter with a certain energy (number)
        """
        
        transmission_in_energy_z = self.barrier.transmissionCoefficient(self.workFunction - self.energy)

        cumulative_transmission = np.insert(np.cumsum((transmission_in_energy_z[1:]+transmission_in_energy_z[:-1])*(self.energy[1]-self.energy[0])/2), 0, 0)
        
        electron_number = self.supply.FermiDiracFunction(self.energy, self.kT) * cumulative_transmission

        return  self.energy, electron_number
    
    def normalEnergyDistribution(self):
        """Calculates the perpendicular energy distribution of field emitted electrons from metals

        Returns:
            energy (array): Energy space for which the electron emission has been evalated (eV)
            electron_number (array): Number of electrons being emitter with a certain energy (number)
        """
        self._spectralEnergyPoints = []
        self._currentDensityPerEnergy = []
        
        currentDensity = self.currentDensitySlow(saveSpectrum=True)
        
        spectralPointArray = np.array(self._spectralEnergyPoints)
        sortIndices = np.argsort(spectralPointArray)
    
        return  spectralPointArray[sortIndices], np.array(self._currentDensityPerEnergy)[sortIndices] * self.kT * self.SommerfeldConstant
    
    def TEDdifferentialSystem(self, energy, integrand):
        return self.supply.FermiDiracFunction(energy, self.kT) * \
            (self.barrier.transmissionCoefficient(self.workFunction - energy) - integrand * np.exp(energy/self.kT) / self.kT)
            
    def TEDJacobian(self, energy, integrand):
        return  np.array([[-self.supply.FermiDiracFunction(energy, self.kT) * np.exp(energy/self.kT) / self.kT]])
    
    def totalEnergyDistribution(self, refinementLevel:int = 2):
        """Calculates the total energy distribution of field emitted electrons from metals. 
        The methodology
        
        Returns:
            energy (array): Energy space for which the electron emission has been evalated (eV)
            electron_number (array): Number of electrons being emitter with a certain energy (number)
        """
        solution = ig.solve_ivp(self.TEDdifferentialSystem, [self._lowEnergyLimit, self._highEnergyLimit], \
            [1.e-8], method='LSODA', dense_output=True, rtol = 1.e-8, jac=self.TEDJacobian)
        energypoints = solution.t
        for i in range(refinementLevel):
            energypoints = np.concatenate((energypoints, 0.5 * (solution.t[1:] + solution.t[:-1])))
        energypoints = np.sort(energypoints)
        return energypoints, solution.sol(energypoints)[0]
        
  
    def Nottingham_Heat_from_Metals(self):
        """Calculates the Nottingham heat resulted from field emitted electrons from metals

        Returns:
            float: Resulted Nottigham from the emission of electrons from an infinitesinal area (J)
        """

        args = tuple([self.workFunction] + [self.kT] + [self.barrier.minEnergyDepth] + [self.barrier.maxEnergyDepth] + [self.barrier.gamowDerivativeAtMinBarrierDepth] + [self.barrier._gamowDerivativeAtMaxEnergyDepth] + list(self.barrier._gamowPolynomialCoefficients))
        try:
            integ, abserr = ig.quad(self.fastIntegrator.nottinghamHeatInegrand, self.energy[0], self.energy[-1], args, full_output = 0)
        except:
            return 0.
        
        return -self.SommerfeldConstant * self.kT * integ    
    
    def Current_Density_from_Metals_educational(self):
        """Calculates the field emitted current density from metals on a clear equation to code manner"""
        integ = self.supply.logFermiDiracFunction(self.energy, self.kT) * self.Transmission_Coefficient(self.workFunction - self.energy)
        
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
        return self.SommerfeldConstant * self.kT * np.sum(integ) * (self.energy[1] - self.energy[0])
    
    def Energy_Distribution_of_Metals_educational(self):
        """Calculates the energy distribution of field emitted electrons from metals on a clear equation to code manner"""
        transmission_in_energy_z = self.Transmission_Coefficient(self.workFunction - self.energy)

        cumulative_transmission = np.copy(self.energy)
        cumulative_transmission[0] = 0

        for i in range(len(self.energy)-1):
            cumulative_transmission[i+1] = cumulative_transmission[i]+((self.energy[1]-self.energy[0])*(transmission_in_energy_z[i+1]+transmission_in_energy_z[i])/2)
        
        energy_distribution = self.supply.FermiDiracFunction(self.energy, self.kT) * cumulative_transmission

        return  self.energy, energy_distribution
    
    def Nottingham_Heat_from_Metals_educational(self):
        """Calculates the Nottingham heat resulted from field emitted electrons from metals on a clear equation to code manner"""
        transmission_in_energy_z = self.Transmission_Coefficient(self.workFunction - self.energy)

        cumulative_transmission = np.copy(self.energy)
        cumulative_transmission[0] = 0

        for i in range(len(self.energy)-1):
            cumulative_transmission[i+1] = cumulative_transmission[i]+((self.energy[1]-self.energy[0])*(transmission_in_energy_z[i+1]+transmission_in_energy_z[i])/2)
        
        heat_components = self.supply.FermiDiracFunction(self.energy, self.kT) * cumulative_transmission * self.energy
        
        return -self.SommerfeldConstant * np.sum(heat_components) * (self.energy[1] - self.energy[0])
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
        resolution = 128
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
        
        jc_integ = self.supply.logFermiDiracFunction(self._energy_conduction_band, self._kT) * (self.barrier.transmissionCoefficient(-self._Ef-self._energy_conduction_band) - (b*self.barrier.transmissionCoefficient(-self._Ef-a*self._energy_conduction_band)))
        
        return self.SommerfeldConstant * self._kT * np.sum(jc_integ) * (self._energy_conduction_band[1]-self._energy_conduction_band[0]) 
    
    def Current_from_Valence_Band(self):
        """Calculates the conduction band component of the emitted current - Andreas' equation

        Returns:
            float: Emitted current density from the valence band (electrons/area*time)
        """
        
        a = self._mp/self._m
        b = 1+a
        
        jv_integ = self.supply.logFermiDiracFunction(self._energy_valence_band, self._kT) * (self.barrier.transmissionCoefficient(-self._Ef-self._energy_valence_band) - (b*self.barrier.transmissionCoefficient(-self._Ef-b*self._energy_valence_band+a*self._Evhigh)))
        
        return self.SommerfeldConstant * self._kT * ((np.sum(jv_integ) * (self._energy_valence_band[1]-self._energy_valence_band[0])))
    
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
        
        c_number_of_electrons = np.sum(self.supply.FermiDiracFunction(self._energy_conduction_band, self._kT) * c_cumulative_transmission )
        
        total_c = np.sum(self.SommerfeldConstant * c_number_of_electrons * (self._energy_conduction_band[1]-self._energy_conduction_band[0]))
        
        #valence band
        c = self._mp/self._m
        d = 1+c
        
        v_transmission_in_z = self.Transmission_Coefficient(-self._Ef-self._energy_valence_band) - (d*self.Transmission_Coefficient(-self._Ef-d*self._energy_valence_band+c*self._Evhigh))
        
        v_cumulative_transmission = np.insert(np.cumsum((v_transmission_in_z[1:]+v_transmission_in_z[:-1])*(self._energy_valence_band[1]-self._energy_valence_band[0])/2), 0, 0)
            
        v_number_of_electrons = self.supply.FermiDiracFunction(self._energy_valence_band, self._kT) * v_cumulative_transmission 
        
        total_v = np.sum(self.SommerfeldConstant * v_number_of_electrons * (self._energy_valence_band[1]-self._energy_valence_band[0]))
        
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
        
        number_of_electrons = self.supply.FermiDiracFunction(self._energy_conduction_band, self._kT) * cumulative_transmission 

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
            
        number_of_electrons = self.supply.FermiDiracFunction(self._energy_valence_band, self._kT) * cumulative_transmission 
      
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
        
        return -self.SommerfeldConstant* self._kT * np.sum(heat*distribution)
    
    def Nottingham_Heat_from_Valence_Band(self): #NEED MOFICATIONS
        """Calculates the contribution of the valence band to the Nottingham heat resulted from field emitted electrons from semicoductors"""
        """the fuction will be replaced by the full code from the function, if we agree the function works well"""
        energy, distribution = self.Energy_Distribution_from_Valence_Band()#*self._energy_valence_band

        heat = energy

        return -self.SommerfeldConstant * self._kT * np.sum(heat*distribution)

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
        
        energy_distribution = self.supply.FermiDiracFunction(self._energy_conduction_band, self._kT) * cumulative_transmission 

        return self._energy_conduction_band , energy_distribution
    # endregiondata


def current_metal_emitter(Field:array, Radius:array, Gamma:array, Workfunction:array, Temperature:array):

    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature
    
    metal_emitter = MetalEmitter(Barrier(tabulationFolder="tabulated/2D_256x128"),Supply())

    j_metal = np.copy(Field)
    
    for i in range(len(Field)):

        metal_emitter.barrier.setBarrierParameters(Field[i], Radius[i], Gamma[i])
        metal_emitter.barrier._calculateParameters()
        metal_emitter.setParameters(Workfunction[i], kT[i])
        j_metal[i] = metal_emitter.currentDensityFast()
        
    return j_metal

def heat_metal_emitter(Field:array, Radius:array, Gamma:array, Workfunction:array, Temperature:array):
    
    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature
    
    metal_emitter = MetalEmitter(Barrier(),Supply())
    
    nh_metal = np.copy(Field)

    for i in range(len(Field)):
        metal_emitter.barrier.setBarrierParameters(Field[i], Radius[i], Gamma[i])
        metal_emitter.barrier._calculateParameters()
    
        metal_emitter.setParameters(Workfunction[i], kT[i])
    
        nh_metal[i] = metal_emitter.Nottingham_Heat_from_Metals()

    return nh_metal

def heat_metal_emitter2(Field:array, Radius:array, Gamma:array, Workfunction:array, Temperature:array):

    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature
    
    metal_emitter = MetalEmitter(Barrier(),Supply())
    
    nh_metal = np.copy(Field)

    for i in range(len(Field)):
        metal_emitter.barrier.setBarrierParameters(Field[i], Radius[i], Gamma[i])
        metal_emitter.barrier._calculateParameters()
    
        metal_emitter.setParameters(Workfunction[i], kT[i])
    
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
    # print("test for normal barrier:")
    # bar = Barrier(5, 10, 10., tabulationFolder="tabulated/2D_512x256")
    # print("Gamow(E=4.5) = ", bar.calculateGamow(4.5))
    # bar.plotGamow(0., 10., show = True)

    # print("test for for field over the tabulation limits")
    # bar = Barrier(25, 10, 10., preloadedGamowTable=bar._gamowTable, preloadedLimits=bar._limits)
    # print("Gamow(E=4.5) = ", bar.calculateGamow(4.5))
    # bar.plotGamow(0., 10., show = True)

    # print("test for for field below the tabulation limits")
    # bar = Barrier(0.1, 10, 10., preloadedGamowTable=bar._gamowTable, preloadedLimits=bar._limits)
    # print("Gamow(E=4.5) = ", bar.calculateGamow(4.5))
    # bar.plotGamow(0., 10., show = True)

    # print(current_metal_emitter(np.array([5., 6.]), np.array([10.,12.]), np.array([10.,10.]), np.array([4.5, 4.5]), np.array([396., 159])))
    
    bar = Barrier(5, 10, 10., tabulationFolder="tabulated")
    sup = Supply()
    em = MetalEmitter(bar, sup)
    
    kT = BoltzmannConstant * 3000.
    # kT = 1.
    for i in range(8):
        kT *= 0.5
        em.setParameters(4., kT)
        Energy, electronCount = em.totalEnergyDistribution()
        print(em.currentDensityFast(), em.currentDensitySlow())


        plt.semilogy(Energy, electronCount, '.', markersize = 1.5)
    plt.grid()
    # plt.xlim(-2, 4)
    plt.savefig("spectrum.png")




