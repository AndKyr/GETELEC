import getelec as gt
import numpy as np
import copy
import scipy.optimize as opt
import matplotlib.pyplot as plt

class IVDataFitter:
    
    emitterModel:gt.GETELECModel
    fittingParameters: dict
    initialParameters: dict
    minParameters: dict
    maxParameters: dict
    voltageData:np.ndarray
    currentData:np.ndarray

    def __init__(self, emitterModel:gt.GETELECModel = gt.GETELECModel()) -> None:
        self.emitterModel = emitterModel
        self.fittingParameters = {"fieldConversionFactor": 1.}
        self.minParameters = {"fieldConversionFactor": 1.}
        self.initialParameters = {"fieldConversionFactor": 1.}
        self.maxParameters = {"fieldConversionFactor": 1.}


    def currentDensityforVoltages(self, voltageArray:np.ndarray) -> np.ndarray:
        """ Calculates and returns the current density for a given array of voltages, assuming a field conversion factor 
        (F=fieldConversionfactor * Voltage). 
        
        Parameters:
            voltageArray: array of Voltages (Volts)
        Returns:
            currentDensity: Array of the resulting current densities for each voltage element
        """

        #set the fields
        self.emitterModel.setParameters(field=self.fittingParameters["fieldConversionFactor"] * voltageArray)
        
        #set all other parameters that need to be fitted
        self.emitterModel.setParameters(**{key: self.fittingParameters[key] for key in self.fittingParameters if key !="fieldConversionFactor"})
    
        self.emitterModel.calculateCurrentDensity()
        
        return np.array(self.emitterModel.getCurrentDensity())
    
    def logCurrentDensityError(self, parameterList) -> float:
        """ Calculates the error between a given I-V curve and the one calculated by GETELEC for a given set of parameters.
        The main usage of this function is to be called by external minimization function in order to find the parameterList that
        gives the best fit to the given I-V

        Parameters:
            parameterList: list of emitter parameters (fieldConversionFactor, radius, gamma, workFunction, temperature)
            voltageArray: input array of voltages (measured)
            curentDensityArray: input array of current densities
        Returns:
            array of relative errors (in logarithmic scale) between the calculated and the given current densities, after normalizing with a multiplication
            factor that forces them to match at their maximum
        """

        assert len(parameterList) == len(self.fittingParameters), "parameterList has a length of %d while it should have %d"%(len(parameterList), len(self.fittingParameters))
        i = 0
        for paramName in self.fittingParameters:
            self.fittingParameters[paramName] = parameterList[i]
            i += 1


        calculatedCurrentDensities = self.currentDensityforVoltages(self.voltageData)
        
        #normalize by forcing the max values to match
        maxValueRatio = np.max(self.currentData) / np.max(calculatedCurrentDensities)
        error = np.log(maxValueRatio * calculatedCurrentDensities / self.currentData)
        return error

    def setIVdata(self, voltageData:np.ndarray, currentData:np.ndarray) -> None:
        self.voltageData = voltageData
        self.currentData = currentData

    def setParameterRange(self, field:list = [1., 6., 15.], radius = 2000., gamma = 10., workFunction = 4.5, temperature = 300.) ->None:

        """Sets the range of the parameters that will be fitted (variables) and the value of those that will be fixed
        
        Args: 
            field: list of 3 floats. Minimum, initial, and maximum value of the field that is allowed for fitting
            radius, gamma, workFunction, temperature: float or list of 3 floats. If single float, it means that the parameter is considered fixed. If list, it is variable and it works as for field. Default values are 2000. [nm], 10., 4.5 [eV], 300 [K] correspondingly
        """
        self.minParameters ["fieldConversionFactor"] = field[0] / np.min(self.voltageData)
        self.initialParameters["fieldConversionFactor"] = field[1] / np.mean(self.voltageData)
        self.maxParameters["fieldConversionFactor"] = field[2] / np.max(self.voltageData)

        #parameter arguments that are to be fixed
        fixedParameterArgs = {key:value for key,value in locals().items() if type(value) in [float, int]}
        
        #parameter arguments that are to be 
        fittingParameterArgs = {key:value for key,value in locals().items() if type(value) == list and key != "field"}

        # set range for all parameters to be fitted
        for key,value in fittingParameterArgs.items(): 
            self.fittingParameters[key] = value[1]
            self.minParameters[key] = value[0]
            self.initialParameters[key] = value[1]
            self.maxParameters[key] = value[2]
        
        # set the value of the parameters that are fixed
        self.emitterModel.setParameters(**fixedParameterArgs)

        #remove all fixed parameters from the fittingParameter dictionary
        for key in fixedParameterArgs:
            try:
                self.fittingParameters.pop(key)
            except KeyError:
                pass

        #sanity checks
        if "radius" in self.fittingParameters:
            assert self.emitterModel.emitter.barrier.getDimension() >= 2, "you cannot fit for Radius of Gamma if the tabulation dimension is < 2"

        if "gamma" in self.fittingParameters:
            assert self.emitterModel.emitter.barrier.getDimension() == 3, "You cannot fit for Gamma if the tabulation dimension is < 3"

    def fitIVCurve(self) -> None:
        """Performs the fitting"""
        self.optimizationData = opt.least_squares(fun = self.logCurrentDensityError, x0 = list(self.initialParameters.values()), \
            bounds =(list(self.minParameters.values()), list(self.maxParameters.values())), method = 'trf', jac = '3-point')
        
        i = 0
        for parameterName in self.fittingParameters:
            self.fittingParameters[parameterName] = self.optimizationData.x[i]
            i += 1
        
        unshiftedCurrentCurve = self.currentDensityforVoltages(self.voltageData)
        self.prefactor = np.max(self.currentData) / np.max(unshiftedCurrentCurve)
        self.fittedCurrent = self.prefactor * unshiftedCurrentCurve

    def getOptCurrentCurve(self, voltageData:np.ndarray = None) -> np.ndarray:
        """Calculates and returns the current curve for a given input voltage data array
        Args:
            voltageData[optional]. array of input voltage values [V]. If not given, the internal voltageData array is used
        Returns:
            Array of currents calculated by the optimized system 
        """
        if (voltageData is None):
            return self.fittedCurrent
        else:
            return self.prefactor * self.currentDensityforVoltages(voltageData)



if (__name__ == "__main__"): #some testing operations
    

    fitter = IVDataFitter()


    voltageData = np.array([6666.66666667, 6528.49740933, 6395.93908629, 6268.65671642, \
       6146.34146341, 6028.70813397, 5915.49295775, 5806.4516129, 5701.35746606, 5600., 5502.18340611, 5407.72532189,\
       5316.4556962 , 5228.21576763, 5142.85714286, 5060.24096386,\
       4980.23715415, 4902.72373541, 4827.5862069 , 4754.71698113,\
       4684.01486989, 4615.38461538, 4548.73646209, 4483.98576512,\
       4421.05263158, 4359.8615917 , 4300.34129693, 4242.42424242,\
       4186.04651163, 4131.14754098, 4077.66990291, 4025.55910543,\
       3974.76340694, 3925.23364486, 3876.92307692, 3829.78723404,\
       3783.78378378, 3738.87240356, 3695.01466276, 3652.17391304,\
       3610.31518625, 3569.40509915, 3529.41176471, 3490.30470914,\
       3452.05479452, 3414.63414634, 3378.01608579, 3342.17506631,\
       3307.08661417, 3272.72727273, 3239.07455013, 3206.10687023,\
       3173.80352645, 3142.1446384 , 3111.11111111, 3080.68459658,\
       3050.84745763, 3021.58273381, 2992.87410926, 2964.70588235,\
       2937.06293706, 2909.93071594, 2883.29519451, 2857.14285714])
    
    currentData = np.array([5.65354016e-05, 4.28464283e-05, 3.24541251e-05, 2.45683389e-05, \
       1.85875533e-05, 1.40539879e-05, 1.06193664e-05, 8.01882052e-06,\
       6.05098032e-06, 4.56287757e-06, 3.43825950e-06, 2.58892706e-06,\
       1.94793392e-06, 1.46451644e-06, 1.10020423e-06, 8.25856631e-07,\
       6.19415049e-07, 4.64192603e-07, 3.47574447e-07, 2.60030670e-07,\
       1.94366991e-07, 1.45156296e-07, 1.08307623e-07, 8.07395611e-08,\
       6.01329743e-08, 4.47437594e-08, 3.32614033e-08, 2.47019853e-08,\
       1.83274191e-08, 1.35845100e-08, 1.00590041e-08, 7.44096972e-09,\
       5.49873640e-09, 4.05929377e-09, 2.99355710e-09, 2.20530977e-09,\
       1.62290148e-09, 1.19302882e-09, 8.76074733e-10, 6.42626597e-10,\
       4.70867981e-10, 3.44633799e-10, 2.51960496e-10, 1.83999483e-10,\
       1.34216889e-10, 9.77911792e-11, 7.11688901e-11, 5.17338317e-11,\
       3.75620267e-11, 2.72401384e-11, 1.97311118e-11, 1.42748585e-11,\
       1.03149269e-11, 7.44441693e-12, 5.36613434e-12, 3.86326821e-12,\
       2.77783297e-12, 1.99486455e-12, 1.43077440e-12, 1.02488643e-12,\
       7.33204336e-13, 5.23857148e-13, 3.73798141e-13, 2.66375053e-13])


    fitter.setIVdata(voltageData=voltageData, currentData=currentData)
    fitter.setParameterRange()
    fitter.fitIVCurve()
    fittedCurrent = fitter.getOptCurrentCurve(voltageData)
    plt.semilogy(1./voltageData, currentData, '.', label="data" )
    plt.semilogy(1./voltageData, fittedCurrent, label = "beta = %g"%fitter.fittingParameters["fieldConversionFactor"])

    fitter.setParameterRange(radius=[1., 10., 2000.])

    fitter.fitIVCurve()
    fittedCurrent = fitter.getOptCurrentCurve()
    plt.semilogy(1./voltageData, fittedCurrent, label = "beta = %g, R = %g"%(fitter.fittingParameters["fieldConversionFactor"], fitter.fittingParameters["radius"]))

    plt.legend()

    plt.savefig("fittedcurve.png")