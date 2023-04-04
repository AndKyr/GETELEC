import getelec as gt
import numpy as np
import copy
import scipy.optimize as opt
import matplotlib.pyplot as plt

class IVDataFitter:
    
    parameters: dict
    initialParameters: dict
    minParameters: dict
    maxParameters: dict

    def __init__(self, emitterModel:gt.GETELECModel = gt.GETELECModel(), fitRadius:bool = False, fitGamma: bool = False ) -> None:
        self.emitterModel = emitterModel
        self.fitRadius = fitRadius
        self.fitGamma = fitGamma

        self.parameters = {"fieldConversionFactor": 1., "workFunction": 4.5, "temperature": 300.}

        if fitRadius:
            assert self.emitterModel.emitter.barrier.getDimension() >= 2, "you cannot fit for Radius of Gamma if the tabulation dimension is < 2"

        if (fitGamma):
            assert self.emitterModel.emitter.barrier.getDimension() == 3 and fitRadius, "You cannot fit for Gamma if the tabulation dimension is < 3 or without fitting for Radius as well"

        if self.fitRadius:
            self.parameters["radius"] = 10.
        
        if self.fitGamma:
            self.parameters["gamma"] = 10.

    

    def currentDensityforVoltages(self, voltageArray:np.ndarray ):
        """ Calculates and returns the current density for a given array of voltages, assuming a field conversion factor 
        (F=fieldConversionfactor * Voltage). 
        
        Parameters:
            voltageArray: array of Voltages (Volts)
            fieldConversionFactor: ratio between local field and voltage (1/nm). Default 1.
            radius: radius of curvature of the emitter. Default 20.
            gamma: gamma barrier parameter. Default 10.
            workFunction: Work function of the emitter (eV). Default 4.5
            temperature: local temperature of the emitter (K). Default 300.

        Returns:
            currentDensity: Array of the resulting current densities for each voltage element
        """

        self.emitterModel.setParameters(field=self.parameters["fieldConversionFactor"] * voltageArray, temperature=self.parameters["temperature"], workFunction=self.parameters["workFunction"])

        if (self.fitRadius):
            self.emitterModel.setParameters(radius=self.parameters["radius"])
        if (self.fitGamma):
            self.emitterModel.setParameters(gamma=self.parameters["Gamma"])
        
    
        self.emitterModel.calculateCurrentDensity()
        
        return np.array(self.emitterModel.getCurrentDensity())
    
    def logCurrentDensityError(self, parameterList, voltageData, currentDensityData):
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

        assert len(parameterList) == len(self.parameters), "parameterList has a length of %d while it should have %d"%(len(parameterList), len(self.parameters))
        i = 0
        for paramName in self.parameters:
            self.parameters[paramName] = parameterList[i]
            i += 1


        calculatedCurrentDensities = self.currentDensityforVoltages(voltageData)
        
        #normalize by forcing the max values to match
        maxValueRatio = np.max(currentDensityData) / np.max(calculatedCurrentDensities)
        error = np.log(maxValueRatio * calculatedCurrentDensities / currentDensityData)
        return error

    def setIVcurve(self, voltageData:np.ndarray, currentData:np.ndarray):
        self.voltageData = voltageData
        self.currentData = currentData
        return

    def setParameterRange(self, fieldRange = [1., 6., 12.], radiusRange = [1., 20., 1000.], gammaRange = [1., 10., 1000.], \
        workFunctionRange = [1., 4.5, 7.5], temperatureRange = [30., 300., 5000.]):
        self.minParameters = copy.deepcopy(self.parameters)
        self.maxParameters = copy.deepcopy(self.parameters)
        self.initialParameters = copy.deepcopy(self.parameters)

        self.minParameters["fieldConversionFactor"] = fieldRange[0] / np.min(self.voltageData)
        self.initialParameters["fieldConversionFactor"] = fieldRange[1] / np.mean(self.voltageData)
        self.maxParameters["fieldConversionFactor"] = fieldRange[2] / np.max(self.voltageData)
        
        self.minParameters["workFunction"] = workFunctionRange[0]
        self.initialParameters["workFunction"] = workFunctionRange[1]
        self.maxParameters["workFunction"] = workFunctionRange[2]

        self.minParameters["temperature"] = temperatureRange[0]
        self.initialParameters["temperature"] = temperatureRange[1]
        self.maxParameters["temperature"] = temperatureRange[2]

        if (self.fitRadius):
            self.minParameters["radius"] = radiusRange[0]
            self.initialParameters["radius"] = radiusRange[1]
            self.maxParameters["radius"] = radiusRange[2]
        if (self.fitGamma == 3):
            self.minParameters["gamma"] = gammaRange[0]
            self.initialParameters["gamma"] = gammaRange[1]
            self.maxParameters["gamma"] = gammaRange[2]

    def fitIVCurve(self) -> None:
        self.optimizationData = opt.least_squares(fun = self.logCurrentDensityError, args = (self.voltageData, self.currentData), x0 = list(self.initialParameters.values()), \
            bounds =(list(self.minParameters.values()), list(self.maxParameters.values())), method = 'trf', jac = '3-point')
        
        i = 0
        for parameterName in self.parameters:
            self.parameters[parameterName] = self.optimizationData.x[i]
            i += 1
        
        unshiftedCurrentCurve = self.currentDensityforVoltages(self.voltageData)
        self.prefactor = np.max(self.currentData) / np.max(unshiftedCurrentCurve)
        self.fittedCurrent = self.prefactor * unshiftedCurrentCurve

    def getOptCurrentCurve(self, voltageData = None) -> np.ndarray:
        if (len(voltageData) == 0):
            return self.fittedCurrent
        else:
            return self.prefactor * self.currentDensityforVoltages(voltageData)



if (__name__ == "__main__"): #some testing operations
    

    fitter = IVDataFitter()




    xFN = np.linspace(0.15, 0.3, 32)

    voltage = 1./xFN
    
    currentDensity = fitter.currentDensityforVoltages(voltage)
    


    fitter.setIVcurve(voltageData=voltage, currentData=currentDensity)
    fitter.setParameterRange()
    fitter.fitIVCurve()
    fittedCurrent = fitter.getOptCurrentCurve(voltage)
    plt.semilogy(voltage, currentDensity, '.')
    plt.semilogy(voltage, fittedCurrent)
    plt.savefig("fittedcurve.png")