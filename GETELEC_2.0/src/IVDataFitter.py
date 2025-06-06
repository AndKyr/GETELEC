import getelec as gt
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt


class generalFitter:
    emitterModel:gt.GETELECModel
    fittingParameters: dict
    initialParameters: dict
    minParameters: dict
    maxParameters: dict
    voltageData:np.ndarray
    currentData:np.ndarray
    optimizationData:opt.OptimizeResult

    def __init__(self, emitterModel:gt.GETELECModel = gt.GETELECModel()) -> None:
        self.emitterModel = emitterModel
        self.fittingParameters = {"fieldConversionFactor": 1.}
        self.minParameters = {"fieldConversionFactor": 1.}
        self.initialParameters = {"fieldConversionFactor": 1.}
        self.maxParameters = {"fieldConversionFactor": 1.}
        self.preFactor = 1.
     
    def setIVdata(self, voltageData:np.ndarray, currentData:np.ndarray) -> None:
        self.voltageData = voltageData
        self.currentData = currentData

    def getFittingError(self) -> None:
        return self.optimizationData.cost

class IVDataFitter(generalFitter):
    
    def __init__(self, emitterModel:gt.GETELECModel = gt.GETELECModel()) -> None:
        super().__init__(emitterModel)

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
        self.emitterModel.setParameters(**{key: self.fittingParameters[key] for key in self.fittingParameters if key in ["radius", "gamma", "temperature", "workFunction"]})
    
        self.emitterModel.calculateCurrentDensity()
                
        return np.array(self.emitterModel.getCurrentDensity())

    
    def logCurrentDensityError(self, parameterList) -> float:
        """ Calculates the error between a given I-V curve and the one calculated by GETELEC for a given set of parameters.
        The main usage of this function is to be called by external minimization function in order to find the parameterList that
        gives the best fit to the given I-V

        Parameters:
            parameterList: list of emitter parameters (fieldConversionFactor, radius, gamma, workFunction, temperature)
        Returns:
            array of relative errors (in logarithmic scale) between the calculated and the given current densities, after normalizing with a multiplication
            factor that forces them to match at their maximum
        """

        assert len(parameterList) == len(self.fittingParameters), "parameterList has a length of %d while it should have %d"%(len(parameterList), len(self.fittingParameters))
        
        #copy parameter list into the parameter dictionary
        i = 0
        for paramName in self.fittingParameters:
            self.fittingParameters[paramName] = parameterList[i]
            i += 1

        self.calculatedCurrentDensities = self.currentDensityforVoltages(self.voltageData)
        
        #normalize by forcing the max values to match
        logOfPrefactor = np.mean(np.log(self.currentData)) - np.mean(np.log(self.calculatedCurrentDensities))
        self.preFactor = max(np.exp(logOfPrefactor), self.prefactorBounds[0])
        self.prefactor = min(self.preFactor, self.prefactorBounds[1])
        error = np.log(self.calculatedCurrentDensities / self.currentData) + logOfPrefactor
        return error

    
    def fitIVCurve(self) -> None:
        """Performs the fitting"""
        self.optimizationData = opt.least_squares(fun = self.logCurrentDensityError, x0 = list(self.initialParameters.values()), \
            bounds =(list(self.minParameters.values()), list(self.maxParameters.values())), method = 'trf', jac = '3-point', xtol=1.e-12)
        
        i = 0
        for parameterName in self.fittingParameters:
            self.fittingParameters[parameterName] = self.optimizationData.x[i]
            i += 1
        
        self.fittedCurrent = self.preFactor * self.currentDensityforVoltages(self.voltageData)

    def setParameterRange(self, field:list = [1., 6., 15.], radius = 2000., gamma = 10., workFunction = 4.5, temperature = 300., prefactorBounds = None) ->None:

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
        
        #parameter arguments that are to be fitted
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
                self.initialParameters.pop(key)
                self.minParameters.pop(key)
                self.maxParameters.pop(key)
            except KeyError:
                pass

        if prefactorBounds is None:
            self.prefactorBounds = (np.finfo(np.float64).tiny, np.finfo(np.float64).max)
        else:
            self.prefactorBounds = prefactorBounds


        #sanity checks
        if "radius" in self.fittingParameters:
            assert self.emitterModel.emitter.barrier.getDimension() >= 2, "you cannot fit for Radius of Gamma if the tabulation dimension is < 2"

        if "gamma" in self.fittingParameters:
            assert self.emitterModel.emitter.barrier.getDimension() == 3, "You cannot fit for Gamma if the tabulation dimension is < 3"

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
            return self.preFactor * self.currentDensityforVoltages(voltageData)
    
    def printFittingData(self) -> None:
        print("Optimum parameters: ", self.fittingParameters)
        print("sigma * Aeff: ", self.preFactor)
        print("Fitting Error (relative): ", self.getFittingError())
        print("orthodoxy status: ", self.orthodoxyStatus())

    def orthodoxyStatus(self) -> str:
        maxCurrentDensity = np.max(self.calculatedCurrentDensities)
        if maxCurrentDensity > 1.e-5:
            self.orthodoxyMessage = """
            The local field and current density predicted by the fit is too high. Your emitter would blow up if such a current density were true. The input data do not correspond to orthodox field emission"""
            return "unorthodox (high)"
        if maxCurrentDensity > 1.e-7:
            self.orthodoxyMessage = """
            The local field and current density predicted by the fit quite high. Your emitter would be close to blowing up if such a current density were true. The input are likely to not correspond to orthodox field emission"""
            return "apparently orthodox (high)"
        
        if self.preFactor > 1.e10:
            self.orthodoxyMessage = """The local field and current density predicted by the fit are too low and the corresponding effective emission area is too high. Do you think it is reasonable your effective emission area (total current / max current density) to be in the order of %g nm^2? If not (e.g. you are not measuring a huge amount of emitters or your emitter is not 0.1mm in radius), the input data do not correspond to orthodox field emission"""%self.preFactor
            return "unorthodox (low)"
        if self.preFactor > 1.e5:
            self.orthodoxyMessage = """The local field and current density predicted by the fit are quite low and the corresponding effective emission area is pretty high. Do you think it is reasonable your effective emission area (total current / max current density) to be in the order of %g nm^2? If not (e.g. you are not measuring a big amount of emitters or your emitter is not micron-scale in radius), the input data are likely to not correspond to orthodox field emission"""%self.preFactor
            return "apparently orthodox (low)"

        self.orthodoxyMessage = "The data are within field emission orthodoxy limits"
        return "orthodox"
    

class twoEmitterFitter(generalFitter):
    def __init__(self, emitterModel: gt.GETELECModel = gt.GETELECModel()) -> None:
        super().__init__(emitterModel)

        self.fittingParameters = {"fieldConversionFactor": 1., "secondFieldConversionFactor": 2., "prefactor": 1., "secondPrefactor": 1.}
        self.minParameters ={"fieldConversionFactor": 1., "secondFieldConversionFactor": 2., "prefactor": 1., "secondPrefactor": 1.}
        self.initialParameters = {"fieldConversionFactor": 1., "secondFieldConversionFactor": 2., "prefactor": 1., "secondPrefactor": 1.}
        self.maxParameters = {"fieldConversionFactor": 1., "secondFieldConversionFactor": 2., "prefactor": 1., "secondPrefactor": 1.}

    def setParameterRange(self, field:list = [1.e-5, 5., 15.], radius = 2000., gamma = 10., workFunction = 4.5, temperature = 300., prefactor = [0.01, 100., 1.e10], secondPrefactor = [0.01, 100., 1.e10]) ->None:

        """Sets the range of the parameters that will be fitted (variables) and the value of those that will be fixed
        
        Args: 
            field: list of 3 floats. Minimum, initial, and maximum value of the field that is allowed for fitting
            radius, gamma, workFunction, temperature: float or list of 3 floats. If single float, it means that the parameter is considered fixed. If list, it is variable and it works as for field. Default values are 2000. [nm], 10., 4.5 [eV], 300 [K] correspondingly
        """
        self.minParameters ["fieldConversionFactor"] = field[0] / np.min(self.voltageData)
        self.initialParameters["fieldConversionFactor"] = field[1] / np.mean(self.voltageData)
        self.maxParameters["fieldConversionFactor"] = field[2] / np.max(self.voltageData)

        self.minParameters ["secondFieldConversionFactor"] = field[0] / np.min(self.voltageData)
        self.initialParameters["secondFieldConversionFactor"] = field[1] / np.mean(self.voltageData)
        self.maxParameters["secondFieldConversionFactor"] = field[2] / np.max(self.voltageData)

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
                self.initialParameters.pop(key)
                self.minParameters.pop(key)
                self.maxParameters.pop(key)
            except KeyError:
                pass

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
        self.emitterModel.calculateCurrentDensity()
        self.currentDensityFromFristEmitter = np.array(self.emitterModel.getCurrentDensity())

        self.emitterModel.setParameters(field=self.fittingParameters["secondFieldConversionFactor"] * voltageArray)
        self.emitterModel.calculateCurrentDensity()
        self.currentDensityFromSecondEmitter = np.array(self.emitterModel.getCurrentDensity())

        
        return self.fittingParameters["prefactor"] * self.currentDensityFromFristEmitter + self.fittingParameters["secondPrefactor"] * self.currentDensityFromSecondEmitter
    
    def logCurrentDensityError(self, parameterList) -> float:
        """ Calculates the error between a given I-V curve and the one calculated by GETELEC for a given set of parameters.
        The main usage of this function is to be called by external minimization function in order to find the parameterList that
        gives the best fit to the given I-V

        Parameters:
            parameterList: list of emitter parameters (fieldConversionFactor, radius, gamma, workFunction, temperature)
        Returns:
            array of relative errors (in logarithmic scale) between the calculated and the given current densities, after normalizing with a multiplication
            factor that forces them to match at their maximum
        """

        assert len(parameterList) == len(self.fittingParameters), "parameterList has a length of %d while it should have %d"%(len(parameterList), len(self.fittingParameters))
        
        #copy parameter list into the parameter dictionary
        i = 0
        for paramName in self.fittingParameters:
            self.fittingParameters[paramName] = parameterList[i]
            i += 1

        self.calculatedCurrentDensities = self.currentDensityforVoltages(self.voltageData)
        return np.log(self.calculatedCurrentDensities / self.currentData)

    def fitIVCurve(self) -> None:
        """Performs the fitting"""
        self.optimizationData = opt.least_squares(fun = self.logCurrentDensityError, x0 = list(self.initialParameters.values()), \
            bounds =(list(self.minParameters.values()), list(self.maxParameters.values())), method = 'trf', jac = '3-point')
        
        i = 0
        for parameterName in self.fittingParameters:
            self.fittingParameters[parameterName] = self.optimizationData.x[i]
            i += 1
        
        self.fittedCurrent = self.currentDensityforVoltages(self.voltageData)

    def printFittingData(self) -> None:
        betas = sorted([self.fittingParameters["fieldConversionFactor"], self.fittingParameters["secondFieldConversionFactor"]])
        prefactors = sorted([self.fittingParameters["prefactor"], self.fittingParameters["secondPrefactor"]])
        
        print("Optimum parameters: beta_1: %.4g, beta_2: %.4g, Aeff_1: %.4g, Aeff_1: %.4g"%(betas[1], betas[0], prefactors[0], prefactors[1]))
        print("Fitting Error (relative): ", self.getFittingError())


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
            return self.preFactor * self.currentDensityforVoltages(voltageData)



def performFullEmissionFitting(voltageData:np.ndarray, currentData:np.ndarray, field:list = [1., 6., 15.], radius = 2000., gamma = 10., workFunction = 4.5, temperature = 300.):
    
    ivFitter = IVDataFitter()
    ivFitter.setIVdata(voltageData=voltageData, currentData=currentData)
    ivFitter.setParameterRange(field, radius, gamma, workFunction, temperature)
    ivFitter.fitIVCurve()
    ivFitter.printFittingData()
    plotVoltages = 1. / np.linspace(1./min(voltageData), 1./max(voltageData), 128)
    plotFields = ivFitter.fittingParameters["fieldConversionFactor"] * plotVoltages

    fittedCurrentCurve = ivFitter.getOptCurrentCurve(plotVoltages)
    return {"plotFields": plotFields, "fittedCurrents": fittedCurrentCurve} | ivFitter.fittingParameters | {"fittingError": ivFitter.getFittingError(), "orthodoxyStatus": ivFitter.orthodoxyStatus(), "orthodoxyMessage": ivFitter.orthodoxyMessage}


if (__name__ == "__main__"): #some testing operations
    

    ivFitter = IVDataFitter()


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


    outData = performFullEmissionFitting(voltageData, currentData)

    fittedCurrent = outData["fittedCurrents"]
    plt.semilogy(1./voltageData, currentData, '.', label="data" )
    plt.semilogy(outData["fieldConversionFactor"] / outData["plotFields"] , fittedCurrent, label = "beta=%.2g, err=%.2g"%(outData["fieldConversionFactor"], outData["fittingError"]))

    # print(outData)
    plt.savefig("fittedCurveOnlyBeta.png")
    plt.close()

    voltageData = np.array([0.0729348 , 0.07171391, 0.07087798, 0.07008738, 0.06931421, \
       0.06769618, 0.06645514, 0.06525878, 0.06399593, 0.06234496, \
       0.06105263, 0.05937943, 0.05763645])
    
    currentData = np.array([4.24256972e-06, 3.96317508e-06, 3.75021449e-06, 3.51572519e-06, \
       3.30859717e-06, 2.96829892e-06, 2.61283693e-06, 2.40444228e-06, \
       2.17155982e-06, 1.69466134e-06, 1.54505191e-06, 1.34140291e-06, 1.11327359e-06])

    plt.figure()
    ivFitter.setIVdata(voltageData=voltageData, currentData=currentData)
    ivFitter.setParameterRange(radius=[3., 10., 2000.])
    ivFitter.fitIVCurve()
    fittedCurrent = ivFitter.getOptCurrentCurve()

    ivFitter.printFittingData()
    plt.semilogy(1 / ivFitter.fittingParameters["fieldConversionFactor"] / voltageData, fittedCurrent / ivFitter.preFactor, label = "beta = %.2g, R = %.2g, err=%.2g"%(ivFitter.fittingParameters["fieldConversionFactor"], ivFitter.fittingParameters["radius"], ivFitter.getFittingError()))
    plt.semilogy(1 / ivFitter.fittingParameters["fieldConversionFactor"]/voltageData, currentData / ivFitter.preFactor, ".")
    plt.legend()
    plt.grid()
    plt.savefig("fittedCurveChen.png")

    plt.close()
    plt.figure()


    ivFitter.setParameterRange(workFunction=4.6, radius=[1., 10., 200.])

    ivFitter.setIVdata(voltageData, currentData)
    ivFitter.fitIVCurve()
    fittedCurrent = ivFitter.getOptCurrentCurve()
    plt.semilogy(1./voltageData, fittedCurrent, label = "beta = %.2g, R = %.2g, err = %.2g"%(ivFitter.fittingParameters["fieldConversionFactor"], ivFitter.fittingParameters["radius"], ivFitter.getFittingError()))

    plt.legend()

    plt.savefig("fittedcurve2Emitter.png")
    plt.close()

    ivFitter.printFittingData()


    voltageData = 1./np.array([0.05401709, 0.05757835, 0.05907407, 0.0605698 , 0.06185185, \
       0.06349003, 0.06519943, 0.06655271, 0.06733618, 0.06840456, \
       0.06990028, 0.07153846, 0.07353276, 0.07566952, 0.07766382,\
       0.07958689, 0.0810114 , 0.0824359 , 0.08407407, 0.0855698 ,\
       0.08777778, 0.08934473, 0.09112536, 0.09247863, 0.09447293,\
       0.09603989, 0.09817664, 0.09988604, 0.10173789, 0.10330484,\
       0.10529915, 0.10665242, 0.10814815, 0.10935897, 0.11064103,\
       0.1122792 , 0.1139886 ])
    
    currentData = 1.e-9 * np.array([1.32987103e+06, 6.30957344e+05, 4.64158883e+05, 3.34048498e+05,\
       2.68269580e+05, 1.93069773e+05, 1.35935639e+05, 1.04483478e+05,\
       8.76712387e+04, 7.04075201e+04, 5.17947468e+04, 3.64673967e+04,\
       2.40409918e+04, 1.55051578e+04, 9.78309319e+03, 6.59246176e+03,\
       4.84969343e+03, 3.56763941e+03, 2.74217545e+03, 1.93069773e+03,\
       1.30102522e+03, 9.16019591e+02, 6.73862717e+02, 5.06712835e+02,\
       3.56763941e+02, 2.56757897e+02, 1.65595152e+02, 1.19176459e+02,\
       8.76712387e+01, 6.59246176e+01, 4.64158883e+01, 3.56763941e+01,\
       2.92864456e+01, 2.35195264e+01, 1.88881958e+01, 1.55051578e+01,\
       1.14062492e+01])
    
    ivFitter = twoEmitterFitter()

    ivFitter.setIVdata(voltageData=voltageData, currentData=currentData)
    ivFitter.setParameterRange()
    ivFitter.fitIVCurve()
    fittedCurrent = ivFitter.getOptCurrentCurve(voltageData)
    plt.semilogy(1./voltageData, currentData, '.', label="data" )
    plt.semilogy(1./voltageData, fittedCurrent, label = "beta1=%.2g, beta2:%.2g, pref1:%.2g, pref2:%.2g, err=%.2g"%(ivFitter.fittingParameters["fieldConversionFactor"], ivFitter.fittingParameters["secondFieldConversionFactor"], ivFitter.fittingParameters["prefactor"], ivFitter.fittingParameters["secondPrefactor"] , ivFitter.getFittingError()))

    ivFitter.printFittingData()
    plt.savefig("twoEmitterFitting.png")