from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Type, TypeVar, Union, cast

import numpy as np
import copy
import json
import inspect

from getelec import MetalEmitter, SemiconductorEmitter, BandEmitter, Globals

def getArgument(arg, idx):

    if isinstance(arg, (np.ndarray, list)):
        if idx >= len(arg):
            print(f"WARNING, one of the arrays is shorter than the others. Last value of the array was copied to match length.")
            return arg[-1]
        else:
            return arg[idx]
    else:
        return arg

class GETELECModel():

    def __init__(
        
        #AK: set the parameters below to match the definitions given in the refactored code

        self,
        emitterType: Optional[str] = None,
        field: Optional[np.ndarray] = None,
        radius: Optional[np.ndarray] = None,
        gamma: Optional[np.ndarray] = None,
        workFunction: Optional[np.ndarray] = None,
        temperature: Optional[np.ndarray] = None,
        emitter: Optional[BandEmitter] = None,
        conductionBandBottom: Optional[float] = None,
        fermiLevel: Optional[float] = None,
        bandGap: Optional[float] = None,
        effectiveMassConduction: Optional[float] = None,
        effectiveMassValence: Optional[float] = None,
        numberOfSpectrumPoints: Optional[int] = None,

        **kwargs: Any

    ) -> None:

        """Create GETELEC model.

        Parameters
        ----------

        emitterType:
            Type of emitter, "metal" or "semiconductor"
            Defining this overrides previous Emitter instance
            and creates a new, empty one.

        field:
            NueffectiveMassValencey array of field float values, [V/nm]

        radius:
            NueffectiveMassValencey array of radius float values, [nm]

        gamma:
            NueffectiveMassValencey array of gamma float values, dimensionless

        workFunction:
            NueffectiveMassValencey array of work function float values, [eV]

        temperature:
            NueffectiveMassValencey array of temperature float values, [K]

        emitter:
            Object of the emitter, MetalEmitter or SemiconductorEmitter

        conductionBandBottom:
            Energy of the bottom of conduction band, [eV]

        fermiLevel:
            Fermi level, [eV]

        bandGap:
            Band gap, [eV]

        effectiveMassConduction:
            Effective electron mass, fraction m/me

        effectiveMassValence:
            Effective hole mass, fraction m/me
            
        numberOfSpectrumPoints:
            Number of spectras to create, int

        """

        self.emitterType = emitterType
        self.field = field
        self.radius = radius
        self.gamma = gamma
        self.workFunction = workFunction
        self.temperature = temperature
        self.conductionBandBottom = conductionBandBottom
        self.fermiLevel = fermiLevel
        self.bandGap = bandGap
        self.effectiveMassConduction = effectiveMassConduction
        self.effectiveMassValence = effectiveMassValence

        self.numberOfSpectrumPoints = numberOfSpectrumPoints

        self.currentDensity = None
        self.nottinghamHeat = None
        self.electronSpectrum = None
        

        if(emitter == None):
            if(emitterType != None):
                if(emitterType == 'metal'): self.emitter = MetalEmitter()
                elif(emitterType == 'semiconductor'): self.emitter = SemiconductorEmitter()
                else: raise ValueError("emitterType has to be 'metal' or 'semiconductor'")
            else: 
                if(emitterType != None):
                    raise ValueError("If emitter object is not passed, emitterType has to be set to 'metal' or 'semiconductor'")
        else:
            if(type(emitter) == MetalEmitter):
                self.emitterType = 'metal'
            elif(type(emitter) == SemiconductorEmitter):
                self.emitterType = 'semiconductor'
            else: raise ValueError("Passed object as 'emitter' has to be of type 'MetalEmitter' or 'SemiconductorEmitter'")
            self.emitter = emitter


        if(emitterType != None):
            if(emitterType == 'metal'): self.emitter = MetalEmitter()
            elif(emitterType == 'semiconductor'): self.emitter = SemiconductorEmitter()
            else: raise ValueError("emitterType has to be 'metal' or 'semiconductor'")
        else: self.emitter = emitter

        if kwargs:
            self.kwargs = kwargs

    @classmethod
    def _getParameterNames(cls):
        """
        Get parameter names for the estimator

        Returns
        -------
        params: dict
            Sorted parameter names
        """

        init = getattr(cls.__init__, "deprecated_original", cls.__init__)
        if init is object.__init__:
            return []

        init_signature = inspect.signature(init)
        parameters = [
            p
            for p in init_signature.parameters.values()
            if p.name != "self" and p.kind != p.VAR_KEYWORD
        ]

        return sorted([p.name for p in parameters])

    def getParameters(self):
        """
        Get parameters for this estimator.

        Returns
        -------
        params : dict
            Parameter names mapped to their values.
        """
        out = dict()
        for key in self._getParameterNames():
            value = getattr(self, key)
            if hasattr(value, "getParameters") and not isinstance(value, type):
                deep_items = value.getParameters().items()
                out.update((key + "__" + k, val) for k, val in deep_items)
            out[key] = value
        return out

    def setParameters(self, **params: Any) -> "GETELECModel":
        """Set the parameters of this estimator.

        Returns
        -------
        self

        """
        if not params:
            return self

        for key, value in params.items():
            if hasattr(self, key):
                setattr(self, key, value)
                if(key == 'emitterType'):
                    if(value == 'metal'):
                        self.emitter = MetalEmitter()
                    elif(value == 'semiconductor'):
                        self.emitter = SemiconductorEmitter()
                    else:
                        raise ValueError("emitterType has to be 'metal' or 'semiconductor'")
            else:
                if not hasattr(self, "kwargs"):
                    self.kwargs = {}
                self.kwargs[key] = value

        self.currentDensity = None 
        self.nottinghamHeat = None
        self.electronSpectrum = None

        return self

    def saveModel():
        return
    
    def loadModel():
        return

    def run(self, calculateCurrent: Optional[bool] = False, calculateNottinghamHeat: Optional[bool] = False, calculateSpectrum: Optional[bool] = False):

        """Runs the model by specifying what properties to compute.
        Allows running multiple calculations at the same time

        Parameters
        ----------

        calculateCurrent:
            Select whether model should calculate emitted current from the emitter, boolean [true/false]
        
        calculateNottinghamHeat:
            Select whether model should calculate emitted Nottingham heat from the emitter, boolean [true/false]

        calculateSpectrum:
            Select whether model should evaluate electron spectrum, boolean [true/false]

        -------
        Returns
        -------

            self

        """

        if(calculateCurrent): _currentDensity = np.array(self.field, dtype=float)
        if(calculateNottinghamHeat): _nottinghamHeat = np.array(self.field, dtype=float)
        
        if(calculateSpectrum): 

            energy = []
            electronCount = []

        nPoints = 256 if self.numberOfSpectrumPoints is None else self.numberOfSpectrumPoints

        kT = Globals.BoltzmannConstant * np.array(self.temperature, dtype=float)

        for i in range(len(self.field)):

            self.emitter.barrier.setParameters(getArgument(self.field, i), getArgument(self.radius, i), getArgument(self.gamma, i))

            if(self.emitterType == 'metal'):

                self.emitter.setParameters(getArgument(self.workFunction, i), getArgument(kT, i))

                if(calculateCurrent): _currentDensity[i] = self.emitter.currentDensity()
                if(calculateNottinghamHeat): _nottinghamHeat[i] = self.emitter.nottinghamHeat()

                if(calculateSpectrum): 

                    _energy, _electronCount = self.emitter.totalEnergySpectrumArrays(nPoints)

                    energy.append(_energy)
                    electronCount.append(_electronCount)

            elif(self.emitterType == 'semiconductor'):

                self.emitter.setParameters(getArgument(self.field, i), getArgument(self.radius, i), getArgument(self.gamma, i), self.conductionBandBottom, self.fermiLevel, self.bandGap, getArgument(kT, i), self.effectiveMassConduction, self.effectiveMassValence)

                if(calculateCurrent): _currentDensity[i] = self.emitter.currentDensity()

                if(calculateNottinghamHeat): _nottinghamHeat[i] = self.emitter.nottinghamHeat()

                if(calculateSpectrum):
                    
                    _energy, _electronCount = self.emitter.totalEnergySpectrumArrays(nPoints)

                    energy.append(_energy)
                    electronCount.append(_electronCount)

        if(calculateCurrent): self.currentDensity = _currentDensity
        if(calculateNottinghamHeat): self.nottinghamHeat = _nottinghamHeat
        if(calculateSpectrum): self.electronSpectrum = {'energy': energy, 'electronCount': electronCount}

        return self

    def getNottinghamHeat(self):

        """Get calculated Nottingham heat from the emitter
        
        Returns numpy array of Nottigham heat values
        """
        #if(self.nottinghamHeat == None): print("WARNING, you have asked Nottingham heat without calculating it first. You might want to run getelecModel.run(calculateNottinghamHeat=True) first")

        return self.nottinghamHeat
    
    def getCurrentDensity(self):

        """Get calculated emitted current density from the emitter.

        Returns numpy array of emitted currents (electrons / area * time)
        
        """
        #if(self.currentDensity == None): print("WARNING, you have asked currentDensity without calculating it first. You might want to run getelecModel.run(calculateCurrentDensity=True) first")

        return self.currentDensity

    def getElectronSpectrum(self):

        """Get calculated electron spectrum from the emitter
        
        Returns dictionary of numpy arrays {energy, electron_count}
        
        """
        #if(self.nottinghamHeat == None): print("WARNING, you have asked Electron Spectrum without calculating it first. You might want to run getelecModel.run(calculateElectronSpectrum=True) first")

        return self.electronSpectrum

    def calculateNottinghamHeat(self):
        """Calculate Nottingham heat for metals

        Returns self
        """
        self.run(calculateNottinghamHeat=True)

        return self
  
    def calculateCurrentDensity(self):
        """Calculate emitted current for metals
        
        Returns self
        """

        self.run(calculateCurrent=True)

        return self

    def calculateElectronSpectrum(self):
        """Calculate electron spectrum for metals
        
        Returns self
        """

        self.run(calculateSpectrum=True)

        return self