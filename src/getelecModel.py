from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Type, TypeVar, Union, cast

import numpy as np
import copy
import json
import inspect

from getelec import MetalEmitter, SemiconductorEmitter, BandEmitter, Globals

def getArgument(arg, idx):

    if isinstance(arg, (np.ndarray, list)):
        if idx >= len(arg):
            return arg[-1]
        else:
            return arg[idx]
    else:
        return arg

class GETELECModel():

    def __init__(

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
            and creates a new, eeffectiveMassValencety one.

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
    
    def setParameters(self, **params) -> "GETELECModel":
        """Set the parameters of the model. Allows unknown kwargs.
        Returns self.
        """

        if not params: return self

        for key, value in params.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                if not hasattr(self, "kwargs"):
                    self.kwargs = {}
                self.kwargs[key] = value

        return self

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
                if not hasattr(self, "kwargs"):
                    self.kwargs = {}
                self.kwargs[key] = value


        return self

    def saveModel():
        return
    
    def loadModel():
        return

    def getNottinghamHeat(self):

        """Calculate Nottigham heat from the emitter
        
        Returns nueffectiveMassValencey array of Nottigham heat values
        """

        if(self.emitterType == 'metal'):

            return self._nottinghamHeatMetal().nottinghamHeat

        if(self.emitterType == 'semiconductor'):

            return self._nottinghamHeatSemiconductor().nottinghamHeat

        raise ValueError("emitterType has to be 'metal' or 'semiconductor'")
    
    def getCurrentDensity(self):

        """Calculate emitted current density from the emitter.

        Returns nueffectiveMassValencey array of emitted currents (electrons / area * time)
        
        """

        if(self.emitterType == 'metal'): 

            return (self._currentDensityMetal()).currentDensity

        if(self.emitterType == 'semiconductor'):

            return (self._currentDensitySemiconductor()).currentDensity

        raise ValueError("emitterType has to be 'metal' or 'semiconductor'")

    def getElectronSpectrum(self):

        """Calculate electron spectrum from the emitter
        
        Returns dictionary of nueffectiveMassValencey arrays (energy, electron_count)
        
        """

        if(self.emitterType == 'metal'):

            _energy, _electronCount = self._electronSpectrumMetal().electronSpectrum

            return {'energy': _energy, 'electron_count': _electronCount}

        if(self.emitterType == 'semiconductor'):

            _energy, _electronCount = self._electronSpectrumSemiconductor().electronSpectrum

            return {'energy': _energy, 'electron_count': _electronCount}

        raise ValueError("emitterType has to be 'metal' or 'semiconductor'")

    def _nottinghamHeatMetal(self):
        """Calculate Nottingham heat for metals

        Returns self
        """
        kT = Globals.BoltzmannConstant * self.temperature

        _nottinghamHeat = np.array(self.field, dtype=float)

        for i in range(len(_nottinghamHeat)):

            self.emitter.barrier.setParameters(getArgument(self.field, i), getArgument(self.radius, i), getArgument(self.gamma, i))
            self.emitter.setParameters(getArgument(self.workFunction, i), getArgument(kT, i))

            _nottinghamHeat[i] = self.emitter.nottinghamHeat()

        self.nottinghamHeat = _nottinghamHeat

        return self
    
    def _nottinghamHeatSemiconductor(self):
        """Calculate Nottingham heat for semiconductors
        
        Returns self
        """
        kT = Globals.BoltzmannConstant * self.temperature

        _nottinghamHeat = np.array(self.field, dtype = float)

        for i in range(len(_nottinghamHeat)):

            self.emitter.barrier.setParameters(getArgument(self.field, i), getArgument(self.radius, i), getArgument(self.gamma, i))
            self.emitter.setParameters(getArgument(self.field, i), getArgument(self.radius, i), getArgument(self.gamma, i), self.conductionBandBottom, self.fermiLevel, self.bandGap, getArgument(kT, i), self.effectiveMassConduction, self.effectiveMassValence)

            _nottinghamHeat[i] = self.emitter.nottinghamHeat()

        self.nottinghamHeat = _nottinghamHeat

        return self
    
    def _currentDensityMetal(self):
        """Calculate emitted current for metals
        
        Returns self
        """

        
        kT = Globals.BoltzmannConstant * self.temperature

        _currentDensity = np.array(self.field, dtype=float)

        for i in range(len(_currentDensity)):

            self.emitter.barrier.setParameters(getArgument(self.field, i), getArgument(self.radius, i), getArgument(self.gamma, i))
            self.emitter.setParameters(getArgument(self.workFunction, i), getArgument(kT, i))

            _currentDensity[i] = self.emitter.currentDensity()

        self.currentDensity = _currentDensity

        return self

    def _currentDensitySemiconductor(self):
        """Calculate emitted current for semiconductors
        
        Returns self
        """
        kT = Globals.BoltzmannConstant * self.temperature
        
        _currentDensity = np.array(self.field, dtype=float)

        for i in range(len(_currentDensity)):

            self.emitter.barrier.setParameters(getArgument(self.field, i), getArgument(self.radius, i), getArgument(self.gamma, i))
            self.emitter.setParameters(getArgument(self.field, i), getArgument(self.radius, i), getArgument(self.gamma, i), self.conductionBandBottom, self.fermiLevel, self.bandGap, getArgument(kT, i), self.effectiveMassConduction, self.effectiveMassValence)

            _currentDensity[i] = self.emitter.currentDensity()

        self.currentDensity = _currentDensity

        return self

    def _electronSpectrumMetal(self):
        """Calculate electron spectrum for metals
        
        Returns self
        """
        kT = Globals.BoltzmannConstant * np.array(self.temperature, dtype=float)

        energy = []
        electron_count = []

        n_points = 256 if self.numberOfSpectrumPoints is None else self.numberOfSpectrumPoints

        for i in range(len(self.field)):

            self.emitter.barrier.setParameters(getArgument(self.field, i), getArgument(self.radius, i), getArgument(self.gamma, i))
            self.emitter.setParameters(getArgument(self.workFunction, i), getArgument(kT, i))

            _energy, _electron_count = self.emitter.totalEnergySpectrumArrays(n_points)

            energy.append(_energy)
            electron_count.append(_electron_count)

        self.electronSpectrum = (energy, electron_count)

        return self

    def _electronSpectrumSemiconductor(self):
        """Calculate electron spectrum for semiconductors
        
        Returns self
        """

        kT = Globals.BoltzmannConstant * np.array(self.temperature, dtype=float)

        energy = []
        electron_count = []

        n_points = 256 if self.numberOfSpectrumPoints is None else self.numberOfSpectrumPoints

        for i in range(len(self.field)):

            self.emitter.barrier.setParameters(getArgument(self.field, i), getArgument(self.radius, i), getArgument(self.gamma, i))
            self.emitter.setParameters(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i), self.conductionBandBottom, self.fermiLevel, self.bandGap, getArgument(kT, i), self.effectiveMassConduction, self.effectiveMassValence)

            _energy, _electron_count = self.emitter.totalEnergyDistribution(n_points)

            energy.append(_energy)
            electron_count.append(_electron_count)

        self.electronSpectrum = (energy, electron_count)

        return self

#########################

if(__name__ == '__main__'):

    model = GETELECModel()

    field = np.array([1, 2, 3, 4, 5])
    radius = np.full(5, 50)
    gamma = np.full(5, 10)
    workFunction = np.full(5, 4)
    temperature = np.full(5, 300)

    conductionBandBottom = 4.05
    fermiLevel = 4.61
    bandGap = 1.12
    effectiveMassConduction = 0.7
    effectiveMassValence = 0.5


    model.setParameters(emitterType='semiconductor', field = field, radius = radius, gamma = gamma, temperature = temperature, conductionBandBottom = conductionBandBottom, fermiLevel = fermiLevel, bandGap = bandGap, effectiveMassConduction = effectiveMassConduction, effectiveMassValence = effectiveMassValence, numberOfSpectrumPoints = 4)

    print(model.getParameters())

    res = model.getElectronSpectrum()

    print(res)


