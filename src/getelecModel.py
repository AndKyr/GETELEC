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
        emitter_type: Optional[str] = None,
        field: Optional[np.ndarray] = None,
        radius: Optional[np.ndarray] = None,
        gamma: Optional[np.ndarray] = None,
        work_function: Optional[np.ndarray] = None,
        temperature: Optional[np.ndarray] = None,
        emitter: Optional[BandEmitter] = None,
        ec: Optional[float] = None,
        ef: Optional[float] = None,
        eg: Optional[float] = None,
        me: Optional[float] = None,
        mp: Optional[float] = None,

        **kwargs: Any

    ) -> None:

        """Create GETELEC model.

        Parameters
        ----------

        emitter_type:
            Type of emitter, "metal" or "semiconductor"
            Defining this overrides previous Emitter instance
            and creates a new, empty one.
        field:
            NumPy array of field float values, [V/nm]
        radius:
            NumPy array of radius float values, [nm]
        gamma:
            NumPy array of gamma float values, dimensionless
        work_function:
            NumPy array of work function float values, [eV]
        temperature:
            NumPy array of temperature float values, [K]
        emitter:
            Object of the emitter, MetalEmitter or SemiconductorEmitter
        Ec:
            Energy of the bottom of conduction band, [eV]
        Ef:
            Fermi level, [eV]
        Eg:
            Ban gap, [eV]
        me:
            effective electron mass, fraction m/me
        mp:
            effective hole mass, fraction m/me
        """

        self.emitter_type = emitter_type
        self.field = field
        self.radius = radius
        self.gamma = gamma
        self.work_function = work_function
        self.temperature = temperature
        self.ec = ec
        self.ef = ef
        self.eg = eg
        self.me = me
        self.mp = mp

        self.emitted_current_density = None
        self.nottingham_heat = None
        self.electron_spectrum = None

        if(emitter == None):
            if(emitter_type != None):
                if(emitter_type == 'metal'): self.emitter = MetalEmitter()
                elif(emitter_type == 'semiconductor'): self.emitter = SemiconductorEmitter()
                else: raise ValueError("emitter_type has to be 'metal' or 'semiconductor'")
            else: 
                if(emitter_type != None):
                    raise ValueError("If emitter object is not passed, emitter_type has to be set to 'metal' or 'semiconductor'")
        else:
            if(type(emitter) == MetalEmitter):
                self.emitter_type = 'metal'
            elif(type(emitter) == SemiconductorEmitter):
                self.emitter_type = 'semiconductor'
            else: raise ValueError("Passed object as 'emitter' has to be of type 'MetalEmitter' or 'SemiconductorEmitter'")
            self.emitter = emitter


        if(emitter_type != None):
            if(emitter_type == 'metal'): self.emitter = MetalEmitter()
            elif(emitter_type == 'semiconductor'): self.emitter = SemiconductorEmitter()
            else: raise ValueError("emitter_type has to be 'metal' or 'semiconductor'")
        else: self.emitter = emitter

        if kwargs:
            self.kwargs = kwargs
    
    def set_params(self, **params) -> "GETELECModel":
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
    def _get_param_names(cls):
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

    def get_params(self):
        """
        Get parameters for this estimator.

        Returns
        -------
        params : dict
            Parameter names mapped to their values.
        """
        out = dict()
        for key in self._get_param_names():
            value = getattr(self, key)
            if hasattr(value, "get_params") and not isinstance(value, type):
                deep_items = value.get_params().items()
                out.update((key + "__" + k, val) for k, val in deep_items)
            out[key] = value
        return out

    def set_params(self, **params: Any) -> "GETELECModel":
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
                if(key == 'emitter_type'):
                    if(value == 'metal'):
                        self.emitter = MetalEmitter()
                    elif(value == 'semiconductor'):
                        self.emitter = SemiconductorEmitter()
            else:
                if not hasattr(self, "kwargs"):
                    self.kwargs = {}
                self.kwargs[key] = value


        return self

    def save_model():
        return
    
    def load_model():
        return

    def get_nottingham_heat(self):

        """Calculate Nottigham heat from the emitter
        
        Returns numpy array of Nottigham heat values
        """

        if(self.emitter_type == 'metal'):

            return self._nottingham_heat_metal().nottingham_heat

        if(self.emitter_type == 'semiconductor'):

            return self._nottingham_heat_semiconductor().nottingham_heat

        raise ValueError("emitter_type has to be 'metal' or 'semiconductor'")
    
    def get_emitted_current_density(self):

        """Calculate emitted current density from the emitter.

        Returns numpy array of emitted currents (electrons / area * time)
        
        """

        if(self.emitter_type == 'metal'): 

            return (self._emitted_current_density_metal()).emitted_current_density

        if(self.emitter_type == 'semiconductor'):

            return (self._emitted_current_density_semiconductor()).emitted_current_density

        raise ValueError("emitter_type has to be 'metal' or 'semiconductor'")

    def get_electron_spectrum(self):

        """Calculate electron spectrum from the emitter
        
        Returns
        
        """

        if(self.emitter_type == 'metal'):

            return self._electron_spectrum_metal()

        if(self.emitter_type == 'semiconductor'):

            return self._electron_spectrum_semiconductor()

        raise ValueError("emitter_type has to be 'metal' or 'semiconductor'")

    def _nottingham_heat_metal(self):
        """Calculate Nottingham heat for metals

        Returns self
        """
        kT = Globals.BoltzmannConstant * self.temperature

        _nottingham_heat = np.array(self.field, dtype=float)

        for i in range(len(_nottingham_heat)):

            self.emitter.barrier.setParameters(getArgument(self.field, i), getArgument(self.radius, i), getArgument(self.gamma, i))
            self.emitter.setParameters(getArgument(self.work_function, i), getArgument(kT, i))

            _nottingham_heat[i] = self.emitter.nottinghamHeat()

        self.nottingham_heat = _nottingham_heat

        return self
    
    def _nottingham_heat_semiconductor(self):
        """Calculate Nottingham heat for semiconductors
        
        Returns self
        """
        kT = Globals.BoltzmannConstant * self.temperature

        _nottingham_heat = np.array(self.field, dtype = float)

        for i in range(len(_nottingham_heat)):

            self.emitter.barrier.setParameters(getArgument(self.field, i), getArgument(self.radius, i), getArgument(self.gamma, i))
            self.emitter.setParameters(getArgument(self.field, i), getArgument(self.radius, i), getArgument(self.gamma, i), self.ec, self.ef, self.eg, getArgument(kT, i), self.me, self.mp)

            _nottingham_heat[i] = self.emitter.nottinghamHeat()


        return self
    
    def _emitted_current_density_metal(self):
        """Calculate emitted current for metals
        
        Returns self
        """
        kT = Globals.BoltzmannConstant * self.temperature

        _emitted_current_density = np.array(self.field, dtype=float)

        for i in range(len(_emitted_current_density)):

            self.emitter.barrier.setParameters(getArgument(self.field, i), getArgument(self.radius, i), getArgument(self.gamma, i))
            self.emitter.setParameters(getArgument(self.work_function, i), getArgument(kT, i))

            _emitted_current_density[i] = self.emitter.currentDensity()

        self.emitted_current_density = _emitted_current_density

        return self

    def _emitted_current_density_semiconductor(self):
        """Calculate emitted current for semiconductors
        
        Returns self
        """
        kT = Globals.BoltzmannConstant * self.temperature
        
        _emitted_current_density = np.array(self.field, dtype=float)

        for i in range(len(_emitted_current_density)):

            self.emitter.barrier.setParameters(getArgument(self.field, i), getArgument(self.radius, i), getArgument(self.gamma, i))
            self.emitter.setParameters(getArgument(self.field, i), getArgument(self.radius, i), getArgument(self.gamma, i), self.ec, self.ef, self.eg, getArgument(kT, i), self.me, self.mp)

            _emitted_current_density[i] = self.emitter.currentDensity()

        self.emitted_current_density = _emitted_current_density

        return self

    def _electron_spectrum_metal(self):
        """Calculate electron spectrum for metals
        
        Returns self
        """

        return

    def _electron_spectrum_semiconductor(self):
        """Calculate electron spectrum for semiconductors
        
        Returns self
        """

        return

#########################


if(__name__ == "__main__"):

    model = GETELECModel()

    field = np.array([1, 2, 3, 4, 5])
    radius = np.full(5, 50)
    gamma = np.full(5, 10)
    work_function = np.full(5, 4)
    temperature = np.full(5, 300)

    model.set_params(emitter_type='metal', field = field, radius = radius, gamma = gamma, work_function = work_function, temperature = temperature)

    print(model.get_params())

    current = model.get_emitted_current_density()

    print(current)

    print(dir(model.get_params()['emitter']))


