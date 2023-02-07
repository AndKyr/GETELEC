from typing import (
    Any,
    Callable,
    Dict,
    List,
    Optional,
    Sequence,
    Tuple,
    Type,
    TypeVar,
    Union,
    cast,
)

import numpy as np
import copy
import json
import inspect
import xgboost as xgb

from getelec import MetalEmitter, SemiconductorEmitter, IVDataFitter, BandEmitter, Globals

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
        field: Optional[Union[np.ndarray, List[float], Tuple[float]]] = None,
        radius: Optional[Union[np.ndarray, List[float], Tuple[float]]] = None,
        gamma: Optional[Union[np.ndarray, List[float], Tuple[float]]] = None,
        work_function: Optional[Union[np.ndarray, List[float], Tuple[float]]] = None,
        temperature: Optional[Union[np.ndarray, List[float], Tuple[float]]] = None,
        emitter: Optional[BandEmitter] = None,
        **kwargs: Any

    ) -> None:

        """Create GETELEC model.

        Parameters
        ----------

        emitter_type:
            Type of emitter, "metal" or "semiconductor"
        field:
            NumPy array of field float values, V/nm
        radius:
            NumPy array of radius float values, nm
        gamma:
            NumPy array of gamma float values,  dimensionless
        work_function:
            NumPy array of work function float values, eV
        temperature:
            NumPy array of temperature float values, K
        emitter:
            Object of the emitter, MetalEmitter or SemiconductorEmitter or BandEmitter
        """

        self.emitter_type = emitter_type
        self.field = field
        self.radius = radius
        self.gamma = gamma
        self.work_function = work_function
        self.temperature = temperature
        self.emitter = emitter

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

            return self._nottingham_heat_metal()

        if(self.emitter_type == 'semiconductor'):

            return self._nottingham_heat_semiconductor()

        raise ValueError("emitter_type has to be 'metal' or 'semiconductor'")
    
    def get_emitted_current_density(self):

        """Calculate emitted current density from the emitter.

        Returns numpy array of emitted currents (electrons / area * time)
        
        """

        if(self.emitter_type == 'metal'): 

            return self._emitted_current_density_metal()

        if(self.emitter_type == 'semiconductor'):

            return self._emitted_current_density_semiconductor()

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
    
    def _nottingham_heat_semiconductor(self):
        """Calculate Nottingham heat for semiconductors
        
        Returns self
        """

        return
    
    def _emitted_current_density_metal(self):
        """Calculate emitted current for metals
        
        Returns emitted current density array
        """

        self.emitter = MetalEmitter()
        kT = Globals.BoltzmannConstant * self.temperature

        _emitted_current_density = np.copy(self.field, )

        for i in range(len(_emitted_current_density)):

            self.emitter.barrier.setParameters(getArgument(self.field, i), getArgument(self.radius, i), getArgument(self.gamma, i))
            self.emitter.setParameters(getArgument(self.work_function, i), getArgument(kT, i))

            _emitted_current_density[i] = self.emitter.currentDensity()

        return _emitted_current_density

    def _emitted_current_density_semiconductor(self):
        """Calculate emitted current for semiconductors
        
        Returns self
        """

        return

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


