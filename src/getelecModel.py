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

from getelec import MetalEmitter, SemiconductorEmitter, IVDataFitter, BandEmitter

class GETELECModel():

    def __init__(

        self,
        emitter_type: Optional[str] = None,
        field: Optional[np.ndarray] = None,
        radius: Optional[np.ndarray] = None,
        gamma: Optional[np.ndarray] = None,
        work_function: Optional[np.ndarray] = None,
        temperature: Optional[np.ndarray] = None,
        emitter_object: Optional[BandEmitter] = None,
        **kwargs: Any

    ) -> None:

        """
        Param emitter_type: Type of emitter, "metal" or "semiconductor"
        Param field: NumPy array of field values, V/nm
        Param radius: NumPy array of radius values, nm
        Param gamma: NumPy array of gamma values,  dimensionless
        Param work_function: NumPy array of work function values, eV
        Param temperature: NumPy array of temperature values, K
        Param emitter_object: Object of the emitter, MetalEmitter or SemiconductorEmitter or BandEmitter
        """

        self.emitter_type = emitter_type
        self.field = field
        self.radius = radius
        self.gamma = gamma
        self.work_function = work_function
        self.temperature = temperature
        self.emitter_object = emitter_object

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
        """Get parameter names for the estimator"""
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

    def nottingham_heat(self):

        if(self.emitter_type == 'metal'):

            MetalEmitter

            return

        if(self.emitter_type == 'semiconductor'):

            SemiconductorEmitter

            return

        raise ValueError("emitter_type has to be 'metal' or 'semiconductor'")
    
    def emitted_current(self):

        if(self.emitter_type == 'metal'):

            MetalEmitter

            return

        if(self.emitter_type == 'semiconductor'):

            SemiconductorEmitter

            return

        raise ValueError("emitter_type has to be 'metal' or 'semiconductor'")

    def electron_spectrum(self):

        if(self.emitter_type == 'metal'):

            MetalEmitter

            return

        if(self.emitter_type == 'semiconductor'):

            SemiconductorEmitter

            return

        raise ValueError("emitter_type has to be 'metal' or 'semiconductor'")

#########################








model = GETELECModel()

model.set_params(emitter='metal')

model.get_params()

current = model.emitted_current()

