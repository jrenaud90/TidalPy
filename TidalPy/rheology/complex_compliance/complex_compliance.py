from typing import List, Dict, TYPE_CHECKING, Callable, Tuple

import numpy as np

from . import known_models, known_model_live_args, known_model_const_args
from .defaults import complex_compliance_defaults
from ...exceptions import InitiatedPropertyChangeError, IncorrectMethodToSetStateProperty, OuterscopePropertySetError
from ...tides.mode_manipulation import FreqSig
from ...utilities.classes.model import LayerModelHolder
from ...utilities.performance import njit
from ...utilities.types import FloatArray, ComplexArray

if TYPE_CHECKING:
    from ..rheology import Rheology
    from ...structures.layers import PhysicalLayerType


@njit(cacheable=False)
def compliance_dict_helper(tidal_frequencies: Dict[FreqSig, FloatArray], compliance_func: Callable,
                           live_inputs: Tuple[FloatArray, ...], inputs: Tuple[float, ...]):
    """ Njit-safe Complex compliance calculator helper - Array version

    Parameters
    ----------
    tidal_frequencies : Dict[FreqSig, FloatArray]
        Njit-safe TypedDict of tidal frequencies.
    compliance_func : Callable
        Complex compliance function.
    live_inputs : Tuple[FloatArray, ...]
        Live inputs to the complex compliance function.
    inputs : Tuple[float, ...]
        Static inputs to the complex compliance function.

    Returns
    -------
    complex_compliance_by_tidal_freq : Dict[FreqSig, ComplexArray]
    """

    # Build fake dictionary so that njit can compile the function
    fake_index = list(tidal_frequencies.keys())[0]
    fake_freq = tidal_frequencies[fake_index]
    compliance = live_inputs[0]
    complex_compliance_by_tidal_freq = {(-100, -100): fake_freq * compliance * (1. + 1.j)}

    # Find the complex compliance for each frequency
    for freq_sig, freq in tidal_frequencies.items():
        complex_compliance = compliance_func(freq, *live_inputs, *inputs)
        complex_compliance_by_tidal_freq[freq_sig] = complex_compliance

    # Delete the fake complex compliance we added earlier
    del complex_compliance_by_tidal_freq[(-100, -100)]

    return complex_compliance_by_tidal_freq

# @njit(cacheable=True)
# def compliance_dict_helper_float(tidal_frequencies: Dict[FreqSig, FloatArray], compliance_func_float: Callable,
#                                  live_inputs: Tuple[FloatArray, ...], inputs: Tuple[float, ...]):
#     """ Njit-safe Complex compliance calculator helper - Float version
#     # TODO: Do we need a separte func for this? Can we use the one above?
#
#     Parameters
#     ----------
#     tidal_frequencies : Dict[FreqSig, FloatArray]
#         Njit-safe TypedDict of tidal frequencies.
#     compliance_func_float : Callable
#         Complex compliance function - float version.
#     live_inputs : Tuple[FloatArray, ...]
#         Live inputs to the complex compliance function.
#     inputs : Tuple[float, ...]
#         Static inputs to the complex compliance function.
#
#     Returns
#     -------
#     complex_compliance_by_tidal_freq : Dict[FreqSig, ComplexArray]
#     """
#
#     # Build fake dictionary so that njit can compile the function
#     fake_index = list(tidal_frequencies.keys())[0]
#     fake_freq = tidal_frequencies[fake_index]
#     complex_compliance_by_tidal_freq = {(-100, -100): fake_freq * (1. + 1.j)}
#
#     # Find the complex compliance for each frequency
#     for freq_sig, freq in tidal_frequencies.items():
#         complex_compliance_by_tidal_freq[freq_sig] = compliance_func_float(freq, *live_inputs, *inputs)
#
#     # Delete the fake complex compliance we added earlier
#     del complex_compliance_by_tidal_freq[(-100, -100)]
#
#     return complex_compliance_by_tidal_freq


class ComplexCompliance(LayerModelHolder):

    """ ComplexCompliance
    Complex compliance provides the functionality to calculate the complex compliance once provided a viscosity, shear
        modulus, and forcing frequency. Different rheological models provide different functional forms of the
        complex compliance function.

    See Also
    --------
    TidalPy.utilities.methods.model.LayerModelHolder
    TidalPy.rheology.Rheology
    """

    default_config = complex_compliance_defaults
    known_models = known_models
    known_model_const_args = known_model_const_args
    known_model_live_args = known_model_live_args
    model_config_key = ('rheology', 'complex_compliance')

    def __init__(self, layer: 'PhysicalLayerType', rheology_class: 'Rheology', model_name: str = None,
                 store_config_in_layer: bool = True, initialize: bool = True):
        """ Constructor for ComplexCompliance class

        Parameters
        ----------
        layer : PhysicalLayerType
            The layer instance which the complex compliance should perform calculations on.
        rheology_class : Rheology
            Rheology class instance where the complex compliance model is stored.
        model_name : str = None
            The user-provided complex compliance name.
        store_config_in_layer: bool = True
            Flag that determines if the final complex compliance model's configuration dictionary should be copied
                into the `layer.config` dictionary.
        initialize : bool = True
            Determines if initial reinit should be performed on the model (loading in data from its `self.config`).
        """

        super().__init__(layer, model_name, store_config_in_layer, initialize=False)

        # Initialized properties
        self._rheology_class = rheology_class

        # State properties
        self._complex_compliances = None

        if initialize:
            self.reinit(initial_init=True)

    def reinit(self, initial_init: bool = False):
        """ Reinit method for the ComplexCompliance class

        Model will look at the user-provided configurations and pull out model information including constants

        Parameters
        ----------
        initial_init : bool = False
            Must be set to `True` if this is the first time this method has been called (additional steps may be
                preformed during the first reinit call).
        """

        super().reinit(initial_init)

    def clear_state(self):
        """ Clears the current state of the class without destroying data set by initialization.
        """

        super().clear_state()

        self._complex_compliances = None

    def _calculate(self) -> Dict[FreqSig, FloatArray]:
        """ Calculate the complex compliance for forcing frequency mode the planet is experiencing.

        Returns
        -------
        complex_compliance : List[np.ndarray]
            Complex compliance of the layer/material [Pa-1]
        """

        if self.tidal_freqs is not None and self.compliance is not None:

            # Check if we need to use the float or array version of the complex compliance calculator
            # OPT: Put this check in the setter for the live args / freqs?
            use_float = True
            for input_ in self.live_inputs:
                if type(input_) == np.ndarray:
                    use_float = False
                    break
            if use_float:
                test_index = list(self.tidal_freqs.keys())[0]
                if type(self.tidal_freqs[test_index]) == np.ndarray:
                    use_float = False
            if use_float:
                comp_func = self.func
            else:
                comp_func = self.func_array

            complex_compliances = compliance_dict_helper(self.tidal_freqs, comp_func,
                                                         self.live_inputs, self.inputs)
            self._complex_compliances = complex_compliances

        return self.complex_compliances


    # # Initialized properties
    @property
    def rheology_class(self) -> 'Rheology':
        """ The rheology class instance where the complex compliance model is stored """
        return self._rheology_class

    @rheology_class.setter
    def rheology_class(self, value):
        raise InitiatedPropertyChangeError


    # # State properties
    @property
    def complex_compliances(self) -> Dict[FreqSig, ComplexArray]:
        """ Complex compliances stored as a dictionary for each unique frequency signature """
        return self._complex_compliances

    @complex_compliances.setter
    def complex_compliances(self, value):
        raise IncorrectMethodToSetStateProperty


    # # Outer-scope properties
    # Rheology class
    @property
    def tidal_freqs(self):
        """ Outer-scope wrapper for rheology.unique_tidal_freqs """
        return self.rheology_class.unique_tidal_frequencies

    @tidal_freqs.setter
    def tidal_freqs(self, value):
        raise OuterscopePropertySetError

    @property
    def compliance(self):
        """ Outer-scope wrapper for rheology.unique_tidal_freqs """

        return self.rheology_class.compliance

    @compliance.setter
    def compliance(self, value):
        raise OuterscopePropertySetError

    @property
    def viscosity(self):
        """ Outer-scope wrapper for rheology.viscosity """
        return self.rheology_class.viscosity

    @viscosity.setter
    def viscosity(self, value):
        raise OuterscopePropertySetError

    @property
    def quality_factor(self):
        """ Outer-scope wrapper for rheology.quality_factor """
        return self.layer.world.tides.fixed_q

    @quality_factor.setter
    def quality_factor(self, value):
        raise OuterscopePropertySetError

    @property
    def beta(self):
        """ Outer-scope wrapper for rheology.beta """
        return self.layer.world.beta

    @beta.setter
    def beta(self, value):
        raise OuterscopePropertySetError