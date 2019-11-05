import copy
from functools import partial
from typing import Tuple

import numpy as np
from . import calculate_tides
from ..exceptions import (AttributeNotSetError, ImplementationException, ParameterMissingError, UnknownModelError)
from ..rheology import andrade_frequency_models, compliance_models
from ..rheology.compliance import ComplianceModelSearcher
from ..rheology.defaults import rheology_param_defaults
from ..types import FloatArray
from ..utilities.model import LayerModel


FAKE_FREQS = (np.asarray([0.]),)


class Tides(LayerModel):
    default_config = copy.deepcopy(rheology_param_defaults)
    config_key = 'rheology'

    def __init__(self, layer):

        super().__init__(layer=layer, function_searcher=None, call_reinit=True)

        # Override model if layer has tidal off
        if not self.layer.config['is_tidal']:
            self.model = 'off'

        # Pull out information about the planet/layer
        self.config['planet_beta'] = self.layer.world.beta
        self.config['quality_factor'] = self.layer.world.fixed_q

        # Setup Love number calculator
        self.calc_love = None
        self.calc_effective_rigidity = None
        self.order_l = None
        self.full_calculation = True
        if self.model == '2order':
            self.order_l = 2
        elif self.model == 'general':
            self.order_l = self.config['order_l']
        elif self.model == 'off':
            self.full_calculation = False
        elif self.model == 'multilayer':
            # TODO: Multilayer code
            raise ImplementationException
        else:
            raise UnknownModelError

        if self.order_l > 2:
            raise ImplementationException

        # Setup complex compliance calculator
        model_searcher = ComplianceModelSearcher(compliance_models, andrade_frequency_models, self.config,
                                                 defaults_require_key=False)

        if not self.full_calculation:
            self.rheology_name = 'off'
        else:
            self.rheology_name = self.config['compliance_model']

        comp_func, comp_input, comp_live_args = model_searcher.find_model(self.rheology_name)
        # TODO: As of 0.1.0a no compliance functions have live args, so comp_live_args is never used.
        self.compliance = comp_func
        self.compliance_inputs = comp_input

        # Switches
        self.is_spin_sync = self.layer.is_spin_sync

    def _calculate(self) -> Tuple[FloatArray, Tuple[FloatArray], Tuple[FloatArray]]:
        """ Calculates the Complex Compliance, Effective Rigidity, and Love Number

        """

        # Physical and Frequency Data
        try:
            shear = self.layer.shear_modulus
            visco = self.layer.viscosity
            orbital_freq = self.layer.orbital_freq
            spin_freq = self.layer.spin_freq
        except AttributeNotSetError:
            raise ParameterMissingError
        if shear is None:
            raise ParameterMissingError
        if visco is None:
            raise ParameterMissingError
        if orbital_freq is None:
            raise ParameterMissingError
        if spin_freq is None:
            if self.is_spin_sync:
                spin_freq = orbital_freq
            else:
                raise ParameterMissingError

        # Make a call to the calculate tides function
        calculate_tides


        eff_rigidity = self.calc_effective_rigidity(shear, self.layer.gravity,
                                                    self.layer.radius, self.layer.density)  # type: FloatArray
        compliance = shear**(-1)
        if self.full_calculation:
            tidal_freqs = self.layer.tidal_freqs
            if tidal_freqs is None:
                raise ParameterMissingError
        else:
            tidal_freqs = FAKE_FREQS

        # Tides Functions
        complex_compliance_tupl = list()
        complex_love_tupl = list()
        for freq in tidal_freqs:
            complex_compliance = self.compliance(compliance, visco, freq, *self.compliance_inputs)  # type: FloatArray
            complex_love = self.calc_love(complex_compliance, shear, eff_rigidity)  # type: FloatArray
            complex_love_tupl.append(np.nan_to_num(complex_love))
            complex_compliance_tupl.append(np.nan_to_num(complex_compliance))

        return eff_rigidity, tuple(complex_compliance_tupl), tuple(complex_love_tupl)
