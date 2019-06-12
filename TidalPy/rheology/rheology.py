from functools import partial
import numpy as np
from TidalPy.rheology.compliance import ComplianceModelSearcher
from TidalPy.utilities.classes import LayerModel
from . import andrade_frequency_models, compliance_models
from .love_1d import complex_love, complex_love_general, effective_rigidity, effective_rigidity_general
from ..exceptions import ImplementationError, UnknownModelError, IncompatibleModelError, ParameterMissingError
from ..types import FloatArray
from .defaults import rheology_param_defaults
from typing import List, Tuple

FAKE_MODES = (np.asarray(0.), )

class Rheology(LayerModel):

    default_config = rheology_param_defaults
    config_key = 'rheology'

    def __init__(self, layer):

        super().__init__(layer=layer, function_searcher=None, automate=True)

        # Override model if layer has tidal off
        if not self.layer.config['is_tidal']:
            self.model = 'off'

        # Setup Love number calculator
        self.calc_love = None
        self.calc_effective_rigidity = None
        self.order_l = self.config['order_l']
        self.full_calculation = True
        if self.model == '2order':
            self.calc_love = complex_love
            self.calc_effective_rigidity = effective_rigidity
        elif self.model == 'general':
            self.calc_love = partial(complex_love_general, order_l=self.order_l)
            self.calc_effective_rigidity = partial(effective_rigidity_general, order_l=self.order_l)
        elif self.model == 'off':
            self.calc_love = partial(complex_love_general, order_l=self.order_l)
            self.calc_effective_rigidity = partial(effective_rigidity_general, order_l=self.order_l)
            self.full_calculation = False
        elif self.model == 'multilayer':
            # TODO: Multilayer code
            raise ImplementationError
        else:
            raise UnknownModelError

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

    def _calculate(self) -> Tuple[np.ndarray, Tuple[np.ndarray], Tuple[np.ndarray]]:
        """ Calculates the Complex Compliance, Effective Rigidity, and Love Number

        """

        # Physical and Frequency Data
        shear = self.layer.shear_modulus
        visco = self.layer.viscosity
        if shear is None:
            raise ParameterMissingError
        if visco is None:
            raise ParameterMissingError

        eff_rigidity = self.calc_effective_rigidity(shear, self.layer.gravity,
                                                    self.layer.radius, self.layer.density)
        compliance = shear**(-1)
        if self.full_calculation:
            tidal_modes = self.layer.tidal_modes
            if tidal_modes is None:
                raise ParameterMissingError
        else:
            tidal_modes = FAKE_MODES

        # Rheology Functions
        complex_compliance_tupl = list()
        complex_love_tupl = list()
        for mode in tidal_modes:
            complex_compliance = self.compliance(compliance, visco, mode, *self.compliance_inputs)
            complex_love_tupl.append(self.calc_love(complex_compliance, shear, eff_rigidity))
            complex_compliance_tupl.append(complex_compliance)

        return effective_rigidity, tuple(complex_compliance_tupl), tuple(complex_love_tupl)
