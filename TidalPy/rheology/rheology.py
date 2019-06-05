from functools import partial

from TidalPy.rheology.compliance import ComplianceModelSearcher
from TidalPy.utilities.classes import LayerModel
from . import andrade_frequency_models, compliance_models
from .love_1d import complex_love, complex_love_general, effective_rigidity, effective_rigidity_general
from ..exceptions import ImplementationError, UnknownModelError
from ..types import FloatArray
from .defaults import rheology_param_defaults
from typing import List

class Rheology(LayerModel):

    default_config = rheology_param_defaults
    config_key = 'rheology'

    def __init__(self, layer):

        super().__init__(layer=layer, function_searcher=None, automate=True)

        # Setup Love number calculator
        self.calc_love = None
        self.calc_effective_rigidity = None
        self.order_l = None
        if self.model == '2order':
            self.order_l = 2
            self.calc_love = complex_love
            self.calc_effective_rigidity = effective_rigidity
        elif self.model == 'general':
            self.order_l = self.config['order_l']
            self.calc_love = partial(complex_love_general, order_l=self.order_l)
            self.calc_effective_rigidity = partial(effective_rigidity_general, order_l=self.order_l)
        elif self.model == 'multilayer':
            # TODO: Multilayer code
            raise ImplementationError
        else:
            raise UnknownModelError

        # Setup complex compliance calculator
        model_searcher = ComplianceModelSearcher(compliance_models, andrade_frequency_models, self.config)

        self.compliances = list()
        self.compliance_inputs = list()
        for rheology_name in self.config['compliances']:
            comp_func, comp_input = model_searcher.find_model(rheology_name)
            self.compliances.append(comp_func)
            self.compliance_inputs.append(comp_input)
        self.compliances = tuple(self.compliances)
        self.compliance_inputs = tuple(self.compliance_inputs)
        self.compliances_byname = {rheology_name: (comp_func, comp_input) for rheology_name, comp_func, comp_input in
                                   zip(self.config['compliances'], self.compliances, self.compliance_inputs)}


    def _calculate(self):
            """ Calculates the Complex Compliance, Effective Rigidity, and Love Number

            """


            shear = self.layer.shear_modulus
            visco = self.layer.viscosity
            eff_rigidity = self.calc_effective_rigidity(shear, self.layer.gravity,
                                                        self.layer.radius, self.layer.density)
            compliance = shear**(-1)

            rheology_compliances = dict()
            rheology_loves = dict()
            for rheology_name, (compliance_func, compliance_input) in self.compliances_byname.items():
                compliance_mode_results = list()
                love_mode_results = list()
                for mode in self.layer.tidal_modes:
                    comp_compliance = compliance_func(compliance, visco, mode, *compliance_input)
                    love_mode_results.append(self.calc_love(comp_compliance, shear, eff_rigidity))
                    compliance_mode_results.append(comp_compliance)
                rheology_compliances[rheology_name] = tuple(compliance_mode_results)
                rheology_loves[rheology_name] = tuple(love_mode_results)

            return effective_rigidity, rheology_compliances, rheology_loves
