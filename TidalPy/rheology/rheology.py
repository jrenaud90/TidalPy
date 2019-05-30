from functools import partial

from TidalPy.rheology.compliance import ComplianceModelSearcher
from TidalPy.utilities.classes import LayerModel
from . import andrade_frequency_models, compliance_models
from .love_1d import complex_love, complex_love_general, effective_rigidity, effective_rigidity_general
from ..exceptions import ImplementationError, UnknownModelError
from ..types import FloatArray
from .defaults import rheology_param_defaults

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


    def _calculate(self, frequency: FloatArray):
            """ Calculates the Complex Compliance, Effective Rigidity, and Love Number

            :param frequency:     <FloatArray> Tidal Frequency [Rad s-1]
            :return:              <Tuple[FloatArray, dict]> effective_rigidity, {Rheology Name: Complex Compliance,
                                                                                                Complex Love Number
            """


            shear = self.layer.shear_modulus
            visco = self.layer.viscosity
            eff_rigidity = self.calc_effective_rigidity(shear, self.layer.gravity,
                                                        self.layer.radius, self.layer.density)
            compliance = shear**(-1)

            rheology_results = dict()
            for rheology_name, (compliance_func, compliance_input) in self.compliances_byname.items():
                comp_compliance = compliance_func(compliance, visco, frequency, *compliance_input)
                comp_love = self.calc_love(comp_compliance, shear, eff_rigidity)
                rheology_results[rheology_name] = comp_compliance, comp_love

            return effective_rigidity, rheology_results
