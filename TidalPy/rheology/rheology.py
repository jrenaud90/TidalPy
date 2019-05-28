from numba import njit

from TidalPy.rheology.compliance import ComplianceModelSearcher
from TidalPy.structures.layers import LayerType
from TidalPy.utilities.classes import ModelHolder, LayerModel
from ..types import FloatArray
from . import compliance_models, andrade_frequency_models
from .love_1d import effective_rigidity, complex_love, complex_love_general, effective_rigidity_general
from ..exceptions import ImplementationError, UnknownModelError
from functools import partial

# Rheology (compliance) Model Finder
rheology_param_defaults = {
    'ice': {
            'model': '2order',
            'solid_viscosity': {
                'model': 'arrhenius',
                # Matching Moore2006 for Volume Diffusion
                'arrhenius_coeff': 9.06e-8**(-1),
                'additional_temp_dependence': True,
                'stress': 1.0,
                'stress_expo': 1.0,
                'grain_size': 5.0e-4,
                'grain_size_expo': 2.0,
                'molar_activation_energy': 59.4e3,
                'molar_activation_volume': -1.3e-5
            },
            'liquid_viscosity': {
                'model': 'reference',
                'reference_viscosity': 0.89e-3,
                'reference_temperature': 25.0 + 273.15,
                'molar_activation_energy': 1.62e4,
                'molar_activation_volume': 0.0
            },
            'order_l': 2,
            'compliances': ['Maxwell', 'Andrade'],
            'voigt_compliance_offset': .2,
            'voigt_viscosity_offset': .02,
            'alpha': .3,
            'zeta': 1.,
            'andrade_frequency_model': 'exponential',
            'andrade_critical_freq': 2.0e-7
        },
    'rock': {
            'model': '2order',
            'solid_viscosity': {
                'model': 'reference',
                'reference_viscosity': 1.0e22,
                'reference_temperature': 1000.0,
                'molar_activation_energy': 300000.0,
                'molar_activation_volume': 0.0
            },
            'liquid_viscosity': {
                'model': 'reference',
                'reference_viscosity': 0.2,
                'reference_temperature': 2000.0,
                'molar_activation_energy': 6.64e-20,
                'molar_activation_volume': 0.0
            },
            'order_l': 2,
            'compliances': ['Maxwell', 'Andrade'],
            'voigt_compliance_offset': .2,
            'voigt_viscosity_offset': .02,
            'alpha': .3,
            'zeta': 1.,
            'andrade_frequency_model': 'exponential',
            'andrade_critical_freq': 2.0e-7
        },
    'iron': {
            'model': 'off',
            # These are just placeholders. Right now we do not calculate the tides in iron layer, so no need to have
            #    correct values at the moment.
            'solid_viscosity': {
                'model': 'constant',
                'reference_viscosity': 1.0e20,
            },
            'liquid_viscosity': {
                # These values match Wijs et al 1998 (their work actually does not show much change in the liquid visc
                #    at Earth's core pressure, so a constant model may not be too incorrect).
                'model': 'constant',
                'reference_viscosity': 1.3e-2,
            },
        }
    }


class Rheology(LayerModel):

    default_config = rheology_param_defaults
    config_key = 'rheology'

    def __init__(self, layer: LayerType):

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
