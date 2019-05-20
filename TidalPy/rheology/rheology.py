from numba import njit

from TidalPy.rheology.compliance import ComplianceModelSearcher
from TidalPy.utilities.classes import ModelHolder
from TidalPy.utilities.search import ModelSearcher
from ..types import FloatArray
from . import compliance_models, andrade_frequency_models

# Rheology (compliance) Model Finder
rheology_param_defaults = {
    'ice': {
            'model': '2order',
            'order_l': 2,
            'compliances': ['Maxwell', 'Andrade'],
            'voigt_compliance_offset': .2,
            'voigt_viscosity_offset': .02,
            'alpha': .3,
            'zeta': 1.,
            'andrade_frequency_model': 'exponential',
            'andrade_critical_freq': 2.0e-7
        },
    'rocky': {
            'model': '2order',
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
            'model': 'off'
        }
    }

@njit
def complex_love(complex_compliance: FloatArray, shear_modulus: FloatArray, eff_rigidity: FloatArray):
    """ Calculates the 2nd order complex Love number

    :param complex_compliance: <FloatArray> Complex compliance (rheology based) [Pa-1]
    :param shear_modulus:      <FloatArray> Temperature modulated rigidity [Pa]
    :param eff_rigidity:       <FloatArray> 2nd order effective rigidity
    :return:                   <FloatArray> Complex Love Number
    """

    return (3. / 2.) * (1. + eff_rigidity / (shear_modulus * complex_compliance))**(-1)


@njit
def complex_love_general(complex_compliance: FloatArray, shear_modulus: FloatArray, eff_rigidity_general: FloatArray,
                         order_l: int = 2):
    """ Calculates the l-th order complex Love number

    :param complex_compliance:      <FloatArray> Complex compliance (rheology based) [Pa-1]
    :param shear_modulus:           <FloatArray> Temperature modulated rigidity [Pa]
    :param eff_rigidity_general:    <FloatArray> l-th order effective rigidity
    :param order_l:                 <int> (optional) Outer-most fourier summation index
    :return:                        <FloatArray> Complex Love Number
    """

    return (3. / (2. * (order_l - 1.))) * (1. + eff_rigidity_general / (shear_modulus * complex_compliance))**(-1)

@njit
def effective_rigidity(shear_modulus: FloatArray, gravity: float, radius: float, density: float):
    """ Calculates the 2nd order effective rigidity

    :param shear_modulus: <FloatArray> Temperature modulated rigidity
    :param gravity:       <float> Surface gravity [m s-2]
    :param radius:        <float> Surface radius [m]
    :param density:       <float> Bulk density [kg m-3]
    :return:              <FloatArray> 2nd order Effective Rigidity
    """

    return (19. / 2.) * shear_modulus / (gravity * radius * density)

@njit
def effective_rigidity_general(shear_modulus: FloatArray, gravity: float, radius: float, density: float,
                               order_l: int = 2):
    """ Calculates the l-th order effective rigidity

    :param shear_modulus: <FloatArray> Temperature modulated rigidity
    :param gravity:       <float> Surface gravity [m s-2]
    :param radius:        <float> Surface radius [m]
    :param density:       <float> Bulk density [kg m-3]
    :param order_l:       <int> (optional) Outer-most fourier summation index
    :return:              <FloatArray> 2nd order Effective Rigidity
    """

    return (2. * order_l**2 + 4. * order_l + 3. / order_l) * shear_modulus / (gravity * radius * density)


class Rheology(ModelHolder):

    name = 'Rheology'

    def __init__(self, layer_type: str, rheology_config: dict = None):

        model_searcher = ComplianceModelSearcher(compliance_models, andrade_frequency_models,
                                                 rheology_param_defaults[layer_type])
        model_searcher.user_parameters = rheology_config

        super().__init__(user_config=rheology_config, default_config=rheology_param_defaults[layer_type],
                         function_searcher=model_searcher, automate=True)

        if self.model == 'off':
            # No tidal heating
