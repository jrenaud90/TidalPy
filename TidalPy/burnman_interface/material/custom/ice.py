""" Low and High Pressure Custom Ice Minerals for BurnMan

    Most values are from Sotin et. al 2007, Icarus 191, 337 (For Ice-I and Ice-VII)
"""

low_pressure_ice_eos = 'bm3'
high_pressure_ice_eos = 'mgd3'

from burnman import Mineral
from burnman.processchemistry import dictionarize_formula, formula_mass

from .constant import ConstantMaterial


class LowPressureIceConst(ConstantMaterial):

    def __init__(self):
        self.params = {
            'name'               : 'ice_lp',
            'equation_of_state'  : low_pressure_ice_eos,
            # Following are at room pressure/temperature
            'density'            : 917.,
            'shear_wave_velocity': 1.89e3,
            'p_wave_velocity'    : 3.98e3,
            'thermal_expansivity': 5.0e-5,
            'molar_mass'         : 18.01528 / 1000,  # molar mass in units of [kg/mol]
            'n'                  : 3,  # number of atoms per formula uni
            'latent_heat'        : 3.3e5,
        }

        # Calculated
        self.params['shear_modulus'] = self.params['shear_wave_velocity']**2 * self.params['density']
        self.params['adiabatic_bulk_modulus'] = self.params['p_wave_velocity']**2 * self.params['density'] - \
                                                (4. / 3.) * self.params['shear_modulus']
        self.params['molar_heat_capacity_p'] = 1960 * self.params['molar_mass']

        super().__init__()


class HighPressureIceConst(ConstantMaterial):

    def __init__(self):
        self.params = {
            'name'               : 'ice_hp',
            'equation_of_state'  : high_pressure_ice_eos,
            # Following are at room pressure/temperature
            'density'            : 1600.,
            'shear_wave_velocity': 1.89e3,
            'p_wave_velocity'    : 3.98e3,
            'thermal_expansivity': 1.5e-4,
            'molar_mass'         : 18.01528 / 1000,  # molar mass in units of [kg/mol]
            'n'                  : 3,  # number of atoms per formula unit
            'latent_heat'        : 3.3e5,
        }

        # Calculated
        self.params['shear_modulus'] = self.params['shear_wave_velocity']**2 * self.params['density']
        self.params['adiabatic_bulk_modulus'] = self.params['p_wave_velocity']**2 * self.params['density'] - \
                                                (4. / 3.) * self.params['shear_modulus']
        self.params['molar_heat_capacity_p'] = 1960 * self.params['molar_mass']
        super().__init__()


class UnterbornIce(Mineral):
    """ Water used in Exoplex """
    def __init__(self):
        self.params = {
            'equation_of_state': 'bm4',
            'V_0': 18.797e-6,
            'K_0': 2.06e9,
            'Kprime_0': 6.29,
            'molar_mass': 0.01801528,
            'Kprime_prime_0': (-1.89/2.06e9),
            'n': 1}
        super().__init__()


class Water(Mineral):
    """ Water data from Lide+ (2005) via Sotin+ 2007

    """
    def __init__(self):
        formula = 'H2O'
        formula = dictionarize_formula(formula)
        molar_mass = formula_mass(formula)
        self.params = {
            'equation_of_state': 'bm3',
            'V_0': (molar_mass / 1000.),
            'K_0': 2.2e9,
            'Kprime_0': 4.0,
            'molar_mass': molar_mass,
            'n': sum(formula.values())}
        super().__init__()


class HighPressureIce(Mineral):
    """ High-Pressure Ice (Ice VII) data from Fei+ (1993) via Sotin+ 2007

    """
    def __init__(self):
        formula = 'H2O'
        formula = dictionarize_formula(formula)
        molar_mass = formula_mass(formula)
        self.params = {
            'equation_of_state': 'mgd3',
            'V_0': (molar_mass / 1460.),
            'K_0': 23.9e9,
            'Kprime_0': 4.2,
            'molar_mass': molar_mass,
            'n': sum(formula.values()),
            'Debye_0': 1470,
            'grueneisen_0': 1.2,
            'q_0': 1.0}
        super().__init__()


class IceX_Fu2010(Mineral):
    """ Very High-Pressure Ice (Ice X) data from Loubeyre+ (1999) via Fu+ (2010)

    """

    def __init__(self):
        formula = 'H2O'
        formula = dictionarize_formula(formula)
        molar_mass = formula_mass(formula)
        self.params = {
            'equation_of_state': 'bm3',
            'V_0': (molar_mass / 1239.),
            'K_0': 4.26e9,
            'Kprime_0': 7.75,
            'molar_mass': molar_mass,
            'n': sum(formula.values())}
        super().__init__()


class IceVII_Fu2010(Mineral):
    """ High-Pressure Ice (Ice VII) data from Hemley+ (1987) via Fu+ (2010)

    """

    def __init__(self):
        formula = 'H2O'
        formula = dictionarize_formula(formula)
        molar_mass = formula_mass(formula)
        self.params = {
            'equation_of_state': 'bm3',
            'V_0'              : (molar_mass / 1463.),
            'K_0'              : 23.7e9,
            'Kprime_0'         : 4.15,
            'molar_mass'       : molar_mass,
            'n'                : sum(formula.values())
        }
        super().__init__()

class IceIh_Fu2010(Mineral):
    """ Low-Pressure Ice (Ice Ih) data from Weast (1969) and Strassle+ (2005) via Fu+ (2010)

    """

    def __init__(self):
        formula = 'H2O'
        formula = dictionarize_formula(formula)
        molar_mass = formula_mass(formula)
        self.params = {
            'equation_of_state': 'bm3',
            'V_0'              : (molar_mass / 917.),
            'K_0'              : 9.86e9,
            'Kprime_0'         : 6.6,
            'molar_mass'       : molar_mass,
            'n'                : sum(formula.values())
        }
        super().__init__()
