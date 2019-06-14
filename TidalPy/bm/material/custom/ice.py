""" Low and High Pressure Custom Ice Minerals for BurnMan

    Most values are from Sotin et. al 2007, Icarus 191, 337 (For Ice-I and Ice-VII)
"""

low_pressure_ice_eos = 'bm3'
high_pressure_ice_eos = 'mgd3'



import burnman

from .constant import ConstantMaterial


class LowPressureIceConst(ConstantMaterial):

    def __init__(self):
        self.params = {
            'name'             : 'ice_lp',
            'equation_of_state': low_pressure_ice_eos,
            # Following are at room pressure/temperature
            'density'          : 917.,
            'shear_wave_velocity': 1.89e3,
            'p_wave_velocity': 3.98e3,
            'thermal_expansivity': 5.0e-5,
            'molar_mass'       : 18.01528/1000,  # molar mass in units of [kg/mol]
            'n'                : 3,  # number of atoms per formula uni
            'latent_heat'          : 3.3e5,
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



# TODO: Actual equation of state
class LowPressureIce(burnman.Mineral):

    def __init__(self):
        self.params = {
            'name'             : 'ice_lp',
            'equation_of_state': low_pressure_ice_eos,
            # Following are at room pressure/temperature
            'V_0'              : 10.844e-6,  # Molar volume [m^3/(mole molecules)]
            'K_0'              : 135.19e9,  # Reference bulk modulus [Pa]
            'Kprime_0'         : 6.04,  # pressure derivative of bulk modulus
            'G_0'              : 175.0e9,  # reference shear modulus
            'Gprime_0'         : 1.7,  # pressure derivative of shear modulus
            'molar_mass'       : .055845,  # molar mass in units of [kg/mol]
            'n'                : 1,  # number of atoms per formula unit
            'Debye_0'          : 998.85,  # Debye temperature for material.
            'grueneisen_0'     : 1.368,  # Gruneisen parameter for material.
            'q_0'              : 0.917,  # isotropic strain derivative of gruneisen
            'eta_s_0'          : 3.0  # full strain derivative of gruneisen parameter
        }
        burnman.Mineral.__init__(self)


class HighPressureIce(burnman.Mineral):

    def __init__(self):
        self.params = {
            'name'             : 'ice_hp',
            'equation_of_state': high_pressure_ice_eos,
            # Following are at room pressure/temperature
            'V_0'              : 10.844e-6,  # Molar volume [m^3/(mole molecules)]
            ##'K_0'              : 23.9e9,  # Reference bulk modulus [Pa]
            'Kprime_0'         : 6.04,  # pressure derivative of bulk modulus
            'G_0'              : 175.0e9,  # reference shear modulus
            'Gprime_0'         : 1.7,  # pressure derivative of shear modulus
            'molar_mass'       : .055845,  # molar mass in units of [kg/mol]
            'n'                : 1,  # number of atoms per formula unit
            ##'Debye_0'          : 1470.0,  # Debye temperature for material.
            'grueneisen_0'     : 1.368,  # Gruneisen parameter for material.
            'q_0'              : 0.917,  # isotropic strain derivative of gruneisen
            'eta_s_0'          : 3.0  # full strain derivative of gruneisen parameter
        }
        burnman.Mineral.__init__(self)