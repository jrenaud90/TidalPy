import numpy as np

from TidalPy.Extending.burnman import burnman_installed, Material, material_property


class ConstantMaterial(Material):
    """
    Base class for a burnman-like material that only returns constant values (set by a child class).
    """

    def __init__(self):

        if not burnman_installed:
            raise ImportError('Burnman package not found.')

        super().__init__()

    @material_property
    def molar_internal_energy(self):
        return self.params['molar_internal_energy'] * np.ones_like(self.temperature)

    @material_property
    def molar_gibbs(self):
        return self.params['molar_gibbs'] * np.ones_like(self.temperature)

    @material_property
    def molar_helmholtz(self):
        return self.params['molar_helmholtz'] * np.ones_like(self.temperature)

    @material_property
    def molar_mass(self):
        return self.params['molar_mass']

    @material_property
    def molar_volume(self):
        return self.params['molar_volume'] * np.ones_like(self.temperature)

    @material_property
    def density(self):
        return self.params['density'] * np.ones_like(self.temperature)

    @material_property
    def molar_entropy(self):
        return self.params['molar_entropy'] * np.ones_like(self.temperature)

    @material_property
    def molar_enthalpy(self):
        return self.params['molar_enthalpy'] * np.ones_like(self.temperature)

    @material_property
    def isothermal_bulk_modulus(self):
        return self.params['isothermal_bulk_modulus'] * np.ones_like(self.temperature)

    @material_property
    def isothermal_bulk_modulus_reuss(self):
        return self.params['isothermal_bulk_modulus_reuss'] * np.ones_like(self.temperature)

    @material_property
    def isothermal_compressibility(self):
        return self.params['isothermal_compressibility'] * np.ones_like(self.temperature)

    @material_property
    def adiabatic_compressibility(self):
        return self.params['adiabatic_compressibility'] * np.ones_like(self.temperature)

    @material_property
    def shear_modulus(self):
        return self.params['shear_modulus'] * np.ones_like(self.temperature)

    @material_property
    def p_wave_velocity(self):
        return self.params['p_wave_velocity'] * np.ones_like(self.temperature)

    @material_property
    def bulk_sound_velocity(self):
        return self.params['bulk_sound_velocity'] * np.ones_like(self.temperature)

    @material_property
    def shear_wave_velocity(self):
        return self.params['shear_wave_velocity'] * np.ones_like(self.temperature)

    @material_property
    def grueneisen_parameter(self):
        return self.params['grueneisen_parameter'] * np.ones_like(self.temperature)

    @material_property
    def thermal_expansivity(self):
        return self.params['thermal_expansivity'] * np.ones_like(self.temperature)

    @material_property
    def molar_heat_capacity_v(self):
        return self.params['molar_heat_capacity_v'] * np.ones_like(self.temperature)

    @material_property
    def molar_heat_capacity_p(self):
        return self.params['molar_heat_capacity_p'] * np.ones_like(self.temperature)
