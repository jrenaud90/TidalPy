import numpy as np

from TidalPy.constants import R_pluto

# General Inputs
andrade_alpha = 0.3
andrade_zeta = 1.

rock_ref_shear = 5.e10
rock_ref_visco = 1.e20

# Host Inputs
host_radius = R_pluto
host_core_density = 3252.941227949796
host_crust_density = 1000.
host_rock_frac_vol = 0.4486265828438465

# Target Inputs
target_radius = .5093 * R_pluto
target_core_density = 3240.397685371033
target_crust_density = 1000.
target_rock_frac_vol = 0.3759615385144888

# Derivations
host_volume = (4. / 3.) * np.pi * host_radius**3
target_volume = (4. / 3.) * np.pi * target_radius**3
host_core_volume = host_rock_frac_vol * host_volume
target_core_volume = target_rock_frac_vol * target_volume
host_crust_volume = (1. - host_rock_frac_vol) * host_volume
target_crust_volume = (1. - target_rock_frac_vol) * target_volume
host_core_mass = host_core_volume * host_core_density
host_crust_mass = host_crust_volume * host_crust_density
host_mass = host_core_mass + host_crust_mass
target_core_mass = target_core_volume * target_core_density
target_crust_mass = target_crust_volume * target_crust_density
target_mass = target_core_mass + target_crust_mass

def diffeq(time, variables, diff_loop: bool = False):


    semi_major_axis = variables[0, :]
    eccentricity = variables[1, :]
    target_spin_freq = variables[2, :]
    target_mantle_temp = variables[3, :]
    target_crust_visco_dx = variables[4, :]
    target_crust_elas_dx = variables[5, :]
    host_spin_freq = variables[6, :]
    host_mantle_temp = variables[7, :]
    host_crust_visco_dx = variables[8, :]
    host_crust_elas_dx = variables[9, :]

    # Update Rocky Strength and Cooling
    #     - Host
    #     - Target

    # Update Ice Shell Thickness
    #     - Host
    #     - Target

    # Update Radiogenics
    #     - Host
    #     - Target

    # Update Rheology and Tides
    #     - Host
    #     - Target

    # Update Change in Temperature
    #     - Host
    #     - Target

    # Update Change in Thickness
    #     - Host
    #     - Target

    # Update Change in Spin-Rate
    #     - Host
    #     - Target

    # Update Change in Orbit

    # Return
