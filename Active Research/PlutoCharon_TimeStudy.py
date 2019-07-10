import time
from functools import partial
from numba import njit

import numpy as np
from scipy.constants import R
from scipy.integrate import solve_ivp

import TidalPy
from TidalPy import log
from TidalPy.dynamics import spin_rate_derivative
from TidalPy.dynamics.duel_dissipation import eccentricity_derivative, semi_major_axis_derivative
from TidalPy.orbit.modes import nsr_modes
from TidalPy.orbit.orbit import OrbitBase as Orbit
from TidalPy.radiogenics import radiogenic_isotope, standard_isotope_input
from TidalPy.rheology.compliance_models import andrade, maxwell, off
from TidalPy.thermal import convection
from TidalPy.thermal.viscosity_models import reference as calc_reference_viscosity
from TidalPy.tides import calc_tides
from TidalPy.utilities import progress_bar
from TidalPy.utilities.conversions import convert_to_hms, myr2sec, sec2myr, semi_a2orbital_motion
from TidalPy.utilities.numpy_help import neg_array_for_log_plot

# Planet Building
target = TidalPy.build_planet('charon', force_build=True)
host = TidalPy.build_planet('pluto', force_build=True)
star = TidalPy.build_planet('sol', force_build=True)
orbit = Orbit(star, host, target, duel_dissipation=True, time_study=True)

# Initial Conditions
time_span = myr2sec(np.asarray((0., 2000)))
modern_semi_major_axis = 1.9599e7
modern_orb_freq = semi_a2orbital_motion(modern_semi_major_axis, host.mass, target.mass)

semi_major_axis_init = 1.5 * modern_semi_major_axis
eccentricity_init = 0.2
spin_freq_host_init = modern_orb_freq * 30
temperature_core_host_init = 1000.
thickness_visco_host_init = 0.1 * host.crust.thickness
thickness_elast_host_init = 0.1 * host.crust.thickness
spin_freq_targ_init = modern_orb_freq * 30
temperature_core_targ_init = 1000.
thickness_visco_targ_init = 0.1 * target.crust.thickness
thickness_elast_targ_init = 0.1 * target.crust.thickness

# Other Conditions
melting_temp_ice = 273.15
ocean_temp = melting_temp_ice

# Switches
lock_to_1to1 = True
lock_zero_eccen = True
lock_ocean_freeze = True

# Integration Parameters
rtol = 1.0e-5
method = 'LSODA'
use_nocore = True

if use_nocore:
    initial_conditions = [
        semi_major_axis_init,
        eccentricity_init,
        spin_freq_host_init,
        thickness_visco_host_init,
        thickness_elast_host_init,
        spin_freq_targ_init,
        thickness_visco_targ_init,
        thickness_elast_targ_init
    ]
else:
    initial_conditions = [
        semi_major_axis_init,
        eccentricity_init,
        spin_freq_host_init,
        temperature_core_host_init,
        thickness_visco_host_init,
        thickness_elast_host_init,
        spin_freq_targ_init,
        temperature_core_targ_init,
        thickness_visco_targ_init,
        thickness_elast_targ_init
    ]
initial_conditions = tuple([np.asarray(init_cond) for init_cond in initial_conditions])


# Rocky Material Inputs
specific_heat_rock = 1.2255e3
thermal_expansion_rock = 5.0e-5
thermal_cond_rock = 4.0
convection_alpha_rock = 1.0
convection_beta_rock = 1. / 3.
critical_rayleigh_rock = 1000.
viscosity_rock_ref = 1.e15
viscosity_ref_temp = 1600.
activation_eng_rock = 360.0e3
andrade_alpha_rock = 0.3
andrade_zeta_rock = 1.

# Icy Material Inputs (These are for the dissipating zone!)
density_ice = 1000.
specific_heat_ice = 1925.0
thermal_expansion_ice = 1.56e-4
thermal_cond_ice = 2.27
thermal_diff_ice = thermal_cond_ice / (density_ice * specific_heat_ice)
latent_heat_ice = 284.0e3
activation_eng_ice = 58.4e3
convection_alpha_ice = 1.0
convection_beta_ice = 1. / 3.
critical_rayleigh_ice = 900.
shear_modulus_crust_host = 3.3e9
shear_modulus_crust_targ = shear_modulus_crust_host
viscosity_crust_host = np.asarray([5.e13])
viscosity_crust_targ = viscosity_crust_host.copy()
crust_viscosity_drop = 10.
andrade_alpha_ice = 0.3
andrade_zeta_ice = 1.

# Orbital Inputs
static_inclination = np.asarray([np.deg2rad(0.)])

# Calculated Inputs
density_rock_host = host.core.density_bulk
density_rock_targ = target.core.density_bulk
thermal_diff_rock_host = thermal_cond_rock / (density_rock_host * specific_heat_rock)
thermal_diff_rock_targ = thermal_cond_rock / (density_rock_targ * specific_heat_rock)
crust_viscosity_drop_log = np.log(crust_viscosity_drop)
crust_viscosity_l = activation_eng_ice / (R * melting_temp_ice)
ice_top_temp_host = melting_temp_ice / ((crust_viscosity_drop_log / crust_viscosity_l) + 1.)
ice_top_temp_targ = ice_top_temp_host

# Overshoot Values
MIN_SEMI_A = (target.radius + host.radius) * 1.1
MIN_ECCEN = 1.e-10
log(f'Min Semi-major Axis [frac of modern]: {MIN_SEMI_A/modern_semi_major_axis:.03f}', level='info')

# Store values for njit optimizations
thickness_crust_host = host.crust.thickness
thickness_crust_target = target.crust.thickness
radius_host = host.radius
radius_target = target.radius
radius_core_host = host.core.radius
radius_core_target = target.core.radius
mass_host = host.mass
mass_target = target.mass
mass_core_host = host.core.mass
mass_core_target = target.core.mass
volume_host = host.volume
volume_target = target.volume
gravity_surf_host = host.gravity_surf
gravity_surf_target = target.gravity_surf
tidal_susceptibility_inflated_host = host.tidal_susceptibility_inflated
tidal_susceptibility_inflated_target = target.tidal_susceptibility_inflated
moi_host = host.moi
moi_target = target.moi
surface_temp_host = host.surface_temperature
surface_temp_targ = target.surface_temperature

@progress_bar(time_span, verbose=True)
def diffeq_withcore(time, variables, complex_compliance_func, complex_compliance_input, diff_loop):
    semi_major_axis = variables[0, :]
    eccentricity = variables[1, :]
    spin_freq_host = variables[2, :]
    temperature_core_host = variables[3, :]
    thickness_visco_host = variables[4, :]
    thickness_elast_host = variables[5, :]
    spin_freq_targ = variables[6, :]
    temperature_core_targ = variables[7, :]
    thickness_visco_targ = variables[8, :]
    thickness_elast_targ = variables[9, :]
    array_shape = semi_major_axis.shape

    # Helper Values
    thickness_ice_host = thickness_visco_host + thickness_elast_host
    thickness_ice_targ = thickness_visco_targ + thickness_elast_targ
    overshoot_index_host = np.copy(thickness_ice_host > 0.99 * host.crust.thickness)
    overshoot_index_targ = np.copy(thickness_ice_targ > 0.99 * target.crust.thickness)

    # Check for value logic
    eccentricity[eccentricity < 1.e-10] = 0.
    semi_major_axis[semi_major_axis < MIN_SEMI_A] = MIN_SEMI_A
    temperature_core_host[temperature_core_host < 273.15] = 273.15
    temperature_core_targ[temperature_core_targ < 273.15] = 273.15
    # For now, any time the ice becomes too thick break it up into 75% visco, 24% elastic, and 1% ocean
    thickness_elast_host[overshoot_index_host] = 0.24 * host.crust.thickness
    thickness_visco_host[overshoot_index_host] = 0.75 * host.crust.thickness
    thickness_elast_targ[overshoot_index_targ] = 0.24 * target.crust.thickness
    thickness_visco_targ[overshoot_index_targ] = 0.75 * target.crust.thickness

    # Conversions
    time_myr = sec2myr(time)
    orbital_freq = semi_a2orbital_motion(semi_major_axis, host.mass, target.mass)

    # Find Tidal Modes
    #     - Host
    modes_host, freqs_host, heating_coeffs_host, torque_coeffs_host = nsr_modes(orbital_freq, spin_freq_host,
                                                                                eccentricity, static_inclination)
    #     - Target
    modes_targ, freqs_targ, heating_coeffs_targ, torque_coeffs_targ = nsr_modes(orbital_freq, spin_freq_targ,
                                                                                eccentricity, static_inclination)

    # Calc Tidal Susceptibility
    #     - Host
    tidal_suscept_host = host.tidal_susceptibility_inflated / semi_major_axis**6
    #     - Target
    tidal_suscept_targ = target.tidal_susceptibility_inflated / semi_major_axis**6

    # Update Rocky Strength
    #     - Host
    viscosity_core_host = calc_reference_viscosity(temperature_core_host, 0., viscosity_rock_ref, viscosity_ref_temp,
                                                   activation_eng_rock, 0.)
    #     - Target
    viscosity_core_targ = calc_reference_viscosity(temperature_core_targ, 0., viscosity_rock_ref, viscosity_ref_temp,
                                                   activation_eng_rock, 0.)

    # Update Ice Shell Thickness
    #     - Host
    radius_elast_host = host.radius
    radius_visco_host = radius_elast_host - thickness_elast_host
    radius_ocean_host = radius_visco_host - thickness_visco_host
    volume_visco_host = (4. / 3.) * np.pi * (radius_visco_host**3 - radius_ocean_host**3)
    volume_ocean_host = (4. / 3.) * np.pi * (radius_ocean_host**3 - host.core.radius**3)

    #     - Target
    radius_elast_targ = target.radius
    radius_visco_targ = radius_elast_targ - thickness_elast_targ
    radius_ocean_targ = radius_visco_targ - thickness_visco_targ
    volume_visco_targ = (4. / 3.) * np.pi * (radius_visco_targ**3 - radius_ocean_targ**3)
    volume_ocean_targ = (4. / 3.) * np.pi * (radius_ocean_targ**3 - target.core.radius**3)

    # Where do tides concentrate?
    tidalvolfrac_crust_host = volume_visco_host / host.volume
    tidalvolfrac_crust_targ = volume_visco_targ / target.volume

    # Calculate Core Cooling
    #     - Host
    delta_temp_core_host = temperature_core_host - ocean_temp
    q_core_conv_host, blt_core_host, rayleigh_core_host, neusselt_core_host = convection(
            delta_temp_core_host, host.core.thickness, thermal_cond_rock, viscosity_core_host,
            thermal_diff_rock_host, thermal_expansion_rock, host.core.gravity_surf, host.core.density_bulk,
            convection_alpha_rock, convection_beta_rock, critical_rayleigh_rock)
    core_cooling_host = q_core_conv_host * 4. * np.pi * host.core.radius**2
    #     - Target
    delta_temp_core_targ = temperature_core_targ - ocean_temp
    q_core_conv_targ, blt_core_targ, rayleigh_core_targ, neusselt_core_targ = convection(
            delta_temp_core_targ, target.core.thickness, thermal_cond_rock, viscosity_core_targ,
            thermal_diff_rock_targ, thermal_expansion_rock, target.core.gravity_surf, target.core.density_bulk,
            convection_alpha_rock, convection_beta_rock, critical_rayleigh_rock)
    core_cooling_targ = q_core_conv_targ * 4. * np.pi * target.core.radius**2

    # Calculate Ice Cooling
    #     - Host
    delta_temp_crust_host = melting_temp_ice - ice_top_temp_host
    q_crust_conv_host, blt_crust_host, rayleigh_crust_host, neusselt_crust_host = convection(
            delta_temp_crust_host, thickness_visco_host, thermal_cond_ice, viscosity_crust_host,
            thermal_diff_ice, thermal_expansion_ice, host.gravity_surf, density_ice,
            convection_alpha_ice, convection_beta_ice, critical_rayleigh_ice)
    crust_cooling_host = q_crust_conv_host * 4. * np.pi * radius_visco_host**2
    #     - Target
    delta_temp_crust_targ = melting_temp_ice - ice_top_temp_targ
    q_crust_conv_targ, blt_crust_targ, rayleigh_crust_targ, neusselt_crust_targ = convection(
            delta_temp_crust_targ, thickness_visco_targ, thermal_cond_ice, viscosity_crust_targ,
            thermal_diff_ice, thermal_expansion_ice, target.gravity_surf, density_ice,
            convection_alpha_ice, convection_beta_ice, critical_rayleigh_ice)
    crust_cooling_targ = q_crust_conv_targ * 4. * np.pi * radius_visco_targ**2

    # Update Radiogenics
    #     - Host
    radio_heating_core_host = radiogenic_isotope(time_myr, host.core.mass, *standard_isotope_input)
    #     - Target
    radio_heating_core_targ = radiogenic_isotope(time_myr, target.core.mass, *standard_isotope_input)

    # Calculate Compliance and Tides
    #     - Host
    cmplx_comp_crust_bymode_host = [
        complex_compliance_func(shear_modulus_crust_host**-1, viscosity_crust_host, tidal_freq,
                                *complex_compliance_input)
        for tidal_freq in freqs_host
    ]

    effective_rigid_crust_host, tidal_heating_crust_host, ztorque_crust_host, _, _ = \
        calc_tides(host.gravity_surf, host.radius, host.crust.density_bulk, shear_modulus_crust_host,
                   tidal_suscept_host,
                   cmplx_comp_crust_bymode_host, heating_coeffs_host, torque_coeffs_host)

    # Scale down by the volume factor
    tidal_heating_crust_host *= tidalvolfrac_crust_host
    ztorque_crust_host *= tidalvolfrac_crust_host

    # TODO: In this version, tides are restricted to the ice layer only
    tidal_heating_core_host = 0.
    ztorque_core_host = 0.

    tidal_heating_host = tidal_heating_crust_host + tidal_heating_core_host
    ztorque_host = ztorque_crust_host + ztorque_core_host

    #     - Target
    cmplx_comp_crust_bymode_targ = [
        complex_compliance_func(shear_modulus_crust_targ**-1, viscosity_crust_targ, tidal_freq,
                                *complex_compliance_input)
        for tidal_freq in freqs_targ
    ]

    effective_rigid_crust_targ, tidal_heating_crust_targ, ztorque_crust_targ, _, _ = \
        calc_tides(target.gravity_surf, target.radius, target.crust.density_bulk, shear_modulus_crust_targ,
                   tidal_suscept_targ,
                   cmplx_comp_crust_bymode_targ, heating_coeffs_targ, torque_coeffs_targ)

    # Scale down by the volume factor
    tidal_heating_crust_targ *= tidalvolfrac_crust_targ
    ztorque_crust_targ *= tidalvolfrac_crust_targ

    # TODO: In this version, tides are restricted to the ice layer only
    tidal_heating_core_targ = 0.
    ztorque_core_targ = 0.

    tidal_heating_targ = tidal_heating_crust_targ + tidal_heating_core_targ
    ztorque_targ = ztorque_crust_targ + ztorque_core_targ

    # Update Change in Core Temperature
    #     - Host
    change_core_temp_host = (tidal_heating_core_host + radio_heating_core_host - core_cooling_host) / \
                            (density_rock_host * specific_heat_rock * host.core.volume)
    #     - Target
    change_core_temp_targ = (tidal_heating_core_targ + radio_heating_core_targ - core_cooling_targ) / \
                            (density_rock_targ * specific_heat_rock * target.core.volume)

    # Find Heat Fluxes
    #     - Host
    delta_temp_elast_host = melting_temp_ice - ice_top_temp_host
    q_visco_in_host = (core_cooling_host + tidal_heating_crust_host) / (4. * np.pi * radius_visco_host**2)
    q_visco_out_host = thermal_cond_ice * delta_temp_elast_host / blt_crust_host
    q_elast_in_host = q_visco_out_host
    q_elast_out_host = thermal_cond_ice * (ice_top_temp_host - host.surface_temperature) / thickness_elast_host
    #     - Targ
    delta_temp_elast_targ = melting_temp_ice - ice_top_temp_targ
    q_visco_in_targ = (core_cooling_targ + tidal_heating_crust_targ) / (4. * np.pi * radius_visco_targ**2)
    q_visco_out_targ = thermal_cond_ice * delta_temp_elast_targ / blt_crust_targ
    q_elast_in_targ = q_visco_out_targ
    q_elast_out_targ = thermal_cond_ice * (ice_top_temp_targ - target.surface_temperature) / thickness_elast_targ

    # Update Change in Thickness
    #     - Host
    delta_temp_elast_host = melting_temp_ice - ice_top_temp_host
    change_visco_ice_host = (q_visco_in_host - q_visco_out_host) / (density_ice * latent_heat_ice)
    change_elast_ice_host = (q_elast_in_host - q_elast_out_host) / \
                            (density_ice * specific_heat_ice * delta_temp_elast_host)
    #     - Target
    delta_temp_elast_targ = melting_temp_ice - ice_top_temp_targ
    change_visco_ice_targ = (q_visco_in_targ - q_visco_out_targ) / (density_ice * latent_heat_ice)
    change_elast_ice_targ = (q_elast_in_targ - q_elast_out_targ) / \
                            (density_ice * specific_heat_ice * delta_temp_elast_targ)

    # Update Change in Spin-Rate
    #     - Host
    change_spin_freq_host = spin_rate_derivative(ztorque_host, host.moi)
    #     - Target
    change_spin_freq_targ = spin_rate_derivative(ztorque_targ, target.moi)

    # Update Change in Orbit
    change_semi_major_axis = semi_major_axis_derivative(semi_major_axis, host.mass, target.mass,
                                                        spin_freq_host, ztorque_host, tidal_heating_host,
                                                        spin_freq_targ, ztorque_targ, tidal_heating_targ)
    change_eccentricity = eccentricity_derivative(semi_major_axis, eccentricity, host.mass, target.mass,
                                                  spin_freq_host, ztorque_host, tidal_heating_host,
                                                  spin_freq_targ, ztorque_targ, tidal_heating_targ)

    # Return

    if diff_loop:
        return (
            change_semi_major_axis,
            change_eccentricity,
            change_spin_freq_host,
            change_core_temp_host,
            change_visco_ice_host,
            change_elast_ice_host,
            change_spin_freq_targ,
            change_core_temp_targ,
            change_visco_ice_targ,
            change_elast_ice_targ
        )
    else:
        return {
            'host'  : {
                'radiogenics'   : radio_heating_core_host,
                'tidal_heating' : tidal_heating_host,
                'crust_cooling' : crust_cooling_host,
                'core_cooling'  : core_cooling_host,
                'core_rayleigh' : rayleigh_core_host,
                'crust_rayleigh': rayleigh_crust_host,
                'ocean_volume'  : volume_ocean_host
            },
            'target': {
                'radiogenics'   : radio_heating_core_targ,
                'tidal_heating' : tidal_heating_targ,
                'crust_cooling' : crust_cooling_targ,
                'core_cooling'  : core_cooling_targ,
                'core_rayleigh' : rayleigh_core_targ,
                'crust_rayleigh': rayleigh_crust_targ,
                'ocean_volume'  : volume_ocean_targ
            }
        }

@njit
def diffeq_test(time, variables, complex_compliance_func, complex_compliance_input):

    semi_major_axis = variables[0, :]
    eccentricity = variables[1, :]
    spin_freq_host = variables[2, :]
    thickness_visco_host = variables[3, :]
    thickness_elast_host = variables[4, :]
    spin_freq_targ = variables[5, :]
    thickness_visco_targ = variables[6, :]
    thickness_elast_targ = variables[7, :]

    # Helper Values
    thickness_ice_host = thickness_visco_host + thickness_elast_host
    thickness_ice_targ = thickness_visco_targ + thickness_elast_targ
    overshoot_index_host = thickness_ice_host > 1. * thickness_crust_host
    overshoot_index_targ = thickness_ice_targ > 1. * thickness_crust_target

    # Check for value logic
    low_eccen = eccentricity < MIN_ECCEN
    if lock_zero_eccen:
        eccentricity[low_eccen] = 0.
    low_semi_a = semi_major_axis < MIN_SEMI_A
    semi_major_axis[low_semi_a] = MIN_SEMI_A

    thickness_visco_host[thickness_visco_host < 10.] = 10.
    thickness_visco_targ[thickness_visco_targ < 10.] = 10.
    thickness_elast_host[thickness_elast_host < 10.] = 10.
    thickness_elast_targ[thickness_elast_targ < 10.] = 10.
    # For now, any time the ice becomes too thick break it up into 75% visco, 24% elastic, and 1% ocean
    thickness_elast_host[overshoot_index_host] = 0.24 * thickness_crust_host
    thickness_visco_host[overshoot_index_host] = 0.75 * thickness_crust_host
    thickness_elast_targ[overshoot_index_targ] = 0.24 * thickness_crust_target
    thickness_visco_targ[overshoot_index_targ] = 0.75 * thickness_crust_target

    ocean_freezeout_index_host = (thickness_visco_host + thickness_elast_host) >= 0.98 * thickness_crust_host
    ocean_freezeout_index_targ = (thickness_visco_targ + thickness_elast_targ) >= 0.98 * thickness_crust_target
    ice_freezeout_index_host = thickness_elast_host >= 0.99 * thickness_crust_host
    ice_freezeout_index_targ = thickness_elast_targ >= 0.99 * thickness_crust_target

    # Conversions
    time_myr = sec2myr(time)
    orbital_freq = semi_a2orbital_motion(semi_major_axis, mass_host, mass_target)

    # Lock into 1:1 when either body is close to it. Ensures integrator is stable.
    lock_in_index_host = np.abs(spin_freq_host / orbital_freq) - 1. <= 0.1
    lock_in_index_targ = np.abs(spin_freq_targ / orbital_freq) - 1. <= 0.1
    if lock_to_1to1:
        spin_freq_host[lock_in_index_host] = orbital_freq[lock_in_index_host]
        spin_freq_targ[lock_in_index_targ] = orbital_freq[lock_in_index_targ]

    # Find Tidal Modes
    #     - Host
    modes_host, freqs_host, heating_coeffs_host, torque_coeffs_host = nsr_modes(orbital_freq, spin_freq_host,
                                                                                eccentricity, static_inclination)
    #     - Target
    modes_targ, freqs_targ, heating_coeffs_targ, torque_coeffs_targ = nsr_modes(orbital_freq, spin_freq_targ,
                                                                                eccentricity, static_inclination)

    # Calc Tidal Susceptibility
    #     - Host
    tidal_suscept_host = tidal_susceptibility_inflated_host / semi_major_axis**6
    #     - Target
    tidal_suscept_targ = tidal_susceptibility_inflated_target / semi_major_axis**6

    # Update Ice Shell Thickness
    #     - Host
    radius_elast_host = radius_host
    radius_visco_host = radius_elast_host - thickness_elast_host
    radius_ocean_host = radius_visco_host - thickness_visco_host
    volume_visco_host = (4. / 3.) * np.pi * (radius_visco_host**3 - radius_ocean_host**3)
    volume_ocean_host = (4. / 3.) * np.pi * (radius_ocean_host**3 - radius_core_host**3)
    #     - Target
    radius_elast_targ = radius_target
    radius_visco_targ = radius_elast_targ - thickness_elast_targ
    radius_ocean_targ = radius_visco_targ - thickness_visco_targ
    volume_visco_targ = (4. / 3.) * np.pi * (radius_visco_targ**3 - radius_ocean_targ**3)
    volume_ocean_targ = (4. / 3.) * np.pi * (radius_ocean_targ**3 - radius_core_target**3)

    # Where do tides concentrate?
    tidalvolfrac_crust_host = volume_visco_host / volume_host
    tidalvolfrac_crust_targ = volume_visco_targ / volume_target

    # Calculate Ice Cooling
    #     - Host
    delta_temp_crust_host = melting_temp_ice - ice_top_temp_host
    q_crust_conv_host, blt_crust_host, rayleigh_crust_host, neusselt_crust_host = convection(
            delta_temp_crust_host, thickness_visco_host, thermal_cond_ice, viscosity_crust_host,
            thermal_diff_ice, thermal_expansion_ice, gravity_surf_host, density_ice,
            convection_alpha_ice, convection_beta_ice, critical_rayleigh_ice)
    crust_cooling_host = q_crust_conv_host * 4. * np.pi * radius_visco_host**2
    #     - Target
    delta_temp_crust_targ = melting_temp_ice - ice_top_temp_targ
    q_crust_conv_targ, blt_crust_targ, rayleigh_crust_targ, neusselt_crust_targ = convection(
            delta_temp_crust_targ, thickness_visco_targ, thermal_cond_ice, viscosity_crust_targ,
            thermal_diff_ice, thermal_expansion_ice, gravity_surf_target, density_ice,
            convection_alpha_ice, convection_beta_ice, critical_rayleigh_ice)
    crust_cooling_targ = q_crust_conv_targ * 4. * np.pi * radius_visco_targ**2

    # Update Radiogenics
    #     - Host
    radio_heating_core_host = radiogenic_isotope(time_myr, mass_core_host, *standard_isotope_input)
    #     - Target
    radio_heating_core_targ = radiogenic_isotope(time_myr, mass_core_target, *standard_isotope_input)

    # Calculate Compliance
    #     - Host
    cmplx_comp_crust_bymode_host = list()
    for tidal_freq_h in freqs_host:
        cmplx_comp_crust_bymode_host.append(complex_compliance_func(shear_modulus_crust_host**-1, viscosity_crust_host,
                                                                    tidal_freq_h, *complex_compliance_input))
    #     - Target
    cmplx_comp_crust_bymode_targ = list()
    for tidal_freq_t in freqs_targ:
        cmplx_comp_crust_bymode_targ.append(complex_compliance_func(shear_modulus_crust_targ**-1, viscosity_crust_targ,
                                                                    tidal_freq_t, *complex_compliance_input))

    # Calculate Tides
    #     - Host
    effective_rigid_crust_host, tidal_heating_crust_host, ztorque_crust_host, _, _ = \
        calc_tides(gravity_surf_host, radius_host, density_ice, shear_modulus_crust_host,
                   tidal_suscept_host,
                   cmplx_comp_crust_bymode_host, heating_coeffs_host, torque_coeffs_host)
    #     - Target
    effective_rigid_crust_targ, tidal_heating_crust_targ, ztorque_crust_targ, _, _ = \
        calc_tides(gravity_surf_target, radius_target, density_ice, shear_modulus_crust_targ,
                   tidal_suscept_targ,
                   cmplx_comp_crust_bymode_targ, heating_coeffs_targ, torque_coeffs_targ)

    # Scale down by the volume factor
    tidal_heating_crust_host *= tidalvolfrac_crust_host
    tidal_heating_crust_targ *= tidalvolfrac_crust_targ
    ztorque_crust_host *= tidalvolfrac_crust_host
    ztorque_crust_targ *= tidalvolfrac_crust_targ

    # TODO: In this version, tides are restricted to the ice layer only
    tidal_heating_core_host = 0.
    tidal_heating_core_targ = 0.
    ztorque_core_host = 0.
    ztorque_core_targ = 0.

    tidal_heating_host = tidal_heating_crust_host + tidal_heating_core_host
    tidal_heating_targ = tidal_heating_crust_targ + tidal_heating_core_targ
    ztorque_host = ztorque_crust_host + ztorque_core_host
    ztorque_targ = ztorque_crust_targ + ztorque_core_targ

    # Find Heat Fluxes
    #     - Host
    q_visco_in_host = (radio_heating_core_host + tidal_heating_crust_host) / (4. * np.pi * radius_visco_host**2)
    q_visco_out_host = q_crust_conv_host
    q_elast_in_host = q_visco_out_host
    q_elast_out_host = thermal_cond_ice * (ice_top_temp_host - surface_temp_host) / thickness_elast_host
    #     - Target
    q_visco_in_targ = (radio_heating_core_targ + tidal_heating_crust_targ) / (4. * np.pi * radius_visco_targ**2)
    q_visco_out_targ = q_crust_conv_targ
    q_elast_in_targ = q_visco_out_targ
    q_elast_out_targ = thermal_cond_ice * (ice_top_temp_targ - surface_temp_targ) / thickness_elast_targ

    # Update Change in Thickness
    #     - Host
    delta_temp_elast_host = melting_temp_ice - ice_top_temp_host
    change_visco_ice_host = (q_visco_out_host - q_visco_in_host) / (density_ice * latent_heat_ice)
    change_elast_ice_host = (q_elast_out_host - q_elast_in_host) / \
                            (density_ice * specific_heat_ice * delta_temp_elast_host)
    #     - Target
    delta_temp_elast_targ = melting_temp_ice - ice_top_temp_targ
    change_visco_ice_targ = (q_visco_out_targ - q_visco_in_targ) / (density_ice * latent_heat_ice)
    change_elast_ice_targ = (q_elast_out_targ - q_elast_in_targ) / \
                            (density_ice * specific_heat_ice * delta_temp_elast_targ)

    if lock_ocean_freeze:
        # Once ocean is frozen (or mostly frozen) then the viscoelastic layer only grows as the inverse of the elastic.
        change_elast_ice_host[ocean_freezeout_index_host] = (q_elast_out_host[ocean_freezeout_index_host]) / \
                                (density_ice * specific_heat_ice * delta_temp_elast_host)
        change_elast_ice_targ[ocean_freezeout_index_targ] = (q_elast_out_targ[ocean_freezeout_index_targ]) / \
                                (density_ice * specific_heat_ice * delta_temp_elast_targ)

        change_visco_ice_host[ocean_freezeout_index_host] = -change_elast_ice_host[ocean_freezeout_index_host]
        change_visco_ice_targ[ocean_freezeout_index_targ] = -change_elast_ice_targ[ocean_freezeout_index_targ]

        # Once the elastic ice shell is takes over then all growth stops.
        change_visco_ice_host[ice_freezeout_index_host] = 0.
        change_elast_ice_host[ice_freezeout_index_host] = 0.
        change_visco_ice_targ[ice_freezeout_index_targ] = 0.
        change_elast_ice_targ[ice_freezeout_index_targ] = 0.

    surface_heat_flux_host = q_elast_out_host
    surface_heat_flux_target = q_elast_out_targ

    # Update Change in Spin-Rate
    #     - Host
    change_spin_freq_host = spin_rate_derivative(ztorque_host, moi_host)
    #     - Target
    change_spin_freq_targ = spin_rate_derivative(ztorque_targ, moi_target)

    if lock_to_1to1:
        change_spin_freq_host[lock_in_index_host] = 0.
        change_spin_freq_targ[lock_in_index_targ] = 0.

    # Update Change in Orbit
    change_semi_major_axis = semi_major_axis_derivative(semi_major_axis, mass_host, mass_target,
                                                        spin_freq_host, ztorque_host, tidal_heating_host,
                                                        spin_freq_targ, ztorque_targ, tidal_heating_targ)
    change_eccentricity = eccentricity_derivative(semi_major_axis, eccentricity, mass_host, mass_target,
                                                  spin_freq_host, ztorque_host, tidal_heating_host,
                                                  spin_freq_targ, ztorque_targ, tidal_heating_targ)

    # de_dt is prop to 1/e so the derivative can go crazy if e is too low.
    if lock_zero_eccen:
        # de_dt is prop to 1/e so the derivative can go crazy if e is too low.
        change_eccentricity[low_eccen] = 0.

    change_semi_major_axis[low_semi_a] = 0.

    return change_semi_major_axis, change_eccentricity, change_spin_freq_host, change_visco_ice_host,\
           change_elast_ice_host, change_spin_freq_targ, change_visco_ice_targ, change_elast_ice_targ, \
           radio_heating_core_host, tidal_heating_host, ztorque_host, surface_heat_flux_host, rayleigh_crust_host, volume_ocean_host, blt_crust_host, tidalvolfrac_crust_host, \
           radio_heating_core_targ, tidal_heating_targ, ztorque_targ, surface_heat_flux_target, rayleigh_crust_targ, volume_ocean_targ, blt_crust_targ, tidalvolfrac_crust_targ

@progress_bar(time_span, verbose=True)
def diffeq_nocore(time, variables, complex_compliance_func, complex_compliance_input, diff_loop):
    # semi_major_axis = variables[0, :]
    # eccentricity = variables[1, :]
    # spin_freq_host = variables[2, :]
    # thickness_visco_host = variables[3, :]
    # thickness_elast_host = variables[4, :]
    # spin_freq_targ = variables[5, :]
    # thickness_visco_targ = variables[6, :]
    # thickness_elast_targ = variables[7, :]
    # array_shape = semi_major_axis.shape
    #
    # # Helper Values
    # thickness_ice_host = thickness_visco_host + thickness_elast_host
    # thickness_ice_targ = thickness_visco_targ + thickness_elast_targ
    # overshoot_index_host = np.copy(thickness_ice_host > 0.99 * thickness_crust_host)
    # overshoot_index_targ = np.copy(thickness_ice_targ > 0.99 * thickness_crust_target)
    #
    # if diff_loop:
    #     # Check for value logic
    #     eccentricity[eccentricity < 1.e-10] = 0.
    #     semi_major_axis[semi_major_axis < MIN_SEMI_A] = MIN_SEMI_A
    #     # For now, any time the ice becomes too thick break it up into 75% visco, 24% elastic, and 1% ocean
    #     thickness_elast_host[overshoot_index_host] = 0.24 * thickness_crust_host
    #     thickness_visco_host[overshoot_index_host] = 0.75 * thickness_crust_host
    #     thickness_elast_targ[overshoot_index_targ] = 0.24 * thickness_crust_target
    #     thickness_visco_targ[overshoot_index_targ] = 0.75 * thickness_crust_target
    #
    # # Conversions
    # time_myr = sec2myr(time)
    # orbital_freq = semi_a2orbital_motion(semi_major_axis, mass_host, mass_target)
    #
    # if diff_loop:
    #     spin_freq_targ[np.abs(spin_freq_targ/orbital_freq) < .1] = orbital_freq[np.abs(spin_freq_host/orbital_freq) < .1]
    #     spin_freq_host[np.abs(spin_freq_host/orbital_freq) < .1] = orbital_freq[np.abs(spin_freq_host/orbital_freq) < .1]
    #
    # variables[0, :] = semi_major_axis
    # variables[1, :] = eccentricity
    # variables[2, :] = spin_freq_host
    # variables[3, :] = thickness_visco_host
    # variables[4, :] = thickness_elast_host
    # variables[5, :] = spin_freq_targ
    # variables[6, :] = thickness_visco_targ
    # variables[7, :] = thickness_elast_targ

    # # Find Tidal Modes
    # #     - Host
    # modes_host, freqs_host, heating_coeffs_host, torque_coeffs_host = nsr_modes(orbital_freq, spin_freq_host,
    #                                                                             eccentricity, static_inclination)
    # #     - Target
    # modes_targ, freqs_targ, heating_coeffs_targ, torque_coeffs_targ = nsr_modes(orbital_freq, spin_freq_targ,
    #                                                                             eccentricity, static_inclination)
    #
    # # Calc Tidal Susceptibility
    # #     - Host
    # tidal_suscept_host = tidal_susceptibility_inflated_host / semi_major_axis**6
    # #     - Target
    # tidal_suscept_targ = tidal_susceptibility_inflated_target / semi_major_axis**6
    #
    # # Update Ice Shell Thickness
    # #     - Host
    # radius_elast_host = radius_host
    # radius_visco_host = radius_elast_host - thickness_elast_host
    # radius_ocean_host = radius_visco_host - thickness_visco_host
    # volume_visco_host = (4. / 3.) * np.pi * (radius_visco_host**3 - radius_ocean_host**3)
    # volume_ocean_host = (4. / 3.) * np.pi * (radius_ocean_host**3 - radius_core_host**3)
    # #     - Target
    # radius_elast_targ = radius_target
    # radius_visco_targ = radius_elast_targ - thickness_elast_targ
    # radius_ocean_targ = radius_visco_targ - thickness_visco_targ
    # volume_visco_targ = (4. / 3.) * np.pi * (radius_visco_targ**3 - radius_ocean_targ**3)
    # volume_ocean_targ = (4. / 3.) * np.pi * (radius_ocean_targ**3 - radius_core_target**3)
    #
    # # Where do tides concentrate?
    # tidalvolfrac_crust_host = volume_visco_host / volume_host
    # tidalvolfrac_crust_targ = volume_visco_targ / volume_target
    #
    # # Calculate Ice Cooling
    # #     - Host
    # delta_temp_crust_host = melting_temp_ice - ice_top_temp_host
    # q_crust_conv_host, blt_crust_host, rayleigh_crust_host, neusselt_crust_host = convection(
    #         delta_temp_crust_host, thickness_visco_host, thermal_cond_ice, viscosity_crust_host,
    #         thermal_diff_ice, thermal_expansion_ice, gravity_surf_host, density_ice,
    #         convection_alpha_ice, convection_beta_ice, critical_rayleigh_ice)
    # crust_cooling_host = q_crust_conv_host * 4. * np.pi * radius_visco_host**2
    # #     - Target
    # delta_temp_crust_targ = melting_temp_ice - ice_top_temp_targ
    # q_crust_conv_targ, blt_crust_targ, rayleigh_crust_targ, neusselt_crust_targ = convection(
    #         delta_temp_crust_targ, thickness_visco_targ, thermal_cond_ice, viscosity_crust_targ,
    #         thermal_diff_ice, thermal_expansion_ice, gravity_surf_target, density_ice,
    #         convection_alpha_ice, convection_beta_ice, critical_rayleigh_ice)
    # crust_cooling_targ = q_crust_conv_targ * 4. * np.pi * radius_visco_targ**2
    #
    # # Update Radiogenics
    # #     - Host
    # radio_heating_core_host = radiogenic_isotope(time_myr, mass_core_host, *standard_isotope_input)
    # #     - Target
    # radio_heating_core_targ = radiogenic_isotope(time_myr, mass_core_target, *standard_isotope_input)
    #
    # # Calculate Compliance and Tides
    # #     - Host
    # cmplx_comp_crust_bymode_host = [
    #     complex_compliance_func(shear_modulus_crust_host**-1, viscosity_crust_host, tidal_freq,
    #                             *complex_compliance_input)
    #     for tidal_freq in freqs_host
    # ]
    #
    # effective_rigid_crust_host, tidal_heating_crust_host, ztorque_crust_host, _, _ = \
    #     calc_tides(gravity_surf_host, radius_host, density_ice, shear_modulus_crust_host,
    #                tidal_suscept_host,
    #                cmplx_comp_crust_bymode_host, heating_coeffs_host, torque_coeffs_host)
    #
    # # Scale down by the volume factor
    # tidal_heating_crust_host *= tidalvolfrac_crust_host
    # ztorque_crust_host *= tidalvolfrac_crust_host
    #
    # # TODO: In this version, tides are restricted to the ice layer only
    # tidal_heating_core_host = 0.
    # ztorque_core_host = 0.
    #
    # tidal_heating_host = tidal_heating_crust_host + tidal_heating_core_host
    # ztorque_host = ztorque_crust_host + ztorque_core_host
    #
    # #     - Target
    # cmplx_comp_crust_bymode_targ = [
    #     complex_compliance_func(shear_modulus_crust_targ**-1, viscosity_crust_targ, tidal_freq,
    #                             *complex_compliance_input)
    #     for tidal_freq in freqs_targ
    # ]
    #
    # effective_rigid_crust_targ, tidal_heating_crust_targ, ztorque_crust_targ, _, _ = \
    #     calc_tides(gravity_surf_target, radius_target, density_ice, shear_modulus_crust_targ,
    #                tidal_suscept_targ,
    #                cmplx_comp_crust_bymode_targ, heating_coeffs_targ, torque_coeffs_targ)
    #
    # # Scale down by the volume factor
    # tidal_heating_crust_targ *= tidalvolfrac_crust_targ
    # ztorque_crust_targ *= tidalvolfrac_crust_targ
    #
    # # TODO: In this version, tides are restricted to the ice layer only
    # tidal_heating_core_targ = 0.
    # ztorque_core_targ = 0.
    #
    # tidal_heating_targ = tidal_heating_crust_targ + tidal_heating_core_targ
    # ztorque_targ = ztorque_crust_targ + ztorque_core_targ
    #
    # # Find Heat Fluxes
    # #     - Host
    # delta_temp_elast_host = melting_temp_ice - ice_top_temp_host
    # q_visco_in_host = (radio_heating_core_host + tidal_heating_crust_host) / (4. * np.pi * radius_visco_host**2)
    # q_visco_out_host = thermal_cond_ice * delta_temp_elast_host / blt_crust_host
    # q_elast_in_host = q_visco_out_host
    # q_elast_out_host = thermal_cond_ice * (ice_top_temp_host - host.surface_temperature) / thickness_elast_host
    # #     - Target
    # delta_temp_elast_targ = melting_temp_ice - ice_top_temp_targ
    # q_visco_in_targ = (radio_heating_core_targ + tidal_heating_crust_targ) / (4. * np.pi * radius_visco_targ**2)
    # q_visco_out_targ = thermal_cond_ice * delta_temp_elast_targ / blt_crust_targ
    # q_elast_in_targ = q_visco_out_targ
    # q_elast_out_targ = thermal_cond_ice * (ice_top_temp_targ - target.surface_temperature) / thickness_elast_targ
    #
    # # Update Change in Thickness
    # #     - Host
    # delta_temp_elast_host = melting_temp_ice - ice_top_temp_host
    # change_visco_ice_host = (q_visco_in_host - q_visco_out_host) / (density_ice * latent_heat_ice)
    # change_elast_ice_host = (q_elast_in_host - q_elast_out_host) / \
    #                         (density_ice * specific_heat_ice * delta_temp_elast_host)
    # #     - Target
    # delta_temp_elast_targ = melting_temp_ice - ice_top_temp_targ
    # change_visco_ice_targ = (q_visco_in_targ - q_visco_out_targ) / (density_ice * latent_heat_ice)
    # change_elast_ice_targ = (q_elast_in_targ - q_elast_out_targ) / \
    #                         (density_ice * specific_heat_ice * delta_temp_elast_targ)
    #
    # # Update Change in Spin-Rate
    # #     - Host
    # change_spin_freq_host = spin_rate_derivative(ztorque_host, moi_host)
    # #     - Target
    # change_spin_freq_targ = spin_rate_derivative(ztorque_targ, moi_target)
    #
    # # Update Change in Orbit
    # change_semi_major_axis = semi_major_axis_derivative(semi_major_axis, mass_host, mass_target,
    #                                                     spin_freq_host, ztorque_host, tidal_heating_host,
    #                                                     spin_freq_targ, ztorque_targ, tidal_heating_targ)
    # change_eccentricity = eccentricity_derivative(semi_major_axis, eccentricity, mass_host, mass_target,
    #                                               spin_freq_host, ztorque_host, tidal_heating_host,
    #                                               spin_freq_targ, ztorque_targ, tidal_heating_targ)

    # Return

    change_semi_major_axis, change_eccentricity, change_spin_freq_host, change_visco_ice_host, \
    change_elast_ice_host, change_spin_freq_targ, change_visco_ice_targ, change_elast_ice_targ, \
    radio_heating_core_host, tidal_heating_host, ztorque_host, surface_heat_flux_host, rayleigh_crust_host, volume_ocean_host, blt_crust_host, tidalvolfrac_crust_host,\
    radio_heating_core_targ, tidal_heating_targ, ztorque_targ, surface_heat_flux_target, rayleigh_crust_targ, volume_ocean_targ, blt_crust_targ, tidalvolfrac_crust_targ = \
        diffeq_test(time, variables, complex_compliance_func, complex_compliance_input)

    if diff_loop:
        return (
            change_semi_major_axis,
            change_eccentricity,
            change_spin_freq_host,
            change_visco_ice_host,
            change_elast_ice_host,
            change_spin_freq_targ,
            change_visco_ice_targ,
            change_elast_ice_targ
        )
    else:
        return {
            'derivatives': {
                'semi_major_axis': change_semi_major_axis,
                'eccentricity'   : change_eccentricity,
                'spin_freq_host' : change_spin_freq_host,
                'visco_ice_host' : change_visco_ice_host,
                'elast_ice_host' : change_elast_ice_host,
                'spin_freq_targ' : change_spin_freq_targ,
                'visco_ice_targ' : change_visco_ice_targ,
                'elast_ice_targ' : change_elast_ice_targ,
            },
            'host'  : {
                'crust_blt'     : blt_crust_host,
                'radiogenics'   : radio_heating_core_host,
                'tidal_heating' : tidal_heating_host,
                'tidal_torque'  : ztorque_host,
                'crust_cooling_flux' : surface_heat_flux_host,
                'crust_rayleigh': rayleigh_crust_host,
                'ocean_volume'  : volume_ocean_host,
                'tvf': tidalvolfrac_crust_host
            },
            'target': {
                'crust_blt'     : blt_crust_targ,
                'radiogenics'   : radio_heating_core_targ,
                'tidal_heating' : tidal_heating_targ,
                'tidal_torque'  : ztorque_targ,
                'crust_cooling_flux' : surface_heat_flux_target,
                'crust_rayleigh': rayleigh_crust_targ,
                'ocean_volume'  : volume_ocean_targ,
                'tvf': tidalvolfrac_crust_targ
            }
        }


def integrate(rheology_name, compliance_func, compliance_input, time_span, initial_conds):

    if use_nocore:
        diffeq = diffeq_nocore
    else:
        diffeq = diffeq_withcore

    log('Starting integration on rheology: {}'.format(rheology_name))
    diff_start_clock = time.time()
    diffeq_partial = partial(diffeq, complex_compliance_func=compliance_func, complex_compliance_input=compliance_input,
                             diff_loop=True)

    solution = solve_ivp(diffeq_partial, time_span, initial_conds,
                         method=method, vectorized=True, rtol=rtol)

    if not solution.success:
        log('\nWARNING - Integration may not have been preformed correctly.')
        log('ODEINT MSG:\n' + solution.message)
    else:
        log('\nNo obvious issues with integration found!')

    # Grab other data
    log('Grabbing auxiliary data')

    if len(solution.t) > 50000:
        log('Solution data size is very large. Reducing to avoid memory errors in auxiliary grab.')
        solution_time_domain = np.linspace(solution.t[0], solution.t[-1], 50000)
        solution_y = np.asarray([np.interp(solution_time_domain, solution.t, solution.y[i, :])
                                 for i in range(solution.y.shape[0])])
    else:
        solution_time_domain = solution.t
        solution_y = solution.y
    del solution

    # Call the differential function in non-diff mode to gather auxiliary data.
    all_data = diffeq(solution_time_domain, solution_y,
                      complex_compliance_func=compliance_func, complex_compliance_input=compliance_input,
                      diff_loop=False)
    # Put all data, including independent variables, into one dictionary.
    all_data['time_domain'] = solution_time_domain
    all_data['time_domain_myr'] = sec2myr(solution_time_domain)
    if use_nocore:
        all_data['semi_major_axis'] = solution_y[0, :]
        all_data['eccentricity'] = solution_y[1, :]
        all_data['host']['spin_freq'] = solution_y[2, :]
        all_data['host']['thickness_visco_ice'] = solution_y[3, :]
        all_data['host']['thickness_elast_ice'] = solution_y[4, :]
        all_data['target']['spin_freq'] = solution_y[5, :]
        all_data['target']['thickness_visco_ice'] = solution_y[6, :]
        all_data['target']['thickness_elast_ice'] = solution_y[7, :]
    else:
        all_data['semi_major_axis'] = solution_y[0, :]
        all_data['eccentricity'] = solution_y[1, :]
        all_data['host']['spin_freq'] = solution_y[2, :]
        all_data['host']['temperature_core'] = solution_y[3, :]
        all_data['host']['thickness_visco_ice'] = solution_y[4, :]
        all_data['host']['thickness_elast_ice'] = solution_y[5, :]
        all_data['target']['spin_freq'] = solution_y[6, :]
        all_data['target']['temperature_core'] = solution_y[7, :]
        all_data['target']['thickness_visco_ice'] = solution_y[8, :]
        all_data['target']['thickness_elast_ice'] = solution_y[9, :]


    log('Integration took {:0>2.0f} days {:0>2.0f}:{:0>2.0f}::{:0>4.1f}'.format(*convert_to_hms(time.time() -
                                                                                                diff_start_clock)))
    log('Finished integration on {}'.format(rheology_name))

    return all_data


def main():
    all_data_by_rheo = dict()
    rheo_inputs = {
        'Andrade': (andrade_alpha_ice, andrade_zeta_ice),
        'Maxwell': tuple(),
        'off': tuple()
    }
    rheo_funcs = {
        'Andrade': andrade,
        'Maxwell': maxwell,
        'off': off
    }
    rheo_colors = {
        'Andrade': 'b',
        'Maxwell': 'g',
        'off': 'k'
    }

    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    # Figure 1:
    fig_1, axes_1 = plt.subplots(2, 2, figsize=(12, 8))
    fig_1.subplots_adjust(wspace=4)
    for axis in axes_1.flatten():
        axis.set_xlabel('Time [Myr]')
    axes_1_axis_1duel = axes_1[0, 0].twinx()

    axes_1_axis_1duel.set_ylabel('Eccentricity')
    axes_1_axis_1duel.set_yscale('log')
    axes_1[0, 0].set_ylabel('Semi-Major Axis [modern frac]')
    axes_1[0, 1].set_ylabel('Spin / n')
    axes_1[1, 0].set_ylabel('Radiogenics [Watts]')
    axes_1[1, 0].set_yscale('log')
    axes_1[1, 1].set_ylabel('Tidal Heating [Watts]')
    axes_1[1, 1].set_yscale('log')

    # Figure 2:
    fig_2, axes_2 = plt.subplots(2, 2, figsize=(12, 8))
    fig_2.subplots_adjust(wspace=4)
    for axis in axes_2.flatten():
        axis.set_xlabel('Time [Myr]')
    # axes_2_axis_0duel = axes_2[0, 0].twinx()
    # axes_2_axis_1duel = axes_2[0, 1].twinx()

    axes_2[0, 0].set_ylabel('H2O Thickness [frac] (Host)')
    axes_2[0, 1].set_ylabel('H2O Thickness [frac] (Target)')
    # axes_2_axis_0duel.set_ylabel('Visco. Boundary Layer [perc] (dot-dash)')
    # axes_2_axis_1duel.set_ylabel('Visco. Boundary Layer [perc] (dot-dash)')
    axes_2[1, 0].set_ylabel('Surface Flux [W m-2]')
    axes_2[1, 0].set_yscale('log')
    axes_2[1, 1].set_ylabel('$\Delta{}D$ [km]')

    # Figure 3:
    fig_3, axes_3 = plt.subplots(2, 2, figsize=(12, 8))
    fig_3.subplots_adjust(wspace=4)
    for axis in axes_3.flatten():
        axis.set_xlabel('Time [Myr]')
    axes_3_axis_2duel = axes_3[1, 0].twinx()

    axes_3[0, 0].set_ylabel('Tidal Torque (Host) [N m]')
    axes_3[0, 0].set_yscale('log')
    axes_3[0, 1].set_ylabel('Tidal Torque (Target) [N m]')
    axes_3[0, 1].set_yscale('log')
    axes_3[1, 0].set_ylabel('Circularization Timescale [Myr]')
    axes_3_axis_2duel.set_ylabel('da/dt [km Myr-1]')
    axes_3[1, 0].set_yscale('log')
    axes_3_axis_2duel.set_yscale('log')
    axes_3[1, 1].set_ylabel('Spin Dampening Timescale [Myr]')
    axes_3[1, 1].set_yscale('log')

    # Calculate No Tides
    all_data_notide = integrate('off', off, tuple(), time_span, initial_conditions)
    notide_x = all_data_notide['time_domain_myr']
    host_notide_visco_dx = all_data_notide['host']['thickness_visco_ice']
    host_notide_elast_dx = all_data_notide['host']['thickness_elast_ice']
    host_notide_ocean_dx = host.crust.thickness - (host_notide_visco_dx + host_notide_elast_dx)
    targ_notide_visco_dx = all_data_notide['target']['thickness_visco_ice']
    targ_notide_elast_dx = all_data_notide['target']['thickness_elast_ice']
    targ_notide_ocean_dx = target.crust.thickness - (targ_notide_visco_dx + targ_notide_elast_dx)

    for r_i, (rheo_name, comp_func) in enumerate(rheo_funcs.items()):
        comp_input = rheo_inputs[rheo_name]
        rheo_color = rheo_colors[rheo_name]

        # Calculate
        all_data = integrate(rheo_name, comp_func, comp_input, time_span, initial_conditions)

        # Store
        all_data_by_rheo[rheo_name] = all_data

        # Plot Figure 1
        # Pane 1: Orbit
        x = all_data['time_domain_myr']
        axes_1[0, 0].plot(x, all_data['semi_major_axis'] / modern_semi_major_axis, c=rheo_color)
        axes_1_axis_1duel.plot(x, all_data['eccentricity'], c=rheo_color, ls=':')

        # Pane 2: Spin
        orb_freq = semi_a2orbital_motion(all_data['semi_major_axis'], host.mass, target.mass)
        host_spin = all_data['host']['spin_freq']
        target_spin = all_data['target']['spin_freq']
        axes_1[0, 1].plot(x, host_spin / orb_freq, c=rheo_color, ls='-')
        axes_1[0, 1].plot(x, target_spin / orb_freq, c=rheo_color, ls='--')

        # Pane 3: Radiogenics
        host_radio_heating = all_data['host']['radiogenics']
        target_radio_heating = all_data['target']['radiogenics']
        axes_1[1, 0].plot(x, host_radio_heating, c=rheo_color, ls='-')
        axes_1[1, 0].plot(x, target_radio_heating, c=rheo_color, ls='--')

        # Pane 4: Heating
        host_tidal_heating = all_data['host']['tidal_heating']
        target_tidal_heating = all_data['target']['tidal_heating']
        axes_1[1, 1].plot(x, host_tidal_heating, c=rheo_color, ls='-')
        axes_1[1, 1].plot(x, target_tidal_heating, c=rheo_color, ls='--')

        # Plot Figure 2:
        # Pane 1: Host Thickness
        host_visco_dx = all_data['host']['thickness_visco_ice']
        host_elast_dx = all_data['host']['thickness_elast_ice']
        host_blt_frac = all_data['host']['crust_blt'] / host_visco_dx
        host_ocean_dx = (host.crust.thickness - host_visco_dx - host_elast_dx)
        axes_2[0, 0].plot(x, host_elast_dx / host.crust.thickness, c=rheo_color, ls='-')
        axes_2[0, 0].plot(x, host_visco_dx / host.crust.thickness, c=rheo_color, ls='--')
        axes_2[0, 0].plot(x, host_ocean_dx / host.crust.thickness, c=rheo_color, ls=':')
        # axes_2_axis_0duel.plot(x, host_blt_frac, c=rheo_color, ls='-.', alpha = 0.3)

        # Pane 2: Target Thickness
        target_visco_dx = all_data['target']['thickness_visco_ice']
        target_elast_dx = all_data['target']['thickness_elast_ice']
        target_blt_frac = all_data['target']['crust_blt'] / target_visco_dx
        target_ocean_dx = (target.crust.thickness - target_visco_dx - target_elast_dx)
        axes_2[0, 1].plot(x, target_elast_dx / target.crust.thickness, c=rheo_color, ls='-')
        axes_2[0, 1].plot(x, target_visco_dx / target.crust.thickness, c=rheo_color, ls='--')
        axes_2[0, 1].plot(x, target_ocean_dx / target.crust.thickness, c=rheo_color, ls=':')
        # axes_2_axis_1duel.plot(x, target_blt_frac, c=rheo_color, ls='-.', alpha = 0.3)

        # Pane 3: Crust Cooling
        host_cooling = all_data['host']['crust_cooling_flux']
        target_cooling = all_data['target']['crust_cooling_flux']
        axes_2[1, 0].plot(x, host_cooling, c=rheo_color, ls='-')
        axes_2[1, 0].plot(x, target_cooling, c=rheo_color, ls='--')

        # Pane 4: Delta_D and Tidal Volume Fraction
        delta_d_host = host_ocean_dx - np.interp(x, notide_x, host_notide_ocean_dx)
        delta_d_targ = target_ocean_dx - np.interp(x, notide_x, targ_notide_ocean_dx)
        axes_2[1, 1].plot(x, delta_d_host / 1000, c=rheo_color, ls='-')
        axes_2[1, 1].plot(x, delta_d_targ / 1000, c=rheo_color, ls='--')

        # Plot Figure 3
        # Pane 1+2: Tidal Torque
        host_ztorque_pos, host_ztorque_neg = neg_array_for_log_plot(all_data['host']['tidal_torque'])
        targ_ztorque_pos, targ_ztorque_neg = neg_array_for_log_plot(all_data['target']['tidal_torque'])
        axes_3[0, 0].plot(x, host_ztorque_pos, c=rheo_color, ls='-')
        axes_3[0, 0].plot(x, host_ztorque_neg, c=rheo_color, ls='-', alpha=0.3)
        axes_3[0, 1].plot(x, targ_ztorque_pos, c=rheo_color, ls='-')
        axes_3[0, 1].plot(x, targ_ztorque_neg, c=rheo_color, ls='-', alpha=0.3)

        # Pane 3: Semi-A and Eccen Derivatives
        semi_a_deriv = all_data['derivatives']['semi_major_axis']
        eccen_deriv = all_data['derivatives']['eccentricity']
        semi_a_deriv *= 3.154e13 / 1000. # Convert from [m s-1] to [km myr-1]
        circ_timescale = sec2myr(-all_data['eccentricity'] / eccen_deriv)
        semi_a_deriv_pos, semi_a_deriv_neg = neg_array_for_log_plot(semi_a_deriv)
        circ_timescale_pos, circ_timescale_neg = neg_array_for_log_plot(circ_timescale)
        axes_3[1, 0].plot(x, semi_a_deriv_pos, c=rheo_color)
        axes_3[1, 0].plot(x, semi_a_deriv_neg, c=rheo_color, alpha=0.3)
        axes_3_axis_2duel.plot(x, circ_timescale_pos, c=rheo_color, ls=':')
        axes_3_axis_2duel.plot(x, circ_timescale_neg, c=rheo_color, ls=':', alpha=0.3)

        # Pane 4: Spin-Freq Timescale
        orbital_motion = semi_a2orbital_motion(all_data['semi_major_axis'], host.mass, target.mass)
        spin_freq_deriv_host = all_data['derivatives']['spin_freq_host']
        spin_freq_deriv_targ = all_data['derivatives']['spin_freq_targ']
        spin_timescale_host = -np.abs(orbital_motion - all_data['host']['spin_freq']) / spin_freq_deriv_host
        spin_timescale_targ = -np.abs(orbital_motion - all_data['target']['spin_freq']) / spin_freq_deriv_targ
        spin_timescale_host = sec2myr(spin_timescale_host)
        spin_timescale_targ = sec2myr(spin_timescale_targ)
        spin_timescale_host_pos, spin_timescale_host_neg = neg_array_for_log_plot(spin_timescale_host)
        spin_timescale_targ_pos, spin_timescale_targ_neg = neg_array_for_log_plot(spin_timescale_targ)
        axes_3[1, 1].plot(x, spin_timescale_host_pos, c=rheo_color, ls='-')
        axes_3[1, 1].plot(x, spin_timescale_targ_pos, c=rheo_color, ls='--')
        axes_3[1, 1].plot(x, spin_timescale_host_neg, c=rheo_color, ls='-', alpha=0.3)
        axes_3[1, 1].plot(x, spin_timescale_targ_neg, c=rheo_color, ls='--', alpha=0.3)

    # Legends
    solid  = Line2D([0], [0], c='k', ls='-', lw=2)
    dashed = Line2D([0], [0], c='k', ls='-.', lw=2)
    dotted = Line2D([0], [0], c='k', ls=':', lw=2)
    axes_1[0, 0].legend([solid, dotted], ['Semi-a', 'Eccen.'], loc='best')
    axes_1[0, 1].legend([solid, dashed], [f'{host.name}', f'{target.name}'], loc='best')
    axes_2[0, 0].legend([solid, dashed, dotted], ['Elastic', 'Visco', 'Ocean'], loc='best')
    axes_2[0, 1].legend([solid, dashed, dotted], ['Elastic', 'Visco', 'Ocean'], loc='best')
    axes_2[1, 1].legend([solid, dashed], [f'{host.name}', f'{target.name}'], loc='best')
    axes_3[1, 0].legend([solid, dotted], [f'da/dt', f'Circ.'], loc='best')
    axes_3[1, 1].legend([solid, dashed], [f'{host.name}', f'{target.name}'], loc='best')

    rheo_legend_names, rheo_legend_lines = list(), list()
    for rheo_name, rheo_color in rheo_colors.items():
        rheo_legend_names.append(rheo_name)
        rheo_legend_lines.append(Line2D([0], [0], c=rheo_color, ls='-', lw=2))

    for axis in [axes_1[1, 0], axes_2[1, 0], axes_3[1, 0]]:

        axis.legend(rheo_legend_lines, rheo_legend_names,
                    loc='upper center', bbox_to_anchor=(.5, -0.15), ncol=3)

    plt.tight_layout()
    plt.show()
    return all_data_by_rheo, fig_1, fig_2, fig_3


if __name__ == '__main__':
    main()
