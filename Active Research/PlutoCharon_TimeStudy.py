import time
from functools import partial

import numpy as np
from scipy.constants import R
from scipy.integrate import solve_ivp

import TidalPy
from TidalPy import log
from TidalPy.orbit.orbit import OrbitBase as Orbit
from TidalPy.constants import R_pluto
from TidalPy.utilities import progress_bar
from TidalPy.utilities.conversions import sec2myr, myr2sec, semi_a2orbital_motion, convert_to_hms
from TidalPy.radiogenics import standard_isotope_input, radiogenic_isotope
from TidalPy.thermal import convection
from TidalPy.thermal.viscosity_models import reference as calc_reference_viscosity
from TidalPy.rheology import complex_love, effective_rigidity
from TidalPy.rheology.compliance_models import maxwell, andrade
from TidalPy.orbit.modes import nsr_modes
from TidalPy.dynamics import spin_rate_derivative
from TidalPy.dynamics.duel_dissipation import semi_major_axis_derivative, eccentricity_derivative

# Planet Building
target = TidalPy.build_planet('charon', force_build=True)
host = TidalPy.build_planet('pluto', force_build=True)
star = TidalPy.build_planet('sol', force_build=True)
orbit = Orbit(star, host, target, duel_dissipation=True, time_study=True)

# Initial Conditions
time_span = myr2sec(np.asarray((0, 5)))
modern_semi_major_axis = 1.9599e7
modern_orb_freq = semi_a2orbital_motion(modern_semi_major_axis, host.mass, target.mass)

semi_major_axis_init = 1.0 * modern_semi_major_axis
eccentricity_init = 0.4
spin_freq_host_init = modern_orb_freq * 20
temperature_core_host_init = 1000.
thickness_visco_host_init = 0.2 * host.crust.thickness
thickness_elast_host_init = 0.1 * host.crust.thickness
spin_freq_targ_init = modern_orb_freq * 20
temperature_core_targ_init = 1000.
thickness_visco_targ_init = 0.2 * target.crust.thickness
thickness_elast_targ_init = 0.1 * target.crust.thickness

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

# Other Conditions
melting_temp_ice = 273.15
ocean_temp = melting_temp_ice

# Integration Parameters
rtol = 1.e-8
method = 'LSODA'

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
static_inclination = np.asarray([0.])

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


@progress_bar(time_span, verbose=True)
def diffeq(time, variables, complex_compliance_func, complex_compliance_input, diff_loop):

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
    crust_cooling_host = q_crust_conv_host * 4. * np.pi * radius_visco_targ**2
    #     - Target
    delta_temp_crust_targ = melting_temp_ice - ice_top_temp_targ
    q_crust_conv_targ, blt_crust_targ, rayleigh_crust_targ, neusselt_crust_targ = convection(
            delta_temp_crust_targ, thickness_visco_targ, thermal_cond_ice, viscosity_crust_targ,
            thermal_diff_ice, thermal_expansion_ice, target.gravity_surf, density_ice,
            convection_alpha_ice, convection_beta_ice, critical_rayleigh_ice)
    crust_cooling_targ = q_crust_conv_targ * 4. * np.pi * radius_visco_targ**2

    # Update Radiogenics
    #     - Host
    radio_heating_core_host = 0.*radiogenic_isotope(time_myr, host.core.mass, *standard_isotope_input)
    #     - Target
    radio_heating_core_targ = 0.*radiogenic_isotope(time_myr, target.core.mass, *standard_isotope_input)

    # Update Rheology and Tides
    #     - Host
    effective_rigid_crust_host = effective_rigidity(shear_modulus_crust_host,
                                                    host.gravity_surf, host.radius, host.density_bulk)
    cmplx_comp_crust_bymode_host = [
        complex_compliance_func(shear_modulus_crust_host**-1, viscosity_crust_host, tidal_freq, *complex_compliance_input)
        for tidal_freq in freqs_host
    ]

    for mode_i, cmplx_comp in enumerate(cmplx_comp_crust_bymode_host):
        new_comp = cmplx_comp.copy()
        new_comp[np.abs(np.imag(cmplx_comp)) == np.inf] = shear_modulus_crust_host**-1 + 0.j
        cmplx_comp_crust_bymode_host[mode_i] = new_comp

    cmplx_love_crust_bymode_host = [
        complex_love(cmplx_comp, shear_modulus_crust_host, effective_rigid_crust_host)
        for cmplx_comp in cmplx_comp_crust_bymode_host]
    tidal_hating_crust_bymode_host = [
        -np.imag(love_num) * tidal_suscept_host * heating_coeff
        for love_num, heating_coeff in zip(cmplx_love_crust_bymode_host, heating_coeffs_host)
    ]
    ztorque_crust_bymode_host = [
        -np.imag(love_num) * tidal_suscept_host * torque_coeff
        for love_num, torque_coeff in zip(cmplx_love_crust_bymode_host, torque_coeffs_host)
    ]
    tidal_heating_crust_host = sum(tidal_hating_crust_bymode_host)
    ztorque_crust_host = sum(ztorque_crust_bymode_host)
    # Scale down by the volume factor
    tidal_heating_crust_host *= tidalvolfrac_crust_host
    ztorque_crust_host *= tidalvolfrac_crust_host

    # TODO: In this version, tides are restricted to the ice layer only
    tidal_heating_core_host = 0.
    ztorque_core_host = 0.

    tidal_heating_host = tidal_heating_crust_host + tidal_heating_core_host
    ztorque_host = ztorque_crust_host + ztorque_core_host

    #     - Target
    effective_rigid_crust_targ = effective_rigidity(shear_modulus_crust_targ,
                                                    target.gravity_surf, target.radius, target.density_bulk)
    cmplx_comp_crust_bymode_targ = [
        complex_compliance_func(shear_modulus_crust_targ**-1, viscosity_crust_targ, tidal_freq, *complex_compliance_input)
        for tidal_freq in freqs_targ
    ]

    for mode_i, cmplx_comp in enumerate(cmplx_comp_crust_bymode_targ):
        new_comp = cmplx_comp.copy()
        new_comp[np.abs(np.imag(cmplx_comp)) == np.inf] = shear_modulus_crust_targ**-1 + 0.j
        cmplx_comp_crust_bymode_targ[mode_i] = new_comp

    cmplx_love_crust_bymode_targ = [
        complex_love(cmplx_comp, shear_modulus_crust_targ, effective_rigid_crust_targ)
        for cmplx_comp in cmplx_comp_crust_bymode_targ]
    tidal_hating_crust_bymode_targ = [
        -np.imag(love_num) * tidal_suscept_targ * heating_coeff
        for love_num, heating_coeff in zip(cmplx_love_crust_bymode_targ, heating_coeffs_targ)
    ]
    ztorque_crust_bymode_targ = [
        -np.imag(love_num) * tidal_suscept_targ * torque_coeff
        for love_num, torque_coeff in zip(cmplx_love_crust_bymode_targ, torque_coeffs_targ)
    ]
    tidal_heating_crust_targ = sum(tidal_hating_crust_bymode_targ)
    ztorque_crust_targ = sum(ztorque_crust_bymode_targ)
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
            'host': {
                'radiogenics': radio_heating_core_host,
                'tidal_heating': tidal_heating_host,
                'crust_cooling': crust_cooling_host,
                'core_cooling': core_cooling_host,
                'core_rayleigh': rayleigh_core_host,
                'crust_rayleigh': rayleigh_crust_host,
                'ocean_volume': volume_ocean_host
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

def integrate(rheology_name, compliance_func, compliance_input, time_span, initial_conds):

    log('Starting integration on {}'.format(rheology_name))
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
        'Maxwell': tuple()
    }
    rheo_funcs = {
        'Andrade': andrade,
        'Maxwell': maxwell
    }
    rheo_colors = {
        'Andrade': 'b',
        'Maxwell': 'g'
    }

    import matplotlib.pyplot as plt

    # Figure 1:
    fig_1, axes_1 = plt.subplots(2, 2, figsize=(12, 8))
    fig_1.subplots_adjust(wspace=2.8)
    axes_1_axis_1duel = axes_1[0, 0].twinx()
    axes_1_axis_1duel.set_xlabel('Time [Myr]')
    for axis in axes_1.flatten():
        axis.set_xlabel('Time [Myr]')

    axes_1_axis_1duel.set_ylabel('Eccentricity (dashed)')
    axes_1_axis_1duel.set_yscale('log')
    axes_1[0, 0].set_ylabel('Semi-Major Axis [modern frac]')
    axes_1[0, 1].set_ylabel('Spin / n')
    axes_1[1, 0].set_ylabel('Thickness [% of H2O]')
    axes_1[1, 1].set_ylabel('Heating [Watts]')
    axes_1[1, 1].set_yscale('log')

    for rheo_name, comp_func in rheo_funcs.items():
        comp_input = rheo_inputs[rheo_name]
        rheo_color = rheo_colors[rheo_name]

        # Calculate
        all_data = integrate(rheo_name, comp_func, comp_input, time_span, initial_conditions)

        # Store
        all_data_by_rheo[rheo_name] = all_data

        # Plot
        # Pane 1: Orbit
        x = all_data['time_domain_myr']
        axes_1[0, 0].plot(x, all_data['semi_major_axis'] / modern_semi_major_axis, c=rheo_color)
        axes_1_axis_1duel.plot(x, all_data['eccentricity'] / modern_semi_major_axis, c=rheo_color, ls='--')

        # Pane 2: Spin
        orb_freq = semi_a2orbital_motion(all_data['semi_major_axis'], host.mass, target.mass)
        host_spin = all_data['host']['spin_freq']
        target_spin = all_data['target']['spin_freq']
        axes_1[0, 1].plot(x, host_spin / orb_freq, c=rheo_color, ls='-')
        axes_1[0, 1].plot(x, target_spin / orb_freq, c=rheo_color, ls='--')

        # Pane 3: Ocean
        host_visco_dx = all_data['host']['thickness_visco_ice']
        target_visco_dx = all_data['target']['thickness_visco_ice']
        host_elast_dx = all_data['host']['thickness_elast_ice']
        target_elast_dx = all_data['target']['thickness_elast_ice']
        host_ocean_dx = (host.crust.thickness - host_visco_dx - host_elast_dx)
        target_ocean_dx = (target.crust.thickness - target_visco_dx - target_elast_dx)

        axes_1[1, 0].plot(x, host_visco_dx / host.crust.thickness, c=rheo_color, ls='-')
        axes_1[1, 0].plot(x, target_visco_dx / host.crust.thickness, c=rheo_color, ls='--')
        axes_1[1, 0].plot(x, host_ocean_dx / host.crust.thickness, c=rheo_color, ls='-.')
        axes_1[1, 0].plot(x, target_ocean_dx / host.crust.thickness, c=rheo_color, ls=':')

        # Pane 4: Heating
        host_tidal_heating = all_data['host']['tidal_heating']
        target_tidal_heating  = all_data['target']['tidal_heating']
        host_radiogenics = all_data['host']['radiogenics']
        target_radiogenics = all_data['target']['radiogenics']
        axes_1[1, 1].plot(x, host_tidal_heating, c=rheo_color, ls='-')
        axes_1[1, 1].plot(x, target_tidal_heating, c=rheo_color, ls='--')
        axes_1[1, 1].plot(x, host_radiogenics, c=rheo_color, ls='-.')
        axes_1[1, 1].plot(x, target_radiogenics, c=rheo_color, ls=':')

    plt.tight_layout()
    plt.show()
    return all_data_by_rheo, fig_1


if __name__ == '__main__':
    main()