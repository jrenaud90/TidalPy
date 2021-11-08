import numpy as np

from TidalPy.utilities.types import float_eps
from TidalPy.constants import G
from TidalPy.utilities.performance import njit
from TidalPy.toolbox.conversions import orbital_motion2semi_a
from TidalPy.cooling.cooling_models import convection, conduction
from TidalPy.rheology.partial_melt.partialmelt import calculate_melt_fraction
from TidalPy.rheology.viscosity.viscosity_models import reference, arrhenius
from TidalPy.tides.mode_manipulation import find_mode_manipulators
from TidalPy.tides import calc_tidal_susceptibility
from TidalPy.rheology.complex_compliance.complex_compliance import compliance_dict_helper

# Hussmann Functions
@njit(cacheable=True)
def shear_falloff(temperature, premelt_shear, solidus_temperature, param_0, param_1):

    if temperature < solidus_temperature:
        return premelt_shear
    else:
        return param_0 * np.exp(param_1 / temperature)


# Model switches
model_use_obliquity = False
model_eccentricity_truncation = 10
model_max_tidal_order_l = 2
model_use_tidal_scale = True
model_restrict_tidal_scale_to_visco_layer_volume = True
model_use_planetary_params_for_tide_calc = False

# Planet Properties
host_mass = 5.683e26
target_name = 'Enceladus'
target_radius = 252.1e3
target_mass = 1.08e20
target_volume = (4. / 3.) * np.pi * target_radius**3
target_density_bulk = target_mass / target_volume
target_gravity = G * target_mass / target_radius**2
surface_temperature = 72.15

# Layer Properties
# # Crust
crust_visco_solidus = 272.0
crust_visco_liquidus = 273.15
# # Ice Properties
ice_specific_heat = 1925.
ice_latent_heat = 284000.0
ice_thermal_conductivity = 2.27
ice_thermal_expansion = 1.56e-4
ice_density = 1000.
water_density = 950.
ice_thermal_diffusivity = ice_thermal_conductivity / (ice_density * ice_specific_heat)
ice_visco_tempdrop = 1. / (1. + (np.log(10) / 27.))
ice_convection_alpha = 1.0
ice_convection_beta = 1. / 3.
ice_critical_rayleigh = 900.
ice_elastic_threshold_temp = ice_visco_tempdrop * crust_visco_solidus
# # Geometry and Structure
crust_radius = target_radius
crust_total_thickness = (35. + 30.) * 1.e3
crust_sa = 4. * np.pi * crust_radius**2
crust_volume = (4. / 3.) * np.pi * (crust_radius**3 - (crust_radius - crust_total_thickness)**3)
crust_pressure = 0.0
crust_estimated_mass = np.average((water_density, ice_density)) * crust_volume

# # Core
# # Silicate Properties
silicate_convection_alpha = 1.0
silicate_convection_beta = 1. / 3.
silicate_critical_rayleigh = 1000.
silicate_density = 3400.
silicate_thermal_expansion = 5.0e-5
silicate_thermal_conductivity = 4.0
silicate_specific_heat = 1225.5
# silicate_latent_heat = NA
# # Geometry and Structure
core_radius = target_radius - crust_total_thickness
core_thickness = core_radius
core_sa = 4. * np.pi * core_radius**2
core_volume = (4. / 3.) * np.pi * core_radius**3
# TODO: Check how this density compares with standard silicate densities.
silicate_density = (1. / core_volume) * (target_mass - crust_estimated_mass)
silicate_thermal_diffusivity = silicate_thermal_conductivity / (silicate_density * silicate_specific_heat)
core_solidus = 1600.0
core_liquidus = 2000.0
core_pressure = 0.0
core_mass = core_volume * silicate_density
core_gravity = G * core_mass / core_radius**2

# # Debug Flags (set all to false for normal operation)
crust_viscoelastic_layer_passes_all_heat = False
crust_elastic_layer_passes_all_heat = False
core_passes_all_heat = False
dont_calculate_tides = False

# Prep and setup various functions
calculate_tidal_terms, collapse_modes, eccentricity_func, inclination_func = \
    find_mode_manipulators(model_max_tidal_order_l, model_eccentricity_truncation, use_obliquity=model_use_obliquity)

@njit(cacheable=True)
def diffeq_scipy(time, variables, extra_params):

    # Variables - eccentricity, orbital_motion, rotation_rate, core_temperature, crust_temperature, crust_visco_dx, crust_elastic_dx
    eccentricity = variables[0]
    obliquity = variables[1]
    orbital_motion = variables[2]
    rotation_rate = variables[3]
    core_temperature = variables[4]
    crust_visco_temp_from_dt = variables[5]
    crust_visco_dx = variables[6]
    crust_elastic_dx = variables[7]

    # Reset Flags
    is_fully_stagnant_lid = False
    is_stagnant_lid_gone = False
    is_ocean_present = False
    is_visco_grounded_out = False

    # Fix unexpected values
    if crust_elastic_dx <= float_eps:
        crust_elastic_dx = float_eps
    if crust_visco_dx <= float_eps:
        crust_visco_dx = float_eps
    if core_temperature <= 10.:
        core_temperature = 10.
    if crust_visco_temp_from_dt <= 10.:
        crust_visco_temp_from_dt = 10.
    if eccentricity >= 1.:
        eccentricity = 0.99
    elif eccentricity <= float_eps:
        eccentricity = 0.

    # Other Parameters
    fixed_eccentricity = extra_params[0]
    static_core_heating = extra_params[1]
    static_crust_heating = extra_params[2]

    # Derived Parameters
    semi_major_axis = orbital_motion2semi_a(orbital_motion, host_mass, target_mass)

    # Calculate viscosity and shear of layers
    # # Core
    core_viscosity = 1.0e22
    core_shear = 50.0e9
    # core_melt_fraction = calculate_melt_fraction(core_temperature, core_solidus, core_liquidus)
    # core_premelt_viscosity = core_viscosity_model(core_temperature, core_pressure, *core_viscosity_inputs)
    # core_viscosity, core_shear_modulus = core_partialmelt_model(core_melt_fraction, core_temperature)

    # # Viscoelastic Ice
    crust_viscosity = 1.0e18
    crust_shear = 3.3e9
    # crust_melt_fraction = calculate_melt_fraction(crust_visco_temperature, crust_visco_solidus, crust_visco_liquidus)

    # Target Geometry
    # Check if ocean is present
    if (crust_elastic_dx + crust_visco_dx) < crust_total_thickness:
        # The icy layers do not add up to the total crust thickness, the ocean is the only remaining portion.
        is_ocean_present = True

    # Check if the layer is totally elastic or if the elastic layer is totally gone
    if (crust_elastic_dx >= crust_total_thickness) or \
            abs(crust_visco_temp_from_dt - ice_elastic_threshold_temp) <= float_eps:
        # Layer is all elastic. No viscoelastic portion, make it small so we don't have divide by zero issues
        is_fully_stagnant_lid = True
        crust_visco_dx = 0.0
        # Correct overshoot
        crust_elastic_dx = crust_total_thickness

    elif crust_elastic_dx <= float_eps:
        # Elastic layer is non-existent. Make it small so we don't have divide by zero issues
        is_stagnant_lid_gone = True
        crust_elastic_dx = 0.0

    if is_fully_stagnant_lid:
        # There is no viscoelastic ice layer. The viscoelastic temperature won't matter so set it to nan.
        crust_viscoelastic_temperature = 0.0
        crust_viscoelastic_bottom_temperature = 0.0
        crust_viscoelastic_top_temperature = 0.0
    else:
        if is_ocean_present:
            # There is a viscoelastic layer and an underlying ocean.
            # # Determine temperature within viscoelastic ice
            crust_viscoelastic_bottom_temperature = crust_visco_solidus
            crust_viscoelastic_top_temperature = crust_viscoelastic_bottom_temperature * ice_visco_tempdrop
            crust_viscoelastic_temperature = \
                0.5 * (crust_viscoelastic_top_temperature + crust_viscoelastic_bottom_temperature)
        else:
            # Grounded Out: There is a viscoelastic layer, but no ocean. Viscoelastic layer is touching the top of
            #    the silicate core.
            is_visco_grounded_out = True
            # The temperature of the viscoelastic layer is determined by a differential equation.
            # Assume the bottom of the viscoelastic layer is what the differential equation is tracking
            crust_viscoelastic_bottom_temperature = crust_visco_temp_from_dt
            crust_viscoelastic_top_temperature = crust_viscoelastic_bottom_temperature * ice_visco_tempdrop
            crust_viscoelastic_temperature = \
                0.5 * (crust_viscoelastic_top_temperature + crust_viscoelastic_bottom_temperature)

    # Determine icy layer geometry and structure
    crust_visco_radius = target_radius - crust_elastic_dx
    crust_ocean_radius = crust_visco_radius - crust_visco_dx
    crust_ocean_volume = (4. / 3.) * np.pi * (crust_ocean_radius**3 - core_radius**3)
    crust_elastic_volume = (4. / 3.) * np.pi * (target_radius**3 - crust_visco_radius**3)
    crust_visco_volume = (4. / 3.) * np.pi * (crust_visco_radius**3 - crust_ocean_radius**3)
    crust_elastic_sa = 4. * np.pi * target_radius**2
    crust_visco_sa = 4. * np.pi * crust_visco_radius**2
    crust_elastic_gravity = target_gravity
    mass_below_visco = core_mass + (crust_ocean_volume * water_density) + (crust_visco_volume * ice_density)
    crust_visco_gravity = G * mass_below_visco / crust_visco_radius**2

    # Update tidal scales
    if model_use_tidal_scale:
        core_tidal_scale = core_volume / target_volume
        crust_tidal_scale = crust_visco_volume / target_volume
    else:
        core_tidal_scale = 1.
        crust_tidal_scale = 1.

    # Calculate heat sources
    # # Radiogenic
    core_radio_heating = 0.0
    crust_radio_heating = 0.0
    # # Tidal
    core_tidal_heating = 0.0
    crust_tidal_heating = 0.0
    # # Static
    core_heating = core_radio_heating + core_tidal_heating + static_core_heating
    crust_visco_invitro_heating = crust_radio_heating + crust_tidal_heating + static_crust_heating

    # # Calculate orbital and rotational changes
    # Calculate eccentricity and obliquity functions
    eccentricity_results = eccentricity_func(eccentricity)
    obliquity_results = inclination_func(obliquity)

    # Calculate tidal modes and susceptibility
    unique_frequencies, tidal_results_by_frequency = \
        calculate_tidal_terms(rotation_rate, orbital_motion, semi_major_axis, target_radius,
                              eccentricity_results, obliquity_results)

    # # Calculate tidal dissipation
    tidal_susceptibility = \
        calc_tidal_susceptibility(host_mass, target_radius, semi_major_axis)

    if dont_calculate_tides:
        dUdM_total = 0.
        dUdw_total = 0.
        dUdO_total = 0.
        tidal_heating_global = 0.
        core_tidal_heating = 0.
        crust_tidal_heating = 0.
    else:
        # Calculate tides!

        # Calculate the rheological response of the core and crust


        complex_compliance_func = complex_compliance_funcs[object_i][layer_i]
        rheology_input = rheology_inputs[object_i][layer_i]
        # Calculate complex compliance based on layer's strength and the unique tidal forcing frequencies

        unique_complex_compliances = \
            compliance_dict_helper(unique_frequencies, core_compliance_func,
                                   (core_shear**(-1), core_viscosity), core_compliance_input)

        # Calculate tidal dissipation
        if model_use_planetary_params_for_tide_calc:
            core_tidal_radius = target_radius
            crust_tidal_radius = target_radius
            core_tidal_gravity = target_gravity
            crust_tidal_gravity = target_gravity
            core_tidal_density = target_density_bulk
            crust_tidal_density = target_density_bulk
        else:
            core_tidal_radius = core_radius
            crust_tidal_radius = crust_visco_radius
            core_tidal_gravity = core_gravity
            crust_tidal_gravity = crust_visco_gravity
            core_tidal_density = silicate_density
            crust_tidal_density = ice_density

        # Calculate tidal dissipation for the core and crust
        core_tidal_heating, core_dUdM, core_dUdw, core_dUdO, core_love_number_by_order_l,\
            core_negative_imk_by_order_l, core_effective_q_by_order_l = \
                collapse_modes(core_tidal_gravity, core_tidal_radius, core_tidal_density,
                               core_shear, core_tidal_scale,
                               host_mass, tidal_susceptibility,
                               core_unique_complex_compliances,
                               tidal_results_by_frequency,
                               model_max_tidal_order_l, cpl_ctl_method=False)

        crust_tidal_heating, crust_dUdM, crust_dUdw, crust_dUdO, crust_love_number_by_order_l, \
            crust_negative_imk_by_order_l, crust_effective_q_by_order_l = \
                collapse_modes(crust_tidal_gravity, crust_tidal_radius, crust_tidal_density,
                               crust_shear, crust_tidal_scale,
                               host_mass, tidal_susceptibility,
                               crust_unique_complex_compliances,
                               tidal_results_by_frequency,
                               model_max_tidal_order_l, cpl_ctl_method=False)

        # Combine the tidal results linearly (note the tidal_scale scales each layer's values if enabled)
        dUdM_total = core_dUdM + crust_dUdM
        dUdw_total = core_dUdw + crust_dUdw
        dUdO_total = core_dUdO + crust_dUdO
        tidal_heating_global = core_tidal_heating + crust_tidal_heating

    # Calculate heat balance and cooling (Start from the bottom to the top)
    # # Core
    # TODO: Right now the model does not allow for a stagnant lid on the silicate core
    # Determine the temperature at the top
    if is_ocean_present:
        core_top_temperature = crust_visco_liquidus
    else:
        core_top_temperature = crust_viscoelastic_bottom_temperature

    core_heat_flux = core_heating / core_sa
    if core_passes_all_heat:
        core_cooling_flux = core_heat_flux
        core_delta_temp = 0.
        core_boundary_dx = float_eps
        core_rayleigh = 0.
        core_nusselt = 1.
        core_cooling = core_heating
    else:
        # Perform convection calculations
        core_delta_temp = core_temperature - core_top_temperature
        core_cooling_flux, core_boundary_dx, core_rayleigh, core_nusselt = \
            convection(
                core_delta_temp, core_viscosity,
                silicate_thermal_conductivity, silicate_thermal_diffusivity, silicate_thermal_expansion,
                core_thickness, core_gravity, silicate_density,
                silicate_convection_alpha, silicate_convection_beta, silicate_critical_rayleigh
            )
        core_cooling = core_cooling_flux * core_sa

    # Update core temperature
    if core_passes_all_heat:
        # If the debug flag is on then ignore the change in core temperature
        dt_core_temperature = 0.0
    else:
        dt_core_temperature = \
            (core_heating - core_cooling) / (core_volume * silicate_density * silicate_specific_heat)

    # # Crust
    crust_visco_heating = crust_visco_invitro_heating + core_cooling
    # Grow or shrink ice layers based on temperature balance
    crust_visco_base_heat_flux = crust_visco_heating / crust_visco_sa
    if is_fully_stagnant_lid or crust_viscoelastic_layer_passes_all_heat:
        # Skip convection calculation if the viscoelastic layer is gone or a debug switch is on
        crust_visco_cooling_flux = crust_visco_base_heat_flux
        crust_visco_delta_temp = 0.
        crust_visco_boundary_dx = float_eps
        crust_visco_rayleigh = 0.
        crust_visco_nusselt = 1.
        crust_visco_cooling = crust_visco_heating
    else:
        # Perform convection calculations
        crust_visco_delta_temp = crust_viscoelastic_bottom_temperature - crust_viscoelastic_top_temperature
        crust_visco_cooling_flux, crust_visco_boundary_dx, crust_visco_rayleigh, crust_visco_nusselt = \
            convection(
                crust_visco_delta_temp, crust_viscosity,
                ice_thermal_conductivity, ice_thermal_diffusivity, ice_thermal_expansion,
                crust_visco_dx, crust_visco_gravity, ice_density,
                ice_convection_alpha, ice_convection_beta, ice_critical_rayleigh
            )
        crust_visco_cooling = crust_visco_cooling_flux * crust_visco_sa

    # Heating entering the base of the elastic layer is equal to that exiting the viscoelastic one
    crust_elastic_heating = crust_visco_cooling_flux * crust_visco_sa
    crust_elastic_base_heat_flux = crust_visco_cooling_flux
    if is_stagnant_lid_gone or crust_elastic_layer_passes_all_heat:
        # There is either no elastic layer or all of the heat is leaving it instantly due to debug switch
        crust_elastic_cooling_flux = crust_elastic_heating / crust_elastic_sa
        crust_elastic_delta_temp = 0.
        crust_elastic_boundary_dx = float_eps
        crust_elastic_rayleigh = 0.
        crust_elastic_nusselt = 1.
        crust_elastic_cooling = crust_elastic_heating
    else:
        # There is an elastic layer, calculate conductive flux
        crust_elastic_delta_temp = crust_viscoelastic_top_temperature - surface_temperature
        crust_elastic_cooling_flux, crust_elastic_boundary_dx, crust_elastic_rayleigh, crust_elastic_nusselt = \
                conduction(crust_elastic_delta_temp, ice_thermal_conductivity, crust_elastic_dx)
        crust_elastic_cooling = crust_elastic_cooling_flux * crust_elastic_sa

    # Calculate the change in the elastic layer's thickness
    crust_visco_energy_to_freeze = ice_specific_heat * crust_visco_delta_temp
    if crust_elastic_layer_passes_all_heat:
        # If the debug flag is on then ignore elastic layer growth.
        dt_crust_elastic_dx = 0.0
    else:
        dt_crust_elastic_dx = \
            (crust_elastic_cooling_flux - crust_elastic_base_heat_flux) / (ice_density * crust_visco_energy_to_freeze)
    # Check for over/undershoots
    if is_stagnant_lid_gone and dt_crust_elastic_dx < 0.:
        # Elastic layer is already gone, make sure it is not growing negative
        #    (this shouldn't happen, just here as a double check).
        dt_crust_elastic_dx = 0.0

    if is_fully_stagnant_lid:
        # There is no viscoelastic layer. No change in its thickness, no change in its temperature
        dt_crust_visco_dx = 0.0
        dt_crust_visco_temperature = 0.0
    elif is_ocean_present:
        # There is an ocean present. We will need to calculate the change in the viscoelastic layer's thickness,
        #    but do not need to calculate the change in temperature
        dt_crust_visco_temperature = 0.0
        dt_crust_visco_dx = \
                (crust_visco_cooling_flux - crust_visco_base_heat_flux) / (ice_density * ice_latent_heat)
    elif is_visco_grounded_out:
        # There is no more ocean so the only way that the viscoelastic layer can change is be stealing or giving up
        #   ground to the elastic layer
        dt_crust_visco_dx = -dt_crust_elastic_dx
        # The temperature of this layer is now allowed to decrease via a differential equation until it too becomes
        #    totally elastic.
        dt_crust_visco_temperature = \
            (crust_visco_heating - crust_visco_cooling) / (crust_visco_volume * ice_density * ice_specific_heat)
    else:
        # How did you end up here.
        raise Exception

    # Build output array
    output_variables = (
        dt_eccentricity,
        dt_obliquity,
        dt_orbital_motion,
        dt_rotation_rate,
        dt_core_temperature,
        dt_crust_visco_temperature,
        dt_crust_visco_dx,
        dt_crust_elastic_dx
    )

    return output_variables




    if is_visco_grounded_out:
        # There is no more ocean so the only way that the viscoelastic layer can change is be stealing or giving up
        #   ground to the elastic layer
        dt_crust_visco_layer_dx = -dt_crust_elastic_layer_dx


        if layer_freeze_out:
            # The viscoelastic layer's temperature is allowed to decrease if the layer is frozen out.
            viscoelastic_layer_change = -elastic_layer_change
            viscoelastic_temperature_change = \
                (total_incoming_heating - viscoelastic_cooling)/ \
                (viscoelastic_volume*material_density*specific_heat)
        else:
            viscoelastic_temperature_change = 0.
            # If there is an ocean then the viscoelastic layer will change depending on the
            #     energy balance
            viscoelastic_layer_change = \
                (viscoelastic_cooling_flux - viscoelastic_heating_flux)/ \
                (material_density*latent_heat)

    wth_layer:

    # Determine heating fluxes
    viscoelastic_heating_flux = total_incoming_heating/viscoelastic_surf_area

    if viscoelastic_passes_all_heat:
        viscoelastic_cooling_flux = viscoelastic_heating_flux

    elastic_heating_flux = viscoelastic_cooling_flux

    if elastic_passes_all_heat:
        elastic_cooling_flux = viscoelastic_cooling_flux
        elastic_cooling = viscoelastic_cooling_flux*elastic_surf_area
        elastic_boundary_layer_thickness = 0.
        elastic_rayleigh = 0.
        elastic_nusselt = 1.
    else:
        # Determine Cooling for the elastic layer
        if top_layer:
            elastic_delta_temperature = visco_top_temp - surface_temperatures[object_i]
        else:
            elastic_delta_temperature = visco_top_temp - temperatures[object_i][layer_i + 1]
        #     conduction function
        elastic_cooling_flux, elastic_boundary_layer_thickness, elastic_rayleigh, elastic_nusselt = \
            conduction(elastic_delta_temperature, thermal_conductivity, elastic_thickness)
        elastic_cooling = elastic_surf_area*elastic_cooling_flux

    if layer_stagnant:
        # If the layer is stagnant (all elastic) then there is no change in temperature.
        elastic_layer_change = 0.
        viscoelastic_layer_change = 0.
        viscoelastic_temperature_change = 0.
    else:
        energy_to_freeze = specific_heat*viscoelastic_temperature_delta
        elastic_layer_change = \
            (elastic_cooling_flux - elastic_heating_flux)/(material_density*energy_to_freeze)
        if layer_freeze_out:
            # The viscoelastic layer's temperature is allowed to decrease if the layer is frozen out.
            viscoelastic_layer_change = -elastic_layer_change
            viscoelastic_temperature_change = \
                (total_incoming_heating - viscoelastic_cooling)/ \
                (viscoelastic_volume*material_density*specific_heat)
        else:
            viscoelastic_temperature_change = 0.
            # If there is an ocean then the viscoelastic layer will change depending on the
            #     energy balance
            viscoelastic_layer_change = \
                (viscoelastic_cooling_flux - viscoelastic_heating_flux)/ \
                (material_density*latent_heat)

    # Heat coming out of the layer is equal to the elastic cooling
    layer_cooling = elastic_cooling


    # Target Planet Cooling - Work from button to top
    # # Core
    if is_ocean_present:
        core_top_temp = ocean_temperature
    else:
        if is_fully_stagnant_lid:
            # TODO
        else:
            # There is a viscoelastic layer, use its temperature for core convection calculations
            core_top_temp = crust_visco_temperature

    core_delta_temp = core_temperature - core_top_temp
    if core_delta_temp > 0:

    core_cooling = convection()

    # # Crust - Assume all heat flux bypasses any present liquid ocean and enters the viscoelastic region


    core_heating = core_tidal_heating + core_radiogenic_heating + core_other_heating

    crust_visco_heat = core_cooling + crust_visco_tidal_heating + crust_visco_other_heating
    crust_visco_heat_flux = crust_visco_heating / crust_visco_sa








    # # Viscoelastic Crust
    crust_visco_temperature_delta = crust_visco_temperature - visco_top_temp
    crust_visco_cooling_flux, crust_visco_boundary_layer_thickness, crust_visco_rayleigh, crust_visco_nusselt = \
        convection(crust_visco_temperature_delta, crust_viscosity, ice_thermal_conductivity,
                   ice_thermal_diffusivity, ice_thermal_expansion,
                   crust_visco_dx, viscoelastic_gravity, material_density, alpha_conv,
                   beta_conv, critical_rayleigh)

    # Calculate viscoelastic cooling
    crust_elastic_top_temp = planet_surface_temp
    crust_visco_top_temp =


    if is_growth_layer:
        if layer_lost_lid:
            if top_layer:
                visco_top_temp = surface_temperatures[object_i]
            else:
                visco_top_temp = bottom_temperatures[layer_i + 1]
        else:
            visco_top_temp = viscoelastic_top_temperatures[object_i][layer_i]
    else:
        if top_layer:
            visco_top_temp = surface_temperatures[object_i]
        else:
            visco_top_temp = bottom_temperatures[layer_i + 1]


    # # Elastic Crust




    # # Progress bar
    # percent_done = round(100000. * (time / time_interval)) / 100000.
    # print('Percent Done:', 100. * percent_done, '%')
    # print('\rPercent Done: {:0>5.2f}%'.format(100. * percent_done), flush=True, end='')

    # Pull out independent variables
    temperatures = list()
    elastic_dxs = list()
    viscoelastic_dxs = list()
    spin_rates = list()
    index = 0
    for object_i in range(2):
        n_layers = num_layers[object_i]
        temperatures_by_layer = list()
        elastic_dx_by_layer = list()
        viscoelastic_dx_by_layer = list()
        for layer_i in range(n_layers):
            temperatures_by_layer.append(variables[index])
            elastic_dx_by_layer.append(variables[index + 1])
            viscoelastic_dx_by_layer.append(variables[index + 2])
            index += 3
        temperatures.append(temperatures_by_layer)
        elastic_dxs.append(elastic_dx_by_layer)
        viscoelastic_dxs.append(viscoelastic_dx_by_layer)
        # Next item is the spin-rate
        spin_rates.append(variables[index])
        index += 1
    orbital_motion = variables[index]
    eccentricity = variables[index + 1]

    # Calculate eccentricity functions
    eccentricity_results = eccentricity_func(eccentricity)

    # Calculate parameters that only depend on orbital properties
    semi_major_axis = orbital_motion2semi_a(orbital_motion, object_masses[0], object_masses[1])

    # Derivative storage
    tidal_derivative_storage = list()
    derivative_storage = list()

    # Perform thermal evolution on each object's layers
    for object_i in range(2):
        number_of_layers = num_layers[object_i]

        # Need opposite object index for tidal calculations
        opposite_object_i = 1
        if object_i == 1:
            opposite_object_i = 0

        # Pull out parameters referenced often
        object_radius = object_radii[object_i]
        object_mass = object_masses[object_i]
        object_gravity = object_gravities[object_i]
        object_density_bulk = object_densities_bulk[object_i]
        tidal_host_mass = object_masses[opposite_object_i]
        spin_rate = spin_rates[object_i]
        spin_locked = False
        if lock_at_1to1:
            if (abs((spin_rate / orbital_motion) - 1.) < 0.01) and eccentricity <= 0.05:
                spin_locked = True
                spin_rate = orbital_motion

        # Calculate obliquity results
        # TODO: For now obliquity is not tracked so it can be calculated outside the main function
        obliquity_results = inclination_func(object_obliquities[object_i])

        # Calculate tidal modes and susceptibility
        unique_frequencies, tidal_results_by_frequency = \
            calculate_tidal_terms(spin_rate, orbital_motion, semi_major_axis, object_radius,
                                  eccentricity_results, obliquity_results)

        tidal_susceptibility = \
            calc_tidal_susceptibility(tidal_host_mass, object_radius, semi_major_axis)

        # Tidal parameter storage
        dUdM_total = 0.
        dUdw_total = 0.
        dUdO_total = 0.

        # First loop through layers to determine geometry and find bottom temperatures which are used for next loop
        bottom_temperatures = list()
        layer_geometries = list()
        for layer_i in range(number_of_layers):
            bottom_layer = layer_i == 0
            top_layer = layer_i == num_layers[object_i] - 1
            # Pull out often used parameters
            #    State Properties
            viscoelastic_temperature = temperatures[object_i][layer_i]
            layer_radius_upper = layer_radii_upper[object_i][layer_i]
            layer_radius_lower = layer_radii_lower[object_i][layer_i]
            #    Layer Properties
            layer_mass_below = layer_masses_below[object_i][layer_i]
            #    Material Properties
            material_density = material_densities[object_i][layer_i]

            # Determine layer geometry
            layer_stagnant = False
            layer_freeze_out = False
            layer_lost_lid = False
            is_growth_layer = growth_layer_flags[object_i][layer_i]

            if is_growth_layer:
                elastic_radius_upper = layer_radius_upper
                elastic_radius_lower = layer_radius_upper - elastic_dxs[object_i][layer_i]

                # Check if the layer is totally elastic or if the elastic layer is totally gone
                if (elastic_radius_lower <= layer_radius_lower) or \
                        abs(viscoelastic_temperature - viscoelastic_top_temperatures[object_i][layer_i]) < 1.:
                    elastic_radius_lower = layer_radius_lower
                    # Layer is all elastic. No viscoelastic portion
                    layer_stagnant = True

                elif elastic_radius_lower >= layer_radius_upper:
                    # Elastic layer is basically non-existent. Make it small so we don't have divide by zero issues
                    elastic_radius_lower = layer_radius_upper
                    elastic_radius_upper = layer_radius_upper
                    layer_lost_lid = True

                if layer_stagnant:
                    # Viscoelastic layer is non-existent
                    viscoelastic_radius_lower = layer_radius_lower
                    viscoelastic_radius_upper = layer_radius_lower
                    if top_layer:
                        # TODO: this is not right, but we don't really care about the thermal evolution after freezeout.
                        viscoelastic_temperature = 0.5 * (surface_temperatures[object_i] +
                                                          constant_viscoelastic_temperatures[object_i][layer_i])
                    else:
                        # TODO: And this should really be a based off the bottom temperature of the above layer.
                        viscoelastic_temperature = constant_viscoelastic_temperatures[object_i][layer_i]
                    bottom_temperature = viscoelastic_temperature
                else:
                    # Viscoelastic layer is present
                    viscoelastic_radius_upper = elastic_radius_lower
                    viscoelastic_radius_lower = elastic_radius_lower - viscoelastic_dxs[object_i][layer_i]
                    if viscoelastic_radius_lower < layer_radius_lower:
                        viscoelastic_radius_lower = layer_radius_lower
                        # Layer is all either elastic or viscoelastic. No ocean
                        ocean_radius_lower = layer_radius_lower
                        ocean_radius_upper = layer_radius_lower
                        layer_freeze_out = True
                        # Bottom temperature is now equal to the viscoelastic temperature
                        bottom_temperature = viscoelastic_temperature
                    else:
                        # Ocean layer is still present
                        ocean_radius_lower = layer_radius_lower
                        ocean_radius_upper = viscoelastic_radius_lower
                        # If the ocean layer is present then the viscoelastic temperature is a constant
                        viscoelastic_temperature = constant_viscoelastic_temperatures[object_i][layer_i]
                        # And bottom temperature is the ocean temperature
                        bottom_temperature = constant_ocean_temperatures[object_i][layer_i]

                # Calculate derivec properties
                elastic_thickness = elastic_radius_upper - elastic_radius_lower
                elastic_volume = (4. / 3.) * np.pi * \
                    (elastic_radius_upper * elastic_radius_upper * elastic_radius_upper -
                     elastic_radius_lower * elastic_radius_lower * elastic_radius_lower)
                elastic_surf_area = 4. * np.pi * (elastic_radius_upper * elastic_radius_upper)
                elastic_mass = elastic_volume * material_density
                elastic_mass_below = layer_mass_below
                elastic_gravity = G * elastic_mass_below / \
                                  (elastic_radius_upper * elastic_radius_upper)

                viscoelastic_thickness = viscoelastic_radius_upper - viscoelastic_radius_lower
                viscoelastic_volume = (4. / 3.) * np.pi * \
                    (viscoelastic_radius_upper * viscoelastic_radius_upper * viscoelastic_radius_upper -
                     viscoelastic_radius_lower * viscoelastic_radius_lower * viscoelastic_radius_lower)
                viscoelastic_surf_area = 4. * np.pi * (viscoelastic_radius_upper * viscoelastic_radius_upper)
                viscoelastic_mass = viscoelastic_volume * material_density
                viscoelastic_mass_below = layer_mass_below - elastic_mass
                viscoelastic_gravity = G * viscoelastic_mass_below / \
                                       (viscoelastic_radius_upper * viscoelastic_radius_upper)

            else:
                # For the non-growth model, the entire layer is assumed to be viscoelastic (no stagnant lid)
                viscoelastic_radius_lower = layer_radii_lower[object_i][layer_i]
                viscoelastic_radius_upper = layer_radii_upper[object_i][layer_i]
                viscoelastic_thickness = layer_thicknesses[object_i][layer_i]
                viscoelastic_volume = layer_volumes[object_i][layer_i]
                viscoelastic_surf_area = layer_surf_areas[object_i][layer_i]
                viscoelastic_mass = layer_masses[object_i][layer_i]
                viscoelastic_gravity = layer_gravities[object_i][layer_i]

                # Elastic parameters will not be used.
                elastic_radius_lower = 0.
                elastic_radius_upper = 0.
                elastic_thickness = 0.
                elastic_volume = 0.
                elastic_surf_area = 0.
                elastic_mass = 0.
                elastic_gravity = 0.

                # Determine the temperature at the base of the layer (used for other layer's thermal evolution)
                bottom_temperature = viscoelastic_temperature

            if use_tidal_scale:
                if use_visco_volume_for_tidal_scale:
                    layer_tidal_scale = viscoelastic_volume / object_volumes[object_i]
                else:
                    layer_tidal_scale = tidal_scales[object_i][layer_i]
            else:
                layer_tidal_scale = 1.

            bottom_temperatures.append(bottom_temperature)
            layer_geometries.append(
                    (layer_tidal_scale, layer_freeze_out, layer_stagnant, layer_lost_lid,
                     viscoelastic_radius_lower, viscoelastic_radius_upper, viscoelastic_thickness,
                     viscoelastic_volume, viscoelastic_surf_area, viscoelastic_mass, viscoelastic_gravity,
                     viscoelastic_temperature,
                     elastic_radius_lower, elastic_radius_upper, elastic_thickness,
                     elastic_volume, elastic_surf_area, elastic_mass, elastic_gravity)
            )

        # Storage for various parameters that need to be accessed during the layer looping
        layer_coolings = list()

        # Now that layer's bottom temperatures are known we can loop through again and calculate everything else.
        for layer_i in range(number_of_layers):
            bottom_layer = layer_i == 0
            top_layer = layer_i == num_layers[object_i] - 1
            is_growth_layer = growth_layer_flags[object_i][layer_i]

            # Pull out often used parameters
            #    Material Properties
            material_density = material_densities[object_i][layer_i]
            thermal_conductivity = thermal_conductivities[object_i][layer_i]
            thermal_expansion = thermal_expansions[object_i][layer_i]
            thermal_diffusivity = thermal_diffusivities[object_i][layer_i]
            specific_heat = specific_heats[object_i][layer_i]
            latent_heat = latent_heats[object_i][layer_i]
            solidus = solidus_temperature[object_i][layer_i]
            liquidus = liquidus_temperature[object_i][layer_i]
            #    Convection Properties
            alpha_conv = convection_alphas[object_i][layer_i]
            beta_conv = convection_betas[object_i][layer_i]
            critical_rayleigh = critical_rayleighs[object_i][layer_i]

            # Pull out geometry information that was just calculated in the previous loop
            layer_tidal_scale, layer_freeze_out, layer_stagnant, layer_lost_lid, \
            viscoelastic_radius_lower, viscoelastic_radius_upper, viscoelastic_thickness, \
            viscoelastic_volume, viscoelastic_surf_area, viscoelastic_mass, viscoelastic_gravity, \
            viscoelastic_temperature, \
            elastic_radius_lower, elastic_radius_upper, elastic_thickness, \
            elastic_volume, elastic_surf_area, elastic_mass, elastic_gravity = \
                layer_geometries[layer_i]

            elastic_passes_all_heat = False
            if layer_lost_lid:
                elastic_passes_all_heat = True

            viscoelastic_passes_all_heat = False
            # Calculate viscoelastic portion's strength
            if viscoelastic_volume != 0.:
            #    TODO: Assume 0. pressure dependence for now (the zero in the input below)
                pre_melt_viscosity = \
                    viscosity_funcs[object_i][layer_i](viscoelastic_temperature, 0., *viscosity_inputs[object_i][layer_i])
                pre_melt_shear = static_shears[object_i][layer_i]
                #    Apply partial melting
                melt_fraction = (viscoelastic_temperature - solidus) / (liquidus - solidus)
                if melt_fraction < 0.:
                    melt_fraction = 0.
                elif melt_fraction > 1.:
                    melt_fraction = 1.
                viscosity, shear_modulus = \
                    partial_melt_funcs[object_i][layer_i](melt_fraction, pre_melt_viscosity, pre_melt_shear,
                                                          *partial_melt_inputs[object_i][layer_i])
                #    Check for over/undershoots
                if viscosity < 0.1:
                    viscosity = 0.1
                elif viscosity > 1.e100:
                    viscosity = 1.e100
                if shear_modulus < 0.1:
                    shear_modulus = 0.1
                elif shear_modulus > 1.e100:
                    shear_modulus = 1.e100
                compliance = 1. / shear_modulus

                # Calculate viscoelastic cooling
                if is_growth_layer:
                    if layer_lost_lid:
                        if top_layer:
                            visco_top_temp = surface_temperatures[object_i]
                        else:
                            visco_top_temp = bottom_temperatures[layer_i + 1]
                    else:
                        visco_top_temp = viscoelastic_top_temperatures[object_i][layer_i]
                else:
                    if top_layer:
                        visco_top_temp = surface_temperatures[object_i]
                    else:
                        visco_top_temp = bottom_temperatures[layer_i + 1]

                # Convective cooling for viscoelastic layer
                viscoelastic_temperature_delta = viscoelastic_temperature - visco_top_temp
                viscoelastic_cooling_flux, viscoelastic_boundary_layer_thickness, viscoelastic_rayleigh, \
                    viscoelastic_nusselt = \
                        convection(viscoelastic_temperature_delta, viscosity, thermal_conductivity,
                                   thermal_diffusivity, thermal_expansion,
                                   viscoelastic_thickness, viscoelastic_gravity, material_density, alpha_conv,
                                   beta_conv, critical_rayleigh)
                viscoelastic_cooling = viscoelastic_surf_area * viscoelastic_cooling_flux


                # Tidal Calculations
                calculate_tides = (not force_tides_off_flags[object_i][layer_i]) and tides_on_flags[object_i] \
                                  and (layer_tidal_scale != 0.)

                if not calculate_tides:
                    tidal_heating = 0.
                    dUdM, dUdw, dUdO = 0., 0., 0.0
                    love_number, negative_imk = 0. + 0.j, 0.
                else:
                    # Calculate tides!
                    complex_compliance_func = complex_compliance_funcs[object_i][layer_i]
                    rheology_input = rheology_inputs[object_i][layer_i]
                    # Calculate complex compliance based on layer's strength and the unique tidal forcing frequencies

                    unique_complex_compliances = \
                        compliance_dict_helper(unique_frequencies, complex_compliance_func,
                                               (compliance, viscosity), rheology_input)

                    # Calculate tidal dissipation
                    if use_planetary_params_for_tide_calc:
                        _radius = object_radius
                        _gravity = object_gravity
                        _density = object_density_bulk
                    else:
                        _radius = viscoelastic_radius_upper
                        _gravity = viscoelastic_gravity
                        _density = viscoelastic_mass / viscoelastic_volume
                    tidal_heating, dUdM, dUdw, dUdO, love_number_by_order_l, negative_imk_by_order_l, \
                        effective_q_by_order_l = \
                            collapse_modes(_gravity, _radius, _density, shear_modulus, layer_tidal_scale,
                                           tidal_host_mass, tidal_susceptibility,
                                           unique_complex_compliances,
                                           tidal_results_by_frequency,
                                           max_tidal_order_l, cpl_ctl_method=False)

            else:
                # No viscoelastic layer present. No tides or convection
                viscosity = 0.
                shear_modulus = 0.
                melt_fraction = 0.
                compliance = 0.
                tidal_heating = 0.
                dUdM, dUdw, dUdO = 0., 0., 0.
                visco_top_temp = bottom_temperatures[layer_i]

                viscoelastic_temperature_delta = 0.
                viscoelastic_cooling = 0.
                viscoelastic_cooling_flux = 0.
                viscoelastic_passes_all_heat = True

            # Add this layer's tidal dissipation to the total for the object
            dUdM_total += dUdM
            dUdw_total += dUdw
            dUdO_total += dUdO

            # Calculate radiogenic heating
            if use_radiogenic_flags[object_i][layer_i]:
                radiogenic_mass = layer_masses[object_i][layer_i]
                radiogenic_time = sec2myr(time)
                radiogenic_input = radiogenic_model_inputs[object_i][layer_i]
                radiogenic_heating = \
                    radiogenic_funcs[object_i][layer_i](radiogenic_time, radiogenic_mass, *radiogenic_input)
            else:
                radiogenic_heating = 0.

            # Determine Heating
            total_incoming_heating = tidal_heating + radiogenic_heating
            print(layer_i, tidal_heating/1e12)
            if not bottom_layer:
                total_incoming_heating += layer_coolings[layer_i - 1]

            # Determine thermal evolution model: growing layer or static
            if is_growth_layer:

                # Determine heating fluxes
                viscoelastic_heating_flux = total_incoming_heating / viscoelastic_surf_area

                if viscoelastic_passes_all_heat:
                    viscoelastic_cooling_flux = viscoelastic_heating_flux

                elastic_heating_flux = viscoelastic_cooling_flux

                if elastic_passes_all_heat:
                    elastic_cooling_flux = viscoelastic_cooling_flux
                    elastic_cooling = viscoelastic_cooling_flux * elastic_surf_area
                    elastic_boundary_layer_thickness = 0.
                    elastic_rayleigh = 0.
                    elastic_nusselt = 1.
                else:
                    # Determine Cooling for the elastic layer
                    if top_layer:
                        elastic_delta_temperature = visco_top_temp - surface_temperatures[object_i]
                    else:
                        elastic_delta_temperature = visco_top_temp - temperatures[object_i][layer_i + 1]
                    #     conduction function
                    elastic_cooling_flux, elastic_boundary_layer_thickness, elastic_rayleigh, elastic_nusselt = \
                        conduction(elastic_delta_temperature, thermal_conductivity, elastic_thickness)
                    elastic_cooling = elastic_surf_area * elastic_cooling_flux

                if layer_stagnant:
                    # If the layer is stagnant (all elastic) then there is no change in temperature.
                    elastic_layer_change = 0.
                    viscoelastic_layer_change = 0.
                    viscoelastic_temperature_change = 0.
                else:
                    energy_to_freeze = specific_heat * viscoelastic_temperature_delta
                    elastic_layer_change = \
                        (elastic_cooling_flux - elastic_heating_flux) / (material_density * energy_to_freeze)
                    if layer_freeze_out:
                        # The viscoelastic layer's temperature is allowed to decrease if the layer is frozen out.
                        viscoelastic_layer_change = -elastic_layer_change
                        viscoelastic_temperature_change = \
                            (total_incoming_heating - viscoelastic_cooling) / \
                            (viscoelastic_volume * material_density * specific_heat)
                    else:
                        viscoelastic_temperature_change = 0.
                        # If there is an ocean then the viscoelastic layer will change depending on the
                        #     energy balance
                        viscoelastic_layer_change = \
                            (viscoelastic_cooling_flux - viscoelastic_heating_flux) / \
                            (material_density * latent_heat)

                # Heat coming out of the layer is equal to the elastic cooling
                layer_cooling = elastic_cooling

            else:
                # Static layer has no viscoelastic/elastic layer growth. Only the average temperature is tracked.
                viscoelastic_temperature_change = \
                    (total_incoming_heating - viscoelastic_cooling) / \
                    (viscoelastic_volume * material_density * specific_heat)

                elastic_layer_change = 0.
                viscoelastic_layer_change = 0.
                layer_cooling = viscoelastic_cooling

            # Store results and derivatives
            layer_coolings.append(layer_cooling)
            derivative_storage.append(viscoelastic_temperature_change)
            derivative_storage.append(elastic_layer_change)
            derivative_storage.append(viscoelastic_layer_change)


        # Determine change in this object's spin-rate
        if spin_locked or dUdO_total == 0.:
            spin_rate_change = 0.
        else:
            spin_rate_change = spin_rate_derivative(dUdO_total, object_mois[object_i], tidal_host_mass)
        derivative_storage.append(spin_rate_change)

        # Store tidal partial derivatives for orbital calculations
        tidal_derivative_storage.append((dUdM_total, dUdw_total, dUdO_total))

    # Determine orbital changes
    obj1_dUdM, obj1_dUdw, obj1_dUdO = tidal_derivative_storage[0]
    obj2_dUdM, obj2_dUdw, obj2_dUdO = tidal_derivative_storage[1]
    if obj1_dUdM == 0. and obj1_dUdw == 0. and obj2_dUdM == 0. and obj2_dUdw == 0.:
        # No tides. No change.
        eccentricity_change = 0.
        orbital_motion_change = 0.
    else:
        da_dt, de_dt = \
            semia_eccen_derivatives_dual(semi_major_axis, orbital_motion, eccentricity,
                                         object_masses[0], obj1_dUdM, obj1_dUdw,
                                         object_masses[1], obj2_dUdM, obj2_dUdw)
        eccentricity_change = de_dt
        orbital_motion_change = (-3. / 2.) * (orbital_motion / semi_major_axis) * da_dt

    if eccentricity < 0.:
        eccentricity_change = 0.

    derivative_storage.append(orbital_motion_change)
    derivative_storage.append(eccentricity_change)

    return derivative_storage


def diffeq_julia(variables, parameters, time):

    output = diffeq_scipy(time, variables)

    return output

if use_julia:
    diffeq = diffeq_julia
else:
    diffeq = diffeq_scipy