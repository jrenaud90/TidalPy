import numpy as np
from scipy.integrate import solve_ivp

from TidalPy.utilities.types import float_eps
from TidalPy.constants import G
from TidalPy.utilities.performance import njit
from TidalPy.toolbox.conversions import orbital_motion2semi_a, semi_a2orbital_motion
from TidalPy.cooling.cooling_models import convection, conduction
from TidalPy.tides.mode_manipulation import find_mode_manipulators
from TidalPy.tides import calc_tidal_susceptibility
from TidalPy.toolbox.conversions import sec2myr, myr2sec, days2rads
from TidalPy.rheology.complex_compliance.complex_compliance import compliance_dict_helper
from TidalPy.rheology.complex_compliance.compliance_models import sundberg

# Over/Undershoots
MAX_VISCOSITY = 1.0e30
MIN_VISCOSITY = 1.
MAX_SHEAR = 1.0e30
MIN_SHEAR = 0.1
MIN_ELASTIC_THICKNESS = 5.
MIN_VISCOELASTIC_THICKNESS = 2.5
MAX_ITERATIONS = 10

# Model switches
model_use_obliquity = False
model_eccentricity_truncation = 10
model_max_tidal_order_l = 2
model_use_tidal_scale = True
model_restrict_tidal_scale_to_visco_layer_volume = True
model_use_planetary_params_for_tide_calc = False
model_use_partialmelt_ice = False
model_use_partialmelt_silicate = True
model_use_static_heating = True
model_force_tidal_lock = False
rheo_model = sundberg
rheology_input = (0.2, .02, 0.33333, 1.)

# Planet Properties
host_mass = 5.683e26
host_radius = 5.83e7
host_volume = (4. / 3.) * np.pi * host_radius**3
host_density_bulk = host_mass / host_volume
target_name = 'Enceladus'
target_radius = 252.1e3
target_mass = 1.08e20
target_volume = (4. / 3.) * np.pi * target_radius**3
target_density_bulk = target_mass / target_volume
target_gravity = G * target_mass / target_radius**2
surface_temperature = 72.15
semi_a_roche = host_radius * (2. * (host_density_bulk / target_density_bulk))**(1 / 3)
period_roche = semi_a2orbital_motion(semi_a_roche, host_mass, target_mass)

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
ice_elastic_threshold_temp = ice_visco_tempdrop * crust_visco_liquidus
ice_visco_mid_temp = 0.5 * (crust_visco_liquidus + ice_elastic_threshold_temp)
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
dont_calculate_core_tides = False
dont_calculate_crust_tides = False

# Prep and setup various functions
calculate_tidal_terms, collapse_modes, eccentricity_func, inclination_func = \
    find_mode_manipulators(model_max_tidal_order_l, model_eccentricity_truncation, use_obliquity=model_use_obliquity)

# Viscosity and Shear Functions
# # Ice
# Pre-Melt Parameters
ice_melting_visco = 5.0e13
ice_premelt_slope = 27.
ice_premelt_shear = 3.3e9

# Partial Melt Parameters
ice_melt_critical_mf = 0.9
ice_melt_mf_width = 0.1
ice_melt_visco_critical_temp = \
    ice_melt_critical_mf * (crust_visco_liquidus - crust_visco_solidus) + crust_visco_solidus
ice_liquid_visco = 0.1
ice_liquid_shear = 0.1
ice_premelt_visco = \
    ice_melting_visco * np.exp(ice_premelt_slope * ((crust_visco_liquidus / ice_melt_visco_critical_temp) - 1.))
ice_melt_visco_slope = -np.log(ice_liquid_visco / ice_premelt_visco) * (1. / ice_melt_mf_width)
ice_melt_shear_slope = -np.log(ice_liquid_shear / ice_premelt_shear) * (1. / ice_melt_mf_width)

@njit(cacheable=True)
def calc_viscoshear_ice(temperature):
    # Critical Melt Fraction
    melt_fraction = (temperature - crust_visco_solidus) / (crust_visco_liquidus - crust_visco_solidus)
    if melt_fraction < 0.:
        melt_fraction = 0.
    elif melt_fraction > 1.:
        melt_fraction = 1.

    if model_use_partialmelt_ice:
        if melt_fraction < ice_melt_critical_mf:
            # No partial melt - follow HS2004
            viscosity = \
                ice_melting_visco * np.exp(ice_premelt_slope * ((crust_visco_liquidus / temperature) - 1.))
            shear = ice_premelt_shear
        elif ice_melt_critical_mf <= melt_fraction < ice_melt_critical_mf + ice_melt_mf_width:
            # Partial melt regime
            viscosity = ice_premelt_visco * np.exp(-ice_melt_visco_slope * (melt_fraction - ice_melt_critical_mf))
            shear = ice_premelt_shear * np.exp(-ice_melt_shear_slope * (melt_fraction - ice_melt_critical_mf))
        else:
            # Acts as a liquid even if there are some ice crystals
            viscosity = ice_liquid_visco
            shear = ice_liquid_shear
    else:
        if melt_fraction < 1.:
            # No partial melt - follow HS2004
            viscosity = \
                ice_melting_visco * np.exp(ice_premelt_slope * ((crust_visco_liquidus / temperature) - 1.))
            shear = ice_premelt_shear
        else:
            # Acts as a liquid even if there are some ice crystals
            viscosity = ice_liquid_visco
            shear = ice_liquid_shear

    return viscosity, shear

# # Silicate
# Pre-Melt Parameters
silicate_melting_visco = 1.5e14
silicate_premelt_slope = 27.
silicate_premelt_shear = 50.e9

# Partial Melt Parameters
silicate_melt_critical_mf = 0.5
silicate_melt_mf_width = 0.1
silicate_melt_visco_critical_temp = \
    silicate_melt_critical_mf * (crust_visco_liquidus - crust_visco_solidus) + crust_visco_solidus
silicate_liquid_visco = 1.e3
silicate_liquid_shear = 0.1
silicate_premelt_visco = \
    silicate_melting_visco * np.exp(silicate_premelt_slope * ((crust_visco_liquidus / silicate_melt_visco_critical_temp) - 1.))
silicate_melt_visco_slope = -np.log(silicate_liquid_visco / silicate_premelt_visco) * (1. / silicate_melt_mf_width)
silicate_melt_shear_slope = -np.log(silicate_liquid_shear / silicate_premelt_shear) * (1. / silicate_melt_mf_width)

@njit(cacheable=True)
def calc_viscoshear_silicate(temperature):
    # Critical Melt Fraction
    melt_fraction = (temperature - core_solidus) / (core_liquidus - core_solidus)
    if melt_fraction < 0.:
        melt_fraction = 0.
    elif melt_fraction > 1.:
        melt_fraction = 1.

    if model_use_partialmelt_silicate:
        if melt_fraction < silicate_melt_critical_mf:
            # No partial melt - follow HS2004
            viscosity = \
                silicate_melting_visco * np.exp(silicate_premelt_slope * ((core_liquidus / temperature) - 1.))
            shear = silicate_premelt_shear
        elif silicate_melt_critical_mf <= melt_fraction < silicate_melt_critical_mf + silicate_melt_mf_width:
            # Partial melt regime
            viscosity = silicate_premelt_visco * np.exp(-silicate_melt_visco_slope * (melt_fraction - silicate_melt_critical_mf))
            shear = silicate_premelt_shear * np.exp(-silicate_melt_shear_slope * (melt_fraction - silicate_melt_critical_mf))
        else:
            # Acts as a liquid even if there are some ice crystals
            viscosity = silicate_liquid_visco
            shear = silicate_liquid_shear
    else:
        if melt_fraction < 1.:
            # No partial melt - follow HS2004
            viscosity = \
                silicate_melting_visco * np.exp(silicate_premelt_slope * ((core_liquidus / temperature) - 1.))
            shear = silicate_premelt_shear
        else:
            # Acts as a liquid even if there are some ice crystals
            viscosity = silicate_liquid_visco
            shear = silicate_liquid_shear

    return viscosity, shear

# @njit(cacheable=True)
def diffeq_scipy(time, variables, extra_params):

    # Variables - eccentricity, orbital_motion, rotation_rate, core_temperature, crust_temperature, crust_visco_dx, crust_elastic_dx
    eccentricity = variables[0]
    obliquity = variables[1]
    orbital_motion = variables[2]
    rotation_rate = variables[3]
    core_temperature = variables[4]
    crust_viscoelastic_temperature = variables[5]
    crust_visco_dx = variables[6]
    crust_elastic_dx = variables[7]

    # Reset Flags
    is_ocean_present = False
    is_fully_ocean = False
    is_fully_stagnant_lid = False
    is_stagnant_lid_gone = False
    is_visco_grounded_out = False
    is_visco_layer_gone = False
    is_visco_temp_diffeq_needed = False
    is_core_cooling_loop_required = False
    is_max_iteration_hit = False
    is_ocean_reemergence = False
    is_visco_reemergence = False
    is_roche_hit = False

    # Derivative corrections (for any changes made in this loop that are not directly due to the diffeq)
    dt_correct_eccentricity = 0.
    dt_correct_obliqutiy = 0.
    dt_correct_orb_motion = 0.
    dt_correct_rotation = 0.
    dt_correct_core_temp = 0.
    dt_correct_crust_visco_temp = 0.
    dt_correct_crust_visco_dx = 0.
    dt_correct_crust_elast_dx = 0.

    # Fix unexpected values
    if crust_elastic_dx <= float_eps:
        dt_correct_crust_elast_dx += np.abs(crust_elastic_dx)
        crust_elastic_dx = float_eps
    if crust_visco_dx <= float_eps:
        dt_correct_crust_visco_dx += np.abs(crust_visco_dx)
        crust_visco_dx = float_eps
    if core_temperature < surface_temperature:
        dt_correct_core_temp += surface_temperature - core_temperature
        core_temperature = surface_temperature
    if crust_viscoelastic_temperature < ice_elastic_threshold_temp:
        dt_correct_crust_visco_temp += ice_elastic_threshold_temp - crust_viscoelastic_temperature
        crust_visco_temp_from_dt = ice_elastic_threshold_temp
    elif crust_viscoelastic_temperature > crust_visco_liquidus:
        dt_correct_crust_visco_temp += crust_visco_liquidus - crust_viscoelastic_temperature
        crust_visco_temp_from_dt = crust_visco_liquidus
    if eccentricity > .99:
        dt_correct_eccentricity += .99 - eccentricity
        eccentricity = .99
    elif eccentricity <= float_eps:
        dt_correct_eccentricity += np.abs(eccentricity)
        eccentricity = 0.
    if np.abs(orbital_motion) <= period_roche:
        # Planet is interior to the roche limit... it is gone :(
        dt_correct_orb_motion = period_roche - np.abs(orbital_motion)
        orbital_motion = np.sign(orbital_motion) * period_roche
        is_roche_hit = True
    if model_force_tidal_lock:
        if np.abs(orbital_motion - rotation_rate) > float_eps:
            # Model forced tidal lock, but the planet is not tidally locked. Make a correction
            dt_correct_rotation = orbital_motion - rotation_rate
            rotation_rate = orbital_motion

    # Other Parameters
    static_core_heating_t0 = extra_params[0]
    static_crust_heating_t0 = extra_params[1]
    static_core_heating_delta_time = extra_params[2]
    static_crust_heating_delta_time = extra_params[3]
    static_core_heating_ref = extra_params[4]
    static_crust_heating_ref = extra_params[5]
    show_progress = extra_params[6]
    time_i = extra_params[7]
    time_f = extra_params[8]

    # Progress bar
    if show_progress:
        percent_done = time / (time_f - time_i)
        percent_done = round(percent_done * 1000., 0) / 10.
        percent_done = int(percent_done)
        if percent_done % 2 == 0:
            if percent_done < 100:
                if percent_done < 10:
                    spaces = '  '
                else:
                    spaces = ' '
            else:
                spaces = ''
            if is_roche_hit:
                print_stmt = 'Time Study: ' + spaces + str(percent_done) +  ' % - ROCHE HIT\r'
            else:
                print_stmt = 'Time Study: ' + spaces + str(percent_done) +  ' %            \r'
            print(print_stmt)

    # Derived Parameters
    semi_major_axis = orbital_motion2semi_a(orbital_motion, host_mass, target_mass)

    # If roche limit is hit then we really don't care about the planet any more. So set all derivatives to zero which
    #    should end an adaptive time step integrator quickly.
    if is_roche_hit:
        output_variables = (
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
        )
        return output_variables

    # Find static heating values (if used)
    static_crust_heating = 0.
    static_core_heating = 0.
    if model_use_static_heating:
        if static_core_heating_t0 <= time <= static_core_heating_t0 + static_core_heating_delta_time:
            static_core_heating = static_core_heating_ref
        if static_crust_heating_t0 <= time <= static_crust_heating_t0 + static_crust_heating_delta_time:
            static_crust_heating = static_crust_heating_ref

    # Target Geometry
    _total_crust_thickness_allowed_ocean = crust_total_thickness - MIN_ELASTIC_THICKNESS - MIN_VISCOELASTIC_THICKNESS
    _half_visco_dx = 0.5 * crust_visco_dx

    # First check if ocean is present
    if (crust_elastic_dx + crust_visco_dx) < _total_crust_thickness_allowed_ocean:
        # The icy layers do not add up to the total crust thickness, the ocean is the only remaining portion.
        is_ocean_present = True
    if (crust_elastic_dx + crust_visco_dx) <= MIN_ELASTIC_THICKNESS:
        # There is no solid ice. Ocean goes to the top. This is unrealistic so always provide a thin layer of ice.
        # Also our model assumes that if there is a ocean then there is also a viscoelastic layer at the base of the
        #    elastic ice.
        is_ocean_present = True
        is_fully_ocean = True
        is_visco_layer_gone = False
        is_stagnant_lid_gone = False
        dt_correct_crust_visco_dx = MIN_VISCOELASTIC_THICKNESS - crust_visco_dx
        dt_correct_crust_elast_dx = MIN_ELASTIC_THICKNESS - crust_elastic_dx
        crust_elastic_dx = MIN_ELASTIC_THICKNESS
        crust_visco_dx = MIN_VISCOELASTIC_THICKNESS
    elif (crust_visco_liquidus - crust_viscoelastic_temperature) <= float_eps:
        # Viscoelastic layer is very warm. Give up half of it to the ocean
        is_visco_layer_gone = False
        crust_visco_dx = _half_visco_dx
        if not is_ocean_present:
            # There was no ocean present, but based on this temperature there should be. Assume ocean forms
            #    with half the visco layer's thickness
            is_ocean_reemergence = True
            is_ocean_present = True

        # Correct temperature (jump to mid-temp)
        dt_correct_crust_visco_temp = ice_visco_mid_temp - crust_viscoelastic_temperature
        crust_visco_temperature = ice_visco_mid_temp

    if is_ocean_present:
        # There is at least a thin layer of ocean. Set the temperatures based off the ocean being adiabatic
        core_top_temp = crust_visco_liquidus
        visco_bottom_temp = crust_visco_liquidus
        elastic_bottom_temp = ice_elastic_threshold_temp

    # Now check if the viscoelastic layer is present
    if crust_elastic_dx >= crust_total_thickness:
        # Elastic layer comprises full water layer
        is_fully_stagnant_lid = True
        is_visco_layer_gone = True
        dt_correct_crust_visco_dx = -crust_visco_dx
        crust_visco_dx = 0.

        # We will need to calculate the temperature at the base of this layer later using an interactive method.
        is_core_cooling_loop_required = True
    if crust_visco_dx <= float_eps:
        # No or nearly no viscoelastic layer found.
        # If there is no ocean then assume elastic layer is grounded out.
        if is_ocean_present or is_fully_ocean:
            # If there is an ocean, assume a think viscoelastic layer separates it from the elastic layer
            is_visco_layer_gone = False
            dt_correct_crust_visco_dx = MIN_VISCOELASTIC_THICKNESS - crust_visco_dx
            crust_visco_dx = MIN_VISCOELASTIC_THICKNESS

            # Update temperatures
            core_top_temp = crust_visco_liquidus
            visco_bottom_temp = crust_visco_liquidus
            elastic_bottom_temp = ice_elastic_threshold_temp
        else:
            # If there is no ocean then assume elastic layer is grounded out.
            is_visco_layer_gone = True
            dt_correct_crust_visco_dx = -crust_visco_dx
            crust_visco_dx = 0.

            # We will need to calculate the temperature at the base of this layer later using an interactive method.
            is_core_cooling_loop_required = True
    elif (crust_viscoelastic_temperature - ice_elastic_threshold_temp) <= float_eps:
        # Viscoelastic layer is very cold. Start transitioning it into elastic.
        is_visco_layer_gone = False
        dt_correct_crust_visco_dx = -_half_visco_dx
        crust_visco_dx = _half_visco_dx
        dt_correct_crust_elast_dx = _half_visco_dx
        crust_elastic_dx = crust_elastic_dx + _half_visco_dx

        # Correct temperature
        dt_correct_crust_visco_temp = ice_visco_mid_temp - crust_viscoelastic_temperature
        crust_visco_temperature = ice_visco_mid_temp

    # Check if elastic layer is gone
    if crust_elastic_dx <= float_eps:
        # The elastic layer is gone but we assume that there is always at least a thin layer over the ocean or visco.
        dt_correct_crust_visco_dx = MIN_ELASTIC_THICKNESS - crust_elastic_dx
        crust_elastic_dx = MIN_ELASTIC_THICKNESS
        elastic_bottom_temp = ice_elastic_threshold_temp


    """
    elif crust_elastic_dx <= MIN_ELASTIC_THICKNESS:
        # Elastic layer is non-existent. Make it small so we don't have divide by zero issues
        is_stagnant_lid_gone = True
        crust_elastic_dx = MIN_ELASTIC_THICKNESS

    if is_fully_stagnant_lid:
        # There is no viscoelastic ice layer. The viscoelastic temperature won't matter so set it to nan.
        is_visco_layer_gone = True
        crust_viscoelastic_temperature = 0.0
        crust_viscoelastic_top_temperature = 0.0

        # The bottom temperature is still required for core cooling.
        # An iterative loop will be required to solve for this bottom temperature of the elastic layer
        is_core_cooling_loop_required = True
        crust_bottom_temperature = 0.0
        crust_elastic_bottom_temperature = 0.0
    else:
        if is_fully_ocean:
            # There is no viscoelastic ice layer due to the large ocean
            is_visco_layer_gone = True
            crust_viscoelastic_temperature = 0.0
            crust_viscoelastic_top_temperature = 0.0
            crust_bottom_temperature = crust_visco_liquidus
            crust_elastic_bottom_temperature = crust_bottom_temperature
        elif is_ocean_present:
            # There is a viscoelastic layer and an underlying ocean.
            # # Determine temperature within viscoelastic ice
            is_visco_layer_gone = True
            crust_bottom_temperature = crust_visco_liquidus
            crust_viscoelastic_top_temperature = crust_bottom_temperature * ice_visco_tempdrop
            crust_elastic_bottom_temperature = crust_viscoelastic_top_temperature
            crust_viscoelastic_temperature = \
                0.5 * (crust_viscoelastic_top_temperature + crust_bottom_temperature)
        else:
            # Grounded Out: There is a viscoelastic layer, but no ocean. Viscoelastic layer is touching the top of
            #    the silicate core.
            is_visco_layer_gone = True
            is_visco_grounded_out = True
            # The temperature of the viscoelastic layer is determined by a differential equation.
            # Assume the bottom of the viscoelastic layer is what the differential equation is tracking
            # TODO: This is not 100% accurate
            crust_bottom_temperature = crust_visco_temp_from_dt
            crust_viscoelastic_top_temperature = crust_bottom_temperature * ice_visco_tempdrop
            crust_elastic_bottom_temperature = crust_viscoelastic_top_temperature
            crust_viscoelastic_temperature = \
                0.5 * (crust_viscoelastic_top_temperature + crust_bottom_temperature)
    """




    # Calculate viscosity and shear modulus of layers
    # # Core
    core_viscosity, core_shear = calc_viscoshear_silicate(core_temperature)
    # # Viscoelastic Ice
    if is_visco_layer_gone:
        crust_viscosity = MIN_VISCOSITY
        crust_shear = MIN_SHEAR
    else:
        crust_viscosity, crust_shear = calc_viscoshear_ice(crust_viscoelastic_temperature)

    # Check for over/undershoots in viscosity and shear
    if crust_viscosity < ice_liquid_visco:
        crust_viscosity = ice_liquid_visco
    elif crust_viscosity > MAX_VISCOSITY:
        crust_viscosity = MAX_VISCOSITY
    if crust_shear < ice_liquid_shear:
        crust_shear = ice_liquid_shear
    elif crust_shear > MAX_SHEAR:
        crust_shear = MAX_SHEAR
    if core_viscosity < silicate_liquid_visco:
        core_viscosity = silicate_liquid_visco
    elif core_viscosity > MAX_VISCOSITY:
        core_viscosity = MAX_VISCOSITY
    if core_shear < silicate_liquid_shear:
        core_shear = silicate_liquid_shear
    elif core_shear > MAX_SHEAR:
        core_shear = MAX_SHEAR

    # Determine icy layer geometry and structure
    crust_elastic_radius = target_radius
    crust_visco_radius = target_radius - crust_elastic_dx
    crust_ocean_radius = crust_visco_radius - crust_visco_dx
    crust_elastic_volume = (4. / 3.) * np.pi * (crust_elastic_radius**3 - crust_visco_radius**3)
    crust_visco_volume = (4. / 3.) * np.pi * (crust_visco_radius**3 - crust_ocean_radius**3)
    crust_ocean_volume = (4. / 3.) * np.pi * (crust_ocean_radius**3 - core_radius**3)
    crust_elastic_sa = 4. * np.pi * target_radius**2
    crust_visco_sa = 4. * np.pi * crust_visco_radius**2
    crust_ocean_sa = 4. * np.pi * crust_ocean_radius**2
    crust_elastic_gravity = target_gravity
    _mass_below_visco = core_mass + (crust_ocean_volume * water_density) + (crust_visco_volume * ice_density)
    crust_visco_gravity = G * _mass_below_visco / crust_visco_radius**2

    # Update tidal scales
    if model_use_tidal_scale:
        core_tidal_scale = core_volume / target_volume
        if is_visco_layer_gone:
            crust_tidal_scale = 0.
        else:
            crust_tidal_scale = crust_visco_volume / target_volume
    else:
        core_tidal_scale = 1.
        crust_tidal_scale = 1.

    # # Calculate Tidal Response and Dissipation
    # Calculate eccentricity and obliquity functions
    eccentricity_results = eccentricity_func(eccentricity)
    obliquity_results = inclination_func(obliquity)
    tidal_susceptibility = \
        calc_tidal_susceptibility(host_mass, target_radius, semi_major_axis)

    # Calculate tidal modes and susceptibility
    unique_frequencies, tidal_results_by_frequency = \
        calculate_tidal_terms(rotation_rate, orbital_motion, semi_major_axis, target_radius,
                              eccentricity_results, obliquity_results)

    dUdM_total = 0.
    dUdw_total = 0.
    dUdO_total = 0.
    core_dUdM, crust_dUdM = 0., 0.
    core_dUdw, crust_dUdw = 0., 0.
    core_dUdO, crust_dUdO = 0., 0.
    tidal_heating_global = 0.
    core_tidal_heating = 0.
    crust_tidal_heating = 0.
    if not dont_calculate_tides:
        # Calculate tides!
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

        # Calculate the rheological response of the core and crust
        # Calculate complex compliance based on layer's strength and the unique tidal forcing frequencies
        if not dont_calculate_core_tides:
            core_unique_complex_compliances = \
                compliance_dict_helper(
                    unique_frequencies, rheo_model,
                    (core_shear**(-1), core_viscosity), rheology_input
                    )

            # Calculate tidal dissipation for the core
            core_tidal_heating, core_dUdM, core_dUdw, core_dUdO, core_love_number_by_order_l, \
                core_negative_imk_by_order_l, core_effective_q_by_order_l = \
                    collapse_modes(
                        core_tidal_gravity, core_tidal_radius, core_tidal_density,
                        core_shear, core_tidal_scale,
                        host_mass, tidal_susceptibility,
                        core_unique_complex_compliances,
                        tidal_results_by_frequency,
                        model_max_tidal_order_l, cpl_ctl_method=False
                        )

        # Now calculate for the viscoelastic lid, if present.
        if not is_visco_layer_gone and not dont_calculate_crust_tides:
            crust_unique_complex_compliances = \
                compliance_dict_helper(
                    unique_frequencies, rheo_model,
                    (crust_shear**(-1), crust_viscosity), rheology_input
                    )

            crust_tidal_heating, crust_dUdM, crust_dUdw, crust_dUdO, crust_love_number_by_order_l, \
                crust_negative_imk_by_order_l, crust_effective_q_by_order_l = \
                    collapse_modes(
                        crust_tidal_gravity, crust_tidal_radius, crust_tidal_density,
                        crust_shear, crust_tidal_scale,
                        host_mass, tidal_susceptibility,
                        crust_unique_complex_compliances,
                        tidal_results_by_frequency,
                        model_max_tidal_order_l, cpl_ctl_method=False
                        )

        # Combine the tidal results linearly (note the tidal_scale scales each layer's values if enabled)
        dUdM_total = core_dUdM + crust_dUdM
        dUdw_total = core_dUdw + crust_dUdw
        dUdO_total = core_dUdO + crust_dUdO
        tidal_heating_global = core_tidal_heating + crust_tidal_heating

    # Update orbital and rotation parameters
    # TODO: For now this is all off
    dt_eccentricity = 0.
    dt_obliquity = 0.
    dt_orbital_motion = 0.
    dt_rotation_rate = 0.

    # Calculate heat balance and cooling (Start from the bottom to the top)
    # Calculate heat sources
    # # Radiogenic
    # TODO: No radiogenic for now
    core_radio_heating = .0
    crust_radio_heating = .0
    # # Static
    core_heating = core_radio_heating + core_tidal_heating + static_core_heating
    crust_visco_invitro_heating = crust_radio_heating + crust_tidal_heating + static_crust_heating

    # # Core
    # Determine the temperature at the top of the core
    core_heat_flux = core_heating / core_sa
    _new_visco_dx = 0.
    _new_ocean_dx = 0.
    if not is_core_cooling_loop_required:
        # There is either ocean or viscoelastic ice over top of the core. The top temperature is known and does not
        #    need to be solved for. We can jump straight to cooling

    else:
        # There is elastic ice on top of the core. We don't know what the temperature is at the core-ice boundary.
        #    we will use an iterative loop to balance the heat leaving the core and that leaving the icy shell
        #    to solve for this temperature.















    if not is_core_cooling_loop_required:
        # An iterative loop will be required to find the temperature at the top of the core.
        # Start the temperature to be equal to the min for elastic ice.
        core_top_temperature = ice_elastic_threshold_temp
        crust_elastic_cooling = 1000.
        core_cooling = 1.
        current_iteration = 0
        while abs(core_cooling - crust_elastic_cooling) / core_cooling > 0.05:
            if current_iteration > MAX_ITERATIONS:
                # Max iterations hit. Go with whatever is calculated now.
                is_max_iteration_hit = True
                # TODO: Any sort of warning here other than the flag?
                break

            # Perform convection calculations
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

            # Update new bottom temperature based on heat flux
            # TODO: Is this right? It assumes equilibrium which may not be reached in a given timestep.
            core_top_temperature = \
                surface_temperature + core_cooling_flux * crust_total_thickness / ice_thermal_conductivity

            # Check for if the temperature gets hot enough for ocean or viscoelastic reemergence.
            if core_top_temperature > ice_elastic_threshold_temp:
                # Viscoelastic reemergence
                is_visco_reemergence = True
            elif core_top_temperature >= crust_visco_liquidus:
                is_visco_reemergence = True
                is_ocean_reemergence = True
                core_top_temperature = crust_visco_liquidus
            else:
                is_visco_reemergence = False
                is_ocean_reemergence = False

            # Calculate elastic cooling
            crust_elastic_delta_temp = core_top_temperature - surface_temperature
            crust_elastic_cooling_flux, crust_elastic_boundary_dx, crust_elastic_rayleigh, crust_elastic_nusselt = \
                conduction(crust_elastic_delta_temp, ice_thermal_conductivity, crust_elastic_dx)
            crust_elastic_cooling = crust_elastic_cooling_flux * crust_elastic_sa

            # Next loop
            current_iteration += 1

    else:
        # Perform convection calculation just once.
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
    if is_visco_layer_gone or crust_viscoelastic_layer_passes_all_heat:

        # Pass all heat through the layer
        crust_visco_cooling_flux = crust_visco_base_heat_flux
        crust_visco_delta_temp = 0.
        crust_visco_boundary_dx = float_eps
        crust_visco_rayleigh = 0.
        crust_visco_nusselt = 1.
        crust_visco_cooling = crust_visco_heating

        if is_visco_reemergence:
            # Viscoelastic layer may have reemerged during this loop.
            # Update elastic bottom temperature
            crust_elastic_bottom_temperature = ice_elastic_threshold_temp

            # Assume at this early stage the layer is purely conductive and thin.
            # The cooling will be a large fraction of the heat flux.
            dt_crust_visco_dx = \
                (0.99 * crust_visco_cooling_flux - crust_visco_base_heat_flux) / (ice_density * ice_latent_heat)
        else:
            dt_crust_visco_dx = 0.

    else:
        # Solve viscoelastic layer cooling using a convection model
        crust_visco_delta_temp = crust_bottom_temperature - crust_viscoelastic_top_temperature
        crust_visco_cooling_flux, crust_visco_boundary_dx, crust_visco_rayleigh, crust_visco_nusselt = \
            convection(
                crust_visco_delta_temp, crust_viscosity,
                ice_thermal_conductivity, ice_thermal_diffusivity, ice_thermal_expansion,
                crust_visco_dx, crust_visco_gravity, ice_density,
                ice_convection_alpha, ice_convection_beta, ice_critical_rayleigh
            )
        crust_visco_cooling = crust_visco_cooling_flux*crust_visco_sa

        dt_crust_visco_dx = \
            (crust_visco_cooling_flux - crust_visco_base_heat_flux)/(ice_density*ice_latent_heat)

    if crust_viscoelastic_layer_passes_all_heat:
        # Pass all heat through the layer
        crust_visco_cooling_flux = crust_visco_base_heat_flux
        crust_visco_delta_temp = 0.
        crust_visco_boundary_dx = float_eps
        crust_visco_rayleigh = 0.
        crust_visco_nusselt = 1.
        crust_visco_cooling = crust_visco_heating


    crust_visco_heating = crust_visco_invitro_heating + core_cooling
    # Grow or shrink ice layers based on temperature balance
    crust_visco_base_heat_flux = crust_visco_heating / crust_visco_sa
    if is_visco_layer_gone or crust_viscoelastic_layer_passes_all_heat:
        # Skip convection calculation if the viscoelastic layer is gone or a debug switch is on
        crust_visco_cooling_flux = crust_visco_base_heat_flux
        crust_visco_delta_temp = 0.
        crust_visco_boundary_dx = float_eps
        crust_visco_rayleigh = 0.
        crust_visco_nusselt = 1.
        crust_visco_cooling = crust_visco_heating
    else:
        # Perform convection calculations
        crust_visco_delta_temp = crust_bottom_temperature - crust_viscoelastic_top_temperature
        crust_visco_cooling_flux, crust_visco_boundary_dx, crust_visco_rayleigh, crust_visco_nusselt = \
            convection(
                crust_visco_delta_temp, crust_viscosity,
                ice_thermal_conductivity, ice_thermal_diffusivity, ice_thermal_expansion,
                crust_visco_dx, crust_visco_gravity, ice_density,
                ice_convection_alpha, ice_convection_beta, ice_critical_rayleigh
            )
        crust_visco_cooling = crust_visco_cooling_flux * crust_visco_sa

    # Heating entering the base of the elastic layer is equal to that exiting the viscoelastic one
    crust_elastic_heating = crust_visco_cooling
    crust_elastic_base_heat_flux = crust_visco_cooling_flux
    if is_stagnant_lid_gone or crust_elastic_layer_passes_all_heat:
        # There is either no elastic layer or all of the heat is leaving it instantly due to debug switch
        crust_elastic_cooling_flux = crust_elastic_heating / crust_elastic_sa
        crust_elastic_delta_temp = 0.
        crust_elastic_boundary_dx = MIN_ELASTIC_THICKNESS
        crust_elastic_rayleigh = 0.
        crust_elastic_nusselt = 1.
        crust_elastic_cooling = crust_elastic_heating
    else:
        # There is an elastic layer, calculate conductive flux
        crust_elastic_delta_temp = crust_elastic_bottom_temperature - surface_temperature
        crust_elastic_cooling_flux, crust_elastic_boundary_dx, crust_elastic_rayleigh, crust_elastic_nusselt = \
                conduction(crust_elastic_delta_temp, ice_thermal_conductivity, crust_elastic_dx)
        crust_elastic_cooling = crust_elastic_cooling_flux * crust_elastic_sa

    # Calculate the change in the elastic layer's thickness
    crust_visco_energy_to_freeze = ice_specific_heat * (crust_visco_liquidus - crust_viscoelastic_top_temperature)
    if crust_elastic_layer_passes_all_heat or is_fully_stagnant_lid:
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

    # Make corrections for the viscoelastic temperature
    if dt_crust_visco_temperature != 0.:
        if crust_visco_temp_from_dt <= surface_temperature:
            dt_crust_visco_temperature = 0.
        elif crust_visco_temp_from_dt >= crust_visco_liquidus:
            dt_crust_visco_temperature = 0.

    # Combine diffeq variables with there corrections and correct any additional problems
    dt_eccentricity += dt_correct_eccentricity
    dt_obliquity += dt_correct_obliqutiy
    dt_orbital_motion += dt_correct_orb_motion
    dt_rotation_rate += dt_correct_rotation
    dt_core_temperature += dt_correct_core_temp
    dt_crust_visco_temperature += dt_correct_crust_visco_temp
    dt_crust_visco_dx += dt_correct_crust_visco_dx
    dt_crust_elastic_dx += dt_correct_crust_elast_dx

    # Store diffeq results in output array
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

def integrate_model(input_dict):

    # Get model parameters
    initial_time = myr2sec(input_dict['initial_time'])
    final_time = myr2sec(input_dict['final_time'])
    rtol = input_dict['integration_rtol']
    atol = input_dict['integration_atol']
    int_method = input_dict['integration_method']

    # Get initial conditions
    i_eccentricity = input_dict['initial_eccentricity']
    i_obliquity = input_dict['initial_obliquity']
    i_orbital_motion = input_dict['initial_orbital_motion']
    i_rotation_rate = input_dict['initial_rotation_rate']
    i_core_temp = input_dict['initial_core_temp']
    i_visco_temp = input_dict['initial_crust_visco_temp']
    i_crust_visco_dx = input_dict['initial_crust_visco_dx']
    i_crust_elastic_dx = input_dict['initial_crust_elastic_dx']

    # Build initial condition array
    if i_visco_temp is None:
        i_visco_temp = 0.5 * (ice_elastic_threshold_temp + crust_visco_liquidus)
    initial_conditions = (
        i_eccentricity,
        i_obliquity,
        i_orbital_motion,
        i_rotation_rate,
        i_core_temp,
        i_visco_temp,
        i_crust_visco_dx,
        i_crust_elastic_dx
    )

    # Build other parameters
    other_args = (
        input_dict['core_heating_t0'],
        input_dict['crust_heating_t0'],
        input_dict['core_heating_delta_time'],
        input_dict['crust_heating_delta_time'],
        input_dict['static_core_heating'],
        input_dict['static_crust_heating'],
        True,
        initial_time,
        final_time
    )
    other_args = (other_args,)

    # Perform integration
    print('Integration Starting...\n\n')
    solution = solve_ivp(diffeq_scipy, (initial_time, final_time), initial_conditions, args=other_args,
                         method=int_method, rtol=rtol, atol=atol)

    if not solution.success:
        print('\nIntegration failed: ', solution.message)
    else:
        print('\nIntegration finished!')
        return solution

def solution_plotter(integration_solution):
    import matplotlib.pyplot as plt

    # Setup figure
    fig, axes = plt.subplots(nrows=2, ncols=2)
    ax_orbit = axes[0, 0]
    ax_spin = axes[0, 1]
    ax_core_temp = axes[1, 0]
    ax_dxs = axes[1, 1]
    ax_eccen = ax_orbit.twinx()
    ax_crust_temp = ax_core_temp.twinx()
    ax_eccen.set(ylabel='Eccentricity (Dotted)')
    ax_orbit.set(ylabel='Semi-major Axis [R$_{S}$]')
    ax_spin.set(ylabel='Spin / Orbital Motion')
    ax_core_temp.set(ylabel='Core Temperature [K]')
    ax_crust_temp.set(ylabel='Ice Temperature (dotted) [K]')
    ax_dxs.set(ylabel='Ice Crust Thicknesses [km]')

    # Pull out results
    time_domain = integration_solution.t
    y = integration_solution.y
    eccentricity = y[0, :]
    obliquity = y[1, :]
    orbital_motion = y[2, :]
    rotation_rate = y[3, :]
    core_temp = y[4, :]
    visco_temp = y[5, :]
    crust_visco_dx = y[6, :]
    crust_elastic_dx = y[7, :]

    # Make any conversions
    semi_a = orbital_motion2semi_a(orbital_motion, host_mass, target_mass)
    semi_a = semi_a / host_radius
    time_domain = sec2myr(time_domain)
    ocean_dx = crust_total_thickness - crust_visco_dx - crust_elastic_dx
    ocean_dx = ocean_dx / 1000.
    crust_visco_dx = crust_visco_dx / 1000.
    crust_elastic_dx = crust_elastic_dx / 1000.
    spin_frac = rotation_rate / orbital_motion

    # Plot results
    ax_orbit.plot(time_domain, semi_a, c='k', ls='-')
    ax_orbit.plot(time_domain, eccentricity, c='k', ls=':')
    ax_spin.plot(time_domain, spin_frac, c='k', ls='-')
    ax_core_temp.plot(time_domain, core_temp, c='k', ls='-')
    ax_crust_temp.plot(time_domain, visco_temp, c='k', ls=':')
    ax_dxs.plot(time_domain, crust_elastic_dx, c='k', ls='-', label='Elastic')
    ax_dxs.plot(time_domain, crust_visco_dx, c='k', ls='--', label='Viscoelastic')
    ax_dxs.plot(time_domain, ocean_dx, c='k', ls=':', label='Ocean')

    # Show and save plot
    fig.tight_layout()
    plt.show()


default_dict = {
    'initial_time' : 0.0,
    'final_time' : 1000.,
    'integration_rtol': 1.0e-4,
    'integration_atol': 1.0e-6,
    'integration_method': 'RK45',
    'initial_eccentricity': 0.0047,
    'initial_obliquity': 0.0,
    'initial_orbital_motion': days2rads(1.370218),
    'initial_rotation_rate': days2rads(1.370218),
    'initial_core_temp': 800.,
    'initial_crust_visco_temp': None,
    'initial_crust_visco_dx': crust_total_thickness * 0.35,
    'initial_crust_elastic_dx': crust_total_thickness * 0.6,
    'core_heating_t0': 0.,
    'crust_heating_t0': 0.,
    'core_heating_delta_time': 50.,
    'crust_heating_delta_time': 50.,
    'static_core_heating': 0.,
    'static_crust_heating': 0.
    }

if __name__ == '__main__':

    sol = integrate_model(default_dict)

    if sol is not None:
        solution_plotter(sol)
