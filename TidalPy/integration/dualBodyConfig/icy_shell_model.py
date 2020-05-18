import os

import numpy as np
import matplotlib.pyplot as plt

from ...constants import G
from ...tools.conversions import orbital_motion2semi_a, sec2myr
from ...utilities.performance import njit
from ...tides.mode_manipulation import build_mode_manipulators
from ...tides.dissipation import calc_tidal_susceptibility
from ...cooling.cooling_models import conduction, convection

plt.rcParams.update({'font.size': 14})

MIN_INTERVAL_SCALE = 0.005
MAX_DATA_SIZE = 2000

def build_2layer_icy_shell_diffeq(obj0_config: dict, obj1_config: dict, orbital_config: dict, integration_config: dict):

    # Load planetary data
    object_configs = (obj0_config, obj1_config)
    object_names = tuple([object_config['name'] for object_config in object_configs])
    object_masses = tuple([object_config['mass'] for object_config in object_configs])
    object_radii = tuple([object_config['radius'] for object_config in object_configs])
    object_obliquities = tuple([object_config['constant_obliquity'] for object_config in object_configs])
    surface_temperatures = tuple([object_config['surface_temperature'] for object_config in object_configs])
    num_layers = tuple([len(object_config['layers']) for object_config in object_configs])

    # Layer data are stored as the following: ( (obj0_layer_0, obj0_layer_1), (obj1_layer_0, obj1_layer1) )
    layer_names = tuple(
        [tuple([layer_config['name'] for layer_config in object_config['layers']])
        for object_config in object_configs]
    )
    layer_masses = tuple(
        [tuple([layer_config['mass'] for layer_config in object_config['layers']])
        for object_config in object_configs]
    )
    radii_upper = tuple(
        [tuple([layer_config['radius_upper'] for layer_config in object_config['layers']])
        for object_config in object_configs]
    )
    radii_lower = tuple(
        [tuple([layer_config['radius_lower'] for layer_config in object_config['layers']])
        for object_config in object_configs]
    )
    thermal_conductivities = tuple(
        [tuple([layer_config['thermal_conductivity'] for layer_config in object_config['layers']])
        for object_config in object_configs]
    )
    thermal_expansions = tuple(
        [tuple([layer_config['thermal_expansion'] for layer_config in object_config['layers']])
        for object_config in object_configs]
    )
    specific_heats = tuple(
        [tuple([layer_config['specific_heat'] for layer_config in object_config['layers']])
        for object_config in object_configs]
    )
    convection_betas = tuple(
        [tuple([layer_config['convection_beta'] for layer_config in object_config['layers']])
        for object_config in object_configs]
    )
    convection_alphas = tuple(
        [tuple([layer_config['convection_alpha'] for layer_config in object_config['layers']])
        for object_config in object_configs]
    )
    critical_rayleighs = tuple(
        [tuple([layer_config['critical_rayleigh'] for layer_config in object_config['layers']])
         for object_config in object_configs]
    )
    static_shears = tuple(
        [tuple([layer_config['static_shear'] for layer_config in object_config['layers']])
        for object_config in object_configs]
    )
    viscosity_model_names = tuple(
        [tuple([layer_config['viscosity_model'] for layer_config in object_config['layers']])
        for object_config in object_configs]
    )
    viscosity_inputs = tuple(
        [tuple([layer_config['viscosity_input'] for layer_config in object_config['layers']])
        for object_config in object_configs]
    )
    partial_melt_model_names = tuple(
        [tuple([layer_config['partial_melt_model'] for layer_config in object_config['layers']])
        for object_config in object_configs]
    )
    partial_melt_inputs = tuple(
        [tuple([layer_config['partial_melt_input'] for layer_config in object_config['layers']])
        for object_config in object_configs]
    )
    rheology_model_names = tuple(
        [tuple([layer_config['rheology_model'] for layer_config in object_config['layers']])
        for object_config in object_configs]
    )
    rheology_inputs = tuple(
        [tuple([layer_config['rheology_input'] for layer_config in object_config['layers']])
        for object_config in object_configs]
    )
    tidal_scales = tuple(
        [tuple([layer_config.getfloat('tidal_scale', 1.) for layer_config in object_config['layers']])
        for object_config in object_configs]
    )
    growth_layer_flags = tuple(
        [tuple([layer_config['growth_layer'] for layer_config in object_config['layers']])
         for object_config in object_configs]
    )
    constant_layer_temperatures = tuple(
        [tuple([layer_config.getfloat('constant_layer_temperature', 100.) for layer_config in object_config['layers']])
         for object_config in object_configs]
    )
    viscoelastic_top_temperatures = tuple(
        [tuple([layer_config.getfloat('viscoelastic_top_temperature', 100.) for layer_config in object_config['layers']])
         for object_config in object_configs]
    )
    material_densities = tuple(
        [tuple([layer_config['material_densities'] for layer_config in object_config['layers']])
         for object_config in object_configs]
    )

    # Load orbital configurations
    eccentricity_truncation = orbital_config['eccentricity_truncation']
    max_tidal_order_l = orbital_config['max_tidal_order_l']
    use_obliquity = orbital_config['use_obliquity']
    #    Build tidal mode manipulators
    calculate_tidal_terms, collapse_modes = \
        build_mode_manipulators(max_tidal_order_l, eccentricity_truncation, use_obliquity=use_obliquity)

    # Load integration configurations
    use_planetary_params_for_tide_calc = integration_config['use_planetary_params_for_tides']
    use_tidal_scale = integration_config['use_tidal_scale']
    use_visco_volume_for_tidal_scale = integration_config['use_visco_volume_for_tidal_scale']
    use_julia = integration_config['use_julia']
    time_span = integration_config['time_span']
    time_interval = time_span[-1] - time_span[0]

    # Calculate derived properties
    object_volumes = tuple([(4. / 3.) * np.pi * _radius**3 for _radius in object_radii])
    object_densities_bulk = tuple([_mass / _volume for _mass, _volume in zip(object_masses, object_volumes)])
    object_gravities = tuple([G * _mass / _radius**2 for _mass, _radius in zip(object_masses, object_radii)])
    layer_volumes = tuple(
        [tuple([(4. / 3.) * np.pi * (_radius_upper**3 - _radius_lower**3) for _radius_upper, _radius_lower
                in zip(_obj_radius_upper, _obj_radius_lower)])
         for _obj_radius_upper, _obj_radius_lower in zip(radii_upper, radii_lower)]
    )
    layer_thicknesses = tuple(
        [tuple([_radius_upper - _radius_lower for _radius_upper, _radius_lower
                in zip(_obj_radii_upper, _obj_radii_lower)])
         for _obj_radii_upper, _obj_radii_lower in zip(radii_upper, radii_lower)]
    )
    surf_areas = tuple(
        [tuple([4. * np.pi * _radius_upper**2 for _radius_upper in _obj_radius_upper])
         for _obj_radius_upper in radii_upper]
    )
    layer_masses_below = list()
    for object_i in range(2):
        n_layer = num_layers[object_i]
        _mass = 0.
        _layer_masses_below_for_object = list()
        for layer_i in range(n_layer):
            _mass += layer_masses[object_i][layer_i]
            _layer_masses_below_for_object.append(_mass)
        layer_masses_below.append(tuple(_layer_masses_below_for_object))
    layer_masses_below = tuple(layer_masses_below)
    layer_densities_bulk = tuple(
        [tuple([_mass / _volume for _mass, _volume in zip(_obj_masses, _obj_volumes)])
         for _obj_masses, _obj_volumes in zip(layer_masses, layer_volumes)]
    )
    layer_gravities = tuple(
        [tuple([G * _mass_below / _radius_upper**2 for _mass_below, _radius_upper
                in zip(_obj_masses_below, _obj_upper_radii)])
         for _obj_masses_below, _obj_upper_radii in zip(layer_masses_below, radii_upper)]
    )
    layer_betas = tuple(
        [tuple([_radius * _gravity * _density_bulk for _radius, _gravity, _density_bulk
                in zip(_obj_radii, _obj_gravities, _obj_densities)])
         for _obj_radii, _obj_gravities, _obj_densities in zip(radii_upper, layer_gravities, layer_densities_bulk)]
    )
    thermal_diffusivities = tuple(
        [tuple([_k / (_rho * _c_p) for _k, _rho, _c_p in zip(_ks, _rhos, _c_ps)])
         for _ks, _rhos, _c_ps in zip(thermal_conductivities, material_densities, specific_heats)]
    )
    del _ks, _rhos, _c_ps, _k, _rho, _c_p, _radius, _gravity, _density_bulk, _obj_radii, _obj_gravities, \
        _obj_densities, _obj_masses_below, _obj_upper_radii, _obj_masses_below, _obj_upper_radii, _obj_masses, \
        _obj_volumes, _obj_masses, _obj_volumes, _mass, _volume

    @njit
    def diffeq_scipy(time, variables, parameters):

        # Progress bar
        percent_done = time / time_interval
        print('\rPercent Done: {:0>5.2f}%'.format(100. * percent_done), flush=True, end='')

        # Pull out independent variables
        temperature = list()
        elastic_dx = list()
        viscoelastic_dx = list()
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
            temperature.append(temperatures_by_layer)
            elastic_dx.append(elastic_dx_by_layer)
            viscoelastic_dx.append(viscoelastic_dx_by_layer)
        spin_rate = (variables[index], variables[index + 1])
        orbital_motion = variables[index + 2]
        eccentricity = variables[index + 3]

        # Calculate parameters that only depend on orbital properties
        semi_major_axis = orbital_motion2semi_a(orbital_motion, object_masses[0], object_masses[1])

        # Derivative storage
        derivative_storage = list()

        # Perform thermal evolution on each object's layers
        for object_i in range(2):
            n_layer = num_layers[object_i]

            # Need opposite object index for tidal calculations
            opposite_object_i = 1
            if object_i == 1:
                opposite_object_i = 0

            # Pull out parameters referenced often
            obj_radius = object_radii[object_i]
            tidal_host_mass = object_masses[opposite_object_i]

            # # Calculate tidal modes and susceptibility
            # unique_frequencies, tidal_results_by_frequency = \
            #     calculate_tidal_terms(orbital_motion, spin_rate[object_i], eccentricity, object_obliquity[object_i],
            #                           semi_major_axis, obj_radius)
            # tidal_susceptibility = \
            #     calc_tidal_susceptibility(tidal_host_mass, obj_radius, semi_major_axis)

            # Storage for parameters used by other layers or for later collection
            viscosity_by_layer = list()
            shear_by_layer = list()
            heat_flux_by_layer = list()

            for layer_i in range(n_layer):
                bottom_layer = layer_i == 0
                top_layer = layer_i == num_layers[object_i] - 1

                # Pull out often used parameters
                #    State Properties
                layer_temp = temperature[object_i][layer_i]
                layer_radius_upper = radii_upper[object_i][layer_i]
                layer_radius_lower = radii_lower[object_i][layer_i]
                #    Material Properties
                material_density = material_densities[object_i][layer_i]
                thermal_conductivity = thermal_conductivities[object_i][layer_i],
                thermal_expansion = thermal_expansion[object_i][layer_i]
                thermal_diffusivity = thermal_diffusivities[object_i][layer_i]
                specific_heat = specific_heats[object_i][layer_i]

                # Determine layer geometry
                layer_stagnant = False
                layer_freeze_out = False
                is_growth_layer = growth_layer_flags[object_i][layer_i]
                if is_growth_layer:
                    elastic_radius_upper = layer_radius_upper
                    elastic_radius_lower = layer_radius_upper - elastic_dx[object_i][layer_i]
                    if elastic_radius_lower < layer_radius_lower + 1.:
                        elastic_radius_lower = layer_radius_lower
                        # Layer is all elastic. No viscoelastic portion
                        layer_stagnant = True
                    elif elastic_radius_lower > layer_radius_upper - 0.5:
                        elastic_radius_lower = 0.5

                    if layer_stagnant:
                        viscoelastic_radius_lower = layer_radius_lower
                        viscoelastic_radius_upper = layer_radius_lower
                        # If no viscoelastic layer, then viscoelastic temperature will not matter
                        viscoelastic_temperature = 0.
                    else:
                        viscoelastic_radius_upper = elastic_radius_lower
                        viscoelastic_radius_lower = elastic_radius_lower - viscoelastic_dx[object_i][layer_i]
                        # If the ocean layer is present then the viscoelastic temperature is a constant
                        viscoelastic_temperature = constant_layer_temperatures[object_i][layer_i]
                        if viscoelastic_radius_lower < layer_radius_lower:
                            viscoelastic_radius_lower = layer_radius_lower
                            # Layer is all either elastic or viscoelastic. No ocean
                            layer_freeze_out = True
                            # If the ocean layer is gone, then the viscoelastic temperature will start to decrease
                            viscoelastic_temperature = layer_temp
                    ocean_radius_lower = layer_radius_lower
                    if layer_freeze_out or layer_stagnant:
                        ocean_radius_upper = layer_radius_upper
                    else:
                        ocean_radius_upper = viscoelastic_radius_lower
                    elastic_thickness = elastic_radius_upper - elastic_radius_lower
                    viscoelastic_thickness = viscoelastic_radius_upper - viscoelastic_radius_lower
                    elastic_volume = (4. / 3.) * np.pi * \
                        (elastic_radius_upper * elastic_radius_upper * elastic_radius_upper -
                         elastic_radius_lower * elastic_radius_lower * elastic_radius_lower)
                    viscoelastic_volume = (4. / 3.) * np.pi * \
                        (viscoelastic_radius_upper * viscoelastic_radius_upper * viscoelastic_radius_upper -
                         viscoelastic_radius_lower * viscoelastic_radius_lower * viscoelastic_radius_lower)
                    elastic_surf_area = 4. * np.pi * (elastic_radius_upper * elastic_radius_upper)
                    viscoelastic_surf_area = 4. * np.pi * (viscoelastic_radius_upper * viscoelastic_radius_upper)
                    elastic_mass = elastic_volume * material_density
                    viscoelastic_mass = layer_masses_below[object_i][layer_i] - elastic_mass
                    viscoelastic_gravity = G * viscoelastic_mass / \
                                           (viscoelastic_radius_upper * viscoelastic_radius_upper)

                else:
                    viscoelastic_temperature = 0.
                    viscoelastic_gravity = 0.
                    viscoelastic_thickness = layer_thicknesses[object_i][layer_i]
                    viscoelastic_volume = layer_volumes[object_i][layer_i]
                    elastic_thickness = 0.
                    elastic_volume = 0.

                if use_tidal_scale:
                    if use_visco_volume_for_tidal_scale:
                        layer_tidal_scale = viscoelastic_volume / object_volumes[object_i]
                    else:
                        layer_tidal_scale = tidal_scales[object_i][layer_i]
                else:
                    layer_tidal_scale = 1.


                #
                # # Determine layer's strength (0. pressure dependence for now)
                # viscosity = viscosity_func[object_i][layer_i](
                #     layer_temp, 0., *viscosity_input[object_i][layer_i]
                # )
                # shear_modulus = static_shear[object_i][layer_i]
                # #    Apply partial melting
                # melt_fraction = \
                #     (layer_temp - solidus_temperature[object_i][layer_i]) / \
                #     (liquidus_temperature[object_i][layer_i] - solidus_temperature[object_i][layer_i])
                # if melt_fraction < 0.:
                #     melt_fraction = 0.
                # elif melt_fraction > 1.:
                #     melt_fraction = 1.
                # viscosity, shear_modulus = partial_melt_func[object_i][layer_i](
                #     melt_fraction, viscosity, shear_modulus, *partial_melt_input[object_i][layer_i]
                # )
                # #    Check for over/undershoots
                # if viscosity < 0.1:
                #     viscosity = 0.1
                # elif viscosity > 1.e100:
                #     viscosity = 1.e100
                # if shear_modulus < 0.1:
                #     shear_modulus = 0.1
                # elif shear_modulus > 1.e100:
                #     shear_modulus = 1.e100
                # #    Calculate complex compliance based on layer's strength and the unique tidal forcing frequencies
                # compliance = 1. / shear_modulus
                # unique_complex_compliances = \
                #     {freq_sig: complex_compliance_func(freq, compliance, viscosity, *rheology_input[object_i][layer_i])
                #      for freq_sig, freq in unique_frequencies.items()}
                #
                # # Calculate tidal dissipation
                # if use_planetary_params_for_tide_calc:
                #     _radius = object_radius[object_i]
                #     _gravity = object_gravity[object_i]
                #     _density = object_density_bulk[object_i]
                # else:
                #     _radius = radius_upper[object_i][layer_i]
                #     _gravity = layer_gravity[object_i][layer_i]
                #     _density = layer_density_bulk[object_i][layer_i]
                # tidal_heating, dUdM, dUdw, dUdO, love_number, negative_imk = \
                #     collapse_modes(_gravity, _radius, _density, shear_modulus, unique_complex_compliances,
                #                    tidal_results_by_frequency, tidal_susceptibility, tidal_host_mass, tidal_scale=,
                #                    cpl_ctl_method=False)


                tidal_heating = 0.
                radiogenic_heating = 0.

                # Determine Cooling
                if is_growth_layer:

                    viscoelastic_temperature_delta = viscoelastic_temperature - \
                                                     viscoelastic_top_temperatures[object_i][layer_i]

                    visco_cooling_flux, visco_boundary_layer_thickness, visco_rayleigh, visco_nusselt = \
                        convection(viscoelastic_temperature_delta, viscosity, thermal_conductivities[object_i][layer_i],
                                   thermal_diffusivity[object_i][layer_i], thermal_expansion[object_i][layer_i],
                                   viscoelastic_thickness, viscoelastic_gravity, material_density, convection_alphas[object_i][layer_i],
                                   convection_betas[object_i][layer_i], critical_rayleighs[object_i][layer_i])

                    viscosity: float, thermal_conductivity: float, thermal_diffusivity: float, thermal_expansion: float,
                    layer_thickness: float, gravity: float, density: float,
                    convection_alpha: float, convection_beta: float, critical_rayleigh: float

                conduction, convection



                if top_layer:
                    delta_temperature = layer_temp - surface_temperatures[object_i]
                else:
                    delta_temperature = layer_temp - temperature[object_i][layer_i + 1]

                delta_temperature = layer_temp - surf_temp
                pluto_delta_temp_core = pluto_temperature_core - pluto_temperature_visco
                pluto_q_core_conv, pluto_blt_core, pluto_rayleigh_core, pluto_nusselt_core = \
                    convection_j(pluto_delta_temp_core, pluto_core_thickness, thermal_cond_rock, pluto_viscosity_core,
                                 pluto_thermal_diff_rock, thermal_expansion_rock, pluto_core_surf_gravity,
                                 pluto_core_density,
                                 convection_alpha_rock, convection_beta_rock, critical_rayleigh_rock)
                pluto_core_cooling = pluto_q_core_conv * 4. * np.pi * pluto_core_radius**2

                # Determine Temperature Derivative
                if is_growth_layer:
                    # TODO: add freezeout condition
                    diff_temperature = 0.

                    diff_elastic_dx = (elastic_cooling - elastic_heaitng) / \
                                      ((4. * np.pi * elastic_radius_upper**2) * material_densities[object_i][layer_i] * )

                    pluto_change_visco_ice = (pluto_visco_cooling - pluto_visco_heating) / \
                                             (4. * np.pi * pluto_radius_visco**2 * density_ice * latent_heat_ice)
                    pluto_energy_to_freeze = specific_heat_ice * (pluto_temperature_visco - ice_top_temp)
                    pluto_change_elast_ice = (pluto_elast_cooling - pluto_elast_heating) / \
                                             (4. * np.pi * pluto_radius_elast**2 * density_ice * pluto_energy_to_freeze)
                    if pluto_ice_freezeout:
                        pluto_change_elast_ice = 0.
                        pluto_change_visco_ice = 0.
                    elif pluto_ocean_freezeout:
                        pluto_change_visco_ice = -pluto_change_elast_ice
                    pluto_surface_heat_flux = pluto_elast_cooling / (4. * np.pi * pluto_radius_elast**2)


                else:
                    (pluto_core_heating - pluto_core_cooling) / \
                    (pluto_core_volume * pluto_core_density * specific_heat_rock)

                    pluto_visco_heating = pluto_core_cooling + pluto_tidal_heating
                    pluto_elast_heating = pluto_visco_cooling
                    pluto_elast_cooling = (4. * np.pi * pluto_radius_elast**2) * \
                                          thermal_cond_ice * (ice_top_temp - pluto_surface_temp) / pluto_thickness_elast
                    # :Target
                    charon_visco_heating = charon_core_cooling + charon_tidal_heating
                    charon_elast_heating = charon_visco_cooling
                    charon_elast_cooling = (4. * np.pi * charon_radius_elast**2) * \
                                           thermal_cond_ice * (
                                                       ice_top_temp - charon_surface_temp) / charon_thickness_elast

                    heating_into_layer = 0.
                    if not bottom_layer:
                        heat_flow_bottom = heat_flux_by_layer[layer_i - 1]
                        heating_into_layer = heat_flow_bottom * surf_areas[object_i][layer_i]

                # Store Derivatives
                derivative_storage.append(diff_temperature)
                derivative_storage.append(diff_elastic_dx)
                derivative_storage.append(diff_viscoelastic_dx)



    @njit
    def diffeq_julia(variables, parameters, time):

        output = diffeq_scipy(time, variables, parameters)

        return output

    if use_julia:
        diffeq = diffeq_julia
    else:
        diffeq = diffeq_scipy


    def plotter(variables, time_domain, logtime: bool = False, save_locale: str = None):

        # Plot Styles
        object_colors = ('r', 'b')

        # Saving information
        run_save_name = object_names[0].lower() + '_' + object_names[1].lower()
        if save_locale is None:
            save_locale = os.getcwd()
        else:
            if not os.path.exists(save_locale):
                os.makedirs(save_locale)

        # Setup x-axis
        if logtime:
            x_label = 'Time [yr]'
            x = time_domain / 3.154e7
        else:
            x_label = 'Time [Myr]'
            x = time_domain / 3.154e13

        # Figure 1: Orbital Motion & Spin Plot
        orbspin_fig, orbspin_axes = plt.subplots(ncols=2, figsize=(6.4*1.75, 4.8), tight_layout=True)
        ax_semia = orbspin_axes[0]
        ax_spin = orbspin_axes[1]
        for ax in [ax_semia, ax_spin]:
            ax.set(xscale=logtime, xlabel=x_label)
        ax_eccen = ax_semia.twinx()
        ax_semia.set(ylabel='Semi-major Axis [km]')
        ax_eccen.set(ylabel='Eccentricity (dotted)')
        ax_spin.set(ylabel='Spin Period [days]')

        # Plot non-object dependent parameters
        semi_major_axis = variables['semi_major_axis'] / 1000.
        eccentricity = variables['eccentricity']
        ax_semia.plot(x, semi_major_axis, ls='-')
        ax_eccen.plot(x, eccentricity, ls=':')

        # Plot object dependent parameters
        planet_figures = list()
        for object_i, object_name in enumerate(object_names):
            object_variables = variables[object_name.lower()]
            color = object_colors[object_i]

            # Plot spin period
            spin_period = (2. * np.pi / (object_variables['spin_rate'])) / 86400.
            ax_spin.plot(x, spin_period, label=object_name, ls='-', c=color)

            # Plot layer properties
            n_layers = num_layers[object_i]
            width_factor = (1.75 / 2.) * n_layers
            layer_fig, layer_axes = plt.subplots(ncols=n_layers, figsize=(6.4*width_factor, 4.8), tight_layout=True)
            layer_fig.suptitle(object_name)
            growth_model_encountered = False
            for layer_i, layer_name in enumerate(layer_names[object_i]):
                ax = layer_axes[layer_i]
                ax.title(layer_name)
                layer_variables = object_variables[layer_name.lower()]

                # Determine if this is a growth or regular layer
                use_growth_model = growth_layer_flags[object_i][layer_i]
                if use_growth_model:
                    ax.set(ylabel='Thickness [layer %]')
                    elastic_thickness = layer_variables['elastic_thickness']
                    viscoelastic_thickness = layer_variables['viscoelastic_thickness']
                    ocean_thickness = layer_thicknesses[object_i][layer_i] - \
                                      (elastic_thickness + viscoelastic_thickness)
                    ax.plot(x, elastic_thickness, c=color, ls='-', label='Elastic')
                    ax.plot(x, viscoelastic_thickness, c=color, ls='--', label='Viscoelastic')
                    ax.plot(x, ocean_thickness, c=color, ls=':', label='Ocean')

                    if not growth_model_encountered:
                        ax.legend(loc='best')
                    growth_model_encountered = True
                else:
                    ax.set(ylabel='Temperature [K]')
                    temperature = layer_variables['temperature']
                    ax.plot(x, temperature, c=color, ls='-')

            # Finish up this object's layer plot
            layer_fig.savefig(os.path.join(save_locale, f'{run_save_name}_{object_name}_LayerPlot.pdf'))
            planet_figures.append(layer_fig)

        # Finish up Orbital & Spin Figure
        ax_spin.legend(loc='best')
        orbspin_fig.savefig(os.path.join(save_locale, f'{run_save_name}_OrbitalSpin.pdf'))

        plt.show()

        return orbspin_fig, planet_figures


    def integrator(initial_conditions, diffeq=diffeq,
                   integration_rtol: float = 1.e-8,
                   auto_plot: bool = True, save_data: bool = False,
                   save_locale: str = None, **plotter_kwargs):

        # Determine save location
        run_save_name = object_names[0].lower() + '_' + object_names[1].lower()
        if save_locale is None:
            save_locale = os.getcwd()
        else:
            if not os.path.exists(save_locale):
                os.makedirs(save_locale)

        result_dict = dict()
        print('Integrating Dual Body System:')
        print(f'\t{object_names[0]} and {object_names[1]}')
        if use_julia:
            print('f\tUsing Julia Diffeq...')

            # Import Julia's Diffeqpy and setup the problem
            from diffeqpy import de
            min_interval = MIN_INTERVAL_SCALE * time_interval
            problem = de.ODEProblem(diffeq, initial_conditions, time_span)
            print(f'Solving...')
            solution = de.solve(problem, de.BS3(), saveat=min_interval, abstol=1e-8, reltol=integration_rtol)
            print('\nIntegration Done!')

            y = np.transpose(solution.u)
            t = solution.t

            del solution
        else:
            from scipy.integrate import solve_ivp

            solution = solve_ivp(diffeq, time_span, initial_conditions,
                                 method='LSODA', vectorized=False, rtol=integration_rtol)

            # Pull out dependent variables
            if len(solution.t) > MAX_DATA_SIZE:
                print('Solution data size is very large. Reducing to avoid memory errors in auxiliary grab.')
                t = np.linspace(solution.t[0], solution.t[-1], MAX_DATA_SIZE)
                y = np.asarray([np.interp(t, solution.t, solution.y[i, :])
                                         for i in range(solution.y.shape[0])])
            else:
                t = solution.t
                y = solution.y
            del solution

        # Pull out independent variables
        starting_index = 0
        for object_i, object_name in enumerate(object_names):
            result_dict[object_name] = dict()
            for layer_i, layer_name in enumerate(layer_names):
                result_dict[object_name][layer_name] = dict()
                result_dict[object_name][layer_name]['temperature'] = \
                    y[starting_index + 0, :]
                result_dict[object_name][layer_name]['viscoelastic_thickness'] = \
                    y[starting_index + 1, :]
                result_dict[object_name][layer_name]['elastic_thickness'] = \
                    y[starting_index + 2, :]
                starting_index += 3
        for object_i, object_name in enumerate(object_names):
            result_dict[object_name]['spin_rate'] = y[starting_index, :]
            starting_index += 1
        result_dict['orbital_motion'] = y[starting_index, :]
        result_dict['eccentricity'] = y[starting_index + 1, :]
        time_domain = sec2myr(t)

        del y

        if save_data:
            data_save_dir = save_locale
            print(f'Saving Data to:\n\t{data_save_dir}')
            np.save(os.path.join(data_save_dir, f'{run_save_name}_TimeDomainMyr.npy'), time_domain)
            for object_i, object_name in enumerate(object_names):
                object_results = result_dict[object_name]
                np.save(os.path.join(data_save_dir, f'{run_save_name}_{object_name}_TimeDomainMyr.npy'),
                        object_results['spin_rate'])
                for layer_name in layer_names[object_i]:
                    layer_results = result_dict[object_name][layer_name]
                    for result_name, result in layer_results.items():
                        np.save(os.path.join(data_save_dir,
                                             f'{run_save_name}_{object_name}_{layer_name}_{result_name}.npy'),
                                result)
            print('Data Saved.')

        if auto_plot:
            print('Calling Plotter...')
            plotter(result_dict, time_domain, save_locale=save_locale, **plotter_kwargs)
            print('Plotter Finished.')