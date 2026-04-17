from typing import List, Tuple

import numpy as np

from TidalPy.RadialSolver import build_rs_input_homogeneous_layers
from TidalPy.utilities.spherical_helper.volume import calculate_voxel_volumes

from structures import Planet


def build_rs_inputs(planet: Planet):

    density_tuple                 = tuple(layer_data.density for layer_data in planet.layers)
    static_bulk_modulus_tuple     = tuple(layer_data.bulk_modulus for layer_data in planet.layers)
    static_shear_modulus_tuple    = tuple(layer_data.shear_modulus for layer_data in planet.layers)
    bulk_viscosity_tuple          = tuple(layer_data.bulk_viscosity for layer_data in planet.layers)
    shear_viscosity_tuple         = tuple(layer_data.shear_viscosity for layer_data in planet.layers)
    layer_type_tuple              = tuple(layer_data.layer_type for layer_data in planet.layers)
    layer_is_static_tuple         = tuple(layer_data.is_static for layer_data in planet.layers)
    layer_is_incompressible_tuple = tuple(layer_data.is_incompressible for layer_data in planet.layers)
    shear_rheology_model_tuple    = tuple(layer_data.shear_rheology_model for layer_data in planet.layers)
    bulk_rheology_model_tuple     = tuple(layer_data.bulk_rheology_model for layer_data in planet.layers)
    radius_fraction_tuple         = tuple(layer_data.radius / planet.radius for layer_data in planet.layers)

    rs_inputs = build_rs_input_homogeneous_layers(
        planet.radius,
        # For now assume the orbital frequency, we can adjust later if needed.
        planet.orbital_freq,
        density_tuple,
        static_bulk_modulus_tuple,
        static_shear_modulus_tuple,
        bulk_viscosity_tuple,
        shear_viscosity_tuple,
        layer_type_tuple,
        layer_is_static_tuple,
        layer_is_incompressible_tuple,
        shear_rheology_model_tuple,
        bulk_rheology_model_tuple,
        radius_fraction_tuple,
        thickness_fraction_tuple = None,
        volume_fraction_tuple = None,
        slices_tuple = None,
        slice_per_layer = 25,
        perform_checks = True)
    
    # Find layer indices
    num_layers = len(density_tuple)
    static_shear_array = np.zeros_like(rs_inputs.radius_array)
    static_bulk_array = np.zeros_like(rs_inputs.radius_array)
    shear_visc_array = np.zeros_like(rs_inputs.radius_array)
    bulk_visc_array = np.zeros_like(rs_inputs.radius_array)
    layer_indices = list()
    last_upper_r = 0.0
    radius_array = rs_inputs.radius_array
    for layer_i in range(num_layers):
        upper_r = rs_inputs.upper_radius_bylayer_array[layer_i]
        layer_index = np.logical_and(radius_array <= upper_r, radius_array >= last_upper_r)
        # The above will over count at layer interfaces. Fix lower one:
        if last_upper_r != 0.0:
            loc = np.where(radius_array==last_upper_r)[0][0]
            layer_index[loc] = False
        # Fix upper one
        if upper_r != planet.radius:
            loc = np.where(radius_array==upper_r)[0][1]
            layer_index[loc] = False
        layer_indices.append(layer_index)

        # Update arrays
        static_shear_array[layer_index] = static_shear_modulus_tuple[layer_i]
        static_bulk_array[layer_index] = static_bulk_modulus_tuple[layer_i]
        shear_visc_array[layer_index] = shear_viscosity_tuple[layer_i]
        bulk_visc_array[layer_index] = bulk_viscosity_tuple[layer_i]

        last_upper_r = upper_r
    
    # Determine obliquity that will be used in heating
    use_obliquity = None
    if planet.obliquity != 0.0:
        use_obliquity = planet.obliquity

    # Build other domains and the final multidimensional arrays
    N = 20
    colatitude = np.radians(np.linspace(0.1, 179.9, N))
    longitude = np.radians(np.linspace(0., 360., N+1))
    time = np.linspace(0., 2. * np.pi / planet.orbital_freq, N+2)
    voxel_volumes = calculate_voxel_volumes(radius_array, longitude, colatitude)
    longitude_matrix, colatitude_matrix, time_matrix = \
        np.meshgrid(longitude, colatitude, time, indexing='ij')

    cmm_inputs = (
        planet.orbital_freq,
        planet.spin_freq,
        planet.semi_major_axis,
        planet.eccentricity,
        planet.host_mass,
        rs_inputs.planet_bulk_density,
        rs_inputs.radius_array,
        rs_inputs.density_array,
        static_bulk_array,
        static_shear_array,
        bulk_visc_array,
        shear_visc_array,
        shear_rheology_model_tuple,
        bulk_rheology_model_tuple,
        rs_inputs.upper_radius_bylayer_array,
        longitude_matrix,
        colatitude_matrix,
        time_matrix,
        voxel_volumes,
        layer_type_tuple,
        layer_is_static_tuple,
        layer_is_incompressible_tuple,
        use_obliquity
    )

    return rs_inputs, cmm_inputs, (colatitude, longitude, time, layer_indices)

rs_optional_args = dict(
    use_kamata = False,
    integration_method = 'DOP853',
    integration_rtol = 1.0e-5,
    integration_atol = 1.0e-8,
    scale_rtols_bylayer_type = False,
    max_num_steps = 500_000,
    expected_size = 1000,
    max_ram_MB = 500,
    max_step = 0,
    # Propagation matrix method parameters
    use_prop_matrix = False,
    core_model = 0,
    # Equation of State solver parameters
    eos_method_bylayer = None,
    surface_pressure = 0.0,
    eos_integration_method = 'DOP853',
    eos_rtol = 1.0e-3,
    eos_atol = 1.0e-5,
    eos_pressure_tol = 1.0e-3,
    eos_max_iters = 40,
    # Error and log reporting
    verbose = False,
    warnings = True,
    raise_on_fail = False,
    perform_checks = True,
    log_info = False
)

cmm_optional_args = dict(
    use_general_obliquity = False,
    solve_load_numbers = False,
    force_mode_calculation = False,
    degree_l = 2,
    use_modes = True,
    use_static_potential = False,
    use_simple_potential = False,
    orbit_average_results = True,
    use_kamata = False,
    integration_method = 'DOP853',
    integration_rtol = 1.0e-6,
    integration_atol = 1.0e-8,
    verbose = False,
    nondimensionalize = True
)