import numpy as np

from ...utilities.types import FloatArray
from ...utilities.performance import njit
from ...tides.potential import tidal_potential_nsr_modes, tidal_potential_nsr, tidal_potential_simple
from ...tides.multilayer.stress_strain import calculate_strain_stress_heating
from . import KNOWN_INTERIOR_MODELS
from typing import Tuple, Callable



def calculate_tidal_y(
    complex_compliance_func: Callable, forcing_frequency: FloatArray,
    radius_array: np.ndarray, gravity_array: np.ndarray, density_array,
    shear_array: np.ndarray, bulk_array: np.ndarray, viscosity_array: np.ndarray,
    interior_layer_model: str,
    layer_index_tuple: Tuple[np.ndarray], layer_is_static_tuple: Tuple[bool, ...],
    order_l: int = 2, complex_compliance_input: Tuple[float, ...] = None,
    **interior_integration_kwargs):

    # Calculate Complex Compliance
    if complex_compliance_input is None:
        complex_compliance_input = tuple()
    complex_compliances = \
        complex_compliance_func(forcing_frequency, shear_array**(-1), viscosity_array, *complex_compliance_input)
    complex_shears = complex_compliances**(-1)

    # Clean up complex shears based on any zero viscosities (short hand to indicate zero dissipation)
    complex_shears[viscosity_array == 0.] = shear_array[viscosity_array == 0.] * (1. + 0.j)
    # Make zero shears very small instead to avoid 1/0 issues.
    complex_shears[shear_array == 0.] = (1. + 1.j) * 1.e-50
    complex_shears[np.imag(complex_shears) == 0.] = \
        np.real(complex_shears[np.imag(complex_shears) == 0.]) + 1.0j * 1e-40

    # Find interface radii for each layer
    num_layers = len(layer_index_tuple)
    radius_of_interfaces = [radius_array[layer_index][-1] for layer_index in layer_index_tuple]

    # Determine input based on the model used.
    interior_calc_input = None
    if interior_layer_model.lower() not in KNOWN_INTERIOR_MODELS:
        raise KeyError('Unknown interior model provided to find_tidal_y.')
    calculate_tidal_y = KNOWN_INTERIOR_MODELS[interior_layer_model.lower()]

    # Build input to the model based on which model is being used.
    if num_layers > 1:
        if interior_layer_model.lower() not in ['homogeneous']:
            raise ValueError('Mismatch between number of layers and interior model method.')
        else:
            interior_calc_input = {
                'use_static': layer_is_static_tuple[0]
                }

    elif num_layers == 2:
        if interior_layer_model.lower() not in ['liquid-solid']:
            raise ValueError('Mismatch between number of layers and interior model method.')
        else:
            interior_calc_input = {
                'interface_1_radius': radius_of_interfaces[0],
                'layer_0_static': layer_is_static_tuple[0],
                'layer_1_static': layer_is_static_tuple[1]
                 }

    elif num_layers == 3:
        if interior_layer_model.lower() not in ['solid-liquid-solid']:
            raise ValueError('Mismatch between number of layers and interior model method.')
        else:
            interior_calc_input = {
                'interface_1_radius': radius_of_interfaces[0],
                'interface_2_radius': radius_of_interfaces[1],
                'layer_0_static'    : layer_is_static_tuple[0],
                'layer_1_static'    : layer_is_static_tuple[1],
                'layer_2_static'    : layer_is_static_tuple[2]
                }
    elif num_layers == 4:
        if interior_layer_model.lower() not in ['solid-solid-liquid-solid']:
            raise ValueError('Mismatch between number of layers and interior model method.')
        else:
            interior_calc_input = {
                'interface_1_radius': radius_of_interfaces[0],
                'interface_2_radius': radius_of_interfaces[1],
                'interface_3_radius': radius_of_interfaces[2],
                'layer_0_static'    : layer_is_static_tuple[0],
                'layer_1_static'    : layer_is_static_tuple[1],
                'layer_2_static'    : layer_is_static_tuple[2],
                'layer_3_static'    : layer_is_static_tuple[3]
                }
    else:
        raise NotImplementedError('Only 1 to 4 layers are currently supported by interior integrator.')

    # Calculate the radial solution to the viscoelastic-gravitational problem
    tidal_y, tidal_y_deriv = \
        calculate_tidal_y(
            radius_array, complex_shears, bulk_array, density_array, gravity_array, forcing_frequency,
            **interior_calc_input, order_l=order_l,
            **interior_integration_kwargs)

    return complex_shears, tidal_y, tidal_y_deriv


def calculate_mode_response_coupled(
    radius_matrix: np.ndarray, colatitude_matrix: np.ndarray,
    time_domain, longitude_domain, colatitude_domain,
    radius_array: np.ndarray, gravity_array: np.ndarray, density_array,
    shear_array: np.ndarray, bulk_array: np.ndarray, viscosity_array: np.ndarray,
    mode_frequency: FloatArray,
    tidal_potential_tuple: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    complex_compliance_function: Callable, interior_layer_model: str,
    layer_index_tuple: Tuple[np.ndarray, ...], layer_is_static_tuple: Tuple[bool, ...],
    complex_compliance_input: Tuple[float, ...] = None,
    order_l: int = 2, tidal_y_integration_kwargs = None):

    # Setup flags
    mode_skipped = False

    # Clean input
    if tidal_y_integration_kwargs is None:
        tidal_y_integration_kwargs = dict()

    # Calculate rheology and radial response
    if mode_frequency < 1.0e-15:
        # If frequency is ~ 0.0 then there will be no tidal response. Skip the calculation of tidal y, etc.
        tidal_y_at_mode = np.zeros((6, *radius_array.shape), dtype=np.complex128)
        tidal_y_deriv_at_mode = np.zeros((6, *radius_array.shape), dtype=np.complex128)
        complex_shears_at_mode = shear_array + 0.0j
        strains_at_mode = np.zeros((6, *radius_matrix.shape), dtype=np.complex128)
        stresses_at_mode = np.zeros((6, *radius_matrix.shape), dtype=np.complex128)
        mode_skipped = True
    else:
        # Calculate the radial functions using a shooting integration method.
        complex_shears_at_mode, tidal_y_at_mode, tidal_y_deriv_at_mode = \
            calculate_tidal_y(
                complex_compliance_function, mode_frequency,
                radius_array, gravity_array, density_array, shear_array, bulk_array, viscosity_array,
                interior_layer_model=interior_layer_model,
                layer_index_tuple=layer_index_tuple, layer_is_static_tuple=layer_is_static_tuple,
                order_l=order_l, complex_compliance_input=complex_compliance_input,
                **tidal_y_integration_kwargs)

        # We need to recast the inputs into the correct dimensions
        complex_shears_higher_dim, _, _, _ = \
            np.meshgrid(complex_shears_at_mode, time_domain, longitude_domain, colatitude_domain, indexing='ij')

        # TODO: Bulk dissipation has not been implemented.
        complex_bulks_higher_dim, _, _, _ = \
            np.meshgrid(bulk_array, time_domain, longitude_domain, colatitude_domain, indexing='ij')

        # We need to recast the tidal y solution into the correct dimensions
        #    (it does not care about long/lat or time - at least in this context).
        tidal_y_higher_dim = np.zeros((6, *radius_matrix.shape), dtype=np.complex128)
        tidal_y_deriv_higher_dim = np.zeros((6, *radius_matrix.shape), dtype=np.complex128)
        for i in range(6):
            tidal_y_higher_dim[i, :, :, :, :], _, _, _ = \
                np.meshgrid(tidal_y_at_mode[i, :], time_domain, longitude_domain, colatitude_domain, indexing='ij')
            tidal_y_deriv_higher_dim[i, :, :, :, :], _, _, _ = \
                np.meshgrid(tidal_y_deriv_at_mode[i, :], time_domain, longitude_domain, colatitude_domain, indexing='ij')

        # Unpack tidal potential terms
        potential, potential_dtheta, potential_dphi, potential_d2theta, potential_d2phi, potential_dtheta_dphi = \
            tidal_potential_tuple
        # Calculate stresses and heating
        strains_at_mode, stresses_at_mode, OLD_volumetric_heating_at_mode = calculate_strain_stress_heating(
            tidal_potential=potential,
            tidal_potential_partial_theta=potential_dtheta, tidal_potential_partial_phi=potential_dphi,
            tidal_potential_partial2_theta2=potential_d2theta,
            tidal_potential_partial2_phi2=potential_d2phi,
            tidal_potential_partial2_theta_phi=potential_dtheta_dphi,
            tidal_solution_y=tidal_y_higher_dim, tidal_solution_y_derivative=tidal_y_deriv_higher_dim,
            colatitude=colatitude_matrix,
            radius=radius_matrix, shear_moduli=complex_shears_higher_dim,
            bulk_moduli=complex_bulks_higher_dim,
            frequency=mode_frequency
            )

    return mode_skipped, strains_at_mode, stresses_at_mode, complex_shears_at_mode, tidal_y_at_mode


def calculate_multimode_response_coupled(
    radius_array, volume_array, voxel_volumes, gravity_array, density_array, shear_array, bulk_array, viscosity_array,
    radius_matrix, time_matrix, longitude_matrix, colatitude_matrix,
    orbital_frequency, spin_frequency, eccentricity, obliquity, host_mass, semi_major_axis,
    complex_compliance_function: Callable, interior_layer_model: str,
    layer_is_static_tuple: Tuple[bool, ...], layer_index_tuple: Tuple[np.ndarray, ...],
    complex_compliance_input: Tuple[float, ...] = None,
    tidally_locked: bool = False, use_modes: bool = True, use_simple_potential: bool = False, use_static_potential: bool = False,
    scale_average_by_normalized_freq: bool = True,
    order_l: int = 2, tidal_y_integration_kwargs: dict = None):

    # Pull out arrays from matrices
    longitude_domain = longitude_matrix[0, 0, :, 0]
    colatitude_domain = colatitude_matrix[0, 0, 0, :]
    time_domain = time_matrix[0, :, 0, 0]
    shape = time_matrix.shape

    # For the voxel volumes we need to insert a new axis for time (since it already has ones for radius, lat, long).
    voxel_volume_array_higher_dim = np.repeat(voxel_volumes[:, np.newaxis, :, :], len(time_domain), axis=1)

    if use_modes and (not tidally_locked):
        # First calculate the modes and potential
        tidal_frequencies, tidal_modes, potential_dict, potential_dtheta_dict, potential_dphi_dict, potential_d2theta_dict,\
            potential_d2phi_dict, potential_dtheta_dphi_dict = \
            tidal_potential_nsr_modes(
                radius_matrix, longitude_matrix, colatitude_matrix, orbital_frequency, eccentricity, time_matrix,
                spin_frequency, world_radius=radius_array[-1], host_mass=host_mass, semi_major_axis=semi_major_axis,
                use_static=use_static_potential
                )
    else:
        # Calculate the non-mode version of the tidal potential
        if use_simple_potential:
            # Calculate the simplified tidal potential
            potential, potential_dtheta, potential_dphi, potential_d2theta, potential_d2phi, potential_dtheta_dphi = \
                tidal_potential_simple(
                    radius_matrix, longitude_matrix, colatitude_matrix, orbital_frequency, eccentricity, time_matrix)
        else:
            potential, potential_dtheta, potential_dphi, potential_d2theta, potential_d2phi, potential_dtheta_dphi = \
                tidal_potential_nsr(
                    radius_matrix, longitude_matrix, colatitude_matrix, orbital_frequency, eccentricity, time_matrix,
                    spin_frequency, world_radius=radius_array[-1], host_mass=host_mass,
                    semi_major_axis=semi_major_axis, use_static=use_static_potential
                )
        # Add the results to a mode dictionary so that the mode vs. non-mode calculation steps are identical.
        tidal_frequencies = {'n': np.abs(orbital_frequency)}
        tidal_modes = {'n': orbital_frequency}
        potential_dict = {'n': potential}
        potential_dtheta_dict = {'n': potential_dtheta}
        potential_dphi_dict = {'n': potential_dphi}
        potential_d2theta_dict = {'n': potential_d2theta}
        potential_d2phi_dict = {'n': potential_d2phi}
        potential_dtheta_dphi_dict = {'n': potential_dtheta_dphi}

    # Record how many modes are skipped (used in average)
    num_modes_skipped = 0

    # Build storages for all modes
    modes_skipped = dict()
    love_k_by_mode = dict()
    love_h_by_mode = dict()
    love_l_by_mode = dict()

    # Large arrays must be added to continuously to avoid memory overload
    complex_shears_avg = np.zeros(shape, dtype=np.complex128)
    tidal_y_avg = np.zeros(shape, dtype=np.complex128)
    stresses = np.zeros(shape, dtype=np.complex128)
    strains = np.zeros(shape, dtype=np.complex128)
    # TODO: Tobie 2005 shows the tidal potential without any radius dependence. Not sure why that is the case.
    tidal_potential = np.zeros(shape, dtype=np.complex128)
    total_potential = np.zeros(shape, dtype=np.complex128)

    # Since it is unclear how to collapse some of the items that depend on tidal modes, the user can choose
    #    to average by the normalized frequency or just by the average. The max and min frequency are required to
    #    perform the scaling
    max_freq = np.max(list(tidal_frequencies.values()))
    min_freq = np.min(list(tidal_frequencies.values()))

    # TODO: multiprocessor
    for mode_name, mode_frequency in tidal_frequencies.items():
        potential_tuple = (
            potential_dict[mode_name], potential_dtheta_dict[mode_name], potential_dphi_dict[mode_name],
            potential_d2theta_dict[mode_name], potential_d2phi_dict[mode_name],
            potential_dtheta_dphi_dict[mode_name])

        # Calculate response at mode
        mode_skipped, strains_at_mode, stresses_at_mode, complex_shears_at_mode, tidal_y_at_mode = \
            calculate_mode_response_coupled(
                radius_matrix, colatitude_matrix,
                time_domain, longitude_domain, colatitude_domain,
                radius_array, gravity_array, density_array, shear_array, bulk_array, viscosity_array,
                mode_frequency, potential_tuple, complex_compliance_function,
                interior_layer_model=interior_layer_model, layer_index_tuple=layer_index_tuple,
                layer_is_static_tuple=layer_is_static_tuple, complex_compliance_input=complex_compliance_input,
                order_l=order_l, tidal_y_integration_kwargs=tidal_y_integration_kwargs)

        if mode_skipped:
            num_modes_skipped += 1

        # Collapse Modes
        # Add items defined for each mode to lists
        modes_skipped[mode_name] = mode_skipped
        love_k_by_mode[mode_name] = tidal_y_at_mode[4, -1] - 1.
        love_h_by_mode[mode_name] = tidal_y_at_mode[0, -1] * gravity_array[-1]
        love_l_by_mode[mode_name] = tidal_y_at_mode[2, -1] * gravity_array[-1]

        # Stresses, strains, and potentials are added together.
        stresses += stresses_at_mode
        strains += strains_at_mode
        tidal_potential += potential_dict[mode_name]
        total_potential += potential_dict[mode_name] * tidal_y_at_mode[4, :]

        # The other parameters it is not clear what category they fall into. TODO: For now let's take the average of them.
        complex_shears_avg += complex_shears_at_mode
        tidal_y_avg += tidal_y_at_mode

    # Finish taking the average of the avg parameters
    complex_shears_avg = complex_shears_avg / max((len(tidal_modes) - num_modes_skipped), 1)
    tidal_y_avg = tidal_y_avg / max((len(tidal_modes) - num_modes_skipped), 1)

    # Calculate final heating as a function of strains and stresses
    # TODO: which frequency... Using orbital freq for now.
    frequency = orbital_frequency
    e_rr, e_thth, e_phph, e_rth, e_rph, e_thph = strains
    s_rr, s_thth, s_phph, s_rth, s_rph, s_thph = stresses
    volumetric_heating = (frequency / 2.) * (
            np.imag(s_rr) * np.real(e_rr) - np.real(s_rr) * np.imag(e_rr) +
            np.imag(s_thth) * np.real(e_thth) - np.real(s_thth) * np.imag(e_thth) +
            np.imag(s_phph) * np.real(e_phph) - np.real(s_phph) * np.imag(e_phph) +
            2. * (np.imag(s_rth) * np.real(e_rth) - np.real(s_rth) * np.imag(e_rth)) +
            2. * (np.imag(s_rph) * np.real(e_rph) - np.real(s_rph) * np.imag(e_rph)) +
            2. * (np.imag(s_thph) * np.real(e_thph) - np.real(s_thph) * np.imag(e_thph))
    )
    heating = volumetric_heating * voxel_volume_array_higher_dim

    return complex_shears_avg, tidal_y_avg, strains, stresses, heating, time_domain