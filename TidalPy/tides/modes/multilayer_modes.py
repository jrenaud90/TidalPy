""" Multilayer calculator for multiple tidal modes.

Each tidal mode imparts a potentially different frequency which needs to be propagated throughout the planet's interior.

This module contains functions to assist with calculating the response at each of these frequencies and then collapse
    the findings into a final value.

"""

from typing import Callable, Dict, Tuple

import numpy as np

from ..multilayer.stress_strain import calculate_strain_stress
from ...utilities.types import FloatArray
from ..potential import tidal_potential_nsr_modes, tidal_potential_simple, tidal_potential_nsr
from ...toolbox.multilayer import KNOWN_INTERIOR_MODELS


def calculate_tidal_y(
    complex_compliance_func: Callable, forcing_frequency: FloatArray,
    radius_array: np.ndarray, gravity_array: np.ndarray, density_array,
    shear_array: np.ndarray, bulk_array: np.ndarray, viscosity_array: np.ndarray,
    interior_layer_model: str,
    layer_index_tuple: Tuple[np.ndarray], layer_is_static_tuple: Tuple[bool, ...],
    order_l: int = 2, complex_compliance_input: Tuple[float, ...] = None,
    **interior_integration_kwargs):
    """ Wrapper for the tidal-y integrator that allows for different interior models.

    Parameters
    ----------
    complex_compliance_func : Callable
        Complex compliance function
    forcing_frequency : FloatArray
        Tidal forcing frequency [rad s-1]
    radius_array : np.ndarray
        Radius array through the planet [m]
    gravity_array : np.ndarray
        Acceleration due to gravity array through the planet [m s-2]
    density_array : np.ndarray
        Local density array through the planet [kg m-3]
    shear_array : np.ndarray
        Shear modulus array through the planet [Pa]
    bulk_array : np.ndarray
        Bulk modulus array through the planet [Pa]
    viscosity_array : np.ndarray
        Viscosity array through the planet [Pa s]
    interior_layer_model : str
        Interior model used for radial solution
    layer_index_tuple : Tuple[np.ndarray, ...]
        Tuple of indices for each internal layer.
    layer_is_static_tuple : Tuple[bool, ...]
        Tuple of flags if each internal layer is static or not.
    order_l : int = 2
        Tidal harmonic order.
    complex_compliance_input : Tuple[float, ...]
        Tuple of floats used in the complex compliance function.
    interior_integration_kwargs : dict
        Keyword arguments for the interior integration function.

    Returns
    -------
    complex_shears : np.ndarray
        Complex shear modulus throughout the planet [Pa]
    tidal_y : np.ndarray
        Viscoelastic-gravitational solution as a function of radius.
    tidal_y_deriv : np.ndarray
        Radial derivative of the viscoelastic-gravitational solution as a function of radius.

    """

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
        np.real(complex_shears[np.imag(complex_shears) == 0.]) + 1.0j * 1e-50

    # Find interface radii for each layer
    num_layers = len(layer_index_tuple)
    radius_of_interfaces = [radius_array[layer_index][-1] for layer_index in layer_index_tuple]

    # Determine input based on the model used.
    interior_calc_input = None
    if interior_layer_model.lower() not in KNOWN_INTERIOR_MODELS:
        raise KeyError('Unknown interior model provided to find_tidal_y.')
    tidal_y_function = KNOWN_INTERIOR_MODELS[interior_layer_model.lower()]

    # Build input to the model based on which model is being used.
    if num_layers == 1:
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
        tidal_y_function(
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
    order_l: int = 2, tidal_y_integration_kwargs: dict = None,
    force_mode_calculation: bool = False):
    """ Given a tidal frequency, this function will call on the interior integration routine with the proper inputs and
        collect the results as well as calculate tidal heating as a function of radius, longitude, latitude, and time.

    Parameters
    ----------
    radius_matrix : np.ndarray
        Radius matrix with increased dimensions [m].
        Shape: (N_r, N_long, N_colat, N_time)
    colatitude_matrix : np.ndarray
        Colatitude matrix with increased dimensions [rads].
        Shape: (N_r, N_long, N_colat, N_time)
    time_domain : np.ndarray
        Time domain [s]
    longitude_domain : np.ndarray
        Longitude domain [rad]
    colatitude_domain : np.ndarray
        Colatitude domain [rad]
    radius_array : np.ndarray
        Radius array [m]
    gravity_array : np.ndarray
        Acceleration due to gravity as a function of radius [m s-2]
    density_array : np.ndarray
        Density as a function of radius [kg m-3]
    shear_array : np.ndarray
        Shear modulus as a function of radius [Pa]
    bulk_array : np.ndarray
        Bulk modulus as a function of radius [Pa]
    viscosity_array : np.ndarray
        Viscosity as a function of radius [Pa s]
    mode_frequency : float
        Tidal forcing frequency at requested mode [rad s-1]
    tidal_potential_tuple : Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        Inputs used to calculate the tidal potential
    complex_compliance_function : Callable
        Complex compliance function (set by rheology)
    interior_layer_model : str
        Interior model type used to propagate viscoelastic-gravitational solution.
    layer_index_tuple : Tuple[np.ndarray, ...]
        Tuple of boolean arrays indicating which radius values are associated with which layer.
    layer_is_static_tuple : Tuple[bool, ...]
        Tuple of booleans to indicate which layers should be treated with the static assumption.
    complex_compliance_input : Tuple[float, ...]
        Inputs used to calculate the complex compliance.
    order_l : int = 2
        Tidal harmonic order.
    tidal_y_integration_kwargs : dict = None
        Keyword arguments passed to the interior integrator function.
    force_mode_calculation : bool = False
        If True, then interior propagation will be performed for all modes regardless of frequency value.

    Returns
    -------
    mode_skipped : bool
        If True, then the mode's frequency was close to zero so no propagation was calculated.
    strains_at_mode : np.ndarray
        Tidal strains as a function of radius, longitude, colatitude, and time.
    stresses_at_mode : np.ndarray
        Tidal stresses as a function of radius, longitude, colatitude, and time.
    complex_shears_at_mode : np.ndarray
        Complex shear modulus as a function of radius.
    tidal_y_at_mode : np.ndarray
        Viscoelastic-gravitational radial solutions as a function of radius.

    """

    # Setup flags
    mode_skipped = False

    # Clean input
    if tidal_y_integration_kwargs is None:
        tidal_y_integration_kwargs = dict()

    # Calculate rheology and radial response
    if (not force_mode_calculation) and (mode_frequency < 1.0e-15):
        # If frequency is ~ 0.0 then there will be no tidal response. Skip the calculation of tidal y, etc.
        tidal_y_at_mode = np.zeros((6, *radius_array.shape), dtype=np.complex128)
        tidal_y_derivative_at_mode = np.zeros((6, *radius_array.shape), dtype=np.complex128)
        complex_shears_at_mode = shear_array + 0.0j
        strains_at_mode = np.zeros((6, *radius_matrix.shape), dtype=np.complex128)
        stresses_at_mode = np.zeros((6, *radius_matrix.shape), dtype=np.complex128)
        mode_skipped = True
    else:
        # Calculate the radial functions using a shooting integration method.
        complex_shears_at_mode, tidal_y_at_mode, tidal_y_derivative_at_mode = \
            calculate_tidal_y(
                complex_compliance_function, mode_frequency,
                radius_array, gravity_array, density_array, shear_array, bulk_array, viscosity_array,
                interior_layer_model=interior_layer_model,
                layer_index_tuple=layer_index_tuple, layer_is_static_tuple=layer_is_static_tuple,
                order_l=order_l, complex_compliance_input=complex_compliance_input,
                **tidal_y_integration_kwargs)

        # We need to recast the inputs into the correct dimensions
        complex_shears_higher_dim, _, _, _ = \
            np.meshgrid(complex_shears_at_mode, longitude_domain, colatitude_domain, time_domain, indexing='ij')

        # TODO: Bulk dissipation has not been implemented.
        complex_bulks_higher_dim, _, _, _ = \
            np.meshgrid(bulk_array, longitude_domain, colatitude_domain, time_domain, indexing='ij')

        # We need to recast the tidal y solution into the correct dimensions
        #    (it does not care about long/lat or time - at least in this context).
        # OPT: It will probably be better from a memory conservation point of view if we don't make these higher dimension variables and instead use a non-numpy array technique during the potential calculations.
        tidal_y_higher_dim = np.repeat(tidal_y_at_mode[:, :, np.newaxis],
                                       len(longitude_domain), axis=2)
        tidal_y_higher_dim = np.repeat(tidal_y_higher_dim[:, :, :, np.newaxis],
                                       len(colatitude_domain), axis=3)
        tidal_y_higher_dim = np.repeat(tidal_y_higher_dim[:, :, :, :, np.newaxis],
                                       len(time_domain), axis=4)

        tidal_y_derivative_higher_dim = np.repeat(tidal_y_at_mode[:, :, np.newaxis],
                                                  len(longitude_domain), axis=2)
        tidal_y_derivative_higher_dim = np.repeat(tidal_y_derivative_higher_dim[:, :, :, np.newaxis],
                                                  len(colatitude_domain), axis=3)
        tidal_y_derivative_higher_dim = np.repeat(tidal_y_derivative_higher_dim[:, :, :, :, np.newaxis],
                                                  len(time_domain), axis=4)

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
            tidal_solution_y=tidal_y_higher_dim, tidal_solution_y_derivative=tidal_y_derivative_higher_dim,
            colatitude=colatitude_matrix,
            radius=radius_matrix, shear_moduli=complex_shears_higher_dim,
            bulk_moduli=complex_bulks_higher_dim,
            frequency=mode_frequency
            )

    return mode_skipped, strains_at_mode, stresses_at_mode, complex_shears_at_mode, tidal_y_at_mode

def collapse_multilayer_modes(
    shear_modulus: np.ndarray, viscosity: np.ndarray, bulk_modulus: np.ndarray,
    density: np.ndarray, gravity: np.ndarray,
    complex_compliance_func: Callable,
    radius_matrix, longitude_matrix, colatitude_matrix, time_matrix, voxel_volume,
    orbital_frequency, spin_frequency, semi_major_axis, eccentricity, host_mass,
    propagation_function: Callable,
    interface_properties: dict,
    order_l: int = 2, complex_compliance_input: tuple = None,
    use_modes: bool = True, use_static_potential: bool = False, use_simple_potential: bool = False,
    orbit_average_results: bool = True,
    integration_parameters: dict = None
    ):
    """ Calculate the multilayer tidal response of a planet over a range of applicable tidal modes. Collapse
    individual modal results into final heating distribution.

    Some of these variables are ndarrays the shape of the planet's radius array, N_r
    Others are multidimensional arrays the shape of (N_r, N_long, N_colat, N_time).
    These multidim arrays MUST be in this order [radius_domain, longitude_domain, colatitude_domain, time_domain]

    Parameters
    ----------
    shear_modulus : np.ndarray
        Non-complex shear modulus as a function of radius [Pa].
        Shape: N_r
    viscosity : np.ndarray
        Viscosity as a function of radius [Pa].
        Shape: N_r
    bulk_modulus : np.ndarray
        Non-complex bulk modulus as a function of radius [Pa].
        Shape: N_r
    density : np.ndarray
        Density as a function of radius [Pa].
        Shape: N_r
    gravity : np.ndarray
        Gravity as a function of radius [Pa].
        Shape: N_r
    complex_compliance_func : Callable
        Complex compliance function to use based on the rheology. See TidalPy.rheology.complex_compliance.
    radius_matrix : np.ndarray
        Radius matrix with increased dimensions [m].
        Shape: (N_r, N_long, N_colat, N_time)
    longitude_matrix : np.ndarray
        Longitude matrix with increased dimensions [rads].
        Shape: (N_r, N_long, N_colat, N_time)
    colatitude_matrix : np.ndarray
        Colatitude matrix with increased dimensions [rads].
        Shape: (N_r, N_long, N_colat, N_time)
    time_matrix : np.ndarray
        Time matrix with increased dimensions [s].
        It is important that the time domain be for one orbital period.
        Shape: (N_r, N_long, N_colat, N_time)
    voxel_volume : np.ndarray
        Volume per voxel matrix with increased dimensions [m3].
        Shape: (N_r, N_long, N_colat, N_time)
    orbital_frequency : float
        Orbital mean motion [rad s-1]
    spin_frequency : float
        Rotation frequency [rad s-1]
    semi_major_axis : float
        Orbital semi-major axis [m]
    eccentricity : float
        Orbital eccentricity
    host_mass : float
        Mass of tidal host [kg]
    propagation_function : Callable
        The interior propagation function to be used.
        See TidalPy.toolbox.multilayer for more information.
    interface_properties : dict
        Dictionary of properties for the planet's internal layer interfaces.
        See TidalPy.toolbox.multilayer for more information.
    order_l : int = 2
        Tidal harmonic order.
    complex_compliance_input : tuple = None
        Additional inputs to the complex compliance function. If None, then defaults will be used.
    use_modes : bool = True
        If True, the interior integration will occur across multiple tidal modes.
        This can be set to False if the planet is tidally locked AND the eccentricity is low (e <~ 0.1) and is
           not expected to increase.
    use_static_potential : bool = False
        If True, then static terms within the tidal potential (usually phase terms like sin(2*phi) will be included.
        These terms should not be used to calculate tidal heating since it is a time derivative of the potential.
    use_simple_potential : bool = False
        If True, then a simplified version of the tidal potential will be used when use_modes is set to False.
    orbit_average_results : bool = True
        If True, then the function will orbit average the heating, stress, and strain results. This will reduce the
        final output's dimension by one.
    integration_parameters : dict = None
        Dictionary of integration parameters used to integrate through the planet's interior.
        See TidalPy.toolbox.multilayer for more information.

    Returns
    -------
    heating : np.ndarray
        Heating within each voxel [W].
        Shape: (N_r, N_long, N_colat)
        If orbit_average_results is false then the results will be given at each time with a new shape of:
            (N_r, N_long, N_colat, N_time)
    volumetric_heating : np.ndarray
        Volumetric Heating within each voxel [W].
        Shape: (N_r, N_long, N_colat)
        If orbit_average_results is false then the results will be given at each time with a new shape of:
            (N_r, N_long, N_colat, N_time)
    volumetric_heating_by_mode : Dict[str, np.ndarray]
        Volumetric heating within each voxel broken up by each tidal mode [W].
        Shape of each stored array: (N_r, N_long, N_colat)
        If orbit_average_results is false then the results will be given at each time with a new shape of:
            (N_r, N_long, N_colat, N_time)
    strains : np.ndarray
        Tidal strains within each voxel [m m-1].
        Shape: (N_r, N_long, N_colat)
        If orbit_average_results is false then the results will be given at each time with a new shape of:
            (N_r, N_long, N_colat, N_time)
    stresses : np.ndarray
        Tidal stresses within each voxel [m m-1].
        Shape: (N_r, N_long, N_colat)
        If orbit_average_results is false then the results will be given at each time with a new shape of:
            (N_r, N_long, N_colat, N_time)

    """

    # If no integration parameters were provided then just use an empty dictionary which will tell the function
    #     to use its default values.
    if integration_parameters is None:
        integration_parameters = dict()

    # If no inputs to the complex compliance function were provided then set it equal to an empty tuple which
    #    will cause the complex compliance function to resort to defaults.
    if complex_compliance_input is None:
        complex_compliance_input = tuple()

    # Variable may be a matrix calculated from a mesh grid over [radius_matrix, longitude_matrix, latitude, time_matrix] <- must be in this order

    # Pull out individual arrays
    radius_array = radius_matrix[:, 0, 0, 0]
    longitude_domain = longitude_matrix[0, :, 0, 0]
    colatitude_domain = colatitude_matrix[0, 0, :, 0]
    time_domain = time_matrix[0, 0, 0, :]
    shape = time_matrix.shape

    # Check that dimensions make sense
    assert radius_array.shape == shear_modulus.shape

    # Check that the time domain has the correct end points
    orbital_period = 2. * np.pi / orbital_frequency
    if orbit_average_results:
        # In order for the orbit average routine to work correctly, the time domain must start at zero and end
        #    after 1 orbital period.
        assert time_domain[0] == 0.
        assert time_domain[-1] == orbital_period

    # TODO: Currently this function (and other multilayer code) does not work for l>2. An additional loop will be
    #   required to loop over each l and sum the results.
    if order_l != 2:
        raise Exception

    # Expand voxel volume dimensions to include a time domain
    voxel_volume_higher_dim = np.repeat(voxel_volume[:, :, :, np.newaxis], len(time_domain), axis=3)

    # Calculate the tidal modes and the tidal potential and its partial derivatives.
    planet_radius = radius_array[-1]

    if use_modes:
        # First calculate the modes and potential
        tidal_frequencies, tidal_modes, potential_dict, potential_dtheta_dict, potential_dphi_dict, potential_d2theta_dict, \
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
                    radius_matrix, longitude_matrix, colatitude_matrix, orbital_frequency, eccentricity, time_matrix
                    )
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
                radius_array, gravity, density, shear_modulus, bulk_modulus, viscosity,
                mode_frequency, potential_tuple, complex_compliance_function,
                interior_layer_model=interior_layer_model, layer_index_tuple=layer_index_tuple,
                layer_is_static_tuple=layer_is_static_tuple, complex_compliance_input=complex_compliance_input,
                order_l=order_l, tidal_y_integration_kwargs=tidal_y_integration_kwargs
                )

        if mode_skipped:
            num_modes_skipped += 1

        # Collapse Modes
        # Add items defined for each mode to lists
        modes_skipped[mode_name] = mode_skipped
        love_k_by_mode[mode_name] = tidal_y_at_mode[4, -1] - 1.
        love_h_by_mode[mode_name] = tidal_y_at_mode[0, -1] * gravity[-1]
        love_l_by_mode[mode_name] = tidal_y_at_mode[2, -1] * gravity[-1]

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

    # Pull out individual stresses and strains
    e_rr, e_thth, e_phph, e_rth, e_rph, e_thph = strains
    s_rr, s_thth, s_phph, s_rth, s_rph, s_thph = stresses

    # Perform orbital averaging
    if orbit_average_results:
        e_rr = (1. / orbital_period) * np.trapz(e_rr, time_domain, axis=3)
        e_thth = (1. / orbital_period) * np.trapz(e_thth, time_domain, axis=3)
        e_phph = (1. / orbital_period) * np.trapz(e_phph, time_domain, axis=3)
        e_rth = (1. / orbital_period) * np.trapz(e_rth, time_domain, axis=3)
        e_rph = (1. / orbital_period) * np.trapz(e_rph, time_domain, axis=3)
        e_thph = (1. / orbital_period) * np.trapz(e_thph, time_domain, axis=3)
        stresses = np.zeros((6, *radius_matrix.shape[:-1]), dtype=np.complex128)
        stresses[0, :, :, :] = e_rr
        stresses[1, :, :, :] = e_thth
        stresses[2, :, :, :] = e_phph
        stresses[3, :, :, :] = e_rth
        stresses[4, :, :, :] = e_rph
        stresses[5, :, :, :] = e_thph

        s_rr = (1. / orbital_period) * np.trapz(s_rr, time_domain, axis=3)
        s_thth = (1. / orbital_period) * np.trapz(s_thth, time_domain, axis=3)
        s_phph = (1. / orbital_period) * np.trapz(s_phph, time_domain, axis=3)
        s_rth = (1. / orbital_period) * np.trapz(s_rth, time_domain, axis=3)
        s_rph = (1. / orbital_period) * np.trapz(s_rph, time_domain, axis=3)
        s_thph = (1. / orbital_period) * np.trapz(s_thph, time_domain, axis=3)
        strains = np.zeros((6, *radius_matrix.shape[:-1]), dtype=np.complex128)
        strains[0, :, :, :] = s_rr
        strains[1, :, :, :] = s_thth
        strains[2, :, :, :] = s_phph
        strains[3, :, :, :] = s_rth
        strains[4, :, :, :] = s_rph
        strains[5, :, :, :] = s_thph

    # Calculate final heating as a function of strains and stresses
    # TODO: which frequency... Using orbital freq for now.
    frequency = orbital_frequency

    # Calculate Tidal Heating
    volumetric_heating = (frequency / 2.) * (
            np.imag(s_rr) * np.real(e_rr) - np.real(s_rr) * np.imag(e_rr) +
            np.imag(s_thth) * np.real(e_thth) - np.real(s_thth) * np.imag(e_thth) +
            np.imag(s_phph) * np.real(e_phph) - np.real(s_phph) * np.imag(e_phph) +
            2. * (np.imag(s_rth) * np.real(e_rth) - np.real(s_rth) * np.imag(e_rth)) +
            2. * (np.imag(s_rph) * np.real(e_rph) - np.real(s_rph) * np.imag(e_rph)) +
            2. * (np.imag(s_thph) * np.real(e_thph) - np.real(s_thph) * np.imag(e_thph))
    )

    if orbit_average_results:
        heating = volumetric_heating * voxel_volume
    else:
        heating = volumetric_heating * voxel_volume_higher_dim

    return heating, volumetric_heating, strains, stresses
