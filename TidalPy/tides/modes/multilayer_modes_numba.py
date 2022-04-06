""" Multilayer calculator for multiple tidal modes.

Each tidal mode imparts a potentially different frequency which needs to be propagated throughout the planet's interior.

This module contains functions to assist with calculating the response at each of these frequencies and then collapse
    the findings into a final value.

"""
from typing import Dict, List, Tuple

import numpy as np

from ..multilayer.numerical_int.solver_numba import tidal_y_solver
from ..multilayer.stress_strain import calculate_strain_stress
from ..potential import (TidalPotentialOutput, tidal_potential_nsr, tidal_potential_nsr_modes,
                         tidal_potential_obliquity_nsr, tidal_potential_obliquity_nsr_modes, tidal_potential_simple)
from ...utilities.performance import njit, nbList, nbDict, bool_, complex128, float64, nbUnicode


@njit(cacheable=True)
def calculate_mode_response_coupled(
    interior_model_name: str, mode_frequency: float,
    radius_array: np.ndarray, shear_array: np.ndarray, bulk_array: np.ndarray, viscosity_array: np.ndarray,
    density_array: np.ndarray, gravity_array: np.ndarray, colatitude_matrix: np.ndarray,
    tidal_potential_tuple: TidalPotentialOutput, complex_compliance_function: callable,
    is_solid_by_layer: List[bool], is_static_by_layer: List[bool], indices_by_layer: List[np.ndarray],
    surface_boundary_conditions: np.ndarray = None, solve_load_numbers: bool = False,
    complex_compliance_input: Tuple[float, ...] = tuple(), force_mode_calculation: bool = False,
    order_l: int = 2, use_kamata: bool = False,
    int_rtol: float = 1.0e-8, int_atol: float = 1.0e-12, rk_method: int = 1,
    verbose: bool = False, nondimensionalize: bool = True, planet_bulk_density: float = None
    ) -> Tuple[bool, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """ Given a tidal frequency, this function will call on the interior integration routine with the proper inputs and
        collect the results as well as calculate tidal stress, strain, and heating as a
        function of radius, longitude, latitude, and time.

    Parameters
    ----------
    interior_model_name : str
        Interior model used in the calculation of the radial functions.
        See options in TidalPy.tides.multilayer.numerical_int.collapse.
    mode_frequency : float
        Tidal forcing frequency at this mode [rad]
    radius_array : np.ndarray
        Radius array defined throughout the planet [m]
        Shape: (n_r,) Type: float
    shear_array : np.ndarray
        Static shear modulus array defined throughout the planet at each radius r [Pa]
        Shape: (n_r,) Type: float
    bulk_array : np.ndarray
        Static bulk modulus array defined throughout the planet at each radius r [Pa]
        Shape: (n_r,) Type: float
        # TODO: Bulk dissipation has not been implemented yet.
    viscosity_array : np.ndarray
        Viscosity array defined throughout the planet at each radius r [Pa s]
        Shape: (n_r,) Type: float
    density_array : np.ndarray
        Density array defined throughout the planet at each radius r [kg m-3]
        Shape: (n_r,) Type: float
    gravity_array : np.ndarray
        Acceleration due to gravity defined throughout the planet at the top of each radius r [m s-2]
        Shape: (n_r,) Type: float
    colatitude_matrix : np.ndarray
        Ndarray of colatitude defined at longitude, colatitude, and time [rad]
        Shape: (n_long, n_colat, n_time) Type: float
    tidal_potential_tuple : TidalPotentialOutput
        Results of the tidal potential calculated for this mode.
    complex_compliance_function : callable
        Complex compliance function used to find the complex shear modulus.
    is_solid_by_layer : List[bool]
        List or tuple of flags indicating which layers are solid (True) or liquid (False).
    is_static_by_layer : List[bool]
        List or tuple of flags indicating which layers are static (True) or dynamic (False).
    indices_by_layer : List[np.ndarray]
        List or tuple of np.ndarrays flagging which radius indices below to which layer.
        Each away: Shape: (n_r,) Type: np.boolean
    surface_boundary_conditions : np.ndarray = None
        Surface conditions applied to the radial function solutions at the top of the planet.
        If equal to None then the function will calculate the surface boundary conditions of either the tidal or load
        number.
    solve_load_numbers : bool = False
        If True, then the surface conditions will be calculated to find the load numbers rather than the tidal ones.
    complex_compliance_input  : Tuple[float, ...] = None
        Optional additional inputs to the complex compliance function.
    force_mode_calculation : bool = False
        If True, then this mode's tidal response will be calculated even if its frequency value falls below the
        minimum threshold.
    order_l: int = 2
        Tidal harmonic order.
    use_kamata : bool = False
        If True, the Kamata+2015 initial conditions will be used at the base of layer 0.
        Otherwise, the Takeuchi & Saito 1972 initial conditions will be used.
    int_rtol : float = 1.0e-6
        Integration relative error.
    int_atol : float = 1.0e-4
        Integration absolute error.
    rk_method : int = 1
        The type of RK method used for integration
            0 = RK23
            1 = RK45
            # TODO NotImplemented 2 = DOP
    verbose: bool = False
        If True, the function will print some information to console during calculation (may cause a slow down).
    nondimensionalize : bool = False
        If True, integration will use dimensionless variables. These will be converted back before output is given to
        the user.
    planet_bulk_density : float = None
        Must be provided if non_dimensionalize is True. Bulk density of the planet.

    Returns
    -------
    mode_skipped : bool
        If True, then the mode's frequency was close to zero so no propagation was calculated.
    strains_at_mode : np.ndarray
        Tidal strains as a function of radius, longitude, colatitude, and time. [m m-1]
    stresses_at_mode : np.ndarray
        Tidal stresses as a function of radius, longitude, colatitude, and time. [Pa]
    complex_shears_at_mode : np.ndarray
        Complex shear modulus as a function of radius. [Pa]
    tidal_y_at_mode : np.ndarray
        Viscoelastic-gravitational radial solutions as a function of radius.

    """

    # Setup flags
    mode_skipped = False

    # Calculate rheology and radial response
    if (not force_mode_calculation) and (mode_frequency < 1.0e-15):
        # If frequency is ~ 0.0 then there will be no tidal response. Skip the calculation of tidal y, etc.
        tidal_y_at_mode = np.zeros((6, *radius_array.shape), dtype=np.complex128)
        complex_shears_at_mode = shear_array + 0.0j
        strains_at_mode = np.zeros((6, *radius_array.shape, *colatitude_matrix.shape), dtype=np.complex128)
        stresses_at_mode = np.zeros((6, *radius_array.shape, *colatitude_matrix.shape), dtype=np.complex128)
        mode_skipped = True
    else:
        # Calculate Complex Compliance
        complex_compliances = \
            complex_compliance_function(mode_frequency, shear_array**(-1), viscosity_array, *complex_compliance_input)
        complex_shears_at_mode = complex_compliances**(-1)

        # TODO: Calculate complex bulk modulus. (When bulk dissipation is added)

        # Calculate the radial functions using a shooting integration method.
        # OPT: the option of scipy or julia integrators, rather than custom numba ones, prevents this function from
        #   being njited.
        tidal_y_at_mode = \
            tidal_y_solver(
                interior_model_name, radius_array, complex_shears_at_mode, bulk_array, density_array, gravity_array,
                mode_frequency,
                is_solid_by_layer, is_static_by_layer, indices_by_layer,
                order_l, surface_boundary_conditions,
                solve_load_numbers, use_kamata,
                int_rtol, int_atol, rk_method, verbose,
                nondimensionalize, planet_bulk_density
                )

        tidal_potential, \
        tidal_potential_partial_theta, tidal_potential_partial_phi, \
        tidal_potential_partial2_theta2, tidal_potential_partial2_phi2, \
        tidal_potential_partial2_theta_phi = \
            tidal_potential_tuple

        # Calculate stresses and heating
        strains_at_mode, stresses_at_mode = calculate_strain_stress(
            tidal_potential,
            tidal_potential_partial_theta, tidal_potential_partial_phi,
            tidal_potential_partial2_theta2, tidal_potential_partial2_phi2,
            tidal_potential_partial2_theta_phi,
            tidal_y_at_mode,
            colatitude_matrix, radius_array, shear_array, bulk_array,
            mode_frequency, order_l
            )

    return mode_skipped, strains_at_mode, stresses_at_mode, complex_shears_at_mode, tidal_y_at_mode


@njit(cacheable=True)
def collapse_multilayer_modes(
    interior_model_name: str,
    orbital_frequency: float, spin_frequency: float, semi_major_axis: float,
    eccentricity: float, host_mass: float,
    radius_array: np.ndarray, shear_array: np.ndarray, bulk_array: np.ndarray, viscosity_array: np.ndarray,
    density_array: np.ndarray, gravity_array: np.ndarray,
    longitude_matrix: np.ndarray, colatitude_matrix: np.ndarray, time_matrix: np.ndarray, voxel_volume: np.ndarray,
    complex_compliance_function: callable,
    is_solid_by_layer: List[bool], is_static_by_layer: List[bool], indices_by_layer: List[np.ndarray],
    obliquity: float = None,
    surface_boundary_conditions: np.ndarray = None, solve_load_numbers: bool = False,
    complex_compliance_input: Tuple[float, ...] = tuple(), force_mode_calculation: bool = False,
    order_l: int = 2,
    use_modes: bool = True, use_static_potential: bool = False, use_simple_potential: bool = False,
    orbit_average_results: bool = True, use_kamata: bool = False,
    int_rtol: float = 1.0e-8, int_atol: float = 1.0e-12, rk_method: int = 1,
    verbose: bool = False, nondimensionalize: bool = True, planet_bulk_density: float = None
    ):
    """ Calculate the multilayer tidal response of a planet over a range of applicable tidal modes. Collapse
    individual modal results into final heating distribution.

    This is a numba-safe version of the function found in multilayer_modes.py

    Some of these variables are ndarrays the shape of the planet's radius array, N_r
    Others are multidimensional arrays the shape of (N_long, N_colat, N_time).
    These multidim arrays MUST be in this order [longitude_domain, colatitude_domain, time_domain]

    Parameters
    ----------
    interior_model_name : str
        Interior model used in the calculation of the radial functions.
        See options in TidalPy.tides.multilayer.numerical_int.collapse.
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
    radius_array : np.ndarray
        Radius array defined throughout the planet [m]
        Shape: (n_r,) Type: float
    shear_array : np.ndarray
        Static shear modulus array defined throughout the planet at each radius r [Pa]
        Shape: (n_r,) Type: float
    bulk_array : np.ndarray
        Static bulk modulus array defined throughout the planet at each radius r [Pa]
        Shape: (n_r,) Type: float
        # TODO: Bulk dissipation has not been implemented yet.
    viscosity_array : np.ndarray
        Viscosity array defined throughout the planet at each radius r [Pa s]
        Shape: (n_r,) Type: float
    density_array : np.ndarray
        Density array defined throughout the planet at each radius r [kg m-3]
        Shape: (n_r,) Type: float
    gravity_array : np.ndarray
        Acceleration due to gravity defined throughout the planet at the top of each radius r [m s-2]
        Shape: (n_r,) Type: float
    longitude_matrix : np.ndarray
        Ndarray of longitude defined at longitude, colatitude, and time [rad]
        Shape: (n_long, n_colat, n_time) Type: float
    colatitude_matrix : np.ndarray
        Ndarray of colatitude defined at longitude, colatitude, and time [rad]
        Shape: (n_long, n_colat, n_time) Type: float
    time_matrix : np.ndarray
        Ndarray of time defined at longitude, colatitude, and time [s]
        Shape: (n_long, n_colat, n_time) Type: float
    voxel_volume : np.ndarray
        Ndarray of volume defined at each voxel defined at each radius, longitude, colatitude [m-3]
        Shape: (n_r, n_long, n_colat) Type: float
    complex_compliance_function : callable
        Complex compliance function used to find the complex shear modulus.
    is_solid_by_layer : List[bool]
        List or tuple of flags indicating which layers are solid (True) or liquid (False).
    is_static_by_layer : List[bool]
        List or tuple of flags indicating which layers are static (True) or dynamic (False).
    indices_by_layer : List[np.ndarray]
        List or tuple of np.ndarrays flagging which radius indices below to which layer.
        Each away: Shape: (n_r,) Type: np.boolean
    obliquity : float = None
        If not None then the tidal potential that accounts for obliquity will be used.
        Obliquity is relative to the orbital plane [rad]
    surface_boundary_conditions : np.ndarray = None
        Surface conditions applied to the radial function solutions at the top of the planet.
        If equal to None then the function will calculate the surface boundary conditions of either the tidal or load
        number.
    solve_load_numbers : bool = False
        If True, then the surface conditions will be calculated to find the load numbers rather than the tidal ones.
    complex_compliance_input  : Tuple[float, ...] = None
        Optional additional inputs to the complex compliance function.
    force_mode_calculation : bool = False
        If True, then this mode's tidal response will be calculated even if its frequency value falls below the
        minimum threshold.
    order_l: int = 2
        Tidal harmonic order.
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
    use_kamata : bool = False
        If True, the Kamata+2015 initial conditions will be used at the base of layer 0.
        Otherwise, the Takeuchi & Saito 1972 initial conditions will be used.
    int_rtol : float = 1.0e-6
        Integration relative error.
    int_atol : float = 1.0e-4
        Integration absolute error.
    rk_method : int = 1
        The type of RK method used for integration
            0 = RK23
            1 = RK45
            # TODO NotImplemented 2 = DOP
    verbose: bool = False
        If True, the function will print some information to console during calculation (may cause a slow down).
    nondimensionalize : bool = False
        If True, integration will use dimensionless variables. These will be converted back before output is given to
        the user.
    planet_bulk_density : float = None
        Must be provided if non_dimensionalize is True. Bulk density of the planet.

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

    # Clean up inputs
    is_solid_by_layer = nbList(is_solid_by_layer)
    is_static_by_layer = nbList(is_static_by_layer)
    indices_by_layer = nbList(indices_by_layer)

    # Certain variables are calculated across the radius, longitude, colatitude, and time domains.
    #   longitude, colatitude, and time are provided as matrices that must be in this order:
    #   [longitude_N, latitude_N, time_N]
    # Check that dimensions make sense
    assert radius_array.shape == shear_array.shape

    # The shear array may have zero values for liquid layers. This will cause an issue with complex compliance calc.
    #     make it small instead
    shear_array[shear_array <= 0.] = 1.e-50

    # Pull out individual arrays
    longitude_domain = longitude_matrix[:, 0, 0]
    colatitude_domain = colatitude_matrix[0, :, 0]
    time_domain = time_matrix[0, 0, :]
    r_shape = radius_array.shape
    colat_shape = colatitude_matrix.shape
    mixed_shape = (*r_shape, *colat_shape)

    # Check that the time domain has the correct end points
    orbital_period = 2. * np.pi / orbital_frequency
    if orbit_average_results:
        # In order for the orbit average routine to work correctly, the time domain must start at zero and end
        #    after 1 orbital period.
        assert time_domain[0] == 0.
        assert time_domain[-1] == np.abs(orbital_period)

    planet_radius = radius_array[-1]

    # TODO: Currently this function (and other multilayer code) does not work for l>2. An additional loop will be
    #   required to loop over each l and sum the results.
    #   Implementation for this function is straight forward, basically another loop around everything from
    #       l = range(2, max_l+1)
    #   The tricker part will be updating the tidal potential which is hardcoded for l=2 at the moment.
    if order_l != 2:
        raise Exception(
            'Not Implemented:'
            'Multilayer tides (specifically the tidal potential) only works for l=2 at the moment.'
            )

    # Setup tidal potential functions
    # Check if obliquity is being used.
    # Calculate the tidal modes and the tidal potential and its partial derivatives.
    # Tidal potential is only calculated at the surface of the planet.
    if obliquity is not None:
        if use_simple_potential:
            raise Exception('The simple version of the tidal potential does not account for a non-zero obliquity')
        elif use_modes:
            potential_output = tidal_potential_obliquity_nsr_modes(
                planet_radius, longitude_matrix, colatitude_matrix, time_matrix, orbital_frequency,
                spin_frequency, eccentricity, obliquity, host_mass, semi_major_axis, use_static_potential
                )
        else:
            potential_output = tidal_potential_obliquity_nsr(
                planet_radius, longitude_matrix, colatitude_matrix, time_matrix, orbital_frequency,
                spin_frequency, eccentricity, obliquity, host_mass, semi_major_axis, use_static_potential
                )
    else:
        # Obliquity is not used. pick the appropriate potentials
        if use_simple_potential:
            # Simple potential assumes very low eccentricity, no obliquity, and synchronous rotation.
            potential_output = tidal_potential_simple(
                planet_radius, longitude_matrix, colatitude_matrix, time_matrix, orbital_frequency,
                eccentricity, host_mass, semi_major_axis
                )
        elif use_modes:
            potential_output = tidal_potential_nsr_modes(
                planet_radius, longitude_matrix, colatitude_matrix, time_matrix, orbital_frequency,
                spin_frequency,
                eccentricity, host_mass, semi_major_axis, use_static_potential
                )
        else:
            potential_output = tidal_potential_nsr(
                planet_radius, longitude_matrix, colatitude_matrix, time_matrix, orbital_frequency,
                spin_frequency,
                eccentricity, host_mass, semi_major_axis, use_static_potential
                )

    # Pull out results from potential calculation
    tidal_frequencies, tidal_modes, tidal_potential_tuple_by_mode = potential_output

    # Record how many modes are skipped (used in average)
    num_modes_skipped = 0

    # Build storages for all modes
    modes_skipped = nbDict.empty(
        key_type=nbUnicode,
        value_type=bool_
        )
    love_k_by_mode = nbDict.empty(
        key_type=nbUnicode,
        value_type=complex128
        )
    love_h_by_mode = nbDict.empty(
        key_type=nbUnicode,
        value_type=complex128
        )
    love_l_by_mode = nbDict.empty(
        key_type=nbUnicode,
        value_type=complex128
        )

    # Large arrays must be added to continuously to avoid memory overload
    complex_shears_avg = np.zeros(r_shape, dtype=np.complex128)
    tidal_y_avg = np.zeros((6, *r_shape), dtype=np.complex128)
    stresses = np.zeros((6, *mixed_shape), dtype=np.complex128)
    strains = np.zeros((6, *mixed_shape), dtype=np.complex128)
    tidal_potential = np.zeros(colat_shape, dtype=np.complex128)
    total_potential = np.zeros(mixed_shape, dtype=np.complex128)

    # Opt: multiprocessor could speed this up when there are several modes to calculate. The continuous summation
    #   would just need to be done at the end. The issue may be that several large arrays would have to be stored in
    #   memory at once. This could negate some/all of the benefits of multiprocessing, especially on machines with
    #   a low amount of ram.
    for mode_name, mode_frequency in tidal_frequencies.items():

        tidal_potential_tuple = tidal_potential_tuple_by_mode[mode_name]
        tidal_potential_at_mode = tidal_potential_tuple[0]

        # Calculate response at mode
        mode_skipped, strains_at_mode, stresses_at_mode, complex_shears_at_mode, tidal_y_at_mode = \
            calculate_mode_response_coupled(
                interior_model_name, mode_frequency,
                radius_array, shear_array, bulk_array, viscosity_array,
                density_array, gravity_array, colatitude_matrix,
                tidal_potential_tuple, complex_compliance_function,
                is_solid_by_layer, is_static_by_layer, indices_by_layer,
                surface_boundary_conditions, solve_load_numbers,
                complex_compliance_input, force_mode_calculation,
                order_l, use_kamata,
                int_rtol, int_atol, rk_method, verbose,
                nondimensionalize, planet_bulk_density
                )

        if mode_skipped:
            num_modes_skipped += 1
            modes_skipped[mode_name] = mode_skipped

        # Collapse Modes
        # Add items defined for each mode to lists
        # The max here is to not lead to negative k2s when y5 = 0 (when frequency = 0)
        love_k = tidal_y_at_mode[4, -1] - 1.
        love_k_by_mode[mode_name] = max(0., np.real(love_k)) + (1.j * np.imag(love_k))
        love_h_by_mode[mode_name] = tidal_y_at_mode[0, -1] * gravity_array[-1]
        love_l_by_mode[mode_name] = tidal_y_at_mode[2, -1] * gravity_array[-1]

        if not mode_skipped:

            # Stresses and strains used in heat calculation are added together.
            # TODO: should this scale be w or w/2pi or w/2? w/2 seems to give the best comparison to homogen equation.
            freq_half = mode_frequency / 2.
            stresses_at_mode *= freq_half
            strains_at_mode *= freq_half

            # Estimate the "orbit" averaged response.
            # TODO: I believe to calculate the orbit averaged response then we need to find Int_0^T 1/T f(t)dt
            #   However, using the multi-mode approach, there are many different periods, T. So we will multiple all
            #   relevant functions by 1/T_i where T_i is the period at this mode. Proceed with the summation and then
            #   at the very end perform the integration for the average.
            # Now we can scale other values by the mode frequency
            if orbit_average_results and (mode_frequency > 1.0e-10):
                period_inv = freq_half / np.pi
                tidal_potential_at_mode *= period_inv
                # # TODO: The following two optimizations are not scaled by the same frequency. w/2 instead of w/2pi.
                # #     How does that affect things? Is it worth the optimization?
                # stresses_at_mode = stress_scaled_at_mode
                # strains_at_mode = strain_scaled_at_mode

        # Stresses, strains, and potentials are added together.
        # TODO: I was tracking two versions of stress and two of strain, one that got this multiplier and one that
        #   didn't. but these arrays are large and really impact performance. We will try to track just one
        #   and see how the results look.
        stresses += stresses_at_mode
        strains += strains_at_mode

        tidal_potential += tidal_potential_at_mode
        total_potential_at_mode = np.expand_dims(tidal_potential_at_mode, axis=0) * \
                                  np.expand_dims(
                                      np.expand_dims(
                                          np.expand_dims(
                                              tidal_y_at_mode[4, :],
                                              axis=-1
                                              ),
                                          axis=-1
                                          ),
                                      axis=-1
                                      )

        total_potential += np.ascontiguousarray(total_potential_at_mode)

        # The other parameters it is not clear what category they fall into.
        # TODO: For now let's take the average of them. Should they also receive a similar 1/T treatment like the above
        #     properties?
        complex_shears_avg += complex_shears_at_mode
        tidal_y_avg += tidal_y_at_mode

    # Finish taking the average of the avg parameters
    # TODO: see previous todo
    complex_shears_avg = complex_shears_avg / max((len(tidal_modes) - num_modes_skipped), 1)
    tidal_y_avg = tidal_y_avg / max((len(tidal_modes) - num_modes_skipped), 1)

    # Calculate tidal heating
    # This stress/strain term is used in the final calculation of tidal heating. It tracks a different coefficient
    #   for frequency averaging (if stress and strain are used on their own then heating would be prop to T^{-2}
    #   rather than T^{-1}
    # Heating is equal to imag[o] * real[s] - real[o] * imag[s] but we need to multiply by two for the cross terms
    #    since it is part of a symmetric matrix but only one side of the matrix is calculated in the previous steps.
    # First calculate the trace terms which are not multiply by two.
    volumetric_heating = (
        # Im[\sigma_ij] * Re[\epsilon_ij]
            (
                    np.sum(np.imag(stresses[:3]) * np.real(strains[:3]), axis=0) +
                    # Now add the cross terms where we do multiply by two
                    2. * np.sum(np.imag(stresses[:3]) * np.real(strains[:3]), axis=0)
            ) -
            # minus Re[\sigma_ij] * Im[\epsilon_ij]
            (
                    np.sum(np.real(stresses[3:]) * np.imag(strains[3:]), axis=0) +
                    2. * np.sum(np.real(stresses[3:]) * np.imag(strains[3:]), axis=0)
            )
    )

    # TODO: Without this abs term the resulting heating maps are very blotchy around
    #    Europa book does have an abs at Equation 42, Page 102
    volumetric_heating = np.abs(volumetric_heating)

    # Perform orbital averaging
    if orbit_average_results:
        # These trapz functions should have an axis=-1 as an additional argument but numba does not support it
        #    and by default uses the last axis anyways so we should be fine.
        strains_int = np.trapz(strains, time_domain)
        stresses_int = np.trapz(stresses, time_domain)
        tidal_potential_int = np.trapz(tidal_potential, time_domain)
        total_potential_int = np.trapz(total_potential, time_domain)
        volumetric_heating_int = np.trapz(volumetric_heating, time_domain)

        # Numba does not like that orbit_average_results can change the dimensions of the returned values. So even
        #   though we just collapsed a dimension with the trapz functions, let's add it back.
        strains = np.expand_dims(strains_int, axis=-1)
        stresses = np.expand_dims(stresses_int, axis=-1)
        tidal_potential = np.expand_dims(tidal_potential_int, axis=-1)
        total_potential = np.expand_dims(total_potential_int, axis=-1)
        volumetric_heating = np.expand_dims(volumetric_heating_int, axis=-1)

        # Now lets make sure that these arrays are ordered correctly
        strains = np.ascontiguousarray(strains)
        stresses = np.ascontiguousarray(stresses)
        tidal_potential = np.ascontiguousarray(tidal_potential)
        total_potential = np.ascontiguousarray(total_potential)
        volumetric_heating = np.ascontiguousarray(volumetric_heating)

    # To find the total heating (rather than volumetric) we need to multiply by the volume in each voxel.
    # The voxel_volume has a shape of (r_N, long_N, colat_N).
    voxel_volume_higher_dim = np.expand_dims(voxel_volume, axis=-1)
    heating = volumetric_heating * voxel_volume_higher_dim

    # Put the love results into a tuple
    love_results = (love_k_by_mode, love_h_by_mode, love_l_by_mode)

    return heating, volumetric_heating, strains, stresses, total_potential, tidal_potential, \
           complex_shears_avg, tidal_y_avg, love_results, tidal_modes, modes_skipped
