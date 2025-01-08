""" Multilayer calculator for multiple tidal modes.

Each tidal mode imparts a potentially different frequency which needs to be propagated throughout the planet's interior.

This module contains functions to assist with calculating the response at each of these frequencies and then collapse
    the findings into a final value.

"""
from typing import Dict, List, Tuple, TYPE_CHECKING

import numpy as np

from TidalPy.utilities.performance import nbList
from TidalPy.rheology.base import RheologyModelBase
from TidalPy.RadialSolver import radial_solver

from ..multilayer.stress_strain import calculate_strain_stress
from ..potential import (tidal_potential_nsr, tidal_potential_nsr_modes,
                         tidal_potential_gen_obliquity_nsr_modes, tidal_potential_gen_obliquity_nsr,
                         tidal_potential_obliquity_nsr, tidal_potential_obliquity_nsr_modes, tidal_potential_simple)

from ..heating import calculate_volumetric_heating

if TYPE_CHECKING:
    from ..potential import TidalPotentialOutput


def calculate_mode_response_coupled(
        mode_frequency: float,
        radius_array: np.ndarray,
        density_array: np.ndarray,
        shear_array: np.ndarray,
        bulk_array: np.ndarray,
        shear_viscosity_array: np.ndarray,
        bulk_viscosity_array: np.ndarray,
        upper_radius_bylayer_array: np.ndarray,
        longitude_domain: np.ndarray,
        colatitude_domain: np.ndarray,
        time_domain: np.ndarray,
        tidal_potential_tuple: 'TidalPotentialOutput',
        shear_rheology_inst: 'RheologyModelBase',
        bulk_rheology_inst: 'RheologyModelBase',
        layer_types: Tuple[str],
        is_static_bylayer: Tuple[bool],
        is_incompressible_bylayer: Tuple[bool],
        solve_load_numbers: bool = False,
        force_mode_calculation: bool = False,
        degree_l: int = 2,
        use_kamata: bool = False,
        integration_method: str = 'DOP853',
        integration_rtol: float = 1.0e-6,
        integration_atol: float = 1.0e-10,
        verbose: bool = False,
        nondimensionalize: bool = True,
        planet_bulk_density: float = None
        ):
    """ Given a tidal frequency, this function will call on the interior integration routine with the proper inputs and
        collect the results as well as calculate tidal stress, strain, and heating as a
        function of radius, longitude, latitude, and time.

    TODO: update input param docs

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
    # Get size of arrays
    n_radius     = len(radius_array)
    n_longitude  = len(longitude_domain)
    n_colatitude = len(colatitude_domain)
    n_time       = len(time_domain)

    # Setup flags
    mode_skipped = False
    
    if solve_load_numbers:
        solve_for_tuple = ('tidal', 'loading')
    else:
        solve_for_tuple = ('tidal',)

    # Calculate rheology and radial response
    if (not force_mode_calculation) and (mode_frequency < 1.0e-15):
        # If frequency is ~ 0.0 then there will be no tidal response. Skip the calculation of tidal y, etc.
        tidal_y_at_mode = np.zeros((6, n_radius), dtype=np.complex128)
        love_numbers = (0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j)
        complex_shears_at_mode = shear_array + 0.0j
        strains_at_mode = np.zeros((6, n_radius, n_longitude, n_colatitude, n_time), dtype=np.complex128)
        stresses_at_mode = np.zeros((6, n_radius, n_longitude, n_colatitude, n_time), dtype=np.complex128)
        mode_skipped = True

    else:
        # Calculate Complex Shear Modulus at this mode's frequency
        complex_shears_at_mode = np.empty(radius_array.size, dtype=np.complex128)

        shear_rheology_inst.vectorize_modulus_viscosity(
            mode_frequency,
            shear_array,
            shear_viscosity_array,
            complex_shears_at_mode)

        # Calculate Complex Bulk Modulus at this mode's frequency
        complex_bulk_modulus = np.empty(radius_array.size, dtype=np.complex128)

        bulk_rheology_inst.vectorize_modulus_viscosity(
            mode_frequency,
            bulk_array,
            bulk_viscosity_array,
            complex_bulk_modulus)
        
        # Calculate the radial functions using a shooting integration method.
        radial_solution_at_mode = \
            radial_solver(
                radius_array,
                density_array,
                complex_bulk_modulus,
                complex_shears_at_mode,
                mode_frequency,
                planet_bulk_density,
                layer_types,
                is_static_bylayer,
                is_incompressible_bylayer,
                upper_radius_bylayer_array,
                degree_l=degree_l,
                solve_for=solve_for_tuple,
                core_model=0,
                use_kamata=use_kamata,
                starting_radius = 0.0,
                start_radius_tolerance = 1.0e-5,
                integration_method=integration_method,
                integration_rtol=integration_rtol,
                integration_atol=integration_atol,
                scale_rtols_bylayer_type = False,
                max_num_steps = 50_000_000,
                expected_size = 500,
                max_ram_MB = 2500,
                max_step = 0,
                nondimensionalize=nondimensionalize,
                use_prop_matrix = False,
                eos_method_bylayer = None,
                eos_integration_method = 'RK45',
                eos_rtol = 1.0e-3,
                eos_atol = 1.0e-5,
                eos_pressure_tol = 1.0e-3,
                eos_max_iters = 40,
                verbose = verbose,
                warnings = verbose,
                raise_on_fail = False,
                perform_checks = True)
        if not radial_solution_at_mode.success:
            raise Exception(f"RadialSolver failed to find a solution at this mode: {radial_solution_at_mode.message}")

        tidal_y_at_mode = radial_solution_at_mode.result
        if solve_load_numbers:
            love_numbers = (radial_solution_at_mode.k[0], radial_solution_at_mode.h[0], radial_solution_at_mode.l[0])
        else:
            love_numbers = (radial_solution_at_mode.k, radial_solution_at_mode.h, radial_solution_at_mode.l)

        # Calculate stresses and heating
        strains_at_mode, stresses_at_mode = calculate_strain_stress(
            *tidal_potential_tuple,
            tidal_y_at_mode,
            longitude_domain, colatitude_domain, time_domain, radius_array, complex_shears_at_mode, complex_bulk_modulus,
            frequency=mode_frequency, degree_l=degree_l
            )

    return mode_skipped, strains_at_mode, stresses_at_mode, complex_shears_at_mode, tidal_y_at_mode, love_numbers


def collapse_multilayer_modes(
    orbital_frequency: float,
    spin_frequency: float,
    semi_major_axis: float,
    eccentricity: float,
    host_mass: float,
    planet_bulk_density: float,
    radius_array: np.ndarray,
    density_array: np.ndarray,
    bulk_array: np.ndarray,
    shear_array: np.ndarray,
    bulk_viscosity_array: np.ndarray,
    shear_viscosity_array: np.ndarray,
    shear_rheology_inst: 'RheologyModelBase',
    bulk_rheology_inst: 'RheologyModelBase',
    upper_radius_bylayer_array: np.ndarray,
    longitude_matrix: np.ndarray,
    colatitude_matrix: np.ndarray,
    time_matrix: np.ndarray,
    voxel_volume: np.ndarray,
    layer_types: tuple[str],
    is_static_bylayer: tuple[bool],
    is_incompressible_bylayer: tuple[bool],
    obliquity: float = None,
    use_general_obliquity: bool = False,
    solve_load_numbers: bool = False,
    force_mode_calculation: bool = False,
    degree_l: int = 2,
    use_modes: bool = True,
    use_static_potential: bool = False,
    use_simple_potential: bool = False,
    orbit_average_results: bool = True,
    use_kamata: bool = False,
    integration_method: str = 'DOP853',
    integration_rtol: float = 1.0e-6,
    integration_atol: float = 1.0e-8,
    verbose: bool = False,
    nondimensionalize: bool = True
    ):
    """ Calculate the multilayer tidal response of a planet over a range of applicable tidal modes. Collapse
    individual modal results into final heating distribution.

    Some of these variables are ndarrays the shape of the planet's radius array, N_r
    Others are multidimensional arrays the shape of (N_long, N_colat, N_time).
    These multidim arrays MUST be in this order [longitude_domain, colatitude_domain, time_domain]

    Parameters
    ----------
    TODO: Update input parameter docs

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

    # Certain variables are calculated across the radius, longitude, colatitude, and time domains.
    #   longitude, colatitude, and time are provided as matrices that must be in this order:
    #   [longitude_N, latitude_N, time_N]
    # Check that dimensions make sense
    assert radius_array.shape == shear_array.shape

    # Cleanup inputs
    layer_types               = tuple(layer_types)
    is_static_bylayer         = tuple(is_static_bylayer)
    is_incompressible_bylayer = tuple(is_incompressible_bylayer)

    # The shear array may have zero values for liquid layers. This will cause an issue with complex compliance calc.
    #     make it small instead
    shear_array[shear_array < 0.] = 0.0

    # Pull out individual arrays
    longitude_domain  = longitude_matrix[:, 0, 0]
    colatitude_domain = colatitude_matrix[0, :, 0]
    time_domain       = time_matrix[0, 0, :]
    r_shape           = radius_array.shape
    colat_shape       = colatitude_matrix.shape
    mixed_shape       = (*r_shape, *colat_shape)

    # Check that the time domain has the correct end points
    orbital_period = 2. * np.pi / orbital_frequency
    if orbit_average_results:
        # In order for the orbit average routine to work correctly, the time domain must start at zero and end
        #    after 1 orbital period.
        assert time_domain[0] == 0.
        assert np.isclose(time_domain[-1], np.abs(orbital_period))

    planet_radius = radius_array[-1]

    # TODO: Currently this function (and other multilayer code) does not work for l>2. An additional loop will be
    #   required to loop over each l and sum the results.
    #   Implementation for this function is straight forward, basically another loop around everything from
    #       l = range(2, max_l+1)
    #   The tricker part will be updating the tidal potential which is hardcoded for l=2 at the moment.
    if degree_l != 2:
        raise NotImplementedError(
            'Multilayer tides (specifically the tidal potential) only works for '
            'l=2 at the moment.'
            )

    # Setup tidal potential functions
    # Check if obliquity is being used.
    # Calculate the tidal modes and the tidal potential and its partial derivatives.
    # Tidal potential is only calculated at the surface of the planet.
    if obliquity is not None:
        if use_simple_potential:
            raise Exception('The simple version of the tidal potential does not account for a non-zero obliquity')
        elif use_modes:
            if use_general_obliquity:
                potential_output = tidal_potential_gen_obliquity_nsr_modes(
                    planet_radius, longitude_matrix, colatitude_matrix, time_matrix, orbital_frequency,
                    spin_frequency, eccentricity, obliquity, host_mass, semi_major_axis, use_static_potential
                    )
            else:
                potential_output = tidal_potential_obliquity_nsr_modes(
                    planet_radius, longitude_matrix, colatitude_matrix, time_matrix, orbital_frequency,
                    spin_frequency, eccentricity, obliquity, host_mass, semi_major_axis, use_static_potential
                    )
        else:
            if use_general_obliquity:
                potential_output = tidal_potential_gen_obliquity_nsr(
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
    modes_skipped  = dict()
    love_k_by_mode = dict()
    love_h_by_mode = dict()
    love_l_by_mode = dict()

    # Large arrays must be added to continuously to avoid memory overload
    complex_shears_avg = np.zeros(r_shape, dtype=np.complex128)
    tidal_y_avg        = np.zeros((6, *r_shape), dtype=np.complex128)
    stresses           = np.zeros((6, *mixed_shape), dtype=np.complex128)
    strains            = np.zeros((6, *mixed_shape), dtype=np.complex128)
    tidal_potential    = np.zeros(colat_shape, dtype=np.complex128)
    total_potential    = np.zeros(mixed_shape, dtype=np.complex128)

    # Opt: multiprocessor could speed this up when there are several modes to calculate. The continuous summation
    #   would just need to be done at the end. The issue may be that several large arrays would have to be stored in
    #   memory at once. This could negate some/all of the benefits of multiprocessing, especially on machines with
    #   a low amount of ram.
    for mode_name, mode_frequency in tidal_frequencies.items():

        tidal_potential_tuple   = tidal_potential_tuple_by_mode[mode_name]
        tidal_potential_at_mode = tidal_potential_tuple[0]

        # Calculate response at mode
        mode_skipped, strains_at_mode, stresses_at_mode, complex_shears_at_mode, tidal_y_at_mode, love_numbers = \
            calculate_mode_response_coupled(
                mode_frequency,
                radius_array,
                density_array,
                shear_array,
                bulk_array,
                shear_viscosity_array,
                bulk_viscosity_array,
                upper_radius_bylayer_array,
                longitude_domain,
                colatitude_domain,
                time_domain,
                tidal_potential_tuple,
                shear_rheology_inst,
                bulk_rheology_inst,
                layer_types,
                is_static_bylayer,
                is_incompressible_bylayer,
                solve_load_numbers,
                force_mode_calculation,
                degree_l,
                use_kamata,
                integration_method,
                integration_rtol,
                integration_atol,
                verbose,
                nondimensionalize,
                planet_bulk_density
                )

        if mode_skipped:
            num_modes_skipped += 1
            modes_skipped[mode_name] = mode_skipped

        # Collapse Modes
        # Add items defined for each mode to lists
        love_k_by_mode[mode_name] = love_numbers[0]
        love_h_by_mode[mode_name] = love_numbers[1]
        love_l_by_mode[mode_name] = love_numbers[2]
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
        total_potential += tidal_potential_at_mode[np.newaxis, :, :, :] * \
                           tidal_y_at_mode[4, :, np.newaxis, np.newaxis, np.newaxis]

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
    volumetric_heating = calculate_volumetric_heating(stresses, strains)

    # Perform orbital averaging
    if orbit_average_results:
        strains            = np.trapz(strains, time_domain, axis=-1)
        stresses           = np.trapz(stresses, time_domain, axis=-1)
        tidal_potential    = np.trapz(tidal_potential, time_domain, axis=-1)
        total_potential    = np.trapz(total_potential, time_domain, axis=-1)
        volumetric_heating = np.trapz(volumetric_heating, time_domain, axis=-1)

    # To find the total heating (rather than volumetric) we need to multiply by the volume in each voxel.
    # The voxel_volume has a shape of (r_N, long_N, colat_N).
    if orbit_average_results:
        heating = volumetric_heating * voxel_volume
    else:
        # If we do not orbit average then we need to expand the voxel volume dimensions to allow for
        #   ndarray multiplication
        voxel_volume_higher_dim = voxel_volume[:, :, :, np.newaxis]
        heating = volumetric_heating * voxel_volume_higher_dim

    # Put the love results into a tuple
    love_results = (love_k_by_mode, love_h_by_mode, love_l_by_mode)

    return heating, volumetric_heating, strains, stresses, total_potential, tidal_potential, \
           complex_shears_avg, tidal_y_avg, love_results, tidal_modes, modes_skipped
