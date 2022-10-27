""" Function to collapse multilayer solutions for multiple tidal modes. """

# If no inputs to the complex compliance function were provided then set it equal to an empty tuple which
#    will cause the complex compliance function to resort to default.
if complex_compliance_input is None:
    complex_compliance_input = tuple()

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
modes_skipped = dict()
love_k_by_mode = dict()
love_h_by_mode = dict()
love_l_by_mode = dict()

# Large arrays must be added to continuously to avoid memory overload
complex_shears_avg = np.zeros(r_shape, dtype=np.complex128)
tidal_y_avg = np.zeros((6, *r_shape), dtype=np.complex128)
stresses = np.zeros((6, *mixed_shape), dtype=np.complex128)
strains = np.zeros((6, *mixed_shape), dtype=np.complex128)
tidal_potential = np.zeros(colat_shape, dtype=np.complex128)
total_potential = np.zeros(mixed_shape, dtype=np.complex128)


from typing import Dict
import numpy as np


MIN_MODE_FREQ = 10. * np.finfo(np.float64).eps


def collapse_multilayer_modes(tidal_modes: Dict[str, float], collapse_like_modes: bool = False):

    # Opt: multiprocessor could speed this up when there are several modes to calculate. The continuous summation
    #   would just need to be done at the end. The issue may be that several large arrays would have to be stored in
    #   memory at once. This could negate some/all of the benefits of multiprocessing, especially on machines with
    #   a low amount of ram.

    # Pre-build storage to allow for numba
    modes_skipped = {'x': False}
    collapsed_modes = {'x': ['x']}
    delete_collapse_mode_fake_value = False

    if collapse_like_modes:
        tidal_modes_to_use = {'x': 0.}
        modes_to_skip = []
        modes_done = []
        # Loop through the provided tidal modes twice, merging like-valued modes.
        for mode, value in tidal_modes.items():
            if mode in modes_to_skip:
                continue
            tidal_modes_to_use[mode] = value
            modes_done.append(mode)
            for mode_check, value_check in tidal_modes.items():
                if mode_check in modes_done:
                    continue
                if abs(value - value_check) < MIN_MODE_FREQ:
                    # These two modes are very similar in value, merge them together
                    modes_to_skip.append(mode_check)



    else:
        # Otherwise use the tidal modes provided by the user
        tidal_modes_to_use = tidal_modes

    if delete_collapse_mode_fake_value:
        del tidal_modes_to_use['x']

    for mode_name, tidal_mode in tidal_modes.items():

        if np.abs(tidal_mode) < MIN_MODE_FREQ:
            modes_skipped[mode_name] = True

        else:
            modes_skipped[mode_name] = False

    # Delete fake values used for dictionary building
    del modes_skipped['x']
    del collapsed_modes['x']



    love_h_by_mode
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
                surface_boundary_conditions=surface_boundary_conditions, solve_load_numbers=solve_load_numbers,
                complex_compliance_input=complex_compliance_input, force_mode_calculation=force_mode_calculation,
                order_l=order_l, use_kamata=use_kamata,
                integrator=integrator, integration_method=integration_method,
                integration_rtol=integration_rtol, integration_atol=integration_atol,
                verbose=verbose, nondimensionalize=nondimensionalize, planet_bulk_density=planet_bulk_density,
                incompressible=incompressible
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
    # Heating is equal to imag[o] * real[s] - real[o] * imag[s] but we need to multiply by two for the cross terms
    #    since it is part of a symmetric matrix but only one side of the matrix is calculated in the previous steps.
    volumetric_heating = (
        # Im[s_rr] Re[e_rr] - Re[s_rr] Im[e_rr]
            np.imag(stresses[0]) * np.real(strains[0]) - np.real(stresses[0]) * np.imag(strains[0]) +
            # Im[s_thth] Re[e_thth] - Re[s_thth] Im[e_thth]
            np.imag(stresses[1]) * np.real(strains[1]) - np.real(stresses[1]) * np.imag(strains[1]) +
            # Im[s_phiphi] Re[e_phiphi] - Re[s_phiphi] Im[e_phiphi]
            np.imag(stresses[2]) * np.real(strains[2]) - np.real(stresses[2]) * np.imag(strains[2]) +
            # Im[s_rth] Re[e_rth] - Re[s_rth] Im[e_rth]
            2. * (np.imag(stresses[3]) * np.real(strains[3]) - np.real(stresses[3]) * np.imag(strains[3])) +
            # Im[s_rphi] Re[e_rphi] - Re[s_rphi] Im[e_rphi]
            2. * (np.imag(stresses[4]) * np.real(strains[4]) - np.real(stresses[4]) * np.imag(strains[4])) +
            # Im[s_thphi] Re[e_thphi] - Re[s_thphi] Im[e_thphi]
            2. * (np.imag(stresses[5]) * np.real(strains[5]) - np.real(stresses[5]) * np.imag(strains[5]))
    )

    # TODO: Without this abs term the resulting heating maps are very blotchy around
    #    Europa book does have an abs at Equation 42, Page 102
    volumetric_heating = np.abs(volumetric_heating)

    # Perform orbital averaging
    if orbit_average_results:
        strains = np.trapz(strains, time_domain, axis=-1)
        stresses = np.trapz(stresses, time_domain, axis=-1)
        tidal_potential = np.trapz(tidal_potential, time_domain, axis=-1)
        total_potential = np.trapz(total_potential, time_domain, axis=-1)
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