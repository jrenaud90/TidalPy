""" Functionality to non-dimensionalize common variables used for multi-layer tidal calculations.

The scheme is based on that proposed by Martens16 (around page 99)

References
----------
Martens16 : H. Martens, PhD Thesis (CalTech), 2016, DOI: 10.7907/Z9N29TX7

"""

from scipy.constants import G as G_
from libc.math cimport pi, sqrt

cdef double G = G_


cdef void cf_non_dimensionalize_physicals(
        size_t num_radius,
        double frequency,
        double mean_radius,
        double bulk_density,
        double* radius_array_ptr,
        double* density_array_ptr,
        double* gravity_array_ptr,
        double* bulk_array_ptr,
        double_numeric* shear_array_ptr,
        double* frequency_to_use,
        double* G_to_use
        ) noexcept nogil:

    # Setup loop variables
    cdef size_t i

    # Setup conversions
    cdef double second2_conversion, second_conversion, length_conversion
    cdef double density_conversion, mass_conversion, pascal_conversion
    second2_conversion = 1. / (pi * G * bulk_density)
    second_conversion  = sqrt(second2_conversion)
    length_conversion  = mean_radius
    density_conversion = bulk_density
    mass_conversion    = bulk_density * mean_radius**3
    pascal_conversion  = mass_conversion / (length_conversion * second2_conversion)

    # Convert array pointers
    for i in range(num_radius):
        radius_array_ptr[i]    /= length_conversion
        density_array_ptr[i]   /= density_conversion
        gravity_array_ptr[i]   /= (length_conversion / second2_conversion)
        bulk_array_ptr[i]      /= pascal_conversion
        shear_array_ptr[i]     /= pascal_conversion

    # Convert non-array pointers
    G_to_use[0]         = G / (length_conversion**3 / (mass_conversion * second2_conversion))
    frequency_to_use[0] = frequency / (1. / second_conversion)

cdef void cf_redimensionalize_physicals(
        size_t num_radius,
        double frequency,
        double mean_radius,
        double bulk_density,
        double* radius_array_ptr,
        double* density_array_ptr,
        double* gravity_array_ptr,
        double* bulk_array_ptr,
        double_numeric* shear_array_ptr,
        double* frequency_to_use,
        double* G_to_use
        ) noexcept nogil:

    # Setup loop variables
    cdef size_t i

    # Setup conversions
    cdef double second2_conversion, second_conversion, length_conversion
    cdef double density_conversion, mass_conversion, pascal_conversion
    second2_conversion = 1. / (pi * G * bulk_density)
    second_conversion  = sqrt(second2_conversion)
    length_conversion  = mean_radius
    density_conversion = bulk_density
    mass_conversion    = bulk_density * mean_radius**3
    pascal_conversion  = mass_conversion / (length_conversion * second2_conversion)

    # Convert array pointers
    for i in range(num_radius):
        radius_array_ptr[i]    *= length_conversion
        density_array_ptr[i]   *= density_conversion
        gravity_array_ptr[i]   *= (length_conversion / second2_conversion)
        bulk_array_ptr[i]      *= pascal_conversion
        shear_array_ptr[i]     *= pascal_conversion

    # Convert non-array pointers
    G_to_use[0]         = G * (length_conversion**3 / (mass_conversion * second2_conversion))
    frequency_to_use[0] = frequency * (1. / second_conversion)



cdef void cf_redimensionalize_radial_functions(
    double complex* radial_function_ptr,
    double mean_radius,
    double bulk_density,
    size_t num_slices,
    size_t num_solutions = 1) noexcept nogil:
    """ A function to re-dimensionalize physical parameters that have been previously non-dimensionalized.

    Parameters
    ----------
    radial_function_ptr : complex128*
        Non-dimensionalized radial solutions as a function of radius.
    mean_radius : float64
        Mean radius of the planet, used in scaling [m]
    bulk_density : float64
        Bulk density of the planet, used in scaling [m]
    num_slices : uint32
        Number of radial slices, used for looping.
    num_solutions : uint32, default=1
        Number of solutions to loop through (size of radial_function_ptr is 6 * num_solutions * num_slices)

    """
    # Loop variables
    cdef size_t slice_i, solver_i
    # Setup conversions
    cdef double second2_conversion = 1. / (pi * G * bulk_density)
    cdef double mass_conversion = bulk_density * mean_radius**3
    cdef double length_conversion = mean_radius

    for solver_i in range(num_solutions):
        for slice_i in range(num_slices):
            # Convert displacements
            #    y1, y3 are the radial and tangential displacements with units of [s2 m-1]
            #    y2, y4 are the radial and tangential stresses with units of [kg m-3]
            #    y5 is the tidal potential which is unitless and thus needs no conversion.
            #    y6 is a "potential stress" with units of [m-1]

            radial_function_ptr[slice_i * 6 * num_solutions + solver_i * 6 + 0] *= \
                (second2_conversion / length_conversion)
            radial_function_ptr[slice_i * 6 * num_solutions + solver_i * 6 + 2] *= \
                (second2_conversion / length_conversion)

            radial_function_ptr[slice_i * 6 * num_solutions + solver_i * 6 + 1] *= \
                (mass_conversion / length_conversion**3)
            radial_function_ptr[slice_i * 6 * num_solutions + solver_i * 6 + 3] *= \
                (mass_conversion / length_conversion**3)

            radial_function_ptr[slice_i * 6 * num_solutions + solver_i * 6 + 4] *= \
                1
            radial_function_ptr[slice_i * 6 * num_solutions + solver_i * 6 + 5] *= \
                (1. / length_conversion)


def non_dimensionalize_physicals(
        double frequency,
        double mean_radius,
        double bulk_density,
        double[:] radius_array_view,
        double[:] density_array_view,
        double[:] gravity_array_view,
        double[:] bulk_array_view,
        double_numeric[:] shear_array_view,
        ):

    cdef size_t num_radius = radius_array_view.size
    cdef double frequency_to_use, G_to_use
    
    cf_non_dimensionalize_physicals(
        num_radius, frequency, mean_radius, bulk_density,
        &radius_array_view[0], &density_array_view[0], &gravity_array_view[0],
        &bulk_array_view[0], &shear_array_view[0],
        &frequency_to_use, &G_to_use
        )
    
    return frequency_to_use, G_to_use

def redimensionalize_physicals(
        double frequency,
        double mean_radius,
        double bulk_density,
        double[:] radius_array_view,
        double[:] density_array_view,
        double[:] gravity_array_view,
        double[:] bulk_array_view,
        double_numeric[:] shear_array_view,
        ):

    cdef size_t num_radius = radius_array_view.size
    cdef double frequency_to_use, G_to_use
    
    cf_redimensionalize_physicals(
        num_radius, frequency, mean_radius, bulk_density,
        &radius_array_view[0], &density_array_view[0], &gravity_array_view[0],
        &bulk_array_view[0], &shear_array_view[0],
        &frequency_to_use, &G_to_use
        )
    
    return frequency_to_use, G_to_use

def redimensionalize_radial_functions(
    double complex[:, ::1] radial_function_view,
    double mean_radius,
    double bulk_density):
    """ A function to re-dimensionalize physical parameters that have been previously non-dimensionalized.

    Parameters
    ----------
    radial_function_ptr : complex128*
        Non-dimensionalized radial solutions as a function of radius.
    mean_radius : float64
        Mean radius of the planet, used in scaling [m]
    bulk_density : float64
        Bulk density of the planet, used in scaling [m]
    num_slices : uint32
        Number of radial slices, used for looping.
    num_solutions : uint32, default=1
        Number of solutions to loop through (size of radial_function_ptr is 6 * num_solutions * num_slices)

    """
    # Size of arrays
    cdef size_t num_solutions, num_slices
    num_solutions = <size_t> (radial_function_view.shape[0] / 6)
    num_slices = radial_function_view.shape[1]

    # Call cython function
    cf_redimensionalize_radial_functions(
        &radial_function_view.T[0, 0],
        mean_radius,
        bulk_density,
        num_slices,
        num_solutions
    )
