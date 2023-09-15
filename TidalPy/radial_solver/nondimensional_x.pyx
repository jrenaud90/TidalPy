""" Functionality to non-dimensionalize common variables used for multi-layer tidal calculations.

The scheme is based on that proposed by Martens16 (around page 99)

References
----------
Martens16 : H. Martens, PhD Thesis (CalTech), 2016, DOI: 10.7907/Z9N29TX7

"""

from TidalPy.constants import G as G_
from libc.math cimport pi, sqrt

cdef double G = G_

cdef void non_dimensionalize_physicals_x(
        Py_ssize_t num_radius,
        double* radius_array_ptr,
        double* density_array_ptr,
        double* gravity_array_ptr,
        double* bulk_array_ptr,
        double* shear_array_ptr,
        double* viscosity_array_ptr,
        double* frequency,
        double* G_to_use,
        double mean_radius,
        double bulk_density
        ):

    # Setup loop variables
    cdef Py_ssize_t i

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
        viscosity_array_ptr[i] /= (pascal_conversion * second_conversion)

    # Convert non-array pointers
    G_to_use[0]   = G / (length_conversion**3 / (mass_conversion * second2_conversion))
    frequency[0] /= (1. / second_conversion)
