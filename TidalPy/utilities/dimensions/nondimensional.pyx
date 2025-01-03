# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
""" Functionality to non-dimensionalize common variables used for multi-layer tidal calculations.

The scheme is based on that proposed by Martens16 (around page 99)

References
----------
Martens16 : H. Martens, PhD Thesis (CalTech), 2016, DOI: 10.7907/Z9N29TX7

"""

from libc.math cimport sqrt

from TidalPy.constants cimport d_G, d_PI_DBL, d_NAN_DBL


cdef class NonDimensionalScalesClass:
    """ Python wrapper for the `NonDimensionalScales` struct. """

    cdef NonDimensionalScalesCC nondim_scales

    def __init__(self):
        # Initialize everything to nan
        self.nondim_scales.second2_conversion = d_NAN_DBL
        self.nondim_scales.second_conversion  = d_NAN_DBL
        self.nondim_scales.length_conversion  = d_NAN_DBL
        self.nondim_scales.length3_conversion = d_NAN_DBL
        self.nondim_scales.density_conversion = d_NAN_DBL
        self.nondim_scales.mass_conversion    = d_NAN_DBL
        self.nondim_scales.pascal_conversion  = d_NAN_DBL
    
    @property
    def second2_conversion(self):
        return self.nondim_scales.second2_conversion

    @property
    def second_conversion(self):
        return self.nondim_scales.second_conversion
    
    @property
    def length_conversion(self):
        return self.nondim_scales.length_conversion
    
    @property
    def length3_conversion(self):
        return self.nondim_scales.length3_conversion
    
    @property
    def density_conversion(self):
        return self.nondim_scales.density_conversion

    @property
    def mass_conversion(self):
        return self.nondim_scales.mass_conversion
    
    @property
    def pascal_conversion(self):
        return self.nondim_scales.pascal_conversion


cdef void cf_build_nondimensional_scales(
        NonDimensionalScalesCC* non_dim_scales_ptr,
        double frequency,
        double mean_radius,
        double bulk_density
        ) noexcept nogil:

    non_dim_scales_ptr.second2_conversion = 1. / (d_PI_DBL * d_G * bulk_density)
    non_dim_scales_ptr.second_conversion  = sqrt(non_dim_scales_ptr.second2_conversion)
    non_dim_scales_ptr.length_conversion  = mean_radius
    non_dim_scales_ptr.length3_conversion = mean_radius * mean_radius * mean_radius
    non_dim_scales_ptr.density_conversion = bulk_density
    non_dim_scales_ptr.mass_conversion    = bulk_density * non_dim_scales_ptr.length3_conversion
    non_dim_scales_ptr.pascal_conversion  = \
        non_dim_scales_ptr.mass_conversion / (non_dim_scales_ptr.length_conversion * non_dim_scales_ptr.second2_conversion)


def build_nondimensional_scales(
        double frequency,
        double mean_radius,
        double bulk_density
        ):

    cdef NonDimensionalScalesClass non_dim_scales = NonDimensionalScalesClass()

    cf_build_nondimensional_scales(
        &non_dim_scales.nondim_scales,
        frequency,
        mean_radius,
        bulk_density
        )

    return non_dim_scales
