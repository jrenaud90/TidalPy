# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Non-dimensionalization scales for multi-layer tidal calculations.

The scheme is based on that proposed by Martens (2016, around page 99): a frequency-independent
time scale from 1/(pi G rho_bulk), a length scale from the mean radius, and a density scale from
the bulk density; the mass and pascal scales follow.

References
----------
Martens16 : H. Martens, PhD Thesis (CalTech), 2016, DOI: 10.7907/Z9N29TX7
"""

from TidalPy.constants cimport get_shared_config_address, set_tidalpy_config_ptr

# Wire this DLL's shared pointer to the process-wide TidalPy config singleton: the scales read G through it.
set_tidalpy_config_ptr(get_shared_config_address())


cdef class NonDimensionalScalesClass:
    """Python wrapper for the ``c_NonDimensionalScales`` conversion-scale struct; every scale is NaN until built."""

    cdef c_NonDimensionalScales nondim_scales

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


def build_nondimensional_scales(
        double mean_radius,
        double bulk_density
        ):
    """Build a populated :class:`NonDimensionalScalesClass` from a planet's scales.

    Parameters
    ----------
    mean_radius : float
        Planet mean radius [m]; sets the length scale.
    bulk_density : float
        Planet bulk density [kg m-3]; sets the density scale and, with G, the time scale.

    Returns
    -------
    NonDimensionalScalesClass
        Conversion factors from non-dimensional solve units back to MKS.

    Notes
    -----
    The time scale ``sqrt(1 / (pi G rho_bulk))`` carries no forcing frequency, so one set of scales serves
    every frequency in a sweep.
    """

    cdef NonDimensionalScalesClass non_dim_scales = NonDimensionalScalesClass()
    non_dim_scales.nondim_scales = c_NonDimensionalScales(mean_radius, bulk_density)
    return non_dim_scales
