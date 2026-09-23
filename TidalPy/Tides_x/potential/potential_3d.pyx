# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Python and Cython wrapper for the dynamic Kaula 3D tidal-potential engine.

``tidal_potential_3d_modes`` returns every active tidal mode's degree, signed forcing frequency, and
complex potential angular-factor amplitude row, the mode's time factor pulled out so
``U(t) = Re[U_c e^{i omega t}]``. Modes and coefficients come from the eccentricity and obliquity
functions the global 1D path uses, plus the associated Legendre functions, following Kaula and Efroimsky
and Williams (2009) Eq. 18.

All quantities MKS. The potential's ``r^2`` coefficient uses the supplied ``planet_radius``; pass the
surface radius for the 3D kernel.
"""

import numpy as np
cimport numpy as cnp
cnp.import_array()

from libcpp.vector cimport vector

from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Tides_x.potential.potential_3d cimport (
    c_TidalPotential3DMode,
    c_tidal_potential_3d_modes,
)

# Wire this extension's config singleton pointer to the process-wide TidalPy config.
set_tidalpy_config_ptr(get_shared_config_address())


def tidal_potential_3d_modes(
        double planet_radius,
        double orbital_frequency,
        double spin_frequency,
        double eccentricity,
        double obliquity,
        double semi_major_axis,
        double host_mass,
        double G_to_use,
        double colatitude,
        double longitude,
        int min_degree_l=2,
        int max_degree_l=2,
        int eccentricity_truncation=3,
        int obliquity_truncation=10):
    """Active tidal modes with complex potential angular-factor amplitudes at one point.

    The body radius comes first, then the orbital state in the same order as the world's ``calc_tides``,
    then Newton's constant, the point's colatitude and longitude, and the degree range and truncations. The
    obliquity truncation defaults to the general form (10), as the world's does; 0 ignores the obliquity.

    Returns
    -------
    degrees : numpy.ndarray of int32, shape (num_modes,)
        Spherical-harmonic degree l of each active mode.
    modes : numpy.ndarray of float64, shape (num_modes,)
        Signed forcing frequency omega_lmpq [rad s-1] of each active mode.
    potentials : numpy.ndarray of complex128, shape (num_modes, 6)
        Row i: the complex amplitudes ``(U, dU/dtheta, dU/dphi, d2U/dtheta2, d2U/dphi2,
        d2U/dtheta_dphi)`` for mode i (with the mode's e^{i omega t} pulled out).
    """
    if eccentricity_truncation not in (1, 2, 3, 4, 5, 10, 15, 20):
        raise NotImplementedError(
            f'Eccentricity truncation {eccentricity_truncation} is not tabulated. '
            'Supported levels: 1, 2, 3, 4, 5, 10, 15, 20.')
    if obliquity_truncation not in (0, 1, 2, 10):
        raise NotImplementedError(
            f'Obliquity truncation {obliquity_truncation} is not tabulated. '
            'Supported levels: 0 (off), 1, 2, 10 (fully general).')

    cdef int error_code = 0
    cdef vector[c_TidalPotential3DMode] modes = c_tidal_potential_3d_modes(
        planet_radius,
        semi_major_axis,
        orbital_frequency,
        spin_frequency,
        obliquity,
        eccentricity,
        host_mass,
        G_to_use,
        min_degree_l,
        max_degree_l,
        obliquity_truncation,
        eccentricity_truncation,
        colatitude,
        longitude,
        &error_code)
    if error_code != 0:
        raise ValueError(
            f"TidalPy: tidal potential engine failed (error {error_code}); check degree/truncation levels")

    cdef Py_ssize_t num = <Py_ssize_t>modes.size()
    cdef cnp.ndarray degrees = np.empty(num, dtype=np.int32)
    cdef cnp.ndarray freqs = np.empty(num, dtype=np.float64)
    cdef cnp.ndarray pots = np.empty((num, 6), dtype=np.complex128)
    cdef cnp.int32_t[::1] degrees_mv = degrees
    cdef double[::1] freqs_mv = freqs
    cdef double complex[:, ::1] pots_mv = pots
    cdef Py_ssize_t i
    cdef c_TidalPotential3DMode* mode_ptr = NULL

    # One element pointer and C stores, rather than eight indexings and six boxed complex objects per mode.
    with nogil:
        for i in range(num):
            mode_ptr = &modes[i]
            degrees_mv[i] = mode_ptr.degree_l
            freqs_mv[i]   = mode_ptr.mode_frequency
            pots_mv[i, 0] = mode_ptr.potential.U.real() + 1j * mode_ptr.potential.U.imag()
            pots_mv[i, 1] = mode_ptr.potential.dU_dtheta.real() + 1j * mode_ptr.potential.dU_dtheta.imag()
            pots_mv[i, 2] = mode_ptr.potential.dU_dphi.real() + 1j * mode_ptr.potential.dU_dphi.imag()
            pots_mv[i, 3] = mode_ptr.potential.d2U_dtheta2.real() + 1j * mode_ptr.potential.d2U_dtheta2.imag()
            pots_mv[i, 4] = mode_ptr.potential.d2U_dphi2.real() + 1j * mode_ptr.potential.d2U_dphi2.imag()
            pots_mv[i, 5] = (mode_ptr.potential.d2U_dtheta_dphi.real() +
                             1j * mode_ptr.potential.d2U_dtheta_dphi.imag())
    return degrees, freqs, pots
