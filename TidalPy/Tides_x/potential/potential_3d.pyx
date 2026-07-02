# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
potential_3d.pyx
Python/Cython wrapper for the dynamic Kaula 3D tidal-potential engine (potential_3d_.hpp).

``tidal_potential_3d_modes(...)`` returns every active tidal mode's degree ``l``, signed forcing
frequency ``omega_lmpq`` [rad s-1], and potential angular factor row
``(U, dU/dtheta, dU/dphi, d2U/dtheta2, d2U/dphi2, d2U/dtheta_dphi)`` evaluated at a point. The modes and
coefficients are built from the eccentricity/obliquity functions (the same used by the global 1D path)
plus the associated Legendre functions, following Kaula / Efroimsky & Williams (2009) Eq. 18. This is
the class-free replacement for the old per-scenario tidal-potential models.

All quantities MKS; frequencies rad s-1; angles radians. The potential's ``r^2`` coefficient uses the
supplied ``planet_radius`` (pass the surface radius for the 3D kernel).
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
        double orbital_frequency,
        double spin_frequency,
        double eccentricity,
        double obliquity,
        double host_mass,
        double semi_major_axis,
        double planet_radius,
        double colatitude,
        double longitude,
        double time,
        double G_to_use,
        int max_degree_l=2,
        int eccentricity_truncation=6,
        int obliquity_truncation=0,
        int min_degree_l=2):
    """Active tidal modes and their potential angular factors at one point.

    Returns
    -------
    degrees : numpy.ndarray of int32, shape (num_modes,)
        Spherical-harmonic degree l of each active mode.
    modes : numpy.ndarray of float64, shape (num_modes,)
        Signed forcing frequency omega_lmpq [rad s-1] of each active mode.
    potentials : numpy.ndarray of float64, shape (num_modes, 6)
        Row i: ``(U, dU/dtheta, dU/dphi, d2U/dtheta2, d2U/dphi2, d2U/dtheta_dphi)`` for mode i.
    """
    cdef int error_code = 0
    cdef vector[c_TidalPotential3DMode] modes = c_tidal_potential_3d_modes(
        planet_radius, semi_major_axis, orbital_frequency, spin_frequency,
        obliquity, eccentricity, host_mass, G_to_use,
        min_degree_l, max_degree_l, obliquity_truncation, eccentricity_truncation,
        colatitude, longitude, time, &error_code)
    if error_code != 0:
        raise ValueError(
            f"TidalPy: tidal potential engine failed (error {error_code}); check degree/truncation levels")

    cdef Py_ssize_t num = <Py_ssize_t>modes.size()
    cdef cnp.ndarray[cnp.int32_t, ndim=1] degrees = np.empty(num, dtype=np.int32)
    cdef cnp.ndarray[cnp.float64_t, ndim=1] freqs = np.empty(num, dtype=np.float64)
    cdef cnp.ndarray[cnp.float64_t, ndim=2] pots = np.empty((num, 6), dtype=np.float64)
    cdef Py_ssize_t i
    for i in range(num):
        degrees[i] = modes[i].degree_l
        freqs[i]   = modes[i].mode_frequency
        pots[i, 0] = modes[i].potential.U
        pots[i, 1] = modes[i].potential.dU_dtheta
        pots[i, 2] = modes[i].potential.dU_dphi
        pots[i, 3] = modes[i].potential.d2U_dtheta2
        pots[i, 4] = modes[i].potential.d2U_dphi2
        pots[i, 5] = modes[i].potential.d2U_dtheta_dphi
    return degrees, freqs, pots
