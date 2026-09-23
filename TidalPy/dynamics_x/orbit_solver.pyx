# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrappers for TidalPy's orbital rate calculator. Rates only; the System class integrates them."""

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from libc.math cimport NAN

from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef c_OrbitState cy_make_state(
        double orbital_frequency,
        double semi_major_axis,
        double eccentricity,
        double target_mass,
        double host_mass):
    cdef c_OrbitState state
    state.orbital_frequency = orbital_frequency
    state.semi_major_axis   = semi_major_axis
    state.eccentricity      = eccentricity
    state.target_mass       = target_mass
    state.host_mass         = host_mass
    return state


cdef class OrbitSolver:
    """Orbital rate calculator from tidal dissipation.

    ``dU_dM`` (wrt mean anomaly) and ``dU_dw`` (wrt argument of pericenter) [J kg-1 rad-1] come from the
    global tidal solve. For a dual-body system the two bodies' rates are additive.
    """

    def calc_da_dt(
            self,
            double orbital_frequency,
            double semi_major_axis,
            double eccentricity,
            double target_mass,
            double host_mass,
            double dU_dM) -> float:
        """Semi-major-axis rate [m s-1]: ``da/dt = (2 / (n a)) dR/dM``."""
        cdef c_OrbitState state = cy_make_state(
            orbital_frequency,
            semi_major_axis,
            eccentricity,
            target_mass,
            host_mass)
        return self._solver.calc_da_dt(state, dU_dM)

    def calc_de_dt(
            self,
            double orbital_frequency,
            double semi_major_axis,
            double eccentricity,
            double target_mass,
            double host_mass,
            double dU_dM,
            double dU_dw,
            dU_dM_minus_dw=None) -> float:
        """Eccentricity rate [s-1]: ``de/dt = (sqrt(1-e^2)/(n a^2 e))(sqrt(1-e^2) dR/dM - dR/dw)``.

        Zero for a circular orbit, where the ``1/e`` term is indeterminate. ``dU_dM_minus_dw``, the collapse's
        per-mode sum of ``dU_dM - dU_dw`` (``tide_result["dUdM_minus_dw"]``), keeps the rate exact at small
        eccentricity, where the two separate sums nearly cancel; without it the difference is formed from them.
        """
        cdef c_OrbitState state = cy_make_state(
            orbital_frequency,
            semi_major_axis,
            eccentricity,
            target_mass,
            host_mass)
        return self._solver.calc_de_dt(state, dU_dM, dU_dw, NAN if dU_dM_minus_dw is None else dU_dM_minus_dw)

    def calc_dn_dt(
            self,
            double orbital_frequency,
            double semi_major_axis,
            double da_dt) -> float:
        """Mean-motion rate [rad s-2]: ``dn/dt = -(3/2)(n / a) da/dt`` (Kepler's third law)."""
        return self._solver.calc_dn_dt(orbital_frequency, semi_major_axis, da_dt)

    def calc_derivatives(
            self,
            double orbital_frequency,
            double semi_major_axis,
            double eccentricity,
            double target_mass,
            double host_mass,
            double dU_dM,
            double dU_dw,
            dU_dM_minus_dw=None) -> dict:
        """All three rates as a dict with keys ``da_dt``, ``de_dt``, ``dn_dt``; see ``calc_de_dt`` for
        ``dU_dM_minus_dw``."""
        cdef c_OrbitState state = cy_make_state(
            orbital_frequency,
            semi_major_axis,
            eccentricity,
            target_mass,
            host_mass)
        cdef c_OrbitDerivatives out = self._solver.calc_derivatives(
            state, dU_dM, dU_dw, NAN if dU_dM_minus_dw is None else dU_dM_minus_dw)
        return {'da_dt': out.da_dt, 'de_dt': out.de_dt, 'dn_dt': out.dn_dt}
