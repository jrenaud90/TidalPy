# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
orbit_solver.pyx
Cython/Python wrapper for TidalPy's orbital rate calculator (Dynamics_x).

OrbitSolver turns the tidal-potential derivatives of a dissipating body (``dU_dM``, ``dU_dw`` from
``world.calc_tides``) plus the orbital state into the instantaneous rates of change of the semi-major
axis, eccentricity, and mean motion. It computes rates only; the System class integrates them.
"""

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef c_OrbitState _make_state(
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


# =====================================================================================================================
# OrbitSolver
# =====================================================================================================================
cdef class OrbitSolver:
    """Orbital rate calculator from tidal dissipation (rates only).

    The tidal-potential derivatives ``dU_dM`` (wrt mean anomaly) and ``dU_dw`` (wrt argument of
    pericenter) of the dissipating (target) body come from the global tidal solve
    (``world.calc_tides`` / ``world.get_tidal_potential_derivatives``). All quantities are MKS:
    ``orbital_frequency`` [rad s-1], ``semi_major_axis`` [m], ``eccentricity`` [-], masses [kg],
    ``dU_d*`` [J kg-1 rad-1]; the rates are ``da/dt`` [m s-1], ``de/dt`` [s-1], ``dn/dt`` [rad s-2].

    For a dual-body dissipation system the two bodies' rates are additive in the disturbing-function
    derivatives; sum the per-body results.
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
        cdef c_OrbitState state = _make_state(
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
            double dU_dw) -> float:
        """Eccentricity rate [s-1]: ``de/dt = (sqrt(1-e^2)/(n a^2 e))(sqrt(1-e^2) dR/dM - dR/dw)``.

        Returns ``0.0`` for a circular (or degenerate) orbit, where the ``1/e`` term is indeterminate.
        """
        cdef c_OrbitState state = _make_state(
            orbital_frequency,
            semi_major_axis,
            eccentricity,
            target_mass,
            host_mass)
        return self._solver.calc_de_dt(state, dU_dM, dU_dw)

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
            double dU_dw) -> dict:
        """All three rates as a dict with keys ``da_dt``, ``de_dt``, ``dn_dt``."""
        cdef c_OrbitState state = _make_state(
            orbital_frequency,
            semi_major_axis,
            eccentricity,
            target_mass,
            host_mass)
        cdef c_OrbitDerivatives out = self._solver.calc_derivatives(state, dU_dM, dU_dw)
        return {'da_dt': out.da_dt, 'de_dt': out.de_dt, 'dn_dt': out.dn_dt}
