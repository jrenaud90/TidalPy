# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrappers for TidalPy's spin-dynamics calculator.

``Spin`` provides the moment of inertia, the tidal spin-rate change dspin/dt from the tidal
potential derivative dU/dO, and the synchronous spin rate. It computes rates only; the System class
integrates them.
"""

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# Spin
# =====================================================================================================================
cdef class Spin:
    """Spin-dynamics calculator for a tidally interacting body (rates only).

    Parameters
    ----------
    moment_of_inertia_factor : float, optional
        Conventional dimensionless moment-of-inertia factor ``C / (M R^2)``: ``0.4`` for a uniform sphere
        (default), smaller for a centrally condensed body (``0.3307`` for the Earth). Must lie within
        ``(0, 2/3]``, where ``2/3`` is a thin hollow shell.

    Raises
    ------
    ValueError
        If ``moment_of_inertia_factor`` is not finite or lies outside ``(0, 2/3]``.
    """

    def __init__(self, double moment_of_inertia_factor=0.4):
        cdef c_SpinConfig config
        config.moment_of_inertia_factor = moment_of_inertia_factor
        self._spin = c_Spin(config)

    @property
    def moment_of_inertia_factor(self) -> float:
        """Conventional dimensionless moment-of-inertia factor ``C / (M R^2)`` (0.4 for a uniform sphere)."""
        return self._spin.get_config().moment_of_inertia_factor

    def calc_moment_of_inertia(
            self,
            double mass,
            double radius) -> float:
        """Moment of inertia [kg m2]: ``I = moment_of_inertia_factor * M R^2``.

        Parameters
        ----------
        mass : float
            Body mass [kg].
        radius : float
            Body (outer) radius [m].

        Returns
        -------
        float
            Moment of inertia [kg m2].
        """
        return self._spin.calc_moment_of_inertia(mass, radius)

    def calc_dspin_dt(
            self,
            double host_mass,
            double dU_dO,
            double moment_of_inertia) -> float:
        """Tidal spin-rate change [rad s-2]: ``dspin/dt = M_host * dU_dO / I``.

        ``dU_dO`` is the tidal potential derivative with respect to the longitude of the node
        [J kg-1 rad-1] from the global tidal solve (``world.calc_tides`` /
        ``get_tidal_potential_derivatives``). Returns NaN for a non-positive moment of inertia.
        """
        return self._spin.calc_dspin_dt(host_mass, dU_dO, moment_of_inertia)

    def calc_synchronous_spin(self, double orbital_frequency) -> float:
        """Synchronous spin rate [rad s-1]: equal to the orbital mean motion ``orbital_frequency``."""
        return self._spin.calc_synchronous_spin(orbital_frequency)
