# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
base.pyx
Cython/Python wrapper for TidalPy's base world class.

BaseWorld: world-level identity and orbital/thermal scalars (albedo, emissivity,
obliquity, spin frequency) plus bulk geometry and equilibrium-temperature
calculations. Layered worlds (terrestrial, gas giant) and stars subclass this.
"""

from libcpp cimport bool as cpp_bool

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport StructureBase, c_TidalPyBaseClass

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# BaseWorld
# =====================================================================================================================

cdef class BaseWorld(StructureBase):
    """Base world: identity, orbital/thermal scalars, and bulk geometry.

    Parameters
    ----------
    name : str
        Human-readable world name.
    radius_m : float
        World radius [m].
    mass_kg : float
        World mass [kg].
    world_type : str, optional
        Free-form type label (e.g. ``"terrestrial"``). Default ``"world"``.
    albedo : float, optional
        Bond albedo [dimensionless]. Default ``0.3``.
    emissivity : float, optional
        Surface emissivity [dimensionless]. Default ``1.0``.
    obliquity_rad : float, optional
        Axial obliquity [rad]. Default ``0.0``.
    spin_frequency_rad_s : float, optional
        Rotation rate [rad/s]. Default ``0.0``.

    Assumptions
    -----------
    - Spherically symmetric world.
    - All values are in MKS units.
    """

    def __cinit__(self, *args, **kwargs):
        pass  # unique_ptr<c_BaseWorld> auto-inits to nullptr; _ptr set in __init__

    def __init__(
            self,
            str    name,
            double radius_m,
            double mass_kg,
            str    world_type           = "world",
            double albedo               = 0.3,
            double emissivity           = 1.0,
            double obliquity_rad        = 0.0,
            double spin_frequency_rad_s = 0.0):
        cdef c_WorldConfig config
        config.name                 = name.encode("utf-8")
        config.world_type_str       = world_type.encode("utf-8")
        config.radius_m             = radius_m
        config.mass_kg              = mass_kg
        config.albedo               = albedo
        config.emissivity           = emissivity
        config.obliquity_rad        = obliquity_rad
        config.spin_frequency_rad_s = spin_frequency_rad_s
        self._world_ptr.reset(new c_BaseWorld(config))
        self._ptr = <c_TidalPyBaseClass*>self._world_ptr.get()

    def __dealloc__(self):
        self._world_ptr.reset()
        self._ptr = NULL

    # ------------------------------------------------------------------------------------------------------------------
    # Base class property overrides (read from the most-derived C++ world)
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def radius(self) -> float:
        """World radius [m]."""
        return self._world_ptr.get().get_radius()

    @property
    def mass(self) -> float:
        """World mass [kg]."""
        return self._world_ptr.get().get_mass()

    # ------------------------------------------------------------------------------------------------------------------
    # World properties
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def name(self) -> str:
        """World name."""
        return self._world_ptr.get().get_name().decode("utf-8")

    @property
    def world_type(self) -> str:
        """World type label."""
        return self._world_ptr.get().get_world_type().decode("utf-8")

    @property
    def albedo(self) -> float:
        """Bond albedo [dimensionless]."""
        return self._world_ptr.get().get_albedo()

    @property
    def emissivity(self) -> float:
        """Surface emissivity [dimensionless]."""
        return self._world_ptr.get().get_emissivity()

    @property
    def obliquity(self) -> float:
        """Axial obliquity [rad]."""
        return self._world_ptr.get().get_obliquity()

    @property
    def spin_frequency(self) -> float:
        """Rotation rate [rad/s]."""
        return self._world_ptr.get().get_spin_frequency()

    # ------------------------------------------------------------------------------------------------------------------
    # Calculations
    # ------------------------------------------------------------------------------------------------------------------
    def calc_surface_gravity(self) -> float:
        """Surface gravitational acceleration [m/s^2] = G·M/R²."""
        return self._world_ptr.get().calc_surface_gravity()

    def calc_escape_velocity(self) -> float:
        """Escape velocity [m/s] = sqrt(2·G·M/R)."""
        return self._world_ptr.get().calc_escape_velocity()

    def calc_mean_density(self) -> float:
        """Mean density [kg/m^3] = M / V_sphere(R)."""
        return self._world_ptr.get().calc_mean_density()

    def calc_equilibrium_temperature(self, double insolation_flux_w_m2) -> float:
        """Radiative-equilibrium temperature [K] for a given insolation flux.

        T_eq = [ (1 − A)·F / (4·ε·σ) ]^(1/4)

        Parameters
        ----------
        insolation_flux_w_m2 : float
            Incident stellar flux [W/m^2].

        Returns
        -------
        float
            Equilibrium temperature [K]; 0.0 for non-positive flux.

        Assumptions
        -----------
        - Fast-rotator, uniform-temperature surface.
        """
        return self._world_ptr.get().calc_equilibrium_temperature(insolation_flux_w_m2)

    # ------------------------------------------------------------------------------------------------------------------
    # Mutators
    # ------------------------------------------------------------------------------------------------------------------
    def set_spin_frequency(self, double freq_rad_s):
        """Set the rotation rate [rad/s]."""
        self._world_ptr.get().set_spin_frequency(freq_rad_s)

    def set_obliquity(self, double obliq_rad):
        """Set the axial obliquity [rad]."""
        self._world_ptr.get().set_obliquity(obliq_rad)

    # ------------------------------------------------------------------------------------------------------------------
    # Config
    # ------------------------------------------------------------------------------------------------------------------
    cpdef dict get_config_dict(self):
        """Return all world-level configuration values as a Python dict (MKS).

        Returns
        -------
        dict
            Keys: ``name``, ``world_type``, ``radius_m``, ``mass_kg``,
            ``albedo``, ``emissivity``, ``obliquity_rad``,
            ``spin_frequency_rad_s``.
        """
        cdef c_BaseWorld* p = self._world_ptr.get()
        return {
            "name":                 p.get_name().decode("utf-8"),
            "world_type":           p.get_world_type().decode("utf-8"),
            "radius_m":             p.get_radius(),
            "mass_kg":              p.get_mass(),
            "albedo":               p.get_albedo(),
            "emissivity":           p.get_emissivity(),
            "obliquity_rad":        p.get_obliquity(),
            "spin_frequency_rad_s": p.get_spin_frequency(),
        }
