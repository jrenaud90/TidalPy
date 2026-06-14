# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
stellar.pyx
Cython/Python wrapper for TidalPy's star world class.

StarWorld: a star with no internal layers and no equation of state. Carries an
effective temperature and luminosity kept consistent via the Stefan-Boltzmann
law (L = 4·pi·R²·sigma·T⁴). A richer luminosity-model hierarchy arrives in a
later phase.
"""

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport c_TidalPyBaseClass
from TidalPy.structures_x.worlds.base cimport BaseWorld, c_BaseWorld

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# StarWorld
# =====================================================================================================================

cdef class StarWorld(BaseWorld):
    """A star: no layers, no EOS; effective temperature and luminosity.

    Parameters
    ----------
    name : str
        Star name.
    radius_m : float
        Stellar radius [m].
    mass_kg : float
        Stellar mass [kg].
    effective_temperature_k : float, optional
        Effective temperature [K]. Default ``5772.0`` (solar).
    luminosity_w : float, optional
        Luminosity [W]. Default ``0.0`` => derived from the effective
        temperature via the Stefan-Boltzmann law.
    world_type : str, optional
        Type label. Default ``"star"``.
    albedo, emissivity, obliquity_rad, spin_frequency_rad_s : float, optional
        See :class:`BaseWorld`.
    """

    def __cinit__(self, *args, **kwargs):
        self._star_ptr = NULL

    def __init__(
            self,
            str    name,
            double radius_m,
            double mass_kg,
            double effective_temperature_k = 5772.0,
            double luminosity_w            = 0.0,
            str    world_type              = "star",
            double albedo                  = 0.0,
            double emissivity              = 1.0,
            double obliquity_rad           = 0.0,
            double spin_frequency_rad_s    = 0.0):
        cdef c_StarConfig config
        config.name                    = name.encode("utf-8")
        config.world_type_str          = world_type.encode("utf-8")
        config.radius_m                = radius_m
        config.mass_kg                 = mass_kg
        config.albedo                  = albedo
        config.emissivity              = emissivity
        config.obliquity_rad           = obliquity_rad
        config.spin_frequency_rad_s    = spin_frequency_rad_s
        config.effective_temperature_k = effective_temperature_k
        config.luminosity_w            = luminosity_w
        cdef c_StarWorld* raw = new c_StarWorld(config)
        self._world_ptr.reset(<c_BaseWorld*>raw)
        self._star_ptr = raw
        self._ptr      = <c_TidalPyBaseClass*>raw

    def __dealloc__(self):
        self._star_ptr = NULL  # base's unique_ptr owns the C++ object

    # ------------------------------------------------------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def effective_temperature(self) -> float:
        """Effective temperature [K]."""
        return self._star_ptr.get_effective_temperature()

    @property
    def luminosity(self) -> float:
        """Luminosity [W]."""
        return self._star_ptr.get_luminosity()

    # ------------------------------------------------------------------------------------------------------------------
    # Calculations
    # ------------------------------------------------------------------------------------------------------------------
    def calc_luminosity_from_temperature(self, double temperature_k) -> float:
        """Stefan-Boltzmann luminosity [W] = 4·pi·R²·sigma·T⁴."""
        return self._star_ptr.calc_luminosity_from_temperature(temperature_k)

    def calc_temperature_from_luminosity(self, double luminosity_w) -> float:
        """Effective temperature [K] from luminosity via Stefan-Boltzmann."""
        return self._star_ptr.calc_temperature_from_luminosity(luminosity_w)

    # ------------------------------------------------------------------------------------------------------------------
    # Mutators (keep T and L consistent)
    # ------------------------------------------------------------------------------------------------------------------
    def set_effective_temperature(self, double temperature_k):
        """Set effective temperature [K]; recomputes luminosity."""
        self._star_ptr.set_effective_temperature(temperature_k)

    def set_luminosity(self, double luminosity_w):
        """Set luminosity [W]; recomputes effective temperature."""
        self._star_ptr.set_luminosity(luminosity_w)

    # ------------------------------------------------------------------------------------------------------------------
    # Config
    # ------------------------------------------------------------------------------------------------------------------
    cpdef dict get_config_dict(self):
        """Return the world config plus stellar fields.

        Returns
        -------
        dict
            All :class:`BaseWorld` keys plus ``effective_temperature_k`` and
            ``luminosity_w``.
        """
        d = BaseWorld.get_config_dict(self)
        d["effective_temperature_k"] = self._star_ptr.get_effective_temperature()
        d["luminosity_w"]            = self._star_ptr.get_luminosity()
        return d
