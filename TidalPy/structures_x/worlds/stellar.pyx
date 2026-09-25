# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrapper for TidalPy's star world class.

StarWorld has no internal layers and no equation of state. Its effective temperature and luminosity are kept
consistent through the Stefan-Boltzmann law (L = 4·pi·R²·sigma·T⁴), and an optional ``LuminosityBase`` model
derives both from the star's mass.
"""

from libcpp.utility cimport move
from libcpp.memory cimport make_shared, shared_ptr, static_pointer_cast

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport c_TidalPyBaseClass, c_PhysicsBase, cy_physics_model_config
from TidalPy.structures_x.worlds.base cimport BaseWorld, c_BaseWorld, cy_fill_world_config
from TidalPy.stellar_x.luminosity cimport LuminosityBase

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef class StarWorld(BaseWorld):
    """A star: no layers, no EOS; effective temperature and luminosity.

    Parameters
    ----------
    name : str
        Star name.
    radius : float
        Stellar radius [m].
    mass : float
        Stellar mass [kg].
    effective_temperature : float, optional
        Effective temperature [K]. Default ``5772.0`` (solar).
    luminosity : float, optional
        Luminosity [W]. Default ``0.0`` => derived from the effective
        temperature via the Stefan-Boltzmann law.
    world_type : str, optional
        Type label. Default ``"star"``.
    albedo, emissivity, obliquity, spin_frequency : float, optional
        See :class:`BaseWorld`.
    """

    def __cinit__(self, *args, **kwargs):
        self._star_ptr = NULL

    def __init__(
            self,
            str    name,
            double radius,
            double mass,
            double effective_temperature = 5772.0,
            double luminosity = 0.0,
            str    world_type = "star",
            double albedo     = 0.0,
            double emissivity = 1.0,
            double obliquity  = 0.0,
            double spin_frequency = 0.0):
        cdef c_StarConfig config
        cy_fill_world_config(
            &config,
            name,
            radius,
            mass,
            world_type,
            albedo,
            emissivity,
            obliquity,
            spin_frequency)
        config.effective_temperature = effective_temperature
        config.luminosity            = luminosity
        self._bind(static_pointer_cast[c_BaseWorld, c_StarWorld](make_shared[c_StarWorld](config)))

    def __dealloc__(self):
        self._star_ptr = NULL  # BaseWorld._world_ptr owns the C++ object

    cdef void _bind(self, shared_ptr[c_BaseWorld] ptr):
        BaseWorld._bind(self, ptr)
        self._star_ptr = <c_StarWorld*>ptr.get()

    @staticmethod
    cdef StarWorld _wrap(shared_ptr[c_BaseWorld] ptr):
        """Wrap an already-constructed C++ star world (no new C++ object is built)."""
        cdef StarWorld world = StarWorld.__new__(StarWorld)
        world._bind(ptr)
        return world

    @property
    def effective_temperature(self) -> float:
        """Effective temperature [K]."""
        return self._star_ptr.get_effective_temperature()

    @property
    def luminosity(self) -> float:
        """Luminosity [W]."""
        return self._star_ptr.get_luminosity()

    def calc_luminosity_from_temperature(self, double temperature) -> float:
        """Stefan-Boltzmann luminosity [W] = 4·pi·R²·sigma·T⁴."""
        return self._star_ptr.calc_luminosity_from_temperature(temperature)

    def calc_temperature_from_luminosity(self, double luminosity) -> float:
        """Effective temperature [K] from luminosity via Stefan-Boltzmann."""
        return self._star_ptr.calc_temperature_from_luminosity(luminosity)

    def set_effective_temperature(self, double temperature):
        """Set effective temperature [K]; recomputes luminosity."""
        self._star_ptr.set_effective_temperature(temperature)

    def set_luminosity(self, double luminosity):
        """Set luminosity [W]; recomputes effective temperature."""
        self._star_ptr.set_luminosity(luminosity)

    # Luminosity model (mass -> luminosity, using the star's own mass and radius)
    def set_luminosity_model(self, LuminosityBase model not None):
        """Attach a :class:`~TidalPy.stellar_x.LuminosityBase` model (transfers ownership).

        Ownership of the C++ model moves out of ``model``, which is left an empty shell and must not be reused.
        Once attached, the star can derive its luminosity and effective temperature from its own mass.
        """
        if model._luminosity_ptr.get() == NULL:
            raise ValueError("This luminosity model holds no C++ object (already attached or moved).")
        self._star_ptr.set_luminosity_model(move(model._luminosity_ptr))
        model._ptr = NULL

    @property
    def luminosity_model_set(self) -> bool:
        """Whether a luminosity model has been attached."""
        return True if self._star_ptr.has_luminosity_model() else False

    def calc_luminosity_from_mass(self) -> float:
        """Luminosity [W] derived from the star's mass; raises RuntimeError when no model is attached."""
        return self._star_ptr.calc_luminosity_from_mass()

    def calc_effective_temperature_from_mass(self) -> float:
        """Effective temperature [K] from the star's mass (mass -> L -> T); raises RuntimeError without a model."""
        return self._star_ptr.calc_effective_temperature_from_mass()

    def update_luminosity_from_mass(self):
        """Update the stored luminosity and effective temperature from the star's mass. Raises RuntimeError when
        no luminosity model is attached."""
        self._star_ptr.update_luminosity_from_mass()

    def family_world_type(self) -> str:
        """Builder world ``type`` for stars."""
        return "star"

    cpdef dict get_config_dict(self):
        """Return the BaseWorld config dict plus the stellar values and the attached luminosity model.

        Returns
        -------
        dict
            The :class:`BaseWorld` keys, ``effective_temperature_k``, ``luminosity_w``, and a ``luminosity`` table
            (``model`` plus the model's parameters) when a luminosity model is attached.
        """
        cdef dict d = BaseWorld.get_config_dict(self)
        d["effective_temperature_k"] = self._star_ptr.get_effective_temperature()
        d["luminosity_w"]            = self._star_ptr.get_luminosity()
        cdef const c_PhysicsBase* model_ptr = <const c_PhysicsBase*>self._star_ptr.get_luminosity_model()
        if model_ptr != NULL:
            d["luminosity"] = cy_physics_model_config(model_ptr)
        return d
