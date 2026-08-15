# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
gasgiant.pyx
Cython/Python wrapper for TidalPy's gas-giant world class.

GasGiantWorld: a layered world representing a gas giant. Behaves like a
LayeredWorld (owns layers, supports the whole-planet EOS solve) but carries a
distinct world type and binary class id.
"""

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport c_TidalPyBaseClass
from TidalPy.structures_x.worlds.base cimport c_BaseWorld, c_WorldConfig
from TidalPy.structures_x.worlds.layered cimport LayeredWorld, c_LayeredWorld

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# GasGiantWorld
# =====================================================================================================================

cdef class GasGiantWorld(LayeredWorld):
    """A layered world representing a gas giant.

    Identical construction and API to :class:`LayeredWorld`; the ``world_type``
    defaults to ``"gasgiant"`` and the binary records use a dedicated class id.
    """

    def __cinit__(self, *args, **kwargs):
        self._gasgiant_ptr = NULL

    def __init__(
            self,
            str    name,
            double radius_m,
            double mass_kg,
            str    world_type           = "gasgiant",
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
        cdef c_GasGiantWorld* raw = new c_GasGiantWorld(config)
        self._world_ptr.reset(<c_BaseWorld*>raw)
        self._layered_ptr  = <c_LayeredWorld*>raw
        self._gasgiant_ptr = raw
        self._ptr          = <c_TidalPyBaseClass*>raw

    def __dealloc__(self):
        self._gasgiant_ptr = NULL  # base's unique_ptr owns the C++ object
        self._layered_ptr  = NULL

    @staticmethod
    cdef GasGiantWorld _wrap(shared_ptr[c_BaseWorld] ptr):
        """Wrap an already-constructed C++ gas-giant world (no new C++ object is built)."""
        cdef GasGiantWorld world = GasGiantWorld.__new__(GasGiantWorld)
        world._world_ptr = ptr
        world._ptr = <c_TidalPyBaseClass*>ptr.get()
        world._layered_ptr = <c_LayeredWorld*>ptr.get()
        world._gasgiant_ptr = <c_GasGiantWorld*>ptr.get()
        world._layer_views = None
        world._layer_view_by_name = None
        return world
