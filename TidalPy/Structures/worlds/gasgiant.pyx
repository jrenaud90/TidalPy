# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrapper for TidalPy's gas-giant world class.

GasGiantWorld behaves like a BaseWorld (it owns layers and supports the whole-planet EOS solve) but carries its own
world type and binary class id.
"""

from libcpp.memory cimport make_shared, shared_ptr, static_pointer_cast

from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Structures.worlds.base cimport (
    BaseWorld, c_BaseWorld, c_WorldConfig, cy_fill_world_config, cy_set_default_spin)

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef class GasGiantWorld(BaseWorld):
    """A world representing a gas giant.

    Identical construction and API to :class:`BaseWorld`; the ``world_type`` defaults to ``"gasgiant"`` and the
    binary records use a dedicated class id.
    """

    def __init__(
            self,
            str name,
            double radius,
            double mass,
            str world_type = "gasgiant",
            albedo = None,
            emissivity = None,
            obliquity = None,
            spin_frequency = None):
        cdef c_WorldConfig config
        cdef dict defaults = cy_fill_world_config(
            &config,
            name,
            radius,
            mass,
            world_type,
            albedo,
            emissivity,
            obliquity,
            spin_frequency)
        self._bind(static_pointer_cast[c_BaseWorld, c_GasGiantWorld](make_shared[c_GasGiantWorld](config)))
        cy_set_default_spin(self, defaults)

    @staticmethod
    cdef GasGiantWorld _wrap(shared_ptr[c_BaseWorld] ptr):
        """Wrap an already-constructed C++ gas-giant world (no new C++ object is built)."""
        cdef GasGiantWorld world = GasGiantWorld.__new__(GasGiantWorld)
        world._bind(ptr)
        return world

    def family_world_type(self) -> str:
        """Builder world ``type`` for gas giants."""
        return "gasgiant"
