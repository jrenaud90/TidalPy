# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrapper for TidalPy's terrestrial world class.

TerrestrialWorld behaves like a BaseWorld (it owns layers and supports the whole-planet EOS, Love, and tidal solves)
but carries its own world type and binary class id. Rocky and icy planets and moons are built as terrestrial worlds.
"""

from libcpp.memory cimport make_shared, shared_ptr, static_pointer_cast

from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Structures.worlds.base cimport BaseWorld, c_BaseWorld, c_WorldConfig, cy_fill_world_config

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef class TerrestrialWorld(BaseWorld):
    """A world representing a rocky or icy planet or moon.

    Identical construction and API to :class:`BaseWorld`; the ``world_type`` defaults to ``"terrestrial"`` and the
    binary records use a dedicated class id.
    """

    def __init__(
            self,
            str name,
            double radius,
            double mass,
            str world_type = "terrestrial",
            double albedo = 0.3,
            double emissivity = 1.0,
            double obliquity = 0.0,
            double spin_frequency = 0.0):
        cdef c_WorldConfig config
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
        self._bind(static_pointer_cast[c_BaseWorld, c_TerrestrialWorld](make_shared[c_TerrestrialWorld](config)))

    @staticmethod
    cdef TerrestrialWorld _wrap(shared_ptr[c_BaseWorld] ptr):
        """Wrap an already-constructed C++ terrestrial world (no new C++ object is built)."""
        cdef TerrestrialWorld world = TerrestrialWorld.__new__(TerrestrialWorld)
        world._bind(ptr)
        return world

    def family_world_type(self) -> str:
        """Builder world ``type`` for terrestrial worlds."""
        return "terrestrial"
