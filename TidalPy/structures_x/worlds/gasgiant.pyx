# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrapper for TidalPy's gas-giant world class.

GasGiantWorld behaves like a LayeredWorld (it owns layers and supports the whole-planet EOS solve) but carries
its own world type and binary class id.
"""

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport c_TidalPyBaseClass
from TidalPy.structures_x.worlds.base cimport c_BaseWorld, c_WorldConfig, cy_fill_world_config
from TidalPy.structures_x.worlds.layered cimport LayeredWorld, c_LayeredWorld
from libcpp.memory cimport make_shared, shared_ptr, static_pointer_cast

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


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
            double radius,
            double mass,
            str    world_type = "gasgiant",
            double albedo     = 0.3,
            double emissivity = 1.0,
            double obliquity  = 0.0,
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
        self._bind(static_pointer_cast[c_BaseWorld, c_GasGiantWorld](make_shared[c_GasGiantWorld](config)))

    def family_world_type(self) -> str:
        """Builder world ``type`` for gas giants."""
        return "gasgiant"

    def __dealloc__(self):
        self._gasgiant_ptr = NULL  # BaseWorld._world_ptr owns the C++ object
        self._layered_ptr  = NULL

    cdef void _bind(self, shared_ptr[c_BaseWorld] ptr):
        LayeredWorld._bind(self, ptr)
        self._gasgiant_ptr = <c_GasGiantWorld*>ptr.get()

    @staticmethod
    cdef GasGiantWorld _wrap(shared_ptr[c_BaseWorld] ptr):
        """Wrap an already-constructed C++ gas-giant world (no new C++ object is built)."""
        cdef GasGiantWorld world = GasGiantWorld.__new__(GasGiantWorld)
        world._bind(ptr)
        return world
