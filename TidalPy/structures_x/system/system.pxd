# distutils: language = c++
"""
system.pxd
Cython declarations for TidalPy's system class.

Exports c_System and the Python wrapper System so other extensions can cimport and use C-speed access.
"""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr, shared_ptr

from TidalPy.structures_x.worlds.base cimport BaseWorld, c_BaseWorld


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "system_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_System:
        c_System() except +
        c_System(const string& name) except +
        const string& get_name() const
        void   set_name(const string& name)
        size_t add_world(
            shared_ptr[c_BaseWorld] world,
            cpp_bool is_host,
            cpp_bool is_star,
            double semi_major_axis,
            double eccentricity) except +
        size_t   get_num_worlds() const
        int      find_world_index(const string& name)
        cpp_bool has_host() const
        int      get_host_index() const
        void     set_host(size_t index) except +
        double   get_host_mass() except +
        cpp_bool has_star() const
        int      get_star_index() const
        void     set_star(size_t index) except +
        double   get_star_mass() except +
        double   get_star_luminosity() except +
        void     set_semi_major_axis(size_t index, double semi_major_axis) except +
        void     set_eccentricity(size_t index, double eccentricity) except +
        double   get_semi_major_axis(size_t index) except +
        double   get_eccentricity(size_t index) except +
        double   calc_gravitational_parameter(size_t index) except +
        double   calc_orbital_frequency(size_t index) except +
        double   calc_semi_major_axis_from_frequency(size_t index, double orbital_frequency) except +
        void     set_stellar_semi_major_axis(size_t index, double semi_major_axis) except +
        void     set_stellar_eccentricity(size_t index, double eccentricity) except +
        double   get_stellar_semi_major_axis(size_t index) except +
        double   get_stellar_eccentricity(size_t index) except +
        double   calc_stellar_gravitational_parameter(size_t index) except +
        double   calc_stellar_orbital_frequency(size_t index) except +
        double   calc_insolation_flux(size_t index) except +
        double   calc_equilibrium_temperature(size_t index) except +


# =====================================================================================================================
# Cython wrapper class declaration
# =====================================================================================================================
cdef class System:
    cdef unique_ptr[c_System] _system
    cdef list _world_wrappers   # Python list of the added BaseWorld wrappers (co-own the C++ worlds)
    cdef Py_ssize_t _resolve_index(self, object world) except *
