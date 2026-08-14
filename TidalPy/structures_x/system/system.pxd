# distutils: language = c++
"""
system.pxd
Cython declarations for TidalPy's system class.

Exports c_System and the Python wrapper System so other extensions can cimport and use C-speed access.
"""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr, shared_ptr
from libcpp.vector cimport vector

# Cimported so the CyRK solver implementation (cysolve.cpp / cysolution.cpp, force-included by CyRK's
# own pxd) compiles into this extension: the world tidal solve the system drives runs the CyRK-backed
# radial solver, whose symbols must be linked here (the same cimport the world extensions use).
from CyRK cimport ODEMethod

from TidalPy.structures_x.worlds.base cimport BaseWorld, c_BaseWorld


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "system_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_WorldEvolution:
        size_t   world_index
        cpp_bool evolved
        double   orbital_frequency
        double   semi_major_axis
        double   eccentricity
        double   spin_frequency
        double   host_mass
        double   target_mass
        double   tidal_heating
        double   dU_dM
        double   dU_dw
        double   dU_dO
        double   da_dt
        double   de_dt
        double   dn_dt
        double   dspin_dt
        double   moment_of_inertia
        cpp_bool has_spin
        double   dE_orbit_dt
        double   dE_spin_dt
        double   energy_residual

    cdef cppclass c_PairEvolution:
        size_t           world_index
        size_t           host_index
        cpp_bool         evolved
        double           orbital_frequency
        double           semi_major_axis
        double           eccentricity
        double           da_dt
        double           de_dt
        double           dn_dt
        c_WorldEvolution world
        c_WorldEvolution host
        double           tidal_heating_total
        double           dE_orbit_dt
        double           dE_spin_dt_total
        double           energy_residual

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
        c_WorldEvolution         calc_world_evolution(size_t index) except +
        vector[c_WorldEvolution] calc_system_evolution() except +
        c_PairEvolution          calc_pair_evolution(size_t index) except +
        double   calc_orbital_energy_derivative(const c_WorldEvolution& evolution) except +
        double   calc_spin_energy_derivative(const c_WorldEvolution& evolution)
        double   calc_energy_residual(const c_WorldEvolution& evolution) except +


# =====================================================================================================================
# Cython wrapper class declaration
# =====================================================================================================================
cdef class System:
    cdef unique_ptr[c_System] _system
    cdef list _world_wrappers   # Python list of the added BaseWorld wrappers (co-own the C++ worlds)
    cdef public dict source_config   # system config the system was built from (or None if built directly)
    cdef Py_ssize_t _resolve_index(self, object world) except *
    cpdef dict get_config_dict(self)
