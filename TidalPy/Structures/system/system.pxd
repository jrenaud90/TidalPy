# distutils: language = c++
"""Cython declarations for TidalPy's system class: c_System and the Python wrapper System."""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr, shared_ptr
from libcpp.vector cimport vector

# Cimported so the CyRK solver implementation (cysolve.cpp / cysolution.cpp, force-included by CyRK's
# own pxd) compiles into this extension: the world tidal solve the system drives runs the CyRK-backed
# radial solver, whose symbols must be linked here (the same cimport the world extensions use).
from CyRK cimport ODEMethod

from TidalPy.Utilities.classes.classes cimport TidalPyBaseClass, c_TidalPyBaseClass
from TidalPy.Structures.worlds.base cimport BaseWorld, c_BaseWorld


# World binary-dispatch helpers (worlds/factory_.hpp): reconstruct + type-discriminate a loaded world.
cdef extern from "factory_.hpp" namespace "tidalpy" nogil:
    int c_world_kind(const c_BaseWorld* world)


cdef extern from "system_.hpp" namespace "tidalpy" nogil:
    # std::invalid_argument for a semi-major axis that is not positive or an eccentricity outside [0, 1); NaN skips.
    void c_check_orbit(double semi_major_axis, double eccentricity, const string& world_name) except +

    cdef cppclass c_TidalDissipation:
        size_t   world_index
        size_t   companion_index
        cpp_bool solved
        cpp_bool has_tide_model
        double   orbital_frequency
        double   semi_major_axis
        double   eccentricity
        double   spin_frequency
        double   obliquity
        double   companion_mass
        double   target_mass
        double   tidal_heating
        double   dU_dM
        double   dU_dw
        double   dU_dO
        double   dU_dM_minus_dw
        double   moment_of_inertia

    cdef cppclass c_WorldEvolution:
        size_t   world_index
        cpp_bool evolved
        cpp_bool has_tide_model
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
        double dU_dM_minus_dw
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
        size_t           first_index
        size_t           second_index
        cpp_bool         evolved
        cpp_bool         has_tide_model
        double           orbital_frequency
        double           semi_major_axis
        double           eccentricity
        double           da_dt
        double           de_dt
        double           dn_dt
        c_WorldEvolution first
        c_WorldEvolution second
        double           tidal_heating_total
        double           dE_orbit_dt
        double           dE_spin_dt_total
        double           energy_residual

    cdef cppclass c_System(c_TidalPyBaseClass):
        c_System() except +
        c_System(const string& name) except +
        const string& get_name() const
        void   set_name(const string& name)
        const shared_ptr[c_BaseWorld]& get_world(size_t index) except +
        size_t add_world(
            shared_ptr[c_BaseWorld] world,
            cpp_bool is_star,
            double semi_major_axis,
            double eccentricity) except +
        size_t   get_num_worlds() const
        int      find_world_index(const string& name)
        cpp_bool has_tidal_host(size_t index) except +
        int      get_tidal_host_index(size_t index) except +
        void     set_tidal_host(size_t index, size_t host_index) except +
        void     clear_tidal_host(size_t index) except +
        double   get_tidal_host_mass(size_t index) except +
        cpp_bool is_mutual_pair(size_t index) except +
        cpp_bool is_hosted_by_star(size_t index) except +
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
        c_TidalDissipation       calc_dissipation(size_t index) except +
        c_WorldEvolution         calc_world_evolution(size_t index) except +
        vector[c_WorldEvolution] calc_system_evolution() except +
        c_PairEvolution          calc_pair_evolution(size_t index) except +
        c_PairEvolution          calc_pair_evolution(size_t first_index, size_t second_index) except +
        double   calc_orbital_energy_derivative(const c_WorldEvolution& evolution) except +
        double   calc_spin_energy_derivative(const c_WorldEvolution& evolution)
        double   calc_energy_residual(const c_WorldEvolution& evolution) except +


# The pair evolution driver (evolution_.hpp): System.evolve and the spin-root tools its tests drive with synthetic
# balances.
cdef extern from "evolution_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_PairEvolveSettings:
        cpp_bool evolve_thermal
        double   orbit_rtol
        double   thermal_rtol
        double   spin_rtol
        double   atol
        double   capture_band
        double   capture_margin
        double   root_tolerance
        double   resolution
        double   max_wall_time

    cdef cppclass c_BodySegment:
        int    mode
        double spin_ratio
        double root
        double window_lower
        double window_upper
        double target
        string reason

    cdef cppclass c_PairSegment:
        double                time
        vector[c_BodySegment] bodies
        string                ended

    cdef cppclass c_BodyEvolutionRecord:
        size_t                world_index
        cpp_bool              rigid
        cpp_bool              thermal
        size_t                num_layers
        vector[double]        spin_ratio
        vector[double]        spin_frequency
        vector[double]        tidal_heating
        vector[double]        da_dt
        vector[double]        de_dt
        vector[double]        dspin_dt
        vector[unsigned char] tracked
        vector[double]        temperature
        size_t                num_tide_solves
        size_t                num_eos_solves

    cdef cppclass c_PairEvolutionRecord:
        vector[double]                time
        vector[double]                semi_major_axis
        vector[double]                eccentricity
        vector[double]                da_dt
        vector[double]                de_dt
        vector[double]                dn_dt
        vector[c_BodyEvolutionRecord] bodies
        vector[c_PairSegment]         segments
        cpp_bool                      success
        string                        message
        double                        elapsed

    shared_ptr[c_PairEvolutionRecord] c_evolve_pair(
        c_System* system_ptr,
        size_t world_index,
        double t_start,
        double t_end,
        const c_PairEvolveSettings& settings) except +
    shared_ptr[c_PairEvolutionRecord] c_new_pair_evolution_record() except +

    cdef cppclass c_SpinRootResult:
        cpp_bool found
        double   lower
        double   upper
        double   spin_ratio
        double   balance
        double   noise
        size_t   num_probes

    ctypedef double (*c_SpinBalanceFunction)(void*, double) noexcept

    vector[double] c_spin_search_points(double start, double center, cpp_bool moving_up, double finest) except +
    c_SpinRootResult c_locate_callback_root(
        c_SpinBalanceFunction function,
        void* context,
        double spin_ratio,
        double tolerance,
        double finest) except +
    c_SpinRootResult c_refine_callback_root(
        c_SpinBalanceFunction function,
        void* context,
        double lower,
        double upper,
        double tolerance) except +
    c_SpinRootResult c_build_callback_window(
        c_SpinBalanceFunction function,
        void* context,
        double root_ratio,
        double resolution) except +

cdef class System(TidalPyBaseClass):
    cdef unique_ptr[c_System] _system
    cdef list _world_wrappers   # Python list of the added BaseWorld wrappers (co-own the C++ worlds)
    cdef public dict source_config   # system config the system was built from (or None if built directly)
    cdef public object source_dir    # folder of the system file its relative paths are relative to, or None
    cdef Py_ssize_t _resolve_index(self, object world) except *
    cdef void _rebuild_world_wrappers(self)
    cpdef dict get_config_dict(self)
