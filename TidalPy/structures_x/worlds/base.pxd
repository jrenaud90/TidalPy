# distutils: language = c++
"""Cython declarations for TidalPy's base world class: c_WorldConfig, c_BaseWorld, and the Python wrapper
BaseWorld."""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr, shared_ptr
from libcpp.complex cimport complex as cpp_complex
from libcpp.vector cimport vector

from TidalPy.Utilities_x.classes_x.classes cimport (
    TidalPyBaseClass,
    StructureBase,
    c_StructureBase,
    c_TidalPyBaseClass,
)
from TidalPy.Tides_x.classes.tide cimport c_TideBase


# =====================================================================================================================
# Global (1D) tidal solve config structs (global namespace, tide_result_.hpp). Declared here on
# the base world (the analytic tide path is common to all world types); cimported by subclasses.
# =====================================================================================================================
cdef extern from "tide_result_.hpp" nogil:
    cdef cppclass c_TideConfig:
        int min_degree_l
        int max_degree_l
        int eccentricity_truncation
        int obliquity_truncation
        double tidal_timescale_width_decades
        int love_method
        double love_fixed_q
        double love_fixed_dt

    cdef cppclass c_TideSolveConfig:
        double orbital_frequency
        double spin_frequency
        double eccentricity
        double obliquity
        double semi_major_axis
        double host_mass

    cdef cppclass c_Heating3DCollapseConfig:
        cpp_bool orbit_averaged
        cpp_bool latitude_summed
        cpp_bool longitude_summed
        cpp_bool radial_summed
        int      latitude_nodes
        int      longitude_nodes
        int      radial_slices
        int      num_threads
        cpp_bool latitude_analytic
        double   colatitude_min
        double   colatitude_max

    cdef cppclass c_Heating3DCollapsed:
        vector[double] values
        vector[size_t] shape
        vector[double] radii
        vector[double] colatitudes
        vector[double] longitudes
        vector[double] times
        vector[double] layer_totals
        size_t n_layers
        size_t n_times
        cpp_bool all_spatial_summed


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "base_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_WorldConfig:
        string   name
        string   world_type_str
        double   radius
        double   mass
        double   albedo
        double   emissivity
        double   obliquity
        double   spin_frequency

    cdef cppclass c_BaseWorld(c_StructureBase):
        c_BaseWorld()
        c_BaseWorld(const c_WorldConfig& cfg) except +
        const string& get_name()         const
        const string& get_world_type()   const
        double   get_albedo()            const
        double   get_emissivity()        const
        double   get_obliquity()         const
        double   get_spin_frequency()    const
        double   calc_surface_gravity()  const
        double   calc_escape_velocity()  const
        double   calc_mean_density()     const
        double   calc_equilibrium_temperature(double insolation_flux) const
        void     set_name(const string& name) except +
        void     set_spin_frequency(double freq)
        void     set_obliquity(double obliq)
        # Global (1D) tidal dissipation (analytic path; calc_tides defined in world_tides_base_.hpp).
        void                 set_tide_model(unique_ptr[c_TideBase] tide)
        cpp_bool             get_tide_model_set() const
        const c_TideBase*    get_tide_model() const
        void                 set_tide_config(const c_TideConfig& cfg)
        const c_TideConfig&  get_tide_config() const
        void                 calc_tides(const c_TideSolveConfig& state) except +
        cpp_bool             get_tide_state(c_TideSolveConfig& state_out) except +
        cpp_bool             get_tides_solved() const
        double               get_tidal_heating() const
        double               get_tidal_dU_dM() const
        double               get_tidal_dU_dw() const
        double               get_tidal_dU_dO() const
        int                  get_num_tidal_modes() const
        cpp_complex[double]  get_tidal_love_k(int degree_l, int m, int p, int q) const


# =====================================================================================================================
# Cython wrapper class declaration
# =====================================================================================================================
cdef class BaseWorld(StructureBase):
    cdef shared_ptr[c_BaseWorld] _world_ptr   # owns the most-derived C++ world object (shared so a System can co-own it)
    cdef public dict source_config            # normalized config the world was built from (or None)
    cdef public dict portable_config          # a data-file world's config as given, for save_to_toml (or None)
    cpdef dict get_config_dict(self)
    # Wrap an already-constructed C++ world without building a new one; each subclass returns its own type.
    @staticmethod
    cdef BaseWorld _wrap(shared_ptr[c_BaseWorld] ptr)
