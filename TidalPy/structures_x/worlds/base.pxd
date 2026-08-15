# distutils: language = c++
"""
base.pxd
Cython declarations for TidalPy's base world class.

Exports c_WorldConfig, c_BaseWorld, and the Python wrapper BaseWorld so other
extensions can cimport and use C-speed access.
"""

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
        cpp_bool latitude_analytic

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
        double   radius_m
        double   mass_kg
        double   albedo
        double   emissivity
        double   obliquity_rad
        double   spin_frequency_rad_s

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
        double   calc_equilibrium_temperature(double insolation_flux_w_m2) const
        void     set_name(const string& name) except +
        void     set_spin_frequency(double freq_rad_s)
        void     set_obliquity(double obliq_rad)
        # Global (1D) tidal dissipation (analytic path; calc_tides defined in world_tides_base_.hpp).
        void                 set_tide_model(unique_ptr[c_TideBase] tide)
        cpp_bool             get_tide_model_set() const
        void                 set_tide_config(const c_TideConfig& cfg)
        void                 calc_tides(const c_TideSolveConfig& state) except +
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
    cpdef dict get_config_dict(self)
    # Wrap an already-constructed C++ world (e.g. one loaded by c_System::read_binary) as a Python
    # wrapper without building a new C++ object. Each subclass overrides to return its own type.
    @staticmethod
    cdef BaseWorld _wrap(shared_ptr[c_BaseWorld] ptr)
