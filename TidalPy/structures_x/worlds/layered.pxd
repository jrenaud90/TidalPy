# distutils: language = c++
"""
layered.pxd
Cython declarations for TidalPy's layered world class (Phase 8b).
"""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.complex cimport complex as cpp_complex

from CyRK cimport ODEMethod

from TidalPy.structures_x.worlds.base cimport BaseWorld, c_BaseWorld, c_WorldConfig
from TidalPy.structures_x.layers.base cimport BaseLayer, c_BaseLayer
from TidalPy.Material_x.eos.eos_solution cimport c_EOSSolution


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "layered_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_WorldEOSSolveConfig:
        double    surface_pressure
        size_t    slices_per_layer
        double    G_to_use
        ODEMethod integration_method
        double    rtol
        double    atol
        double    pressure_tol
        size_t    max_iters
        double    temperature
        cpp_bool  verbose

    cdef cppclass c_LayeredWorld(c_BaseWorld):
        c_LayeredWorld()
        c_LayeredWorld(const c_WorldConfig& cfg) except +
        void         add_layer(unique_ptr[c_BaseLayer] layer) except +
        cpp_bool     accepts_layer(const c_BaseLayer& layer) const
        c_BaseLayer* get_layer(size_t index) except +
        size_t       get_num_layers() const
        double       calc_total_mass() const
        double       calc_internal_heating(double time_s) const
        cpp_bool     validate_layers() const
        void         solve_eos(const c_WorldEOSSolveConfig& cfg) except +
        double       get_density(double radius_m) const
        double       get_gravity(double radius_m) const
        double       get_pressure(double radius_m) const
        double       get_shear_modulus(double radius_m) const
        double       get_bulk_modulus(double radius_m) const
        double       get_shear_viscosity(double radius_m) const
        double       get_bulk_viscosity(double radius_m) const
        double       get_premelt_shear_modulus(double radius_m) const
        double       get_premelt_bulk_modulus(double radius_m) const
        double       get_premelt_shear_viscosity(double radius_m) const
        double       get_premelt_bulk_viscosity(double radius_m) const
        cpp_complex[double] calc_complex_shear_modulus(double radius_m, double frequency_rad_s) const
        cpp_complex[double] calc_complex_bulk_modulus(double radius_m, double frequency_rad_s) const
        cpp_bool     get_eos_solved() const
        cpp_bool     get_all_eos_set() const
        cpp_bool     get_eos_success() const
        const string& get_eos_message() const
        int          get_eos_iterations() const
        double       get_eos_pressure_error() const
        double       get_surface_gravity_eos() const
        double       get_surface_pressure_eos() const
        double       get_central_pressure() const
        double       get_planet_mass_eos() const
        double       get_planet_moi_eos() const
        const c_EOSSolution* get_eos_solution() const


# =====================================================================================================================
# Cython wrapper class declaration
# =====================================================================================================================
cdef class LayeredWorld(BaseWorld):
    cdef c_LayeredWorld* _layered_ptr   # non-owning; ownership via BaseWorld._world_ptr
    cpdef dict get_config_dict(self)
