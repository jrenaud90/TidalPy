# distutils: language = c++
"""Cython declarations for the material EOS models, config struct, factory, and wrapper classes."""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_PhysicsBase


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "material_eos_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_MaterialEOSBase(c_PhysicsBase):
        double calc_density(double pressure, double temperature, double radius) const
        double calc_bulk_modulus(double pressure, double temperature, double radius) const
        double get_thermal_expansion() const
        double get_reference_temperature() const
        double calc_static_shear_modulus(double radius) const
        double calc_static_bulk_modulus(double radius) const
        double calc_shear_viscosity(double radius) const
        double calc_bulk_viscosity(double radius) const

    cdef cppclass c_MaterialEOSConfig:
        double reference_density
        double reference_bulk_modulus
        double bulk_modulus_derivative
        double thermal_expansion
        double reference_temperature
        double invert_rtol
        int    invert_max_iters
        vector[double] radius
        vector[double] density
        vector[double] shear_modulus
        vector[double] bulk_modulus
        vector[double] shear_viscosity
        vector[double] bulk_viscosity

    cdef cppclass c_ConstantDensityEOS(c_MaterialEOSBase):
        c_ConstantDensityEOS() except +
        c_ConstantDensityEOS(const c_MaterialEOSConfig& cfg) except +
        double get_reference_density() const

    cdef cppclass c_BirchMurnaghanEOS(c_MaterialEOSBase):
        c_BirchMurnaghanEOS() except +
        c_BirchMurnaghanEOS(const c_MaterialEOSConfig& cfg) except +
        double get_reference_density()       const
        double get_reference_bulk_modulus()  const
        double get_bulk_modulus_derivative() const
        double get_invert_rtol()             const
        int    get_invert_max_iters()        const

    cdef cppclass c_VinetEOS(c_MaterialEOSBase):
        c_VinetEOS() except +
        c_VinetEOS(const c_MaterialEOSConfig& cfg) except +
        double get_reference_density()       const
        double get_reference_bulk_modulus()  const
        double get_bulk_modulus_derivative() const
        double get_invert_rtol()             const
        int    get_invert_max_iters()        const

    cdef cppclass c_InterpolatedEOS(c_MaterialEOSBase):
        c_InterpolatedEOS() except +
        c_InterpolatedEOS(const c_MaterialEOSConfig& cfg) except +
        size_t get_num_points() const
        cpp_bool has_shear_modulus() const
        cpp_bool has_bulk_modulus() const
        cpp_bool has_shear_viscosity() const
        cpp_bool has_bulk_viscosity() const

    double eos_bm_pressure(double eta, double K0, double K0_prime)
    double eos_vinet_pressure(double eta, double K0, double K0_prime)

    cdef enum class c_MaterialEOSModel:
        Constant
        BirchMurnaghan
        Vinet
        Interpolated

    c_MaterialEOSModel c_material_eos_model_from_name(const string& model_name) except +
    unique_ptr[c_MaterialEOSBase] c_find_material_eos(
        c_MaterialEOSModel model, const c_MaterialEOSConfig& cfg) except +


# =====================================================================================================================
# Cython wrapper class declarations
# =====================================================================================================================
cdef class MaterialEOSBase(PhysicsBase):
    cdef unique_ptr[c_MaterialEOSBase] _eos_ptr   # owns the most-derived C++ model; the typed pointers below do not


cdef class ConstantDensityEOS(MaterialEOSBase):
    cdef c_ConstantDensityEOS* _constant_ptr


cdef class BirchMurnaghanEOS(MaterialEOSBase):
    cdef c_BirchMurnaghanEOS* _bm_ptr


cdef class VinetEOS(MaterialEOSBase):
    cdef c_VinetEOS* _vinet_ptr


cdef class InterpolatedEOS(MaterialEOSBase):
    cdef c_InterpolatedEOS* _interp_ptr
