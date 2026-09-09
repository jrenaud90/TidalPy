# distutils: language = c++
"""
viscosity.pxd
Cython declarations for TidalPy's viscosity model hierarchy.

Exports the C++ models, the combined config struct, the enum factory, and the
Python wrapper classes so other extensions (layers, worlds) can cimport and build
or attach viscosity models at C speed.
"""

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp cimport bool as cpp_bool

from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_PhysicsBase


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "viscosity_base_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_ViscosityBase(c_PhysicsBase):
        double calc_viscosity(double temperature_k, double pressure_pa) const


cdef extern from "viscosity_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_ViscosityConfig:
        double reference_viscosity
        double reference_temperature
        double molar_activation_energy
        double molar_activation_volume
        double arrhenius_coeff
        double stress
        double stress_expo
        double grain_size
        double grain_size_expo
        cpp_bool additional_temp_dependence

    cdef cppclass c_ConstantViscosity(c_ViscosityBase):
        c_ConstantViscosity() except +
        c_ConstantViscosity(const c_ViscosityConfig& cfg) except +
        double get_reference_viscosity() const

    cdef cppclass c_ReferenceViscosity(c_ViscosityBase):
        c_ReferenceViscosity() except +
        c_ReferenceViscosity(const c_ViscosityConfig& cfg) except +
        double get_reference_viscosity() const
        double get_reference_temperature() const
        double get_molar_activation_energy() const
        double get_molar_activation_volume() const

    cdef cppclass c_ArrheniusViscosity(c_ViscosityBase):
        c_ArrheniusViscosity() except +
        c_ArrheniusViscosity(const c_ViscosityConfig& cfg) except +
        double get_arrhenius_coeff() const
        double get_stress() const
        double get_stress_expo() const
        double get_grain_size() const
        double get_grain_size_expo() const
        double get_molar_activation_energy() const
        double get_molar_activation_volume() const
        cpp_bool get_additional_temp_dependence() const

    cdef enum class c_ViscosityModel:
        Arrhenius
        Reference
        Constant

    c_ViscosityModel c_viscosity_model_from_name(const string& model_name) except +
    unique_ptr[c_ViscosityBase] c_find_viscosity(
        c_ViscosityModel model, const c_ViscosityConfig& cfg) except +


# =====================================================================================================================
# Cython wrapper class declarations
# =====================================================================================================================
cdef class ViscosityBase(PhysicsBase):
    cdef unique_ptr[c_ViscosityBase] _visc_ptr   # owns the most-derived C++ model object
    cpdef dict get_config_dict(self)


cdef class ConstantViscosity(ViscosityBase):
    cdef c_ConstantViscosity* _constant_ptr      # non-owning; ownership via ViscosityBase._visc_ptr
    cpdef dict get_config_dict(self)


cdef class ReferenceViscosity(ViscosityBase):
    cdef c_ReferenceViscosity* _ref_ptr          # non-owning; ownership via ViscosityBase._visc_ptr
    cpdef dict get_config_dict(self)


cdef class ArrheniusViscosity(ViscosityBase):
    cdef c_ArrheniusViscosity* _arr_ptr          # non-owning; ownership via ViscosityBase._visc_ptr
    cpdef dict get_config_dict(self)
