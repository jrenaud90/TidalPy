# distutils: language = c++

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

from TidalPy.Utilities.classes.classes cimport PhysicsBase, c_PhysicsBase, c_ParamMap, c_ThermoPoint


cdef extern from "eos_law_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_EOSPoint:
        double density
        double bulk_modulus
        double adiabatic_bulk_modulus
        double thermal_expansion

    cdef cppclass c_EOSBase(c_PhysicsBase):
        void calc_eos(const c_ThermoPoint& point, cpp_bool thermal, c_EOSPoint& out) const
        double calc_density(const c_ThermoPoint& point, cpp_bool thermal) const
        void calc_eos_vectorize(
            const vector[double]& pressure,
            const vector[double]& temperature,
            const vector[double]& radius,
            cpp_bool thermal,
            vector[double]& out_density,
            vector[double]& out_bulk_modulus,
            vector[double]& out_adiabatic_bulk_modulus,
            vector[double]& out_thermal_expansion) except +

    unique_ptr[c_EOSBase] c_find_eos(const string& model_name, const c_ParamMap& params) except +
    string c_eos_canonical_name(const string& model_name) except +


cdef extern from "shear_modulus_law_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_ShearModulusBase(c_PhysicsBase):
        double calc_shear_modulus(const c_ThermoPoint& point) const
        void calc_shear_modulus_vectorize(
            const vector[double]& pressure,
            const vector[double]& temperature,
            const vector[double]& radius,
            vector[double]& out_shear_modulus) except +

    unique_ptr[c_ShearModulusBase] c_find_shear_modulus(const string& model_name, const c_ParamMap& params) except +
    string c_shear_modulus_canonical_name(const string& model_name) except +


cdef class EOSBase(PhysicsBase):
    cdef c_EOSBase* _eos(self) except NULL


cdef class ShearModulusBase(PhysicsBase):
    cdef c_ShearModulusBase* _shear_modulus(self) except NULL
