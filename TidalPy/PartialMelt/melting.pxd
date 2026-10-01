# distutils: language = c++

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

from TidalPy.Utilities.classes.classes cimport PhysicsBase, c_PhysicsBase, c_ParamMap


cdef extern from "melting_curve_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_MeltingCurveBase(c_PhysicsBase):
        double calc_melting_temperature(double pressure) const
        void calc_melting_temperature_vectorize(
            const vector[double]& pressure,
            vector[double]& out_temperature) except +

    unique_ptr[c_MeltingCurveBase] c_find_melting_curve(const string& model_name, const c_ParamMap& params) except +
    string c_melting_curve_canonical_name(const string& model_name) except +


cdef extern from "melt_weakening_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_MeltWeakeningInputs:
        double temperature
        double solidus
        double liquidus
        double melt_fraction
        double solid_shear
        double solid_viscosity
        double liquid_shear
        double liquid_viscosity

    cdef cppclass c_MeltWeakeningResult:
        double shear_modulus
        double viscosity

    cdef cppclass c_MeltWeakeningBase(c_PhysicsBase):
        c_MeltWeakeningResult calc_weakening(const c_MeltWeakeningInputs& inputs) const

    unique_ptr[c_MeltWeakeningBase] c_find_melt_weakening(const string& model_name, const c_ParamMap& params) except +
    string c_melt_weakening_canonical_name(const string& model_name) except +


cdef extern from "melt_mixing_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_BulkModulusMixingBase(c_PhysicsBase):
        double calc_bulk_modulus(
            double solid_bulk_modulus,
            double liquid_bulk_modulus,
            double framework_shear_modulus,
            double melt_fraction) const

    cdef cppclass c_BulkViscosityMixingBase(c_PhysicsBase):
        double calc_bulk_viscosity(
            double solid_bulk_viscosity,
            double postmelt_shear_viscosity,
            double melt_fraction) const

    unique_ptr[c_BulkModulusMixingBase] c_find_bulk_modulus_mixing(
        const string& model_name, const c_ParamMap& params) except +
    string c_bulk_modulus_mixing_canonical_name(const string& model_name) except +
    unique_ptr[c_BulkViscosityMixingBase] c_find_bulk_viscosity_mixing(
        const string& model_name, const c_ParamMap& params) except +
    string c_bulk_viscosity_mixing_canonical_name(const string& model_name) except +


cdef class MeltingCurveBase(PhysicsBase):
    cdef c_MeltingCurveBase* _curve(self) except NULL


cdef class MeltWeakeningBase(PhysicsBase):
    cdef c_MeltWeakeningBase* _weakening(self) except NULL


cdef class BulkModulusMixingBase(PhysicsBase):
    cdef c_BulkModulusMixingBase* _mixing(self) except NULL


cdef class BulkViscosityMixingBase(PhysicsBase):
    cdef c_BulkViscosityMixingBase* _mixing(self) except NULL
