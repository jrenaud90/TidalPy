# distutils: language = c++

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

from TidalPy.Utilities.classes.classes cimport PhysicsBase, c_PhysicsBase, c_ParamMap


cdef extern from "luminosity_base_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_LuminosityBase(c_PhysicsBase):
        double calc_luminosity(double mass) const
        double calc_luminosity_from_temperature(double temperature, double radius) const
        double calc_temperature_from_luminosity(double luminosity, double radius) const
        double calc_effective_temperature(double mass, double radius) const
        void calc_luminosity_vectorize_mass(
            const vector[double]& mass,
            vector[double]& out_luminosity) except +


cdef extern from "luminosity_.hpp" namespace "tidalpy" nogil:

    unique_ptr[c_LuminosityBase] c_find_luminosity(const string& model_name, const c_ParamMap& params) except +
    string c_luminosity_canonical_name(const string& model_name) except +
    vector[string] c_luminosity_model_names() except +


cdef class LuminosityBase(PhysicsBase):
    cdef c_LuminosityBase* _luminosity(self) except NULL


cdef class FixedLuminosity(LuminosityBase):
    pass


cdef class MassToLuminosity(LuminosityBase):
    pass


cdef class PowerLawLuminosity(LuminosityBase):
    pass
