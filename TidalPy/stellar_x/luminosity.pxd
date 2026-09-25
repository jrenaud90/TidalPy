# distutils: language = c++

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_PhysicsBase


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

    cdef cppclass c_LuminosityConfig:
        double luminosity
        double power_law_coeff
        double power_law_exponent

    cdef cppclass c_FixedLuminosity(c_LuminosityBase):
        c_FixedLuminosity() except +
        c_FixedLuminosity(const c_LuminosityConfig& config) except +
        double get_luminosity() const

    cdef cppclass c_MassToLuminosity(c_LuminosityBase):
        c_MassToLuminosity() except +
        c_MassToLuminosity(const c_LuminosityConfig& config) except +

    cdef cppclass c_PowerLawLuminosity(c_LuminosityBase):
        c_PowerLawLuminosity() except +
        c_PowerLawLuminosity(const c_LuminosityConfig& config) except +
        double get_coeff()    const
        double get_exponent() const

    cdef enum class c_LuminosityModel:
        Fixed
        MassToLuminosity
        PowerLaw

    # Raises ValueError on an unknown name.
    c_LuminosityModel c_luminosity_model_from_name(const string& model_name) except +

    unique_ptr[c_LuminosityBase] c_find_luminosity(
        c_LuminosityModel model, const c_LuminosityConfig& config) except +


cdef class LuminosityBase(PhysicsBase):
    cdef unique_ptr[c_LuminosityBase] _luminosity_ptr   # owns the most-derived C++ model object
    cdef void _adopt(self, unique_ptr[c_LuminosityBase]& model) noexcept


cdef class FixedLuminosity(LuminosityBase):
    pass


cdef class MassToLuminosity(LuminosityBase):
    pass


cdef class PowerLawLuminosity(LuminosityBase):
    pass
