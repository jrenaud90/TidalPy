# distutils: language = c++
"""
luminosity.pxd
Cython declarations for TidalPy's stellar luminosity model hierarchy.

Exports the C++ luminosity models, the config struct, the enum factory, and the Python wrapper
classes so other extensions (e.g. the star world) can cimport and build/attach luminosity models.

Usage::

    from TidalPy.stellar_x.luminosity cimport (
        LuminosityBase, MassToLuminosity,
        c_LuminosityBase, c_MassToLuminosity, c_LuminosityConfig)
"""

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_PhysicsBase


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "luminosity_base_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_LuminosityBase(c_PhysicsBase):
        double calc_luminosity(double mass_kg) const
        double calc_luminosity_from_temperature(double temperature_k, double radius_m) const
        double calc_temperature_from_luminosity(double luminosity_w, double radius_m) const
        double calc_effective_temperature(double mass_kg, double radius_m) const
        void calc_luminosity_vectorize_mass(
            const vector[double]& mass_kg,
            vector[double]& out_luminosity) except +


cdef extern from "luminosity_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_LuminosityConfig:
        double luminosity_w
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

    # Enum naming each luminosity model.
    cdef enum class c_LuminosityModel:
        Fixed
        MassToLuminosity
        PowerLaw

    # Map a name/alias to the enum (raises ValueError on unknown name).
    c_LuminosityModel c_luminosity_model_from_name(const string& model_name) except +

    # Build the model named by the enum; returns an owning unique_ptr.
    unique_ptr[c_LuminosityBase] c_find_luminosity(
        c_LuminosityModel model, const c_LuminosityConfig& config) except +


# =====================================================================================================================
# Cython wrapper class declarations
# =====================================================================================================================
cdef class LuminosityBase(PhysicsBase):
    cdef unique_ptr[c_LuminosityBase] _luminosity_ptr   # owns the most-derived C++ model object
    cpdef dict get_config_dict(self)


cdef class FixedLuminosity(LuminosityBase):
    cdef c_FixedLuminosity* _fixed_ptr   # non-owning; ownership via LuminosityBase._luminosity_ptr
    cpdef dict get_config_dict(self)


cdef class MassToLuminosity(LuminosityBase):
    pass


cdef class PowerLawLuminosity(LuminosityBase):
    cdef c_PowerLawLuminosity* _power_law_ptr   # non-owning; ownership via LuminosityBase._luminosity_ptr
    cpdef dict get_config_dict(self)
