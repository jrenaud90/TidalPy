# distutils: language = c++

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector
from libcpp cimport bool as cpp_bool

from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_PhysicsBase


cdef extern from "partial_melt_base_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_PartialMeltInputs:
        double temperature
        double premelt_viscosity
        double premelt_shear

    cdef cppclass c_PartialMeltResult:
        double melt_fraction
        double postmelt_viscosity
        double postmelt_shear_modulus

    cdef cppclass c_PartialMeltBase(c_PhysicsBase):
        double get_solidus() const
        double get_liquidus() const
        double get_liquid_shear() const
        double get_liquid_viscosity() const
        cpp_bool get_bulk_melt_weakening() const
        double get_liquid_bulk_modulus() const
        double calc_melt_fraction(double temperature) const
        c_PartialMeltResult calc_partial_melt(const c_PartialMeltInputs& inputs) const
        double calc_bulk_modulus_melt(double temperature, double premelt_bulk, double framework_shear) const


cdef extern from "partial_melt_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_PartialMeltConfig:
        double solidus
        double liquidus
        double liquid_shear
        double liquid_viscosity
        cpp_bool bulk_melt_weakening
        double liquid_bulk_modulus
        double fs_visc_power_slope
        double fs_visc_log10_at_solidus
        double fs_shear_power_slope
        double fs_shear_log10_at_solidus
        double crit_melt_frac
        double crit_melt_frac_width
        double hn_visc_slope_1
        double hn_visc_falloff_slope
        double hn_shear_param_1
        double hn_shear_falloff_slope

    cdef cppclass c_OffPartialMelt(c_PartialMeltBase):
        c_OffPartialMelt() except +
        c_OffPartialMelt(const c_PartialMeltConfig& cfg) except +

    cdef cppclass c_SpohnPartialMelt(c_PartialMeltBase):
        c_SpohnPartialMelt() except +
        c_SpohnPartialMelt(const c_PartialMeltConfig& cfg) except +
        double get_visc_power_slope() const
        double get_visc_log10_at_solidus() const
        double get_shear_power_slope() const
        double get_shear_log10_at_solidus() const

    cdef cppclass c_HenningPartialMelt(c_PartialMeltBase):
        c_HenningPartialMelt() except +
        c_HenningPartialMelt(const c_PartialMeltConfig& cfg) except +
        double get_crit_melt_frac() const
        double get_crit_melt_frac_width() const
        double get_visc_slope_1() const
        double get_visc_falloff_slope() const
        double get_shear_param_1() const
        double get_shear_falloff_slope() const

    cdef enum class c_PartialMeltModel:
        Off
        Spohn
        Henning

    c_PartialMeltModel c_partial_melt_model_from_name(const string& model_name) except +
    unique_ptr[c_PartialMeltBase] c_find_partial_melt(
        c_PartialMeltModel model, const c_PartialMeltConfig& cfg) except +


cdef class PartialMeltBase(PhysicsBase):
    cdef unique_ptr[c_PartialMeltBase] _melt_ptr   # owns the most-derived C++ model object
    cdef void _adopt(self, unique_ptr[c_PartialMeltBase]& model) noexcept


cdef class OffPartialMelt(PartialMeltBase):
    pass


cdef class SpohnPartialMelt(PartialMeltBase):
    pass


cdef class HenningPartialMelt(PartialMeltBase):
    pass
