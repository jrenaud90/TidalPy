# distutils: language = c++
"""
partial_melt.pxd
Cython declarations for TidalPy's partial-melt model hierarchy.

Exports the C++ models, the combined config + input/result structs, the enum
factory, and the Python wrapper classes so other extensions (layers, worlds) can
cimport and build or attach partial-melt models at C speed.
"""

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector
from libcpp cimport bool as cpp_bool

from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_PhysicsBase


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "partial_melt_base_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_PartialMeltInputs:
        double temperature_k
        double premelt_viscosity
        double premelt_shear
        double liquid_viscosity

    cdef cppclass c_PartialMeltResult:
        double melt_fraction
        double postmelt_viscosity
        double postmelt_shear_modulus

    cdef cppclass c_PartialMeltBase(c_PhysicsBase):
        double get_solidus() const
        double get_liquidus() const
        double get_liquid_shear() const
        double calc_melt_fraction(double temperature_k) const
        c_PartialMeltResult calc_partial_melt(const c_PartialMeltInputs& inputs) const


cdef extern from "partial_melt_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_PartialMeltConfig:
        double solidus_k
        double liquidus_k
        double liquid_shear_pa
        double fs_visc_power_slope
        double fs_visc_power_phase
        double fs_shear_power_slope
        double fs_shear_power_phase
        double crit_melt_frac
        double crit_melt_frac_width
        double hn_visc_slope_1
        double hn_visc_falloff_slope
        double hn_shear_param_1
        double hn_shear_param_2
        double hn_shear_falloff_slope

    cdef cppclass c_OffPartialMelt(c_PartialMeltBase):
        c_OffPartialMelt() except +
        c_OffPartialMelt(const c_PartialMeltConfig& cfg) except +

    cdef cppclass c_SpohnPartialMelt(c_PartialMeltBase):
        c_SpohnPartialMelt() except +
        c_SpohnPartialMelt(const c_PartialMeltConfig& cfg) except +
        double get_visc_power_slope() const
        double get_visc_power_phase() const
        double get_shear_power_slope() const
        double get_shear_power_phase() const

    cdef cppclass c_HenningPartialMelt(c_PartialMeltBase):
        c_HenningPartialMelt() except +
        c_HenningPartialMelt(const c_PartialMeltConfig& cfg) except +
        double get_crit_melt_frac() const
        double get_crit_melt_frac_width() const
        double get_visc_slope_1() const
        double get_visc_falloff_slope() const
        double get_shear_param_1() const
        double get_shear_param_2() const
        double get_shear_falloff_slope() const

    cdef enum class c_PartialMeltModel:
        Off
        Spohn
        Henning

    c_PartialMeltModel c_partial_melt_model_from_name(const string& model_name) except +
    unique_ptr[c_PartialMeltBase] c_find_partial_melt(
        c_PartialMeltModel model, const c_PartialMeltConfig& cfg) except +


# =====================================================================================================================
# Cython wrapper class declarations
# =====================================================================================================================
cdef class PartialMeltBase(PhysicsBase):
    cdef unique_ptr[c_PartialMeltBase] _melt_ptr   # owns the most-derived C++ model object
    cpdef dict get_config_dict(self)


cdef class OffPartialMelt(PartialMeltBase):
    cdef c_OffPartialMelt* _off_ptr                # non-owning; ownership via PartialMeltBase._melt_ptr
    cpdef dict get_config_dict(self)


cdef class SpohnPartialMelt(PartialMeltBase):
    cdef c_SpohnPartialMelt* _spohn_ptr            # non-owning; ownership via PartialMeltBase._melt_ptr
    cpdef dict get_config_dict(self)


cdef class HenningPartialMelt(PartialMeltBase):
    cdef c_HenningPartialMelt* _henning_ptr        # non-owning; ownership via PartialMeltBase._melt_ptr
    cpdef dict get_config_dict(self)
