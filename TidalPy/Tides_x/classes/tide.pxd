# distutils: language = c++
"""
tide.pxd
Cython declarations for TidalPy's global (1D) tidal dissipation model hierarchy.

Exports the C++ models, the per-degree config struct, the enum factory, and the Python
wrapper classes so other extensions (the global-potential collapse, the world) can
cimport and build or attach a tide model at C speed.
"""

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector
from libcpp cimport bool as cpp_bool

from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_PhysicsBase
from TidalPy.Tides_x.love.love cimport c_LoveNumbers


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "tide_base_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_TideBase(c_PhysicsBase):
        c_LoveNumbers calc_love_numbers(int degree_l, double frequency, const c_LoveNumbers& solver_love) const
        double calc_neg_imk(int degree_l, double frequency, const c_LoveNumbers& solver_love) const
        cpp_bool needs_radial_solve() const


cdef extern from "tide_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_TideModelConfig:
        vector[double] fixed_k
        vector[double] fixed_q
        vector[double] fixed_dt

    cdef cppclass c_RheologyTide(c_TideBase):
        c_RheologyTide() except +
        c_RheologyTide(const c_TideModelConfig& cfg) except +

    cdef cppclass c_FixedQTide(c_TideBase):
        c_FixedQTide() except +
        c_FixedQTide(const c_TideModelConfig& cfg) except +
        double get_fixed_k(int degree_l) const
        double get_fixed_q(int degree_l) const

    cdef cppclass c_FixedLagTide(c_TideBase):
        c_FixedLagTide() except +
        c_FixedLagTide(const c_TideModelConfig& cfg) except +
        double get_fixed_k(int degree_l) const
        double get_fixed_dt(int degree_l) const

    cdef cppclass c_CTLQTide(c_TideBase):
        c_CTLQTide() except +
        c_CTLQTide(const c_TideModelConfig& cfg) except +
        double get_fixed_k(int degree_l) const
        double get_fixed_dt(int degree_l) const
        double get_fixed_q(int degree_l) const

    cdef enum class c_TideModel:
        Rheology
        FixedQ
        FixedLag
        CTLQ

    c_TideModel c_tide_model_from_name(const string& model_name) except +
    unique_ptr[c_TideBase] c_find_tide(c_TideModel model, const c_TideModelConfig& cfg) except +


# =====================================================================================================================
# Cython wrapper class declarations
# =====================================================================================================================
cdef class TideBase(PhysicsBase):
    cdef unique_ptr[c_TideBase] _tide_ptr   # owns the most-derived C++ model object
    cpdef dict get_config_dict(self)


cdef class RheologyTide(TideBase):
    cdef c_RheologyTide* _rheology_ptr       # non-owning; ownership via TideBase._tide_ptr


cdef class FixedQTide(TideBase):
    cdef c_FixedQTide* _fixedq_ptr           # non-owning; ownership via TideBase._tide_ptr
    cpdef dict get_config_dict(self)


cdef class FixedLagTide(TideBase):
    cdef c_FixedLagTide* _fixedlag_ptr       # non-owning; ownership via TideBase._tide_ptr
    cpdef dict get_config_dict(self)


cdef class CTLQTide(TideBase):
    cdef c_CTLQTide* _ctlq_ptr               # non-owning; ownership via TideBase._tide_ptr
    cpdef dict get_config_dict(self)
