# distutils: language = c++
"""Cython declarations for TidalPy's global (1D) tide models: the C++ base class, the name-based factory, and the
Python wrapper classes.
"""

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector
from libcpp cimport bool as cpp_bool

from TidalPy.Utilities.classes.classes cimport PhysicsBase, c_PhysicsBase, c_ParamMap
from TidalPy.Tides.love.love cimport c_LoveNumbers


cdef extern from "tide_base_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_TideBase(c_PhysicsBase):
        c_LoveNumbers calc_love_numbers(int degree_l, double frequency, const c_LoveNumbers& solver_love) const
        double calc_neg_imk(int degree_l, double frequency, const c_LoveNumbers& solver_love) const
        cpp_bool needs_radial_solve() const
        double get_fixed_k(int degree_l) const
        double get_fixed_q(int degree_l) const
        double get_fixed_dt(int degree_l) const


cdef extern from "tide_.hpp" namespace "tidalpy" nogil:

    unique_ptr[c_TideBase] c_find_tide(const string& model_name, const c_ParamMap& params) except +
    string c_tide_canonical_name(const string& model_name) except +
    vector[string] c_tide_model_names() except +


cdef class TideBase(PhysicsBase):
    cdef c_TideBase* _tide(self) except NULL


cdef class RheologyTide(TideBase):
    pass


cdef class FixedQTide(TideBase):
    pass


cdef class FixedLagTide(TideBase):
    pass


cdef class CTLQTide(TideBase):
    pass
