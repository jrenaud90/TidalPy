# distutils: language = c++
"""
love.pxd
Cython declarations for TidalPy's Love numbers container.

Exports c_LoveNumbers and the Python wrapper LoveNumbers so other extensions
can cimport and work with Love numbers at C speed.

Usage::

    from TidalPy.Tides_x.love.love cimport LoveNumbers, c_LoveNumbers
"""

from libcpp.complex cimport complex as cpp_complex
from libcpp.string cimport string


# =====================================================================================================================
# C++ struct declaration
# =====================================================================================================================
cdef extern from "love_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_LoveNumbers:
        cpp_complex[double] k   # potential Love number           [dimensionless]
        cpp_complex[double] h   # radial displacement Love number [dimensionless]
        cpp_complex[double] l   # tangential displacement Love number [dimensionless]
        c_LoveNumbers() except +
        c_LoveNumbers(cpp_complex[double], cpp_complex[double], cpp_complex[double]) except +


# =====================================================================================================================
# Love-number methods and the homogeneous-sphere formulas (love_method_.hpp)
# =====================================================================================================================
cdef extern from "love_method_.hpp" namespace "tidalpy" nogil:

    # Method names <-> c_LoveMethod index (0 radial_solver, 1 propagation_matrix, 2 homogeneous, 3 cpl, 4 ctl,
    # 5 laterally_inhomogeneous). Unknown names / indices throw std::invalid_argument (ValueError).
    int    c_parse_love_method_int(const string& name) except +
    string c_love_method_name_int(int value) except +

    double c_calc_effective_rigidity_real(
        double shear_modulus, double density, double gravity, double radius, int degree_l) except +
    cpp_complex[double] c_calc_effective_rigidity_complex(
        cpp_complex[double] shear_modulus, double density, double gravity, double radius,
        int degree_l) except +
    c_LoveNumbers c_calc_homogeneous_love_numbers(
        cpp_complex[double] complex_shear_modulus, double density, double gravity, double radius,
        int degree_l) except +
    c_LoveNumbers c_apply_fixed_q(const c_LoveNumbers& love, double fixed_q) except +
    c_LoveNumbers c_apply_fixed_dt(const c_LoveNumbers& love, double frequency, double fixed_dt) except +


# =====================================================================================================================
# Cython wrapper class declaration
# =====================================================================================================================
cdef class LoveNumbers:
    cdef c_LoveNumbers _love
    cpdef dict to_dict(self)
