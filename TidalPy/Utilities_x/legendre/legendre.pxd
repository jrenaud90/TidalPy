# distutils: language = c++
"""
legendre.pxd
Cython declarations for TidalPy's associated-Legendre utilities (Utilities_x/legendre).

Exposes the precomputed table driver (l = 2..10) and the generic xsf-backed evaluator so other
extensions can cimport and evaluate P_lm(cos theta) with its first/second colatitude derivatives.
"""


cdef extern from "legendre_common_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_LegendreValue:
        double p
        double dp_dtheta
        double d2p_dtheta2


cdef extern from "legendre_driver_.hpp" namespace "tidalpy" nogil:
    c_LegendreValue c_legendre(int degree_l, int order_m, double colatitude)


cdef extern from "legendre_generic_.hpp" namespace "tidalpy" nogil:
    c_LegendreValue c_legendre_generic(int degree_l, int order_m, double colatitude)
