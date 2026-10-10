# distutils: language = c++

cdef extern from "legendre_common_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_LegendreValue:
        double p
        double dp_dtheta
        double d2p_dtheta2


cdef extern from "legendre_driver_.hpp" namespace "tidalpy" nogil:
    c_LegendreValue c_legendre(int degree_l, int order_m, double colatitude)


cdef extern from "legendre_generic_.hpp" namespace "tidalpy" nogil:
    c_LegendreValue c_legendre_generic(int degree_l, int order_m, double colatitude)
