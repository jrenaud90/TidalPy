# distutils: language = c++

cdef extern from "conversions_.hpp" namespace "tidalpy" nogil:
    double c_semi_a2orbital_motion(double semi_major_axis, double gravitational_parameter)
    double c_orbital_motion2semi_a(double orbital_motion, double gravitational_parameter)
