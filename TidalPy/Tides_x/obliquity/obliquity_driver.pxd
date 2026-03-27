from libcpp.pair cimport pair

from TidalPy.Tides_x.obliquity.obliquity_common cimport ObliquityFuncOutput

cdef extern from "obliquity_driver_.hpp" nogil:
    ObliquityFuncOutput c_obliquity_func(
        int* error_code_ptr,
        double obliquity,
        int degree_l,
        int truncation)
