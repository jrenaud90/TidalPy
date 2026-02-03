from libcpp.pair cimport pair

from TidalPy.Tides_x.eccentricity.eccentricity_common cimport EccentricityFuncOutput

cdef extern from "eccentricity_driver_.hpp" nogil:
    EccentricityFuncOutput c_eccentricity_func(
        int* error_code_ptr,
        double eccentricity,
        int degree_l,
        int truncation)
