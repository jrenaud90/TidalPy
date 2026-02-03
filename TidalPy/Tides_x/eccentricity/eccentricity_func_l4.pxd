from TidalPy.Tides_x.eccentricity.eccentricity_common cimport EccentricityFuncOutput

cdef extern from "eccentricity_func_l4_.hpp" nogil:

    EccentricityFuncOutput c_eccentricity_function_l4_e1(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l4_e2(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l4_e3(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l4_e4(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l4_e5(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l4_e10(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l4_e15(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l4_e20(double eccentricity)
