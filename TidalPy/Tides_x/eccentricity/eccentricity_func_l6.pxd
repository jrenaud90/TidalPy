from TidalPy.Tides_x.eccentricity.eccentricity_common cimport EccentricityFuncOutput

cdef extern from "eccentricity_func_l6_.hpp" nogil:

    EccentricityFuncOutput c_eccentricity_function_l6_e1(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l6_e2(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l6_e3(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l6_e4(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l6_e5(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l6_e10(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l6_e15(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l6_e20(double eccentricity)
