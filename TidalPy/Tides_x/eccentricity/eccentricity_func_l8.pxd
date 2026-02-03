from TidalPy.Tides_x.eccentricity.eccentricity_common cimport EccentricityFuncOutput

cdef extern from "eccentricity_func_l8_.hpp" nogil:

    EccentricityFuncOutput c_eccentricity_function_l8_e1(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l8_e2(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l8_e3(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l8_e4(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l8_e5(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l8_e10(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l8_e15(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l8_e20(double eccentricity)
