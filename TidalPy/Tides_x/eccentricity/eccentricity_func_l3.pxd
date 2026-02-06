from TidalPy.Tides_x.eccentricity.eccentricity_common cimport EccentricityFuncOutput

cdef extern from "eccentricity_func_l3_.hpp" nogil:

    EccentricityFuncOutput c_eccentricity_function_l3_e1(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l3_e2(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l3_e3(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l3_e4(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l3_e5(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l3_e10(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l3_e15(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l3_e20(double eccentricity)
