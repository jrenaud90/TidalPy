from TidalPy.Tides_x.eccentricity.eccentricity_common cimport EccentricityFuncOutput

cdef extern from "eccentricity_func_l7_.hpp" nogil:

    EccentricityFuncOutput c_eccentricity_function_l7_e1(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l7_e2(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l7_e3(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l7_e4(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l7_e5(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l7_e10(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l7_e15(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l7_e20(double eccentricity)
