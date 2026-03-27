from libcpp.pair cimport pair

from TidalPy.Tides_x.eccentricity.eccentricity_common cimport EccentricityFuncOutput
from TidalPy.utilities.lookups cimport c_IntMap, c_Key3, c_Key2, c_Key1

cdef extern from "eccentricity_func_l9_.hpp" nogil:
    ctypedef pair[c_IntMap[c_Key3, double], c_IntMap[c_Key2, c_IntMap[c_Key1, double]]] EccentricityFuncOutput

    EccentricityFuncOutput c_eccentricity_function_l9_e1(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l9_e2(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l9_e3(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l9_e4(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l9_e5(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l9_e10(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l9_e15(double eccentricity)
    EccentricityFuncOutput c_eccentricity_function_l9_e20(double eccentricity)
