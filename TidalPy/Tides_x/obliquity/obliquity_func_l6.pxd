from libcpp.pair cimport pair

from TidalPy.utilities.lookups cimport c_IntMap, c_Key3, c_Key2, c_Key1
from TidalPy.Tides_x.obliquity.obliquity_common cimport ObliquityFuncOutput

cdef extern from "obliquity_func_l6_.hpp" nogil:
    ctypedef pair[c_IntMap[c_Key3, double], c_IntMap[c_Key2, c_IntMap[c_Key1, double]]] ObliquityFuncOutput

    ObliquityFuncOutput c_obliquity_function_l6_gen(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l6_2(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l6_4(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l6_off(double obliquity)
