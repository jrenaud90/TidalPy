from libcpp.pair cimport pair

from TidalPy.utilities.lookups cimport c_IntMap, c_Key3, c_Key2, c_Key1
from TidalPy.Tides_x.obliquity.obliquity_common cimport ObliquityFuncOutput

cdef extern from "obliquity_func_l4_.hpp" nogil:
    ctypedef pair[c_IntMap[c_Key3, double], c_IntMap[c_Key2, c_IntMap[c_Key1, double]]] ObliquityFuncOutput

    ObliquityFuncOutput c_obliquity_function_l4_gen(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l4_2(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l4_4(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l4_off(double obliquity)
