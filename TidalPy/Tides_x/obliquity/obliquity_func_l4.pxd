from TidalPy.Tides_x.obliquity.obliquity_common cimport ObliquityFuncOutput

cdef extern from "obliquity_func_l4_.hpp" nogil:

    ObliquityFuncOutput c_obliquity_function_l4_gen(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l4_2(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l4_4(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l4_off(double obliquity)
