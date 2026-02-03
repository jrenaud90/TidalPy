from TidalPy.Tides_x.obliquity.obliquity_common cimport ObliquityFuncOutput

cdef extern from "obliquity_func_l6_.hpp" nogil:

    ObliquityFuncOutput c_obliquity_function_l6_gen(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l6_2(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l6_4(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l6_off(double obliquity)
