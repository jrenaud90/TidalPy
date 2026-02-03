from TidalPy.Tides_x.obliquity.obliquity_common cimport ObliquityFuncOutput

cdef extern from "obliquity_func_l8_.hpp" nogil:

    ObliquityFuncOutput c_obliquity_function_l8_gen(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l8_2(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l8_4(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l8_off(double obliquity)
