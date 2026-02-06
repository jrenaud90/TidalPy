from TidalPy.Tides_x.obliquity.obliquity_common cimport ObliquityFuncOutput

cdef extern from "obliquity_func_l5_.hpp" nogil:

    ObliquityFuncOutput c_obliquity_function_l5_gen(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l5_2(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l5_4(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l5_off(double obliquity)
