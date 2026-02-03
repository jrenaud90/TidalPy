from TidalPy.Tides_x.obliquity.obliquity_common cimport ObliquityFuncOutput

cdef extern from "obliquity_func_l2_.hpp" nogil:

    ObliquityFuncOutput c_obliquity_function_l2_gen(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l2_2(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l2_4(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l2_off(double obliquity)
