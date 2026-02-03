from TidalPy.Tides_x.obliquity.obliquity_common cimport ObliquityFuncOutput

cdef extern from "obliquity_func_l3_.hpp" nogil:

    ObliquityFuncOutput c_obliquity_function_l3_gen(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l3_2(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l3_4(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l3_off(double obliquity)
