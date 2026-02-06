from TidalPy.Tides_x.obliquity.obliquity_common cimport ObliquityFuncOutput

cdef extern from "obliquity_func_l7_.hpp" nogil:

    ObliquityFuncOutput c_obliquity_function_l7_gen(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l7_2(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l7_4(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l7_off(double obliquity)
