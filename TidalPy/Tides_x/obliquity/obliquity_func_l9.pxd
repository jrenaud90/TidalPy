from TidalPy.Tides_x.obliquity.obliquity_common cimport ObliquityFuncOutput

cdef extern from "obliquity_func_l9_.hpp" nogil:

    ObliquityFuncOutput c_obliquity_function_l9_gen(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l9_2(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l9_4(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l9_off(double obliquity)
