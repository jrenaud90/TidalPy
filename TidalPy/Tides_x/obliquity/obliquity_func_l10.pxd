from TidalPy.Tides_x.obliquity.obliquity_common cimport ObliquityFuncOutput

cdef extern from "obliquity_func_l10_.hpp" nogil:

    ObliquityFuncOutput c_obliquity_function_l10_gen(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l10_2(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l10_4(double obliquity)
    ObliquityFuncOutput c_obliquity_function_l10_off(double obliquity)
