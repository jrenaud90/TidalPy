from TidalPy.Tides_x.mode_func_common cimport c_ModeFuncOutput
from TidalPy.Tides_x.obliquity.obliquity_common cimport c_ObliquitySeriesTable, c_ObliquityGeneralTable

cdef extern from "obliquity_driver_.hpp" nogil:
    c_ObliquitySeriesTable c_obliquity_series_table(
        int* error_code_ptr,
        int degree_l,
        int truncation)
    c_ObliquityGeneralTable c_obliquity_general_table(
        int* error_code_ptr,
        int degree_l)
    c_ModeFuncOutput c_obliquity_func(
        int* error_code_ptr,
        double obliquity,
        int degree_l,
        int truncation) except +
    c_ModeFuncOutput c_obliquity_squared_func(
        int* error_code_ptr,
        double obliquity,
        int degree_l,
        int truncation) except +
