from TidalPy.Tides_x.obliquity.obliquity_common cimport c_ObliquitySeriesTable, c_ObliquityGeneralTable

cdef extern from "obliquity_func_l9_.hpp" nogil:
    c_ObliquitySeriesTable c_obliquity_series_l9(int truncation)
    c_ObliquityGeneralTable c_obliquity_general_l9()
