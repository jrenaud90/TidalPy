from TidalPy.Tides_x.obliquity.obliquity_common cimport c_ObliquitySeriesTable, c_ObliquityGeneralTable

cdef extern from "obliquity_func_l10_.hpp" nogil:
    c_ObliquitySeriesTable c_obliquity_series_l10(int truncation)
    c_ObliquityGeneralTable c_obliquity_general_l10()
