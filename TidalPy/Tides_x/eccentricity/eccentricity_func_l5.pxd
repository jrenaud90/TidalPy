from TidalPy.Tides_x.eccentricity.eccentricity_common cimport c_EccentricitySeriesTable

cdef extern from "eccentricity_func_l5_.hpp" nogil:
    c_EccentricitySeriesTable c_eccentricity_series_l5(int truncation)
