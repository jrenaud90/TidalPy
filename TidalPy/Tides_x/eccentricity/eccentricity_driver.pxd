from TidalPy.Tides_x.eccentricity.eccentricity_common cimport EccentricityFuncOutput, c_EccentricitySeriesTable

cdef extern from "eccentricity_driver_.hpp" nogil:
    c_EccentricitySeriesTable c_eccentricity_series_table(
        int* error_code_ptr,
        int degree_l,
        int truncation)
    EccentricityFuncOutput c_eccentricity_func(
        int* error_code_ptr,
        double eccentricity,
        int degree_l,
        int truncation)
    EccentricityFuncOutput c_eccentricity_squared_func(
        int* error_code_ptr,
        double eccentricity,
        int degree_l,
        int truncation)
