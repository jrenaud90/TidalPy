from libcpp cimport bool as cpp_bool


cdef extern from "eccentricity_accuracy_.hpp" nogil:
    const int C_ECCENTRICITY_EXACT
    const double C_ECCENTRICITY_EXACT_TOLERANCE
    double c_eccentricity_accuracy_limit(
        int eccentricity_truncation,
        double tolerance,
        int max_degree_l,
        double spin_ratio)
    double c_eccentricity_truncation_limit(int eccentricity_truncation, int max_degree_l, double spin_ratio)
    int c_recommend_eccentricity_truncation(
        double eccentricity,
        double tolerance,
        int max_degree_l,
        double spin_ratio)


cdef extern from "eccentricity_common_.hpp" nogil:
    const int C_NUM_ECCENTRICITY_TRUNCATIONS
    const int* C_ECCENTRICITY_TRUNCATIONS

    cdef cppclass c_EccentricitySeriesTable:
        int degree_l
        int truncation
        cpp_bool valid()
