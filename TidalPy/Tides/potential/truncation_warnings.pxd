from libcpp.string cimport string

cdef extern from "truncation_warnings_.hpp" nogil:
    void c_warn_standalone_tide_truncations(
        const string& function_name,
        double eccentricity,
        double obliquity,
        int eccentricity_truncation,
        int obliquity_truncation,
        int max_degree_l,
        double spin_ratio) except +
