from libcpp cimport bool as cpp_bool
from libcpp.pair cimport pair

from TidalPy.Utilities_x.lookups cimport c_IntMap, c_Key1, c_Key2, c_Key3


cdef extern from "obliquity_accuracy_.hpp" nogil:
    const int C_OBLIQUITY_OFF
    const int C_OBLIQUITY_GENERAL
    double c_obliquity_accuracy_limit(int obliquity_truncation, double tolerance, int max_degree_l)
    double c_obliquity_truncation_limit(int obliquity_truncation, int max_degree_l)
    int c_recommend_obliquity_truncation(double obliquity, double tolerance, int max_degree_l)


cdef extern from "obliquity_common_.hpp" nogil:
    ctypedef pair[c_IntMap[c_Key3, double], c_IntMap[c_Key2, c_IntMap[c_Key1, double]]] ObliquityFuncOutput

    const int C_NUM_OBLIQUITY_TRUNCATIONS
    const int* C_OBLIQUITY_TRUNCATIONS
    cpp_bool c_obliquity_truncation_tabulated(int truncation)

    cdef cppclass c_ObliquitySeriesTable:
        int degree_l
        int truncation
        cpp_bool valid()

    cdef cppclass c_ObliquityGeneralTable:
        int degree_l
        cpp_bool valid()
