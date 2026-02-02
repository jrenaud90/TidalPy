from libcpp.pair cimport pair

from TidalPy.utilities.lookups cimport c_IntMap, c_Key1, c_Key2, c_Key3


cdef extern from "obliquity_common_.hpp" nogil:
    ctypedef pair[c_IntMap[c_Key3, double], c_IntMap[c_Key2, c_IntMap[c_Key1, double]]] ObliquityFuncOutput
