from libcpp.pair cimport pair

from TidalPy.Utilities_x.lookups cimport c_IntMap, c_Key1, c_Key2, c_Key3


cdef extern from "eccentricity_common_.hpp" nogil:
    ctypedef pair[c_IntMap[c_Key3, double], c_IntMap[c_Key2, c_IntMap[c_Key1, double]]] EccentricityFuncOutput
