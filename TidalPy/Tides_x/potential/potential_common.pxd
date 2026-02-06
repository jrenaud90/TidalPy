from libc.stdint cimport int16_t

from TidalPy.utilities.lookups cimport c_IntMap, c_Key2, c_Key4


cdef extern from "potential_common_.hpp" nogil:

    struct UniqueFrequencyIndex:
        int16_t l
        int16_t k
    
    ctypedef c_IntMap[c_Key4, double] ModeMap
    ctypedef c_IntMap[c_Key4, UniqueFrequencyIndex] UniqueFreqIndexMap
    ctypedef c_IntMap[c_Key2, double] UniqueFreqMap
