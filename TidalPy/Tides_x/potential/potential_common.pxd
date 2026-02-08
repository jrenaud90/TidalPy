from libc.stdint cimport int16_t

from libcpp.vector cimport vector

from TidalPy.utilities.lookups cimport c_IntMap, c_Key2, c_Key4


cdef extern from "potential_common_.hpp" nogil:

    struct c_FrequencyStorage:
        size_t num_instances
        double frequency
        c_FrequencyStorage(double freq)
    
    struct c_ModeStorage:
        double mode
        double mode_strength
        int n_coeff
        int o_coeff
    
    ctypedef c_IntMap[c_Key4, c_ModeStorage] c_ModeMap
    ctypedef c_IntMap[c_Key4, size_t] c_UniqueFreqIndexMap
    ctypedef vector[c_FrequencyStorage] c_UniqueFreqMap

    inline c_IntMap[c_Key2, double]& c_get_lm_coeff_map()
