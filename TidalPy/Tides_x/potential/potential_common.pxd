from libc.stdint cimport int16_t

from libcpp cimport bool as cpp_bool
from libcpp.vector cimport vector

from TidalPy.utilities.lookups cimport c_IntMap, c_Key2, c_Key4


cdef extern from "potential_common_.hpp" nogil:

    cdef cppclass c_FrequencyStorage:
        size_t num_instances
        double frequency
        c_FrequencyStorage() except +
        c_FrequencyStorage(double freq) except +
    
    cdef cppclass c_ModeStorage:
        double mode
        double mode_strength
        int n_coeff
        int o_coeff
        c_ModeStorage() except +
        c_ModeStorage(double mode, double mode_strength, int n_coeff, int o_coeff) except +
        c_ModeStorage(int n_coeff, int o_coeff) except +
    
    ctypedef c_IntMap[c_Key4, c_ModeStorage] c_ModeMap
    ctypedef c_IntMap[c_Key4, size_t] c_UniqueFreqIndexMap
    ctypedef vector[c_FrequencyStorage] c_UniqueFreqMap

    inline c_IntMap[c_Key2, double]& c_get_lm_coeff_map()

cdef class ModeMap:

    cdef c_ModeMap mode_map_cinst

    cdef void c_reserve(self, size_t n) noexcept nogil
    cdef void c_clear(self) noexcept nogil
    cdef void c_set(self, c_Key4& key, c_ModeStorage& value) noexcept nogil
    cdef size_t c_size(self) noexcept nogil
    cdef cpp_bool c_get(self, c_ModeStorage& result, c_Key4& key) noexcept nogil


cdef tuple c_convert_from_mode_storage(c_ModeStorage mode_storage_inst)
cdef c_ModeStorage c_convert_to_mode_storage(tuple mode_storage_tuple)