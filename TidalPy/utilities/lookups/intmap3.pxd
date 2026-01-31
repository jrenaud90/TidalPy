from libcpp cimport bool as cpp_bool
from libc.stdint cimport uint32_t


cdef extern from "intmap3_.cpp" nogil:
    ctypedef uint32_t KeyType
    
    cdef KeyType convet_3key(size_t l, size_t m, size_t p)

    cdef cppclass c_IntMap3:
        c_IntMap3() except +
        void reserve(size_t n)
        void clear()
        void set(size_t l, size_t m, size_t p, double value)
        double get(cpp_bool& o_found, size_t l, size_t m, size_t p) const
        size_t size() const


cdef class IntMap3:
    cdef c_IntMap3 intmap3_cinst

    cdef void c_reserve(self, size_t n) noexcept nogil
    cdef void c_clear(self) noexcept nogil
    cdef void c_set(self, size_t l, size_t m, size_t p, double value) noexcept nogil
    cdef size_t c_size(self) noexcept nogil    
    cdef cpp_bool c_get(self, double& result, size_t l, size_t m, size_t p) noexcept nogil
