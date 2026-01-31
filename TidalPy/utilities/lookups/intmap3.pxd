from libcpp cimport bool as cpp_bool
from libc.stdint cimport uint64_t, int16_t


cdef extern from "intmap3_.hpp" nogil:
    ctypedef uint64_t KeyType3
    
    cdef KeyType3 convet_3key(int16_t l, int16_t m, int16_t p)

    cdef cppclass c_IntMap3[T]:
        c_IntMap3() except +
        void reserve(size_t n)
        void clear()
        void set(int16_t l, int16_t m, int16_t p, T value)
        T get(cpp_bool& o_found, int16_t l, int16_t m, int16_t p) const
        size_t size() const


cdef class IntMap3:
    cdef c_IntMap3[double] intmap3_cinst

    cdef void c_reserve(self, size_t n) noexcept nogil
    cdef void c_clear(self) noexcept nogil
    cdef void c_set(self, int16_t l, int16_t m, int16_t p, double value) noexcept nogil
    cdef size_t c_size(self) noexcept nogil    
    cdef cpp_bool c_get(self, double& result, int16_t l, int16_t m, int16_t p) noexcept nogil

cdef class IntMap3Complex:
    cdef c_IntMap3[double complex] intmap3_cinst

    cdef void c_reserve(self, size_t n) noexcept nogil
    cdef void c_clear(self) noexcept nogil
    cdef void c_set(self, int16_t l, int16_t m, int16_t p, double complex value) noexcept nogil
    cdef size_t c_size(self) noexcept nogil    
    cdef cpp_bool c_get(self, double complex& result, int16_t l, int16_t m, int16_t p) noexcept nogil
