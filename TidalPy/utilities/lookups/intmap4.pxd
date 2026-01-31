from libcpp cimport bool as cpp_bool
from libc.stdint cimport uint64_t, int16_t


cdef extern from "intmap4_.hpp" nogil:
    ctypedef uint64_t KeyType4
    
    cdef KeyType4 convet_3key(int16_t l, int16_t m, int16_t p)

    cdef cppclass c_IntMap4[T]:
        c_IntMap3() except +
        void reserve(size_t n)
        void clear()
        void set(int16_t l, int16_t m, int16_t p, int16_t q, T value)
        T get(cpp_bool& o_found, int16_t l, int16_t m, int16_t p, int16_t q) const
        size_t size() const


cdef class IntMap4:
    cdef c_IntMap4[double] intmap4_cinst

    cdef void c_reserve(self, size_t n) noexcept nogil
    cdef void c_clear(self) noexcept nogil
    cdef void c_set(self, int16_t l, int16_t m, int16_t p, int16_t q, double value) noexcept nogil
    cdef size_t c_size(self) noexcept nogil    
    cdef cpp_bool c_get(self, double& result, int16_t l, int16_t m, int16_t p, int16_t q) noexcept nogil

cdef class IntMap4Complex:
    cdef c_IntMap4[double complex] intmap4_cinst

    cdef void c_reserve(self, size_t n) noexcept nogil
    cdef void c_clear(self) noexcept nogil
    cdef void c_set(self, int16_t l, int16_t m, int16_t p, int16_t q, double complex value) noexcept nogil
    cdef size_t c_size(self) noexcept nogil    
    cdef cpp_bool c_get(self, double complex& result, int16_t l, int16_t m, int16_t p, int16_t q) noexcept nogil
