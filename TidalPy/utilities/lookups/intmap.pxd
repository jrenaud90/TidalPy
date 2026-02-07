from libc.stdint cimport uint64_t, int16_t

from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp cimport bool as cpp_bool


cdef extern from "keys_.hpp" nogil:
    ctypedef uint64_t RefKeyType

    cdef RefKeyType convert_4key(int16_t a, int16_t b, int16_t c, int16_t d)
    cdef RefKeyType convert_3key(int16_t a, int16_t b, int16_t c)
    cdef RefKeyType convert_2key(int16_t a, int16_t b)
    cdef RefKeyType convert_1key(int16_t a)

    cdef cppclass c_Key4:
        int16_t a
        int16_t b
        int16_t c
        int16_t d
        RefKeyType reference
        c_Key4() except +
        c_Key4(int16_t a, int16_t b, int16_t c, int16_t d) except +
        void rebuild_reference()
    
    cdef cppclass c_Key3:
        int16_t a
        int16_t b
        int16_t c
        RefKeyType reference
        c_Key3() except +
        c_Key3(int16_t a, int16_t b, int16_t c) except +
        void rebuild_reference()
    
    cdef cppclass c_Key2:
        int16_t a
        int16_t b
        RefKeyType reference
        c_Key2() except +
        c_Key2(int16_t a, int16_t b) except +
        void rebuild_reference()
    
    cdef cppclass c_Key1:
        int16_t a
        RefKeyType reference
        c_Key1() except +
        c_Key1(int16_t a) except +
        void rebuild_reference()

cdef extern from "intmap_.hpp" nogil:

    cdef cppclass c_IntMap[KeyType, ValueType]:
        vector[pair[KeyType, ValueType]] data
        c_IntMap() except +  # For simplicity, only define the regular no-argument constructor.
        void reserve(size_t n)
        void clear()
        size_t size() const
        void set(const KeyType& key, const ValueType& value)
        ValueType get(cpp_bool& o_found, const KeyType& key) const
        const ValueType* get_ptr(bool& o_found, const KeyType& key) const

cdef class IntMap4:
    cdef c_IntMap[c_Key4, double] intmap_cinst

    cdef void c_reserve(self, size_t n) noexcept nogil
    cdef void c_clear(self) noexcept nogil
    cdef void c_set(self, c_Key4& key, double value) noexcept nogil
    cdef size_t c_size(self) noexcept nogil    
    cdef cpp_bool c_get(self, double& result, c_Key4& key) noexcept nogil

cdef class IntMap3:
    cdef c_IntMap[c_Key3, double] intmap_cinst

    cdef void c_reserve(self, size_t n) noexcept nogil
    cdef void c_clear(self) noexcept nogil
    cdef void c_set(self, c_Key3& key, double value) noexcept nogil
    cdef size_t c_size(self) noexcept nogil    
    cdef cpp_bool c_get(self, double& result, c_Key3& key) noexcept nogil

cdef class IntMap2:
    cdef c_IntMap[c_Key2, double] intmap_cinst

    cdef void c_reserve(self, size_t n) noexcept nogil
    cdef void c_clear(self) noexcept nogil
    cdef void c_set(self, c_Key2& key, double value) noexcept nogil
    cdef size_t c_size(self) noexcept nogil    
    cdef cpp_bool c_get(self, double& result, c_Key2& key) noexcept nogil

cdef class IntMap1:
    cdef c_IntMap[c_Key1, double] intmap_cinst

    cdef void c_reserve(self, size_t n) noexcept nogil
    cdef void c_clear(self) noexcept nogil
    cdef void c_set(self, c_Key1& key, double value) noexcept nogil
    cdef size_t c_size(self) noexcept nogil    
    cdef cpp_bool c_get(self, double& result, c_Key1& key) noexcept nogil

cdef class IntMap4Complex:
    cdef c_IntMap[c_Key4, double complex] intmap_cinst

    cdef void c_reserve(self, size_t n) noexcept nogil
    cdef void c_clear(self) noexcept nogil
    cdef void c_set(self, c_Key4& key, double complex value) noexcept nogil
    cdef size_t c_size(self) noexcept nogil    
    cdef cpp_bool c_get(self, double complex& result, c_Key4& key) noexcept nogil

cdef class IntMap3Complex:
    cdef c_IntMap[c_Key3, double complex] intmap_cinst

    cdef void c_reserve(self, size_t n) noexcept nogil
    cdef void c_clear(self) noexcept nogil
    cdef void c_set(self, c_Key3& key, double complex value) noexcept nogil
    cdef size_t c_size(self) noexcept nogil    
    cdef cpp_bool c_get(self, double complex& result, c_Key3& key) noexcept nogil

cdef class IntMap2Complex:
    cdef c_IntMap[c_Key2, double complex] intmap_cinst

    cdef void c_reserve(self, size_t n) noexcept nogil
    cdef void c_clear(self) noexcept nogil
    cdef void c_set(self, c_Key2& key, double complex value) noexcept nogil
    cdef size_t c_size(self) noexcept nogil    
    cdef cpp_bool c_get(self, double complex& result, c_Key2& key) noexcept nogil

cdef class IntMap1Complex:
    cdef c_IntMap[c_Key1, double complex] intmap_cinst

    cdef void c_reserve(self, size_t n) noexcept nogil
    cdef void c_clear(self) noexcept nogil
    cdef void c_set(self, c_Key1& key, double complex value) noexcept nogil
    cdef size_t c_size(self) noexcept nogil    
    cdef cpp_bool c_get(self, double complex& result, c_Key1& key) noexcept nogil