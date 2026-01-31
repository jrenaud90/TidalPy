# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

cdef class IntMap4:
    cdef void c_reserve(self, size_t n) noexcept nogil:
        self.intmap4_cinst.reserve(n)

    cdef void c_clear(self) noexcept nogil:
        self.intmap4_cinst.clear()

    cdef void c_set(self, int16_t l, int16_t m, int16_t p, int16_t q, double value) noexcept nogil:
        self.intmap4_cinst.set(l, m, p, q, value)

    cdef size_t c_size(self) noexcept nogil:
        return self.intmap4_cinst.size()
    
    cdef cpp_bool c_get(self, double& result, int16_t l, int16_t m, int16_t p, int16_t q) noexcept nogil:
        cdef cpp_bool found = False
        result = self.intmap4_cinst.get(found, l, m, p, q)
        return found

    ## Python wrappers
    def reserve(self, size_t n):
        self.c_reserve(n)
    
    def clear(self):
        self.c_clear()
    
    def set(self, int16_t l, int16_t m, int16_t p, int16_t q, double value):
        self.c_set(l, m, p, q, value)

    def size(self):
        return self.c_size()
    
    def get(self, int16_t l, int16_t m, int16_t p, int16_t q):
        cdef double result = 0.0
        cdef cpp_bool found = self.c_get(result, l, m, p, q)
        if not found:
            raise KeyError(f"Can not find entry for key: ({l}, {m}, {p}, {q}).")
        return result
    
    def __setitem__(self, tuple key, double value):
        if len(key) != 4:
            raise ValueError("Key must be a tuple of 4 integers (l, m, p, q)")
        # Unpack tuple
        cdef int16_t l = key[0]
        cdef int16_t m = key[1]
        cdef int16_t p = key[2]
        cdef int16_t q = key[3]
        self.c_set(l, m, p, q, value)

    def __getitem__(self, tuple key):
        if len(key) != 4:
            raise ValueError("Key must be a tuple of 4 integers (l, m, p, q)")
        cdef int16_t l = key[0]
        cdef int16_t m = key[1]
        cdef int16_t p = key[2]
        cdef int16_t q = key[3]
        cdef double result = 0.0
        cdef cpp_bool found = self.c_get(result, l, m, p, q)
        if not found:
            raise KeyError(f"Can not find result for ({l}, {m}, {p}, {q}).")
        return result


cdef class IntMap4Complex:
    cdef void c_reserve(self, size_t n) noexcept nogil:
        self.intmap4_cinst.reserve(n)

    cdef void c_clear(self) noexcept nogil:
        self.intmap4_cinst.clear()

    cdef void c_set(self, int16_t l, int16_t m, int16_t p, int16_t q, double complex value) noexcept nogil:
        self.intmap4_cinst.set(l, m, p, q, value)

    cdef size_t c_size(self) noexcept nogil:
        return self.intmap4_cinst.size()
    
    cdef cpp_bool c_get(self, double complex& result, int16_t l, int16_t m, int16_t p, int16_t q) noexcept nogil:
        cdef cpp_bool found = False
        result = self.intmap4_cinst.get(found, l, m, p, q)
        return found

    ## Python wrappers
    def reserve(self, size_t n):
        self.c_reserve(n)
    
    def clear(self):
        self.c_clear()
    
    def set(self, int16_t l, int16_t m, int16_t p, int16_t q, double complex value):
        self.c_set(l, m, p, q, value)

    def size(self):
        return self.c_size()
    
    def get(self, int16_t l, int16_t m, int16_t p, int16_t q):
        cdef double complex result = 0.0
        cdef cpp_bool found = self.c_get(result, l, m, p, q)
        if not found:
            raise KeyError(f"Can not find entry for key: ({l}, {m}, {p}, {q}).")
        return result
    
    def __setitem__(self, tuple key, double complex value):
        if len(key) != 4:
            raise ValueError("Key must be a tuple of 4 integers (l, m, p, q)")
        # Unpack tuple
        cdef int16_t l = key[0]
        cdef int16_t m = key[1]
        cdef int16_t p = key[2]
        cdef int16_t q = key[3]
        self.c_set(l, m, p, q, value)

    def __getitem__(self, tuple key):
        if len(key) != 4:
            raise ValueError("Key must be a tuple of 4 integers (l, m, p, q)")
        cdef int16_t l = key[0]
        cdef int16_t m = key[1]
        cdef int16_t p = key[2]
        cdef int16_t q = key[3]
        cdef double complex result = 0.0
        cdef cpp_bool found = self.c_get(result, l, m, p, q)
        if not found:
            raise KeyError(f"Can not find result for ({l}, {m}, {p}, {q}).")
        return result
