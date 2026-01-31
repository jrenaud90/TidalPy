# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

cdef class IntMap3:
    cdef void c_reserve(self, size_t n) noexcept nogil:
        self.intmap3_cinst.reserve(n)

    cdef void c_clear(self) noexcept nogil:
        self.intmap3_cinst.clear()

    cdef void c_set(self, int16_t l, int16_t m, int16_t p, double value) noexcept nogil:
        self.intmap3_cinst.set(l, m, p, value)

    cdef size_t c_size(self) noexcept nogil:
        return self.intmap3_cinst.size()
    
    cdef cpp_bool c_get(self, double& result, int16_t l, int16_t m, int16_t p) noexcept nogil:
        cdef cpp_bool found = False
        result = self.intmap3_cinst.get(found, l, m, p)
        return found
    
    ## Python wrappers
    def reserve(self, size_t n):
        self.c_reserve(n)
    
    def clear(self):
        self.c_clear()
    
    def set(self, int16_t l, int16_t m, int16_t p, double value):
        self.c_set(l, m, p, value)

    def size(self):
        return self.c_size()
    
    def get(self, int16_t l, int16_t m, int16_t p):
        cdef double result = 0.0
        cdef cpp_bool found = self.c_get(result, l, m, p)
        if not found:
            raise KeyError(f"Can not find entry for key: ({l}, {m}, {p}).")
        return result
    
    def __setitem__(self, tuple key, double value):
        if len(key) != 3:
            raise ValueError("Key must be a tuple of 3 integers (l, m, p)")
        # Unpack tuple
        cdef int16_t l = key[0]
        cdef int16_t m = key[1]
        cdef int16_t p = key[2]
        self.c_set(l, m, p, value)

    def __getitem__(self, tuple key):
        if len(key) != 3:
            raise ValueError("Key must be a tuple of 3 integers (l, m, p)")
        cdef int16_t l = key[0]
        cdef int16_t m = key[1]
        cdef int16_t p = key[2]
        cdef double result = 0.0
        cdef cpp_bool found = self.c_get(result, l, m, p)
        if not found:
            raise KeyError(f"Can not find result for ({l}, {m}, {p}).")
        return result
