# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

cdef class IntMap4:
    cdef void c_reserve(self, size_t n) noexcept nogil:
        self.intmap_cinst.reserve(n)

    cdef void c_clear(self) noexcept nogil:
        self.intmap_cinst.clear()

    cdef void c_set(self, c_Key4& key, double value) noexcept nogil:
        self.intmap_cinst.set(key, value)

    cdef size_t c_size(self) noexcept nogil:
        return self.intmap_cinst.size()
    
    cdef cpp_bool c_get(self, double& result, c_Key4& key) noexcept nogil:
        cdef cpp_bool found = False
        result = self.intmap_cinst.get(found, key)
        return found
    
    ## Python wrappers
    def reserve(self, size_t n):
        self.c_reserve(n)
    
    def clear(self):
        self.c_clear()
    
    def size(self):
        return self.c_size()
    
    def set(self, tuple key, double value):
        if len(key) != 4:
            raise ValueError("Key must be a tuple of 3 integers (l, m, p, q)")

        cdef int16_t a = key[0]
        cdef int16_t b = key[1]
        cdef int16_t c = key[2]
        cdef int16_t d = key[3]
        cdef c_Key4 c_key = c_Key4(a, b, c, d)
        self.c_set(c_key, value)
    
    def get(self, tuple key):
        if len(key) != 4:
            raise ValueError("Key must be a tuple of 3 integers (l, m, p, q)")

        cdef int16_t a = key[0]
        cdef int16_t b = key[1]
        cdef int16_t c = key[2]
        cdef int16_t d = key[3]
        cdef c_Key4 c_key = c_Key4(a, b, c, d)

        cdef double result = 0.0
        cdef cpp_bool found = self.c_get(result, c_key)
        if not found:
            raise KeyError(f"Can not find entry for key: ({key}).")
        return result
    
    def __setitem__(self, tuple key, double value):
        self.set(key, value)

    def __getitem__(self, tuple key):
        return self.get(key)


cdef class IntMap3:
    cdef void c_reserve(self, size_t n) noexcept nogil:
        self.intmap_cinst.reserve(n)

    cdef void c_clear(self) noexcept nogil:
        self.intmap_cinst.clear()

    cdef void c_set(self, c_Key3& key, double value) noexcept nogil:
        self.intmap_cinst.set(key, value)

    cdef size_t c_size(self) noexcept nogil:
        return self.intmap_cinst.size()
    
    cdef cpp_bool c_get(self, double& result, c_Key3& key) noexcept nogil:
        cdef cpp_bool found = False
        result = self.intmap_cinst.get(found, key)
        return found
    
    ## Python wrappers
    def reserve(self, size_t n):
        self.c_reserve(n)
    
    def clear(self):
        self.c_clear()
    
    def size(self):
        return self.c_size()
    
    def set(self, tuple key, double value):
        if len(key) != 3:
            raise ValueError("Key must be a tuple of 3 integers (l, m, p)")

        cdef int16_t a = key[0]
        cdef int16_t b = key[1]
        cdef int16_t c = key[2]
        cdef c_Key3 c_key = c_Key3(a, b, c)
        self.c_set(c_key, value)
    
    def get(self, tuple key):
        if len(key) != 3:
            raise ValueError("Key must be a tuple of 3 integers (l, m, p)")

        cdef int16_t a = key[0]
        cdef int16_t b = key[1]
        cdef int16_t c = key[2]
        cdef c_Key3 c_key = c_Key3(a, b, c)

        cdef double result = 0.0
        cdef cpp_bool found = self.c_get(result, c_key)
        if not found:
            raise KeyError(f"Can not find entry for key: ({key}).")
        return result
    
    def __setitem__(self, tuple key, double value):
        self.set(key, value)

    def __getitem__(self, tuple key):
        return self.get(key)


cdef class IntMap2:
    cdef void c_reserve(self, size_t n) noexcept nogil:
        self.intmap_cinst.reserve(n)

    cdef void c_clear(self) noexcept nogil:
        self.intmap_cinst.clear()

    cdef void c_set(self, c_Key2& key, double value) noexcept nogil:
        self.intmap_cinst.set(key, value)

    cdef size_t c_size(self) noexcept nogil:
        return self.intmap_cinst.size()
    
    cdef cpp_bool c_get(self, double& result, c_Key2& key) noexcept nogil:
        cdef cpp_bool found = False
        result = self.intmap_cinst.get(found, key)
        return found
    
    ## Python wrappers
    def reserve(self, size_t n):
        self.c_reserve(n)
    
    def clear(self):
        self.c_clear()
    
    def size(self):
        return self.c_size()
    
    def set(self, tuple key, double value):
        if len(key) != 2:
            raise ValueError("Key must be a tuple of 2 integers (l, m)")

        cdef int16_t a = key[0]
        cdef int16_t b = key[1]
        cdef c_Key2 c_key = c_Key2(a, b)
        self.c_set(c_key, value)
    
    def get(self, tuple key):
        if len(key) != 2:
            raise ValueError("Key must be a tuple of 2 integers (l, m)")

        cdef int16_t a = key[0]
        cdef int16_t b = key[1]
        cdef c_Key2 c_key = c_Key2(a, b)

        cdef double result = 0.0
        cdef cpp_bool found = self.c_get(result, c_key)
        if not found:
            raise KeyError(f"Can not find entry for key: ({key}).")
        return result
    
    def __setitem__(self, tuple key, double value):
        self.set(key, value)

    def __getitem__(self, tuple key):
        return self.get(key)


cdef class IntMap1:
    cdef void c_reserve(self, size_t n) noexcept nogil:
        self.intmap_cinst.reserve(n)

    cdef void c_clear(self) noexcept nogil:
        self.intmap_cinst.clear()

    cdef void c_set(self, c_Key1& key, double value) noexcept nogil:
        self.intmap_cinst.set(key, value)

    cdef size_t c_size(self) noexcept nogil:
        return self.intmap_cinst.size()
    
    cdef cpp_bool c_get(self, double& result, c_Key1& key) noexcept nogil:
        cdef cpp_bool found = False
        result = self.intmap_cinst.get(found, key)
        return found
    
    ## Python wrappers
    def reserve(self, size_t n):
        self.c_reserve(n)
    
    def clear(self):
        self.c_clear()
    
    def size(self):
        return self.c_size()
    
    def set(self, tuple key, double value):
        if len(key) != 1:
            raise ValueError("Key must be a tuple of 1 integers (l,)")

        cdef int16_t a = key[0]
        cdef c_Key1 c_key = c_Key1(a)
        self.c_set(c_key, value)
    
    def get(self, tuple key):
        if len(key) != 1:
            raise ValueError("Key must be a tuple of 1 integers (l)")

        cdef int16_t a = key[0]
        cdef c_Key1 c_key = c_Key1(a)

        cdef double result = 0.0
        cdef cpp_bool found = self.c_get(result, c_key)
        if not found:
            raise KeyError(f"Can not find entry for key: ({key}).")
        return result
    
    def __setitem__(self, tuple key, double value):
        self.set(key, value)

    def __getitem__(self, tuple key):
        return self.get(key)


cdef class IntMap4Complex:
    cdef void c_reserve(self, size_t n) noexcept nogil:
        self.intmap_cinst.reserve(n)

    cdef void c_clear(self) noexcept nogil:
        self.intmap_cinst.clear()

    cdef void c_set(self, c_Key4& key, double complex value) noexcept nogil:
        self.intmap_cinst.set(key, value)

    cdef size_t c_size(self) noexcept nogil:
        return self.intmap_cinst.size()
    
    cdef cpp_bool c_get(self, double complex& result, c_Key4& key) noexcept nogil:
        cdef cpp_bool found = False
        result = self.intmap_cinst.get(found, key)
        return found
    
    ## Python wrappers
    def reserve(self, size_t n):
        self.c_reserve(n)
    
    def clear(self):
        self.c_clear()
    
    def size(self):
        return self.c_size()
    
    def set(self, tuple key, double complex value):
        if len(key) != 4:
            raise ValueError("Key must be a tuple of 3 integers (l, m, p, q)")

        cdef int16_t a = key[0]
        cdef int16_t b = key[1]
        cdef int16_t c = key[2]
        cdef int16_t d = key[3]
        cdef c_Key4 c_key = c_Key4(a, b, c, d)
        self.c_set(c_key, value)
    
    def get(self, tuple key):
        if len(key) != 4:
            raise ValueError("Key must be a tuple of 3 integers (l, m, p, q)")

        cdef int16_t a = key[0]
        cdef int16_t b = key[1]
        cdef int16_t c = key[2]
        cdef int16_t d = key[3]
        cdef c_Key4 c_key = c_Key4(a, b, c, d)

        cdef double complex result = 0.0
        cdef cpp_bool found = self.c_get(result, c_key)
        if not found:
            raise KeyError(f"Can not find entry for key: ({key}).")
        return result
    
    def __setitem__(self, tuple key, double complex value):
        self.set(key, value)

    def __getitem__(self, tuple key):
        return self.get(key)


cdef class IntMap3Complex:
    cdef void c_reserve(self, size_t n) noexcept nogil:
        self.intmap_cinst.reserve(n)

    cdef void c_clear(self) noexcept nogil:
        self.intmap_cinst.clear()

    cdef void c_set(self, c_Key3& key, double complex value) noexcept nogil:
        self.intmap_cinst.set(key, value)

    cdef size_t c_size(self) noexcept nogil:
        return self.intmap_cinst.size()
    
    cdef cpp_bool c_get(self, double complex& result, c_Key3& key) noexcept nogil:
        cdef cpp_bool found = False
        result = self.intmap_cinst.get(found, key)
        return found
    
    ## Python wrappers
    def reserve(self, size_t n):
        self.c_reserve(n)
    
    def clear(self):
        self.c_clear()
    
    def size(self):
        return self.c_size()
    
    def set(self, tuple key, double complex value):
        if len(key) != 3:
            raise ValueError("Key must be a tuple of 3 integers (l, m, p)")

        cdef int16_t a = key[0]
        cdef int16_t b = key[1]
        cdef int16_t c = key[2]
        cdef c_Key3 c_key = c_Key3(a, b, c)
        self.c_set(c_key, value)
    
    def get(self, tuple key):
        if len(key) != 3:
            raise ValueError("Key must be a tuple of 3 integers (l, m, p)")

        cdef int16_t a = key[0]
        cdef int16_t b = key[1]
        cdef int16_t c = key[2]
        cdef c_Key3 c_key = c_Key3(a, b, c)

        cdef double complex result = 0.0
        cdef cpp_bool found = self.c_get(result, c_key)
        if not found:
            raise KeyError(f"Can not find entry for key: ({key}).")
        return result
    
    def __setitem__(self, tuple key, double complex value):
        self.set(key, value)

    def __getitem__(self, tuple key):
        return self.get(key)


cdef class IntMap2Complex:
    cdef void c_reserve(self, size_t n) noexcept nogil:
        self.intmap_cinst.reserve(n)

    cdef void c_clear(self) noexcept nogil:
        self.intmap_cinst.clear()

    cdef void c_set(self, c_Key2& key, double complex value) noexcept nogil:
        self.intmap_cinst.set(key, value)

    cdef size_t c_size(self) noexcept nogil:
        return self.intmap_cinst.size()
    
    cdef cpp_bool c_get(self, double complex& result, c_Key2& key) noexcept nogil:
        cdef cpp_bool found = False
        result = self.intmap_cinst.get(found, key)
        return found
    
    ## Python wrappers
    def reserve(self, size_t n):
        self.c_reserve(n)
    
    def clear(self):
        self.c_clear()
    
    def size(self):
        return self.c_size()
    
    def set(self, tuple key, double complex value):
        if len(key) != 2:
            raise ValueError("Key must be a tuple of 2 integers (l, m)")

        cdef int16_t a = key[0]
        cdef int16_t b = key[1]
        cdef c_Key2 c_key = c_Key2(a, b)
        self.c_set(c_key, value)
    
    def get(self, tuple key):
        if len(key) != 2:
            raise ValueError("Key must be a tuple of 2 integers (l, m)")

        cdef int16_t a = key[0]
        cdef int16_t b = key[1]
        cdef c_Key2 c_key = c_Key2(a, b)

        cdef double complex result = 0.0
        cdef cpp_bool found = self.c_get(result, c_key)
        if not found:
            raise KeyError(f"Can not find entry for key: ({key}).")
        return result
    
    def __setitem__(self, tuple key, double complex value):
        self.set(key, value)

    def __getitem__(self, tuple key):
        return self.get(key)


cdef class IntMap1Complex:
    cdef void c_reserve(self, size_t n) noexcept nogil:
        self.intmap_cinst.reserve(n)

    cdef void c_clear(self) noexcept nogil:
        self.intmap_cinst.clear()

    cdef void c_set(self, c_Key1& key, double complex value) noexcept nogil:
        self.intmap_cinst.set(key, value)

    cdef size_t c_size(self) noexcept nogil:
        return self.intmap_cinst.size()
    
    cdef cpp_bool c_get(self, double complex& result, c_Key1& key) noexcept nogil:
        cdef cpp_bool found = False
        result = self.intmap_cinst.get(found, key)
        return found
    
    ## Python wrappers
    def reserve(self, size_t n):
        self.c_reserve(n)
    
    def clear(self):
        self.c_clear()
    
    def size(self):
        return self.c_size()
    
    def set(self, tuple key, double complex value):
        if len(key) != 1:
            raise ValueError("Key must be a tuple of 1 integers (l,)")

        cdef int16_t a = key[0]
        cdef c_Key1 c_key = c_Key1(a)
        self.c_set(c_key, value)
    
    def get(self, tuple key):
        if len(key) != 1:
            raise ValueError("Key must be a tuple of 1 integers (l)")

        cdef int16_t a = key[0]
        cdef c_Key1 c_key = c_Key1(a)

        cdef double complex result = 0.0
        cdef cpp_bool found = self.c_get(result, c_key)
        if not found:
            raise KeyError(f"Can not find entry for key: ({key}).")
        return result
    
    def __setitem__(self, tuple key, double complex value):
        self.set(key, value)

    def __getitem__(self, tuple key):
        return self.get(key)