# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

def convert_from_mode_storage(c_ModeStorage mode_storage_inst):
    """Converts C++ struct `c_ModeStorage` to Python type."""
    cdef double mode = mode_storage_inst.mode
    cdef double mode_strength = mode_storage_inst.mode_strength
    cdef int n_coeff = mode_storage_inst.n_coeff
    cdef int o_coeff = mode_storage_inst.o_coeff
    cdef tuple result = (mode, mode_strength, n_coeff, o_coeff)
    return result

cdef c_ModeStorage c_convert_to_mode_storage(tuple mode_storage_tuple):
    """Converts Python tuple to C++ struct `c_ModeStorage`."""
    if len(mode_storage_tuple) != 4:
        raise ValueError("`mode_storage_tuple` must be a tuple of 4 values: (double, double, int, int).")
    if not isinstance(mode_storage_tuple[0], float):
        raise TypeError("`mode_storage_tuple[0] must be a floating point value (the 'mode').")
    if not isinstance(mode_storage_tuple[1], float):
        raise TypeError("`mode_storage_tuple[1] must be a floating point value (the 'mode strength').")
    if not isinstance(mode_storage_tuple[2], int):
        raise TypeError("`mode_storage_tuple[2] must be a integer value (the 'n coefficient').")
    if not isinstance(mode_storage_tuple[3], int):
        raise TypeError("`mode_storage_tuple[3] must be a integer value (the 'O coefficient').")
    cdef double mode = mode_storage_tuple[0]
    cdef double mode_strength = mode_storage_tuple[1]
    cdef int n_coeff = mode_storage_tuple[2]
    cdef int o_coeff = mode_storage_tuple[3]
    return c_ModeStorage(mode, mode_strength, n_coeff, o_coeff)


cdef class ModeMap:

    def __cinit__(self, c_ModeMap mode_map_cinst_):

        # Store the C++ mode map instance
        self.mode_map_cinst = mode_map_cinst_
    
    cdef void c_reserve(self, size_t n) noexcept nogil:
        self.mode_map_cinst.reserve(n)

    cdef void c_clear(self) noexcept nogil:
        self.mode_map_cinst.clear()

    cdef void c_set(self, c_Key4& key, c_ModeStorage& value) noexcept nogil:
        self.mode_map_cinst.set(key, value)

    cdef size_t c_size(self) noexcept nogil:
        return self.mode_map_cinst.size()
    
    cdef cpp_bool c_get(self, c_ModeStorage& result, c_Key4& key) noexcept nogil:
        cdef cpp_bool found = False
        result = self.mode_map_cinst.get(found, key)
        return found
    
    ## Python wrappers
    def reserve(self, size_t n):
        self.c_reserve(n)
    
    def clear(self):
        self.c_clear()
    
    def size(self):
        return self.c_size()
    
    def set(self, tuple key, tuple mode_storage_tuple):
        if len(key) != 4:
            raise ValueError("Key must be a tuple of 4 integers (l, m, p, q)")
        if len(mode_storage_tuple) != 4:
            raise ValueError("`mode_storage_tuple` must be a tuple of 4 values: (double, double, int, int).")
        if not isinstance(mode_storage_tuple[0], float):
            raise TypeError("`mode_storage_tuple[0] must be a floating point value (the 'mode').")
        if not isinstance(mode_storage_tuple[1], float):
            raise TypeError("`mode_storage_tuple[1] must be a floating point value (the 'mode strength').")
        if not isinstance(mode_storage_tuple[2], int):
            raise TypeError("`mode_storage_tuple[2] must be a integer value (the 'n coefficient').")
        if not isinstance(mode_storage_tuple[3], int):
            raise TypeError("`mode_storage_tuple[3] must be a integer value (the 'O coefficient').")

        # Build Key
        cdef int16_t l = key[0]
        cdef int16_t m = key[1]
        cdef int16_t p = key[2]
        cdef int16_t q = key[4]
        cdef c_Key4 c_key = c_Key4(l, m, p, q)

        # Build value
        cdef c_ModeStorage mode_storage = c_convert_to_mode_storage(mode_storage_tuple)

        self.c_set(c_key, mode_storage)
    
    def get(self, tuple key):
        if len(key) != 4:
            raise ValueError("Key must be a tuple of 4 integers (l, m, p, q)")

        cdef int16_t l = key[0]
        cdef int16_t m = key[1]
        cdef int16_t p = key[2]
        cdef int16_t q = key[4]
        cdef c_Key4 c_key = c_Key4(l, m, p, q)

        cdef c_ModeStorage result_cinst
        cdef cpp_bool found = self.c_get(result_cinst, c_key)
        if not found:
            raise KeyError(f"Can not find entry for key: ({key}).")
        
        # Convert result to python readable.
        return convert_from_mode_storage(result_cinst)

    def __setitem__(self, tuple key, tuple value):
        self.set(key, value)

    def __getitem__(self, tuple key):
        return self.get(key)
    
    def __len__(self):
        return self.mode_map_cinst.size()

    def __iter__(self):
        """
        Yields pairs of ((a, b, c), value).
        """

        cdef size_t i
        cdef c_Key4 key
        cdef c_ModeStorage value

        for i in range(self.mode_map_cinst.size()):
            key = self.mode_map_cinst.data[i].first
            value = self.mode_map_cinst.data[i].second

            yield ((key.a, key.b, key.c, key.d), convert_from_mode_storage(value))


def test_mode_map():

    cdef c_ModeMap mode_map_cinst

    # Fill with test data
    cdef size_t l, m, p, q
    cdef size_t total_size = 0
    cdef c_Key4 key
    cdef c_ModeStorage storage
    for l in range(2, 4):
        key.a = l
        for m in range(0, l + 1):
            key.b = m
            for p in range(0, l + 1):
                key.c = p
                for q in range(-1, 2):
                    key.d = q
                    key.rebuild_reference()

                    storage.mode = (l - 2 * p + q) * 10.0 - m * 5.0
                    storage.mode_strength = l + p + q + m - 4.5
                    storage.n_coeff = (l - 2 * p + q)
                    storage.o_coeff = -m

                    mode_map_cinst.set(key, storage)
                    total_size += 1

    cdef ModeMap mode_map = ModeMap(mode_map_cinst)

    return mode_map, total_size