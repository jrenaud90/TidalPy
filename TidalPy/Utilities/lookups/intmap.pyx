# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Python wrappers of the C++ integer-keyed maps (``intmap_.hpp``): one class per key length (1 to 4 small integers)
and value type (``double`` or ``double complex``).

Cython has no class templates, so each class holds its own typed ``c_IntMap`` and the typed ``c_*`` calls other
extensions use; everything Python sees (``set``, ``get``, indexing, ``len``, iteration) lives once in
``IntMapBase``, which reaches the map through each class's ``_set_item``, ``_get_item``, and ``_entry``.
"""


cdef class IntMapBase:
    """Shared Python interface of the IntMap classes. Not for direct use; instantiate ``IntMap1`` to ``IntMap4`` or
    their ``Complex`` variants."""

    def __init__(self):
        if type(self) is IntMapBase:
            raise TypeError("IntMapBase is abstract; instantiate IntMap1 to IntMap4 or their Complex variants.")

    cdef size_t c_size(self) noexcept nogil:
        return 0

    cdef void c_reserve(self, size_t n) noexcept nogil:
        pass

    cdef void c_clear(self) noexcept nogil:
        pass

    cdef void _set_item(self, tuple key, object value) except *:
        raise NotImplementedError

    cdef object _get_item(self, tuple key, cpp_bool* found):
        raise NotImplementedError

    cdef tuple _entry(self, size_t i):
        raise NotImplementedError

    cdef void _check_key(self, tuple key) except *:
        if len(key) != self.key_length:
            raise ValueError(f"Key must be a tuple of {self.key_length} integers.")

    def reserve(self, size_t n):
        self.c_reserve(n)

    def clear(self):
        self.c_clear()

    def size(self):
        return self.c_size()

    def set(self, tuple key, value):
        self._check_key(key)
        self._set_item(key, value)

    def get(self, tuple key):
        self._check_key(key)
        cdef cpp_bool found = False
        value = self._get_item(key, &found)
        if not found:
            raise KeyError(f"Can not find entry for key: ({key}).")
        return value

    def __setitem__(self, tuple key, value):
        self.set(key, value)

    def __getitem__(self, tuple key):
        return self.get(key)

    def __len__(self):
        return self.c_size()

    def __iter__(self):
        """Yields pairs of (key tuple, value)."""
        cdef size_t i
        for i in range(self.c_size()):
            yield self._entry(i)


cdef class IntMap4(IntMapBase):
    """4-integer keys to float values."""

    def __cinit__(self):
        self.key_length = 4

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

    cdef void _set_item(self, tuple key, object value) except *:
        cdef c_Key4 c_key = c_Key4(key[0], key[1], key[2], key[3])
        self.c_set(c_key, value)

    cdef object _get_item(self, tuple key, cpp_bool* found):
        cdef c_Key4 c_key = c_Key4(key[0], key[1], key[2], key[3])
        cdef double result = 0.0
        found[0] = self.c_get(result, c_key)
        return result

    cdef tuple _entry(self, size_t i):
        cdef c_Key4 key = self.intmap_cinst.data[i].first
        return ((key.a, key.b, key.c, key.d), self.intmap_cinst.data[i].second)


cdef class IntMap3(IntMapBase):
    """3-integer keys to float values."""

    def __cinit__(self):
        self.key_length = 3

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

    cdef void _set_item(self, tuple key, object value) except *:
        cdef c_Key3 c_key = c_Key3(key[0], key[1], key[2])
        self.c_set(c_key, value)

    cdef object _get_item(self, tuple key, cpp_bool* found):
        cdef c_Key3 c_key = c_Key3(key[0], key[1], key[2])
        cdef double result = 0.0
        found[0] = self.c_get(result, c_key)
        return result

    cdef tuple _entry(self, size_t i):
        cdef c_Key3 key = self.intmap_cinst.data[i].first
        return ((key.a, key.b, key.c), self.intmap_cinst.data[i].second)


cdef class IntMap2(IntMapBase):
    """2-integer keys to float values."""

    def __cinit__(self):
        self.key_length = 2

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

    cdef void _set_item(self, tuple key, object value) except *:
        cdef c_Key2 c_key = c_Key2(key[0], key[1])
        self.c_set(c_key, value)

    cdef object _get_item(self, tuple key, cpp_bool* found):
        cdef c_Key2 c_key = c_Key2(key[0], key[1])
        cdef double result = 0.0
        found[0] = self.c_get(result, c_key)
        return result

    cdef tuple _entry(self, size_t i):
        cdef c_Key2 key = self.intmap_cinst.data[i].first
        return ((key.a, key.b), self.intmap_cinst.data[i].second)


cdef class IntMap1(IntMapBase):
    """1-integer keys to float values."""

    def __cinit__(self):
        self.key_length = 1

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

    cdef void _set_item(self, tuple key, object value) except *:
        cdef c_Key1 c_key = c_Key1(key[0])
        self.c_set(c_key, value)

    cdef object _get_item(self, tuple key, cpp_bool* found):
        cdef c_Key1 c_key = c_Key1(key[0])
        cdef double result = 0.0
        found[0] = self.c_get(result, c_key)
        return result

    cdef tuple _entry(self, size_t i):
        cdef c_Key1 key = self.intmap_cinst.data[i].first
        return ((key.a,), self.intmap_cinst.data[i].second)


cdef class IntMap4Complex(IntMapBase):
    """4-integer keys to complex values."""

    def __cinit__(self):
        self.key_length = 4

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

    cdef void _set_item(self, tuple key, object value) except *:
        cdef c_Key4 c_key = c_Key4(key[0], key[1], key[2], key[3])
        self.c_set(c_key, value)

    cdef object _get_item(self, tuple key, cpp_bool* found):
        cdef c_Key4 c_key = c_Key4(key[0], key[1], key[2], key[3])
        cdef double complex result = 0.0
        found[0] = self.c_get(result, c_key)
        return result

    cdef tuple _entry(self, size_t i):
        cdef c_Key4 key = self.intmap_cinst.data[i].first
        return ((key.a, key.b, key.c, key.d), self.intmap_cinst.data[i].second)


cdef class IntMap3Complex(IntMapBase):
    """3-integer keys to complex values."""

    def __cinit__(self):
        self.key_length = 3

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

    cdef void _set_item(self, tuple key, object value) except *:
        cdef c_Key3 c_key = c_Key3(key[0], key[1], key[2])
        self.c_set(c_key, value)

    cdef object _get_item(self, tuple key, cpp_bool* found):
        cdef c_Key3 c_key = c_Key3(key[0], key[1], key[2])
        cdef double complex result = 0.0
        found[0] = self.c_get(result, c_key)
        return result

    cdef tuple _entry(self, size_t i):
        cdef c_Key3 key = self.intmap_cinst.data[i].first
        return ((key.a, key.b, key.c), self.intmap_cinst.data[i].second)


cdef class IntMap2Complex(IntMapBase):
    """2-integer keys to complex values."""

    def __cinit__(self):
        self.key_length = 2

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

    cdef void _set_item(self, tuple key, object value) except *:
        cdef c_Key2 c_key = c_Key2(key[0], key[1])
        self.c_set(c_key, value)

    cdef object _get_item(self, tuple key, cpp_bool* found):
        cdef c_Key2 c_key = c_Key2(key[0], key[1])
        cdef double complex result = 0.0
        found[0] = self.c_get(result, c_key)
        return result

    cdef tuple _entry(self, size_t i):
        cdef c_Key2 key = self.intmap_cinst.data[i].first
        return ((key.a, key.b), self.intmap_cinst.data[i].second)


cdef class IntMap1Complex(IntMapBase):
    """1-integer keys to complex values."""

    def __cinit__(self):
        self.key_length = 1

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

    cdef void _set_item(self, tuple key, object value) except *:
        cdef c_Key1 c_key = c_Key1(key[0])
        self.c_set(c_key, value)

    cdef object _get_item(self, tuple key, cpp_bool* found):
        cdef c_Key1 c_key = c_Key1(key[0])
        cdef double complex result = 0.0
        found[0] = self.c_get(result, c_key)
        return result

    cdef tuple _entry(self, size_t i):
        cdef c_Key1 key = self.intmap_cinst.data[i].first
        return ((key.a,), self.intmap_cinst.data[i].second)
