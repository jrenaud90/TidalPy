# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

# Maximum size for array building
cdef size_t MAX_NUM_Y      = 6
cdef size_t MAX_NUM_Y_REAL = 2 * MAX_NUM_Y
cdef size_t MAX_NUM_SOL    = 3
