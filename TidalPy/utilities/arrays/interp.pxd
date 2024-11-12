cdef extern from "interp_.cpp" nogil:
    size_t cf_binary_search_with_guess(double key, double* array, size_t length, size_t guess)
