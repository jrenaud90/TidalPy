# distutils: language = c++

# Functions defined in `TidalPy.utilities.performance.special.interp`
cdef int binary_search_with_guess(double key, double[:] array, int length, int guess) nogil

cpdef double interp(double desired_x, double[:] x_domain, double[:] dependent_values) nogil

cpdef double complex interp_complex(double desired_x, double[:] x_domain, double complex[:] dependent_variables) nogil
