# distutils: language = c++
import cython
import numpy as np
cimport numpy as np
np.import_array()
from libc.math cimport isnan

# Get machine precision.
cdef double EPS
EPS = np.finfo(dtype=np.float64).eps

# Determine cache limits.
cdef int LIKELY_IN_CACHE_SIZE = 8

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cdef int binary_search_with_guess(double key, double[:] array, int length, int guess) nogil:
    """ Binary search with guess.
    
    Based on `numpy`'s `binary_search_with_guess` function.
    
    Parameters
    ----------
    key : float
        Key index to search for.
    array : np.ndarray
        Array to search in.
    length : int
        Length of array.
    guess : int 
        Initial guess of where key might be.
    
    Returns
    -------
    guess : int
        Corrected guess after search.
    """

    cdef int imin = 0
    cdef int imax = length

    if key > array[length - 1]:
        return length
    elif key < array[0]:
        return -1

    # If len <= 4 use linear search.
#     if length <= 4:
#         raise NotImplemented

    if guess > (length - 3):
        guess = length - 3
    if guess < 1:
        guess = 1

    # check most likely values: guess - 1, guess, guess + 1
    if key < array[guess]:
        if key < array[guess - 1]:
            imax = guess - 1
            # last attempt to restrict search to items in cache
            if guess > LIKELY_IN_CACHE_SIZE and key >= array[guess - LIKELY_IN_CACHE_SIZE]:
                imin = guess - LIKELY_IN_CACHE_SIZE

        else:
            return guess - 1
    else:
        if key < array[guess + 1]:
            return guess
        else:
            if key < array[guess + 2]:
                return guess + 1
            else:
                imin = guess + 2
                # last attempt to restrict search to items in cache
                if guess < (length - LIKELY_IN_CACHE_SIZE - 1) and key < array[guess + LIKELY_IN_CACHE_SIZE]:
                    imax = guess + LIKELY_IN_CACHE_SIZE

    # Finally, find index by bisection
    cdef int imid
    while imin < imax:
        imid = imin + ((imax - imin) >> 1)
        if key >= array[imid]:
            imin = imid + 1
        else:
            imax = imid

    return imin - 1

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef double interp(double desired_x, double[:] x_domain, double[:] dependent_values) nogil:
    """ Interpolation function for floats.
    
    Provided a domain, `x_domain` and a dependent array `dependent_values` search domain for value closest to 
    `desired_x` and return the value of `dependent_values` at that location if it is defined. Otherwise, use local 
    slopes of `x_domain` and `dependent_values` to interpolate a value of `dependent_values` at `desired_x`.

    Based on `numpy`'s `interp` function.

    Parameters
    ----------
    desired_x : float
        Location where `dependent_variables` is desired.
    x_domain : np.ndarray[float]
        Domain to search for the correct location.
    dependent_values : np.ndarray[float]
        Dependent values that are to be returned after search and interpolation.

    Returns
    -------
    result : float
        Desired value of `dependent_values`.
    
    """

    cdef int lenx
    lenx = len(x_domain)
    # TODO: Needs to be at least 3 item long array. Add exception here?

    cdef double left_value
    left_value = dependent_values[0]
    cdef double right_value
    right_value = dependent_values[lenx - 1]

    # Binary Search with Guess
    cdef int i, j
    j = 0
    cdef double slope

    cdef double result
    cdef double fp_at_j
    cdef double xp_at_j
    cdef double fp_at_jp1
    cdef double xp_at_jp1

    # Perform binary search with guess
    j = binary_search_with_guess(desired_x, x_domain, lenx, j)

    if j == -1:
        result = left_value
    elif j == lenx:
        result = right_value
    else:
        fp_at_j = dependent_values[j]
        xp_at_j = x_domain[j]
        if j == lenx - 1:
            result = fp_at_j
        elif xp_at_j == desired_x:
            result = fp_at_j
        else:
            fp_at_jp1 = dependent_values[j + 1]
            xp_at_jp1 = x_domain[j + 1]
            slope = (fp_at_jp1 - fp_at_j) / (xp_at_jp1 - xp_at_j)

            # If we get nan in one direction, try the other
            result = slope * (desired_x - xp_at_j) + fp_at_j
            if isnan(result):
                result = slope * (desired_x - xp_at_jp1) + fp_at_jp1
                if isnan(result) and (fp_at_jp1 == fp_at_j):
                    result = fp_at_j

    return result


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef double complex interp_complex(double desired_x, double[:] x_domain,
                                    double complex[:] dependent_variables) nogil:
    """ Interpolation function for complex numbers.

    Provided a domain, `desired_x` and a dependent array `dependent_values` search domain for value closest to 
    `desired_x` and return the value of `dependent_values` at that location if it is defined. Otherwise, use local 
    slopes of `desired_x` and `dependent_values` to interpolate a value of `dependent_values` at `desired_x`.

    Based on `numpy`'s `interp` function.

    Parameters
    ----------
    desired_x : float
        Location where `dependent_variables` is desired.
    desired_x : np.ndarray[float]
        Domain to search for the correct location.
    dependent_values : np.ndarray[complex]
        Dependent values that are to be returned after search and interpolation.

    Returns
    -------
    result : complex
        Desired value of `dependent_values`.

    """
    
    cdef int lenx
    lenx = len(x_domain)
    # Note: Needs to be at least 3 item long array. Add exception here?
    
    cdef double complex left_value
    left_value = dependent_variables[0]
    cdef double complex right_value
    right_value = dependent_variables[lenx - 1]
    
    # Binary Search with Guess
    cdef int i, j
    j = 0
    cdef double slope_real
    cdef double slope_imag
    cdef double x_slope_inverse
    
    cdef double result_real
    cdef double result_imag
    cdef double fp_at_j_real
    cdef double fp_at_j_imag
    cdef double xp_at_j
    cdef double fp_at_jp1_real
    cdef double fp_at_jp1_imag
    cdef double xp_at_jp1

    # Perform binary search with guess
    j = binary_search_with_guess(desired_x, x_domain, lenx, j)
    
    if j == -1:
        result_real = left_value.real
        result_imag = left_value.imag
    elif j == lenx:
        result_real = right_value.real
        result_imag = right_value.imag
    else:
        fp_at_j_real = dependent_variables[j].real
        fp_at_j_imag = dependent_variables[j].imag
        xp_at_j = x_domain[j]
        if j == lenx - 1:
            result_real = fp_at_j_real
            result_imag = fp_at_j_imag
        elif xp_at_j == desired_x:
            result_real = fp_at_j_real
            result_imag = fp_at_j_imag
        else:
            fp_at_jp1_real = dependent_variables[j + 1].real
            fp_at_jp1_imag = dependent_variables[j + 1].imag
            xp_at_jp1 = x_domain[j + 1]
            x_slope_inverse = 1.0 / (xp_at_jp1 - xp_at_j)
            slope_real = (fp_at_jp1_real - fp_at_j_real) * x_slope_inverse
            slope_imag = (fp_at_jp1_imag - fp_at_j_imag) * x_slope_inverse
    
            # If we get nan in one direction try the other
            # Real Part
            result_real = slope_real * (desired_x - xp_at_j) + fp_at_j_real
            if isnan(result_real):
                result_real = slope_real * (desired_x - xp_at_jp1) + fp_at_jp1_real
                if isnan(result_real) and (fp_at_jp1_real == fp_at_j_real):
                    result_real = fp_at_j_real
            
            # Imaginary Part
            result_imag = slope_imag * (desired_x - xp_at_j) + fp_at_j_imag
            if isnan(result_imag):
                result_imag = slope_imag * (desired_x - xp_at_jp1) + fp_at_jp1_imag
                if isnan(result_imag) and (fp_at_jp1_imag == fp_at_j_imag):
                    result_imag = fp_at_j_imag
                    
    cdef double complex result
    result = result_real + 1.0j * result_imag
                    
    return result
