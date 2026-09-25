# distutils: language = c++
# cython: boundscheck=False, wraparound=False
"""Inline helpers shared by the physics-model wrappers: NumPy inputs into std::vector, std::vector results back into
NumPy, and the broadcasting the vectorized model calls accept.

The including module must call ``cnp.import_array()``. Inputs are read through const memoryviews, so read-only arrays
are accepted.
"""

from libc.string cimport memcpy
from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex
from libcpp.vector cimport vector

cimport numpy as cnp


cdef inline void cy_fill_vector(const double[::1] source, vector[double]& destination) noexcept nogil:
    cdef Py_ssize_t num_values = source.shape[0]
    destination.resize(num_values)
    if num_values > 0:
        memcpy(destination.data(), &source[0], num_values * sizeof(double))


cdef inline object cy_vector_to_ndarray(vector[double]& source, tuple shape):
    """A float64 array of ``shape`` holding a copy of ``source``."""
    cdef cnp.npy_intp num_values = <cnp.npy_intp>source.size()
    cdef cnp.ndarray result = cnp.PyArray_EMPTY(1, &num_values, cnp.NPY_FLOAT64, 0)
    if num_values > 0:
        memcpy(cnp.PyArray_DATA(result), source.data(), num_values * sizeof(double))
    return result.reshape(shape)


cdef inline object cy_complex_vector_to_ndarray(vector[cpp_complex[double]]& source, tuple shape):
    """A complex128 array of ``shape`` holding a copy of ``source``; std::complex<double> is laid out as two
    doubles, as complex128 is."""
    cdef cnp.npy_intp num_values = <cnp.npy_intp>source.size()
    cdef cnp.ndarray result = cnp.PyArray_EMPTY(1, &num_values, cnp.NPY_COMPLEX128, 0)
    if num_values > 0:
        memcpy(cnp.PyArray_DATA(result), source.data(), num_values * 2 * sizeof(double))
    return result.reshape(shape)


cdef inline object cy_broadcast_inputs(tuple values, vector[vector[double]]& destinations, cpp_bool flatten):
    """Prepare the inputs of an element-wise model call, each a float or an array, for a broadcast-aware C++ sweep.

    Returns None, leaving ``destinations`` untouched, when ``flatten`` is off and no input is an ndarray: the caller
    then makes one scalar call. Otherwise ``destinations`` gets one vector per input, in order, holding a single value
    for an input of one element and one value per point of the broadcast shape for any other, and the output shape is
    returned. ``flatten`` ravels every input first, so the output is one-dimensional. The shape has at least one
    dimension, so an all-0-d call gives a one-element array.

    Raises
    ------
    ValueError
        The input shapes do not broadcast together.
    """
    cdef object value
    cdef cpp_bool any_array = flatten
    if not any_array:
        for value in values:
            if isinstance(value, cnp.ndarray):
                any_array = True
                break
        if not any_array:
            return None

    # Fast path for the common call, floats and 1-D contiguous float64 arrays of one length: no NumPy broadcasting.
    cdef Py_ssize_t num_points = 1
    cdef Py_ssize_t array_length
    cdef cpp_bool simple = True
    cdef cnp.ndarray simple_array
    cdef const double[::1] simple_view
    for value in values:
        if isinstance(value, cnp.ndarray):
            simple_array = <cnp.ndarray>value
            if (cnp.PyArray_NDIM(simple_array) != 1 or cnp.PyArray_TYPE(simple_array) != cnp.NPY_FLOAT64
                    or not cnp.PyArray_IS_C_CONTIGUOUS(simple_array)):
                simple = False
                break
            array_length = cnp.PyArray_DIM(simple_array, 0)
            if array_length != 1:
                if num_points != 1 and array_length != num_points:
                    simple = False
                    break
                num_points = array_length
        elif not isinstance(value, float):
            simple = False
            break
    cdef Py_ssize_t value_i
    if simple:
        destinations.resize(len(values))
        for value_i in range(len(values)):
            value = values[value_i]
            if isinstance(value, cnp.ndarray):
                simple_view = value
                cy_fill_vector(simple_view, destinations[value_i])
            else:
                destinations[value_i].assign(1, <double>value)
        return (num_points,)

    # Deferred: a .pxd cannot import at module level, and the paths above need no NumPy call.
    import numpy as np
    cdef list arrays = [np.asarray(value, dtype=np.float64) for value in values]
    if flatten:
        arrays = [array.reshape(-1) for array in arrays]
    cdef tuple shape = np.broadcast_shapes(*[array.shape for array in arrays])
    if len(shape) == 0:
        shape = (1,)

    cdef Py_ssize_t num_inputs = len(arrays)
    cdef Py_ssize_t input_i
    cdef object array
    cdef const double[::1] view
    destinations.resize(num_inputs)
    for input_i in range(num_inputs):
        array = arrays[input_i]
        if array.size == 1:
            view = np.ascontiguousarray(array).reshape(-1)
        else:
            view = np.ascontiguousarray(np.broadcast_to(array, shape)).reshape(-1)
        cy_fill_vector(view, destinations[input_i])
    return shape
