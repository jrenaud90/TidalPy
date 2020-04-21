from typing import Union

import numpy as np

# General helpers
NoneType = type(None)

# Numpy and Array type lists
float_like = [float, np.float, np.float64]
floatarray_like = [float, np.float, np.float64, np.ndarray]
int_like = [int, np.int, np.int64]

numpy_float_info = np.finfo(dtype=np.float)
float_eps = numpy_float_info.eps
float_max = numpy_float_info.max
float_min = numpy_float_info.min
float_log10_max = np.log10(float_max)
float_lognat_max = np.log(float_max)

# Other type lists
list_like = [list, tuple, set]

# Type-annotation helpers
NumericalType = Union[float, int, complex, np.float, np.int, np.complex]
FloatArray = Union[float, np.float, np.float64, np.ndarray]
ComplexArray = Union[complex, np.complex, np.complex128, np.ndarray]
NumArray = Union[FloatArray, ComplexArray]
TupleNone = Union[tuple, NoneType]
ArrayNone = Union[np.ndarray, NoneType]