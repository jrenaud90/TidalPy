""" Numba-safe implementation of np.trapz

"""
import numpy as np
from ..performance import njit, bool_
from numba import types

# TODO: Development Stopped.
"""
@njit(cacheable=True)
def diff(a, n=1, axis=-1, is_bool: bool = False):
    if n == 0:
        return a
    if n < 0:
        raise Exception("order must be non-negative")

    nd = a.ndim
    if nd == 0:
        raise Exception("diff requires input that is at least one dimensional")

    # Check is axis is last axis
    if axis == -1:
        axis = nd - 1

    # Numba does not like a lot of the slicing operations that the normal numpy code uses. So we have to get
    #    creative and use some flattening.
    dim_length = a.shape

    a_flat = a.flatten()

    # Make slices that cover the whole flattened array
    slice1 = np.ones_like(a_flat, dtype=np.bool)
    slice2 = np.ones_like(a_flat, dtype=np.bool)

    # Find the 1 x nd-1 value of each array that we need to remove.
    value_1_index = (axis + 1) * 2





    if axis == 0 or (nd == 1 and axis == -1):
        slice1_0 = (slice(1, None),)
        slice2_0 = (slice(None, -1),)
    else:
        slice1_0 = (slice(None),)
        slice2_0 = (slice(None),)

    for i in range(1, nd):

        if ((i == nd - 1) and axis == -1) or (i == axis):
            slice1_1 = slice1_0 + (slice(1, None),)
            slice2_1 = slice2_0 + (slice(None, -1),)
        else:
            slice1_1 = slice1_0 + (slice(None),)
            slice2_1 = slice2_0 + (slice(None),)

        slice1_0 = slice1_1
        slice2_0 = slice2_1

    # slice1_tuple = types.Tuple()
    # slice1 = tuple(slice1)
    # slice2 = tuple(slice2)

    a1 = np.ascontiguousarray(a[slice1_0])
    a2 = np.ascontiguousarray(a[slice2_0])

    if is_bool:
        result = np.not_equal(a1, a2)
    else:
        result = np.subtract(a1, a2)

    return result

@njit(cacheable=True)
def trapz(y, x: np.ndarray = None, dx: float = 1., axis: int = -1):


    # Determine dx
    if x is None:
        d = dx
    else:
        if x.ndim == 1:
            d = diff(x)
            # reshape to correct shape
            shape = [1] * y.ndim
            shape[axis] = d.shape[0]
            d = d.reshape(shape)
        else:
            d = diff(x, axis=axis)

    nd = y.ndim
    slice1 = [slice(None)] * nd
    slice2 = [slice(None)] * nd
    slice1[axis] = slice(1, None)
    slice2[axis] = slice(None, -1)

    # Complete integration
    result = (d * (y[tuple(slice1)] + y[tuple(slice2)]) / 2.0).sum(axis)

    return result
"""