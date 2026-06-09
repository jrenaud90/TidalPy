# distutils: language = c++
"""
love.pxd
Cython declarations for TidalPy's Love numbers container.

Exports c_LoveNumbers and the Python wrapper LoveNumbers so other extensions
can cimport and work with Love numbers at C speed.

Usage::

    from TidalPy.Tides_x.love.love cimport LoveNumbers, c_LoveNumbers
"""

from libcpp.complex cimport complex as cpp_complex


# =====================================================================================================================
# C++ struct declaration
# =====================================================================================================================
cdef extern from "love_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_LoveNumbers:
        cpp_complex[double] k   # potential Love number           [dimensionless]
        cpp_complex[double] h   # radial displacement Love number [dimensionless]
        cpp_complex[double] l   # tangential displacement Love number [dimensionless]
        c_LoveNumbers() except +
        c_LoveNumbers(cpp_complex[double], cpp_complex[double], cpp_complex[double]) except +


# =====================================================================================================================
# Cython wrapper class declaration
# =====================================================================================================================
cdef class LoveNumbers:
    cdef c_LoveNumbers _love
    cpdef dict to_dict(self)
