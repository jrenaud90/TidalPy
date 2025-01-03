# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.math cimport tgamma

from TidalPy.constants cimport d_NAN_DBL

cdef double[51] pre_calculated_doubles
cdef double* pre_calculated_doubles_ptr = &pre_calculated_doubles[0]
pre_calculated_doubles_ptr[  0] = 1.00
pre_calculated_doubles_ptr[  1] = 1.00
pre_calculated_doubles_ptr[  2] = 2.00
pre_calculated_doubles_ptr[  3] = 3.00
pre_calculated_doubles_ptr[  4] = 8.00
pre_calculated_doubles_ptr[  5] = 15.00
pre_calculated_doubles_ptr[  6] = 48.00
pre_calculated_doubles_ptr[  7] = 105.00
pre_calculated_doubles_ptr[  8] = 384.00
pre_calculated_doubles_ptr[  9] = 945.00
pre_calculated_doubles_ptr[ 10] = 3840.00
pre_calculated_doubles_ptr[ 11] = 10395.00
pre_calculated_doubles_ptr[ 12] = 46080.00
pre_calculated_doubles_ptr[ 13] = 135135.00
pre_calculated_doubles_ptr[ 14] = 645120.00
pre_calculated_doubles_ptr[ 15] = 2027025.00
pre_calculated_doubles_ptr[ 16] = 10321920.00
pre_calculated_doubles_ptr[ 17] = 34459425.00
pre_calculated_doubles_ptr[ 18] = 185794560.00
pre_calculated_doubles_ptr[ 19] = 654729075.00
pre_calculated_doubles_ptr[ 20] = 3715891200.00
pre_calculated_doubles_ptr[ 21] = 13749310575.00
pre_calculated_doubles_ptr[ 22] = 81749606400.00
pre_calculated_doubles_ptr[ 23] = 316234143225.00
pre_calculated_doubles_ptr[ 24] = 1961990553600.00
pre_calculated_doubles_ptr[ 25] = 7905853580625.01
pre_calculated_doubles_ptr[ 26] = 51011754393599.96
pre_calculated_doubles_ptr[ 27] = 213458046676875.16
pre_calculated_doubles_ptr[ 28] = 1428329123020799.25
pre_calculated_doubles_ptr[ 29] = 6190283353629379.00
pre_calculated_doubles_ptr[ 30] = 42849873690623960.00
pre_calculated_doubles_ptr[ 31] = 191898783962510816.00
pre_calculated_doubles_ptr[ 32] = 1371195958099966720.00
pre_calculated_doubles_ptr[ 33] = 6332659870762856448.00
pre_calculated_doubles_ptr[ 34] = 46620662575398879232.00
pre_calculated_doubles_ptr[ 35] = 221643095476699824128.00
pre_calculated_doubles_ptr[ 36] = 1678343852714360832000.00
pre_calculated_doubles_ptr[ 37] = 8200794532637892935680.00
pre_calculated_doubles_ptr[ 38] = 63777066403145720004608.00
pre_calculated_doubles_ptr[ 39] = 319830986772877752139776.00
pre_calculated_doubles_ptr[ 40] = 2551082656125828464640000.00
pre_calculated_doubles_ptr[ 41] = 13113070457687983475654656.00
pre_calculated_doubles_ptr[ 42] = 107145471557284812694749184.00
pre_calculated_doubles_ptr[ 43] = 563862029680583787669356544.00
pre_calculated_doubles_ptr[ 44] = 4714400748520528253875650560.00
pre_calculated_doubles_ptr[ 45] = 25373791335626273400058544128.00
pre_calculated_doubles_ptr[ 46] = 216862434431944368947512475648.00
pre_calculated_doubles_ptr[ 47] = 1192568192774434172503588864000.00
pre_calculated_doubles_ptr[ 48] = 10409396852733329709480598831104.00
pre_calculated_doubles_ptr[ 49] = 58435841445947288807899666579456.00
pre_calculated_doubles_ptr[ 50] = 520469842636666553028024352112640.00


cdef double cf_double_factorial(int n) noexcept nogil:

    if n < 51:
        # Use precalculated doubles
        return pre_calculated_doubles_ptr[n]
    elif n >= 51 and n < 171:
        # Calculate using recursion
        return tgamma(<double>n + 1.) / cf_double_factorial(n - 1)
    else:
        # ValueError('C function `tgamma` experiences overflow for l > 170.')
        return d_NAN_DBL


def double_factorial(int n):
    
    if n >= 171:
        raise ValueError('C function `tgamma` experiences overflow for l > 170.')

    return cf_double_factorial(n)
