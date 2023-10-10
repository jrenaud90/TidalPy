from libc.math cimport tgamma

cdef double[51] pre_calculated_doubles
cdef double* pre_calculated_doubles_ptr = &pre_calculated_doubles[0]
pre_calculated_doubles[  0] = 1.00
pre_calculated_doubles[  1] = 1.00
pre_calculated_doubles[  2] = 2.00
pre_calculated_doubles[  3] = 3.00
pre_calculated_doubles[  4] = 8.00
pre_calculated_doubles[  5] = 15.00
pre_calculated_doubles[  6] = 48.00
pre_calculated_doubles[  7] = 105.00
pre_calculated_doubles[  8] = 384.00
pre_calculated_doubles[  9] = 945.00
pre_calculated_doubles[ 10] = 3840.00
pre_calculated_doubles[ 11] = 10395.00
pre_calculated_doubles[ 12] = 46080.00
pre_calculated_doubles[ 13] = 135135.00
pre_calculated_doubles[ 14] = 645120.00
pre_calculated_doubles[ 15] = 2027025.00
pre_calculated_doubles[ 16] = 10321920.00
pre_calculated_doubles[ 17] = 34459425.00
pre_calculated_doubles[ 18] = 185794560.00
pre_calculated_doubles[ 19] = 654729075.00
pre_calculated_doubles[ 20] = 3715891200.00
pre_calculated_doubles[ 21] = 13749310575.00
pre_calculated_doubles[ 22] = 81749606400.00
pre_calculated_doubles[ 23] = 316234143225.00
pre_calculated_doubles[ 24] = 1961990553600.00
pre_calculated_doubles[ 25] = 7905853580625.01
pre_calculated_doubles[ 26] = 51011754393599.96
pre_calculated_doubles[ 27] = 213458046676875.16
pre_calculated_doubles[ 28] = 1428329123020799.25
pre_calculated_doubles[ 29] = 6190283353629379.00
pre_calculated_doubles[ 30] = 42849873690623960.00
pre_calculated_doubles[ 31] = 191898783962510816.00
pre_calculated_doubles[ 32] = 1371195958099966720.00
pre_calculated_doubles[ 33] = 6332659870762856448.00
pre_calculated_doubles[ 34] = 46620662575398879232.00
pre_calculated_doubles[ 35] = 221643095476699824128.00
pre_calculated_doubles[ 36] = 1678343852714360832000.00
pre_calculated_doubles[ 37] = 8200794532637892935680.00
pre_calculated_doubles[ 38] = 63777066403145720004608.00
pre_calculated_doubles[ 39] = 319830986772877752139776.00
pre_calculated_doubles[ 40] = 2551082656125828464640000.00
pre_calculated_doubles[ 41] = 13113070457687983475654656.00
pre_calculated_doubles[ 42] = 107145471557284812694749184.00
pre_calculated_doubles[ 43] = 563862029680583787669356544.00
pre_calculated_doubles[ 44] = 4714400748520528253875650560.00
pre_calculated_doubles[ 45] = 25373791335626273400058544128.00
pre_calculated_doubles[ 46] = 216862434431944368947512475648.00
pre_calculated_doubles[ 47] = 1192568192774434172503588864000.00
pre_calculated_doubles[ 48] = 10409396852733329709480598831104.00
pre_calculated_doubles[ 49] = 58435841445947288807899666579456.00
pre_calculated_doubles[ 50] = 520469842636666553028024352112640.00


cdef double double_factorial(unsigned char n) nogil:

    if n < 51:
        # Use precalculated doubles
        return pre_calculated_doubles_ptr[n]
    elif n >= 51 and n < 171:
        # Calculate using recursion
        return tgamma(<double>n + 1.) / double_factorial(n - 1)
    else:
        raise ValueError('C function `tgamma` experiences overflow for l > 170.')


def double_factorial_(unsigned char n):
    return double_factorial(n)
