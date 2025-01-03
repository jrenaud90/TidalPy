# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
Modified from https://github.com/numpy/numpy/blob/main/numpy/core/src/npymath/npy_math_complex.c.src */
/*==========================================================
 * Helper functions
 *
 * These are necessary because we do not count on using a
 * C99 compiler.
 *=========================================================*/
"""

from libc.math cimport isfinite, isinf, isnan, copysign, \
    sqrt, fabs, signbit, exp, cos, sin, log, log1p, ldexp, atan2, frexp, ceil

cdef int DBL_MANT_DIG_INT = <int>d_DBL_MANT_DIG

SQRT2 = 1.414213562373095048801688724209698079  # sqrt 2
LOGE2 = 0.693147180559945309417232121458176568  # log_e 2

# We risk spurious overflow for components >= d_DBL_MAX / (1 + sqrt(2)).
SQRT2_INV = 1. / (1.0 + SQRT2)
THRESH    = SQRT2_INV * d_DBL_MAX
DBL_MAX_4 = 0.25 * d_DBL_MAX

# scaled_cexp precison constant
#if @precision@ == 1
# precision for float
SCALED_CEXP_K_F  = 235
# precision for double
SCALED_CEXP_K_D  = 1799
# precision for long double
SCALED_CEXP_K_LD = 19547

SCALED_K_LOGE2_D = SCALED_CEXP_K_D * LOGE2

SCALED_CEXP_LOWERF = 88.722839
SCALED_CEXP_UPPERF = 192.69492
SCALED_CEXP_LOWER  = 710.47586007394386
SCALED_CEXP_UPPER  = 1454.9159319953251
SCALED_CEXP_LOWERL = 11357.216553474703895
SCALED_CEXP_UPPERL = 22756.021937783004509


cdef inline double complex cf_build_dblcmplx(const double a, const double b) noexcept nogil:
    cdef double complex result
    cdef double* result_dbl_ptr = <double *> &result
    result_dbl_ptr[0] = a
    result_dbl_ptr[1] = b

    return result

cdef double cf_cabs(double complex z) noexcept nogil:
    cdef double z_real = z.real
    cdef double z_imag = z.imag
    return sqrt((z_real * z_real) + (z_imag * z_imag))

cdef double cf_cabs2(double complex z) noexcept nogil:
    cdef double z_real = z.real
    cdef double z_imag = z.imag
    return (z_real * z_real) + (z_imag * z_imag)

cdef double cf_carg(double complex z) noexcept nogil:
    return atan2(z.imag, z.real)

cdef double complex cf_cinv(double complex z) noexcept nogil:

    cdef double z_real = z.real
    cdef double z_imag = z.imag
    cdef double denom = ((z_real * z_real) + (z_imag * z_imag))
    cdef double complex inv

    # Check for extreme values
    if denom == 0.:
        # This is to match python's behavior, it is not strictly accurate depending on your definition.
        inv = cf_build_dblcmplx(d_INF_DBL, d_NAN_DBL)
    else:
        inv = cf_build_dblcmplx(z_real / denom, -z_imag / denom)

    return inv

cdef double cf_hypot(const double x, const double y) noexcept nogil:

    cdef double yx
    cdef double temp
    cdef double x_abs, y_abs

    if isinf(x) or isinf(y):
        return d_INF_DBL
    
    if isnan(x) or isnan(y):
        return d_NAN_DBL

    x_abs = fabs(x)
    y_abs = fabs(y)
    if x_abs < y_abs:
        temp  = x_abs
        x_abs = y_abs
        y_abs = temp

    if x_abs == 0.:
        return 0.
    else:
        yx = y_abs / x_abs
        return x_abs * sqrt(1. + yx * yx)


cdef double complex cf_csqrt(const double complex z) noexcept nogil:
    cdef double complex result
    cdef double t
    cdef char scale

    cdef double z_real = z.real
    cdef double z_imag = z.imag

    # Handle special cases.
    if z_imag == 0.0:
        if z_real == 0.0:
            return cf_build_dblcmplx(0., 0.)
        elif z_real > 0.0:
            return cf_build_dblcmplx(sqrt(z_real), 0.0)

    if isinf(z_imag):
        # Return inf +- inf (where sign is the same as input)
        return cf_build_dblcmplx(d_INF_DBL, z_imag)

    if isnan(z_real):
        # Return d_NAN_DBL + d_NAN_DBL i
        return cf_build_dblcmplx(d_NAN_DBL, d_NAN_DBL)

    if isinf(z_real):
        # csqrt(-inf + d_NAN_DBL i) = d_NAN_DBL +- inf i
        # csqrt(-inf + y i)   = 0   +  inf i
        # csqrt(inf + d_NAN_DBL i)  = inf +  d_NAN_DBL i
        # csqrt(inf + y i)    = inf +  0 i
        
        if signbit(z_real):
            # Negative z_real
            if isnan(z_imag):
                return cf_build_dblcmplx(d_NAN_DBL, d_INF_DBL)
            else:
                return cf_build_dblcmplx(0., d_INF_DBL)
        else:
            if isnan(z_imag):
                return cf_build_dblcmplx(d_INF_DBL, d_NAN_DBL)
            else:
                return cf_build_dblcmplx(d_INF_DBL, 0.)

    # The remaining special case (b is d_NAN_DBL) is handled just fine by the normal code path below.
    # Scale to avoid overflow.
    if (fabs(z_real) >= THRESH) or (fabs(z_imag) >= THRESH):
        z_real *= 0.25
        z_imag *= 0.25
        scale = 1
    else:
        scale = 0

    # Algorithm 312, CACM vol 10, Oct 1967.
    if z_real >= 0.0:
        t = sqrt((z_real + cf_hypot(z_real, z_imag)) * 0.5)
        result = cf_build_dblcmplx(t, (z_imag / (2.0 * t)))
    else:
        t = sqrt((-z_real + cf_hypot(z_real, z_imag)) * 0.5)
        result = cf_build_dblcmplx((fabs(z_imag) / (2.0 * t)), copysign(t, z_imag))

    # Rescale.
    if scale == 1:
        return cf_build_dblcmplx(result.real * 2.0, result.imag)
    else:
        return result


cdef double complex cf_scaled_cexp(const double x, const double y, const int expt) noexcept nogil:

    cdef double mant, mantcos, mantsin
    cdef int ex, excos, exsin, expt_scaled
    
    mant = frexp(exp(x - SCALED_K_LOGE2_D), &ex)
    mantcos = frexp(cos(y), &excos)
    mantsin = frexp(sin(y), &exsin)

    expt_scaled = expt + ex + SCALED_CEXP_K_D
    return cf_build_dblcmplx(ldexp(mant * mantcos, expt_scaled + excos), ldexp(mant * mantsin, expt_scaled + exsin))


cdef double complex cf_cexp(const double complex z) noexcept nogil:
    cdef double x, c, s
    cdef double complex ret
    
    cdef double z_real = z.real
    cdef double z_imag = z.imag

    if isfinite(z_real):
        if (z_real >= SCALED_CEXP_LOWER) and (z_real <= SCALED_CEXP_UPPER):
            ret = cf_scaled_cexp(z_real, z_imag, 0)
        else:
            x = exp(z_real)
            c = cos(z_imag)
            s = sin(z_imag)

            if isfinite(z_imag):
                ret = cf_build_dblcmplx((x * c), (x * s))
            else:
                ret = cf_build_dblcmplx(d_NAN_DBL, copysign(d_NAN_DBL, z_imag))
    elif isnan(z_real):
        # z_real is d_NAN_DBL
        if z_imag == 0:
            ret = z
        else:
            ret = cf_build_dblcmplx(z_real, copysign(d_NAN_DBL, z_imag))
    else:
        # z_real is +- inf
        if z_real > 0:
            if z_imag == 0:
                ret = z
            elif isfinite(z_imag):
                c = cos(z_imag)
                s = sin(z_imag)

                ret = cf_build_dblcmplx((z_real * c), (z_real * s))
            else:
                # x = +inf, y = +-inf | d_NAN_DBL
                ret = cf_build_dblcmplx(z_real, d_NAN_DBL)
        else:
            if isfinite(z_imag):
                x = exp(z_real)
                c = cos(z_imag)
                s = sin(z_imag)

                ret = cf_build_dblcmplx((x * c), (x * s))
            else:
                # x = -inf, y = d_NAN_DBL | +i inf
                ret = cf_build_dblcmplx(0.0, 0.0)
    return ret


cdef double complex cf_clog(const double complex z) noexcept nogil:
    # algorithm from cpython, rev. d86f5686cef9
    # 
    # The usual formula for the real part is log(hypot(z.real, z.imag)).
    # There are four situations where this formula is potentially
    # problematic:
    # 
    # (1) the absolute value of z is subnormal.  Then hypot is subnormal,
    #  so has fewer than the usual number of bits of accuracy, hence may
    #  have large relative error.  This then gives a large absolute error
    #  in the log.  This can be solved by rescaling z by a suitable power
    #  of 2.
    # 
    # (2) the absolute value of z is greater than d_DBL_MAX (e.g. when both
    #  z.real and z.imag are within a factor of 1/sqrt(2) of d_DBL_MAX)
    #  Again, rescaling solves this.
    #  
    # (3) the absolute value of z is close to 1.  In this case it's
    #  difficult to achieve good accuracy, at least in part because a
    #  change of 1ulp in the real or imaginary part of z can result in a
    #  change of billions of ulps in the correctly rounded answer.
    #  
    # (4) z = 0.  The simplest thing to do here is to call the
    #  floating-point log with an argument of 0, and let its behaviour
    #  (returning -d_INF_DBL, signaling a floating-point exception, setting
    #  errno, or whatever) determine that of c_log.  So the usual formula
    #  is fine here.

    cdef double z_real = z.real
    cdef double z_imag = z.imag
    cdef double z_real_abs = fabs(z_real)
    cdef double z_imag_abs = fabs(z_imag)

    cdef double h, a_max, a_min
    cdef double r_real
    cdef double r_imag

    if (z_real_abs > DBL_MAX_4) or (z_imag_abs > DBL_MAX_4):
        r_real = log(cf_hypot(z_real_abs / 2., z_imag_abs / 2.)) + LOGE2
    
    elif (z_real_abs < d_DBL_MIN) and (z_imag_abs < d_DBL_MIN):
        if (z_real_abs > 0) or (z_imag_abs > 0):
            # catch cases where hypot(z_real_abs, z_imag_abs) is subnormal
            r_real = log(cf_hypot(
                ldexp(z_real_abs, DBL_MANT_DIG_INT),
                ldexp(z_imag_abs, DBL_MANT_DIG_INT)
                )) - d_DBL_MANT_DIG * LOGE2
        else:
            # log(+/-0 +/- 0i)
            # raise divide-by-zero floating point exception
            r_real = -1.0 / z_real
            r_real = copysign(r_real, -1)
            r_imag = cf_carg(z)
            return cf_build_dblcmplx(r_real, r_imag)
    else:
        h = cf_hypot(z_real_abs, z_imag_abs)
        if (0.71 <= h) and (h <= 1.73):
            if z_real_abs > z_imag_abs:
                a_max = z_real_abs
                a_min = z_imag_abs
            else:
                a_max = z_imag_abs
                a_min = z_real_abs
            r_real = log1p((a_max - 1) * (a_max + 1) + a_min * a_min) / 2.
        else:
            r_real = log(h)
    r_imag = cf_carg(z)
    return cf_build_dblcmplx(r_real, r_imag)


cdef double complex cf_cpow(const double complex a, const double complex b) noexcept nogil:
    cdef double a_real = a.real
    cdef double b_real = b.real
    cdef double a_imag = a.imag
    cdef double b_imag = b.imag
    cdef double complex loga
    cdef int mask, b_asint
    
    cdef double complex intermediate_value1, intermediate_value2, result

    # Checking if in a^b, if b is zero.
    # If a is not zero then by definition of logarithm a^0 is 1.
    # If a is also zero then 0^0 is best defined as 1.
    if b_real == 0. and b_imag == 0.:
        return 1.0 + 0.0j
    elif a_real == 0. and a_imag == 0.:
        # case 0^b
        #  If a is a complex zero (ai=ar=0), then the result depends 
        #  upon values of br and bi. The result is either:
        #  0 (in magnitude), undefined or 1.
        #  The later case is for br=bi=0 and independent of ar and ai
        #  but is handled above).
        if b_real > 0.:
            # If the real part of b is positive (br>0) then this is
            #  the zero complex with positive sign on both the
            #  real and imaginary part.
            return cf_build_dblcmplx(0.0, 0.0)
        else:
            # else we are in the case where the
            # real part of b is negative (br<0).
            # Here we should return a complex d_NAN_DBL
            # and raise FloatingPointError: invalid value...
            return cf_build_dblcmplx(d_NAN_DBL, d_NAN_DBL)
    if (b_imag == 0.0) and (b_real > -100) and (b_real < 100) and (ceil(b_real) == b_real) \
            and isfinite(b_real) and isfinite(b_imag):
        # b_real can be cast as an integer that is between -100 and 100.
        b_asint = <int> b_real

        if (b_asint == 1):
            return a
        elif (b_asint == 2):
            return a * a
        elif (b_asint == 3):
            return a * a * a
        elif (b_asint > -100 and b_asint < 100):
            mask = 1
            if (b_asint < 0):
                b_asint = -b_asint
            intermediate_value1 = 1.0
            intermediate_value2 = a
            while True:
                if b_asint & mask:
                    intermediate_value1 = intermediate_value1 * intermediate_value2
                mask <<= 1
                if (b_asint < mask) or (mask <= 0):
                    break
                intermediate_value2 = intermediate_value2 * intermediate_value2
            result = intermediate_value1
            if (b_real < 0):
                result = 1. / result
            return result
    else:
        loga = cf_clog(a)
        a_real = loga.real
        a_imag = loga.imag
        return cf_cexp(
            cf_build_dblcmplx(
                (a_real * b_real - a_imag * b_imag), 
                (a_real * b_imag + a_imag * b_real)
                )
            )


cdef double complex cf_cipow(const double complex a, const int b) noexcept nogil:
    cdef double a_real = a.real
    cdef double a_imag = a.imag
    cdef double complex loga
    cdef int mask
    cdef char negative_pow = 0
    cdef int b_abs = b

    # This used to be "signbit(<double> b)" but it does not really make sesne to convert to double here just to check sign.
    if b < 0:
        negative_pow = 1
        b_abs = -b
    
    cdef double complex intermediate_value1, intermediate_value2, result


    # Checking if in a^b, if b is zero.
    # If a is not zero then by definition of logarithm a^0 is 1.
    # If a is also zero then 0^0 is best defined as 1.
    if b == 0:
        return cf_build_dblcmplx(1.0, 0.0)
    elif a_real == 0. and a_imag == 0.:
        # case 0^b
        #  If a is a complex zero (ai=ar=0), then the result depends 
        #  upon values of br and bi. The result is either:
        #  0 (in magnitude), undefined or 1.
        #  The later case is for br=bi=0 and independent of ar and ai
        #  but is handled above).
        if b > 0:
            # If the real part of b is positive (br>0) then this is
            #  the zero complex with positive sign on both the
            #  real and imaginary part.
            return cf_build_dblcmplx(0.0, 0.0)
        else:
            # else we are in the case where the
            # real part of b is negative (br<0).
            # Here we should return a complex d_NAN_DBL
            # and raise FloatingPointError: invalid value...
            return cf_build_dblcmplx(d_NAN_DBL, d_NAN_DBL)

    if (b_abs < 100):
        # b is between -100 and 100.
        if (b_abs == 1):
            if negative_pow:
                return 1. / a
            else:
                return a
        elif (b_abs == 2):
            if negative_pow:
                return 1. / (a * a)
            else:
                return a * a
        elif (b_abs == 3):
            if negative_pow:
                return 1. / (a * a * a)
            else:
                return a * a * a
        elif (b_abs < 100):
            mask = 1
            intermediate_value1 = cf_build_dblcmplx(1., 0.)
            intermediate_value2 = a
            while True:
                if b_abs & mask:
                    intermediate_value1 *= intermediate_value2
                mask <<= 1
                if (b_abs < mask) or (mask <= 0):
                    break
                intermediate_value2 *= intermediate_value2
            result = intermediate_value1
            if negative_pow:
                result = 1. / result
            return result
    else:
        loga   = cf_clog(a)
        a_real = loga.real
        a_imag = loga.imag
        return cf_cexp(cf_build_dblcmplx((a_real * b), (a_imag * b)))

########################################################################################################################
# Constants
########################################################################################################################
cdef double complex cmplx_NAN  = cf_build_dblcmplx(d_NAN_DBL, d_NAN_DBL)
cdef double complex cmplx_zero = cf_build_dblcmplx(0.0, 0.0)
cdef double complex cmplx_one  = cf_build_dblcmplx(1.0, 0.0)

########################################################################################################################
# Python wrappers
########################################################################################################################
def cinv(double complex z):
    return cf_cinv(z)

def cabs(double complex z):
    return cf_cabs(z)

def cabs2(double complex z):
    return cf_cabs2(z)

def hypot(const double x, const double y):
    return cf_hypot(x, y)

def csqrt(const double complex z):
    return cf_csqrt(z)

def scaled_cexp(const double x, const double y, const int expt):
    return cf_scaled_cexp(x, y, expt)

def cexp(const double complex z):
    return cf_cexp(z)

def clog(const double complex z):
    return cf_clog(z)

def cpow(const double complex a, const double complex b):
    return cf_cpow(a, b)

def cipow(const double complex a, const int b):
    return cf_cipow (a, b)
