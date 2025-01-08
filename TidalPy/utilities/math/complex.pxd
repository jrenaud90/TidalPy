from TidalPy.constants cimport d_DBL_MAX, d_DBL_MIN, d_DBL_MANT_DIG, d_NAN_DBL, d_INF_DBL

cdef double complex cmplx_NAN
cdef double complex cmplx_zero
cdef double complex cmplx_one
cdef double SQRT2
cdef double LOGE2
cdef double SQRT2_INV
cdef double THRESH
cdef double DBL_MAX_4
cdef int SCALED_CEXP_K_F
cdef int SCALED_CEXP_K_D
cdef int SCALED_CEXP_K_LD
cdef double SCALED_K_LOGE2_D
cdef float SCALED_CEXP_LOWERF
cdef float SCALED_CEXP_UPPERF
cdef double SCALED_CEXP_LOWER
cdef double SCALED_CEXP_UPPER
cdef long double SCALED_CEXP_LOWERL
cdef long double SCALED_CEXP_UPPERL


cdef double complex cf_build_dblcmplx(double a, double b) noexcept nogil

cdef double cf_cabs(double complex z) noexcept nogil

cdef double cf_cabs2(double complex z) noexcept nogil

cdef double complex cf_cinv(double complex z) noexcept nogil

cdef double cf_hypot(const double x, const double y) noexcept nogil

cdef double complex cf_csqrt(const double complex z) noexcept nogil

cdef double complex cf_scaled_cexp(const double x, const double y, const int expt) noexcept nogil

cdef double complex cf_cexp(const double complex z) noexcept nogil

cdef double complex cf_clog(const double complex z) noexcept nogil

cdef double complex cf_cpow(const double complex a, const double complex b) noexcept nogil

cdef double complex cf_cipow(const double complex a, const int b) noexcept nogil
