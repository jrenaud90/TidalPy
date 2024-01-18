#include <math.h>
#include <complex.h>

#ifdef _MSC_VER // using Microsoft's implementation:
    // there is also _Fcomplex for float complex and _Lcomplex for long double complex
    typedef _Dcomplex double_complex;
    #define mk_dbl_cmplx _Cbuild
    // we need this on MSVC since it does not support the '+, *, -, /' operators
    #define cMult(a,b) _Cmulcc(a,b)
    #define cAdd(a,b) mk_dbl_cmplx(creal(a)+creal(b),cimag(a)+cimag(b))
#else // using standard C:
    typedef double complex double_complex;
    #define mk_dbl_cmplx CMPLX
    #define cMult(a,b) a*b
    #define cAdd(a,b) a+b
#endif

void find_love_cf(
        double_complex* complex_love_numbers_ptr,
        double_complex* surface_solutions_ptr,
        double surface_gravity) {

    // Extract the required radial solution values, isolating the real and imaginary portions.
    double_complex y1 = surface_solutions_ptr[0];
    double_complex y3 = surface_solutions_ptr[2];
    double_complex y5 = surface_solutions_ptr[4];
    double_complex c_surface_gravity = mk_dbl_cmplx(surface_gravity, 0.);

    // Calculate Love and Shida numbers
    // Note Im(k2) = -Im(y5) (Henning & Hurford 2014 eq. A9), opposite convention of Tobie et al. (2005, eqs. 9 & 36)
    // And k2 = |-y5-1| (Roberts & Nimmo 2008 equation A8), not 1-y5 as in Henning & Hurford (2014) equation A9.
    // 
    // Okay, some clarification on this. It looks like VS04 that HH14 is based on used a different convention for y5,
    //     Tobie05's y5 = -y5 of SV; we follow that format here.
    // Love k
    complex_love_numbers_ptr[0] = mk_dbl_cmplx(creal(y5) - 1.0, cimag(y5));
    // Love h
    complex_love_numbers_ptr[1] = cMult(y1, c_surface_gravity);
    // Shida l
    complex_love_numbers_ptr[2] = cMult(y3, c_surface_gravity);
}
