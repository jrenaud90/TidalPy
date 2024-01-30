#include <complex>


void find_love_cf(
        std::complex<double>* complex_love_numbers_ptr,
        std::complex<double>* surface_solutions_ptr,
        double surface_gravity) {

    // Extract the required radial solution values, isolating the real and imaginary portions.
    std::complex<double> y1 = surface_solutions_ptr[0];
    std::complex<double> y3 = surface_solutions_ptr[2];
    std::complex<double> y5 = surface_solutions_ptr[4];
    std::complex<double> c_surface_gravity = std::complex<double>(surface_gravity, 0.0);

    // Calculate Love and Shida numbers
    // Note Im(k2) = -Im(y5) (Henning & Hurford 2014 eq. A9), opposite convention of Tobie et al. (2005, eqs. 9 & 36)
    // And k2 = |-y5-1| (Roberts & Nimmo 2008 equation A8), not 1-y5 as in Henning & Hurford (2014) equation A9.
    // 
    // Okay, some clarification on this. It looks like VS04 that HH14 is based on used a different convention for y5,
    //     Tobie05's y5 = -y5 of SV; we follow that format here.
    // Love k
    complex_love_numbers_ptr[0] = std::complex<double>(std::real<double>(y5) - 1.0, std::imag<double>(y5));
    // Love h
    complex_love_numbers_ptr[1] = y1 * c_surface_gravity;
    // Shida l
    complex_love_numbers_ptr[2] = y3 * c_surface_gravity;
}
