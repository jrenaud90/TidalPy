#include "love_.hpp"

void find_love_cf(
        double* complex_love_numbers_ptr,  // These are double pointers that pointer to a double complex array.
        double* surface_solutions_ptr,     // Same as above
        double surface_gravity)
{
    const double y1_real = surface_solutions_ptr[0];
    const double y1_imag = surface_solutions_ptr[1];
    const double y3_real = surface_solutions_ptr[4];
    const double y3_imag = surface_solutions_ptr[5];
    const double y5_real = surface_solutions_ptr[8];
    const double y5_imag = surface_solutions_ptr[9];

    // Calculate Love and Shida numbers
    // Note Im(k2) = -Im(y5) (Henning & Hurford 2014 eq. A9), opposite convention of Tobie et al. (2005, eqs. 9 & 36)
    // And k2 = |-y5-1| (Roberts & Nimmo 2008 equation A8), not 1-y5 as in Henning & Hurford (2014) equation A9.
     
    // Okay, some clarification on this. It looks like VS04 that HH14 is based on used a different convention for y5,
    //      Tobie05's y5 = -y5 of SV; we follow that format here.

    // Love k
    complex_love_numbers_ptr[0] = y5_real - 1.0;
    complex_love_numbers_ptr[1] = y5_imag;

    // Love h
    complex_love_numbers_ptr[2] = y1_real * surface_gravity;
    complex_love_numbers_ptr[3] = y1_imag * surface_gravity;

    // Shida l
    complex_love_numbers_ptr[4] = y3_real * surface_gravity;
    complex_love_numbers_ptr[5] = y3_imag * surface_gravity;
}
