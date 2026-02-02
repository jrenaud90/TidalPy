#pragma once

#include "obliquity_common_.hpp"
#include "obliquity_func_l2_.hpp"


ObliquityFuncOutput c_obliquity_func(
        int* error_code_ptr,
        double obliquity,
        int degree_l,
        int truncation)
{
    error_code_ptr[0] = 0;

    switch (degree_l)
    {
    case 2:
        switch (truncation)
        {
        case 0:
            // Obliquity is off
            return c_obliquity_function_l2_off(obliquity);
        case 2:
            // 2nd order truncation in obliquity
            return c_obliquity_function_l2_2(obliquity);
        case 4:
            // 4th order truncation in obliquity
            return c_obliquity_function_l2_4(obliquity);
        case 10:
            // 10 == the request for a general obliquity.
            return c_obliquity_function_l2_gen(obliquity);
        default:
            // Unsupported / Not implemented truncation provided.
            error_code_ptr[0] = -1;
            return ObliquityFuncOutput();
        }
    default:
        // Unsupported / Not implemented degree l provided.
        error_code_ptr[0] = -2;
        return ObliquityFuncOutput();
    }
}


