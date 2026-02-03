#pragma once

#include "obliquity_common_.hpp"
#include "obliquity_func_l2_.hpp"
#include "obliquity_func_l3_.hpp"
#include "obliquity_func_l4_.hpp"
#include "obliquity_func_l5_.hpp"
#include "obliquity_func_l6_.hpp"
#include "obliquity_func_l7_.hpp"
#include "obliquity_func_l8_.hpp"
#include "obliquity_func_l9_.hpp"
#include "obliquity_func_l10_.hpp"


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
    case 3:
        switch (truncation)
        {
        case 0:
            // Obliquity is off
            return c_obliquity_function_l3_off(obliquity);
        case 2:
            // 2nd order truncation in obliquity
            return c_obliquity_function_l3_2(obliquity);
        case 4:
            // 4th order truncation in obliquity
            return c_obliquity_function_l3_4(obliquity);
        case 10:
            // 10 == the request for a general obliquity.
            return c_obliquity_function_l3_gen(obliquity);
        default:
            // Unsupported / Not implemented truncation provided.
            error_code_ptr[0] = -1;
            return ObliquityFuncOutput();
        }
    case 4:
        switch (truncation)
        {
        case 0:
            // Obliquity is off
            return c_obliquity_function_l4_off(obliquity);
        case 2:
            // 2nd order truncation in obliquity
            return c_obliquity_function_l4_2(obliquity);
        case 4:
            // 4th order truncation in obliquity
            return c_obliquity_function_l4_4(obliquity);
        case 10:
            // 10 == the request for a general obliquity.
            return c_obliquity_function_l4_gen(obliquity);
        default:
            // Unsupported / Not implemented truncation provided.
            error_code_ptr[0] = -1;
            return ObliquityFuncOutput();
        }
    case 5:
        switch (truncation)
        {
        case 0:
            // Obliquity is off
            return c_obliquity_function_l5_off(obliquity);
        case 2:
            // 2nd order truncation in obliquity
            return c_obliquity_function_l5_2(obliquity);
        case 4:
            // 4th order truncation in obliquity
            return c_obliquity_function_l5_4(obliquity);
        case 10:
            // 10 == the request for a general obliquity.
            return c_obliquity_function_l5_gen(obliquity);
        default:
            // Unsupported / Not implemented truncation provided.
            error_code_ptr[0] = -1;
            return ObliquityFuncOutput();
        }
    case 6:
        switch (truncation)
        {
        case 0:
            // Obliquity is off
            return c_obliquity_function_l6_off(obliquity);
        case 2:
            // 2nd order truncation in obliquity
            return c_obliquity_function_l6_2(obliquity);
        case 4:
            // 4th order truncation in obliquity
            return c_obliquity_function_l6_4(obliquity);
        case 10:
            // 10 == the request for a general obliquity.
            return c_obliquity_function_l6_gen(obliquity);
        default:
            // Unsupported / Not implemented truncation provided.
            error_code_ptr[0] = -1;
            return ObliquityFuncOutput();
        }
    case 7:
        switch (truncation)
        {
        case 0:
            // Obliquity is off
            return c_obliquity_function_l7_off(obliquity);
        case 2:
            // 2nd order truncation in obliquity
            return c_obliquity_function_l7_2(obliquity);
        case 4:
            // 4th order truncation in obliquity
            return c_obliquity_function_l7_4(obliquity);
        case 10:
            // 10 == the request for a general obliquity.
            return c_obliquity_function_l7_gen(obliquity);
        default:
            // Unsupported / Not implemented truncation provided.
            error_code_ptr[0] = -1;
            return ObliquityFuncOutput();
        }
    case 8:
        switch (truncation)
        {
        case 0:
            // Obliquity is off
            return c_obliquity_function_l8_off(obliquity);
        case 2:
            // 2nd order truncation in obliquity
            return c_obliquity_function_l8_2(obliquity);
        case 4:
            // 4th order truncation in obliquity
            return c_obliquity_function_l8_4(obliquity);
        case 10:
            // 10 == the request for a general obliquity.
            return c_obliquity_function_l8_gen(obliquity);
        default:
            // Unsupported / Not implemented truncation provided.
            error_code_ptr[0] = -1;
            return ObliquityFuncOutput();
        }
    case 9:
        switch (truncation)
        {
        case 0:
            // Obliquity is off
            return c_obliquity_function_l9_off(obliquity);
        case 2:
            // 2nd order truncation in obliquity
            return c_obliquity_function_l9_2(obliquity);
        case 4:
            // 4th order truncation in obliquity
            return c_obliquity_function_l9_4(obliquity);
        case 10:
            // 10 == the request for a general obliquity.
            return c_obliquity_function_l9_gen(obliquity);
        default:
            // Unsupported / Not implemented truncation provided.
            error_code_ptr[0] = -1;
            return ObliquityFuncOutput();
        }
    case 10:
        switch (truncation)
        {
        case 0:
            // Obliquity is off
            return c_obliquity_function_l10_off(obliquity);
        case 2:
            // 2nd order truncation in obliquity
            return c_obliquity_function_l10_2(obliquity);
        case 4:
            // 4th order truncation in obliquity
            return c_obliquity_function_l10_4(obliquity);
        case 10:
            // 10 == the request for a general obliquity.
            return c_obliquity_function_l10_gen(obliquity);
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


