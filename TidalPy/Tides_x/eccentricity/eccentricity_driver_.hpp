#pragma once

#include "eccentricity_common_.hpp"
#include "eccentricity_func_l2_.hpp"
#include "eccentricity_func_l3_.hpp"
#include "eccentricity_func_l4_.hpp"
#include "eccentricity_func_l5_.hpp"
#include "eccentricity_func_l6_.hpp"
#include "eccentricity_func_l7_.hpp"
#include "eccentricity_func_l8_.hpp"
#include "eccentricity_func_l9_.hpp"
#include "eccentricity_func_l10_.hpp"


EccentricityFuncOutput c_eccentricity_func(
        int* error_code_ptr,
        double eccentricity,
        int degree_l,
        int truncation)
{
    error_code_ptr[0] = 0;

    switch (degree_l)
    {
    case 2:
        switch (truncation)
        {
        case 1:
            return c_eccentricity_function_l2_e1(eccentricity);
        case 2:
            return c_eccentricity_function_l2_e2(eccentricity);
        case 3:
            return c_eccentricity_function_l2_e3(eccentricity);
        case 4:
            return c_eccentricity_function_l2_e4(eccentricity);
        case 5:
            return c_eccentricity_function_l2_e5(eccentricity);
        case 10:
            return c_eccentricity_function_l2_e10(eccentricity);
        case 15:
            return c_eccentricity_function_l2_e15(eccentricity);
        case 20:
            return c_eccentricity_function_l2_e20(eccentricity);
        default:
            // Unsupported / Not implemented truncation provided.
            error_code_ptr[0] = -1;
            return EccentricityFuncOutput();
        }
    case 3:
        switch (truncation)
        {
        case 1:
            return c_eccentricity_function_l3_e1(eccentricity);
        case 2:
            return c_eccentricity_function_l3_e2(eccentricity);
        case 3:
            return c_eccentricity_function_l3_e3(eccentricity);
        case 4:
            return c_eccentricity_function_l3_e4(eccentricity);
        case 5:
            return c_eccentricity_function_l3_e5(eccentricity);
        case 10:
            return c_eccentricity_function_l3_e10(eccentricity);
        case 15:
            return c_eccentricity_function_l3_e15(eccentricity);
        case 20:
            return c_eccentricity_function_l3_e20(eccentricity);
        default:
            // Unsupported / Not implemented truncation provided.
            error_code_ptr[0] = -1;
            return EccentricityFuncOutput();
        }
    case 4:
        switch (truncation)
        {
        case 1:
            return c_eccentricity_function_l4_e1(eccentricity);
        case 2:
            return c_eccentricity_function_l4_e2(eccentricity);
        case 3:
            return c_eccentricity_function_l4_e3(eccentricity);
        case 4:
            return c_eccentricity_function_l4_e4(eccentricity);
        case 5:
            return c_eccentricity_function_l4_e5(eccentricity);
        case 10:
            return c_eccentricity_function_l4_e10(eccentricity);
        case 15:
            return c_eccentricity_function_l4_e15(eccentricity);
        case 20:
            return c_eccentricity_function_l4_e20(eccentricity);
        default:
            // Unsupported / Not implemented truncation provided.
            error_code_ptr[0] = -1;
            return EccentricityFuncOutput();
        }
    case 5:
        switch (truncation)
        {
        case 1:
            return c_eccentricity_function_l5_e1(eccentricity);
        case 2:
            return c_eccentricity_function_l5_e2(eccentricity);
        case 3:
            return c_eccentricity_function_l5_e3(eccentricity);
        case 4:
            return c_eccentricity_function_l5_e4(eccentricity);
        case 5:
            return c_eccentricity_function_l5_e5(eccentricity);
        case 10:
            return c_eccentricity_function_l5_e10(eccentricity);
        case 15:
            return c_eccentricity_function_l5_e15(eccentricity);
        case 20:
            return c_eccentricity_function_l5_e20(eccentricity);
        default:
            // Unsupported / Not implemented truncation provided.
            error_code_ptr[0] = -1;
            return EccentricityFuncOutput();
        }
    case 6:
        switch (truncation)
        {
        case 1:
            return c_eccentricity_function_l6_e1(eccentricity);
        case 2:
            return c_eccentricity_function_l6_e2(eccentricity);
        case 3:
            return c_eccentricity_function_l6_e3(eccentricity);
        case 4:
            return c_eccentricity_function_l6_e4(eccentricity);
        case 5:
            return c_eccentricity_function_l6_e5(eccentricity);
        case 10:
            return c_eccentricity_function_l6_e10(eccentricity);
        case 15:
            return c_eccentricity_function_l6_e15(eccentricity);
        case 20:
            return c_eccentricity_function_l6_e20(eccentricity);
        default:
            // Unsupported / Not implemented truncation provided.
            error_code_ptr[0] = -1;
            return EccentricityFuncOutput();
        }
    case 7:
        switch (truncation)
        {
        case 1:
            return c_eccentricity_function_l7_e1(eccentricity);
        case 2:
            return c_eccentricity_function_l7_e2(eccentricity);
        case 3:
            return c_eccentricity_function_l7_e3(eccentricity);
        case 4:
            return c_eccentricity_function_l7_e4(eccentricity);
        case 5:
            return c_eccentricity_function_l7_e5(eccentricity);
        case 10:
            return c_eccentricity_function_l7_e10(eccentricity);
        case 15:
            return c_eccentricity_function_l7_e15(eccentricity);
        case 20:
            return c_eccentricity_function_l7_e20(eccentricity);
        default:
            // Unsupported / Not implemented truncation provided.
            error_code_ptr[0] = -1;
            return EccentricityFuncOutput();
        }
    case 8:
        switch (truncation)
        {
        case 1:
            return c_eccentricity_function_l8_e1(eccentricity);
        case 2:
            return c_eccentricity_function_l8_e2(eccentricity);
        case 3:
            return c_eccentricity_function_l8_e3(eccentricity);
        case 4:
            return c_eccentricity_function_l8_e4(eccentricity);
        case 5:
            return c_eccentricity_function_l8_e5(eccentricity);
        case 10:
            return c_eccentricity_function_l8_e10(eccentricity);
        case 15:
            return c_eccentricity_function_l8_e15(eccentricity);
        case 20:
            return c_eccentricity_function_l8_e20(eccentricity);
        default:
            // Unsupported / Not implemented truncation provided.
            error_code_ptr[0] = -1;
            return EccentricityFuncOutput();
        }
    case 9:
        switch (truncation)
        {
        case 1:
            return c_eccentricity_function_l9_e1(eccentricity);
        case 2:
            return c_eccentricity_function_l9_e2(eccentricity);
        case 3:
            return c_eccentricity_function_l9_e3(eccentricity);
        case 4:
            return c_eccentricity_function_l9_e4(eccentricity);
        case 5:
            return c_eccentricity_function_l9_e5(eccentricity);
        case 10:
            return c_eccentricity_function_l9_e10(eccentricity);
        case 15:
            return c_eccentricity_function_l9_e15(eccentricity);
        case 20:
            return c_eccentricity_function_l9_e20(eccentricity);
        default:
            // Unsupported / Not implemented truncation provided.
            error_code_ptr[0] = -1;
            return EccentricityFuncOutput();
        }
    case 10:
        switch (truncation)
        {
        case 1:
            return c_eccentricity_function_l10_e1(eccentricity);
        case 2:
            return c_eccentricity_function_l10_e2(eccentricity);
        case 3:
            return c_eccentricity_function_l10_e3(eccentricity);
        case 4:
            return c_eccentricity_function_l10_e4(eccentricity);
        case 5:
            return c_eccentricity_function_l10_e5(eccentricity);
        case 10:
            return c_eccentricity_function_l10_e10(eccentricity);
        case 15:
            return c_eccentricity_function_l10_e15(eccentricity);
        case 20:
            return c_eccentricity_function_l10_e20(eccentricity);
        default:
            // Unsupported / Not implemented truncation provided.
            error_code_ptr[0] = -1;
            return EccentricityFuncOutput();
        }
    default:
        // Unsupported / Not implemented degree l provided.
        error_code_ptr[0] = -2;
        return EccentricityFuncOutput();
    }
}


