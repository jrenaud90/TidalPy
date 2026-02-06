#pragma once

#include "obliquity_common_.hpp"

ObliquityFuncOutput c_obliquity_function_l3_gen(double obliquity)
{
    // Inclination Functions Calculated for l = 3.
    // Functions are fully general with no assumptions on the magnitude of obliquity.

    //  Functions Calculated for l = 3.

    c_IntMap<c_Key3, double> result_by_lmp(16);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(4);
    // Optimizations
    double sin_i = std::sin(obliquity);
    double sin_i_2 = sin_i * sin_i;
    double sin_i_3 = sin_i_2 * sin_i;
    double cos_i_half = std::cos((1.0/2.0)*obliquity);
    double sin_i_half = std::sin((1.0/2.0)*obliquity);
    double cos_i_half_2 = cos_i_half * cos_i_half;
    double cos_i_half_4 = cos_i_half_2 * cos_i_half_2;
    double sin_i_half_2 = sin_i_half * sin_i_half;
    double cos_i = std::cos(obliquity);
    double cos_i_2 = cos_i * cos_i;
    double sin_i_half_4 = sin_i_half_2 * sin_i_half_2;
    double sin_i_half_5 = sin_i_half_4 * sin_i_half;
    double cos_i_half_3 = cos_i_half_2 * cos_i_half;
    double cos_i_half_6 = cos_i_half_3 * cos_i_half_3;
    double sin_i_half_3 = sin_i_half_2 * sin_i_half;
    double sin_i_half_6 = sin_i_half_3 * sin_i_half_3;

    c_IntMap<c_Key1, double> result_by_p(3);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (3, 0).
    // p = 0
    tmp_double = -0.3125*sin_i_3;
    result_by_lmp.set(c_Key3(3, 0, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 0.1875*(5.0*sin_i_2 - 4.0)*sin_i;
    result_by_lmp.set(c_Key3(3, 0, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 0.1875*(4.0 - 5.0*sin_i_2)*sin_i;
    result_by_lmp.set(c_Key3(3, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 0.3125*sin_i_3;
    result_by_lmp.set(c_Key3(3, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(3, 0), result_by_p);
    result_by_p.clear();

    // l , m = (3, 1).
    // p = 0
    tmp_double = -7.5*sin_i_half_2*cos_i_half_4;
    result_by_lmp.set(c_Key3(3, 1, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 0.375*(-15.0*cos_i_2 + 10.0*cos_i + 1.0)*cos_i_half_2;
    result_by_lmp.set(c_Key3(3, 1, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 0.375*(15.0*sin_i_2 - 10.0*cos_i - 14.0)*sin_i_half_2;
    result_by_lmp.set(c_Key3(3, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -7.5*sin_i_half_4*cos_i_half_2;
    result_by_lmp.set(c_Key3(3, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(3, 1), result_by_p);
    result_by_p.clear();

    // l , m = (3, 2).
    // p = 0
    tmp_double = 1.875*std::pow(cos_i + 1.0, 2)*sin_i;
    result_by_lmp.set(c_Key3(3, 2, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 1.875*(3.0*sin_i_2 - 2.0*cos_i - 2.0)*sin_i;
    result_by_lmp.set(c_Key3(3, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 15.0*sin_i_half_5*cos_i_half - 3.75*sin_i_3;
    result_by_lmp.set(c_Key3(3, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -15.0*sin_i_half_5*cos_i_half;
    result_by_lmp.set(c_Key3(3, 2, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(3, 2), result_by_p);
    result_by_p.clear();

    // l , m = (3, 3).
    // p = 0
    tmp_double = 15.0*cos_i_half_6;
    result_by_lmp.set(c_Key3(3, 3, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 45.0*sin_i_half_2*cos_i_half_4;
    result_by_lmp.set(c_Key3(3, 3, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 45.0*sin_i_half_4*cos_i_half_2;
    result_by_lmp.set(c_Key3(3, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 15.0*sin_i_half_6;
    result_by_lmp.set(c_Key3(3, 3, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(3, 3), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l3_2(double obliquity)
{
    // Inclination Functions Calculated for l = 3.
    // Functions are taylor expanded and truncated at obliquity^2.

    //  Functions Calculated for l = 3.

    c_IntMap<c_Key3, double> result_by_lmp(6);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(4);
    // Optimizations

    c_IntMap<c_Key1, double> result_by_p(3);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (3, 0).
    // p = 1
    tmp_double = -0.75*obliquity;
    result_by_lmp.set(c_Key3(3, 0, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 0.75*obliquity;
    result_by_lmp.set(c_Key3(3, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(3, 0), result_by_p);
    result_by_p.clear();

    // l , m = (3, 1).
    // p = 1
    tmp_double = -1.5;
    result_by_lmp.set(c_Key3(3, 1, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(3, 1), result_by_p);
    result_by_p.clear();

    // l , m = (3, 2).
    // p = 0
    tmp_double = 7.5*obliquity;
    result_by_lmp.set(c_Key3(3, 2, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -7.5*obliquity;
    result_by_lmp.set(c_Key3(3, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(3, 2), result_by_p);
    result_by_p.clear();

    // l , m = (3, 3).
    // p = 0
    tmp_double = 15.0;
    result_by_lmp.set(c_Key3(3, 3, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(3, 3), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l3_4(double obliquity)
{
    // Inclination Functions Calculated for l = 3.
    // Functions are taylor expanded and truncated at obliquity^4.

    //  Functions Calculated for l = 3.

    c_IntMap<c_Key3, double> result_by_lmp(12);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(4);
    // Optimizations
    double obliquity_2 = obliquity * obliquity;
    double obliquity_3 = obliquity_2 * obliquity;

    c_IntMap<c_Key1, double> result_by_p(3);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (3, 0).
    // p = 0
    tmp_double = -0.3125*obliquity_3;
    result_by_lmp.set(c_Key3(3, 0, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 1.0625*obliquity_3 - 0.75*obliquity;
    result_by_lmp.set(c_Key3(3, 0, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -1.0625*obliquity_3 + 0.75*obliquity;
    result_by_lmp.set(c_Key3(3, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 0.3125*obliquity_3;
    result_by_lmp.set(c_Key3(3, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(3, 0), result_by_p);
    result_by_p.clear();

    // l , m = (3, 1).
    // p = 0
    tmp_double = -1.875*obliquity_2;
    result_by_lmp.set(c_Key3(3, 1, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 4.125*obliquity_2 - 1.5;
    result_by_lmp.set(c_Key3(3, 1, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -2.25*obliquity_2;
    result_by_lmp.set(c_Key3(3, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(3, 1), result_by_p);
    result_by_p.clear();

    // l , m = (3, 2).
    // p = 0
    tmp_double = -5.0*obliquity_3 + 7.5*obliquity;
    result_by_lmp.set(c_Key3(3, 2, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 8.75*obliquity_3 - 7.5*obliquity;
    result_by_lmp.set(c_Key3(3, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -3.75*obliquity_3;
    result_by_lmp.set(c_Key3(3, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(3, 2), result_by_p);
    result_by_p.clear();

    // l , m = (3, 3).
    // p = 0
    tmp_double = 15.0 - 11.25*obliquity_2;
    result_by_lmp.set(c_Key3(3, 3, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 11.25*obliquity_2;
    result_by_lmp.set(c_Key3(3, 3, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(3, 3), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l3_off(double obliquity)
{
    // Inclination Functions Calculated for l = 3.
    // Functions assume `obliquity=0.0`. However some modes are still activated (non-zero `F(l,m,p)`).

    //  Functions Calculated for l = 3.

    c_IntMap<c_Key3, double> result_by_lmp(2);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(2);
    // Optimizations

    c_IntMap<c_Key1, double> result_by_p(3);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (3, 1).
    // p = 1
    tmp_double = -1.5;
    result_by_lmp.set(c_Key3(3, 1, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(3, 1), result_by_p);
    result_by_p.clear();

    // l , m = (3, 3).
    // p = 0
    tmp_double = 15.0;
    result_by_lmp.set(c_Key3(3, 3, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(3, 3), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}
