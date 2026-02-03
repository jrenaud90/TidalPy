#pragma once

#include "obliquity_common_.hpp"

ObliquityFuncOutput c_obliquity_function_l2_gen(double obliquity)
{
    // Inclination Functions Calculated for l = 2.
    // Functions are fully general with no assumptions on the magnitude of obliquity.

    //  Functions Calculated for l = 2.

    c_IntMap<c_Key3, double> result_by_lmp(9);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(3);
    // Optimizations
    double sin_i = std::sin(obliquity);
    double sin_i_2 = sin_i * sin_i;
    double sin_i_half = std::sin((1.0/2.0)*obliquity);
    double sin_i_half_2 = sin_i_half * sin_i_half;
    double sin_i_half_4 = sin_i_half_2 * sin_i_half_2;
    double cos_i_half = std::cos((1.0/2.0)*obliquity);
    double cos_i_half_2 = cos_i_half * cos_i_half;
    double cos_i_half_3 = cos_i_half_2 * cos_i_half;
    double sin_i_double = std::sin(2*obliquity);
    double sin_i_half_3 = sin_i_half_2 * sin_i_half;
    double cos_i_half_4 = cos_i_half_2 * cos_i_half_2;

    c_IntMap<c_Key1, double> result_by_p(2);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (2, 0).
    // p = 0
    tmp_double = -0.375*sin_i_2;
    result_by_lmp.set(c_Key3(2, 0, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -sin_i_half_4 + sin_i_half_2 + 0.5*sin_i_2 - 0.5;
    result_by_lmp.set(c_Key3(2, 0, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -0.375*sin_i_2;
    result_by_lmp.set(c_Key3(2, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(2, 0), result_by_p);
    result_by_p.clear();

    // l , m = (2, 1).
    // p = 0
    tmp_double = 3.0*sin_i_half*cos_i_half_3;
    result_by_lmp.set(c_Key3(2, 1, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -0.75*sin_i_double;
    result_by_lmp.set(c_Key3(2, 1, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -3.0*sin_i_half_3*cos_i_half;
    result_by_lmp.set(c_Key3(2, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(2, 1), result_by_p);
    result_by_p.clear();

    // l , m = (2, 2).
    // p = 0
    tmp_double = 3.0*cos_i_half_4;
    result_by_lmp.set(c_Key3(2, 2, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 1.5*sin_i_2;
    result_by_lmp.set(c_Key3(2, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 3.0*sin_i_half_4;
    result_by_lmp.set(c_Key3(2, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(2, 2), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l2_2(double obliquity)
{
    // Inclination Functions Calculated for l = 2.
    // Functions are taylor expanded and truncated at obliquity^2.

    //  Functions Calculated for l = 2.

    c_IntMap<c_Key3, double> result_by_lmp(4);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(3);
    // Optimizations

    c_IntMap<c_Key1, double> result_by_p(2);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (2, 0).
    // p = 1
    tmp_double = -0.5;
    result_by_lmp.set(c_Key3(2, 0, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(2, 0), result_by_p);
    result_by_p.clear();

    // l , m = (2, 1).
    // p = 0
    tmp_double = 1.5*obliquity;
    result_by_lmp.set(c_Key3(2, 1, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -1.5*obliquity;
    result_by_lmp.set(c_Key3(2, 1, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(2, 1), result_by_p);
    result_by_p.clear();

    // l , m = (2, 2).
    // p = 0
    tmp_double = 3.0;
    result_by_lmp.set(c_Key3(2, 2, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(2, 2), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l2_4(double obliquity)
{
    // Inclination Functions Calculated for l = 2.
    // Functions are taylor expanded and truncated at obliquity^4.

    //  Functions Calculated for l = 2.

    c_IntMap<c_Key3, double> result_by_lmp(8);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(3);
    // Optimizations
    double obliquity_2 = obliquity * obliquity;
    double obliquity_3 = obliquity_2 * obliquity;

    c_IntMap<c_Key1, double> result_by_p(2);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (2, 0).
    // p = 0
    tmp_double = -0.375*obliquity_2;
    result_by_lmp.set(c_Key3(2, 0, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 0.75*obliquity_2 - 0.5;
    result_by_lmp.set(c_Key3(2, 0, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -0.375*obliquity_2;
    result_by_lmp.set(c_Key3(2, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(2, 0), result_by_p);
    result_by_p.clear();

    // l , m = (2, 1).
    // p = 0
    tmp_double = -0.625*obliquity_3 + 1.5*obliquity;
    result_by_lmp.set(c_Key3(2, 1, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = obliquity_3 - 1.5*obliquity;
    result_by_lmp.set(c_Key3(2, 1, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -0.375*obliquity_3;
    result_by_lmp.set(c_Key3(2, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(2, 1), result_by_p);
    result_by_p.clear();

    // l , m = (2, 2).
    // p = 0
    tmp_double = 3.0 - 1.5*obliquity_2;
    result_by_lmp.set(c_Key3(2, 2, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 1.5*obliquity_2;
    result_by_lmp.set(c_Key3(2, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(2, 2), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l2_off(double obliquity)
{
    // Inclination Functions Calculated for l = 2.
    // Functions assume `obliquity=0.0`. However some modes are still activated (non-zero `F(l,m,p)`).

    //  Functions Calculated for l = 2.

    c_IntMap<c_Key3, double> result_by_lmp(2);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(2);
    // Optimizations

    c_IntMap<c_Key1, double> result_by_p(2);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (2, 0).
    // p = 1
    tmp_double = -0.5;
    result_by_lmp.set(c_Key3(2, 0, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(2, 0), result_by_p);
    result_by_p.clear();

    // l , m = (2, 2).
    // p = 0
    tmp_double = 3.0;
    result_by_lmp.set(c_Key3(2, 2, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(2, 2), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}
