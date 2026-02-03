#pragma once

#include "obliquity_common_.hpp"

ObliquityFuncOutput c_obliquity_function_l4_gen(double obliquity)
{
    // Inclination Functions Calculated for l = 4.
    // Functions are fully general with no assumptions on the magnitude of obliquity.

    //  Functions Calculated for l = 4.

    c_IntMap<c_Key3, double> result_by_lmp(25);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(5);
    // Optimizations
    double cos_i_half = std::cos((1.0/2.0)*obliquity);
    double sin_i_half = std::sin((1.0/2.0)*obliquity);
    double cos_i_half_2 = cos_i_half * cos_i_half;
    double cos_i_half_4 = cos_i_half_2 * cos_i_half_2;
    double sin_i_half_2 = sin_i_half * sin_i_half;
    double sin_i_half_4 = sin_i_half_2 * sin_i_half_2;
    double sin_i_half_3 = sin_i_half_2 * sin_i_half;
    double sin_i_half_6 = sin_i_half_3 * sin_i_half_3;
    double cos_i_half_3 = cos_i_half_2 * cos_i_half;
    double cos_i_half_6 = cos_i_half_3 * cos_i_half_3;
    double cos_i_half_8 = cos_i_half_4 * cos_i_half_4;
    double sin_i_half_8 = sin_i_half_4 * sin_i_half_4;
    double cos_i_half_5 = cos_i_half_4 * cos_i_half;
    double cos_i = std::cos(obliquity);
    double sin_i = std::sin(obliquity);
    double sin_i_half_5 = sin_i_half_4 * sin_i_half;
    double sin_i_half_7 = sin_i_half_6 * sin_i_half;
    double sin_i_2 = sin_i * sin_i;
    double sin_i_3 = sin_i_2 * sin_i;
    double sin_i_double = std::sin(2*obliquity);

    c_IntMap<c_Key1, double> result_by_p(4);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (4, 0).
    // p = 0
    tmp_double = 4.375*sin_i_half_4*cos_i_half_4;
    result_by_lmp.set(c_Key3(4, 0, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 1.25*(-14.0*sin_i_half_6 + 28.0*sin_i_half_4 - 17.0*sin_i_half_2 + 3.0)*sin_i_half_2;
    result_by_lmp.set(c_Key3(4, 0, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 19.875*sin_i_half_8 - 33.0*sin_i_half_6 + 13.5*sin_i_half_4 + 6.375*cos_i_half_8 - 6.0*cos_i_half_6;
    result_by_lmp.set(c_Key3(4, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 1.25*(-14.0*sin_i_half_6 + 28.0*sin_i_half_4 - 17.0*sin_i_half_2 + 3.0)*sin_i_half_2;
    result_by_lmp.set(c_Key3(4, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 4.375*sin_i_half_4*cos_i_half_4;
    result_by_lmp.set(c_Key3(4, 0, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(4, 0), result_by_p);
    result_by_p.clear();

    // l , m = (4, 1).
    // p = 0
    tmp_double = -17.5*sin_i_half_3*cos_i_half_5;
    result_by_lmp.set(c_Key3(4, 1, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -0.46875*std::pow(cos_i + 1.0, 3)*sin_i - 25.0*sin_i_half_5*cos_i_half_3 + 37.5*sin_i_half_3*cos_i_half_5;
    result_by_lmp.set(c_Key3(4, 1, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 0.46875*std::pow(cos_i + 1.0, 3)*sin_i - 7.5*sin_i_half_7*cos_i_half + 45.0*sin_i_half_5*cos_i_half_3 - 45.0*sin_i_half_3*cos_i_half_5;
    result_by_lmp.set(c_Key3(4, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 2.5*(-7.0*sin_i_2 + 3.5*cos_i + 6.5)*sin_i_half_3*cos_i_half;
    result_by_lmp.set(c_Key3(4, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 17.5*sin_i_half_5*cos_i_half_3;
    result_by_lmp.set(c_Key3(4, 1, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(4, 1), result_by_p);
    result_by_p.clear();

    // l , m = (4, 2).
    // p = 0
    tmp_double = -52.5*sin_i_half_2*cos_i_half_6;
    result_by_lmp.set(c_Key3(4, 2, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 1.875*std::pow(cos_i + 1.0, 2)*(7.0*sin_i_2 + 7.0*cos_i - 8.0);
    result_by_lmp.set(c_Key3(4, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 22.5*(14.0*sin_i_half_6 - 28.0*sin_i_half_4 + 17.0*sin_i_half_2 - 3.0)*sin_i_half_2;
    result_by_lmp.set(c_Key3(4, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = (52.5*sin_i_2 - 52.5*cos_i - 60.0)*sin_i_half_4;
    result_by_lmp.set(c_Key3(4, 2, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -52.5*sin_i_half_6*cos_i_half_2;
    result_by_lmp.set(c_Key3(4, 2, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(4, 2), result_by_p);
    result_by_p.clear();

    // l , m = (4, 3).
    // p = 0
    tmp_double = 6.5625*std::pow(cos_i + 1.0, 3)*sin_i;
    result_by_lmp.set(c_Key3(4, 3, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -6.5625*std::pow(cos_i + 1.0, 3)*sin_i + 315.0*sin_i_half_3*cos_i_half_5;
    result_by_lmp.set(c_Key3(4, 3, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -39.375*sin_i_3*cos_i;
    result_by_lmp.set(c_Key3(4, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -13.125*(sin_i + sin_i_double)*std::pow(cos_i - 1.0, 2);
    result_by_lmp.set(c_Key3(4, 3, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -105.0*sin_i_half_7*cos_i_half;
    result_by_lmp.set(c_Key3(4, 3, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(4, 3), result_by_p);
    result_by_p.clear();

    // l , m = (4, 4).
    // p = 0
    tmp_double = 105.0*cos_i_half_8;
    result_by_lmp.set(c_Key3(4, 4, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 420.0*sin_i_half_2*cos_i_half_6;
    result_by_lmp.set(c_Key3(4, 4, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 630.0*sin_i_half_4*cos_i_half_4;
    result_by_lmp.set(c_Key3(4, 4, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 420.0*sin_i_half_6*cos_i_half_2;
    result_by_lmp.set(c_Key3(4, 4, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 105.0*sin_i_half_8;
    result_by_lmp.set(c_Key3(4, 4, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(4, 4), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l4_2(double obliquity)
{
    // Inclination Functions Calculated for l = 4.
    // Functions are taylor expanded and truncated at obliquity^2.

    //  Functions Calculated for l = 4.

    c_IntMap<c_Key3, double> result_by_lmp(7);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(5);
    // Optimizations

    c_IntMap<c_Key1, double> result_by_p(4);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (4, 0).
    // p = 2
    tmp_double = 0.375;
    result_by_lmp.set(c_Key3(4, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(4, 0), result_by_p);
    result_by_p.clear();

    // l , m = (4, 1).
    // p = 1
    tmp_double = -3.75*obliquity;
    result_by_lmp.set(c_Key3(4, 1, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 3.75*obliquity;
    result_by_lmp.set(c_Key3(4, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(4, 1), result_by_p);
    result_by_p.clear();

    // l , m = (4, 2).
    // p = 1
    tmp_double = -7.5;
    result_by_lmp.set(c_Key3(4, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(4, 2), result_by_p);
    result_by_p.clear();

    // l , m = (4, 3).
    // p = 0
    tmp_double = 52.5*obliquity;
    result_by_lmp.set(c_Key3(4, 3, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -52.5*obliquity;
    result_by_lmp.set(c_Key3(4, 3, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(4, 3), result_by_p);
    result_by_p.clear();

    // l , m = (4, 4).
    // p = 0
    tmp_double = 105.0;
    result_by_lmp.set(c_Key3(4, 4, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(4, 4), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l4_4(double obliquity)
{
    // Inclination Functions Calculated for l = 4.
    // Functions are taylor expanded and truncated at obliquity^4.

    //  Functions Calculated for l = 4.

    c_IntMap<c_Key3, double> result_by_lmp(15);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(5);
    // Optimizations
    double obliquity_2 = obliquity * obliquity;
    double obliquity_3 = obliquity_2 * obliquity;

    c_IntMap<c_Key1, double> result_by_p(4);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (4, 0).
    // p = 1
    tmp_double = 0.9375*obliquity_2;
    result_by_lmp.set(c_Key3(4, 0, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 0.375 - 1.875*obliquity_2;
    result_by_lmp.set(c_Key3(4, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 0.9375*obliquity_2;
    result_by_lmp.set(c_Key3(4, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(4, 0), result_by_p);
    result_by_p.clear();

    // l , m = (4, 1).
    // p = 0
    tmp_double = -2.1875*obliquity_3;
    result_by_lmp.set(c_Key3(4, 1, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 8.125*obliquity_3 - 3.75*obliquity;
    result_by_lmp.set(c_Key3(4, 1, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -9.0625*obliquity_3 + 3.75*obliquity;
    result_by_lmp.set(c_Key3(4, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 3.125*obliquity_3;
    result_by_lmp.set(c_Key3(4, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(4, 1), result_by_p);
    result_by_p.clear();

    // l , m = (4, 2).
    // p = 0
    tmp_double = -13.125*obliquity_2;
    result_by_lmp.set(c_Key3(4, 2, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 30.0*obliquity_2 - 7.5;
    result_by_lmp.set(c_Key3(4, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -16.875*obliquity_2;
    result_by_lmp.set(c_Key3(4, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(4, 2), result_by_p);
    result_by_p.clear();

    // l , m = (4, 3).
    // p = 0
    tmp_double = -48.125*obliquity_3 + 52.5*obliquity;
    result_by_lmp.set(c_Key3(4, 3, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 87.5*obliquity_3 - 52.5*obliquity;
    result_by_lmp.set(c_Key3(4, 3, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -39.375*obliquity_3;
    result_by_lmp.set(c_Key3(4, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(4, 3), result_by_p);
    result_by_p.clear();

    // l , m = (4, 4).
    // p = 0
    tmp_double = 105.0 - 105.0*obliquity_2;
    result_by_lmp.set(c_Key3(4, 4, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 105.0*obliquity_2;
    result_by_lmp.set(c_Key3(4, 4, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(4, 4), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l4_off(double obliquity)
{
    // Inclination Functions Calculated for l = 4.
    // Functions assume `obliquity=0.0`. However some modes are still activated (non-zero `F(l,m,p)`).

    //  Functions Calculated for l = 4.

    c_IntMap<c_Key3, double> result_by_lmp(3);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(3);
    // Optimizations

    c_IntMap<c_Key1, double> result_by_p(4);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (4, 0).
    // p = 2
    tmp_double = 0.375;
    result_by_lmp.set(c_Key3(4, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(4, 0), result_by_p);
    result_by_p.clear();

    // l , m = (4, 2).
    // p = 1
    tmp_double = -7.5;
    result_by_lmp.set(c_Key3(4, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(4, 2), result_by_p);
    result_by_p.clear();

    // l , m = (4, 4).
    // p = 0
    tmp_double = 105.0;
    result_by_lmp.set(c_Key3(4, 4, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(4, 4), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}
