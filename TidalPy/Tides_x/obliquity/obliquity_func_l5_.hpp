#pragma once

#include "obliquity_common_.hpp"

ObliquityFuncOutput c_obliquity_function_l5_gen(double obliquity)
{
    // Inclination Functions Calculated for l = 5.
    // Functions are fully general with no assumptions on the magnitude of obliquity.

    //  Functions Calculated for l = 5.

    c_IntMap<c_Key3, double> result_by_lmp(36);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(6);
    // Optimizations
    double cos_i_half = std::cos((1.0/2.0)*obliquity);
    double sin_i_half = std::sin((1.0/2.0)*obliquity);
    double cos_i_half_2 = cos_i_half * cos_i_half;
    double cos_i_half_4 = cos_i_half_2 * cos_i_half_2;
    double cos_i_half_5 = cos_i_half_4 * cos_i_half;
    double sin_i_half_2 = sin_i_half * sin_i_half;
    double sin_i_half_4 = sin_i_half_2 * sin_i_half_2;
    double sin_i_half_5 = sin_i_half_4 * sin_i_half;
    double sin_i = std::sin(obliquity);
    double sin_i_2 = sin_i * sin_i;
    double sin_i_3 = sin_i_2 * sin_i;
    double cos_i = std::cos(obliquity);
    double cos_i_half_3 = cos_i_half_2 * cos_i_half;
    double cos_i_half_6 = cos_i_half_3 * cos_i_half_3;
    double cos_i_half_7 = cos_i_half_6 * cos_i_half;
    double sin_i_half_3 = sin_i_half_2 * sin_i_half;
    double sin_i_half_6 = sin_i_half_3 * sin_i_half_3;
    double sin_i_half_7 = sin_i_half_6 * sin_i_half;
    double sin_i_half_8 = sin_i_half_4 * sin_i_half_4;
    double sin_i_half_9 = sin_i_half_8 * sin_i_half;
    double cos_i_double = std::cos(2*obliquity);
    double cos_i_triple = std::cos(3*obliquity);
    double cos_i_half_8 = cos_i_half_4 * cos_i_half_4;
    double cos_i_half_10 = cos_i_half_5 * cos_i_half_5;
    double sin_i_half_10 = sin_i_half_5 * sin_i_half_5;
    double cos_i_2 = cos_i * cos_i;

    c_IntMap<c_Key1, double> result_by_p(5);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (5, 0).
    // p = 0
    tmp_double = 7.875*sin_i_half_5*cos_i_half_5;
    result_by_lmp.set(c_Key3(5, 0, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 0.13671875*(8.0 - 9.0*sin_i_2)*sin_i_3;
    result_by_lmp.set(c_Key3(5, 0, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 0.05859375*std::pow(cos_i + 1.0, 4)*sin_i + 1.875*sin_i_half_9*cos_i_half - 18.75*sin_i_half_7*cos_i_half_3 + 37.5*sin_i_half_5*cos_i_half_5 - 18.75*sin_i_half_3*cos_i_half_7;
    result_by_lmp.set(c_Key3(5, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -0.05859375*std::pow(cos_i + 1.0, 4)*sin_i - 1.875*sin_i_half_9*cos_i_half + 18.75*sin_i_half_7*cos_i_half_3 - 37.5*sin_i_half_5*cos_i_half_5 + 18.75*sin_i_half_3*cos_i_half_7;
    result_by_lmp.set(c_Key3(5, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 0.13671875*(9.0*sin_i_2 - 8.0)*sin_i_3;
    result_by_lmp.set(c_Key3(5, 0, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = -7.875*sin_i_half_5*cos_i_half_5;
    result_by_lmp.set(c_Key3(5, 0, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 0), result_by_p);
    result_by_p.clear();

    // l , m = (5, 1).
    // p = 0
    tmp_double = 39.375*sin_i_half_4*cos_i_half_6;
    result_by_lmp.set(c_Key3(5, 1, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 0.1025390625*std::pow(cos_i + 1.0, 2)*(-65.0*cos_i + 42.0*cos_i_double - 15.0*cos_i_triple + 38.0);
    result_by_lmp.set(c_Key3(5, 1, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -346.875*sin_i_half_10 + 834.375*sin_i_half_8 - 656.25*sin_i_half_6 + 168.75*sin_i_half_4 + 46.875*cos_i_half_10 - 45.0*cos_i_half_8;
    result_by_lmp.set(c_Key3(5, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 1.875*(115.0*sin_i_half_8 - 204.0*sin_i_half_6 + 90.0*sin_i_half_4 + 95.0*cos_i_half_8 - 80.0*cos_i_half_6)*sin_i_half_2;
    result_by_lmp.set(c_Key3(5, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 13.125*(-15.0*sin_i_half_6 + 33.0*sin_i_half_4 - 23.0*sin_i_half_2 + 5.0)*sin_i_half_4;
    result_by_lmp.set(c_Key3(5, 1, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 39.375*sin_i_half_6*cos_i_half_4;
    result_by_lmp.set(c_Key3(5, 1, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 1), result_by_p);
    result_by_p.clear();

    // l , m = (5, 2).
    // p = 0
    tmp_double = -157.5*sin_i_half_3*cos_i_half_7;
    result_by_lmp.set(c_Key3(5, 2, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -1.640625*std::pow(cos_i + 1.0, 4)*sin_i - 367.5*sin_i_half_5*cos_i_half_5 + 367.5*sin_i_half_3*cos_i_half_7;
    result_by_lmp.set(c_Key3(5, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 1.640625*std::pow(cos_i + 1.0, 4)*sin_i - 262.5*sin_i_half_7*cos_i_half_3 + 787.5*sin_i_half_5*cos_i_half_5 - 472.5*sin_i_half_3*cos_i_half_7;
    result_by_lmp.set(c_Key3(5, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 52.5*(-25.0*sin_i_half_6 + 39.0*sin_i_half_4 - 15.0*sin_i_half_2 + 5.0*cos_i_half_6)*sin_i_half_3*cos_i_half;
    result_by_lmp.set(c_Key3(5, 2, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 52.5*(15.0*sin_i_half_4 - 21.0*sin_i_half_2 + 7.0)*sin_i_half_5*cos_i_half;
    result_by_lmp.set(c_Key3(5, 2, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 157.5*sin_i_half_7*cos_i_half_3;
    result_by_lmp.set(c_Key3(5, 2, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 2), result_by_p);
    result_by_p.clear();

    // l , m = (5, 3).
    // p = 0
    tmp_double = -472.5*sin_i_half_2*cos_i_half_8;
    result_by_lmp.set(c_Key3(5, 3, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 1.640625*std::pow(cos_i + 1.0, 3)*(-45.0*cos_i_2 + 54.0*cos_i - 13.0);
    result_by_lmp.set(c_Key3(5, 3, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 2.4609375*std::pow(cos_i + 1.0, 2)*(65.0*cos_i - 42.0*cos_i_double + 15.0*cos_i_triple - 38.0);
    result_by_lmp.set(c_Key3(5, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = (4725.0*sin_i_half_6 - 10395.0*sin_i_half_4 + 7245.0*sin_i_half_2 - 1575.0)*sin_i_half_4;
    result_by_lmp.set(c_Key3(5, 3, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 13.125*(45.0*sin_i_2 - 54.0*cos_i - 58.0)*sin_i_half_6;
    result_by_lmp.set(c_Key3(5, 3, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = -472.5*sin_i_half_8*cos_i_half_2;
    result_by_lmp.set(c_Key3(5, 3, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 3), result_by_p);
    result_by_p.clear();

    // l , m = (5, 4).
    // p = 0
    tmp_double = 29.53125*std::pow(cos_i + 1.0, 4)*sin_i;
    result_by_lmp.set(c_Key3(5, 4, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -29.53125*std::pow(cos_i + 1.0, 4)*sin_i + 3780.0*sin_i_half_3*cos_i_half_7;
    result_by_lmp.set(c_Key3(5, 4, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = (945.0 - 4725.0*cos_i)*sin_i_half_3*cos_i_half_5;
    result_by_lmp.set(c_Key3(5, 4, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = (-4725.0*cos_i - 945.0)*sin_i_half_5*cos_i_half_3;
    result_by_lmp.set(c_Key3(5, 4, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 945.0*(5.0*sin_i_half_2 - 4.0)*sin_i_half_7*cos_i_half;
    result_by_lmp.set(c_Key3(5, 4, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = -945.0*sin_i_half_9*cos_i_half;
    result_by_lmp.set(c_Key3(5, 4, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 4), result_by_p);
    result_by_p.clear();

    // l , m = (5, 5).
    // p = 0
    tmp_double = 945.0*cos_i_half_10;
    result_by_lmp.set(c_Key3(5, 5, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 4725.0*sin_i_half_2*cos_i_half_8;
    result_by_lmp.set(c_Key3(5, 5, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 9450.0*sin_i_half_4*cos_i_half_6;
    result_by_lmp.set(c_Key3(5, 5, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 9450.0*sin_i_half_6*cos_i_half_4;
    result_by_lmp.set(c_Key3(5, 5, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 4725.0*sin_i_half_8*cos_i_half_2;
    result_by_lmp.set(c_Key3(5, 5, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 945.0*sin_i_half_10;
    result_by_lmp.set(c_Key3(5, 5, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 5), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l5_2(double obliquity)
{
    // Inclination Functions Calculated for l = 5.
    // Functions are taylor expanded and truncated at obliquity^2.

    //  Functions Calculated for l = 5.

    c_IntMap<c_Key3, double> result_by_lmp(9);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(6);
    // Optimizations

    c_IntMap<c_Key1, double> result_by_p(5);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (5, 0).
    // p = 2
    tmp_double = 0.9375*obliquity;
    result_by_lmp.set(c_Key3(5, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -0.9375*obliquity;
    result_by_lmp.set(c_Key3(5, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 0), result_by_p);
    result_by_p.clear();

    // l , m = (5, 1).
    // p = 2
    tmp_double = 1.875;
    result_by_lmp.set(c_Key3(5, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 1), result_by_p);
    result_by_p.clear();

    // l , m = (5, 2).
    // p = 1
    tmp_double = -26.25*obliquity;
    result_by_lmp.set(c_Key3(5, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 26.25*obliquity;
    result_by_lmp.set(c_Key3(5, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 2), result_by_p);
    result_by_p.clear();

    // l , m = (5, 3).
    // p = 1
    tmp_double = -52.5;
    result_by_lmp.set(c_Key3(5, 3, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 3), result_by_p);
    result_by_p.clear();

    // l , m = (5, 4).
    // p = 0
    tmp_double = 472.5*obliquity;
    result_by_lmp.set(c_Key3(5, 4, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -472.5*obliquity;
    result_by_lmp.set(c_Key3(5, 4, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 4), result_by_p);
    result_by_p.clear();

    // l , m = (5, 5).
    // p = 0
    tmp_double = 945.0;
    result_by_lmp.set(c_Key3(5, 5, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 5), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l5_4(double obliquity)
{
    // Inclination Functions Calculated for l = 5.
    // Functions are taylor expanded and truncated at obliquity^4.

    //  Functions Calculated for l = 5.

    c_IntMap<c_Key3, double> result_by_lmp(19);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(6);
    // Optimizations
    double obliquity_2 = obliquity * obliquity;
    double obliquity_3 = obliquity_2 * obliquity;

    c_IntMap<c_Key1, double> result_by_p(5);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (5, 0).
    // p = 1
    tmp_double = 1.09375*obliquity_3;
    result_by_lmp.set(c_Key3(5, 0, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -3.4375*obliquity_3 + 0.9375*obliquity;
    result_by_lmp.set(c_Key3(5, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 3.4375*obliquity_3 - 0.9375*obliquity;
    result_by_lmp.set(c_Key3(5, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -1.09375*obliquity_3;
    result_by_lmp.set(c_Key3(5, 0, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 0), result_by_p);
    result_by_p.clear();

    // l , m = (5, 1).
    // p = 1
    tmp_double = 6.5625*obliquity_2;
    result_by_lmp.set(c_Key3(5, 1, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 1.875 - 13.59375*obliquity_2;
    result_by_lmp.set(c_Key3(5, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 7.03125*obliquity_2;
    result_by_lmp.set(c_Key3(5, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 1), result_by_p);
    result_by_p.clear();

    // l , m = (5, 2).
    // p = 0
    tmp_double = -19.6875*obliquity_3;
    result_by_lmp.set(c_Key3(5, 2, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 76.5625*obliquity_3 - 26.25*obliquity;
    result_by_lmp.set(c_Key3(5, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -89.6875*obliquity_3 + 26.25*obliquity;
    result_by_lmp.set(c_Key3(5, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 32.8125*obliquity_3;
    result_by_lmp.set(c_Key3(5, 2, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 2), result_by_p);
    result_by_p.clear();

    // l , m = (5, 3).
    // p = 0
    tmp_double = -118.125*obliquity_2;
    result_by_lmp.set(c_Key3(5, 3, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 275.625*obliquity_2 - 52.5;
    result_by_lmp.set(c_Key3(5, 3, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -157.5*obliquity_2;
    result_by_lmp.set(c_Key3(5, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 3), result_by_p);
    result_by_p.clear();

    // l , m = (5, 4).
    // p = 0
    tmp_double = -551.25*obliquity_3 + 472.5*obliquity;
    result_by_lmp.set(c_Key3(5, 4, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 1023.75*obliquity_3 - 472.5*obliquity;
    result_by_lmp.set(c_Key3(5, 4, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -472.5*obliquity_3;
    result_by_lmp.set(c_Key3(5, 4, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 4), result_by_p);
    result_by_p.clear();

    // l , m = (5, 5).
    // p = 0
    tmp_double = 945.0 - 1181.25*obliquity_2;
    result_by_lmp.set(c_Key3(5, 5, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 1181.25*obliquity_2;
    result_by_lmp.set(c_Key3(5, 5, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 5), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l5_off(double obliquity)
{
    // Inclination Functions Calculated for l = 5.
    // Functions assume `obliquity=0.0`. However some modes are still activated (non-zero `F(l,m,p)`).

    //  Functions Calculated for l = 5.

    c_IntMap<c_Key3, double> result_by_lmp(3);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(3);
    // Optimizations

    c_IntMap<c_Key1, double> result_by_p(5);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (5, 1).
    // p = 2
    tmp_double = 1.875;
    result_by_lmp.set(c_Key3(5, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 1), result_by_p);
    result_by_p.clear();

    // l , m = (5, 3).
    // p = 1
    tmp_double = -52.5;
    result_by_lmp.set(c_Key3(5, 3, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 3), result_by_p);
    result_by_p.clear();

    // l , m = (5, 5).
    // p = 0
    tmp_double = 945.0;
    result_by_lmp.set(c_Key3(5, 5, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(5, 5), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}
