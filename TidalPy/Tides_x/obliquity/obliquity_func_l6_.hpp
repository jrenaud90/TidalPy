#pragma once

#include "obliquity_common_.hpp"

ObliquityFuncOutput c_obliquity_function_l6_gen(double obliquity)
{
    // Inclination Functions Calculated for l = 6.
    // Functions are fully general with no assumptions on the magnitude of obliquity.

    //  Functions Calculated for l = 6.

    c_IntMap<c_Key3, double> result_by_lmp(49);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(7);
    // Optimizations
    double cos_i_half = std::cos((1.0/2.0)*obliquity);
    double sin_i_half = std::sin((1.0/2.0)*obliquity);
    double cos_i_half_2 = cos_i_half * cos_i_half;
    double cos_i_half_3 = cos_i_half_2 * cos_i_half;
    double cos_i_half_6 = cos_i_half_3 * cos_i_half_3;
    double sin_i_half_2 = sin_i_half * sin_i_half;
    double sin_i_half_3 = sin_i_half_2 * sin_i_half;
    double sin_i_half_6 = sin_i_half_3 * sin_i_half_3;
    double sin_i = std::sin(obliquity);
    double sin_i_2 = sin_i * sin_i;
    double sin_i_4 = sin_i_2 * sin_i_2;
    double cos_i_half_4 = cos_i_half_2 * cos_i_half_2;
    double cos_i_half_8 = cos_i_half_4 * cos_i_half_4;
    double cos_i_half_5 = cos_i_half_4 * cos_i_half;
    double cos_i_half_10 = cos_i_half_5 * cos_i_half_5;
    double sin_i_half_4 = sin_i_half_2 * sin_i_half_2;
    double sin_i_half_8 = sin_i_half_4 * sin_i_half_4;
    double sin_i_half_5 = sin_i_half_4 * sin_i_half;
    double sin_i_half_10 = sin_i_half_5 * sin_i_half_5;
    double cos_i_half_12 = cos_i_half_6 * cos_i_half_6;
    double sin_i_half_12 = sin_i_half_6 * sin_i_half_6;
    double cos_i_half_7 = cos_i_half_6 * cos_i_half;
    double cos_i = std::cos(obliquity);
    double cos_i_half_9 = cos_i_half_8 * cos_i_half;
    double sin_i_half_7 = sin_i_half_6 * sin_i_half;
    double sin_i_half_9 = sin_i_half_8 * sin_i_half;
    double sin_i_half_11 = sin_i_half_10 * sin_i_half;
    double cos_i_double = std::cos(2*obliquity);
    double cos_i_triple = std::cos(3*obliquity);

    c_IntMap<c_Key1, double> result_by_p(6);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (6, 0).
    // p = 0
    tmp_double = -14.4375*sin_i_half_6*cos_i_half_6;
    result_by_lmp.set(c_Key3(6, 0, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 0.123046875*(11.0*sin_i_2 - 10.0)*sin_i_4;
    result_by_lmp.set(c_Key3(6, 0, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 6.5625*(24.0*sin_i_half_10 - 62.0*sin_i_half_8 + 53.0*sin_i_half_6 - 15.0*sin_i_half_4 - 9.0*cos_i_half_10 + 8.0*cos_i_half_8)*sin_i_half_2;
    result_by_lmp.set(c_Key3(6, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -81.875*sin_i_half_12 + 151.875*sin_i_half_10 - 70.3125*sin_i_half_8 + 125.0*sin_i_half_6*cos_i_half_6 - 70.3125*sin_i_half_4*cos_i_half_8 - 11.5625*cos_i_half_12 + 11.25*cos_i_half_10;
    result_by_lmp.set(c_Key3(6, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 6.5625*(24.0*sin_i_half_10 - 62.0*sin_i_half_8 + 53.0*sin_i_half_6 - 15.0*sin_i_half_4 - 9.0*cos_i_half_10 + 8.0*cos_i_half_8)*sin_i_half_2;
    result_by_lmp.set(c_Key3(6, 0, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 0.123046875*(11.0*sin_i_2 - 10.0)*sin_i_4;
    result_by_lmp.set(c_Key3(6, 0, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = -14.4375*sin_i_half_6*cos_i_half_6;
    result_by_lmp.set(c_Key3(6, 0, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 0), result_by_p);
    result_by_p.clear();

    // l , m = (6, 1).
    // p = 0
    tmp_double = 86.625*sin_i_half_5*cos_i_half_7;
    result_by_lmp.set(c_Key3(6, 1, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 7.875*(66.0*sin_i_half_4 - 55.0*sin_i_half_2 + 10.0)*sin_i_half_3*cos_i_half_5;
    result_by_lmp.set(c_Key3(6, 1, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 0.205078125*std::pow(cos_i + 1.0, 5)*sin_i + 91.875*sin_i_half_9*cos_i_half_3 - 459.375*sin_i_half_7*cos_i_half_5 + 551.25*sin_i_half_5*cos_i_half_7 - 183.75*sin_i_half_3*cos_i_half_9;
    result_by_lmp.set(c_Key3(6, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -0.205078125*std::pow(cos_i + 1.0, 5)*sin_i + 13.125*sin_i_half_11*cos_i_half - 196.875*sin_i_half_9*cos_i_half_3 + 656.25*sin_i_half_7*cos_i_half_5 - 656.25*sin_i_half_5*cos_i_half_7 + 196.875*sin_i_half_3*cos_i_half_9;
    result_by_lmp.set(c_Key3(6, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 13.125*(-57.0*sin_i_half_8 + 98.0*sin_i_half_6 - 42.0*sin_i_half_4 - 42.0*cos_i_half_8 + 35.0*cos_i_half_6)*sin_i_half_3*cos_i_half;
    result_by_lmp.set(c_Key3(6, 1, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 7.875*(-66.0*sin_i_half_4 + 77.0*sin_i_half_2 - 21.0)*sin_i_half_5*cos_i_half_3;
    result_by_lmp.set(c_Key3(6, 1, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = -86.625*sin_i_half_7*cos_i_half_5;
    result_by_lmp.set(c_Key3(6, 1, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 1), result_by_p);
    result_by_p.clear();

    // l , m = (6, 2).
    // p = 0
    tmp_double = 433.125*sin_i_half_4*cos_i_half_8;
    result_by_lmp.set(c_Key3(6, 2, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 0.3076171875*std::pow(cos_i + 1.0, 3)*(-191.0*cos_i + 110.0*cos_i_double - 33.0*cos_i_triple + 114.0);
    result_by_lmp.set(c_Key3(6, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 13.125*(42.0*std::pow(cos_i - 1.0, 2) + 462.0*sin_i_half_8 - 560.0*sin_i_half_6 + 33.0*cos_i_half_8 - 32.0*cos_i_half_6)*cos_i_half_4;
    result_by_lmp.set(c_Key3(6, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 262.5*(-24.0*sin_i_half_10 + 62.0*sin_i_half_8 - 53.0*sin_i_half_6 + 15.0*sin_i_half_4 + 9.0*cos_i_half_10 - 8.0*cos_i_half_8)*sin_i_half_2;
    result_by_lmp.set(c_Key3(6, 2, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 13.125*(201.0*sin_i_half_8 - 368.0*sin_i_half_6 + 168.0*sin_i_half_4 + 294.0*cos_i_half_8 - 224.0*cos_i_half_6)*sin_i_half_4;
    result_by_lmp.set(c_Key3(6, 2, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 78.75*(-33.0*sin_i_half_6 + 77.0*sin_i_half_4 - 58.0*sin_i_half_2 + 14.0)*sin_i_half_6;
    result_by_lmp.set(c_Key3(6, 2, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 433.125*sin_i_half_8*cos_i_half_4;
    result_by_lmp.set(c_Key3(6, 2, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 2), result_by_p);
    result_by_p.clear();

    // l , m = (6, 3).
    // p = 0
    tmp_double = -1732.5*sin_i_half_3*cos_i_half_9;
    result_by_lmp.set(c_Key3(6, 3, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -7.3828125*std::pow(cos_i + 1.0, 5)*sin_i - 5670.0*sin_i_half_5*cos_i_half_7 + 4252.5*sin_i_half_3*cos_i_half_9;
    result_by_lmp.set(c_Key3(6, 3, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 7.3828125*std::pow(cos_i + 1.0, 5)*sin_i - 6615.0*sin_i_half_7*cos_i_half_5 + 13230.0*sin_i_half_5*cos_i_half_7 - 5670.0*sin_i_half_3*cos_i_half_9;
    result_by_lmp.set(c_Key3(6, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 3150.0*(11.0*sin_i_half_4 - 11.0*sin_i_half_2 + 2.0)*sin_i_half_3*std::sin((1.0/2.0)*obliquity + M_PI_4)*cos_i_half_3*std::cos((1.0/2.0)*obliquity + M_PI_4);
    result_by_lmp.set(c_Key3(6, 3, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 472.5*(-41.0*sin_i_half_6 + 68.0*sin_i_half_4 - 28.0*sin_i_half_2 + 14.0*cos_i_half_6)*sin_i_half_5*cos_i_half;
    result_by_lmp.set(c_Key3(6, 3, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 472.5*(22.0*sin_i_half_4 - 33.0*sin_i_half_2 + 12.0)*sin_i_half_7*cos_i_half;
    result_by_lmp.set(c_Key3(6, 3, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 1732.5*sin_i_half_9*cos_i_half_3;
    result_by_lmp.set(c_Key3(6, 3, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 3), result_by_p);
    result_by_p.clear();

    // l , m = (6, 4).
    // p = 0
    tmp_double = -5197.5*sin_i_half_2*cos_i_half_10;
    result_by_lmp.set(c_Key3(6, 4, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 14.765625*std::pow(cos_i + 1.0, 4)*(33.0*sin_i_2 + 44.0*cos_i - 46.0);
    result_by_lmp.set(c_Key3(6, 4, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 9.228515625*std::pow(cos_i + 1.0, 3)*(191.0*cos_i - 110.0*cos_i_double + 33.0*cos_i_triple - 114.0);
    result_by_lmp.set(c_Key3(6, 4, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = (25987.5*sin_i_2 - 23625.0)*sin_i_half_4*cos_i_half_4;
    result_by_lmp.set(c_Key3(6, 4, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 2362.5*(33.0*sin_i_half_6 - 77.0*sin_i_half_4 + 58.0*sin_i_half_2 - 14.0)*sin_i_half_6;
    result_by_lmp.set(c_Key3(6, 4, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = (7796.25*sin_i_2 - 10395.0*cos_i - 10867.5)*sin_i_half_8;
    result_by_lmp.set(c_Key3(6, 4, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = -5197.5*sin_i_half_10*cos_i_half_2;
    result_by_lmp.set(c_Key3(6, 4, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 4), result_by_p);
    result_by_p.clear();

    // l , m = (6, 5).
    // p = 0
    tmp_double = 162.421875*std::pow(cos_i + 1.0, 5)*sin_i;
    result_by_lmp.set(c_Key3(6, 5, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -162.421875*std::pow(cos_i + 1.0, 5)*sin_i + 51975.0*sin_i_half_3*cos_i_half_9;
    result_by_lmp.set(c_Key3(6, 5, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 25987.5*(1.0 - 3.0*cos_i)*sin_i_half_3*cos_i_half_7;
    result_by_lmp.set(c_Key3(6, 5, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -103950.0*sin_i_half_5*cos_i_half_5*cos_i;
    result_by_lmp.set(c_Key3(6, 5, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = (51975.0 - 155925.0*cos_i_half_2)*sin_i_half_7*cos_i_half_3;
    result_by_lmp.set(c_Key3(6, 5, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = -324.84375*std::pow(cos_i - 1.0, 4)*(3.0*cos_i + 2.0)*sin_i;
    result_by_lmp.set(c_Key3(6, 5, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = -10395.0*sin_i_half_11*cos_i_half;
    result_by_lmp.set(c_Key3(6, 5, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 5), result_by_p);
    result_by_p.clear();

    // l , m = (6, 6).
    // p = 0
    tmp_double = 10395.0*cos_i_half_12;
    result_by_lmp.set(c_Key3(6, 6, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 62370.0*sin_i_half_2*cos_i_half_10;
    result_by_lmp.set(c_Key3(6, 6, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 155925.0*sin_i_half_4*cos_i_half_8;
    result_by_lmp.set(c_Key3(6, 6, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 207900.0*sin_i_half_6*cos_i_half_6;
    result_by_lmp.set(c_Key3(6, 6, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 155925.0*sin_i_half_8*cos_i_half_4;
    result_by_lmp.set(c_Key3(6, 6, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 62370.0*sin_i_half_10*cos_i_half_2;
    result_by_lmp.set(c_Key3(6, 6, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 10395.0*sin_i_half_12;
    result_by_lmp.set(c_Key3(6, 6, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 6), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l6_2(double obliquity)
{
    // Inclination Functions Calculated for l = 6.
    // Functions are taylor expanded and truncated at obliquity^2.

    //  Functions Calculated for l = 6.

    c_IntMap<c_Key3, double> result_by_lmp(10);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(7);
    // Optimizations

    c_IntMap<c_Key1, double> result_by_p(6);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (6, 0).
    // p = 3
    tmp_double = -0.3125;
    result_by_lmp.set(c_Key3(6, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 0), result_by_p);
    result_by_p.clear();

    // l , m = (6, 1).
    // p = 2
    tmp_double = 6.5625*obliquity;
    result_by_lmp.set(c_Key3(6, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -6.5625*obliquity;
    result_by_lmp.set(c_Key3(6, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 1), result_by_p);
    result_by_p.clear();

    // l , m = (6, 2).
    // p = 2
    tmp_double = 13.125;
    result_by_lmp.set(c_Key3(6, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 2), result_by_p);
    result_by_p.clear();

    // l , m = (6, 3).
    // p = 1
    tmp_double = -236.25*obliquity;
    result_by_lmp.set(c_Key3(6, 3, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 236.25*obliquity;
    result_by_lmp.set(c_Key3(6, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 3), result_by_p);
    result_by_p.clear();

    // l , m = (6, 4).
    // p = 1
    tmp_double = -472.5;
    result_by_lmp.set(c_Key3(6, 4, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 4), result_by_p);
    result_by_p.clear();

    // l , m = (6, 5).
    // p = 0
    tmp_double = 5197.5*obliquity;
    result_by_lmp.set(c_Key3(6, 5, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -5197.5*obliquity;
    result_by_lmp.set(c_Key3(6, 5, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 5), result_by_p);
    result_by_p.clear();

    // l , m = (6, 6).
    // p = 0
    tmp_double = 10395.0;
    result_by_lmp.set(c_Key3(6, 6, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 6), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l6_4(double obliquity)
{
    // Inclination Functions Calculated for l = 6.
    // Functions are taylor expanded and truncated at obliquity^4.

    //  Functions Calculated for l = 6.

    c_IntMap<c_Key3, double> result_by_lmp(22);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(7);
    // Optimizations
    double obliquity_2 = obliquity * obliquity;
    double obliquity_3 = obliquity_2 * obliquity;

    c_IntMap<c_Key1, double> result_by_p(6);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (6, 0).
    // p = 2
    tmp_double = -1.640625*obliquity_2;
    result_by_lmp.set(c_Key3(6, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 3.28125*obliquity_2 - 0.3125;
    result_by_lmp.set(c_Key3(6, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -1.640625*obliquity_2;
    result_by_lmp.set(c_Key3(6, 0, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 0), result_by_p);
    result_by_p.clear();

    // l , m = (6, 1).
    // p = 1
    tmp_double = 9.84375*obliquity_3;
    result_by_lmp.set(c_Key3(6, 1, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -32.265625*obliquity_3 + 6.5625*obliquity;
    result_by_lmp.set(c_Key3(6, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 33.90625*obliquity_3 - 6.5625*obliquity;
    result_by_lmp.set(c_Key3(6, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -11.484375*obliquity_3;
    result_by_lmp.set(c_Key3(6, 1, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 1), result_by_p);
    result_by_p.clear();

    // l , m = (6, 2).
    // p = 1
    tmp_double = 59.0625*obliquity_2;
    result_by_lmp.set(c_Key3(6, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 13.125 - 124.6875*obliquity_2;
    result_by_lmp.set(c_Key3(6, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 65.625*obliquity_2;
    result_by_lmp.set(c_Key3(6, 2, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 2), result_by_p);
    result_by_p.clear();

    // l , m = (6, 3).
    // p = 0
    tmp_double = -216.5625*obliquity_3;
    result_by_lmp.set(c_Key3(6, 3, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 866.25*obliquity_3 - 236.25*obliquity;
    result_by_lmp.set(c_Key3(6, 3, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -1043.4375*obliquity_3 + 236.25*obliquity;
    result_by_lmp.set(c_Key3(6, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 393.75*obliquity_3;
    result_by_lmp.set(c_Key3(6, 3, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 3), result_by_p);
    result_by_p.clear();

    // l , m = (6, 4).
    // p = 0
    tmp_double = -1299.375*obliquity_2;
    result_by_lmp.set(c_Key3(6, 4, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 3071.25*obliquity_2 - 472.5;
    result_by_lmp.set(c_Key3(6, 4, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -1771.875*obliquity_2;
    result_by_lmp.set(c_Key3(6, 4, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 4), result_by_p);
    result_by_p.clear();

    // l , m = (6, 5).
    // p = 0
    tmp_double = -7363.125*obliquity_3 + 5197.5*obliquity;
    result_by_lmp.set(c_Key3(6, 5, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 13860.0*obliquity_3 - 5197.5*obliquity;
    result_by_lmp.set(c_Key3(6, 5, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -6496.875*obliquity_3;
    result_by_lmp.set(c_Key3(6, 5, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 5), result_by_p);
    result_by_p.clear();

    // l , m = (6, 6).
    // p = 0
    tmp_double = 10395.0 - 15592.5*obliquity_2;
    result_by_lmp.set(c_Key3(6, 6, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 15592.5*obliquity_2;
    result_by_lmp.set(c_Key3(6, 6, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 6), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l6_off(double obliquity)
{
    // Inclination Functions Calculated for l = 6.
    // Functions assume `obliquity=0.0`. However some modes are still activated (non-zero `F(l,m,p)`).

    //  Functions Calculated for l = 6.

    c_IntMap<c_Key3, double> result_by_lmp(4);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(4);
    // Optimizations

    c_IntMap<c_Key1, double> result_by_p(6);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (6, 0).
    // p = 3
    tmp_double = -0.3125;
    result_by_lmp.set(c_Key3(6, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 0), result_by_p);
    result_by_p.clear();

    // l , m = (6, 2).
    // p = 2
    tmp_double = 13.125;
    result_by_lmp.set(c_Key3(6, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 2), result_by_p);
    result_by_p.clear();

    // l , m = (6, 4).
    // p = 1
    tmp_double = -472.5;
    result_by_lmp.set(c_Key3(6, 4, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 4), result_by_p);
    result_by_p.clear();

    // l , m = (6, 6).
    // p = 0
    tmp_double = 10395.0;
    result_by_lmp.set(c_Key3(6, 6, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(6, 6), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}
