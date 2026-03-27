#pragma once

#include "obliquity_common_.hpp"

ObliquityFuncOutput c_obliquity_function_l7_gen(double obliquity)
{
    // Inclination Functions Calculated for l = 7.
    // Functions are fully general with no assumptions on the magnitude of obliquity.

    //  Functions Calculated for l = 7.

    c_IntMap<c_Key3, double> result_by_lmp(64);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(8);
    // Optimizations
    double cos_i_half = std::cos((1.0/2.0)*obliquity);
    double sin_i_half = std::sin((1.0/2.0)*obliquity);
    double cos_i_half_2 = cos_i_half * cos_i_half;
    double cos_i_half_3 = cos_i_half_2 * cos_i_half;
    double cos_i_half_6 = cos_i_half_3 * cos_i_half_3;
    double cos_i_half_7 = cos_i_half_6 * cos_i_half;
    double sin_i_half_2 = sin_i_half * sin_i_half;
    double sin_i_half_3 = sin_i_half_2 * sin_i_half;
    double sin_i_half_6 = sin_i_half_3 * sin_i_half_3;
    double sin_i_half_7 = sin_i_half_6 * sin_i_half;
    double cos_i_half_4 = cos_i_half_2 * cos_i_half_2;
    double cos_i_half_5 = cos_i_half_4 * cos_i_half;
    double sin_i_half_4 = sin_i_half_2 * sin_i_half_2;
    double sin_i_half_5 = sin_i_half_4 * sin_i_half;
    double cos_i_half_8 = cos_i_half_4 * cos_i_half_4;
    double sin_i_half_8 = sin_i_half_4 * sin_i_half_4;
    double cos_i = std::cos(obliquity);
    double sin_i = std::sin(obliquity);
    double cos_i_half_9 = cos_i_half_8 * cos_i_half;
    double cos_i_half_10 = cos_i_half_5 * cos_i_half_5;
    double cos_i_half_11 = cos_i_half_10 * cos_i_half;
    double sin_i_half_9 = sin_i_half_8 * sin_i_half;
    double sin_i_half_10 = sin_i_half_5 * sin_i_half_5;
    double sin_i_half_11 = sin_i_half_10 * sin_i_half;
    double sin_i_half_12 = sin_i_half_6 * sin_i_half_6;
    double sin_i_half_13 = sin_i_half_12 * sin_i_half;
    double sin_i_2 = sin_i * sin_i;
    double cos_i_half_12 = cos_i_half_6 * cos_i_half_6;
    double cos_i_half_14 = cos_i_half_7 * cos_i_half_7;
    double sin_i_half_14 = sin_i_half_7 * sin_i_half_7;
    double cos_i_double = std::cos(2*obliquity);
    double cos_i_triple = std::cos(3*obliquity);

    c_IntMap<c_Key1, double> result_by_p(7);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (7, 0).
    // p = 0
    tmp_double = -26.8125*sin_i_half_7*cos_i_half_7;
    result_by_lmp.set(c_Key3(7, 0, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 14.4375*(-13.0*sin_i_half_4 + 13.0*sin_i_half_2 - 3.0)*sin_i_half_5*cos_i_half_5;
    result_by_lmp.set(c_Key3(7, 0, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 3.9375*(-103.0*sin_i_half_8 + 161.0*sin_i_half_6 - 63.0*sin_i_half_4 - 40.0*cos_i_half_8 + 35.0*cos_i_half_6)*sin_i_half_3*cos_i_half_3;
    result_by_lmp.set(c_Key3(7, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -0.01708984375*std::pow(cos_i + 1.0, 6)*sin_i - 2.1875*sin_i_half_13*cos_i_half + 45.9375*sin_i_half_11*cos_i_half_3 - 229.6875*sin_i_half_9*cos_i_half_5 + 382.8125*sin_i_half_7*cos_i_half_7 - 229.6875*sin_i_half_5*cos_i_half_9 + 45.9375*sin_i_half_3*cos_i_half_11;
    result_by_lmp.set(c_Key3(7, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 0.01708984375*std::pow(cos_i + 1.0, 6)*sin_i + 2.1875*sin_i_half_13*cos_i_half - 45.9375*sin_i_half_11*cos_i_half_3 + 229.6875*sin_i_half_9*cos_i_half_5 - 382.8125*sin_i_half_7*cos_i_half_7 + 229.6875*sin_i_half_5*cos_i_half_9 - 45.9375*sin_i_half_3*cos_i_half_11;
    result_by_lmp.set(c_Key3(7, 0, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 3.9375*(103.0*sin_i_half_8 - 161.0*sin_i_half_6 + 63.0*sin_i_half_4 + 40.0*cos_i_half_8 - 35.0*cos_i_half_6)*sin_i_half_3*cos_i_half_3;
    result_by_lmp.set(c_Key3(7, 0, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 14.4375*(13.0*sin_i_half_4 - 13.0*sin_i_half_2 + 3.0)*sin_i_half_5*cos_i_half_5;
    result_by_lmp.set(c_Key3(7, 0, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 26.8125*sin_i_half_7*cos_i_half_7;
    result_by_lmp.set(c_Key3(7, 0, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 0), result_by_p);
    result_by_p.clear();

    // l , m = (7, 1).
    // p = 0
    tmp_double = -187.6875*sin_i_half_6*cos_i_half_8;
    result_by_lmp.set(c_Key3(7, 1, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 0.11279296875*std::pow(cos_i - 1.0, 2)*std::pow(cos_i + 1.0, 3)*(91.0*sin_i_2 + 26.0*cos_i - 86.0);
    result_by_lmp.set(c_Key3(7, 1, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 3.9375*(-826.0*sin_i_half_8 + 1176.0*sin_i_half_6 - 420.0*sin_i_half_4 - 175.0*cos_i_half_8 + 160.0*cos_i_half_6)*sin_i_half_2*cos_i_half_4;
    result_by_lmp.set(c_Key3(7, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 3093.125*sin_i_half_14 - 8421.875*sin_i_half_12 + 7625.625*sin_i_half_10 - 2296.875*sin_i_half_8 + 2450.0*sin_i_half_6*cos_i_half_8 - 918.75*sin_i_half_4*cos_i_half_10 - 107.1875*cos_i_half_14 + 105.0*cos_i_half_12;
    result_by_lmp.set(c_Key3(7, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 2.1875*(-469.0*sin_i_half_12 + 888.0*sin_i_half_10 - 420.0*sin_i_half_8 + 1120.0*sin_i_half_6*cos_i_half_6 - 1050.0*sin_i_half_4*cos_i_half_8 - 364.0*cos_i_half_12 + 336.0*cos_i_half_10)*sin_i_half_2;
    result_by_lmp.set(c_Key3(7, 1, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 3.9375*(595.0*sin_i_half_10 - 1595.0*sin_i_half_8 + 1420.0*sin_i_half_6 - 420.0*sin_i_half_4 - 406.0*cos_i_half_10 + 336.0*cos_i_half_8)*sin_i_half_4;
    result_by_lmp.set(c_Key3(7, 1, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 0.11279296875*std::pow(cos_i - 1.0, 3)*std::pow(cos_i + 1.0, 2)*(-91.0*sin_i_2 + 26.0*cos_i + 86.0);
    result_by_lmp.set(c_Key3(7, 1, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = -187.6875*sin_i_half_8*cos_i_half_6;
    result_by_lmp.set(c_Key3(7, 1, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 1), result_by_p);
    result_by_p.clear();

    // l , m = (7, 2).
    // p = 0
    tmp_double = 1126.125*sin_i_half_5*cos_i_half_9;
    result_by_lmp.set(c_Key3(7, 2, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 86.625*(91.0*sin_i_half_4 - 65.0*sin_i_half_2 + 10.0)*sin_i_half_3*cos_i_half_7;
    result_by_lmp.set(c_Key3(7, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 0.9228515625*std::pow(cos_i + 1.0, 6)*sin_i + 2976.75*sin_i_half_9*cos_i_half_5 - 9922.5*sin_i_half_7*cos_i_half_7 + 8505.0*sin_i_half_5*cos_i_half_9 - 2126.25*sin_i_half_3*cos_i_half_11;
    result_by_lmp.set(c_Key3(7, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -0.9228515625*std::pow(cos_i + 1.0, 6)*sin_i + 1102.5*sin_i_half_11*cos_i_half_3 - 8268.75*sin_i_half_9*cos_i_half_5 + 16537.5*sin_i_half_7*cos_i_half_7 - 11025.0*sin_i_half_5*cos_i_half_9 + 2362.5*sin_i_half_3*cos_i_half_11;
    result_by_lmp.set(c_Key3(7, 2, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 39.375*(343.0*sin_i_half_10 - 620.0*sin_i_half_8 + 280.0*sin_i_half_6 - 420.0*sin_i_half_4*cos_i_half_6 - 238.0*cos_i_half_10 + 210.0*cos_i_half_8)*sin_i_half_3*cos_i_half;
    result_by_lmp.set(c_Key3(7, 2, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 23.625*(-455.0*sin_i_half_8 + 810.0*sin_i_half_6 - 360.0*sin_i_half_4 - 546.0*cos_i_half_8 + 420.0*cos_i_half_6)*sin_i_half_5*cos_i_half;
    result_by_lmp.set(c_Key3(7, 2, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 86.625*(-91.0*sin_i_half_4 + 117.0*sin_i_half_2 - 36.0)*sin_i_half_7*cos_i_half_3;
    result_by_lmp.set(c_Key3(7, 2, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = -1126.125*sin_i_half_9*cos_i_half_5;
    result_by_lmp.set(c_Key3(7, 2, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 2), result_by_p);
    result_by_p.clear();

    // l , m = (7, 3).
    // p = 0
    tmp_double = 5630.625*sin_i_half_4*cos_i_half_10;
    result_by_lmp.set(c_Key3(7, 3, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 0.845947265625*std::pow(cos_i + 1.0, 4)*(-629.0*cos_i + 338.0*cos_i_double - 91.0*cos_i_triple + 382.0);
    result_by_lmp.set(c_Key3(7, 3, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 118.125*(960.0*sin_i_half_8 - 1020.0*sin_i_half_6 + 270.0*sin_i_half_4 + 41.0*cos_i_half_8 - 40.0*cos_i_half_6)*cos_i_half_6;
    result_by_lmp.set(c_Key3(7, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 196.875*(826.0*sin_i_half_8 - 1176.0*sin_i_half_6 + 420.0*sin_i_half_4 + 175.0*cos_i_half_8 - 160.0*cos_i_half_6)*sin_i_half_2*cos_i_half_4;
    result_by_lmp.set(c_Key3(7, 3, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 196.875*(-595.0*sin_i_half_10 + 1595.0*sin_i_half_8 - 1420.0*sin_i_half_6 + 420.0*sin_i_half_4 + 406.0*cos_i_half_10 - 336.0*cos_i_half_8)*sin_i_half_4;
    result_by_lmp.set(c_Key3(7, 3, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 118.125*(311.0*sin_i_half_8 - 580.0*sin_i_half_6 + 270.0*sin_i_half_4 + 690.0*cos_i_half_8 - 480.0*cos_i_half_6)*sin_i_half_6;
    result_by_lmp.set(c_Key3(7, 3, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 433.125*(-91.0*sin_i_half_6 + 221.0*sin_i_half_4 - 175.0*sin_i_half_2 + 45.0)*sin_i_half_8;
    result_by_lmp.set(c_Key3(7, 3, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 5630.625*sin_i_half_10*cos_i_half_4;
    result_by_lmp.set(c_Key3(7, 3, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 3), result_by_p);
    result_by_p.clear();

    // l , m = (7, 4).
    // p = 0
    tmp_double = -22522.5*sin_i_half_3*cos_i_half_11;
    result_by_lmp.set(c_Key3(7, 4, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -40.60546875*std::pow(cos_i + 1.0, 6)*sin_i - 95287.5*sin_i_half_5*cos_i_half_9 + 57172.5*sin_i_half_3*cos_i_half_11;
    result_by_lmp.set(c_Key3(7, 4, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 40.60546875*std::pow(cos_i + 1.0, 6)*sin_i - 155925.0*sin_i_half_7*cos_i_half_7 + 233887.5*sin_i_half_5*cos_i_half_9 - 77962.5*sin_i_half_3*cos_i_half_11;
    result_by_lmp.set(c_Key3(7, 4, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = (-744975.0*sin_i_half_6 + 883575.0*sin_i_half_4 - 259875.0*sin_i_half_2 + 43312.5*cos_i_half_6)*sin_i_half_3*cos_i_half_5;
    result_by_lmp.set(c_Key3(7, 4, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = (-667012.5*sin_i_half_6 + 987525.0*sin_i_half_4 - 363825.0*sin_i_half_2 + 121275.0*cos_i_half_6)*sin_i_half_5*cos_i_half_3;
    result_by_lmp.set(c_Key3(7, 4, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 5197.5*(-61.0*sin_i_half_6 + 105.0*sin_i_half_4 - 45.0*sin_i_half_2 + 30.0*cos_i_half_6)*sin_i_half_7*cos_i_half;
    result_by_lmp.set(c_Key3(7, 4, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 1732.5*(91.0*sin_i_half_4 - 143.0*sin_i_half_2 + 55.0)*sin_i_half_9*cos_i_half;
    result_by_lmp.set(c_Key3(7, 4, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 22522.5*sin_i_half_11*cos_i_half_3;
    result_by_lmp.set(c_Key3(7, 4, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 4), result_by_p);
    result_by_p.clear();

    // l , m = (7, 5).
    // p = 0
    tmp_double = -67567.5*sin_i_half_2*cos_i_half_12;
    result_by_lmp.set(c_Key3(7, 5, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 40.60546875*std::pow(cos_i + 1.0, 5)*(91.0*sin_i_2 + 130.0*cos_i - 134.0);
    result_by_lmp.set(c_Key3(7, 5, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 30.4541015625*std::pow(cos_i + 1.0, 4)*(629.0*cos_i - 338.0*cos_i_double + 91.0*cos_i_triple - 382.0);
    result_by_lmp.set(c_Key3(7, 5, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 203.02734375*std::pow(cos_i - 1.0, 2)*std::pow(cos_i + 1.0, 3)*(91.0*sin_i_2 + 26.0*cos_i - 86.0);
    result_by_lmp.set(c_Key3(7, 5, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 203.02734375*std::pow(cos_i - 1.0, 3)*std::pow(cos_i + 1.0, 2)*(-91.0*sin_i_2 + 26.0*cos_i + 86.0);
    result_by_lmp.set(c_Key3(7, 5, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 15592.5*(91.0*sin_i_half_6 - 221.0*sin_i_half_4 + 175.0*sin_i_half_2 - 45.0)*sin_i_half_8;
    result_by_lmp.set(c_Key3(7, 5, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 1299.375*(91.0*sin_i_2 - 130.0*cos_i - 134.0)*sin_i_half_10;
    result_by_lmp.set(c_Key3(7, 5, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = -67567.5*sin_i_half_12*cos_i_half_2;
    result_by_lmp.set(c_Key3(7, 5, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 5), result_by_p);
    result_by_p.clear();

    // l , m = (7, 6).
    // p = 0
    tmp_double = 1055.7421875*std::pow(cos_i + 1.0, 6)*sin_i;
    result_by_lmp.set(c_Key3(7, 6, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -1055.7421875*std::pow(cos_i + 1.0, 6)*sin_i + 810810.0*sin_i_half_3*cos_i_half_11;
    result_by_lmp.set(c_Key3(7, 6, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 202702.5*(3.0 - 7.0*cos_i)*sin_i_half_3*cos_i_half_9;
    result_by_lmp.set(c_Key3(7, 6, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 337837.5*(1.0 - 7.0*cos_i)*sin_i_half_5*cos_i_half_7;
    result_by_lmp.set(c_Key3(7, 6, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = (2027025.0 - 4729725.0*cos_i_half_2)*sin_i_half_7*cos_i_half_5;
    result_by_lmp.set(c_Key3(7, 6, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = (810810.0 - 2837835.0*cos_i_half_2)*sin_i_half_9*cos_i_half_3;
    result_by_lmp.set(c_Key3(7, 6, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 1055.7421875*std::pow(cos_i - 1.0, 5)*(7.0*cos_i + 5.0)*sin_i;
    result_by_lmp.set(c_Key3(7, 6, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = -135135.0*sin_i_half_13*cos_i_half;
    result_by_lmp.set(c_Key3(7, 6, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 6), result_by_p);
    result_by_p.clear();

    // l , m = (7, 7).
    // p = 0
    tmp_double = 135135.0*cos_i_half_14;
    result_by_lmp.set(c_Key3(7, 7, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 945945.0*sin_i_half_2*cos_i_half_12;
    result_by_lmp.set(c_Key3(7, 7, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 2837835.0*sin_i_half_4*cos_i_half_10;
    result_by_lmp.set(c_Key3(7, 7, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 4729725.0*sin_i_half_6*cos_i_half_8;
    result_by_lmp.set(c_Key3(7, 7, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 4729725.0*sin_i_half_8*cos_i_half_6;
    result_by_lmp.set(c_Key3(7, 7, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 2837835.0*sin_i_half_10*cos_i_half_4;
    result_by_lmp.set(c_Key3(7, 7, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 945945.0*sin_i_half_12*cos_i_half_2;
    result_by_lmp.set(c_Key3(7, 7, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 135135.0*sin_i_half_14;
    result_by_lmp.set(c_Key3(7, 7, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 7), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l7_2(double obliquity)
{
    // Inclination Functions Calculated for l = 7.
    // Functions are taylor expanded and truncated at obliquity^2.

    //  Functions Calculated for l = 7.

    c_IntMap<c_Key3, double> result_by_lmp(12);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(8);
    // Optimizations

    c_IntMap<c_Key1, double> result_by_p(7);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (7, 0).
    // p = 3
    tmp_double = -1.09375*obliquity;
    result_by_lmp.set(c_Key3(7, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 1.09375*obliquity;
    result_by_lmp.set(c_Key3(7, 0, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 0), result_by_p);
    result_by_p.clear();

    // l , m = (7, 1).
    // p = 3
    tmp_double = -2.1875;
    result_by_lmp.set(c_Key3(7, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 1), result_by_p);
    result_by_p.clear();

    // l , m = (7, 2).
    // p = 2
    tmp_double = 59.0625*obliquity;
    result_by_lmp.set(c_Key3(7, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -59.0625*obliquity;
    result_by_lmp.set(c_Key3(7, 2, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 2), result_by_p);
    result_by_p.clear();

    // l , m = (7, 3).
    // p = 2
    tmp_double = 118.125;
    result_by_lmp.set(c_Key3(7, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 3), result_by_p);
    result_by_p.clear();

    // l , m = (7, 4).
    // p = 1
    tmp_double = -2598.75*obliquity;
    result_by_lmp.set(c_Key3(7, 4, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 2598.75*obliquity;
    result_by_lmp.set(c_Key3(7, 4, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 4), result_by_p);
    result_by_p.clear();

    // l , m = (7, 5).
    // p = 1
    tmp_double = -5197.5;
    result_by_lmp.set(c_Key3(7, 5, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 5), result_by_p);
    result_by_p.clear();

    // l , m = (7, 6).
    // p = 0
    tmp_double = 67567.5*obliquity;
    result_by_lmp.set(c_Key3(7, 6, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -67567.5*obliquity;
    result_by_lmp.set(c_Key3(7, 6, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 6), result_by_p);
    result_by_p.clear();

    // l , m = (7, 7).
    // p = 0
    tmp_double = 135135.0;
    result_by_lmp.set(c_Key3(7, 7, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 7), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l7_4(double obliquity)
{
    // Inclination Functions Calculated for l = 7.
    // Functions are taylor expanded and truncated at obliquity^4.

    //  Functions Calculated for l = 7.

    c_IntMap<c_Key3, double> result_by_lmp(26);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(8);
    // Optimizations
    double obliquity_2 = obliquity * obliquity;
    double obliquity_3 = obliquity_2 * obliquity;

    c_IntMap<c_Key1, double> result_by_p(7);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (7, 0).
    // p = 2
    tmp_double = -2.4609375*obliquity_3;
    result_by_lmp.set(c_Key3(7, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 7.5651041666666667*obliquity_3 - 1.09375*obliquity;
    result_by_lmp.set(c_Key3(7, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -7.5651041666666667*obliquity_3 + 1.09375*obliquity;
    result_by_lmp.set(c_Key3(7, 0, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 2.4609375*obliquity_3;
    result_by_lmp.set(c_Key3(7, 0, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 0), result_by_p);
    result_by_p.clear();

    // l , m = (7, 1).
    // p = 2
    tmp_double = -14.765625*obliquity_2;
    result_by_lmp.set(c_Key3(7, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 30.078125*obliquity_2 - 2.1875;
    result_by_lmp.set(c_Key3(7, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -15.3125*obliquity_2;
    result_by_lmp.set(c_Key3(7, 1, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 1), result_by_p);
    result_by_p.clear();

    // l , m = (7, 2).
    // p = 1
    tmp_double = 108.28125*obliquity_3;
    result_by_lmp.set(c_Key3(7, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -364.21875*obliquity_3 + 59.0625*obliquity;
    result_by_lmp.set(c_Key3(7, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 393.75*obliquity_3 - 59.0625*obliquity;
    result_by_lmp.set(c_Key3(7, 2, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -137.8125*obliquity_3;
    result_by_lmp.set(c_Key3(7, 2, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 2), result_by_p);
    result_by_p.clear();

    // l , m = (7, 3).
    // p = 1
    tmp_double = 649.6875*obliquity_2;
    result_by_lmp.set(c_Key3(7, 3, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 118.125 - 1387.96875*obliquity_2;
    result_by_lmp.set(c_Key3(7, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 738.28125*obliquity_2;
    result_by_lmp.set(c_Key3(7, 3, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 3), result_by_p);
    result_by_p.clear();

    // l , m = (7, 4).
    // p = 0
    tmp_double = -2815.3125*obliquity_3;
    result_by_lmp.set(c_Key3(7, 4, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 11477.8125*obliquity_3 - 2598.75*obliquity;
    result_by_lmp.set(c_Key3(7, 4, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -14076.5625*obliquity_3 + 2598.75*obliquity;
    result_by_lmp.set(c_Key3(7, 4, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 5414.0625*obliquity_3;
    result_by_lmp.set(c_Key3(7, 4, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 4), result_by_p);
    result_by_p.clear();

    // l , m = (7, 5).
    // p = 0
    tmp_double = -16891.875*obliquity_2;
    result_by_lmp.set(c_Key3(7, 5, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 40280.625*obliquity_2 - 5197.5;
    result_by_lmp.set(c_Key3(7, 5, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -23388.75*obliquity_2;
    result_by_lmp.set(c_Key3(7, 5, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 5), result_by_p);
    result_by_p.clear();

    // l , m = (7, 6).
    // p = 0
    tmp_double = -112612.5*obliquity_3 + 67567.5*obliquity;
    result_by_lmp.set(c_Key3(7, 6, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 213963.75*obliquity_3 - 67567.5*obliquity;
    result_by_lmp.set(c_Key3(7, 6, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -101351.25*obliquity_3;
    result_by_lmp.set(c_Key3(7, 6, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 6), result_by_p);
    result_by_p.clear();

    // l , m = (7, 7).
    // p = 0
    tmp_double = 135135.0 - 236486.25*obliquity_2;
    result_by_lmp.set(c_Key3(7, 7, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 236486.25*obliquity_2;
    result_by_lmp.set(c_Key3(7, 7, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 7), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l7_off(double obliquity)
{
    // Inclination Functions Calculated for l = 7.
    // Functions assume `obliquity=0.0`. However some modes are still activated (non-zero `F(l,m,p)`).

    //  Functions Calculated for l = 7.

    c_IntMap<c_Key3, double> result_by_lmp(4);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(4);
    // Optimizations

    c_IntMap<c_Key1, double> result_by_p(7);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (7, 1).
    // p = 3
    tmp_double = -2.1875;
    result_by_lmp.set(c_Key3(7, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 1), result_by_p);
    result_by_p.clear();

    // l , m = (7, 3).
    // p = 2
    tmp_double = 118.125;
    result_by_lmp.set(c_Key3(7, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 3), result_by_p);
    result_by_p.clear();

    // l , m = (7, 5).
    // p = 1
    tmp_double = -5197.5;
    result_by_lmp.set(c_Key3(7, 5, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 5), result_by_p);
    result_by_p.clear();

    // l , m = (7, 7).
    // p = 0
    tmp_double = 135135.0;
    result_by_lmp.set(c_Key3(7, 7, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(7, 7), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}
