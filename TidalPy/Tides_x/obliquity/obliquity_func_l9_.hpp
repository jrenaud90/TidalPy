#pragma once

#include "obliquity_common_.hpp"

ObliquityFuncOutput c_obliquity_function_l9_gen(double obliquity)
{
    // Inclination Functions Calculated for l = 9.
    // Functions are fully general with no assumptions on the magnitude of obliquity.

    //  Functions Calculated for l = 9.

    c_IntMap<c_Key3, double> result_by_lmp(100);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(10);
    // Optimizations
    double cos_i_half = std::cos((1.0/2.0)*obliquity);
    double sin_i_half = std::sin((1.0/2.0)*obliquity);
    double cos_i_half_2 = cos_i_half * cos_i_half;
    double cos_i_half_4 = cos_i_half_2 * cos_i_half_2;
    double cos_i_half_8 = cos_i_half_4 * cos_i_half_4;
    double cos_i_half_9 = cos_i_half_8 * cos_i_half;
    double sin_i_half_2 = sin_i_half * sin_i_half;
    double sin_i_half_4 = sin_i_half_2 * sin_i_half_2;
    double sin_i_half_8 = sin_i_half_4 * sin_i_half_4;
    double sin_i_half_9 = sin_i_half_8 * sin_i_half;
    double sin_i = std::sin(obliquity);
    double sin_i_2 = sin_i * sin_i;
    double sin_i_3 = sin_i_2 * sin_i;
    double sin_i_6 = sin_i_3 * sin_i_3;
    double sin_i_7 = sin_i_6 * sin_i;
    double cos_i_half_5 = cos_i_half_4 * cos_i_half;
    double cos_i_half_3 = cos_i_half_2 * cos_i_half;
    double cos_i_half_6 = cos_i_half_3 * cos_i_half_3;
    double sin_i_half_5 = sin_i_half_4 * sin_i_half;
    double sin_i_half_3 = sin_i_half_2 * sin_i_half;
    double sin_i_half_6 = sin_i_half_3 * sin_i_half_3;
    double cos_i_half_10 = cos_i_half_5 * cos_i_half_5;
    double cos_i_half_12 = cos_i_half_6 * cos_i_half_6;
    double sin_i_half_10 = sin_i_half_5 * sin_i_half_5;
    double sin_i_half_12 = sin_i_half_6 * sin_i_half_6;
    double cos_i = std::cos(obliquity);
    double cos_i_half_7 = cos_i_half_6 * cos_i_half;
    double cos_i_half_11 = cos_i_half_10 * cos_i_half;
    double cos_i_half_13 = cos_i_half_12 * cos_i_half;
    double cos_i_half_14 = cos_i_half_7 * cos_i_half_7;
    double cos_i_half_15 = cos_i_half_14 * cos_i_half;
    double sin_i_half_7 = sin_i_half_6 * sin_i_half;
    double sin_i_half_11 = sin_i_half_10 * sin_i_half;
    double sin_i_half_13 = sin_i_half_12 * sin_i_half;
    double sin_i_half_14 = sin_i_half_7 * sin_i_half_7;
    double sin_i_half_15 = sin_i_half_14 * sin_i_half;
    double cos_i_half_16 = cos_i_half_8 * cos_i_half_8;
    double sin_i_half_16 = sin_i_half_8 * sin_i_half_8;
    double sin_i_4 = sin_i_2 * sin_i_2;
    double cos_i_double = std::cos(2*obliquity);
    double cos_i_triple = std::cos(3*obliquity);

    c_IntMap<c_Key1, double> result_by_p(9);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (9, 0).
    // p = 0
    tmp_double = 94.9609375*sin_i_half_9*cos_i_half_9;
    result_by_lmp.set(c_Key3(9, 0, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 0.0981903076171875*(16.0 - 17.0*sin_i_2)*sin_i_7;
    result_by_lmp.set(c_Key3(9, 0, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 20.109375*(121.0*sin_i_half_8 - 186.0*sin_i_half_6 + 72.0*sin_i_half_4 + 49.0*cos_i_half_8 - 42.0*cos_i_half_6)*sin_i_half_5*cos_i_half_5;
    result_by_lmp.set(c_Key3(9, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 18.046875*(137.0*sin_i_half_12 - 243.0*sin_i_half_10 + 108.0*sin_i_half_8 - 168.0*sin_i_half_6*cos_i_half_6 + 108.0*sin_i_half_4*cos_i_half_8 + 29.0*cos_i_half_12 - 27.0*cos_i_half_10)*sin_i_half_3*cos_i_half_3;
    result_by_lmp.set(c_Key3(9, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 0.0048065185546875*std::pow(cos_i + 1.0, 8)*sin_i + 2.4609375*std::pow(sin_i_half, 17)*cos_i_half - 88.59375*sin_i_half_15*cos_i_half_3 + 826.875*sin_i_half_13*cos_i_half_5 - 2894.0625*sin_i_half_11*cos_i_half_7 + 4341.09375*sin_i_half_9*cos_i_half_9 - 2894.0625*sin_i_half_7*cos_i_half_11 + 826.875*sin_i_half_5*cos_i_half_13 - 88.59375*sin_i_half_3*cos_i_half_15;
    result_by_lmp.set(c_Key3(9, 0, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = -0.0048065185546875*std::pow(cos_i + 1.0, 8)*sin_i - 2.4609375*std::pow(sin_i_half, 17)*cos_i_half + 88.59375*sin_i_half_15*cos_i_half_3 - 826.875*sin_i_half_13*cos_i_half_5 + 2894.0625*sin_i_half_11*cos_i_half_7 - 4341.09375*sin_i_half_9*cos_i_half_9 + 2894.0625*sin_i_half_7*cos_i_half_11 - 826.875*sin_i_half_5*cos_i_half_13 + 88.59375*sin_i_half_3*cos_i_half_15;
    result_by_lmp.set(c_Key3(9, 0, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 18.046875*(-137.0*sin_i_half_12 + 243.0*sin_i_half_10 - 108.0*sin_i_half_8 + 168.0*sin_i_half_6*cos_i_half_6 - 108.0*sin_i_half_4*cos_i_half_8 - 29.0*cos_i_half_12 + 27.0*cos_i_half_10)*sin_i_half_3*cos_i_half_3;
    result_by_lmp.set(c_Key3(9, 0, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 20.109375*(-121.0*sin_i_half_8 + 186.0*sin_i_half_6 - 72.0*sin_i_half_4 - 49.0*cos_i_half_8 + 42.0*cos_i_half_6)*sin_i_half_5*cos_i_half_5;
    result_by_lmp.set(c_Key3(9, 0, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 0.0981903076171875*(17.0*sin_i_2 - 16.0)*sin_i_7;
    result_by_lmp.set(c_Key3(9, 0, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = -94.9609375*sin_i_half_9*cos_i_half_9;
    result_by_lmp.set(c_Key3(9, 0, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 0), result_by_p);
    result_by_p.clear();

    // l , m = (9, 1).
    // p = 0
    tmp_double = 854.6484375*sin_i_half_8*cos_i_half_10;
    result_by_lmp.set(c_Key3(9, 1, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 0.0981903076171875*std::pow(cos_i - 1.0, 3)*std::pow(cos_i + 1.0, 4)*(153.0*sin_i_2 + 34.0*cos_i - 146.0);
    result_by_lmp.set(c_Key3(9, 1, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 100.546875*(243.0*sin_i_half_8 - 348.0*sin_i_half_6 + 126.0*sin_i_half_4 + 63.0*cos_i_half_8 - 56.0*cos_i_half_6)*sin_i_half_4*cos_i_half_6;
    result_by_lmp.set(c_Key3(9, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 54.140625*(579.0*sin_i_half_12 - 984.0*sin_i_half_10 + 420.0*sin_i_half_8 - 480.0*sin_i_half_6*cos_i_half_6 + 225.0*sin_i_half_4*cos_i_half_8 + 42.0*cos_i_half_12 - 40.0*cos_i_half_10)*sin_i_half_2*cos_i_half_4;
    result_by_lmp.set(c_Key3(9, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -16943.5546875*std::pow(sin_i_half, 18) + 48246.6796875*sin_i_half_16 - 45773.4375*sin_i_half_14 + 14470.3125*sin_i_half_12 - 34728.75*sin_i_half_10*cos_i_half_8 + 36175.78125*sin_i_half_8*cos_i_half_10 - 16537.5*sin_i_half_6*cos_i_half_12 + 3100.78125*sin_i_half_4*cos_i_half_14 + 199.3359375*std::pow(cos_i_half, 18) - 196.875*cos_i_half_16;
    result_by_lmp.set(c_Key3(9, 1, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 2.4609375*(1341.0*sin_i_half_16 - 2600.0*sin_i_half_14 + 1260.0*sin_i_half_12 - 6720.0*sin_i_half_10*cos_i_half_6 + 14700.0*sin_i_half_8*cos_i_half_8 - 14112.0*sin_i_half_6*cos_i_half_10 + 5880.0*sin_i_half_4*cos_i_half_12 + 1005.0*cos_i_half_16 - 960.0*cos_i_half_14)*sin_i_half_2;
    result_by_lmp.set(c_Key3(9, 1, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 54.140625*(-267.0*sin_i_half_14 + 757.0*sin_i_half_12 - 715.0*sin_i_half_10 + 225.0*sin_i_half_8 - 480.0*sin_i_half_6*cos_i_half_8 + 420.0*sin_i_half_4*cos_i_half_10 + 159.0*cos_i_half_14 - 144.0*cos_i_half_12)*sin_i_half_4;
    result_by_lmp.set(c_Key3(9, 1, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 100.546875*(189.0*sin_i_half_8 - 308.0*sin_i_half_6 + 126.0*sin_i_half_4 + 117.0*cos_i_half_8 - 96.0*cos_i_half_6)*sin_i_half_6*cos_i_half_4;
    result_by_lmp.set(c_Key3(9, 1, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 0.0981903076171875*std::pow(cos_i - 1.0, 4)*std::pow(cos_i + 1.0, 3)*(-153.0*sin_i_2 + 34.0*cos_i + 146.0);
    result_by_lmp.set(c_Key3(9, 1, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = 854.6484375*sin_i_half_10*cos_i_half_8;
    result_by_lmp.set(c_Key3(9, 1, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 1), result_by_p);
    result_by_p.clear();

    // l , m = (9, 2).
    // p = 0
    tmp_double = -6837.1875*sin_i_half_7*cos_i_half_11;
    result_by_lmp.set(c_Key3(9, 2, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 402.1875*(-153.0*sin_i_half_4 + 119.0*sin_i_half_2 - 21.0)*sin_i_half_5*cos_i_half_9;
    result_by_lmp.set(c_Key3(9, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 402.1875*(-528.0*sin_i_half_8 + 693.0*sin_i_half_6 - 231.0*sin_i_half_4 - 84.0*cos_i_half_8 + 77.0*cos_i_half_6)*sin_i_half_3*cos_i_half_7;
    result_by_lmp.set(c_Key3(9, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -0.4229736328125*std::pow(cos_i + 1.0, 8)*sin_i - 14293.125*sin_i_half_13*cos_i_half_5 + 100051.875*sin_i_half_11*cos_i_half_7 - 214396.875*sin_i_half_9*cos_i_half_9 + 178664.0625*sin_i_half_7*cos_i_half_11 - 59554.6875*sin_i_half_5*cos_i_half_13 + 7146.5625*sin_i_half_3*cos_i_half_15;
    result_by_lmp.set(c_Key3(9, 2, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 0.4229736328125*std::pow(cos_i + 1.0, 8)*sin_i - 3248.4375*sin_i_half_15*cos_i_half_3 + 45478.125*sin_i_half_13*cos_i_half_5 - 191008.125*sin_i_half_11*cos_i_half_7 + 318346.875*sin_i_half_9*cos_i_half_9 - 227390.625*sin_i_half_7*cos_i_half_11 + 68217.1875*sin_i_half_5*cos_i_half_13 - 7579.6875*sin_i_half_3*cos_i_half_15;
    result_by_lmp.set(c_Key3(9, 2, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 216.5625*(-351.0*sin_i_half_14 + 665.0*sin_i_half_12 - 315.0*sin_i_half_10 + 1050.0*sin_i_half_8*cos_i_half_6 - 1470.0*sin_i_half_6*cos_i_half_8 + 882.0*sin_i_half_4*cos_i_half_10 + 225.0*cos_i_half_14 - 210.0*cos_i_half_12)*sin_i_half_3*cos_i_half;
    result_by_lmp.set(c_Key3(9, 2, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 216.5625*(309.0*sin_i_half_12 - 583.0*sin_i_half_10 + 275.0*sin_i_half_8 - 825.0*sin_i_half_6*cos_i_half_6 + 990.0*sin_i_half_4*cos_i_half_8 + 528.0*cos_i_half_12 - 462.0*cos_i_half_10)*sin_i_half_5*cos_i_half;
    result_by_lmp.set(c_Key3(9, 2, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 402.1875*(315.0*sin_i_half_8 - 539.0*sin_i_half_6 + 231.0*sin_i_half_4 + 297.0*cos_i_half_8 - 231.0*cos_i_half_6)*sin_i_half_7*cos_i_half_3;
    result_by_lmp.set(c_Key3(9, 2, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 402.1875*(153.0*sin_i_half_4 - 187.0*sin_i_half_2 + 55.0)*sin_i_half_9*cos_i_half_5;
    result_by_lmp.set(c_Key3(9, 2, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = 6837.1875*sin_i_half_11*cos_i_half_7;
    result_by_lmp.set(c_Key3(9, 2, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 2), result_by_p);
    result_by_p.clear();

    // l , m = (9, 3).
    // p = 0
    tmp_double = -47860.3125*sin_i_half_6*cos_i_half_12;
    result_by_lmp.set(c_Key3(9, 3, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 16.4959716796875*std::pow(cos_i + 1.0, 3)*(51.0*sin_i_2 + 34.0*cos_i - 54.0)*sin_i_4;
    result_by_lmp.set(c_Key3(9, 3, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 8445.9375*(-187.0*sin_i_half_8 + 220.0*sin_i_half_6 - 66.0*sin_i_half_4 - 17.0*cos_i_half_8 + 16.0*cos_i_half_6)*sin_i_half_2*cos_i_half_8;
    result_by_lmp.set(c_Key3(9, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 216.5625*(-13101.0*sin_i_half_12 + 19602.0*sin_i_half_10 - 7425.0*sin_i_half_8 + 4400.0*sin_i_half_6*cos_i_half_6 - 990.0*sin_i_half_4*cos_i_half_8 - 73.0*cos_i_half_12 + 72.0*cos_i_half_10)*cos_i_half_6;
    result_by_lmp.set(c_Key3(9, 3, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 4547.8125*(-579.0*sin_i_half_12 + 984.0*sin_i_half_10 - 420.0*sin_i_half_8 + 480.0*sin_i_half_6*cos_i_half_6 - 225.0*sin_i_half_4*cos_i_half_8 - 42.0*cos_i_half_12 + 40.0*cos_i_half_10)*sin_i_half_2*cos_i_half_4;
    result_by_lmp.set(c_Key3(9, 3, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 4547.8125*(267.0*sin_i_half_14 - 757.0*sin_i_half_12 + 715.0*sin_i_half_10 - 225.0*sin_i_half_8 + 480.0*sin_i_half_6*cos_i_half_8 - 420.0*sin_i_half_4*cos_i_half_10 - 159.0*cos_i_half_14 + 144.0*cos_i_half_12)*sin_i_half_4;
    result_by_lmp.set(c_Key3(9, 3, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 216.5625*(-1063.0*sin_i_half_12 + 2052.0*sin_i_half_10 - 990.0*sin_i_half_8 + 4400.0*sin_i_half_6*cos_i_half_6 - 7425.0*sin_i_half_4*cos_i_half_8 - 5676.0*cos_i_half_12 + 4752.0*cos_i_half_10)*sin_i_half_6;
    result_by_lmp.set(c_Key3(9, 3, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 8445.9375*(83.0*sin_i_half_10 - 231.0*sin_i_half_8 + 214.0*sin_i_half_6 - 66.0*sin_i_half_4 - 121.0*cos_i_half_10 + 88.0*cos_i_half_8)*sin_i_half_8;
    result_by_lmp.set(c_Key3(9, 3, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 16.4959716796875*std::pow(cos_i - 1.0, 3)*(-51.0*sin_i_2 + 34.0*cos_i + 54.0)*sin_i_4;
    result_by_lmp.set(c_Key3(9, 3, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = -47860.3125*sin_i_half_12*cos_i_half_6;
    result_by_lmp.set(c_Key3(9, 3, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 3), result_by_p);
    result_by_p.clear();

    // l , m = (9, 4).
    // p = 0
    tmp_double = 287161.875*sin_i_half_5*cos_i_half_13;
    result_by_lmp.set(c_Key3(9, 4, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 16891.875*(153.0*sin_i_half_4 - 85.0*sin_i_half_2 + 10.0)*sin_i_half_3*cos_i_half_11;
    result_by_lmp.set(c_Key3(9, 4, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 32.991943359375*std::pow(cos_i + 1.0, 8)*sin_i + 2415538.125*sin_i_half_9*cos_i_half_9 - 4831076.25*sin_i_half_7*cos_i_half_11 + 2635132.5*sin_i_half_5*cos_i_half_13 - 439188.75*sin_i_half_3*cos_i_half_15;
    result_by_lmp.set(c_Key3(9, 4, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -32.991943359375*std::pow(cos_i + 1.0, 8)*sin_i + 2229727.5*sin_i_half_11*cos_i_half_7 - 8361478.125*sin_i_half_9*cos_i_half_9 + 9290531.25*sin_i_half_7*cos_i_half_11 - 3716212.5*sin_i_half_5*cos_i_half_13 + 506756.25*sin_i_half_3*cos_i_half_15;
    result_by_lmp.set(c_Key3(9, 4, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 118243.125*(189.0*sin_i_half_10 - 300.0*sin_i_half_8 + 120.0*sin_i_half_6 - 90.0*sin_i_half_4*cos_i_half_6 - 27.0*cos_i_half_10 + 25.0*cos_i_half_8)*sin_i_half_3*cos_i_half_5;
    result_by_lmp.set(c_Key3(9, 4, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 118243.125*(117.0*sin_i_half_10 - 205.0*sin_i_half_8 + 90.0*sin_i_half_6 - 120.0*sin_i_half_4*cos_i_half_6 - 69.0*cos_i_half_10 + 60.0*cos_i_half_8)*sin_i_half_5*cos_i_half_3;
    result_by_lmp.set(c_Key3(9, 4, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 16891.875*(251.0*sin_i_half_10 - 470.0*sin_i_half_8 + 220.0*sin_i_half_6 - 550.0*sin_i_half_4*cos_i_half_6 - 627.0*cos_i_half_10 + 495.0*cos_i_half_8)*sin_i_half_7*cos_i_half;
    result_by_lmp.set(c_Key3(9, 4, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 16891.875*(-183.0*sin_i_half_8 + 338.0*sin_i_half_6 - 156.0*sin_i_half_4 - 429.0*cos_i_half_8 + 286.0*cos_i_half_6)*sin_i_half_9*cos_i_half;
    result_by_lmp.set(c_Key3(9, 4, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 16891.875*(-153.0*sin_i_half_4 + 221.0*sin_i_half_2 - 78.0)*sin_i_half_11*cos_i_half_3;
    result_by_lmp.set(c_Key3(9, 4, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = -287161.875*sin_i_half_13*cos_i_half_5;
    result_by_lmp.set(c_Key3(9, 4, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 4), result_by_p);
    result_by_p.clear();

    // l , m = (9, 5).
    // p = 0
    tmp_double = 1435809.375*sin_i_half_4*cos_i_half_14;
    result_by_lmp.set(c_Key3(9, 5, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 41.23992919921875*std::pow(cos_i + 1.0, 6)*(-1303.0*cos_i + 646.0*cos_i_double - 153.0*cos_i_triple + 810.0);
    result_by_lmp.set(c_Key3(9, 5, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 16891.875*(3003.0*sin_i_half_8 - 2548.0*sin_i_half_6 + 546.0*sin_i_half_4 + 57.0*cos_i_half_8 - 56.0*cos_i_half_6)*cos_i_half_10;
    result_by_lmp.set(c_Key3(9, 5, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 591215.625*(187.0*sin_i_half_8 - 220.0*sin_i_half_6 + 66.0*sin_i_half_4 + 17.0*cos_i_half_8 - 16.0*cos_i_half_6)*sin_i_half_2*cos_i_half_8;
    result_by_lmp.set(c_Key3(9, 5, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 591215.625*(243.0*sin_i_half_8 - 348.0*sin_i_half_6 + 126.0*sin_i_half_4 + 63.0*cos_i_half_8 - 56.0*cos_i_half_6)*sin_i_half_4*cos_i_half_6;
    result_by_lmp.set(c_Key3(9, 5, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 591215.625*(189.0*sin_i_half_8 - 308.0*sin_i_half_6 + 126.0*sin_i_half_4 + 117.0*cos_i_half_8 - 96.0*cos_i_half_6)*sin_i_half_6*cos_i_half_4;
    result_by_lmp.set(c_Key3(9, 5, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 591215.625*(-83.0*sin_i_half_10 + 231.0*sin_i_half_8 - 214.0*sin_i_half_6 + 66.0*sin_i_half_4 + 121.0*cos_i_half_10 - 88.0*cos_i_half_8)*sin_i_half_8;
    result_by_lmp.set(c_Key3(9, 5, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 16891.875*(603.0*sin_i_half_8 - 1148.0*sin_i_half_6 + 546.0*sin_i_half_4 + 2457.0*cos_i_half_8 - 1456.0*cos_i_half_6)*sin_i_half_10;
    result_by_lmp.set(c_Key3(9, 5, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 84459.375*(-153.0*sin_i_half_6 + 391.0*sin_i_half_4 - 329.0*sin_i_half_2 + 91.0)*sin_i_half_12;
    result_by_lmp.set(c_Key3(9, 5, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = 1435809.375*sin_i_half_14*cos_i_half_4;
    result_by_lmp.set(c_Key3(9, 5, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 5), result_by_p);
    result_by_p.clear();

    // l , m = (9, 6).
    // p = 0
    tmp_double = -5743237.5*sin_i_half_3*cos_i_half_15;
    result_by_lmp.set(c_Key3(9, 6, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -1979.5166015625*std::pow(cos_i + 1.0, 8)*sin_i - 35472937.5*sin_i_half_5*cos_i_half_13 + 15202687.5*sin_i_half_3*cos_i_half_15;
    result_by_lmp.set(c_Key3(9, 6, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 1979.5166015625*std::pow(cos_i + 1.0, 8)*sin_i - 92229637.5*sin_i_half_7*cos_i_half_11 + 92229637.5*sin_i_half_5*cos_i_half_13 - 21283762.5*sin_i_half_3*cos_i_half_15;
    result_by_lmp.set(c_Key3(9, 6, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 2364862.5*(-199.0*sin_i_half_6 + 189.0*sin_i_half_4 - 45.0*sin_i_half_2 + 5.0*cos_i_half_6)*sin_i_half_3*cos_i_half_9;
    result_by_lmp.set(c_Key3(9, 6, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 7094587.5*(-95.0*sin_i_half_6 + 115.0*sin_i_half_4 - 35.0*sin_i_half_2 + 7.0*cos_i_half_6)*sin_i_half_5*cos_i_half_7;
    result_by_lmp.set(c_Key3(9, 6, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 7094587.5*(-87.0*sin_i_half_6 + 125.0*sin_i_half_4 - 45.0*sin_i_half_2 + 15.0*cos_i_half_6)*sin_i_half_7*cos_i_half_5;
    result_by_lmp.set(c_Key3(9, 6, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 2364862.5*(-149.0*sin_i_half_6 + 243.0*sin_i_half_4 - 99.0*sin_i_half_2 + 55.0*cos_i_half_6)*sin_i_half_9*cos_i_half_3;
    result_by_lmp.set(c_Key3(9, 6, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 1013512.5*(-113.0*sin_i_half_6 + 203.0*sin_i_half_4 - 91.0*sin_i_half_2 + 91.0*cos_i_half_6)*sin_i_half_11*cos_i_half;
    result_by_lmp.set(c_Key3(9, 6, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 1013512.5*(51.0*sin_i_half_4 - 85.0*sin_i_half_2 + 35.0)*sin_i_half_13*cos_i_half;
    result_by_lmp.set(c_Key3(9, 6, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = 5743237.5*sin_i_half_15*cos_i_half_3;
    result_by_lmp.set(c_Key3(9, 6, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 6), result_by_p);
    result_by_p.clear();

    // l , m = (9, 7).
    // p = 0
    tmp_double = -17229712.5*sin_i_half_2*cos_i_half_16;
    result_by_lmp.set(c_Key3(9, 7, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 1979.5166015625*std::pow(cos_i + 1.0, 7)*(153.0*sin_i_2 + 238.0*cos_i - 242.0);
    result_by_lmp.set(c_Key3(9, 7, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 1979.5166015625*std::pow(cos_i + 1.0, 6)*(1303.0*cos_i - 646.0*cos_i_double + 153.0*cos_i_triple - 810.0);
    result_by_lmp.set(c_Key3(9, 7, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = (361823962.5*sin_i_2 + 241215975.0*cos_i - 383107725.0)*sin_i_half_4*cos_i_half_10;
    result_by_lmp.set(c_Key3(9, 7, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 3547293.75*(153.0*sin_i_2 + 34.0*cos_i - 146.0)*sin_i_half_6*cos_i_half_8;
    result_by_lmp.set(c_Key3(9, 7, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 3547293.75*(153.0*sin_i_2 - 34.0*cos_i - 146.0)*sin_i_half_8*cos_i_half_6;
    result_by_lmp.set(c_Key3(9, 7, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = (361823962.5*sin_i_2 - 241215975.0*cos_i - 383107725.0)*sin_i_half_10*cos_i_half_4;
    result_by_lmp.set(c_Key3(9, 7, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = (620269650.0*sin_i_half_6 - 1585133550.0*sin_i_half_4 + 1333782450.0*sin_i_half_2 - 368918550.0)*sin_i_half_12;
    result_by_lmp.set(c_Key3(9, 7, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 253378.125*(153.0*sin_i_2 - 238.0*cos_i - 242.0)*sin_i_half_14;
    result_by_lmp.set(c_Key3(9, 7, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = -17229712.5*sin_i_half_16*cos_i_half_2;
    result_by_lmp.set(c_Key3(9, 7, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 7), result_by_p);
    result_by_p.clear();

    // l , m = (9, 8).
    // p = 0
    tmp_double = 67303.564453125*std::pow(cos_i + 1.0, 8)*sin_i;
    result_by_lmp.set(c_Key3(9, 8, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -67303.564453125*std::pow(cos_i + 1.0, 8)*sin_i + 275675400.0*sin_i_half_3*cos_i_half_15;
    result_by_lmp.set(c_Key3(9, 8, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = (344594250.0 - 620269650.0*cos_i)*sin_i_half_3*cos_i_half_13;
    result_by_lmp.set(c_Key3(9, 8, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = (482431950.0 - 1447295850.0*cos_i)*sin_i_half_5*cos_i_half_11;
    result_by_lmp.set(c_Key3(9, 8, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = (241215975.0 - 2170943775.0*cos_i)*sin_i_half_7*cos_i_half_9;
    result_by_lmp.set(c_Key3(9, 8, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = (-2170943775.0*cos_i - 241215975.0)*sin_i_half_9*cos_i_half_7;
    result_by_lmp.set(c_Key3(9, 8, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = (-1447295850.0*cos_i - 482431950.0)*sin_i_half_11*cos_i_half_5;
    result_by_lmp.set(c_Key3(9, 8, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = (-620269650.0*cos_i - 344594250.0)*sin_i_half_13*cos_i_half_3;
    result_by_lmp.set(c_Key3(9, 8, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 67303.564453125*std::pow(cos_i - 1.0, 7)*(9.0*cos_i + 7.0)*sin_i;
    result_by_lmp.set(c_Key3(9, 8, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = -34459425.0*std::pow(sin_i_half, 17)*cos_i_half;
    result_by_lmp.set(c_Key3(9, 8, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 8), result_by_p);
    result_by_p.clear();

    // l , m = (9, 9).
    // p = 0
    tmp_double = 34459425.0*std::pow(cos_i_half, 18);
    result_by_lmp.set(c_Key3(9, 9, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 310134825.0*sin_i_half_2*cos_i_half_16;
    result_by_lmp.set(c_Key3(9, 9, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 1240539300.0*sin_i_half_4*cos_i_half_14;
    result_by_lmp.set(c_Key3(9, 9, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 2894591700.0*sin_i_half_6*cos_i_half_12;
    result_by_lmp.set(c_Key3(9, 9, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 4341887550.0*sin_i_half_8*cos_i_half_10;
    result_by_lmp.set(c_Key3(9, 9, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 4341887550.0*sin_i_half_10*cos_i_half_8;
    result_by_lmp.set(c_Key3(9, 9, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 2894591700.0*sin_i_half_12*cos_i_half_6;
    result_by_lmp.set(c_Key3(9, 9, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 1240539300.0*sin_i_half_14*cos_i_half_4;
    result_by_lmp.set(c_Key3(9, 9, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 310134825.0*sin_i_half_16*cos_i_half_2;
    result_by_lmp.set(c_Key3(9, 9, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = 34459425.0*std::pow(sin_i_half, 18);
    result_by_lmp.set(c_Key3(9, 9, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 9), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l9_2(double obliquity)
{
    // Inclination Functions Calculated for l = 9.
    // Functions are taylor expanded and truncated at obliquity^2.

    //  Functions Calculated for l = 9.

    c_IntMap<c_Key3, double> result_by_lmp(15);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(10);
    // Optimizations

    c_IntMap<c_Key1, double> result_by_p(9);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (9, 0).
    // p = 4
    tmp_double = 1.23046875*obliquity;
    result_by_lmp.set(c_Key3(9, 0, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = -1.23046875*obliquity;
    result_by_lmp.set(c_Key3(9, 0, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 0), result_by_p);
    result_by_p.clear();

    // l , m = (9, 1).
    // p = 4
    tmp_double = 2.4609375;
    result_by_lmp.set(c_Key3(9, 1, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 1), result_by_p);
    result_by_p.clear();

    // l , m = (9, 2).
    // p = 3
    tmp_double = -108.28125*obliquity;
    result_by_lmp.set(c_Key3(9, 2, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 108.28125*obliquity;
    result_by_lmp.set(c_Key3(9, 2, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 2), result_by_p);
    result_by_p.clear();

    // l , m = (9, 3).
    // p = 3
    tmp_double = -216.5625;
    result_by_lmp.set(c_Key3(9, 3, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 3), result_by_p);
    result_by_p.clear();

    // l , m = (9, 4).
    // p = 2
    tmp_double = 8445.9375*obliquity;
    result_by_lmp.set(c_Key3(9, 4, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -8445.9375*obliquity;
    result_by_lmp.set(c_Key3(9, 4, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 4), result_by_p);
    result_by_p.clear();

    // l , m = (9, 5).
    // p = 2
    tmp_double = 16891.875;
    result_by_lmp.set(c_Key3(9, 5, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 5), result_by_p);
    result_by_p.clear();

    // l , m = (9, 6).
    // p = 1
    tmp_double = -506756.25*obliquity;
    result_by_lmp.set(c_Key3(9, 6, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 506756.25*obliquity;
    result_by_lmp.set(c_Key3(9, 6, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 6), result_by_p);
    result_by_p.clear();

    // l , m = (9, 7).
    // p = 1
    tmp_double = -1013512.5;
    result_by_lmp.set(c_Key3(9, 7, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 7), result_by_p);
    result_by_p.clear();

    // l , m = (9, 8).
    // p = 0
    tmp_double = 17229712.5*obliquity;
    result_by_lmp.set(c_Key3(9, 8, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -17229712.5*obliquity;
    result_by_lmp.set(c_Key3(9, 8, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 8), result_by_p);
    result_by_p.clear();

    // l , m = (9, 9).
    // p = 0
    tmp_double = 34459425.0;
    result_by_lmp.set(c_Key3(9, 9, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 9), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l9_4(double obliquity)
{
    // Inclination Functions Calculated for l = 9.
    // Functions are taylor expanded and truncated at obliquity^4.

    //  Functions Calculated for l = 9.

    c_IntMap<c_Key3, double> result_by_lmp(33);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(10);
    // Optimizations
    double obliquity_2 = obliquity * obliquity;
    double obliquity_3 = obliquity_2 * obliquity;

    c_IntMap<c_Key1, double> result_by_p(9);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (9, 0).
    // p = 3
    tmp_double = 4.51171875*obliquity_3;
    result_by_lmp.set(c_Key3(9, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -13.740234375*obliquity_3 + 1.23046875*obliquity;
    result_by_lmp.set(c_Key3(9, 0, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 13.740234375*obliquity_3 - 1.23046875*obliquity;
    result_by_lmp.set(c_Key3(9, 0, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = -4.51171875*obliquity_3;
    result_by_lmp.set(c_Key3(9, 0, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 0), result_by_p);
    result_by_p.clear();

    // l , m = (9, 1).
    // p = 3
    tmp_double = 27.0703125*obliquity_2;
    result_by_lmp.set(c_Key3(9, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 2.4609375 - 54.755859375*obliquity_2;
    result_by_lmp.set(c_Key3(9, 1, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 27.685546875*obliquity_2;
    result_by_lmp.set(c_Key3(9, 1, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 1), result_by_p);
    result_by_p.clear();

    // l , m = (9, 2).
    // p = 2
    tmp_double = -351.9140625*obliquity_3;
    result_by_lmp.set(c_Key3(9, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 1127.9296875*obliquity_3 - 108.28125*obliquity;
    result_by_lmp.set(c_Key3(9, 2, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -1182.0703125*obliquity_3 + 108.28125*obliquity;
    result_by_lmp.set(c_Key3(9, 2, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 406.0546875*obliquity_3;
    result_by_lmp.set(c_Key3(9, 2, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 2), result_by_p);
    result_by_p.clear();

    // l , m = (9, 3).
    // p = 2
    tmp_double = -2111.484375*obliquity_2;
    result_by_lmp.set(c_Key3(9, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 4385.390625*obliquity_2 - 216.5625;
    result_by_lmp.set(c_Key3(9, 3, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -2273.90625*obliquity_2;
    result_by_lmp.set(c_Key3(9, 3, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 3), result_by_p);
    result_by_p.clear();

    // l , m = (9, 4).
    // p = 1
    tmp_double = 21114.84375*obliquity_3;
    result_by_lmp.set(c_Key3(9, 4, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -73198.125*obliquity_3 + 8445.9375*obliquity;
    result_by_lmp.set(c_Key3(9, 4, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 81644.0625*obliquity_3 - 8445.9375*obliquity;
    result_by_lmp.set(c_Key3(9, 4, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -29560.78125*obliquity_3;
    result_by_lmp.set(c_Key3(9, 4, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 4), result_by_p);
    result_by_p.clear();

    // l , m = (9, 5).
    // p = 1
    tmp_double = 126689.0625*obliquity_2;
    result_by_lmp.set(c_Key3(9, 5, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 16891.875 - 274492.96875*obliquity_2;
    result_by_lmp.set(c_Key3(9, 5, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 147803.90625*obliquity_2;
    result_by_lmp.set(c_Key3(9, 5, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 5), result_by_p);
    result_by_p.clear();

    // l , m = (9, 6).
    // p = 0
    tmp_double = -717904.6875*obliquity_3;
    result_by_lmp.set(c_Key3(9, 6, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 2998307.8125*obliquity_3 - 506756.25*obliquity;
    result_by_lmp.set(c_Key3(9, 6, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -3758442.1875*obliquity_3 + 506756.25*obliquity;
    result_by_lmp.set(c_Key3(9, 6, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 1478039.0625*obliquity_3;
    result_by_lmp.set(c_Key3(9, 6, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 6), result_by_p);
    result_by_p.clear();

    // l , m = (9, 7).
    // p = 0
    tmp_double = -4307428.125*obliquity_2;
    result_by_lmp.set(c_Key3(9, 7, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 10388503.125*obliquity_2 - 1013512.5;
    result_by_lmp.set(c_Key3(9, 7, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -6081075.0*obliquity_2;
    result_by_lmp.set(c_Key3(9, 7, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 7), result_by_p);
    result_by_p.clear();

    // l , m = (9, 8).
    // p = 0
    tmp_double = -37331043.75*obliquity_3 + 17229712.5*obliquity;
    result_by_lmp.set(c_Key3(9, 8, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 71790468.75*obliquity_3 - 17229712.5*obliquity;
    result_by_lmp.set(c_Key3(9, 8, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -34459425.0*obliquity_3;
    result_by_lmp.set(c_Key3(9, 8, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 8), result_by_p);
    result_by_p.clear();

    // l , m = (9, 9).
    // p = 0
    tmp_double = 34459425.0 - 77533706.25*obliquity_2;
    result_by_lmp.set(c_Key3(9, 9, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 77533706.25*obliquity_2;
    result_by_lmp.set(c_Key3(9, 9, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 9), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l9_off(double obliquity)
{
    // Inclination Functions Calculated for l = 9.
    // Functions assume `obliquity=0.0`. However some modes are still activated (non-zero `F(l,m,p)`).

    //  Functions Calculated for l = 9.

    c_IntMap<c_Key3, double> result_by_lmp(5);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(5);
    // Optimizations

    c_IntMap<c_Key1, double> result_by_p(9);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (9, 1).
    // p = 4
    tmp_double = 2.4609375;
    result_by_lmp.set(c_Key3(9, 1, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 1), result_by_p);
    result_by_p.clear();

    // l , m = (9, 3).
    // p = 3
    tmp_double = -216.5625;
    result_by_lmp.set(c_Key3(9, 3, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 3), result_by_p);
    result_by_p.clear();

    // l , m = (9, 5).
    // p = 2
    tmp_double = 16891.875;
    result_by_lmp.set(c_Key3(9, 5, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 5), result_by_p);
    result_by_p.clear();

    // l , m = (9, 7).
    // p = 1
    tmp_double = -1013512.5;
    result_by_lmp.set(c_Key3(9, 7, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 7), result_by_p);
    result_by_p.clear();

    // l , m = (9, 9).
    // p = 0
    tmp_double = 34459425.0;
    result_by_lmp.set(c_Key3(9, 9, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(9, 9), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}
