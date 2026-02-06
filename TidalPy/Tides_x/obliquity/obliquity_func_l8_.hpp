#pragma once

#include "obliquity_common_.hpp"

ObliquityFuncOutput c_obliquity_function_l8_gen(double obliquity)
{
    // Inclination Functions Calculated for l = 8.
    // Functions are fully general with no assumptions on the magnitude of obliquity.

    //  Functions Calculated for l = 8.

    c_IntMap<c_Key3, double> result_by_lmp(81);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(9);
    // Optimizations
    double cos_i_half = std::cos((1.0/2.0)*obliquity);
    double sin_i_half = std::sin((1.0/2.0)*obliquity);
    double cos_i_half_2 = cos_i_half * cos_i_half;
    double cos_i_half_4 = cos_i_half_2 * cos_i_half_2;
    double cos_i_half_8 = cos_i_half_4 * cos_i_half_4;
    double sin_i_half_2 = sin_i_half * sin_i_half;
    double sin_i_half_4 = sin_i_half_2 * sin_i_half_2;
    double sin_i_half_8 = sin_i_half_4 * sin_i_half_4;
    double sin_i = std::sin(obliquity);
    double sin_i_2 = sin_i * sin_i;
    double sin_i_3 = sin_i_2 * sin_i;
    double sin_i_6 = sin_i_3 * sin_i_3;
    double cos_i_half_3 = cos_i_half_2 * cos_i_half;
    double cos_i_half_6 = cos_i_half_3 * cos_i_half_3;
    double sin_i_half_3 = sin_i_half_2 * sin_i_half;
    double sin_i_half_6 = sin_i_half_3 * sin_i_half_3;
    double cos_i_half_5 = cos_i_half_4 * cos_i_half;
    double cos_i_half_10 = cos_i_half_5 * cos_i_half_5;
    double cos_i_half_12 = cos_i_half_6 * cos_i_half_6;
    double cos_i_half_7 = cos_i_half_6 * cos_i_half;
    double cos_i_half_14 = cos_i_half_7 * cos_i_half_7;
    double sin_i_half_5 = sin_i_half_4 * sin_i_half;
    double sin_i_half_10 = sin_i_half_5 * sin_i_half_5;
    double sin_i_half_12 = sin_i_half_6 * sin_i_half_6;
    double sin_i_half_7 = sin_i_half_6 * sin_i_half;
    double sin_i_half_14 = sin_i_half_7 * sin_i_half_7;
    double cos_i_half_16 = cos_i_half_8 * cos_i_half_8;
    double sin_i_half_16 = sin_i_half_8 * sin_i_half_8;
    double cos_i_half_9 = cos_i_half_8 * cos_i_half;
    double cos_i = std::cos(obliquity);
    double cos_i_2 = cos_i * cos_i;
    double cos_i_half_11 = cos_i_half_10 * cos_i_half;
    double cos_i_half_13 = cos_i_half_12 * cos_i_half;
    double sin_i_half_9 = sin_i_half_8 * sin_i_half;
    double sin_i_half_11 = sin_i_half_10 * sin_i_half;
    double sin_i_half_13 = sin_i_half_12 * sin_i_half;
    double sin_i_half_15 = sin_i_half_14 * sin_i_half;
    double sin_i_4 = sin_i_2 * sin_i_2;
    double cos_i_double = std::cos(2*obliquity);
    double cos_i_triple = std::cos(3*obliquity);
    double sin_i_5 = sin_i_4 * sin_i;

    c_IntMap<c_Key1, double> result_by_p(8);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (8, 0).
    // p = 0
    tmp_double = 50.2734375*sin_i_half_8*cos_i_half_8;
    result_by_lmp.set(c_Key3(8, 0, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 0.104736328125*(14.0 - 15.0*sin_i_2)*sin_i_6;
    result_by_lmp.set(c_Key3(8, 0, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 10.828125*(93.0*sin_i_half_8 - 144.0*sin_i_half_6 + 56.0*sin_i_half_4 + 37.0*cos_i_half_8 - 32.0*cos_i_half_6)*sin_i_half_4*cos_i_half_4;
    result_by_lmp.set(c_Key3(8, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 9.84375*(-87.0*sin_i_half_14 + 243.0*sin_i_half_12 - 226.0*sin_i_half_10 + 70.0*sin_i_half_8 - 112.0*sin_i_half_6*cos_i_half_8 + 70.0*sin_i_half_4*cos_i_half_10 + 17.0*cos_i_half_14 - 16.0*cos_i_half_12)*sin_i_half_2;
    result_by_lmp.set(c_Key3(8, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 232.1484375*sin_i_half_16 - 446.25*sin_i_half_14 + 214.375*sin_i_half_12 - 857.5*sin_i_half_10*cos_i_half_6 + 1339.84375*sin_i_half_8*cos_i_half_8 - 857.5*sin_i_half_6*cos_i_half_10 + 214.375*sin_i_half_4*cos_i_half_12 + 17.7734375*cos_i_half_16 - 17.5*cos_i_half_14;
    result_by_lmp.set(c_Key3(8, 0, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 9.84375*(-87.0*sin_i_half_14 + 243.0*sin_i_half_12 - 226.0*sin_i_half_10 + 70.0*sin_i_half_8 - 112.0*sin_i_half_6*cos_i_half_8 + 70.0*sin_i_half_4*cos_i_half_10 + 17.0*cos_i_half_14 - 16.0*cos_i_half_12)*sin_i_half_2;
    result_by_lmp.set(c_Key3(8, 0, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 10.828125*(93.0*sin_i_half_8 - 144.0*sin_i_half_6 + 56.0*sin_i_half_4 + 37.0*cos_i_half_8 - 32.0*cos_i_half_6)*sin_i_half_4*cos_i_half_4;
    result_by_lmp.set(c_Key3(8, 0, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 0.104736328125*(14.0 - 15.0*sin_i_2)*sin_i_6;
    result_by_lmp.set(c_Key3(8, 0, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 50.2734375*sin_i_half_8*cos_i_half_8;
    result_by_lmp.set(c_Key3(8, 0, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 0), result_by_p);
    result_by_p.clear();

    // l , m = (8, 1).
    // p = 0
    tmp_double = -402.1875*sin_i_half_7*cos_i_half_9;
    result_by_lmp.set(c_Key3(8, 1, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 40.21875*(-20.0*cos_i_2 + 5.0*cos_i + 1.0)*sin_i_half_5*cos_i_half_7;
    result_by_lmp.set(c_Key3(8, 1, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 43.3125*(-210.0*sin_i_half_8 + 300.0*sin_i_half_6 - 108.0*sin_i_half_4 - 50.0*cos_i_half_8 + 45.0*cos_i_half_6)*sin_i_half_3*cos_i_half_5;
    result_by_lmp.set(c_Key3(8, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -0.076904296875*std::pow(cos_i + 1.0, 7)*sin_i - 236.25*sin_i_half_13*cos_i_half_3 + 2480.625*sin_i_half_11*cos_i_half_5 - 7441.875*sin_i_half_9*cos_i_half_7 + 8268.75*sin_i_half_7*cos_i_half_9 - 3543.75*sin_i_half_5*cos_i_half_11 + 531.5625*sin_i_half_3*cos_i_half_13;
    result_by_lmp.set(c_Key3(8, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 0.076904296875*std::pow(cos_i + 1.0, 7)*sin_i - 19.6875*sin_i_half_15*cos_i_half + 551.25*sin_i_half_13*cos_i_half_3 - 3858.75*sin_i_half_11*cos_i_half_5 + 9646.875*sin_i_half_9*cos_i_half_7 - 9646.875*sin_i_half_7*cos_i_half_9 + 3858.75*sin_i_half_5*cos_i_half_11 - 551.25*sin_i_half_3*cos_i_half_13;
    result_by_lmp.set(c_Key3(8, 1, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 19.6875*(208.0*sin_i_half_12 - 387.0*sin_i_half_10 + 180.0*sin_i_half_8 - 420.0*sin_i_half_6*cos_i_half_6 + 378.0*sin_i_half_4*cos_i_half_8 + 138.0*cos_i_half_12 - 126.0*cos_i_half_10)*sin_i_half_3*cos_i_half;
    result_by_lmp.set(c_Key3(8, 1, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 43.3125*(158.0*sin_i_half_8 - 261.0*sin_i_half_6 + 108.0*sin_i_half_4 + 102.0*cos_i_half_8 - 84.0*cos_i_half_6)*sin_i_half_5*cos_i_half_3;
    result_by_lmp.set(c_Key3(8, 1, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 40.21875*(20.0*cos_i_2 + 5.0*cos_i - 1.0)*sin_i_half_7*cos_i_half_5;
    result_by_lmp.set(c_Key3(8, 1, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 402.1875*sin_i_half_9*cos_i_half_7;
    result_by_lmp.set(c_Key3(8, 1, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 1), result_by_p);
    result_by_p.clear();

    // l , m = (8, 2).
    // p = 0
    tmp_double = -2815.3125*sin_i_half_6*cos_i_half_10;
    result_by_lmp.set(c_Key3(8, 2, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 43.9892578125*(1.0 - 2.0*cos_i)*std::pow(cos_i + 1.0, 2)*sin_i_4*cos_i;
    result_by_lmp.set(c_Key3(8, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 216.5625*(-321.0*sin_i_half_8 + 414.0*sin_i_half_6 - 135.0*sin_i_half_4 - 43.0*cos_i_half_8 + 40.0*cos_i_half_6)*sin_i_half_2*cos_i_half_6;
    result_by_lmp.set(c_Key3(8, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 19.6875*(-4872.0*sin_i_half_12 + 7812.0*sin_i_half_10 - 3150.0*sin_i_half_8 + 2400.0*sin_i_half_6*cos_i_half_6 - 736.0*cos_i_half_12 + 1410.0*cos_i_half_10 - 675.0*cos_i_half_8)*cos_i_half_4;
    result_by_lmp.set(c_Key3(8, 2, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 689.0625*(87.0*sin_i_half_14 - 243.0*sin_i_half_12 + 226.0*sin_i_half_10 - 70.0*sin_i_half_8 + 112.0*sin_i_half_6*cos_i_half_8 - 70.0*sin_i_half_4*cos_i_half_10 - 17.0*cos_i_half_14 + 16.0*cos_i_half_12)*sin_i_half_2;
    result_by_lmp.set(c_Key3(8, 2, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 19.6875*(-736.0*sin_i_half_12 + 1410.0*sin_i_half_10 - 675.0*sin_i_half_8 + 2400.0*sin_i_half_6*cos_i_half_6 - 3150.0*sin_i_half_4*cos_i_half_8 - 1722.0*cos_i_half_12 + 1512.0*cos_i_half_10)*sin_i_half_4;
    result_by_lmp.set(c_Key3(8, 2, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 216.5625*(178.0*sin_i_half_10 - 488.0*sin_i_half_8 + 445.0*sin_i_half_6 - 135.0*sin_i_half_4 - 186.0*cos_i_half_10 + 144.0*cos_i_half_8)*sin_i_half_6;
    result_by_lmp.set(c_Key3(8, 2, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = -43.9892578125*std::pow(cos_i - 1.0, 2)*(2.0*cos_i + 1.0)*sin_i_4*cos_i;
    result_by_lmp.set(c_Key3(8, 2, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = -2815.3125*sin_i_half_10*cos_i_half_6;
    result_by_lmp.set(c_Key3(8, 2, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 2), result_by_p);
    result_by_p.clear();

    // l , m = (8, 3).
    // p = 0
    tmp_double = 16891.875*sin_i_half_5*cos_i_half_11;
    result_by_lmp.set(c_Key3(8, 3, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 2815.3125*(-12.0*sin_i_2 - 9.0*cos_i + 13.0)*sin_i_half_3*cos_i_half_9;
    result_by_lmp.set(c_Key3(8, 3, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 5.07568359375*std::pow(cos_i + 1.0, 7)*sin_i + 85758.75*sin_i_half_9*cos_i_half_7 - 214396.875*sin_i_half_7*cos_i_half_9 + 142931.25*sin_i_half_5*cos_i_half_11 - 28586.25*sin_i_half_3*cos_i_half_13;
    result_by_lmp.set(c_Key3(8, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -5.07568359375*std::pow(cos_i + 1.0, 7)*sin_i + 54573.75*sin_i_half_11*cos_i_half_5 - 272868.75*sin_i_half_9*cos_i_half_7 + 389812.5*sin_i_half_7*cos_i_half_9 - 194906.25*sin_i_half_5*cos_i_half_11 + 32484.375*sin_i_half_3*cos_i_half_13;
    result_by_lmp.set(c_Key3(8, 3, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 30318.75*(-29.0*sin_i_half_8 + 47.0*sin_i_half_6 - 19.0*sin_i_half_4 - 10.0*cos_i_half_8 + 9.0*cos_i_half_6)*sin_i_half_3*std::sin((1.0/2.0)*obliquity + M_PI_4)*cos_i_half_3*std::cos((1.0/2.0)*obliquity + M_PI_4);
    result_by_lmp.set(c_Key3(8, 3, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 1299.375*(176.0*sin_i_half_10 - 325.0*sin_i_half_8 + 150.0*sin_i_half_6 - 300.0*sin_i_half_4*cos_i_half_6 - 252.0*cos_i_half_10 + 210.0*cos_i_half_8)*sin_i_half_5*cos_i_half;
    result_by_lmp.set(c_Key3(8, 3, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 1299.375*(-133.0*sin_i_half_8 + 242.0*sin_i_half_6 - 110.0*sin_i_half_4 - 231.0*cos_i_half_8 + 165.0*cos_i_half_6)*sin_i_half_7*cos_i_half;
    result_by_lmp.set(c_Key3(8, 3, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = -2815.3125*(12.0*cos_i_2 + 9.0*cos_i + 1.0)*sin_i_half_9*cos_i_half_3;
    result_by_lmp.set(c_Key3(8, 3, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = -16891.875*sin_i_half_11*cos_i_half_5;
    result_by_lmp.set(c_Key3(8, 3, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 3), result_by_p);
    result_by_p.clear();

    // l , m = (8, 4).
    // p = 0
    tmp_double = 84459.375*sin_i_half_4*cos_i_half_12;
    result_by_lmp.set(c_Key3(8, 4, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 131.9677734375*std::pow(cos_i + 1.0, 5)*(-39.0*cos_i + 20.0*cos_i_double - 5.0*cos_i_triple + 24.0);
    result_by_lmp.set(c_Key3(8, 4, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 1299.375*(1771.0*sin_i_half_8 - 1672.0*sin_i_half_6 + 396.0*sin_i_half_4 + 49.0*cos_i_half_8 - 48.0*cos_i_half_6)*cos_i_half_8;
    result_by_lmp.set(c_Key3(8, 4, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 12993.75*(321.0*sin_i_half_8 - 414.0*sin_i_half_6 + 135.0*sin_i_half_4 + 43.0*cos_i_half_8 - 40.0*cos_i_half_6)*sin_i_half_2*cos_i_half_6;
    result_by_lmp.set(c_Key3(8, 4, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 45478.125*(93.0*sin_i_half_8 - 144.0*sin_i_half_6 + 56.0*sin_i_half_4 + 37.0*cos_i_half_8 - 32.0*cos_i_half_6)*sin_i_half_4*cos_i_half_4;
    result_by_lmp.set(c_Key3(8, 4, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 12993.75*(-178.0*sin_i_half_10 + 488.0*sin_i_half_8 - 445.0*sin_i_half_6 + 135.0*sin_i_half_4 + 186.0*cos_i_half_10 - 144.0*cos_i_half_8)*sin_i_half_6;
    result_by_lmp.set(c_Key3(8, 4, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 1299.375*(445.0*sin_i_half_8 - 840.0*sin_i_half_6 + 396.0*sin_i_half_4 + 1375.0*cos_i_half_8 - 880.0*cos_i_half_6)*sin_i_half_8;
    result_by_lmp.set(c_Key3(8, 4, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 33783.75*(-20.0*sin_i_half_6 + 50.0*sin_i_half_4 - 41.0*sin_i_half_2 + 11.0)*sin_i_half_10;
    result_by_lmp.set(c_Key3(8, 4, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 84459.375*sin_i_half_12*cos_i_half_4;
    result_by_lmp.set(c_Key3(8, 4, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 4), result_by_p);
    result_by_p.clear();

    // l , m = (8, 5).
    // p = 0
    tmp_double = -337837.5*sin_i_half_3*cos_i_half_13;
    result_by_lmp.set(c_Key3(8, 5, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -263.935546875*std::pow(cos_i + 1.0, 7)*sin_i - 1756755.0*sin_i_half_5*cos_i_half_11 + 878377.5*sin_i_half_3*cos_i_half_13;
    result_by_lmp.set(c_Key3(8, 5, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 263.935546875*std::pow(cos_i + 1.0, 7)*sin_i - 3716212.5*sin_i_half_7*cos_i_half_9 + 4459455.0*sin_i_half_5*cos_i_half_11 - 1216215.0*sin_i_half_3*cos_i_half_13;
    result_by_lmp.set(c_Key3(8, 5, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 337837.5*(-54.0*sin_i_half_6 + 57.0*sin_i_half_4 - 15.0*sin_i_half_2 + 2.0*cos_i_half_6)*sin_i_half_3*cos_i_half_7;
    result_by_lmp.set(c_Key3(8, 5, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 18475.48828125*(4.0 - 5.0*sin_i_2)*sin_i_5*cos_i;
    result_by_lmp.set(c_Key3(8, 5, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 337837.5*(-44.0*sin_i_half_6 + 69.0*sin_i_half_4 - 27.0*sin_i_half_2 + 12.0*cos_i_half_6)*sin_i_half_7*cos_i_half_3;
    result_by_lmp.set(c_Key3(8, 5, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 67567.5*(-85.0*sin_i_half_6 + 150.0*sin_i_half_4 - 66.0*sin_i_half_2 + 55.0*cos_i_half_6)*sin_i_half_9*cos_i_half;
    result_by_lmp.set(c_Key3(8, 5, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 67567.5*(-10.0*sin_i_2 + 12.5*cos_i + 13.5)*sin_i_half_11*cos_i_half;
    result_by_lmp.set(c_Key3(8, 5, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 337837.5*sin_i_half_13*cos_i_half_3;
    result_by_lmp.set(c_Key3(8, 5, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 5), result_by_p);
    result_by_p.clear();

    // l , m = (8, 6).
    // p = 0
    tmp_double = -1013512.5*sin_i_half_2*cos_i_half_14;
    result_by_lmp.set(c_Key3(8, 6, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 1055.7421875*std::pow(cos_i + 1.0, 6)*(30.0*sin_i_2 + 45.0*cos_i - 46.0);
    result_by_lmp.set(c_Key3(8, 6, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 1418917.5*(5.0*sin_i_2 + 5.0*cos_i - 6.0)*sin_i_half_2*cos_i_half_10;
    result_by_lmp.set(c_Key3(8, 6, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 110852.9296875*(1.0 - 2.0*cos_i)*std::pow(cos_i + 1.0, 2)*sin_i_4*cos_i;
    result_by_lmp.set(c_Key3(8, 6, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 18475.48828125*(15.0*sin_i_2 - 14.0)*sin_i_6;
    result_by_lmp.set(c_Key3(8, 6, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = -110852.9296875*std::pow(cos_i - 1.0, 2)*(2.0*cos_i + 1.0)*sin_i_4*cos_i;
    result_by_lmp.set(c_Key3(8, 6, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 1418917.5*(20.0*sin_i_half_6 - 50.0*sin_i_half_4 + 41.0*sin_i_half_2 - 11.0)*sin_i_half_10;
    result_by_lmp.set(c_Key3(8, 6, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = (2027025.0*sin_i_2 - 3040537.5*cos_i - 3108105.0)*sin_i_half_12;
    result_by_lmp.set(c_Key3(8, 6, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = -1013512.5*sin_i_half_14*cos_i_half_2;
    result_by_lmp.set(c_Key3(8, 6, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 6), result_by_p);
    result_by_p.clear();

    // l , m = (8, 7).
    // p = 0
    tmp_double = 7918.06640625*std::pow(cos_i + 1.0, 7)*sin_i;
    result_by_lmp.set(c_Key3(8, 7, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -7918.06640625*std::pow(cos_i + 1.0, 7)*sin_i + 14189175.0*sin_i_half_3*cos_i_half_13;
    result_by_lmp.set(c_Key3(8, 7, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = (14189175.0 - 28378350.0*cos_i)*sin_i_half_3*cos_i_half_11;
    result_by_lmp.set(c_Key3(8, 7, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = (14189175.0 - 56756700.0*cos_i)*sin_i_half_5*cos_i_half_9;
    result_by_lmp.set(c_Key3(8, 7, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -70945875.0*sin_i_half_7*cos_i_half_7*cos_i;
    result_by_lmp.set(c_Key3(8, 7, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = (-56756700.0*cos_i - 14189175.0)*sin_i_half_9*cos_i_half_5;
    result_by_lmp.set(c_Key3(8, 7, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = (-28378350.0*cos_i - 14189175.0)*sin_i_half_11*cos_i_half_3;
    result_by_lmp.set(c_Key3(8, 7, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = -15836.1328125*std::pow(cos_i - 1.0, 6)*(4.0*cos_i + 3.0)*sin_i;
    result_by_lmp.set(c_Key3(8, 7, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = -2027025.0*sin_i_half_15*cos_i_half;
    result_by_lmp.set(c_Key3(8, 7, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 7), result_by_p);
    result_by_p.clear();

    // l , m = (8, 8).
    // p = 0
    tmp_double = 2027025.0*cos_i_half_16;
    result_by_lmp.set(c_Key3(8, 8, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 16216200.0*sin_i_half_2*cos_i_half_14;
    result_by_lmp.set(c_Key3(8, 8, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 56756700.0*sin_i_half_4*cos_i_half_12;
    result_by_lmp.set(c_Key3(8, 8, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 113513400.0*sin_i_half_6*cos_i_half_10;
    result_by_lmp.set(c_Key3(8, 8, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 141891750.0*sin_i_half_8*cos_i_half_8;
    result_by_lmp.set(c_Key3(8, 8, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 113513400.0*sin_i_half_10*cos_i_half_6;
    result_by_lmp.set(c_Key3(8, 8, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 56756700.0*sin_i_half_12*cos_i_half_4;
    result_by_lmp.set(c_Key3(8, 8, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 16216200.0*sin_i_half_14*cos_i_half_2;
    result_by_lmp.set(c_Key3(8, 8, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 2027025.0*sin_i_half_16;
    result_by_lmp.set(c_Key3(8, 8, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 8), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l8_2(double obliquity)
{
    // Inclination Functions Calculated for l = 8.
    // Functions are taylor expanded and truncated at obliquity^2.

    //  Functions Calculated for l = 8.

    c_IntMap<c_Key3, double> result_by_lmp(13);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(9);
    // Optimizations

    c_IntMap<c_Key1, double> result_by_p(8);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (8, 0).
    // p = 4
    tmp_double = 0.2734375;
    result_by_lmp.set(c_Key3(8, 0, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 0), result_by_p);
    result_by_p.clear();

    // l , m = (8, 1).
    // p = 3
    tmp_double = -9.84375*obliquity;
    result_by_lmp.set(c_Key3(8, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 9.84375*obliquity;
    result_by_lmp.set(c_Key3(8, 1, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 1), result_by_p);
    result_by_p.clear();

    // l , m = (8, 2).
    // p = 3
    tmp_double = -19.6875;
    result_by_lmp.set(c_Key3(8, 2, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 2), result_by_p);
    result_by_p.clear();

    // l , m = (8, 3).
    // p = 2
    tmp_double = 649.6875*obliquity;
    result_by_lmp.set(c_Key3(8, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -649.6875*obliquity;
    result_by_lmp.set(c_Key3(8, 3, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 3), result_by_p);
    result_by_p.clear();

    // l , m = (8, 4).
    // p = 2
    tmp_double = 1299.375;
    result_by_lmp.set(c_Key3(8, 4, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 4), result_by_p);
    result_by_p.clear();

    // l , m = (8, 5).
    // p = 1
    tmp_double = -33783.75*obliquity;
    result_by_lmp.set(c_Key3(8, 5, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 33783.75*obliquity;
    result_by_lmp.set(c_Key3(8, 5, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 5), result_by_p);
    result_by_p.clear();

    // l , m = (8, 6).
    // p = 1
    tmp_double = -67567.5;
    result_by_lmp.set(c_Key3(8, 6, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 6), result_by_p);
    result_by_p.clear();

    // l , m = (8, 7).
    // p = 0
    tmp_double = 1013512.5*obliquity;
    result_by_lmp.set(c_Key3(8, 7, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -1013512.5*obliquity;
    result_by_lmp.set(c_Key3(8, 7, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 7), result_by_p);
    result_by_p.clear();

    // l , m = (8, 8).
    // p = 0
    tmp_double = 2027025.0;
    result_by_lmp.set(c_Key3(8, 8, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 8), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l8_4(double obliquity)
{
    // Inclination Functions Calculated for l = 8.
    // Functions are taylor expanded and truncated at obliquity^4.

    //  Functions Calculated for l = 8.

    c_IntMap<c_Key3, double> result_by_lmp(29);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(9);
    // Optimizations
    double obliquity_2 = obliquity * obliquity;
    double obliquity_3 = obliquity_2 * obliquity;

    c_IntMap<c_Key1, double> result_by_p(8);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (8, 0).
    // p = 3
    tmp_double = 2.4609375*obliquity_2;
    result_by_lmp.set(c_Key3(8, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 0.2734375 - 4.921875*obliquity_2;
    result_by_lmp.set(c_Key3(8, 0, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 2.4609375*obliquity_2;
    result_by_lmp.set(c_Key3(8, 0, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 0), result_by_p);
    result_by_p.clear();

    // l , m = (8, 1).
    // p = 2
    tmp_double = -27.0703125*obliquity_3;
    result_by_lmp.set(c_Key3(8, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 85.3125*obliquity_3 - 9.84375*obliquity;
    result_by_lmp.set(c_Key3(8, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -87.7734375*obliquity_3 + 9.84375*obliquity;
    result_by_lmp.set(c_Key3(8, 1, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 29.53125*obliquity_3;
    result_by_lmp.set(c_Key3(8, 1, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 1), result_by_p);
    result_by_p.clear();

    // l , m = (8, 2).
    // p = 2
    tmp_double = -162.421875*obliquity_2;
    result_by_lmp.set(c_Key3(8, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 334.6875*obliquity_2 - 19.6875;
    result_by_lmp.set(c_Key3(8, 2, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -172.265625*obliquity_2;
    result_by_lmp.set(c_Key3(8, 2, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 2), result_by_p);
    result_by_p.clear();

    // l , m = (8, 3).
    // p = 1
    tmp_double = 1407.65625*obliquity_3;
    result_by_lmp.set(c_Key3(8, 3, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -4818.515625*obliquity_3 + 649.6875*obliquity;
    result_by_lmp.set(c_Key3(8, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 5305.78125*obliquity_3 - 649.6875*obliquity;
    result_by_lmp.set(c_Key3(8, 3, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -1894.921875*obliquity_3;
    result_by_lmp.set(c_Key3(8, 3, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 3), result_by_p);
    result_by_p.clear();

    // l , m = (8, 4).
    // p = 1
    tmp_double = 8445.9375*obliquity_2;
    result_by_lmp.set(c_Key3(8, 4, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 1299.375 - 18191.25*obliquity_2;
    result_by_lmp.set(c_Key3(8, 4, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 9745.3125*obliquity_2;
    result_by_lmp.set(c_Key3(8, 4, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 4), result_by_p);
    result_by_p.clear();

    // l , m = (8, 5).
    // p = 0
    tmp_double = -42229.6875*obliquity_3;
    result_by_lmp.set(c_Key3(8, 5, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 174549.375*obliquity_3 - 33783.75*obliquity;
    result_by_lmp.set(c_Key3(8, 5, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -216779.0625*obliquity_3 + 33783.75*obliquity;
    result_by_lmp.set(c_Key3(8, 5, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 84459.375*obliquity_3;
    result_by_lmp.set(c_Key3(8, 5, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 5), result_by_p);
    result_by_p.clear();

    // l , m = (8, 6).
    // p = 0
    tmp_double = -253378.125*obliquity_2;
    result_by_lmp.set(c_Key3(8, 6, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 608107.5*obliquity_2 - 67567.5;
    result_by_lmp.set(c_Key3(8, 6, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -354729.375*obliquity_2;
    result_by_lmp.set(c_Key3(8, 6, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 6), result_by_p);
    result_by_p.clear();

    // l , m = (8, 7).
    // p = 0
    tmp_double = -1942565.625*obliquity_3 + 1013512.5*obliquity;
    result_by_lmp.set(c_Key3(8, 7, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 3716212.5*obliquity_3 - 1013512.5*obliquity;
    result_by_lmp.set(c_Key3(8, 7, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -1773646.875*obliquity_3;
    result_by_lmp.set(c_Key3(8, 7, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 7), result_by_p);
    result_by_p.clear();

    // l , m = (8, 8).
    // p = 0
    tmp_double = 2027025.0 - 4054050.0*obliquity_2;
    result_by_lmp.set(c_Key3(8, 8, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 4054050.0*obliquity_2;
    result_by_lmp.set(c_Key3(8, 8, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 8), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l8_off(double obliquity)
{
    // Inclination Functions Calculated for l = 8.
    // Functions assume `obliquity=0.0`. However some modes are still activated (non-zero `F(l,m,p)`).

    //  Functions Calculated for l = 8.

    c_IntMap<c_Key3, double> result_by_lmp(5);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(5);
    // Optimizations

    c_IntMap<c_Key1, double> result_by_p(8);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (8, 0).
    // p = 4
    tmp_double = 0.2734375;
    result_by_lmp.set(c_Key3(8, 0, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 0), result_by_p);
    result_by_p.clear();

    // l , m = (8, 2).
    // p = 3
    tmp_double = -19.6875;
    result_by_lmp.set(c_Key3(8, 2, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 2), result_by_p);
    result_by_p.clear();

    // l , m = (8, 4).
    // p = 2
    tmp_double = 1299.375;
    result_by_lmp.set(c_Key3(8, 4, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 4), result_by_p);
    result_by_p.clear();

    // l , m = (8, 6).
    // p = 1
    tmp_double = -67567.5;
    result_by_lmp.set(c_Key3(8, 6, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 6), result_by_p);
    result_by_p.clear();

    // l , m = (8, 8).
    // p = 0
    tmp_double = 2027025.0;
    result_by_lmp.set(c_Key3(8, 8, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(8, 8), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}
