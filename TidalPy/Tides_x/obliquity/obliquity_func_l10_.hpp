#pragma once

#include "obliquity_common_.hpp"

ObliquityFuncOutput c_obliquity_function_l10_gen(double obliquity)
{
    // Inclination Functions Calculated for l = 10.
    // Functions are fully general with no assumptions on the magnitude of obliquity.

    //  Functions Calculated for l = 10.

    c_IntMap<c_Key3, double> result_by_lmp(121);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(11);
    // Optimizations
    double cos_i_half = std::cos((1.0/2.0)*obliquity);
    double sin_i_half = std::sin((1.0/2.0)*obliquity);
    double cos_i_half_2 = cos_i_half * cos_i_half;
    double cos_i_half_4 = cos_i_half_2 * cos_i_half_2;
    double cos_i_half_5 = cos_i_half_4 * cos_i_half;
    double cos_i_half_10 = cos_i_half_5 * cos_i_half_5;
    double sin_i_half_2 = sin_i_half * sin_i_half;
    double sin_i_half_4 = sin_i_half_2 * sin_i_half_2;
    double sin_i_half_5 = sin_i_half_4 * sin_i_half;
    double sin_i_half_10 = sin_i_half_5 * sin_i_half_5;
    double sin_i = std::sin(obliquity);
    double sin_i_2 = sin_i * sin_i;
    double sin_i_4 = sin_i_2 * sin_i_2;
    double sin_i_8 = sin_i_4 * sin_i_4;
    double cos_i_half_3 = cos_i_half_2 * cos_i_half;
    double cos_i_half_6 = cos_i_half_3 * cos_i_half_3;
    double cos_i_half_8 = cos_i_half_4 * cos_i_half_4;
    double sin_i_half_3 = sin_i_half_2 * sin_i_half;
    double sin_i_half_6 = sin_i_half_3 * sin_i_half_3;
    double sin_i_half_8 = sin_i_half_4 * sin_i_half_4;
    double cos_i_half_12 = cos_i_half_6 * cos_i_half_6;
    double sin_i_half_12 = sin_i_half_6 * sin_i_half_6;
    double cos_i_half_7 = cos_i_half_6 * cos_i_half;
    double cos_i_half_14 = cos_i_half_7 * cos_i_half_7;
    double cos_i_half_16 = cos_i_half_8 * cos_i_half_8;
    double sin_i_half_7 = sin_i_half_6 * sin_i_half;
    double sin_i_half_14 = sin_i_half_7 * sin_i_half_7;
    double sin_i_half_16 = sin_i_half_8 * sin_i_half_8;
    double cos_i_half_11 = cos_i_half_10 * cos_i_half;
    double sin_i_half_9 = sin_i_half_8 * sin_i_half;
    double cos_i_half_9 = cos_i_half_8 * cos_i_half;
    double cos_i = std::cos(obliquity);
    double cos_i_half_13 = cos_i_half_12 * cos_i_half;
    double cos_i_half_15 = cos_i_half_14 * cos_i_half;
    double sin_i_half_11 = sin_i_half_10 * sin_i_half;
    double sin_i_half_13 = sin_i_half_12 * sin_i_half;
    double sin_i_half_15 = sin_i_half_14 * sin_i_half;
    double sin_i_3 = sin_i_2 * sin_i;
    double sin_i_6 = sin_i_3 * sin_i_3;
    double sin_i_5 = sin_i_4 * sin_i;
    double cos_i_double = std::cos(2*obliquity);
    double cos_i_triple = std::cos(3*obliquity);
    double cos_i_2 = cos_i * cos_i;

    c_IntMap<c_Key1, double> result_by_p(10);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (10, 0).
    // p = 0
    tmp_double = -180.42578125*sin_i_half_10*cos_i_half_10;
    result_by_lmp.set(c_Key3(10, 0, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 0.09273529052734375*(19.0*sin_i_2 - 18.0)*sin_i_8;
    result_by_lmp.set(c_Key3(10, 0, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 25.13671875*(-229.0*sin_i_half_8 + 350.0*sin_i_half_6 - 135.0*sin_i_half_4 - 94.0*cos_i_half_8 + 80.0*cos_i_half_6)*sin_i_half_6*cos_i_half_6;
    result_by_lmp.set(c_Key3(10, 0, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 16.7578125*(-406.0*sin_i_half_12 + 714.0*sin_i_half_10 - 315.0*sin_i_half_8 + 480.0*sin_i_half_6*cos_i_half_6 - 315.0*sin_i_half_4*cos_i_half_8 - 91.0*cos_i_half_12 + 84.0*cos_i_half_10)*sin_i_half_4*cos_i_half_4;
    result_by_lmp.set(c_Key3(10, 0, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 4.51171875*(713.0*std::pow(sin_i_half, 18) - 2053.0*sin_i_half_16 + 1970.0*sin_i_half_14 - 630.0*sin_i_half_12 + 2016.0*sin_i_half_10*cos_i_half_8 - 2940.0*sin_i_half_8*cos_i_half_10 + 2016.0*sin_i_half_6*cos_i_half_12 - 630.0*sin_i_half_4*cos_i_half_14 - 83.0*std::pow(cos_i_half, 18) + 80.0*cos_i_half_16)*sin_i_half_2;
    result_by_lmp.set(c_Key3(10, 0, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = -523.1953125*std::pow(sin_i_half, 20) + 1021.2890625*std::pow(sin_i_half, 18) - 498.33984375*sin_i_half_16 + 3543.75*sin_i_half_14*cos_i_half_6 - 10852.734375*sin_i_half_12*cos_i_half_8 + 15627.9375*sin_i_half_10*cos_i_half_10 - 10852.734375*sin_i_half_8*cos_i_half_12 + 3543.75*sin_i_half_6*cos_i_half_14 - 498.33984375*sin_i_half_4*cos_i_half_16 - 24.85546875*std::pow(cos_i_half, 20) + 24.609375*std::pow(cos_i_half, 18);
    result_by_lmp.set(c_Key3(10, 0, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 4.51171875*(713.0*std::pow(sin_i_half, 18) - 2053.0*sin_i_half_16 + 1970.0*sin_i_half_14 - 630.0*sin_i_half_12 + 2016.0*sin_i_half_10*cos_i_half_8 - 2940.0*sin_i_half_8*cos_i_half_10 + 2016.0*sin_i_half_6*cos_i_half_12 - 630.0*sin_i_half_4*cos_i_half_14 - 83.0*std::pow(cos_i_half, 18) + 80.0*cos_i_half_16)*sin_i_half_2;
    result_by_lmp.set(c_Key3(10, 0, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 16.7578125*(-406.0*sin_i_half_12 + 714.0*sin_i_half_10 - 315.0*sin_i_half_8 + 480.0*sin_i_half_6*cos_i_half_6 - 315.0*sin_i_half_4*cos_i_half_8 - 91.0*cos_i_half_12 + 84.0*cos_i_half_10)*sin_i_half_4*cos_i_half_4;
    result_by_lmp.set(c_Key3(10, 0, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 25.13671875*(-229.0*sin_i_half_8 + 350.0*sin_i_half_6 - 135.0*sin_i_half_4 - 94.0*cos_i_half_8 + 80.0*cos_i_half_6)*sin_i_half_6*cos_i_half_6;
    result_by_lmp.set(c_Key3(10, 0, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = 0.09273529052734375*(19.0*sin_i_2 - 18.0)*sin_i_8;
    result_by_lmp.set(c_Key3(10, 0, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // p = 10
    tmp_double = -180.42578125*sin_i_half_10*cos_i_half_10;
    result_by_lmp.set(c_Key3(10, 0, 10), tmp_double);
    result_by_p.set(c_Key1(10), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 0), result_by_p);
    result_by_p.clear();

    // l , m = (10, 1).
    // p = 0
    tmp_double = 1804.2578125*sin_i_half_9*cos_i_half_11;
    result_by_lmp.set(c_Key3(10, 1, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 94.9609375*(190.0*sin_i_half_4 - 171.0*sin_i_half_2 + 36.0)*sin_i_half_7*cos_i_half_9;
    result_by_lmp.set(c_Key3(10, 1, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 50.2734375*(1265.0*sin_i_half_8 - 1815.0*sin_i_half_6 + 660.0*sin_i_half_4 + 350.0*cos_i_half_8 - 308.0*cos_i_half_6)*sin_i_half_5*cos_i_half_7;
    result_by_lmp.set(c_Key3(10, 1, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 33.515625*(2750.0*sin_i_half_12 - 4653.0*sin_i_half_10 + 1980.0*sin_i_half_8 - 2310.0*sin_i_half_6*cos_i_half_6 + 1155.0*sin_i_half_4*cos_i_half_8 + 245.0*cos_i_half_12 - 231.0*cos_i_half_10)*sin_i_half_3*cos_i_half_5;
    result_by_lmp.set(c_Key3(10, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 0.02643585205078125*std::pow(cos_i + 1.0, 9)*sin_i + 496.2890625*std::pow(sin_i_half, 17)*cos_i_half_3 - 8933.203125*sin_i_half_15*cos_i_half_5 + 50025.9375*sin_i_half_13*cos_i_half_7 - 116727.1875*sin_i_half_11*cos_i_half_9 + 125064.84375*sin_i_half_9*cos_i_half_11 - 62532.421875*sin_i_half_7*cos_i_half_13 + 13896.09375*sin_i_half_5*cos_i_half_15 - 1191.09375*sin_i_half_3*std::pow(cos_i_half, 17);
    result_by_lmp.set(c_Key3(10, 1, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = -0.02643585205078125*std::pow(cos_i + 1.0, 9)*sin_i + 27.0703125*std::pow(sin_i_half, 19)*cos_i_half - 1218.1640625*std::pow(sin_i_half, 17)*cos_i_half_3 + 14617.96875*sin_i_half_15*cos_i_half_5 - 68217.1875*sin_i_half_13*cos_i_half_7 + 143256.09375*sin_i_half_11*cos_i_half_9 - 143256.09375*sin_i_half_9*cos_i_half_11 + 68217.1875*sin_i_half_7*cos_i_half_13 - 14617.96875*sin_i_half_5*cos_i_half_15 + 1218.1640625*sin_i_half_3*std::pow(cos_i_half, 17);
    result_by_lmp.set(c_Key3(10, 1, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 9.0234375*(-1675.0*sin_i_half_16 + 3212.0*sin_i_half_14 - 1540.0*sin_i_half_12 + 6930.0*sin_i_half_10*cos_i_half_6 - 13860.0*sin_i_half_8*cos_i_half_8 + 12936.0*sin_i_half_6*cos_i_half_10 - 5544.0*sin_i_half_4*cos_i_half_12 - 1045.0*cos_i_half_16 + 990.0*cos_i_half_14)*sin_i_half_3*cos_i_half;
    result_by_lmp.set(c_Key3(10, 1, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 33.515625*(-1400.0*sin_i_half_12 + 2541.0*sin_i_half_10 - 1155.0*sin_i_half_8 + 2310.0*sin_i_half_6*cos_i_half_6 - 1980.0*sin_i_half_4*cos_i_half_8 - 770.0*cos_i_half_12 + 693.0*cos_i_half_10)*sin_i_half_5*cos_i_half_3;
    result_by_lmp.set(c_Key3(10, 1, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 50.2734375*(-1010.0*sin_i_half_8 + 1628.0*sin_i_half_6 - 660.0*sin_i_half_4 - 605.0*cos_i_half_8 + 495.0*cos_i_half_6)*sin_i_half_7*cos_i_half_5;
    result_by_lmp.set(c_Key3(10, 1, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = 94.9609375*(-190.0*sin_i_half_4 + 209.0*sin_i_half_2 - 55.0)*sin_i_half_9*cos_i_half_7;
    result_by_lmp.set(c_Key3(10, 1, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // p = 10
    tmp_double = -1804.2578125*sin_i_half_11*cos_i_half_9;
    result_by_lmp.set(c_Key3(10, 1, 10), tmp_double);
    result_by_p.set(c_Key1(10), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 1), result_by_p);
    result_by_p.clear();

    // l , m = (10, 2).
    // p = 0
    tmp_double = 16238.3203125*sin_i_half_8*cos_i_half_12;
    result_by_lmp.set(c_Key3(10, 2, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 1.6692352294921875*std::pow(cos_i + 1.0, 2)*(-95.0*sin_i_2 - 38.0*cos_i + 94.0)*sin_i_6;
    result_by_lmp.set(c_Key3(10, 2, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 150.8203125*(4103.0*sin_i_half_8 - 5456.0*sin_i_half_6 + 1848.0*sin_i_half_4 + 742.0*cos_i_half_8 - 672.0*cos_i_half_6)*sin_i_half_4*cos_i_half_8;
    result_by_lmp.set(c_Key3(10, 2, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 201.09375*(5280.0*sin_i_half_12 - 8514.0*sin_i_half_10 + 3465.0*sin_i_half_8 - 3080.0*sin_i_half_6*cos_i_half_6 + 1155.0*sin_i_half_4*cos_i_half_8 + 175.0*cos_i_half_12 - 168.0*cos_i_half_10)*sin_i_half_2*cos_i_half_6;
    result_by_lmp.set(c_Key3(10, 2, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 27.0703125*(32703.0*sin_i_half_16 - 58080.0*sin_i_half_14 + 25872.0*sin_i_half_12 - 44352.0*sin_i_half_10*cos_i_half_6 + 34650.0*sin_i_half_8*cos_i_half_8 - 12320.0*sin_i_half_6*cos_i_half_10 + 1945.0*cos_i_half_16 - 3792.0*cos_i_half_14 + 1848.0*cos_i_half_12)*cos_i_half_4;
    result_by_lmp.set(c_Key3(10, 2, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 487.265625*(-713.0*std::pow(sin_i_half, 18) + 2053.0*sin_i_half_16 - 1970.0*sin_i_half_14 + 630.0*sin_i_half_12 - 2016.0*sin_i_half_10*cos_i_half_8 + 2940.0*sin_i_half_8*cos_i_half_10 - 2016.0*sin_i_half_6*cos_i_half_12 + 630.0*sin_i_half_4*cos_i_half_14 + 83.0*std::pow(cos_i_half, 18) - 80.0*cos_i_half_16)*sin_i_half_2;
    result_by_lmp.set(c_Key3(10, 2, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 27.0703125*(1945.0*sin_i_half_16 - 3792.0*sin_i_half_14 + 1848.0*sin_i_half_12 - 12320.0*sin_i_half_10*cos_i_half_6 + 34650.0*sin_i_half_8*cos_i_half_8 - 44352.0*sin_i_half_6*cos_i_half_10 + 25872.0*sin_i_half_4*cos_i_half_12 + 6831.0*cos_i_half_16 - 6336.0*cos_i_half_14)*sin_i_half_4;
    result_by_lmp.set(c_Key3(10, 2, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 201.09375*(-1330.0*sin_i_half_14 + 3808.0*sin_i_half_12 - 3633.0*sin_i_half_10 + 1155.0*sin_i_half_8 - 3080.0*sin_i_half_6*cos_i_half_8 + 3465.0*sin_i_half_4*cos_i_half_10 + 1815.0*cos_i_half_14 - 1584.0*cos_i_half_12)*sin_i_half_6;
    result_by_lmp.set(c_Key3(10, 2, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 150.8203125*(2590.0*sin_i_half_8 - 4368.0*sin_i_half_6 + 1848.0*sin_i_half_4 + 2255.0*cos_i_half_8 - 1760.0*cos_i_half_6)*sin_i_half_8*cos_i_half_4;
    result_by_lmp.set(c_Key3(10, 2, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = 1.6692352294921875*std::pow(cos_i - 1.0, 2)*(-95.0*sin_i_2 + 38.0*cos_i + 94.0)*sin_i_6;
    result_by_lmp.set(c_Key3(10, 2, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // p = 10
    tmp_double = 16238.3203125*sin_i_half_12*cos_i_half_8;
    result_by_lmp.set(c_Key3(10, 2, 10), tmp_double);
    result_by_p.set(c_Key1(10), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 2), result_by_p);
    result_by_p.clear();

    // l , m = (10, 3).
    // p = 0
    tmp_double = -129906.5625*sin_i_half_7*cos_i_half_13;
    result_by_lmp.set(c_Key3(10, 3, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 6837.1875*(-190.0*sin_i_half_4 + 133.0*sin_i_half_2 - 21.0)*sin_i_half_5*cos_i_half_11;
    result_by_lmp.set(c_Key3(10, 3, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 1206.5625*(-4355.0*sin_i_half_8 + 5278.0*sin_i_half_6 - 1638.0*sin_i_half_4 - 490.0*cos_i_half_8 + 455.0*cos_i_half_6)*sin_i_half_3*cos_i_half_9;
    result_by_lmp.set(c_Key3(10, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -2.74932861328125*std::pow(cos_i + 1.0, 9)*sin_i - 690153.75*sin_i_half_13*cos_i_half_7 + 3623307.1875*sin_i_half_11*cos_i_half_9 - 6038845.3125*sin_i_half_9*cos_i_half_11 + 4025896.875*sin_i_half_7*cos_i_half_13 - 1097971.875*sin_i_half_5*cos_i_half_15 + 109797.1875*sin_i_half_3*std::pow(cos_i_half, 17);
    result_by_lmp.set(c_Key3(10, 3, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 2.74932861328125*std::pow(cos_i + 1.0, 9)*sin_i - 278715.9375*sin_i_half_15*cos_i_half_5 + 2601348.75*sin_i_half_13*cos_i_half_7 - 7804046.25*sin_i_half_11*cos_i_half_9 + 9755057.8125*sin_i_half_9*cos_i_half_11 - 5419476.5625*sin_i_half_7*cos_i_half_13 + 1300674.375*sin_i_half_5*cos_i_half_15 - 118243.125*sin_i_half_3*std::pow(cos_i_half, 17);
    result_by_lmp.set(c_Key3(10, 3, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 50675.625*(191.0*sin_i_half_12 - 345.0*sin_i_half_10 + 156.0*sin_i_half_8 - 264.0*sin_i_half_6*cos_i_half_6 + 156.0*sin_i_half_4*cos_i_half_8 + 35.0*cos_i_half_12 - 33.0*cos_i_half_10)*sin_i_half_3*std::sin((1.0/2.0)*obliquity + M_PI_4)*cos_i_half_3*std::cos((1.0/2.0)*obliquity + M_PI_4);
    result_by_lmp.set(c_Key3(10, 3, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 2815.3125*(-505.0*sin_i_half_14 + 966.0*sin_i_half_12 - 462.0*sin_i_half_10 + 1925.0*sin_i_half_8*cos_i_half_6 - 3465.0*sin_i_half_6*cos_i_half_8 + 2772.0*sin_i_half_4*cos_i_half_10 + 1023.0*cos_i_half_14 - 924.0*cos_i_half_12)*sin_i_half_5*cos_i_half;
    result_by_lmp.set(c_Key3(10, 3, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 402.1875*(3010.0*sin_i_half_12 - 5733.0*sin_i_half_10 + 2730.0*sin_i_half_8 - 10010.0*sin_i_half_6*cos_i_half_6 + 15015.0*sin_i_half_4*cos_i_half_8 + 10725.0*cos_i_half_12 - 9009.0*cos_i_half_10)*sin_i_half_7*cos_i_half;
    result_by_lmp.set(c_Key3(10, 3, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 1206.5625*(2128.0*sin_i_half_8 - 3731.0*sin_i_half_6 + 1638.0*sin_i_half_4 + 2717.0*cos_i_half_8 - 2002.0*cos_i_half_6)*sin_i_half_9*cos_i_half_3;
    result_by_lmp.set(c_Key3(10, 3, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = 6837.1875*(190.0*sin_i_half_4 - 247.0*sin_i_half_2 + 78.0)*sin_i_half_11*cos_i_half_5;
    result_by_lmp.set(c_Key3(10, 3, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // p = 10
    tmp_double = 129906.5625*sin_i_half_13*cos_i_half_7;
    result_by_lmp.set(c_Key3(10, 3, 10), tmp_double);
    result_by_p.set(c_Key1(10), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 3), result_by_p);
    result_by_p.clear();

    // l , m = (10, 4).
    // p = 0
    tmp_double = -909345.9375*sin_i_half_6*cos_i_half_14;
    result_by_lmp.set(c_Key3(10, 4, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 93.4771728515625*std::pow(cos_i + 1.0, 4)*(95.0*sin_i_2 + 76.0*cos_i - 106.0)*sin_i_4;
    result_by_lmp.set(c_Key3(10, 4, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 8445.9375*(-4550.0*sin_i_half_8 + 4914.0*sin_i_half_6 - 1365.0*sin_i_half_4 - 295.0*cos_i_half_8 + 280.0*cos_i_half_6)*sin_i_half_2*cos_i_half_10;
    result_by_lmp.set(c_Key3(10, 4, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 2815.3125*(-30030.0*sin_i_half_12 + 42042.0*sin_i_half_10 - 15015.0*sin_i_half_8 + 7280.0*sin_i_half_6*cos_i_half_6 - 1365.0*sin_i_half_4*cos_i_half_8 - 85.0*cos_i_half_12 + 84.0*cos_i_half_10)*cos_i_half_8;
    result_by_lmp.set(c_Key3(10, 4, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 19707.1875*(-5280.0*sin_i_half_12 + 8514.0*sin_i_half_10 - 3465.0*sin_i_half_8 + 3080.0*sin_i_half_6*cos_i_half_6 - 1155.0*sin_i_half_4*cos_i_half_8 - 175.0*cos_i_half_12 + 168.0*cos_i_half_10)*sin_i_half_2*cos_i_half_6;
    result_by_lmp.set(c_Key3(10, 4, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 177364.6875*(-406.0*sin_i_half_12 + 714.0*sin_i_half_10 - 315.0*sin_i_half_8 + 480.0*sin_i_half_6*cos_i_half_6 - 315.0*sin_i_half_4*cos_i_half_8 - 91.0*cos_i_half_12 + 84.0*cos_i_half_10)*sin_i_half_4*cos_i_half_4;
    result_by_lmp.set(c_Key3(10, 4, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 19707.1875*(1330.0*sin_i_half_14 - 3808.0*sin_i_half_12 + 3633.0*sin_i_half_10 - 1155.0*sin_i_half_8 + 3080.0*sin_i_half_6*cos_i_half_8 - 3465.0*sin_i_half_4*cos_i_half_10 - 1815.0*cos_i_half_14 + 1584.0*cos_i_half_12)*sin_i_half_6;
    result_by_lmp.set(c_Key3(10, 4, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 2815.3125*(-1450.0*sin_i_half_12 + 2814.0*sin_i_half_10 - 1365.0*sin_i_half_8 + 7280.0*sin_i_half_6*cos_i_half_6 - 15015.0*sin_i_half_4*cos_i_half_8 - 15015.0*cos_i_half_12 + 12012.0*cos_i_half_10)*sin_i_half_8;
    result_by_lmp.set(c_Key3(10, 4, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 8445.9375*(1660.0*sin_i_half_10 - 4670.0*sin_i_half_8 + 4375.0*sin_i_half_6 - 1365.0*sin_i_half_4 - 3185.0*cos_i_half_10 + 2184.0*cos_i_half_8)*sin_i_half_10;
    result_by_lmp.set(c_Key3(10, 4, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = 93.4771728515625*std::pow(cos_i - 1.0, 4)*(95.0*sin_i_2 - 76.0*cos_i - 106.0)*sin_i_4;
    result_by_lmp.set(c_Key3(10, 4, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // p = 10
    tmp_double = -909345.9375*sin_i_half_14*cos_i_half_6;
    result_by_lmp.set(c_Key3(10, 4, 10), tmp_double);
    result_by_p.set(c_Key1(10), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 4), result_by_p);
    result_by_p.clear();

    // l , m = (10, 5).
    // p = 0
    tmp_double = 5456075.625*sin_i_half_5*cos_i_half_15;
    result_by_lmp.set(c_Key3(10, 5, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 1435809.375*(38.0*sin_i_half_4 - 19.0*sin_i_half_2 + 2.0)*sin_i_half_3*cos_i_half_13;
    result_by_lmp.set(c_Key3(10, 5, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 247.4395751953125*std::pow(cos_i + 1.0, 9)*sin_i + 69172228.125*sin_i_half_9*cos_i_half_11 - 115287046.875*sin_i_half_7*cos_i_half_13 + 53209406.25*sin_i_half_5*cos_i_half_15 - 7601343.75*sin_i_half_3*std::pow(cos_i_half, 17);
    result_by_lmp.set(c_Key3(10, 5, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -247.4395751953125*std::pow(cos_i + 1.0, 9)*sin_i + 84543834.375*sin_i_half_11*cos_i_half_9 - 253631503.125*sin_i_half_9*cos_i_half_11 + 230574093.75*sin_i_half_7*cos_i_half_13 - 76858031.25*sin_i_half_5*cos_i_half_15 + 8868234.375*sin_i_half_3*std::pow(cos_i_half, 17);
    result_by_lmp.set(c_Key3(10, 5, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 591215.625*(1364.0*sin_i_half_10 - 2035.0*sin_i_half_8 + 770.0*sin_i_half_6 - 462.0*sin_i_half_4*cos_i_half_6 - 112.0*cos_i_half_10 + 105.0*cos_i_half_8)*sin_i_half_3*cos_i_half_7;
    result_by_lmp.set(c_Key3(10, 5, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 4156.98486328125*(-323.0*sin_i_4 + 476.0*sin_i_2 - 168.0)*sin_i_5*cos_i;
    result_by_lmp.set(c_Key3(10, 5, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 591215.625*(574.0*sin_i_half_10 - 1029.0*sin_i_half_8 + 462.0*sin_i_half_6 - 770.0*sin_i_half_4*cos_i_half_6 - 594.0*cos_i_half_10 + 495.0*cos_i_half_8)*sin_i_half_7*cos_i_half_3;
    result_by_lmp.set(c_Key3(10, 5, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 84459.375*(1018.0*sin_i_half_10 - 1925.0*sin_i_half_8 + 910.0*sin_i_half_6 - 2730.0*sin_i_half_4*cos_i_half_6 - 4004.0*cos_i_half_10 + 3003.0*cos_i_half_8)*sin_i_half_9*cos_i_half;
    result_by_lmp.set(c_Key3(10, 5, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 253378.125*(-241.0*sin_i_half_8 + 450.0*sin_i_half_6 - 210.0*sin_i_half_4 - 728.0*cos_i_half_8 + 455.0*cos_i_half_6)*sin_i_half_11*cos_i_half;
    result_by_lmp.set(c_Key3(10, 5, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = 1435809.375*(-38.0*sin_i_half_4 + 57.0*sin_i_half_2 - 21.0)*sin_i_half_13*cos_i_half_3;
    result_by_lmp.set(c_Key3(10, 5, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // p = 10
    tmp_double = -5456075.625*sin_i_half_15*cos_i_half_5;
    result_by_lmp.set(c_Key3(10, 5, 10), tmp_double);
    result_by_p.set(c_Key1(10), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 5), result_by_p);
    result_by_p.clear();

    // l , m = (10, 6).
    // p = 0
    tmp_double = 27280378.125*sin_i_half_4*cos_i_half_16;
    result_by_lmp.set(c_Key3(10, 6, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 701.07879638671875*std::pow(cos_i + 1.0, 7)*(-865.0*cos_i + 418.0*cos_i_double - 95.0*cos_i_triple + 542.0);
    result_by_lmp.set(c_Key3(10, 6, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 253378.125*(4780.0*sin_i_half_8 - 3680.0*sin_i_half_6 + 720.0*sin_i_half_4 + 65.0*cos_i_half_8 - 64.0*cos_i_half_6)*cos_i_half_12;
    result_by_lmp.set(c_Key3(10, 6, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = (3074321250.0*sin_i_half_8 - 3320266950.0*sin_i_half_6 + 922296375.0*sin_i_half_4 + 199324125.0*cos_i_half_8 - 189189000.0*cos_i_half_6)*sin_i_half_2*cos_i_half_10;
    result_by_lmp.set(c_Key3(10, 6, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 1182431.25*(4103.0*sin_i_half_8 - 5456.0*sin_i_half_6 + 1848.0*sin_i_half_4 + 742.0*cos_i_half_8 - 672.0*cos_i_half_6)*sin_i_half_4*cos_i_half_8;
    result_by_lmp.set(c_Key3(10, 6, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 21283762.5*(229.0*sin_i_half_8 - 350.0*sin_i_half_6 + 135.0*sin_i_half_4 + 94.0*cos_i_half_8 - 80.0*cos_i_half_6)*sin_i_half_6*cos_i_half_6;
    result_by_lmp.set(c_Key3(10, 6, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 1182431.25*(2590.0*sin_i_half_8 - 4368.0*sin_i_half_6 + 1848.0*sin_i_half_4 + 2255.0*cos_i_half_8 - 1760.0*cos_i_half_6)*sin_i_half_8*cos_i_half_4;
    result_by_lmp.set(c_Key3(10, 6, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = (-1121620500.0*sin_i_half_10 + 3155402250.0*sin_i_half_8 - 2956078125.0*sin_i_half_6 + 922296375.0*sin_i_half_4 + 2152024875.0*cos_i_half_10 - 1475674200.0*cos_i_half_8)*sin_i_half_10;
    result_by_lmp.set(c_Key3(10, 6, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 253378.125*(785.0*sin_i_half_8 - 1504.0*sin_i_half_6 + 720.0*sin_i_half_4 + 4060.0*cos_i_half_8 - 2240.0*cos_i_half_6)*sin_i_half_12;
    result_by_lmp.set(c_Key3(10, 6, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = 2871618.75*(-95.0*sin_i_half_6 + 247.0*sin_i_half_4 - 212.0*sin_i_half_2 + 60.0)*sin_i_half_14;
    result_by_lmp.set(c_Key3(10, 6, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // p = 10
    tmp_double = 27280378.125*sin_i_half_16*cos_i_half_4;
    result_by_lmp.set(c_Key3(10, 6, 10), tmp_double);
    result_by_p.set(c_Key1(10), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 6), result_by_p);
    result_by_p.clear();

    // l , m = (10, 7).
    // p = 0
    tmp_double = -109121512.5*sin_i_half_3*std::pow(cos_i_half, 17);
    result_by_lmp.set(c_Key3(10, 7, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -16825.89111328125*std::pow(cos_i + 1.0, 9)*sin_i - 781080300.0*sin_i_half_5*cos_i_half_15 + 292905112.5*sin_i_half_3*std::pow(cos_i_half, 17);
    result_by_lmp.set(c_Key3(10, 7, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 16825.89111328125*std::pow(cos_i + 1.0, 9)*sin_i - 2412159750.0*sin_i_half_7*cos_i_half_13 + 2067565500.0*sin_i_half_5*cos_i_half_15 - 413513100.0*sin_i_half_3*std::pow(cos_i_half, 17);
    result_by_lmp.set(c_Key3(10, 7, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = (-12864852000.0*sin_i_half_6 + 11095934850.0*sin_i_half_4 - 2412159750.0*sin_i_half_2 + 229729500.0*cos_i_half_6)*sin_i_half_3*cos_i_half_11;
    result_by_lmp.set(c_Key3(10, 7, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = (-21789843075.0*sin_i_half_6 + 24121597500.0*sin_i_half_4 - 6754047300.0*sin_i_half_2 + 1125674550.0*cos_i_half_6)*sin_i_half_5*cos_i_half_9;
    result_by_lmp.set(c_Key3(10, 7, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 1447295850.0*(19.0*sin_i_half_4 - 19.0*sin_i_half_2 + 4.0)*sin_i_half_7*std::sin((1.0/2.0)*obliquity + M_PI_4)*cos_i_half_7*std::cos((1.0/2.0)*obliquity + M_PI_4);
    result_by_lmp.set(c_Key3(10, 7, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = (-18493224750.0*sin_i_half_6 + 27981053100.0*sin_i_half_4 - 10613502900.0*sin_i_half_2 + 4422292875.0*cos_i_half_6)*sin_i_half_9*cos_i_half_5;
    result_by_lmp.set(c_Key3(10, 7, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = (-8913504600.0*sin_i_half_6 + 14955390450.0*sin_i_half_4 - 6271615350.0*sin_i_half_2 + 4181076900.0*cos_i_half_6)*sin_i_half_11*cos_i_half_3;
    result_by_lmp.set(c_Key3(10, 7, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 17229712.5*(-145.0*sin_i_half_6 + 264.0*sin_i_half_4 - 120.0*sin_i_half_2 + 140.0*cos_i_half_6)*sin_i_half_13*cos_i_half;
    result_by_lmp.set(c_Key3(10, 7, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = 5743237.5*(190.0*sin_i_half_4 - 323.0*sin_i_half_2 + 136.0)*sin_i_half_15*cos_i_half;
    result_by_lmp.set(c_Key3(10, 7, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // p = 10
    tmp_double = 109121512.5*std::pow(sin_i_half, 17)*cos_i_half_3;
    result_by_lmp.set(c_Key3(10, 7, 10), tmp_double);
    result_by_p.set(c_Key1(10), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 7), result_by_p);
    result_by_p.clear();

    // l , m = (10, 8).
    // p = 0
    tmp_double = -327364537.5*sin_i_half_2*std::pow(cos_i_half, 18);
    result_by_lmp.set(c_Key3(10, 8, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 33651.7822265625*std::pow(cos_i + 1.0, 8)*(95.0*sin_i_2 + 152.0*cos_i - 154.0);
    result_by_lmp.set(c_Key3(10, 8, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 37858.255004882813*std::pow(cos_i + 1.0, 7)*(865.0*cos_i - 418.0*cos_i_double + 95.0*cos_i_triple - 542.0);
    result_by_lmp.set(c_Key3(10, 8, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = (9820936125.0*sin_i_2 + 7856748900.0*cos_i - 10958097150.0)*sin_i_half_4*cos_i_half_12;
    result_by_lmp.set(c_Key3(10, 8, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 180911981.25*(-95.0*cos_i_2 + 38.0*cos_i + 1.0)*sin_i_half_6*cos_i_half_10;
    result_by_lmp.set(c_Key3(10, 8, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = (20623965862.5*sin_i_2 - 19538493975.0)*sin_i_half_8*cos_i_half_8;
    result_by_lmp.set(c_Key3(10, 8, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 180911981.25*(-95.0*cos_i_2 - 38.0*cos_i + 1.0)*sin_i_half_10*cos_i_half_6;
    result_by_lmp.set(c_Key3(10, 8, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = (9820936125.0*sin_i_2 - 7856748900.0*cos_i - 10958097150.0)*sin_i_half_12*cos_i_half_4;
    result_by_lmp.set(c_Key3(10, 8, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 155067412.5*(95.0*sin_i_half_6 - 247.0*sin_i_half_4 + 212.0*sin_i_half_2 - 60.0)*sin_i_half_14;
    result_by_lmp.set(c_Key3(10, 8, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = (818411343.75*sin_i_2 - 1309458150.0*cos_i - 1326687862.5)*sin_i_half_16;
    result_by_lmp.set(c_Key3(10, 8, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // p = 10
    tmp_double = -327364537.5*std::pow(sin_i_half, 18)*cos_i_half_2;
    result_by_lmp.set(c_Key3(10, 8, 10), tmp_double);
    result_by_p.set(c_Key1(10), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 8), result_by_p);
    result_by_p.clear();

    // l , m = (10, 9).
    // p = 0
    tmp_double = 639383.8623046875*std::pow(cos_i + 1.0, 9)*sin_i;
    result_by_lmp.set(c_Key3(10, 9, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -639383.8623046875*std::pow(cos_i + 1.0, 9)*sin_i + 5892561675.0*sin_i_half_3*std::pow(cos_i_half, 17);
    result_by_lmp.set(c_Key3(10, 9, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 2946280837.5*(3.0 - 5.0*cos_i)*sin_i_half_3*cos_i_half_15;
    result_by_lmp.set(c_Key3(10, 9, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = (15713497800.0 - 39283744500.0*cos_i)*sin_i_half_5*cos_i_half_13;
    result_by_lmp.set(c_Key3(10, 9, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = (13749310575.0 - 68746552875.0*cos_i)*sin_i_half_7*cos_i_half_11;
    result_by_lmp.set(c_Key3(10, 9, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = -82495863450.0*sin_i_half_9*cos_i_half_9*cos_i;
    result_by_lmp.set(c_Key3(10, 9, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = (-68746552875.0*cos_i - 13749310575.0)*sin_i_half_11*cos_i_half_7;
    result_by_lmp.set(c_Key3(10, 9, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = (-39283744500.0*cos_i - 15713497800.0)*sin_i_half_13*cos_i_half_5;
    result_by_lmp.set(c_Key3(10, 9, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = -2946280837.5*(5.0*cos_i + 3.0)*sin_i_half_15*cos_i_half_3;
    result_by_lmp.set(c_Key3(10, 9, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = -654729075.0*(5.0*cos_i + 4.0)*std::pow(sin_i_half, 17)*cos_i_half;
    result_by_lmp.set(c_Key3(10, 9, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // p = 10
    tmp_double = -654729075.0*std::pow(sin_i_half, 19)*cos_i_half;
    result_by_lmp.set(c_Key3(10, 9, 10), tmp_double);
    result_by_p.set(c_Key1(10), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 9), result_by_p);
    result_by_p.clear();

    // l , m = (10, 10).
    // p = 0
    tmp_double = 654729075.0*std::pow(cos_i_half, 20);
    result_by_lmp.set(c_Key3(10, 10, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 6547290750.0*sin_i_half_2*std::pow(cos_i_half, 18);
    result_by_lmp.set(c_Key3(10, 10, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 29462808375.0*sin_i_half_4*cos_i_half_16;
    result_by_lmp.set(c_Key3(10, 10, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 78567489000.0*sin_i_half_6*cos_i_half_14;
    result_by_lmp.set(c_Key3(10, 10, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 137493105750.0*sin_i_half_8*cos_i_half_12;
    result_by_lmp.set(c_Key3(10, 10, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 164991726900.0*sin_i_half_10*cos_i_half_10;
    result_by_lmp.set(c_Key3(10, 10, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = 137493105750.0*sin_i_half_12*cos_i_half_8;
    result_by_lmp.set(c_Key3(10, 10, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // p = 7
    tmp_double = 78567489000.0*sin_i_half_14*cos_i_half_6;
    result_by_lmp.set(c_Key3(10, 10, 7), tmp_double);
    result_by_p.set(c_Key1(7), tmp_double);
    // p = 8
    tmp_double = 29462808375.0*sin_i_half_16*cos_i_half_4;
    result_by_lmp.set(c_Key3(10, 10, 8), tmp_double);
    result_by_p.set(c_Key1(8), tmp_double);
    // p = 9
    tmp_double = 6547290750.0*std::pow(sin_i_half, 18)*cos_i_half_2;
    result_by_lmp.set(c_Key3(10, 10, 9), tmp_double);
    result_by_p.set(c_Key1(9), tmp_double);
    // p = 10
    tmp_double = 654729075.0*std::pow(sin_i_half, 20);
    result_by_lmp.set(c_Key3(10, 10, 10), tmp_double);
    result_by_p.set(c_Key1(10), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 10), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l10_2(double obliquity)
{
    // Inclination Functions Calculated for l = 10.
    // Functions are taylor expanded and truncated at obliquity^2.

    //  Functions Calculated for l = 10.

    c_IntMap<c_Key3, double> result_by_lmp(16);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(11);
    // Optimizations

    c_IntMap<c_Key1, double> result_by_p(10);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (10, 0).
    // p = 5
    tmp_double = -0.24609375;
    result_by_lmp.set(c_Key3(10, 0, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 0), result_by_p);
    result_by_p.clear();

    // l , m = (10, 1).
    // p = 4
    tmp_double = 13.53515625*obliquity;
    result_by_lmp.set(c_Key3(10, 1, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = -13.53515625*obliquity;
    result_by_lmp.set(c_Key3(10, 1, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 1), result_by_p);
    result_by_p.clear();

    // l , m = (10, 2).
    // p = 4
    tmp_double = 27.0703125;
    result_by_lmp.set(c_Key3(10, 2, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 2), result_by_p);
    result_by_p.clear();

    // l , m = (10, 3).
    // p = 3
    tmp_double = -1407.65625*obliquity;
    result_by_lmp.set(c_Key3(10, 3, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 1407.65625*obliquity;
    result_by_lmp.set(c_Key3(10, 3, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 3), result_by_p);
    result_by_p.clear();

    // l , m = (10, 4).
    // p = 3
    tmp_double = -2815.3125;
    result_by_lmp.set(c_Key3(10, 4, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 4), result_by_p);
    result_by_p.clear();

    // l , m = (10, 5).
    // p = 2
    tmp_double = 126689.0625*obliquity;
    result_by_lmp.set(c_Key3(10, 5, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = -126689.0625*obliquity;
    result_by_lmp.set(c_Key3(10, 5, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 5), result_by_p);
    result_by_p.clear();

    // l , m = (10, 6).
    // p = 2
    tmp_double = 253378.125;
    result_by_lmp.set(c_Key3(10, 6, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 6), result_by_p);
    result_by_p.clear();

    // l , m = (10, 7).
    // p = 1
    tmp_double = -8614856.25*obliquity;
    result_by_lmp.set(c_Key3(10, 7, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 8614856.25*obliquity;
    result_by_lmp.set(c_Key3(10, 7, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 7), result_by_p);
    result_by_p.clear();

    // l , m = (10, 8).
    // p = 1
    tmp_double = -17229712.5;
    result_by_lmp.set(c_Key3(10, 8, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 8), result_by_p);
    result_by_p.clear();

    // l , m = (10, 9).
    // p = 0
    tmp_double = 327364537.5*obliquity;
    result_by_lmp.set(c_Key3(10, 9, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = -327364537.5*obliquity;
    result_by_lmp.set(c_Key3(10, 9, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 9), result_by_p);
    result_by_p.clear();

    // l , m = (10, 10).
    // p = 0
    tmp_double = 654729075.0;
    result_by_lmp.set(c_Key3(10, 10, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 10), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l10_4(double obliquity)
{
    // Inclination Functions Calculated for l = 10.
    // Functions are taylor expanded and truncated at obliquity^4.

    //  Functions Calculated for l = 10.

    c_IntMap<c_Key3, double> result_by_lmp(36);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(11);
    // Optimizations
    double obliquity_2 = obliquity * obliquity;
    double obliquity_3 = obliquity_2 * obliquity;

    c_IntMap<c_Key1, double> result_by_p(10);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (10, 0).
    // p = 4
    tmp_double = -3.3837890625*obliquity_2;
    result_by_lmp.set(c_Key3(10, 0, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 6.767578125*obliquity_2 - 0.24609375;
    result_by_lmp.set(c_Key3(10, 0, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = -3.3837890625*obliquity_2;
    result_by_lmp.set(c_Key3(10, 0, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 0), result_by_p);
    result_by_p.clear();

    // l , m = (10, 1).
    // p = 3
    tmp_double = 58.65234375*obliquity_3;
    result_by_lmp.set(c_Key3(10, 1, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -181.5966796875*obliquity_3 + 13.53515625*obliquity;
    result_by_lmp.set(c_Key3(10, 1, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 184.98046875*obliquity_3 - 13.53515625*obliquity;
    result_by_lmp.set(c_Key3(10, 1, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // p = 6
    tmp_double = -62.0361328125*obliquity_3;
    result_by_lmp.set(c_Key3(10, 1, 6), tmp_double);
    result_by_p.set(c_Key1(6), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 1), result_by_p);
    result_by_p.clear();

    // l , m = (10, 2).
    // p = 3
    tmp_double = 351.9140625*obliquity_2;
    result_by_lmp.set(c_Key3(10, 2, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = 27.0703125 - 717.36328125*obliquity_2;
    result_by_lmp.set(c_Key3(10, 2, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 365.44921875*obliquity_2;
    result_by_lmp.set(c_Key3(10, 2, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 2), result_by_p);
    result_by_p.clear();

    // l , m = (10, 3).
    // p = 2
    tmp_double = -5278.7109375*obliquity_3;
    result_by_lmp.set(c_Key3(10, 3, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 17126.484375*obliquity_3 - 1407.65625*obliquity;
    result_by_lmp.set(c_Key3(10, 3, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -18182.2265625*obliquity_3 + 1407.65625*obliquity;
    result_by_lmp.set(c_Key3(10, 3, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // p = 5
    tmp_double = 6334.453125*obliquity_3;
    result_by_lmp.set(c_Key3(10, 3, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 3), result_by_p);
    result_by_p.clear();

    // l , m = (10, 4).
    // p = 2
    tmp_double = -31672.265625*obliquity_2;
    result_by_lmp.set(c_Key3(10, 4, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 66159.84375*obliquity_2 - 2815.3125;
    result_by_lmp.set(c_Key3(10, 4, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -34487.578125*obliquity_2;
    result_by_lmp.set(c_Key3(10, 4, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 4), result_by_p);
    result_by_p.clear();

    // l , m = (10, 5).
    // p = 1
    tmp_double = 358952.34375*obliquity_3;
    result_by_lmp.set(c_Key3(10, 5, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -1256333.203125*obliquity_3 + 126689.0625*obliquity;
    result_by_lmp.set(c_Key3(10, 5, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 1414694.53125*obliquity_3 - 126689.0625*obliquity;
    result_by_lmp.set(c_Key3(10, 5, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // p = 4
    tmp_double = -517313.671875*obliquity_3;
    result_by_lmp.set(c_Key3(10, 5, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 5), result_by_p);
    result_by_p.clear();

    // l , m = (10, 6).
    // p = 1
    tmp_double = 2153714.0625*obliquity_2;
    result_by_lmp.set(c_Key3(10, 6, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = 253378.125 - 4687495.3125*obliquity_2;
    result_by_lmp.set(c_Key3(10, 6, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 2533781.25*obliquity_2;
    result_by_lmp.set(c_Key3(10, 6, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 6), result_by_p);
    result_by_p.clear();

    // l , m = (10, 7).
    // p = 0
    tmp_double = -13640189.0625*obliquity_3;
    result_by_lmp.set(c_Key3(10, 7, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 57432375.0*obliquity_3 - 8614856.25*obliquity;
    result_by_lmp.set(c_Key3(10, 7, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -72508373.4375*obliquity_3 + 8614856.25*obliquity;
    result_by_lmp.set(c_Key3(10, 7, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // p = 3
    tmp_double = 28716187.5*obliquity_3;
    result_by_lmp.set(c_Key3(10, 7, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 7), result_by_p);
    result_by_p.clear();

    // l , m = (10, 8).
    // p = 0
    tmp_double = -81841134.375*obliquity_2;
    result_by_lmp.set(c_Key3(10, 8, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 198141693.75*obliquity_2 - 17229712.5;
    result_by_lmp.set(c_Key3(10, 8, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -116300559.375*obliquity_2;
    result_by_lmp.set(c_Key3(10, 8, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 8), result_by_p);
    result_by_p.clear();

    // l , m = (10, 9).
    // p = 0
    tmp_double = -791130965.625*obliquity_3 + 327364537.5*obliquity;
    result_by_lmp.set(c_Key3(10, 9, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 1527701175.0*obliquity_3 - 327364537.5*obliquity;
    result_by_lmp.set(c_Key3(10, 9, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // p = 2
    tmp_double = -736570209.375*obliquity_3;
    result_by_lmp.set(c_Key3(10, 9, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 9), result_by_p);
    result_by_p.clear();

    // l , m = (10, 10).
    // p = 0
    tmp_double = 654729075.0 - 1636822687.5*obliquity_2;
    result_by_lmp.set(c_Key3(10, 10, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // p = 1
    tmp_double = 1636822687.5*obliquity_2;
    result_by_lmp.set(c_Key3(10, 10, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 10), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}


ObliquityFuncOutput c_obliquity_function_l10_off(double obliquity)
{
    // Inclination Functions Calculated for l = 10.
    // Functions assume `obliquity=0.0`. However some modes are still activated (non-zero `F(l,m,p)`).

    //  Functions Calculated for l = 10.

    c_IntMap<c_Key3, double> result_by_lmp(6);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lm(6);
    // Optimizations

    c_IntMap<c_Key1, double> result_by_p(10);  // We don't know what size the inner loop will be but it should not be larger than l_
    double tmp_double;

    // Obliquity function by mode:
    // l , m = (10, 0).
    // p = 5
    tmp_double = -0.24609375;
    result_by_lmp.set(c_Key3(10, 0, 5), tmp_double);
    result_by_p.set(c_Key1(5), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 0), result_by_p);
    result_by_p.clear();

    // l , m = (10, 2).
    // p = 4
    tmp_double = 27.0703125;
    result_by_lmp.set(c_Key3(10, 2, 4), tmp_double);
    result_by_p.set(c_Key1(4), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 2), result_by_p);
    result_by_p.clear();

    // l , m = (10, 4).
    // p = 3
    tmp_double = -2815.3125;
    result_by_lmp.set(c_Key3(10, 4, 3), tmp_double);
    result_by_p.set(c_Key1(3), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 4), result_by_p);
    result_by_p.clear();

    // l , m = (10, 6).
    // p = 2
    tmp_double = 253378.125;
    result_by_lmp.set(c_Key3(10, 6, 2), tmp_double);
    result_by_p.set(c_Key1(2), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 6), result_by_p);
    result_by_p.clear();

    // l , m = (10, 8).
    // p = 1
    tmp_double = -17229712.5;
    result_by_lmp.set(c_Key3(10, 8, 1), tmp_double);
    result_by_p.set(c_Key1(1), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 8), result_by_p);
    result_by_p.clear();

    // l , m = (10, 10).
    // p = 0
    tmp_double = 654729075.0;
    result_by_lmp.set(c_Key3(10, 10, 0), tmp_double);
    result_by_p.set(c_Key1(0), tmp_double);
    // Store the p table into the results_lm then reset it.
    
    result_by_lm.set(c_Key2(10, 10), result_by_p);
    result_by_p.clear();

    return ObliquityFuncOutput(result_by_lmp, result_by_lm);
}
