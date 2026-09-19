#pragma once

#include "eccentricity_common_.hpp"

EccentricityFuncOutput c_eccentricity_function_l8_e1(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 8.
    // Functions are taylor expanded and truncated at eccentricity^1.

    // Functions Calculated for l = 8.

    c_IntMap<c_Key3, double> result_by_lpq(27);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(9);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = -3.5*eccentricity;
    double common_term_1 = 12.5*eccentricity;
    double common_term_2 = -1.5*eccentricity;
    double common_term_3 = 10.5*eccentricity;
    double common_term_4 = 0.5*eccentricity;
    double common_term_5 = 8.5*eccentricity;
    double common_term_6 = 2.5*eccentricity;
    double common_term_7 = 6.5*eccentricity;
    double common_term_8 = 4.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(3);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (8, 0).
    // q = -1
    result_by_lpq.set(c_Key3(8, 0, -1), common_term_0);
    result_by_q.set(c_Key1(-1), common_term_0);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(8, 0, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(8, 0, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 0), result_by_q);
    result_by_q.clear();

    // l , p = (8, 1).
    // q = -1
    result_by_lpq.set(c_Key3(8, 1, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(8, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(8, 1, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 1), result_by_q);
    result_by_q.clear();

    // l , p = (8, 2).
    // q = -1
    result_by_lpq.set(c_Key3(8, 2, -1), common_term_4);
    result_by_q.set(c_Key1(-1), common_term_4);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(8, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(8, 2, 1), common_term_5);
    result_by_q.set(c_Key1(1), common_term_5);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 2), result_by_q);
    result_by_q.clear();

    // l , p = (8, 3).
    // q = -1
    result_by_lpq.set(c_Key3(8, 3, -1), common_term_6);
    result_by_q.set(c_Key1(-1), common_term_6);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(8, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(8, 3, 1), common_term_7);
    result_by_q.set(c_Key1(1), common_term_7);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 3), result_by_q);
    result_by_q.clear();

    // l , p = (8, 4).
    // q = -1
    result_by_lpq.set(c_Key3(8, 4, -1), common_term_8);
    result_by_q.set(c_Key1(-1), common_term_8);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -7.5);
    result_by_lpq.set(c_Key3(8, 4, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(8, 4, 1), common_term_8);
    result_by_q.set(c_Key1(1), common_term_8);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 4), result_by_q);
    result_by_q.clear();

    // l , p = (8, 5).
    // q = -1
    result_by_lpq.set(c_Key3(8, 5, -1), common_term_7);
    result_by_q.set(c_Key1(-1), common_term_7);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(8, 5, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(8, 5, 1), common_term_6);
    result_by_q.set(c_Key1(1), common_term_6);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 5), result_by_q);
    result_by_q.clear();

    // l , p = (8, 6).
    // q = -1
    result_by_lpq.set(c_Key3(8, 6, -1), common_term_5);
    result_by_q.set(c_Key1(-1), common_term_5);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(8, 6, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(8, 6, 1), common_term_4);
    result_by_q.set(c_Key1(1), common_term_4);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 6), result_by_q);
    result_by_q.clear();

    // l , p = (8, 7).
    // q = -1
    result_by_lpq.set(c_Key3(8, 7, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(8, 7, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(8, 7, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 7), result_by_q);
    result_by_q.clear();

    // l , p = (8, 8).
    // q = -1
    result_by_lpq.set(c_Key3(8, 8, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(8, 8, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(8, 8, 1), common_term_0);
    result_by_q.set(c_Key1(1), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 8), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l8_e2(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 8.
    // Functions are taylor expanded and truncated at eccentricity^2.

    // Functions Calculated for l = 8.

    c_IntMap<c_Key3, double> result_by_lpq(45);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(9);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = 4.5*eccentricity_2;
    double common_term_1 = -3.5*eccentricity;
    double common_term_2 = 1.0 - 46.0*eccentricity_2;
    double common_term_3 = 12.5*eccentricity;
    double common_term_4 = 86.5*eccentricity_2;
    double common_term_5 = 0.75*eccentricity_2;
    double common_term_6 = -1.5*eccentricity;
    double common_term_7 = 1.0 - 18.0*eccentricity_2;
    double common_term_8 = 10.5*eccentricity;
    double common_term_9 = 62.25*eccentricity_2;
    double common_term_10 = 0.5*eccentricity;
    double common_term_11 = 2.0*eccentricity_2 + 1.0;
    double common_term_12 = 8.5*eccentricity;
    double common_term_13 = 42.0*eccentricity_2;
    double common_term_14 = 5.25*eccentricity_2*std::pow(1.0 - eccentricity_2, -7.5);
    double common_term_15 = 2.5*eccentricity;
    double common_term_16 = 14.0*eccentricity_2 + 1.0;
    double common_term_17 = 6.5*eccentricity;
    double common_term_18 = 25.75*eccentricity_2;
    double common_term_19 = 13.5*eccentricity_2;
    double common_term_20 = 4.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(5);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (8, 0).
    // q = -2
    result_by_lpq.set(c_Key3(8, 0, -2), common_term_0);
    result_by_q.set(c_Key1(-2), common_term_0);
    // q = -1
    result_by_lpq.set(c_Key3(8, 0, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    result_by_lpq.set(c_Key3(8, 0, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(8, 0, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // q = 2
    result_by_lpq.set(c_Key3(8, 0, 2), common_term_4);
    result_by_q.set(c_Key1(2), common_term_4);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 0), result_by_q);
    result_by_q.clear();

    // l , p = (8, 1).
    // q = -2
    result_by_lpq.set(c_Key3(8, 1, -2), common_term_5);
    result_by_q.set(c_Key1(-2), common_term_5);
    // q = -1
    result_by_lpq.set(c_Key3(8, 1, -1), common_term_6);
    result_by_q.set(c_Key1(-1), common_term_6);
    // q = 0
    result_by_lpq.set(c_Key3(8, 1, 0), common_term_7);
    result_by_q.set(c_Key1(0), common_term_7);
    // q = 1
    result_by_lpq.set(c_Key3(8, 1, 1), common_term_8);
    result_by_q.set(c_Key1(1), common_term_8);
    // q = 2
    result_by_lpq.set(c_Key3(8, 1, 2), common_term_9);
    result_by_q.set(c_Key1(2), common_term_9);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 1), result_by_q);
    result_by_q.clear();

    // l , p = (8, 2).
    // q = -2
    result_by_lpq.set(c_Key3(8, 2, -2), eccentricity_2);
    result_by_q.set(c_Key1(-2), eccentricity_2);
    // q = -1
    result_by_lpq.set(c_Key3(8, 2, -1), common_term_10);
    result_by_q.set(c_Key1(-1), common_term_10);
    // q = 0
    result_by_lpq.set(c_Key3(8, 2, 0), common_term_11);
    result_by_q.set(c_Key1(0), common_term_11);
    // q = 1
    result_by_lpq.set(c_Key3(8, 2, 1), common_term_12);
    result_by_q.set(c_Key1(1), common_term_12);
    // q = 2
    result_by_lpq.set(c_Key3(8, 2, 2), common_term_13);
    result_by_q.set(c_Key1(2), common_term_13);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 2), result_by_q);
    result_by_q.clear();

    // l , p = (8, 3).
    // q = -2
    result_by_lpq.set(c_Key3(8, 3, -2), common_term_14);
    result_by_q.set(c_Key1(-2), common_term_14);
    // q = -1
    result_by_lpq.set(c_Key3(8, 3, -1), common_term_15);
    result_by_q.set(c_Key1(-1), common_term_15);
    // q = 0
    result_by_lpq.set(c_Key3(8, 3, 0), common_term_16);
    result_by_q.set(c_Key1(0), common_term_16);
    // q = 1
    result_by_lpq.set(c_Key3(8, 3, 1), common_term_17);
    result_by_q.set(c_Key1(1), common_term_17);
    // q = 2
    result_by_lpq.set(c_Key3(8, 3, 2), common_term_18);
    result_by_q.set(c_Key1(2), common_term_18);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 3), result_by_q);
    result_by_q.clear();

    // l , p = (8, 4).
    // q = -2
    result_by_lpq.set(c_Key3(8, 4, -2), common_term_19);
    result_by_q.set(c_Key1(-2), common_term_19);
    // q = -1
    result_by_lpq.set(c_Key3(8, 4, -1), common_term_20);
    result_by_q.set(c_Key1(-1), common_term_20);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -7.5)*(10.5*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(8, 4, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(8, 4, 1), common_term_20);
    result_by_q.set(c_Key1(1), common_term_20);
    // q = 2
    result_by_lpq.set(c_Key3(8, 4, 2), common_term_19);
    result_by_q.set(c_Key1(2), common_term_19);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 4), result_by_q);
    result_by_q.clear();

    // l , p = (8, 5).
    // q = -2
    result_by_lpq.set(c_Key3(8, 5, -2), common_term_18);
    result_by_q.set(c_Key1(-2), common_term_18);
    // q = -1
    result_by_lpq.set(c_Key3(8, 5, -1), common_term_17);
    result_by_q.set(c_Key1(-1), common_term_17);
    // q = 0
    result_by_lpq.set(c_Key3(8, 5, 0), common_term_16);
    result_by_q.set(c_Key1(0), common_term_16);
    // q = 1
    result_by_lpq.set(c_Key3(8, 5, 1), common_term_15);
    result_by_q.set(c_Key1(1), common_term_15);
    // q = 2
    result_by_lpq.set(c_Key3(8, 5, 2), common_term_14);
    result_by_q.set(c_Key1(2), common_term_14);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 5), result_by_q);
    result_by_q.clear();

    // l , p = (8, 6).
    // q = -2
    result_by_lpq.set(c_Key3(8, 6, -2), common_term_13);
    result_by_q.set(c_Key1(-2), common_term_13);
    // q = -1
    result_by_lpq.set(c_Key3(8, 6, -1), common_term_12);
    result_by_q.set(c_Key1(-1), common_term_12);
    // q = 0
    result_by_lpq.set(c_Key3(8, 6, 0), common_term_11);
    result_by_q.set(c_Key1(0), common_term_11);
    // q = 1
    result_by_lpq.set(c_Key3(8, 6, 1), common_term_10);
    result_by_q.set(c_Key1(1), common_term_10);
    // q = 2
    result_by_lpq.set(c_Key3(8, 6, 2), eccentricity_2);
    result_by_q.set(c_Key1(2), eccentricity_2);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 6), result_by_q);
    result_by_q.clear();

    // l , p = (8, 7).
    // q = -2
    result_by_lpq.set(c_Key3(8, 7, -2), common_term_9);
    result_by_q.set(c_Key1(-2), common_term_9);
    // q = -1
    result_by_lpq.set(c_Key3(8, 7, -1), common_term_8);
    result_by_q.set(c_Key1(-1), common_term_8);
    // q = 0
    result_by_lpq.set(c_Key3(8, 7, 0), common_term_7);
    result_by_q.set(c_Key1(0), common_term_7);
    // q = 1
    result_by_lpq.set(c_Key3(8, 7, 1), common_term_6);
    result_by_q.set(c_Key1(1), common_term_6);
    // q = 2
    result_by_lpq.set(c_Key3(8, 7, 2), common_term_5);
    result_by_q.set(c_Key1(2), common_term_5);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 7), result_by_q);
    result_by_q.clear();

    // l , p = (8, 8).
    // q = -2
    result_by_lpq.set(c_Key3(8, 8, -2), common_term_4);
    result_by_q.set(c_Key1(-2), common_term_4);
    // q = -1
    result_by_lpq.set(c_Key3(8, 8, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    result_by_lpq.set(c_Key3(8, 8, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(8, 8, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // q = 2
    result_by_lpq.set(c_Key3(8, 8, 2), common_term_0);
    result_by_q.set(c_Key1(2), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 8), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l8_e3(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 8.
    // Functions are taylor expanded and truncated at eccentricity^3.

    // Functions Calculated for l = 8.

    c_IntMap<c_Key3, double> result_by_lpq(63);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(9);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double common_term_0 = -2.6041666666666667*eccentricity_3;
    double common_term_1 = 4.5*eccentricity_2;
    double common_term_2 = 63.4375*eccentricity_3 - 3.5*eccentricity;
    double common_term_3 = 1.0 - 46.0*eccentricity_2;
    double common_term_4 = -333.5625*eccentricity_3 + 12.5*eccentricity;
    double common_term_5 = 86.5*eccentricity_2;
    double common_term_6 = 437.72916666666667*eccentricity_3;
    double common_term_7 = -0.0625*eccentricity_3;
    double common_term_8 = 0.75*eccentricity_2;
    double common_term_9 = 10.3125*eccentricity_3 - 1.5*eccentricity;
    double common_term_10 = 1.0 - 18.0*eccentricity_2;
    double common_term_11 = -119.4375*eccentricity_3 + 10.5*eccentricity;
    double common_term_12 = 62.25*eccentricity_2;
    double common_term_13 = 274.1875*eccentricity_3;
    double common_term_14 = 1.4791666666666667*eccentricity_3;
    double common_term_15 = 6.1875*eccentricity_3 + 0.5*eccentricity;
    double common_term_16 = 2.0*eccentricity_2 + 1.0;
    double common_term_17 = -0.3125*eccentricity_3 + 8.5*eccentricity;
    double common_term_18 = 42.0*eccentricity_2;
    double common_term_19 = 157.64583333333333*eccentricity_3;
    double common_term_20 = 10.020833333333333*eccentricity_3;
    double common_term_21 = 5.25*eccentricity_2*std::pow(1.0 - eccentricity_2, -7.5);
    double common_term_22 = 27.0625*eccentricity_3 + 2.5*eccentricity;
    double common_term_23 = 14.0*eccentricity_2 + 1.0;
    double common_term_24 = 47.8125*eccentricity_3 + 6.5*eccentricity;
    double common_term_25 = 25.75*eccentricity_2;
    double common_term_26 = 80.104166666666667*eccentricity_3;
    double common_term_27 = 33.5625*eccentricity_3;
    double common_term_28 = 13.5*eccentricity_2;
    double common_term_29 = 48.9375*eccentricity_3 + 4.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(7);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (8, 0).
    // q = -3
    result_by_lpq.set(c_Key3(8, 0, -3), common_term_0);
    result_by_q.set(c_Key1(-3), common_term_0);
    // q = -2
    result_by_lpq.set(c_Key3(8, 0, -2), common_term_1);
    result_by_q.set(c_Key1(-2), common_term_1);
    // q = -1
    result_by_lpq.set(c_Key3(8, 0, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    result_by_lpq.set(c_Key3(8, 0, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(8, 0, 1), common_term_4);
    result_by_q.set(c_Key1(1), common_term_4);
    // q = 2
    result_by_lpq.set(c_Key3(8, 0, 2), common_term_5);
    result_by_q.set(c_Key1(2), common_term_5);
    // q = 3
    result_by_lpq.set(c_Key3(8, 0, 3), common_term_6);
    result_by_q.set(c_Key1(3), common_term_6);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 0), result_by_q);
    result_by_q.clear();

    // l , p = (8, 1).
    // q = -3
    result_by_lpq.set(c_Key3(8, 1, -3), common_term_7);
    result_by_q.set(c_Key1(-3), common_term_7);
    // q = -2
    result_by_lpq.set(c_Key3(8, 1, -2), common_term_8);
    result_by_q.set(c_Key1(-2), common_term_8);
    // q = -1
    result_by_lpq.set(c_Key3(8, 1, -1), common_term_9);
    result_by_q.set(c_Key1(-1), common_term_9);
    // q = 0
    result_by_lpq.set(c_Key3(8, 1, 0), common_term_10);
    result_by_q.set(c_Key1(0), common_term_10);
    // q = 1
    result_by_lpq.set(c_Key3(8, 1, 1), common_term_11);
    result_by_q.set(c_Key1(1), common_term_11);
    // q = 2
    result_by_lpq.set(c_Key3(8, 1, 2), common_term_12);
    result_by_q.set(c_Key1(2), common_term_12);
    // q = 3
    result_by_lpq.set(c_Key3(8, 1, 3), common_term_13);
    result_by_q.set(c_Key1(3), common_term_13);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 1), result_by_q);
    result_by_q.clear();

    // l , p = (8, 2).
    // q = -3
    result_by_lpq.set(c_Key3(8, 2, -3), common_term_14);
    result_by_q.set(c_Key1(-3), common_term_14);
    // q = -2
    result_by_lpq.set(c_Key3(8, 2, -2), eccentricity_2);
    result_by_q.set(c_Key1(-2), eccentricity_2);
    // q = -1
    result_by_lpq.set(c_Key3(8, 2, -1), common_term_15);
    result_by_q.set(c_Key1(-1), common_term_15);
    // q = 0
    result_by_lpq.set(c_Key3(8, 2, 0), common_term_16);
    result_by_q.set(c_Key1(0), common_term_16);
    // q = 1
    result_by_lpq.set(c_Key3(8, 2, 1), common_term_17);
    result_by_q.set(c_Key1(1), common_term_17);
    // q = 2
    result_by_lpq.set(c_Key3(8, 2, 2), common_term_18);
    result_by_q.set(c_Key1(2), common_term_18);
    // q = 3
    result_by_lpq.set(c_Key3(8, 2, 3), common_term_19);
    result_by_q.set(c_Key1(3), common_term_19);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 2), result_by_q);
    result_by_q.clear();

    // l , p = (8, 3).
    // q = -3
    result_by_lpq.set(c_Key3(8, 3, -3), common_term_20);
    result_by_q.set(c_Key1(-3), common_term_20);
    // q = -2
    result_by_lpq.set(c_Key3(8, 3, -2), common_term_21);
    result_by_q.set(c_Key1(-2), common_term_21);
    // q = -1
    result_by_lpq.set(c_Key3(8, 3, -1), common_term_22);
    result_by_q.set(c_Key1(-1), common_term_22);
    // q = 0
    result_by_lpq.set(c_Key3(8, 3, 0), common_term_23);
    result_by_q.set(c_Key1(0), common_term_23);
    // q = 1
    result_by_lpq.set(c_Key3(8, 3, 1), common_term_24);
    result_by_q.set(c_Key1(1), common_term_24);
    // q = 2
    result_by_lpq.set(c_Key3(8, 3, 2), common_term_25);
    result_by_q.set(c_Key1(2), common_term_25);
    // q = 3
    result_by_lpq.set(c_Key3(8, 3, 3), common_term_26);
    result_by_q.set(c_Key1(3), common_term_26);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 3), result_by_q);
    result_by_q.clear();

    // l , p = (8, 4).
    // q = -3
    result_by_lpq.set(c_Key3(8, 4, -3), common_term_27);
    result_by_q.set(c_Key1(-3), common_term_27);
    // q = -2
    result_by_lpq.set(c_Key3(8, 4, -2), common_term_28);
    result_by_q.set(c_Key1(-2), common_term_28);
    // q = -1
    result_by_lpq.set(c_Key3(8, 4, -1), common_term_29);
    result_by_q.set(c_Key1(-1), common_term_29);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -7.5)*(10.5*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(8, 4, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(8, 4, 1), common_term_29);
    result_by_q.set(c_Key1(1), common_term_29);
    // q = 2
    result_by_lpq.set(c_Key3(8, 4, 2), common_term_28);
    result_by_q.set(c_Key1(2), common_term_28);
    // q = 3
    result_by_lpq.set(c_Key3(8, 4, 3), common_term_27);
    result_by_q.set(c_Key1(3), common_term_27);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 4), result_by_q);
    result_by_q.clear();

    // l , p = (8, 5).
    // q = -3
    result_by_lpq.set(c_Key3(8, 5, -3), common_term_26);
    result_by_q.set(c_Key1(-3), common_term_26);
    // q = -2
    result_by_lpq.set(c_Key3(8, 5, -2), common_term_25);
    result_by_q.set(c_Key1(-2), common_term_25);
    // q = -1
    result_by_lpq.set(c_Key3(8, 5, -1), common_term_24);
    result_by_q.set(c_Key1(-1), common_term_24);
    // q = 0
    result_by_lpq.set(c_Key3(8, 5, 0), common_term_23);
    result_by_q.set(c_Key1(0), common_term_23);
    // q = 1
    result_by_lpq.set(c_Key3(8, 5, 1), common_term_22);
    result_by_q.set(c_Key1(1), common_term_22);
    // q = 2
    result_by_lpq.set(c_Key3(8, 5, 2), common_term_21);
    result_by_q.set(c_Key1(2), common_term_21);
    // q = 3
    result_by_lpq.set(c_Key3(8, 5, 3), common_term_20);
    result_by_q.set(c_Key1(3), common_term_20);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 5), result_by_q);
    result_by_q.clear();

    // l , p = (8, 6).
    // q = -3
    result_by_lpq.set(c_Key3(8, 6, -3), common_term_19);
    result_by_q.set(c_Key1(-3), common_term_19);
    // q = -2
    result_by_lpq.set(c_Key3(8, 6, -2), common_term_18);
    result_by_q.set(c_Key1(-2), common_term_18);
    // q = -1
    result_by_lpq.set(c_Key3(8, 6, -1), common_term_17);
    result_by_q.set(c_Key1(-1), common_term_17);
    // q = 0
    result_by_lpq.set(c_Key3(8, 6, 0), common_term_16);
    result_by_q.set(c_Key1(0), common_term_16);
    // q = 1
    result_by_lpq.set(c_Key3(8, 6, 1), common_term_15);
    result_by_q.set(c_Key1(1), common_term_15);
    // q = 2
    result_by_lpq.set(c_Key3(8, 6, 2), eccentricity_2);
    result_by_q.set(c_Key1(2), eccentricity_2);
    // q = 3
    result_by_lpq.set(c_Key3(8, 6, 3), common_term_14);
    result_by_q.set(c_Key1(3), common_term_14);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 6), result_by_q);
    result_by_q.clear();

    // l , p = (8, 7).
    // q = -3
    result_by_lpq.set(c_Key3(8, 7, -3), common_term_13);
    result_by_q.set(c_Key1(-3), common_term_13);
    // q = -2
    result_by_lpq.set(c_Key3(8, 7, -2), common_term_12);
    result_by_q.set(c_Key1(-2), common_term_12);
    // q = -1
    result_by_lpq.set(c_Key3(8, 7, -1), common_term_11);
    result_by_q.set(c_Key1(-1), common_term_11);
    // q = 0
    result_by_lpq.set(c_Key3(8, 7, 0), common_term_10);
    result_by_q.set(c_Key1(0), common_term_10);
    // q = 1
    result_by_lpq.set(c_Key3(8, 7, 1), common_term_9);
    result_by_q.set(c_Key1(1), common_term_9);
    // q = 2
    result_by_lpq.set(c_Key3(8, 7, 2), common_term_8);
    result_by_q.set(c_Key1(2), common_term_8);
    // q = 3
    result_by_lpq.set(c_Key3(8, 7, 3), common_term_7);
    result_by_q.set(c_Key1(3), common_term_7);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 7), result_by_q);
    result_by_q.clear();

    // l , p = (8, 8).
    // q = -3
    result_by_lpq.set(c_Key3(8, 8, -3), common_term_6);
    result_by_q.set(c_Key1(-3), common_term_6);
    // q = -2
    result_by_lpq.set(c_Key3(8, 8, -2), common_term_5);
    result_by_q.set(c_Key1(-2), common_term_5);
    // q = -1
    result_by_lpq.set(c_Key3(8, 8, -1), common_term_4);
    result_by_q.set(c_Key1(-1), common_term_4);
    // q = 0
    result_by_lpq.set(c_Key3(8, 8, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(8, 8, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // q = 2
    result_by_lpq.set(c_Key3(8, 8, 2), common_term_1);
    result_by_q.set(c_Key1(2), common_term_1);
    // q = 3
    result_by_lpq.set(c_Key3(8, 8, 3), common_term_0);
    result_by_q.set(c_Key1(3), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 8), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l8_e4(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 8.
    // Functions are taylor expanded and truncated at eccentricity^4.

    // Functions Calculated for l = 8.

    c_IntMap<c_Key3, double> result_by_lpq(81);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(9);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double common_term_0 = 0.66666666666666667*eccentricity_4;
    double common_term_1 = -2.6041666666666667*eccentricity_3;
    double common_term_2 = -40.5*eccentricity_4 + 4.5*eccentricity_2;
    double common_term_3 = 63.4375*eccentricity_3 - 3.5*eccentricity;
    double common_term_4 = 490.75*eccentricity_4 - 46.0*eccentricity_2 + 1.0;
    double common_term_5 = -333.5625*eccentricity_3 + 12.5*eccentricity;
    double common_term_6 = -1764.1666666666667*eccentricity_4 + 86.5*eccentricity_2;
    double common_term_7 = 437.72916666666667*eccentricity_3;
    double common_term_8 = 1808.25*eccentricity_4;
    double common_term_9 = 0.0625*eccentricity_4;
    double common_term_10 = -0.0625*eccentricity_3;
    double common_term_11 = -1.625*eccentricity_4 + 0.75*eccentricity_2;
    double common_term_12 = 10.3125*eccentricity_3 - 1.5*eccentricity;
    double common_term_13 = 78.1875*eccentricity_4 - 18.0*eccentricity_2 + 1.0;
    double common_term_14 = -119.4375*eccentricity_3 + 10.5*eccentricity;
    double common_term_15 = -580.375*eccentricity_4 + 62.25*eccentricity_2;
    double common_term_16 = 274.1875*eccentricity_3;
    double common_term_17 = 998.75*eccentricity_4;
    double common_term_18 = 2.1875*eccentricity_4*std::pow(1.0 - eccentricity_2, -7.5);
    double common_term_19 = 1.4791666666666667*eccentricity_3;
    double common_term_20 = 8.5833333333333333*eccentricity_4 + eccentricity_2;
    double common_term_21 = 6.1875*eccentricity_3 + 0.5*eccentricity;
    double common_term_22 = 23.5*eccentricity_4 + 2.0*eccentricity_2 + 1.0;
    double common_term_23 = -0.3125*eccentricity_3 + 8.5*eccentricity;
    double common_term_24 = -37.25*eccentricity_4 + 42.0*eccentricity_2;
    double common_term_25 = 157.64583333333333*eccentricity_3;
    double common_term_26 = 497.97916666666667*eccentricity_4;
    double common_term_27 = 18.041666666666667*eccentricity_4;
    double common_term_28 = 10.020833333333333*eccentricity_3;
    double common_term_29 = std::pow(1.0 - eccentricity_2, -7.5)*(8.75*eccentricity_4 + 5.25*eccentricity_2);
    double common_term_30 = 27.0625*eccentricity_3 + 2.5*eccentricity;
    double common_term_31 = 86.6875*eccentricity_4 + 14.0*eccentricity_2 + 1.0;
    double common_term_32 = 47.8125*eccentricity_3 + 6.5*eccentricity;
    double common_term_33 = 127.20833333333333*eccentricity_4 + 25.75*eccentricity_2;
    double common_term_34 = 80.104166666666667*eccentricity_3;
    double common_term_35 = 214.9375*eccentricity_4;
    double common_term_36 = 74.625*eccentricity_4;
    double common_term_37 = 33.5625*eccentricity_3;
    double common_term_38 = 111.0*eccentricity_4 + 13.5*eccentricity_2;
    double common_term_39 = 48.9375*eccentricity_3 + 4.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(9);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (8, 0).
    // q = -4
    result_by_lpq.set(c_Key3(8, 0, -4), common_term_0);
    result_by_q.set(c_Key1(-4), common_term_0);
    // q = -3
    result_by_lpq.set(c_Key3(8, 0, -3), common_term_1);
    result_by_q.set(c_Key1(-3), common_term_1);
    // q = -2
    result_by_lpq.set(c_Key3(8, 0, -2), common_term_2);
    result_by_q.set(c_Key1(-2), common_term_2);
    // q = -1
    result_by_lpq.set(c_Key3(8, 0, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    result_by_lpq.set(c_Key3(8, 0, 0), common_term_4);
    result_by_q.set(c_Key1(0), common_term_4);
    // q = 1
    result_by_lpq.set(c_Key3(8, 0, 1), common_term_5);
    result_by_q.set(c_Key1(1), common_term_5);
    // q = 2
    result_by_lpq.set(c_Key3(8, 0, 2), common_term_6);
    result_by_q.set(c_Key1(2), common_term_6);
    // q = 3
    result_by_lpq.set(c_Key3(8, 0, 3), common_term_7);
    result_by_q.set(c_Key1(3), common_term_7);
    // q = 4
    result_by_lpq.set(c_Key3(8, 0, 4), common_term_8);
    result_by_q.set(c_Key1(4), common_term_8);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 0), result_by_q);
    result_by_q.clear();

    // l , p = (8, 1).
    // q = -4
    result_by_lpq.set(c_Key3(8, 1, -4), common_term_9);
    result_by_q.set(c_Key1(-4), common_term_9);
    // q = -3
    result_by_lpq.set(c_Key3(8, 1, -3), common_term_10);
    result_by_q.set(c_Key1(-3), common_term_10);
    // q = -2
    result_by_lpq.set(c_Key3(8, 1, -2), common_term_11);
    result_by_q.set(c_Key1(-2), common_term_11);
    // q = -1
    result_by_lpq.set(c_Key3(8, 1, -1), common_term_12);
    result_by_q.set(c_Key1(-1), common_term_12);
    // q = 0
    result_by_lpq.set(c_Key3(8, 1, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(8, 1, 1), common_term_14);
    result_by_q.set(c_Key1(1), common_term_14);
    // q = 2
    result_by_lpq.set(c_Key3(8, 1, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // q = 3
    result_by_lpq.set(c_Key3(8, 1, 3), common_term_16);
    result_by_q.set(c_Key1(3), common_term_16);
    // q = 4
    result_by_lpq.set(c_Key3(8, 1, 4), common_term_17);
    result_by_q.set(c_Key1(4), common_term_17);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 1), result_by_q);
    result_by_q.clear();

    // l , p = (8, 2).
    // q = -4
    result_by_lpq.set(c_Key3(8, 2, -4), common_term_18);
    result_by_q.set(c_Key1(-4), common_term_18);
    // q = -3
    result_by_lpq.set(c_Key3(8, 2, -3), common_term_19);
    result_by_q.set(c_Key1(-3), common_term_19);
    // q = -2
    result_by_lpq.set(c_Key3(8, 2, -2), common_term_20);
    result_by_q.set(c_Key1(-2), common_term_20);
    // q = -1
    result_by_lpq.set(c_Key3(8, 2, -1), common_term_21);
    result_by_q.set(c_Key1(-1), common_term_21);
    // q = 0
    result_by_lpq.set(c_Key3(8, 2, 0), common_term_22);
    result_by_q.set(c_Key1(0), common_term_22);
    // q = 1
    result_by_lpq.set(c_Key3(8, 2, 1), common_term_23);
    result_by_q.set(c_Key1(1), common_term_23);
    // q = 2
    result_by_lpq.set(c_Key3(8, 2, 2), common_term_24);
    result_by_q.set(c_Key1(2), common_term_24);
    // q = 3
    result_by_lpq.set(c_Key3(8, 2, 3), common_term_25);
    result_by_q.set(c_Key1(3), common_term_25);
    // q = 4
    result_by_lpq.set(c_Key3(8, 2, 4), common_term_26);
    result_by_q.set(c_Key1(4), common_term_26);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 2), result_by_q);
    result_by_q.clear();

    // l , p = (8, 3).
    // q = -4
    result_by_lpq.set(c_Key3(8, 3, -4), common_term_27);
    result_by_q.set(c_Key1(-4), common_term_27);
    // q = -3
    result_by_lpq.set(c_Key3(8, 3, -3), common_term_28);
    result_by_q.set(c_Key1(-3), common_term_28);
    // q = -2
    result_by_lpq.set(c_Key3(8, 3, -2), common_term_29);
    result_by_q.set(c_Key1(-2), common_term_29);
    // q = -1
    result_by_lpq.set(c_Key3(8, 3, -1), common_term_30);
    result_by_q.set(c_Key1(-1), common_term_30);
    // q = 0
    result_by_lpq.set(c_Key3(8, 3, 0), common_term_31);
    result_by_q.set(c_Key1(0), common_term_31);
    // q = 1
    result_by_lpq.set(c_Key3(8, 3, 1), common_term_32);
    result_by_q.set(c_Key1(1), common_term_32);
    // q = 2
    result_by_lpq.set(c_Key3(8, 3, 2), common_term_33);
    result_by_q.set(c_Key1(2), common_term_33);
    // q = 3
    result_by_lpq.set(c_Key3(8, 3, 3), common_term_34);
    result_by_q.set(c_Key1(3), common_term_34);
    // q = 4
    result_by_lpq.set(c_Key3(8, 3, 4), common_term_35);
    result_by_q.set(c_Key1(4), common_term_35);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 3), result_by_q);
    result_by_q.clear();

    // l , p = (8, 4).
    // q = -4
    result_by_lpq.set(c_Key3(8, 4, -4), common_term_36);
    result_by_q.set(c_Key1(-4), common_term_36);
    // q = -3
    result_by_lpq.set(c_Key3(8, 4, -3), common_term_37);
    result_by_q.set(c_Key1(-3), common_term_37);
    // q = -2
    result_by_lpq.set(c_Key3(8, 4, -2), common_term_38);
    result_by_q.set(c_Key1(-2), common_term_38);
    // q = -1
    result_by_lpq.set(c_Key3(8, 4, -1), common_term_39);
    result_by_q.set(c_Key1(-1), common_term_39);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -7.5)*(13.125*eccentricity_4 + 10.5*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(8, 4, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(8, 4, 1), common_term_39);
    result_by_q.set(c_Key1(1), common_term_39);
    // q = 2
    result_by_lpq.set(c_Key3(8, 4, 2), common_term_38);
    result_by_q.set(c_Key1(2), common_term_38);
    // q = 3
    result_by_lpq.set(c_Key3(8, 4, 3), common_term_37);
    result_by_q.set(c_Key1(3), common_term_37);
    // q = 4
    result_by_lpq.set(c_Key3(8, 4, 4), common_term_36);
    result_by_q.set(c_Key1(4), common_term_36);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 4), result_by_q);
    result_by_q.clear();

    // l , p = (8, 5).
    // q = -4
    result_by_lpq.set(c_Key3(8, 5, -4), common_term_35);
    result_by_q.set(c_Key1(-4), common_term_35);
    // q = -3
    result_by_lpq.set(c_Key3(8, 5, -3), common_term_34);
    result_by_q.set(c_Key1(-3), common_term_34);
    // q = -2
    result_by_lpq.set(c_Key3(8, 5, -2), common_term_33);
    result_by_q.set(c_Key1(-2), common_term_33);
    // q = -1
    result_by_lpq.set(c_Key3(8, 5, -1), common_term_32);
    result_by_q.set(c_Key1(-1), common_term_32);
    // q = 0
    result_by_lpq.set(c_Key3(8, 5, 0), common_term_31);
    result_by_q.set(c_Key1(0), common_term_31);
    // q = 1
    result_by_lpq.set(c_Key3(8, 5, 1), common_term_30);
    result_by_q.set(c_Key1(1), common_term_30);
    // q = 2
    result_by_lpq.set(c_Key3(8, 5, 2), common_term_29);
    result_by_q.set(c_Key1(2), common_term_29);
    // q = 3
    result_by_lpq.set(c_Key3(8, 5, 3), common_term_28);
    result_by_q.set(c_Key1(3), common_term_28);
    // q = 4
    result_by_lpq.set(c_Key3(8, 5, 4), common_term_27);
    result_by_q.set(c_Key1(4), common_term_27);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 5), result_by_q);
    result_by_q.clear();

    // l , p = (8, 6).
    // q = -4
    result_by_lpq.set(c_Key3(8, 6, -4), common_term_26);
    result_by_q.set(c_Key1(-4), common_term_26);
    // q = -3
    result_by_lpq.set(c_Key3(8, 6, -3), common_term_25);
    result_by_q.set(c_Key1(-3), common_term_25);
    // q = -2
    result_by_lpq.set(c_Key3(8, 6, -2), common_term_24);
    result_by_q.set(c_Key1(-2), common_term_24);
    // q = -1
    result_by_lpq.set(c_Key3(8, 6, -1), common_term_23);
    result_by_q.set(c_Key1(-1), common_term_23);
    // q = 0
    result_by_lpq.set(c_Key3(8, 6, 0), common_term_22);
    result_by_q.set(c_Key1(0), common_term_22);
    // q = 1
    result_by_lpq.set(c_Key3(8, 6, 1), common_term_21);
    result_by_q.set(c_Key1(1), common_term_21);
    // q = 2
    result_by_lpq.set(c_Key3(8, 6, 2), common_term_20);
    result_by_q.set(c_Key1(2), common_term_20);
    // q = 3
    result_by_lpq.set(c_Key3(8, 6, 3), common_term_19);
    result_by_q.set(c_Key1(3), common_term_19);
    // q = 4
    result_by_lpq.set(c_Key3(8, 6, 4), common_term_18);
    result_by_q.set(c_Key1(4), common_term_18);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 6), result_by_q);
    result_by_q.clear();

    // l , p = (8, 7).
    // q = -4
    result_by_lpq.set(c_Key3(8, 7, -4), common_term_17);
    result_by_q.set(c_Key1(-4), common_term_17);
    // q = -3
    result_by_lpq.set(c_Key3(8, 7, -3), common_term_16);
    result_by_q.set(c_Key1(-3), common_term_16);
    // q = -2
    result_by_lpq.set(c_Key3(8, 7, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(8, 7, -1), common_term_14);
    result_by_q.set(c_Key1(-1), common_term_14);
    // q = 0
    result_by_lpq.set(c_Key3(8, 7, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(8, 7, 1), common_term_12);
    result_by_q.set(c_Key1(1), common_term_12);
    // q = 2
    result_by_lpq.set(c_Key3(8, 7, 2), common_term_11);
    result_by_q.set(c_Key1(2), common_term_11);
    // q = 3
    result_by_lpq.set(c_Key3(8, 7, 3), common_term_10);
    result_by_q.set(c_Key1(3), common_term_10);
    // q = 4
    result_by_lpq.set(c_Key3(8, 7, 4), common_term_9);
    result_by_q.set(c_Key1(4), common_term_9);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 7), result_by_q);
    result_by_q.clear();

    // l , p = (8, 8).
    // q = -4
    result_by_lpq.set(c_Key3(8, 8, -4), common_term_8);
    result_by_q.set(c_Key1(-4), common_term_8);
    // q = -3
    result_by_lpq.set(c_Key3(8, 8, -3), common_term_7);
    result_by_q.set(c_Key1(-3), common_term_7);
    // q = -2
    result_by_lpq.set(c_Key3(8, 8, -2), common_term_6);
    result_by_q.set(c_Key1(-2), common_term_6);
    // q = -1
    result_by_lpq.set(c_Key3(8, 8, -1), common_term_5);
    result_by_q.set(c_Key1(-1), common_term_5);
    // q = 0
    result_by_lpq.set(c_Key3(8, 8, 0), common_term_4);
    result_by_q.set(c_Key1(0), common_term_4);
    // q = 1
    result_by_lpq.set(c_Key3(8, 8, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // q = 2
    result_by_lpq.set(c_Key3(8, 8, 2), common_term_2);
    result_by_q.set(c_Key1(2), common_term_2);
    // q = 3
    result_by_lpq.set(c_Key3(8, 8, 3), common_term_1);
    result_by_q.set(c_Key1(3), common_term_1);
    // q = 4
    result_by_lpq.set(c_Key3(8, 8, 4), common_term_0);
    result_by_q.set(c_Key1(4), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 8), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l8_e5(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 8.
    // Functions are taylor expanded and truncated at eccentricity^5.

    // Functions Calculated for l = 8.

    c_IntMap<c_Key3, double> result_by_lpq(99);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(9);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_5 = eccentricity_4 * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double common_term_0 = -0.06328125*eccentricity_5;
    double common_term_1 = 0.66666666666666667*eccentricity_4;
    double common_term_2 = 11.881510416666667*eccentricity_5 - 2.6041666666666667*eccentricity_3;
    double common_term_3 = -40.5*eccentricity_4 + 4.5*eccentricity_2;
    double common_term_4 = -343.30989583333333*eccentricity_5 + 63.4375*eccentricity_3 - 3.5*eccentricity;
    double common_term_5 = 490.75*eccentricity_4 - 46.0*eccentricity_2 + 1.0;
    double common_term_6 = 2757.0234375*eccentricity_5 - 333.5625*eccentricity_3 + 12.5*eccentricity;
    double common_term_7 = -1764.1666666666667*eccentricity_4 + 86.5*eccentricity_2;
    double common_term_8 = -7600.3268229166667*eccentricity_5 + 437.72916666666667*eccentricity_3;
    double common_term_9 = 1808.25*eccentricity_4;
    double common_term_10 = 6461.7950520833333*eccentricity_5;
    double common_term_11 = 0.08046875*eccentricity_5;
    double common_term_12 = 0.0625*eccentricity_4;
    double common_term_13 = 0.43359375*eccentricity_5 - 0.0625*eccentricity_3;
    double common_term_14 = -1.625*eccentricity_4 + 0.75*eccentricity_2;
    double common_term_15 = -17.3828125*eccentricity_5 + 10.3125*eccentricity_3 - 1.5*eccentricity;
    double common_term_16 = 78.1875*eccentricity_4 - 18.0*eccentricity_2 + 1.0;
    double common_term_17 = 430.1171875*eccentricity_5 - 119.4375*eccentricity_3 + 10.5*eccentricity;
    double common_term_18 = -580.375*eccentricity_4 + 62.25*eccentricity_2;
    double common_term_19 = -2307.59765625*eccentricity_5 + 274.1875*eccentricity_3;
    double common_term_20 = 998.75*eccentricity_4;
    double common_term_21 = 3181.34921875*eccentricity_5;
    double common_term_22 = 3.2138020833333333*eccentricity_5;
    double common_term_23 = 2.1875*eccentricity_4*std::pow(1.0 - eccentricity_2, -7.5);
    double common_term_24 = 12.141927083333333*eccentricity_5 + 1.4791666666666667*eccentricity_3;
    double common_term_25 = 8.5833333333333333*eccentricity_4 + eccentricity_2;
    double common_term_26 = 29.6484375*eccentricity_5 + 6.1875*eccentricity_3 + 0.5*eccentricity;
    double common_term_27 = 23.5*eccentricity_4 + 2.0*eccentricity_2 + 1.0;
    double common_term_28 = 73.815104166666667*eccentricity_5 - 0.3125*eccentricity_3 + 8.5*eccentricity;
    double common_term_29 = -37.25*eccentricity_4 + 42.0*eccentricity_2;
    double common_term_30 = -226.96223958333333*eccentricity_5 + 157.64583333333333*eccentricity_3;
    double common_term_31 = 497.97916666666667*eccentricity_4;
    double common_term_32 = 1395.14296875*eccentricity_5;
    double common_term_33 = 31.21171875*eccentricity_5;
    double common_term_34 = 18.041666666666667*eccentricity_4;
    double common_term_35 = 81.131510416666667*eccentricity_5 + 10.020833333333333*eccentricity_3;
    double common_term_36 = std::pow(1.0 - eccentricity_2, -7.5)*(8.75*eccentricity_4 + 5.25*eccentricity_2);
    double common_term_37 = 146.03385416666667*eccentricity_5 + 27.0625*eccentricity_3 + 2.5*eccentricity;
    double common_term_38 = 86.6875*eccentricity_4 + 14.0*eccentricity_2 + 1.0;
    double common_term_39 = 219.8671875*eccentricity_5 + 47.8125*eccentricity_3 + 6.5*eccentricity;
    double common_term_40 = 127.20833333333333*eccentricity_4 + 25.75*eccentricity_2;
    double common_term_41 = 287.45442708333333*eccentricity_5 + 80.104166666666667*eccentricity_3;
    double common_term_42 = 214.9375*eccentricity_4;
    double common_term_43 = 521.30130208333333*eccentricity_5;
    double common_term_44 = 153.94921875*eccentricity_5;
    double common_term_45 = 74.625*eccentricity_4;
    double common_term_46 = 225.52734375*eccentricity_5 + 33.5625*eccentricity_3;
    double common_term_47 = 111.0*eccentricity_4 + 13.5*eccentricity_2;
    double common_term_48 = 264.0234375*eccentricity_5 + 48.9375*eccentricity_3 + 4.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(11);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (8, 0).
    // q = -5
    result_by_lpq.set(c_Key3(8, 0, -5), common_term_0);
    result_by_q.set(c_Key1(-5), common_term_0);
    // q = -4
    result_by_lpq.set(c_Key3(8, 0, -4), common_term_1);
    result_by_q.set(c_Key1(-4), common_term_1);
    // q = -3
    result_by_lpq.set(c_Key3(8, 0, -3), common_term_2);
    result_by_q.set(c_Key1(-3), common_term_2);
    // q = -2
    result_by_lpq.set(c_Key3(8, 0, -2), common_term_3);
    result_by_q.set(c_Key1(-2), common_term_3);
    // q = -1
    result_by_lpq.set(c_Key3(8, 0, -1), common_term_4);
    result_by_q.set(c_Key1(-1), common_term_4);
    // q = 0
    result_by_lpq.set(c_Key3(8, 0, 0), common_term_5);
    result_by_q.set(c_Key1(0), common_term_5);
    // q = 1
    result_by_lpq.set(c_Key3(8, 0, 1), common_term_6);
    result_by_q.set(c_Key1(1), common_term_6);
    // q = 2
    result_by_lpq.set(c_Key3(8, 0, 2), common_term_7);
    result_by_q.set(c_Key1(2), common_term_7);
    // q = 3
    result_by_lpq.set(c_Key3(8, 0, 3), common_term_8);
    result_by_q.set(c_Key1(3), common_term_8);
    // q = 4
    result_by_lpq.set(c_Key3(8, 0, 4), common_term_9);
    result_by_q.set(c_Key1(4), common_term_9);
    // q = 5
    result_by_lpq.set(c_Key3(8, 0, 5), common_term_10);
    result_by_q.set(c_Key1(5), common_term_10);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 0), result_by_q);
    result_by_q.clear();

    // l , p = (8, 1).
    // q = -5
    result_by_lpq.set(c_Key3(8, 1, -5), common_term_11);
    result_by_q.set(c_Key1(-5), common_term_11);
    // q = -4
    result_by_lpq.set(c_Key3(8, 1, -4), common_term_12);
    result_by_q.set(c_Key1(-4), common_term_12);
    // q = -3
    result_by_lpq.set(c_Key3(8, 1, -3), common_term_13);
    result_by_q.set(c_Key1(-3), common_term_13);
    // q = -2
    result_by_lpq.set(c_Key3(8, 1, -2), common_term_14);
    result_by_q.set(c_Key1(-2), common_term_14);
    // q = -1
    result_by_lpq.set(c_Key3(8, 1, -1), common_term_15);
    result_by_q.set(c_Key1(-1), common_term_15);
    // q = 0
    result_by_lpq.set(c_Key3(8, 1, 0), common_term_16);
    result_by_q.set(c_Key1(0), common_term_16);
    // q = 1
    result_by_lpq.set(c_Key3(8, 1, 1), common_term_17);
    result_by_q.set(c_Key1(1), common_term_17);
    // q = 2
    result_by_lpq.set(c_Key3(8, 1, 2), common_term_18);
    result_by_q.set(c_Key1(2), common_term_18);
    // q = 3
    result_by_lpq.set(c_Key3(8, 1, 3), common_term_19);
    result_by_q.set(c_Key1(3), common_term_19);
    // q = 4
    result_by_lpq.set(c_Key3(8, 1, 4), common_term_20);
    result_by_q.set(c_Key1(4), common_term_20);
    // q = 5
    result_by_lpq.set(c_Key3(8, 1, 5), common_term_21);
    result_by_q.set(c_Key1(5), common_term_21);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 1), result_by_q);
    result_by_q.clear();

    // l , p = (8, 2).
    // q = -5
    result_by_lpq.set(c_Key3(8, 2, -5), common_term_22);
    result_by_q.set(c_Key1(-5), common_term_22);
    // q = -4
    result_by_lpq.set(c_Key3(8, 2, -4), common_term_23);
    result_by_q.set(c_Key1(-4), common_term_23);
    // q = -3
    result_by_lpq.set(c_Key3(8, 2, -3), common_term_24);
    result_by_q.set(c_Key1(-3), common_term_24);
    // q = -2
    result_by_lpq.set(c_Key3(8, 2, -2), common_term_25);
    result_by_q.set(c_Key1(-2), common_term_25);
    // q = -1
    result_by_lpq.set(c_Key3(8, 2, -1), common_term_26);
    result_by_q.set(c_Key1(-1), common_term_26);
    // q = 0
    result_by_lpq.set(c_Key3(8, 2, 0), common_term_27);
    result_by_q.set(c_Key1(0), common_term_27);
    // q = 1
    result_by_lpq.set(c_Key3(8, 2, 1), common_term_28);
    result_by_q.set(c_Key1(1), common_term_28);
    // q = 2
    result_by_lpq.set(c_Key3(8, 2, 2), common_term_29);
    result_by_q.set(c_Key1(2), common_term_29);
    // q = 3
    result_by_lpq.set(c_Key3(8, 2, 3), common_term_30);
    result_by_q.set(c_Key1(3), common_term_30);
    // q = 4
    result_by_lpq.set(c_Key3(8, 2, 4), common_term_31);
    result_by_q.set(c_Key1(4), common_term_31);
    // q = 5
    result_by_lpq.set(c_Key3(8, 2, 5), common_term_32);
    result_by_q.set(c_Key1(5), common_term_32);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 2), result_by_q);
    result_by_q.clear();

    // l , p = (8, 3).
    // q = -5
    result_by_lpq.set(c_Key3(8, 3, -5), common_term_33);
    result_by_q.set(c_Key1(-5), common_term_33);
    // q = -4
    result_by_lpq.set(c_Key3(8, 3, -4), common_term_34);
    result_by_q.set(c_Key1(-4), common_term_34);
    // q = -3
    result_by_lpq.set(c_Key3(8, 3, -3), common_term_35);
    result_by_q.set(c_Key1(-3), common_term_35);
    // q = -2
    result_by_lpq.set(c_Key3(8, 3, -2), common_term_36);
    result_by_q.set(c_Key1(-2), common_term_36);
    // q = -1
    result_by_lpq.set(c_Key3(8, 3, -1), common_term_37);
    result_by_q.set(c_Key1(-1), common_term_37);
    // q = 0
    result_by_lpq.set(c_Key3(8, 3, 0), common_term_38);
    result_by_q.set(c_Key1(0), common_term_38);
    // q = 1
    result_by_lpq.set(c_Key3(8, 3, 1), common_term_39);
    result_by_q.set(c_Key1(1), common_term_39);
    // q = 2
    result_by_lpq.set(c_Key3(8, 3, 2), common_term_40);
    result_by_q.set(c_Key1(2), common_term_40);
    // q = 3
    result_by_lpq.set(c_Key3(8, 3, 3), common_term_41);
    result_by_q.set(c_Key1(3), common_term_41);
    // q = 4
    result_by_lpq.set(c_Key3(8, 3, 4), common_term_42);
    result_by_q.set(c_Key1(4), common_term_42);
    // q = 5
    result_by_lpq.set(c_Key3(8, 3, 5), common_term_43);
    result_by_q.set(c_Key1(5), common_term_43);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 3), result_by_q);
    result_by_q.clear();

    // l , p = (8, 4).
    // q = -5
    result_by_lpq.set(c_Key3(8, 4, -5), common_term_44);
    result_by_q.set(c_Key1(-5), common_term_44);
    // q = -4
    result_by_lpq.set(c_Key3(8, 4, -4), common_term_45);
    result_by_q.set(c_Key1(-4), common_term_45);
    // q = -3
    result_by_lpq.set(c_Key3(8, 4, -3), common_term_46);
    result_by_q.set(c_Key1(-3), common_term_46);
    // q = -2
    result_by_lpq.set(c_Key3(8, 4, -2), common_term_47);
    result_by_q.set(c_Key1(-2), common_term_47);
    // q = -1
    result_by_lpq.set(c_Key3(8, 4, -1), common_term_48);
    result_by_q.set(c_Key1(-1), common_term_48);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -7.5)*(13.125*eccentricity_4 + 10.5*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(8, 4, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(8, 4, 1), common_term_48);
    result_by_q.set(c_Key1(1), common_term_48);
    // q = 2
    result_by_lpq.set(c_Key3(8, 4, 2), common_term_47);
    result_by_q.set(c_Key1(2), common_term_47);
    // q = 3
    result_by_lpq.set(c_Key3(8, 4, 3), common_term_46);
    result_by_q.set(c_Key1(3), common_term_46);
    // q = 4
    result_by_lpq.set(c_Key3(8, 4, 4), common_term_45);
    result_by_q.set(c_Key1(4), common_term_45);
    // q = 5
    result_by_lpq.set(c_Key3(8, 4, 5), common_term_44);
    result_by_q.set(c_Key1(5), common_term_44);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 4), result_by_q);
    result_by_q.clear();

    // l , p = (8, 5).
    // q = -5
    result_by_lpq.set(c_Key3(8, 5, -5), common_term_43);
    result_by_q.set(c_Key1(-5), common_term_43);
    // q = -4
    result_by_lpq.set(c_Key3(8, 5, -4), common_term_42);
    result_by_q.set(c_Key1(-4), common_term_42);
    // q = -3
    result_by_lpq.set(c_Key3(8, 5, -3), common_term_41);
    result_by_q.set(c_Key1(-3), common_term_41);
    // q = -2
    result_by_lpq.set(c_Key3(8, 5, -2), common_term_40);
    result_by_q.set(c_Key1(-2), common_term_40);
    // q = -1
    result_by_lpq.set(c_Key3(8, 5, -1), common_term_39);
    result_by_q.set(c_Key1(-1), common_term_39);
    // q = 0
    result_by_lpq.set(c_Key3(8, 5, 0), common_term_38);
    result_by_q.set(c_Key1(0), common_term_38);
    // q = 1
    result_by_lpq.set(c_Key3(8, 5, 1), common_term_37);
    result_by_q.set(c_Key1(1), common_term_37);
    // q = 2
    result_by_lpq.set(c_Key3(8, 5, 2), common_term_36);
    result_by_q.set(c_Key1(2), common_term_36);
    // q = 3
    result_by_lpq.set(c_Key3(8, 5, 3), common_term_35);
    result_by_q.set(c_Key1(3), common_term_35);
    // q = 4
    result_by_lpq.set(c_Key3(8, 5, 4), common_term_34);
    result_by_q.set(c_Key1(4), common_term_34);
    // q = 5
    result_by_lpq.set(c_Key3(8, 5, 5), common_term_33);
    result_by_q.set(c_Key1(5), common_term_33);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 5), result_by_q);
    result_by_q.clear();

    // l , p = (8, 6).
    // q = -5
    result_by_lpq.set(c_Key3(8, 6, -5), common_term_32);
    result_by_q.set(c_Key1(-5), common_term_32);
    // q = -4
    result_by_lpq.set(c_Key3(8, 6, -4), common_term_31);
    result_by_q.set(c_Key1(-4), common_term_31);
    // q = -3
    result_by_lpq.set(c_Key3(8, 6, -3), common_term_30);
    result_by_q.set(c_Key1(-3), common_term_30);
    // q = -2
    result_by_lpq.set(c_Key3(8, 6, -2), common_term_29);
    result_by_q.set(c_Key1(-2), common_term_29);
    // q = -1
    result_by_lpq.set(c_Key3(8, 6, -1), common_term_28);
    result_by_q.set(c_Key1(-1), common_term_28);
    // q = 0
    result_by_lpq.set(c_Key3(8, 6, 0), common_term_27);
    result_by_q.set(c_Key1(0), common_term_27);
    // q = 1
    result_by_lpq.set(c_Key3(8, 6, 1), common_term_26);
    result_by_q.set(c_Key1(1), common_term_26);
    // q = 2
    result_by_lpq.set(c_Key3(8, 6, 2), common_term_25);
    result_by_q.set(c_Key1(2), common_term_25);
    // q = 3
    result_by_lpq.set(c_Key3(8, 6, 3), common_term_24);
    result_by_q.set(c_Key1(3), common_term_24);
    // q = 4
    result_by_lpq.set(c_Key3(8, 6, 4), common_term_23);
    result_by_q.set(c_Key1(4), common_term_23);
    // q = 5
    result_by_lpq.set(c_Key3(8, 6, 5), common_term_22);
    result_by_q.set(c_Key1(5), common_term_22);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 6), result_by_q);
    result_by_q.clear();

    // l , p = (8, 7).
    // q = -5
    result_by_lpq.set(c_Key3(8, 7, -5), common_term_21);
    result_by_q.set(c_Key1(-5), common_term_21);
    // q = -4
    result_by_lpq.set(c_Key3(8, 7, -4), common_term_20);
    result_by_q.set(c_Key1(-4), common_term_20);
    // q = -3
    result_by_lpq.set(c_Key3(8, 7, -3), common_term_19);
    result_by_q.set(c_Key1(-3), common_term_19);
    // q = -2
    result_by_lpq.set(c_Key3(8, 7, -2), common_term_18);
    result_by_q.set(c_Key1(-2), common_term_18);
    // q = -1
    result_by_lpq.set(c_Key3(8, 7, -1), common_term_17);
    result_by_q.set(c_Key1(-1), common_term_17);
    // q = 0
    result_by_lpq.set(c_Key3(8, 7, 0), common_term_16);
    result_by_q.set(c_Key1(0), common_term_16);
    // q = 1
    result_by_lpq.set(c_Key3(8, 7, 1), common_term_15);
    result_by_q.set(c_Key1(1), common_term_15);
    // q = 2
    result_by_lpq.set(c_Key3(8, 7, 2), common_term_14);
    result_by_q.set(c_Key1(2), common_term_14);
    // q = 3
    result_by_lpq.set(c_Key3(8, 7, 3), common_term_13);
    result_by_q.set(c_Key1(3), common_term_13);
    // q = 4
    result_by_lpq.set(c_Key3(8, 7, 4), common_term_12);
    result_by_q.set(c_Key1(4), common_term_12);
    // q = 5
    result_by_lpq.set(c_Key3(8, 7, 5), common_term_11);
    result_by_q.set(c_Key1(5), common_term_11);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 7), result_by_q);
    result_by_q.clear();

    // l , p = (8, 8).
    // q = -5
    result_by_lpq.set(c_Key3(8, 8, -5), common_term_10);
    result_by_q.set(c_Key1(-5), common_term_10);
    // q = -4
    result_by_lpq.set(c_Key3(8, 8, -4), common_term_9);
    result_by_q.set(c_Key1(-4), common_term_9);
    // q = -3
    result_by_lpq.set(c_Key3(8, 8, -3), common_term_8);
    result_by_q.set(c_Key1(-3), common_term_8);
    // q = -2
    result_by_lpq.set(c_Key3(8, 8, -2), common_term_7);
    result_by_q.set(c_Key1(-2), common_term_7);
    // q = -1
    result_by_lpq.set(c_Key3(8, 8, -1), common_term_6);
    result_by_q.set(c_Key1(-1), common_term_6);
    // q = 0
    result_by_lpq.set(c_Key3(8, 8, 0), common_term_5);
    result_by_q.set(c_Key1(0), common_term_5);
    // q = 1
    result_by_lpq.set(c_Key3(8, 8, 1), common_term_4);
    result_by_q.set(c_Key1(1), common_term_4);
    // q = 2
    result_by_lpq.set(c_Key3(8, 8, 2), common_term_3);
    result_by_q.set(c_Key1(2), common_term_3);
    // q = 3
    result_by_lpq.set(c_Key3(8, 8, 3), common_term_2);
    result_by_q.set(c_Key1(3), common_term_2);
    // q = 4
    result_by_lpq.set(c_Key3(8, 8, 4), common_term_1);
    result_by_q.set(c_Key1(4), common_term_1);
    // q = 5
    result_by_lpq.set(c_Key3(8, 8, 5), common_term_0);
    result_by_q.set(c_Key1(5), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 8), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l8_e10(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 8.
    // Functions are taylor expanded and truncated at eccentricity^10.

    // Functions Calculated for l = 8.

    c_IntMap<c_Key3, double> result_by_lpq(187);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(9);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_5 = eccentricity_4 * eccentricity;
    double eccentricity_10 = eccentricity_5 * eccentricity_5;
    double eccentricity_8 = eccentricity_4 * eccentricity_4;
    double eccentricity_9 = eccentricity_8 * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double eccentricity_6 = eccentricity_3 * eccentricity_3;
    double eccentricity_7 = eccentricity_6 * eccentricity;
    double common_term_0 = 2.7557319223985891e-7*eccentricity_10;
    double common_term_1 = 5.3822889109347443e-9*eccentricity_9;
    double common_term_2 = -2.2767082093253968e-6*eccentricity_9 - 1.5500992063492063e-6*eccentricity_7;
    double common_term_3 = 0.00090525793650793651*eccentricity_10 + 0.00099206349206349206*eccentricity_8 + 0.0013888888888888889*eccentricity_6;
    double common_term_4 = -0.0015537806919642857*eccentricity_9 + 0.02373046875*eccentricity_7 - 0.06328125*eccentricity_5;
    double common_term_5 = -0.14656084656084656*eccentricity_10 + 0.58888888888888889*eccentricity_8 - 1.3333333333333333*eccentricity_6 + 0.66666666666666667*eccentricity_4;
    double common_term_6 = 8.5598980938946759*eccentricity_9 - 16.215006510416667*eccentricity_7 + 11.881510416666667*eccentricity_5 - 2.6041666666666667*eccentricity_3;
    double common_term_7 = 83.11640625*eccentricity_10 - 139.6125*eccentricity_8 + 115.3125*eccentricity_6 - 40.5*eccentricity_4 + 4.5*eccentricity_2;
    double common_term_8 = -935.56048651801215*eccentricity_9 + 800.71652560763889*eccentricity_7 - 343.30989583333333*eccentricity_5 + 63.4375*eccentricity_3 - 3.5*eccentricity;
    double common_term_9 = -5182.6406944444444*eccentricity_10 + 4439.1657986111111*eccentricity_8 - 2099.8611111111111*eccentricity_6 + 490.75*eccentricity_4 - 46.0*eccentricity_2 + 1.0;
    double common_term_10 = 20854.591314697266*eccentricity_9 - 10362.46142578125*eccentricity_7 + 2757.0234375*eccentricity_5 - 333.5625*eccentricity_3 + 12.5*eccentricity;
    double common_term_11 = 86174.851707175926*eccentricity_10 - 43779.236111111111*eccentricity_8 + 12571.145833333333*eccentricity_6 - 1764.1666666666667*eccentricity_4 + 86.5*eccentricity_2;
    double common_term_12 = -164205.81790861907*eccentricity_9 + 49339.254264322917*eccentricity_7 - 7600.3268229166667*eccentricity_5 + 437.72916666666667*eccentricity_3;
    double common_term_13 = -560229.34285714286*eccentricity_10 + 172758.65625*eccentricity_8 - 28273.05*eccentricity_6 + 1808.25*eccentricity_4;
    double common_term_14 = 552633.89492439391*eccentricity_9 - 94090.807693142361*eccentricity_7 + 6461.7950520833333*eccentricity_5;
    double common_term_15 = 1642304.3425347222*eccentricity_10 - 286780.58888888889*eccentricity_8 + 20690.784722222222*eccentricity_6;
    double common_term_16 = -813848.52869306292*eccentricity_9 + 60764.489606584821*eccentricity_7;
    double common_term_17 = -2176625.8677358907*eccentricity_10 + 166371.02557043651*eccentricity_8;
    double common_term_18 = 429802.86250706695*eccentricity_9;
    double common_term_19 = 1057233.6862946429*eccentricity_10;
    double common_term_20 = 0.37329840443121693*eccentricity_10;
    double common_term_21 = 0.27462201799665179*eccentricity_9;
    double common_term_22 = 1.4143146494708995*eccentricity_10 + 0.20204613095238095*eccentricity_8;
    double common_term_23 = 1.0777563185918899*eccentricity_9 + 0.14865606398809524*eccentricity_7;
    double common_term_24 = 0.109375*eccentricity_6*std::pow(1.0 - eccentricity_2, -7.5);
    double common_term_25 = 2.7210094633556548*eccentricity_9 + 0.62366536458333333*eccentricity_7 + 0.08046875*eccentricity_5;
    double common_term_26 = 7.0086185515873016*eccentricity_10 + 2.1216145833333333*eccentricity_8 + 0.475*eccentricity_6 + 0.0625*eccentricity_4;
    double common_term_27 = 5.57567138671875*eccentricity_9 + 1.64326171875*eccentricity_7 + 0.43359375*eccentricity_5 - 0.0625*eccentricity_3;
    double common_term_28 = 12.522894965277778*eccentricity_10 + 4.0916666666666667*eccentricity_8 + 2.5*eccentricity_6 - 1.625*eccentricity_4 + 0.75*eccentricity_2;
    double common_term_29 = 5.5236307779947917*eccentricity_9 + 16.171061197916667*eccentricity_7 - 17.3828125*eccentricity_5 + 10.3125*eccentricity_3 - 1.5*eccentricity;
    double common_term_30 = -19.72640625*eccentricity_10 + 101.90625*eccentricity_8 - 122.3125*eccentricity_6 + 78.1875*eccentricity_4 - 18.0*eccentricity_2 + 1.0;
    double common_term_31 = 560.69659627278646*eccentricity_9 - 662.57438151041667*eccentricity_7 + 430.1171875*eccentricity_5 - 119.4375*eccentricity_3 + 10.5*eccentricity;
    double common_term_32 = 2661.0693576388889*eccentricity_10 - 2985.2182291666667*eccentricity_8 + 1915.953125*eccentricity_6 - 580.375*eccentricity_4 + 62.25*eccentricity_2;
    double common_term_33 = -11708.047570800781*eccentricity_9 + 7329.28740234375*eccentricity_7 - 2307.59765625*eccentricity_5 + 274.1875*eccentricity_3;
    double common_term_34 = -41181.945374503968*eccentricity_10 + 24965.729166666667*eccentricity_8 - 7957.5625*eccentricity_6 + 998.75*eccentricity_4;
    double common_term_35 = 77582.160425385975*eccentricity_9 - 24651.222298177083*eccentricity_7 + 3181.34921875*eccentricity_5;
    double common_term_36 = 223753.60797991071*eccentricity_10 - 70208.108035714286*eccentricity_8 + 9163.8375*eccentricity_6;
    double common_term_37 = -186830.13835812523*eccentricity_9 + 24400.922632998512*eccentricity_7;
    double common_term_38 = -470044.07146990741*eccentricity_10 + 60988.455208333333*eccentricity_8;
    double common_term_39 = 144690.15621730259*eccentricity_9;
    double common_term_40 = 328564.26045841601*eccentricity_10;
    double common_term_41 = 20.492165178571429*eccentricity_10;
    double common_term_42 = 14.253836199509824*eccentricity_9;
    double common_term_43 = 62.845375881834215*eccentricity_10 + 9.8822296626984127*eccentricity_8;
    double common_term_44 = 45.722182355608259*eccentricity_9 + 6.8260463169642857*eccentricity_7;
    double common_term_45 = 134.20566716269841*eccentricity_10 + 33.096825396825397*eccentricity_8 + 4.6951388888888889*eccentricity_6;
    double common_term_46 = 100.35344373914931*eccentricity_9 + 23.830978732638889*eccentricity_7 + 3.2138020833333333*eccentricity_5;
    double common_term_47 = std::pow(1.0 - eccentricity_2, -7.5)*(0.65625*eccentricity_6 + 2.1875*eccentricity_4);
    double common_term_48 = 184.36379711009838*eccentricity_9 + 55.222623697916667*eccentricity_7 + 12.141927083333333*eccentricity_5 + 1.4791666666666667*eccentricity_3;
    double common_term_49 = 392.9662037037037*eccentricity_10 + 139.89722222222222*eccentricity_8 + 40.614583333333333*eccentricity_6 + 8.5833333333333333*eccentricity_4 + eccentricity_2;
    double common_term_50 = 304.38755493164062*eccentricity_9 + 105.59130859375*eccentricity_7 + 29.6484375*eccentricity_5 + 6.1875*eccentricity_3 + 0.5*eccentricity;
    double common_term_51 = 595.95076388888889*eccentricity_10 + 234.75564236111111*eccentricity_8 + 78.555555555555556*eccentricity_6 + 23.5*eccentricity_4 + 2.0*eccentricity_2 + 1.0;
    double common_term_52 = 469.22595553927951*eccentricity_9 + 173.39816623263889*eccentricity_7 + 73.815104166666667*eccentricity_5 - 0.3125*eccentricity_3 + 8.5*eccentricity;
    double common_term_53 = 873.31640625*eccentricity_10 + 321.95625*eccentricity_8 + 218.90625*eccentricity_6 - 37.25*eccentricity_4 + 42.0*eccentricity_2;
    double common_term_54 = 458.83704653139468*eccentricity_9 + 648.71975911458333*eccentricity_7 - 226.96223958333333*eccentricity_5 + 157.64583333333333*eccentricity_3;
    double common_term_55 = 232.20453042328042*eccentricity_10 + 1927.3425347222222*eccentricity_8 - 932.79583333333333*eccentricity_6 + 497.97916666666667*eccentricity_4;
    double common_term_56 = 5629.1753435407366*eccentricity_9 - 3134.50556640625*eccentricity_7 + 1395.14296875*eccentricity_5;
    double common_term_57 = 15875.237475198413*eccentricity_10 - 9247.1106150793651*eccentricity_8 + 3575.9618055555556*eccentricity_6;
    double common_term_58 = -24842.062921966068*eccentricity_9 + 8555.9166837177579*eccentricity_7;
    double common_term_59 = -62129.775892857143*eccentricity_10 + 19375.531473214286*eccentricity_8;
    double common_term_60 = 41945.743762018651*eccentricity_9;
    double common_term_61 = 87458.838555169753*eccentricity_10;
    double common_term_62 = 347.53703627921076*eccentricity_10;
    double common_term_63 = 221.16175382637683*eccentricity_9;
    double common_term_64 = 717.35178571428571*eccentricity_10 + 139.04497767857143*eccentricity_8;
    double common_term_65 = 483.05924309624566*eccentricity_9 + 86.155606708829365*eccentricity_7;
    double common_term_66 = 1181.0311042906746*eccentricity_10 + 319.79875992063492*eccentricity_8 + 52.440972222222222*eccentricity_6;
    double common_term_67 = 813.49096854073661*eccentricity_9 + 207.61494140625*eccentricity_7 + 31.21171875*eccentricity_5;
    double common_term_68 = 1721.8277571097884*eccentricity_10 + 550.18836805555556*eccentricity_8 + 131.67916666666667*eccentricity_6 + 18.041666666666667*eccentricity_4;
    double common_term_69 = 1200.9678082501447*eccentricity_9 + 364.11526692708333*eccentricity_7 + 81.131510416666667*eccentricity_5 + 10.020833333333333*eccentricity_3;
    double common_term_70 = std::pow(1.0 - eccentricity_2, -7.5)*(1.640625*eccentricity_6 + 8.75*eccentricity_4 + 5.25*eccentricity_2);
    double common_term_71 = 1633.5173536512587*eccentricity_9 + 547.90695529513889*eccentricity_7 + 146.03385416666667*eccentricity_5 + 27.0625*eccentricity_3 + 2.5*eccentricity;
    double common_term_72 = 2979.5158159722222*eccentricity_10 + 1121.8220486111111*eccentricity_8 + 354.40972222222222*eccentricity_6 + 86.6875*eccentricity_4 + 14.0*eccentricity_2 + 1.0;
    double common_term_73 = 2096.8638977050781*eccentricity_9 + 749.23486328125*eccentricity_7 + 219.8671875*eccentricity_5 + 47.8125*eccentricity_3 + 6.5*eccentricity;
    double common_term_74 = 3662.7114945023148*eccentricity_10 + 1439.2086805555556*eccentricity_8 + 482.55208333333333*eccentricity_6 + 127.20833333333333*eccentricity_4 + 25.75*eccentricity_2;
    double common_term_75 = 2573.966889558015*eccentricity_9 + 958.58414713541667*eccentricity_7 + 287.45442708333333*eccentricity_5 + 80.104166666666667*eccentricity_3;
    double common_term_76 = 4352.0655133928571*eccentricity_10 + 1772.05546875*eccentricity_8 + 573.0*eccentricity_6 + 214.9375*eccentricity_4;
    double common_term_77 = 3111.831438530816*eccentricity_9 + 1021.5815646701389*eccentricity_7 + 521.30130208333333*eccentricity_5;
    double common_term_78 = 5286.5708922371032*eccentricity_10 + 1618.5184771825397*eccentricity_8 + 1174.3086805555556*eccentricity_6;
    double common_term_79 = 2188.9004241943359*eccentricity_9 + 2499.8066545758929*eccentricity_7;
    double common_term_80 = 2167.5672054122575*eccentricity_10 + 5088.2304067460317*eccentricity_8;
    double common_term_81 = 9986.2402226469925*eccentricity_9;
    double common_term_82 = 19015.278426339286*eccentricity_10;
    double common_term_83 = 3137.8065166170635*eccentricity_10;
    double common_term_84 = 1812.6737448556083*eccentricity_9;
    double common_term_85 = 3283.9923115079365*eccentricity_10 + 1024.7152901785714*eccentricity_8;
    double common_term_86 = 2087.4189422607422*eccentricity_9 + 564.30244140625*eccentricity_7;
    double common_term_87 = 3957.3959263392857*eccentricity_10 + 1280.1776785714286*eccentricity_8 + 300.8625*eccentricity_6;
    double common_term_88 = 2515.3772408621652*eccentricity_9 + 754.56982421875*eccentricity_7 + 153.94921875*eccentricity_5;
    double common_term_89 = 4414.2023809523809*eccentricity_10 + 1544.878125*eccentricity_8 + 424.575*eccentricity_6 + 74.625*eccentricity_4;
    double common_term_90 = 2798.1171752929687*eccentricity_9 + 909.78837890625*eccentricity_7 + 225.52734375*eccentricity_5 + 33.5625*eccentricity_3;
    double common_term_91 = 4692.5872395833333*eccentricity_10 + 1705.44375*eccentricity_8 + 507.9375*eccentricity_6 + 111.0*eccentricity_4 + 13.5*eccentricity_2;
    double common_term_92 = 2941.4128967285156*eccentricity_9 + 988.83935546875*eccentricity_7 + 264.0234375*eccentricity_5 + 48.9375*eccentricity_3 + 4.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(21);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (8, 0).
    // q = -10
    result_by_lpq.set(c_Key3(8, 0, -10), common_term_0);
    result_by_q.set(c_Key1(-10), common_term_0);
    // q = -9
    result_by_lpq.set(c_Key3(8, 0, -9), common_term_1);
    result_by_q.set(c_Key1(-9), common_term_1);
    // q = -7
    result_by_lpq.set(c_Key3(8, 0, -7), common_term_2);
    result_by_q.set(c_Key1(-7), common_term_2);
    // q = -6
    result_by_lpq.set(c_Key3(8, 0, -6), common_term_3);
    result_by_q.set(c_Key1(-6), common_term_3);
    // q = -5
    result_by_lpq.set(c_Key3(8, 0, -5), common_term_4);
    result_by_q.set(c_Key1(-5), common_term_4);
    // q = -4
    result_by_lpq.set(c_Key3(8, 0, -4), common_term_5);
    result_by_q.set(c_Key1(-4), common_term_5);
    // q = -3
    result_by_lpq.set(c_Key3(8, 0, -3), common_term_6);
    result_by_q.set(c_Key1(-3), common_term_6);
    // q = -2
    result_by_lpq.set(c_Key3(8, 0, -2), common_term_7);
    result_by_q.set(c_Key1(-2), common_term_7);
    // q = -1
    result_by_lpq.set(c_Key3(8, 0, -1), common_term_8);
    result_by_q.set(c_Key1(-1), common_term_8);
    // q = 0
    result_by_lpq.set(c_Key3(8, 0, 0), common_term_9);
    result_by_q.set(c_Key1(0), common_term_9);
    // q = 1
    result_by_lpq.set(c_Key3(8, 0, 1), common_term_10);
    result_by_q.set(c_Key1(1), common_term_10);
    // q = 2
    result_by_lpq.set(c_Key3(8, 0, 2), common_term_11);
    result_by_q.set(c_Key1(2), common_term_11);
    // q = 3
    result_by_lpq.set(c_Key3(8, 0, 3), common_term_12);
    result_by_q.set(c_Key1(3), common_term_12);
    // q = 4
    result_by_lpq.set(c_Key3(8, 0, 4), common_term_13);
    result_by_q.set(c_Key1(4), common_term_13);
    // q = 5
    result_by_lpq.set(c_Key3(8, 0, 5), common_term_14);
    result_by_q.set(c_Key1(5), common_term_14);
    // q = 6
    result_by_lpq.set(c_Key3(8, 0, 6), common_term_15);
    result_by_q.set(c_Key1(6), common_term_15);
    // q = 7
    result_by_lpq.set(c_Key3(8, 0, 7), common_term_16);
    result_by_q.set(c_Key1(7), common_term_16);
    // q = 8
    result_by_lpq.set(c_Key3(8, 0, 8), common_term_17);
    result_by_q.set(c_Key1(8), common_term_17);
    // q = 9
    result_by_lpq.set(c_Key3(8, 0, 9), common_term_18);
    result_by_q.set(c_Key1(9), common_term_18);
    // q = 10
    result_by_lpq.set(c_Key3(8, 0, 10), common_term_19);
    result_by_q.set(c_Key1(10), common_term_19);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 0), result_by_q);
    result_by_q.clear();

    // l , p = (8, 1).
    // q = -10
    result_by_lpq.set(c_Key3(8, 1, -10), common_term_20);
    result_by_q.set(c_Key1(-10), common_term_20);
    // q = -9
    result_by_lpq.set(c_Key3(8, 1, -9), common_term_21);
    result_by_q.set(c_Key1(-9), common_term_21);
    // q = -8
    result_by_lpq.set(c_Key3(8, 1, -8), common_term_22);
    result_by_q.set(c_Key1(-8), common_term_22);
    // q = -7
    result_by_lpq.set(c_Key3(8, 1, -7), common_term_23);
    result_by_q.set(c_Key1(-7), common_term_23);
    // q = -6
    result_by_lpq.set(c_Key3(8, 1, -6), common_term_24);
    result_by_q.set(c_Key1(-6), common_term_24);
    // q = -5
    result_by_lpq.set(c_Key3(8, 1, -5), common_term_25);
    result_by_q.set(c_Key1(-5), common_term_25);
    // q = -4
    result_by_lpq.set(c_Key3(8, 1, -4), common_term_26);
    result_by_q.set(c_Key1(-4), common_term_26);
    // q = -3
    result_by_lpq.set(c_Key3(8, 1, -3), common_term_27);
    result_by_q.set(c_Key1(-3), common_term_27);
    // q = -2
    result_by_lpq.set(c_Key3(8, 1, -2), common_term_28);
    result_by_q.set(c_Key1(-2), common_term_28);
    // q = -1
    result_by_lpq.set(c_Key3(8, 1, -1), common_term_29);
    result_by_q.set(c_Key1(-1), common_term_29);
    // q = 0
    result_by_lpq.set(c_Key3(8, 1, 0), common_term_30);
    result_by_q.set(c_Key1(0), common_term_30);
    // q = 1
    result_by_lpq.set(c_Key3(8, 1, 1), common_term_31);
    result_by_q.set(c_Key1(1), common_term_31);
    // q = 2
    result_by_lpq.set(c_Key3(8, 1, 2), common_term_32);
    result_by_q.set(c_Key1(2), common_term_32);
    // q = 3
    result_by_lpq.set(c_Key3(8, 1, 3), common_term_33);
    result_by_q.set(c_Key1(3), common_term_33);
    // q = 4
    result_by_lpq.set(c_Key3(8, 1, 4), common_term_34);
    result_by_q.set(c_Key1(4), common_term_34);
    // q = 5
    result_by_lpq.set(c_Key3(8, 1, 5), common_term_35);
    result_by_q.set(c_Key1(5), common_term_35);
    // q = 6
    result_by_lpq.set(c_Key3(8, 1, 6), common_term_36);
    result_by_q.set(c_Key1(6), common_term_36);
    // q = 7
    result_by_lpq.set(c_Key3(8, 1, 7), common_term_37);
    result_by_q.set(c_Key1(7), common_term_37);
    // q = 8
    result_by_lpq.set(c_Key3(8, 1, 8), common_term_38);
    result_by_q.set(c_Key1(8), common_term_38);
    // q = 9
    result_by_lpq.set(c_Key3(8, 1, 9), common_term_39);
    result_by_q.set(c_Key1(9), common_term_39);
    // q = 10
    result_by_lpq.set(c_Key3(8, 1, 10), common_term_40);
    result_by_q.set(c_Key1(10), common_term_40);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 1), result_by_q);
    result_by_q.clear();

    // l , p = (8, 2).
    // q = -10
    result_by_lpq.set(c_Key3(8, 2, -10), common_term_41);
    result_by_q.set(c_Key1(-10), common_term_41);
    // q = -9
    result_by_lpq.set(c_Key3(8, 2, -9), common_term_42);
    result_by_q.set(c_Key1(-9), common_term_42);
    // q = -8
    result_by_lpq.set(c_Key3(8, 2, -8), common_term_43);
    result_by_q.set(c_Key1(-8), common_term_43);
    // q = -7
    result_by_lpq.set(c_Key3(8, 2, -7), common_term_44);
    result_by_q.set(c_Key1(-7), common_term_44);
    // q = -6
    result_by_lpq.set(c_Key3(8, 2, -6), common_term_45);
    result_by_q.set(c_Key1(-6), common_term_45);
    // q = -5
    result_by_lpq.set(c_Key3(8, 2, -5), common_term_46);
    result_by_q.set(c_Key1(-5), common_term_46);
    // q = -4
    result_by_lpq.set(c_Key3(8, 2, -4), common_term_47);
    result_by_q.set(c_Key1(-4), common_term_47);
    // q = -3
    result_by_lpq.set(c_Key3(8, 2, -3), common_term_48);
    result_by_q.set(c_Key1(-3), common_term_48);
    // q = -2
    result_by_lpq.set(c_Key3(8, 2, -2), common_term_49);
    result_by_q.set(c_Key1(-2), common_term_49);
    // q = -1
    result_by_lpq.set(c_Key3(8, 2, -1), common_term_50);
    result_by_q.set(c_Key1(-1), common_term_50);
    // q = 0
    result_by_lpq.set(c_Key3(8, 2, 0), common_term_51);
    result_by_q.set(c_Key1(0), common_term_51);
    // q = 1
    result_by_lpq.set(c_Key3(8, 2, 1), common_term_52);
    result_by_q.set(c_Key1(1), common_term_52);
    // q = 2
    result_by_lpq.set(c_Key3(8, 2, 2), common_term_53);
    result_by_q.set(c_Key1(2), common_term_53);
    // q = 3
    result_by_lpq.set(c_Key3(8, 2, 3), common_term_54);
    result_by_q.set(c_Key1(3), common_term_54);
    // q = 4
    result_by_lpq.set(c_Key3(8, 2, 4), common_term_55);
    result_by_q.set(c_Key1(4), common_term_55);
    // q = 5
    result_by_lpq.set(c_Key3(8, 2, 5), common_term_56);
    result_by_q.set(c_Key1(5), common_term_56);
    // q = 6
    result_by_lpq.set(c_Key3(8, 2, 6), common_term_57);
    result_by_q.set(c_Key1(6), common_term_57);
    // q = 7
    result_by_lpq.set(c_Key3(8, 2, 7), common_term_58);
    result_by_q.set(c_Key1(7), common_term_58);
    // q = 8
    result_by_lpq.set(c_Key3(8, 2, 8), common_term_59);
    result_by_q.set(c_Key1(8), common_term_59);
    // q = 9
    result_by_lpq.set(c_Key3(8, 2, 9), common_term_60);
    result_by_q.set(c_Key1(9), common_term_60);
    // q = 10
    result_by_lpq.set(c_Key3(8, 2, 10), common_term_61);
    result_by_q.set(c_Key1(10), common_term_61);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 2), result_by_q);
    result_by_q.clear();

    // l , p = (8, 3).
    // q = -10
    result_by_lpq.set(c_Key3(8, 3, -10), common_term_62);
    result_by_q.set(c_Key1(-10), common_term_62);
    // q = -9
    result_by_lpq.set(c_Key3(8, 3, -9), common_term_63);
    result_by_q.set(c_Key1(-9), common_term_63);
    // q = -8
    result_by_lpq.set(c_Key3(8, 3, -8), common_term_64);
    result_by_q.set(c_Key1(-8), common_term_64);
    // q = -7
    result_by_lpq.set(c_Key3(8, 3, -7), common_term_65);
    result_by_q.set(c_Key1(-7), common_term_65);
    // q = -6
    result_by_lpq.set(c_Key3(8, 3, -6), common_term_66);
    result_by_q.set(c_Key1(-6), common_term_66);
    // q = -5
    result_by_lpq.set(c_Key3(8, 3, -5), common_term_67);
    result_by_q.set(c_Key1(-5), common_term_67);
    // q = -4
    result_by_lpq.set(c_Key3(8, 3, -4), common_term_68);
    result_by_q.set(c_Key1(-4), common_term_68);
    // q = -3
    result_by_lpq.set(c_Key3(8, 3, -3), common_term_69);
    result_by_q.set(c_Key1(-3), common_term_69);
    // q = -2
    result_by_lpq.set(c_Key3(8, 3, -2), common_term_70);
    result_by_q.set(c_Key1(-2), common_term_70);
    // q = -1
    result_by_lpq.set(c_Key3(8, 3, -1), common_term_71);
    result_by_q.set(c_Key1(-1), common_term_71);
    // q = 0
    result_by_lpq.set(c_Key3(8, 3, 0), common_term_72);
    result_by_q.set(c_Key1(0), common_term_72);
    // q = 1
    result_by_lpq.set(c_Key3(8, 3, 1), common_term_73);
    result_by_q.set(c_Key1(1), common_term_73);
    // q = 2
    result_by_lpq.set(c_Key3(8, 3, 2), common_term_74);
    result_by_q.set(c_Key1(2), common_term_74);
    // q = 3
    result_by_lpq.set(c_Key3(8, 3, 3), common_term_75);
    result_by_q.set(c_Key1(3), common_term_75);
    // q = 4
    result_by_lpq.set(c_Key3(8, 3, 4), common_term_76);
    result_by_q.set(c_Key1(4), common_term_76);
    // q = 5
    result_by_lpq.set(c_Key3(8, 3, 5), common_term_77);
    result_by_q.set(c_Key1(5), common_term_77);
    // q = 6
    result_by_lpq.set(c_Key3(8, 3, 6), common_term_78);
    result_by_q.set(c_Key1(6), common_term_78);
    // q = 7
    result_by_lpq.set(c_Key3(8, 3, 7), common_term_79);
    result_by_q.set(c_Key1(7), common_term_79);
    // q = 8
    result_by_lpq.set(c_Key3(8, 3, 8), common_term_80);
    result_by_q.set(c_Key1(8), common_term_80);
    // q = 9
    result_by_lpq.set(c_Key3(8, 3, 9), common_term_81);
    result_by_q.set(c_Key1(9), common_term_81);
    // q = 10
    result_by_lpq.set(c_Key3(8, 3, 10), common_term_82);
    result_by_q.set(c_Key1(10), common_term_82);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 3), result_by_q);
    result_by_q.clear();

    // l , p = (8, 4).
    // q = -10
    result_by_lpq.set(c_Key3(8, 4, -10), common_term_83);
    result_by_q.set(c_Key1(-10), common_term_83);
    // q = -9
    result_by_lpq.set(c_Key3(8, 4, -9), common_term_84);
    result_by_q.set(c_Key1(-9), common_term_84);
    // q = -8
    result_by_lpq.set(c_Key3(8, 4, -8), common_term_85);
    result_by_q.set(c_Key1(-8), common_term_85);
    // q = -7
    result_by_lpq.set(c_Key3(8, 4, -7), common_term_86);
    result_by_q.set(c_Key1(-7), common_term_86);
    // q = -6
    result_by_lpq.set(c_Key3(8, 4, -6), common_term_87);
    result_by_q.set(c_Key1(-6), common_term_87);
    // q = -5
    result_by_lpq.set(c_Key3(8, 4, -5), common_term_88);
    result_by_q.set(c_Key1(-5), common_term_88);
    // q = -4
    result_by_lpq.set(c_Key3(8, 4, -4), common_term_89);
    result_by_q.set(c_Key1(-4), common_term_89);
    // q = -3
    result_by_lpq.set(c_Key3(8, 4, -3), common_term_90);
    result_by_q.set(c_Key1(-3), common_term_90);
    // q = -2
    result_by_lpq.set(c_Key3(8, 4, -2), common_term_91);
    result_by_q.set(c_Key1(-2), common_term_91);
    // q = -1
    result_by_lpq.set(c_Key3(8, 4, -1), common_term_92);
    result_by_q.set(c_Key1(-1), common_term_92);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -7.5)*(2.1875*eccentricity_6 + 13.125*eccentricity_4 + 10.5*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(8, 4, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(8, 4, 1), common_term_92);
    result_by_q.set(c_Key1(1), common_term_92);
    // q = 2
    result_by_lpq.set(c_Key3(8, 4, 2), common_term_91);
    result_by_q.set(c_Key1(2), common_term_91);
    // q = 3
    result_by_lpq.set(c_Key3(8, 4, 3), common_term_90);
    result_by_q.set(c_Key1(3), common_term_90);
    // q = 4
    result_by_lpq.set(c_Key3(8, 4, 4), common_term_89);
    result_by_q.set(c_Key1(4), common_term_89);
    // q = 5
    result_by_lpq.set(c_Key3(8, 4, 5), common_term_88);
    result_by_q.set(c_Key1(5), common_term_88);
    // q = 6
    result_by_lpq.set(c_Key3(8, 4, 6), common_term_87);
    result_by_q.set(c_Key1(6), common_term_87);
    // q = 7
    result_by_lpq.set(c_Key3(8, 4, 7), common_term_86);
    result_by_q.set(c_Key1(7), common_term_86);
    // q = 8
    result_by_lpq.set(c_Key3(8, 4, 8), common_term_85);
    result_by_q.set(c_Key1(8), common_term_85);
    // q = 9
    result_by_lpq.set(c_Key3(8, 4, 9), common_term_84);
    result_by_q.set(c_Key1(9), common_term_84);
    // q = 10
    result_by_lpq.set(c_Key3(8, 4, 10), common_term_83);
    result_by_q.set(c_Key1(10), common_term_83);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 4), result_by_q);
    result_by_q.clear();

    // l , p = (8, 5).
    // q = -10
    result_by_lpq.set(c_Key3(8, 5, -10), common_term_82);
    result_by_q.set(c_Key1(-10), common_term_82);
    // q = -9
    result_by_lpq.set(c_Key3(8, 5, -9), common_term_81);
    result_by_q.set(c_Key1(-9), common_term_81);
    // q = -8
    result_by_lpq.set(c_Key3(8, 5, -8), common_term_80);
    result_by_q.set(c_Key1(-8), common_term_80);
    // q = -7
    result_by_lpq.set(c_Key3(8, 5, -7), common_term_79);
    result_by_q.set(c_Key1(-7), common_term_79);
    // q = -6
    result_by_lpq.set(c_Key3(8, 5, -6), common_term_78);
    result_by_q.set(c_Key1(-6), common_term_78);
    // q = -5
    result_by_lpq.set(c_Key3(8, 5, -5), common_term_77);
    result_by_q.set(c_Key1(-5), common_term_77);
    // q = -4
    result_by_lpq.set(c_Key3(8, 5, -4), common_term_76);
    result_by_q.set(c_Key1(-4), common_term_76);
    // q = -3
    result_by_lpq.set(c_Key3(8, 5, -3), common_term_75);
    result_by_q.set(c_Key1(-3), common_term_75);
    // q = -2
    result_by_lpq.set(c_Key3(8, 5, -2), common_term_74);
    result_by_q.set(c_Key1(-2), common_term_74);
    // q = -1
    result_by_lpq.set(c_Key3(8, 5, -1), common_term_73);
    result_by_q.set(c_Key1(-1), common_term_73);
    // q = 0
    result_by_lpq.set(c_Key3(8, 5, 0), common_term_72);
    result_by_q.set(c_Key1(0), common_term_72);
    // q = 1
    result_by_lpq.set(c_Key3(8, 5, 1), common_term_71);
    result_by_q.set(c_Key1(1), common_term_71);
    // q = 2
    result_by_lpq.set(c_Key3(8, 5, 2), common_term_70);
    result_by_q.set(c_Key1(2), common_term_70);
    // q = 3
    result_by_lpq.set(c_Key3(8, 5, 3), common_term_69);
    result_by_q.set(c_Key1(3), common_term_69);
    // q = 4
    result_by_lpq.set(c_Key3(8, 5, 4), common_term_68);
    result_by_q.set(c_Key1(4), common_term_68);
    // q = 5
    result_by_lpq.set(c_Key3(8, 5, 5), common_term_67);
    result_by_q.set(c_Key1(5), common_term_67);
    // q = 6
    result_by_lpq.set(c_Key3(8, 5, 6), common_term_66);
    result_by_q.set(c_Key1(6), common_term_66);
    // q = 7
    result_by_lpq.set(c_Key3(8, 5, 7), common_term_65);
    result_by_q.set(c_Key1(7), common_term_65);
    // q = 8
    result_by_lpq.set(c_Key3(8, 5, 8), common_term_64);
    result_by_q.set(c_Key1(8), common_term_64);
    // q = 9
    result_by_lpq.set(c_Key3(8, 5, 9), common_term_63);
    result_by_q.set(c_Key1(9), common_term_63);
    // q = 10
    result_by_lpq.set(c_Key3(8, 5, 10), common_term_62);
    result_by_q.set(c_Key1(10), common_term_62);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 5), result_by_q);
    result_by_q.clear();

    // l , p = (8, 6).
    // q = -10
    result_by_lpq.set(c_Key3(8, 6, -10), common_term_61);
    result_by_q.set(c_Key1(-10), common_term_61);
    // q = -9
    result_by_lpq.set(c_Key3(8, 6, -9), common_term_60);
    result_by_q.set(c_Key1(-9), common_term_60);
    // q = -8
    result_by_lpq.set(c_Key3(8, 6, -8), common_term_59);
    result_by_q.set(c_Key1(-8), common_term_59);
    // q = -7
    result_by_lpq.set(c_Key3(8, 6, -7), common_term_58);
    result_by_q.set(c_Key1(-7), common_term_58);
    // q = -6
    result_by_lpq.set(c_Key3(8, 6, -6), common_term_57);
    result_by_q.set(c_Key1(-6), common_term_57);
    // q = -5
    result_by_lpq.set(c_Key3(8, 6, -5), common_term_56);
    result_by_q.set(c_Key1(-5), common_term_56);
    // q = -4
    result_by_lpq.set(c_Key3(8, 6, -4), common_term_55);
    result_by_q.set(c_Key1(-4), common_term_55);
    // q = -3
    result_by_lpq.set(c_Key3(8, 6, -3), common_term_54);
    result_by_q.set(c_Key1(-3), common_term_54);
    // q = -2
    result_by_lpq.set(c_Key3(8, 6, -2), common_term_53);
    result_by_q.set(c_Key1(-2), common_term_53);
    // q = -1
    result_by_lpq.set(c_Key3(8, 6, -1), common_term_52);
    result_by_q.set(c_Key1(-1), common_term_52);
    // q = 0
    result_by_lpq.set(c_Key3(8, 6, 0), common_term_51);
    result_by_q.set(c_Key1(0), common_term_51);
    // q = 1
    result_by_lpq.set(c_Key3(8, 6, 1), common_term_50);
    result_by_q.set(c_Key1(1), common_term_50);
    // q = 2
    result_by_lpq.set(c_Key3(8, 6, 2), common_term_49);
    result_by_q.set(c_Key1(2), common_term_49);
    // q = 3
    result_by_lpq.set(c_Key3(8, 6, 3), common_term_48);
    result_by_q.set(c_Key1(3), common_term_48);
    // q = 4
    result_by_lpq.set(c_Key3(8, 6, 4), common_term_47);
    result_by_q.set(c_Key1(4), common_term_47);
    // q = 5
    result_by_lpq.set(c_Key3(8, 6, 5), common_term_46);
    result_by_q.set(c_Key1(5), common_term_46);
    // q = 6
    result_by_lpq.set(c_Key3(8, 6, 6), common_term_45);
    result_by_q.set(c_Key1(6), common_term_45);
    // q = 7
    result_by_lpq.set(c_Key3(8, 6, 7), common_term_44);
    result_by_q.set(c_Key1(7), common_term_44);
    // q = 8
    result_by_lpq.set(c_Key3(8, 6, 8), common_term_43);
    result_by_q.set(c_Key1(8), common_term_43);
    // q = 9
    result_by_lpq.set(c_Key3(8, 6, 9), common_term_42);
    result_by_q.set(c_Key1(9), common_term_42);
    // q = 10
    result_by_lpq.set(c_Key3(8, 6, 10), common_term_41);
    result_by_q.set(c_Key1(10), common_term_41);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 6), result_by_q);
    result_by_q.clear();

    // l , p = (8, 7).
    // q = -10
    result_by_lpq.set(c_Key3(8, 7, -10), common_term_40);
    result_by_q.set(c_Key1(-10), common_term_40);
    // q = -9
    result_by_lpq.set(c_Key3(8, 7, -9), common_term_39);
    result_by_q.set(c_Key1(-9), common_term_39);
    // q = -8
    result_by_lpq.set(c_Key3(8, 7, -8), common_term_38);
    result_by_q.set(c_Key1(-8), common_term_38);
    // q = -7
    result_by_lpq.set(c_Key3(8, 7, -7), common_term_37);
    result_by_q.set(c_Key1(-7), common_term_37);
    // q = -6
    result_by_lpq.set(c_Key3(8, 7, -6), common_term_36);
    result_by_q.set(c_Key1(-6), common_term_36);
    // q = -5
    result_by_lpq.set(c_Key3(8, 7, -5), common_term_35);
    result_by_q.set(c_Key1(-5), common_term_35);
    // q = -4
    result_by_lpq.set(c_Key3(8, 7, -4), common_term_34);
    result_by_q.set(c_Key1(-4), common_term_34);
    // q = -3
    result_by_lpq.set(c_Key3(8, 7, -3), common_term_33);
    result_by_q.set(c_Key1(-3), common_term_33);
    // q = -2
    result_by_lpq.set(c_Key3(8, 7, -2), common_term_32);
    result_by_q.set(c_Key1(-2), common_term_32);
    // q = -1
    result_by_lpq.set(c_Key3(8, 7, -1), common_term_31);
    result_by_q.set(c_Key1(-1), common_term_31);
    // q = 0
    result_by_lpq.set(c_Key3(8, 7, 0), common_term_30);
    result_by_q.set(c_Key1(0), common_term_30);
    // q = 1
    result_by_lpq.set(c_Key3(8, 7, 1), common_term_29);
    result_by_q.set(c_Key1(1), common_term_29);
    // q = 2
    result_by_lpq.set(c_Key3(8, 7, 2), common_term_28);
    result_by_q.set(c_Key1(2), common_term_28);
    // q = 3
    result_by_lpq.set(c_Key3(8, 7, 3), common_term_27);
    result_by_q.set(c_Key1(3), common_term_27);
    // q = 4
    result_by_lpq.set(c_Key3(8, 7, 4), common_term_26);
    result_by_q.set(c_Key1(4), common_term_26);
    // q = 5
    result_by_lpq.set(c_Key3(8, 7, 5), common_term_25);
    result_by_q.set(c_Key1(5), common_term_25);
    // q = 6
    result_by_lpq.set(c_Key3(8, 7, 6), common_term_24);
    result_by_q.set(c_Key1(6), common_term_24);
    // q = 7
    result_by_lpq.set(c_Key3(8, 7, 7), common_term_23);
    result_by_q.set(c_Key1(7), common_term_23);
    // q = 8
    result_by_lpq.set(c_Key3(8, 7, 8), common_term_22);
    result_by_q.set(c_Key1(8), common_term_22);
    // q = 9
    result_by_lpq.set(c_Key3(8, 7, 9), common_term_21);
    result_by_q.set(c_Key1(9), common_term_21);
    // q = 10
    result_by_lpq.set(c_Key3(8, 7, 10), common_term_20);
    result_by_q.set(c_Key1(10), common_term_20);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 7), result_by_q);
    result_by_q.clear();

    // l , p = (8, 8).
    // q = -10
    result_by_lpq.set(c_Key3(8, 8, -10), common_term_19);
    result_by_q.set(c_Key1(-10), common_term_19);
    // q = -9
    result_by_lpq.set(c_Key3(8, 8, -9), common_term_18);
    result_by_q.set(c_Key1(-9), common_term_18);
    // q = -8
    result_by_lpq.set(c_Key3(8, 8, -8), common_term_17);
    result_by_q.set(c_Key1(-8), common_term_17);
    // q = -7
    result_by_lpq.set(c_Key3(8, 8, -7), common_term_16);
    result_by_q.set(c_Key1(-7), common_term_16);
    // q = -6
    result_by_lpq.set(c_Key3(8, 8, -6), common_term_15);
    result_by_q.set(c_Key1(-6), common_term_15);
    // q = -5
    result_by_lpq.set(c_Key3(8, 8, -5), common_term_14);
    result_by_q.set(c_Key1(-5), common_term_14);
    // q = -4
    result_by_lpq.set(c_Key3(8, 8, -4), common_term_13);
    result_by_q.set(c_Key1(-4), common_term_13);
    // q = -3
    result_by_lpq.set(c_Key3(8, 8, -3), common_term_12);
    result_by_q.set(c_Key1(-3), common_term_12);
    // q = -2
    result_by_lpq.set(c_Key3(8, 8, -2), common_term_11);
    result_by_q.set(c_Key1(-2), common_term_11);
    // q = -1
    result_by_lpq.set(c_Key3(8, 8, -1), common_term_10);
    result_by_q.set(c_Key1(-1), common_term_10);
    // q = 0
    result_by_lpq.set(c_Key3(8, 8, 0), common_term_9);
    result_by_q.set(c_Key1(0), common_term_9);
    // q = 1
    result_by_lpq.set(c_Key3(8, 8, 1), common_term_8);
    result_by_q.set(c_Key1(1), common_term_8);
    // q = 2
    result_by_lpq.set(c_Key3(8, 8, 2), common_term_7);
    result_by_q.set(c_Key1(2), common_term_7);
    // q = 3
    result_by_lpq.set(c_Key3(8, 8, 3), common_term_6);
    result_by_q.set(c_Key1(3), common_term_6);
    // q = 4
    result_by_lpq.set(c_Key3(8, 8, 4), common_term_5);
    result_by_q.set(c_Key1(4), common_term_5);
    // q = 5
    result_by_lpq.set(c_Key3(8, 8, 5), common_term_4);
    result_by_q.set(c_Key1(5), common_term_4);
    // q = 6
    result_by_lpq.set(c_Key3(8, 8, 6), common_term_3);
    result_by_q.set(c_Key1(6), common_term_3);
    // q = 7
    result_by_lpq.set(c_Key3(8, 8, 7), common_term_2);
    result_by_q.set(c_Key1(7), common_term_2);
    // q = 9
    result_by_lpq.set(c_Key3(8, 8, 9), common_term_1);
    result_by_q.set(c_Key1(9), common_term_1);
    // q = 10
    result_by_lpq.set(c_Key3(8, 8, 10), common_term_0);
    result_by_q.set(c_Key1(10), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 8), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l8_e15(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 8.
    // Functions are taylor expanded and truncated at eccentricity^15.

    // Functions Calculated for l = 8.

    c_IntMap<c_Key3, double> result_by_lpq(277);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(9);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double eccentricity_6 = eccentricity_3 * eccentricity_3;
    double eccentricity_7 = eccentricity_6 * eccentricity;
    double eccentricity_14 = eccentricity_7 * eccentricity_7;
    double eccentricity_15 = eccentricity_14 * eccentricity;
    double eccentricity_12 = eccentricity_6 * eccentricity_6;
    double eccentricity_13 = eccentricity_12 * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_5 = eccentricity_4 * eccentricity;
    double eccentricity_10 = eccentricity_5 * eccentricity_5;
    double eccentricity_11 = eccentricity_10 * eccentricity;
    double eccentricity_8 = eccentricity_4 * eccentricity_4;
    double eccentricity_9 = eccentricity_8 * eccentricity;
    double common_term_0 = 0.00011079522764105174*eccentricity_15;
    double common_term_1 = 5.4864220600827744e-5*eccentricity_14;
    double common_term_2 = 7.1362201676549913e-5*eccentricity_15 + 2.3929840083154462e-5*eccentricity_13;
    double common_term_3 = 2.4995580551136107e-5*eccentricity_14 + 8.5511196622307733e-6*eccentricity_12;
    double common_term_4 = 1.074487572306996e-5*eccentricity_15 + 6.0945362239689022e-6*eccentricity_13 + 2.1669462129667208e-6*eccentricity_11;
    double common_term_5 = 1.21398341884453e-6*eccentricity_14 + 7.2651114317780984e-7*eccentricity_12 + 2.7557319223985891e-7*eccentricity_10;
    double common_term_6 = 2.5365045247212074e-8*eccentricity_15 + 1.9744743382642146e-8*eccentricity_13 + 1.2782936163470018e-8*eccentricity_11 + 5.3822889109347443e-9*eccentricity_9;
    double common_term_7 = -2.7806556216540593e-6*eccentricity_15 - 2.7637997492140423e-6*eccentricity_13 - 2.627229774650022e-6*eccentricity_11 - 2.2767082093253968e-6*eccentricity_9 - 1.5500992063492063e-6*eccentricity_7;
    double common_term_8 = 0.0007178452013521458*eccentricity_14 + 0.00080421443268665491*eccentricity_12 + 0.00090525793650793651*eccentricity_10 + 0.00099206349206349206*eccentricity_8 + 0.0013888888888888889*eccentricity_6;
    double common_term_9 = 0.0025082426241465977*eccentricity_15 + 0.0025113752910069057*eccentricity_13 + 0.0027588435581752232*eccentricity_11 - 0.0015537806919642857*eccentricity_9 + 0.02373046875*eccentricity_7 - 0.06328125*eccentricity_5;
    double common_term_10 = -0.0020223398001175779*eccentricity_14 + 0.016881613756613757*eccentricity_12 - 0.14656084656084656*eccentricity_10 + 0.58888888888888889*eccentricity_8 - 1.3333333333333333*eccentricity_6 + 0.66666666666666667*eccentricity_4;
    double common_term_11 = -0.088458628195326818*eccentricity_15 + 0.4550276136902905*eccentricity_13 - 2.6191353167175616*eccentricity_11 + 8.5598980938946759*eccentricity_9 - 16.215006510416667*eccentricity_7 + 11.881510416666667*eccentricity_5 - 2.6041666666666667*eccentricity_3;
    double common_term_12 = 7.017890625*eccentricity_14 - 30.029263392857143*eccentricity_12 + 83.11640625*eccentricity_10 - 139.6125*eccentricity_8 + 115.3125*eccentricity_6 - 40.5*eccentricity_4 + 4.5*eccentricity_2;
    double common_term_13 = 72.291540367743115*eccentricity_15 - 253.89152100763203*eccentricity_13 + 610.67264824196144*eccentricity_11 - 935.56048651801215*eccentricity_9 + 800.71652560763889*eccentricity_7 - 343.30989583333333*eccentricity_5 + 63.4375*eccentricity_3 - 3.5*eccentricity;
    double common_term_14 = -1710.8899369724742*eccentricity_14 + 3655.4365760030864*eccentricity_12 - 5182.6406944444444*eccentricity_10 + 4439.1657986111111*eccentricity_8 - 2099.8611111111111*eccentricity_6 + 490.75*eccentricity_4 - 46.0*eccentricity_2 + 1.0;
    double common_term_15 = -9678.2896001151265*eccentricity_15 + 18682.353862372807*eccentricity_13 - 24732.193187713623*eccentricity_11 + 20854.591314697266*eccentricity_9 - 10362.46142578125*eccentricity_7 + 2757.0234375*eccentricity_5 - 333.5625*eccentricity_3 + 12.5*eccentricity;
    double common_term_16 = 84165.626218550209*eccentricity_14 - 104690.64519262566*eccentricity_12 + 86174.851707175926*eccentricity_10 - 43779.236111111111*eccentricity_8 + 12571.145833333333*eccentricity_6 - 1764.1666666666667*eccentricity_4 + 86.5*eccentricity_2;
    double common_term_17 = 341959.49825736394*eccentricity_15 - 401594.35779903985*eccentricity_13 + 321267.57119522902*eccentricity_11 - 164205.81790861907*eccentricity_9 + 49339.254264322917*eccentricity_7 - 7600.3268229166667*eccentricity_5 + 437.72916666666667*eccentricity_3;
    double common_term_18 = -1418742.2822042411*eccentricity_14 + 1100497.1265904018*eccentricity_12 - 560229.34285714286*eccentricity_10 + 172758.65625*eccentricity_8 - 28273.05*eccentricity_6 + 1808.25*eccentricity_4;
    double common_term_19 = -4673634.5473623291*eccentricity_15 + 3511333.9634982762*eccentricity_13 - 1768754.423578666*eccentricity_11 + 552633.89492439391*eccentricity_9 - 94090.807693142361*eccentricity_7 + 6461.7950520833333*eccentricity_5;
    double common_term_20 = 10545772.889832256*eccentricity_14 - 5233760.1698913323*eccentricity_12 + 1642304.3425347222*eccentricity_10 - 286780.58888888889*eccentricity_8 + 20690.784722222222*eccentricity_6;
    double common_term_21 = 30060994.488160975*eccentricity_15 - 14656166.780386141*eccentricity_13 + 4590201.1259119851*eccentricity_11 - 813848.52869306292*eccentricity_9 + 60764.489606584821*eccentricity_7;
    double common_term_22 = -39137278.192853393*eccentricity_14 + 12179584.482504685*eccentricity_12 - 2176625.8677358907*eccentricity_10 + 166371.02557043651*eccentricity_8;
    double common_term_23 = -100269202.88532857*eccentricity_15 + 30904655.242274058*eccentricity_13 - 5536752.3833186226*eccentricity_11 + 429802.86250706695*eccentricity_9;
    double common_term_24 = 75427724.874900061*eccentricity_14 - 13491473.585716315*eccentricity_12 + 1057233.6862946429*eccentricity_10;
    double common_term_25 = 177910875.05159583*eccentricity_15 - 31671129.351751119*eccentricity_13 + 2493744.2467064668*eccentricity_11;
    double common_term_26 = -71955647.267619634*eccentricity_14 + 5672187.3666974703*eccentricity_12;
    double common_term_27 = -158819986.31638283*eccentricity_15 + 12497962.364256332*eccentricity_13;
    double common_term_28 = 26775478.224996169*eccentricity_14;
    double common_term_29 = 55948914.794879542*eccentricity_15;
    double common_term_30 = 1.7361061035642716*eccentricity_15;
    double common_term_31 = 1.2762573453994817*eccentricity_14;
    double common_term_32 = 5.3905542732752184*eccentricity_15 + 0.93835479278592447*eccentricity_13;
    double common_term_33 = 4.1378410456730769*eccentricity_14 + 0.69002422382305195*eccentricity_12;
    double common_term_34 = 11.73128761832799*eccentricity_15 + 3.1708903963106726*eccentricity_13 + 0.50749178561182651*eccentricity_11;
    double common_term_35 = 9.2376283748321509*eccentricity_14 + 2.4261317884950697*eccentricity_12 + 0.37329840443121693*eccentricity_10;
    double common_term_36 = 21.530950284507368*eccentricity_15 + 7.2599628672661719*eccentricity_13 + 1.8536271122523717*eccentricity_11 + 0.27462201799665179*eccentricity_9;
    double common_term_37 = 17.26639119368286*eccentricity_14 + 5.6951335668816138*eccentricity_12 + 1.4143146494708995*eccentricity_10 + 0.20204613095238095*eccentricity_8;
    double common_term_38 = 35.630606926419767*eccentricity_15 + 13.819591798849627*eccentricity_13 + 4.4596811829420625*eccentricity_11 + 1.0777563185918899*eccentricity_9 + 0.14865606398809524*eccentricity_7;
    double common_term_39 = 0.109375*eccentricity_6*std::pow(1.0 - eccentricity_2, -7.5);
    double common_term_40 = 54.939037858091544*eccentricity_15 + 23.526691280302758*eccentricity_13 + 8.803474856179858*eccentricity_11 + 2.7210094633556548*eccentricity_9 + 0.62366536458333333*eccentricity_7 + 0.08046875*eccentricity_5;
    double common_term_41 = 45.201081590332892*eccentricity_14 + 19.065321955605159*eccentricity_12 + 7.0086185515873016*eccentricity_10 + 2.1216145833333333*eccentricity_8 + 0.475*eccentricity_6 + 0.0625*eccentricity_4;
    double common_term_42 = 80.436430373958179*eccentricity_15 + 37.125392362730844*eccentricity_13 + 15.425517817905971*eccentricity_11 + 5.57567138671875*eccentricity_9 + 1.64326171875*eccentricity_7 + 0.43359375*eccentricity_5 - 0.0625*eccentricity_3;
    double common_term_43 = 66.824699332423942*eccentricity_14 + 30.437646639384921*eccentricity_12 + 12.522894965277778*eccentricity_10 + 4.0916666666666667*eccentricity_8 + 2.5*eccentricity_6 - 1.625*eccentricity_4 + 0.75*eccentricity_2;
    double common_term_44 = 113.19083135138609*eccentricity_15 + 55.269214627604005*eccentricity_13 + 25.94196531507704*eccentricity_11 + 5.5236307779947917*eccentricity_9 + 16.171061197916667*eccentricity_7 - 17.3828125*eccentricity_5 + 10.3125*eccentricity_3 - 1.5*eccentricity;
    double common_term_45 = 92.612688137755102*eccentricity_14 + 57.024970703125*eccentricity_12 - 19.72640625*eccentricity_10 + 101.90625*eccentricity_8 - 122.3125*eccentricity_6 + 78.1875*eccentricity_4 - 18.0*eccentricity_2 + 1.0;
    double common_term_46 = 133.14591825483022*eccentricity_15 + 167.7827899007444*eccentricity_13 - 229.92991219414605*eccentricity_11 + 560.69659627278646*eccentricity_9 - 662.57438151041667*eccentricity_7 + 430.1171875*eccentricity_5 - 119.4375*eccentricity_3 + 10.5*eccentricity;
    double common_term_47 = 687.22760082103588*eccentricity_14 - 1391.7037338789683*eccentricity_12 + 2661.0693576388889*eccentricity_10 - 2985.2182291666667*eccentricity_8 + 1915.953125*eccentricity_6 - 580.375*eccentricity_4 + 62.25*eccentricity_2;
    double common_term_48 = 3140.3596028706857*eccentricity_15 - 6674.072648970059*eccentricity_13 + 11096.136155809675*eccentricity_11 - 11708.047570800781*eccentricity_9 + 7329.28740234375*eccentricity_7 - 2307.59765625*eccentricity_5 + 274.1875*eccentricity_3;
    double common_term_49 = -27709.08882154569*eccentricity_14 + 41527.90535094246*eccentricity_12 - 41181.945374503968*eccentricity_10 + 24965.729166666667*eccentricity_8 - 7957.5625*eccentricity_6 + 998.75*eccentricity_4;
    double common_term_50 = -103358.24492683926*eccentricity_15 + 142048.9446050667*eccentricity_13 - 132645.21026563493*eccentricity_11 + 77582.160425385975*eccentricity_9 - 24651.222298177083*eccentricity_7 + 3181.34921875*eccentricity_5;
    double common_term_51 = 450549.36138532366*eccentricity_14 - 397243.92134486607*eccentricity_12 + 223753.60797991071*eccentricity_10 - 70208.108035714286*eccentricity_8 + 9163.8375*eccentricity_6;
    double common_term_52 = 1340315.2745059372*eccentricity_15 - 1118951.471470651*eccentricity_13 + 606550.27029326706*eccentricity_11 - 186830.13835812523*eccentricity_9 + 24400.922632998512*eccentricity_7;
    double common_term_53 = -2991179.6441564297*eccentricity_14 + 1560397.2811733218*eccentricity_12 - 470044.07146990741*eccentricity_10 + 60988.455208333333*eccentricity_8;
    double common_term_54 = -7642510.1706927773*eccentricity_15 + 3838292.1718349859*eccentricity_13 - 1128041.6384100233*eccentricity_11 + 144690.15621730259*eccentricity_9;
    double common_term_55 = 9081820.6651140959*eccentricity_14 - 2600155.3861269127*eccentricity_12 + 328564.26045841601*eccentricity_10;
    double common_term_56 = 20770338.744247634*eccentricity_15 - 5788030.6720844356*eccentricity_13 + 718795.38038070565*eccentricity_11;
    double common_term_57 = -12497715.918680811*eccentricity_14 + 1522709.4454134537*eccentricity_12;
    double common_term_58 = -26270225.293139121*eccentricity_15 + 3136496.2068759774*eccentricity_13;
    double common_term_59 = 6303051.860868987*eccentricity_14;
    double common_term_60 = 12392155.598679251*eccentricity_15;
    double common_term_61 = 121.14309121499761*eccentricity_15;
    double common_term_62 = 85.273143039149751*eccentricity_14;
    double common_term_63 = 286.76737477258496*eccentricity_15 + 59.908468079385938*eccentricity_13;
    double common_term_64 = 213.74585329958268*eccentricity_14 + 41.99941369350332*eccentricity_12;
    double common_term_65 = 533.78907190873898*eccentricity_15 + 158.5322733228732*eccentricity_13 + 29.375192652658448*eccentricity_11;
    double common_term_66 = 408.8027098341112*eccentricity_14 + 117.0156742086039*eccentricity_12 + 20.492165178571429*eccentricity_10;
    double common_term_67 = 872.73445321116142*eccentricity_15 + 311.64210511511558*eccentricity_13 + 85.960494009485531*eccentricity_11 + 14.253836199509824*eccentricity_9;
    double common_term_68 = 681.15639702873644*eccentricity_14 + 236.46823033785273*eccentricity_12 + 62.845375881834215*eccentricity_10 + 9.8822296626984127*eccentricity_8;
    double common_term_69 = 1317.1940303954908*eccentricity_15 + 529.32218106269836*eccentricity_13 + 178.57800832475935*eccentricity_11 + 45.722182355608259*eccentricity_9 + 6.8260463169642857*eccentricity_7;
    double common_term_70 = 1042.9932195675338*eccentricity_14 + 409.50493126745297*eccentricity_12 + 134.20566716269841*eccentricity_10 + 33.096825396825397*eccentricity_8 + 4.6951388888888889*eccentricity_6;
    double common_term_71 = 1881.1768521114715*eccentricity_15 + 822.4778344982025*eccentricity_13 + 315.36323912080634*eccentricity_11 + 100.35344373914931*eccentricity_9 + 23.830978732638889*eccentricity_7 + 3.2138020833333333*eccentricity_5;
    double common_term_72 = std::pow(1.0 - eccentricity_2, -7.5)*(0.65625*eccentricity_6 + 2.1875*eccentricity_4);
    double common_term_73 = 2579.1932129021996*eccentricity_15 + 1202.4828373346127*eccentricity_13 + 504.9256145134174*eccentricity_11 + 184.36379711009838*eccentricity_9 + 55.222623697916667*eccentricity_7 + 12.141927083333333*eccentricity_5 + 1.4791666666666667*eccentricity_3;
    double common_term_74 = 2086.2112137965719*eccentricity_14 + 955.65467282572751*eccentricity_12 + 392.9662037037037*eccentricity_10 + 139.89722222222222*eccentricity_8 + 40.614583333333333*eccentricity_6 + 8.5833333333333333*eccentricity_4 + eccentricity_2;
    double common_term_75 = 3426.2534545571828*eccentricity_15 + 1681.1850924055917*eccentricity_13 + 756.32119033813477*eccentricity_11 + 304.38755493164062*eccentricity_9 + 105.59130859375*eccentricity_7 + 29.6484375*eccentricity_5 + 6.1875*eccentricity_3 + 0.5*eccentricity;
    double common_term_76 = 2794.3606267026487*eccentricity_14 + 1349.584135320216*eccentricity_12 + 595.95076388888889*eccentricity_10 + 234.75564236111111*eccentricity_8 + 78.555555555555556*eccentricity_6 + 23.5*eccentricity_4 + 2.0*eccentricity_2 + 1.0;
    double common_term_77 = 4437.8718071490992*eccentricity_15 + 2270.9555916096588*eccentricity_13 + 1078.7587194089536*eccentricity_11 + 469.22595553927951*eccentricity_9 + 173.39816623263889*eccentricity_7 + 73.815104166666667*eccentricity_5 - 0.3125*eccentricity_3 + 8.5*eccentricity;
    double common_term_78 = 3646.0642445591518*eccentricity_14 + 1835.5665234375*eccentricity_12 + 873.31640625*eccentricity_10 + 321.95625*eccentricity_8 + 218.90625*eccentricity_6 - 37.25*eccentricity_4 + 42.0*eccentricity_2;
    double common_term_79 = 5634.8871911270048*eccentricity_15 + 2960.8069926811148*eccentricity_13 + 1569.906435168231*eccentricity_11 + 458.83704653139468*eccentricity_9 + 648.71975911458333*eccentricity_7 - 226.96223958333333*eccentricity_5 + 157.64583333333333*eccentricity_3;
    double common_term_80 = 4514.8397351799034*eccentricity_14 + 2868.7049109726356*eccentricity_12 + 232.20453042328042*eccentricity_10 + 1927.3425347222222*eccentricity_8 - 932.79583333333333*eccentricity_6 + 497.97916666666667*eccentricity_4;
    double common_term_81 = 6332.2276472097976*eccentricity_15 + 5690.9572262048721*eccentricity_13 - 1628.573877116612*eccentricity_11 + 5629.1753435407366*eccentricity_9 - 3134.50556640625*eccentricity_7 + 1395.14296875*eccentricity_5;
    double common_term_82 = 12851.250283826035*eccentricity_14 - 9021.6337471524104*eccentricity_12 + 15875.237475198413*eccentricity_10 - 9247.1106150793651*eccentricity_8 + 3575.9618055555556*eccentricity_6;
    double common_term_83 = 32633.328711873985*eccentricity_15 - 32777.046731996615*eccentricity_13 + 42844.59543250015*eccentricity_11 - 24842.062921966068*eccentricity_9 + 8555.9166837177579*eccentricity_7;
    double common_term_84 = -100913.54372666396*eccentricity_14 + 110448.49248883929*eccentricity_12 - 62129.775892857143*eccentricity_10 + 19375.531473214286*eccentricity_8;
    double common_term_85 = -281991.78143195508*eccentricity_15 + 272490.68950632819*eccentricity_13 - 146790.81474965427*eccentricity_11 + 41945.743762018651*eccentricity_9;
    double common_term_86 = 645651.00100172004*eccentricity_14 - 331033.60035002806*eccentricity_12 + 87458.838555169753*eccentricity_10;
    double common_term_87 = 1475032.0988348779*eccentricity_15 - 718011.57625047417*eccentricity_13 + 176637.60430073428*eccentricity_11;
    double common_term_88 = -1506624.8069187101*eccentricity_14 + 347116.40083112144*eccentricity_12;
    double common_term_89 = -3072374.4261272348*eccentricity_15 + 666099.66287382847*eccentricity_13;
    double common_term_90 = 1251825.6868847893*eccentricity_14;
    double common_term_91 = 2309591.5418258783*eccentricity_15;
    double common_term_92 = 2926.6658339357941*eccentricity_15;
    double common_term_93 = 1937.4461366158953*eccentricity_14;
    double common_term_94 = 4179.9295842065951*eccentricity_15 + 1274.9490624733074*eccentricity_13;
    double common_term_95 = 3019.4252945214588*eccentricity_14 + 833.35244924025306*eccentricity_12;
    double common_term_96 = 6196.3667029721152*eccentricity_15 + 2151.49817300716*eccentricity_13 + 540.54460558705516*eccentricity_11;
    double common_term_97 = 4552.5318770684952*eccentricity_14 + 1512.621510137962*eccentricity_12 + 347.53703627921076*eccentricity_10;
    double common_term_98 = 8413.136035539822*eccentricity_15 + 3309.5092659818069*eccentricity_13 + 1049.0975573042881*eccentricity_11 + 221.16175382637683*eccentricity_9;
    double common_term_99 = 6271.2194009993912*eccentricity_14 + 2378.333583984375*eccentricity_12 + 717.35178571428571*eccentricity_10 + 139.04497767857143*eccentricity_8;
    double common_term_100 = 10856.895523331472*eccentricity_15 + 4625.5713758638226*eccentricity_13 + 1687.7395159773515*eccentricity_11 + 483.05924309624566*eccentricity_9 + 86.155606708829365*eccentricity_7;
    double common_term_101 = 8173.0600789632018*eccentricity_14 + 3372.4788363807136*eccentricity_12 + 1181.0311042906746*eccentricity_10 + 319.79875992063492*eccentricity_8 + 52.440972222222222*eccentricity_6;
    double common_term_102 = 13501.639044812935*eccentricity_15 + 6087.4061027492796*eccentricity_13 + 2427.3691937582833*eccentricity_11 + 813.49096854073661*eccentricity_9 + 207.61494140625*eccentricity_7 + 31.21171875*eccentricity_5;
    double common_term_103 = 10235.202401781121*eccentricity_14 + 4480.4448167266038*eccentricity_12 + 1721.8277571097884*eccentricity_10 + 550.18836805555556*eccentricity_8 + 131.67916666666667*eccentricity_6 + 18.041666666666667*eccentricity_4;
    double common_term_104 = 16320.07322272787*eccentricity_15 + 7674.3587744977487*eccentricity_13 + 3253.691183524157*eccentricity_11 + 1200.9678082501447*eccentricity_9 + 364.11526692708333*eccentricity_7 + 81.131510416666667*eccentricity_5 + 10.020833333333333*eccentricity_3;
    double common_term_105 = std::pow(1.0 - eccentricity_2, -7.5)*(1.640625*eccentricity_6 + 8.75*eccentricity_4 + 5.25*eccentricity_2);
    double common_term_106 = 19282.115231111371*eccentricity_15 + 9363.2274901158637*eccentricity_13 + 4149.5350677772805*eccentricity_11 + 1633.5173536512587*eccentricity_9 + 547.90695529513889*eccentricity_7 + 146.03385416666667*eccentricity_5 + 27.0625*eccentricity_3 + 2.5*eccentricity;
    double common_term_107 = 14736.952994648959*eccentricity_14 + 6960.0260164689429*eccentricity_12 + 2979.5158159722222*eccentricity_10 + 1121.8220486111111*eccentricity_8 + 354.40972222222222*eccentricity_6 + 86.6875*eccentricity_4 + 14.0*eccentricity_2 + 1.0;
    double common_term_108 = 22354.895010375599*eccentricity_15 + 11128.178120109013*eccentricity_13 + 5095.2647703552246*eccentricity_11 + 2096.8638977050781*eccentricity_9 + 749.23486328125*eccentricity_7 + 219.8671875*eccentricity_5 + 47.8125*eccentricity_3 + 6.5*eccentricity;
    double common_term_109 = 17118.068152247299*eccentricity_14 + 8286.3820745494378*eccentricity_12 + 3662.7114945023148*eccentricity_10 + 1439.2086805555556*eccentricity_8 + 482.55208333333333*eccentricity_6 + 127.20833333333333*eccentricity_4 + 25.75*eccentricity_2;
    double common_term_110 = 25502.760370540711*eccentricity_15 + 12940.742162879182*eccentricity_13 + 6068.8542543643366*eccentricity_11 + 2573.966889558015*eccentricity_9 + 958.58414713541667*eccentricity_7 + 287.45442708333333*eccentricity_5 + 80.104166666666667*eccentricity_3;
    double common_term_111 = 19542.408988560268*eccentricity_14 + 9636.7225655691964*eccentricity_12 + 4352.0655133928571*eccentricity_10 + 1772.05546875*eccentricity_8 + 573.0*eccentricity_6 + 214.9375*eccentricity_4;
    double common_term_112 = 28686.436395698423*eccentricity_15 + 14774.660635731053*eccentricity_13 + 7024.8494779092294*eccentricity_11 + 3111.831438530816*eccentricity_9 + 1021.5815646701389*eccentricity_7 + 521.30130208333333*eccentricity_5;
    double common_term_113 = 22000.731539307358*eccentricity_14 + 10884.473990081662*eccentricity_12 + 5286.5708922371032*eccentricity_10 + 1618.5184771825397*eccentricity_8 + 1174.3086805555556*eccentricity_6;
    double common_term_114 = 31989.224256097365*eccentricity_15 + 16204.755631538119*eccentricity_13 + 8851.5450169699533*eccentricity_11 + 2188.9004241943359*eccentricity_9 + 2499.8066545758929*eccentricity_7;
    double common_term_115 = 23058.242806039166*eccentricity_14 + 14894.877307925485*eccentricity_12 + 2167.5672054122575*eccentricity_10 + 5088.2304067460317*eccentricity_8;
    double common_term_116 = 30834.009431963243*eccentricity_15 + 25660.809450528702*eccentricity_13 + 138.15669965323527*eccentricity_11 + 9986.2402226469925*eccentricity_9;
    double common_term_117 = 45857.294896129261*eccentricity_14 - 7051.1112312297078*eccentricity_12 + 19015.278426339286*eccentricity_10;
    double common_term_118 = 85307.598335647122*eccentricity_15 - 25928.956411922032*eccentricity_13 + 35295.352656819605*eccentricity_11;
    double common_term_119 = -69371.145960718407*eccentricity_14 + 64097.985970790077*eccentricity_12;
    double common_term_120 = -161863.38631548219*eccentricity_15 + 114223.29016644697*eccentricity_13;
    double common_term_121 = 200206.09740457872*eccentricity_14;
    double common_term_122 = 345825.64133652176*eccentricity_15;
    double common_term_123 = 38830.271627030094*eccentricity_15;
    double common_term_124 = 24045.531161456816*eccentricity_14;
    double common_term_125 = 17864.457369950933*eccentricity_15 + 14736.269914774821*eccentricity_13;
    double common_term_126 = 14014.760135060252*eccentricity_14 + 8925.3544891436688*eccentricity_12;
    double common_term_127 = 26932.785682546405*eccentricity_15 + 10366.636800894334*eccentricity_13 + 5333.4715783509754*eccentricity_11;
    double common_term_128 = 18970.483798292636*eccentricity_14 + 7331.8707470914502*eccentricity_12 + 3137.8065166170635*eccentricity_10;
    double common_term_129 = 29049.230429099558*eccentricity_15 + 13195.548930184624*eccentricity_13 + 4992.8608943285261*eccentricity_11 + 1812.6737448556083*eccentricity_9;
    double common_term_130 = 20765.734099138709*eccentricity_14 + 9026.688396577381*eccentricity_12 + 3283.9923115079365*eccentricity_10 + 1024.7152901785714*eccentricity_8;
    double common_term_131 = 31215.17029686814*eccentricity_15 + 14579.106566331298*eccentricity_13 + 6049.8636864980062*eccentricity_11 + 2087.4189422607422*eccentricity_9 + 564.30244140625*eccentricity_7;
    double common_term_132 = 22327.635171595982*eccentricity_14 + 10030.963895089286*eccentricity_12 + 3957.3959263392857*eccentricity_10 + 1280.1776785714286*eccentricity_8 + 300.8625*eccentricity_6;
    double common_term_133 = 32843.366350184594*eccentricity_15 + 15674.387087534975*eccentricity_13 + 6744.2472852979388*eccentricity_11 + 2515.3772408621652*eccentricity_9 + 754.56982421875*eccentricity_7 + 153.94921875*eccentricity_5;
    double common_term_134 = 23451.645246155754*eccentricity_14 + 10769.689983258929*eccentricity_12 + 4414.2023809523809*eccentricity_10 + 1544.878125*eccentricity_8 + 424.575*eccentricity_6 + 74.625*eccentricity_4;
    double common_term_135 = 33936.454313233495*eccentricity_15 + 16411.23981387479*eccentricity_13 + 7216.2684065682547*eccentricity_11 + 2798.1171752929687*eccentricity_9 + 909.78837890625*eccentricity_7 + 225.52734375*eccentricity_5 + 33.5625*eccentricity_3;
    double common_term_136 = 24130.470153924851*eccentricity_14 + 11217.020379464286*eccentricity_12 + 4692.5872395833333*eccentricity_10 + 1705.44375*eccentricity_8 + 507.9375*eccentricity_6 + 111.0*eccentricity_4 + 13.5*eccentricity_2;
    double common_term_137 = 34485.263931086781*eccentricity_15 + 16781.810886405158*eccentricity_13 + 7454.2881489562988*eccentricity_11 + 2941.4128967285156*eccentricity_9 + 988.83935546875*eccentricity_7 + 264.0234375*eccentricity_5 + 48.9375*eccentricity_3 + 4.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(31);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (8, 0).
    // q = -15
    result_by_lpq.set(c_Key3(8, 0, -15), common_term_0);
    result_by_q.set(c_Key1(-15), common_term_0);
    // q = -14
    result_by_lpq.set(c_Key3(8, 0, -14), common_term_1);
    result_by_q.set(c_Key1(-14), common_term_1);
    // q = -13
    result_by_lpq.set(c_Key3(8, 0, -13), common_term_2);
    result_by_q.set(c_Key1(-13), common_term_2);
    // q = -12
    result_by_lpq.set(c_Key3(8, 0, -12), common_term_3);
    result_by_q.set(c_Key1(-12), common_term_3);
    // q = -11
    result_by_lpq.set(c_Key3(8, 0, -11), common_term_4);
    result_by_q.set(c_Key1(-11), common_term_4);
    // q = -10
    result_by_lpq.set(c_Key3(8, 0, -10), common_term_5);
    result_by_q.set(c_Key1(-10), common_term_5);
    // q = -9
    result_by_lpq.set(c_Key3(8, 0, -9), common_term_6);
    result_by_q.set(c_Key1(-9), common_term_6);
    // q = -7
    result_by_lpq.set(c_Key3(8, 0, -7), common_term_7);
    result_by_q.set(c_Key1(-7), common_term_7);
    // q = -6
    result_by_lpq.set(c_Key3(8, 0, -6), common_term_8);
    result_by_q.set(c_Key1(-6), common_term_8);
    // q = -5
    result_by_lpq.set(c_Key3(8, 0, -5), common_term_9);
    result_by_q.set(c_Key1(-5), common_term_9);
    // q = -4
    result_by_lpq.set(c_Key3(8, 0, -4), common_term_10);
    result_by_q.set(c_Key1(-4), common_term_10);
    // q = -3
    result_by_lpq.set(c_Key3(8, 0, -3), common_term_11);
    result_by_q.set(c_Key1(-3), common_term_11);
    // q = -2
    result_by_lpq.set(c_Key3(8, 0, -2), common_term_12);
    result_by_q.set(c_Key1(-2), common_term_12);
    // q = -1
    result_by_lpq.set(c_Key3(8, 0, -1), common_term_13);
    result_by_q.set(c_Key1(-1), common_term_13);
    // q = 0
    result_by_lpq.set(c_Key3(8, 0, 0), common_term_14);
    result_by_q.set(c_Key1(0), common_term_14);
    // q = 1
    result_by_lpq.set(c_Key3(8, 0, 1), common_term_15);
    result_by_q.set(c_Key1(1), common_term_15);
    // q = 2
    result_by_lpq.set(c_Key3(8, 0, 2), common_term_16);
    result_by_q.set(c_Key1(2), common_term_16);
    // q = 3
    result_by_lpq.set(c_Key3(8, 0, 3), common_term_17);
    result_by_q.set(c_Key1(3), common_term_17);
    // q = 4
    result_by_lpq.set(c_Key3(8, 0, 4), common_term_18);
    result_by_q.set(c_Key1(4), common_term_18);
    // q = 5
    result_by_lpq.set(c_Key3(8, 0, 5), common_term_19);
    result_by_q.set(c_Key1(5), common_term_19);
    // q = 6
    result_by_lpq.set(c_Key3(8, 0, 6), common_term_20);
    result_by_q.set(c_Key1(6), common_term_20);
    // q = 7
    result_by_lpq.set(c_Key3(8, 0, 7), common_term_21);
    result_by_q.set(c_Key1(7), common_term_21);
    // q = 8
    result_by_lpq.set(c_Key3(8, 0, 8), common_term_22);
    result_by_q.set(c_Key1(8), common_term_22);
    // q = 9
    result_by_lpq.set(c_Key3(8, 0, 9), common_term_23);
    result_by_q.set(c_Key1(9), common_term_23);
    // q = 10
    result_by_lpq.set(c_Key3(8, 0, 10), common_term_24);
    result_by_q.set(c_Key1(10), common_term_24);
    // q = 11
    result_by_lpq.set(c_Key3(8, 0, 11), common_term_25);
    result_by_q.set(c_Key1(11), common_term_25);
    // q = 12
    result_by_lpq.set(c_Key3(8, 0, 12), common_term_26);
    result_by_q.set(c_Key1(12), common_term_26);
    // q = 13
    result_by_lpq.set(c_Key3(8, 0, 13), common_term_27);
    result_by_q.set(c_Key1(13), common_term_27);
    // q = 14
    result_by_lpq.set(c_Key3(8, 0, 14), common_term_28);
    result_by_q.set(c_Key1(14), common_term_28);
    // q = 15
    result_by_lpq.set(c_Key3(8, 0, 15), common_term_29);
    result_by_q.set(c_Key1(15), common_term_29);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 0), result_by_q);
    result_by_q.clear();

    // l , p = (8, 1).
    // q = -15
    result_by_lpq.set(c_Key3(8, 1, -15), common_term_30);
    result_by_q.set(c_Key1(-15), common_term_30);
    // q = -14
    result_by_lpq.set(c_Key3(8, 1, -14), common_term_31);
    result_by_q.set(c_Key1(-14), common_term_31);
    // q = -13
    result_by_lpq.set(c_Key3(8, 1, -13), common_term_32);
    result_by_q.set(c_Key1(-13), common_term_32);
    // q = -12
    result_by_lpq.set(c_Key3(8, 1, -12), common_term_33);
    result_by_q.set(c_Key1(-12), common_term_33);
    // q = -11
    result_by_lpq.set(c_Key3(8, 1, -11), common_term_34);
    result_by_q.set(c_Key1(-11), common_term_34);
    // q = -10
    result_by_lpq.set(c_Key3(8, 1, -10), common_term_35);
    result_by_q.set(c_Key1(-10), common_term_35);
    // q = -9
    result_by_lpq.set(c_Key3(8, 1, -9), common_term_36);
    result_by_q.set(c_Key1(-9), common_term_36);
    // q = -8
    result_by_lpq.set(c_Key3(8, 1, -8), common_term_37);
    result_by_q.set(c_Key1(-8), common_term_37);
    // q = -7
    result_by_lpq.set(c_Key3(8, 1, -7), common_term_38);
    result_by_q.set(c_Key1(-7), common_term_38);
    // q = -6
    result_by_lpq.set(c_Key3(8, 1, -6), common_term_39);
    result_by_q.set(c_Key1(-6), common_term_39);
    // q = -5
    result_by_lpq.set(c_Key3(8, 1, -5), common_term_40);
    result_by_q.set(c_Key1(-5), common_term_40);
    // q = -4
    result_by_lpq.set(c_Key3(8, 1, -4), common_term_41);
    result_by_q.set(c_Key1(-4), common_term_41);
    // q = -3
    result_by_lpq.set(c_Key3(8, 1, -3), common_term_42);
    result_by_q.set(c_Key1(-3), common_term_42);
    // q = -2
    result_by_lpq.set(c_Key3(8, 1, -2), common_term_43);
    result_by_q.set(c_Key1(-2), common_term_43);
    // q = -1
    result_by_lpq.set(c_Key3(8, 1, -1), common_term_44);
    result_by_q.set(c_Key1(-1), common_term_44);
    // q = 0
    result_by_lpq.set(c_Key3(8, 1, 0), common_term_45);
    result_by_q.set(c_Key1(0), common_term_45);
    // q = 1
    result_by_lpq.set(c_Key3(8, 1, 1), common_term_46);
    result_by_q.set(c_Key1(1), common_term_46);
    // q = 2
    result_by_lpq.set(c_Key3(8, 1, 2), common_term_47);
    result_by_q.set(c_Key1(2), common_term_47);
    // q = 3
    result_by_lpq.set(c_Key3(8, 1, 3), common_term_48);
    result_by_q.set(c_Key1(3), common_term_48);
    // q = 4
    result_by_lpq.set(c_Key3(8, 1, 4), common_term_49);
    result_by_q.set(c_Key1(4), common_term_49);
    // q = 5
    result_by_lpq.set(c_Key3(8, 1, 5), common_term_50);
    result_by_q.set(c_Key1(5), common_term_50);
    // q = 6
    result_by_lpq.set(c_Key3(8, 1, 6), common_term_51);
    result_by_q.set(c_Key1(6), common_term_51);
    // q = 7
    result_by_lpq.set(c_Key3(8, 1, 7), common_term_52);
    result_by_q.set(c_Key1(7), common_term_52);
    // q = 8
    result_by_lpq.set(c_Key3(8, 1, 8), common_term_53);
    result_by_q.set(c_Key1(8), common_term_53);
    // q = 9
    result_by_lpq.set(c_Key3(8, 1, 9), common_term_54);
    result_by_q.set(c_Key1(9), common_term_54);
    // q = 10
    result_by_lpq.set(c_Key3(8, 1, 10), common_term_55);
    result_by_q.set(c_Key1(10), common_term_55);
    // q = 11
    result_by_lpq.set(c_Key3(8, 1, 11), common_term_56);
    result_by_q.set(c_Key1(11), common_term_56);
    // q = 12
    result_by_lpq.set(c_Key3(8, 1, 12), common_term_57);
    result_by_q.set(c_Key1(12), common_term_57);
    // q = 13
    result_by_lpq.set(c_Key3(8, 1, 13), common_term_58);
    result_by_q.set(c_Key1(13), common_term_58);
    // q = 14
    result_by_lpq.set(c_Key3(8, 1, 14), common_term_59);
    result_by_q.set(c_Key1(14), common_term_59);
    // q = 15
    result_by_lpq.set(c_Key3(8, 1, 15), common_term_60);
    result_by_q.set(c_Key1(15), common_term_60);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 1), result_by_q);
    result_by_q.clear();

    // l , p = (8, 2).
    // q = -15
    result_by_lpq.set(c_Key3(8, 2, -15), common_term_61);
    result_by_q.set(c_Key1(-15), common_term_61);
    // q = -14
    result_by_lpq.set(c_Key3(8, 2, -14), common_term_62);
    result_by_q.set(c_Key1(-14), common_term_62);
    // q = -13
    result_by_lpq.set(c_Key3(8, 2, -13), common_term_63);
    result_by_q.set(c_Key1(-13), common_term_63);
    // q = -12
    result_by_lpq.set(c_Key3(8, 2, -12), common_term_64);
    result_by_q.set(c_Key1(-12), common_term_64);
    // q = -11
    result_by_lpq.set(c_Key3(8, 2, -11), common_term_65);
    result_by_q.set(c_Key1(-11), common_term_65);
    // q = -10
    result_by_lpq.set(c_Key3(8, 2, -10), common_term_66);
    result_by_q.set(c_Key1(-10), common_term_66);
    // q = -9
    result_by_lpq.set(c_Key3(8, 2, -9), common_term_67);
    result_by_q.set(c_Key1(-9), common_term_67);
    // q = -8
    result_by_lpq.set(c_Key3(8, 2, -8), common_term_68);
    result_by_q.set(c_Key1(-8), common_term_68);
    // q = -7
    result_by_lpq.set(c_Key3(8, 2, -7), common_term_69);
    result_by_q.set(c_Key1(-7), common_term_69);
    // q = -6
    result_by_lpq.set(c_Key3(8, 2, -6), common_term_70);
    result_by_q.set(c_Key1(-6), common_term_70);
    // q = -5
    result_by_lpq.set(c_Key3(8, 2, -5), common_term_71);
    result_by_q.set(c_Key1(-5), common_term_71);
    // q = -4
    result_by_lpq.set(c_Key3(8, 2, -4), common_term_72);
    result_by_q.set(c_Key1(-4), common_term_72);
    // q = -3
    result_by_lpq.set(c_Key3(8, 2, -3), common_term_73);
    result_by_q.set(c_Key1(-3), common_term_73);
    // q = -2
    result_by_lpq.set(c_Key3(8, 2, -2), common_term_74);
    result_by_q.set(c_Key1(-2), common_term_74);
    // q = -1
    result_by_lpq.set(c_Key3(8, 2, -1), common_term_75);
    result_by_q.set(c_Key1(-1), common_term_75);
    // q = 0
    result_by_lpq.set(c_Key3(8, 2, 0), common_term_76);
    result_by_q.set(c_Key1(0), common_term_76);
    // q = 1
    result_by_lpq.set(c_Key3(8, 2, 1), common_term_77);
    result_by_q.set(c_Key1(1), common_term_77);
    // q = 2
    result_by_lpq.set(c_Key3(8, 2, 2), common_term_78);
    result_by_q.set(c_Key1(2), common_term_78);
    // q = 3
    result_by_lpq.set(c_Key3(8, 2, 3), common_term_79);
    result_by_q.set(c_Key1(3), common_term_79);
    // q = 4
    result_by_lpq.set(c_Key3(8, 2, 4), common_term_80);
    result_by_q.set(c_Key1(4), common_term_80);
    // q = 5
    result_by_lpq.set(c_Key3(8, 2, 5), common_term_81);
    result_by_q.set(c_Key1(5), common_term_81);
    // q = 6
    result_by_lpq.set(c_Key3(8, 2, 6), common_term_82);
    result_by_q.set(c_Key1(6), common_term_82);
    // q = 7
    result_by_lpq.set(c_Key3(8, 2, 7), common_term_83);
    result_by_q.set(c_Key1(7), common_term_83);
    // q = 8
    result_by_lpq.set(c_Key3(8, 2, 8), common_term_84);
    result_by_q.set(c_Key1(8), common_term_84);
    // q = 9
    result_by_lpq.set(c_Key3(8, 2, 9), common_term_85);
    result_by_q.set(c_Key1(9), common_term_85);
    // q = 10
    result_by_lpq.set(c_Key3(8, 2, 10), common_term_86);
    result_by_q.set(c_Key1(10), common_term_86);
    // q = 11
    result_by_lpq.set(c_Key3(8, 2, 11), common_term_87);
    result_by_q.set(c_Key1(11), common_term_87);
    // q = 12
    result_by_lpq.set(c_Key3(8, 2, 12), common_term_88);
    result_by_q.set(c_Key1(12), common_term_88);
    // q = 13
    result_by_lpq.set(c_Key3(8, 2, 13), common_term_89);
    result_by_q.set(c_Key1(13), common_term_89);
    // q = 14
    result_by_lpq.set(c_Key3(8, 2, 14), common_term_90);
    result_by_q.set(c_Key1(14), common_term_90);
    // q = 15
    result_by_lpq.set(c_Key3(8, 2, 15), common_term_91);
    result_by_q.set(c_Key1(15), common_term_91);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 2), result_by_q);
    result_by_q.clear();

    // l , p = (8, 3).
    // q = -15
    result_by_lpq.set(c_Key3(8, 3, -15), common_term_92);
    result_by_q.set(c_Key1(-15), common_term_92);
    // q = -14
    result_by_lpq.set(c_Key3(8, 3, -14), common_term_93);
    result_by_q.set(c_Key1(-14), common_term_93);
    // q = -13
    result_by_lpq.set(c_Key3(8, 3, -13), common_term_94);
    result_by_q.set(c_Key1(-13), common_term_94);
    // q = -12
    result_by_lpq.set(c_Key3(8, 3, -12), common_term_95);
    result_by_q.set(c_Key1(-12), common_term_95);
    // q = -11
    result_by_lpq.set(c_Key3(8, 3, -11), common_term_96);
    result_by_q.set(c_Key1(-11), common_term_96);
    // q = -10
    result_by_lpq.set(c_Key3(8, 3, -10), common_term_97);
    result_by_q.set(c_Key1(-10), common_term_97);
    // q = -9
    result_by_lpq.set(c_Key3(8, 3, -9), common_term_98);
    result_by_q.set(c_Key1(-9), common_term_98);
    // q = -8
    result_by_lpq.set(c_Key3(8, 3, -8), common_term_99);
    result_by_q.set(c_Key1(-8), common_term_99);
    // q = -7
    result_by_lpq.set(c_Key3(8, 3, -7), common_term_100);
    result_by_q.set(c_Key1(-7), common_term_100);
    // q = -6
    result_by_lpq.set(c_Key3(8, 3, -6), common_term_101);
    result_by_q.set(c_Key1(-6), common_term_101);
    // q = -5
    result_by_lpq.set(c_Key3(8, 3, -5), common_term_102);
    result_by_q.set(c_Key1(-5), common_term_102);
    // q = -4
    result_by_lpq.set(c_Key3(8, 3, -4), common_term_103);
    result_by_q.set(c_Key1(-4), common_term_103);
    // q = -3
    result_by_lpq.set(c_Key3(8, 3, -3), common_term_104);
    result_by_q.set(c_Key1(-3), common_term_104);
    // q = -2
    result_by_lpq.set(c_Key3(8, 3, -2), common_term_105);
    result_by_q.set(c_Key1(-2), common_term_105);
    // q = -1
    result_by_lpq.set(c_Key3(8, 3, -1), common_term_106);
    result_by_q.set(c_Key1(-1), common_term_106);
    // q = 0
    result_by_lpq.set(c_Key3(8, 3, 0), common_term_107);
    result_by_q.set(c_Key1(0), common_term_107);
    // q = 1
    result_by_lpq.set(c_Key3(8, 3, 1), common_term_108);
    result_by_q.set(c_Key1(1), common_term_108);
    // q = 2
    result_by_lpq.set(c_Key3(8, 3, 2), common_term_109);
    result_by_q.set(c_Key1(2), common_term_109);
    // q = 3
    result_by_lpq.set(c_Key3(8, 3, 3), common_term_110);
    result_by_q.set(c_Key1(3), common_term_110);
    // q = 4
    result_by_lpq.set(c_Key3(8, 3, 4), common_term_111);
    result_by_q.set(c_Key1(4), common_term_111);
    // q = 5
    result_by_lpq.set(c_Key3(8, 3, 5), common_term_112);
    result_by_q.set(c_Key1(5), common_term_112);
    // q = 6
    result_by_lpq.set(c_Key3(8, 3, 6), common_term_113);
    result_by_q.set(c_Key1(6), common_term_113);
    // q = 7
    result_by_lpq.set(c_Key3(8, 3, 7), common_term_114);
    result_by_q.set(c_Key1(7), common_term_114);
    // q = 8
    result_by_lpq.set(c_Key3(8, 3, 8), common_term_115);
    result_by_q.set(c_Key1(8), common_term_115);
    // q = 9
    result_by_lpq.set(c_Key3(8, 3, 9), common_term_116);
    result_by_q.set(c_Key1(9), common_term_116);
    // q = 10
    result_by_lpq.set(c_Key3(8, 3, 10), common_term_117);
    result_by_q.set(c_Key1(10), common_term_117);
    // q = 11
    result_by_lpq.set(c_Key3(8, 3, 11), common_term_118);
    result_by_q.set(c_Key1(11), common_term_118);
    // q = 12
    result_by_lpq.set(c_Key3(8, 3, 12), common_term_119);
    result_by_q.set(c_Key1(12), common_term_119);
    // q = 13
    result_by_lpq.set(c_Key3(8, 3, 13), common_term_120);
    result_by_q.set(c_Key1(13), common_term_120);
    // q = 14
    result_by_lpq.set(c_Key3(8, 3, 14), common_term_121);
    result_by_q.set(c_Key1(14), common_term_121);
    // q = 15
    result_by_lpq.set(c_Key3(8, 3, 15), common_term_122);
    result_by_q.set(c_Key1(15), common_term_122);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 3), result_by_q);
    result_by_q.clear();

    // l , p = (8, 4).
    // q = -15
    result_by_lpq.set(c_Key3(8, 4, -15), common_term_123);
    result_by_q.set(c_Key1(-15), common_term_123);
    // q = -14
    result_by_lpq.set(c_Key3(8, 4, -14), common_term_124);
    result_by_q.set(c_Key1(-14), common_term_124);
    // q = -13
    result_by_lpq.set(c_Key3(8, 4, -13), common_term_125);
    result_by_q.set(c_Key1(-13), common_term_125);
    // q = -12
    result_by_lpq.set(c_Key3(8, 4, -12), common_term_126);
    result_by_q.set(c_Key1(-12), common_term_126);
    // q = -11
    result_by_lpq.set(c_Key3(8, 4, -11), common_term_127);
    result_by_q.set(c_Key1(-11), common_term_127);
    // q = -10
    result_by_lpq.set(c_Key3(8, 4, -10), common_term_128);
    result_by_q.set(c_Key1(-10), common_term_128);
    // q = -9
    result_by_lpq.set(c_Key3(8, 4, -9), common_term_129);
    result_by_q.set(c_Key1(-9), common_term_129);
    // q = -8
    result_by_lpq.set(c_Key3(8, 4, -8), common_term_130);
    result_by_q.set(c_Key1(-8), common_term_130);
    // q = -7
    result_by_lpq.set(c_Key3(8, 4, -7), common_term_131);
    result_by_q.set(c_Key1(-7), common_term_131);
    // q = -6
    result_by_lpq.set(c_Key3(8, 4, -6), common_term_132);
    result_by_q.set(c_Key1(-6), common_term_132);
    // q = -5
    result_by_lpq.set(c_Key3(8, 4, -5), common_term_133);
    result_by_q.set(c_Key1(-5), common_term_133);
    // q = -4
    result_by_lpq.set(c_Key3(8, 4, -4), common_term_134);
    result_by_q.set(c_Key1(-4), common_term_134);
    // q = -3
    result_by_lpq.set(c_Key3(8, 4, -3), common_term_135);
    result_by_q.set(c_Key1(-3), common_term_135);
    // q = -2
    result_by_lpq.set(c_Key3(8, 4, -2), common_term_136);
    result_by_q.set(c_Key1(-2), common_term_136);
    // q = -1
    result_by_lpq.set(c_Key3(8, 4, -1), common_term_137);
    result_by_q.set(c_Key1(-1), common_term_137);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -7.5)*(2.1875*eccentricity_6 + 13.125*eccentricity_4 + 10.5*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(8, 4, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(8, 4, 1), common_term_137);
    result_by_q.set(c_Key1(1), common_term_137);
    // q = 2
    result_by_lpq.set(c_Key3(8, 4, 2), common_term_136);
    result_by_q.set(c_Key1(2), common_term_136);
    // q = 3
    result_by_lpq.set(c_Key3(8, 4, 3), common_term_135);
    result_by_q.set(c_Key1(3), common_term_135);
    // q = 4
    result_by_lpq.set(c_Key3(8, 4, 4), common_term_134);
    result_by_q.set(c_Key1(4), common_term_134);
    // q = 5
    result_by_lpq.set(c_Key3(8, 4, 5), common_term_133);
    result_by_q.set(c_Key1(5), common_term_133);
    // q = 6
    result_by_lpq.set(c_Key3(8, 4, 6), common_term_132);
    result_by_q.set(c_Key1(6), common_term_132);
    // q = 7
    result_by_lpq.set(c_Key3(8, 4, 7), common_term_131);
    result_by_q.set(c_Key1(7), common_term_131);
    // q = 8
    result_by_lpq.set(c_Key3(8, 4, 8), common_term_130);
    result_by_q.set(c_Key1(8), common_term_130);
    // q = 9
    result_by_lpq.set(c_Key3(8, 4, 9), common_term_129);
    result_by_q.set(c_Key1(9), common_term_129);
    // q = 10
    result_by_lpq.set(c_Key3(8, 4, 10), common_term_128);
    result_by_q.set(c_Key1(10), common_term_128);
    // q = 11
    result_by_lpq.set(c_Key3(8, 4, 11), common_term_127);
    result_by_q.set(c_Key1(11), common_term_127);
    // q = 12
    result_by_lpq.set(c_Key3(8, 4, 12), common_term_126);
    result_by_q.set(c_Key1(12), common_term_126);
    // q = 13
    result_by_lpq.set(c_Key3(8, 4, 13), common_term_125);
    result_by_q.set(c_Key1(13), common_term_125);
    // q = 14
    result_by_lpq.set(c_Key3(8, 4, 14), common_term_124);
    result_by_q.set(c_Key1(14), common_term_124);
    // q = 15
    result_by_lpq.set(c_Key3(8, 4, 15), common_term_123);
    result_by_q.set(c_Key1(15), common_term_123);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 4), result_by_q);
    result_by_q.clear();

    // l , p = (8, 5).
    // q = -15
    result_by_lpq.set(c_Key3(8, 5, -15), common_term_122);
    result_by_q.set(c_Key1(-15), common_term_122);
    // q = -14
    result_by_lpq.set(c_Key3(8, 5, -14), common_term_121);
    result_by_q.set(c_Key1(-14), common_term_121);
    // q = -13
    result_by_lpq.set(c_Key3(8, 5, -13), common_term_120);
    result_by_q.set(c_Key1(-13), common_term_120);
    // q = -12
    result_by_lpq.set(c_Key3(8, 5, -12), common_term_119);
    result_by_q.set(c_Key1(-12), common_term_119);
    // q = -11
    result_by_lpq.set(c_Key3(8, 5, -11), common_term_118);
    result_by_q.set(c_Key1(-11), common_term_118);
    // q = -10
    result_by_lpq.set(c_Key3(8, 5, -10), common_term_117);
    result_by_q.set(c_Key1(-10), common_term_117);
    // q = -9
    result_by_lpq.set(c_Key3(8, 5, -9), common_term_116);
    result_by_q.set(c_Key1(-9), common_term_116);
    // q = -8
    result_by_lpq.set(c_Key3(8, 5, -8), common_term_115);
    result_by_q.set(c_Key1(-8), common_term_115);
    // q = -7
    result_by_lpq.set(c_Key3(8, 5, -7), common_term_114);
    result_by_q.set(c_Key1(-7), common_term_114);
    // q = -6
    result_by_lpq.set(c_Key3(8, 5, -6), common_term_113);
    result_by_q.set(c_Key1(-6), common_term_113);
    // q = -5
    result_by_lpq.set(c_Key3(8, 5, -5), common_term_112);
    result_by_q.set(c_Key1(-5), common_term_112);
    // q = -4
    result_by_lpq.set(c_Key3(8, 5, -4), common_term_111);
    result_by_q.set(c_Key1(-4), common_term_111);
    // q = -3
    result_by_lpq.set(c_Key3(8, 5, -3), common_term_110);
    result_by_q.set(c_Key1(-3), common_term_110);
    // q = -2
    result_by_lpq.set(c_Key3(8, 5, -2), common_term_109);
    result_by_q.set(c_Key1(-2), common_term_109);
    // q = -1
    result_by_lpq.set(c_Key3(8, 5, -1), common_term_108);
    result_by_q.set(c_Key1(-1), common_term_108);
    // q = 0
    result_by_lpq.set(c_Key3(8, 5, 0), common_term_107);
    result_by_q.set(c_Key1(0), common_term_107);
    // q = 1
    result_by_lpq.set(c_Key3(8, 5, 1), common_term_106);
    result_by_q.set(c_Key1(1), common_term_106);
    // q = 2
    result_by_lpq.set(c_Key3(8, 5, 2), common_term_105);
    result_by_q.set(c_Key1(2), common_term_105);
    // q = 3
    result_by_lpq.set(c_Key3(8, 5, 3), common_term_104);
    result_by_q.set(c_Key1(3), common_term_104);
    // q = 4
    result_by_lpq.set(c_Key3(8, 5, 4), common_term_103);
    result_by_q.set(c_Key1(4), common_term_103);
    // q = 5
    result_by_lpq.set(c_Key3(8, 5, 5), common_term_102);
    result_by_q.set(c_Key1(5), common_term_102);
    // q = 6
    result_by_lpq.set(c_Key3(8, 5, 6), common_term_101);
    result_by_q.set(c_Key1(6), common_term_101);
    // q = 7
    result_by_lpq.set(c_Key3(8, 5, 7), common_term_100);
    result_by_q.set(c_Key1(7), common_term_100);
    // q = 8
    result_by_lpq.set(c_Key3(8, 5, 8), common_term_99);
    result_by_q.set(c_Key1(8), common_term_99);
    // q = 9
    result_by_lpq.set(c_Key3(8, 5, 9), common_term_98);
    result_by_q.set(c_Key1(9), common_term_98);
    // q = 10
    result_by_lpq.set(c_Key3(8, 5, 10), common_term_97);
    result_by_q.set(c_Key1(10), common_term_97);
    // q = 11
    result_by_lpq.set(c_Key3(8, 5, 11), common_term_96);
    result_by_q.set(c_Key1(11), common_term_96);
    // q = 12
    result_by_lpq.set(c_Key3(8, 5, 12), common_term_95);
    result_by_q.set(c_Key1(12), common_term_95);
    // q = 13
    result_by_lpq.set(c_Key3(8, 5, 13), common_term_94);
    result_by_q.set(c_Key1(13), common_term_94);
    // q = 14
    result_by_lpq.set(c_Key3(8, 5, 14), common_term_93);
    result_by_q.set(c_Key1(14), common_term_93);
    // q = 15
    result_by_lpq.set(c_Key3(8, 5, 15), common_term_92);
    result_by_q.set(c_Key1(15), common_term_92);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 5), result_by_q);
    result_by_q.clear();

    // l , p = (8, 6).
    // q = -15
    result_by_lpq.set(c_Key3(8, 6, -15), common_term_91);
    result_by_q.set(c_Key1(-15), common_term_91);
    // q = -14
    result_by_lpq.set(c_Key3(8, 6, -14), common_term_90);
    result_by_q.set(c_Key1(-14), common_term_90);
    // q = -13
    result_by_lpq.set(c_Key3(8, 6, -13), common_term_89);
    result_by_q.set(c_Key1(-13), common_term_89);
    // q = -12
    result_by_lpq.set(c_Key3(8, 6, -12), common_term_88);
    result_by_q.set(c_Key1(-12), common_term_88);
    // q = -11
    result_by_lpq.set(c_Key3(8, 6, -11), common_term_87);
    result_by_q.set(c_Key1(-11), common_term_87);
    // q = -10
    result_by_lpq.set(c_Key3(8, 6, -10), common_term_86);
    result_by_q.set(c_Key1(-10), common_term_86);
    // q = -9
    result_by_lpq.set(c_Key3(8, 6, -9), common_term_85);
    result_by_q.set(c_Key1(-9), common_term_85);
    // q = -8
    result_by_lpq.set(c_Key3(8, 6, -8), common_term_84);
    result_by_q.set(c_Key1(-8), common_term_84);
    // q = -7
    result_by_lpq.set(c_Key3(8, 6, -7), common_term_83);
    result_by_q.set(c_Key1(-7), common_term_83);
    // q = -6
    result_by_lpq.set(c_Key3(8, 6, -6), common_term_82);
    result_by_q.set(c_Key1(-6), common_term_82);
    // q = -5
    result_by_lpq.set(c_Key3(8, 6, -5), common_term_81);
    result_by_q.set(c_Key1(-5), common_term_81);
    // q = -4
    result_by_lpq.set(c_Key3(8, 6, -4), common_term_80);
    result_by_q.set(c_Key1(-4), common_term_80);
    // q = -3
    result_by_lpq.set(c_Key3(8, 6, -3), common_term_79);
    result_by_q.set(c_Key1(-3), common_term_79);
    // q = -2
    result_by_lpq.set(c_Key3(8, 6, -2), common_term_78);
    result_by_q.set(c_Key1(-2), common_term_78);
    // q = -1
    result_by_lpq.set(c_Key3(8, 6, -1), common_term_77);
    result_by_q.set(c_Key1(-1), common_term_77);
    // q = 0
    result_by_lpq.set(c_Key3(8, 6, 0), common_term_76);
    result_by_q.set(c_Key1(0), common_term_76);
    // q = 1
    result_by_lpq.set(c_Key3(8, 6, 1), common_term_75);
    result_by_q.set(c_Key1(1), common_term_75);
    // q = 2
    result_by_lpq.set(c_Key3(8, 6, 2), common_term_74);
    result_by_q.set(c_Key1(2), common_term_74);
    // q = 3
    result_by_lpq.set(c_Key3(8, 6, 3), common_term_73);
    result_by_q.set(c_Key1(3), common_term_73);
    // q = 4
    result_by_lpq.set(c_Key3(8, 6, 4), common_term_72);
    result_by_q.set(c_Key1(4), common_term_72);
    // q = 5
    result_by_lpq.set(c_Key3(8, 6, 5), common_term_71);
    result_by_q.set(c_Key1(5), common_term_71);
    // q = 6
    result_by_lpq.set(c_Key3(8, 6, 6), common_term_70);
    result_by_q.set(c_Key1(6), common_term_70);
    // q = 7
    result_by_lpq.set(c_Key3(8, 6, 7), common_term_69);
    result_by_q.set(c_Key1(7), common_term_69);
    // q = 8
    result_by_lpq.set(c_Key3(8, 6, 8), common_term_68);
    result_by_q.set(c_Key1(8), common_term_68);
    // q = 9
    result_by_lpq.set(c_Key3(8, 6, 9), common_term_67);
    result_by_q.set(c_Key1(9), common_term_67);
    // q = 10
    result_by_lpq.set(c_Key3(8, 6, 10), common_term_66);
    result_by_q.set(c_Key1(10), common_term_66);
    // q = 11
    result_by_lpq.set(c_Key3(8, 6, 11), common_term_65);
    result_by_q.set(c_Key1(11), common_term_65);
    // q = 12
    result_by_lpq.set(c_Key3(8, 6, 12), common_term_64);
    result_by_q.set(c_Key1(12), common_term_64);
    // q = 13
    result_by_lpq.set(c_Key3(8, 6, 13), common_term_63);
    result_by_q.set(c_Key1(13), common_term_63);
    // q = 14
    result_by_lpq.set(c_Key3(8, 6, 14), common_term_62);
    result_by_q.set(c_Key1(14), common_term_62);
    // q = 15
    result_by_lpq.set(c_Key3(8, 6, 15), common_term_61);
    result_by_q.set(c_Key1(15), common_term_61);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 6), result_by_q);
    result_by_q.clear();

    // l , p = (8, 7).
    // q = -15
    result_by_lpq.set(c_Key3(8, 7, -15), common_term_60);
    result_by_q.set(c_Key1(-15), common_term_60);
    // q = -14
    result_by_lpq.set(c_Key3(8, 7, -14), common_term_59);
    result_by_q.set(c_Key1(-14), common_term_59);
    // q = -13
    result_by_lpq.set(c_Key3(8, 7, -13), common_term_58);
    result_by_q.set(c_Key1(-13), common_term_58);
    // q = -12
    result_by_lpq.set(c_Key3(8, 7, -12), common_term_57);
    result_by_q.set(c_Key1(-12), common_term_57);
    // q = -11
    result_by_lpq.set(c_Key3(8, 7, -11), common_term_56);
    result_by_q.set(c_Key1(-11), common_term_56);
    // q = -10
    result_by_lpq.set(c_Key3(8, 7, -10), common_term_55);
    result_by_q.set(c_Key1(-10), common_term_55);
    // q = -9
    result_by_lpq.set(c_Key3(8, 7, -9), common_term_54);
    result_by_q.set(c_Key1(-9), common_term_54);
    // q = -8
    result_by_lpq.set(c_Key3(8, 7, -8), common_term_53);
    result_by_q.set(c_Key1(-8), common_term_53);
    // q = -7
    result_by_lpq.set(c_Key3(8, 7, -7), common_term_52);
    result_by_q.set(c_Key1(-7), common_term_52);
    // q = -6
    result_by_lpq.set(c_Key3(8, 7, -6), common_term_51);
    result_by_q.set(c_Key1(-6), common_term_51);
    // q = -5
    result_by_lpq.set(c_Key3(8, 7, -5), common_term_50);
    result_by_q.set(c_Key1(-5), common_term_50);
    // q = -4
    result_by_lpq.set(c_Key3(8, 7, -4), common_term_49);
    result_by_q.set(c_Key1(-4), common_term_49);
    // q = -3
    result_by_lpq.set(c_Key3(8, 7, -3), common_term_48);
    result_by_q.set(c_Key1(-3), common_term_48);
    // q = -2
    result_by_lpq.set(c_Key3(8, 7, -2), common_term_47);
    result_by_q.set(c_Key1(-2), common_term_47);
    // q = -1
    result_by_lpq.set(c_Key3(8, 7, -1), common_term_46);
    result_by_q.set(c_Key1(-1), common_term_46);
    // q = 0
    result_by_lpq.set(c_Key3(8, 7, 0), common_term_45);
    result_by_q.set(c_Key1(0), common_term_45);
    // q = 1
    result_by_lpq.set(c_Key3(8, 7, 1), common_term_44);
    result_by_q.set(c_Key1(1), common_term_44);
    // q = 2
    result_by_lpq.set(c_Key3(8, 7, 2), common_term_43);
    result_by_q.set(c_Key1(2), common_term_43);
    // q = 3
    result_by_lpq.set(c_Key3(8, 7, 3), common_term_42);
    result_by_q.set(c_Key1(3), common_term_42);
    // q = 4
    result_by_lpq.set(c_Key3(8, 7, 4), common_term_41);
    result_by_q.set(c_Key1(4), common_term_41);
    // q = 5
    result_by_lpq.set(c_Key3(8, 7, 5), common_term_40);
    result_by_q.set(c_Key1(5), common_term_40);
    // q = 6
    result_by_lpq.set(c_Key3(8, 7, 6), common_term_39);
    result_by_q.set(c_Key1(6), common_term_39);
    // q = 7
    result_by_lpq.set(c_Key3(8, 7, 7), common_term_38);
    result_by_q.set(c_Key1(7), common_term_38);
    // q = 8
    result_by_lpq.set(c_Key3(8, 7, 8), common_term_37);
    result_by_q.set(c_Key1(8), common_term_37);
    // q = 9
    result_by_lpq.set(c_Key3(8, 7, 9), common_term_36);
    result_by_q.set(c_Key1(9), common_term_36);
    // q = 10
    result_by_lpq.set(c_Key3(8, 7, 10), common_term_35);
    result_by_q.set(c_Key1(10), common_term_35);
    // q = 11
    result_by_lpq.set(c_Key3(8, 7, 11), common_term_34);
    result_by_q.set(c_Key1(11), common_term_34);
    // q = 12
    result_by_lpq.set(c_Key3(8, 7, 12), common_term_33);
    result_by_q.set(c_Key1(12), common_term_33);
    // q = 13
    result_by_lpq.set(c_Key3(8, 7, 13), common_term_32);
    result_by_q.set(c_Key1(13), common_term_32);
    // q = 14
    result_by_lpq.set(c_Key3(8, 7, 14), common_term_31);
    result_by_q.set(c_Key1(14), common_term_31);
    // q = 15
    result_by_lpq.set(c_Key3(8, 7, 15), common_term_30);
    result_by_q.set(c_Key1(15), common_term_30);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 7), result_by_q);
    result_by_q.clear();

    // l , p = (8, 8).
    // q = -15
    result_by_lpq.set(c_Key3(8, 8, -15), common_term_29);
    result_by_q.set(c_Key1(-15), common_term_29);
    // q = -14
    result_by_lpq.set(c_Key3(8, 8, -14), common_term_28);
    result_by_q.set(c_Key1(-14), common_term_28);
    // q = -13
    result_by_lpq.set(c_Key3(8, 8, -13), common_term_27);
    result_by_q.set(c_Key1(-13), common_term_27);
    // q = -12
    result_by_lpq.set(c_Key3(8, 8, -12), common_term_26);
    result_by_q.set(c_Key1(-12), common_term_26);
    // q = -11
    result_by_lpq.set(c_Key3(8, 8, -11), common_term_25);
    result_by_q.set(c_Key1(-11), common_term_25);
    // q = -10
    result_by_lpq.set(c_Key3(8, 8, -10), common_term_24);
    result_by_q.set(c_Key1(-10), common_term_24);
    // q = -9
    result_by_lpq.set(c_Key3(8, 8, -9), common_term_23);
    result_by_q.set(c_Key1(-9), common_term_23);
    // q = -8
    result_by_lpq.set(c_Key3(8, 8, -8), common_term_22);
    result_by_q.set(c_Key1(-8), common_term_22);
    // q = -7
    result_by_lpq.set(c_Key3(8, 8, -7), common_term_21);
    result_by_q.set(c_Key1(-7), common_term_21);
    // q = -6
    result_by_lpq.set(c_Key3(8, 8, -6), common_term_20);
    result_by_q.set(c_Key1(-6), common_term_20);
    // q = -5
    result_by_lpq.set(c_Key3(8, 8, -5), common_term_19);
    result_by_q.set(c_Key1(-5), common_term_19);
    // q = -4
    result_by_lpq.set(c_Key3(8, 8, -4), common_term_18);
    result_by_q.set(c_Key1(-4), common_term_18);
    // q = -3
    result_by_lpq.set(c_Key3(8, 8, -3), common_term_17);
    result_by_q.set(c_Key1(-3), common_term_17);
    // q = -2
    result_by_lpq.set(c_Key3(8, 8, -2), common_term_16);
    result_by_q.set(c_Key1(-2), common_term_16);
    // q = -1
    result_by_lpq.set(c_Key3(8, 8, -1), common_term_15);
    result_by_q.set(c_Key1(-1), common_term_15);
    // q = 0
    result_by_lpq.set(c_Key3(8, 8, 0), common_term_14);
    result_by_q.set(c_Key1(0), common_term_14);
    // q = 1
    result_by_lpq.set(c_Key3(8, 8, 1), common_term_13);
    result_by_q.set(c_Key1(1), common_term_13);
    // q = 2
    result_by_lpq.set(c_Key3(8, 8, 2), common_term_12);
    result_by_q.set(c_Key1(2), common_term_12);
    // q = 3
    result_by_lpq.set(c_Key3(8, 8, 3), common_term_11);
    result_by_q.set(c_Key1(3), common_term_11);
    // q = 4
    result_by_lpq.set(c_Key3(8, 8, 4), common_term_10);
    result_by_q.set(c_Key1(4), common_term_10);
    // q = 5
    result_by_lpq.set(c_Key3(8, 8, 5), common_term_9);
    result_by_q.set(c_Key1(5), common_term_9);
    // q = 6
    result_by_lpq.set(c_Key3(8, 8, 6), common_term_8);
    result_by_q.set(c_Key1(6), common_term_8);
    // q = 7
    result_by_lpq.set(c_Key3(8, 8, 7), common_term_7);
    result_by_q.set(c_Key1(7), common_term_7);
    // q = 9
    result_by_lpq.set(c_Key3(8, 8, 9), common_term_6);
    result_by_q.set(c_Key1(9), common_term_6);
    // q = 10
    result_by_lpq.set(c_Key3(8, 8, 10), common_term_5);
    result_by_q.set(c_Key1(10), common_term_5);
    // q = 11
    result_by_lpq.set(c_Key3(8, 8, 11), common_term_4);
    result_by_q.set(c_Key1(11), common_term_4);
    // q = 12
    result_by_lpq.set(c_Key3(8, 8, 12), common_term_3);
    result_by_q.set(c_Key1(12), common_term_3);
    // q = 13
    result_by_lpq.set(c_Key3(8, 8, 13), common_term_2);
    result_by_q.set(c_Key1(13), common_term_2);
    // q = 14
    result_by_lpq.set(c_Key3(8, 8, 14), common_term_1);
    result_by_q.set(c_Key1(14), common_term_1);
    // q = 15
    result_by_lpq.set(c_Key3(8, 8, 15), common_term_0);
    result_by_q.set(c_Key1(15), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 8), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l8_e20(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 8.
    // Functions are taylor expanded and truncated at eccentricity^20.

    // Functions Calculated for l = 8.

    c_IntMap<c_Key3, double> result_by_lpq(367);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(9);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_5 = eccentricity_4 * eccentricity;
    double eccentricity_10 = eccentricity_5 * eccentricity_5;
    double eccentricity_20 = eccentricity_10 * eccentricity_10;
    double eccentricity_8 = eccentricity_4 * eccentricity_4;
    double eccentricity_9 = eccentricity_8 * eccentricity;
    double eccentricity_18 = eccentricity_9 * eccentricity_9;
    double eccentricity_19 = eccentricity_18 * eccentricity;
    double eccentricity_16 = eccentricity_8 * eccentricity_8;
    double eccentricity_17 = eccentricity_16 * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double eccentricity_6 = eccentricity_3 * eccentricity_3;
    double eccentricity_7 = eccentricity_6 * eccentricity;
    double eccentricity_14 = eccentricity_7 * eccentricity_7;
    double eccentricity_15 = eccentricity_14 * eccentricity;
    double eccentricity_12 = eccentricity_6 * eccentricity_6;
    double eccentricity_13 = eccentricity_12 * eccentricity;
    double eccentricity_11 = eccentricity_10 * eccentricity;
    double common_term_0 = 0.0015027972469812364*eccentricity_20;
    double common_term_1 = 0.00095895113867134476*eccentricity_19;
    double common_term_2 = 0.0016620394442409114*eccentricity_20 + 0.00059582546114296823*eccentricity_18;
    double common_term_3 = 0.0010284473608455816*eccentricity_19 + 0.00035772082116368055*eccentricity_17;
    double common_term_4 = 0.001154850769402763*eccentricity_20 + 0.00060375584807579828*eccentricity_18 + 0.00020527698834577141*eccentricity_16;
    double common_term_5 = 0.00063178992802255527*eccentricity_19 + 0.00033065450749126379*eccentricity_17 + 0.00011079522764105174*eccentricity_15;
    double common_term_6 = 0.00047778264462198778*eccentricity_20 + 0.00031238315604596297*eccentricity_18 + 0.00016459266180248323*eccentricity_16 + 5.4864220600827744e-5*eccentricity_14;
    double common_term_7 = 0.00020109298600167295*eccentricity_19 + 0.00013355485897599823*eccentricity_17 + 7.1362201676549913e-5*eccentricity_15 + 2.3929840083154462e-5*eccentricity_13;
    double common_term_8 = 8.7488333079162003e-5*eccentricity_20 + 6.7255756850113111e-5*eccentricity_18 + 4.5692109184172676e-5*eccentricity_16 + 2.4995580551136107e-5*eccentricity_14 + 8.5511196622307733e-6*eccentricity_12;
    double common_term_9 = 1.9391225459910339e-5*eccentricity_19 + 1.5320777605973961e-5*eccentricity_17 + 1.074487572306996e-5*eccentricity_15 + 6.0945362239689022e-6*eccentricity_13 + 2.1669462129667208e-6*eccentricity_11;
    double common_term_10 = 2.3070255827124489e-6*eccentricity_20 + 2.0225001171698427e-6*eccentricity_18 + 1.6566241971334564e-6*eccentricity_16 + 1.21398341884453e-6*eccentricity_14 + 7.2651114317780984e-7*eccentricity_12 + 2.7557319223985891e-7*eccentricity_10;
    double common_term_11 = 3.2432832271015019e-8*eccentricity_19 + 2.9533511844736341e-8*eccentricity_17 + 2.5365045247212074e-8*eccentricity_15 + 1.9744743382642146e-8*eccentricity_13 + 1.2782936163470018e-8*eccentricity_11 + 5.3822889109347443e-9*eccentricity_9;
    double common_term_12 = -2.6458880099528992e-6*eccentricity_19 - 2.7311478440815041e-6*eccentricity_17 - 2.7806556216540593e-6*eccentricity_15 - 2.7637997492140423e-6*eccentricity_13 - 2.627229774650022e-6*eccentricity_11 - 2.2767082093253968e-6*eccentricity_9 - 1.5500992063492063e-6*eccentricity_7;
    double common_term_13 = 0.00052702803953777072*eccentricity_20 + 0.0005809993713428561*eccentricity_18 + 0.00064401434149698039*eccentricity_16 + 0.0007178452013521458*eccentricity_14 + 0.00080421443268665491*eccentricity_12 + 0.00090525793650793651*eccentricity_10 + 0.00099206349206349206*eccentricity_8 + 0.0013888888888888889*eccentricity_6;
    double common_term_14 = 0.0023067077318973604*eccentricity_19 + 0.0024204935448897349*eccentricity_17 + 0.0025082426241465977*eccentricity_15 + 0.0025113752910069057*eccentricity_13 + 0.0027588435581752232*eccentricity_11 - 0.0015537806919642857*eccentricity_9 + 0.02373046875*eccentricity_7 - 0.06328125*eccentricity_5;
    double common_term_15 = 0.0030530492484801588*eccentricity_20 + 0.0024712432017326197*eccentricity_18 + 0.0016454169116206153*eccentricity_16 - 0.0020223398001175779*eccentricity_14 + 0.016881613756613757*eccentricity_12 - 0.14656084656084656*eccentricity_10 + 0.58888888888888889*eccentricity_8 - 1.3333333333333333*eccentricity_6 + 0.66666666666666667*eccentricity_4;
    double common_term_16 = -0.009379638490624503*eccentricity_19 - 0.007832995989643554*eccentricity_17 - 0.088458628195326818*eccentricity_15 + 0.4550276136902905*eccentricity_13 - 2.6191353167175616*eccentricity_11 + 8.5598980938946759*eccentricity_9 - 16.215006510416667*eccentricity_7 + 11.881510416666667*eccentricity_5 - 2.6041666666666667*eccentricity_3;
    double common_term_17 = -0.067968305994716779*eccentricity_20 + 0.09161074144013074*eccentricity_18 - 1.3278838488520408*eccentricity_16 + 7.017890625*eccentricity_14 - 30.029263392857143*eccentricity_12 + 83.11640625*eccentricity_10 - 139.6125*eccentricity_8 + 115.3125*eccentricity_6 - 40.5*eccentricity_4 + 4.5*eccentricity_2;
    double common_term_18 = 2.2728819790508349*eccentricity_19 - 15.566054155525454*eccentricity_17 + 72.291540367743115*eccentricity_15 - 253.89152100763203*eccentricity_13 + 610.67264824196144*eccentricity_11 - 935.56048651801215*eccentricity_9 + 800.71652560763889*eccentricity_7 - 343.30989583333333*eccentricity_5 + 63.4375*eccentricity_3 - 3.5*eccentricity;
    double common_term_19 = 26.806480209488636*eccentricity_20 - 141.73836034686102*eccentricity_18 + 567.20383346376926*eccentricity_16 - 1710.8899369724742*eccentricity_14 + 3655.4365760030864*eccentricity_12 - 5182.6406944444444*eccentricity_10 + 4439.1657986111111*eccentricity_8 - 2099.8611111111111*eccentricity_6 + 490.75*eccentricity_4 - 46.0*eccentricity_2 + 1.0;
    double common_term_20 = -1041.687894251133*eccentricity_19 + 3638.8731086200601*eccentricity_17 - 9678.2896001151265*eccentricity_15 + 18682.353862372807*eccentricity_13 - 24732.193187713623*eccentricity_11 + 20854.591314697266*eccentricity_9 - 10362.46142578125*eccentricity_7 + 2757.0234375*eccentricity_5 - 333.5625*eccentricity_3 + 12.5*eccentricity;
    double common_term_21 = -6439.3943107803421*eccentricity_20 + 19954.730070049046*eccentricity_18 - 47644.057747626845*eccentricity_16 + 84165.626218550209*eccentricity_14 - 104690.64519262566*eccentricity_12 + 86174.851707175926*eccentricity_10 - 43779.236111111111*eccentricity_8 + 12571.145833333333*eccentricity_6 - 1764.1666666666667*eccentricity_4 + 86.5*eccentricity_2;
    double common_term_22 = 96429.735024918023*eccentricity_19 - 209443.63208524453*eccentricity_17 + 341959.49825736394*eccentricity_15 - 401594.35779903985*eccentricity_13 + 321267.57119522902*eccentricity_11 - 164205.81790861907*eccentricity_9 + 49339.254264322917*eccentricity_7 - 7600.3268229166667*eccentricity_5 + 437.72916666666667*eccentricity_3;
    double common_term_23 = 419817.76288135779*eccentricity_20 - 838067.51362725026*eccentricity_18 + 1274506.6427682757*eccentricity_16 - 1418742.2822042411*eccentricity_14 + 1100497.1265904018*eccentricity_12 - 560229.34285714286*eccentricity_10 + 172758.65625*eccentricity_8 - 28273.05*eccentricity_6 + 1808.25*eccentricity_4;
    double common_term_24 = -3097192.3517843532*eccentricity_19 + 4414728.9869249735*eccentricity_17 - 4673634.5473623291*eccentricity_15 + 3511333.9634982762*eccentricity_13 - 1768754.423578666*eccentricity_11 + 552633.89492439391*eccentricity_9 - 94090.807693142361*eccentricity_7 + 6461.7950520833333*eccentricity_5;
    double common_term_25 = -10692166.213784524*eccentricity_20 + 14358498.473378641*eccentricity_18 - 14497460.482039825*eccentricity_16 + 10545772.889832256*eccentricity_14 - 5233760.1698913323*eccentricity_12 + 1642304.3425347222*eccentricity_10 - 286780.58888888889*eccentricity_8 + 20690.784722222222*eccentricity_6;
    double common_term_26 = 44209683.159485378*eccentricity_19 - 42679756.898512769*eccentricity_17 + 30060994.488160975*eccentricity_15 - 14656166.780386141*eccentricity_13 + 4590201.1259119851*eccentricity_11 - 813848.52869306292*eccentricity_9 + 60764.489606584821*eccentricity_7;
    double common_term_27 = 129725595.59089078*eccentricity_20 - 120010190.61575461*eccentricity_18 + 81872132.247217504*eccentricity_16 - 39137278.192853393*eccentricity_14 + 12179584.482504685*eccentricity_12 - 2176625.8677358907*eccentricity_10 + 166371.02557043651*eccentricity_8;
    double common_term_28 = -324017885.75938178*eccentricity_19 + 214208189.52493535*eccentricity_17 - 100269202.88532857*eccentricity_15 + 30904655.242274058*eccentricity_13 - 5536752.3833186226*eccentricity_11 + 429802.86250706695*eccentricity_9;
    double common_term_29 = -843695113.25735776*eccentricity_20 + 540825026.83907737*eccentricity_18 - 247686910.06394587*eccentricity_16 + 75427724.874900061*eccentricity_14 - 13491473.585716315*eccentricity_12 + 1057233.6862946429*eccentricity_10;
    double common_term_30 = 1322620254.9455328*eccentricity_19 - 592345218.95010126*eccentricity_17 + 177910875.05159583*eccentricity_15 - 31671129.351751119*eccentricity_13 + 2493744.2467064668*eccentricity_11;
    double common_term_31 = 3143093702.8017533*eccentricity_20 - 1376164819.4890721*eccentricity_18 + 407123249.08456135*eccentricity_16 - 71955647.267619634*eccentricity_14 + 5672187.3666974703*eccentricity_12;
    double common_term_32 = -3114927026.0251179*eccentricity_19 + 906795599.64628042*eccentricity_17 - 158819986.31638283*eccentricity_15 + 12497962.364256332*eccentricity_13;
    double common_term_33 = -6886231975.0692121*eccentricity_20 + 1971248404.6502638*eccentricity_18 - 341626516.23594254*eccentricity_16 + 26775478.224996169*eccentricity_14;
    double common_term_34 = 4192128368.6438094*eccentricity_19 - 718051649.36163461*eccentricity_17 + 55948914.794879542*eccentricity_15;
    double common_term_35 = 8738961908.7768649*eccentricity_20 - 1478077112.9483162*eccentricity_18 + 114323547.98647246*eccentricity_16;
    double common_term_36 = -2985496129.0127606*eccentricity_19 + 228945502.38171663*eccentricity_17;
    double common_term_37 = -5927085183.6951043*eccentricity_20 + 450200722.51609635*eccentricity_18;
    double common_term_38 = 870704611.08444028*eccentricity_19;
    double common_term_39 = 1658617930.9503999*eccentricity_20;
    double common_term_40 = 8.1027892369674124*eccentricity_20;
    double common_term_41 = 5.9528087483326965*eccentricity_19;
    double common_term_42 = 19.591766756782454*eccentricity_20 + 4.3737632425780729*eccentricity_18;
    double common_term_43 = 15.211600409391111*eccentricity_19 + 3.2139496051148163*eccentricity_17;
    double common_term_44 = 37.516397349130769*eccentricity_20 + 11.777771122524456*eccentricity_18 + 2.3619891394852352*eccentricity_16;
    double common_term_45 = 29.872038187499631*eccentricity_19 + 9.0963132818330784*eccentricity_17 + 1.7361061035642716*eccentricity_15;
    double common_term_46 = 62.919119482414333*eccentricity_20 + 23.727646763736514*eccentricity_18 + 7.0095621042565526*eccentricity_16 + 1.2762573453994817*eccentricity_14;
    double common_term_47 = 50.996262906200171*eccentricity_19 + 18.803348544880662*eccentricity_17 + 5.3905542732752184*eccentricity_15 + 0.93835479278592447*eccentricity_13;
    double common_term_48 = 97.174808081350195*eccentricity_20 + 41.239924825749725*eccentricity_18 + 14.867948803167396*eccentricity_16 + 4.1378410456730769*eccentricity_14 + 0.69002422382305195*eccentricity_12;
    double common_term_49 = 79.833392743658076*eccentricity_19 + 33.277260174044665*eccentricity_17 + 11.73128761832799*eccentricity_15 + 3.1708903963106726*eccentricity_13 + 0.50749178561182651*eccentricity_11;
    double common_term_50 = 141.71947253977302*eccentricity_20 + 65.450049277739603*eccentricity_18 + 26.795008788974175*eccentricity_16 + 9.2376283748321509*eccentricity_14 + 2.4261317884950697*eccentricity_12 + 0.37329840443121693*eccentricity_10;
    double common_term_51 = 117.70136781515984*eccentricity_19 + 53.548408720447653*eccentricity_17 + 21.530950284507368*eccentricity_15 + 7.2599628672661719*eccentricity_13 + 1.8536271122523717*eccentricity_11 + 0.27462201799665179*eccentricity_9;
    double common_term_52 = 198.06627759859368*eccentricity_20 + 97.564243743701643*eccentricity_18 + 43.723172022089801*eccentricity_16 + 17.26639119368286*eccentricity_14 + 5.6951335668816138*eccentricity_12 + 1.4143146494708995*eccentricity_10 + 0.20204613095238095*eccentricity_8;
    double common_term_53 = 165.99386411744743*eccentricity_19 + 80.717673397516724*eccentricity_17 + 35.630606926419767*eccentricity_15 + 13.819591798849627*eccentricity_13 + 4.4596811829420625*eccentricity_11 + 1.0777563185918899*eccentricity_9 + 0.14865606398809524*eccentricity_7;
    double common_term_54 = 0.109375*eccentricity_6*std::pow(1.0 - eccentricity_2, -7.5);
    double common_term_55 = 226.18153481709673*eccentricity_19 + 115.95843015327581*eccentricity_17 + 54.939037858091544*eccentricity_15 + 23.526691280302758*eccentricity_13 + 8.803474856179858*eccentricity_11 + 2.7210094633556548*eccentricity_9 + 0.62366536458333333*eccentricity_7 + 0.08046875*eccentricity_5;
    double common_term_56 = 352.61358299181287*eccentricity_20 + 190.70208210702517*eccentricity_18 + 96.660970437168164*eccentricity_16 + 45.201081590332892*eccentricity_14 + 19.065321955605159*eccentricity_12 + 7.0086185515873016*eccentricity_10 + 2.1216145833333333*eccentricity_8 + 0.475*eccentricity_6 + 0.0625*eccentricity_4;
    double common_term_57 = 299.81420056994565*eccentricity_19 + 160.51944121754565*eccentricity_17 + 80.436430373958179*eccentricity_15 + 37.125392362730844*eccentricity_13 + 15.425517817905971*eccentricity_11 + 5.57567138671875*eccentricity_9 + 1.64326171875*eccentricity_7 + 0.43359375*eccentricity_5 - 0.0625*eccentricity_3;
    double common_term_58 = 454.24002084994062*eccentricity_20 + 254.51711451680568*eccentricity_18 + 134.89235375634146*eccentricity_16 + 66.824699332423942*eccentricity_14 + 30.437646639384921*eccentricity_12 + 12.522894965277778*eccentricity_10 + 4.0916666666666667*eccentricity_8 + 2.5*eccentricity_6 - 1.625*eccentricity_4 + 0.75*eccentricity_2;
    double common_term_59 = 388.52071842664467*eccentricity_19 + 215.72270339543276*eccentricity_17 + 113.19083135138609*eccentricity_15 + 55.269214627604005*eccentricity_13 + 25.94196531507704*eccentricity_11 + 5.5236307779947917*eccentricity_9 + 16.171061197916667*eccentricity_7 - 17.3828125*eccentricity_5 + 10.3125*eccentricity_3 - 1.5*eccentricity;
    double common_term_60 = 574.51357144635259*eccentricity_20 + 331.77071928063217*eccentricity_18 + 182.87232126664142*eccentricity_16 + 92.612688137755102*eccentricity_14 + 57.024970703125*eccentricity_12 - 19.72640625*eccentricity_10 + 101.90625*eccentricity_8 - 122.3125*eccentricity_6 + 78.1875*eccentricity_4 - 18.0*eccentricity_2 + 1.0;
    double common_term_61 = 493.43664018153295*eccentricity_19 + 286.7277149523058*eccentricity_17 + 133.14591825483022*eccentricity_15 + 167.7827899007444*eccentricity_13 - 229.92991219414605*eccentricity_11 + 560.69659627278646*eccentricity_9 - 662.57438151041667*eccentricity_7 + 430.1171875*eccentricity_5 - 119.4375*eccentricity_3 + 10.5*eccentricity;
    double common_term_62 = 709.5685140066866*eccentricity_20 + 457.5956742505296*eccentricity_18 + 84.734203898675838*eccentricity_16 + 687.22760082103588*eccentricity_14 - 1391.7037338789683*eccentricity_12 + 2661.0693576388889*eccentricity_10 - 2985.2182291666667*eccentricity_8 + 1915.953125*eccentricity_6 - 580.375*eccentricity_4 + 62.25*eccentricity_2;
    double common_term_63 = 853.62013814635036*eccentricity_19 - 581.81354296275882*eccentricity_17 + 3140.3596028706857*eccentricity_15 - 6674.072648970059*eccentricity_13 + 11096.136155809675*eccentricity_11 - 11708.047570800781*eccentricity_9 + 7329.28740234375*eccentricity_7 - 2307.59765625*eccentricity_5 + 274.1875*eccentricity_3;
    double common_term_64 = 2267.7785591974907*eccentricity_20 - 4355.5196817547584*eccentricity_18 + 13760.787842945097*eccentricity_16 - 27709.08882154569*eccentricity_14 + 41527.90535094246*eccentricity_12 - 41181.945374503968*eccentricity_10 + 24965.729166666667*eccentricity_8 - 7957.5625*eccentricity_6 + 998.75*eccentricity_4;
    double common_term_65 = -21502.829243771265*eccentricity_19 + 55496.840763838927*eccentricity_17 - 103358.24492683926*eccentricity_15 + 142048.9446050667*eccentricity_13 - 132645.21026563493*eccentricity_11 + 77582.160425385975*eccentricity_9 - 24651.222298177083*eccentricity_7 + 3181.34921875*eccentricity_5;
    double common_term_66 = -90285.265770808062*eccentricity_20 + 205830.54338031957*eccentricity_18 - 353853.23002890371*eccentricity_16 + 450549.36138532366*eccentricity_14 - 397243.92134486607*eccentricity_12 + 223753.60797991071*eccentricity_10 - 70208.108035714286*eccentricity_8 + 9163.8375*eccentricity_6;
    double common_term_67 = 708137.22483336192*eccentricity_19 - 1128088.2209649158*eccentricity_17 + 1340315.2745059372*eccentricity_15 - 1118951.471470651*eccentricity_13 + 606550.27029326706*eccentricity_11 - 186830.13835812523*eccentricity_9 + 24400.922632998512*eccentricity_7;
    double common_term_68 = 2281516.9481211117*eccentricity_20 - 3384969.9592361881*eccentricity_18 + 3773709.8652757867*eccentricity_16 - 2991179.6441564297*eccentricity_14 + 1560397.2811733218*eccentricity_12 - 470044.07146990741*eccentricity_10 + 60988.455208333333*eccentricity_8;
    double common_term_69 = -9639760.8470058049*eccentricity_19 + 10129601.942911905*eccentricity_17 - 7642510.1706927773*eccentricity_15 + 3838292.1718349859*eccentricity_13 - 1128041.6384100233*eccentricity_11 + 144690.15621730259*eccentricity_9;
    double common_term_70 = -26227730.516120504*eccentricity_20 + 26076913.958303522*eccentricity_18 - 18770791.61849961*eccentricity_16 + 9081820.6651140959*eccentricity_14 - 2600155.3861269127*eccentricity_12 + 328564.26045841601*eccentricity_10;
    double common_term_71 = 64697083.086877482*eccentricity_19 - 44526994.941952148*eccentricity_17 + 20770338.744247634*eccentricity_15 - 5788030.6720844356*eccentricity_13 + 718795.38038070565*eccentricity_11;
    double common_term_72 = 155327809.93409297*eccentricity_20 - 102412479.29229473*eccentricity_18 + 46097656.206926124*eccentricity_16 - 12497715.918680811*eccentricity_14 + 1522709.4454134537*eccentricity_12;
    double common_term_73 = -229136140.40606587*eccentricity_19 + 99613352.919932352*eccentricity_17 - 26270225.293139121*eccentricity_15 + 3136496.2068759774*eccentricity_13;
    double common_term_74 = -500094015.66830428*eccentricity_20 + 210169901.52995076*eccentricity_18 - 53917426.040027371*eccentricity_16 + 6303051.860868987*eccentricity_14;
    double common_term_75 = 433976913.65316237*eccentricity_19 - 108322148.5424877*eccentricity_17 + 12392155.598679251*eccentricity_15;
    double common_term_76 = 878795499.97859888*eccentricity_20 - 213478009.68628576*eccentricity_18 + 23891754.456744972*eccentricity_16;
    double common_term_77 = -413455221.46545222*eccentricity_19 + 45259959.230096186*eccentricity_17;
    double common_term_78 = -788181684.18023007*eccentricity_20 + 84387825.843883373*eccentricity_18;
    double common_term_79 = 155087839.44832738*eccentricity_19;
    double common_term_80 = 281291802.15753558*eccentricity_20;
    double common_term_81 = 684.54629028972052*eccentricity_20;
    double common_term_82 = 485.49372998260163*eccentricity_19;
    double common_term_83 = 1144.9780867795499*eccentricity_20 + 343.88236671708722*eccentricity_18;
    double common_term_84 = 879.34219777805822*eccentricity_19 + 243.2418036824339*eccentricity_17;
    double common_term_85 = 1892.6731095066038*eccentricity_20 + 670.56700405485279*eccentricity_18 + 171.79859304581788*eccentricity_16;
    double common_term_86 = 1482.6265695480704*eccentricity_19 + 508.10538707517321*eccentricity_17 + 121.14309121499761*eccentricity_15;
    double common_term_87 = 2829.509041078752*eccentricity_20 + 1156.2690220226189*eccentricity_18 + 382.7583192940536*eccentricity_16 + 85.273143039149751*eccentricity_14;
    double common_term_88 = 2254.8466000770337*eccentricity_19 + 897.70267058890213*eccentricity_17 + 286.76737477258496*eccentricity_15 + 59.908468079385938*eccentricity_13;
    double common_term_89 = 3990.275347808563*eccentricity_20 + 1789.6005942428417*eccentricity_18 + 693.80782462991271*eccentricity_16 + 213.74585329958268*eccentricity_14 + 41.99941369350332*eccentricity_12;
    double common_term_90 = 3221.2591713049364*eccentricity_19 + 1414.5253154822474*eccentricity_17 + 533.78907190873898*eccentricity_15 + 158.5322733228732*eccentricity_13 + 29.375192652658448*eccentricity_11;
    double common_term_91 = 5396.3367098348751*eccentricity_20 + 2590.7749708475467*eccentricity_18 + 1113.4240953236217*eccentricity_16 + 408.8027098341112*eccentricity_14 + 117.0156742086039*eccentricity_12 + 20.492165178571429*eccentricity_10;
    double common_term_92 = 4401.9002061049085*eccentricity_19 + 2075.8142892893831*eccentricity_17 + 872.73445321116142*eccentricity_15 + 311.64210511511558*eccentricity_13 + 85.960494009485531*eccentricity_11 + 14.253836199509824*eccentricity_9;
    double common_term_93 = 7070.2707579604941*eccentricity_20 + 3578.2942370745069*eccentricity_18 + 1656.810500052729*eccentricity_16 + 681.15639702873644*eccentricity_14 + 236.46823033785273*eccentricity_12 + 62.845375881834215*eccentricity_10 + 9.8822296626984127*eccentricity_8;
    double common_term_94 = 5817.5666102846243*eccentricity_19 + 2898.5180971909654*eccentricity_17 + 1317.1940303954908*eccentricity_15 + 529.32218106269836*eccentricity_13 + 178.57800832475935*eccentricity_11 + 45.722182355608259*eccentricity_9 + 6.8260463169642857*eccentricity_7;
    double common_term_95 = 9035.210792432971*eccentricity_20 + 4771.2620644388882*eccentricity_18 + 2339.410327046298*eccentricity_16 + 1042.9932195675338*eccentricity_14 + 409.50493126745297*eccentricity_12 + 134.20566716269841*eccentricity_10 + 33.096825396825397*eccentricity_8 + 4.6951388888888889*eccentricity_6;
    double common_term_96 = 7489.6065989639505*eccentricity_19 + 3900.1292586073831*eccentricity_17 + 1881.1768521114715*eccentricity_15 + 822.4778344982025*eccentricity_13 + 315.36323912080634*eccentricity_11 + 100.35344373914931*eccentricity_9 + 23.830978732638889*eccentricity_7 + 3.2138020833333333*eccentricity_5;
    double common_term_97 = std::pow(1.0 - eccentricity_2, -7.5)*(0.65625*eccentricity_6 + 2.1875*eccentricity_4);
    double common_term_98 = 9439.9185340714445*eccentricity_19 + 5098.6663329202123*eccentricity_17 + 2579.1932129021996*eccentricity_15 + 1202.4828373346127*eccentricity_13 + 504.9256145134174*eccentricity_11 + 184.36379711009838*eccentricity_9 + 55.222623697916667*eccentricity_7 + 12.141927083333333*eccentricity_5 + 1.4791666666666667*eccentricity_3;
    double common_term_99 = 13933.457732411072*eccentricity_20 + 7852.653663726652*eccentricity_18 + 4186.6089429039684*eccentricity_16 + 2086.2112137965719*eccentricity_14 + 955.65467282572751*eccentricity_12 + 392.9662037037037*eccentricity_10 + 139.89722222222222*eccentricity_8 + 40.614583333333333*eccentricity_6 + 8.5833333333333333*eccentricity_4 + eccentricity_2;
    double common_term_100 = 11690.9490795956*eccentricity_19 + 6512.6725946163192*eccentricity_17 + 3426.2534545571828*eccentricity_15 + 1681.1850924055917*eccentricity_13 + 756.32119033813477*eccentricity_11 + 304.38755493164062*eccentricity_9 + 105.59130859375*eccentricity_7 + 29.6484375*eccentricity_5 + 6.1875*eccentricity_3 + 0.5*eccentricity;
    double common_term_101 = 16915.840809042576*eccentricity_20 + 9781.9766586569588*eccentricity_18 + 5384.6709878829914*eccentricity_16 + 2794.3606267026487*eccentricity_14 + 1349.584135320216*eccentricity_12 + 595.95076388888889*eccentricity_10 + 234.75564236111111*eccentricity_8 + 78.555555555555556*eccentricity_6 + 23.5*eccentricity_4 + 2.0*eccentricity_2 + 1.0;
    double common_term_102 = 14265.698283182097*eccentricity_19 + 8161.2224713156346*eccentricity_17 + 4437.8718071490992*eccentricity_15 + 2270.9555916096588*eccentricity_13 + 1078.7587194089536*eccentricity_11 + 469.22595553927951*eccentricity_9 + 173.39816623263889*eccentricity_7 + 73.815104166666667*eccentricity_5 - 0.3125*eccentricity_3 + 8.5*eccentricity;
    double common_term_103 = 20287.386428669096*eccentricity_20 + 11998.559434286312*eccentricity_18 + 6788.8154052734375*eccentricity_16 + 3646.0642445591518*eccentricity_14 + 1835.5665234375*eccentricity_12 + 873.31640625*eccentricity_10 + 321.95625*eccentricity_8 + 218.90625*eccentricity_6 - 37.25*eccentricity_4 + 42.0*eccentricity_2;
    double common_term_104 = 17187.812675214463*eccentricity_19 + 10063.167107271199*eccentricity_17 + 5634.8871911270048*eccentricity_15 + 2960.8069926811148*eccentricity_13 + 1569.906435168231*eccentricity_11 + 458.83704653139468*eccentricity_9 + 648.71975911458333*eccentricity_7 - 226.96223958333333*eccentricity_5 + 157.64583333333333*eccentricity_3;
    double common_term_105 = 24075.003333778354*eccentricity_20 + 14517.835033201334*eccentricity_18 + 8450.9566388897022*eccentricity_16 + 4514.8397351799034*eccentricity_14 + 2868.7049109726356*eccentricity_12 + 232.20453042328042*eccentricity_10 + 1927.3425347222222*eccentricity_8 - 932.79583333333333*eccentricity_6 + 497.97916666666667*eccentricity_4;
    double common_term_106 = 20438.506083660771*eccentricity_19 + 12433.812172157926*eccentricity_17 + 6332.2276472097976*eccentricity_15 + 5690.9572262048721*eccentricity_13 - 1628.573877116612*eccentricity_11 + 5629.1753435407366*eccentricity_9 - 3134.50556640625*eccentricity_7 + 1395.14296875*eccentricity_5;
    double common_term_107 = 28063.145189358618*eccentricity_20 + 18323.69426441576*eccentricity_18 + 7347.0933989638552*eccentricity_16 + 12851.250283826035*eccentricity_14 - 9021.6337471524104*eccentricity_12 + 15875.237475198413*eccentricity_10 - 9247.1106150793651*eccentricity_8 + 3575.9618055555556*eccentricity_6;
    double common_term_108 = 28223.992889681186*eccentricity_19 + 3467.047735114143*eccentricity_17 + 32633.328711873985*eccentricity_15 - 32777.046731996615*eccentricity_13 + 42844.59543250015*eccentricity_11 - 24842.062921966068*eccentricity_9 + 8555.9166837177579*eccentricity_7;
    double common_term_109 = 48702.277797475021*eccentricity_20 - 18597.077919061895*eccentricity_18 + 88137.731547360111*eccentricity_16 - 100913.54372666396*eccentricity_14 + 110448.49248883929*eccentricity_12 - 62129.775892857143*eccentricity_10 + 19375.531473214286*eccentricity_8;
    double common_term_110 = -98098.07201527014*eccentricity_19 + 240434.83814053971*eccentricity_17 - 281991.78143195508*eccentricity_15 + 272490.68950632819*eccentricity_13 - 146790.81474965427*eccentricity_11 + 41945.743762018651*eccentricity_9;
    double common_term_111 = -343403.35026978921*eccentricity_20 + 643373.0744003425*eccentricity_18 - 736303.14226992921*eccentricity_16 + 645651.00100172004*eccentricity_14 - 331033.60035002806*eccentricity_12 + 87458.838555169753*eccentricity_10;
    double common_term_112 = 1668417.5994084481*eccentricity_19 - 1824503.4811775299*eccentricity_17 + 1475032.0988348779*eccentricity_15 - 718011.57625047417*eccentricity_13 + 176637.60430073428*eccentricity_11;
    double common_term_113 = 4179720.2338346356*eccentricity_20 - 4331823.4154171438*eccentricity_18 + 3261446.0175656087*eccentricity_16 - 1506624.8069187101*eccentricity_14 + 347116.40083112144*eccentricity_12;
    double common_term_114 = -9919713.8674406614*eccentricity_19 + 7003657.568782102*eccentricity_17 - 3072374.4261272348*eccentricity_15 + 666099.66287382847*eccentricity_13;
    double common_term_115 = -22016384.411622141*eccentricity_20 + 14651272.502469821*eccentricity_18 - 6111166.6236146133*eccentricity_16 + 1251825.6868847893*eccentricity_14;
    double common_term_116 = 29938134.655573779*eccentricity_19 - 11891770.972113935*eccentricity_17 + 2309591.5418258783*eccentricity_15;
    double common_term_117 = 59894316.096227556*eccentricity_20 - 22693887.384820599*eccentricity_18 + 4191663.6764965681*eccentricity_16;
    double common_term_118 = -42560311.713931959*eccentricity_19 + 7496082.4211335116*eccentricity_17;
    double common_term_119 = -78576101.974312903*eccentricity_20 + 13228372.610073542*eccentricity_18;
    double common_term_120 = 23064366.879176214*eccentricity_19;
    double common_term_121 = 39774622.640618849*eccentricity_20;
    double common_term_122 = 21457.541315085186*eccentricity_20;
    double common_term_123 = 14520.426142632637*eccentricity_19;
    double common_term_124 = 16433.491927399273*eccentricity_20 + 9790.8208210865492*eccentricity_18;
    double common_term_125 = 13059.732286428915*eccentricity_19 + 6575.8288627245346*eccentricity_17;
    double common_term_126 = 25449.720230708691*eccentricity_20 + 10109.125707205264*eccentricity_18 + 4397.4592565370715*eccentricity_16;
    double common_term_127 = 19450.263484654546*eccentricity_19 + 7661.4517293291343*eccentricity_17 + 2926.6658339357941*eccentricity_15;
    double common_term_128 = 31982.138955127009*eccentricity_20 + 14782.127106249834*eccentricity_18 + 5703.6922359837428*eccentricity_16 + 1937.4461366158953*eccentricity_14;
    double common_term_129 = 24883.813909130604*eccentricity_19 + 11156.013360185008*eccentricity_17 + 4179.9295842065951*eccentricity_15 + 1274.9490624733074*eccentricity_13;
    double common_term_130 = 39377.235774997017*eccentricity_20 + 19213.327264401882*eccentricity_18 + 8351.639543728627*eccentricity_16 + 3019.4252945214588*eccentricity_14 + 833.35244924025306*eccentricity_12;
    double common_term_131 = 30906.554058116621*eccentricity_19 + 14716.750973884218*eccentricity_17 + 6196.3667029721152*eccentricity_15 + 2151.49817300716*eccentricity_13 + 540.54460558705516*eccentricity_11;
    double common_term_132 = 47260.616346966225*eccentricity_20 + 24083.818416646534*eccentricity_18 + 11177.520827804715*eccentricity_16 + 4552.5318770684952*eccentricity_14 + 1512.621510137962*eccentricity_12 + 347.53703627921076*eccentricity_10;
    double common_term_133 = 37352.689398178225*eccentricity_19 + 18623.497918340524*eccentricity_17 + 8413.136035539822*eccentricity_15 + 3309.5092659818069*eccentricity_13 + 1049.0975573042881*eccentricity_11 + 221.16175382637683*eccentricity_9;
    double common_term_134 = 55607.01013687221*eccentricity_20 + 29313.684256228025*eccentricity_18 + 14282.906556364778*eccentricity_16 + 6271.2194009993912*eccentricity_14 + 2378.333583984375*eccentricity_12 + 717.35178571428571*eccentricity_10 + 139.04497767857143*eccentricity_8;
    double common_term_135 = 44189.475958295735*eccentricity_19 + 22830.633653517989*eccentricity_17 + 10856.895523331472*eccentricity_15 + 4625.5713758638226*eccentricity_13 + 1687.7395159773515*eccentricity_11 + 483.05924309624566*eccentricity_9 + 86.155606708829365*eccentricity_7;
    double common_term_136 = 64373.015351742061*eccentricity_20 + 34869.628497130946*eccentricity_18 + 17635.868583446799*eccentricity_16 + 8173.0600789632018*eccentricity_14 + 3372.4788363807136*eccentricity_12 + 1181.0311042906746*eccentricity_10 + 319.79875992063492*eccentricity_8 + 52.440972222222222*eccentricity_6;
    double common_term_137 = 51376.607452654685*eccentricity_19 + 27306.48426836369*eccentricity_17 + 13501.639044812935*eccentricity_15 + 6087.4061027492796*eccentricity_13 + 2427.3691937582833*eccentricity_11 + 813.49096854073661*eccentricity_9 + 207.61494140625*eccentricity_7 + 31.21171875*eccentricity_5;
    double common_term_138 = 73512.349191304375*eccentricity_20 + 40714.145849651768*eccentricity_18 + 21206.903160576236*eccentricity_16 + 10235.202401781121*eccentricity_14 + 4480.4448167266038*eccentricity_12 + 1721.8277571097884*eccentricity_10 + 550.18836805555556*eccentricity_8 + 131.67916666666667*eccentricity_6 + 18.041666666666667*eccentricity_4;
    double common_term_139 = 58870.770642880133*eccentricity_19 + 32016.166273716718*eccentricity_17 + 16320.07322272787*eccentricity_15 + 7674.3587744977487*eccentricity_13 + 3253.691183524157*eccentricity_11 + 1200.9678082501447*eccentricity_9 + 364.11526692708333*eccentricity_7 + 81.131510416666667*eccentricity_5 + 10.020833333333333*eccentricity_3;
    double common_term_140 = std::pow(1.0 - eccentricity_2, -7.5)*(1.640625*eccentricity_6 + 8.75*eccentricity_4 + 5.25*eccentricity_2);
    double common_term_141 = 66625.583787931284*eccentricity_19 + 36921.869880090993*eccentricity_17 + 19282.115231111371*eccentricity_15 + 9363.2274901158637*eccentricity_13 + 4149.5350677772805*eccentricity_11 + 1633.5173536512587*eccentricity_9 + 547.90695529513889*eccentricity_7 + 146.03385416666667*eccentricity_5 + 27.0625*eccentricity_3 + 2.5*eccentricity;
    double common_term_142 = 92710.20293008311*eccentricity_20 + 53103.938964984605*eccentricity_18 + 28870.682726496405*eccentricity_16 + 14736.952994648959*eccentricity_14 + 6960.0260164689429*eccentricity_12 + 2979.5158159722222*eccentricity_10 + 1121.8220486111111*eccentricity_8 + 354.40972222222222*eccentricity_6 + 86.6875*eccentricity_4 + 14.0*eccentricity_2 + 1.0;
    double common_term_143 = 74591.591253809641*eccentricity_19 + 41982.85119931195*eccentricity_17 + 22354.895010375599*eccentricity_15 + 11128.178120109013*eccentricity_13 + 5095.2647703552246*eccentricity_11 + 2096.8638977050781*eccentricity_9 + 749.23486328125*eccentricity_7 + 219.8671875*eccentricity_5 + 47.8125*eccentricity_3 + 6.5*eccentricity;
    double common_term_144 = 102660.47235567459*eccentricity_20 + 59559.225825986113*eccentricity_18 + 32890.024577477676*eccentricity_16 + 17118.068152247299*eccentricity_14 + 8286.3820745494378*eccentricity_12 + 3662.7114945023148*eccentricity_10 + 1439.2086805555556*eccentricity_8 + 482.55208333333333*eccentricity_6 + 127.20833333333333*eccentricity_4 + 25.75*eccentricity_2;
    double common_term_145 = 82716.266730041409*eccentricity_19 + 47155.435965436727*eccentricity_17 + 25502.760370540711*eccentricity_15 + 12940.742162879182*eccentricity_13 + 6068.8542543643366*eccentricity_11 + 2573.966889558015*eccentricity_9 + 958.58414713541667*eccentricity_7 + 287.45442708333333*eccentricity_5 + 80.104166666666667*eccentricity_3;
    double common_term_146 = 112767.57294417368*eccentricity_20 + 66123.111110768957*eccentricity_18 + 36980.648179234096*eccentricity_16 + 19542.408988560268*eccentricity_14 + 9636.7225655691964*eccentricity_12 + 4352.0655133928571*eccentricity_10 + 1772.05546875*eccentricity_8 + 573.0*eccentricity_6 + 214.9375*eccentricity_4;
    double common_term_147 = 90944.004022405646*eccentricity_19 + 52393.141663500394*eccentricity_17 + 28686.436395698423*eccentricity_15 + 14774.660635731053*eccentricity_13 + 7024.8494779092294*eccentricity_11 + 3111.831438530816*eccentricity_9 + 1021.5815646701389*eccentricity_7 + 521.30130208333333*eccentricity_5;
    double common_term_148 = 122969.4109833409*eccentricity_20 + 72744.048051908649*eccentricity_18 + 41093.082797605334*eccentricity_16 + 22000.731539307358*eccentricity_14 + 10884.473990081662*eccentricity_12 + 5286.5708922371032*eccentricity_10 + 1618.5184771825397*eccentricity_8 + 1174.3086805555556*eccentricity_6;
    double common_term_149 = 99222.195850010797*eccentricity_19 + 57615.763986263487*eccentricity_17 + 31989.224256097365*eccentricity_15 + 16204.755631538119*eccentricity_13 + 8851.5450169699533*eccentricity_11 + 2188.9004241943359*eccentricity_9 + 2499.8066545758929*eccentricity_7;
    double common_term_150 = 133233.37892880022*eccentricity_20 + 79223.710053646716*eccentricity_18 + 45685.97610933302*eccentricity_16 + 23058.242806039166*eccentricity_14 + 14894.877307925485*eccentricity_12 + 2167.5672054122575*eccentricity_10 + 5088.2304067460317*eccentricity_8;
    double common_term_151 = 106902.16979188773*eccentricity_19 + 64610.243249226481*eccentricity_17 + 30834.009431963243*eccentricity_15 + 25660.809450528702*eccentricity_13 + 138.15669965323527*eccentricity_11 + 9986.2402226469925*eccentricity_9;
    double common_term_152 = 141312.49332766574*eccentricity_20 + 91630.502552760306*eccentricity_18 + 37074.123898830662*eccentricity_16 + 45857.294896129261*eccentricity_14 - 7051.1112312297078*eccentricity_12 + 19015.278426339286*eccentricity_10;
    double common_term_153 = 132879.10032182419*eccentricity_19 + 34830.4032642377*eccentricity_17 + 85307.598335647122*eccentricity_15 - 25928.956411922032*eccentricity_13 + 35295.352656819605*eccentricity_11;
    double common_term_154 = 202278.12442513929*eccentricity_20 + 6928.1237253255422*eccentricity_18 + 164169.08244098456*eccentricity_16 - 69371.145960718407*eccentricity_14 + 64097.985970790077*eccentricity_12;
    double common_term_155 = -86032.923749765168*eccentricity_19 + 322941.51944711517*eccentricity_17 - 161863.38631548219*eccentricity_15 + 114223.29016644697*eccentricity_13;
    double common_term_156 = -330019.70183264538*eccentricity_20 + 641232.70388495405*eccentricity_18 - 348665.8130864197*eccentricity_16 + 200206.09740457872*eccentricity_14;
    double common_term_157 = 1272258.805030937*eccentricity_19 - 711536.5464652625*eccentricity_17 + 345825.64133652176*eccentricity_15;
    double common_term_158 = 2505187.0198386038*eccentricity_20 - 1395364.84038749*eccentricity_18 + 589653.60386437926*eccentricity_16;
    double common_term_159 = -2652779.1382360331*eccentricity_19 + 993773.49077001778*eccentricity_17;
    double common_term_160 = -4918121.3688505491*eccentricity_20 + 1657411.0514342354*eccentricity_18;
    double common_term_161 = 2738137.0118639785*eccentricity_19;
    double common_term_162 = 4484698.6790054913*eccentricity_20;
    double common_term_163 = 377787.54726339713*eccentricity_20;
    double common_term_164 = 242955.63362904255*eccentricity_19;
    double common_term_165 = -66096.229281512686*eccentricity_20 + 155290.31789344062*eccentricity_18;
    double common_term_166 = -11144.000417105628*eccentricity_19 + 98590.791164105207*eccentricity_17;
    double common_term_167 = 159817.14962769454*eccentricity_20 + 12788.373692859968*eccentricity_18 + 62129.432838845761*eccentricity_16;
    double common_term_168 = 108624.48124853734*eccentricity_19 + 20659.600133680822*eccentricity_17 + 38830.271627030094*eccentricity_15;
    double common_term_169 = 121745.7485319993*eccentricity_20 + 75593.797570510655*eccentricity_18 + 20846.48947785767*eccentricity_16 + 24045.531161456816*eccentricity_14;
    double common_term_170 = 94725.791656940672*eccentricity_19 + 53394.726806174315*eccentricity_17 + 17864.457369950933*eccentricity_15 + 14736.269914774821*eccentricity_13;
    double common_term_171 = 134426.43187566241*eccentricity_20 + 72227.962053432338*eccentricity_18 + 37948.391661483745*eccentricity_16 + 14014.760135060252*eccentricity_14 + 8925.3544891436688*eccentricity_12;
    double common_term_172 = 102791.78618185255*eccentricity_19 + 54152.713887361873*eccentricity_17 + 26932.785682546405*eccentricity_15 + 10366.636800894334*eccentricity_13 + 5333.4715783509754*eccentricity_11;
    double common_term_173 = 140957.5406282885*eccentricity_20 + 77782.844443627939*eccentricity_18 + 39972.847645232492*eccentricity_16 + 18970.483798292636*eccentricity_14 + 7331.8707470914502*eccentricity_12 + 3137.8065166170635*eccentricity_10;
    double common_term_174 = 108013.04578354859*eccentricity_19 + 58165.735450234899*eccentricity_17 + 29049.230429099558*eccentricity_15 + 13195.548930184624*eccentricity_13 + 4992.8608943285261*eccentricity_11 + 1812.6737448556083*eccentricity_9;
    double common_term_175 = 146512.9742400296*eccentricity_20 + 81837.843034865561*eccentricity_18 + 42925.4009106688*eccentricity_16 + 20765.734099138709*eccentricity_14 + 9026.688396577381*eccentricity_12 + 3283.9923115079365*eccentricity_10 + 1024.7152901785714*eccentricity_8;
    double common_term_176 = 112268.54521804222*eccentricity_19 + 61238.870862752279*eccentricity_17 + 31215.17029686814*eccentricity_15 + 14579.106566331298*eccentricity_13 + 6049.8636864980062*eccentricity_11 + 2087.4189422607422*eccentricity_9 + 564.30244140625*eccentricity_7;
    double common_term_177 = 150855.0546051197*eccentricity_20 + 85026.680068939605*eccentricity_18 + 45195.984850826907*eccentricity_16 + 22327.635171595982*eccentricity_14 + 10030.963895089286*eccentricity_12 + 3957.3959263392857*eccentricity_10 + 1280.1776785714286*eccentricity_8 + 300.8625*eccentricity_6;
    double common_term_178 = 115477.11687331004*eccentricity_19 + 63563.068261970348*eccentricity_17 + 32843.366350184594*eccentricity_15 + 15674.387087534975*eccentricity_13 + 6744.2472852979388*eccentricity_11 + 2515.3772408621652*eccentricity_9 + 754.56982421875*eccentricity_7 + 153.94921875*eccentricity_5;
    double common_term_179 = 153969.26674670476*eccentricity_20 + 87316.524018284058*eccentricity_18 + 46829.94445196587*eccentricity_16 + 23451.645246155754*eccentricity_14 + 10769.689983258929*eccentricity_12 + 4414.2023809523809*eccentricity_10 + 1544.878125*eccentricity_8 + 424.575*eccentricity_6 + 74.625*eccentricity_4;
    double common_term_180 = 117624.4669626307*eccentricity_19 + 65120.458644287094*eccentricity_17 + 33936.454313233495*eccentricity_15 + 16411.23981387479*eccentricity_13 + 7216.2684065682547*eccentricity_11 + 2798.1171752929687*eccentricity_9 + 909.78837890625*eccentricity_7 + 225.52734375*eccentricity_5 + 33.5625*eccentricity_3;
    double common_term_181 = 155842.88171796954*eccentricity_20 + 88695.303067545112*eccentricity_18 + 47814.966918611918*eccentricity_16 + 24130.470153924851*eccentricity_14 + 11217.020379464286*eccentricity_12 + 4692.5872395833333*eccentricity_10 + 1705.44375*eccentricity_8 + 507.9375*eccentricity_6 + 111.0*eccentricity_4 + 13.5*eccentricity_2;
    double common_term_182 = 118700.63286320097*eccentricity_19 + 65901.534789433163*eccentricity_17 + 34485.263931086781*eccentricity_15 + 16781.810886405158*eccentricity_13 + 7454.2881489562988*eccentricity_11 + 2941.4128967285156*eccentricity_9 + 988.83935546875*eccentricity_7 + 264.0234375*eccentricity_5 + 48.9375*eccentricity_3 + 4.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(41);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (8, 0).
    // q = -20
    result_by_lpq.set(c_Key3(8, 0, -20), common_term_0);
    result_by_q.set(c_Key1(-20), common_term_0);
    // q = -19
    result_by_lpq.set(c_Key3(8, 0, -19), common_term_1);
    result_by_q.set(c_Key1(-19), common_term_1);
    // q = -18
    result_by_lpq.set(c_Key3(8, 0, -18), common_term_2);
    result_by_q.set(c_Key1(-18), common_term_2);
    // q = -17
    result_by_lpq.set(c_Key3(8, 0, -17), common_term_3);
    result_by_q.set(c_Key1(-17), common_term_3);
    // q = -16
    result_by_lpq.set(c_Key3(8, 0, -16), common_term_4);
    result_by_q.set(c_Key1(-16), common_term_4);
    // q = -15
    result_by_lpq.set(c_Key3(8, 0, -15), common_term_5);
    result_by_q.set(c_Key1(-15), common_term_5);
    // q = -14
    result_by_lpq.set(c_Key3(8, 0, -14), common_term_6);
    result_by_q.set(c_Key1(-14), common_term_6);
    // q = -13
    result_by_lpq.set(c_Key3(8, 0, -13), common_term_7);
    result_by_q.set(c_Key1(-13), common_term_7);
    // q = -12
    result_by_lpq.set(c_Key3(8, 0, -12), common_term_8);
    result_by_q.set(c_Key1(-12), common_term_8);
    // q = -11
    result_by_lpq.set(c_Key3(8, 0, -11), common_term_9);
    result_by_q.set(c_Key1(-11), common_term_9);
    // q = -10
    result_by_lpq.set(c_Key3(8, 0, -10), common_term_10);
    result_by_q.set(c_Key1(-10), common_term_10);
    // q = -9
    result_by_lpq.set(c_Key3(8, 0, -9), common_term_11);
    result_by_q.set(c_Key1(-9), common_term_11);
    // q = -7
    result_by_lpq.set(c_Key3(8, 0, -7), common_term_12);
    result_by_q.set(c_Key1(-7), common_term_12);
    // q = -6
    result_by_lpq.set(c_Key3(8, 0, -6), common_term_13);
    result_by_q.set(c_Key1(-6), common_term_13);
    // q = -5
    result_by_lpq.set(c_Key3(8, 0, -5), common_term_14);
    result_by_q.set(c_Key1(-5), common_term_14);
    // q = -4
    result_by_lpq.set(c_Key3(8, 0, -4), common_term_15);
    result_by_q.set(c_Key1(-4), common_term_15);
    // q = -3
    result_by_lpq.set(c_Key3(8, 0, -3), common_term_16);
    result_by_q.set(c_Key1(-3), common_term_16);
    // q = -2
    result_by_lpq.set(c_Key3(8, 0, -2), common_term_17);
    result_by_q.set(c_Key1(-2), common_term_17);
    // q = -1
    result_by_lpq.set(c_Key3(8, 0, -1), common_term_18);
    result_by_q.set(c_Key1(-1), common_term_18);
    // q = 0
    result_by_lpq.set(c_Key3(8, 0, 0), common_term_19);
    result_by_q.set(c_Key1(0), common_term_19);
    // q = 1
    result_by_lpq.set(c_Key3(8, 0, 1), common_term_20);
    result_by_q.set(c_Key1(1), common_term_20);
    // q = 2
    result_by_lpq.set(c_Key3(8, 0, 2), common_term_21);
    result_by_q.set(c_Key1(2), common_term_21);
    // q = 3
    result_by_lpq.set(c_Key3(8, 0, 3), common_term_22);
    result_by_q.set(c_Key1(3), common_term_22);
    // q = 4
    result_by_lpq.set(c_Key3(8, 0, 4), common_term_23);
    result_by_q.set(c_Key1(4), common_term_23);
    // q = 5
    result_by_lpq.set(c_Key3(8, 0, 5), common_term_24);
    result_by_q.set(c_Key1(5), common_term_24);
    // q = 6
    result_by_lpq.set(c_Key3(8, 0, 6), common_term_25);
    result_by_q.set(c_Key1(6), common_term_25);
    // q = 7
    result_by_lpq.set(c_Key3(8, 0, 7), common_term_26);
    result_by_q.set(c_Key1(7), common_term_26);
    // q = 8
    result_by_lpq.set(c_Key3(8, 0, 8), common_term_27);
    result_by_q.set(c_Key1(8), common_term_27);
    // q = 9
    result_by_lpq.set(c_Key3(8, 0, 9), common_term_28);
    result_by_q.set(c_Key1(9), common_term_28);
    // q = 10
    result_by_lpq.set(c_Key3(8, 0, 10), common_term_29);
    result_by_q.set(c_Key1(10), common_term_29);
    // q = 11
    result_by_lpq.set(c_Key3(8, 0, 11), common_term_30);
    result_by_q.set(c_Key1(11), common_term_30);
    // q = 12
    result_by_lpq.set(c_Key3(8, 0, 12), common_term_31);
    result_by_q.set(c_Key1(12), common_term_31);
    // q = 13
    result_by_lpq.set(c_Key3(8, 0, 13), common_term_32);
    result_by_q.set(c_Key1(13), common_term_32);
    // q = 14
    result_by_lpq.set(c_Key3(8, 0, 14), common_term_33);
    result_by_q.set(c_Key1(14), common_term_33);
    // q = 15
    result_by_lpq.set(c_Key3(8, 0, 15), common_term_34);
    result_by_q.set(c_Key1(15), common_term_34);
    // q = 16
    result_by_lpq.set(c_Key3(8, 0, 16), common_term_35);
    result_by_q.set(c_Key1(16), common_term_35);
    // q = 17
    result_by_lpq.set(c_Key3(8, 0, 17), common_term_36);
    result_by_q.set(c_Key1(17), common_term_36);
    // q = 18
    result_by_lpq.set(c_Key3(8, 0, 18), common_term_37);
    result_by_q.set(c_Key1(18), common_term_37);
    // q = 19
    result_by_lpq.set(c_Key3(8, 0, 19), common_term_38);
    result_by_q.set(c_Key1(19), common_term_38);
    // q = 20
    result_by_lpq.set(c_Key3(8, 0, 20), common_term_39);
    result_by_q.set(c_Key1(20), common_term_39);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 0), result_by_q);
    result_by_q.clear();

    // l , p = (8, 1).
    // q = -20
    result_by_lpq.set(c_Key3(8, 1, -20), common_term_40);
    result_by_q.set(c_Key1(-20), common_term_40);
    // q = -19
    result_by_lpq.set(c_Key3(8, 1, -19), common_term_41);
    result_by_q.set(c_Key1(-19), common_term_41);
    // q = -18
    result_by_lpq.set(c_Key3(8, 1, -18), common_term_42);
    result_by_q.set(c_Key1(-18), common_term_42);
    // q = -17
    result_by_lpq.set(c_Key3(8, 1, -17), common_term_43);
    result_by_q.set(c_Key1(-17), common_term_43);
    // q = -16
    result_by_lpq.set(c_Key3(8, 1, -16), common_term_44);
    result_by_q.set(c_Key1(-16), common_term_44);
    // q = -15
    result_by_lpq.set(c_Key3(8, 1, -15), common_term_45);
    result_by_q.set(c_Key1(-15), common_term_45);
    // q = -14
    result_by_lpq.set(c_Key3(8, 1, -14), common_term_46);
    result_by_q.set(c_Key1(-14), common_term_46);
    // q = -13
    result_by_lpq.set(c_Key3(8, 1, -13), common_term_47);
    result_by_q.set(c_Key1(-13), common_term_47);
    // q = -12
    result_by_lpq.set(c_Key3(8, 1, -12), common_term_48);
    result_by_q.set(c_Key1(-12), common_term_48);
    // q = -11
    result_by_lpq.set(c_Key3(8, 1, -11), common_term_49);
    result_by_q.set(c_Key1(-11), common_term_49);
    // q = -10
    result_by_lpq.set(c_Key3(8, 1, -10), common_term_50);
    result_by_q.set(c_Key1(-10), common_term_50);
    // q = -9
    result_by_lpq.set(c_Key3(8, 1, -9), common_term_51);
    result_by_q.set(c_Key1(-9), common_term_51);
    // q = -8
    result_by_lpq.set(c_Key3(8, 1, -8), common_term_52);
    result_by_q.set(c_Key1(-8), common_term_52);
    // q = -7
    result_by_lpq.set(c_Key3(8, 1, -7), common_term_53);
    result_by_q.set(c_Key1(-7), common_term_53);
    // q = -6
    result_by_lpq.set(c_Key3(8, 1, -6), common_term_54);
    result_by_q.set(c_Key1(-6), common_term_54);
    // q = -5
    result_by_lpq.set(c_Key3(8, 1, -5), common_term_55);
    result_by_q.set(c_Key1(-5), common_term_55);
    // q = -4
    result_by_lpq.set(c_Key3(8, 1, -4), common_term_56);
    result_by_q.set(c_Key1(-4), common_term_56);
    // q = -3
    result_by_lpq.set(c_Key3(8, 1, -3), common_term_57);
    result_by_q.set(c_Key1(-3), common_term_57);
    // q = -2
    result_by_lpq.set(c_Key3(8, 1, -2), common_term_58);
    result_by_q.set(c_Key1(-2), common_term_58);
    // q = -1
    result_by_lpq.set(c_Key3(8, 1, -1), common_term_59);
    result_by_q.set(c_Key1(-1), common_term_59);
    // q = 0
    result_by_lpq.set(c_Key3(8, 1, 0), common_term_60);
    result_by_q.set(c_Key1(0), common_term_60);
    // q = 1
    result_by_lpq.set(c_Key3(8, 1, 1), common_term_61);
    result_by_q.set(c_Key1(1), common_term_61);
    // q = 2
    result_by_lpq.set(c_Key3(8, 1, 2), common_term_62);
    result_by_q.set(c_Key1(2), common_term_62);
    // q = 3
    result_by_lpq.set(c_Key3(8, 1, 3), common_term_63);
    result_by_q.set(c_Key1(3), common_term_63);
    // q = 4
    result_by_lpq.set(c_Key3(8, 1, 4), common_term_64);
    result_by_q.set(c_Key1(4), common_term_64);
    // q = 5
    result_by_lpq.set(c_Key3(8, 1, 5), common_term_65);
    result_by_q.set(c_Key1(5), common_term_65);
    // q = 6
    result_by_lpq.set(c_Key3(8, 1, 6), common_term_66);
    result_by_q.set(c_Key1(6), common_term_66);
    // q = 7
    result_by_lpq.set(c_Key3(8, 1, 7), common_term_67);
    result_by_q.set(c_Key1(7), common_term_67);
    // q = 8
    result_by_lpq.set(c_Key3(8, 1, 8), common_term_68);
    result_by_q.set(c_Key1(8), common_term_68);
    // q = 9
    result_by_lpq.set(c_Key3(8, 1, 9), common_term_69);
    result_by_q.set(c_Key1(9), common_term_69);
    // q = 10
    result_by_lpq.set(c_Key3(8, 1, 10), common_term_70);
    result_by_q.set(c_Key1(10), common_term_70);
    // q = 11
    result_by_lpq.set(c_Key3(8, 1, 11), common_term_71);
    result_by_q.set(c_Key1(11), common_term_71);
    // q = 12
    result_by_lpq.set(c_Key3(8, 1, 12), common_term_72);
    result_by_q.set(c_Key1(12), common_term_72);
    // q = 13
    result_by_lpq.set(c_Key3(8, 1, 13), common_term_73);
    result_by_q.set(c_Key1(13), common_term_73);
    // q = 14
    result_by_lpq.set(c_Key3(8, 1, 14), common_term_74);
    result_by_q.set(c_Key1(14), common_term_74);
    // q = 15
    result_by_lpq.set(c_Key3(8, 1, 15), common_term_75);
    result_by_q.set(c_Key1(15), common_term_75);
    // q = 16
    result_by_lpq.set(c_Key3(8, 1, 16), common_term_76);
    result_by_q.set(c_Key1(16), common_term_76);
    // q = 17
    result_by_lpq.set(c_Key3(8, 1, 17), common_term_77);
    result_by_q.set(c_Key1(17), common_term_77);
    // q = 18
    result_by_lpq.set(c_Key3(8, 1, 18), common_term_78);
    result_by_q.set(c_Key1(18), common_term_78);
    // q = 19
    result_by_lpq.set(c_Key3(8, 1, 19), common_term_79);
    result_by_q.set(c_Key1(19), common_term_79);
    // q = 20
    result_by_lpq.set(c_Key3(8, 1, 20), common_term_80);
    result_by_q.set(c_Key1(20), common_term_80);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 1), result_by_q);
    result_by_q.clear();

    // l , p = (8, 2).
    // q = -20
    result_by_lpq.set(c_Key3(8, 2, -20), common_term_81);
    result_by_q.set(c_Key1(-20), common_term_81);
    // q = -19
    result_by_lpq.set(c_Key3(8, 2, -19), common_term_82);
    result_by_q.set(c_Key1(-19), common_term_82);
    // q = -18
    result_by_lpq.set(c_Key3(8, 2, -18), common_term_83);
    result_by_q.set(c_Key1(-18), common_term_83);
    // q = -17
    result_by_lpq.set(c_Key3(8, 2, -17), common_term_84);
    result_by_q.set(c_Key1(-17), common_term_84);
    // q = -16
    result_by_lpq.set(c_Key3(8, 2, -16), common_term_85);
    result_by_q.set(c_Key1(-16), common_term_85);
    // q = -15
    result_by_lpq.set(c_Key3(8, 2, -15), common_term_86);
    result_by_q.set(c_Key1(-15), common_term_86);
    // q = -14
    result_by_lpq.set(c_Key3(8, 2, -14), common_term_87);
    result_by_q.set(c_Key1(-14), common_term_87);
    // q = -13
    result_by_lpq.set(c_Key3(8, 2, -13), common_term_88);
    result_by_q.set(c_Key1(-13), common_term_88);
    // q = -12
    result_by_lpq.set(c_Key3(8, 2, -12), common_term_89);
    result_by_q.set(c_Key1(-12), common_term_89);
    // q = -11
    result_by_lpq.set(c_Key3(8, 2, -11), common_term_90);
    result_by_q.set(c_Key1(-11), common_term_90);
    // q = -10
    result_by_lpq.set(c_Key3(8, 2, -10), common_term_91);
    result_by_q.set(c_Key1(-10), common_term_91);
    // q = -9
    result_by_lpq.set(c_Key3(8, 2, -9), common_term_92);
    result_by_q.set(c_Key1(-9), common_term_92);
    // q = -8
    result_by_lpq.set(c_Key3(8, 2, -8), common_term_93);
    result_by_q.set(c_Key1(-8), common_term_93);
    // q = -7
    result_by_lpq.set(c_Key3(8, 2, -7), common_term_94);
    result_by_q.set(c_Key1(-7), common_term_94);
    // q = -6
    result_by_lpq.set(c_Key3(8, 2, -6), common_term_95);
    result_by_q.set(c_Key1(-6), common_term_95);
    // q = -5
    result_by_lpq.set(c_Key3(8, 2, -5), common_term_96);
    result_by_q.set(c_Key1(-5), common_term_96);
    // q = -4
    result_by_lpq.set(c_Key3(8, 2, -4), common_term_97);
    result_by_q.set(c_Key1(-4), common_term_97);
    // q = -3
    result_by_lpq.set(c_Key3(8, 2, -3), common_term_98);
    result_by_q.set(c_Key1(-3), common_term_98);
    // q = -2
    result_by_lpq.set(c_Key3(8, 2, -2), common_term_99);
    result_by_q.set(c_Key1(-2), common_term_99);
    // q = -1
    result_by_lpq.set(c_Key3(8, 2, -1), common_term_100);
    result_by_q.set(c_Key1(-1), common_term_100);
    // q = 0
    result_by_lpq.set(c_Key3(8, 2, 0), common_term_101);
    result_by_q.set(c_Key1(0), common_term_101);
    // q = 1
    result_by_lpq.set(c_Key3(8, 2, 1), common_term_102);
    result_by_q.set(c_Key1(1), common_term_102);
    // q = 2
    result_by_lpq.set(c_Key3(8, 2, 2), common_term_103);
    result_by_q.set(c_Key1(2), common_term_103);
    // q = 3
    result_by_lpq.set(c_Key3(8, 2, 3), common_term_104);
    result_by_q.set(c_Key1(3), common_term_104);
    // q = 4
    result_by_lpq.set(c_Key3(8, 2, 4), common_term_105);
    result_by_q.set(c_Key1(4), common_term_105);
    // q = 5
    result_by_lpq.set(c_Key3(8, 2, 5), common_term_106);
    result_by_q.set(c_Key1(5), common_term_106);
    // q = 6
    result_by_lpq.set(c_Key3(8, 2, 6), common_term_107);
    result_by_q.set(c_Key1(6), common_term_107);
    // q = 7
    result_by_lpq.set(c_Key3(8, 2, 7), common_term_108);
    result_by_q.set(c_Key1(7), common_term_108);
    // q = 8
    result_by_lpq.set(c_Key3(8, 2, 8), common_term_109);
    result_by_q.set(c_Key1(8), common_term_109);
    // q = 9
    result_by_lpq.set(c_Key3(8, 2, 9), common_term_110);
    result_by_q.set(c_Key1(9), common_term_110);
    // q = 10
    result_by_lpq.set(c_Key3(8, 2, 10), common_term_111);
    result_by_q.set(c_Key1(10), common_term_111);
    // q = 11
    result_by_lpq.set(c_Key3(8, 2, 11), common_term_112);
    result_by_q.set(c_Key1(11), common_term_112);
    // q = 12
    result_by_lpq.set(c_Key3(8, 2, 12), common_term_113);
    result_by_q.set(c_Key1(12), common_term_113);
    // q = 13
    result_by_lpq.set(c_Key3(8, 2, 13), common_term_114);
    result_by_q.set(c_Key1(13), common_term_114);
    // q = 14
    result_by_lpq.set(c_Key3(8, 2, 14), common_term_115);
    result_by_q.set(c_Key1(14), common_term_115);
    // q = 15
    result_by_lpq.set(c_Key3(8, 2, 15), common_term_116);
    result_by_q.set(c_Key1(15), common_term_116);
    // q = 16
    result_by_lpq.set(c_Key3(8, 2, 16), common_term_117);
    result_by_q.set(c_Key1(16), common_term_117);
    // q = 17
    result_by_lpq.set(c_Key3(8, 2, 17), common_term_118);
    result_by_q.set(c_Key1(17), common_term_118);
    // q = 18
    result_by_lpq.set(c_Key3(8, 2, 18), common_term_119);
    result_by_q.set(c_Key1(18), common_term_119);
    // q = 19
    result_by_lpq.set(c_Key3(8, 2, 19), common_term_120);
    result_by_q.set(c_Key1(19), common_term_120);
    // q = 20
    result_by_lpq.set(c_Key3(8, 2, 20), common_term_121);
    result_by_q.set(c_Key1(20), common_term_121);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 2), result_by_q);
    result_by_q.clear();

    // l , p = (8, 3).
    // q = -20
    result_by_lpq.set(c_Key3(8, 3, -20), common_term_122);
    result_by_q.set(c_Key1(-20), common_term_122);
    // q = -19
    result_by_lpq.set(c_Key3(8, 3, -19), common_term_123);
    result_by_q.set(c_Key1(-19), common_term_123);
    // q = -18
    result_by_lpq.set(c_Key3(8, 3, -18), common_term_124);
    result_by_q.set(c_Key1(-18), common_term_124);
    // q = -17
    result_by_lpq.set(c_Key3(8, 3, -17), common_term_125);
    result_by_q.set(c_Key1(-17), common_term_125);
    // q = -16
    result_by_lpq.set(c_Key3(8, 3, -16), common_term_126);
    result_by_q.set(c_Key1(-16), common_term_126);
    // q = -15
    result_by_lpq.set(c_Key3(8, 3, -15), common_term_127);
    result_by_q.set(c_Key1(-15), common_term_127);
    // q = -14
    result_by_lpq.set(c_Key3(8, 3, -14), common_term_128);
    result_by_q.set(c_Key1(-14), common_term_128);
    // q = -13
    result_by_lpq.set(c_Key3(8, 3, -13), common_term_129);
    result_by_q.set(c_Key1(-13), common_term_129);
    // q = -12
    result_by_lpq.set(c_Key3(8, 3, -12), common_term_130);
    result_by_q.set(c_Key1(-12), common_term_130);
    // q = -11
    result_by_lpq.set(c_Key3(8, 3, -11), common_term_131);
    result_by_q.set(c_Key1(-11), common_term_131);
    // q = -10
    result_by_lpq.set(c_Key3(8, 3, -10), common_term_132);
    result_by_q.set(c_Key1(-10), common_term_132);
    // q = -9
    result_by_lpq.set(c_Key3(8, 3, -9), common_term_133);
    result_by_q.set(c_Key1(-9), common_term_133);
    // q = -8
    result_by_lpq.set(c_Key3(8, 3, -8), common_term_134);
    result_by_q.set(c_Key1(-8), common_term_134);
    // q = -7
    result_by_lpq.set(c_Key3(8, 3, -7), common_term_135);
    result_by_q.set(c_Key1(-7), common_term_135);
    // q = -6
    result_by_lpq.set(c_Key3(8, 3, -6), common_term_136);
    result_by_q.set(c_Key1(-6), common_term_136);
    // q = -5
    result_by_lpq.set(c_Key3(8, 3, -5), common_term_137);
    result_by_q.set(c_Key1(-5), common_term_137);
    // q = -4
    result_by_lpq.set(c_Key3(8, 3, -4), common_term_138);
    result_by_q.set(c_Key1(-4), common_term_138);
    // q = -3
    result_by_lpq.set(c_Key3(8, 3, -3), common_term_139);
    result_by_q.set(c_Key1(-3), common_term_139);
    // q = -2
    result_by_lpq.set(c_Key3(8, 3, -2), common_term_140);
    result_by_q.set(c_Key1(-2), common_term_140);
    // q = -1
    result_by_lpq.set(c_Key3(8, 3, -1), common_term_141);
    result_by_q.set(c_Key1(-1), common_term_141);
    // q = 0
    result_by_lpq.set(c_Key3(8, 3, 0), common_term_142);
    result_by_q.set(c_Key1(0), common_term_142);
    // q = 1
    result_by_lpq.set(c_Key3(8, 3, 1), common_term_143);
    result_by_q.set(c_Key1(1), common_term_143);
    // q = 2
    result_by_lpq.set(c_Key3(8, 3, 2), common_term_144);
    result_by_q.set(c_Key1(2), common_term_144);
    // q = 3
    result_by_lpq.set(c_Key3(8, 3, 3), common_term_145);
    result_by_q.set(c_Key1(3), common_term_145);
    // q = 4
    result_by_lpq.set(c_Key3(8, 3, 4), common_term_146);
    result_by_q.set(c_Key1(4), common_term_146);
    // q = 5
    result_by_lpq.set(c_Key3(8, 3, 5), common_term_147);
    result_by_q.set(c_Key1(5), common_term_147);
    // q = 6
    result_by_lpq.set(c_Key3(8, 3, 6), common_term_148);
    result_by_q.set(c_Key1(6), common_term_148);
    // q = 7
    result_by_lpq.set(c_Key3(8, 3, 7), common_term_149);
    result_by_q.set(c_Key1(7), common_term_149);
    // q = 8
    result_by_lpq.set(c_Key3(8, 3, 8), common_term_150);
    result_by_q.set(c_Key1(8), common_term_150);
    // q = 9
    result_by_lpq.set(c_Key3(8, 3, 9), common_term_151);
    result_by_q.set(c_Key1(9), common_term_151);
    // q = 10
    result_by_lpq.set(c_Key3(8, 3, 10), common_term_152);
    result_by_q.set(c_Key1(10), common_term_152);
    // q = 11
    result_by_lpq.set(c_Key3(8, 3, 11), common_term_153);
    result_by_q.set(c_Key1(11), common_term_153);
    // q = 12
    result_by_lpq.set(c_Key3(8, 3, 12), common_term_154);
    result_by_q.set(c_Key1(12), common_term_154);
    // q = 13
    result_by_lpq.set(c_Key3(8, 3, 13), common_term_155);
    result_by_q.set(c_Key1(13), common_term_155);
    // q = 14
    result_by_lpq.set(c_Key3(8, 3, 14), common_term_156);
    result_by_q.set(c_Key1(14), common_term_156);
    // q = 15
    result_by_lpq.set(c_Key3(8, 3, 15), common_term_157);
    result_by_q.set(c_Key1(15), common_term_157);
    // q = 16
    result_by_lpq.set(c_Key3(8, 3, 16), common_term_158);
    result_by_q.set(c_Key1(16), common_term_158);
    // q = 17
    result_by_lpq.set(c_Key3(8, 3, 17), common_term_159);
    result_by_q.set(c_Key1(17), common_term_159);
    // q = 18
    result_by_lpq.set(c_Key3(8, 3, 18), common_term_160);
    result_by_q.set(c_Key1(18), common_term_160);
    // q = 19
    result_by_lpq.set(c_Key3(8, 3, 19), common_term_161);
    result_by_q.set(c_Key1(19), common_term_161);
    // q = 20
    result_by_lpq.set(c_Key3(8, 3, 20), common_term_162);
    result_by_q.set(c_Key1(20), common_term_162);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 3), result_by_q);
    result_by_q.clear();

    // l , p = (8, 4).
    // q = -20
    result_by_lpq.set(c_Key3(8, 4, -20), common_term_163);
    result_by_q.set(c_Key1(-20), common_term_163);
    // q = -19
    result_by_lpq.set(c_Key3(8, 4, -19), common_term_164);
    result_by_q.set(c_Key1(-19), common_term_164);
    // q = -18
    result_by_lpq.set(c_Key3(8, 4, -18), common_term_165);
    result_by_q.set(c_Key1(-18), common_term_165);
    // q = -17
    result_by_lpq.set(c_Key3(8, 4, -17), common_term_166);
    result_by_q.set(c_Key1(-17), common_term_166);
    // q = -16
    result_by_lpq.set(c_Key3(8, 4, -16), common_term_167);
    result_by_q.set(c_Key1(-16), common_term_167);
    // q = -15
    result_by_lpq.set(c_Key3(8, 4, -15), common_term_168);
    result_by_q.set(c_Key1(-15), common_term_168);
    // q = -14
    result_by_lpq.set(c_Key3(8, 4, -14), common_term_169);
    result_by_q.set(c_Key1(-14), common_term_169);
    // q = -13
    result_by_lpq.set(c_Key3(8, 4, -13), common_term_170);
    result_by_q.set(c_Key1(-13), common_term_170);
    // q = -12
    result_by_lpq.set(c_Key3(8, 4, -12), common_term_171);
    result_by_q.set(c_Key1(-12), common_term_171);
    // q = -11
    result_by_lpq.set(c_Key3(8, 4, -11), common_term_172);
    result_by_q.set(c_Key1(-11), common_term_172);
    // q = -10
    result_by_lpq.set(c_Key3(8, 4, -10), common_term_173);
    result_by_q.set(c_Key1(-10), common_term_173);
    // q = -9
    result_by_lpq.set(c_Key3(8, 4, -9), common_term_174);
    result_by_q.set(c_Key1(-9), common_term_174);
    // q = -8
    result_by_lpq.set(c_Key3(8, 4, -8), common_term_175);
    result_by_q.set(c_Key1(-8), common_term_175);
    // q = -7
    result_by_lpq.set(c_Key3(8, 4, -7), common_term_176);
    result_by_q.set(c_Key1(-7), common_term_176);
    // q = -6
    result_by_lpq.set(c_Key3(8, 4, -6), common_term_177);
    result_by_q.set(c_Key1(-6), common_term_177);
    // q = -5
    result_by_lpq.set(c_Key3(8, 4, -5), common_term_178);
    result_by_q.set(c_Key1(-5), common_term_178);
    // q = -4
    result_by_lpq.set(c_Key3(8, 4, -4), common_term_179);
    result_by_q.set(c_Key1(-4), common_term_179);
    // q = -3
    result_by_lpq.set(c_Key3(8, 4, -3), common_term_180);
    result_by_q.set(c_Key1(-3), common_term_180);
    // q = -2
    result_by_lpq.set(c_Key3(8, 4, -2), common_term_181);
    result_by_q.set(c_Key1(-2), common_term_181);
    // q = -1
    result_by_lpq.set(c_Key3(8, 4, -1), common_term_182);
    result_by_q.set(c_Key1(-1), common_term_182);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -7.5)*(2.1875*eccentricity_6 + 13.125*eccentricity_4 + 10.5*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(8, 4, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(8, 4, 1), common_term_182);
    result_by_q.set(c_Key1(1), common_term_182);
    // q = 2
    result_by_lpq.set(c_Key3(8, 4, 2), common_term_181);
    result_by_q.set(c_Key1(2), common_term_181);
    // q = 3
    result_by_lpq.set(c_Key3(8, 4, 3), common_term_180);
    result_by_q.set(c_Key1(3), common_term_180);
    // q = 4
    result_by_lpq.set(c_Key3(8, 4, 4), common_term_179);
    result_by_q.set(c_Key1(4), common_term_179);
    // q = 5
    result_by_lpq.set(c_Key3(8, 4, 5), common_term_178);
    result_by_q.set(c_Key1(5), common_term_178);
    // q = 6
    result_by_lpq.set(c_Key3(8, 4, 6), common_term_177);
    result_by_q.set(c_Key1(6), common_term_177);
    // q = 7
    result_by_lpq.set(c_Key3(8, 4, 7), common_term_176);
    result_by_q.set(c_Key1(7), common_term_176);
    // q = 8
    result_by_lpq.set(c_Key3(8, 4, 8), common_term_175);
    result_by_q.set(c_Key1(8), common_term_175);
    // q = 9
    result_by_lpq.set(c_Key3(8, 4, 9), common_term_174);
    result_by_q.set(c_Key1(9), common_term_174);
    // q = 10
    result_by_lpq.set(c_Key3(8, 4, 10), common_term_173);
    result_by_q.set(c_Key1(10), common_term_173);
    // q = 11
    result_by_lpq.set(c_Key3(8, 4, 11), common_term_172);
    result_by_q.set(c_Key1(11), common_term_172);
    // q = 12
    result_by_lpq.set(c_Key3(8, 4, 12), common_term_171);
    result_by_q.set(c_Key1(12), common_term_171);
    // q = 13
    result_by_lpq.set(c_Key3(8, 4, 13), common_term_170);
    result_by_q.set(c_Key1(13), common_term_170);
    // q = 14
    result_by_lpq.set(c_Key3(8, 4, 14), common_term_169);
    result_by_q.set(c_Key1(14), common_term_169);
    // q = 15
    result_by_lpq.set(c_Key3(8, 4, 15), common_term_168);
    result_by_q.set(c_Key1(15), common_term_168);
    // q = 16
    result_by_lpq.set(c_Key3(8, 4, 16), common_term_167);
    result_by_q.set(c_Key1(16), common_term_167);
    // q = 17
    result_by_lpq.set(c_Key3(8, 4, 17), common_term_166);
    result_by_q.set(c_Key1(17), common_term_166);
    // q = 18
    result_by_lpq.set(c_Key3(8, 4, 18), common_term_165);
    result_by_q.set(c_Key1(18), common_term_165);
    // q = 19
    result_by_lpq.set(c_Key3(8, 4, 19), common_term_164);
    result_by_q.set(c_Key1(19), common_term_164);
    // q = 20
    result_by_lpq.set(c_Key3(8, 4, 20), common_term_163);
    result_by_q.set(c_Key1(20), common_term_163);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 4), result_by_q);
    result_by_q.clear();

    // l , p = (8, 5).
    // q = -20
    result_by_lpq.set(c_Key3(8, 5, -20), common_term_162);
    result_by_q.set(c_Key1(-20), common_term_162);
    // q = -19
    result_by_lpq.set(c_Key3(8, 5, -19), common_term_161);
    result_by_q.set(c_Key1(-19), common_term_161);
    // q = -18
    result_by_lpq.set(c_Key3(8, 5, -18), common_term_160);
    result_by_q.set(c_Key1(-18), common_term_160);
    // q = -17
    result_by_lpq.set(c_Key3(8, 5, -17), common_term_159);
    result_by_q.set(c_Key1(-17), common_term_159);
    // q = -16
    result_by_lpq.set(c_Key3(8, 5, -16), common_term_158);
    result_by_q.set(c_Key1(-16), common_term_158);
    // q = -15
    result_by_lpq.set(c_Key3(8, 5, -15), common_term_157);
    result_by_q.set(c_Key1(-15), common_term_157);
    // q = -14
    result_by_lpq.set(c_Key3(8, 5, -14), common_term_156);
    result_by_q.set(c_Key1(-14), common_term_156);
    // q = -13
    result_by_lpq.set(c_Key3(8, 5, -13), common_term_155);
    result_by_q.set(c_Key1(-13), common_term_155);
    // q = -12
    result_by_lpq.set(c_Key3(8, 5, -12), common_term_154);
    result_by_q.set(c_Key1(-12), common_term_154);
    // q = -11
    result_by_lpq.set(c_Key3(8, 5, -11), common_term_153);
    result_by_q.set(c_Key1(-11), common_term_153);
    // q = -10
    result_by_lpq.set(c_Key3(8, 5, -10), common_term_152);
    result_by_q.set(c_Key1(-10), common_term_152);
    // q = -9
    result_by_lpq.set(c_Key3(8, 5, -9), common_term_151);
    result_by_q.set(c_Key1(-9), common_term_151);
    // q = -8
    result_by_lpq.set(c_Key3(8, 5, -8), common_term_150);
    result_by_q.set(c_Key1(-8), common_term_150);
    // q = -7
    result_by_lpq.set(c_Key3(8, 5, -7), common_term_149);
    result_by_q.set(c_Key1(-7), common_term_149);
    // q = -6
    result_by_lpq.set(c_Key3(8, 5, -6), common_term_148);
    result_by_q.set(c_Key1(-6), common_term_148);
    // q = -5
    result_by_lpq.set(c_Key3(8, 5, -5), common_term_147);
    result_by_q.set(c_Key1(-5), common_term_147);
    // q = -4
    result_by_lpq.set(c_Key3(8, 5, -4), common_term_146);
    result_by_q.set(c_Key1(-4), common_term_146);
    // q = -3
    result_by_lpq.set(c_Key3(8, 5, -3), common_term_145);
    result_by_q.set(c_Key1(-3), common_term_145);
    // q = -2
    result_by_lpq.set(c_Key3(8, 5, -2), common_term_144);
    result_by_q.set(c_Key1(-2), common_term_144);
    // q = -1
    result_by_lpq.set(c_Key3(8, 5, -1), common_term_143);
    result_by_q.set(c_Key1(-1), common_term_143);
    // q = 0
    result_by_lpq.set(c_Key3(8, 5, 0), common_term_142);
    result_by_q.set(c_Key1(0), common_term_142);
    // q = 1
    result_by_lpq.set(c_Key3(8, 5, 1), common_term_141);
    result_by_q.set(c_Key1(1), common_term_141);
    // q = 2
    result_by_lpq.set(c_Key3(8, 5, 2), common_term_140);
    result_by_q.set(c_Key1(2), common_term_140);
    // q = 3
    result_by_lpq.set(c_Key3(8, 5, 3), common_term_139);
    result_by_q.set(c_Key1(3), common_term_139);
    // q = 4
    result_by_lpq.set(c_Key3(8, 5, 4), common_term_138);
    result_by_q.set(c_Key1(4), common_term_138);
    // q = 5
    result_by_lpq.set(c_Key3(8, 5, 5), common_term_137);
    result_by_q.set(c_Key1(5), common_term_137);
    // q = 6
    result_by_lpq.set(c_Key3(8, 5, 6), common_term_136);
    result_by_q.set(c_Key1(6), common_term_136);
    // q = 7
    result_by_lpq.set(c_Key3(8, 5, 7), common_term_135);
    result_by_q.set(c_Key1(7), common_term_135);
    // q = 8
    result_by_lpq.set(c_Key3(8, 5, 8), common_term_134);
    result_by_q.set(c_Key1(8), common_term_134);
    // q = 9
    result_by_lpq.set(c_Key3(8, 5, 9), common_term_133);
    result_by_q.set(c_Key1(9), common_term_133);
    // q = 10
    result_by_lpq.set(c_Key3(8, 5, 10), common_term_132);
    result_by_q.set(c_Key1(10), common_term_132);
    // q = 11
    result_by_lpq.set(c_Key3(8, 5, 11), common_term_131);
    result_by_q.set(c_Key1(11), common_term_131);
    // q = 12
    result_by_lpq.set(c_Key3(8, 5, 12), common_term_130);
    result_by_q.set(c_Key1(12), common_term_130);
    // q = 13
    result_by_lpq.set(c_Key3(8, 5, 13), common_term_129);
    result_by_q.set(c_Key1(13), common_term_129);
    // q = 14
    result_by_lpq.set(c_Key3(8, 5, 14), common_term_128);
    result_by_q.set(c_Key1(14), common_term_128);
    // q = 15
    result_by_lpq.set(c_Key3(8, 5, 15), common_term_127);
    result_by_q.set(c_Key1(15), common_term_127);
    // q = 16
    result_by_lpq.set(c_Key3(8, 5, 16), common_term_126);
    result_by_q.set(c_Key1(16), common_term_126);
    // q = 17
    result_by_lpq.set(c_Key3(8, 5, 17), common_term_125);
    result_by_q.set(c_Key1(17), common_term_125);
    // q = 18
    result_by_lpq.set(c_Key3(8, 5, 18), common_term_124);
    result_by_q.set(c_Key1(18), common_term_124);
    // q = 19
    result_by_lpq.set(c_Key3(8, 5, 19), common_term_123);
    result_by_q.set(c_Key1(19), common_term_123);
    // q = 20
    result_by_lpq.set(c_Key3(8, 5, 20), common_term_122);
    result_by_q.set(c_Key1(20), common_term_122);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 5), result_by_q);
    result_by_q.clear();

    // l , p = (8, 6).
    // q = -20
    result_by_lpq.set(c_Key3(8, 6, -20), common_term_121);
    result_by_q.set(c_Key1(-20), common_term_121);
    // q = -19
    result_by_lpq.set(c_Key3(8, 6, -19), common_term_120);
    result_by_q.set(c_Key1(-19), common_term_120);
    // q = -18
    result_by_lpq.set(c_Key3(8, 6, -18), common_term_119);
    result_by_q.set(c_Key1(-18), common_term_119);
    // q = -17
    result_by_lpq.set(c_Key3(8, 6, -17), common_term_118);
    result_by_q.set(c_Key1(-17), common_term_118);
    // q = -16
    result_by_lpq.set(c_Key3(8, 6, -16), common_term_117);
    result_by_q.set(c_Key1(-16), common_term_117);
    // q = -15
    result_by_lpq.set(c_Key3(8, 6, -15), common_term_116);
    result_by_q.set(c_Key1(-15), common_term_116);
    // q = -14
    result_by_lpq.set(c_Key3(8, 6, -14), common_term_115);
    result_by_q.set(c_Key1(-14), common_term_115);
    // q = -13
    result_by_lpq.set(c_Key3(8, 6, -13), common_term_114);
    result_by_q.set(c_Key1(-13), common_term_114);
    // q = -12
    result_by_lpq.set(c_Key3(8, 6, -12), common_term_113);
    result_by_q.set(c_Key1(-12), common_term_113);
    // q = -11
    result_by_lpq.set(c_Key3(8, 6, -11), common_term_112);
    result_by_q.set(c_Key1(-11), common_term_112);
    // q = -10
    result_by_lpq.set(c_Key3(8, 6, -10), common_term_111);
    result_by_q.set(c_Key1(-10), common_term_111);
    // q = -9
    result_by_lpq.set(c_Key3(8, 6, -9), common_term_110);
    result_by_q.set(c_Key1(-9), common_term_110);
    // q = -8
    result_by_lpq.set(c_Key3(8, 6, -8), common_term_109);
    result_by_q.set(c_Key1(-8), common_term_109);
    // q = -7
    result_by_lpq.set(c_Key3(8, 6, -7), common_term_108);
    result_by_q.set(c_Key1(-7), common_term_108);
    // q = -6
    result_by_lpq.set(c_Key3(8, 6, -6), common_term_107);
    result_by_q.set(c_Key1(-6), common_term_107);
    // q = -5
    result_by_lpq.set(c_Key3(8, 6, -5), common_term_106);
    result_by_q.set(c_Key1(-5), common_term_106);
    // q = -4
    result_by_lpq.set(c_Key3(8, 6, -4), common_term_105);
    result_by_q.set(c_Key1(-4), common_term_105);
    // q = -3
    result_by_lpq.set(c_Key3(8, 6, -3), common_term_104);
    result_by_q.set(c_Key1(-3), common_term_104);
    // q = -2
    result_by_lpq.set(c_Key3(8, 6, -2), common_term_103);
    result_by_q.set(c_Key1(-2), common_term_103);
    // q = -1
    result_by_lpq.set(c_Key3(8, 6, -1), common_term_102);
    result_by_q.set(c_Key1(-1), common_term_102);
    // q = 0
    result_by_lpq.set(c_Key3(8, 6, 0), common_term_101);
    result_by_q.set(c_Key1(0), common_term_101);
    // q = 1
    result_by_lpq.set(c_Key3(8, 6, 1), common_term_100);
    result_by_q.set(c_Key1(1), common_term_100);
    // q = 2
    result_by_lpq.set(c_Key3(8, 6, 2), common_term_99);
    result_by_q.set(c_Key1(2), common_term_99);
    // q = 3
    result_by_lpq.set(c_Key3(8, 6, 3), common_term_98);
    result_by_q.set(c_Key1(3), common_term_98);
    // q = 4
    result_by_lpq.set(c_Key3(8, 6, 4), common_term_97);
    result_by_q.set(c_Key1(4), common_term_97);
    // q = 5
    result_by_lpq.set(c_Key3(8, 6, 5), common_term_96);
    result_by_q.set(c_Key1(5), common_term_96);
    // q = 6
    result_by_lpq.set(c_Key3(8, 6, 6), common_term_95);
    result_by_q.set(c_Key1(6), common_term_95);
    // q = 7
    result_by_lpq.set(c_Key3(8, 6, 7), common_term_94);
    result_by_q.set(c_Key1(7), common_term_94);
    // q = 8
    result_by_lpq.set(c_Key3(8, 6, 8), common_term_93);
    result_by_q.set(c_Key1(8), common_term_93);
    // q = 9
    result_by_lpq.set(c_Key3(8, 6, 9), common_term_92);
    result_by_q.set(c_Key1(9), common_term_92);
    // q = 10
    result_by_lpq.set(c_Key3(8, 6, 10), common_term_91);
    result_by_q.set(c_Key1(10), common_term_91);
    // q = 11
    result_by_lpq.set(c_Key3(8, 6, 11), common_term_90);
    result_by_q.set(c_Key1(11), common_term_90);
    // q = 12
    result_by_lpq.set(c_Key3(8, 6, 12), common_term_89);
    result_by_q.set(c_Key1(12), common_term_89);
    // q = 13
    result_by_lpq.set(c_Key3(8, 6, 13), common_term_88);
    result_by_q.set(c_Key1(13), common_term_88);
    // q = 14
    result_by_lpq.set(c_Key3(8, 6, 14), common_term_87);
    result_by_q.set(c_Key1(14), common_term_87);
    // q = 15
    result_by_lpq.set(c_Key3(8, 6, 15), common_term_86);
    result_by_q.set(c_Key1(15), common_term_86);
    // q = 16
    result_by_lpq.set(c_Key3(8, 6, 16), common_term_85);
    result_by_q.set(c_Key1(16), common_term_85);
    // q = 17
    result_by_lpq.set(c_Key3(8, 6, 17), common_term_84);
    result_by_q.set(c_Key1(17), common_term_84);
    // q = 18
    result_by_lpq.set(c_Key3(8, 6, 18), common_term_83);
    result_by_q.set(c_Key1(18), common_term_83);
    // q = 19
    result_by_lpq.set(c_Key3(8, 6, 19), common_term_82);
    result_by_q.set(c_Key1(19), common_term_82);
    // q = 20
    result_by_lpq.set(c_Key3(8, 6, 20), common_term_81);
    result_by_q.set(c_Key1(20), common_term_81);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 6), result_by_q);
    result_by_q.clear();

    // l , p = (8, 7).
    // q = -20
    result_by_lpq.set(c_Key3(8, 7, -20), common_term_80);
    result_by_q.set(c_Key1(-20), common_term_80);
    // q = -19
    result_by_lpq.set(c_Key3(8, 7, -19), common_term_79);
    result_by_q.set(c_Key1(-19), common_term_79);
    // q = -18
    result_by_lpq.set(c_Key3(8, 7, -18), common_term_78);
    result_by_q.set(c_Key1(-18), common_term_78);
    // q = -17
    result_by_lpq.set(c_Key3(8, 7, -17), common_term_77);
    result_by_q.set(c_Key1(-17), common_term_77);
    // q = -16
    result_by_lpq.set(c_Key3(8, 7, -16), common_term_76);
    result_by_q.set(c_Key1(-16), common_term_76);
    // q = -15
    result_by_lpq.set(c_Key3(8, 7, -15), common_term_75);
    result_by_q.set(c_Key1(-15), common_term_75);
    // q = -14
    result_by_lpq.set(c_Key3(8, 7, -14), common_term_74);
    result_by_q.set(c_Key1(-14), common_term_74);
    // q = -13
    result_by_lpq.set(c_Key3(8, 7, -13), common_term_73);
    result_by_q.set(c_Key1(-13), common_term_73);
    // q = -12
    result_by_lpq.set(c_Key3(8, 7, -12), common_term_72);
    result_by_q.set(c_Key1(-12), common_term_72);
    // q = -11
    result_by_lpq.set(c_Key3(8, 7, -11), common_term_71);
    result_by_q.set(c_Key1(-11), common_term_71);
    // q = -10
    result_by_lpq.set(c_Key3(8, 7, -10), common_term_70);
    result_by_q.set(c_Key1(-10), common_term_70);
    // q = -9
    result_by_lpq.set(c_Key3(8, 7, -9), common_term_69);
    result_by_q.set(c_Key1(-9), common_term_69);
    // q = -8
    result_by_lpq.set(c_Key3(8, 7, -8), common_term_68);
    result_by_q.set(c_Key1(-8), common_term_68);
    // q = -7
    result_by_lpq.set(c_Key3(8, 7, -7), common_term_67);
    result_by_q.set(c_Key1(-7), common_term_67);
    // q = -6
    result_by_lpq.set(c_Key3(8, 7, -6), common_term_66);
    result_by_q.set(c_Key1(-6), common_term_66);
    // q = -5
    result_by_lpq.set(c_Key3(8, 7, -5), common_term_65);
    result_by_q.set(c_Key1(-5), common_term_65);
    // q = -4
    result_by_lpq.set(c_Key3(8, 7, -4), common_term_64);
    result_by_q.set(c_Key1(-4), common_term_64);
    // q = -3
    result_by_lpq.set(c_Key3(8, 7, -3), common_term_63);
    result_by_q.set(c_Key1(-3), common_term_63);
    // q = -2
    result_by_lpq.set(c_Key3(8, 7, -2), common_term_62);
    result_by_q.set(c_Key1(-2), common_term_62);
    // q = -1
    result_by_lpq.set(c_Key3(8, 7, -1), common_term_61);
    result_by_q.set(c_Key1(-1), common_term_61);
    // q = 0
    result_by_lpq.set(c_Key3(8, 7, 0), common_term_60);
    result_by_q.set(c_Key1(0), common_term_60);
    // q = 1
    result_by_lpq.set(c_Key3(8, 7, 1), common_term_59);
    result_by_q.set(c_Key1(1), common_term_59);
    // q = 2
    result_by_lpq.set(c_Key3(8, 7, 2), common_term_58);
    result_by_q.set(c_Key1(2), common_term_58);
    // q = 3
    result_by_lpq.set(c_Key3(8, 7, 3), common_term_57);
    result_by_q.set(c_Key1(3), common_term_57);
    // q = 4
    result_by_lpq.set(c_Key3(8, 7, 4), common_term_56);
    result_by_q.set(c_Key1(4), common_term_56);
    // q = 5
    result_by_lpq.set(c_Key3(8, 7, 5), common_term_55);
    result_by_q.set(c_Key1(5), common_term_55);
    // q = 6
    result_by_lpq.set(c_Key3(8, 7, 6), common_term_54);
    result_by_q.set(c_Key1(6), common_term_54);
    // q = 7
    result_by_lpq.set(c_Key3(8, 7, 7), common_term_53);
    result_by_q.set(c_Key1(7), common_term_53);
    // q = 8
    result_by_lpq.set(c_Key3(8, 7, 8), common_term_52);
    result_by_q.set(c_Key1(8), common_term_52);
    // q = 9
    result_by_lpq.set(c_Key3(8, 7, 9), common_term_51);
    result_by_q.set(c_Key1(9), common_term_51);
    // q = 10
    result_by_lpq.set(c_Key3(8, 7, 10), common_term_50);
    result_by_q.set(c_Key1(10), common_term_50);
    // q = 11
    result_by_lpq.set(c_Key3(8, 7, 11), common_term_49);
    result_by_q.set(c_Key1(11), common_term_49);
    // q = 12
    result_by_lpq.set(c_Key3(8, 7, 12), common_term_48);
    result_by_q.set(c_Key1(12), common_term_48);
    // q = 13
    result_by_lpq.set(c_Key3(8, 7, 13), common_term_47);
    result_by_q.set(c_Key1(13), common_term_47);
    // q = 14
    result_by_lpq.set(c_Key3(8, 7, 14), common_term_46);
    result_by_q.set(c_Key1(14), common_term_46);
    // q = 15
    result_by_lpq.set(c_Key3(8, 7, 15), common_term_45);
    result_by_q.set(c_Key1(15), common_term_45);
    // q = 16
    result_by_lpq.set(c_Key3(8, 7, 16), common_term_44);
    result_by_q.set(c_Key1(16), common_term_44);
    // q = 17
    result_by_lpq.set(c_Key3(8, 7, 17), common_term_43);
    result_by_q.set(c_Key1(17), common_term_43);
    // q = 18
    result_by_lpq.set(c_Key3(8, 7, 18), common_term_42);
    result_by_q.set(c_Key1(18), common_term_42);
    // q = 19
    result_by_lpq.set(c_Key3(8, 7, 19), common_term_41);
    result_by_q.set(c_Key1(19), common_term_41);
    // q = 20
    result_by_lpq.set(c_Key3(8, 7, 20), common_term_40);
    result_by_q.set(c_Key1(20), common_term_40);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 7), result_by_q);
    result_by_q.clear();

    // l , p = (8, 8).
    // q = -20
    result_by_lpq.set(c_Key3(8, 8, -20), common_term_39);
    result_by_q.set(c_Key1(-20), common_term_39);
    // q = -19
    result_by_lpq.set(c_Key3(8, 8, -19), common_term_38);
    result_by_q.set(c_Key1(-19), common_term_38);
    // q = -18
    result_by_lpq.set(c_Key3(8, 8, -18), common_term_37);
    result_by_q.set(c_Key1(-18), common_term_37);
    // q = -17
    result_by_lpq.set(c_Key3(8, 8, -17), common_term_36);
    result_by_q.set(c_Key1(-17), common_term_36);
    // q = -16
    result_by_lpq.set(c_Key3(8, 8, -16), common_term_35);
    result_by_q.set(c_Key1(-16), common_term_35);
    // q = -15
    result_by_lpq.set(c_Key3(8, 8, -15), common_term_34);
    result_by_q.set(c_Key1(-15), common_term_34);
    // q = -14
    result_by_lpq.set(c_Key3(8, 8, -14), common_term_33);
    result_by_q.set(c_Key1(-14), common_term_33);
    // q = -13
    result_by_lpq.set(c_Key3(8, 8, -13), common_term_32);
    result_by_q.set(c_Key1(-13), common_term_32);
    // q = -12
    result_by_lpq.set(c_Key3(8, 8, -12), common_term_31);
    result_by_q.set(c_Key1(-12), common_term_31);
    // q = -11
    result_by_lpq.set(c_Key3(8, 8, -11), common_term_30);
    result_by_q.set(c_Key1(-11), common_term_30);
    // q = -10
    result_by_lpq.set(c_Key3(8, 8, -10), common_term_29);
    result_by_q.set(c_Key1(-10), common_term_29);
    // q = -9
    result_by_lpq.set(c_Key3(8, 8, -9), common_term_28);
    result_by_q.set(c_Key1(-9), common_term_28);
    // q = -8
    result_by_lpq.set(c_Key3(8, 8, -8), common_term_27);
    result_by_q.set(c_Key1(-8), common_term_27);
    // q = -7
    result_by_lpq.set(c_Key3(8, 8, -7), common_term_26);
    result_by_q.set(c_Key1(-7), common_term_26);
    // q = -6
    result_by_lpq.set(c_Key3(8, 8, -6), common_term_25);
    result_by_q.set(c_Key1(-6), common_term_25);
    // q = -5
    result_by_lpq.set(c_Key3(8, 8, -5), common_term_24);
    result_by_q.set(c_Key1(-5), common_term_24);
    // q = -4
    result_by_lpq.set(c_Key3(8, 8, -4), common_term_23);
    result_by_q.set(c_Key1(-4), common_term_23);
    // q = -3
    result_by_lpq.set(c_Key3(8, 8, -3), common_term_22);
    result_by_q.set(c_Key1(-3), common_term_22);
    // q = -2
    result_by_lpq.set(c_Key3(8, 8, -2), common_term_21);
    result_by_q.set(c_Key1(-2), common_term_21);
    // q = -1
    result_by_lpq.set(c_Key3(8, 8, -1), common_term_20);
    result_by_q.set(c_Key1(-1), common_term_20);
    // q = 0
    result_by_lpq.set(c_Key3(8, 8, 0), common_term_19);
    result_by_q.set(c_Key1(0), common_term_19);
    // q = 1
    result_by_lpq.set(c_Key3(8, 8, 1), common_term_18);
    result_by_q.set(c_Key1(1), common_term_18);
    // q = 2
    result_by_lpq.set(c_Key3(8, 8, 2), common_term_17);
    result_by_q.set(c_Key1(2), common_term_17);
    // q = 3
    result_by_lpq.set(c_Key3(8, 8, 3), common_term_16);
    result_by_q.set(c_Key1(3), common_term_16);
    // q = 4
    result_by_lpq.set(c_Key3(8, 8, 4), common_term_15);
    result_by_q.set(c_Key1(4), common_term_15);
    // q = 5
    result_by_lpq.set(c_Key3(8, 8, 5), common_term_14);
    result_by_q.set(c_Key1(5), common_term_14);
    // q = 6
    result_by_lpq.set(c_Key3(8, 8, 6), common_term_13);
    result_by_q.set(c_Key1(6), common_term_13);
    // q = 7
    result_by_lpq.set(c_Key3(8, 8, 7), common_term_12);
    result_by_q.set(c_Key1(7), common_term_12);
    // q = 9
    result_by_lpq.set(c_Key3(8, 8, 9), common_term_11);
    result_by_q.set(c_Key1(9), common_term_11);
    // q = 10
    result_by_lpq.set(c_Key3(8, 8, 10), common_term_10);
    result_by_q.set(c_Key1(10), common_term_10);
    // q = 11
    result_by_lpq.set(c_Key3(8, 8, 11), common_term_9);
    result_by_q.set(c_Key1(11), common_term_9);
    // q = 12
    result_by_lpq.set(c_Key3(8, 8, 12), common_term_8);
    result_by_q.set(c_Key1(12), common_term_8);
    // q = 13
    result_by_lpq.set(c_Key3(8, 8, 13), common_term_7);
    result_by_q.set(c_Key1(13), common_term_7);
    // q = 14
    result_by_lpq.set(c_Key3(8, 8, 14), common_term_6);
    result_by_q.set(c_Key1(14), common_term_6);
    // q = 15
    result_by_lpq.set(c_Key3(8, 8, 15), common_term_5);
    result_by_q.set(c_Key1(15), common_term_5);
    // q = 16
    result_by_lpq.set(c_Key3(8, 8, 16), common_term_4);
    result_by_q.set(c_Key1(16), common_term_4);
    // q = 17
    result_by_lpq.set(c_Key3(8, 8, 17), common_term_3);
    result_by_q.set(c_Key1(17), common_term_3);
    // q = 18
    result_by_lpq.set(c_Key3(8, 8, 18), common_term_2);
    result_by_q.set(c_Key1(18), common_term_2);
    // q = 19
    result_by_lpq.set(c_Key3(8, 8, 19), common_term_1);
    result_by_q.set(c_Key1(19), common_term_1);
    // q = 20
    result_by_lpq.set(c_Key3(8, 8, 20), common_term_0);
    result_by_q.set(c_Key1(20), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(8, 8), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}
