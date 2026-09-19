#pragma once

#include "eccentricity_common_.hpp"

EccentricityFuncOutput c_eccentricity_function_l7_e1(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 7.
    // Functions are taylor expanded and truncated at eccentricity^1.

    // Functions Calculated for l = 7.

    c_IntMap<c_Key3, double> result_by_lpq(24);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(8);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = -3.0*eccentricity;
    double common_term_1 = 11.0*eccentricity;
    double common_term_2 = -eccentricity;
    double common_term_3 = 9.0*eccentricity;
    double common_term_4 = 7.0*eccentricity;
    double common_term_5 = 3.0*eccentricity*std::pow(1.0 - eccentricity_2, -6.5);
    double common_term_6 = 5.0*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(3);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (7, 0).
    // q = -1
    result_by_lpq.set(c_Key3(7, 0, -1), common_term_0);
    result_by_q.set(c_Key1(-1), common_term_0);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(7, 0, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(7, 0, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 0), result_by_q);
    result_by_q.clear();

    // l , p = (7, 1).
    // q = -1
    result_by_lpq.set(c_Key3(7, 1, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(7, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(7, 1, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 1), result_by_q);
    result_by_q.clear();

    // l , p = (7, 2).
    // q = -1
    result_by_lpq.set(c_Key3(7, 2, -1), eccentricity);
    result_by_q.set(c_Key1(-1), eccentricity);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(7, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(7, 2, 1), common_term_4);
    result_by_q.set(c_Key1(1), common_term_4);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 2), result_by_q);
    result_by_q.clear();

    // l , p = (7, 3).
    // q = -1
    result_by_lpq.set(c_Key3(7, 3, -1), common_term_5);
    result_by_q.set(c_Key1(-1), common_term_5);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(7, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(7, 3, 1), common_term_6);
    result_by_q.set(c_Key1(1), common_term_6);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 3), result_by_q);
    result_by_q.clear();

    // l , p = (7, 4).
    // q = -1
    result_by_lpq.set(c_Key3(7, 4, -1), common_term_6);
    result_by_q.set(c_Key1(-1), common_term_6);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(7, 4, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(7, 4, 1), common_term_5);
    result_by_q.set(c_Key1(1), common_term_5);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 4), result_by_q);
    result_by_q.clear();

    // l , p = (7, 5).
    // q = -1
    result_by_lpq.set(c_Key3(7, 5, -1), common_term_4);
    result_by_q.set(c_Key1(-1), common_term_4);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(7, 5, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(7, 5, 1), eccentricity);
    result_by_q.set(c_Key1(1), eccentricity);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 5), result_by_q);
    result_by_q.clear();

    // l , p = (7, 6).
    // q = -1
    result_by_lpq.set(c_Key3(7, 6, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(7, 6, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(7, 6, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 6), result_by_q);
    result_by_q.clear();

    // l , p = (7, 7).
    // q = -1
    result_by_lpq.set(c_Key3(7, 7, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(7, 7, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(7, 7, 1), common_term_0);
    result_by_q.set(c_Key1(1), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 7), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l7_e2(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 7.
    // Functions are taylor expanded and truncated at eccentricity^2.

    // Functions Calculated for l = 7.

    c_IntMap<c_Key3, double> result_by_lpq(40);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(8);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = 3.125*eccentricity_2;
    double common_term_1 = -3.0*eccentricity;
    double common_term_2 = 1.0 - 35.0*eccentricity_2;
    double common_term_3 = 11.0*eccentricity;
    double common_term_4 = 67.875*eccentricity_2;
    double common_term_5 = 0.375*eccentricity_2;
    double common_term_6 = -eccentricity;
    double common_term_7 = 1.0 - 11.0*eccentricity_2;
    double common_term_8 = 9.0*eccentricity;
    double common_term_9 = 46.625*eccentricity_2;
    double common_term_10 = 1.625*eccentricity_2;
    double common_term_11 = 5.0*eccentricity_2 + 1.0;
    double common_term_12 = 7.0*eccentricity;
    double common_term_13 = 29.375*eccentricity_2;
    double common_term_14 = 6.875*eccentricity_2;
    double common_term_15 = 3.0*eccentricity*std::pow(1.0 - eccentricity_2, -6.5);
    double common_term_16 = 13.0*eccentricity_2 + 1.0;
    double common_term_17 = 5.0*eccentricity;
    double common_term_18 = 16.125*eccentricity_2;

    c_IntMap<c_Key1, double> result_by_q(5);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (7, 0).
    // q = -2
    result_by_lpq.set(c_Key3(7, 0, -2), common_term_0);
    result_by_q.set(c_Key1(-2), common_term_0);
    // q = -1
    result_by_lpq.set(c_Key3(7, 0, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    result_by_lpq.set(c_Key3(7, 0, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(7, 0, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // q = 2
    result_by_lpq.set(c_Key3(7, 0, 2), common_term_4);
    result_by_q.set(c_Key1(2), common_term_4);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 0), result_by_q);
    result_by_q.clear();

    // l , p = (7, 1).
    // q = -2
    result_by_lpq.set(c_Key3(7, 1, -2), common_term_5);
    result_by_q.set(c_Key1(-2), common_term_5);
    // q = -1
    result_by_lpq.set(c_Key3(7, 1, -1), common_term_6);
    result_by_q.set(c_Key1(-1), common_term_6);
    // q = 0
    result_by_lpq.set(c_Key3(7, 1, 0), common_term_7);
    result_by_q.set(c_Key1(0), common_term_7);
    // q = 1
    result_by_lpq.set(c_Key3(7, 1, 1), common_term_8);
    result_by_q.set(c_Key1(1), common_term_8);
    // q = 2
    result_by_lpq.set(c_Key3(7, 1, 2), common_term_9);
    result_by_q.set(c_Key1(2), common_term_9);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 1), result_by_q);
    result_by_q.clear();

    // l , p = (7, 2).
    // q = -2
    result_by_lpq.set(c_Key3(7, 2, -2), common_term_10);
    result_by_q.set(c_Key1(-2), common_term_10);
    // q = -1
    result_by_lpq.set(c_Key3(7, 2, -1), eccentricity);
    result_by_q.set(c_Key1(-1), eccentricity);
    // q = 0
    result_by_lpq.set(c_Key3(7, 2, 0), common_term_11);
    result_by_q.set(c_Key1(0), common_term_11);
    // q = 1
    result_by_lpq.set(c_Key3(7, 2, 1), common_term_12);
    result_by_q.set(c_Key1(1), common_term_12);
    // q = 2
    result_by_lpq.set(c_Key3(7, 2, 2), common_term_13);
    result_by_q.set(c_Key1(2), common_term_13);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 2), result_by_q);
    result_by_q.clear();

    // l , p = (7, 3).
    // q = -2
    result_by_lpq.set(c_Key3(7, 3, -2), common_term_14);
    result_by_q.set(c_Key1(-2), common_term_14);
    // q = -1
    result_by_lpq.set(c_Key3(7, 3, -1), common_term_15);
    result_by_q.set(c_Key1(-1), common_term_15);
    // q = 0
    result_by_lpq.set(c_Key3(7, 3, 0), common_term_16);
    result_by_q.set(c_Key1(0), common_term_16);
    // q = 1
    result_by_lpq.set(c_Key3(7, 3, 1), common_term_17);
    result_by_q.set(c_Key1(1), common_term_17);
    // q = 2
    result_by_lpq.set(c_Key3(7, 3, 2), common_term_18);
    result_by_q.set(c_Key1(2), common_term_18);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 3), result_by_q);
    result_by_q.clear();

    // l , p = (7, 4).
    // q = -2
    result_by_lpq.set(c_Key3(7, 4, -2), common_term_18);
    result_by_q.set(c_Key1(-2), common_term_18);
    // q = -1
    result_by_lpq.set(c_Key3(7, 4, -1), common_term_17);
    result_by_q.set(c_Key1(-1), common_term_17);
    // q = 0
    result_by_lpq.set(c_Key3(7, 4, 0), common_term_16);
    result_by_q.set(c_Key1(0), common_term_16);
    // q = 1
    result_by_lpq.set(c_Key3(7, 4, 1), common_term_15);
    result_by_q.set(c_Key1(1), common_term_15);
    // q = 2
    result_by_lpq.set(c_Key3(7, 4, 2), common_term_14);
    result_by_q.set(c_Key1(2), common_term_14);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 4), result_by_q);
    result_by_q.clear();

    // l , p = (7, 5).
    // q = -2
    result_by_lpq.set(c_Key3(7, 5, -2), common_term_13);
    result_by_q.set(c_Key1(-2), common_term_13);
    // q = -1
    result_by_lpq.set(c_Key3(7, 5, -1), common_term_12);
    result_by_q.set(c_Key1(-1), common_term_12);
    // q = 0
    result_by_lpq.set(c_Key3(7, 5, 0), common_term_11);
    result_by_q.set(c_Key1(0), common_term_11);
    // q = 1
    result_by_lpq.set(c_Key3(7, 5, 1), eccentricity);
    result_by_q.set(c_Key1(1), eccentricity);
    // q = 2
    result_by_lpq.set(c_Key3(7, 5, 2), common_term_10);
    result_by_q.set(c_Key1(2), common_term_10);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 5), result_by_q);
    result_by_q.clear();

    // l , p = (7, 6).
    // q = -2
    result_by_lpq.set(c_Key3(7, 6, -2), common_term_9);
    result_by_q.set(c_Key1(-2), common_term_9);
    // q = -1
    result_by_lpq.set(c_Key3(7, 6, -1), common_term_8);
    result_by_q.set(c_Key1(-1), common_term_8);
    // q = 0
    result_by_lpq.set(c_Key3(7, 6, 0), common_term_7);
    result_by_q.set(c_Key1(0), common_term_7);
    // q = 1
    result_by_lpq.set(c_Key3(7, 6, 1), common_term_6);
    result_by_q.set(c_Key1(1), common_term_6);
    // q = 2
    result_by_lpq.set(c_Key3(7, 6, 2), common_term_5);
    result_by_q.set(c_Key1(2), common_term_5);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 6), result_by_q);
    result_by_q.clear();

    // l , p = (7, 7).
    // q = -2
    result_by_lpq.set(c_Key3(7, 7, -2), common_term_4);
    result_by_q.set(c_Key1(-2), common_term_4);
    // q = -1
    result_by_lpq.set(c_Key3(7, 7, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    result_by_lpq.set(c_Key3(7, 7, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(7, 7, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // q = 2
    result_by_lpq.set(c_Key3(7, 7, 2), common_term_0);
    result_by_q.set(c_Key1(2), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 7), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l7_e3(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 7.
    // Functions are taylor expanded and truncated at eccentricity^3.

    // Functions Calculated for l = 7.

    c_IntMap<c_Key3, double> result_by_lpq(56);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(8);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double common_term_0 = -1.3333333333333333*eccentricity_3;
    double common_term_1 = 3.125*eccentricity_2;
    double common_term_2 = 39.75*eccentricity_3 - 3.0*eccentricity;
    double common_term_3 = 1.0 - 35.0*eccentricity_2;
    double common_term_4 = -228.0*eccentricity_3 + 11.0*eccentricity;
    double common_term_5 = 67.875*eccentricity_2;
    double common_term_6 = 309.58333333333333*eccentricity_3;
    double common_term_7 = 0.083333333333333333*eccentricity_3;
    double common_term_8 = 0.375*eccentricity_2;
    double common_term_9 = 4.5*eccentricity_3 - eccentricity;
    double common_term_10 = 1.0 - 11.0*eccentricity_2;
    double common_term_11 = -66.75*eccentricity_3 + 9.0*eccentricity;
    double common_term_12 = 46.625*eccentricity_2;
    double common_term_13 = 182.16666666666667*eccentricity_3;
    double common_term_14 = 2.5*eccentricity_3*std::pow(1.0 - eccentricity_2, -6.5);
    double common_term_15 = 1.625*eccentricity_2;
    double common_term_16 = 8.25*eccentricity_3 + eccentricity;
    double common_term_17 = 5.0*eccentricity_2 + 1.0;
    double common_term_18 = 13.5*eccentricity_3 + 7.0*eccentricity;
    double common_term_19 = 29.375*eccentricity_2;
    double common_term_20 = 95.75*eccentricity_3;
    double common_term_21 = 13.916666666666667*eccentricity_3;
    double common_term_22 = 6.875*eccentricity_2;
    double common_term_23 = std::pow(1.0 - eccentricity_2, -6.5)*(7.5*eccentricity_3 + 3.0*eccentricity);
    double common_term_24 = 13.0*eccentricity_2 + 1.0;
    double common_term_25 = 36.75*eccentricity_3 + 5.0*eccentricity;
    double common_term_26 = 16.125*eccentricity_2;
    double common_term_27 = 42.333333333333333*eccentricity_3;

    c_IntMap<c_Key1, double> result_by_q(7);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (7, 0).
    // q = -3
    result_by_lpq.set(c_Key3(7, 0, -3), common_term_0);
    result_by_q.set(c_Key1(-3), common_term_0);
    // q = -2
    result_by_lpq.set(c_Key3(7, 0, -2), common_term_1);
    result_by_q.set(c_Key1(-2), common_term_1);
    // q = -1
    result_by_lpq.set(c_Key3(7, 0, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    result_by_lpq.set(c_Key3(7, 0, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(7, 0, 1), common_term_4);
    result_by_q.set(c_Key1(1), common_term_4);
    // q = 2
    result_by_lpq.set(c_Key3(7, 0, 2), common_term_5);
    result_by_q.set(c_Key1(2), common_term_5);
    // q = 3
    result_by_lpq.set(c_Key3(7, 0, 3), common_term_6);
    result_by_q.set(c_Key1(3), common_term_6);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 0), result_by_q);
    result_by_q.clear();

    // l , p = (7, 1).
    // q = -3
    result_by_lpq.set(c_Key3(7, 1, -3), common_term_7);
    result_by_q.set(c_Key1(-3), common_term_7);
    // q = -2
    result_by_lpq.set(c_Key3(7, 1, -2), common_term_8);
    result_by_q.set(c_Key1(-2), common_term_8);
    // q = -1
    result_by_lpq.set(c_Key3(7, 1, -1), common_term_9);
    result_by_q.set(c_Key1(-1), common_term_9);
    // q = 0
    result_by_lpq.set(c_Key3(7, 1, 0), common_term_10);
    result_by_q.set(c_Key1(0), common_term_10);
    // q = 1
    result_by_lpq.set(c_Key3(7, 1, 1), common_term_11);
    result_by_q.set(c_Key1(1), common_term_11);
    // q = 2
    result_by_lpq.set(c_Key3(7, 1, 2), common_term_12);
    result_by_q.set(c_Key1(2), common_term_12);
    // q = 3
    result_by_lpq.set(c_Key3(7, 1, 3), common_term_13);
    result_by_q.set(c_Key1(3), common_term_13);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 1), result_by_q);
    result_by_q.clear();

    // l , p = (7, 2).
    // q = -3
    result_by_lpq.set(c_Key3(7, 2, -3), common_term_14);
    result_by_q.set(c_Key1(-3), common_term_14);
    // q = -2
    result_by_lpq.set(c_Key3(7, 2, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(7, 2, -1), common_term_16);
    result_by_q.set(c_Key1(-1), common_term_16);
    // q = 0
    result_by_lpq.set(c_Key3(7, 2, 0), common_term_17);
    result_by_q.set(c_Key1(0), common_term_17);
    // q = 1
    result_by_lpq.set(c_Key3(7, 2, 1), common_term_18);
    result_by_q.set(c_Key1(1), common_term_18);
    // q = 2
    result_by_lpq.set(c_Key3(7, 2, 2), common_term_19);
    result_by_q.set(c_Key1(2), common_term_19);
    // q = 3
    result_by_lpq.set(c_Key3(7, 2, 3), common_term_20);
    result_by_q.set(c_Key1(3), common_term_20);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 2), result_by_q);
    result_by_q.clear();

    // l , p = (7, 3).
    // q = -3
    result_by_lpq.set(c_Key3(7, 3, -3), common_term_21);
    result_by_q.set(c_Key1(-3), common_term_21);
    // q = -2
    result_by_lpq.set(c_Key3(7, 3, -2), common_term_22);
    result_by_q.set(c_Key1(-2), common_term_22);
    // q = -1
    result_by_lpq.set(c_Key3(7, 3, -1), common_term_23);
    result_by_q.set(c_Key1(-1), common_term_23);
    // q = 0
    result_by_lpq.set(c_Key3(7, 3, 0), common_term_24);
    result_by_q.set(c_Key1(0), common_term_24);
    // q = 1
    result_by_lpq.set(c_Key3(7, 3, 1), common_term_25);
    result_by_q.set(c_Key1(1), common_term_25);
    // q = 2
    result_by_lpq.set(c_Key3(7, 3, 2), common_term_26);
    result_by_q.set(c_Key1(2), common_term_26);
    // q = 3
    result_by_lpq.set(c_Key3(7, 3, 3), common_term_27);
    result_by_q.set(c_Key1(3), common_term_27);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 3), result_by_q);
    result_by_q.clear();

    // l , p = (7, 4).
    // q = -3
    result_by_lpq.set(c_Key3(7, 4, -3), common_term_27);
    result_by_q.set(c_Key1(-3), common_term_27);
    // q = -2
    result_by_lpq.set(c_Key3(7, 4, -2), common_term_26);
    result_by_q.set(c_Key1(-2), common_term_26);
    // q = -1
    result_by_lpq.set(c_Key3(7, 4, -1), common_term_25);
    result_by_q.set(c_Key1(-1), common_term_25);
    // q = 0
    result_by_lpq.set(c_Key3(7, 4, 0), common_term_24);
    result_by_q.set(c_Key1(0), common_term_24);
    // q = 1
    result_by_lpq.set(c_Key3(7, 4, 1), common_term_23);
    result_by_q.set(c_Key1(1), common_term_23);
    // q = 2
    result_by_lpq.set(c_Key3(7, 4, 2), common_term_22);
    result_by_q.set(c_Key1(2), common_term_22);
    // q = 3
    result_by_lpq.set(c_Key3(7, 4, 3), common_term_21);
    result_by_q.set(c_Key1(3), common_term_21);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 4), result_by_q);
    result_by_q.clear();

    // l , p = (7, 5).
    // q = -3
    result_by_lpq.set(c_Key3(7, 5, -3), common_term_20);
    result_by_q.set(c_Key1(-3), common_term_20);
    // q = -2
    result_by_lpq.set(c_Key3(7, 5, -2), common_term_19);
    result_by_q.set(c_Key1(-2), common_term_19);
    // q = -1
    result_by_lpq.set(c_Key3(7, 5, -1), common_term_18);
    result_by_q.set(c_Key1(-1), common_term_18);
    // q = 0
    result_by_lpq.set(c_Key3(7, 5, 0), common_term_17);
    result_by_q.set(c_Key1(0), common_term_17);
    // q = 1
    result_by_lpq.set(c_Key3(7, 5, 1), common_term_16);
    result_by_q.set(c_Key1(1), common_term_16);
    // q = 2
    result_by_lpq.set(c_Key3(7, 5, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // q = 3
    result_by_lpq.set(c_Key3(7, 5, 3), common_term_14);
    result_by_q.set(c_Key1(3), common_term_14);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 5), result_by_q);
    result_by_q.clear();

    // l , p = (7, 6).
    // q = -3
    result_by_lpq.set(c_Key3(7, 6, -3), common_term_13);
    result_by_q.set(c_Key1(-3), common_term_13);
    // q = -2
    result_by_lpq.set(c_Key3(7, 6, -2), common_term_12);
    result_by_q.set(c_Key1(-2), common_term_12);
    // q = -1
    result_by_lpq.set(c_Key3(7, 6, -1), common_term_11);
    result_by_q.set(c_Key1(-1), common_term_11);
    // q = 0
    result_by_lpq.set(c_Key3(7, 6, 0), common_term_10);
    result_by_q.set(c_Key1(0), common_term_10);
    // q = 1
    result_by_lpq.set(c_Key3(7, 6, 1), common_term_9);
    result_by_q.set(c_Key1(1), common_term_9);
    // q = 2
    result_by_lpq.set(c_Key3(7, 6, 2), common_term_8);
    result_by_q.set(c_Key1(2), common_term_8);
    // q = 3
    result_by_lpq.set(c_Key3(7, 6, 3), common_term_7);
    result_by_q.set(c_Key1(3), common_term_7);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 6), result_by_q);
    result_by_q.clear();

    // l , p = (7, 7).
    // q = -3
    result_by_lpq.set(c_Key3(7, 7, -3), common_term_6);
    result_by_q.set(c_Key1(-3), common_term_6);
    // q = -2
    result_by_lpq.set(c_Key3(7, 7, -2), common_term_5);
    result_by_q.set(c_Key1(-2), common_term_5);
    // q = -1
    result_by_lpq.set(c_Key3(7, 7, -1), common_term_4);
    result_by_q.set(c_Key1(-1), common_term_4);
    // q = 0
    result_by_lpq.set(c_Key3(7, 7, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(7, 7, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // q = 2
    result_by_lpq.set(c_Key3(7, 7, 2), common_term_1);
    result_by_q.set(c_Key1(2), common_term_1);
    // q = 3
    result_by_lpq.set(c_Key3(7, 7, 3), common_term_0);
    result_by_q.set(c_Key1(3), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 7), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l7_e4(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 7.
    // Functions are taylor expanded and truncated at eccentricity^4.

    // Functions Calculated for l = 7.

    c_IntMap<c_Key3, double> result_by_lpq(72);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(8);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double common_term_0 = 0.2109375*eccentricity_4;
    double common_term_1 = -1.3333333333333333*eccentricity_3;
    double common_term_2 = -19.270833333333333*eccentricity_4 + 3.125*eccentricity_2;
    double common_term_3 = 39.75*eccentricity_3 - 3.0*eccentricity;
    double common_term_4 = 280.109375*eccentricity_4 - 35.0*eccentricity_2 + 1.0;
    double common_term_5 = -228.0*eccentricity_3 + 11.0*eccentricity;
    double common_term_6 = -1094.0625*eccentricity_4 + 67.875*eccentricity_2;
    double common_term_7 = 309.58333333333333*eccentricity_3;
    double common_term_8 = 1163.0130208333333*eccentricity_4;
    double common_term_9 = 0.13802083333333333*eccentricity_4;
    double common_term_10 = 0.083333333333333333*eccentricity_3;
    double common_term_11 = 0.1875*eccentricity_4 + 0.375*eccentricity_2;
    double common_term_12 = 4.5*eccentricity_3 - eccentricity;
    double common_term_13 = 31.484375*eccentricity_4 - 11.0*eccentricity_2 + 1.0;
    double common_term_14 = -66.75*eccentricity_3 + 9.0*eccentricity;
    double common_term_15 = -297.52083333333333*eccentricity_4 + 46.625*eccentricity_2;
    double common_term_16 = 182.16666666666667*eccentricity_3;
    double common_term_17 = 595.7109375*eccentricity_4;
    double common_term_18 = 3.7942708333333333*eccentricity_4;
    double common_term_19 = 2.5*eccentricity_3*std::pow(1.0 - eccentricity_2, -6.5);
    double common_term_20 = 11.979166666666667*eccentricity_4 + 1.625*eccentricity_2;
    double common_term_21 = 8.25*eccentricity_3 + eccentricity;
    double common_term_22 = 25.734375*eccentricity_4 + 5.0*eccentricity_2 + 1.0;
    double common_term_23 = 13.5*eccentricity_3 + 7.0*eccentricity;
    double common_term_24 = 21.354166666666667*eccentricity_4 + 29.375*eccentricity_2;
    double common_term_25 = 95.75*eccentricity_3;
    double common_term_26 = 267.13802083333333*eccentricity_4;
    double common_term_27 = 26.1796875*eccentricity_4;
    double common_term_28 = 13.916666666666667*eccentricity_3;
    double common_term_29 = 50.104166666666667*eccentricity_4 + 6.875*eccentricity_2;
    double common_term_30 = std::pow(1.0 - eccentricity_2, -6.5)*(7.5*eccentricity_3 + 3.0*eccentricity);
    double common_term_31 = 70.859375*eccentricity_4 + 13.0*eccentricity_2 + 1.0;
    double common_term_32 = 36.75*eccentricity_3 + 5.0*eccentricity;
    double common_term_33 = 84.5625*eccentricity_4 + 16.125*eccentricity_2;
    double common_term_34 = 42.333333333333333*eccentricity_3;
    double common_term_35 = 98.294270833333333*eccentricity_4;

    c_IntMap<c_Key1, double> result_by_q(9);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (7, 0).
    // q = -4
    result_by_lpq.set(c_Key3(7, 0, -4), common_term_0);
    result_by_q.set(c_Key1(-4), common_term_0);
    // q = -3
    result_by_lpq.set(c_Key3(7, 0, -3), common_term_1);
    result_by_q.set(c_Key1(-3), common_term_1);
    // q = -2
    result_by_lpq.set(c_Key3(7, 0, -2), common_term_2);
    result_by_q.set(c_Key1(-2), common_term_2);
    // q = -1
    result_by_lpq.set(c_Key3(7, 0, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    result_by_lpq.set(c_Key3(7, 0, 0), common_term_4);
    result_by_q.set(c_Key1(0), common_term_4);
    // q = 1
    result_by_lpq.set(c_Key3(7, 0, 1), common_term_5);
    result_by_q.set(c_Key1(1), common_term_5);
    // q = 2
    result_by_lpq.set(c_Key3(7, 0, 2), common_term_6);
    result_by_q.set(c_Key1(2), common_term_6);
    // q = 3
    result_by_lpq.set(c_Key3(7, 0, 3), common_term_7);
    result_by_q.set(c_Key1(3), common_term_7);
    // q = 4
    result_by_lpq.set(c_Key3(7, 0, 4), common_term_8);
    result_by_q.set(c_Key1(4), common_term_8);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 0), result_by_q);
    result_by_q.clear();

    // l , p = (7, 1).
    // q = -4
    result_by_lpq.set(c_Key3(7, 1, -4), common_term_9);
    result_by_q.set(c_Key1(-4), common_term_9);
    // q = -3
    result_by_lpq.set(c_Key3(7, 1, -3), common_term_10);
    result_by_q.set(c_Key1(-3), common_term_10);
    // q = -2
    result_by_lpq.set(c_Key3(7, 1, -2), common_term_11);
    result_by_q.set(c_Key1(-2), common_term_11);
    // q = -1
    result_by_lpq.set(c_Key3(7, 1, -1), common_term_12);
    result_by_q.set(c_Key1(-1), common_term_12);
    // q = 0
    result_by_lpq.set(c_Key3(7, 1, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(7, 1, 1), common_term_14);
    result_by_q.set(c_Key1(1), common_term_14);
    // q = 2
    result_by_lpq.set(c_Key3(7, 1, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // q = 3
    result_by_lpq.set(c_Key3(7, 1, 3), common_term_16);
    result_by_q.set(c_Key1(3), common_term_16);
    // q = 4
    result_by_lpq.set(c_Key3(7, 1, 4), common_term_17);
    result_by_q.set(c_Key1(4), common_term_17);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 1), result_by_q);
    result_by_q.clear();

    // l , p = (7, 2).
    // q = -4
    result_by_lpq.set(c_Key3(7, 2, -4), common_term_18);
    result_by_q.set(c_Key1(-4), common_term_18);
    // q = -3
    result_by_lpq.set(c_Key3(7, 2, -3), common_term_19);
    result_by_q.set(c_Key1(-3), common_term_19);
    // q = -2
    result_by_lpq.set(c_Key3(7, 2, -2), common_term_20);
    result_by_q.set(c_Key1(-2), common_term_20);
    // q = -1
    result_by_lpq.set(c_Key3(7, 2, -1), common_term_21);
    result_by_q.set(c_Key1(-1), common_term_21);
    // q = 0
    result_by_lpq.set(c_Key3(7, 2, 0), common_term_22);
    result_by_q.set(c_Key1(0), common_term_22);
    // q = 1
    result_by_lpq.set(c_Key3(7, 2, 1), common_term_23);
    result_by_q.set(c_Key1(1), common_term_23);
    // q = 2
    result_by_lpq.set(c_Key3(7, 2, 2), common_term_24);
    result_by_q.set(c_Key1(2), common_term_24);
    // q = 3
    result_by_lpq.set(c_Key3(7, 2, 3), common_term_25);
    result_by_q.set(c_Key1(3), common_term_25);
    // q = 4
    result_by_lpq.set(c_Key3(7, 2, 4), common_term_26);
    result_by_q.set(c_Key1(4), common_term_26);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 2), result_by_q);
    result_by_q.clear();

    // l , p = (7, 3).
    // q = -4
    result_by_lpq.set(c_Key3(7, 3, -4), common_term_27);
    result_by_q.set(c_Key1(-4), common_term_27);
    // q = -3
    result_by_lpq.set(c_Key3(7, 3, -3), common_term_28);
    result_by_q.set(c_Key1(-3), common_term_28);
    // q = -2
    result_by_lpq.set(c_Key3(7, 3, -2), common_term_29);
    result_by_q.set(c_Key1(-2), common_term_29);
    // q = -1
    result_by_lpq.set(c_Key3(7, 3, -1), common_term_30);
    result_by_q.set(c_Key1(-1), common_term_30);
    // q = 0
    result_by_lpq.set(c_Key3(7, 3, 0), common_term_31);
    result_by_q.set(c_Key1(0), common_term_31);
    // q = 1
    result_by_lpq.set(c_Key3(7, 3, 1), common_term_32);
    result_by_q.set(c_Key1(1), common_term_32);
    // q = 2
    result_by_lpq.set(c_Key3(7, 3, 2), common_term_33);
    result_by_q.set(c_Key1(2), common_term_33);
    // q = 3
    result_by_lpq.set(c_Key3(7, 3, 3), common_term_34);
    result_by_q.set(c_Key1(3), common_term_34);
    // q = 4
    result_by_lpq.set(c_Key3(7, 3, 4), common_term_35);
    result_by_q.set(c_Key1(4), common_term_35);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 3), result_by_q);
    result_by_q.clear();

    // l , p = (7, 4).
    // q = -4
    result_by_lpq.set(c_Key3(7, 4, -4), common_term_35);
    result_by_q.set(c_Key1(-4), common_term_35);
    // q = -3
    result_by_lpq.set(c_Key3(7, 4, -3), common_term_34);
    result_by_q.set(c_Key1(-3), common_term_34);
    // q = -2
    result_by_lpq.set(c_Key3(7, 4, -2), common_term_33);
    result_by_q.set(c_Key1(-2), common_term_33);
    // q = -1
    result_by_lpq.set(c_Key3(7, 4, -1), common_term_32);
    result_by_q.set(c_Key1(-1), common_term_32);
    // q = 0
    result_by_lpq.set(c_Key3(7, 4, 0), common_term_31);
    result_by_q.set(c_Key1(0), common_term_31);
    // q = 1
    result_by_lpq.set(c_Key3(7, 4, 1), common_term_30);
    result_by_q.set(c_Key1(1), common_term_30);
    // q = 2
    result_by_lpq.set(c_Key3(7, 4, 2), common_term_29);
    result_by_q.set(c_Key1(2), common_term_29);
    // q = 3
    result_by_lpq.set(c_Key3(7, 4, 3), common_term_28);
    result_by_q.set(c_Key1(3), common_term_28);
    // q = 4
    result_by_lpq.set(c_Key3(7, 4, 4), common_term_27);
    result_by_q.set(c_Key1(4), common_term_27);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 4), result_by_q);
    result_by_q.clear();

    // l , p = (7, 5).
    // q = -4
    result_by_lpq.set(c_Key3(7, 5, -4), common_term_26);
    result_by_q.set(c_Key1(-4), common_term_26);
    // q = -3
    result_by_lpq.set(c_Key3(7, 5, -3), common_term_25);
    result_by_q.set(c_Key1(-3), common_term_25);
    // q = -2
    result_by_lpq.set(c_Key3(7, 5, -2), common_term_24);
    result_by_q.set(c_Key1(-2), common_term_24);
    // q = -1
    result_by_lpq.set(c_Key3(7, 5, -1), common_term_23);
    result_by_q.set(c_Key1(-1), common_term_23);
    // q = 0
    result_by_lpq.set(c_Key3(7, 5, 0), common_term_22);
    result_by_q.set(c_Key1(0), common_term_22);
    // q = 1
    result_by_lpq.set(c_Key3(7, 5, 1), common_term_21);
    result_by_q.set(c_Key1(1), common_term_21);
    // q = 2
    result_by_lpq.set(c_Key3(7, 5, 2), common_term_20);
    result_by_q.set(c_Key1(2), common_term_20);
    // q = 3
    result_by_lpq.set(c_Key3(7, 5, 3), common_term_19);
    result_by_q.set(c_Key1(3), common_term_19);
    // q = 4
    result_by_lpq.set(c_Key3(7, 5, 4), common_term_18);
    result_by_q.set(c_Key1(4), common_term_18);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 5), result_by_q);
    result_by_q.clear();

    // l , p = (7, 6).
    // q = -4
    result_by_lpq.set(c_Key3(7, 6, -4), common_term_17);
    result_by_q.set(c_Key1(-4), common_term_17);
    // q = -3
    result_by_lpq.set(c_Key3(7, 6, -3), common_term_16);
    result_by_q.set(c_Key1(-3), common_term_16);
    // q = -2
    result_by_lpq.set(c_Key3(7, 6, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(7, 6, -1), common_term_14);
    result_by_q.set(c_Key1(-1), common_term_14);
    // q = 0
    result_by_lpq.set(c_Key3(7, 6, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(7, 6, 1), common_term_12);
    result_by_q.set(c_Key1(1), common_term_12);
    // q = 2
    result_by_lpq.set(c_Key3(7, 6, 2), common_term_11);
    result_by_q.set(c_Key1(2), common_term_11);
    // q = 3
    result_by_lpq.set(c_Key3(7, 6, 3), common_term_10);
    result_by_q.set(c_Key1(3), common_term_10);
    // q = 4
    result_by_lpq.set(c_Key3(7, 6, 4), common_term_9);
    result_by_q.set(c_Key1(4), common_term_9);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 6), result_by_q);
    result_by_q.clear();

    // l , p = (7, 7).
    // q = -4
    result_by_lpq.set(c_Key3(7, 7, -4), common_term_8);
    result_by_q.set(c_Key1(-4), common_term_8);
    // q = -3
    result_by_lpq.set(c_Key3(7, 7, -3), common_term_7);
    result_by_q.set(c_Key1(-3), common_term_7);
    // q = -2
    result_by_lpq.set(c_Key3(7, 7, -2), common_term_6);
    result_by_q.set(c_Key1(-2), common_term_6);
    // q = -1
    result_by_lpq.set(c_Key3(7, 7, -1), common_term_5);
    result_by_q.set(c_Key1(-1), common_term_5);
    // q = 0
    result_by_lpq.set(c_Key3(7, 7, 0), common_term_4);
    result_by_q.set(c_Key1(0), common_term_4);
    // q = 1
    result_by_lpq.set(c_Key3(7, 7, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // q = 2
    result_by_lpq.set(c_Key3(7, 7, 2), common_term_2);
    result_by_q.set(c_Key1(2), common_term_2);
    // q = 3
    result_by_lpq.set(c_Key3(7, 7, 3), common_term_1);
    result_by_q.set(c_Key1(3), common_term_1);
    // q = 4
    result_by_lpq.set(c_Key3(7, 7, 4), common_term_0);
    result_by_q.set(c_Key1(4), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 7), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l7_e5(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 7.
    // Functions are taylor expanded and truncated at eccentricity^5.

    // Functions Calculated for l = 7.

    c_IntMap<c_Key3, double> result_by_lpq(88);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(8);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_5 = eccentricity_4 * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double common_term_0 = -0.0083333333333333333*eccentricity_5;
    double common_term_1 = 0.2109375*eccentricity_4;
    double common_term_2 = 3.6666666666666667*eccentricity_5 - 1.3333333333333333*eccentricity_3;
    double common_term_3 = -19.270833333333333*eccentricity_4 + 3.125*eccentricity_2;
    double common_term_4 = -152.8125*eccentricity_5 + 39.75*eccentricity_3 - 3.0*eccentricity;
    double common_term_5 = 280.109375*eccentricity_4 - 35.0*eccentricity_2 + 1.0;
    double common_term_6 = 1445.0833333333333*eccentricity_5 - 228.0*eccentricity_3 + 11.0*eccentricity;
    double common_term_7 = -1094.0625*eccentricity_4 + 67.875*eccentricity_2;
    double common_term_8 = -4311.9791666666667*eccentricity_5 + 309.58333333333333*eccentricity_3;
    double common_term_9 = 1163.0130208333333*eccentricity_4;
    double common_term_10 = 3808.05*eccentricity_5;
    double common_term_11 = 0.1875*eccentricity_5*std::pow(1.0 - eccentricity_2, -6.5);
    double common_term_12 = 0.13802083333333333*eccentricity_4;
    double common_term_13 = 0.70833333333333333*eccentricity_5 + 0.083333333333333333*eccentricity_3;
    double common_term_14 = 0.1875*eccentricity_4 + 0.375*eccentricity_2;
    double common_term_15 = -2.1666666666666667*eccentricity_5 + 4.5*eccentricity_3 - eccentricity;
    double common_term_16 = 31.484375*eccentricity_4 - 11.0*eccentricity_2 + 1.0;
    double common_term_17 = 162.9375*eccentricity_5 - 66.75*eccentricity_3 + 9.0*eccentricity;
    double common_term_18 = -297.52083333333333*eccentricity_4 + 46.625*eccentricity_2;
    double common_term_19 = -1090.0208333333333*eccentricity_5 + 182.16666666666667*eccentricity_3;
    double common_term_20 = 595.7109375*eccentricity_4;
    double common_term_21 = 1720.3541666666667*eccentricity_5;
    double common_term_22 = 5.6958333333333333*eccentricity_5;
    double common_term_23 = 3.7942708333333333*eccentricity_4;
    double common_term_24 = std::pow(1.0 - eccentricity_2, -6.5)*(0.9375*eccentricity_5 + 2.5*eccentricity_3);
    double common_term_25 = 11.979166666666667*eccentricity_4 + 1.625*eccentricity_2;
    double common_term_26 = 35.604166666666667*eccentricity_5 + 8.25*eccentricity_3 + eccentricity;
    double common_term_27 = 25.734375*eccentricity_4 + 5.0*eccentricity_2 + 1.0;
    double common_term_28 = 64.166666666666667*eccentricity_5 + 13.5*eccentricity_3 + 7.0*eccentricity;
    double common_term_29 = 21.354166666666667*eccentricity_4 + 29.375*eccentricity_2;
    double common_term_30 = -0.5*eccentricity_5 + 95.75*eccentricity_3;
    double common_term_31 = 267.13802083333333*eccentricity_4;
    double common_term_32 = 669.84583333333333*eccentricity_5;
    double common_term_33 = 46.891666666666667*eccentricity_5;
    double common_term_34 = 26.1796875*eccentricity_4;
    double common_term_35 = 86.729166666666667*eccentricity_5 + 13.916666666666667*eccentricity_3;
    double common_term_36 = 50.104166666666667*eccentricity_4 + 6.875*eccentricity_2;
    double common_term_37 = std::pow(1.0 - eccentricity_2, -6.5)*(1.875*eccentricity_5 + 7.5*eccentricity_3 + 3.0*eccentricity);
    double common_term_38 = 70.859375*eccentricity_4 + 13.0*eccentricity_2 + 1.0;
    double common_term_39 = 153.52083333333333*eccentricity_5 + 36.75*eccentricity_3 + 5.0*eccentricity;
    double common_term_40 = 84.5625*eccentricity_4 + 16.125*eccentricity_2;
    double common_term_41 = 170.95833333333333*eccentricity_5 + 42.333333333333333*eccentricity_3;
    double common_term_42 = 98.294270833333333*eccentricity_4;
    double common_term_43 = 210.15*eccentricity_5;

    c_IntMap<c_Key1, double> result_by_q(11);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (7, 0).
    // q = -5
    result_by_lpq.set(c_Key3(7, 0, -5), common_term_0);
    result_by_q.set(c_Key1(-5), common_term_0);
    // q = -4
    result_by_lpq.set(c_Key3(7, 0, -4), common_term_1);
    result_by_q.set(c_Key1(-4), common_term_1);
    // q = -3
    result_by_lpq.set(c_Key3(7, 0, -3), common_term_2);
    result_by_q.set(c_Key1(-3), common_term_2);
    // q = -2
    result_by_lpq.set(c_Key3(7, 0, -2), common_term_3);
    result_by_q.set(c_Key1(-2), common_term_3);
    // q = -1
    result_by_lpq.set(c_Key3(7, 0, -1), common_term_4);
    result_by_q.set(c_Key1(-1), common_term_4);
    // q = 0
    result_by_lpq.set(c_Key3(7, 0, 0), common_term_5);
    result_by_q.set(c_Key1(0), common_term_5);
    // q = 1
    result_by_lpq.set(c_Key3(7, 0, 1), common_term_6);
    result_by_q.set(c_Key1(1), common_term_6);
    // q = 2
    result_by_lpq.set(c_Key3(7, 0, 2), common_term_7);
    result_by_q.set(c_Key1(2), common_term_7);
    // q = 3
    result_by_lpq.set(c_Key3(7, 0, 3), common_term_8);
    result_by_q.set(c_Key1(3), common_term_8);
    // q = 4
    result_by_lpq.set(c_Key3(7, 0, 4), common_term_9);
    result_by_q.set(c_Key1(4), common_term_9);
    // q = 5
    result_by_lpq.set(c_Key3(7, 0, 5), common_term_10);
    result_by_q.set(c_Key1(5), common_term_10);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 0), result_by_q);
    result_by_q.clear();

    // l , p = (7, 1).
    // q = -5
    result_by_lpq.set(c_Key3(7, 1, -5), common_term_11);
    result_by_q.set(c_Key1(-5), common_term_11);
    // q = -4
    result_by_lpq.set(c_Key3(7, 1, -4), common_term_12);
    result_by_q.set(c_Key1(-4), common_term_12);
    // q = -3
    result_by_lpq.set(c_Key3(7, 1, -3), common_term_13);
    result_by_q.set(c_Key1(-3), common_term_13);
    // q = -2
    result_by_lpq.set(c_Key3(7, 1, -2), common_term_14);
    result_by_q.set(c_Key1(-2), common_term_14);
    // q = -1
    result_by_lpq.set(c_Key3(7, 1, -1), common_term_15);
    result_by_q.set(c_Key1(-1), common_term_15);
    // q = 0
    result_by_lpq.set(c_Key3(7, 1, 0), common_term_16);
    result_by_q.set(c_Key1(0), common_term_16);
    // q = 1
    result_by_lpq.set(c_Key3(7, 1, 1), common_term_17);
    result_by_q.set(c_Key1(1), common_term_17);
    // q = 2
    result_by_lpq.set(c_Key3(7, 1, 2), common_term_18);
    result_by_q.set(c_Key1(2), common_term_18);
    // q = 3
    result_by_lpq.set(c_Key3(7, 1, 3), common_term_19);
    result_by_q.set(c_Key1(3), common_term_19);
    // q = 4
    result_by_lpq.set(c_Key3(7, 1, 4), common_term_20);
    result_by_q.set(c_Key1(4), common_term_20);
    // q = 5
    result_by_lpq.set(c_Key3(7, 1, 5), common_term_21);
    result_by_q.set(c_Key1(5), common_term_21);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 1), result_by_q);
    result_by_q.clear();

    // l , p = (7, 2).
    // q = -5
    result_by_lpq.set(c_Key3(7, 2, -5), common_term_22);
    result_by_q.set(c_Key1(-5), common_term_22);
    // q = -4
    result_by_lpq.set(c_Key3(7, 2, -4), common_term_23);
    result_by_q.set(c_Key1(-4), common_term_23);
    // q = -3
    result_by_lpq.set(c_Key3(7, 2, -3), common_term_24);
    result_by_q.set(c_Key1(-3), common_term_24);
    // q = -2
    result_by_lpq.set(c_Key3(7, 2, -2), common_term_25);
    result_by_q.set(c_Key1(-2), common_term_25);
    // q = -1
    result_by_lpq.set(c_Key3(7, 2, -1), common_term_26);
    result_by_q.set(c_Key1(-1), common_term_26);
    // q = 0
    result_by_lpq.set(c_Key3(7, 2, 0), common_term_27);
    result_by_q.set(c_Key1(0), common_term_27);
    // q = 1
    result_by_lpq.set(c_Key3(7, 2, 1), common_term_28);
    result_by_q.set(c_Key1(1), common_term_28);
    // q = 2
    result_by_lpq.set(c_Key3(7, 2, 2), common_term_29);
    result_by_q.set(c_Key1(2), common_term_29);
    // q = 3
    result_by_lpq.set(c_Key3(7, 2, 3), common_term_30);
    result_by_q.set(c_Key1(3), common_term_30);
    // q = 4
    result_by_lpq.set(c_Key3(7, 2, 4), common_term_31);
    result_by_q.set(c_Key1(4), common_term_31);
    // q = 5
    result_by_lpq.set(c_Key3(7, 2, 5), common_term_32);
    result_by_q.set(c_Key1(5), common_term_32);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 2), result_by_q);
    result_by_q.clear();

    // l , p = (7, 3).
    // q = -5
    result_by_lpq.set(c_Key3(7, 3, -5), common_term_33);
    result_by_q.set(c_Key1(-5), common_term_33);
    // q = -4
    result_by_lpq.set(c_Key3(7, 3, -4), common_term_34);
    result_by_q.set(c_Key1(-4), common_term_34);
    // q = -3
    result_by_lpq.set(c_Key3(7, 3, -3), common_term_35);
    result_by_q.set(c_Key1(-3), common_term_35);
    // q = -2
    result_by_lpq.set(c_Key3(7, 3, -2), common_term_36);
    result_by_q.set(c_Key1(-2), common_term_36);
    // q = -1
    result_by_lpq.set(c_Key3(7, 3, -1), common_term_37);
    result_by_q.set(c_Key1(-1), common_term_37);
    // q = 0
    result_by_lpq.set(c_Key3(7, 3, 0), common_term_38);
    result_by_q.set(c_Key1(0), common_term_38);
    // q = 1
    result_by_lpq.set(c_Key3(7, 3, 1), common_term_39);
    result_by_q.set(c_Key1(1), common_term_39);
    // q = 2
    result_by_lpq.set(c_Key3(7, 3, 2), common_term_40);
    result_by_q.set(c_Key1(2), common_term_40);
    // q = 3
    result_by_lpq.set(c_Key3(7, 3, 3), common_term_41);
    result_by_q.set(c_Key1(3), common_term_41);
    // q = 4
    result_by_lpq.set(c_Key3(7, 3, 4), common_term_42);
    result_by_q.set(c_Key1(4), common_term_42);
    // q = 5
    result_by_lpq.set(c_Key3(7, 3, 5), common_term_43);
    result_by_q.set(c_Key1(5), common_term_43);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 3), result_by_q);
    result_by_q.clear();

    // l , p = (7, 4).
    // q = -5
    result_by_lpq.set(c_Key3(7, 4, -5), common_term_43);
    result_by_q.set(c_Key1(-5), common_term_43);
    // q = -4
    result_by_lpq.set(c_Key3(7, 4, -4), common_term_42);
    result_by_q.set(c_Key1(-4), common_term_42);
    // q = -3
    result_by_lpq.set(c_Key3(7, 4, -3), common_term_41);
    result_by_q.set(c_Key1(-3), common_term_41);
    // q = -2
    result_by_lpq.set(c_Key3(7, 4, -2), common_term_40);
    result_by_q.set(c_Key1(-2), common_term_40);
    // q = -1
    result_by_lpq.set(c_Key3(7, 4, -1), common_term_39);
    result_by_q.set(c_Key1(-1), common_term_39);
    // q = 0
    result_by_lpq.set(c_Key3(7, 4, 0), common_term_38);
    result_by_q.set(c_Key1(0), common_term_38);
    // q = 1
    result_by_lpq.set(c_Key3(7, 4, 1), common_term_37);
    result_by_q.set(c_Key1(1), common_term_37);
    // q = 2
    result_by_lpq.set(c_Key3(7, 4, 2), common_term_36);
    result_by_q.set(c_Key1(2), common_term_36);
    // q = 3
    result_by_lpq.set(c_Key3(7, 4, 3), common_term_35);
    result_by_q.set(c_Key1(3), common_term_35);
    // q = 4
    result_by_lpq.set(c_Key3(7, 4, 4), common_term_34);
    result_by_q.set(c_Key1(4), common_term_34);
    // q = 5
    result_by_lpq.set(c_Key3(7, 4, 5), common_term_33);
    result_by_q.set(c_Key1(5), common_term_33);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 4), result_by_q);
    result_by_q.clear();

    // l , p = (7, 5).
    // q = -5
    result_by_lpq.set(c_Key3(7, 5, -5), common_term_32);
    result_by_q.set(c_Key1(-5), common_term_32);
    // q = -4
    result_by_lpq.set(c_Key3(7, 5, -4), common_term_31);
    result_by_q.set(c_Key1(-4), common_term_31);
    // q = -3
    result_by_lpq.set(c_Key3(7, 5, -3), common_term_30);
    result_by_q.set(c_Key1(-3), common_term_30);
    // q = -2
    result_by_lpq.set(c_Key3(7, 5, -2), common_term_29);
    result_by_q.set(c_Key1(-2), common_term_29);
    // q = -1
    result_by_lpq.set(c_Key3(7, 5, -1), common_term_28);
    result_by_q.set(c_Key1(-1), common_term_28);
    // q = 0
    result_by_lpq.set(c_Key3(7, 5, 0), common_term_27);
    result_by_q.set(c_Key1(0), common_term_27);
    // q = 1
    result_by_lpq.set(c_Key3(7, 5, 1), common_term_26);
    result_by_q.set(c_Key1(1), common_term_26);
    // q = 2
    result_by_lpq.set(c_Key3(7, 5, 2), common_term_25);
    result_by_q.set(c_Key1(2), common_term_25);
    // q = 3
    result_by_lpq.set(c_Key3(7, 5, 3), common_term_24);
    result_by_q.set(c_Key1(3), common_term_24);
    // q = 4
    result_by_lpq.set(c_Key3(7, 5, 4), common_term_23);
    result_by_q.set(c_Key1(4), common_term_23);
    // q = 5
    result_by_lpq.set(c_Key3(7, 5, 5), common_term_22);
    result_by_q.set(c_Key1(5), common_term_22);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 5), result_by_q);
    result_by_q.clear();

    // l , p = (7, 6).
    // q = -5
    result_by_lpq.set(c_Key3(7, 6, -5), common_term_21);
    result_by_q.set(c_Key1(-5), common_term_21);
    // q = -4
    result_by_lpq.set(c_Key3(7, 6, -4), common_term_20);
    result_by_q.set(c_Key1(-4), common_term_20);
    // q = -3
    result_by_lpq.set(c_Key3(7, 6, -3), common_term_19);
    result_by_q.set(c_Key1(-3), common_term_19);
    // q = -2
    result_by_lpq.set(c_Key3(7, 6, -2), common_term_18);
    result_by_q.set(c_Key1(-2), common_term_18);
    // q = -1
    result_by_lpq.set(c_Key3(7, 6, -1), common_term_17);
    result_by_q.set(c_Key1(-1), common_term_17);
    // q = 0
    result_by_lpq.set(c_Key3(7, 6, 0), common_term_16);
    result_by_q.set(c_Key1(0), common_term_16);
    // q = 1
    result_by_lpq.set(c_Key3(7, 6, 1), common_term_15);
    result_by_q.set(c_Key1(1), common_term_15);
    // q = 2
    result_by_lpq.set(c_Key3(7, 6, 2), common_term_14);
    result_by_q.set(c_Key1(2), common_term_14);
    // q = 3
    result_by_lpq.set(c_Key3(7, 6, 3), common_term_13);
    result_by_q.set(c_Key1(3), common_term_13);
    // q = 4
    result_by_lpq.set(c_Key3(7, 6, 4), common_term_12);
    result_by_q.set(c_Key1(4), common_term_12);
    // q = 5
    result_by_lpq.set(c_Key3(7, 6, 5), common_term_11);
    result_by_q.set(c_Key1(5), common_term_11);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 6), result_by_q);
    result_by_q.clear();

    // l , p = (7, 7).
    // q = -5
    result_by_lpq.set(c_Key3(7, 7, -5), common_term_10);
    result_by_q.set(c_Key1(-5), common_term_10);
    // q = -4
    result_by_lpq.set(c_Key3(7, 7, -4), common_term_9);
    result_by_q.set(c_Key1(-4), common_term_9);
    // q = -3
    result_by_lpq.set(c_Key3(7, 7, -3), common_term_8);
    result_by_q.set(c_Key1(-3), common_term_8);
    // q = -2
    result_by_lpq.set(c_Key3(7, 7, -2), common_term_7);
    result_by_q.set(c_Key1(-2), common_term_7);
    // q = -1
    result_by_lpq.set(c_Key3(7, 7, -1), common_term_6);
    result_by_q.set(c_Key1(-1), common_term_6);
    // q = 0
    result_by_lpq.set(c_Key3(7, 7, 0), common_term_5);
    result_by_q.set(c_Key1(0), common_term_5);
    // q = 1
    result_by_lpq.set(c_Key3(7, 7, 1), common_term_4);
    result_by_q.set(c_Key1(1), common_term_4);
    // q = 2
    result_by_lpq.set(c_Key3(7, 7, 2), common_term_3);
    result_by_q.set(c_Key1(2), common_term_3);
    // q = 3
    result_by_lpq.set(c_Key3(7, 7, 3), common_term_2);
    result_by_q.set(c_Key1(3), common_term_2);
    // q = 4
    result_by_lpq.set(c_Key3(7, 7, 4), common_term_1);
    result_by_q.set(c_Key1(4), common_term_1);
    // q = 5
    result_by_lpq.set(c_Key3(7, 7, 5), common_term_0);
    result_by_q.set(c_Key1(5), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 7), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l7_e10(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 7.
    // Functions are taylor expanded and truncated at eccentricity^10.

    // Functions Calculated for l = 7.

    c_IntMap<c_Key3, double> result_by_lpq(166);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(8);
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
    double common_term_0 = 1.5890938895089286e-5*eccentricity_10;
    double common_term_1 = 2.7557319223985891e-6*eccentricity_9;
    double common_term_2 = 2.0452697861552028e-7*eccentricity_10 + 9.6881200396825397e-8*eccentricity_8;
    double common_term_3 = 2.7308388361855159e-5*eccentricity_10 + 2.6351686507936508e-5*eccentricity_8 + 2.1701388888888889e-5*eccentricity_6;
    double common_term_4 = -0.0030505952380952381*eccentricity_9 - 0.0034722222222222222*eccentricity_7 - 0.0083333333333333333*eccentricity_5;
    double common_term_5 = -0.019869559151785714*eccentricity_10 + 0.0085693359375*eccentricity_8 - 0.16875*eccentricity_6 + 0.2109375*eccentricity_4;
    double common_term_6 = 0.64351851851851852*eccentricity_9 - 2.5166666666666667*eccentricity_7 + 3.6666666666666667*eccentricity_5 - 1.3333333333333333*eccentricity_3;
    double common_term_7 = 8.9518299809208623*eccentricity_10 - 25.147840711805556*eccentricity_8 + 35.2783203125*eccentricity_6 - 19.270833333333333*eccentricity_4 + 3.125*eccentricity_2;
    double common_term_8 = -183.90703125*eccentricity_9 + 242.25*eccentricity_7 - 152.8125*eccentricity_5 + 39.75*eccentricity_3 - 3.0*eccentricity;
    double common_term_9 = -1068.9444772677951*eccentricity_10 + 1320.0667656792535*eccentricity_8 - 878.40277777777778*eccentricity_6 + 280.109375*eccentricity_4 - 35.0*eccentricity_2 + 1.0;
    double common_term_10 = 6062.7453125*eccentricity_9 - 4088.2222222222222*eccentricity_7 + 1445.0833333333333*eccentricity_5 - 228.0*eccentricity_3 + 11.0*eccentricity;
    double common_term_11 = 24396.094569396973*eccentricity_10 - 16337.53388671875*eccentricity_8 + 6092.2880859375*eccentricity_6 - 1094.0625*eccentricity_4 + 67.875*eccentricity_2;
    double common_term_12 = -58113.017216435185*eccentricity_9 + 22236.979166666667*eccentricity_7 - 4311.9791666666667*eccentricity_5 + 309.58333333333333*eccentricity_3;
    double common_term_13 = -188463.82743972698*eccentricity_10 + 72772.538658311632*eccentricity_8 - 14777.60625*eccentricity_6 + 1163.0130208333333*eccentricity_4;
    double common_term_14 = 218518.31517857143*eccentricity_9 - 45579.75*eccentricity_7 + 3808.05*eccentricity_5;
    double common_term_15 = 611884.91987441532*eccentricity_10 - 129425.82464812748*eccentricity_8 + 11244.611349826389*eccentricity_6;
    double common_term_16 = -343737.49739583333*eccentricity_9 + 30623.263194444444*eccentricity_7;
    double common_term_17 = -863790.37660871233*eccentricity_10 + 78130.892355782645*eccentricity_8;
    double common_term_18 = 188892.72068176808*eccentricity_9;
    double common_term_19 = 436481.20207806919*eccentricity_10;
    double common_term_20 = 0.8708391879180962*eccentricity_10;
    double common_term_21 = 0.64026262125220459*eccentricity_9;
    double common_term_22 = 2.7070338657924107*eccentricity_10 + 0.47087140764508929*eccentricity_8;
    double common_term_23 = 2.0781994047619048*eccentricity_9 + 0.34637896825396825*eccentricity_7;
    double common_term_24 = 5.7975862654428633*eccentricity_10 + 1.5927439856150794*eccentricity_8 + 0.25483940972222222*eccentricity_6;
    double common_term_25 = 0.1875*eccentricity_5*std::pow(1.0 - eccentricity_2, -6.5);
    double common_term_26 = 10.427123765087632*eccentricity_10 + 3.5954942491319444*eccentricity_8 + 0.93125*eccentricity_6 + 0.13802083333333333*eccentricity_4;
    double common_term_27 = 8.3757523148148148*eccentricity_9 + 2.8197916666666667*eccentricity_7 + 0.70833333333333333*eccentricity_5 + 0.083333333333333333*eccentricity_3;
    double common_term_28 = 16.876936340332031*eccentricity_10 + 6.68759765625*eccentricity_8 + 2.2763671875*eccentricity_6 + 0.1875*eccentricity_4 + 0.375*eccentricity_2;
    double common_term_29 = 13.431510416666667*eccentricity_9 + 6.6944444444444444*eccentricity_7 - 2.1666666666666667*eccentricity_5 + 4.5*eccentricity_3 - eccentricity;
    double common_term_30 = 21.909552680121528*eccentricity_10 + 23.896518283420139*eccentricity_8 - 24.444444444444444*eccentricity_6 + 31.484375*eccentricity_4 - 11.0*eccentricity_2 + 1.0;
    double common_term_31 = 105.17109375*eccentricity_9 - 150.796875*eccentricity_7 + 162.9375*eccentricity_5 - 66.75*eccentricity_3 + 9.0*eccentricity;
    double common_term_32 = 474.60708239520038*eccentricity_10 - 708.84745008680556*eccentricity_8 + 687.3837890625*eccentricity_6 - 297.52083333333333*eccentricity_4 + 46.625*eccentricity_2;
    double common_term_33 = -2801.8372106481481*eccentricity_9 + 2496.9395833333333*eccentricity_7 - 1090.0208333333333*eccentricity_5 + 182.16666666666667*eccentricity_3;
    double common_term_34 = -9767.0560860770089*eccentricity_10 + 8088.3712646484375*eccentricity_8 - 3480.13125*eccentricity_6 + 595.7109375*eccentricity_4;
    double common_term_35 = 23931.665736607143*eccentricity_9 - 10027.289930555556*eccentricity_7 + 1720.3541666666667*eccentricity_5;
    double common_term_36 = 65796.128116704547*eccentricity_10 - 26674.834478856647*eccentricity_8 + 4529.7294487847222*eccentricity_6;
    double common_term_37 = -66559.319866071429*eccentricity_9 + 11102.067857142857*eccentricity_7;
    double common_term_38 = -157569.04723622694*eccentricity_10 + 25694.067438712953*eccentricity_8;
    double common_term_39 = 56735.224209104938*eccentricity_9;
    double common_term_40 = 120454.7790510995*eccentricity_10;
    double common_term_41 = 38.97926045547297*eccentricity_10;
    double common_term_42 = 26.819977678571429*eccentricity_9;
    double common_term_43 = 90.591962165092455*eccentricity_10 + 18.369782462952629*eccentricity_8;
    double common_term_44 = 66.046688988095238*eccentricity_9 + 12.515674603174603*eccentricity_7;
    double common_term_45 = 162.23449042184012*eccentricity_10 + 47.78232421875*eccentricity_8 + 8.4744140625*eccentricity_6;
    double common_term_46 = 121.77518601190476*eccentricity_9 + 34.289930555555556*eccentricity_7 + 5.6958333333333333*eccentricity_5;
    double common_term_47 = 255.31468111875827*eccentricity_10 + 90.715578884548611*eccentricity_8 + 24.39375*eccentricity_6 + 3.7942708333333333*eccentricity_4;
    double common_term_48 = std::pow(1.0 - eccentricity_2, -6.5)*(0.9375*eccentricity_5 + 2.5*eccentricity_3);
    double common_term_49 = 371.54114950674552*eccentricity_10 + 148.61791449652778*eccentricity_8 + 49.0947265625*eccentricity_6 + 11.979166666666667*eccentricity_4 + 1.625*eccentricity_2;
    double common_term_50 = 288.82994791666667*eccentricity_9 + 112.08159722222222*eccentricity_7 + 35.604166666666667*eccentricity_5 + 8.25*eccentricity_3 + eccentricity;
    double common_term_51 = 512.61032958984375*eccentricity_10 + 222.94708251953125*eccentricity_8 + 83.75*eccentricity_6 + 25.734375*eccentricity_4 + 5.0*eccentricity_2 + 1.0;
    double common_term_52 = 403.37942708333333*eccentricity_9 + 170.20138888888889*eccentricity_7 + 64.166666666666667*eccentricity_5 + 13.5*eccentricity_3 + 7.0*eccentricity;
    double common_term_53 = 681.12755528202763*eccentricity_10 + 310.51459418402778*eccentricity_8 + 143.6279296875*eccentricity_6 + 21.354166666666667*eccentricity_4 + 29.375*eccentricity_2;
    double common_term_54 = 513.1875*eccentricity_9 + 312.009375*eccentricity_7 - 0.5*eccentricity_5 + 95.75*eccentricity_3;
    double common_term_55 = 748.01756275318287*eccentricity_10 + 694.08989529079861*eccentricity_8 - 146.54375*eccentricity_6 + 267.13802083333333*eccentricity_4;
    double common_term_56 = 1614.3615327380952*eccentricity_9 - 667.76736111111111*eccentricity_7 + 669.84583333333333*eccentricity_5;
    double common_term_57 = 3882.5545790536063*eccentricity_10 - 2159.8021344866071*eccentricity_8 + 1553.2029296875*eccentricity_6;
    double common_term_58 = -5934.5935639880952*eccentricity_9 + 3391.6381448412698*eccentricity_7;
    double common_term_59 = -14729.937145856154*eccentricity_10 + 7061.7649624294705*eccentricity_8;
    double common_term_60 = 14144.685491071429*eccentricity_9;
    double common_term_61 = 27434.965575509584*eccentricity_10;
    double common_term_62 = 581.64816264561244*eccentricity_10;
    double common_term_63 = 364.20700920414462*eccentricity_9;
    double common_term_64 = 768.54706767517843*eccentricity_10 + 224.84113430447049*eccentricity_8;
    double common_term_65 = 524.66685267857143*eccentricity_9 + 136.44642857142857*eccentricity_7;
    double common_term_66 = 1027.9645262824165*eccentricity_10 + 349.70782955109127*eccentricity_8 + 81.070203993055556*eccentricity_6;
    double common_term_67 = 712.79389880952381*eccentricity_9 + 227.03611111111111*eccentricity_7 + 46.891666666666667*eccentricity_5;
    double common_term_68 = 1274.9632638113839*eccentricity_10 + 483.0539794921875*eccentricity_8 + 142.96875*eccentricity_6 + 26.1796875*eccentricity_4;
    double common_term_69 = 891.58909143518519*eccentricity_9 + 318.51458333333333*eccentricity_7 + 86.729166666666667*eccentricity_5 + 13.916666666666667*eccentricity_3;
    double common_term_70 = 1506.4300543325919*eccentricity_10 + 608.14258897569444*eccentricity_8 + 202.9990234375*eccentricity_6 + 50.104166666666667*eccentricity_4 + 6.875*eccentricity_2;
    double common_term_71 = std::pow(1.0 - eccentricity_2, -6.5)*(1.875*eccentricity_5 + 7.5*eccentricity_3 + 3.0*eccentricity);
    double common_term_72 = 1715.6026451280382*eccentricity_10 + 719.42158338758681*eccentricity_8 + 255.51388888888889*eccentricity_6 + 70.859375*eccentricity_4 + 13.0*eccentricity_2 + 1.0;
    double common_term_73 = 1198.9997395833333*eccentricity_9 + 472.82986111111111*eccentricity_7 + 153.52083333333333*eccentricity_5 + 36.75*eccentricity_3 + 5.0*eccentricity;
    double common_term_74 = 1895.7113754272461*eccentricity_10 + 811.11943359375*eccentricity_8 + 295.7548828125*eccentricity_6 + 84.5625*eccentricity_4 + 16.125*eccentricity_2;
    double common_term_75 = 1315.0491898148148*eccentricity_9 + 525.29166666666667*eccentricity_7 + 170.95833333333333*eccentricity_5 + 42.333333333333333*eccentricity_3;
    double common_term_76 = 2039.7442756128059*eccentricity_10 + 878.53856065538194*eccentricity_8 + 314.40625*eccentricity_6 + 98.294270833333333*eccentricity_4;
    double common_term_77 = 1403.6785714285714*eccentricity_9 + 534.825*eccentricity_7 + 210.15*eccentricity_5;
    double common_term_78 = 2167.3953261481391*eccentricity_10 + 846.33332248263889*eccentricity_8 + 423.28700086805556*eccentricity_6;
    double common_term_79 = 1240.6909598214286*eccentricity_9 + 814.86884920634921*eccentricity_7;
    double common_term_80 = 1653.4173156738281*eccentricity_10 + 1513.8415675571987*eccentricity_8;
    double common_term_81 = 2732.5824997244268*eccentricity_9;
    double common_term_82 = 4816.5759872627595*eccentricity_10;

    c_IntMap<c_Key1, double> result_by_q(21);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (7, 0).
    // q = -10
    result_by_lpq.set(c_Key3(7, 0, -10), common_term_0);
    result_by_q.set(c_Key1(-10), common_term_0);
    // q = -9
    result_by_lpq.set(c_Key3(7, 0, -9), common_term_1);
    result_by_q.set(c_Key1(-9), common_term_1);
    // q = -8
    result_by_lpq.set(c_Key3(7, 0, -8), common_term_2);
    result_by_q.set(c_Key1(-8), common_term_2);
    // q = -6
    result_by_lpq.set(c_Key3(7, 0, -6), common_term_3);
    result_by_q.set(c_Key1(-6), common_term_3);
    // q = -5
    result_by_lpq.set(c_Key3(7, 0, -5), common_term_4);
    result_by_q.set(c_Key1(-5), common_term_4);
    // q = -4
    result_by_lpq.set(c_Key3(7, 0, -4), common_term_5);
    result_by_q.set(c_Key1(-4), common_term_5);
    // q = -3
    result_by_lpq.set(c_Key3(7, 0, -3), common_term_6);
    result_by_q.set(c_Key1(-3), common_term_6);
    // q = -2
    result_by_lpq.set(c_Key3(7, 0, -2), common_term_7);
    result_by_q.set(c_Key1(-2), common_term_7);
    // q = -1
    result_by_lpq.set(c_Key3(7, 0, -1), common_term_8);
    result_by_q.set(c_Key1(-1), common_term_8);
    // q = 0
    result_by_lpq.set(c_Key3(7, 0, 0), common_term_9);
    result_by_q.set(c_Key1(0), common_term_9);
    // q = 1
    result_by_lpq.set(c_Key3(7, 0, 1), common_term_10);
    result_by_q.set(c_Key1(1), common_term_10);
    // q = 2
    result_by_lpq.set(c_Key3(7, 0, 2), common_term_11);
    result_by_q.set(c_Key1(2), common_term_11);
    // q = 3
    result_by_lpq.set(c_Key3(7, 0, 3), common_term_12);
    result_by_q.set(c_Key1(3), common_term_12);
    // q = 4
    result_by_lpq.set(c_Key3(7, 0, 4), common_term_13);
    result_by_q.set(c_Key1(4), common_term_13);
    // q = 5
    result_by_lpq.set(c_Key3(7, 0, 5), common_term_14);
    result_by_q.set(c_Key1(5), common_term_14);
    // q = 6
    result_by_lpq.set(c_Key3(7, 0, 6), common_term_15);
    result_by_q.set(c_Key1(6), common_term_15);
    // q = 7
    result_by_lpq.set(c_Key3(7, 0, 7), common_term_16);
    result_by_q.set(c_Key1(7), common_term_16);
    // q = 8
    result_by_lpq.set(c_Key3(7, 0, 8), common_term_17);
    result_by_q.set(c_Key1(8), common_term_17);
    // q = 9
    result_by_lpq.set(c_Key3(7, 0, 9), common_term_18);
    result_by_q.set(c_Key1(9), common_term_18);
    // q = 10
    result_by_lpq.set(c_Key3(7, 0, 10), common_term_19);
    result_by_q.set(c_Key1(10), common_term_19);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 0), result_by_q);
    result_by_q.clear();

    // l , p = (7, 1).
    // q = -10
    result_by_lpq.set(c_Key3(7, 1, -10), common_term_20);
    result_by_q.set(c_Key1(-10), common_term_20);
    // q = -9
    result_by_lpq.set(c_Key3(7, 1, -9), common_term_21);
    result_by_q.set(c_Key1(-9), common_term_21);
    // q = -8
    result_by_lpq.set(c_Key3(7, 1, -8), common_term_22);
    result_by_q.set(c_Key1(-8), common_term_22);
    // q = -7
    result_by_lpq.set(c_Key3(7, 1, -7), common_term_23);
    result_by_q.set(c_Key1(-7), common_term_23);
    // q = -6
    result_by_lpq.set(c_Key3(7, 1, -6), common_term_24);
    result_by_q.set(c_Key1(-6), common_term_24);
    // q = -5
    result_by_lpq.set(c_Key3(7, 1, -5), common_term_25);
    result_by_q.set(c_Key1(-5), common_term_25);
    // q = -4
    result_by_lpq.set(c_Key3(7, 1, -4), common_term_26);
    result_by_q.set(c_Key1(-4), common_term_26);
    // q = -3
    result_by_lpq.set(c_Key3(7, 1, -3), common_term_27);
    result_by_q.set(c_Key1(-3), common_term_27);
    // q = -2
    result_by_lpq.set(c_Key3(7, 1, -2), common_term_28);
    result_by_q.set(c_Key1(-2), common_term_28);
    // q = -1
    result_by_lpq.set(c_Key3(7, 1, -1), common_term_29);
    result_by_q.set(c_Key1(-1), common_term_29);
    // q = 0
    result_by_lpq.set(c_Key3(7, 1, 0), common_term_30);
    result_by_q.set(c_Key1(0), common_term_30);
    // q = 1
    result_by_lpq.set(c_Key3(7, 1, 1), common_term_31);
    result_by_q.set(c_Key1(1), common_term_31);
    // q = 2
    result_by_lpq.set(c_Key3(7, 1, 2), common_term_32);
    result_by_q.set(c_Key1(2), common_term_32);
    // q = 3
    result_by_lpq.set(c_Key3(7, 1, 3), common_term_33);
    result_by_q.set(c_Key1(3), common_term_33);
    // q = 4
    result_by_lpq.set(c_Key3(7, 1, 4), common_term_34);
    result_by_q.set(c_Key1(4), common_term_34);
    // q = 5
    result_by_lpq.set(c_Key3(7, 1, 5), common_term_35);
    result_by_q.set(c_Key1(5), common_term_35);
    // q = 6
    result_by_lpq.set(c_Key3(7, 1, 6), common_term_36);
    result_by_q.set(c_Key1(6), common_term_36);
    // q = 7
    result_by_lpq.set(c_Key3(7, 1, 7), common_term_37);
    result_by_q.set(c_Key1(7), common_term_37);
    // q = 8
    result_by_lpq.set(c_Key3(7, 1, 8), common_term_38);
    result_by_q.set(c_Key1(8), common_term_38);
    // q = 9
    result_by_lpq.set(c_Key3(7, 1, 9), common_term_39);
    result_by_q.set(c_Key1(9), common_term_39);
    // q = 10
    result_by_lpq.set(c_Key3(7, 1, 10), common_term_40);
    result_by_q.set(c_Key1(10), common_term_40);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 1), result_by_q);
    result_by_q.clear();

    // l , p = (7, 2).
    // q = -10
    result_by_lpq.set(c_Key3(7, 2, -10), common_term_41);
    result_by_q.set(c_Key1(-10), common_term_41);
    // q = -9
    result_by_lpq.set(c_Key3(7, 2, -9), common_term_42);
    result_by_q.set(c_Key1(-9), common_term_42);
    // q = -8
    result_by_lpq.set(c_Key3(7, 2, -8), common_term_43);
    result_by_q.set(c_Key1(-8), common_term_43);
    // q = -7
    result_by_lpq.set(c_Key3(7, 2, -7), common_term_44);
    result_by_q.set(c_Key1(-7), common_term_44);
    // q = -6
    result_by_lpq.set(c_Key3(7, 2, -6), common_term_45);
    result_by_q.set(c_Key1(-6), common_term_45);
    // q = -5
    result_by_lpq.set(c_Key3(7, 2, -5), common_term_46);
    result_by_q.set(c_Key1(-5), common_term_46);
    // q = -4
    result_by_lpq.set(c_Key3(7, 2, -4), common_term_47);
    result_by_q.set(c_Key1(-4), common_term_47);
    // q = -3
    result_by_lpq.set(c_Key3(7, 2, -3), common_term_48);
    result_by_q.set(c_Key1(-3), common_term_48);
    // q = -2
    result_by_lpq.set(c_Key3(7, 2, -2), common_term_49);
    result_by_q.set(c_Key1(-2), common_term_49);
    // q = -1
    result_by_lpq.set(c_Key3(7, 2, -1), common_term_50);
    result_by_q.set(c_Key1(-1), common_term_50);
    // q = 0
    result_by_lpq.set(c_Key3(7, 2, 0), common_term_51);
    result_by_q.set(c_Key1(0), common_term_51);
    // q = 1
    result_by_lpq.set(c_Key3(7, 2, 1), common_term_52);
    result_by_q.set(c_Key1(1), common_term_52);
    // q = 2
    result_by_lpq.set(c_Key3(7, 2, 2), common_term_53);
    result_by_q.set(c_Key1(2), common_term_53);
    // q = 3
    result_by_lpq.set(c_Key3(7, 2, 3), common_term_54);
    result_by_q.set(c_Key1(3), common_term_54);
    // q = 4
    result_by_lpq.set(c_Key3(7, 2, 4), common_term_55);
    result_by_q.set(c_Key1(4), common_term_55);
    // q = 5
    result_by_lpq.set(c_Key3(7, 2, 5), common_term_56);
    result_by_q.set(c_Key1(5), common_term_56);
    // q = 6
    result_by_lpq.set(c_Key3(7, 2, 6), common_term_57);
    result_by_q.set(c_Key1(6), common_term_57);
    // q = 7
    result_by_lpq.set(c_Key3(7, 2, 7), common_term_58);
    result_by_q.set(c_Key1(7), common_term_58);
    // q = 8
    result_by_lpq.set(c_Key3(7, 2, 8), common_term_59);
    result_by_q.set(c_Key1(8), common_term_59);
    // q = 9
    result_by_lpq.set(c_Key3(7, 2, 9), common_term_60);
    result_by_q.set(c_Key1(9), common_term_60);
    // q = 10
    result_by_lpq.set(c_Key3(7, 2, 10), common_term_61);
    result_by_q.set(c_Key1(10), common_term_61);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 2), result_by_q);
    result_by_q.clear();

    // l , p = (7, 3).
    // q = -10
    result_by_lpq.set(c_Key3(7, 3, -10), common_term_62);
    result_by_q.set(c_Key1(-10), common_term_62);
    // q = -9
    result_by_lpq.set(c_Key3(7, 3, -9), common_term_63);
    result_by_q.set(c_Key1(-9), common_term_63);
    // q = -8
    result_by_lpq.set(c_Key3(7, 3, -8), common_term_64);
    result_by_q.set(c_Key1(-8), common_term_64);
    // q = -7
    result_by_lpq.set(c_Key3(7, 3, -7), common_term_65);
    result_by_q.set(c_Key1(-7), common_term_65);
    // q = -6
    result_by_lpq.set(c_Key3(7, 3, -6), common_term_66);
    result_by_q.set(c_Key1(-6), common_term_66);
    // q = -5
    result_by_lpq.set(c_Key3(7, 3, -5), common_term_67);
    result_by_q.set(c_Key1(-5), common_term_67);
    // q = -4
    result_by_lpq.set(c_Key3(7, 3, -4), common_term_68);
    result_by_q.set(c_Key1(-4), common_term_68);
    // q = -3
    result_by_lpq.set(c_Key3(7, 3, -3), common_term_69);
    result_by_q.set(c_Key1(-3), common_term_69);
    // q = -2
    result_by_lpq.set(c_Key3(7, 3, -2), common_term_70);
    result_by_q.set(c_Key1(-2), common_term_70);
    // q = -1
    result_by_lpq.set(c_Key3(7, 3, -1), common_term_71);
    result_by_q.set(c_Key1(-1), common_term_71);
    // q = 0
    result_by_lpq.set(c_Key3(7, 3, 0), common_term_72);
    result_by_q.set(c_Key1(0), common_term_72);
    // q = 1
    result_by_lpq.set(c_Key3(7, 3, 1), common_term_73);
    result_by_q.set(c_Key1(1), common_term_73);
    // q = 2
    result_by_lpq.set(c_Key3(7, 3, 2), common_term_74);
    result_by_q.set(c_Key1(2), common_term_74);
    // q = 3
    result_by_lpq.set(c_Key3(7, 3, 3), common_term_75);
    result_by_q.set(c_Key1(3), common_term_75);
    // q = 4
    result_by_lpq.set(c_Key3(7, 3, 4), common_term_76);
    result_by_q.set(c_Key1(4), common_term_76);
    // q = 5
    result_by_lpq.set(c_Key3(7, 3, 5), common_term_77);
    result_by_q.set(c_Key1(5), common_term_77);
    // q = 6
    result_by_lpq.set(c_Key3(7, 3, 6), common_term_78);
    result_by_q.set(c_Key1(6), common_term_78);
    // q = 7
    result_by_lpq.set(c_Key3(7, 3, 7), common_term_79);
    result_by_q.set(c_Key1(7), common_term_79);
    // q = 8
    result_by_lpq.set(c_Key3(7, 3, 8), common_term_80);
    result_by_q.set(c_Key1(8), common_term_80);
    // q = 9
    result_by_lpq.set(c_Key3(7, 3, 9), common_term_81);
    result_by_q.set(c_Key1(9), common_term_81);
    // q = 10
    result_by_lpq.set(c_Key3(7, 3, 10), common_term_82);
    result_by_q.set(c_Key1(10), common_term_82);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 3), result_by_q);
    result_by_q.clear();

    // l , p = (7, 4).
    // q = -10
    result_by_lpq.set(c_Key3(7, 4, -10), common_term_82);
    result_by_q.set(c_Key1(-10), common_term_82);
    // q = -9
    result_by_lpq.set(c_Key3(7, 4, -9), common_term_81);
    result_by_q.set(c_Key1(-9), common_term_81);
    // q = -8
    result_by_lpq.set(c_Key3(7, 4, -8), common_term_80);
    result_by_q.set(c_Key1(-8), common_term_80);
    // q = -7
    result_by_lpq.set(c_Key3(7, 4, -7), common_term_79);
    result_by_q.set(c_Key1(-7), common_term_79);
    // q = -6
    result_by_lpq.set(c_Key3(7, 4, -6), common_term_78);
    result_by_q.set(c_Key1(-6), common_term_78);
    // q = -5
    result_by_lpq.set(c_Key3(7, 4, -5), common_term_77);
    result_by_q.set(c_Key1(-5), common_term_77);
    // q = -4
    result_by_lpq.set(c_Key3(7, 4, -4), common_term_76);
    result_by_q.set(c_Key1(-4), common_term_76);
    // q = -3
    result_by_lpq.set(c_Key3(7, 4, -3), common_term_75);
    result_by_q.set(c_Key1(-3), common_term_75);
    // q = -2
    result_by_lpq.set(c_Key3(7, 4, -2), common_term_74);
    result_by_q.set(c_Key1(-2), common_term_74);
    // q = -1
    result_by_lpq.set(c_Key3(7, 4, -1), common_term_73);
    result_by_q.set(c_Key1(-1), common_term_73);
    // q = 0
    result_by_lpq.set(c_Key3(7, 4, 0), common_term_72);
    result_by_q.set(c_Key1(0), common_term_72);
    // q = 1
    result_by_lpq.set(c_Key3(7, 4, 1), common_term_71);
    result_by_q.set(c_Key1(1), common_term_71);
    // q = 2
    result_by_lpq.set(c_Key3(7, 4, 2), common_term_70);
    result_by_q.set(c_Key1(2), common_term_70);
    // q = 3
    result_by_lpq.set(c_Key3(7, 4, 3), common_term_69);
    result_by_q.set(c_Key1(3), common_term_69);
    // q = 4
    result_by_lpq.set(c_Key3(7, 4, 4), common_term_68);
    result_by_q.set(c_Key1(4), common_term_68);
    // q = 5
    result_by_lpq.set(c_Key3(7, 4, 5), common_term_67);
    result_by_q.set(c_Key1(5), common_term_67);
    // q = 6
    result_by_lpq.set(c_Key3(7, 4, 6), common_term_66);
    result_by_q.set(c_Key1(6), common_term_66);
    // q = 7
    result_by_lpq.set(c_Key3(7, 4, 7), common_term_65);
    result_by_q.set(c_Key1(7), common_term_65);
    // q = 8
    result_by_lpq.set(c_Key3(7, 4, 8), common_term_64);
    result_by_q.set(c_Key1(8), common_term_64);
    // q = 9
    result_by_lpq.set(c_Key3(7, 4, 9), common_term_63);
    result_by_q.set(c_Key1(9), common_term_63);
    // q = 10
    result_by_lpq.set(c_Key3(7, 4, 10), common_term_62);
    result_by_q.set(c_Key1(10), common_term_62);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 4), result_by_q);
    result_by_q.clear();

    // l , p = (7, 5).
    // q = -10
    result_by_lpq.set(c_Key3(7, 5, -10), common_term_61);
    result_by_q.set(c_Key1(-10), common_term_61);
    // q = -9
    result_by_lpq.set(c_Key3(7, 5, -9), common_term_60);
    result_by_q.set(c_Key1(-9), common_term_60);
    // q = -8
    result_by_lpq.set(c_Key3(7, 5, -8), common_term_59);
    result_by_q.set(c_Key1(-8), common_term_59);
    // q = -7
    result_by_lpq.set(c_Key3(7, 5, -7), common_term_58);
    result_by_q.set(c_Key1(-7), common_term_58);
    // q = -6
    result_by_lpq.set(c_Key3(7, 5, -6), common_term_57);
    result_by_q.set(c_Key1(-6), common_term_57);
    // q = -5
    result_by_lpq.set(c_Key3(7, 5, -5), common_term_56);
    result_by_q.set(c_Key1(-5), common_term_56);
    // q = -4
    result_by_lpq.set(c_Key3(7, 5, -4), common_term_55);
    result_by_q.set(c_Key1(-4), common_term_55);
    // q = -3
    result_by_lpq.set(c_Key3(7, 5, -3), common_term_54);
    result_by_q.set(c_Key1(-3), common_term_54);
    // q = -2
    result_by_lpq.set(c_Key3(7, 5, -2), common_term_53);
    result_by_q.set(c_Key1(-2), common_term_53);
    // q = -1
    result_by_lpq.set(c_Key3(7, 5, -1), common_term_52);
    result_by_q.set(c_Key1(-1), common_term_52);
    // q = 0
    result_by_lpq.set(c_Key3(7, 5, 0), common_term_51);
    result_by_q.set(c_Key1(0), common_term_51);
    // q = 1
    result_by_lpq.set(c_Key3(7, 5, 1), common_term_50);
    result_by_q.set(c_Key1(1), common_term_50);
    // q = 2
    result_by_lpq.set(c_Key3(7, 5, 2), common_term_49);
    result_by_q.set(c_Key1(2), common_term_49);
    // q = 3
    result_by_lpq.set(c_Key3(7, 5, 3), common_term_48);
    result_by_q.set(c_Key1(3), common_term_48);
    // q = 4
    result_by_lpq.set(c_Key3(7, 5, 4), common_term_47);
    result_by_q.set(c_Key1(4), common_term_47);
    // q = 5
    result_by_lpq.set(c_Key3(7, 5, 5), common_term_46);
    result_by_q.set(c_Key1(5), common_term_46);
    // q = 6
    result_by_lpq.set(c_Key3(7, 5, 6), common_term_45);
    result_by_q.set(c_Key1(6), common_term_45);
    // q = 7
    result_by_lpq.set(c_Key3(7, 5, 7), common_term_44);
    result_by_q.set(c_Key1(7), common_term_44);
    // q = 8
    result_by_lpq.set(c_Key3(7, 5, 8), common_term_43);
    result_by_q.set(c_Key1(8), common_term_43);
    // q = 9
    result_by_lpq.set(c_Key3(7, 5, 9), common_term_42);
    result_by_q.set(c_Key1(9), common_term_42);
    // q = 10
    result_by_lpq.set(c_Key3(7, 5, 10), common_term_41);
    result_by_q.set(c_Key1(10), common_term_41);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 5), result_by_q);
    result_by_q.clear();

    // l , p = (7, 6).
    // q = -10
    result_by_lpq.set(c_Key3(7, 6, -10), common_term_40);
    result_by_q.set(c_Key1(-10), common_term_40);
    // q = -9
    result_by_lpq.set(c_Key3(7, 6, -9), common_term_39);
    result_by_q.set(c_Key1(-9), common_term_39);
    // q = -8
    result_by_lpq.set(c_Key3(7, 6, -8), common_term_38);
    result_by_q.set(c_Key1(-8), common_term_38);
    // q = -7
    result_by_lpq.set(c_Key3(7, 6, -7), common_term_37);
    result_by_q.set(c_Key1(-7), common_term_37);
    // q = -6
    result_by_lpq.set(c_Key3(7, 6, -6), common_term_36);
    result_by_q.set(c_Key1(-6), common_term_36);
    // q = -5
    result_by_lpq.set(c_Key3(7, 6, -5), common_term_35);
    result_by_q.set(c_Key1(-5), common_term_35);
    // q = -4
    result_by_lpq.set(c_Key3(7, 6, -4), common_term_34);
    result_by_q.set(c_Key1(-4), common_term_34);
    // q = -3
    result_by_lpq.set(c_Key3(7, 6, -3), common_term_33);
    result_by_q.set(c_Key1(-3), common_term_33);
    // q = -2
    result_by_lpq.set(c_Key3(7, 6, -2), common_term_32);
    result_by_q.set(c_Key1(-2), common_term_32);
    // q = -1
    result_by_lpq.set(c_Key3(7, 6, -1), common_term_31);
    result_by_q.set(c_Key1(-1), common_term_31);
    // q = 0
    result_by_lpq.set(c_Key3(7, 6, 0), common_term_30);
    result_by_q.set(c_Key1(0), common_term_30);
    // q = 1
    result_by_lpq.set(c_Key3(7, 6, 1), common_term_29);
    result_by_q.set(c_Key1(1), common_term_29);
    // q = 2
    result_by_lpq.set(c_Key3(7, 6, 2), common_term_28);
    result_by_q.set(c_Key1(2), common_term_28);
    // q = 3
    result_by_lpq.set(c_Key3(7, 6, 3), common_term_27);
    result_by_q.set(c_Key1(3), common_term_27);
    // q = 4
    result_by_lpq.set(c_Key3(7, 6, 4), common_term_26);
    result_by_q.set(c_Key1(4), common_term_26);
    // q = 5
    result_by_lpq.set(c_Key3(7, 6, 5), common_term_25);
    result_by_q.set(c_Key1(5), common_term_25);
    // q = 6
    result_by_lpq.set(c_Key3(7, 6, 6), common_term_24);
    result_by_q.set(c_Key1(6), common_term_24);
    // q = 7
    result_by_lpq.set(c_Key3(7, 6, 7), common_term_23);
    result_by_q.set(c_Key1(7), common_term_23);
    // q = 8
    result_by_lpq.set(c_Key3(7, 6, 8), common_term_22);
    result_by_q.set(c_Key1(8), common_term_22);
    // q = 9
    result_by_lpq.set(c_Key3(7, 6, 9), common_term_21);
    result_by_q.set(c_Key1(9), common_term_21);
    // q = 10
    result_by_lpq.set(c_Key3(7, 6, 10), common_term_20);
    result_by_q.set(c_Key1(10), common_term_20);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 6), result_by_q);
    result_by_q.clear();

    // l , p = (7, 7).
    // q = -10
    result_by_lpq.set(c_Key3(7, 7, -10), common_term_19);
    result_by_q.set(c_Key1(-10), common_term_19);
    // q = -9
    result_by_lpq.set(c_Key3(7, 7, -9), common_term_18);
    result_by_q.set(c_Key1(-9), common_term_18);
    // q = -8
    result_by_lpq.set(c_Key3(7, 7, -8), common_term_17);
    result_by_q.set(c_Key1(-8), common_term_17);
    // q = -7
    result_by_lpq.set(c_Key3(7, 7, -7), common_term_16);
    result_by_q.set(c_Key1(-7), common_term_16);
    // q = -6
    result_by_lpq.set(c_Key3(7, 7, -6), common_term_15);
    result_by_q.set(c_Key1(-6), common_term_15);
    // q = -5
    result_by_lpq.set(c_Key3(7, 7, -5), common_term_14);
    result_by_q.set(c_Key1(-5), common_term_14);
    // q = -4
    result_by_lpq.set(c_Key3(7, 7, -4), common_term_13);
    result_by_q.set(c_Key1(-4), common_term_13);
    // q = -3
    result_by_lpq.set(c_Key3(7, 7, -3), common_term_12);
    result_by_q.set(c_Key1(-3), common_term_12);
    // q = -2
    result_by_lpq.set(c_Key3(7, 7, -2), common_term_11);
    result_by_q.set(c_Key1(-2), common_term_11);
    // q = -1
    result_by_lpq.set(c_Key3(7, 7, -1), common_term_10);
    result_by_q.set(c_Key1(-1), common_term_10);
    // q = 0
    result_by_lpq.set(c_Key3(7, 7, 0), common_term_9);
    result_by_q.set(c_Key1(0), common_term_9);
    // q = 1
    result_by_lpq.set(c_Key3(7, 7, 1), common_term_8);
    result_by_q.set(c_Key1(1), common_term_8);
    // q = 2
    result_by_lpq.set(c_Key3(7, 7, 2), common_term_7);
    result_by_q.set(c_Key1(2), common_term_7);
    // q = 3
    result_by_lpq.set(c_Key3(7, 7, 3), common_term_6);
    result_by_q.set(c_Key1(3), common_term_6);
    // q = 4
    result_by_lpq.set(c_Key3(7, 7, 4), common_term_5);
    result_by_q.set(c_Key1(4), common_term_5);
    // q = 5
    result_by_lpq.set(c_Key3(7, 7, 5), common_term_4);
    result_by_q.set(c_Key1(5), common_term_4);
    // q = 6
    result_by_lpq.set(c_Key3(7, 7, 6), common_term_3);
    result_by_q.set(c_Key1(6), common_term_3);
    // q = 8
    result_by_lpq.set(c_Key3(7, 7, 8), common_term_2);
    result_by_q.set(c_Key1(8), common_term_2);
    // q = 9
    result_by_lpq.set(c_Key3(7, 7, 9), common_term_1);
    result_by_q.set(c_Key1(9), common_term_1);
    // q = 10
    result_by_lpq.set(c_Key3(7, 7, 10), common_term_0);
    result_by_q.set(c_Key1(10), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 7), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l7_e15(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 7.
    // Functions are taylor expanded and truncated at eccentricity^15.

    // Functions Calculated for l = 7.

    c_IntMap<c_Key3, double> result_by_lpq(246);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(8);
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
    double common_term_0 = 0.00082110795338308566*eccentricity_15;
    double common_term_1 = 0.00047483668989022175*eccentricity_14;
    double common_term_2 = 0.00066751468397673755*eccentricity_15 + 0.00025603302947052947*eccentricity_13;
    double common_term_3 = 0.00032544582513090068*eccentricity_14 + 0.0001244351684324032*eccentricity_12;
    double common_term_4 = 0.0002207668874335541*eccentricity_15 + 0.00013254235476457699*eccentricity_13 + 5.130671797338464e-5*eccentricity_11;
    double common_term_5 = 6.3981344173480938e-5*eccentricity_14 + 3.9727347237723214e-5*eccentricity_12 + 1.5890938895089286e-5*eccentricity_10;
    double common_term_6 = 1.2582073490638305e-5*eccentricity_15 + 9.8987143257976591e-6*eccentricity_13 + 6.4759700176366843e-6*eccentricity_11 + 2.7557319223985891e-6*eccentricity_9;
    double common_term_7 = 3.4561838167659565e-7*eccentricity_14 + 2.8872616076595569e-7*eccentricity_12 + 2.0452697861552028e-7*eccentricity_10 + 9.6881200396825397e-8*eccentricity_8;
    double common_term_8 = 2.5351391617933629e-5*eccentricity_14 + 2.6653318948986809e-5*eccentricity_12 + 2.7308388361855159e-5*eccentricity_10 + 2.6351686507936508e-5*eccentricity_8 + 2.1701388888888889e-5*eccentricity_6;
    double common_term_9 = -0.0017635168650793651*eccentricity_15 - 0.0020673730526161082*eccentricity_13 - 0.0024698247354497354*eccentricity_11 - 0.0030505952380952381*eccentricity_9 - 0.0034722222222222222*eccentricity_7 - 0.0083333333333333333*eccentricity_5;
    double common_term_10 = -0.011672703879220145*eccentricity_14 - 0.013487250464303153*eccentricity_12 - 0.019869559151785714*eccentricity_10 + 0.0085693359375*eccentricity_8 - 0.16875*eccentricity_6 + 0.2109375*eccentricity_4;
    double common_term_11 = -0.037010061483441113*eccentricity_15 - 0.028399470899470899*eccentricity_13 - 0.16982473544973545*eccentricity_11 + 0.64351851851851852*eccentricity_9 - 2.5166666666666667*eccentricity_7 + 3.6666666666666667*eccentricity_5 - 1.3333333333333333*eccentricity_3;
    double common_term_12 = 0.22088288603846568*eccentricity_14 - 2.1656081790015811*eccentricity_12 + 8.9518299809208623*eccentricity_10 - 25.147840711805556*eccentricity_8 + 35.2783203125*eccentricity_6 - 19.270833333333333*eccentricity_4 + 3.125*eccentricity_2;
    double common_term_13 = 4.1779526068239796*eccentricity_15 - 22.384762834821429*eccentricity_13 + 78.94640625*eccentricity_11 - 183.90703125*eccentricity_9 + 242.25*eccentricity_7 - 152.8125*eccentricity_5 + 39.75*eccentricity_3 - 3.0*eccentricity;
    double common_term_14 = -175.43335821740421*eccentricity_14 + 526.83805451334259*eccentricity_12 - 1068.9444772677951*eccentricity_10 + 1320.0667656792535*eccentricity_8 - 878.40277777777778*eccentricity_6 + 280.109375*eccentricity_4 - 35.0*eccentricity_2 + 1.0;
    double common_term_15 = -1101.9652262849584*eccentricity_15 + 2880.6147021053792*eccentricity_13 - 5223.6626851851852*eccentricity_11 + 6062.7453125*eccentricity_9 - 4088.2222222222222*eccentricity_7 + 1445.0833333333333*eccentricity_5 - 228.0*eccentricity_3 + 11.0*eccentricity;
    double common_term_16 = 13530.386257162264*eccentricity_14 - 22294.568359745571*eccentricity_12 + 24396.094569396973*eccentricity_10 - 16337.53388671875*eccentricity_8 + 6092.2880859375*eccentricity_6 - 1094.0625*eccentricity_4 + 67.875*eccentricity_2;
    double common_term_17 = 56343.268666686952*eccentricity_15 - 85355.961694465112*eccentricity_13 + 88336.957981977513*eccentricity_11 - 58113.017216435185*eccentricity_9 + 22236.979166666667*eccentricity_7 - 4311.9791666666667*eccentricity_5 + 309.58333333333333*eccentricity_3;
    double common_term_18 = -298897.11036031678*eccentricity_14 + 293421.13897960829*eccentricity_12 - 188463.82743972698*eccentricity_10 + 72772.538658311632*eccentricity_8 - 14777.60625*eccentricity_6 + 1163.0130208333333*eccentricity_4;
    double common_term_19 = -971396.39525223214*eccentricity_15 + 906983.55929129464*eccentricity_13 - 566805.24776785714*eccentricity_11 + 218518.31517857143*eccentricity_9 - 45579.75*eccentricity_7 + 3808.05*eccentricity_5;
    double common_term_20 = 2637812.1504471274*eccentricity_14 - 1600816.3808199213*eccentricity_12 + 611884.91987441532*eccentricity_10 - 129425.82464812748*eccentricity_8 + 11244.611349826389*eccentricity_6;
    double common_term_21 = 7280654.3524639735*eccentricity_15 - 4286530.2131702514*eccentricity_13 + 1616827.3261091821*eccentricity_11 - 343737.49739583333*eccentricity_9 + 30623.263194444444*eccentricity_7;
    double common_term_22 = -10964172.505874448*eccentricity_14 + 4067895.538355623*eccentricity_12 - 863790.37660871233*eccentricity_10 + 78130.892355782645*eccentricity_8;
    double common_term_23 = -26949151.751089544*eccentricity_15 + 9813445.3178285584*eccentricity_13 - 2071825.9252954145*eccentricity_11 + 188892.72068176808*eccentricity_9;
    double common_term_24 = 22825990.404013062*eccentricity_14 - 4775273.3536105321*eccentricity_12 + 436481.20207806919*eccentricity_10;
    double common_term_25 = 51421560.566231415*eccentricity_15 - 10633351.945539265*eccentricity_13 + 970435.77748782468*eccentricity_11;
    double common_term_26 = -22974539.901608997*eccentricity_14 + 2086930.9663047792*eccentricity_12;
    double common_term_27 = -48335585.125646251*eccentricity_15 + 4359471.0434132282*eccentricity_13;
    double common_term_28 = 8876731.8003844841*eccentricity_14;
    double common_term_29 = 17669216.908951237*eccentricity_15;
    double common_term_30 = 4.0699611066729458*eccentricity_15;
    double common_term_31 = 2.9885374554252938*eccentricity_14;
    double common_term_32 = 9.8401352562385394*eccentricity_15 + 2.1949078031423952*eccentricity_13;
    double common_term_33 = 7.6389896802115383*eccentricity_14 + 1.6124069097235637*eccentricity_12;
    double common_term_34 = 18.399222727955638*eccentricity_15 + 5.9140178571428571*eccentricity_13 + 1.1848006290584416*eccentricity_11;
    double common_term_35 = 14.677134880754364*eccentricity_14 + 4.5674260705064385*eccentricity_12 + 0.8708391879180962*eccentricity_10;
    double common_term_36 = 30.099295869088259*eccentricity_15 + 11.677661204054433*eccentricity_13 + 3.5197512951940035*eccentricity_11 + 0.64026262125220459*eccentricity_9;
    double common_term_37 = 24.464791772520387*eccentricity_14 + 9.2682209750584194*eccentricity_12 + 2.7070338657924107*eccentricity_10 + 0.47087140764508929*eccentricity_8;
    double common_term_38 = 45.354514222395082*eccentricity_15 + 19.836474259167402*eccentricity_13 + 7.3385771467151675*eccentricity_11 + 2.0781994047619048*eccentricity_9 + 0.34637896825396825*eccentricity_7;
    double common_term_39 = 37.390312413537874*eccentricity_14 + 16.045561849559356*eccentricity_12 + 5.7975862654428633*eccentricity_10 + 1.5927439856150794*eccentricity_8 + 0.25483940972222222*eccentricity_6;
    double common_term_40 = 0.1875*eccentricity_5*std::pow(1.0 - eccentricity_2, -6.5);
    double common_term_41 = 53.849767283159869*eccentricity_14 + 25.239834294243464*eccentricity_12 + 10.427123765087632*eccentricity_10 + 3.5954942491319444*eccentricity_8 + 0.93125*eccentricity_6 + 0.13802083333333333*eccentricity_4;
    double common_term_42 = 88.215787991974696*eccentricity_15 + 44.802796792328042*eccentricity_13 + 20.667006655092593*eccentricity_11 + 8.3757523148148148*eccentricity_9 + 2.8197916666666667*eccentricity_7 + 0.70833333333333333*eccentricity_5 + 0.083333333333333333*eccentricity_3;
    double common_term_43 = 74.238951664822442*eccentricity_14 + 37.189341234479632*eccentricity_12 + 16.876936340332031*eccentricity_10 + 6.68759765625*eccentricity_8 + 2.2763671875*eccentricity_6 + 0.1875*eccentricity_4 + 0.375*eccentricity_2;
    double common_term_44 = 116.65986192306784*eccentricity_15 + 62.328202229387125*eccentricity_13 + 30.817219328703704*eccentricity_11 + 13.431510416666667*eccentricity_9 + 6.6944444444444444*eccentricity_7 - 2.1666666666666667*eccentricity_5 + 4.5*eccentricity_3 - eccentricity;
    double common_term_45 = 98.832554347450831*eccentricity_14 + 52.795619199305405*eccentricity_12 + 21.909552680121528*eccentricity_10 + 23.896518283420139*eccentricity_8 - 24.444444444444444*eccentricity_6 + 31.484375*eccentricity_4 - 11.0*eccentricity_2 + 1.0;
    double common_term_46 = 149.22814378388074*eccentricity_15 + 89.920555245535714*eccentricity_13 + 15.6230859375*eccentricity_11 + 105.17109375*eccentricity_9 - 150.796875*eccentricity_7 + 162.9375*eccentricity_5 - 66.75*eccentricity_3 + 9.0*eccentricity;
    double common_term_47 = 174.9804609570827*eccentricity_14 - 100.78621590508355*eccentricity_12 + 474.60708239520038*eccentricity_10 - 708.84745008680556*eccentricity_8 + 687.3837890625*eccentricity_6 - 297.52083333333333*eccentricity_4 + 46.625*eccentricity_2;
    double common_term_48 = 465.78605284556113*eccentricity_15 - 752.55108129753638*eccentricity_13 + 1975.8828786375661*eccentricity_11 - 2801.8372106481481*eccentricity_9 + 2496.9395833333333*eccentricity_7 - 1090.0208333333333*eccentricity_5 + 182.16666666666667*eccentricity_3;
    double common_term_49 = -3561.5581530353001*eccentricity_14 + 7431.9097426210131*eccentricity_12 - 9767.0560860770089*eccentricity_10 + 8088.3712646484375*eccentricity_8 - 3480.13125*eccentricity_6 + 595.7109375*eccentricity_4;
    double common_term_50 = -14047.947527096595*eccentricity_15 + 25445.885708027906*eccentricity_13 - 30894.854342344577*eccentricity_11 + 23931.665736607143*eccentricity_9 - 10027.289930555556*eccentricity_7 + 1720.3541666666667*eccentricity_5;
    double common_term_51 = 80336.397336890789*eccentricity_14 - 90383.69587827253*eccentricity_12 + 65796.128116704547*eccentricity_10 - 26674.834478856647*eccentricity_8 + 4529.7294487847222*eccentricity_6;
    double common_term_52 = 236722.20799994927*eccentricity_15 - 247917.21170758929*eccentricity_13 + 170241.80323660714*eccentricity_11 - 66559.319866071429*eccentricity_9 + 11102.067857142857*eccentricity_7;
    double common_term_53 = -644137.50143032605*eccentricity_14 + 418575.6100948585*eccentricity_12 - 157569.04723622694*eccentricity_10 + 25694.067438712953*eccentricity_8;
    double common_term_54 = -1597867.3251102959*eccentricity_15 + 985361.68525131962*eccentricity_13 - 356951.15354841821*eccentricity_11 + 56735.224209104938*eccentricity_9;
    double common_term_55 = 2234230.984189969*eccentricity_14 - 778911.49955623181*eccentricity_12 + 120454.7790510995*eccentricity_10;
    double common_term_56 = 4903044.5475260969*eccentricity_15 - 1645757.5669389445*eccentricity_13 + 247360.22211577581*eccentricity_11;
    double common_term_57 = -3381063.2515870121*eccentricity_14 + 493629.74340239703*eccentricity_12;
    double common_term_58 = -6776841.5829965194*eccentricity_15 + 960872.32975856175*eccentricity_13;
    double common_term_59 = 1829985.8681430587*eccentricity_14;
    double common_term_60 = 3418541.8629036232*eccentricity_15;
    double common_term_61 = 240.14970281336967*eccentricity_15;
    double common_term_62 = 167.87385594393542*eccentricity_14;
    double common_term_63 = 392.06041777910674*eccentricity_15 + 117.05804771139515*eccentricity_13;
    double common_term_64 = 297.26168560929351*eccentricity_14 + 81.398393964272041*eccentricity_12;
    double common_term_65 = 614.67654218569095*eccentricity_15 + 223.46268072016878*eccentricity_13 + 56.426572107483566*eccentricity_11;
    double common_term_66 = 477.02526942850943*eccentricity_14 + 166.64708811313616*eccentricity_12 + 38.97926045547297*eccentricity_10;
    double common_term_67 = 880.57067982066761*eccentricity_15 + 367.89049043729708*eccentricity_13 + 123.33017299107143*eccentricity_11 + 26.819977678571429*eccentricity_9;
    double common_term_68 = 695.66809893049426*eccentricity_14 + 281.90556371022785*eccentricity_12 + 90.591962165092455*eccentricity_10 + 18.369782462952629*eccentricity_8;
    double common_term_69 = 1194.1469077681578*eccentricity_15 + 546.45176438032775*eccentricity_13 + 214.59179963073192*eccentricity_11 + 66.046688988095238*eccentricity_9 + 12.515674603174603*eccentricity_7;
    double common_term_70 = 956.09849668281419*eccentricity_14 + 426.7005069732666*eccentricity_12 + 162.23449042184012*eccentricity_10 + 47.78232421875*eccentricity_8 + 8.4744140625*eccentricity_6;
    double common_term_71 = 1557.6572694864694*eccentricity_15 + 761.43242399921002*eccentricity_13 + 331.13737495866402*eccentricity_11 + 121.77518601190476*eccentricity_9 + 34.289930555555556*eccentricity_7 + 5.6958333333333333*eccentricity_5;
    double common_term_72 = 1260.4874743123646*eccentricity_14 + 603.034425431963*eccentricity_12 + 255.31468111875827*eccentricity_10 + 90.715578884548611*eccentricity_8 + 24.39375*eccentricity_6 + 3.7942708333333333*eccentricity_4;
    double common_term_73 = std::pow(1.0 - eccentricity_2, -6.5)*(0.9375*eccentricity_5 + 2.5*eccentricity_3);
    double common_term_74 = 1611.0221868718007*eccentricity_14 + 812.85363100324358*eccentricity_12 + 371.54114950674552*eccentricity_10 + 148.61791449652778*eccentricity_8 + 49.0947265625*eccentricity_6 + 11.979166666666667*eccentricity_4 + 1.625*eccentricity_2;
    double common_term_75 = 2443.7033643131732*eccentricity_15 + 1308.9060652970679*eccentricity_13 + 647.41946903935185*eccentricity_11 + 288.82994791666667*eccentricity_9 + 112.08159722222222*eccentricity_7 + 35.604166666666667*eccentricity_5 + 8.25*eccentricity_3 + eccentricity;
    double common_term_76 = 2009.8835076297059*eccentricity_14 + 1058.097120475769*eccentricity_12 + 512.61032958984375*eccentricity_10 + 222.94708251953125*eccentricity_8 + 83.75*eccentricity_6 + 25.734375*eccentricity_4 + 5.0*eccentricity_2 + 1.0;
    double common_term_77 = 2970.8351335152116*eccentricity_15 + 1645.5110541363536*eccentricity_13 + 850.77243344907407*eccentricity_11 + 403.37942708333333*eccentricity_9 + 170.20138888888889*eccentricity_7 + 64.166666666666667*eccentricity_5 + 13.5*eccentricity_3 + 7.0*eccentricity;
    double common_term_78 = 2459.2330386569786*eccentricity_14 + 1340.5223589094858*eccentricity_12 + 681.12755528202763*eccentricity_10 + 310.51459418402778*eccentricity_8 + 143.6279296875*eccentricity_6 + 21.354166666666667*eccentricity_4 + 29.375*eccentricity_2;
    double common_term_79 = 3557.2417344447545*eccentricity_15 + 2025.3457993861607*eccentricity_13 + 1093.8272042410714*eccentricity_11 + 513.1875*eccentricity_9 + 312.009375*eccentricity_7 - 0.5*eccentricity_5 + 95.75*eccentricity_3;
    double common_term_80 = 2951.6210479807756*eccentricity_14 + 1703.3558808485667*eccentricity_12 + 748.01756275318287*eccentricity_10 + 694.08989529079861*eccentricity_8 - 146.54375*eccentricity_6 + 267.13802083333333*eccentricity_4;
    double common_term_81 = 4151.9040680338542*eccentricity_15 + 2645.1808312722755*eccentricity_13 + 851.20330067791005*eccentricity_11 + 1614.3615327380952*eccentricity_9 - 667.76736111111111*eccentricity_7 + 669.84583333333333*eccentricity_5;
    double common_term_82 = 4280.3905432988916*eccentricity_14 + 270.77507416861398*eccentricity_12 + 3882.5545790536063*eccentricity_10 - 2159.8021344866071*eccentricity_8 + 1553.2029296875*eccentricity_6;
    double common_term_83 = 7622.5562947786865*eccentricity_15 - 2563.8781164170433*eccentricity_13 + 9420.4653570050705*eccentricity_11 - 5934.5935639880952*eccentricity_9 + 3391.6381448412698*eccentricity_7;
    double common_term_84 = -11710.969000665865*eccentricity_14 + 22577.04078355224*eccentricity_12 - 14729.937145856154*eccentricity_10 + 7061.7649624294705*eccentricity_8;
    double common_term_85 = -36976.604685217127*eccentricity_15 + 52787.792419845779*eccentricity_13 - 34027.292165178571*eccentricity_11 + 14144.685491071429*eccentricity_9;
    double common_term_86 = 119791.53741267898*eccentricity_14 - 74444.509645218732*eccentricity_12 + 27434.965575509584*eccentricity_10;
    double common_term_87 = 263604.6276901502*eccentricity_15 - 155998.09765191703*eccentricity_13 + 51787.313793578143*eccentricity_11;
    double common_term_88 = -315572.51451485879*eccentricity_14 + 95510.178818067947*eccentricity_12;
    double common_term_89 = -619824.63616060102*eccentricity_15 + 172638.04917161679*eccentricity_13;
    double common_term_90 = 306605.38089750011*eccentricity_14;
    double common_term_91 = 536141.70405830693*eccentricity_15;
    double common_term_92 = 5203.3922597343835*eccentricity_15;
    double common_term_93 = 3410.1798952645757*eccentricity_14;
    double common_term_94 = 3582.9661959357607*eccentricity_15 + 2219.732540779533*eccentricity_13;
    double common_term_95 = 2790.4622432272421*eccentricity_14 + 1433.713227655718*eccentricity_12;
    double common_term_96 = 4997.5802377629363*eccentricity_15 + 2101.5557331776136*eccentricity_13 + 917.8695613438452*eccentricity_11;
    double common_term_97 = 3738.2614052149383*eccentricity_14 + 1539.4267026440509*eccentricity_12 + 581.64816264561244*eccentricity_10;
    double common_term_98 = 5858.9234321358913*eccentricity_15 + 2766.1920268032508*eccentricity_13 + 1100.3229152888007*eccentricity_11 + 364.20700920414462*eccentricity_9;
    double common_term_99 = 4442.8299867050879*eccentricity_14 + 2020.3823286063583*eccentricity_12 + 768.54706767517843*eccentricity_10 + 224.84113430447049*eccentricity_8;
    double common_term_100 = 6734.3571869926948*eccentricity_15 + 3326.6990248325893*eccentricity_13 + 1453.5822823660714*eccentricity_11 + 524.66685267857143*eccentricity_9 + 136.44642857142857*eccentricity_7;
    double common_term_101 = 5138.6539618359076*eccentricity_14 + 2456.4120699082453*eccentricity_12 + 1027.9645262824165*eccentricity_10 + 349.70782955109127*eccentricity_8 + 81.070203993055556*eccentricity_6;
    double common_term_102 = 7577.79759888944*eccentricity_15 + 3871.5873547499633*eccentricity_13 + 1785.5936549272487*eccentricity_11 + 712.79389880952381*eccentricity_9 + 227.03611111111111*eccentricity_7 + 46.891666666666667*eccentricity_5;
    double common_term_103 = 5805.1891895484924*eccentricity_14 + 2875.5705566065652*eccentricity_12 + 1274.9632638113839*eccentricity_10 + 483.0539794921875*eccentricity_8 + 142.96875*eccentricity_6 + 26.1796875*eccentricity_4;
    double common_term_104 = 8380.7718018764621*eccentricity_15 + 4389.2523771288029*eccentricity_13 + 2101.2032159391534*eccentricity_11 + 891.58909143518519*eccentricity_9 + 318.51458333333333*eccentricity_7 + 86.729166666666667*eccentricity_5 + 13.916666666666667*eccentricity_3;
    double common_term_105 = 6433.8736219646334*eccentricity_14 + 3269.3012639484708*eccentricity_12 + 1506.4300543325919*eccentricity_10 + 608.14258897569444*eccentricity_8 + 202.9990234375*eccentricity_6 + 50.104166666666667*eccentricity_4 + 6.875*eccentricity_2;
    double common_term_106 = std::pow(1.0 - eccentricity_2, -6.5)*(1.875*eccentricity_5 + 7.5*eccentricity_3 + 3.0*eccentricity);
    double common_term_107 = 7015.971802648293*eccentricity_14 + 3629.8719170666918*eccentricity_12 + 1715.6026451280382*eccentricity_10 + 719.42158338758681*eccentricity_8 + 255.51388888888889*eccentricity_6 + 70.859375*eccentricity_4 + 13.0*eccentricity_2 + 1.0;
    double common_term_108 = 9828.3969265850872*eccentricity_15 + 5310.0679312720459*eccentricity_13 + 2653.7806134259259*eccentricity_11 + 1198.9997395833333*eccentricity_9 + 472.82986111111111*eccentricity_7 + 153.52083333333333*eccentricity_5 + 36.75*eccentricity_3 + 5.0*eccentricity;
    double common_term_109 = 7542.7149553304059*eccentricity_14 + 3949.5133370753697*eccentricity_12 + 1895.7113754272461*eccentricity_10 + 811.11943359375*eccentricity_8 + 295.7548828125*eccentricity_6 + 84.5625*eccentricity_4 + 16.125*eccentricity_2;
    double common_term_110 = 10454.507961914254*eccentricity_15 + 5696.6772443369709*eccentricity_13 + 2876.2075603505291*eccentricity_11 + 1315.0491898148148*eccentricity_9 + 525.29166666666667*eccentricity_7 + 170.95833333333333*eccentricity_5 + 42.333333333333333*eccentricity_3;
    double common_term_111 = 8005.2878809465793*eccentricity_14 + 4220.4332463514237*eccentricity_12 + 2039.7442756128059*eccentricity_10 + 878.53856065538194*eccentricity_8 + 314.40625*eccentricity_6 + 98.294270833333333*eccentricity_4;
    double common_term_112 = 11003.026275669643*eccentricity_15 + 6023.2007603236607*eccentricity_13 + 3051.6074497767857*eccentricity_11 + 1403.6785714285714*eccentricity_9 + 534.825*eccentricity_7 + 210.15*eccentricity_5;
    double common_term_113 = 8396.2536505617281*eccentricity_14 + 4427.6258389429791*eccentricity_12 + 2167.3953261481391*eccentricity_10 + 846.33332248263889*eccentricity_8 + 423.28700086805556*eccentricity_6;
    double common_term_114 = 11472.022296503289*eccentricity_15 + 6250.12744575801*eccentricity_13 + 3270.1611951609347*eccentricity_11 + 1240.6909598214286*eccentricity_9 + 814.86884920634921*eccentricity_7;
    double common_term_115 = 8589.572218032936*eccentricity_14 + 4880.2275804110936*eccentricity_12 + 1653.4173156738281*eccentricity_10 + 1513.8415675571987*eccentricity_8;
    double common_term_116 = 11458.60758622292*eccentricity_15 + 7306.1606682133212*eccentricity_13 + 1898.8176756090168*eccentricity_11 + 2732.5824997244268*eccentricity_9;
    double common_term_117 = 11145.369200872075*eccentricity_14 + 1550.6867732646478*eccentricity_12 + 4816.5759872627595*eccentricity_10;
    double common_term_118 = 17577.292772305819*eccentricity_15 - 269.91084770698052*eccentricity_13 + 8321.8382508116883*eccentricity_11;
    double common_term_119 = -5263.1133069898864*eccentricity_14 + 14134.835157518386*eccentricity_12;
    double common_term_120 = -16577.762690309124*eccentricity_15 + 23657.161646581219*eccentricity_13;
    double common_term_121 = 39088.493294836992*eccentricity_14;
    double common_term_122 = 63858.064981999293*eccentricity_15;

    c_IntMap<c_Key1, double> result_by_q(31);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (7, 0).
    // q = -15
    result_by_lpq.set(c_Key3(7, 0, -15), common_term_0);
    result_by_q.set(c_Key1(-15), common_term_0);
    // q = -14
    result_by_lpq.set(c_Key3(7, 0, -14), common_term_1);
    result_by_q.set(c_Key1(-14), common_term_1);
    // q = -13
    result_by_lpq.set(c_Key3(7, 0, -13), common_term_2);
    result_by_q.set(c_Key1(-13), common_term_2);
    // q = -12
    result_by_lpq.set(c_Key3(7, 0, -12), common_term_3);
    result_by_q.set(c_Key1(-12), common_term_3);
    // q = -11
    result_by_lpq.set(c_Key3(7, 0, -11), common_term_4);
    result_by_q.set(c_Key1(-11), common_term_4);
    // q = -10
    result_by_lpq.set(c_Key3(7, 0, -10), common_term_5);
    result_by_q.set(c_Key1(-10), common_term_5);
    // q = -9
    result_by_lpq.set(c_Key3(7, 0, -9), common_term_6);
    result_by_q.set(c_Key1(-9), common_term_6);
    // q = -8
    result_by_lpq.set(c_Key3(7, 0, -8), common_term_7);
    result_by_q.set(c_Key1(-8), common_term_7);
    // q = -6
    result_by_lpq.set(c_Key3(7, 0, -6), common_term_8);
    result_by_q.set(c_Key1(-6), common_term_8);
    // q = -5
    result_by_lpq.set(c_Key3(7, 0, -5), common_term_9);
    result_by_q.set(c_Key1(-5), common_term_9);
    // q = -4
    result_by_lpq.set(c_Key3(7, 0, -4), common_term_10);
    result_by_q.set(c_Key1(-4), common_term_10);
    // q = -3
    result_by_lpq.set(c_Key3(7, 0, -3), common_term_11);
    result_by_q.set(c_Key1(-3), common_term_11);
    // q = -2
    result_by_lpq.set(c_Key3(7, 0, -2), common_term_12);
    result_by_q.set(c_Key1(-2), common_term_12);
    // q = -1
    result_by_lpq.set(c_Key3(7, 0, -1), common_term_13);
    result_by_q.set(c_Key1(-1), common_term_13);
    // q = 0
    result_by_lpq.set(c_Key3(7, 0, 0), common_term_14);
    result_by_q.set(c_Key1(0), common_term_14);
    // q = 1
    result_by_lpq.set(c_Key3(7, 0, 1), common_term_15);
    result_by_q.set(c_Key1(1), common_term_15);
    // q = 2
    result_by_lpq.set(c_Key3(7, 0, 2), common_term_16);
    result_by_q.set(c_Key1(2), common_term_16);
    // q = 3
    result_by_lpq.set(c_Key3(7, 0, 3), common_term_17);
    result_by_q.set(c_Key1(3), common_term_17);
    // q = 4
    result_by_lpq.set(c_Key3(7, 0, 4), common_term_18);
    result_by_q.set(c_Key1(4), common_term_18);
    // q = 5
    result_by_lpq.set(c_Key3(7, 0, 5), common_term_19);
    result_by_q.set(c_Key1(5), common_term_19);
    // q = 6
    result_by_lpq.set(c_Key3(7, 0, 6), common_term_20);
    result_by_q.set(c_Key1(6), common_term_20);
    // q = 7
    result_by_lpq.set(c_Key3(7, 0, 7), common_term_21);
    result_by_q.set(c_Key1(7), common_term_21);
    // q = 8
    result_by_lpq.set(c_Key3(7, 0, 8), common_term_22);
    result_by_q.set(c_Key1(8), common_term_22);
    // q = 9
    result_by_lpq.set(c_Key3(7, 0, 9), common_term_23);
    result_by_q.set(c_Key1(9), common_term_23);
    // q = 10
    result_by_lpq.set(c_Key3(7, 0, 10), common_term_24);
    result_by_q.set(c_Key1(10), common_term_24);
    // q = 11
    result_by_lpq.set(c_Key3(7, 0, 11), common_term_25);
    result_by_q.set(c_Key1(11), common_term_25);
    // q = 12
    result_by_lpq.set(c_Key3(7, 0, 12), common_term_26);
    result_by_q.set(c_Key1(12), common_term_26);
    // q = 13
    result_by_lpq.set(c_Key3(7, 0, 13), common_term_27);
    result_by_q.set(c_Key1(13), common_term_27);
    // q = 14
    result_by_lpq.set(c_Key3(7, 0, 14), common_term_28);
    result_by_q.set(c_Key1(14), common_term_28);
    // q = 15
    result_by_lpq.set(c_Key3(7, 0, 15), common_term_29);
    result_by_q.set(c_Key1(15), common_term_29);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 0), result_by_q);
    result_by_q.clear();

    // l , p = (7, 1).
    // q = -15
    result_by_lpq.set(c_Key3(7, 1, -15), common_term_30);
    result_by_q.set(c_Key1(-15), common_term_30);
    // q = -14
    result_by_lpq.set(c_Key3(7, 1, -14), common_term_31);
    result_by_q.set(c_Key1(-14), common_term_31);
    // q = -13
    result_by_lpq.set(c_Key3(7, 1, -13), common_term_32);
    result_by_q.set(c_Key1(-13), common_term_32);
    // q = -12
    result_by_lpq.set(c_Key3(7, 1, -12), common_term_33);
    result_by_q.set(c_Key1(-12), common_term_33);
    // q = -11
    result_by_lpq.set(c_Key3(7, 1, -11), common_term_34);
    result_by_q.set(c_Key1(-11), common_term_34);
    // q = -10
    result_by_lpq.set(c_Key3(7, 1, -10), common_term_35);
    result_by_q.set(c_Key1(-10), common_term_35);
    // q = -9
    result_by_lpq.set(c_Key3(7, 1, -9), common_term_36);
    result_by_q.set(c_Key1(-9), common_term_36);
    // q = -8
    result_by_lpq.set(c_Key3(7, 1, -8), common_term_37);
    result_by_q.set(c_Key1(-8), common_term_37);
    // q = -7
    result_by_lpq.set(c_Key3(7, 1, -7), common_term_38);
    result_by_q.set(c_Key1(-7), common_term_38);
    // q = -6
    result_by_lpq.set(c_Key3(7, 1, -6), common_term_39);
    result_by_q.set(c_Key1(-6), common_term_39);
    // q = -5
    result_by_lpq.set(c_Key3(7, 1, -5), common_term_40);
    result_by_q.set(c_Key1(-5), common_term_40);
    // q = -4
    result_by_lpq.set(c_Key3(7, 1, -4), common_term_41);
    result_by_q.set(c_Key1(-4), common_term_41);
    // q = -3
    result_by_lpq.set(c_Key3(7, 1, -3), common_term_42);
    result_by_q.set(c_Key1(-3), common_term_42);
    // q = -2
    result_by_lpq.set(c_Key3(7, 1, -2), common_term_43);
    result_by_q.set(c_Key1(-2), common_term_43);
    // q = -1
    result_by_lpq.set(c_Key3(7, 1, -1), common_term_44);
    result_by_q.set(c_Key1(-1), common_term_44);
    // q = 0
    result_by_lpq.set(c_Key3(7, 1, 0), common_term_45);
    result_by_q.set(c_Key1(0), common_term_45);
    // q = 1
    result_by_lpq.set(c_Key3(7, 1, 1), common_term_46);
    result_by_q.set(c_Key1(1), common_term_46);
    // q = 2
    result_by_lpq.set(c_Key3(7, 1, 2), common_term_47);
    result_by_q.set(c_Key1(2), common_term_47);
    // q = 3
    result_by_lpq.set(c_Key3(7, 1, 3), common_term_48);
    result_by_q.set(c_Key1(3), common_term_48);
    // q = 4
    result_by_lpq.set(c_Key3(7, 1, 4), common_term_49);
    result_by_q.set(c_Key1(4), common_term_49);
    // q = 5
    result_by_lpq.set(c_Key3(7, 1, 5), common_term_50);
    result_by_q.set(c_Key1(5), common_term_50);
    // q = 6
    result_by_lpq.set(c_Key3(7, 1, 6), common_term_51);
    result_by_q.set(c_Key1(6), common_term_51);
    // q = 7
    result_by_lpq.set(c_Key3(7, 1, 7), common_term_52);
    result_by_q.set(c_Key1(7), common_term_52);
    // q = 8
    result_by_lpq.set(c_Key3(7, 1, 8), common_term_53);
    result_by_q.set(c_Key1(8), common_term_53);
    // q = 9
    result_by_lpq.set(c_Key3(7, 1, 9), common_term_54);
    result_by_q.set(c_Key1(9), common_term_54);
    // q = 10
    result_by_lpq.set(c_Key3(7, 1, 10), common_term_55);
    result_by_q.set(c_Key1(10), common_term_55);
    // q = 11
    result_by_lpq.set(c_Key3(7, 1, 11), common_term_56);
    result_by_q.set(c_Key1(11), common_term_56);
    // q = 12
    result_by_lpq.set(c_Key3(7, 1, 12), common_term_57);
    result_by_q.set(c_Key1(12), common_term_57);
    // q = 13
    result_by_lpq.set(c_Key3(7, 1, 13), common_term_58);
    result_by_q.set(c_Key1(13), common_term_58);
    // q = 14
    result_by_lpq.set(c_Key3(7, 1, 14), common_term_59);
    result_by_q.set(c_Key1(14), common_term_59);
    // q = 15
    result_by_lpq.set(c_Key3(7, 1, 15), common_term_60);
    result_by_q.set(c_Key1(15), common_term_60);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 1), result_by_q);
    result_by_q.clear();

    // l , p = (7, 2).
    // q = -15
    result_by_lpq.set(c_Key3(7, 2, -15), common_term_61);
    result_by_q.set(c_Key1(-15), common_term_61);
    // q = -14
    result_by_lpq.set(c_Key3(7, 2, -14), common_term_62);
    result_by_q.set(c_Key1(-14), common_term_62);
    // q = -13
    result_by_lpq.set(c_Key3(7, 2, -13), common_term_63);
    result_by_q.set(c_Key1(-13), common_term_63);
    // q = -12
    result_by_lpq.set(c_Key3(7, 2, -12), common_term_64);
    result_by_q.set(c_Key1(-12), common_term_64);
    // q = -11
    result_by_lpq.set(c_Key3(7, 2, -11), common_term_65);
    result_by_q.set(c_Key1(-11), common_term_65);
    // q = -10
    result_by_lpq.set(c_Key3(7, 2, -10), common_term_66);
    result_by_q.set(c_Key1(-10), common_term_66);
    // q = -9
    result_by_lpq.set(c_Key3(7, 2, -9), common_term_67);
    result_by_q.set(c_Key1(-9), common_term_67);
    // q = -8
    result_by_lpq.set(c_Key3(7, 2, -8), common_term_68);
    result_by_q.set(c_Key1(-8), common_term_68);
    // q = -7
    result_by_lpq.set(c_Key3(7, 2, -7), common_term_69);
    result_by_q.set(c_Key1(-7), common_term_69);
    // q = -6
    result_by_lpq.set(c_Key3(7, 2, -6), common_term_70);
    result_by_q.set(c_Key1(-6), common_term_70);
    // q = -5
    result_by_lpq.set(c_Key3(7, 2, -5), common_term_71);
    result_by_q.set(c_Key1(-5), common_term_71);
    // q = -4
    result_by_lpq.set(c_Key3(7, 2, -4), common_term_72);
    result_by_q.set(c_Key1(-4), common_term_72);
    // q = -3
    result_by_lpq.set(c_Key3(7, 2, -3), common_term_73);
    result_by_q.set(c_Key1(-3), common_term_73);
    // q = -2
    result_by_lpq.set(c_Key3(7, 2, -2), common_term_74);
    result_by_q.set(c_Key1(-2), common_term_74);
    // q = -1
    result_by_lpq.set(c_Key3(7, 2, -1), common_term_75);
    result_by_q.set(c_Key1(-1), common_term_75);
    // q = 0
    result_by_lpq.set(c_Key3(7, 2, 0), common_term_76);
    result_by_q.set(c_Key1(0), common_term_76);
    // q = 1
    result_by_lpq.set(c_Key3(7, 2, 1), common_term_77);
    result_by_q.set(c_Key1(1), common_term_77);
    // q = 2
    result_by_lpq.set(c_Key3(7, 2, 2), common_term_78);
    result_by_q.set(c_Key1(2), common_term_78);
    // q = 3
    result_by_lpq.set(c_Key3(7, 2, 3), common_term_79);
    result_by_q.set(c_Key1(3), common_term_79);
    // q = 4
    result_by_lpq.set(c_Key3(7, 2, 4), common_term_80);
    result_by_q.set(c_Key1(4), common_term_80);
    // q = 5
    result_by_lpq.set(c_Key3(7, 2, 5), common_term_81);
    result_by_q.set(c_Key1(5), common_term_81);
    // q = 6
    result_by_lpq.set(c_Key3(7, 2, 6), common_term_82);
    result_by_q.set(c_Key1(6), common_term_82);
    // q = 7
    result_by_lpq.set(c_Key3(7, 2, 7), common_term_83);
    result_by_q.set(c_Key1(7), common_term_83);
    // q = 8
    result_by_lpq.set(c_Key3(7, 2, 8), common_term_84);
    result_by_q.set(c_Key1(8), common_term_84);
    // q = 9
    result_by_lpq.set(c_Key3(7, 2, 9), common_term_85);
    result_by_q.set(c_Key1(9), common_term_85);
    // q = 10
    result_by_lpq.set(c_Key3(7, 2, 10), common_term_86);
    result_by_q.set(c_Key1(10), common_term_86);
    // q = 11
    result_by_lpq.set(c_Key3(7, 2, 11), common_term_87);
    result_by_q.set(c_Key1(11), common_term_87);
    // q = 12
    result_by_lpq.set(c_Key3(7, 2, 12), common_term_88);
    result_by_q.set(c_Key1(12), common_term_88);
    // q = 13
    result_by_lpq.set(c_Key3(7, 2, 13), common_term_89);
    result_by_q.set(c_Key1(13), common_term_89);
    // q = 14
    result_by_lpq.set(c_Key3(7, 2, 14), common_term_90);
    result_by_q.set(c_Key1(14), common_term_90);
    // q = 15
    result_by_lpq.set(c_Key3(7, 2, 15), common_term_91);
    result_by_q.set(c_Key1(15), common_term_91);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 2), result_by_q);
    result_by_q.clear();

    // l , p = (7, 3).
    // q = -15
    result_by_lpq.set(c_Key3(7, 3, -15), common_term_92);
    result_by_q.set(c_Key1(-15), common_term_92);
    // q = -14
    result_by_lpq.set(c_Key3(7, 3, -14), common_term_93);
    result_by_q.set(c_Key1(-14), common_term_93);
    // q = -13
    result_by_lpq.set(c_Key3(7, 3, -13), common_term_94);
    result_by_q.set(c_Key1(-13), common_term_94);
    // q = -12
    result_by_lpq.set(c_Key3(7, 3, -12), common_term_95);
    result_by_q.set(c_Key1(-12), common_term_95);
    // q = -11
    result_by_lpq.set(c_Key3(7, 3, -11), common_term_96);
    result_by_q.set(c_Key1(-11), common_term_96);
    // q = -10
    result_by_lpq.set(c_Key3(7, 3, -10), common_term_97);
    result_by_q.set(c_Key1(-10), common_term_97);
    // q = -9
    result_by_lpq.set(c_Key3(7, 3, -9), common_term_98);
    result_by_q.set(c_Key1(-9), common_term_98);
    // q = -8
    result_by_lpq.set(c_Key3(7, 3, -8), common_term_99);
    result_by_q.set(c_Key1(-8), common_term_99);
    // q = -7
    result_by_lpq.set(c_Key3(7, 3, -7), common_term_100);
    result_by_q.set(c_Key1(-7), common_term_100);
    // q = -6
    result_by_lpq.set(c_Key3(7, 3, -6), common_term_101);
    result_by_q.set(c_Key1(-6), common_term_101);
    // q = -5
    result_by_lpq.set(c_Key3(7, 3, -5), common_term_102);
    result_by_q.set(c_Key1(-5), common_term_102);
    // q = -4
    result_by_lpq.set(c_Key3(7, 3, -4), common_term_103);
    result_by_q.set(c_Key1(-4), common_term_103);
    // q = -3
    result_by_lpq.set(c_Key3(7, 3, -3), common_term_104);
    result_by_q.set(c_Key1(-3), common_term_104);
    // q = -2
    result_by_lpq.set(c_Key3(7, 3, -2), common_term_105);
    result_by_q.set(c_Key1(-2), common_term_105);
    // q = -1
    result_by_lpq.set(c_Key3(7, 3, -1), common_term_106);
    result_by_q.set(c_Key1(-1), common_term_106);
    // q = 0
    result_by_lpq.set(c_Key3(7, 3, 0), common_term_107);
    result_by_q.set(c_Key1(0), common_term_107);
    // q = 1
    result_by_lpq.set(c_Key3(7, 3, 1), common_term_108);
    result_by_q.set(c_Key1(1), common_term_108);
    // q = 2
    result_by_lpq.set(c_Key3(7, 3, 2), common_term_109);
    result_by_q.set(c_Key1(2), common_term_109);
    // q = 3
    result_by_lpq.set(c_Key3(7, 3, 3), common_term_110);
    result_by_q.set(c_Key1(3), common_term_110);
    // q = 4
    result_by_lpq.set(c_Key3(7, 3, 4), common_term_111);
    result_by_q.set(c_Key1(4), common_term_111);
    // q = 5
    result_by_lpq.set(c_Key3(7, 3, 5), common_term_112);
    result_by_q.set(c_Key1(5), common_term_112);
    // q = 6
    result_by_lpq.set(c_Key3(7, 3, 6), common_term_113);
    result_by_q.set(c_Key1(6), common_term_113);
    // q = 7
    result_by_lpq.set(c_Key3(7, 3, 7), common_term_114);
    result_by_q.set(c_Key1(7), common_term_114);
    // q = 8
    result_by_lpq.set(c_Key3(7, 3, 8), common_term_115);
    result_by_q.set(c_Key1(8), common_term_115);
    // q = 9
    result_by_lpq.set(c_Key3(7, 3, 9), common_term_116);
    result_by_q.set(c_Key1(9), common_term_116);
    // q = 10
    result_by_lpq.set(c_Key3(7, 3, 10), common_term_117);
    result_by_q.set(c_Key1(10), common_term_117);
    // q = 11
    result_by_lpq.set(c_Key3(7, 3, 11), common_term_118);
    result_by_q.set(c_Key1(11), common_term_118);
    // q = 12
    result_by_lpq.set(c_Key3(7, 3, 12), common_term_119);
    result_by_q.set(c_Key1(12), common_term_119);
    // q = 13
    result_by_lpq.set(c_Key3(7, 3, 13), common_term_120);
    result_by_q.set(c_Key1(13), common_term_120);
    // q = 14
    result_by_lpq.set(c_Key3(7, 3, 14), common_term_121);
    result_by_q.set(c_Key1(14), common_term_121);
    // q = 15
    result_by_lpq.set(c_Key3(7, 3, 15), common_term_122);
    result_by_q.set(c_Key1(15), common_term_122);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 3), result_by_q);
    result_by_q.clear();

    // l , p = (7, 4).
    // q = -15
    result_by_lpq.set(c_Key3(7, 4, -15), common_term_122);
    result_by_q.set(c_Key1(-15), common_term_122);
    // q = -14
    result_by_lpq.set(c_Key3(7, 4, -14), common_term_121);
    result_by_q.set(c_Key1(-14), common_term_121);
    // q = -13
    result_by_lpq.set(c_Key3(7, 4, -13), common_term_120);
    result_by_q.set(c_Key1(-13), common_term_120);
    // q = -12
    result_by_lpq.set(c_Key3(7, 4, -12), common_term_119);
    result_by_q.set(c_Key1(-12), common_term_119);
    // q = -11
    result_by_lpq.set(c_Key3(7, 4, -11), common_term_118);
    result_by_q.set(c_Key1(-11), common_term_118);
    // q = -10
    result_by_lpq.set(c_Key3(7, 4, -10), common_term_117);
    result_by_q.set(c_Key1(-10), common_term_117);
    // q = -9
    result_by_lpq.set(c_Key3(7, 4, -9), common_term_116);
    result_by_q.set(c_Key1(-9), common_term_116);
    // q = -8
    result_by_lpq.set(c_Key3(7, 4, -8), common_term_115);
    result_by_q.set(c_Key1(-8), common_term_115);
    // q = -7
    result_by_lpq.set(c_Key3(7, 4, -7), common_term_114);
    result_by_q.set(c_Key1(-7), common_term_114);
    // q = -6
    result_by_lpq.set(c_Key3(7, 4, -6), common_term_113);
    result_by_q.set(c_Key1(-6), common_term_113);
    // q = -5
    result_by_lpq.set(c_Key3(7, 4, -5), common_term_112);
    result_by_q.set(c_Key1(-5), common_term_112);
    // q = -4
    result_by_lpq.set(c_Key3(7, 4, -4), common_term_111);
    result_by_q.set(c_Key1(-4), common_term_111);
    // q = -3
    result_by_lpq.set(c_Key3(7, 4, -3), common_term_110);
    result_by_q.set(c_Key1(-3), common_term_110);
    // q = -2
    result_by_lpq.set(c_Key3(7, 4, -2), common_term_109);
    result_by_q.set(c_Key1(-2), common_term_109);
    // q = -1
    result_by_lpq.set(c_Key3(7, 4, -1), common_term_108);
    result_by_q.set(c_Key1(-1), common_term_108);
    // q = 0
    result_by_lpq.set(c_Key3(7, 4, 0), common_term_107);
    result_by_q.set(c_Key1(0), common_term_107);
    // q = 1
    result_by_lpq.set(c_Key3(7, 4, 1), common_term_106);
    result_by_q.set(c_Key1(1), common_term_106);
    // q = 2
    result_by_lpq.set(c_Key3(7, 4, 2), common_term_105);
    result_by_q.set(c_Key1(2), common_term_105);
    // q = 3
    result_by_lpq.set(c_Key3(7, 4, 3), common_term_104);
    result_by_q.set(c_Key1(3), common_term_104);
    // q = 4
    result_by_lpq.set(c_Key3(7, 4, 4), common_term_103);
    result_by_q.set(c_Key1(4), common_term_103);
    // q = 5
    result_by_lpq.set(c_Key3(7, 4, 5), common_term_102);
    result_by_q.set(c_Key1(5), common_term_102);
    // q = 6
    result_by_lpq.set(c_Key3(7, 4, 6), common_term_101);
    result_by_q.set(c_Key1(6), common_term_101);
    // q = 7
    result_by_lpq.set(c_Key3(7, 4, 7), common_term_100);
    result_by_q.set(c_Key1(7), common_term_100);
    // q = 8
    result_by_lpq.set(c_Key3(7, 4, 8), common_term_99);
    result_by_q.set(c_Key1(8), common_term_99);
    // q = 9
    result_by_lpq.set(c_Key3(7, 4, 9), common_term_98);
    result_by_q.set(c_Key1(9), common_term_98);
    // q = 10
    result_by_lpq.set(c_Key3(7, 4, 10), common_term_97);
    result_by_q.set(c_Key1(10), common_term_97);
    // q = 11
    result_by_lpq.set(c_Key3(7, 4, 11), common_term_96);
    result_by_q.set(c_Key1(11), common_term_96);
    // q = 12
    result_by_lpq.set(c_Key3(7, 4, 12), common_term_95);
    result_by_q.set(c_Key1(12), common_term_95);
    // q = 13
    result_by_lpq.set(c_Key3(7, 4, 13), common_term_94);
    result_by_q.set(c_Key1(13), common_term_94);
    // q = 14
    result_by_lpq.set(c_Key3(7, 4, 14), common_term_93);
    result_by_q.set(c_Key1(14), common_term_93);
    // q = 15
    result_by_lpq.set(c_Key3(7, 4, 15), common_term_92);
    result_by_q.set(c_Key1(15), common_term_92);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 4), result_by_q);
    result_by_q.clear();

    // l , p = (7, 5).
    // q = -15
    result_by_lpq.set(c_Key3(7, 5, -15), common_term_91);
    result_by_q.set(c_Key1(-15), common_term_91);
    // q = -14
    result_by_lpq.set(c_Key3(7, 5, -14), common_term_90);
    result_by_q.set(c_Key1(-14), common_term_90);
    // q = -13
    result_by_lpq.set(c_Key3(7, 5, -13), common_term_89);
    result_by_q.set(c_Key1(-13), common_term_89);
    // q = -12
    result_by_lpq.set(c_Key3(7, 5, -12), common_term_88);
    result_by_q.set(c_Key1(-12), common_term_88);
    // q = -11
    result_by_lpq.set(c_Key3(7, 5, -11), common_term_87);
    result_by_q.set(c_Key1(-11), common_term_87);
    // q = -10
    result_by_lpq.set(c_Key3(7, 5, -10), common_term_86);
    result_by_q.set(c_Key1(-10), common_term_86);
    // q = -9
    result_by_lpq.set(c_Key3(7, 5, -9), common_term_85);
    result_by_q.set(c_Key1(-9), common_term_85);
    // q = -8
    result_by_lpq.set(c_Key3(7, 5, -8), common_term_84);
    result_by_q.set(c_Key1(-8), common_term_84);
    // q = -7
    result_by_lpq.set(c_Key3(7, 5, -7), common_term_83);
    result_by_q.set(c_Key1(-7), common_term_83);
    // q = -6
    result_by_lpq.set(c_Key3(7, 5, -6), common_term_82);
    result_by_q.set(c_Key1(-6), common_term_82);
    // q = -5
    result_by_lpq.set(c_Key3(7, 5, -5), common_term_81);
    result_by_q.set(c_Key1(-5), common_term_81);
    // q = -4
    result_by_lpq.set(c_Key3(7, 5, -4), common_term_80);
    result_by_q.set(c_Key1(-4), common_term_80);
    // q = -3
    result_by_lpq.set(c_Key3(7, 5, -3), common_term_79);
    result_by_q.set(c_Key1(-3), common_term_79);
    // q = -2
    result_by_lpq.set(c_Key3(7, 5, -2), common_term_78);
    result_by_q.set(c_Key1(-2), common_term_78);
    // q = -1
    result_by_lpq.set(c_Key3(7, 5, -1), common_term_77);
    result_by_q.set(c_Key1(-1), common_term_77);
    // q = 0
    result_by_lpq.set(c_Key3(7, 5, 0), common_term_76);
    result_by_q.set(c_Key1(0), common_term_76);
    // q = 1
    result_by_lpq.set(c_Key3(7, 5, 1), common_term_75);
    result_by_q.set(c_Key1(1), common_term_75);
    // q = 2
    result_by_lpq.set(c_Key3(7, 5, 2), common_term_74);
    result_by_q.set(c_Key1(2), common_term_74);
    // q = 3
    result_by_lpq.set(c_Key3(7, 5, 3), common_term_73);
    result_by_q.set(c_Key1(3), common_term_73);
    // q = 4
    result_by_lpq.set(c_Key3(7, 5, 4), common_term_72);
    result_by_q.set(c_Key1(4), common_term_72);
    // q = 5
    result_by_lpq.set(c_Key3(7, 5, 5), common_term_71);
    result_by_q.set(c_Key1(5), common_term_71);
    // q = 6
    result_by_lpq.set(c_Key3(7, 5, 6), common_term_70);
    result_by_q.set(c_Key1(6), common_term_70);
    // q = 7
    result_by_lpq.set(c_Key3(7, 5, 7), common_term_69);
    result_by_q.set(c_Key1(7), common_term_69);
    // q = 8
    result_by_lpq.set(c_Key3(7, 5, 8), common_term_68);
    result_by_q.set(c_Key1(8), common_term_68);
    // q = 9
    result_by_lpq.set(c_Key3(7, 5, 9), common_term_67);
    result_by_q.set(c_Key1(9), common_term_67);
    // q = 10
    result_by_lpq.set(c_Key3(7, 5, 10), common_term_66);
    result_by_q.set(c_Key1(10), common_term_66);
    // q = 11
    result_by_lpq.set(c_Key3(7, 5, 11), common_term_65);
    result_by_q.set(c_Key1(11), common_term_65);
    // q = 12
    result_by_lpq.set(c_Key3(7, 5, 12), common_term_64);
    result_by_q.set(c_Key1(12), common_term_64);
    // q = 13
    result_by_lpq.set(c_Key3(7, 5, 13), common_term_63);
    result_by_q.set(c_Key1(13), common_term_63);
    // q = 14
    result_by_lpq.set(c_Key3(7, 5, 14), common_term_62);
    result_by_q.set(c_Key1(14), common_term_62);
    // q = 15
    result_by_lpq.set(c_Key3(7, 5, 15), common_term_61);
    result_by_q.set(c_Key1(15), common_term_61);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 5), result_by_q);
    result_by_q.clear();

    // l , p = (7, 6).
    // q = -15
    result_by_lpq.set(c_Key3(7, 6, -15), common_term_60);
    result_by_q.set(c_Key1(-15), common_term_60);
    // q = -14
    result_by_lpq.set(c_Key3(7, 6, -14), common_term_59);
    result_by_q.set(c_Key1(-14), common_term_59);
    // q = -13
    result_by_lpq.set(c_Key3(7, 6, -13), common_term_58);
    result_by_q.set(c_Key1(-13), common_term_58);
    // q = -12
    result_by_lpq.set(c_Key3(7, 6, -12), common_term_57);
    result_by_q.set(c_Key1(-12), common_term_57);
    // q = -11
    result_by_lpq.set(c_Key3(7, 6, -11), common_term_56);
    result_by_q.set(c_Key1(-11), common_term_56);
    // q = -10
    result_by_lpq.set(c_Key3(7, 6, -10), common_term_55);
    result_by_q.set(c_Key1(-10), common_term_55);
    // q = -9
    result_by_lpq.set(c_Key3(7, 6, -9), common_term_54);
    result_by_q.set(c_Key1(-9), common_term_54);
    // q = -8
    result_by_lpq.set(c_Key3(7, 6, -8), common_term_53);
    result_by_q.set(c_Key1(-8), common_term_53);
    // q = -7
    result_by_lpq.set(c_Key3(7, 6, -7), common_term_52);
    result_by_q.set(c_Key1(-7), common_term_52);
    // q = -6
    result_by_lpq.set(c_Key3(7, 6, -6), common_term_51);
    result_by_q.set(c_Key1(-6), common_term_51);
    // q = -5
    result_by_lpq.set(c_Key3(7, 6, -5), common_term_50);
    result_by_q.set(c_Key1(-5), common_term_50);
    // q = -4
    result_by_lpq.set(c_Key3(7, 6, -4), common_term_49);
    result_by_q.set(c_Key1(-4), common_term_49);
    // q = -3
    result_by_lpq.set(c_Key3(7, 6, -3), common_term_48);
    result_by_q.set(c_Key1(-3), common_term_48);
    // q = -2
    result_by_lpq.set(c_Key3(7, 6, -2), common_term_47);
    result_by_q.set(c_Key1(-2), common_term_47);
    // q = -1
    result_by_lpq.set(c_Key3(7, 6, -1), common_term_46);
    result_by_q.set(c_Key1(-1), common_term_46);
    // q = 0
    result_by_lpq.set(c_Key3(7, 6, 0), common_term_45);
    result_by_q.set(c_Key1(0), common_term_45);
    // q = 1
    result_by_lpq.set(c_Key3(7, 6, 1), common_term_44);
    result_by_q.set(c_Key1(1), common_term_44);
    // q = 2
    result_by_lpq.set(c_Key3(7, 6, 2), common_term_43);
    result_by_q.set(c_Key1(2), common_term_43);
    // q = 3
    result_by_lpq.set(c_Key3(7, 6, 3), common_term_42);
    result_by_q.set(c_Key1(3), common_term_42);
    // q = 4
    result_by_lpq.set(c_Key3(7, 6, 4), common_term_41);
    result_by_q.set(c_Key1(4), common_term_41);
    // q = 5
    result_by_lpq.set(c_Key3(7, 6, 5), common_term_40);
    result_by_q.set(c_Key1(5), common_term_40);
    // q = 6
    result_by_lpq.set(c_Key3(7, 6, 6), common_term_39);
    result_by_q.set(c_Key1(6), common_term_39);
    // q = 7
    result_by_lpq.set(c_Key3(7, 6, 7), common_term_38);
    result_by_q.set(c_Key1(7), common_term_38);
    // q = 8
    result_by_lpq.set(c_Key3(7, 6, 8), common_term_37);
    result_by_q.set(c_Key1(8), common_term_37);
    // q = 9
    result_by_lpq.set(c_Key3(7, 6, 9), common_term_36);
    result_by_q.set(c_Key1(9), common_term_36);
    // q = 10
    result_by_lpq.set(c_Key3(7, 6, 10), common_term_35);
    result_by_q.set(c_Key1(10), common_term_35);
    // q = 11
    result_by_lpq.set(c_Key3(7, 6, 11), common_term_34);
    result_by_q.set(c_Key1(11), common_term_34);
    // q = 12
    result_by_lpq.set(c_Key3(7, 6, 12), common_term_33);
    result_by_q.set(c_Key1(12), common_term_33);
    // q = 13
    result_by_lpq.set(c_Key3(7, 6, 13), common_term_32);
    result_by_q.set(c_Key1(13), common_term_32);
    // q = 14
    result_by_lpq.set(c_Key3(7, 6, 14), common_term_31);
    result_by_q.set(c_Key1(14), common_term_31);
    // q = 15
    result_by_lpq.set(c_Key3(7, 6, 15), common_term_30);
    result_by_q.set(c_Key1(15), common_term_30);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 6), result_by_q);
    result_by_q.clear();

    // l , p = (7, 7).
    // q = -15
    result_by_lpq.set(c_Key3(7, 7, -15), common_term_29);
    result_by_q.set(c_Key1(-15), common_term_29);
    // q = -14
    result_by_lpq.set(c_Key3(7, 7, -14), common_term_28);
    result_by_q.set(c_Key1(-14), common_term_28);
    // q = -13
    result_by_lpq.set(c_Key3(7, 7, -13), common_term_27);
    result_by_q.set(c_Key1(-13), common_term_27);
    // q = -12
    result_by_lpq.set(c_Key3(7, 7, -12), common_term_26);
    result_by_q.set(c_Key1(-12), common_term_26);
    // q = -11
    result_by_lpq.set(c_Key3(7, 7, -11), common_term_25);
    result_by_q.set(c_Key1(-11), common_term_25);
    // q = -10
    result_by_lpq.set(c_Key3(7, 7, -10), common_term_24);
    result_by_q.set(c_Key1(-10), common_term_24);
    // q = -9
    result_by_lpq.set(c_Key3(7, 7, -9), common_term_23);
    result_by_q.set(c_Key1(-9), common_term_23);
    // q = -8
    result_by_lpq.set(c_Key3(7, 7, -8), common_term_22);
    result_by_q.set(c_Key1(-8), common_term_22);
    // q = -7
    result_by_lpq.set(c_Key3(7, 7, -7), common_term_21);
    result_by_q.set(c_Key1(-7), common_term_21);
    // q = -6
    result_by_lpq.set(c_Key3(7, 7, -6), common_term_20);
    result_by_q.set(c_Key1(-6), common_term_20);
    // q = -5
    result_by_lpq.set(c_Key3(7, 7, -5), common_term_19);
    result_by_q.set(c_Key1(-5), common_term_19);
    // q = -4
    result_by_lpq.set(c_Key3(7, 7, -4), common_term_18);
    result_by_q.set(c_Key1(-4), common_term_18);
    // q = -3
    result_by_lpq.set(c_Key3(7, 7, -3), common_term_17);
    result_by_q.set(c_Key1(-3), common_term_17);
    // q = -2
    result_by_lpq.set(c_Key3(7, 7, -2), common_term_16);
    result_by_q.set(c_Key1(-2), common_term_16);
    // q = -1
    result_by_lpq.set(c_Key3(7, 7, -1), common_term_15);
    result_by_q.set(c_Key1(-1), common_term_15);
    // q = 0
    result_by_lpq.set(c_Key3(7, 7, 0), common_term_14);
    result_by_q.set(c_Key1(0), common_term_14);
    // q = 1
    result_by_lpq.set(c_Key3(7, 7, 1), common_term_13);
    result_by_q.set(c_Key1(1), common_term_13);
    // q = 2
    result_by_lpq.set(c_Key3(7, 7, 2), common_term_12);
    result_by_q.set(c_Key1(2), common_term_12);
    // q = 3
    result_by_lpq.set(c_Key3(7, 7, 3), common_term_11);
    result_by_q.set(c_Key1(3), common_term_11);
    // q = 4
    result_by_lpq.set(c_Key3(7, 7, 4), common_term_10);
    result_by_q.set(c_Key1(4), common_term_10);
    // q = 5
    result_by_lpq.set(c_Key3(7, 7, 5), common_term_9);
    result_by_q.set(c_Key1(5), common_term_9);
    // q = 6
    result_by_lpq.set(c_Key3(7, 7, 6), common_term_8);
    result_by_q.set(c_Key1(6), common_term_8);
    // q = 8
    result_by_lpq.set(c_Key3(7, 7, 8), common_term_7);
    result_by_q.set(c_Key1(8), common_term_7);
    // q = 9
    result_by_lpq.set(c_Key3(7, 7, 9), common_term_6);
    result_by_q.set(c_Key1(9), common_term_6);
    // q = 10
    result_by_lpq.set(c_Key3(7, 7, 10), common_term_5);
    result_by_q.set(c_Key1(10), common_term_5);
    // q = 11
    result_by_lpq.set(c_Key3(7, 7, 11), common_term_4);
    result_by_q.set(c_Key1(11), common_term_4);
    // q = 12
    result_by_lpq.set(c_Key3(7, 7, 12), common_term_3);
    result_by_q.set(c_Key1(12), common_term_3);
    // q = 13
    result_by_lpq.set(c_Key3(7, 7, 13), common_term_2);
    result_by_q.set(c_Key1(13), common_term_2);
    // q = 14
    result_by_lpq.set(c_Key3(7, 7, 14), common_term_1);
    result_by_q.set(c_Key1(14), common_term_1);
    // q = 15
    result_by_lpq.set(c_Key3(7, 7, 15), common_term_0);
    result_by_q.set(c_Key1(15), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 7), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l7_e20(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 7.
    // Functions are taylor expanded and truncated at eccentricity^20.

    // Functions Calculated for l = 7.

    c_IntMap<c_Key3, double> result_by_lpq(326);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(8);
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
    double common_term_0 = 0.0074497640165616833*eccentricity_20;
    double common_term_1 = 0.0050093241566041212*eccentricity_19;
    double common_term_2 = 0.0072357222281565105*eccentricity_20 + 0.0033127402972282819*eccentricity_18;
    double common_term_3 = 0.0049453513274866363*eccentricity_19 + 0.0021449716601146856*eccentricity_17;
    double common_term_4 = 0.0055831773997595281*eccentricity_20 + 0.0032592341483802006*eccentricity_18 + 0.001351389768840571*eccentricity_16;
    double common_term_5 = 0.0035349904904837989*eccentricity_19 + 0.0020527698834577141*eccentricity_17 + 0.00082110795338308566*eccentricity_15;
    double common_term_6 = 0.0029999021403847899*eccentricity_20 + 0.0021020797681741731*eccentricity_18 + 0.0012187475040515692*eccentricity_16 + 0.00047483668989022175*eccentricity_14;
    double common_term_7 = 0.0016233522772567834*eccentricity_19 + 0.0011471194123956401*eccentricity_17 + 0.00066751468397673755*eccentricity_15 + 0.00025603302947052947*eccentricity_13;
    double common_term_8 = 0.00096396938704716771*eccentricity_20 + 0.00077223825500886885*eccentricity_18 + 0.00055327072227971856*eccentricity_16 + 0.00032544582513090068*eccentricity_14 + 0.0001244351684324032*eccentricity_12;
    double common_term_9 = 0.00036953824475870331*eccentricity_19 + 0.00030179305311580444*eccentricity_17 + 0.0002207668874335541*eccentricity_15 + 0.00013254235476457699*eccentricity_13 + 5.130671797338464e-5*eccentricity_11;
    double common_term_10 = 0.00011329047509849356*eccentricity_20 + 0.00010125527800483074*eccentricity_18 + 8.4869889588980051e-5*eccentricity_16 + 6.3981344173480938e-5*eccentricity_14 + 3.9727347237723214e-5*eccentricity_12 + 1.5890938895089286e-5*eccentricity_10;
    double common_term_11 = 1.5765472710175504e-5*eccentricity_19 + 1.4498454394649418e-5*eccentricity_17 + 1.2582073490638305e-5*eccentricity_15 + 9.8987143257976591e-6*eccentricity_13 + 6.4759700176366843e-6*eccentricity_11 + 2.7557319223985891e-6*eccentricity_9;
    double common_term_12 = 4.05690318420075e-7*eccentricity_20 + 3.9858815411425968e-7*eccentricity_18 + 3.8022437806952005e-7*eccentricity_16 + 3.4561838167659565e-7*eccentricity_14 + 2.8872616076595569e-7*eccentricity_12 + 2.0452697861552028e-7*eccentricity_10 + 9.6881200396825397e-8*eccentricity_8;
    double common_term_13 = 2.0788447505208615e-5*eccentricity_20 + 2.2277968130484441e-5*eccentricity_18 + 2.38281951086454e-5*eccentricity_16 + 2.5351391617933629e-5*eccentricity_14 + 2.6653318948986809e-5*eccentricity_12 + 2.7308388361855159e-5*eccentricity_10 + 2.6351686507936508e-5*eccentricity_8 + 2.1701388888888889e-5*eccentricity_6;
    double common_term_14 = -0.001339985982641182*eccentricity_19 - 0.0015275994337318845*eccentricity_17 - 0.0017635168650793651*eccentricity_15 - 0.0020673730526161082*eccentricity_13 - 0.0024698247354497354*eccentricity_11 - 0.0030505952380952381*eccentricity_9 - 0.0034722222222222222*eccentricity_7 - 0.0083333333333333333*eccentricity_5;
    double common_term_15 = -0.00782272139455873*eccentricity_20 - 0.008829798929600636*eccentricity_18 - 0.010077265454190118*eccentricity_16 - 0.011672703879220145*eccentricity_14 - 0.013487250464303153*eccentricity_12 - 0.019869559151785714*eccentricity_10 + 0.0085693359375*eccentricity_8 - 0.16875*eccentricity_6 + 0.2109375*eccentricity_4;
    double common_term_16 = -0.02765971769146889*eccentricity_19 - 0.031245347187228801*eccentricity_17 - 0.037010061483441113*eccentricity_15 - 0.028399470899470899*eccentricity_13 - 0.16982473544973545*eccentricity_11 + 0.64351851851851852*eccentricity_9 - 2.5166666666666667*eccentricity_7 + 3.6666666666666667*eccentricity_5 - 1.3333333333333333*eccentricity_3;
    double common_term_17 = -0.064944893405450522*eccentricity_20 - 0.069712035045643843*eccentricity_18 - 0.11888498588005268*eccentricity_16 + 0.22088288603846568*eccentricity_14 - 2.1656081790015811*eccentricity_12 + 8.9518299809208623*eccentricity_10 - 25.147840711805556*eccentricity_8 + 35.2783203125*eccentricity_6 - 19.270833333333333*eccentricity_4 + 3.125*eccentricity_2;
    double common_term_18 = -0.069705563292211416*eccentricity_19 - 0.80179276446906888*eccentricity_17 + 4.1779526068239796*eccentricity_15 - 22.384762834821429*eccentricity_13 + 78.94640625*eccentricity_11 - 183.90703125*eccentricity_9 + 242.25*eccentricity_7 - 152.8125*eccentricity_5 + 39.75*eccentricity_3 - 3.0*eccentricity;
    double common_term_19 = 0.82429296836700541*eccentricity_20 - 7.846995057115515*eccentricity_18 + 41.523197229126077*eccentricity_16 - 175.43335821740421*eccentricity_14 + 526.83805451334259*eccentricity_12 - 1068.9444772677951*eccentricity_10 + 1320.0667656792535*eccentricity_8 - 878.40277777777778*eccentricity_6 + 280.109375*eccentricity_4 - 35.0*eccentricity_2 + 1.0;
    double common_term_20 = -66.955996039265668*eccentricity_19 + 308.57640694110298*eccentricity_17 - 1101.9652262849584*eccentricity_15 + 2880.6147021053792*eccentricity_13 - 5223.6626851851852*eccentricity_11 + 6062.7453125*eccentricity_9 - 4088.2222222222222*eccentricity_7 + 1445.0833333333333*eccentricity_5 - 228.0*eccentricity_3 + 11.0*eccentricity;
    double common_term_21 = -468.35558967432224*eccentricity_20 + 1873.6874032978653*eccentricity_18 - 5830.01039704307*eccentricity_16 + 13530.386257162264*eccentricity_14 - 22294.568359745571*eccentricity_12 + 24396.094569396973*eccentricity_10 - 16337.53388671875*eccentricity_8 + 6092.2880859375*eccentricity_6 - 1094.0625*eccentricity_4 + 67.875*eccentricity_2;
    double common_term_22 = 9739.6559420274187*eccentricity_19 - 26925.385392315227*eccentricity_17 + 56343.268666686952*eccentricity_15 - 85355.961694465112*eccentricity_13 + 88336.957981977513*eccentricity_11 - 58113.017216435185*eccentricity_9 + 22236.979166666667*eccentricity_7 - 4311.9791666666667*eccentricity_5 + 309.58333333333333*eccentricity_3;
    double common_term_23 = 44698.440233928717*eccentricity_20 - 111385.53980795084*eccentricity_18 + 212705.07765675281*eccentricity_16 - 298897.11036031678*eccentricity_14 + 293421.13897960829*eccentricity_12 - 188463.82743972698*eccentricity_10 + 72772.538658311632*eccentricity_8 - 14777.60625*eccentricity_6 + 1163.0130208333333*eccentricity_4;
    double common_term_24 = -420648.70110032776*eccentricity_19 + 740077.34689548245*eccentricity_17 - 971396.39525223214*eccentricity_15 + 906983.55929129464*eccentricity_13 - 566805.24776785714*eccentricity_11 + 218518.31517857143*eccentricity_9 - 45579.75*eccentricity_7 + 3808.05*eccentricity_5;
    double common_term_25 = -1471245.339487014*eccentricity_20 + 2403215.8110022119*eccentricity_18 - 2963025.077554623*eccentricity_16 + 2637812.1504471274*eccentricity_14 - 1600816.3808199213*eccentricity_12 + 611884.91987441532*eccentricity_10 - 129425.82464812748*eccentricity_8 + 11244.611349826389*eccentricity_6;
    double common_term_26 = 7354992.78587202*eccentricity_19 - 8558083.3021652723*eccentricity_17 + 7280654.3524639735*eccentricity_15 - 4286530.2131702514*eccentricity_13 + 1616827.3261091821*eccentricity_11 - 343737.49739583333*eccentricity_9 + 30623.263194444444*eccentricity_7;
    double common_term_27 = 21381490.489301462*eccentricity_20 - 23572033.047721056*eccentricity_18 + 19203132.730778407*eccentricity_16 - 10964172.505874448*eccentricity_14 + 4067895.538355623*eccentricity_12 - 863790.37660871233*eccentricity_10 + 78130.892355782645*eccentricity_8;
    double common_term_28 = -62273459.574468785*eccentricity_19 + 48672043.125118458*eccentricity_17 - 26949151.751089544*eccentricity_15 + 9813445.3178285584*eccentricity_13 - 2071825.9252954145*eccentricity_11 + 188892.72068176808*eccentricity_9;
    double common_term_29 = -158549030.95351451*eccentricity_20 + 119095346.42392047*eccentricity_18 - 63961744.577118442*eccentricity_16 + 22825990.404013062*eccentricity_14 - 4775273.3536105321*eccentricity_12 + 436481.20207806919*eccentricity_10;
    double common_term_30 = 282414362.96419294*eccentricity_19 - 147176304.74830212*eccentricity_17 + 51421560.566231415*eccentricity_15 - 10633351.945539265*eccentricity_13 + 970435.77748782468*eccentricity_11;
    double common_term_31 = 651121714.22111634*eccentricity_20 - 329417993.79375457*eccentricity_18 + 112608650.39780295*eccentricity_16 - 22974539.901608997*eccentricity_14 + 2086930.9663047792*eccentricity_12;
    double common_term_32 = -719236687.86747912*eccentricity_19 + 240460756.42401798*eccentricity_17 - 48335585.125646251*eccentricity_15 + 4359471.0434132282*eccentricity_13;
    double common_term_33 = -1535514298.9872321*eccentricity_20 + 501980746.38575283*eccentricity_18 - 99313589.398729141*eccentricity_16 + 8876731.8003844841*eccentricity_14;
    double common_term_34 = 1026731497.3824746*eccentricity_19 - 199777085.96013138*eccentricity_17 + 17669216.908951237*eccentricity_15;
    double common_term_35 = 2061459051.1288307*eccentricity_20 - 394266001.14099111*eccentricity_18 + 34464816.460093673*eccentricity_16;
    double common_term_36 = -764751197.86319018*eccentricity_19 + 66011318.752456782*eccentricity_17;
    double common_term_37 = -1460199031.3862824*eccentricity_20 + 124366939.13016253*eccentricity_18;
    double common_term_38 = 230829268.25695985*eccentricity_19;
    double common_term_39 = 422615149.96673085*eccentricity_20;
    double common_term_40 = 19.108661039323087*eccentricity_20;
    double common_term_41 = 14.021684458148115*eccentricity_19;
    double common_term_42 = 32.995368420332185*eccentricity_20 + 10.289955157331042*eccentricity_18;
    double common_term_43 = 26.147402296492959*eccentricity_19 + 7.5522667924055412*eccentricity_17;
    double common_term_44 = 54.551802193821633*eccentricity_20 + 20.60999375400944*eccentricity_18 + 5.5437017946826046*eccentricity_16;
    double common_term_45 = 44.163960494891565*eccentricity_19 + 16.170450872025934*eccentricity_17 + 4.0699611066729458*eccentricity_15;
    double common_term_46 = 81.594532307041084*eccentricity_20 + 35.641666161619197*eccentricity_18 + 12.636417476653151*eccentricity_16 + 2.9885374554252938*eccentricity_14;
    double common_term_47 = 67.204523595108148*eccentricity_19 + 28.675555981772568*eccentricity_17 + 9.8401352562385394*eccentricity_15 + 2.1949078031423952*eccentricity_13;
    double common_term_48 = 114.9311351797822*eccentricity_20 + 55.198749497531008*eccentricity_18 + 23.002490343829614*eccentricity_16 + 7.6389896802115383*eccentricity_14 + 1.6124069097235637*eccentricity_12;
    double common_term_49 = 95.894543242545854*eccentricity_19 + 45.214731801750147*eccentricity_17 + 18.399222727955638*eccentricity_15 + 5.9140178571428571*eccentricity_13 + 1.1848006290584416*eccentricity_11;
    double common_term_50 = 155.09470589818828*eccentricity_20 + 79.813257484383327*eccentricity_18 + 36.93848950553172*eccentricity_16 + 14.677134880754364*eccentricity_14 + 4.5674260705064385*eccentricity_12 + 0.8708391879180962*eccentricity_10;
    double common_term_51 = 130.74858811596853*eccentricity_19 + 66.266771617129261*eccentricity_17 + 30.099295869088259*eccentricity_15 + 11.677661204054433*eccentricity_13 + 3.5197512951940035*eccentricity_11 + 0.64026262125220459*eccentricity_9;
    double common_term_52 = 202.63243040721222*eccentricity_20 + 109.97675184926306*eccentricity_18 + 54.88732269465633*eccentricity_16 + 24.464791772520387*eccentricity_14 + 9.2682209750584194*eccentricity_12 + 2.7070338657924107*eccentricity_10 + 0.47087140764508929*eccentricity_8;
    double common_term_53 = 172.28923143701832*eccentricity_19 + 92.299012026077609*eccentricity_17 + 45.354514222395082*eccentricity_15 + 19.836474259167402*eccentricity_13 + 7.3385771467151675*eccentricity_11 + 2.0781994047619048*eccentricity_9 + 0.34637896825396825*eccentricity_7;
    double common_term_54 = 258.09839088217321*eccentricity_20 + 146.18775619582771*eccentricity_18 + 77.292114673926577*eccentricity_16 + 37.390312413537874*eccentricity_14 + 16.045561849559356*eccentricity_12 + 5.7975862654428633*eccentricity_10 + 1.5927439856150794*eccentricity_8 + 0.25483940972222222*eccentricity_6;
    double common_term_55 = 0.1875*eccentricity_5*std::pow(1.0 - eccentricity_2, -6.5);
    double common_term_56 = 322.05607583236466*eccentricity_20 + 188.95391748084274*eccentricity_18 + 104.60469921857885*eccentricity_16 + 53.849767283159869*eccentricity_14 + 25.239834294243464*eccentricity_12 + 10.427123765087632*eccentricity_10 + 3.5954942491319444*eccentricity_8 + 0.93125*eccentricity_6 + 0.13802083333333333*eccentricity_4;
    double common_term_57 = 277.56168417608189*eccentricity_19 + 161.21232354007583*eccentricity_17 + 88.215787991974696*eccentricity_15 + 44.802796792328042*eccentricity_13 + 20.667006655092593*eccentricity_11 + 8.3757523148148148*eccentricity_9 + 2.8197916666666667*eccentricity_7 + 0.70833333333333333*eccentricity_5 + 0.083333333333333333*eccentricity_3;
    double common_term_58 = 395.07131042879974*eccentricity_20 + 238.78463626124548*eccentricity_18 + 137.27782991712191*eccentricity_16 + 74.238951664822442*eccentricity_14 + 37.189341234479632*eccentricity_12 + 16.876936340332031*eccentricity_10 + 6.68759765625*eccentricity_8 + 2.2763671875*eccentricity_6 + 0.1875*eccentricity_4 + 0.375*eccentricity_2;
    double common_term_59 = 342.36062329417912*eccentricity_19 + 205.04627276187067*eccentricity_17 + 116.65986192306784*eccentricity_15 + 62.328202229387125*eccentricity_13 + 30.817219328703704*eccentricity_11 + 13.431510416666667*eccentricity_9 + 6.6944444444444444*eccentricity_7 - 2.1666666666666667*eccentricity_5 + 4.5*eccentricity_3 - eccentricity;
    double common_term_60 = 477.68551260075419*eccentricity_20 + 296.15994787141646*eccentricity_18 + 175.73745565953657*eccentricity_16 + 98.832554347450831*eccentricity_14 + 52.795619199305405*eccentricity_12 + 21.909552680121528*eccentricity_10 + 23.896518283420139*eccentricity_8 - 24.444444444444444*eccentricity_6 + 31.484375*eccentricity_4 - 11.0*eccentricity_2 + 1.0;
    double common_term_61 = 415.92003275158941*eccentricity_19 + 255.84337773537149*eccentricity_17 + 149.22814378388074*eccentricity_15 + 89.920555245535714*eccentricity_13 + 15.6230859375*eccentricity_11 + 105.17109375*eccentricity_9 - 150.796875*eccentricity_7 + 162.9375*eccentricity_5 - 66.75*eccentricity_3 + 9.0*eccentricity;
    double common_term_62 = 570.20720781946667*eccentricity_20 + 363.00447291329779*eccentricity_18 + 210.7957761179853*eccentricity_16 + 174.9804609570827*eccentricity_14 - 100.78621590508355*eccentricity_12 + 474.60708239520038*eccentricity_10 - 708.84745008680556*eccentricity_8 + 687.3837890625*eccentricity_6 - 297.52083333333333*eccentricity_4 + 46.625*eccentricity_2;
    double common_term_63 = 511.30440320476695*eccentricity_19 + 246.79823144534667*eccentricity_17 + 465.78605284556113*eccentricity_15 - 752.55108129753638*eccentricity_13 + 1975.8828786375661*eccentricity_11 - 2801.8372106481481*eccentricity_9 + 2496.9395833333333*eccentricity_7 - 1090.0208333333333*eccentricity_5 + 182.16666666666667*eccentricity_3;
    double common_term_64 = 758.27215066420472*eccentricity_20 + 51.566391349264941*eccentricity_18 + 1641.5194412961474*eccentricity_16 - 3561.5581530353001*eccentricity_14 + 7431.9097426210131*eccentricity_12 - 9767.0560860770089*eccentricity_10 + 8088.3712646484375*eccentricity_8 - 3480.13125*eccentricity_6 + 595.7109375*eccentricity_4;
    double common_term_65 = -1295.2141398624953*eccentricity_19 + 6300.7613819967646*eccentricity_17 - 14047.947527096595*eccentricity_15 + 25445.885708027906*eccentricity_13 - 30894.854342344577*eccentricity_11 + 23931.665736607143*eccentricity_9 - 10027.289930555556*eccentricity_7 + 1720.3541666666667*eccentricity_5;
    double common_term_66 = -7396.9541872046814*eccentricity_20 + 23439.460248723672*eccentricity_18 - 49398.92761343685*eccentricity_16 + 80336.397336890789*eccentricity_14 - 90383.69587827253*eccentricity_12 + 65796.128116704547*eccentricity_10 - 26674.834478856647*eccentricity_8 + 4529.7294487847222*eccentricity_6;
    double common_term_67 = 81724.584809720187*eccentricity_19 - 159400.98469376839*eccentricity_17 + 236722.20799994927*eccentricity_15 - 247917.21170758929*eccentricity_13 + 170241.80323660714*eccentricity_11 - 66559.319866071429*eccentricity_9 + 11102.067857142857*eccentricity_7;
    double common_term_68 = 266412.03899170273*eccentricity_20 - 479885.02992708009*eccentricity_18 + 657652.02490062198*eccentricity_16 - 644137.50143032605*eccentricity_14 + 418575.6100948585*eccentricity_12 - 157569.04723622694*eccentricity_10 + 25694.067438712953*eccentricity_8;
    double common_term_69 = -1363144.0420207509*eccentricity_19 + 1736861.0409141105*eccentricity_17 - 1597867.3251102959*eccentricity_15 + 985361.68525131962*eccentricity_13 - 356951.15354841821*eccentricity_11 + 56735.224209104938*eccentricity_9;
    double common_term_70 = -3684104.0423023278*eccentricity_20 + 4389906.6732856117*eccentricity_18 - 3808142.5618735105*eccentricity_16 + 2234230.984189969*eccentricity_14 - 778911.49955623181*eccentricity_12 + 120454.7790510995*eccentricity_10;
    double common_term_71 = 10676926.858129445*eccentricity_19 - 8763779.6629679435*eccentricity_17 + 4903044.5475260969*eccentricity_15 - 1645757.5669389445*eccentricity_13 + 247360.22211577581*eccentricity_11;
    double common_term_72 = 25101654.813911323*eccentricity_20 - 19555787.444076166*eccentricity_18 + 10454980.127179817*eccentricity_16 - 3381063.2515870121*eccentricity_14 + 493629.74340239703*eccentricity_12;
    double common_term_73 = -42458036.657549949*eccentricity_19 + 21733055.53637045*eccentricity_17 - 6776841.5829965194*eccentricity_15 + 960872.32975856175*eccentricity_13;
    double common_term_74 = -89950061.783125983*eccentricity_20 + 44162026.531297567*eccentricity_18 - 13289506.750684231*eccentricity_16 + 1829985.8681430587*eccentricity_14;
    double common_term_75 = 87925586.823635785*eccentricity_19 - 25557550.319644213*eccentricity_17 + 3418541.8629036232*eccentricity_15;
    double common_term_76 = 171861983.53217287*eccentricity_20 - 48297185.961421776*eccentricity_18 + 6277092.8523520407*eccentricity_16;
    double common_term_77 = -89836611.761719483*eccentricity_19 + 11349368.763206806*eccentricity_17;
    double common_term_78 = -164720493.96085622*eccentricity_20 + 20236589.583097144*eccentricity_18;
    double common_term_79 = 35630330.744967336*eccentricity_19;
    double common_term_80 = 62016594.056112162*eccentricity_20;
    double common_term_81 = 1396.005316407883*eccentricity_20;
    double common_term_82 = 985.16583600618937*eccentricity_19;
    double common_term_83 = 1316.0554280324593*eccentricity_20 + 694.1335606152507*eccentricity_18;
    double common_term_84 = 1064.5289498914274*eccentricity_19 + 488.23549713752112*eccentricity_17;
    double common_term_85 = 2005.9141431910298*eccentricity_20 + 845.770632748252*eccentricity_18 + 342.76879721151412*eccentricity_16;
    double common_term_86 = 1599.1174291553588*eccentricity_19 + 662.25331999064342*eccentricity_17 + 240.14970281336967*eccentricity_15;
    double common_term_87 = 2640.1403812077651*eccentricity_20 + 1269.2820979712358*eccentricity_18 + 512.2539266171124*eccentricity_16 + 167.87385594393542*eccentricity_14;
    double common_term_88 = 2141.579596775477*eccentricity_19 + 1002.391144539429*eccentricity_17 + 392.06041777910674*eccentricity_15 + 117.05804771139515*eccentricity_13;
    double common_term_89 = 3371.6397342771172*eccentricity_20 + 1728.4057076891954*eccentricity_18 + 787.24994538544326*eccentricity_16 + 297.26168560929351*eccentricity_14 + 81.398393964272041*eccentricity_12;
    double common_term_90 = 2765.2102233143669*eccentricity_19 + 1387.7872562055079*eccentricity_17 + 614.67654218569095*eccentricity_15 + 223.46268072016878*eccentricity_13 + 56.426572107483566*eccentricity_11;
    double common_term_91 = 4190.6207263127266*eccentricity_20 + 2257.6363734484875*eccentricity_18 + 1108.4482112141791*eccentricity_16 + 477.02526942850943*eccentricity_14 + 166.64708811313616*eccentricity_12 + 38.97926045547297*eccentricity_10;
    double common_term_92 = 3467.8984320632763*eccentricity_19 + 1834.6803249588546*eccentricity_17 + 880.57067982066761*eccentricity_15 + 367.89049043729708*eccentricity_13 + 123.33017299107143*eccentricity_11 + 26.819977678571429*eccentricity_9;
    double common_term_93 = 5100.4052533810418*eccentricity_20 + 2857.8697377853513*eccentricity_18 + 1483.8316618364534*eccentricity_16 + 695.66809893049426*eccentricity_14 + 281.90556371022785*eccentricity_12 + 90.591962165092455*eccentricity_10 + 18.369782462952629*eccentricity_8;
    double common_term_94 = 4252.5387782371965*eccentricity_19 + 2345.021445261128*eccentricity_17 + 1194.1469077681578*eccentricity_15 + 546.45176438032775*eccentricity_13 + 214.59179963073192*eccentricity_11 + 66.046688988095238*eccentricity_9 + 12.515674603174603*eccentricity_7;
    double common_term_95 = 6103.8836830409883*eccentricity_20 + 3531.7912661342848*eccentricity_18 + 1915.6411337052862*eccentricity_16 + 956.09849668281419*eccentricity_14 + 426.7005069732666*eccentricity_12 + 162.23449042184012*eccentricity_10 + 47.78232421875*eccentricity_8 + 8.4744140625*eccentricity_6;
    double common_term_96 = 5121.9097056026083*eccentricity_19 + 2921.3554423049444*eccentricity_17 + 1557.6572694864694*eccentricity_15 + 761.43242399921002*eccentricity_13 + 331.13737495866402*eccentricity_11 + 121.77518601190476*eccentricity_9 + 34.289930555555556*eccentricity_7 + 5.6958333333333333*eccentricity_5;
    double common_term_97 = 7203.9603947628736*eccentricity_20 + 4282.0638643951933*eccentricity_18 + 2406.2996683787144*eccentricity_16 + 1260.4874743123646*eccentricity_14 + 603.034425431963*eccentricity_12 + 255.31468111875827*eccentricity_10 + 90.715578884548611*eccentricity_8 + 24.39375*eccentricity_6 + 3.7942708333333333*eccentricity_4;
    double common_term_98 = std::pow(1.0 - eccentricity_2, -6.5)*(0.9375*eccentricity_5 + 2.5*eccentricity_3);
    double common_term_99 = 8403.5449081973837*eccentricity_20 + 5111.3562312131424*eccentricity_18 + 2958.2351570436047*eccentricity_16 + 1611.0221868718007*eccentricity_14 + 812.85363100324358*eccentricity_12 + 371.54114950674552*eccentricity_10 + 148.61791449652778*eccentricity_8 + 49.0947265625*eccentricity_6 + 11.979166666666667*eccentricity_4 + 1.625*eccentricity_2;
    double common_term_100 = 7125.9934777098144*eccentricity_19 + 4282.1882443671561*eccentricity_17 + 2443.7033643131732*eccentricity_15 + 1308.9060652970679*eccentricity_13 + 647.41946903935185*eccentricity_11 + 288.82994791666667*eccentricity_9 + 112.08159722222222*eccentricity_7 + 35.604166666666667*eccentricity_5 + 8.25*eccentricity_3 + eccentricity;
    double common_term_101 = 9705.5424471876086*eccentricity_20 + 6022.3322955460763*eccentricity_18 + 3573.8700959307997*eccentricity_16 + 2009.8835076297059*eccentricity_14 + 1058.097120475769*eccentricity_12 + 512.61032958984375*eccentricity_10 + 222.94708251953125*eccentricity_8 + 83.75*eccentricity_6 + 25.734375*eccentricity_4 + 5.0*eccentricity_2 + 1.0;
    double common_term_102 = 8266.2699704060532*eccentricity_19 + 5071.7671783067259*eccentricity_17 + 2970.8351335152116*eccentricity_15 + 1645.5110541363536*eccentricity_13 + 850.77243344907407*eccentricity_11 + 403.37942708333333*eccentricity_9 + 170.20138888888889*eccentricity_7 + 64.166666666666667*eccentricity_5 + 13.5*eccentricity_3 + 7.0*eccentricity;
    double common_term_103 = 11112.835385523991*eccentricity_20 + 7017.6301814970689*eccentricity_18 + 4255.5959845365411*eccentricity_16 + 2459.2330386569786*eccentricity_14 + 1340.5223589094858*eccentricity_12 + 681.12755528202763*eccentricity_10 + 310.51459418402778*eccentricity_8 + 143.6279296875*eccentricity_6 + 21.354166666666667*eccentricity_4 + 29.375*eccentricity_2;
    double common_term_104 = 9502.3726549406883*eccentricity_19 + 5937.4395253682039*eccentricity_17 + 3557.2417344447545*eccentricity_15 + 2025.3457993861607*eccentricity_13 + 1093.8272042410714*eccentricity_11 + 513.1875*eccentricity_9 + 312.009375*eccentricity_7 - 0.5*eccentricity_5 + 95.75*eccentricity_3;
    double common_term_105 = 12628.297425678866*eccentricity_20 + 8099.6161128070384*eccentricity_18 + 5007.4366695639339*eccentricity_16 + 2951.6210479807756*eccentricity_14 + 1703.3558808485667*eccentricity_12 + 748.01756275318287*eccentricity_10 + 694.08989529079861*eccentricity_8 - 146.54375*eccentricity_6 + 267.13802083333333*eccentricity_4;
    double common_term_106 = 10835.130271498774*eccentricity_19 + 6892.8475169128045*eccentricity_17 + 4151.9040680338542*eccentricity_15 + 2645.1808312722755*eccentricity_13 + 851.20330067791005*eccentricity_11 + 1614.3615327380952*eccentricity_9 - 667.76736111111111*eccentricity_7 + 669.84583333333333*eccentricity_5;
    double common_term_107 = 14242.682828652588*eccentricity_20 + 9332.1184037124153*eccentricity_18 + 5582.0110597549615*eccentricity_16 + 4280.3905432988916*eccentricity_14 + 270.77507416861398*eccentricity_12 + 3882.5545790536063*eccentricity_10 - 2159.8021344866071*eccentricity_8 + 1553.2029296875*eccentricity_6;
    double common_term_108 = 12554.755334426377*eccentricity_19 + 6915.6997568193253*eccentricity_17 + 7622.5562947786865*eccentricity_15 - 2563.8781164170433*eccentricity_13 + 9420.4653570050705*eccentricity_11 - 5934.5935639880952*eccentricity_9 + 3391.6381448412698*eccentricity_7;
    double common_term_109 = 17152.295725584421*eccentricity_20 + 6938.1457936004047*eccentricity_18 + 15449.936936074537*eccentricity_16 - 11710.969000665865*eccentricity_14 + 22577.04078355224*eccentricity_12 - 14729.937145856154*eccentricity_10 + 7061.7649624294705*eccentricity_8;
    double common_term_110 = 1886.3711924251983*eccentricity_19 + 35033.411718934142*eccentricity_17 - 36976.604685217127*eccentricity_15 + 52787.792419845779*eccentricity_13 - 34027.292165178571*eccentricity_11 + 14144.685491071429*eccentricity_9;
    double common_term_111 = -18820.168298203932*eccentricity_20 + 84503.32794716156*eccentricity_18 - 100868.96418440037*eccentricity_16 + 119791.53741267898*eccentricity_14 - 74444.509645218732*eccentricity_12 + 27434.965575509584*eccentricity_10;
    double common_term_112 = 207106.14854518932*eccentricity_19 - 252947.18749015286*eccentricity_17 + 263604.6276901502*eccentricity_15 - 155998.09765191703*eccentricity_13 + 51787.313793578143*eccentricity_11;
    double common_term_113 = 502026.27534325335*eccentricity_20 - 598604.86799481414*eccentricity_18 + 563177.93610851022*eccentricity_16 - 315572.51451485879*eccentricity_14 + 95510.178818067947*eccentricity_12;
    double common_term_114 = -1355423.4035969743*eccentricity_19 + 1170639.0856813026*eccentricity_17 - 619824.63616060102*eccentricity_15 + 172638.04917161679*eccentricity_13;
    double common_term_115 = -2961191.8926603362*eccentricity_20 + 2373103.6846468836*eccentricity_18 - 1187201.2631509788*eccentricity_16 + 306605.38089750011*eccentricity_14;
    double common_term_116 = 4702625.541307624*eccentricity_19 - 2225090.3424417321*eccentricity_17 + 536141.70405830693*eccentricity_15;
    double common_term_117 = 9129234.4211811879*eccentricity_20 - 4091887.9251178695*eccentricity_18 + 924665.71197961586*eccentricity_16;
    double common_term_118 = -7399765.4894607641*eccentricity_19 + 1575161.1375469564*eccentricity_17;
    double common_term_119 = -13183448.757692872*eccentricity_20 + 2653597.628259475*eccentricity_18;
    double common_term_120 = 4425609.5630584049*eccentricity_19;
    double common_term_121 = 7313656.2983455319*eccentricity_20;
    double common_term_122 = 39733.41520120026*eccentricity_20;
    double common_term_123 = 26697.72247954273*eccentricity_19;
    double common_term_124 = 1137.1697236964025*eccentricity_20 + 17866.240247366696*eccentricity_18;
    double common_term_125 = 4317.4537902775023*eccentricity_19 + 11903.106750708065*eccentricity_17;
    double common_term_126 = 20646.412190686783*eccentricity_20 + 5259.1265813227759*eccentricity_18 + 7891.4336377999706*eccentricity_16;
    double common_term_127 = 15336.783662062483*eccentricity_19 + 5076.2714385486112*eccentricity_17 + 5203.3922597343835*eccentricity_15;
    double common_term_128 = 19767.43838852038*eccentricity_20 + 11547.485902300369*eccentricity_18 + 4402.7894838852435*eccentricity_16 + 3410.1798952645757*eccentricity_14;
    double common_term_129 = 15853.30316603491*eccentricity_19 + 8747.7002225011317*eccentricity_17 + 3582.9661959357607*eccentricity_15 + 2219.732540779533*eccentricity_13;
    double common_term_130 = 22353.836904011141*eccentricity_20 + 12562.318039643584*eccentricity_18 + 6627.1371559397548*eccentricity_16 + 2790.4622432272421*eccentricity_14 + 1433.713227655718*eccentricity_12;
    double common_term_131 = 17894.457712204871*eccentricity_19 + 9846.9707732594357*eccentricity_17 + 4997.5802377629363*eccentricity_15 + 2101.5557331776136*eccentricity_13 + 917.8695613438452*eccentricity_11;
    double common_term_132 = 24614.285891297382*eccentricity_20 + 14211.216808144955*eccentricity_18 + 7637.132778292789*eccentricity_16 + 3738.2614052149383*eccentricity_14 + 1539.4267026440509*eccentricity_12 + 581.64816264561244*eccentricity_10;
    double common_term_133 = 19791.702699385434*eccentricity_19 + 11188.360771111126*eccentricity_17 + 5858.9234321358913*eccentricity_15 + 2766.1920268032508*eccentricity_13 + 1100.3229152888007*eccentricity_11 + 364.20700920414462*eccentricity_9;
    double common_term_134 = 26826.92721848866*eccentricity_20 + 15786.707194277842*eccentricity_18 + 8725.359496784696*eccentricity_16 + 4442.8299867050879*eccentricity_14 + 2020.3823286063583*eccentricity_12 + 768.54706767517843*eccentricity_10 + 224.84113430447049*eccentricity_8;
    double common_term_135 = 21638.247397958792*eccentricity_19 + 12482.909206963081*eccentricity_17 + 6734.3571869926948*eccentricity_15 + 3326.6990248325893*eccentricity_13 + 1453.5822823660714*eccentricity_11 + 524.66685267857143*eccentricity_9 + 136.44642857142857*eccentricity_7;
    double common_term_136 = 28970.755755186163*eccentricity_20 + 17313.154677481836*eccentricity_18 + 9777.0074335308644*eccentricity_16 + 5138.6539618359076*eccentricity_14 + 2456.4120699082453*eccentricity_12 + 1027.9645262824165*eccentricity_10 + 349.70782955109127*eccentricity_8 + 81.070203993055556*eccentricity_6;
    double common_term_137 = 23420.392334131672*eccentricity_19 + 13731.234963925086*eccentricity_17 + 7577.79759888944*eccentricity_15 + 3871.5873547499633*eccentricity_13 + 1785.5936549272487*eccentricity_11 + 712.79389880952381*eccentricity_9 + 227.03611111111111*eccentricity_7 + 46.891666666666667*eccentricity_5;
    double common_term_138 = 31034.180120538619*eccentricity_20 + 18779.267544457304*eccentricity_18 + 10785.381936045949*eccentricity_16 + 5805.1891895484924*eccentricity_14 + 2875.5705566065652*eccentricity_12 + 1274.9632638113839*eccentricity_10 + 483.0539794921875*eccentricity_8 + 142.96875*eccentricity_6 + 26.1796875*eccentricity_4;
    double common_term_139 = 25126.954966084754*eccentricity_19 + 14923.015060615198*eccentricity_17 + 8380.7718018764621*eccentricity_15 + 4389.2523771288029*eccentricity_13 + 2101.2032159391534*eccentricity_11 + 891.58909143518519*eccentricity_9 + 318.51458333333333*eccentricity_7 + 86.729166666666667*eccentricity_5 + 13.916666666666667*eccentricity_3;
    double common_term_140 = 33005.489455747568*eccentricity_20 + 20174.337868841967*eccentricity_18 + 11740.75207279902*eccentricity_16 + 6433.8736219646334*eccentricity_14 + 3269.3012639484708*eccentricity_12 + 1506.4300543325919*eccentricity_10 + 608.14258897569444*eccentricity_8 + 202.9990234375*eccentricity_6 + 50.104166666666667*eccentricity_4 + 6.875*eccentricity_2;
    double common_term_141 = std::pow(1.0 - eccentricity_2, -6.5)*(1.875*eccentricity_5 + 7.5*eccentricity_3 + 3.0*eccentricity);
    double common_term_142 = 34872.947007682995*eccentricity_20 + 21487.629499097101*eccentricity_18 + 12633.382245929482*eccentricity_16 + 7015.971802648293*eccentricity_14 + 3629.8719170666918*eccentricity_12 + 1715.6026451280382*eccentricity_10 + 719.42158338758681*eccentricity_8 + 255.51388888888889*eccentricity_6 + 70.859375*eccentricity_4 + 13.0*eccentricity_2 + 1.0;
    double common_term_143 = 28268.413489724639*eccentricity_19 + 17096.021174203857*eccentricity_17 + 9828.3969265850872*eccentricity_15 + 5310.0679312720459*eccentricity_13 + 2653.7806134259259*eccentricity_11 + 1198.9997395833333*eccentricity_9 + 472.82986111111111*eccentricity_7 + 153.52083333333333*eccentricity_5 + 36.75*eccentricity_3 + 5.0*eccentricity;
    double common_term_144 = 36624.784734106642*eccentricity_20 + 22708.374202354765*eccentricity_18 + 13453.503984032152*eccentricity_16 + 7542.7149553304059*eccentricity_14 + 3949.5133370753697*eccentricity_12 + 1895.7113754272461*eccentricity_10 + 811.11943359375*eccentricity_8 + 295.7548828125*eccentricity_6 + 84.5625*eccentricity_4 + 16.125*eccentricity_2;
    double common_term_145 = 29680.771140870576*eccentricity_19 + 18056.708222861153*eccentricity_17 + 10454.507961914254*eccentricity_15 + 5696.6772443369709*eccentricity_13 + 2876.2075603505291*eccentricity_11 + 1315.0491898148148*eccentricity_9 + 525.29166666666667*eccentricity_7 + 170.95833333333333*eccentricity_5 + 42.333333333333333*eccentricity_3;
    double common_term_146 = 38249.196587171398*eccentricity_20 + 23825.763942354959*eccentricity_18 + 14191.307149393959*eccentricity_16 + 8005.2878809465793*eccentricity_14 + 4220.4332463514237*eccentricity_12 + 2039.7442756128059*eccentricity_10 + 878.53856065538194*eccentricity_8 + 314.40625*eccentricity_6 + 98.294270833333333*eccentricity_4;
    double common_term_147 = 30972.457143244982*eccentricity_19 + 18919.759881720209*eccentricity_17 + 11003.026275669643*eccentricity_15 + 6023.2007603236607*eccentricity_13 + 3051.6074497767857*eccentricity_11 + 1403.6785714285714*eccentricity_9 + 534.825*eccentricity_7 + 210.15*eccentricity_5;
    double common_term_148 = 39734.331985411195*eccentricity_20 + 24828.974004303339*eccentricity_18 + 14836.71395837454*eccentricity_16 + 8396.2536505617281*eccentricity_14 + 4427.6258389429791*eccentricity_12 + 2167.3953261481391*eccentricity_10 + 846.33332248263889*eccentricity_8 + 423.28700086805556*eccentricity_6;
    double common_term_149 = 32132.309314332181*eccentricity_19 + 19673.402706744412*eccentricity_17 + 11472.022296503289*eccentricity_15 + 6250.12744575801*eccentricity_13 + 3270.1611951609347*eccentricity_11 + 1240.6909598214286*eccentricity_9 + 814.86884920634921*eccentricity_7;
    double common_term_150 = 41069.5862725569*eccentricity_20 + 25699.88114063068*eccentricity_18 + 15412.616649217572*eccentricity_16 + 8589.572218032936*eccentricity_14 + 4880.2275804110936*eccentricity_12 + 1653.4173156738281*eccentricity_10 + 1513.8415675571987*eccentricity_8;
    double common_term_151 = 33116.847668802292*eccentricity_19 + 20433.448696257571*eccentricity_17 + 11458.60758622292*eccentricity_15 + 7306.1606682133212*eccentricity_13 + 1898.8176756090168*eccentricity_11 + 2732.5824997244268*eccentricity_9;
    double common_term_152 = 42117.155979466461*eccentricity_20 + 26864.457051724452*eccentricity_18 + 14701.560823456953*eccentricity_16 + 11145.369200872075*eccentricity_14 + 1550.6867732646478*eccentricity_12 + 4816.5759872627595*eccentricity_10;
    double common_term_153 = 35301.763182580986*eccentricity_19 + 17746.136699861941*eccentricity_17 + 17577.292772305819*eccentricity_15 - 269.91084770698052*eccentricity_13 + 8321.8382508116883*eccentricity_11;
    double common_term_154 = 46965.445056193638*eccentricity_20 + 19072.843433047585*eccentricity_18 + 28924.275651635131*eccentricity_16 - 5263.1133069898864*eccentricity_14 + 14134.835157518386*eccentricity_12;
    double common_term_155 = 15125.277009332411*eccentricity_19 + 49696.309867380975*eccentricity_17 - 16577.762690309124*eccentricity_15 + 23657.161646581219*eccentricity_13;
    double common_term_156 = -1856.5609026329713*eccentricity_20 + 88493.353107952137*eccentricity_18 - 39863.158298257686*eccentricity_16 + 39088.493294836992*eccentricity_14;
    double common_term_157 = 161402.17329041829*eccentricity_19 - 85008.477656959714*eccentricity_17 + 63858.064981999293*eccentricity_15;
    double common_term_158 = 297960.08166834292*eccentricity_20 - 168986.74379281262*eccentricity_18 + 103279.80384113961*eccentricity_16;
    double common_term_159 = -320458.741508488*eccentricity_19 + 165543.09372395997*eccentricity_17;
    double common_term_160 = -587161.16185042373*eccentricity_20 + 263205.61254751027*eccentricity_18;
    double common_term_161 = 415434.99856760201*eccentricity_19;
    double common_term_162 = 651364.28261275118*eccentricity_20;

    c_IntMap<c_Key1, double> result_by_q(41);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (7, 0).
    // q = -20
    result_by_lpq.set(c_Key3(7, 0, -20), common_term_0);
    result_by_q.set(c_Key1(-20), common_term_0);
    // q = -19
    result_by_lpq.set(c_Key3(7, 0, -19), common_term_1);
    result_by_q.set(c_Key1(-19), common_term_1);
    // q = -18
    result_by_lpq.set(c_Key3(7, 0, -18), common_term_2);
    result_by_q.set(c_Key1(-18), common_term_2);
    // q = -17
    result_by_lpq.set(c_Key3(7, 0, -17), common_term_3);
    result_by_q.set(c_Key1(-17), common_term_3);
    // q = -16
    result_by_lpq.set(c_Key3(7, 0, -16), common_term_4);
    result_by_q.set(c_Key1(-16), common_term_4);
    // q = -15
    result_by_lpq.set(c_Key3(7, 0, -15), common_term_5);
    result_by_q.set(c_Key1(-15), common_term_5);
    // q = -14
    result_by_lpq.set(c_Key3(7, 0, -14), common_term_6);
    result_by_q.set(c_Key1(-14), common_term_6);
    // q = -13
    result_by_lpq.set(c_Key3(7, 0, -13), common_term_7);
    result_by_q.set(c_Key1(-13), common_term_7);
    // q = -12
    result_by_lpq.set(c_Key3(7, 0, -12), common_term_8);
    result_by_q.set(c_Key1(-12), common_term_8);
    // q = -11
    result_by_lpq.set(c_Key3(7, 0, -11), common_term_9);
    result_by_q.set(c_Key1(-11), common_term_9);
    // q = -10
    result_by_lpq.set(c_Key3(7, 0, -10), common_term_10);
    result_by_q.set(c_Key1(-10), common_term_10);
    // q = -9
    result_by_lpq.set(c_Key3(7, 0, -9), common_term_11);
    result_by_q.set(c_Key1(-9), common_term_11);
    // q = -8
    result_by_lpq.set(c_Key3(7, 0, -8), common_term_12);
    result_by_q.set(c_Key1(-8), common_term_12);
    // q = -6
    result_by_lpq.set(c_Key3(7, 0, -6), common_term_13);
    result_by_q.set(c_Key1(-6), common_term_13);
    // q = -5
    result_by_lpq.set(c_Key3(7, 0, -5), common_term_14);
    result_by_q.set(c_Key1(-5), common_term_14);
    // q = -4
    result_by_lpq.set(c_Key3(7, 0, -4), common_term_15);
    result_by_q.set(c_Key1(-4), common_term_15);
    // q = -3
    result_by_lpq.set(c_Key3(7, 0, -3), common_term_16);
    result_by_q.set(c_Key1(-3), common_term_16);
    // q = -2
    result_by_lpq.set(c_Key3(7, 0, -2), common_term_17);
    result_by_q.set(c_Key1(-2), common_term_17);
    // q = -1
    result_by_lpq.set(c_Key3(7, 0, -1), common_term_18);
    result_by_q.set(c_Key1(-1), common_term_18);
    // q = 0
    result_by_lpq.set(c_Key3(7, 0, 0), common_term_19);
    result_by_q.set(c_Key1(0), common_term_19);
    // q = 1
    result_by_lpq.set(c_Key3(7, 0, 1), common_term_20);
    result_by_q.set(c_Key1(1), common_term_20);
    // q = 2
    result_by_lpq.set(c_Key3(7, 0, 2), common_term_21);
    result_by_q.set(c_Key1(2), common_term_21);
    // q = 3
    result_by_lpq.set(c_Key3(7, 0, 3), common_term_22);
    result_by_q.set(c_Key1(3), common_term_22);
    // q = 4
    result_by_lpq.set(c_Key3(7, 0, 4), common_term_23);
    result_by_q.set(c_Key1(4), common_term_23);
    // q = 5
    result_by_lpq.set(c_Key3(7, 0, 5), common_term_24);
    result_by_q.set(c_Key1(5), common_term_24);
    // q = 6
    result_by_lpq.set(c_Key3(7, 0, 6), common_term_25);
    result_by_q.set(c_Key1(6), common_term_25);
    // q = 7
    result_by_lpq.set(c_Key3(7, 0, 7), common_term_26);
    result_by_q.set(c_Key1(7), common_term_26);
    // q = 8
    result_by_lpq.set(c_Key3(7, 0, 8), common_term_27);
    result_by_q.set(c_Key1(8), common_term_27);
    // q = 9
    result_by_lpq.set(c_Key3(7, 0, 9), common_term_28);
    result_by_q.set(c_Key1(9), common_term_28);
    // q = 10
    result_by_lpq.set(c_Key3(7, 0, 10), common_term_29);
    result_by_q.set(c_Key1(10), common_term_29);
    // q = 11
    result_by_lpq.set(c_Key3(7, 0, 11), common_term_30);
    result_by_q.set(c_Key1(11), common_term_30);
    // q = 12
    result_by_lpq.set(c_Key3(7, 0, 12), common_term_31);
    result_by_q.set(c_Key1(12), common_term_31);
    // q = 13
    result_by_lpq.set(c_Key3(7, 0, 13), common_term_32);
    result_by_q.set(c_Key1(13), common_term_32);
    // q = 14
    result_by_lpq.set(c_Key3(7, 0, 14), common_term_33);
    result_by_q.set(c_Key1(14), common_term_33);
    // q = 15
    result_by_lpq.set(c_Key3(7, 0, 15), common_term_34);
    result_by_q.set(c_Key1(15), common_term_34);
    // q = 16
    result_by_lpq.set(c_Key3(7, 0, 16), common_term_35);
    result_by_q.set(c_Key1(16), common_term_35);
    // q = 17
    result_by_lpq.set(c_Key3(7, 0, 17), common_term_36);
    result_by_q.set(c_Key1(17), common_term_36);
    // q = 18
    result_by_lpq.set(c_Key3(7, 0, 18), common_term_37);
    result_by_q.set(c_Key1(18), common_term_37);
    // q = 19
    result_by_lpq.set(c_Key3(7, 0, 19), common_term_38);
    result_by_q.set(c_Key1(19), common_term_38);
    // q = 20
    result_by_lpq.set(c_Key3(7, 0, 20), common_term_39);
    result_by_q.set(c_Key1(20), common_term_39);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 0), result_by_q);
    result_by_q.clear();

    // l , p = (7, 1).
    // q = -20
    result_by_lpq.set(c_Key3(7, 1, -20), common_term_40);
    result_by_q.set(c_Key1(-20), common_term_40);
    // q = -19
    result_by_lpq.set(c_Key3(7, 1, -19), common_term_41);
    result_by_q.set(c_Key1(-19), common_term_41);
    // q = -18
    result_by_lpq.set(c_Key3(7, 1, -18), common_term_42);
    result_by_q.set(c_Key1(-18), common_term_42);
    // q = -17
    result_by_lpq.set(c_Key3(7, 1, -17), common_term_43);
    result_by_q.set(c_Key1(-17), common_term_43);
    // q = -16
    result_by_lpq.set(c_Key3(7, 1, -16), common_term_44);
    result_by_q.set(c_Key1(-16), common_term_44);
    // q = -15
    result_by_lpq.set(c_Key3(7, 1, -15), common_term_45);
    result_by_q.set(c_Key1(-15), common_term_45);
    // q = -14
    result_by_lpq.set(c_Key3(7, 1, -14), common_term_46);
    result_by_q.set(c_Key1(-14), common_term_46);
    // q = -13
    result_by_lpq.set(c_Key3(7, 1, -13), common_term_47);
    result_by_q.set(c_Key1(-13), common_term_47);
    // q = -12
    result_by_lpq.set(c_Key3(7, 1, -12), common_term_48);
    result_by_q.set(c_Key1(-12), common_term_48);
    // q = -11
    result_by_lpq.set(c_Key3(7, 1, -11), common_term_49);
    result_by_q.set(c_Key1(-11), common_term_49);
    // q = -10
    result_by_lpq.set(c_Key3(7, 1, -10), common_term_50);
    result_by_q.set(c_Key1(-10), common_term_50);
    // q = -9
    result_by_lpq.set(c_Key3(7, 1, -9), common_term_51);
    result_by_q.set(c_Key1(-9), common_term_51);
    // q = -8
    result_by_lpq.set(c_Key3(7, 1, -8), common_term_52);
    result_by_q.set(c_Key1(-8), common_term_52);
    // q = -7
    result_by_lpq.set(c_Key3(7, 1, -7), common_term_53);
    result_by_q.set(c_Key1(-7), common_term_53);
    // q = -6
    result_by_lpq.set(c_Key3(7, 1, -6), common_term_54);
    result_by_q.set(c_Key1(-6), common_term_54);
    // q = -5
    result_by_lpq.set(c_Key3(7, 1, -5), common_term_55);
    result_by_q.set(c_Key1(-5), common_term_55);
    // q = -4
    result_by_lpq.set(c_Key3(7, 1, -4), common_term_56);
    result_by_q.set(c_Key1(-4), common_term_56);
    // q = -3
    result_by_lpq.set(c_Key3(7, 1, -3), common_term_57);
    result_by_q.set(c_Key1(-3), common_term_57);
    // q = -2
    result_by_lpq.set(c_Key3(7, 1, -2), common_term_58);
    result_by_q.set(c_Key1(-2), common_term_58);
    // q = -1
    result_by_lpq.set(c_Key3(7, 1, -1), common_term_59);
    result_by_q.set(c_Key1(-1), common_term_59);
    // q = 0
    result_by_lpq.set(c_Key3(7, 1, 0), common_term_60);
    result_by_q.set(c_Key1(0), common_term_60);
    // q = 1
    result_by_lpq.set(c_Key3(7, 1, 1), common_term_61);
    result_by_q.set(c_Key1(1), common_term_61);
    // q = 2
    result_by_lpq.set(c_Key3(7, 1, 2), common_term_62);
    result_by_q.set(c_Key1(2), common_term_62);
    // q = 3
    result_by_lpq.set(c_Key3(7, 1, 3), common_term_63);
    result_by_q.set(c_Key1(3), common_term_63);
    // q = 4
    result_by_lpq.set(c_Key3(7, 1, 4), common_term_64);
    result_by_q.set(c_Key1(4), common_term_64);
    // q = 5
    result_by_lpq.set(c_Key3(7, 1, 5), common_term_65);
    result_by_q.set(c_Key1(5), common_term_65);
    // q = 6
    result_by_lpq.set(c_Key3(7, 1, 6), common_term_66);
    result_by_q.set(c_Key1(6), common_term_66);
    // q = 7
    result_by_lpq.set(c_Key3(7, 1, 7), common_term_67);
    result_by_q.set(c_Key1(7), common_term_67);
    // q = 8
    result_by_lpq.set(c_Key3(7, 1, 8), common_term_68);
    result_by_q.set(c_Key1(8), common_term_68);
    // q = 9
    result_by_lpq.set(c_Key3(7, 1, 9), common_term_69);
    result_by_q.set(c_Key1(9), common_term_69);
    // q = 10
    result_by_lpq.set(c_Key3(7, 1, 10), common_term_70);
    result_by_q.set(c_Key1(10), common_term_70);
    // q = 11
    result_by_lpq.set(c_Key3(7, 1, 11), common_term_71);
    result_by_q.set(c_Key1(11), common_term_71);
    // q = 12
    result_by_lpq.set(c_Key3(7, 1, 12), common_term_72);
    result_by_q.set(c_Key1(12), common_term_72);
    // q = 13
    result_by_lpq.set(c_Key3(7, 1, 13), common_term_73);
    result_by_q.set(c_Key1(13), common_term_73);
    // q = 14
    result_by_lpq.set(c_Key3(7, 1, 14), common_term_74);
    result_by_q.set(c_Key1(14), common_term_74);
    // q = 15
    result_by_lpq.set(c_Key3(7, 1, 15), common_term_75);
    result_by_q.set(c_Key1(15), common_term_75);
    // q = 16
    result_by_lpq.set(c_Key3(7, 1, 16), common_term_76);
    result_by_q.set(c_Key1(16), common_term_76);
    // q = 17
    result_by_lpq.set(c_Key3(7, 1, 17), common_term_77);
    result_by_q.set(c_Key1(17), common_term_77);
    // q = 18
    result_by_lpq.set(c_Key3(7, 1, 18), common_term_78);
    result_by_q.set(c_Key1(18), common_term_78);
    // q = 19
    result_by_lpq.set(c_Key3(7, 1, 19), common_term_79);
    result_by_q.set(c_Key1(19), common_term_79);
    // q = 20
    result_by_lpq.set(c_Key3(7, 1, 20), common_term_80);
    result_by_q.set(c_Key1(20), common_term_80);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 1), result_by_q);
    result_by_q.clear();

    // l , p = (7, 2).
    // q = -20
    result_by_lpq.set(c_Key3(7, 2, -20), common_term_81);
    result_by_q.set(c_Key1(-20), common_term_81);
    // q = -19
    result_by_lpq.set(c_Key3(7, 2, -19), common_term_82);
    result_by_q.set(c_Key1(-19), common_term_82);
    // q = -18
    result_by_lpq.set(c_Key3(7, 2, -18), common_term_83);
    result_by_q.set(c_Key1(-18), common_term_83);
    // q = -17
    result_by_lpq.set(c_Key3(7, 2, -17), common_term_84);
    result_by_q.set(c_Key1(-17), common_term_84);
    // q = -16
    result_by_lpq.set(c_Key3(7, 2, -16), common_term_85);
    result_by_q.set(c_Key1(-16), common_term_85);
    // q = -15
    result_by_lpq.set(c_Key3(7, 2, -15), common_term_86);
    result_by_q.set(c_Key1(-15), common_term_86);
    // q = -14
    result_by_lpq.set(c_Key3(7, 2, -14), common_term_87);
    result_by_q.set(c_Key1(-14), common_term_87);
    // q = -13
    result_by_lpq.set(c_Key3(7, 2, -13), common_term_88);
    result_by_q.set(c_Key1(-13), common_term_88);
    // q = -12
    result_by_lpq.set(c_Key3(7, 2, -12), common_term_89);
    result_by_q.set(c_Key1(-12), common_term_89);
    // q = -11
    result_by_lpq.set(c_Key3(7, 2, -11), common_term_90);
    result_by_q.set(c_Key1(-11), common_term_90);
    // q = -10
    result_by_lpq.set(c_Key3(7, 2, -10), common_term_91);
    result_by_q.set(c_Key1(-10), common_term_91);
    // q = -9
    result_by_lpq.set(c_Key3(7, 2, -9), common_term_92);
    result_by_q.set(c_Key1(-9), common_term_92);
    // q = -8
    result_by_lpq.set(c_Key3(7, 2, -8), common_term_93);
    result_by_q.set(c_Key1(-8), common_term_93);
    // q = -7
    result_by_lpq.set(c_Key3(7, 2, -7), common_term_94);
    result_by_q.set(c_Key1(-7), common_term_94);
    // q = -6
    result_by_lpq.set(c_Key3(7, 2, -6), common_term_95);
    result_by_q.set(c_Key1(-6), common_term_95);
    // q = -5
    result_by_lpq.set(c_Key3(7, 2, -5), common_term_96);
    result_by_q.set(c_Key1(-5), common_term_96);
    // q = -4
    result_by_lpq.set(c_Key3(7, 2, -4), common_term_97);
    result_by_q.set(c_Key1(-4), common_term_97);
    // q = -3
    result_by_lpq.set(c_Key3(7, 2, -3), common_term_98);
    result_by_q.set(c_Key1(-3), common_term_98);
    // q = -2
    result_by_lpq.set(c_Key3(7, 2, -2), common_term_99);
    result_by_q.set(c_Key1(-2), common_term_99);
    // q = -1
    result_by_lpq.set(c_Key3(7, 2, -1), common_term_100);
    result_by_q.set(c_Key1(-1), common_term_100);
    // q = 0
    result_by_lpq.set(c_Key3(7, 2, 0), common_term_101);
    result_by_q.set(c_Key1(0), common_term_101);
    // q = 1
    result_by_lpq.set(c_Key3(7, 2, 1), common_term_102);
    result_by_q.set(c_Key1(1), common_term_102);
    // q = 2
    result_by_lpq.set(c_Key3(7, 2, 2), common_term_103);
    result_by_q.set(c_Key1(2), common_term_103);
    // q = 3
    result_by_lpq.set(c_Key3(7, 2, 3), common_term_104);
    result_by_q.set(c_Key1(3), common_term_104);
    // q = 4
    result_by_lpq.set(c_Key3(7, 2, 4), common_term_105);
    result_by_q.set(c_Key1(4), common_term_105);
    // q = 5
    result_by_lpq.set(c_Key3(7, 2, 5), common_term_106);
    result_by_q.set(c_Key1(5), common_term_106);
    // q = 6
    result_by_lpq.set(c_Key3(7, 2, 6), common_term_107);
    result_by_q.set(c_Key1(6), common_term_107);
    // q = 7
    result_by_lpq.set(c_Key3(7, 2, 7), common_term_108);
    result_by_q.set(c_Key1(7), common_term_108);
    // q = 8
    result_by_lpq.set(c_Key3(7, 2, 8), common_term_109);
    result_by_q.set(c_Key1(8), common_term_109);
    // q = 9
    result_by_lpq.set(c_Key3(7, 2, 9), common_term_110);
    result_by_q.set(c_Key1(9), common_term_110);
    // q = 10
    result_by_lpq.set(c_Key3(7, 2, 10), common_term_111);
    result_by_q.set(c_Key1(10), common_term_111);
    // q = 11
    result_by_lpq.set(c_Key3(7, 2, 11), common_term_112);
    result_by_q.set(c_Key1(11), common_term_112);
    // q = 12
    result_by_lpq.set(c_Key3(7, 2, 12), common_term_113);
    result_by_q.set(c_Key1(12), common_term_113);
    // q = 13
    result_by_lpq.set(c_Key3(7, 2, 13), common_term_114);
    result_by_q.set(c_Key1(13), common_term_114);
    // q = 14
    result_by_lpq.set(c_Key3(7, 2, 14), common_term_115);
    result_by_q.set(c_Key1(14), common_term_115);
    // q = 15
    result_by_lpq.set(c_Key3(7, 2, 15), common_term_116);
    result_by_q.set(c_Key1(15), common_term_116);
    // q = 16
    result_by_lpq.set(c_Key3(7, 2, 16), common_term_117);
    result_by_q.set(c_Key1(16), common_term_117);
    // q = 17
    result_by_lpq.set(c_Key3(7, 2, 17), common_term_118);
    result_by_q.set(c_Key1(17), common_term_118);
    // q = 18
    result_by_lpq.set(c_Key3(7, 2, 18), common_term_119);
    result_by_q.set(c_Key1(18), common_term_119);
    // q = 19
    result_by_lpq.set(c_Key3(7, 2, 19), common_term_120);
    result_by_q.set(c_Key1(19), common_term_120);
    // q = 20
    result_by_lpq.set(c_Key3(7, 2, 20), common_term_121);
    result_by_q.set(c_Key1(20), common_term_121);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 2), result_by_q);
    result_by_q.clear();

    // l , p = (7, 3).
    // q = -20
    result_by_lpq.set(c_Key3(7, 3, -20), common_term_122);
    result_by_q.set(c_Key1(-20), common_term_122);
    // q = -19
    result_by_lpq.set(c_Key3(7, 3, -19), common_term_123);
    result_by_q.set(c_Key1(-19), common_term_123);
    // q = -18
    result_by_lpq.set(c_Key3(7, 3, -18), common_term_124);
    result_by_q.set(c_Key1(-18), common_term_124);
    // q = -17
    result_by_lpq.set(c_Key3(7, 3, -17), common_term_125);
    result_by_q.set(c_Key1(-17), common_term_125);
    // q = -16
    result_by_lpq.set(c_Key3(7, 3, -16), common_term_126);
    result_by_q.set(c_Key1(-16), common_term_126);
    // q = -15
    result_by_lpq.set(c_Key3(7, 3, -15), common_term_127);
    result_by_q.set(c_Key1(-15), common_term_127);
    // q = -14
    result_by_lpq.set(c_Key3(7, 3, -14), common_term_128);
    result_by_q.set(c_Key1(-14), common_term_128);
    // q = -13
    result_by_lpq.set(c_Key3(7, 3, -13), common_term_129);
    result_by_q.set(c_Key1(-13), common_term_129);
    // q = -12
    result_by_lpq.set(c_Key3(7, 3, -12), common_term_130);
    result_by_q.set(c_Key1(-12), common_term_130);
    // q = -11
    result_by_lpq.set(c_Key3(7, 3, -11), common_term_131);
    result_by_q.set(c_Key1(-11), common_term_131);
    // q = -10
    result_by_lpq.set(c_Key3(7, 3, -10), common_term_132);
    result_by_q.set(c_Key1(-10), common_term_132);
    // q = -9
    result_by_lpq.set(c_Key3(7, 3, -9), common_term_133);
    result_by_q.set(c_Key1(-9), common_term_133);
    // q = -8
    result_by_lpq.set(c_Key3(7, 3, -8), common_term_134);
    result_by_q.set(c_Key1(-8), common_term_134);
    // q = -7
    result_by_lpq.set(c_Key3(7, 3, -7), common_term_135);
    result_by_q.set(c_Key1(-7), common_term_135);
    // q = -6
    result_by_lpq.set(c_Key3(7, 3, -6), common_term_136);
    result_by_q.set(c_Key1(-6), common_term_136);
    // q = -5
    result_by_lpq.set(c_Key3(7, 3, -5), common_term_137);
    result_by_q.set(c_Key1(-5), common_term_137);
    // q = -4
    result_by_lpq.set(c_Key3(7, 3, -4), common_term_138);
    result_by_q.set(c_Key1(-4), common_term_138);
    // q = -3
    result_by_lpq.set(c_Key3(7, 3, -3), common_term_139);
    result_by_q.set(c_Key1(-3), common_term_139);
    // q = -2
    result_by_lpq.set(c_Key3(7, 3, -2), common_term_140);
    result_by_q.set(c_Key1(-2), common_term_140);
    // q = -1
    result_by_lpq.set(c_Key3(7, 3, -1), common_term_141);
    result_by_q.set(c_Key1(-1), common_term_141);
    // q = 0
    result_by_lpq.set(c_Key3(7, 3, 0), common_term_142);
    result_by_q.set(c_Key1(0), common_term_142);
    // q = 1
    result_by_lpq.set(c_Key3(7, 3, 1), common_term_143);
    result_by_q.set(c_Key1(1), common_term_143);
    // q = 2
    result_by_lpq.set(c_Key3(7, 3, 2), common_term_144);
    result_by_q.set(c_Key1(2), common_term_144);
    // q = 3
    result_by_lpq.set(c_Key3(7, 3, 3), common_term_145);
    result_by_q.set(c_Key1(3), common_term_145);
    // q = 4
    result_by_lpq.set(c_Key3(7, 3, 4), common_term_146);
    result_by_q.set(c_Key1(4), common_term_146);
    // q = 5
    result_by_lpq.set(c_Key3(7, 3, 5), common_term_147);
    result_by_q.set(c_Key1(5), common_term_147);
    // q = 6
    result_by_lpq.set(c_Key3(7, 3, 6), common_term_148);
    result_by_q.set(c_Key1(6), common_term_148);
    // q = 7
    result_by_lpq.set(c_Key3(7, 3, 7), common_term_149);
    result_by_q.set(c_Key1(7), common_term_149);
    // q = 8
    result_by_lpq.set(c_Key3(7, 3, 8), common_term_150);
    result_by_q.set(c_Key1(8), common_term_150);
    // q = 9
    result_by_lpq.set(c_Key3(7, 3, 9), common_term_151);
    result_by_q.set(c_Key1(9), common_term_151);
    // q = 10
    result_by_lpq.set(c_Key3(7, 3, 10), common_term_152);
    result_by_q.set(c_Key1(10), common_term_152);
    // q = 11
    result_by_lpq.set(c_Key3(7, 3, 11), common_term_153);
    result_by_q.set(c_Key1(11), common_term_153);
    // q = 12
    result_by_lpq.set(c_Key3(7, 3, 12), common_term_154);
    result_by_q.set(c_Key1(12), common_term_154);
    // q = 13
    result_by_lpq.set(c_Key3(7, 3, 13), common_term_155);
    result_by_q.set(c_Key1(13), common_term_155);
    // q = 14
    result_by_lpq.set(c_Key3(7, 3, 14), common_term_156);
    result_by_q.set(c_Key1(14), common_term_156);
    // q = 15
    result_by_lpq.set(c_Key3(7, 3, 15), common_term_157);
    result_by_q.set(c_Key1(15), common_term_157);
    // q = 16
    result_by_lpq.set(c_Key3(7, 3, 16), common_term_158);
    result_by_q.set(c_Key1(16), common_term_158);
    // q = 17
    result_by_lpq.set(c_Key3(7, 3, 17), common_term_159);
    result_by_q.set(c_Key1(17), common_term_159);
    // q = 18
    result_by_lpq.set(c_Key3(7, 3, 18), common_term_160);
    result_by_q.set(c_Key1(18), common_term_160);
    // q = 19
    result_by_lpq.set(c_Key3(7, 3, 19), common_term_161);
    result_by_q.set(c_Key1(19), common_term_161);
    // q = 20
    result_by_lpq.set(c_Key3(7, 3, 20), common_term_162);
    result_by_q.set(c_Key1(20), common_term_162);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 3), result_by_q);
    result_by_q.clear();

    // l , p = (7, 4).
    // q = -20
    result_by_lpq.set(c_Key3(7, 4, -20), common_term_162);
    result_by_q.set(c_Key1(-20), common_term_162);
    // q = -19
    result_by_lpq.set(c_Key3(7, 4, -19), common_term_161);
    result_by_q.set(c_Key1(-19), common_term_161);
    // q = -18
    result_by_lpq.set(c_Key3(7, 4, -18), common_term_160);
    result_by_q.set(c_Key1(-18), common_term_160);
    // q = -17
    result_by_lpq.set(c_Key3(7, 4, -17), common_term_159);
    result_by_q.set(c_Key1(-17), common_term_159);
    // q = -16
    result_by_lpq.set(c_Key3(7, 4, -16), common_term_158);
    result_by_q.set(c_Key1(-16), common_term_158);
    // q = -15
    result_by_lpq.set(c_Key3(7, 4, -15), common_term_157);
    result_by_q.set(c_Key1(-15), common_term_157);
    // q = -14
    result_by_lpq.set(c_Key3(7, 4, -14), common_term_156);
    result_by_q.set(c_Key1(-14), common_term_156);
    // q = -13
    result_by_lpq.set(c_Key3(7, 4, -13), common_term_155);
    result_by_q.set(c_Key1(-13), common_term_155);
    // q = -12
    result_by_lpq.set(c_Key3(7, 4, -12), common_term_154);
    result_by_q.set(c_Key1(-12), common_term_154);
    // q = -11
    result_by_lpq.set(c_Key3(7, 4, -11), common_term_153);
    result_by_q.set(c_Key1(-11), common_term_153);
    // q = -10
    result_by_lpq.set(c_Key3(7, 4, -10), common_term_152);
    result_by_q.set(c_Key1(-10), common_term_152);
    // q = -9
    result_by_lpq.set(c_Key3(7, 4, -9), common_term_151);
    result_by_q.set(c_Key1(-9), common_term_151);
    // q = -8
    result_by_lpq.set(c_Key3(7, 4, -8), common_term_150);
    result_by_q.set(c_Key1(-8), common_term_150);
    // q = -7
    result_by_lpq.set(c_Key3(7, 4, -7), common_term_149);
    result_by_q.set(c_Key1(-7), common_term_149);
    // q = -6
    result_by_lpq.set(c_Key3(7, 4, -6), common_term_148);
    result_by_q.set(c_Key1(-6), common_term_148);
    // q = -5
    result_by_lpq.set(c_Key3(7, 4, -5), common_term_147);
    result_by_q.set(c_Key1(-5), common_term_147);
    // q = -4
    result_by_lpq.set(c_Key3(7, 4, -4), common_term_146);
    result_by_q.set(c_Key1(-4), common_term_146);
    // q = -3
    result_by_lpq.set(c_Key3(7, 4, -3), common_term_145);
    result_by_q.set(c_Key1(-3), common_term_145);
    // q = -2
    result_by_lpq.set(c_Key3(7, 4, -2), common_term_144);
    result_by_q.set(c_Key1(-2), common_term_144);
    // q = -1
    result_by_lpq.set(c_Key3(7, 4, -1), common_term_143);
    result_by_q.set(c_Key1(-1), common_term_143);
    // q = 0
    result_by_lpq.set(c_Key3(7, 4, 0), common_term_142);
    result_by_q.set(c_Key1(0), common_term_142);
    // q = 1
    result_by_lpq.set(c_Key3(7, 4, 1), common_term_141);
    result_by_q.set(c_Key1(1), common_term_141);
    // q = 2
    result_by_lpq.set(c_Key3(7, 4, 2), common_term_140);
    result_by_q.set(c_Key1(2), common_term_140);
    // q = 3
    result_by_lpq.set(c_Key3(7, 4, 3), common_term_139);
    result_by_q.set(c_Key1(3), common_term_139);
    // q = 4
    result_by_lpq.set(c_Key3(7, 4, 4), common_term_138);
    result_by_q.set(c_Key1(4), common_term_138);
    // q = 5
    result_by_lpq.set(c_Key3(7, 4, 5), common_term_137);
    result_by_q.set(c_Key1(5), common_term_137);
    // q = 6
    result_by_lpq.set(c_Key3(7, 4, 6), common_term_136);
    result_by_q.set(c_Key1(6), common_term_136);
    // q = 7
    result_by_lpq.set(c_Key3(7, 4, 7), common_term_135);
    result_by_q.set(c_Key1(7), common_term_135);
    // q = 8
    result_by_lpq.set(c_Key3(7, 4, 8), common_term_134);
    result_by_q.set(c_Key1(8), common_term_134);
    // q = 9
    result_by_lpq.set(c_Key3(7, 4, 9), common_term_133);
    result_by_q.set(c_Key1(9), common_term_133);
    // q = 10
    result_by_lpq.set(c_Key3(7, 4, 10), common_term_132);
    result_by_q.set(c_Key1(10), common_term_132);
    // q = 11
    result_by_lpq.set(c_Key3(7, 4, 11), common_term_131);
    result_by_q.set(c_Key1(11), common_term_131);
    // q = 12
    result_by_lpq.set(c_Key3(7, 4, 12), common_term_130);
    result_by_q.set(c_Key1(12), common_term_130);
    // q = 13
    result_by_lpq.set(c_Key3(7, 4, 13), common_term_129);
    result_by_q.set(c_Key1(13), common_term_129);
    // q = 14
    result_by_lpq.set(c_Key3(7, 4, 14), common_term_128);
    result_by_q.set(c_Key1(14), common_term_128);
    // q = 15
    result_by_lpq.set(c_Key3(7, 4, 15), common_term_127);
    result_by_q.set(c_Key1(15), common_term_127);
    // q = 16
    result_by_lpq.set(c_Key3(7, 4, 16), common_term_126);
    result_by_q.set(c_Key1(16), common_term_126);
    // q = 17
    result_by_lpq.set(c_Key3(7, 4, 17), common_term_125);
    result_by_q.set(c_Key1(17), common_term_125);
    // q = 18
    result_by_lpq.set(c_Key3(7, 4, 18), common_term_124);
    result_by_q.set(c_Key1(18), common_term_124);
    // q = 19
    result_by_lpq.set(c_Key3(7, 4, 19), common_term_123);
    result_by_q.set(c_Key1(19), common_term_123);
    // q = 20
    result_by_lpq.set(c_Key3(7, 4, 20), common_term_122);
    result_by_q.set(c_Key1(20), common_term_122);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 4), result_by_q);
    result_by_q.clear();

    // l , p = (7, 5).
    // q = -20
    result_by_lpq.set(c_Key3(7, 5, -20), common_term_121);
    result_by_q.set(c_Key1(-20), common_term_121);
    // q = -19
    result_by_lpq.set(c_Key3(7, 5, -19), common_term_120);
    result_by_q.set(c_Key1(-19), common_term_120);
    // q = -18
    result_by_lpq.set(c_Key3(7, 5, -18), common_term_119);
    result_by_q.set(c_Key1(-18), common_term_119);
    // q = -17
    result_by_lpq.set(c_Key3(7, 5, -17), common_term_118);
    result_by_q.set(c_Key1(-17), common_term_118);
    // q = -16
    result_by_lpq.set(c_Key3(7, 5, -16), common_term_117);
    result_by_q.set(c_Key1(-16), common_term_117);
    // q = -15
    result_by_lpq.set(c_Key3(7, 5, -15), common_term_116);
    result_by_q.set(c_Key1(-15), common_term_116);
    // q = -14
    result_by_lpq.set(c_Key3(7, 5, -14), common_term_115);
    result_by_q.set(c_Key1(-14), common_term_115);
    // q = -13
    result_by_lpq.set(c_Key3(7, 5, -13), common_term_114);
    result_by_q.set(c_Key1(-13), common_term_114);
    // q = -12
    result_by_lpq.set(c_Key3(7, 5, -12), common_term_113);
    result_by_q.set(c_Key1(-12), common_term_113);
    // q = -11
    result_by_lpq.set(c_Key3(7, 5, -11), common_term_112);
    result_by_q.set(c_Key1(-11), common_term_112);
    // q = -10
    result_by_lpq.set(c_Key3(7, 5, -10), common_term_111);
    result_by_q.set(c_Key1(-10), common_term_111);
    // q = -9
    result_by_lpq.set(c_Key3(7, 5, -9), common_term_110);
    result_by_q.set(c_Key1(-9), common_term_110);
    // q = -8
    result_by_lpq.set(c_Key3(7, 5, -8), common_term_109);
    result_by_q.set(c_Key1(-8), common_term_109);
    // q = -7
    result_by_lpq.set(c_Key3(7, 5, -7), common_term_108);
    result_by_q.set(c_Key1(-7), common_term_108);
    // q = -6
    result_by_lpq.set(c_Key3(7, 5, -6), common_term_107);
    result_by_q.set(c_Key1(-6), common_term_107);
    // q = -5
    result_by_lpq.set(c_Key3(7, 5, -5), common_term_106);
    result_by_q.set(c_Key1(-5), common_term_106);
    // q = -4
    result_by_lpq.set(c_Key3(7, 5, -4), common_term_105);
    result_by_q.set(c_Key1(-4), common_term_105);
    // q = -3
    result_by_lpq.set(c_Key3(7, 5, -3), common_term_104);
    result_by_q.set(c_Key1(-3), common_term_104);
    // q = -2
    result_by_lpq.set(c_Key3(7, 5, -2), common_term_103);
    result_by_q.set(c_Key1(-2), common_term_103);
    // q = -1
    result_by_lpq.set(c_Key3(7, 5, -1), common_term_102);
    result_by_q.set(c_Key1(-1), common_term_102);
    // q = 0
    result_by_lpq.set(c_Key3(7, 5, 0), common_term_101);
    result_by_q.set(c_Key1(0), common_term_101);
    // q = 1
    result_by_lpq.set(c_Key3(7, 5, 1), common_term_100);
    result_by_q.set(c_Key1(1), common_term_100);
    // q = 2
    result_by_lpq.set(c_Key3(7, 5, 2), common_term_99);
    result_by_q.set(c_Key1(2), common_term_99);
    // q = 3
    result_by_lpq.set(c_Key3(7, 5, 3), common_term_98);
    result_by_q.set(c_Key1(3), common_term_98);
    // q = 4
    result_by_lpq.set(c_Key3(7, 5, 4), common_term_97);
    result_by_q.set(c_Key1(4), common_term_97);
    // q = 5
    result_by_lpq.set(c_Key3(7, 5, 5), common_term_96);
    result_by_q.set(c_Key1(5), common_term_96);
    // q = 6
    result_by_lpq.set(c_Key3(7, 5, 6), common_term_95);
    result_by_q.set(c_Key1(6), common_term_95);
    // q = 7
    result_by_lpq.set(c_Key3(7, 5, 7), common_term_94);
    result_by_q.set(c_Key1(7), common_term_94);
    // q = 8
    result_by_lpq.set(c_Key3(7, 5, 8), common_term_93);
    result_by_q.set(c_Key1(8), common_term_93);
    // q = 9
    result_by_lpq.set(c_Key3(7, 5, 9), common_term_92);
    result_by_q.set(c_Key1(9), common_term_92);
    // q = 10
    result_by_lpq.set(c_Key3(7, 5, 10), common_term_91);
    result_by_q.set(c_Key1(10), common_term_91);
    // q = 11
    result_by_lpq.set(c_Key3(7, 5, 11), common_term_90);
    result_by_q.set(c_Key1(11), common_term_90);
    // q = 12
    result_by_lpq.set(c_Key3(7, 5, 12), common_term_89);
    result_by_q.set(c_Key1(12), common_term_89);
    // q = 13
    result_by_lpq.set(c_Key3(7, 5, 13), common_term_88);
    result_by_q.set(c_Key1(13), common_term_88);
    // q = 14
    result_by_lpq.set(c_Key3(7, 5, 14), common_term_87);
    result_by_q.set(c_Key1(14), common_term_87);
    // q = 15
    result_by_lpq.set(c_Key3(7, 5, 15), common_term_86);
    result_by_q.set(c_Key1(15), common_term_86);
    // q = 16
    result_by_lpq.set(c_Key3(7, 5, 16), common_term_85);
    result_by_q.set(c_Key1(16), common_term_85);
    // q = 17
    result_by_lpq.set(c_Key3(7, 5, 17), common_term_84);
    result_by_q.set(c_Key1(17), common_term_84);
    // q = 18
    result_by_lpq.set(c_Key3(7, 5, 18), common_term_83);
    result_by_q.set(c_Key1(18), common_term_83);
    // q = 19
    result_by_lpq.set(c_Key3(7, 5, 19), common_term_82);
    result_by_q.set(c_Key1(19), common_term_82);
    // q = 20
    result_by_lpq.set(c_Key3(7, 5, 20), common_term_81);
    result_by_q.set(c_Key1(20), common_term_81);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 5), result_by_q);
    result_by_q.clear();

    // l , p = (7, 6).
    // q = -20
    result_by_lpq.set(c_Key3(7, 6, -20), common_term_80);
    result_by_q.set(c_Key1(-20), common_term_80);
    // q = -19
    result_by_lpq.set(c_Key3(7, 6, -19), common_term_79);
    result_by_q.set(c_Key1(-19), common_term_79);
    // q = -18
    result_by_lpq.set(c_Key3(7, 6, -18), common_term_78);
    result_by_q.set(c_Key1(-18), common_term_78);
    // q = -17
    result_by_lpq.set(c_Key3(7, 6, -17), common_term_77);
    result_by_q.set(c_Key1(-17), common_term_77);
    // q = -16
    result_by_lpq.set(c_Key3(7, 6, -16), common_term_76);
    result_by_q.set(c_Key1(-16), common_term_76);
    // q = -15
    result_by_lpq.set(c_Key3(7, 6, -15), common_term_75);
    result_by_q.set(c_Key1(-15), common_term_75);
    // q = -14
    result_by_lpq.set(c_Key3(7, 6, -14), common_term_74);
    result_by_q.set(c_Key1(-14), common_term_74);
    // q = -13
    result_by_lpq.set(c_Key3(7, 6, -13), common_term_73);
    result_by_q.set(c_Key1(-13), common_term_73);
    // q = -12
    result_by_lpq.set(c_Key3(7, 6, -12), common_term_72);
    result_by_q.set(c_Key1(-12), common_term_72);
    // q = -11
    result_by_lpq.set(c_Key3(7, 6, -11), common_term_71);
    result_by_q.set(c_Key1(-11), common_term_71);
    // q = -10
    result_by_lpq.set(c_Key3(7, 6, -10), common_term_70);
    result_by_q.set(c_Key1(-10), common_term_70);
    // q = -9
    result_by_lpq.set(c_Key3(7, 6, -9), common_term_69);
    result_by_q.set(c_Key1(-9), common_term_69);
    // q = -8
    result_by_lpq.set(c_Key3(7, 6, -8), common_term_68);
    result_by_q.set(c_Key1(-8), common_term_68);
    // q = -7
    result_by_lpq.set(c_Key3(7, 6, -7), common_term_67);
    result_by_q.set(c_Key1(-7), common_term_67);
    // q = -6
    result_by_lpq.set(c_Key3(7, 6, -6), common_term_66);
    result_by_q.set(c_Key1(-6), common_term_66);
    // q = -5
    result_by_lpq.set(c_Key3(7, 6, -5), common_term_65);
    result_by_q.set(c_Key1(-5), common_term_65);
    // q = -4
    result_by_lpq.set(c_Key3(7, 6, -4), common_term_64);
    result_by_q.set(c_Key1(-4), common_term_64);
    // q = -3
    result_by_lpq.set(c_Key3(7, 6, -3), common_term_63);
    result_by_q.set(c_Key1(-3), common_term_63);
    // q = -2
    result_by_lpq.set(c_Key3(7, 6, -2), common_term_62);
    result_by_q.set(c_Key1(-2), common_term_62);
    // q = -1
    result_by_lpq.set(c_Key3(7, 6, -1), common_term_61);
    result_by_q.set(c_Key1(-1), common_term_61);
    // q = 0
    result_by_lpq.set(c_Key3(7, 6, 0), common_term_60);
    result_by_q.set(c_Key1(0), common_term_60);
    // q = 1
    result_by_lpq.set(c_Key3(7, 6, 1), common_term_59);
    result_by_q.set(c_Key1(1), common_term_59);
    // q = 2
    result_by_lpq.set(c_Key3(7, 6, 2), common_term_58);
    result_by_q.set(c_Key1(2), common_term_58);
    // q = 3
    result_by_lpq.set(c_Key3(7, 6, 3), common_term_57);
    result_by_q.set(c_Key1(3), common_term_57);
    // q = 4
    result_by_lpq.set(c_Key3(7, 6, 4), common_term_56);
    result_by_q.set(c_Key1(4), common_term_56);
    // q = 5
    result_by_lpq.set(c_Key3(7, 6, 5), common_term_55);
    result_by_q.set(c_Key1(5), common_term_55);
    // q = 6
    result_by_lpq.set(c_Key3(7, 6, 6), common_term_54);
    result_by_q.set(c_Key1(6), common_term_54);
    // q = 7
    result_by_lpq.set(c_Key3(7, 6, 7), common_term_53);
    result_by_q.set(c_Key1(7), common_term_53);
    // q = 8
    result_by_lpq.set(c_Key3(7, 6, 8), common_term_52);
    result_by_q.set(c_Key1(8), common_term_52);
    // q = 9
    result_by_lpq.set(c_Key3(7, 6, 9), common_term_51);
    result_by_q.set(c_Key1(9), common_term_51);
    // q = 10
    result_by_lpq.set(c_Key3(7, 6, 10), common_term_50);
    result_by_q.set(c_Key1(10), common_term_50);
    // q = 11
    result_by_lpq.set(c_Key3(7, 6, 11), common_term_49);
    result_by_q.set(c_Key1(11), common_term_49);
    // q = 12
    result_by_lpq.set(c_Key3(7, 6, 12), common_term_48);
    result_by_q.set(c_Key1(12), common_term_48);
    // q = 13
    result_by_lpq.set(c_Key3(7, 6, 13), common_term_47);
    result_by_q.set(c_Key1(13), common_term_47);
    // q = 14
    result_by_lpq.set(c_Key3(7, 6, 14), common_term_46);
    result_by_q.set(c_Key1(14), common_term_46);
    // q = 15
    result_by_lpq.set(c_Key3(7, 6, 15), common_term_45);
    result_by_q.set(c_Key1(15), common_term_45);
    // q = 16
    result_by_lpq.set(c_Key3(7, 6, 16), common_term_44);
    result_by_q.set(c_Key1(16), common_term_44);
    // q = 17
    result_by_lpq.set(c_Key3(7, 6, 17), common_term_43);
    result_by_q.set(c_Key1(17), common_term_43);
    // q = 18
    result_by_lpq.set(c_Key3(7, 6, 18), common_term_42);
    result_by_q.set(c_Key1(18), common_term_42);
    // q = 19
    result_by_lpq.set(c_Key3(7, 6, 19), common_term_41);
    result_by_q.set(c_Key1(19), common_term_41);
    // q = 20
    result_by_lpq.set(c_Key3(7, 6, 20), common_term_40);
    result_by_q.set(c_Key1(20), common_term_40);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 6), result_by_q);
    result_by_q.clear();

    // l , p = (7, 7).
    // q = -20
    result_by_lpq.set(c_Key3(7, 7, -20), common_term_39);
    result_by_q.set(c_Key1(-20), common_term_39);
    // q = -19
    result_by_lpq.set(c_Key3(7, 7, -19), common_term_38);
    result_by_q.set(c_Key1(-19), common_term_38);
    // q = -18
    result_by_lpq.set(c_Key3(7, 7, -18), common_term_37);
    result_by_q.set(c_Key1(-18), common_term_37);
    // q = -17
    result_by_lpq.set(c_Key3(7, 7, -17), common_term_36);
    result_by_q.set(c_Key1(-17), common_term_36);
    // q = -16
    result_by_lpq.set(c_Key3(7, 7, -16), common_term_35);
    result_by_q.set(c_Key1(-16), common_term_35);
    // q = -15
    result_by_lpq.set(c_Key3(7, 7, -15), common_term_34);
    result_by_q.set(c_Key1(-15), common_term_34);
    // q = -14
    result_by_lpq.set(c_Key3(7, 7, -14), common_term_33);
    result_by_q.set(c_Key1(-14), common_term_33);
    // q = -13
    result_by_lpq.set(c_Key3(7, 7, -13), common_term_32);
    result_by_q.set(c_Key1(-13), common_term_32);
    // q = -12
    result_by_lpq.set(c_Key3(7, 7, -12), common_term_31);
    result_by_q.set(c_Key1(-12), common_term_31);
    // q = -11
    result_by_lpq.set(c_Key3(7, 7, -11), common_term_30);
    result_by_q.set(c_Key1(-11), common_term_30);
    // q = -10
    result_by_lpq.set(c_Key3(7, 7, -10), common_term_29);
    result_by_q.set(c_Key1(-10), common_term_29);
    // q = -9
    result_by_lpq.set(c_Key3(7, 7, -9), common_term_28);
    result_by_q.set(c_Key1(-9), common_term_28);
    // q = -8
    result_by_lpq.set(c_Key3(7, 7, -8), common_term_27);
    result_by_q.set(c_Key1(-8), common_term_27);
    // q = -7
    result_by_lpq.set(c_Key3(7, 7, -7), common_term_26);
    result_by_q.set(c_Key1(-7), common_term_26);
    // q = -6
    result_by_lpq.set(c_Key3(7, 7, -6), common_term_25);
    result_by_q.set(c_Key1(-6), common_term_25);
    // q = -5
    result_by_lpq.set(c_Key3(7, 7, -5), common_term_24);
    result_by_q.set(c_Key1(-5), common_term_24);
    // q = -4
    result_by_lpq.set(c_Key3(7, 7, -4), common_term_23);
    result_by_q.set(c_Key1(-4), common_term_23);
    // q = -3
    result_by_lpq.set(c_Key3(7, 7, -3), common_term_22);
    result_by_q.set(c_Key1(-3), common_term_22);
    // q = -2
    result_by_lpq.set(c_Key3(7, 7, -2), common_term_21);
    result_by_q.set(c_Key1(-2), common_term_21);
    // q = -1
    result_by_lpq.set(c_Key3(7, 7, -1), common_term_20);
    result_by_q.set(c_Key1(-1), common_term_20);
    // q = 0
    result_by_lpq.set(c_Key3(7, 7, 0), common_term_19);
    result_by_q.set(c_Key1(0), common_term_19);
    // q = 1
    result_by_lpq.set(c_Key3(7, 7, 1), common_term_18);
    result_by_q.set(c_Key1(1), common_term_18);
    // q = 2
    result_by_lpq.set(c_Key3(7, 7, 2), common_term_17);
    result_by_q.set(c_Key1(2), common_term_17);
    // q = 3
    result_by_lpq.set(c_Key3(7, 7, 3), common_term_16);
    result_by_q.set(c_Key1(3), common_term_16);
    // q = 4
    result_by_lpq.set(c_Key3(7, 7, 4), common_term_15);
    result_by_q.set(c_Key1(4), common_term_15);
    // q = 5
    result_by_lpq.set(c_Key3(7, 7, 5), common_term_14);
    result_by_q.set(c_Key1(5), common_term_14);
    // q = 6
    result_by_lpq.set(c_Key3(7, 7, 6), common_term_13);
    result_by_q.set(c_Key1(6), common_term_13);
    // q = 8
    result_by_lpq.set(c_Key3(7, 7, 8), common_term_12);
    result_by_q.set(c_Key1(8), common_term_12);
    // q = 9
    result_by_lpq.set(c_Key3(7, 7, 9), common_term_11);
    result_by_q.set(c_Key1(9), common_term_11);
    // q = 10
    result_by_lpq.set(c_Key3(7, 7, 10), common_term_10);
    result_by_q.set(c_Key1(10), common_term_10);
    // q = 11
    result_by_lpq.set(c_Key3(7, 7, 11), common_term_9);
    result_by_q.set(c_Key1(11), common_term_9);
    // q = 12
    result_by_lpq.set(c_Key3(7, 7, 12), common_term_8);
    result_by_q.set(c_Key1(12), common_term_8);
    // q = 13
    result_by_lpq.set(c_Key3(7, 7, 13), common_term_7);
    result_by_q.set(c_Key1(13), common_term_7);
    // q = 14
    result_by_lpq.set(c_Key3(7, 7, 14), common_term_6);
    result_by_q.set(c_Key1(14), common_term_6);
    // q = 15
    result_by_lpq.set(c_Key3(7, 7, 15), common_term_5);
    result_by_q.set(c_Key1(15), common_term_5);
    // q = 16
    result_by_lpq.set(c_Key3(7, 7, 16), common_term_4);
    result_by_q.set(c_Key1(16), common_term_4);
    // q = 17
    result_by_lpq.set(c_Key3(7, 7, 17), common_term_3);
    result_by_q.set(c_Key1(17), common_term_3);
    // q = 18
    result_by_lpq.set(c_Key3(7, 7, 18), common_term_2);
    result_by_q.set(c_Key1(18), common_term_2);
    // q = 19
    result_by_lpq.set(c_Key3(7, 7, 19), common_term_1);
    result_by_q.set(c_Key1(19), common_term_1);
    // q = 20
    result_by_lpq.set(c_Key3(7, 7, 20), common_term_0);
    result_by_q.set(c_Key1(20), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(7, 7), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}
