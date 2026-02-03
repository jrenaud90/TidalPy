#pragma once

#include "eccentricity_common_.hpp"

EccentricityFuncOutput c_eccentricity_function_l9_e1(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 9.
    // Functions are taylor expanded and truncated at eccentricity^1.

    // Functions Calculated for l = 9.

    c_IntMap<c_Key3, double> result_by_lpq(12);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(10);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = 4.0*eccentricity*std::pow(1.0 - eccentricity_2, -8.5);

    c_IntMap<c_Key1, double> result_by_q(2);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (9, 0).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 0, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 0), result_by_q);
    result_by_q.clear();

    // l , p = (9, 1).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 1), result_by_q);
    result_by_q.clear();

    // l , p = (9, 2).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 2), result_by_q);
    result_by_q.clear();

    // l , p = (9, 3).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 3), result_by_q);
    result_by_q.clear();

    // l , p = (9, 4).
    // q = -1
    result_by_lpq.set(c_Key3(9, 4, -1), common_term_0);
    result_by_q.set(c_Key1(-1), common_term_0);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 4, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 4), result_by_q);
    result_by_q.clear();

    // l , p = (9, 5).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 5, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(9, 5, 1), common_term_0);
    result_by_q.set(c_Key1(1), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 5), result_by_q);
    result_by_q.clear();

    // l , p = (9, 6).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 6, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 6), result_by_q);
    result_by_q.clear();

    // l , p = (9, 7).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 7, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 7), result_by_q);
    result_by_q.clear();

    // l , p = (9, 8).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 8, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 8), result_by_q);
    result_by_q.clear();

    // l , p = (9, 9).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 9, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 9), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l9_e2(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 9.
    // Functions are taylor expanded and truncated at eccentricity^2.

    // Functions Calculated for l = 9.

    c_IntMap<c_Key3, double> result_by_lpq(28);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(10);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = -4.0*eccentricity;
    double common_term_1 = 14.0*eccentricity;
    double common_term_2 = -2.0*eccentricity;
    double common_term_3 = 12.0*eccentricity;
    double common_term_4 = 10.0*eccentricity;
    double common_term_5 = 2.0*eccentricity;
    double common_term_6 = 8.0*eccentricity;
    double common_term_7 = 4.0*eccentricity*std::pow(1.0 - eccentricity_2, -8.5);
    double common_term_8 = 6.0*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(3);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (9, 0).
    // q = -1
    result_by_lpq.set(c_Key3(9, 0, -1), common_term_0);
    result_by_q.set(c_Key1(-1), common_term_0);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 0, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(9, 0, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 0), result_by_q);
    result_by_q.clear();

    // l , p = (9, 1).
    // q = -1
    result_by_lpq.set(c_Key3(9, 1, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(9, 1, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 1), result_by_q);
    result_by_q.clear();

    // l , p = (9, 2).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(9, 2, 1), common_term_4);
    result_by_q.set(c_Key1(1), common_term_4);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 2), result_by_q);
    result_by_q.clear();

    // l , p = (9, 3).
    // q = -1
    result_by_lpq.set(c_Key3(9, 3, -1), common_term_5);
    result_by_q.set(c_Key1(-1), common_term_5);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(9, 3, 1), common_term_6);
    result_by_q.set(c_Key1(1), common_term_6);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 3), result_by_q);
    result_by_q.clear();

    // l , p = (9, 4).
    // q = -1
    result_by_lpq.set(c_Key3(9, 4, -1), common_term_7);
    result_by_q.set(c_Key1(-1), common_term_7);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 4, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(9, 4, 1), common_term_8);
    result_by_q.set(c_Key1(1), common_term_8);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 4), result_by_q);
    result_by_q.clear();

    // l , p = (9, 5).
    // q = -1
    result_by_lpq.set(c_Key3(9, 5, -1), common_term_8);
    result_by_q.set(c_Key1(-1), common_term_8);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 5, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(9, 5, 1), common_term_7);
    result_by_q.set(c_Key1(1), common_term_7);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 5), result_by_q);
    result_by_q.clear();

    // l , p = (9, 6).
    // q = -1
    result_by_lpq.set(c_Key3(9, 6, -1), common_term_6);
    result_by_q.set(c_Key1(-1), common_term_6);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 6, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(9, 6, 1), common_term_5);
    result_by_q.set(c_Key1(1), common_term_5);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 6), result_by_q);
    result_by_q.clear();

    // l , p = (9, 7).
    // q = -1
    result_by_lpq.set(c_Key3(9, 7, -1), common_term_4);
    result_by_q.set(c_Key1(-1), common_term_4);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 7, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 7), result_by_q);
    result_by_q.clear();

    // l , p = (9, 8).
    // q = -1
    result_by_lpq.set(c_Key3(9, 8, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 8, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(9, 8, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 8), result_by_q);
    result_by_q.clear();

    // l , p = (9, 9).
    // q = -1
    result_by_lpq.set(c_Key3(9, 9, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(9, 9, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(9, 9, 1), common_term_0);
    result_by_q.set(c_Key1(1), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 9), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l9_e3(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 9.
    // Functions are taylor expanded and truncated at eccentricity^3.

    // Functions Calculated for l = 9.

    c_IntMap<c_Key3, double> result_by_lpq(50);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(10);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double common_term_0 = 6.125*eccentricity_2;
    double common_term_1 = -4.0*eccentricity;
    double common_term_2 = 1.0 - 58.5*eccentricity_2;
    double common_term_3 = 14.0*eccentricity;
    double common_term_4 = 107.375*eccentricity_2;
    double common_term_5 = 1.375*eccentricity_2;
    double common_term_6 = -2.0*eccentricity;
    double common_term_7 = 1.0 - 26.5*eccentricity_2;
    double common_term_8 = 12.0*eccentricity;
    double common_term_9 = 80.125*eccentricity_2;
    double common_term_10 = 0.625*eccentricity_2;
    double common_term_11 = 1.0 - 2.5*eccentricity_2;
    double common_term_12 = 10.0*eccentricity;
    double common_term_13 = 56.875*eccentricity_2;
    double common_term_14 = 7.0*eccentricity_3*std::pow(1.0 - eccentricity_2, -8.5);
    double common_term_15 = 3.875*eccentricity_2;
    double common_term_16 = 2.0*eccentricity;
    double common_term_17 = 13.5*eccentricity_2 + 1.0;
    double common_term_18 = 8.0*eccentricity;
    double common_term_19 = 37.625*eccentricity_2;
    double common_term_20 = 11.125*eccentricity_2;
    double common_term_21 = std::pow(1.0 - eccentricity_2, -8.5)*(21.0*eccentricity_3 + 4.0*eccentricity);
    double common_term_22 = 21.5*eccentricity_2 + 1.0;
    double common_term_23 = 6.0*eccentricity;
    double common_term_24 = 22.375*eccentricity_2;

    c_IntMap<c_Key1, double> result_by_q(6);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (9, 0).
    // q = -2
    result_by_lpq.set(c_Key3(9, 0, -2), common_term_0);
    result_by_q.set(c_Key1(-2), common_term_0);
    // q = -1
    result_by_lpq.set(c_Key3(9, 0, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    result_by_lpq.set(c_Key3(9, 0, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(9, 0, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // q = 2
    result_by_lpq.set(c_Key3(9, 0, 2), common_term_4);
    result_by_q.set(c_Key1(2), common_term_4);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 0), result_by_q);
    result_by_q.clear();

    // l , p = (9, 1).
    // q = -2
    result_by_lpq.set(c_Key3(9, 1, -2), common_term_5);
    result_by_q.set(c_Key1(-2), common_term_5);
    // q = -1
    result_by_lpq.set(c_Key3(9, 1, -1), common_term_6);
    result_by_q.set(c_Key1(-1), common_term_6);
    // q = 0
    result_by_lpq.set(c_Key3(9, 1, 0), common_term_7);
    result_by_q.set(c_Key1(0), common_term_7);
    // q = 1
    result_by_lpq.set(c_Key3(9, 1, 1), common_term_8);
    result_by_q.set(c_Key1(1), common_term_8);
    // q = 2
    result_by_lpq.set(c_Key3(9, 1, 2), common_term_9);
    result_by_q.set(c_Key1(2), common_term_9);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 1), result_by_q);
    result_by_q.clear();

    // l , p = (9, 2).
    // q = -2
    result_by_lpq.set(c_Key3(9, 2, -2), common_term_10);
    result_by_q.set(c_Key1(-2), common_term_10);
    // q = 0
    result_by_lpq.set(c_Key3(9, 2, 0), common_term_11);
    result_by_q.set(c_Key1(0), common_term_11);
    // q = 1
    result_by_lpq.set(c_Key3(9, 2, 1), common_term_12);
    result_by_q.set(c_Key1(1), common_term_12);
    // q = 2
    result_by_lpq.set(c_Key3(9, 2, 2), common_term_13);
    result_by_q.set(c_Key1(2), common_term_13);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 2), result_by_q);
    result_by_q.clear();

    // l , p = (9, 3).
    // q = -3
    result_by_lpq.set(c_Key3(9, 3, -3), common_term_14);
    result_by_q.set(c_Key1(-3), common_term_14);
    // q = -2
    result_by_lpq.set(c_Key3(9, 3, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(9, 3, -1), common_term_16);
    result_by_q.set(c_Key1(-1), common_term_16);
    // q = 0
    result_by_lpq.set(c_Key3(9, 3, 0), common_term_17);
    result_by_q.set(c_Key1(0), common_term_17);
    // q = 1
    result_by_lpq.set(c_Key3(9, 3, 1), common_term_18);
    result_by_q.set(c_Key1(1), common_term_18);
    // q = 2
    result_by_lpq.set(c_Key3(9, 3, 2), common_term_19);
    result_by_q.set(c_Key1(2), common_term_19);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 3), result_by_q);
    result_by_q.clear();

    // l , p = (9, 4).
    // q = -2
    result_by_lpq.set(c_Key3(9, 4, -2), common_term_20);
    result_by_q.set(c_Key1(-2), common_term_20);
    // q = -1
    result_by_lpq.set(c_Key3(9, 4, -1), common_term_21);
    result_by_q.set(c_Key1(-1), common_term_21);
    // q = 0
    result_by_lpq.set(c_Key3(9, 4, 0), common_term_22);
    result_by_q.set(c_Key1(0), common_term_22);
    // q = 1
    result_by_lpq.set(c_Key3(9, 4, 1), common_term_23);
    result_by_q.set(c_Key1(1), common_term_23);
    // q = 2
    result_by_lpq.set(c_Key3(9, 4, 2), common_term_24);
    result_by_q.set(c_Key1(2), common_term_24);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 4), result_by_q);
    result_by_q.clear();

    // l , p = (9, 5).
    // q = -2
    result_by_lpq.set(c_Key3(9, 5, -2), common_term_24);
    result_by_q.set(c_Key1(-2), common_term_24);
    // q = -1
    result_by_lpq.set(c_Key3(9, 5, -1), common_term_23);
    result_by_q.set(c_Key1(-1), common_term_23);
    // q = 0
    result_by_lpq.set(c_Key3(9, 5, 0), common_term_22);
    result_by_q.set(c_Key1(0), common_term_22);
    // q = 1
    result_by_lpq.set(c_Key3(9, 5, 1), common_term_21);
    result_by_q.set(c_Key1(1), common_term_21);
    // q = 2
    result_by_lpq.set(c_Key3(9, 5, 2), common_term_20);
    result_by_q.set(c_Key1(2), common_term_20);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 5), result_by_q);
    result_by_q.clear();

    // l , p = (9, 6).
    // q = -2
    result_by_lpq.set(c_Key3(9, 6, -2), common_term_19);
    result_by_q.set(c_Key1(-2), common_term_19);
    // q = -1
    result_by_lpq.set(c_Key3(9, 6, -1), common_term_18);
    result_by_q.set(c_Key1(-1), common_term_18);
    // q = 0
    result_by_lpq.set(c_Key3(9, 6, 0), common_term_17);
    result_by_q.set(c_Key1(0), common_term_17);
    // q = 1
    result_by_lpq.set(c_Key3(9, 6, 1), common_term_16);
    result_by_q.set(c_Key1(1), common_term_16);
    // q = 2
    result_by_lpq.set(c_Key3(9, 6, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // q = 3
    result_by_lpq.set(c_Key3(9, 6, 3), common_term_14);
    result_by_q.set(c_Key1(3), common_term_14);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 6), result_by_q);
    result_by_q.clear();

    // l , p = (9, 7).
    // q = -2
    result_by_lpq.set(c_Key3(9, 7, -2), common_term_13);
    result_by_q.set(c_Key1(-2), common_term_13);
    // q = -1
    result_by_lpq.set(c_Key3(9, 7, -1), common_term_12);
    result_by_q.set(c_Key1(-1), common_term_12);
    // q = 0
    result_by_lpq.set(c_Key3(9, 7, 0), common_term_11);
    result_by_q.set(c_Key1(0), common_term_11);
    // q = 2
    result_by_lpq.set(c_Key3(9, 7, 2), common_term_10);
    result_by_q.set(c_Key1(2), common_term_10);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 7), result_by_q);
    result_by_q.clear();

    // l , p = (9, 8).
    // q = -2
    result_by_lpq.set(c_Key3(9, 8, -2), common_term_9);
    result_by_q.set(c_Key1(-2), common_term_9);
    // q = -1
    result_by_lpq.set(c_Key3(9, 8, -1), common_term_8);
    result_by_q.set(c_Key1(-1), common_term_8);
    // q = 0
    result_by_lpq.set(c_Key3(9, 8, 0), common_term_7);
    result_by_q.set(c_Key1(0), common_term_7);
    // q = 1
    result_by_lpq.set(c_Key3(9, 8, 1), common_term_6);
    result_by_q.set(c_Key1(1), common_term_6);
    // q = 2
    result_by_lpq.set(c_Key3(9, 8, 2), common_term_5);
    result_by_q.set(c_Key1(2), common_term_5);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 8), result_by_q);
    result_by_q.clear();

    // l , p = (9, 9).
    // q = -2
    result_by_lpq.set(c_Key3(9, 9, -2), common_term_4);
    result_by_q.set(c_Key1(-2), common_term_4);
    // q = -1
    result_by_lpq.set(c_Key3(9, 9, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    result_by_lpq.set(c_Key3(9, 9, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(9, 9, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // q = 2
    result_by_lpq.set(c_Key3(9, 9, 2), common_term_0);
    result_by_q.set(c_Key1(2), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 9), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l9_e4(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 9.
    // Functions are taylor expanded and truncated at eccentricity^4.

    // Functions Calculated for l = 9.

    c_IntMap<c_Key3, double> result_by_lpq(70);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(10);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double common_term_0 = -4.5*eccentricity_3;
    double common_term_1 = 6.125*eccentricity_2;
    double common_term_2 = 95.0*eccentricity_3 - 4.0*eccentricity;
    double common_term_3 = 1.0 - 58.5*eccentricity_2;
    double common_term_4 = -467.5*eccentricity_3 + 14.0*eccentricity;
    double common_term_5 = 107.375*eccentricity_2;
    double common_term_6 = 597.0*eccentricity_3;
    double common_term_7 = -0.33333333333333333*eccentricity_3;
    double common_term_8 = 1.375*eccentricity_2;
    double common_term_9 = 20.5*eccentricity_3 - 2.0*eccentricity;
    double common_term_10 = 1.0 - 26.5*eccentricity_2;
    double common_term_11 = -193.0*eccentricity_3 + 12.0*eccentricity;
    double common_term_12 = 80.125*eccentricity_2;
    double common_term_13 = 392.83333333333333*eccentricity_3;
    double common_term_14 = 0.83333333333333333*eccentricity_3;
    double common_term_15 = 0.625*eccentricity_2;
    double common_term_16 = 5.0*eccentricity_3;
    double common_term_17 = 1.0 - 2.5*eccentricity_2;
    double common_term_18 = -27.5*eccentricity_3 + 10.0*eccentricity;
    double common_term_19 = 56.875*eccentricity_2;
    double common_term_20 = 241.66666666666667*eccentricity_3;
    double common_term_21 = 7.0*eccentricity_3*std::pow(1.0 - eccentricity_2, -8.5);
    double common_term_22 = 3.875*eccentricity_2;
    double common_term_23 = 24.5*eccentricity_3 + 2.0*eccentricity;
    double common_term_24 = 13.5*eccentricity_2 + 1.0;
    double common_term_25 = 53.0*eccentricity_3 + 8.0*eccentricity;
    double common_term_26 = 37.625*eccentricity_2;
    double common_term_27 = 135.5*eccentricity_3;
    double common_term_28 = 26.166666666666667*eccentricity_3;
    double common_term_29 = 11.125*eccentricity_2;
    double common_term_30 = std::pow(1.0 - eccentricity_2, -8.5)*(21.0*eccentricity_3 + 4.0*eccentricity);
    double common_term_31 = 21.5*eccentricity_2 + 1.0;
    double common_term_32 = 72.5*eccentricity_3 + 6.0*eccentricity;
    double common_term_33 = 22.375*eccentricity_2;
    double common_term_34 = 66.333333333333333*eccentricity_3;

    c_IntMap<c_Key1, double> result_by_q(7);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (9, 0).
    // q = -3
    result_by_lpq.set(c_Key3(9, 0, -3), common_term_0);
    result_by_q.set(c_Key1(-3), common_term_0);
    // q = -2
    result_by_lpq.set(c_Key3(9, 0, -2), common_term_1);
    result_by_q.set(c_Key1(-2), common_term_1);
    // q = -1
    result_by_lpq.set(c_Key3(9, 0, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    result_by_lpq.set(c_Key3(9, 0, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(9, 0, 1), common_term_4);
    result_by_q.set(c_Key1(1), common_term_4);
    // q = 2
    result_by_lpq.set(c_Key3(9, 0, 2), common_term_5);
    result_by_q.set(c_Key1(2), common_term_5);
    // q = 3
    result_by_lpq.set(c_Key3(9, 0, 3), common_term_6);
    result_by_q.set(c_Key1(3), common_term_6);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 0), result_by_q);
    result_by_q.clear();

    // l , p = (9, 1).
    // q = -3
    result_by_lpq.set(c_Key3(9, 1, -3), common_term_7);
    result_by_q.set(c_Key1(-3), common_term_7);
    // q = -2
    result_by_lpq.set(c_Key3(9, 1, -2), common_term_8);
    result_by_q.set(c_Key1(-2), common_term_8);
    // q = -1
    result_by_lpq.set(c_Key3(9, 1, -1), common_term_9);
    result_by_q.set(c_Key1(-1), common_term_9);
    // q = 0
    result_by_lpq.set(c_Key3(9, 1, 0), common_term_10);
    result_by_q.set(c_Key1(0), common_term_10);
    // q = 1
    result_by_lpq.set(c_Key3(9, 1, 1), common_term_11);
    result_by_q.set(c_Key1(1), common_term_11);
    // q = 2
    result_by_lpq.set(c_Key3(9, 1, 2), common_term_12);
    result_by_q.set(c_Key1(2), common_term_12);
    // q = 3
    result_by_lpq.set(c_Key3(9, 1, 3), common_term_13);
    result_by_q.set(c_Key1(3), common_term_13);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 1), result_by_q);
    result_by_q.clear();

    // l , p = (9, 2).
    // q = -3
    result_by_lpq.set(c_Key3(9, 2, -3), common_term_14);
    result_by_q.set(c_Key1(-3), common_term_14);
    // q = -2
    result_by_lpq.set(c_Key3(9, 2, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(9, 2, -1), common_term_16);
    result_by_q.set(c_Key1(-1), common_term_16);
    // q = 0
    result_by_lpq.set(c_Key3(9, 2, 0), common_term_17);
    result_by_q.set(c_Key1(0), common_term_17);
    // q = 1
    result_by_lpq.set(c_Key3(9, 2, 1), common_term_18);
    result_by_q.set(c_Key1(1), common_term_18);
    // q = 2
    result_by_lpq.set(c_Key3(9, 2, 2), common_term_19);
    result_by_q.set(c_Key1(2), common_term_19);
    // q = 3
    result_by_lpq.set(c_Key3(9, 2, 3), common_term_20);
    result_by_q.set(c_Key1(3), common_term_20);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 2), result_by_q);
    result_by_q.clear();

    // l , p = (9, 3).
    // q = -3
    result_by_lpq.set(c_Key3(9, 3, -3), common_term_21);
    result_by_q.set(c_Key1(-3), common_term_21);
    // q = -2
    result_by_lpq.set(c_Key3(9, 3, -2), common_term_22);
    result_by_q.set(c_Key1(-2), common_term_22);
    // q = -1
    result_by_lpq.set(c_Key3(9, 3, -1), common_term_23);
    result_by_q.set(c_Key1(-1), common_term_23);
    // q = 0
    result_by_lpq.set(c_Key3(9, 3, 0), common_term_24);
    result_by_q.set(c_Key1(0), common_term_24);
    // q = 1
    result_by_lpq.set(c_Key3(9, 3, 1), common_term_25);
    result_by_q.set(c_Key1(1), common_term_25);
    // q = 2
    result_by_lpq.set(c_Key3(9, 3, 2), common_term_26);
    result_by_q.set(c_Key1(2), common_term_26);
    // q = 3
    result_by_lpq.set(c_Key3(9, 3, 3), common_term_27);
    result_by_q.set(c_Key1(3), common_term_27);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 3), result_by_q);
    result_by_q.clear();

    // l , p = (9, 4).
    // q = -3
    result_by_lpq.set(c_Key3(9, 4, -3), common_term_28);
    result_by_q.set(c_Key1(-3), common_term_28);
    // q = -2
    result_by_lpq.set(c_Key3(9, 4, -2), common_term_29);
    result_by_q.set(c_Key1(-2), common_term_29);
    // q = -1
    result_by_lpq.set(c_Key3(9, 4, -1), common_term_30);
    result_by_q.set(c_Key1(-1), common_term_30);
    // q = 0
    result_by_lpq.set(c_Key3(9, 4, 0), common_term_31);
    result_by_q.set(c_Key1(0), common_term_31);
    // q = 1
    result_by_lpq.set(c_Key3(9, 4, 1), common_term_32);
    result_by_q.set(c_Key1(1), common_term_32);
    // q = 2
    result_by_lpq.set(c_Key3(9, 4, 2), common_term_33);
    result_by_q.set(c_Key1(2), common_term_33);
    // q = 3
    result_by_lpq.set(c_Key3(9, 4, 3), common_term_34);
    result_by_q.set(c_Key1(3), common_term_34);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 4), result_by_q);
    result_by_q.clear();

    // l , p = (9, 5).
    // q = -3
    result_by_lpq.set(c_Key3(9, 5, -3), common_term_34);
    result_by_q.set(c_Key1(-3), common_term_34);
    // q = -2
    result_by_lpq.set(c_Key3(9, 5, -2), common_term_33);
    result_by_q.set(c_Key1(-2), common_term_33);
    // q = -1
    result_by_lpq.set(c_Key3(9, 5, -1), common_term_32);
    result_by_q.set(c_Key1(-1), common_term_32);
    // q = 0
    result_by_lpq.set(c_Key3(9, 5, 0), common_term_31);
    result_by_q.set(c_Key1(0), common_term_31);
    // q = 1
    result_by_lpq.set(c_Key3(9, 5, 1), common_term_30);
    result_by_q.set(c_Key1(1), common_term_30);
    // q = 2
    result_by_lpq.set(c_Key3(9, 5, 2), common_term_29);
    result_by_q.set(c_Key1(2), common_term_29);
    // q = 3
    result_by_lpq.set(c_Key3(9, 5, 3), common_term_28);
    result_by_q.set(c_Key1(3), common_term_28);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 5), result_by_q);
    result_by_q.clear();

    // l , p = (9, 6).
    // q = -3
    result_by_lpq.set(c_Key3(9, 6, -3), common_term_27);
    result_by_q.set(c_Key1(-3), common_term_27);
    // q = -2
    result_by_lpq.set(c_Key3(9, 6, -2), common_term_26);
    result_by_q.set(c_Key1(-2), common_term_26);
    // q = -1
    result_by_lpq.set(c_Key3(9, 6, -1), common_term_25);
    result_by_q.set(c_Key1(-1), common_term_25);
    // q = 0
    result_by_lpq.set(c_Key3(9, 6, 0), common_term_24);
    result_by_q.set(c_Key1(0), common_term_24);
    // q = 1
    result_by_lpq.set(c_Key3(9, 6, 1), common_term_23);
    result_by_q.set(c_Key1(1), common_term_23);
    // q = 2
    result_by_lpq.set(c_Key3(9, 6, 2), common_term_22);
    result_by_q.set(c_Key1(2), common_term_22);
    // q = 3
    result_by_lpq.set(c_Key3(9, 6, 3), common_term_21);
    result_by_q.set(c_Key1(3), common_term_21);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 6), result_by_q);
    result_by_q.clear();

    // l , p = (9, 7).
    // q = -3
    result_by_lpq.set(c_Key3(9, 7, -3), common_term_20);
    result_by_q.set(c_Key1(-3), common_term_20);
    // q = -2
    result_by_lpq.set(c_Key3(9, 7, -2), common_term_19);
    result_by_q.set(c_Key1(-2), common_term_19);
    // q = -1
    result_by_lpq.set(c_Key3(9, 7, -1), common_term_18);
    result_by_q.set(c_Key1(-1), common_term_18);
    // q = 0
    result_by_lpq.set(c_Key3(9, 7, 0), common_term_17);
    result_by_q.set(c_Key1(0), common_term_17);
    // q = 1
    result_by_lpq.set(c_Key3(9, 7, 1), common_term_16);
    result_by_q.set(c_Key1(1), common_term_16);
    // q = 2
    result_by_lpq.set(c_Key3(9, 7, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // q = 3
    result_by_lpq.set(c_Key3(9, 7, 3), common_term_14);
    result_by_q.set(c_Key1(3), common_term_14);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 7), result_by_q);
    result_by_q.clear();

    // l , p = (9, 8).
    // q = -3
    result_by_lpq.set(c_Key3(9, 8, -3), common_term_13);
    result_by_q.set(c_Key1(-3), common_term_13);
    // q = -2
    result_by_lpq.set(c_Key3(9, 8, -2), common_term_12);
    result_by_q.set(c_Key1(-2), common_term_12);
    // q = -1
    result_by_lpq.set(c_Key3(9, 8, -1), common_term_11);
    result_by_q.set(c_Key1(-1), common_term_11);
    // q = 0
    result_by_lpq.set(c_Key3(9, 8, 0), common_term_10);
    result_by_q.set(c_Key1(0), common_term_10);
    // q = 1
    result_by_lpq.set(c_Key3(9, 8, 1), common_term_9);
    result_by_q.set(c_Key1(1), common_term_9);
    // q = 2
    result_by_lpq.set(c_Key3(9, 8, 2), common_term_8);
    result_by_q.set(c_Key1(2), common_term_8);
    // q = 3
    result_by_lpq.set(c_Key3(9, 8, 3), common_term_7);
    result_by_q.set(c_Key1(3), common_term_7);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 8), result_by_q);
    result_by_q.clear();

    // l , p = (9, 9).
    // q = -3
    result_by_lpq.set(c_Key3(9, 9, -3), common_term_6);
    result_by_q.set(c_Key1(-3), common_term_6);
    // q = -2
    result_by_lpq.set(c_Key3(9, 9, -2), common_term_5);
    result_by_q.set(c_Key1(-2), common_term_5);
    // q = -1
    result_by_lpq.set(c_Key3(9, 9, -1), common_term_4);
    result_by_q.set(c_Key1(-1), common_term_4);
    // q = 0
    result_by_lpq.set(c_Key3(9, 9, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(9, 9, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // q = 2
    result_by_lpq.set(c_Key3(9, 9, 2), common_term_1);
    result_by_q.set(c_Key1(2), common_term_1);
    // q = 3
    result_by_lpq.set(c_Key3(9, 9, 3), common_term_0);
    result_by_q.set(c_Key1(3), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 9), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l9_e5(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 9.
    // Functions are taylor expanded and truncated at eccentricity^5.

    // Functions Calculated for l = 9.

    c_IntMap<c_Key3, double> result_by_lpq(92);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(10);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double eccentricity_5 = eccentricity_4 * eccentricity;
    double common_term_0 = 1.6276041666666667*eccentricity_4;
    double common_term_1 = -4.5*eccentricity_3;
    double common_term_2 = -75.541666666666667*eccentricity_4 + 6.125*eccentricity_2;
    double common_term_3 = 95.0*eccentricity_3 - 4.0*eccentricity;
    double common_term_4 = 801.984375*eccentricity_4 - 58.5*eccentricity_2 + 1.0;
    double common_term_5 = -467.5*eccentricity_3 + 14.0*eccentricity;
    double common_term_6 = -2702.7916666666667*eccentricity_4 + 107.375*eccentricity_2;
    double common_term_7 = 597.0*eccentricity_3;
    double common_term_8 = 2689.7213541666667*eccentricity_4;
    double common_term_9 = 0.0546875*eccentricity_4;
    double common_term_10 = -0.33333333333333333*eccentricity_3;
    double common_term_11 = -6.125*eccentricity_4 + 1.375*eccentricity_2;
    double common_term_12 = 20.5*eccentricity_3 - 2.0*eccentricity;
    double common_term_13 = 166.484375*eccentricity_4 - 26.5*eccentricity_2 + 1.0;
    double common_term_14 = -193.0*eccentricity_3 + 12.0*eccentricity;
    double common_term_15 = -1022.875*eccentricity_4 + 80.125*eccentricity_2;
    double common_term_16 = 392.83333333333333*eccentricity_3;
    double common_term_17 = 1577.4609375*eccentricity_4;
    double common_term_18 = 1.75*eccentricity_5*std::pow(1.0 - eccentricity_2, -8.5);
    double common_term_19 = 1.2109375*eccentricity_4;
    double common_term_20 = 0.83333333333333333*eccentricity_3;
    double common_term_21 = 5.625*eccentricity_4 + 0.625*eccentricity_2;
    double common_term_22 = 5.0*eccentricity_3;
    double common_term_23 = 25.859375*eccentricity_4 - 2.5*eccentricity_2 + 1.0;
    double common_term_24 = -27.5*eccentricity_3 + 10.0*eccentricity;
    double common_term_25 = -170.625*eccentricity_4 + 56.875*eccentricity_2;
    double common_term_26 = 241.66666666666667*eccentricity_3;
    double common_term_27 = 852.9296875*eccentricity_4;
    double common_term_28 = 12.096354166666667*eccentricity_4;
    double common_term_29 = std::pow(1.0 - eccentricity_2, -8.5)*(8.75*eccentricity_5 + 7.0*eccentricity_3);
    double common_term_30 = 41.708333333333333*eccentricity_4 + 3.875*eccentricity_2;
    double common_term_31 = 24.5*eccentricity_3 + 2.0*eccentricity;
    double common_term_32 = 92.109375*eccentricity_4 + 13.5*eccentricity_2 + 1.0;
    double common_term_33 = 53.0*eccentricity_3 + 8.0*eccentricity;
    double common_term_34 = 155.95833333333333*eccentricity_4 + 37.625*eccentricity_2;
    double common_term_35 = 135.5*eccentricity_3;
    double common_term_36 = 413.12760416666667*eccentricity_4;
    double common_term_37 = 55.7109375*eccentricity_4;
    double common_term_38 = 26.166666666666667*eccentricity_3;
    double common_term_39 = 120.125*eccentricity_4 + 11.125*eccentricity_2;
    double common_term_40 = std::pow(1.0 - eccentricity_2, -8.5)*(17.5*eccentricity_5 + 21.0*eccentricity_3 + 4.0*eccentricity);
    double common_term_41 = 173.234375*eccentricity_4 + 21.5*eccentricity_2 + 1.0;
    double common_term_42 = 72.5*eccentricity_3 + 6.0*eccentricity;
    double common_term_43 = 194.875*eccentricity_4 + 22.375*eccentricity_2;
    double common_term_44 = 66.333333333333333*eccentricity_3;
    double common_term_45 = 171.0546875*eccentricity_4;

    c_IntMap<c_Key1, double> result_by_q(10);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (9, 0).
    // q = -4
    result_by_lpq.set(c_Key3(9, 0, -4), common_term_0);
    result_by_q.set(c_Key1(-4), common_term_0);
    // q = -3
    result_by_lpq.set(c_Key3(9, 0, -3), common_term_1);
    result_by_q.set(c_Key1(-3), common_term_1);
    // q = -2
    result_by_lpq.set(c_Key3(9, 0, -2), common_term_2);
    result_by_q.set(c_Key1(-2), common_term_2);
    // q = -1
    result_by_lpq.set(c_Key3(9, 0, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    result_by_lpq.set(c_Key3(9, 0, 0), common_term_4);
    result_by_q.set(c_Key1(0), common_term_4);
    // q = 1
    result_by_lpq.set(c_Key3(9, 0, 1), common_term_5);
    result_by_q.set(c_Key1(1), common_term_5);
    // q = 2
    result_by_lpq.set(c_Key3(9, 0, 2), common_term_6);
    result_by_q.set(c_Key1(2), common_term_6);
    // q = 3
    result_by_lpq.set(c_Key3(9, 0, 3), common_term_7);
    result_by_q.set(c_Key1(3), common_term_7);
    // q = 4
    result_by_lpq.set(c_Key3(9, 0, 4), common_term_8);
    result_by_q.set(c_Key1(4), common_term_8);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 0), result_by_q);
    result_by_q.clear();

    // l , p = (9, 1).
    // q = -4
    result_by_lpq.set(c_Key3(9, 1, -4), common_term_9);
    result_by_q.set(c_Key1(-4), common_term_9);
    // q = -3
    result_by_lpq.set(c_Key3(9, 1, -3), common_term_10);
    result_by_q.set(c_Key1(-3), common_term_10);
    // q = -2
    result_by_lpq.set(c_Key3(9, 1, -2), common_term_11);
    result_by_q.set(c_Key1(-2), common_term_11);
    // q = -1
    result_by_lpq.set(c_Key3(9, 1, -1), common_term_12);
    result_by_q.set(c_Key1(-1), common_term_12);
    // q = 0
    result_by_lpq.set(c_Key3(9, 1, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(9, 1, 1), common_term_14);
    result_by_q.set(c_Key1(1), common_term_14);
    // q = 2
    result_by_lpq.set(c_Key3(9, 1, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // q = 3
    result_by_lpq.set(c_Key3(9, 1, 3), common_term_16);
    result_by_q.set(c_Key1(3), common_term_16);
    // q = 4
    result_by_lpq.set(c_Key3(9, 1, 4), common_term_17);
    result_by_q.set(c_Key1(4), common_term_17);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 1), result_by_q);
    result_by_q.clear();

    // l , p = (9, 2).
    // q = -5
    result_by_lpq.set(c_Key3(9, 2, -5), common_term_18);
    result_by_q.set(c_Key1(-5), common_term_18);
    // q = -4
    result_by_lpq.set(c_Key3(9, 2, -4), common_term_19);
    result_by_q.set(c_Key1(-4), common_term_19);
    // q = -3
    result_by_lpq.set(c_Key3(9, 2, -3), common_term_20);
    result_by_q.set(c_Key1(-3), common_term_20);
    // q = -2
    result_by_lpq.set(c_Key3(9, 2, -2), common_term_21);
    result_by_q.set(c_Key1(-2), common_term_21);
    // q = -1
    result_by_lpq.set(c_Key3(9, 2, -1), common_term_22);
    result_by_q.set(c_Key1(-1), common_term_22);
    // q = 0
    result_by_lpq.set(c_Key3(9, 2, 0), common_term_23);
    result_by_q.set(c_Key1(0), common_term_23);
    // q = 1
    result_by_lpq.set(c_Key3(9, 2, 1), common_term_24);
    result_by_q.set(c_Key1(1), common_term_24);
    // q = 2
    result_by_lpq.set(c_Key3(9, 2, 2), common_term_25);
    result_by_q.set(c_Key1(2), common_term_25);
    // q = 3
    result_by_lpq.set(c_Key3(9, 2, 3), common_term_26);
    result_by_q.set(c_Key1(3), common_term_26);
    // q = 4
    result_by_lpq.set(c_Key3(9, 2, 4), common_term_27);
    result_by_q.set(c_Key1(4), common_term_27);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 2), result_by_q);
    result_by_q.clear();

    // l , p = (9, 3).
    // q = -4
    result_by_lpq.set(c_Key3(9, 3, -4), common_term_28);
    result_by_q.set(c_Key1(-4), common_term_28);
    // q = -3
    result_by_lpq.set(c_Key3(9, 3, -3), common_term_29);
    result_by_q.set(c_Key1(-3), common_term_29);
    // q = -2
    result_by_lpq.set(c_Key3(9, 3, -2), common_term_30);
    result_by_q.set(c_Key1(-2), common_term_30);
    // q = -1
    result_by_lpq.set(c_Key3(9, 3, -1), common_term_31);
    result_by_q.set(c_Key1(-1), common_term_31);
    // q = 0
    result_by_lpq.set(c_Key3(9, 3, 0), common_term_32);
    result_by_q.set(c_Key1(0), common_term_32);
    // q = 1
    result_by_lpq.set(c_Key3(9, 3, 1), common_term_33);
    result_by_q.set(c_Key1(1), common_term_33);
    // q = 2
    result_by_lpq.set(c_Key3(9, 3, 2), common_term_34);
    result_by_q.set(c_Key1(2), common_term_34);
    // q = 3
    result_by_lpq.set(c_Key3(9, 3, 3), common_term_35);
    result_by_q.set(c_Key1(3), common_term_35);
    // q = 4
    result_by_lpq.set(c_Key3(9, 3, 4), common_term_36);
    result_by_q.set(c_Key1(4), common_term_36);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 3), result_by_q);
    result_by_q.clear();

    // l , p = (9, 4).
    // q = -4
    result_by_lpq.set(c_Key3(9, 4, -4), common_term_37);
    result_by_q.set(c_Key1(-4), common_term_37);
    // q = -3
    result_by_lpq.set(c_Key3(9, 4, -3), common_term_38);
    result_by_q.set(c_Key1(-3), common_term_38);
    // q = -2
    result_by_lpq.set(c_Key3(9, 4, -2), common_term_39);
    result_by_q.set(c_Key1(-2), common_term_39);
    // q = -1
    result_by_lpq.set(c_Key3(9, 4, -1), common_term_40);
    result_by_q.set(c_Key1(-1), common_term_40);
    // q = 0
    result_by_lpq.set(c_Key3(9, 4, 0), common_term_41);
    result_by_q.set(c_Key1(0), common_term_41);
    // q = 1
    result_by_lpq.set(c_Key3(9, 4, 1), common_term_42);
    result_by_q.set(c_Key1(1), common_term_42);
    // q = 2
    result_by_lpq.set(c_Key3(9, 4, 2), common_term_43);
    result_by_q.set(c_Key1(2), common_term_43);
    // q = 3
    result_by_lpq.set(c_Key3(9, 4, 3), common_term_44);
    result_by_q.set(c_Key1(3), common_term_44);
    // q = 4
    result_by_lpq.set(c_Key3(9, 4, 4), common_term_45);
    result_by_q.set(c_Key1(4), common_term_45);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 4), result_by_q);
    result_by_q.clear();

    // l , p = (9, 5).
    // q = -4
    result_by_lpq.set(c_Key3(9, 5, -4), common_term_45);
    result_by_q.set(c_Key1(-4), common_term_45);
    // q = -3
    result_by_lpq.set(c_Key3(9, 5, -3), common_term_44);
    result_by_q.set(c_Key1(-3), common_term_44);
    // q = -2
    result_by_lpq.set(c_Key3(9, 5, -2), common_term_43);
    result_by_q.set(c_Key1(-2), common_term_43);
    // q = -1
    result_by_lpq.set(c_Key3(9, 5, -1), common_term_42);
    result_by_q.set(c_Key1(-1), common_term_42);
    // q = 0
    result_by_lpq.set(c_Key3(9, 5, 0), common_term_41);
    result_by_q.set(c_Key1(0), common_term_41);
    // q = 1
    result_by_lpq.set(c_Key3(9, 5, 1), common_term_40);
    result_by_q.set(c_Key1(1), common_term_40);
    // q = 2
    result_by_lpq.set(c_Key3(9, 5, 2), common_term_39);
    result_by_q.set(c_Key1(2), common_term_39);
    // q = 3
    result_by_lpq.set(c_Key3(9, 5, 3), common_term_38);
    result_by_q.set(c_Key1(3), common_term_38);
    // q = 4
    result_by_lpq.set(c_Key3(9, 5, 4), common_term_37);
    result_by_q.set(c_Key1(4), common_term_37);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 5), result_by_q);
    result_by_q.clear();

    // l , p = (9, 6).
    // q = -4
    result_by_lpq.set(c_Key3(9, 6, -4), common_term_36);
    result_by_q.set(c_Key1(-4), common_term_36);
    // q = -3
    result_by_lpq.set(c_Key3(9, 6, -3), common_term_35);
    result_by_q.set(c_Key1(-3), common_term_35);
    // q = -2
    result_by_lpq.set(c_Key3(9, 6, -2), common_term_34);
    result_by_q.set(c_Key1(-2), common_term_34);
    // q = -1
    result_by_lpq.set(c_Key3(9, 6, -1), common_term_33);
    result_by_q.set(c_Key1(-1), common_term_33);
    // q = 0
    result_by_lpq.set(c_Key3(9, 6, 0), common_term_32);
    result_by_q.set(c_Key1(0), common_term_32);
    // q = 1
    result_by_lpq.set(c_Key3(9, 6, 1), common_term_31);
    result_by_q.set(c_Key1(1), common_term_31);
    // q = 2
    result_by_lpq.set(c_Key3(9, 6, 2), common_term_30);
    result_by_q.set(c_Key1(2), common_term_30);
    // q = 3
    result_by_lpq.set(c_Key3(9, 6, 3), common_term_29);
    result_by_q.set(c_Key1(3), common_term_29);
    // q = 4
    result_by_lpq.set(c_Key3(9, 6, 4), common_term_28);
    result_by_q.set(c_Key1(4), common_term_28);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 6), result_by_q);
    result_by_q.clear();

    // l , p = (9, 7).
    // q = -4
    result_by_lpq.set(c_Key3(9, 7, -4), common_term_27);
    result_by_q.set(c_Key1(-4), common_term_27);
    // q = -3
    result_by_lpq.set(c_Key3(9, 7, -3), common_term_26);
    result_by_q.set(c_Key1(-3), common_term_26);
    // q = -2
    result_by_lpq.set(c_Key3(9, 7, -2), common_term_25);
    result_by_q.set(c_Key1(-2), common_term_25);
    // q = -1
    result_by_lpq.set(c_Key3(9, 7, -1), common_term_24);
    result_by_q.set(c_Key1(-1), common_term_24);
    // q = 0
    result_by_lpq.set(c_Key3(9, 7, 0), common_term_23);
    result_by_q.set(c_Key1(0), common_term_23);
    // q = 1
    result_by_lpq.set(c_Key3(9, 7, 1), common_term_22);
    result_by_q.set(c_Key1(1), common_term_22);
    // q = 2
    result_by_lpq.set(c_Key3(9, 7, 2), common_term_21);
    result_by_q.set(c_Key1(2), common_term_21);
    // q = 3
    result_by_lpq.set(c_Key3(9, 7, 3), common_term_20);
    result_by_q.set(c_Key1(3), common_term_20);
    // q = 4
    result_by_lpq.set(c_Key3(9, 7, 4), common_term_19);
    result_by_q.set(c_Key1(4), common_term_19);
    // q = 5
    result_by_lpq.set(c_Key3(9, 7, 5), common_term_18);
    result_by_q.set(c_Key1(5), common_term_18);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 7), result_by_q);
    result_by_q.clear();

    // l , p = (9, 8).
    // q = -4
    result_by_lpq.set(c_Key3(9, 8, -4), common_term_17);
    result_by_q.set(c_Key1(-4), common_term_17);
    // q = -3
    result_by_lpq.set(c_Key3(9, 8, -3), common_term_16);
    result_by_q.set(c_Key1(-3), common_term_16);
    // q = -2
    result_by_lpq.set(c_Key3(9, 8, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(9, 8, -1), common_term_14);
    result_by_q.set(c_Key1(-1), common_term_14);
    // q = 0
    result_by_lpq.set(c_Key3(9, 8, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(9, 8, 1), common_term_12);
    result_by_q.set(c_Key1(1), common_term_12);
    // q = 2
    result_by_lpq.set(c_Key3(9, 8, 2), common_term_11);
    result_by_q.set(c_Key1(2), common_term_11);
    // q = 3
    result_by_lpq.set(c_Key3(9, 8, 3), common_term_10);
    result_by_q.set(c_Key1(3), common_term_10);
    // q = 4
    result_by_lpq.set(c_Key3(9, 8, 4), common_term_9);
    result_by_q.set(c_Key1(4), common_term_9);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 8), result_by_q);
    result_by_q.clear();

    // l , p = (9, 9).
    // q = -4
    result_by_lpq.set(c_Key3(9, 9, -4), common_term_8);
    result_by_q.set(c_Key1(-4), common_term_8);
    // q = -3
    result_by_lpq.set(c_Key3(9, 9, -3), common_term_7);
    result_by_q.set(c_Key1(-3), common_term_7);
    // q = -2
    result_by_lpq.set(c_Key3(9, 9, -2), common_term_6);
    result_by_q.set(c_Key1(-2), common_term_6);
    // q = -1
    result_by_lpq.set(c_Key3(9, 9, -1), common_term_5);
    result_by_q.set(c_Key1(-1), common_term_5);
    // q = 0
    result_by_lpq.set(c_Key3(9, 9, 0), common_term_4);
    result_by_q.set(c_Key1(0), common_term_4);
    // q = 1
    result_by_lpq.set(c_Key3(9, 9, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // q = 2
    result_by_lpq.set(c_Key3(9, 9, 2), common_term_2);
    result_by_q.set(c_Key1(2), common_term_2);
    // q = 3
    result_by_lpq.set(c_Key3(9, 9, 3), common_term_1);
    result_by_q.set(c_Key1(3), common_term_1);
    // q = 4
    result_by_lpq.set(c_Key3(9, 9, 4), common_term_0);
    result_by_q.set(c_Key1(4), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 9), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l9_e10(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 9.
    // Functions are taylor expanded and truncated at eccentricity^10.

    // Functions Calculated for l = 9.

    c_IntMap<c_Key3, double> result_by_lpq(188);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(10);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_8 = eccentricity_4 * eccentricity_4;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double eccentricity_6 = eccentricity_3 * eccentricity_3;
    double eccentricity_7 = eccentricity_6 * eccentricity;
    double eccentricity_9 = eccentricity_8 * eccentricity;
    double eccentricity_5 = eccentricity_4 * eccentricity;
    double common_term_0 = 9.6881200396825397e-8*eccentricity_8;
    double common_term_1 = -0.00019841269841269841*eccentricity_9 - 0.00019841269841269841*eccentricity_7;
    double common_term_2 = 0.0158203125*eccentricity_6;
    double common_term_3 = -0.093650793650793651*eccentricity_9 + 0.37777777777777778*eccentricity_7 - 0.26666666666666667*eccentricity_5;
    double common_term_4 = 5.7898627387152778*eccentricity_8 - 5.6966145833333333*eccentricity_6 + 1.6276041666666667*eccentricity_4;
    double common_term_5 = 62.8171875*eccentricity_9 - 66.825*eccentricity_7 + 30.375*eccentricity_5 - 4.5*eccentricity_3;
    double common_term_6 = -557.04748263888889*eccentricity_8 + 306.17024739583333*eccentricity_6 - 75.541666666666667*eccentricity_4 + 6.125*eccentricity_2;
    double common_term_7 = -3664.5597222222222*eccentricity_9 + 2203.2152777777778*eccentricity_7 - 686.83333333333333*eccentricity_5 + 95.0*eccentricity_3 - 4.0*eccentricity;
    double common_term_8 = 12656.252746582031*eccentricity_8 - 4483.0078125*eccentricity_6 + 801.984375*eccentricity_4 - 58.5*eccentricity_2 + 1.0;
    double common_term_9 = 61626.019965277778*eccentricity_9 - 23508.472222222222*eccentricity_7 + 4882.2916666666667*eccentricity_5 - 467.5*eccentricity_3 + 14.0*eccentricity;
    double common_term_10 = -105162.43107638889*eccentricity_8 + 23978.449544270833*eccentricity_6 - 2702.7916666666667*eccentricity_4 + 107.375*eccentricity_2;
    double common_term_11 = -416358.6*eccentricity_9 + 100843.65*eccentricity_7 - 12640.5*eccentricity_5 + 597.0*eccentricity_3;
    double common_term_12 = 376659.80988226997*eccentricity_8 - 50743.262239583333*eccentricity_6 + 2689.7213541666667*eccentricity_4;
    double common_term_13 = 1280231.6090277778*eccentricity_9 - 181294.75972222222*eccentricity_7 + 10416.933333333333*eccentricity_5;
    double common_term_14 = -590554.91350446429*eccentricity_8 + 35952.3310546875*eccentricity_6;
    double common_term_15 = -1784019.4609126984*eccentricity_9 + 113262.81408730159*eccentricity_7;
    double common_term_16 = 331262.53957180447*eccentricity_8;
    double common_term_17 = 910742.26830357143*eccentricity_9;
    double common_term_18 = 0.11545414462081129*eccentricity_9;
    double common_term_19 = 0.084946308438740079*eccentricity_8;
    double common_term_20 = 0.0625*eccentricity_7*std::pow(1.0 - eccentricity_2, -8.5);
    double common_term_21 = 0.40236855158730159*eccentricity_8 + 0.045985243055555556*eccentricity_6;
    double common_term_22 = 1.5158730158730159*eccentricity_9 + 0.30416666666666667*eccentricity_7 + 0.033333333333333333*eccentricity_5;
    double common_term_23 = 1.1750732421875*eccentricity_8 + 0.22109375*eccentricity_6 + 0.0546875*eccentricity_4;
    double common_term_24 = 3.5074074074074074*eccentricity_9 + 0.63333333333333333*eccentricity_7 + 0.83333333333333333*eccentricity_5 - 0.33333333333333333*eccentricity_3;
    double common_term_25 = -1.5199652777777778*eccentricity_8 + 8.9547526041666667*eccentricity_6 - 6.125*eccentricity_4 + 1.375*eccentricity_2;
    double common_term_26 = -35.2265625*eccentricity_9 + 73.875*eccentricity_7 - 59.375*eccentricity_5 + 20.5*eccentricity_3 - 2.0*eccentricity;
    double common_term_27 = 484.53443060980903*eccentricity_8 - 405.70225694444444*eccentricity_6 + 166.484375*eccentricity_4 - 26.5*eccentricity_2 + 1.0;
    double common_term_28 = 2635.5152777777778*eccentricity_9 - 2196.3125*eccentricity_7 + 975.16666666666667*eccentricity_5 - 193.0*eccentricity_3 + 12.0*eccentricity;
    double common_term_29 = -10033.50859375*eccentricity_8 + 4608.4150390625*eccentricity_6 - 1022.875*eccentricity_4 + 80.125*eccentricity_2;
    double common_term_30 = -40225.29369212963*eccentricity_9 + 18651.666666666667*eccentricity_7 - 4409.4583333333333*eccentricity_5 + 392.83333333333333*eccentricity_3;
    double common_term_31 = 67061.771546766493*eccentricity_8 - 16397.329947916667*eccentricity_6 + 1577.4609375*eccentricity_4;
    double common_term_32 = 219503.57142857143*eccentricity_9 - 54513.45*eccentricity_7 + 5494.8*eccentricity_5;
    double common_term_33 = -165903.71507936508*eccentricity_8 + 17190.395334201389*eccentricity_6;
    double common_term_34 = -469949.86909722222*eccentricity_9 + 49423.220833333333*eccentricity_7;
    double common_term_35 = 132700.77527291434*eccentricity_8;
    double common_term_36 = 336685.63266093475*eccentricity_9;
    double common_term_37 = 7.3935736331569665*eccentricity_9;
    double common_term_38 = 5.178314208984375*eccentricity_8;
    double common_term_39 = 29.171378968253968*eccentricity_9 + 3.6175595238095238*eccentricity_7;
    double common_term_40 = 21.172867063492063*eccentricity_8 + 2.5200737847222222*eccentricity_6;
    double common_term_41 = std::pow(1.0 - eccentricity_2, -8.5)*(0.4375*eccentricity_7 + 1.75*eccentricity_5);
    double common_term_42 = 55.334242078993056*eccentricity_8 + 11.032552083333333*eccentricity_6 + 1.2109375*eccentricity_4;
    double common_term_43 = 153.97974537037037*eccentricity_9 + 41.010416666666667*eccentricity_7 + 7.9166666666666667*eccentricity_5 + 0.83333333333333333*eccentricity_3;
    double common_term_44 = 116.73828125*eccentricity_8 + 30.2783203125*eccentricity_6 + 5.625*eccentricity_4 + 0.625*eccentricity_2;
    double common_term_45 = 280.95486111111111*eccentricity_9 + 88.333333333333333*eccentricity_7 + 21.666666666666667*eccentricity_5 + 5.0*eccentricity_3;
    double common_term_46 = 218.55689154730903*eccentricity_8 + 59.943576388888889*eccentricity_6 + 25.859375*eccentricity_4 - 2.5*eccentricity_2 + 1.0;
    double common_term_47 = 490.8359375*eccentricity_9 + 117.96875*eccentricity_7 + 114.375*eccentricity_5 - 27.5*eccentricity_3 + 10.0*eccentricity;
    double common_term_48 = 89.620659722222222*eccentricity_8 + 466.15397135416667*eccentricity_6 - 170.625*eccentricity_4 + 56.875*eccentricity_2;
    double common_term_49 = -584.38425925925926*eccentricity_9 + 1766.2708333333333*eccentricity_7 - 785.41666666666667*eccentricity_5 + 241.66666666666667*eccentricity_3;
    double common_term_50 = 6200.8631591796875*eccentricity_8 - 2985.25390625*eccentricity_6 + 852.9296875*eccentricity_4;
    double common_term_51 = 20183.854786706349*eccentricity_9 - 9906.40625*eccentricity_7 + 2641.7083333333333*eccentricity_5;
    double common_term_52 = -29681.301649305556*eccentricity_8 + 7420.3254123263889*eccentricity_6;
    double common_term_53 = -82084.044642857143*eccentricity_9 + 19313.892857142857*eccentricity_7;
    double common_term_54 = 47283.837234254867*eccentricity_8;
    double common_term_55 = 110067.86361882716*eccentricity_9;
    double common_term_56 = 131.75334821428571*eccentricity_9;
    double common_term_57 = 84.257838173518105*eccentricity_8;
    double common_term_58 = 386.73730158730159*eccentricity_9 + 53.233730158730159*eccentricity_7;
    double common_term_59 = 257.31462053571429*eccentricity_8 + 33.1412109375*eccentricity_6;
    double common_term_60 = 788.78559027777778*eccentricity_9 + 168.60902777777778*eccentricity_7 + 20.258333333333333*eccentricity_5;
    double common_term_61 = 536.18533799913194*eccentricity_8 + 108.48776041666667*eccentricity_6 + 12.096354166666667*eccentricity_4;
    double common_term_62 = std::pow(1.0 - eccentricity_2, -8.5)*(1.3125*eccentricity_7 + 8.75*eccentricity_5 + 7.0*eccentricity_3);
    double common_term_63 = 932.26788194444444*eccentricity_8 + 234.53157552083333*eccentricity_6 + 41.708333333333333*eccentricity_4 + 3.875*eccentricity_2;
    double common_term_64 = 2091.1480902777778*eccentricity_9 + 630.05902777777778*eccentricity_7 + 149.54166666666667*eccentricity_5 + 24.5*eccentricity_3 + 2.0*eccentricity;
    double common_term_65 = 1451.6248168945313*eccentricity_8 + 415.9296875*eccentricity_6 + 92.109375*eccentricity_4 + 13.5*eccentricity_2 + 1.0;
    double common_term_66 = 3004.7777777777778*eccentricity_9 + 986.65277777777778*eccentricity_7 + 266.66666666666667*eccentricity_5 + 53.0*eccentricity_3 + 8.0*eccentricity;
    double common_term_67 = 2093.5160590277778*eccentricity_8 + 656.27571614583333*eccentricity_6 + 155.95833333333333*eccentricity_4 + 37.625*eccentricity_2;
    double common_term_68 = 4083.69375*eccentricity_9 + 1447.70625*eccentricity_7 + 375.0*eccentricity_5 + 135.5*eccentricity_3;
    double common_term_69 = 2966.9605007595486*eccentricity_8 + 757.21588541666667*eccentricity_6 + 413.12760416666667*eccentricity_4;
    double common_term_70 = 5818.1884920634921*eccentricity_9 + 1257.9152777777778*eccentricity_7 + 1122.2833333333333*eccentricity_5;
    double common_term_71 = 1502.6236607142857*eccentricity_8 + 2799.4181640625*eccentricity_6;
    double common_term_72 = 216.5906498015873*eccentricity_9 + 6537.5114087301587*eccentricity_7;
    double common_term_73 = 14485.249283951048*eccentricity_8;
    double common_term_74 = 30744.7*eccentricity_9;
    double common_term_75 = 1185.2050870811287*eccentricity_9;
    double common_term_76 = 683.3218254937066*eccentricity_8;
    double common_term_77 = 2210.6404017857143*eccentricity_9 + 384.77142857142857*eccentricity_7;
    double common_term_78 = 1338.2390873015873*eccentricity_8 + 210.44585503472222*eccentricity_6;
    double common_term_79 = 3311.768253968254*eccentricity_9 + 785.10833333333333*eccentricity_7 + 110.93333333333333*eccentricity_5;
    double common_term_80 = 2027.1390380859375*eccentricity_8 + 443.33671875*eccentricity_6 + 55.7109375*eccentricity_4;
    double common_term_81 = 4376.4506365740741*eccentricity_9 + 1196.6708333333333*eccentricity_7 + 238.45833333333333*eccentricity_5 + 26.166666666666667*eccentricity_3;
    double common_term_82 = 2673.5511284722222*eccentricity_8 + 674.60514322916667*eccentricity_6 + 120.125*eccentricity_4 + 11.125*eccentricity_2;
    double common_term_83 = std::pow(1.0 - eccentricity_2, -8.5)*(2.1875*eccentricity_7 + 17.5*eccentricity_5 + 21.0*eccentricity_3 + 4.0*eccentricity);
    double common_term_84 = 3209.274664984809*eccentricity_8 + 863.58940972222222*eccentricity_6 + 173.234375*eccentricity_4 + 21.5*eccentricity_2 + 1.0;
    double common_term_85 = 6051.5210069444444*eccentricity_9 + 1835.2291666666667*eccentricity_7 + 438.79166666666667*eccentricity_5 + 72.5*eccentricity_3 + 6.0*eccentricity;
    double common_term_86 = 3569.0953125*eccentricity_8 + 971.8896484375*eccentricity_6 + 194.875*eccentricity_4 + 22.375*eccentricity_2;
    double common_term_87 = 6497.5703703703704*eccentricity_9 + 1961.3416666666667*eccentricity_7 + 455.16666666666667*eccentricity_5 + 66.333333333333333*eccentricity_3;
    double common_term_88 = 3692.5153944227431*eccentricity_8 + 963.32942708333333*eccentricity_6 + 171.0546875*eccentricity_4;
    double common_term_89 = 6583.5727678571429*eccentricity_9 + 1891.1*eccentricity_7 + 401.15*eccentricity_5;
    double common_term_90 = 3491.5160590277778*eccentricity_8 + 877.80379774305556*eccentricity_6;
    double common_term_91 = 6112.5556051587302*eccentricity_9 + 1821.7160714285714*eccentricity_7;
    double common_term_92 = 3625.3474897112165*eccentricity_8;
    double common_term_93 = 6972.9033702601411*eccentricity_9;

    c_IntMap<c_Key1, double> result_by_q(19);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (9, 0).
    // q = -8
    result_by_lpq.set(c_Key3(9, 0, -8), common_term_0);
    result_by_q.set(c_Key1(-8), common_term_0);
    // q = -7
    result_by_lpq.set(c_Key3(9, 0, -7), common_term_1);
    result_by_q.set(c_Key1(-7), common_term_1);
    // q = -6
    result_by_lpq.set(c_Key3(9, 0, -6), common_term_2);
    result_by_q.set(c_Key1(-6), common_term_2);
    // q = -5
    result_by_lpq.set(c_Key3(9, 0, -5), common_term_3);
    result_by_q.set(c_Key1(-5), common_term_3);
    // q = -4
    result_by_lpq.set(c_Key3(9, 0, -4), common_term_4);
    result_by_q.set(c_Key1(-4), common_term_4);
    // q = -3
    result_by_lpq.set(c_Key3(9, 0, -3), common_term_5);
    result_by_q.set(c_Key1(-3), common_term_5);
    // q = -2
    result_by_lpq.set(c_Key3(9, 0, -2), common_term_6);
    result_by_q.set(c_Key1(-2), common_term_6);
    // q = -1
    result_by_lpq.set(c_Key3(9, 0, -1), common_term_7);
    result_by_q.set(c_Key1(-1), common_term_7);
    // q = 0
    result_by_lpq.set(c_Key3(9, 0, 0), common_term_8);
    result_by_q.set(c_Key1(0), common_term_8);
    // q = 1
    result_by_lpq.set(c_Key3(9, 0, 1), common_term_9);
    result_by_q.set(c_Key1(1), common_term_9);
    // q = 2
    result_by_lpq.set(c_Key3(9, 0, 2), common_term_10);
    result_by_q.set(c_Key1(2), common_term_10);
    // q = 3
    result_by_lpq.set(c_Key3(9, 0, 3), common_term_11);
    result_by_q.set(c_Key1(3), common_term_11);
    // q = 4
    result_by_lpq.set(c_Key3(9, 0, 4), common_term_12);
    result_by_q.set(c_Key1(4), common_term_12);
    // q = 5
    result_by_lpq.set(c_Key3(9, 0, 5), common_term_13);
    result_by_q.set(c_Key1(5), common_term_13);
    // q = 6
    result_by_lpq.set(c_Key3(9, 0, 6), common_term_14);
    result_by_q.set(c_Key1(6), common_term_14);
    // q = 7
    result_by_lpq.set(c_Key3(9, 0, 7), common_term_15);
    result_by_q.set(c_Key1(7), common_term_15);
    // q = 8
    result_by_lpq.set(c_Key3(9, 0, 8), common_term_16);
    result_by_q.set(c_Key1(8), common_term_16);
    // q = 9
    result_by_lpq.set(c_Key3(9, 0, 9), common_term_17);
    result_by_q.set(c_Key1(9), common_term_17);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 0), result_by_q);
    result_by_q.clear();

    // l , p = (9, 1).
    // q = -9
    result_by_lpq.set(c_Key3(9, 1, -9), common_term_18);
    result_by_q.set(c_Key1(-9), common_term_18);
    // q = -8
    result_by_lpq.set(c_Key3(9, 1, -8), common_term_19);
    result_by_q.set(c_Key1(-8), common_term_19);
    // q = -7
    result_by_lpq.set(c_Key3(9, 1, -7), common_term_20);
    result_by_q.set(c_Key1(-7), common_term_20);
    // q = -6
    result_by_lpq.set(c_Key3(9, 1, -6), common_term_21);
    result_by_q.set(c_Key1(-6), common_term_21);
    // q = -5
    result_by_lpq.set(c_Key3(9, 1, -5), common_term_22);
    result_by_q.set(c_Key1(-5), common_term_22);
    // q = -4
    result_by_lpq.set(c_Key3(9, 1, -4), common_term_23);
    result_by_q.set(c_Key1(-4), common_term_23);
    // q = -3
    result_by_lpq.set(c_Key3(9, 1, -3), common_term_24);
    result_by_q.set(c_Key1(-3), common_term_24);
    // q = -2
    result_by_lpq.set(c_Key3(9, 1, -2), common_term_25);
    result_by_q.set(c_Key1(-2), common_term_25);
    // q = -1
    result_by_lpq.set(c_Key3(9, 1, -1), common_term_26);
    result_by_q.set(c_Key1(-1), common_term_26);
    // q = 0
    result_by_lpq.set(c_Key3(9, 1, 0), common_term_27);
    result_by_q.set(c_Key1(0), common_term_27);
    // q = 1
    result_by_lpq.set(c_Key3(9, 1, 1), common_term_28);
    result_by_q.set(c_Key1(1), common_term_28);
    // q = 2
    result_by_lpq.set(c_Key3(9, 1, 2), common_term_29);
    result_by_q.set(c_Key1(2), common_term_29);
    // q = 3
    result_by_lpq.set(c_Key3(9, 1, 3), common_term_30);
    result_by_q.set(c_Key1(3), common_term_30);
    // q = 4
    result_by_lpq.set(c_Key3(9, 1, 4), common_term_31);
    result_by_q.set(c_Key1(4), common_term_31);
    // q = 5
    result_by_lpq.set(c_Key3(9, 1, 5), common_term_32);
    result_by_q.set(c_Key1(5), common_term_32);
    // q = 6
    result_by_lpq.set(c_Key3(9, 1, 6), common_term_33);
    result_by_q.set(c_Key1(6), common_term_33);
    // q = 7
    result_by_lpq.set(c_Key3(9, 1, 7), common_term_34);
    result_by_q.set(c_Key1(7), common_term_34);
    // q = 8
    result_by_lpq.set(c_Key3(9, 1, 8), common_term_35);
    result_by_q.set(c_Key1(8), common_term_35);
    // q = 9
    result_by_lpq.set(c_Key3(9, 1, 9), common_term_36);
    result_by_q.set(c_Key1(9), common_term_36);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 1), result_by_q);
    result_by_q.clear();

    // l , p = (9, 2).
    // q = -9
    result_by_lpq.set(c_Key3(9, 2, -9), common_term_37);
    result_by_q.set(c_Key1(-9), common_term_37);
    // q = -8
    result_by_lpq.set(c_Key3(9, 2, -8), common_term_38);
    result_by_q.set(c_Key1(-8), common_term_38);
    // q = -7
    result_by_lpq.set(c_Key3(9, 2, -7), common_term_39);
    result_by_q.set(c_Key1(-7), common_term_39);
    // q = -6
    result_by_lpq.set(c_Key3(9, 2, -6), common_term_40);
    result_by_q.set(c_Key1(-6), common_term_40);
    // q = -5
    result_by_lpq.set(c_Key3(9, 2, -5), common_term_41);
    result_by_q.set(c_Key1(-5), common_term_41);
    // q = -4
    result_by_lpq.set(c_Key3(9, 2, -4), common_term_42);
    result_by_q.set(c_Key1(-4), common_term_42);
    // q = -3
    result_by_lpq.set(c_Key3(9, 2, -3), common_term_43);
    result_by_q.set(c_Key1(-3), common_term_43);
    // q = -2
    result_by_lpq.set(c_Key3(9, 2, -2), common_term_44);
    result_by_q.set(c_Key1(-2), common_term_44);
    // q = -1
    result_by_lpq.set(c_Key3(9, 2, -1), common_term_45);
    result_by_q.set(c_Key1(-1), common_term_45);
    // q = 0
    result_by_lpq.set(c_Key3(9, 2, 0), common_term_46);
    result_by_q.set(c_Key1(0), common_term_46);
    // q = 1
    result_by_lpq.set(c_Key3(9, 2, 1), common_term_47);
    result_by_q.set(c_Key1(1), common_term_47);
    // q = 2
    result_by_lpq.set(c_Key3(9, 2, 2), common_term_48);
    result_by_q.set(c_Key1(2), common_term_48);
    // q = 3
    result_by_lpq.set(c_Key3(9, 2, 3), common_term_49);
    result_by_q.set(c_Key1(3), common_term_49);
    // q = 4
    result_by_lpq.set(c_Key3(9, 2, 4), common_term_50);
    result_by_q.set(c_Key1(4), common_term_50);
    // q = 5
    result_by_lpq.set(c_Key3(9, 2, 5), common_term_51);
    result_by_q.set(c_Key1(5), common_term_51);
    // q = 6
    result_by_lpq.set(c_Key3(9, 2, 6), common_term_52);
    result_by_q.set(c_Key1(6), common_term_52);
    // q = 7
    result_by_lpq.set(c_Key3(9, 2, 7), common_term_53);
    result_by_q.set(c_Key1(7), common_term_53);
    // q = 8
    result_by_lpq.set(c_Key3(9, 2, 8), common_term_54);
    result_by_q.set(c_Key1(8), common_term_54);
    // q = 9
    result_by_lpq.set(c_Key3(9, 2, 9), common_term_55);
    result_by_q.set(c_Key1(9), common_term_55);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 2), result_by_q);
    result_by_q.clear();

    // l , p = (9, 3).
    // q = -9
    result_by_lpq.set(c_Key3(9, 3, -9), common_term_56);
    result_by_q.set(c_Key1(-9), common_term_56);
    // q = -8
    result_by_lpq.set(c_Key3(9, 3, -8), common_term_57);
    result_by_q.set(c_Key1(-8), common_term_57);
    // q = -7
    result_by_lpq.set(c_Key3(9, 3, -7), common_term_58);
    result_by_q.set(c_Key1(-7), common_term_58);
    // q = -6
    result_by_lpq.set(c_Key3(9, 3, -6), common_term_59);
    result_by_q.set(c_Key1(-6), common_term_59);
    // q = -5
    result_by_lpq.set(c_Key3(9, 3, -5), common_term_60);
    result_by_q.set(c_Key1(-5), common_term_60);
    // q = -4
    result_by_lpq.set(c_Key3(9, 3, -4), common_term_61);
    result_by_q.set(c_Key1(-4), common_term_61);
    // q = -3
    result_by_lpq.set(c_Key3(9, 3, -3), common_term_62);
    result_by_q.set(c_Key1(-3), common_term_62);
    // q = -2
    result_by_lpq.set(c_Key3(9, 3, -2), common_term_63);
    result_by_q.set(c_Key1(-2), common_term_63);
    // q = -1
    result_by_lpq.set(c_Key3(9, 3, -1), common_term_64);
    result_by_q.set(c_Key1(-1), common_term_64);
    // q = 0
    result_by_lpq.set(c_Key3(9, 3, 0), common_term_65);
    result_by_q.set(c_Key1(0), common_term_65);
    // q = 1
    result_by_lpq.set(c_Key3(9, 3, 1), common_term_66);
    result_by_q.set(c_Key1(1), common_term_66);
    // q = 2
    result_by_lpq.set(c_Key3(9, 3, 2), common_term_67);
    result_by_q.set(c_Key1(2), common_term_67);
    // q = 3
    result_by_lpq.set(c_Key3(9, 3, 3), common_term_68);
    result_by_q.set(c_Key1(3), common_term_68);
    // q = 4
    result_by_lpq.set(c_Key3(9, 3, 4), common_term_69);
    result_by_q.set(c_Key1(4), common_term_69);
    // q = 5
    result_by_lpq.set(c_Key3(9, 3, 5), common_term_70);
    result_by_q.set(c_Key1(5), common_term_70);
    // q = 6
    result_by_lpq.set(c_Key3(9, 3, 6), common_term_71);
    result_by_q.set(c_Key1(6), common_term_71);
    // q = 7
    result_by_lpq.set(c_Key3(9, 3, 7), common_term_72);
    result_by_q.set(c_Key1(7), common_term_72);
    // q = 8
    result_by_lpq.set(c_Key3(9, 3, 8), common_term_73);
    result_by_q.set(c_Key1(8), common_term_73);
    // q = 9
    result_by_lpq.set(c_Key3(9, 3, 9), common_term_74);
    result_by_q.set(c_Key1(9), common_term_74);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 3), result_by_q);
    result_by_q.clear();

    // l , p = (9, 4).
    // q = -9
    result_by_lpq.set(c_Key3(9, 4, -9), common_term_75);
    result_by_q.set(c_Key1(-9), common_term_75);
    // q = -8
    result_by_lpq.set(c_Key3(9, 4, -8), common_term_76);
    result_by_q.set(c_Key1(-8), common_term_76);
    // q = -7
    result_by_lpq.set(c_Key3(9, 4, -7), common_term_77);
    result_by_q.set(c_Key1(-7), common_term_77);
    // q = -6
    result_by_lpq.set(c_Key3(9, 4, -6), common_term_78);
    result_by_q.set(c_Key1(-6), common_term_78);
    // q = -5
    result_by_lpq.set(c_Key3(9, 4, -5), common_term_79);
    result_by_q.set(c_Key1(-5), common_term_79);
    // q = -4
    result_by_lpq.set(c_Key3(9, 4, -4), common_term_80);
    result_by_q.set(c_Key1(-4), common_term_80);
    // q = -3
    result_by_lpq.set(c_Key3(9, 4, -3), common_term_81);
    result_by_q.set(c_Key1(-3), common_term_81);
    // q = -2
    result_by_lpq.set(c_Key3(9, 4, -2), common_term_82);
    result_by_q.set(c_Key1(-2), common_term_82);
    // q = -1
    result_by_lpq.set(c_Key3(9, 4, -1), common_term_83);
    result_by_q.set(c_Key1(-1), common_term_83);
    // q = 0
    result_by_lpq.set(c_Key3(9, 4, 0), common_term_84);
    result_by_q.set(c_Key1(0), common_term_84);
    // q = 1
    result_by_lpq.set(c_Key3(9, 4, 1), common_term_85);
    result_by_q.set(c_Key1(1), common_term_85);
    // q = 2
    result_by_lpq.set(c_Key3(9, 4, 2), common_term_86);
    result_by_q.set(c_Key1(2), common_term_86);
    // q = 3
    result_by_lpq.set(c_Key3(9, 4, 3), common_term_87);
    result_by_q.set(c_Key1(3), common_term_87);
    // q = 4
    result_by_lpq.set(c_Key3(9, 4, 4), common_term_88);
    result_by_q.set(c_Key1(4), common_term_88);
    // q = 5
    result_by_lpq.set(c_Key3(9, 4, 5), common_term_89);
    result_by_q.set(c_Key1(5), common_term_89);
    // q = 6
    result_by_lpq.set(c_Key3(9, 4, 6), common_term_90);
    result_by_q.set(c_Key1(6), common_term_90);
    // q = 7
    result_by_lpq.set(c_Key3(9, 4, 7), common_term_91);
    result_by_q.set(c_Key1(7), common_term_91);
    // q = 8
    result_by_lpq.set(c_Key3(9, 4, 8), common_term_92);
    result_by_q.set(c_Key1(8), common_term_92);
    // q = 9
    result_by_lpq.set(c_Key3(9, 4, 9), common_term_93);
    result_by_q.set(c_Key1(9), common_term_93);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 4), result_by_q);
    result_by_q.clear();

    // l , p = (9, 5).
    // q = -9
    result_by_lpq.set(c_Key3(9, 5, -9), common_term_93);
    result_by_q.set(c_Key1(-9), common_term_93);
    // q = -8
    result_by_lpq.set(c_Key3(9, 5, -8), common_term_92);
    result_by_q.set(c_Key1(-8), common_term_92);
    // q = -7
    result_by_lpq.set(c_Key3(9, 5, -7), common_term_91);
    result_by_q.set(c_Key1(-7), common_term_91);
    // q = -6
    result_by_lpq.set(c_Key3(9, 5, -6), common_term_90);
    result_by_q.set(c_Key1(-6), common_term_90);
    // q = -5
    result_by_lpq.set(c_Key3(9, 5, -5), common_term_89);
    result_by_q.set(c_Key1(-5), common_term_89);
    // q = -4
    result_by_lpq.set(c_Key3(9, 5, -4), common_term_88);
    result_by_q.set(c_Key1(-4), common_term_88);
    // q = -3
    result_by_lpq.set(c_Key3(9, 5, -3), common_term_87);
    result_by_q.set(c_Key1(-3), common_term_87);
    // q = -2
    result_by_lpq.set(c_Key3(9, 5, -2), common_term_86);
    result_by_q.set(c_Key1(-2), common_term_86);
    // q = -1
    result_by_lpq.set(c_Key3(9, 5, -1), common_term_85);
    result_by_q.set(c_Key1(-1), common_term_85);
    // q = 0
    result_by_lpq.set(c_Key3(9, 5, 0), common_term_84);
    result_by_q.set(c_Key1(0), common_term_84);
    // q = 1
    result_by_lpq.set(c_Key3(9, 5, 1), common_term_83);
    result_by_q.set(c_Key1(1), common_term_83);
    // q = 2
    result_by_lpq.set(c_Key3(9, 5, 2), common_term_82);
    result_by_q.set(c_Key1(2), common_term_82);
    // q = 3
    result_by_lpq.set(c_Key3(9, 5, 3), common_term_81);
    result_by_q.set(c_Key1(3), common_term_81);
    // q = 4
    result_by_lpq.set(c_Key3(9, 5, 4), common_term_80);
    result_by_q.set(c_Key1(4), common_term_80);
    // q = 5
    result_by_lpq.set(c_Key3(9, 5, 5), common_term_79);
    result_by_q.set(c_Key1(5), common_term_79);
    // q = 6
    result_by_lpq.set(c_Key3(9, 5, 6), common_term_78);
    result_by_q.set(c_Key1(6), common_term_78);
    // q = 7
    result_by_lpq.set(c_Key3(9, 5, 7), common_term_77);
    result_by_q.set(c_Key1(7), common_term_77);
    // q = 8
    result_by_lpq.set(c_Key3(9, 5, 8), common_term_76);
    result_by_q.set(c_Key1(8), common_term_76);
    // q = 9
    result_by_lpq.set(c_Key3(9, 5, 9), common_term_75);
    result_by_q.set(c_Key1(9), common_term_75);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 5), result_by_q);
    result_by_q.clear();

    // l , p = (9, 6).
    // q = -9
    result_by_lpq.set(c_Key3(9, 6, -9), common_term_74);
    result_by_q.set(c_Key1(-9), common_term_74);
    // q = -8
    result_by_lpq.set(c_Key3(9, 6, -8), common_term_73);
    result_by_q.set(c_Key1(-8), common_term_73);
    // q = -7
    result_by_lpq.set(c_Key3(9, 6, -7), common_term_72);
    result_by_q.set(c_Key1(-7), common_term_72);
    // q = -6
    result_by_lpq.set(c_Key3(9, 6, -6), common_term_71);
    result_by_q.set(c_Key1(-6), common_term_71);
    // q = -5
    result_by_lpq.set(c_Key3(9, 6, -5), common_term_70);
    result_by_q.set(c_Key1(-5), common_term_70);
    // q = -4
    result_by_lpq.set(c_Key3(9, 6, -4), common_term_69);
    result_by_q.set(c_Key1(-4), common_term_69);
    // q = -3
    result_by_lpq.set(c_Key3(9, 6, -3), common_term_68);
    result_by_q.set(c_Key1(-3), common_term_68);
    // q = -2
    result_by_lpq.set(c_Key3(9, 6, -2), common_term_67);
    result_by_q.set(c_Key1(-2), common_term_67);
    // q = -1
    result_by_lpq.set(c_Key3(9, 6, -1), common_term_66);
    result_by_q.set(c_Key1(-1), common_term_66);
    // q = 0
    result_by_lpq.set(c_Key3(9, 6, 0), common_term_65);
    result_by_q.set(c_Key1(0), common_term_65);
    // q = 1
    result_by_lpq.set(c_Key3(9, 6, 1), common_term_64);
    result_by_q.set(c_Key1(1), common_term_64);
    // q = 2
    result_by_lpq.set(c_Key3(9, 6, 2), common_term_63);
    result_by_q.set(c_Key1(2), common_term_63);
    // q = 3
    result_by_lpq.set(c_Key3(9, 6, 3), common_term_62);
    result_by_q.set(c_Key1(3), common_term_62);
    // q = 4
    result_by_lpq.set(c_Key3(9, 6, 4), common_term_61);
    result_by_q.set(c_Key1(4), common_term_61);
    // q = 5
    result_by_lpq.set(c_Key3(9, 6, 5), common_term_60);
    result_by_q.set(c_Key1(5), common_term_60);
    // q = 6
    result_by_lpq.set(c_Key3(9, 6, 6), common_term_59);
    result_by_q.set(c_Key1(6), common_term_59);
    // q = 7
    result_by_lpq.set(c_Key3(9, 6, 7), common_term_58);
    result_by_q.set(c_Key1(7), common_term_58);
    // q = 8
    result_by_lpq.set(c_Key3(9, 6, 8), common_term_57);
    result_by_q.set(c_Key1(8), common_term_57);
    // q = 9
    result_by_lpq.set(c_Key3(9, 6, 9), common_term_56);
    result_by_q.set(c_Key1(9), common_term_56);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 6), result_by_q);
    result_by_q.clear();

    // l , p = (9, 7).
    // q = -9
    result_by_lpq.set(c_Key3(9, 7, -9), common_term_55);
    result_by_q.set(c_Key1(-9), common_term_55);
    // q = -8
    result_by_lpq.set(c_Key3(9, 7, -8), common_term_54);
    result_by_q.set(c_Key1(-8), common_term_54);
    // q = -7
    result_by_lpq.set(c_Key3(9, 7, -7), common_term_53);
    result_by_q.set(c_Key1(-7), common_term_53);
    // q = -6
    result_by_lpq.set(c_Key3(9, 7, -6), common_term_52);
    result_by_q.set(c_Key1(-6), common_term_52);
    // q = -5
    result_by_lpq.set(c_Key3(9, 7, -5), common_term_51);
    result_by_q.set(c_Key1(-5), common_term_51);
    // q = -4
    result_by_lpq.set(c_Key3(9, 7, -4), common_term_50);
    result_by_q.set(c_Key1(-4), common_term_50);
    // q = -3
    result_by_lpq.set(c_Key3(9, 7, -3), common_term_49);
    result_by_q.set(c_Key1(-3), common_term_49);
    // q = -2
    result_by_lpq.set(c_Key3(9, 7, -2), common_term_48);
    result_by_q.set(c_Key1(-2), common_term_48);
    // q = -1
    result_by_lpq.set(c_Key3(9, 7, -1), common_term_47);
    result_by_q.set(c_Key1(-1), common_term_47);
    // q = 0
    result_by_lpq.set(c_Key3(9, 7, 0), common_term_46);
    result_by_q.set(c_Key1(0), common_term_46);
    // q = 1
    result_by_lpq.set(c_Key3(9, 7, 1), common_term_45);
    result_by_q.set(c_Key1(1), common_term_45);
    // q = 2
    result_by_lpq.set(c_Key3(9, 7, 2), common_term_44);
    result_by_q.set(c_Key1(2), common_term_44);
    // q = 3
    result_by_lpq.set(c_Key3(9, 7, 3), common_term_43);
    result_by_q.set(c_Key1(3), common_term_43);
    // q = 4
    result_by_lpq.set(c_Key3(9, 7, 4), common_term_42);
    result_by_q.set(c_Key1(4), common_term_42);
    // q = 5
    result_by_lpq.set(c_Key3(9, 7, 5), common_term_41);
    result_by_q.set(c_Key1(5), common_term_41);
    // q = 6
    result_by_lpq.set(c_Key3(9, 7, 6), common_term_40);
    result_by_q.set(c_Key1(6), common_term_40);
    // q = 7
    result_by_lpq.set(c_Key3(9, 7, 7), common_term_39);
    result_by_q.set(c_Key1(7), common_term_39);
    // q = 8
    result_by_lpq.set(c_Key3(9, 7, 8), common_term_38);
    result_by_q.set(c_Key1(8), common_term_38);
    // q = 9
    result_by_lpq.set(c_Key3(9, 7, 9), common_term_37);
    result_by_q.set(c_Key1(9), common_term_37);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 7), result_by_q);
    result_by_q.clear();

    // l , p = (9, 8).
    // q = -9
    result_by_lpq.set(c_Key3(9, 8, -9), common_term_36);
    result_by_q.set(c_Key1(-9), common_term_36);
    // q = -8
    result_by_lpq.set(c_Key3(9, 8, -8), common_term_35);
    result_by_q.set(c_Key1(-8), common_term_35);
    // q = -7
    result_by_lpq.set(c_Key3(9, 8, -7), common_term_34);
    result_by_q.set(c_Key1(-7), common_term_34);
    // q = -6
    result_by_lpq.set(c_Key3(9, 8, -6), common_term_33);
    result_by_q.set(c_Key1(-6), common_term_33);
    // q = -5
    result_by_lpq.set(c_Key3(9, 8, -5), common_term_32);
    result_by_q.set(c_Key1(-5), common_term_32);
    // q = -4
    result_by_lpq.set(c_Key3(9, 8, -4), common_term_31);
    result_by_q.set(c_Key1(-4), common_term_31);
    // q = -3
    result_by_lpq.set(c_Key3(9, 8, -3), common_term_30);
    result_by_q.set(c_Key1(-3), common_term_30);
    // q = -2
    result_by_lpq.set(c_Key3(9, 8, -2), common_term_29);
    result_by_q.set(c_Key1(-2), common_term_29);
    // q = -1
    result_by_lpq.set(c_Key3(9, 8, -1), common_term_28);
    result_by_q.set(c_Key1(-1), common_term_28);
    // q = 0
    result_by_lpq.set(c_Key3(9, 8, 0), common_term_27);
    result_by_q.set(c_Key1(0), common_term_27);
    // q = 1
    result_by_lpq.set(c_Key3(9, 8, 1), common_term_26);
    result_by_q.set(c_Key1(1), common_term_26);
    // q = 2
    result_by_lpq.set(c_Key3(9, 8, 2), common_term_25);
    result_by_q.set(c_Key1(2), common_term_25);
    // q = 3
    result_by_lpq.set(c_Key3(9, 8, 3), common_term_24);
    result_by_q.set(c_Key1(3), common_term_24);
    // q = 4
    result_by_lpq.set(c_Key3(9, 8, 4), common_term_23);
    result_by_q.set(c_Key1(4), common_term_23);
    // q = 5
    result_by_lpq.set(c_Key3(9, 8, 5), common_term_22);
    result_by_q.set(c_Key1(5), common_term_22);
    // q = 6
    result_by_lpq.set(c_Key3(9, 8, 6), common_term_21);
    result_by_q.set(c_Key1(6), common_term_21);
    // q = 7
    result_by_lpq.set(c_Key3(9, 8, 7), common_term_20);
    result_by_q.set(c_Key1(7), common_term_20);
    // q = 8
    result_by_lpq.set(c_Key3(9, 8, 8), common_term_19);
    result_by_q.set(c_Key1(8), common_term_19);
    // q = 9
    result_by_lpq.set(c_Key3(9, 8, 9), common_term_18);
    result_by_q.set(c_Key1(9), common_term_18);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 8), result_by_q);
    result_by_q.clear();

    // l , p = (9, 9).
    // q = -9
    result_by_lpq.set(c_Key3(9, 9, -9), common_term_17);
    result_by_q.set(c_Key1(-9), common_term_17);
    // q = -8
    result_by_lpq.set(c_Key3(9, 9, -8), common_term_16);
    result_by_q.set(c_Key1(-8), common_term_16);
    // q = -7
    result_by_lpq.set(c_Key3(9, 9, -7), common_term_15);
    result_by_q.set(c_Key1(-7), common_term_15);
    // q = -6
    result_by_lpq.set(c_Key3(9, 9, -6), common_term_14);
    result_by_q.set(c_Key1(-6), common_term_14);
    // q = -5
    result_by_lpq.set(c_Key3(9, 9, -5), common_term_13);
    result_by_q.set(c_Key1(-5), common_term_13);
    // q = -4
    result_by_lpq.set(c_Key3(9, 9, -4), common_term_12);
    result_by_q.set(c_Key1(-4), common_term_12);
    // q = -3
    result_by_lpq.set(c_Key3(9, 9, -3), common_term_11);
    result_by_q.set(c_Key1(-3), common_term_11);
    // q = -2
    result_by_lpq.set(c_Key3(9, 9, -2), common_term_10);
    result_by_q.set(c_Key1(-2), common_term_10);
    // q = -1
    result_by_lpq.set(c_Key3(9, 9, -1), common_term_9);
    result_by_q.set(c_Key1(-1), common_term_9);
    // q = 0
    result_by_lpq.set(c_Key3(9, 9, 0), common_term_8);
    result_by_q.set(c_Key1(0), common_term_8);
    // q = 1
    result_by_lpq.set(c_Key3(9, 9, 1), common_term_7);
    result_by_q.set(c_Key1(1), common_term_7);
    // q = 2
    result_by_lpq.set(c_Key3(9, 9, 2), common_term_6);
    result_by_q.set(c_Key1(2), common_term_6);
    // q = 3
    result_by_lpq.set(c_Key3(9, 9, 3), common_term_5);
    result_by_q.set(c_Key1(3), common_term_5);
    // q = 4
    result_by_lpq.set(c_Key3(9, 9, 4), common_term_4);
    result_by_q.set(c_Key1(4), common_term_4);
    // q = 5
    result_by_lpq.set(c_Key3(9, 9, 5), common_term_3);
    result_by_q.set(c_Key1(5), common_term_3);
    // q = 6
    result_by_lpq.set(c_Key3(9, 9, 6), common_term_2);
    result_by_q.set(c_Key1(6), common_term_2);
    // q = 7
    result_by_lpq.set(c_Key3(9, 9, 7), common_term_1);
    result_by_q.set(c_Key1(7), common_term_1);
    // q = 8
    result_by_lpq.set(c_Key3(9, 9, 8), common_term_0);
    result_by_q.set(c_Key1(8), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 9), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l9_e15(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 9.
    // Functions are taylor expanded and truncated at eccentricity^15.

    // Functions Calculated for l = 9.

    c_IntMap<c_Key3, double> result_by_lpq(288);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(10);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double eccentricity_6 = eccentricity_3 * eccentricity_3;
    double eccentricity_7 = eccentricity_6 * eccentricity;
    double eccentricity_14 = eccentricity_7 * eccentricity_7;
    double eccentricity_12 = eccentricity_6 * eccentricity_6;
    double eccentricity_13 = eccentricity_12 * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_5 = eccentricity_4 * eccentricity;
    double eccentricity_10 = eccentricity_5 * eccentricity_5;
    double eccentricity_11 = eccentricity_10 * eccentricity;
    double eccentricity_8 = eccentricity_4 * eccentricity_4;
    double eccentricity_9 = eccentricity_8 * eccentricity;
    double common_term_0 = 4.2731857291347253e-6*eccentricity_14;
    double common_term_1 = 1.3155568711124267e-6*eccentricity_13;
    double common_term_2 = 8.4385886178030953e-7*eccentricity_14 + 2.708682766208401e-7*eccentricity_12;
    double common_term_3 = 7.3068649457538346e-8*eccentricity_13 + 2.5052108385441719e-8*eccentricity_11;
    double common_term_4 = 1.1892896200333673e-9*eccentricity_14 + 7.0948353825957993e-10*eccentricity_12 + 2.6911444554673721e-10*eccentricity_10;
    double common_term_5 = 2.4010869936942799e-7*eccentricity_14 + 2.1256677267622906e-7*eccentricity_12 + 1.6685095623897707e-7*eccentricity_10 + 9.6881200396825397e-8*eccentricity_8;
    double common_term_6 = -0.00019257973251028807*eccentricity_13 - 0.00019979056437389771*eccentricity_11 - 0.00019841269841269841*eccentricity_9 - 0.00019841269841269841*eccentricity_7;
    double common_term_7 = 0.00076860002108982631*eccentricity_14 + 0.0011388506208147321*eccentricity_12 + 0.0021452767508370536*eccentricity_10 + 0.0158203125*eccentricity_6;
    double common_term_8 = 0.010765726043503821*eccentricity_13 + 0.040013227513227513*eccentricity_11 - 0.093650793650793651*eccentricity_9 + 0.37777777777777778*eccentricity_7 - 0.26666666666666667*eccentricity_5;
    double common_term_9 = -0.035442644045816878*eccentricity_14 + 0.66459573135174141*eccentricity_12 - 2.2781210601645172*eccentricity_10 + 5.7898627387152778*eccentricity_8 - 5.6966145833333333*eccentricity_6 + 1.6276041666666667*eccentricity_4;
    double common_term_10 = 9.3742131696428571*eccentricity_13 - 29.959151785714286*eccentricity_11 + 62.8171875*eccentricity_9 - 66.825*eccentricity_7 + 30.375*eccentricity_5 - 4.5*eccentricity_3;
    double common_term_11 = 99.655787984368295*eccentricity_14 - 279.30555665192781*eccentricity_12 + 519.05630227548105*eccentricity_10 - 557.04748263888889*eccentricity_8 + 306.17024739583333*eccentricity_6 - 75.541666666666667*eccentricity_4 + 6.125*eccentricity_2;
    double common_term_12 = -2049.2777143959439*eccentricity_13 + 3470.4545601851852*eccentricity_11 - 3664.5597222222222*eccentricity_9 + 2203.2152777777778*eccentricity_7 - 686.83333333333333*eccentricity_5 + 95.0*eccentricity_3 - 4.0*eccentricity;
    double common_term_13 = -12537.720637434743*eccentricity_14 + 19631.782016716005*eccentricity_12 - 20179.836617431641*eccentricity_10 + 12656.252746582031*eccentricity_8 - 4483.0078125*eccentricity_6 + 801.984375*eccentricity_4 - 58.5*eccentricity_2 + 1.0;
    double common_term_14 = 97048.000044780642*eccentricity_13 - 96633.478877314817*eccentricity_11 + 61626.019965277778*eccentricity_9 - 23508.472222222222*eccentricity_7 + 4882.2916666666667*eccentricity_5 - 467.5*eccentricity_3 + 14.0*eccentricity;
    double common_term_15 = 429376.97441965474*eccentricity_14 - 413227.00778149963*eccentricity_12 + 263997.47338092946*eccentricity_10 - 105162.43107638889*eccentricity_8 + 23978.449544270833*eccentricity_6 - 2702.7916666666667*eccentricity_4 + 107.375*eccentricity_2;
    double common_term_16 = -1608932.925703125*eccentricity_13 + 1020341.8275669643*eccentricity_11 - 416358.6*eccentricity_9 + 100843.65*eccentricity_7 - 12640.5*eccentricity_5 + 597.0*eccentricity_3;
    double common_term_17 = -5788420.2515178157*eccentricity_14 + 3622721.0416216716*eccentricity_12 - 1495361.3327643824*eccentricity_10 + 376659.80988226997*eccentricity_8 - 50743.262239583333*eccentricity_6 + 2689.7213541666667*eccentricity_4;
    double common_term_18 = 11976103.482719264*eccentricity_13 - 4957815.5905671296*eccentricity_11 + 1280231.6090277778*eccentricity_9 - 181294.75972222222*eccentricity_7 + 10416.933333333333*eccentricity_5;
    double common_term_19 = 37247077.490721536*eccentricity_14 - 15371773.469458989*eccentricity_12 + 4028436.3722256252*eccentricity_10 - 590554.91350446429*eccentricity_8 + 35952.3310546875*eccentricity_6;
    double common_term_20 = -45014191.050091306*eccentricity_13 + 11885168.934584436*eccentricity_11 - 1784019.4609126984*eccentricity_9 + 113262.81408730159*eccentricity_7;
    double common_term_21 = -125471238.23406481*eccentricity_14 + 33196912.336670863*eccentricity_12 - 5061190.124490238*eccentricity_10 + 331262.53957180447*eccentricity_8;
    double common_term_22 = 88450810.42967989*eccentricity_13 - 13613443.384017857*eccentricity_11 + 910742.26830357143*eccentricity_9;
    double common_term_23 = 226175522.16480095*eccentricity_14 - 34977324.016420916*eccentricity_12 + 2376158.6679388382*eccentricity_10;
    double common_term_24 = -86356379.447419904*eccentricity_13 + 5926901.1560890402*eccentricity_11;
    double common_term_25 = -205869261.68385361*eccentricity_14 + 14217267.973173739*eccentricity_12;
    double common_term_26 = 32955001.403697092*eccentricity_13;
    double common_term_27 = 74106571.639012238*eccentricity_14;
    double common_term_28 = 0.53572790204492376*eccentricity_14;
    double common_term_29 = 0.39408009958791209*eccentricity_13;
    double common_term_30 = 2.1016420987547041*eccentricity_14 + 0.28990642380073486*eccentricity_12;
    double common_term_31 = 1.5995791120530704*eccentricity_13 + 0.21328405583613917*eccentricity_11;
    double common_term_32 = 5.3646493616816285*eccentricity_14 + 1.2161206906801694*eccentricity_12 + 0.15692005702427455*eccentricity_10;
    double common_term_33 = 4.177992148669232*eccentricity_13 + 0.92363233024691358*eccentricity_11 + 0.11545414462081129*eccentricity_9;
    double common_term_34 = 11.146105815556053*eccentricity_14 + 3.2491962517295985*eccentricity_12 + 0.70080703654617229*eccentricity_10 + 0.084946308438740079*eccentricity_8;
    double common_term_35 = 0.0625*eccentricity_7*std::pow(1.0 - eccentricity_2, -8.5);
    double common_term_36 = 20.43161983094378*eccentricity_14 + 6.9877937891005348*eccentricity_12 + 1.9572340344625806*eccentricity_10 + 0.40236855158730159*eccentricity_8 + 0.045985243055555556*eccentricity_6;
    double common_term_37 = 16.414319816468254*eccentricity_13 + 5.5202835648148148*eccentricity_11 + 1.5158730158730159*eccentricity_9 + 0.30416666666666667*eccentricity_7 + 0.033333333333333333*eccentricity_5;
    double common_term_38 = 34.379260070664542*eccentricity_14 + 13.167109646115984*eccentricity_12 + 4.3550973074776786*eccentricity_10 + 1.1750732421875*eccentricity_8 + 0.22109375*eccentricity_6 + 0.0546875*eccentricity_4;
    double common_term_39 = 27.94325355489418*eccentricity_13 + 10.544762731481481*eccentricity_11 + 3.5074074074074074*eccentricity_9 + 0.63333333333333333*eccentricity_7 + 0.83333333333333333*eccentricity_5 - 0.33333333333333333*eccentricity_3;
    double common_term_40 = 54.39085740854921*eccentricity_14 + 22.467579957669374*eccentricity_12 + 9.7438227335611979*eccentricity_10 - 1.5199652777777778*eccentricity_8 + 8.9547526041666667*eccentricity_6 - 6.125*eccentricity_4 + 1.375*eccentricity_2;
    double common_term_41 = 41.150580357142857*eccentricity_13 + 33.4621875*eccentricity_11 - 35.2265625*eccentricity_9 + 73.875*eccentricity_7 - 59.375*eccentricity_5 + 20.5*eccentricity_3 - 2.0*eccentricity;
    double common_term_42 = 45.619802982424512*eccentricity_14 + 165.04652462193995*eccentricity_12 - 295.80602715386285*eccentricity_10 + 484.53443060980903*eccentricity_8 - 405.70225694444444*eccentricity_6 + 166.484375*eccentricity_4 - 26.5*eccentricity_2 + 1.0;
    double common_term_43 = 935.22048776455026*eccentricity_13 - 1826.0444675925926*eccentricity_11 + 2635.5152777777778*eccentricity_9 - 2196.3125*eccentricity_7 + 975.16666666666667*eccentricity_5 - 193.0*eccentricity_3 + 12.0*eccentricity;
    double common_term_44 = 5009.9082951394154*eccentricity_14 - 9351.8743494088309*eccentricity_12 + 12333.765953063965*eccentricity_10 - 10033.50859375*eccentricity_8 + 4608.4150390625*eccentricity_6 - 1022.875*eccentricity_4 + 80.125*eccentricity_2;
    double common_term_45 = -41759.86802972057*eccentricity_13 + 51101.179108796296*eccentricity_11 - 40225.29369212963*eccentricity_9 + 18651.666666666667*eccentricity_7 - 4409.4583333333333*eccentricity_5 + 392.83333333333333*eccentricity_3;
    double common_term_46 = -167303.10199640734*eccentricity_14 + 191554.71031180614*eccentricity_12 - 145321.01616462829*eccentricity_10 + 67061.771546766493*eccentricity_8 - 16397.329947916667*eccentricity_6 + 1577.4609375*eccentricity_4;
    double common_term_47 = 660498.54988839286*eccentricity_13 - 482114.6390625*eccentricity_11 + 219503.57142857143*eccentricity_9 - 54513.45*eccentricity_7 + 5494.8*eccentricity_5;
    double common_term_48 = 2121958.5393176698*eccentricity_14 - 1489592.6070817632*eccentricity_12 + 665473.88377873254*eccentricity_10 - 165903.71507936508*eccentricity_8 + 17190.395334201389*eccentricity_6;
    double common_term_49 = -4332675.3715330826*eccentricity_13 + 1892730.2281442901*eccentricity_11 - 469949.86909722222*eccentricity_9 + 49423.220833333333*eccentricity_7;
    double common_term_50 = -11964284.546877746*eccentricity_14 + 5099659.6063139098*eccentricity_12 - 1254253.8821781703*eccentricity_10 + 132700.77527291434*eccentricity_8;
    double common_term_51 = 13115741.563279045*eccentricity_13 - 3183256.661253858*eccentricity_11 + 336685.63266093475*eccentricity_9;
    double common_term_52 = 32395259.567587899*eccentricity_14 - 7738082.4555637731*eccentricity_12 + 814424.14075698691*eccentricity_10;
    double common_term_53 = -18119802.234710836*eccentricity_13 + 1891272.7499431818*eccentricity_11;
    double common_term_54 = -41062126.957570431*eccentricity_14 + 4239530.1961378482*eccentricity_12;
    double common_term_55 = 9214381.6287820975*eccentricity_13;
    double common_term_56 = 19488502.621578406*eccentricity_14;
    double common_term_57 = 42.547407043743222*eccentricity_14;
    double common_term_58 = 30.087678469774824*eccentricity_13;
    double common_term_59 = 137.75201855709393*eccentricity_14 + 21.242784448513769*eccentricity_12;
    double common_term_60 = 101.64655032467532*eccentricity_13 + 14.971971387987013*eccentricity_11;
    double common_term_61 = 309.25502990227388*eccentricity_14 + 74.762625074241336*eccentricity_12 + 10.532229771689763*eccentricity_10;
    double common_term_62 = 234.17501628387045*eccentricity_13 + 54.811539627425044*eccentricity_11 + 7.3935736331569665*eccentricity_9;
    double common_term_63 = 582.69145248524554*eccentricity_14 + 176.71940148217337*eccentricity_12 + 40.053403581891741*eccentricity_10 + 5.178314208984375*eccentricity_8;
    double common_term_64 = 449.59552228009259*eccentricity_13 + 132.90198894951499*eccentricity_11 + 29.171378968253968*eccentricity_9 + 3.6175595238095238*eccentricity_7;
    double common_term_65 = 986.9559695809151*eccentricity_14 + 345.73604182983973*eccentricity_12 + 99.599545130653987*eccentricity_10 + 21.172867063492063*eccentricity_8 + 2.5200737847222222*eccentricity_6;
    double common_term_66 = std::pow(1.0 - eccentricity_2, -8.5)*(0.4375*eccentricity_7 + 1.75*eccentricity_5);
    double common_term_67 = 1554.2469358595731*eccentricity_14 + 603.14970095321615*eccentricity_12 + 202.3493410140749*eccentricity_10 + 55.334242078993056*eccentricity_8 + 11.032552083333333*eccentricity_6 + 1.2109375*eccentricity_4;
    double common_term_68 = 1231.7354021990741*eccentricity_13 + 469.17746569113757*eccentricity_11 + 153.97974537037037*eccentricity_9 + 41.010416666666667*eccentricity_7 + 7.9166666666666667*eccentricity_5 + 0.83333333333333333*eccentricity_3;
    double common_term_69 = 2320.1474348672799*eccentricity_14 + 973.09070412772042*eccentricity_12 + 363.72112274169922*eccentricity_10 + 116.73828125*eccentricity_8 + 30.2783203125*eccentricity_6 + 5.625*eccentricity_4 + 0.625*eccentricity_2;
    double common_term_70 = 1857.4591352513228*eccentricity_13 + 766.29846643518519*eccentricity_11 + 280.95486111111111*eccentricity_9 + 88.333333333333333*eccentricity_7 + 21.666666666666667*eccentricity_5 + 5.0*eccentricity_3;
    double common_term_71 = 3323.7105976727304*eccentricity_14 + 1482.6555449285625*eccentricity_12 + 600.97017076280382*eccentricity_10 + 218.55689154730903*eccentricity_8 + 59.943576388888889*eccentricity_6 + 25.859375*eccentricity_4 - 2.5*eccentricity_2 + 1.0;
    double common_term_72 = 2685.3772600446429*eccentricity_13 + 1174.148046875*eccentricity_11 + 490.8359375*eccentricity_9 + 117.96875*eccentricity_7 + 114.375*eccentricity_5 - 27.5*eccentricity_3 + 10.0*eccentricity;
    double common_term_73 = 4618.2684950833703*eccentricity_14 + 2116.8987029464157*eccentricity_12 + 1071.9672554863824*eccentricity_10 + 89.620659722222222*eccentricity_8 + 466.15397135416667*eccentricity_6 - 170.625*eccentricity_4 + 56.875*eccentricity_2;
    double common_term_74 = 3468.0038804150132*eccentricity_13 + 2473.6350033068783*eccentricity_11 - 584.38425925925926*eccentricity_9 + 1766.2708333333333*eccentricity_7 - 785.41666666666667*eccentricity_5 + 241.66666666666667*eccentricity_3;
    double common_term_75 = 4739.7066311665944*eccentricity_14 + 6444.1837419441768*eccentricity_12 - 4073.581551688058*eccentricity_10 + 6200.8631591796875*eccentricity_8 - 2985.25390625*eccentricity_6 + 852.9296875*eccentricity_4;
    double common_term_76 = 18883.091474798832*eccentricity_13 - 17599.995143022487*eccentricity_11 + 20183.854786706349*eccentricity_9 - 9906.40625*eccentricity_7 + 2641.7083333333333*eccentricity_5;
    double common_term_77 = 58598.085660153467*eccentricity_14 - 62887.841761699499*eccentricity_12 + 61227.443857889327*eccentricity_10 - 29681.301649305556*eccentricity_8 + 7420.3254123263889*eccentricity_6;
    double common_term_78 = -200504.915625*eccentricity_13 + 174330.13950892857*eccentricity_11 - 82084.044642857143*eccentricity_9 + 19313.892857142857*eccentricity_7;
    double common_term_79 = -589381.7735199361*eccentricity_14 + 469335.64014766145*eccentricity_12 - 212777.2675541469*eccentricity_10 + 47283.837234254867*eccentricity_8;
    double common_term_80 = 1202919.9477901848*eccentricity_13 - 522822.35218942901*eccentricity_11 + 110067.86361882716*eccentricity_9;
    double common_term_81 = 2952785.6003944421*eccentricity_14 - 1228061.0608972822*eccentricity_12 + 245612.21217945644*eccentricity_10;
    double common_term_82 = -2775681.8376196865*eccentricity_13 + 528701.3024037498*eccentricity_11;
    double common_term_83 = -6068176.5279497196*eccentricity_14 + 1103304.8232635854*eccentricity_12;
    double common_term_84 = 2240972.8308000593*eccentricity_13;
    double common_term_85 = 4444749.3334428413*eccentricity_14;
    double common_term_86 = 1084.3946115033133*eccentricity_14;
    double common_term_87 = 720.94310661568691*eccentricity_13;
    double common_term_88 = 2495.6297662454838*eccentricity_14 + 476.53095888992409*eccentricity_12;
    double common_term_89 = 1753.0762655677977*eccentricity_13 + 312.91185608565817*eccentricity_11;
    double common_term_90 = 4514.2273971522006*eccentricity_14 + 1220.1931485086392*eccentricity_12 + 203.93073070573218*eccentricity_10;
    double common_term_91 = 3250.4038735288149*eccentricity_13 + 841.059140625*eccentricity_11 + 131.75334821428571*eccentricity_9;
    double common_term_92 = 7167.0033690675334*eccentricity_14 + 2318.8363587372728*eccentricity_12 + 573.65326247980565*eccentricity_10 + 84.257838173518105*eccentricity_8;
    double common_term_93 = 5245.5033547361846*eccentricity_13 + 1637.6462783840388*eccentricity_11 + 386.73730158730159*eccentricity_9 + 53.233730158730159*eccentricity_7;
    double common_term_94 = 10494.656793201481*eccentricity_14 + 3803.0297357831682*eccentricity_12 + 1143.7223088945661*eccentricity_10 + 257.31462053571429*eccentricity_8 + 33.1412109375*eccentricity_6;
    double common_term_95 = 7770.94570737342*eccentricity_13 + 2728.5433428406085*eccentricity_11 + 788.78559027777778*eccentricity_9 + 168.60902777777778*eccentricity_7 + 20.258333333333333*eccentricity_5;
    double common_term_96 = 14527.494697299065*eccentricity_14 + 5698.4299680742637*eccentricity_12 + 1934.7809803473256*eccentricity_10 + 536.18533799913194*eccentricity_8 + 108.48776041666667*eccentricity_6 + 12.096354166666667*eccentricity_4;
    double common_term_97 = std::pow(1.0 - eccentricity_2, -8.5)*(1.3125*eccentricity_7 + 8.75*eccentricity_5 + 7.0*eccentricity_3);
    double common_term_98 = 19285.542654869672*eccentricity_14 + 8022.3546295828037*eccentricity_12 + 2961.233682985659*eccentricity_10 + 932.26788194444444*eccentricity_8 + 234.53157552083333*eccentricity_6 + 41.708333333333333*eccentricity_4 + 3.875*eccentricity_2;
    double common_term_99 = 14496.956052551808*eccentricity_13 + 5864.1605237268518*eccentricity_11 + 2091.1480902777778*eccentricity_9 + 630.05902777777778*eccentricity_7 + 149.54166666666667*eccentricity_5 + 24.5*eccentricity_3 + 2.0*eccentricity;
    double common_term_100 = 24777.936440142807*eccentricity_14 + 10783.215084037781*eccentricity_12 + 4230.4302575683594*eccentricity_10 + 1451.6248168945313*eccentricity_8 + 415.9296875*eccentricity_6 + 92.109375*eccentricity_4 + 13.5*eccentricity_2 + 1.0;
    double common_term_101 = 18715.021332396384*eccentricity_13 + 7924.5632581018519*eccentricity_11 + 3004.7777777777778*eccentricity_9 + 986.65277777777778*eccentricity_7 + 266.66666666666667*eccentricity_5 + 53.0*eccentricity_3 + 8.0*eccentricity;
    double common_term_102 = 31002.333077499843*eccentricity_14 + 13979.906038032007*eccentricity_12 + 5742.2284220942744*eccentricity_10 + 2093.5160590277778*eccentricity_8 + 656.27571614583333*eccentricity_6 + 155.95833333333333*eccentricity_4 + 37.625*eccentricity_2;
    double common_term_103 = 23496.982960379464*eccentricity_13 + 10311.072488839286*eccentricity_11 + 4083.69375*eccentricity_9 + 1447.70625*eccentricity_7 + 375.0*eccentricity_5 + 135.5*eccentricity_3;
    double common_term_104 = 37941.912039481642*eccentricity_14 + 17613.32027336933*eccentricity_12 + 7443.6734243887442*eccentricity_10 + 2966.9605007595486*eccentricity_8 + 757.21588541666667*eccentricity_6 + 413.12760416666667*eccentricity_4;
    double common_term_105 = 28895.653577054857*eccentricity_13 + 12785.68521412037*eccentricity_11 + 5818.1884920634921*eccentricity_9 + 1257.9152777777778*eccentricity_7 + 1122.2833333333333*eccentricity_5;
    double common_term_106 = 45921.587596170391*eccentricity_14 + 20701.191739000593*eccentricity_12 + 11220.506275177002*eccentricity_10 + 1502.6236607142857*eccentricity_8 + 2799.4181640625*eccentricity_6;
    double common_term_107 = 31226.856321591619*eccentricity_13 + 21819.054859733245*eccentricity_11 + 216.5906498015873*eccentricity_9 + 6537.5114087301587*eccentricity_7;
    double common_term_108 = 42294.035758734097*eccentricity_14 + 43575.455233026972*eccentricity_12 - 6033.3518224914658*eccentricity_10 + 14485.249283951048*eccentricity_8;
    double common_term_109 = 90006.839086850649*eccentricity_13 - 25460.433883928571*eccentricity_11 + 30744.7*eccentricity_9;
    double common_term_110 = 191343.08058544473*eccentricity_14 - 76296.55152307011*eccentricity_12 + 62957.192477849997*eccentricity_10;
    double common_term_111 = -196923.72650306596*eccentricity_13 + 125062.34629647166*eccentricity_11;
    double common_term_112 = -464957.44270754158*eccentricity_14 + 242034.08835602736*eccentricity_12;
    double common_term_113 = 457914.39447379427*eccentricity_13;
    double common_term_114 = 849297.60792072932*eccentricity_14;
    double common_term_115 = 14612.741020659801*eccentricity_14;
    double common_term_116 = 9064.6274672202797*eccentricity_13;
    double common_term_117 = 19120.129839481894*eccentricity_14 + 5563.0408091910035*eccentricity_12;
    double common_term_118 = 12940.394119377681*eccentricity_13 + 3372.530453029802*eccentricity_11;
    double common_term_119 = 26701.327805336559*eccentricity_14 + 8587.7030786866027*eccentricity_12 + 2015.861515611921*eccentricity_10;
    double common_term_120 = 18277.500904193222*eccentricity_13 + 5585.2709025022046*eccentricity_11 + 1185.2050870811287*eccentricity_9;
    double common_term_121 = 34123.268146539057*eccentricity_14 + 12302.504021480466*eccentricity_12 + 3555.4650755705657*eccentricity_10 + 683.3218254937066*eccentricity_8;
    double common_term_122 = 23581.400142299107*eccentricity_13 + 8125.6586495535714*eccentricity_11 + 2210.6404017857143*eccentricity_9 + 384.77142857142857*eccentricity_7;
    double common_term_123 = 41421.029536145302*eccentricity_14 + 16010.410773990997*eccentricity_12 + 5252.626898326571*eccentricity_10 + 1338.2390873015873*eccentricity_8 + 210.44585503472222*eccentricity_6;
    double common_term_124 = 28747.072789214065*eccentricity_13 + 10652.907510747354*eccentricity_11 + 3311.768253968254*eccentricity_9 + 785.10833333333333*eccentricity_7 + 110.93333333333333*eccentricity_5;
    double common_term_125 = 48356.915566487312*eccentricity_14 + 19576.206907755988*eccentricity_12 + 6923.8570155552455*eccentricity_10 + 2027.1390380859375*eccentricity_8 + 443.33671875*eccentricity_6 + 55.7109375*eccentricity_4;
    double common_term_126 = 33579.477349537037*eccentricity_13 + 13040.228257275132*eccentricity_11 + 4376.4506365740741*eccentricity_9 + 1196.6708333333333*eccentricity_7 + 238.45833333333333*eccentricity_5 + 26.166666666666667*eccentricity_3;
    double common_term_127 = 54705.23853801273*eccentricity_14 + 22838.89973539968*eccentricity_12 + 8461.7476655748155*eccentricity_10 + 2673.5511284722222*eccentricity_8 + 674.60514322916667*eccentricity_6 + 120.125*eccentricity_4 + 11.125*eccentricity_2;
    double common_term_128 = std::pow(1.0 - eccentricity_2, -8.5)*(2.1875*eccentricity_7 + 17.5*eccentricity_5 + 21.0*eccentricity_3 + 4.0*eccentricity);
    double common_term_129 = 60245.641232276541*eccentricity_14 + 25641.887269770775*eccentricity_12 + 9759.8455874294705*eccentricity_10 + 3209.274664984809*eccentricity_8 + 863.58940972222222*eccentricity_6 + 173.234375*eccentricity_4 + 21.5*eccentricity_2 + 1.0;
    double common_term_130 = 41493.096462673611*eccentricity_13 + 16871.826087962963*eccentricity_11 + 6051.5210069444444*eccentricity_9 + 1835.2291666666667*eccentricity_7 + 438.79166666666667*eccentricity_5 + 72.5*eccentricity_3 + 6.0*eccentricity;
    double common_term_131 = 64764.59881095767*eccentricity_14 + 27834.041443437849*eccentricity_12 + 10715.910200500488*eccentricity_10 + 3569.0953125*eccentricity_8 + 971.8896484375*eccentricity_6 + 194.875*eccentricity_4 + 22.375*eccentricity_2;
    double common_term_132 = 44213.601309937169*eccentricity_13 + 18066.265777943122*eccentricity_11 + 6497.5703703703704*eccentricity_9 + 1961.3416666666667*eccentricity_7 + 455.16666666666667*eccentricity_5 + 66.333333333333333*eccentricity_3;
    double common_term_133 = 68056.858297803344*eccentricity_14 + 29271.151923332265*eccentricity_12 + 11233.310910179501*eccentricity_10 + 3692.5153944227431*eccentricity_8 + 963.32942708333333*eccentricity_6 + 171.0546875*eccentricity_4;
    double common_term_134 = 45882.353724888393*eccentricity_13 + 18623.578404017857*eccentricity_11 + 6583.5727678571429*eccentricity_9 + 1891.1*eccentricity_7 + 401.15*eccentricity_5;
    double common_term_135 = 69927.299509688076*eccentricity_14 + 29814.914947921062*eccentricity_12 + 11233.090602535672*eccentricity_10 + 3491.5160590277778*eccentricity_8 + 877.80379774305556*eccentricity_6;
    double common_term_136 = 46327.005182016093*eccentricity_13 + 18485.898155037478*eccentricity_11 + 6112.5556051587302*eccentricity_9 + 1821.7160714285714*eccentricity_7;
    double common_term_137 = 70129.440291661968*eccentricity_14 + 29532.771608513423*eccentricity_12 + 10188.909559413365*eccentricity_10 + 3625.3474897112165*eccentricity_8;
    double common_term_138 = 46078.888808643854*eccentricity_13 + 16179.124872547399*eccentricity_11 + 6972.9033702601411*eccentricity_9;
    double common_term_139 = 70655.822435672789*eccentricity_14 + 24385.021060059293*eccentricity_12 + 13037.374969302384*eccentricity_10;
    double common_term_140 = 34536.848508522727*eccentricity_13 + 23800.945304383117*eccentricity_11;
    double common_term_141 = 44935.450212509978*eccentricity_14 + 42571.446420460595*eccentricity_12;
    double common_term_142 = 74808.31916897676*eccentricity_13;
    double common_term_143 = 129434.23875676698*eccentricity_14;

    c_IntMap<c_Key1, double> result_by_q(29);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (9, 0).
    // q = -14
    result_by_lpq.set(c_Key3(9, 0, -14), common_term_0);
    result_by_q.set(c_Key1(-14), common_term_0);
    // q = -13
    result_by_lpq.set(c_Key3(9, 0, -13), common_term_1);
    result_by_q.set(c_Key1(-13), common_term_1);
    // q = -12
    result_by_lpq.set(c_Key3(9, 0, -12), common_term_2);
    result_by_q.set(c_Key1(-12), common_term_2);
    // q = -11
    result_by_lpq.set(c_Key3(9, 0, -11), common_term_3);
    result_by_q.set(c_Key1(-11), common_term_3);
    // q = -10
    result_by_lpq.set(c_Key3(9, 0, -10), common_term_4);
    result_by_q.set(c_Key1(-10), common_term_4);
    // q = -8
    result_by_lpq.set(c_Key3(9, 0, -8), common_term_5);
    result_by_q.set(c_Key1(-8), common_term_5);
    // q = -7
    result_by_lpq.set(c_Key3(9, 0, -7), common_term_6);
    result_by_q.set(c_Key1(-7), common_term_6);
    // q = -6
    result_by_lpq.set(c_Key3(9, 0, -6), common_term_7);
    result_by_q.set(c_Key1(-6), common_term_7);
    // q = -5
    result_by_lpq.set(c_Key3(9, 0, -5), common_term_8);
    result_by_q.set(c_Key1(-5), common_term_8);
    // q = -4
    result_by_lpq.set(c_Key3(9, 0, -4), common_term_9);
    result_by_q.set(c_Key1(-4), common_term_9);
    // q = -3
    result_by_lpq.set(c_Key3(9, 0, -3), common_term_10);
    result_by_q.set(c_Key1(-3), common_term_10);
    // q = -2
    result_by_lpq.set(c_Key3(9, 0, -2), common_term_11);
    result_by_q.set(c_Key1(-2), common_term_11);
    // q = -1
    result_by_lpq.set(c_Key3(9, 0, -1), common_term_12);
    result_by_q.set(c_Key1(-1), common_term_12);
    // q = 0
    result_by_lpq.set(c_Key3(9, 0, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(9, 0, 1), common_term_14);
    result_by_q.set(c_Key1(1), common_term_14);
    // q = 2
    result_by_lpq.set(c_Key3(9, 0, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // q = 3
    result_by_lpq.set(c_Key3(9, 0, 3), common_term_16);
    result_by_q.set(c_Key1(3), common_term_16);
    // q = 4
    result_by_lpq.set(c_Key3(9, 0, 4), common_term_17);
    result_by_q.set(c_Key1(4), common_term_17);
    // q = 5
    result_by_lpq.set(c_Key3(9, 0, 5), common_term_18);
    result_by_q.set(c_Key1(5), common_term_18);
    // q = 6
    result_by_lpq.set(c_Key3(9, 0, 6), common_term_19);
    result_by_q.set(c_Key1(6), common_term_19);
    // q = 7
    result_by_lpq.set(c_Key3(9, 0, 7), common_term_20);
    result_by_q.set(c_Key1(7), common_term_20);
    // q = 8
    result_by_lpq.set(c_Key3(9, 0, 8), common_term_21);
    result_by_q.set(c_Key1(8), common_term_21);
    // q = 9
    result_by_lpq.set(c_Key3(9, 0, 9), common_term_22);
    result_by_q.set(c_Key1(9), common_term_22);
    // q = 10
    result_by_lpq.set(c_Key3(9, 0, 10), common_term_23);
    result_by_q.set(c_Key1(10), common_term_23);
    // q = 11
    result_by_lpq.set(c_Key3(9, 0, 11), common_term_24);
    result_by_q.set(c_Key1(11), common_term_24);
    // q = 12
    result_by_lpq.set(c_Key3(9, 0, 12), common_term_25);
    result_by_q.set(c_Key1(12), common_term_25);
    // q = 13
    result_by_lpq.set(c_Key3(9, 0, 13), common_term_26);
    result_by_q.set(c_Key1(13), common_term_26);
    // q = 14
    result_by_lpq.set(c_Key3(9, 0, 14), common_term_27);
    result_by_q.set(c_Key1(14), common_term_27);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 0), result_by_q);
    result_by_q.clear();

    // l , p = (9, 1).
    // q = -14
    result_by_lpq.set(c_Key3(9, 1, -14), common_term_28);
    result_by_q.set(c_Key1(-14), common_term_28);
    // q = -13
    result_by_lpq.set(c_Key3(9, 1, -13), common_term_29);
    result_by_q.set(c_Key1(-13), common_term_29);
    // q = -12
    result_by_lpq.set(c_Key3(9, 1, -12), common_term_30);
    result_by_q.set(c_Key1(-12), common_term_30);
    // q = -11
    result_by_lpq.set(c_Key3(9, 1, -11), common_term_31);
    result_by_q.set(c_Key1(-11), common_term_31);
    // q = -10
    result_by_lpq.set(c_Key3(9, 1, -10), common_term_32);
    result_by_q.set(c_Key1(-10), common_term_32);
    // q = -9
    result_by_lpq.set(c_Key3(9, 1, -9), common_term_33);
    result_by_q.set(c_Key1(-9), common_term_33);
    // q = -8
    result_by_lpq.set(c_Key3(9, 1, -8), common_term_34);
    result_by_q.set(c_Key1(-8), common_term_34);
    // q = -7
    result_by_lpq.set(c_Key3(9, 1, -7), common_term_35);
    result_by_q.set(c_Key1(-7), common_term_35);
    // q = -6
    result_by_lpq.set(c_Key3(9, 1, -6), common_term_36);
    result_by_q.set(c_Key1(-6), common_term_36);
    // q = -5
    result_by_lpq.set(c_Key3(9, 1, -5), common_term_37);
    result_by_q.set(c_Key1(-5), common_term_37);
    // q = -4
    result_by_lpq.set(c_Key3(9, 1, -4), common_term_38);
    result_by_q.set(c_Key1(-4), common_term_38);
    // q = -3
    result_by_lpq.set(c_Key3(9, 1, -3), common_term_39);
    result_by_q.set(c_Key1(-3), common_term_39);
    // q = -2
    result_by_lpq.set(c_Key3(9, 1, -2), common_term_40);
    result_by_q.set(c_Key1(-2), common_term_40);
    // q = -1
    result_by_lpq.set(c_Key3(9, 1, -1), common_term_41);
    result_by_q.set(c_Key1(-1), common_term_41);
    // q = 0
    result_by_lpq.set(c_Key3(9, 1, 0), common_term_42);
    result_by_q.set(c_Key1(0), common_term_42);
    // q = 1
    result_by_lpq.set(c_Key3(9, 1, 1), common_term_43);
    result_by_q.set(c_Key1(1), common_term_43);
    // q = 2
    result_by_lpq.set(c_Key3(9, 1, 2), common_term_44);
    result_by_q.set(c_Key1(2), common_term_44);
    // q = 3
    result_by_lpq.set(c_Key3(9, 1, 3), common_term_45);
    result_by_q.set(c_Key1(3), common_term_45);
    // q = 4
    result_by_lpq.set(c_Key3(9, 1, 4), common_term_46);
    result_by_q.set(c_Key1(4), common_term_46);
    // q = 5
    result_by_lpq.set(c_Key3(9, 1, 5), common_term_47);
    result_by_q.set(c_Key1(5), common_term_47);
    // q = 6
    result_by_lpq.set(c_Key3(9, 1, 6), common_term_48);
    result_by_q.set(c_Key1(6), common_term_48);
    // q = 7
    result_by_lpq.set(c_Key3(9, 1, 7), common_term_49);
    result_by_q.set(c_Key1(7), common_term_49);
    // q = 8
    result_by_lpq.set(c_Key3(9, 1, 8), common_term_50);
    result_by_q.set(c_Key1(8), common_term_50);
    // q = 9
    result_by_lpq.set(c_Key3(9, 1, 9), common_term_51);
    result_by_q.set(c_Key1(9), common_term_51);
    // q = 10
    result_by_lpq.set(c_Key3(9, 1, 10), common_term_52);
    result_by_q.set(c_Key1(10), common_term_52);
    // q = 11
    result_by_lpq.set(c_Key3(9, 1, 11), common_term_53);
    result_by_q.set(c_Key1(11), common_term_53);
    // q = 12
    result_by_lpq.set(c_Key3(9, 1, 12), common_term_54);
    result_by_q.set(c_Key1(12), common_term_54);
    // q = 13
    result_by_lpq.set(c_Key3(9, 1, 13), common_term_55);
    result_by_q.set(c_Key1(13), common_term_55);
    // q = 14
    result_by_lpq.set(c_Key3(9, 1, 14), common_term_56);
    result_by_q.set(c_Key1(14), common_term_56);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 1), result_by_q);
    result_by_q.clear();

    // l , p = (9, 2).
    // q = -14
    result_by_lpq.set(c_Key3(9, 2, -14), common_term_57);
    result_by_q.set(c_Key1(-14), common_term_57);
    // q = -13
    result_by_lpq.set(c_Key3(9, 2, -13), common_term_58);
    result_by_q.set(c_Key1(-13), common_term_58);
    // q = -12
    result_by_lpq.set(c_Key3(9, 2, -12), common_term_59);
    result_by_q.set(c_Key1(-12), common_term_59);
    // q = -11
    result_by_lpq.set(c_Key3(9, 2, -11), common_term_60);
    result_by_q.set(c_Key1(-11), common_term_60);
    // q = -10
    result_by_lpq.set(c_Key3(9, 2, -10), common_term_61);
    result_by_q.set(c_Key1(-10), common_term_61);
    // q = -9
    result_by_lpq.set(c_Key3(9, 2, -9), common_term_62);
    result_by_q.set(c_Key1(-9), common_term_62);
    // q = -8
    result_by_lpq.set(c_Key3(9, 2, -8), common_term_63);
    result_by_q.set(c_Key1(-8), common_term_63);
    // q = -7
    result_by_lpq.set(c_Key3(9, 2, -7), common_term_64);
    result_by_q.set(c_Key1(-7), common_term_64);
    // q = -6
    result_by_lpq.set(c_Key3(9, 2, -6), common_term_65);
    result_by_q.set(c_Key1(-6), common_term_65);
    // q = -5
    result_by_lpq.set(c_Key3(9, 2, -5), common_term_66);
    result_by_q.set(c_Key1(-5), common_term_66);
    // q = -4
    result_by_lpq.set(c_Key3(9, 2, -4), common_term_67);
    result_by_q.set(c_Key1(-4), common_term_67);
    // q = -3
    result_by_lpq.set(c_Key3(9, 2, -3), common_term_68);
    result_by_q.set(c_Key1(-3), common_term_68);
    // q = -2
    result_by_lpq.set(c_Key3(9, 2, -2), common_term_69);
    result_by_q.set(c_Key1(-2), common_term_69);
    // q = -1
    result_by_lpq.set(c_Key3(9, 2, -1), common_term_70);
    result_by_q.set(c_Key1(-1), common_term_70);
    // q = 0
    result_by_lpq.set(c_Key3(9, 2, 0), common_term_71);
    result_by_q.set(c_Key1(0), common_term_71);
    // q = 1
    result_by_lpq.set(c_Key3(9, 2, 1), common_term_72);
    result_by_q.set(c_Key1(1), common_term_72);
    // q = 2
    result_by_lpq.set(c_Key3(9, 2, 2), common_term_73);
    result_by_q.set(c_Key1(2), common_term_73);
    // q = 3
    result_by_lpq.set(c_Key3(9, 2, 3), common_term_74);
    result_by_q.set(c_Key1(3), common_term_74);
    // q = 4
    result_by_lpq.set(c_Key3(9, 2, 4), common_term_75);
    result_by_q.set(c_Key1(4), common_term_75);
    // q = 5
    result_by_lpq.set(c_Key3(9, 2, 5), common_term_76);
    result_by_q.set(c_Key1(5), common_term_76);
    // q = 6
    result_by_lpq.set(c_Key3(9, 2, 6), common_term_77);
    result_by_q.set(c_Key1(6), common_term_77);
    // q = 7
    result_by_lpq.set(c_Key3(9, 2, 7), common_term_78);
    result_by_q.set(c_Key1(7), common_term_78);
    // q = 8
    result_by_lpq.set(c_Key3(9, 2, 8), common_term_79);
    result_by_q.set(c_Key1(8), common_term_79);
    // q = 9
    result_by_lpq.set(c_Key3(9, 2, 9), common_term_80);
    result_by_q.set(c_Key1(9), common_term_80);
    // q = 10
    result_by_lpq.set(c_Key3(9, 2, 10), common_term_81);
    result_by_q.set(c_Key1(10), common_term_81);
    // q = 11
    result_by_lpq.set(c_Key3(9, 2, 11), common_term_82);
    result_by_q.set(c_Key1(11), common_term_82);
    // q = 12
    result_by_lpq.set(c_Key3(9, 2, 12), common_term_83);
    result_by_q.set(c_Key1(12), common_term_83);
    // q = 13
    result_by_lpq.set(c_Key3(9, 2, 13), common_term_84);
    result_by_q.set(c_Key1(13), common_term_84);
    // q = 14
    result_by_lpq.set(c_Key3(9, 2, 14), common_term_85);
    result_by_q.set(c_Key1(14), common_term_85);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 2), result_by_q);
    result_by_q.clear();

    // l , p = (9, 3).
    // q = -14
    result_by_lpq.set(c_Key3(9, 3, -14), common_term_86);
    result_by_q.set(c_Key1(-14), common_term_86);
    // q = -13
    result_by_lpq.set(c_Key3(9, 3, -13), common_term_87);
    result_by_q.set(c_Key1(-13), common_term_87);
    // q = -12
    result_by_lpq.set(c_Key3(9, 3, -12), common_term_88);
    result_by_q.set(c_Key1(-12), common_term_88);
    // q = -11
    result_by_lpq.set(c_Key3(9, 3, -11), common_term_89);
    result_by_q.set(c_Key1(-11), common_term_89);
    // q = -10
    result_by_lpq.set(c_Key3(9, 3, -10), common_term_90);
    result_by_q.set(c_Key1(-10), common_term_90);
    // q = -9
    result_by_lpq.set(c_Key3(9, 3, -9), common_term_91);
    result_by_q.set(c_Key1(-9), common_term_91);
    // q = -8
    result_by_lpq.set(c_Key3(9, 3, -8), common_term_92);
    result_by_q.set(c_Key1(-8), common_term_92);
    // q = -7
    result_by_lpq.set(c_Key3(9, 3, -7), common_term_93);
    result_by_q.set(c_Key1(-7), common_term_93);
    // q = -6
    result_by_lpq.set(c_Key3(9, 3, -6), common_term_94);
    result_by_q.set(c_Key1(-6), common_term_94);
    // q = -5
    result_by_lpq.set(c_Key3(9, 3, -5), common_term_95);
    result_by_q.set(c_Key1(-5), common_term_95);
    // q = -4
    result_by_lpq.set(c_Key3(9, 3, -4), common_term_96);
    result_by_q.set(c_Key1(-4), common_term_96);
    // q = -3
    result_by_lpq.set(c_Key3(9, 3, -3), common_term_97);
    result_by_q.set(c_Key1(-3), common_term_97);
    // q = -2
    result_by_lpq.set(c_Key3(9, 3, -2), common_term_98);
    result_by_q.set(c_Key1(-2), common_term_98);
    // q = -1
    result_by_lpq.set(c_Key3(9, 3, -1), common_term_99);
    result_by_q.set(c_Key1(-1), common_term_99);
    // q = 0
    result_by_lpq.set(c_Key3(9, 3, 0), common_term_100);
    result_by_q.set(c_Key1(0), common_term_100);
    // q = 1
    result_by_lpq.set(c_Key3(9, 3, 1), common_term_101);
    result_by_q.set(c_Key1(1), common_term_101);
    // q = 2
    result_by_lpq.set(c_Key3(9, 3, 2), common_term_102);
    result_by_q.set(c_Key1(2), common_term_102);
    // q = 3
    result_by_lpq.set(c_Key3(9, 3, 3), common_term_103);
    result_by_q.set(c_Key1(3), common_term_103);
    // q = 4
    result_by_lpq.set(c_Key3(9, 3, 4), common_term_104);
    result_by_q.set(c_Key1(4), common_term_104);
    // q = 5
    result_by_lpq.set(c_Key3(9, 3, 5), common_term_105);
    result_by_q.set(c_Key1(5), common_term_105);
    // q = 6
    result_by_lpq.set(c_Key3(9, 3, 6), common_term_106);
    result_by_q.set(c_Key1(6), common_term_106);
    // q = 7
    result_by_lpq.set(c_Key3(9, 3, 7), common_term_107);
    result_by_q.set(c_Key1(7), common_term_107);
    // q = 8
    result_by_lpq.set(c_Key3(9, 3, 8), common_term_108);
    result_by_q.set(c_Key1(8), common_term_108);
    // q = 9
    result_by_lpq.set(c_Key3(9, 3, 9), common_term_109);
    result_by_q.set(c_Key1(9), common_term_109);
    // q = 10
    result_by_lpq.set(c_Key3(9, 3, 10), common_term_110);
    result_by_q.set(c_Key1(10), common_term_110);
    // q = 11
    result_by_lpq.set(c_Key3(9, 3, 11), common_term_111);
    result_by_q.set(c_Key1(11), common_term_111);
    // q = 12
    result_by_lpq.set(c_Key3(9, 3, 12), common_term_112);
    result_by_q.set(c_Key1(12), common_term_112);
    // q = 13
    result_by_lpq.set(c_Key3(9, 3, 13), common_term_113);
    result_by_q.set(c_Key1(13), common_term_113);
    // q = 14
    result_by_lpq.set(c_Key3(9, 3, 14), common_term_114);
    result_by_q.set(c_Key1(14), common_term_114);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 3), result_by_q);
    result_by_q.clear();

    // l , p = (9, 4).
    // q = -14
    result_by_lpq.set(c_Key3(9, 4, -14), common_term_115);
    result_by_q.set(c_Key1(-14), common_term_115);
    // q = -13
    result_by_lpq.set(c_Key3(9, 4, -13), common_term_116);
    result_by_q.set(c_Key1(-13), common_term_116);
    // q = -12
    result_by_lpq.set(c_Key3(9, 4, -12), common_term_117);
    result_by_q.set(c_Key1(-12), common_term_117);
    // q = -11
    result_by_lpq.set(c_Key3(9, 4, -11), common_term_118);
    result_by_q.set(c_Key1(-11), common_term_118);
    // q = -10
    result_by_lpq.set(c_Key3(9, 4, -10), common_term_119);
    result_by_q.set(c_Key1(-10), common_term_119);
    // q = -9
    result_by_lpq.set(c_Key3(9, 4, -9), common_term_120);
    result_by_q.set(c_Key1(-9), common_term_120);
    // q = -8
    result_by_lpq.set(c_Key3(9, 4, -8), common_term_121);
    result_by_q.set(c_Key1(-8), common_term_121);
    // q = -7
    result_by_lpq.set(c_Key3(9, 4, -7), common_term_122);
    result_by_q.set(c_Key1(-7), common_term_122);
    // q = -6
    result_by_lpq.set(c_Key3(9, 4, -6), common_term_123);
    result_by_q.set(c_Key1(-6), common_term_123);
    // q = -5
    result_by_lpq.set(c_Key3(9, 4, -5), common_term_124);
    result_by_q.set(c_Key1(-5), common_term_124);
    // q = -4
    result_by_lpq.set(c_Key3(9, 4, -4), common_term_125);
    result_by_q.set(c_Key1(-4), common_term_125);
    // q = -3
    result_by_lpq.set(c_Key3(9, 4, -3), common_term_126);
    result_by_q.set(c_Key1(-3), common_term_126);
    // q = -2
    result_by_lpq.set(c_Key3(9, 4, -2), common_term_127);
    result_by_q.set(c_Key1(-2), common_term_127);
    // q = -1
    result_by_lpq.set(c_Key3(9, 4, -1), common_term_128);
    result_by_q.set(c_Key1(-1), common_term_128);
    // q = 0
    result_by_lpq.set(c_Key3(9, 4, 0), common_term_129);
    result_by_q.set(c_Key1(0), common_term_129);
    // q = 1
    result_by_lpq.set(c_Key3(9, 4, 1), common_term_130);
    result_by_q.set(c_Key1(1), common_term_130);
    // q = 2
    result_by_lpq.set(c_Key3(9, 4, 2), common_term_131);
    result_by_q.set(c_Key1(2), common_term_131);
    // q = 3
    result_by_lpq.set(c_Key3(9, 4, 3), common_term_132);
    result_by_q.set(c_Key1(3), common_term_132);
    // q = 4
    result_by_lpq.set(c_Key3(9, 4, 4), common_term_133);
    result_by_q.set(c_Key1(4), common_term_133);
    // q = 5
    result_by_lpq.set(c_Key3(9, 4, 5), common_term_134);
    result_by_q.set(c_Key1(5), common_term_134);
    // q = 6
    result_by_lpq.set(c_Key3(9, 4, 6), common_term_135);
    result_by_q.set(c_Key1(6), common_term_135);
    // q = 7
    result_by_lpq.set(c_Key3(9, 4, 7), common_term_136);
    result_by_q.set(c_Key1(7), common_term_136);
    // q = 8
    result_by_lpq.set(c_Key3(9, 4, 8), common_term_137);
    result_by_q.set(c_Key1(8), common_term_137);
    // q = 9
    result_by_lpq.set(c_Key3(9, 4, 9), common_term_138);
    result_by_q.set(c_Key1(9), common_term_138);
    // q = 10
    result_by_lpq.set(c_Key3(9, 4, 10), common_term_139);
    result_by_q.set(c_Key1(10), common_term_139);
    // q = 11
    result_by_lpq.set(c_Key3(9, 4, 11), common_term_140);
    result_by_q.set(c_Key1(11), common_term_140);
    // q = 12
    result_by_lpq.set(c_Key3(9, 4, 12), common_term_141);
    result_by_q.set(c_Key1(12), common_term_141);
    // q = 13
    result_by_lpq.set(c_Key3(9, 4, 13), common_term_142);
    result_by_q.set(c_Key1(13), common_term_142);
    // q = 14
    result_by_lpq.set(c_Key3(9, 4, 14), common_term_143);
    result_by_q.set(c_Key1(14), common_term_143);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 4), result_by_q);
    result_by_q.clear();

    // l , p = (9, 5).
    // q = -14
    result_by_lpq.set(c_Key3(9, 5, -14), common_term_143);
    result_by_q.set(c_Key1(-14), common_term_143);
    // q = -13
    result_by_lpq.set(c_Key3(9, 5, -13), common_term_142);
    result_by_q.set(c_Key1(-13), common_term_142);
    // q = -12
    result_by_lpq.set(c_Key3(9, 5, -12), common_term_141);
    result_by_q.set(c_Key1(-12), common_term_141);
    // q = -11
    result_by_lpq.set(c_Key3(9, 5, -11), common_term_140);
    result_by_q.set(c_Key1(-11), common_term_140);
    // q = -10
    result_by_lpq.set(c_Key3(9, 5, -10), common_term_139);
    result_by_q.set(c_Key1(-10), common_term_139);
    // q = -9
    result_by_lpq.set(c_Key3(9, 5, -9), common_term_138);
    result_by_q.set(c_Key1(-9), common_term_138);
    // q = -8
    result_by_lpq.set(c_Key3(9, 5, -8), common_term_137);
    result_by_q.set(c_Key1(-8), common_term_137);
    // q = -7
    result_by_lpq.set(c_Key3(9, 5, -7), common_term_136);
    result_by_q.set(c_Key1(-7), common_term_136);
    // q = -6
    result_by_lpq.set(c_Key3(9, 5, -6), common_term_135);
    result_by_q.set(c_Key1(-6), common_term_135);
    // q = -5
    result_by_lpq.set(c_Key3(9, 5, -5), common_term_134);
    result_by_q.set(c_Key1(-5), common_term_134);
    // q = -4
    result_by_lpq.set(c_Key3(9, 5, -4), common_term_133);
    result_by_q.set(c_Key1(-4), common_term_133);
    // q = -3
    result_by_lpq.set(c_Key3(9, 5, -3), common_term_132);
    result_by_q.set(c_Key1(-3), common_term_132);
    // q = -2
    result_by_lpq.set(c_Key3(9, 5, -2), common_term_131);
    result_by_q.set(c_Key1(-2), common_term_131);
    // q = -1
    result_by_lpq.set(c_Key3(9, 5, -1), common_term_130);
    result_by_q.set(c_Key1(-1), common_term_130);
    // q = 0
    result_by_lpq.set(c_Key3(9, 5, 0), common_term_129);
    result_by_q.set(c_Key1(0), common_term_129);
    // q = 1
    result_by_lpq.set(c_Key3(9, 5, 1), common_term_128);
    result_by_q.set(c_Key1(1), common_term_128);
    // q = 2
    result_by_lpq.set(c_Key3(9, 5, 2), common_term_127);
    result_by_q.set(c_Key1(2), common_term_127);
    // q = 3
    result_by_lpq.set(c_Key3(9, 5, 3), common_term_126);
    result_by_q.set(c_Key1(3), common_term_126);
    // q = 4
    result_by_lpq.set(c_Key3(9, 5, 4), common_term_125);
    result_by_q.set(c_Key1(4), common_term_125);
    // q = 5
    result_by_lpq.set(c_Key3(9, 5, 5), common_term_124);
    result_by_q.set(c_Key1(5), common_term_124);
    // q = 6
    result_by_lpq.set(c_Key3(9, 5, 6), common_term_123);
    result_by_q.set(c_Key1(6), common_term_123);
    // q = 7
    result_by_lpq.set(c_Key3(9, 5, 7), common_term_122);
    result_by_q.set(c_Key1(7), common_term_122);
    // q = 8
    result_by_lpq.set(c_Key3(9, 5, 8), common_term_121);
    result_by_q.set(c_Key1(8), common_term_121);
    // q = 9
    result_by_lpq.set(c_Key3(9, 5, 9), common_term_120);
    result_by_q.set(c_Key1(9), common_term_120);
    // q = 10
    result_by_lpq.set(c_Key3(9, 5, 10), common_term_119);
    result_by_q.set(c_Key1(10), common_term_119);
    // q = 11
    result_by_lpq.set(c_Key3(9, 5, 11), common_term_118);
    result_by_q.set(c_Key1(11), common_term_118);
    // q = 12
    result_by_lpq.set(c_Key3(9, 5, 12), common_term_117);
    result_by_q.set(c_Key1(12), common_term_117);
    // q = 13
    result_by_lpq.set(c_Key3(9, 5, 13), common_term_116);
    result_by_q.set(c_Key1(13), common_term_116);
    // q = 14
    result_by_lpq.set(c_Key3(9, 5, 14), common_term_115);
    result_by_q.set(c_Key1(14), common_term_115);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 5), result_by_q);
    result_by_q.clear();

    // l , p = (9, 6).
    // q = -14
    result_by_lpq.set(c_Key3(9, 6, -14), common_term_114);
    result_by_q.set(c_Key1(-14), common_term_114);
    // q = -13
    result_by_lpq.set(c_Key3(9, 6, -13), common_term_113);
    result_by_q.set(c_Key1(-13), common_term_113);
    // q = -12
    result_by_lpq.set(c_Key3(9, 6, -12), common_term_112);
    result_by_q.set(c_Key1(-12), common_term_112);
    // q = -11
    result_by_lpq.set(c_Key3(9, 6, -11), common_term_111);
    result_by_q.set(c_Key1(-11), common_term_111);
    // q = -10
    result_by_lpq.set(c_Key3(9, 6, -10), common_term_110);
    result_by_q.set(c_Key1(-10), common_term_110);
    // q = -9
    result_by_lpq.set(c_Key3(9, 6, -9), common_term_109);
    result_by_q.set(c_Key1(-9), common_term_109);
    // q = -8
    result_by_lpq.set(c_Key3(9, 6, -8), common_term_108);
    result_by_q.set(c_Key1(-8), common_term_108);
    // q = -7
    result_by_lpq.set(c_Key3(9, 6, -7), common_term_107);
    result_by_q.set(c_Key1(-7), common_term_107);
    // q = -6
    result_by_lpq.set(c_Key3(9, 6, -6), common_term_106);
    result_by_q.set(c_Key1(-6), common_term_106);
    // q = -5
    result_by_lpq.set(c_Key3(9, 6, -5), common_term_105);
    result_by_q.set(c_Key1(-5), common_term_105);
    // q = -4
    result_by_lpq.set(c_Key3(9, 6, -4), common_term_104);
    result_by_q.set(c_Key1(-4), common_term_104);
    // q = -3
    result_by_lpq.set(c_Key3(9, 6, -3), common_term_103);
    result_by_q.set(c_Key1(-3), common_term_103);
    // q = -2
    result_by_lpq.set(c_Key3(9, 6, -2), common_term_102);
    result_by_q.set(c_Key1(-2), common_term_102);
    // q = -1
    result_by_lpq.set(c_Key3(9, 6, -1), common_term_101);
    result_by_q.set(c_Key1(-1), common_term_101);
    // q = 0
    result_by_lpq.set(c_Key3(9, 6, 0), common_term_100);
    result_by_q.set(c_Key1(0), common_term_100);
    // q = 1
    result_by_lpq.set(c_Key3(9, 6, 1), common_term_99);
    result_by_q.set(c_Key1(1), common_term_99);
    // q = 2
    result_by_lpq.set(c_Key3(9, 6, 2), common_term_98);
    result_by_q.set(c_Key1(2), common_term_98);
    // q = 3
    result_by_lpq.set(c_Key3(9, 6, 3), common_term_97);
    result_by_q.set(c_Key1(3), common_term_97);
    // q = 4
    result_by_lpq.set(c_Key3(9, 6, 4), common_term_96);
    result_by_q.set(c_Key1(4), common_term_96);
    // q = 5
    result_by_lpq.set(c_Key3(9, 6, 5), common_term_95);
    result_by_q.set(c_Key1(5), common_term_95);
    // q = 6
    result_by_lpq.set(c_Key3(9, 6, 6), common_term_94);
    result_by_q.set(c_Key1(6), common_term_94);
    // q = 7
    result_by_lpq.set(c_Key3(9, 6, 7), common_term_93);
    result_by_q.set(c_Key1(7), common_term_93);
    // q = 8
    result_by_lpq.set(c_Key3(9, 6, 8), common_term_92);
    result_by_q.set(c_Key1(8), common_term_92);
    // q = 9
    result_by_lpq.set(c_Key3(9, 6, 9), common_term_91);
    result_by_q.set(c_Key1(9), common_term_91);
    // q = 10
    result_by_lpq.set(c_Key3(9, 6, 10), common_term_90);
    result_by_q.set(c_Key1(10), common_term_90);
    // q = 11
    result_by_lpq.set(c_Key3(9, 6, 11), common_term_89);
    result_by_q.set(c_Key1(11), common_term_89);
    // q = 12
    result_by_lpq.set(c_Key3(9, 6, 12), common_term_88);
    result_by_q.set(c_Key1(12), common_term_88);
    // q = 13
    result_by_lpq.set(c_Key3(9, 6, 13), common_term_87);
    result_by_q.set(c_Key1(13), common_term_87);
    // q = 14
    result_by_lpq.set(c_Key3(9, 6, 14), common_term_86);
    result_by_q.set(c_Key1(14), common_term_86);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 6), result_by_q);
    result_by_q.clear();

    // l , p = (9, 7).
    // q = -14
    result_by_lpq.set(c_Key3(9, 7, -14), common_term_85);
    result_by_q.set(c_Key1(-14), common_term_85);
    // q = -13
    result_by_lpq.set(c_Key3(9, 7, -13), common_term_84);
    result_by_q.set(c_Key1(-13), common_term_84);
    // q = -12
    result_by_lpq.set(c_Key3(9, 7, -12), common_term_83);
    result_by_q.set(c_Key1(-12), common_term_83);
    // q = -11
    result_by_lpq.set(c_Key3(9, 7, -11), common_term_82);
    result_by_q.set(c_Key1(-11), common_term_82);
    // q = -10
    result_by_lpq.set(c_Key3(9, 7, -10), common_term_81);
    result_by_q.set(c_Key1(-10), common_term_81);
    // q = -9
    result_by_lpq.set(c_Key3(9, 7, -9), common_term_80);
    result_by_q.set(c_Key1(-9), common_term_80);
    // q = -8
    result_by_lpq.set(c_Key3(9, 7, -8), common_term_79);
    result_by_q.set(c_Key1(-8), common_term_79);
    // q = -7
    result_by_lpq.set(c_Key3(9, 7, -7), common_term_78);
    result_by_q.set(c_Key1(-7), common_term_78);
    // q = -6
    result_by_lpq.set(c_Key3(9, 7, -6), common_term_77);
    result_by_q.set(c_Key1(-6), common_term_77);
    // q = -5
    result_by_lpq.set(c_Key3(9, 7, -5), common_term_76);
    result_by_q.set(c_Key1(-5), common_term_76);
    // q = -4
    result_by_lpq.set(c_Key3(9, 7, -4), common_term_75);
    result_by_q.set(c_Key1(-4), common_term_75);
    // q = -3
    result_by_lpq.set(c_Key3(9, 7, -3), common_term_74);
    result_by_q.set(c_Key1(-3), common_term_74);
    // q = -2
    result_by_lpq.set(c_Key3(9, 7, -2), common_term_73);
    result_by_q.set(c_Key1(-2), common_term_73);
    // q = -1
    result_by_lpq.set(c_Key3(9, 7, -1), common_term_72);
    result_by_q.set(c_Key1(-1), common_term_72);
    // q = 0
    result_by_lpq.set(c_Key3(9, 7, 0), common_term_71);
    result_by_q.set(c_Key1(0), common_term_71);
    // q = 1
    result_by_lpq.set(c_Key3(9, 7, 1), common_term_70);
    result_by_q.set(c_Key1(1), common_term_70);
    // q = 2
    result_by_lpq.set(c_Key3(9, 7, 2), common_term_69);
    result_by_q.set(c_Key1(2), common_term_69);
    // q = 3
    result_by_lpq.set(c_Key3(9, 7, 3), common_term_68);
    result_by_q.set(c_Key1(3), common_term_68);
    // q = 4
    result_by_lpq.set(c_Key3(9, 7, 4), common_term_67);
    result_by_q.set(c_Key1(4), common_term_67);
    // q = 5
    result_by_lpq.set(c_Key3(9, 7, 5), common_term_66);
    result_by_q.set(c_Key1(5), common_term_66);
    // q = 6
    result_by_lpq.set(c_Key3(9, 7, 6), common_term_65);
    result_by_q.set(c_Key1(6), common_term_65);
    // q = 7
    result_by_lpq.set(c_Key3(9, 7, 7), common_term_64);
    result_by_q.set(c_Key1(7), common_term_64);
    // q = 8
    result_by_lpq.set(c_Key3(9, 7, 8), common_term_63);
    result_by_q.set(c_Key1(8), common_term_63);
    // q = 9
    result_by_lpq.set(c_Key3(9, 7, 9), common_term_62);
    result_by_q.set(c_Key1(9), common_term_62);
    // q = 10
    result_by_lpq.set(c_Key3(9, 7, 10), common_term_61);
    result_by_q.set(c_Key1(10), common_term_61);
    // q = 11
    result_by_lpq.set(c_Key3(9, 7, 11), common_term_60);
    result_by_q.set(c_Key1(11), common_term_60);
    // q = 12
    result_by_lpq.set(c_Key3(9, 7, 12), common_term_59);
    result_by_q.set(c_Key1(12), common_term_59);
    // q = 13
    result_by_lpq.set(c_Key3(9, 7, 13), common_term_58);
    result_by_q.set(c_Key1(13), common_term_58);
    // q = 14
    result_by_lpq.set(c_Key3(9, 7, 14), common_term_57);
    result_by_q.set(c_Key1(14), common_term_57);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 7), result_by_q);
    result_by_q.clear();

    // l , p = (9, 8).
    // q = -14
    result_by_lpq.set(c_Key3(9, 8, -14), common_term_56);
    result_by_q.set(c_Key1(-14), common_term_56);
    // q = -13
    result_by_lpq.set(c_Key3(9, 8, -13), common_term_55);
    result_by_q.set(c_Key1(-13), common_term_55);
    // q = -12
    result_by_lpq.set(c_Key3(9, 8, -12), common_term_54);
    result_by_q.set(c_Key1(-12), common_term_54);
    // q = -11
    result_by_lpq.set(c_Key3(9, 8, -11), common_term_53);
    result_by_q.set(c_Key1(-11), common_term_53);
    // q = -10
    result_by_lpq.set(c_Key3(9, 8, -10), common_term_52);
    result_by_q.set(c_Key1(-10), common_term_52);
    // q = -9
    result_by_lpq.set(c_Key3(9, 8, -9), common_term_51);
    result_by_q.set(c_Key1(-9), common_term_51);
    // q = -8
    result_by_lpq.set(c_Key3(9, 8, -8), common_term_50);
    result_by_q.set(c_Key1(-8), common_term_50);
    // q = -7
    result_by_lpq.set(c_Key3(9, 8, -7), common_term_49);
    result_by_q.set(c_Key1(-7), common_term_49);
    // q = -6
    result_by_lpq.set(c_Key3(9, 8, -6), common_term_48);
    result_by_q.set(c_Key1(-6), common_term_48);
    // q = -5
    result_by_lpq.set(c_Key3(9, 8, -5), common_term_47);
    result_by_q.set(c_Key1(-5), common_term_47);
    // q = -4
    result_by_lpq.set(c_Key3(9, 8, -4), common_term_46);
    result_by_q.set(c_Key1(-4), common_term_46);
    // q = -3
    result_by_lpq.set(c_Key3(9, 8, -3), common_term_45);
    result_by_q.set(c_Key1(-3), common_term_45);
    // q = -2
    result_by_lpq.set(c_Key3(9, 8, -2), common_term_44);
    result_by_q.set(c_Key1(-2), common_term_44);
    // q = -1
    result_by_lpq.set(c_Key3(9, 8, -1), common_term_43);
    result_by_q.set(c_Key1(-1), common_term_43);
    // q = 0
    result_by_lpq.set(c_Key3(9, 8, 0), common_term_42);
    result_by_q.set(c_Key1(0), common_term_42);
    // q = 1
    result_by_lpq.set(c_Key3(9, 8, 1), common_term_41);
    result_by_q.set(c_Key1(1), common_term_41);
    // q = 2
    result_by_lpq.set(c_Key3(9, 8, 2), common_term_40);
    result_by_q.set(c_Key1(2), common_term_40);
    // q = 3
    result_by_lpq.set(c_Key3(9, 8, 3), common_term_39);
    result_by_q.set(c_Key1(3), common_term_39);
    // q = 4
    result_by_lpq.set(c_Key3(9, 8, 4), common_term_38);
    result_by_q.set(c_Key1(4), common_term_38);
    // q = 5
    result_by_lpq.set(c_Key3(9, 8, 5), common_term_37);
    result_by_q.set(c_Key1(5), common_term_37);
    // q = 6
    result_by_lpq.set(c_Key3(9, 8, 6), common_term_36);
    result_by_q.set(c_Key1(6), common_term_36);
    // q = 7
    result_by_lpq.set(c_Key3(9, 8, 7), common_term_35);
    result_by_q.set(c_Key1(7), common_term_35);
    // q = 8
    result_by_lpq.set(c_Key3(9, 8, 8), common_term_34);
    result_by_q.set(c_Key1(8), common_term_34);
    // q = 9
    result_by_lpq.set(c_Key3(9, 8, 9), common_term_33);
    result_by_q.set(c_Key1(9), common_term_33);
    // q = 10
    result_by_lpq.set(c_Key3(9, 8, 10), common_term_32);
    result_by_q.set(c_Key1(10), common_term_32);
    // q = 11
    result_by_lpq.set(c_Key3(9, 8, 11), common_term_31);
    result_by_q.set(c_Key1(11), common_term_31);
    // q = 12
    result_by_lpq.set(c_Key3(9, 8, 12), common_term_30);
    result_by_q.set(c_Key1(12), common_term_30);
    // q = 13
    result_by_lpq.set(c_Key3(9, 8, 13), common_term_29);
    result_by_q.set(c_Key1(13), common_term_29);
    // q = 14
    result_by_lpq.set(c_Key3(9, 8, 14), common_term_28);
    result_by_q.set(c_Key1(14), common_term_28);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 8), result_by_q);
    result_by_q.clear();

    // l , p = (9, 9).
    // q = -14
    result_by_lpq.set(c_Key3(9, 9, -14), common_term_27);
    result_by_q.set(c_Key1(-14), common_term_27);
    // q = -13
    result_by_lpq.set(c_Key3(9, 9, -13), common_term_26);
    result_by_q.set(c_Key1(-13), common_term_26);
    // q = -12
    result_by_lpq.set(c_Key3(9, 9, -12), common_term_25);
    result_by_q.set(c_Key1(-12), common_term_25);
    // q = -11
    result_by_lpq.set(c_Key3(9, 9, -11), common_term_24);
    result_by_q.set(c_Key1(-11), common_term_24);
    // q = -10
    result_by_lpq.set(c_Key3(9, 9, -10), common_term_23);
    result_by_q.set(c_Key1(-10), common_term_23);
    // q = -9
    result_by_lpq.set(c_Key3(9, 9, -9), common_term_22);
    result_by_q.set(c_Key1(-9), common_term_22);
    // q = -8
    result_by_lpq.set(c_Key3(9, 9, -8), common_term_21);
    result_by_q.set(c_Key1(-8), common_term_21);
    // q = -7
    result_by_lpq.set(c_Key3(9, 9, -7), common_term_20);
    result_by_q.set(c_Key1(-7), common_term_20);
    // q = -6
    result_by_lpq.set(c_Key3(9, 9, -6), common_term_19);
    result_by_q.set(c_Key1(-6), common_term_19);
    // q = -5
    result_by_lpq.set(c_Key3(9, 9, -5), common_term_18);
    result_by_q.set(c_Key1(-5), common_term_18);
    // q = -4
    result_by_lpq.set(c_Key3(9, 9, -4), common_term_17);
    result_by_q.set(c_Key1(-4), common_term_17);
    // q = -3
    result_by_lpq.set(c_Key3(9, 9, -3), common_term_16);
    result_by_q.set(c_Key1(-3), common_term_16);
    // q = -2
    result_by_lpq.set(c_Key3(9, 9, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(9, 9, -1), common_term_14);
    result_by_q.set(c_Key1(-1), common_term_14);
    // q = 0
    result_by_lpq.set(c_Key3(9, 9, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(9, 9, 1), common_term_12);
    result_by_q.set(c_Key1(1), common_term_12);
    // q = 2
    result_by_lpq.set(c_Key3(9, 9, 2), common_term_11);
    result_by_q.set(c_Key1(2), common_term_11);
    // q = 3
    result_by_lpq.set(c_Key3(9, 9, 3), common_term_10);
    result_by_q.set(c_Key1(3), common_term_10);
    // q = 4
    result_by_lpq.set(c_Key3(9, 9, 4), common_term_9);
    result_by_q.set(c_Key1(4), common_term_9);
    // q = 5
    result_by_lpq.set(c_Key3(9, 9, 5), common_term_8);
    result_by_q.set(c_Key1(5), common_term_8);
    // q = 6
    result_by_lpq.set(c_Key3(9, 9, 6), common_term_7);
    result_by_q.set(c_Key1(6), common_term_7);
    // q = 7
    result_by_lpq.set(c_Key3(9, 9, 7), common_term_6);
    result_by_q.set(c_Key1(7), common_term_6);
    // q = 8
    result_by_lpq.set(c_Key3(9, 9, 8), common_term_5);
    result_by_q.set(c_Key1(8), common_term_5);
    // q = 10
    result_by_lpq.set(c_Key3(9, 9, 10), common_term_4);
    result_by_q.set(c_Key1(10), common_term_4);
    // q = 11
    result_by_lpq.set(c_Key3(9, 9, 11), common_term_3);
    result_by_q.set(c_Key1(11), common_term_3);
    // q = 12
    result_by_lpq.set(c_Key3(9, 9, 12), common_term_2);
    result_by_q.set(c_Key1(12), common_term_2);
    // q = 13
    result_by_lpq.set(c_Key3(9, 9, 13), common_term_1);
    result_by_q.set(c_Key1(13), common_term_1);
    // q = 14
    result_by_lpq.set(c_Key3(9, 9, 14), common_term_0);
    result_by_q.set(c_Key1(14), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 9), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l9_e20(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 9.
    // Functions are taylor expanded and truncated at eccentricity^20.

    // Functions Calculated for l = 9.

    c_IntMap<c_Key3, double> result_by_lpq(388);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(10);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
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
    double eccentricity_5 = eccentricity_4 * eccentricity;
    double eccentricity_10 = eccentricity_5 * eccentricity_5;
    double eccentricity_11 = eccentricity_10 * eccentricity;
    double common_term_0 = 0.00015679617398499164*eccentricity_19;
    double common_term_1 = 8.9430205290920138e-5*eccentricity_18;
    double common_term_2 = 0.0001623432391492702*eccentricity_19 + 4.8300467846063862e-5*eccentricity_17;
    double common_term_3 = 8.1976248392506114e-5*eccentricity_18 + 2.4236456046480068e-5*eccentricity_16;
    double common_term_4 = 7.6789738171820302e-5*eccentricity_19 + 3.7033348905558727e-5*eccentricity_17 + 1.0972844120165549e-5*eccentricity_15;
    double common_term_5 = 2.9041527230753514e-5*eccentricity_18 + 1.4243952430449084e-5*eccentricity_16 + 4.2731857291347253e-6*eccentricity_14;
    double common_term_6 = 1.3376329072978103e-5*eccentricity_19 + 8.4931725143365355e-6*eccentricity_17 + 4.2755598311153867e-6*eccentricity_15 + 1.3155568711124267e-6*eccentricity_13;
    double common_term_7 = 2.4585063330402758e-6*eccentricity_18 + 1.613628924684074e-6*eccentricity_16 + 8.4385886178030953e-7*eccentricity_14 + 2.708682766208401e-7*eccentricity_12;
    double common_term_8 = 2.4846154334552703e-7*eccentricity_19 + 1.9278834331331024e-7*eccentricity_17 + 1.323265212154101e-7*eccentricity_15 + 7.3068649457538346e-8*eccentricity_13 + 2.5052108385441719e-8*eccentricity_11;
    double common_term_9 = 2.0037804926928512e-9*eccentricity_18 + 1.6312743534479736e-9*eccentricity_16 + 1.1892896200333673e-9*eccentricity_14 + 7.0948353825957993e-10*eccentricity_12 + 2.6911444554673721e-10*eccentricity_10;
    double common_term_10 = 2.6165016474752718e-7*eccentricity_18 + 2.5509397608251903e-7*eccentricity_16 + 2.4010869936942799e-7*eccentricity_14 + 2.1256677267622906e-7*eccentricity_12 + 1.6685095623897707e-7*eccentricity_10 + 9.6881200396825397e-8*eccentricity_8;
    double common_term_11 = -0.00016081069640700235*eccentricity_19 - 0.00017162297404991271*eccentricity_17 - 0.00018252757402062958*eccentricity_15 - 0.00019257973251028807*eccentricity_13 - 0.00019979056437389771*eccentricity_11 - 0.00019841269841269841*eccentricity_9 - 0.00019841269841269841*eccentricity_7;
    double common_term_12 = 0.00031943913803181868*eccentricity_18 + 0.00050155324440497857*eccentricity_16 + 0.00076860002108982631*eccentricity_14 + 0.0011388506208147321*eccentricity_12 + 0.0021452767508370536*eccentricity_10 + 0.0158203125*eccentricity_6;
    double common_term_13 = 0.0066525194809303702*eccentricity_19 + 0.0082167697700105058*eccentricity_17 + 0.010607657260435038*eccentricity_15 + 0.010765726043503821*eccentricity_13 + 0.040013227513227513*eccentricity_11 - 0.093650793650793651*eccentricity_9 + 0.37777777777777778*eccentricity_7 - 0.26666666666666667*eccentricity_5;
    double common_term_14 = 0.038465812085422115*eccentricity_18 + 0.059673072711471593*eccentricity_16 - 0.035442644045816878*eccentricity_14 + 0.66459573135174141*eccentricity_12 - 2.2781210601645172*eccentricity_10 + 5.7898627387152778*eccentricity_8 - 5.6966145833333333*eccentricity_6 + 1.6276041666666667*eccentricity_4;
    double common_term_15 = 0.093066501256453876*eccentricity_19 + 0.44220154555961318*eccentricity_17 - 1.7061077008928623*eccentricity_15 + 9.3742131696428571*eccentricity_13 - 29.959151785714286*eccentricity_11 + 62.8171875*eccentricity_9 - 66.825*eccentricity_7 + 30.375*eccentricity_5 - 4.5*eccentricity_3;
    double common_term_16 = 5.0328579683877018*eccentricity_18 - 24.321003408507072*eccentricity_16 + 99.655787984368295*eccentricity_14 - 279.30555665192781*eccentricity_12 + 519.05630227548105*eccentricity_10 - 557.04748263888889*eccentricity_8 + 306.17024739583333*eccentricity_6 - 75.541666666666667*eccentricity_4 + 6.125*eccentricity_2;
    double common_term_17 = 54.344427907215654*eccentricity_19 - 239.98676288155768*eccentricity_17 + 826.85053773433488*eccentricity_15 - 2049.2777143959439*eccentricity_13 + 3470.4545601851852*eccentricity_11 - 3664.5597222222222*eccentricity_9 + 2203.2152777777778*eccentricity_7 - 686.83333333333333*eccentricity_5 + 95.0*eccentricity_3 - 4.0*eccentricity;
    double common_term_18 = -1868.5133856163701*eccentricity_18 + 5631.221562745031*eccentricity_16 - 12537.720637434743*eccentricity_14 + 19631.782016716005*eccentricity_12 - 20179.836617431641*eccentricity_10 + 12656.252746582031*eccentricity_8 - 4483.0078125*eccentricity_6 + 801.984375*eccentricity_4 - 58.5*eccentricity_2 + 1.0;
    double common_term_19 = -12163.062162179059*eccentricity_19 + 32753.830578230219*eccentricity_17 - 66414.233764433139*eccentricity_15 + 97048.000044780642*eccentricity_13 - 96633.478877314817*eccentricity_11 + 61626.019965277778*eccentricity_9 - 23508.472222222222*eccentricity_7 + 4882.2916666666667*eccentricity_5 - 467.5*eccentricity_3 + 14.0*eccentricity;
    double common_term_20 = 167543.69661696567*eccentricity_18 - 312770.7815736186*eccentricity_16 + 429376.97441965474*eccentricity_14 - 413227.00778149963*eccentricity_12 + 263997.47338092946*eccentricity_10 - 105162.43107638889*eccentricity_8 + 23978.449544270833*eccentricity_6 - 2702.7916666666667*eccentricity_4 + 107.375*eccentricity_2;
    double common_term_21 = 770334.25120671839*eccentricity_19 - 1335383.1991930007*eccentricity_17 + 1731218.8928236607*eccentricity_15 - 1608932.925703125*eccentricity_13 + 1020341.8275669643*eccentricity_11 - 416358.6*eccentricity_9 + 100843.65*eccentricity_7 - 12640.5*eccentricity_5 + 597.0*eccentricity_3;
    double common_term_22 = -5246665.6480784754*eccentricity_18 + 6450317.559779074*eccentricity_16 - 5788420.2515178157*eccentricity_14 + 3622721.0416216716*eccentricity_12 - 1495361.3327643824*eccentricity_10 + 376659.80988226997*eccentricity_8 - 50743.262239583333*eccentricity_6 + 2689.7213541666667*eccentricity_4;
    double common_term_23 = -19192592.348577977*eccentricity_19 + 22454320.432723807*eccentricity_17 - 19463273.081460568*eccentricity_15 + 11976103.482719264*eccentricity_13 - 4957815.5905671296*eccentricity_11 + 1280231.6090277778*eccentricity_9 - 181294.75972222222*eccentricity_7 + 10416.933333333333*eccentricity_5;
    double common_term_24 = 73678225.729192099*eccentricity_18 - 61722424.500421683*eccentricity_16 + 37247077.490721536*eccentricity_14 - 15371773.469458989*eccentricity_12 + 4028436.3722256252*eccentricity_10 - 590554.91350446429*eccentricity_8 + 35952.3310546875*eccentricity_6;
    double common_term_25 = 229519876.91698793*eccentricity_19 - 185961356.43962165*eccentricity_17 + 109880512.00473121*eccentricity_15 - 45014191.050091306*eccentricity_13 + 11885168.934584436*eccentricity_11 - 1784019.4609126984*eccentricity_9 + 113262.81408730159*eccentricity_7;
    double common_term_26 = -535506168.42425663*eccentricity_18 + 309502496.24738192*eccentricity_16 - 125471238.23406481*eccentricity_14 + 33196912.336670863*eccentricity_12 - 5061190.124490238*eccentricity_10 + 331262.53957180447*eccentricity_8;
    double common_term_27 = -1481278187.2168365*eccentricity_19 + 836886426.92802342*eccentricity_17 - 334977290.69402724*eccentricity_15 + 88450810.42967989*eccentricity_13 - 13613443.384017857*eccentricity_11 + 910742.26830357143*eccentricity_9;
    double common_term_28 = 2182094312.7773202*eccentricity_18 - 860936138.22354256*eccentricity_16 + 226175522.16480095*eccentricity_14 - 34977324.016420916*eccentricity_12 + 2376158.6679388382*eccentricity_10;
    double common_term_29 = 5507060858.2080425*eccentricity_19 - 2139135270.8775515*eccentricity_17 + 557784403.3173486*eccentricity_15 - 86356379.447419904*eccentricity_13 + 5926901.1560890402*eccentricity_11;
    double common_term_30 = -5156425810.5370692*eccentricity_18 + 1332074506.7675574*eccentricity_16 - 205869261.68385361*eccentricity_14 + 14217267.973173739*eccentricity_12;
    double common_term_31 = -12094838998.203366*eccentricity_19 + 3091059323.2509553*eccentricity_17 - 475782768.19529186*eccentricity_15 + 32955001.403697092*eccentricity_13;
    double common_term_32 = 6989565013.2417061*eccentricity_18 - 1069525307.0654351*eccentricity_16 + 74106571.639012238*eccentricity_14;
    double common_term_33 = 15439094778.143272*eccentricity_19 - 2345084320.84629*eccentricity_17 + 162200290.83477009*eccentricity_15;
    double common_term_34 = -5027485579.9932091*eccentricity_18 + 346506688.03548021*eccentricity_16;
    double common_term_35 = -10559963319.678437*eccentricity_19 + 724208296.45679759*eccentricity_17;
    double common_term_36 = 1483847817.5993173*eccentricity_18;
    double common_term_37 = 2985746258.5454742*eccentricity_19;
    double common_term_38 = 2.4907058686332918*eccentricity_19;
    double common_term_39 = 1.8313368967413934*eccentricity_18;
    double common_term_40 = 8.0709256249455359*eccentricity_19 + 1.3466438344771134*eccentricity_17;
    double common_term_41 = 6.1846903281442202*eccentricity_18 + 0.99032312947678878*eccentricity_16;
    double common_term_42 = 18.287382056415236*eccentricity_19 + 4.73182144282708*eccentricity_17 + 0.72835150256726604*eccentricity_15;
    double common_term_43 = 14.358148289124904*eccentricity_18 + 3.6149999889130043*eccentricity_16 + 0.53572790204492376*eccentricity_14;
    double common_term_44 = 34.817952113089896*eccentricity_19 + 11.253273875468839*eccentricity_17 + 2.7580669191299772*eccentricity_15 + 0.39408009958791209*eccentricity_13;
    double common_term_45 = 27.820107928660591*eccentricity_18 + 8.8049263926742687*eccentricity_16 + 2.1016420987547041*eccentricity_14 + 0.28990642380073486*eccentricity_12;
    double common_term_46 = 59.589320839363031*eccentricity_19 + 22.189434402682274*eccentricity_17 + 6.878124995985239*eccentricity_15 + 1.5995791120530704*eccentricity_13 + 0.21328405583613917*eccentricity_11;
    double common_term_47 = 48.26708476739048*eccentricity_18 + 17.66802327570739*eccentricity_16 + 5.3646493616816285*eccentricity_14 + 1.2161206906801694*eccentricity_12 + 0.15692005702427455*eccentricity_10;
    double common_term_48 = 94.780159299402032*eccentricity_19 + 39.029165743428891*eccentricity_17 + 14.044498808215533*eccentricity_15 + 4.177992148669232*eccentricity_13 + 0.92363233024691358*eccentricity_11 + 0.11545414462081129*eccentricity_9;
    double common_term_49 = 77.630837356242457*eccentricity_18 + 31.5064525940506*eccentricity_16 + 11.146105815556053*eccentricity_14 + 3.2491962517295985*eccentricity_12 + 0.70080703654617229*eccentricity_10 + 0.084946308438740079*eccentricity_8;
    double common_term_50 = 0.0625*eccentricity_7*std::pow(1.0 - eccentricity_2, -8.5);
    double common_term_51 = 118.08732009962402*eccentricity_18 + 51.825357980413596*eccentricity_16 + 20.43161983094378*eccentricity_14 + 6.9877937891005348*eccentricity_12 + 1.9572340344625806*eccentricity_10 + 0.40236855158730159*eccentricity_8 + 0.045985243055555556*eccentricity_6;
    double common_term_52 = 206.44998450259302*eccentricity_19 + 97.47765277209756*eccentricity_17 + 42.243137116034685*eccentricity_15 + 16.414319816468254*eccentricity_13 + 5.5202835648148148*eccentricity_11 + 1.5158730158730159*eccentricity_9 + 0.30416666666666667*eccentricity_7 + 0.033333333333333333*eccentricity_5;
    double common_term_53 = 172.06455098867397*eccentricity_18 + 80.340918642141989*eccentricity_16 + 34.379260070664542*eccentricity_14 + 13.167109646115984*eccentricity_12 + 4.3550973074776786*eccentricity_10 + 1.1750732421875*eccentricity_8 + 0.22109375*eccentricity_6 + 0.0546875*eccentricity_4;
    double common_term_54 = 288.62921884517223*eccentricity_19 + 143.19605863093051*eccentricity_17 + 66.120573643873702*eccentricity_15 + 27.94325355489418*eccentricity_13 + 10.544762731481481*eccentricity_11 + 3.5074074074074074*eccentricity_9 + 0.63333333333333333*eccentricity_7 + 0.83333333333333333*eccentricity_5 - 0.33333333333333333*eccentricity_3;
    double common_term_55 = 242.27217132009859*eccentricity_18 + 119.01075868919633*eccentricity_16 + 54.39085740854921*eccentricity_14 + 22.467579957669374*eccentricity_12 + 9.7438227335611979*eccentricity_10 - 1.5199652777777778*eccentricity_8 + 8.9547526041666667*eccentricity_6 - 6.125*eccentricity_4 + 1.375*eccentricity_2;
    double common_term_56 = 392.71013272655253*eccentricity_19 + 203.04252242107781*eccentricity_17 + 99.443281648596939*eccentricity_15 + 41.150580357142857*eccentricity_13 + 33.4621875*eccentricity_11 - 35.2265625*eccentricity_9 + 73.875*eccentricity_7 - 59.375*eccentricity_5 + 20.5*eccentricity_3 - 2.0*eccentricity;
    double common_term_57 = 330.58855481583915*eccentricity_18 + 177.84804273777217*eccentricity_16 + 45.619802982424512*eccentricity_14 + 165.04652462193995*eccentricity_12 - 295.80602715386285*eccentricity_10 + 484.53443060980903*eccentricity_8 - 405.70225694444444*eccentricity_6 + 166.484375*eccentricity_4 - 26.5*eccentricity_2 + 1.0;
    double common_term_58 = 508.82817162978021*eccentricity_19 + 351.31774959337892*eccentricity_17 - 144.20446824461855*eccentricity_15 + 935.22048776455026*eccentricity_13 - 1826.0444675925926*eccentricity_11 + 2635.5152777777778*eccentricity_9 - 2196.3125*eccentricity_7 + 975.16666666666667*eccentricity_5 - 193.0*eccentricity_3 + 12.0*eccentricity;
    double common_term_59 = 969.52930284730438*eccentricity_18 - 1602.5445929334106*eccentricity_16 + 5009.9082951394154*eccentricity_14 - 9351.8743494088309*eccentricity_12 + 12333.765953063965*eccentricity_10 - 10033.50859375*eccentricity_8 + 4608.4150390625*eccentricity_6 - 1022.875*eccentricity_4 + 80.125*eccentricity_2;
    double common_term_60 = 3927.0928011800378*eccentricity_19 - 9669.8334979843857*eccentricity_17 + 24103.320728634656*eccentricity_15 - 41759.86802972057*eccentricity_13 + 51101.179108796296*eccentricity_11 - 40225.29369212963*eccentricity_9 + 18651.666666666667*eccentricity_7 - 4409.4583333333333*eccentricity_5 + 392.83333333333333*eccentricity_3;
    double common_term_61 = -47650.850928130657*eccentricity_18 + 104332.15815067429*eccentricity_16 - 167303.10199640734*eccentricity_14 + 191554.71031180614*eccentricity_12 - 145321.01616462829*eccentricity_10 + 67061.771546766493*eccentricity_8 - 16397.329947916667*eccentricity_6 + 1577.4609375*eccentricity_4;
    double common_term_62 = -207258.96951711112*eccentricity_19 + 411502.47848001219*eccentricity_17 - 613080.9355033482*eccentricity_15 + 660498.54988839286*eccentricity_13 - 482114.6390625*eccentricity_11 + 219503.57142857143*eccentricity_9 - 54513.45*eccentricity_7 + 5494.8*eccentricity_5;
    double common_term_63 = 1498143.6594918641*eccentricity_18 - 2084242.6940495565*eccentricity_16 + 2121958.5393176698*eccentricity_14 - 1489592.6070817632*eccentricity_12 + 665473.88377873254*eccentricity_10 - 165903.71507936508*eccentricity_8 + 17190.395334201389*eccentricity_6;
    double common_term_64 = 5090510.9657484841*eccentricity_19 - 6645313.255311166*eccentricity_17 + 6416130.0117816986*eccentricity_15 - 4332675.3715330826*eccentricity_13 + 1892730.2281442901*eccentricity_11 - 469949.86909722222*eccentricity_9 + 49423.220833333333*eccentricity_7;
    double common_term_65 = -20042472.761940465*eccentricity_18 + 18407284.830869912*eccentricity_16 - 11964284.546877746*eccentricity_14 + 5099659.6063139098*eccentricity_12 - 1254253.8821781703*eccentricity_10 + 132700.77527291434*eccentricity_8;
    double common_term_66 = -57579675.748889606*eccentricity_19 + 50435491.547375227*eccentricity_17 - 31578949.293618405*eccentricity_15 + 13115741.563279045*eccentricity_13 - 3183256.661253858*eccentricity_11 + 336685.63266093475*eccentricity_9;
    double common_term_67 = 132695961.37450618*eccentricity_18 - 80108842.593498572*eccentricity_16 + 32395259.567587899*eccentricity_14 - 7738082.4555637731*eccentricity_12 + 814424.14075698691*eccentricity_10;
    double common_term_68 = 336749336.872071*eccentricity_19 - 196204509.35269591*eccentricity_17 + 77223485.483624559*eccentricity_15 - 18119802.234710836*eccentricity_13 + 1891272.7499431818*eccentricity_11;
    double common_term_69 = -465728105.16283728*eccentricity_18 + 178386676.49951927*eccentricity_16 - 41062126.957570431*eccentricity_14 + 4239530.1961378482*eccentricity_12;
    double common_term_70 = -1074839165.2722165*eccentricity_19 + 400676950.94261515*eccentricity_17 - 90395849.290232592*eccentricity_15 + 9214381.6287820975*eccentricity_13;
    double common_term_71 = 877584710.54147598*eccentricity_18 - 193931853.60943599*eccentricity_16 + 19488502.621578406*eccentricity_14;
    double common_term_72 = 1878911649.8394871*eccentricity_19 - 406537619.32195164*eccentricity_17 + 40231551.29956322*eccentricity_15;
    double common_term_73 = -834617408.82347331*eccentricity_18 + 81270731.34323376*eccentricity_16;
    double common_term_74 = -1681332480.0663107*eccentricity_19 + 160996779.73043899*eccentricity_17;
    double common_term_75 = 313341638.38959954*eccentricity_18;
    double common_term_76 = 600108297.80236725*eccentricity_19;
    double common_term_77 = 235.84517689085505*eccentricity_19;
    double common_term_78 = 167.83955953688641*eccentricity_18;
    double common_term_79 = 598.57698142400433*eccentricity_19 + 119.31343413488101*eccentricity_17;
    double common_term_80 = 449.38756471507317*eccentricity_18 + 84.718329023696355*eccentricity_16;
    double common_term_81 = 1180.5777365832363*eccentricity_19 + 336.12471570424968*eccentricity_17 + 60.078592217973909*eccentricity_15;
    double common_term_82 = 909.29183764206608*eccentricity_18 + 250.51655423086687*eccentricity_16 + 42.547407043743222*eccentricity_14;
    double common_term_83 = 2028.2429581290529*eccentricity_19 + 697.94387528096177*eccentricity_17 + 186.07506506565364*eccentricity_15 + 30.087678469774824*eccentricity_13;
    double common_term_84 = 1590.7757864621217*eccentricity_18 + 533.89681016945318*eccentricity_16 + 137.75201855709393*eccentricity_14 + 21.242784448513769*eccentricity_12;
    double common_term_85 = 3198.5706405927098*eccentricity_19 + 1243.628244307311*eccentricity_17 + 407.02550889052354*eccentricity_15 + 101.64655032467532*eccentricity_13 + 14.971971387987013*eccentricity_11;
    double common_term_86 = 2543.9227297253765*eccentricity_18 + 969.08299601944427*eccentricity_16 + 309.25502990227388*eccentricity_14 + 74.762625074241336*eccentricity_12 + 10.532229771689763*eccentricity_10;
    double common_term_87 = 4752.9891022223591*eccentricity_19 + 2017.0881032109421*eccentricity_17 + 752.68564725392706*eccentricity_15 + 234.17501628387045*eccentricity_13 + 54.811539627425044*eccentricity_11 + 7.3935736331569665*eccentricity_9;
    double common_term_88 = 3823.146112365726*eccentricity_18 + 1594.4366268807079*eccentricity_16 + 582.69145248524554*eccentricity_14 + 176.71940148217337*eccentricity_12 + 40.053403581891741*eccentricity_10 + 5.178314208984375*eccentricity_8;
    double common_term_89 = 6757.7455816546254*eccentricity_19 + 3066.339340736588*eccentricity_17 + 1256.4267701558826*eccentricity_15 + 449.59552228009259*eccentricity_13 + 132.90198894951499*eccentricity_11 + 29.171378968253968*eccentricity_9 + 3.6175595238095238*eccentricity_7;
    double common_term_90 = 5487.3723248062063*eccentricity_18 + 2452.1610198754915*eccentricity_16 + 986.9559695809151*eccentricity_14 + 345.73604182983973*eccentricity_12 + 99.599545130653987*eccentricity_10 + 21.172867063492063*eccentricity_8 + 2.5200737847222222*eccentricity_6;
    double common_term_91 = std::pow(1.0 - eccentricity_2, -8.5)*(0.4375*eccentricity_7 + 1.75*eccentricity_5);
    double common_term_92 = 7600.1206599722404*eccentricity_18 + 3588.3911554325066*eccentricity_16 + 1554.2469358595731*eccentricity_14 + 603.14970095321615*eccentricity_12 + 202.3493410140749*eccentricity_10 + 55.334242078993056*eccentricity_8 + 11.032552083333333*eccentricity_6 + 1.2109375*eccentricity_4;
    double common_term_93 = 12407.819254189161*eccentricity_19 + 6205.4493039549483*eccentricity_17 + 2889.5740294273895*eccentricity_15 + 1231.7354021990741*eccentricity_13 + 469.17746569113757*eccentricity_11 + 153.97974537037037*eccentricity_9 + 41.010416666666667*eccentricity_7 + 7.9166666666666667*eccentricity_5 + 0.83333333333333333*eccentricity_3;
    double common_term_94 = 10229.583354917616*eccentricity_18 + 5053.2739611229118*eccentricity_16 + 2320.1474348672799*eccentricity_14 + 973.09070412772042*eccentricity_12 + 363.72112274169922*eccentricity_10 + 116.73828125*eccentricity_8 + 30.2783203125*eccentricity_6 + 5.625*eccentricity_4 + 0.625*eccentricity_2;
    double common_term_95 = 16210.430653642303*eccentricity_19 + 8412.7525818597276*eccentricity_17 + 4103.9084452849427*eccentricity_15 + 1857.4591352513228*eccentricity_13 + 766.29846643518519*eccentricity_11 + 280.95486111111111*eccentricity_9 + 88.333333333333333*eccentricity_7 + 21.666666666666667*eccentricity_5 + 5.0*eccentricity_3;
    double common_term_96 = 13448.71624525219*eccentricity_18 + 6901.0611807801776*eccentricity_16 + 3323.7105976727304*eccentricity_14 + 1482.6555449285625*eccentricity_12 + 600.97017076280382*eccentricity_10 + 218.55689154730903*eccentricity_8 + 59.943576388888889*eccentricity_6 + 25.859375*eccentricity_4 - 2.5*eccentricity_2 + 1.0;
    double common_term_97 = 20778.155729899504*eccentricity_19 + 11130.946955466757*eccentricity_17 + 5646.2195187938457*eccentricity_15 + 2685.3772600446429*eccentricity_13 + 1174.148046875*eccentricity_11 + 490.8359375*eccentricity_9 + 117.96875*eccentricity_7 + 114.375*eccentricity_5 - 27.5*eccentricity_3 + 10.0*eccentricity;
    double common_term_98 = 17335.640310352144*eccentricity_18 + 9188.3402372189019*eccentricity_16 + 4618.2684950833703*eccentricity_14 + 2116.8987029464157*eccentricity_12 + 1071.9672554863824*eccentricity_10 + 89.620659722222222*eccentricity_8 + 466.15397135416667*eccentricity_6 - 170.625*eccentricity_4 + 56.875*eccentricity_2;
    double common_term_99 = 26205.524292838645*eccentricity_19 + 14413.143374118439*eccentricity_17 + 7648.0516265803357*eccentricity_15 + 3468.0038804150132*eccentricity_13 + 2473.6350033068783*eccentricity_11 - 584.38425925925926*eccentricity_9 + 1766.2708333333333*eccentricity_7 - 785.41666666666667*eccentricity_5 + 241.66666666666667*eccentricity_3;
    double common_term_100 = 21853.58427483769*eccentricity_18 + 12460.536510933191*eccentricity_16 + 4739.7066311665944*eccentricity_14 + 6444.1837419441768*eccentricity_12 - 4073.581551688058*eccentricity_10 + 6200.8631591796875*eccentricity_8 - 2985.25390625*eccentricity_6 + 852.9296875*eccentricity_4;
    double common_term_101 = 31882.346906943309*eccentricity_19 + 20838.182255449299*eccentricity_17 + 3211.1375227204353*eccentricity_15 + 18883.091474798832*eccentricity_13 - 17599.995143022487*eccentricity_11 + 20183.854786706349*eccentricity_9 - 9906.40625*eccentricity_7 + 2641.7083333333333*eccentricity_5;
    double common_term_102 = 38535.074565669105*eccentricity_18 - 11827.907447665517*eccentricity_16 + 58598.085660153467*eccentricity_14 - 62887.841761699499*eccentricity_12 + 61227.443857889327*eccentricity_10 - 29681.301649305556*eccentricity_8 + 7420.3254123263889*eccentricity_6;
    double common_term_103 = 85053.62804178366*eccentricity_19 - 76856.735861150566*eccentricity_17 + 181829.25168171672*eccentricity_15 - 200504.915625*eccentricity_13 + 174330.13950892857*eccentricity_11 - 82084.044642857143*eccentricity_9 + 19313.892857142857*eccentricity_7;
    double common_term_104 = -305322.75131784799*eccentricity_18 + 547004.94150160467*eccentricity_16 - 589381.7735199361*eccentricity_14 + 469335.64014766145*eccentricity_12 - 212777.2675541469*eccentricity_10 + 47283.837234254867*eccentricity_8;
    double common_term_105 = -1026768.2930192882*eccentricity_19 + 1577453.0443238732*eccentricity_17 - 1626146.1662386183*eccentricity_15 + 1202919.9477901848*eccentricity_13 - 522822.35218942901*eccentricity_11 + 110067.86361882716*eccentricity_9;
    double common_term_106 = 4353294.8114436442*eccentricity_18 - 4259448.2046274131*eccentricity_16 + 2952785.6003944421*eccentricity_14 - 1228061.0608972822*eccentricity_12 + 245612.21217945644*eccentricity_10;
    double common_term_107 = 11521464.423158283*eccentricity_19 - 10676772.751240997*eccentricity_17 + 6977572.3550935798*eccentricity_15 - 2775681.8376196865*eccentricity_13 + 528701.3024037498*eccentricity_11;
    double common_term_108 = -25763965.701112116*eccentricity_18 + 15942600.254913664*eccentricity_16 - 6068176.5279497196*eccentricity_14 + 1103304.8232635854*eccentricity_12;
    double common_term_109 = -60131777.868279155*eccentricity_19 + 35352270.19040307*eccentricity_17 - 12885593.777100341*eccentricity_15 + 2240972.8308000593*eccentricity_13;
    double common_term_110 = 76323925.243568762*eccentricity_18 - 26668496.000657048*eccentricity_16 + 4444749.3334428413*eccentricity_14;
    double common_term_111 = 160868858.60231702*eccentricity_19 - 53948287.268245645*eccentricity_17 + 8631725.9629193032*eccentricity_15;
    double common_term_112 = -106925425.75607861*eccentricity_18 + 16450065.50093517*eccentricity_16;
    double common_term_113 = -208060778.20081929*eccentricity_19 + 30823818.992713969*eccentricity_17;
    double common_term_114 = 56880344.698169911*eccentricity_18;
    double common_term_115 = 103515326.90073651*eccentricity_19;
    double common_term_116 = 7800.2099694157223*eccentricity_19;
    double common_term_117 = 5297.1040138897926*eccentricity_18;
    double common_term_118 = 12776.585536552788*eccentricity_19 + 3584.8466887409968*eccentricity_17;
    double common_term_119 = 9385.4958769092867*eccentricity_18 + 2416.8896573311643*eccentricity_16;
    double common_term_120 = 20693.332264842291*eccentricity_19 + 6829.8070204824472*eccentricity_17 + 1622.6627042343817*eccentricity_15;
    double common_term_121 = 15484.864928150838*eccentricity_18 + 4925.6142737932461*eccentricity_16 + 1084.3946115033133*eccentricity_14;
    double common_term_122 = 30303.330537974639*eccentricity_19 + 11507.454029140329*eccentricity_17 + 3521.3746621807351*eccentricity_15 + 720.94310661568691*eccentricity_13;
    double common_term_123 = 23037.553490881856*eccentricity_18 + 8489.7626271724033*eccentricity_16 + 2495.6297662454838*eccentricity_14 + 476.53095888992409*eccentricity_12;
    double common_term_124 = 41846.851220701574*eccentricity_19 + 17396.724649284498*eccentricity_17 + 6215.745949760558*eccentricity_15 + 1753.0762655677977*eccentricity_13 + 312.91185608565817*eccentricity_11;
    double common_term_125 = 32182.119592967761*eccentricity_18 + 13044.463546343242*eccentricity_16 + 4514.2273971522006*eccentricity_14 + 1220.1931485086392*eccentricity_12 + 203.93073070573218*eccentricity_10;
    double common_term_126 = 55394.557859176912*eccentricity_19 + 24587.897815949919*eccentricity_17 + 9707.8936330027394*eccentricity_15 + 3250.4038735288149*eccentricity_13 + 841.059140625*eccentricity_11 + 131.75334821428571*eccentricity_9;
    double common_term_127 = 42982.785449612291*eccentricity_18 + 18655.25171315353*eccentricity_16 + 7167.0033690675334*eccentricity_14 + 2318.8363587372728*eccentricity_12 + 573.65326247980565*eccentricity_10 + 84.257838173518105*eccentricity_8;
    double common_term_128 = 71009.966149198223*eccentricity_19 + 33136.402244361205*eccentricity_17 + 14048.734611951192*eccentricity_15 + 5245.5033547361846*eccentricity_13 + 1637.6462783840388*eccentricity_11 + 386.73730158730159*eccentricity_9 + 53.233730158730159*eccentricity_7;
    double common_term_129 = 55492.549945840973*eccentricity_18 + 25368.567032049866*eccentricity_16 + 10494.656793201481*eccentricity_14 + 3803.0297357831682*eccentricity_12 + 1143.7223088945661*eccentricity_10 + 257.31462053571429*eccentricity_8 + 33.1412109375*eccentricity_6;
    double common_term_130 = 88741.621037334954*eccentricity_19 + 43085.692967238858*eccentricity_17 + 19276.262789900702*eccentricity_15 + 7770.94570737342*eccentricity_13 + 2728.5433428406085*eccentricity_11 + 788.78559027777778*eccentricity_9 + 168.60902777777778*eccentricity_7 + 20.258333333333333*eccentricity_5;
    double common_term_131 = 69750.379423513786*eccentricity_18 + 33219.212839995246*eccentricity_16 + 14527.494697299065*eccentricity_14 + 5698.4299680742637*eccentricity_12 + 1934.7809803473256*eccentricity_10 + 536.18533799913194*eccentricity_8 + 108.48776041666667*eccentricity_6 + 12.096354166666667*eccentricity_4;
    double common_term_132 = std::pow(1.0 - eccentricity_2, -8.5)*(1.3125*eccentricity_7 + 8.75*eccentricity_5 + 7.0*eccentricity_3);
    double common_term_133 = 85780.621966221238*eccentricity_18 + 42229.597467844475*eccentricity_16 + 19285.542654869672*eccentricity_14 + 8022.3546295828037*eccentricity_12 + 2961.233682985659*eccentricity_10 + 932.26788194444444*eccentricity_8 + 234.53157552083333*eccentricity_6 + 41.708333333333333*eccentricity_4 + 3.875*eccentricity_2;
    double common_term_134 = 130669.90477486832*eccentricity_19 + 67293.890292321679*eccentricity_17 + 32487.76775017617*eccentricity_15 + 14496.956052551808*eccentricity_13 + 5864.1605237268518*eccentricity_11 + 2091.1480902777778*eccentricity_9 + 630.05902777777778*eccentricity_7 + 149.54166666666667*eccentricity_5 + 24.5*eccentricity_3 + 2.0*eccentricity;
    double common_term_135 = 103592.40651101234*eccentricity_18 + 52409.13291403315*eccentricity_16 + 24777.936440142807*eccentricity_14 + 10783.215084037781*eccentricity_12 + 4230.4302575683594*eccentricity_10 + 1451.6248168945313*eccentricity_8 + 415.9296875*eccentricity_6 + 92.109375*eccentricity_4 + 13.5*eccentricity_2 + 1.0;
    double common_term_136 = 154883.87155756756*eccentricity_19 + 81571.343550757438*eccentricity_17 + 40490.270220576657*eccentricity_15 + 18715.021332396384*eccentricity_13 + 7924.5632581018519*eccentricity_11 + 3004.7777777777778*eccentricity_9 + 986.65277777777778*eccentricity_7 + 266.66666666666667*eccentricity_5 + 53.0*eccentricity_3 + 8.0*eccentricity;
    double common_term_137 = 123179.05201263452*eccentricity_18 + 63753.643896553193*eccentricity_16 + 31002.333077499843*eccentricity_14 + 13979.906038032007*eccentricity_12 + 5742.2284220942744*eccentricity_10 + 2093.5160590277778*eccentricity_8 + 656.27571614583333*eccentricity_6 + 155.95833333333333*eccentricity_4 + 37.625*eccentricity_2;
    double common_term_138 = 181247.58217219784*eccentricity_19 + 97285.494608074578*eccentricity_17 + 49415.366081891741*eccentricity_15 + 23496.982960379464*eccentricity_13 + 10311.072488839286*eccentricity_11 + 4083.69375*eccentricity_9 + 1447.70625*eccentricity_7 + 375.0*eccentricity_5 + 135.5*eccentricity_3;
    double common_term_139 = 144517.44825226858*eccentricity_18 + 76245.177381701057*eccentricity_16 + 37941.912039481642*eccentricity_14 + 17613.32027336933*eccentricity_12 + 7443.6734243887442*eccentricity_10 + 2966.9605007595486*eccentricity_8 + 757.21588541666667*eccentricity_6 + 413.12760416666667*eccentricity_4;
    double common_term_140 = 209725.81751224776*eccentricity_19 + 114411.1257299214*eccentricity_17 + 59222.933874704907*eccentricity_15 + 28895.653577054857*eccentricity_13 + 12785.68521412037*eccentricity_11 + 5818.1884920634921*eccentricity_9 + 1257.9152777777778*eccentricity_7 + 1122.2833333333333*eccentricity_5;
    double common_term_141 = 167589.1850516713*eccentricity_18 + 89752.534770760468*eccentricity_16 + 45921.587596170391*eccentricity_14 + 20701.191739000593*eccentricity_12 + 11220.506275177002*eccentricity_10 + 1502.6236607142857*eccentricity_8 + 2799.4181640625*eccentricity_6;
    double common_term_142 = 240386.8143436027*eccentricity_19 + 132423.2588344265*eccentricity_17 + 71387.39107146902*eccentricity_15 + 31226.856321591619*eccentricity_13 + 21819.054859733245*eccentricity_11 + 216.5906498015873*eccentricity_9 + 6537.5114087301587*eccentricity_7;
    double common_term_143 = 190265.77912825364*eccentricity_18 + 110031.17151171312*eccentricity_16 + 42294.035758734097*eccentricity_14 + 43575.455233026972*eccentricity_12 - 6033.3518224914658*eccentricity_10 + 14485.249283951048*eccentricity_8;
    double common_term_144 = 265111.43977515882*eccentricity_19 + 171601.62748469304*eccentricity_17 + 45696.761713169643*eccentricity_15 + 90006.839086850649*eccentricity_13 - 25460.433883928571*eccentricity_11 + 30744.7*eccentricity_9;
    double common_term_145 = 278544.51341060726*eccentricity_18 + 17322.217296834574*eccentricity_16 + 191343.08058544473*eccentricity_14 - 76296.55152307011*eccentricity_12 + 62957.192477849997*eccentricity_10;
    double common_term_146 = 485082.73177340011*eccentricity_19 - 105319.14793341796*eccentricity_17 + 413654.63490926953*eccentricity_15 - 196923.72650306596*eccentricity_13 + 125062.34629647166*eccentricity_11;
    double common_term_147 = -472680.93989429507*eccentricity_18 + 897642.34882705052*eccentricity_16 - 464957.44270754158*eccentricity_14 + 242034.08835602736*eccentricity_12;
    double common_term_148 = -1428928.0582125809*eccentricity_19 + 1935289.0893354301*eccentricity_17 - 1032410.5904093375*eccentricity_15 + 457914.39447379427*eccentricity_13;
    double common_term_149 = 4117928.7876018397*eccentricity_18 - 2189291.3725260587*eccentricity_16 + 849297.60792072932*eccentricity_14;
    double common_term_150 = 8616610.6516342617*eccentricity_19 - 4476790.8883498675*eccentricity_17 + 1547749.3085034558*eccentricity_15;
    double common_term_151 = -8885911.5430669566*eccentricity_18 + 2776760.7867697294*eccentricity_16;
    double common_term_152 = -17201719.791570884*eccentricity_19 + 4912189.3483176188*eccentricity_17;
    double common_term_153 = 8580406.8363049765*eccentricity_18;
    double common_term_154 = 14816639.982330246*eccentricity_19;
    double common_term_155 = 140602.26694199336*eccentricity_19;
    double common_term_156 = 90637.080451187001*eccentricity_18;
    double common_term_157 = 96203.828192372159*eccentricity_19 + 58068.023344054002*eccentricity_17;
    double common_term_158 = 73583.520588604108*eccentricity_18 + 36949.581941120817*eccentricity_16;
    double common_term_159 = 146350.25228201067*eccentricity_19 + 54486.699573288415*eccentricity_17 + 23334.633009545894*eccentricity_15;
    double common_term_160 = 106269.27074596688*eccentricity_18 + 39293.813193431935*eccentricity_16 + 14612.741020659801*eccentricity_14;
    double common_term_161 = 175599.62720402224*eccentricity_19 + 76530.47537321607*eccentricity_17 + 27696.684354445777*eccentricity_15 + 9064.6274672202797*eccentricity_13;
    double common_term_162 = 129625.06500202773*eccentricity_18 + 54551.943543280274*eccentricity_16 + 19120.129839481894*eccentricity_14 + 5563.0408091910035*eccentricity_12;
    double common_term_163 = 207114.91626794092*eccentricity_19 + 94627.580378948607*eccentricity_17 + 38423.941374400901*eccentricity_15 + 12940.394119377681*eccentricity_13 + 3372.530453029802*eccentricity_11;
    double common_term_164 = 153797.83757834102*eccentricity_18 + 68259.796160641167*eccentricity_16 + 26701.327805336559*eccentricity_14 + 8587.7030786866027*eccentricity_12 + 2015.861515611921*eccentricity_10;
    double common_term_165 = 238058.07353049148*eccentricity_19 + 112956.61292842567*eccentricity_17 + 48606.116862490821*eccentricity_15 + 18277.500904193222*eccentricity_13 + 5585.2709025022046*eccentricity_11 + 1185.2050870811287*eccentricity_9;
    double common_term_166 = 177487.22629588829*eccentricity_18 + 81969.347545127098*eccentricity_16 + 34123.268146539057*eccentricity_14 + 12302.504021480466*eccentricity_12 + 3555.4650755705657*eccentricity_10 + 683.3218254937066*eccentricity_8;
    double common_term_167 = 268056.3925659237*eccentricity_19 + 130843.33524779639*eccentricity_17 + 58699.933590452516*eccentricity_15 + 23581.400142299107*eccentricity_13 + 8125.6586495535714*eccentricity_11 + 2210.6404017857143*eccentricity_9 + 384.77142857142857*eccentricity_7;
    double common_term_168 = 200307.06309585681*eccentricity_18 + 95264.393166060517*eccentricity_16 + 41421.029536145302*eccentricity_14 + 16010.410773990997*eccentricity_12 + 5252.626898326571*eccentricity_10 + 1338.2390873015873*eccentricity_8 + 210.44585503472222*eccentricity_6;
    double common_term_169 = 296638.08537890538*eccentricity_19 + 147934.19985312986*eccentricity_17 + 68404.996242495223*eccentricity_15 + 28747.072789214065*eccentricity_13 + 10652.907510747354*eccentricity_11 + 3311.768253968254*eccentricity_9 + 785.10833333333333*eccentricity_7 + 110.93333333333333*eccentricity_5;
    double common_term_170 = 221846.5355749358*eccentricity_18 + 107835.15564490939*eccentricity_16 + 48356.915566487312*eccentricity_14 + 19576.206907755988*eccentricity_12 + 6923.8570155552455*eccentricity_10 + 2027.1390380859375*eccentricity_8 + 443.33671875*eccentricity_6 + 55.7109375*eccentricity_4;
    double common_term_171 = 323339.96666629194*eccentricity_19 + 163873.24443622077*eccentricity_17 + 77455.204325808272*eccentricity_15 + 33579.477349537037*eccentricity_13 + 13040.228257275132*eccentricity_11 + 4376.4506365740741*eccentricity_9 + 1196.6708333333333*eccentricity_7 + 238.45833333333333*eccentricity_5 + 26.166666666666667*eccentricity_3;
    double common_term_172 = 241702.05727741603*eccentricity_18 + 119375.44527611406*eccentricity_16 + 54705.23853801273*eccentricity_14 + 22838.89973539968*eccentricity_12 + 8461.7476655748155*eccentricity_10 + 2673.5511284722222*eccentricity_8 + 674.60514322916667*eccentricity_6 + 120.125*eccentricity_4 + 11.125*eccentricity_2;
    double common_term_173 = std::pow(1.0 - eccentricity_2, -8.5)*(2.1875*eccentricity_7 + 17.5*eccentricity_5 + 21.0*eccentricity_3 + 4.0*eccentricity);
    double common_term_174 = 259478.39258943249*eccentricity_18 + 129585.90837285566*eccentricity_16 + 60245.641232276541*eccentricity_14 + 25641.887269770775*eccentricity_12 + 9759.8455874294705*eccentricity_10 + 3209.274664984809*eccentricity_8 + 863.58940972222222*eccentricity_6 + 173.234375*eccentricity_4 + 21.5*eccentricity_2 + 1.0;
    double common_term_175 = 369295.06319135044*eccentricity_19 + 190907.96783395968*eccentricity_17 + 92553.566940367437*eccentricity_15 + 41493.096462673611*eccentricity_13 + 16871.826087962963*eccentricity_11 + 6051.5210069444444*eccentricity_9 + 1835.2291666666667*eccentricity_7 + 438.79166666666667*eccentricity_5 + 72.5*eccentricity_3 + 6.0*eccentricity;
    double common_term_176 = 274790.10271172715*eccentricity_18 + 138175.46590017518*eccentricity_16 + 64764.59881095767*eccentricity_14 + 27834.041443437849*eccentricity_12 + 10715.910200500488*eccentricity_10 + 3569.0953125*eccentricity_8 + 971.8896484375*eccentricity_6 + 194.875*eccentricity_4 + 22.375*eccentricity_2;
    double common_term_177 = 387669.76682867132*eccentricity_19 + 201331.89622000637*eccentricity_17 + 98101.472920050093*eccentricity_15 + 44213.601309937169*eccentricity_13 + 18066.265777943122*eccentricity_11 + 6497.5703703703704*eccentricity_9 + 1961.3416666666667*eccentricity_7 + 455.16666666666667*eccentricity_5 + 66.333333333333333*eccentricity_3;
    double common_term_178 = 287262.98592478559*eccentricity_18 + 144862.75360138302*eccentricity_16 + 68056.858297803344*eccentricity_14 + 29271.151923332265*eccentricity_12 + 11233.310910179501*eccentricity_10 + 3692.5153944227431*eccentricity_8 + 963.32942708333333*eccentricity_6 + 171.0546875*eccentricity_4;
    double common_term_179 = 402410.51679016385*eccentricity_19 + 209263.38502796203*eccentricity_17 + 101996.90143470982*eccentricity_15 + 45882.353724888393*eccentricity_13 + 18623.578404017857*eccentricity_11 + 6583.5727678571429*eccentricity_9 + 1891.1*eccentricity_7 + 401.15*eccentricity_5;
    double common_term_180 = 296535.52018028997*eccentricity_18 + 149377.49793638151*eccentricity_16 + 69927.299509688076*eccentricity_14 + 29814.914947921062*eccentricity_12 + 11233.090602535672*eccentricity_10 + 3491.5160590277778*eccentricity_8 + 877.80379774305556*eccentricity_6;
    double common_term_181 = 413110.70086964009*eccentricity_19 + 214394.75921103658*eccentricity_17 + 104017.53081988714*eccentricity_15 + 46327.005182016093*eccentricity_13 + 18485.898155037478*eccentricity_11 + 6112.5556051587302*eccentricity_9 + 1821.7160714285714*eccentricity_7;
    double common_term_182 = 302257.24461886413*eccentricity_18 + 151477.42665135734*eccentricity_16 + 70129.440291661968*eccentricity_14 + 29532.771608513423*eccentricity_12 + 10188.909559413365*eccentricity_10 + 3625.3474897112165*eccentricity_8;
    double common_term_183 = 419363.04208000796*eccentricity_19 + 216503.79739556315*eccentricity_17 + 103701.73914198129*eccentricity_15 + 46078.888808643854*eccentricity_13 + 16179.124872547399*eccentricity_11 + 6972.9033702601411*eccentricity_9;
    double common_term_184 = 304384.69887690526*eccentricity_18 + 150024.97558675347*eccentricity_16 + 70655.822435672789*eccentricity_14 + 24385.021060059293*eccentricity_12 + 13037.374969302384*eccentricity_10;
    double common_term_185 = 421851.18113225513*eccentricity_19 + 212381.37284051328*eccentricity_17 + 107230.93044194087*eccentricity_15 + 34536.848508522727*eccentricity_13 + 23800.945304383117*eccentricity_11;
    double common_term_186 = 293684.07180556041*eccentricity_18 + 162419.59918544881*eccentricity_16 + 44935.450212509978*eccentricity_14 + 42571.446420460595*eccentricity_12;
    double common_term_187 = 394741.96341576399*eccentricity_19 + 247888.4907087542*eccentricity_17 + 50786.218467281376*eccentricity_15 + 74808.31916897676*eccentricity_13;
    double common_term_188 = 385053.17801633135*eccentricity_18 + 41097.796341673401*eccentricity_16 + 129434.23875676698*eccentricity_14;
    double common_term_189 = 614094.61459738136*eccentricity_19 - 6914.5285726141306*eccentricity_17 + 220904.65009367382*eccentricity_15;
    double common_term_190 = -137689.5266567441*eccentricity_18 + 372454.09995328634*eccentricity_16;
    double common_term_191 = -434363.28597093863*eccentricity_19 + 621161.2715737625*eccentricity_17;
    double common_term_192 = 1025812.6753226241*eccentricity_18;
    double common_term_193 = 1679055.7656150983*eccentricity_19;

    c_IntMap<c_Key1, double> result_by_q(39);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (9, 0).
    // q = -19
    result_by_lpq.set(c_Key3(9, 0, -19), common_term_0);
    result_by_q.set(c_Key1(-19), common_term_0);
    // q = -18
    result_by_lpq.set(c_Key3(9, 0, -18), common_term_1);
    result_by_q.set(c_Key1(-18), common_term_1);
    // q = -17
    result_by_lpq.set(c_Key3(9, 0, -17), common_term_2);
    result_by_q.set(c_Key1(-17), common_term_2);
    // q = -16
    result_by_lpq.set(c_Key3(9, 0, -16), common_term_3);
    result_by_q.set(c_Key1(-16), common_term_3);
    // q = -15
    result_by_lpq.set(c_Key3(9, 0, -15), common_term_4);
    result_by_q.set(c_Key1(-15), common_term_4);
    // q = -14
    result_by_lpq.set(c_Key3(9, 0, -14), common_term_5);
    result_by_q.set(c_Key1(-14), common_term_5);
    // q = -13
    result_by_lpq.set(c_Key3(9, 0, -13), common_term_6);
    result_by_q.set(c_Key1(-13), common_term_6);
    // q = -12
    result_by_lpq.set(c_Key3(9, 0, -12), common_term_7);
    result_by_q.set(c_Key1(-12), common_term_7);
    // q = -11
    result_by_lpq.set(c_Key3(9, 0, -11), common_term_8);
    result_by_q.set(c_Key1(-11), common_term_8);
    // q = -10
    result_by_lpq.set(c_Key3(9, 0, -10), common_term_9);
    result_by_q.set(c_Key1(-10), common_term_9);
    // q = -8
    result_by_lpq.set(c_Key3(9, 0, -8), common_term_10);
    result_by_q.set(c_Key1(-8), common_term_10);
    // q = -7
    result_by_lpq.set(c_Key3(9, 0, -7), common_term_11);
    result_by_q.set(c_Key1(-7), common_term_11);
    // q = -6
    result_by_lpq.set(c_Key3(9, 0, -6), common_term_12);
    result_by_q.set(c_Key1(-6), common_term_12);
    // q = -5
    result_by_lpq.set(c_Key3(9, 0, -5), common_term_13);
    result_by_q.set(c_Key1(-5), common_term_13);
    // q = -4
    result_by_lpq.set(c_Key3(9, 0, -4), common_term_14);
    result_by_q.set(c_Key1(-4), common_term_14);
    // q = -3
    result_by_lpq.set(c_Key3(9, 0, -3), common_term_15);
    result_by_q.set(c_Key1(-3), common_term_15);
    // q = -2
    result_by_lpq.set(c_Key3(9, 0, -2), common_term_16);
    result_by_q.set(c_Key1(-2), common_term_16);
    // q = -1
    result_by_lpq.set(c_Key3(9, 0, -1), common_term_17);
    result_by_q.set(c_Key1(-1), common_term_17);
    // q = 0
    result_by_lpq.set(c_Key3(9, 0, 0), common_term_18);
    result_by_q.set(c_Key1(0), common_term_18);
    // q = 1
    result_by_lpq.set(c_Key3(9, 0, 1), common_term_19);
    result_by_q.set(c_Key1(1), common_term_19);
    // q = 2
    result_by_lpq.set(c_Key3(9, 0, 2), common_term_20);
    result_by_q.set(c_Key1(2), common_term_20);
    // q = 3
    result_by_lpq.set(c_Key3(9, 0, 3), common_term_21);
    result_by_q.set(c_Key1(3), common_term_21);
    // q = 4
    result_by_lpq.set(c_Key3(9, 0, 4), common_term_22);
    result_by_q.set(c_Key1(4), common_term_22);
    // q = 5
    result_by_lpq.set(c_Key3(9, 0, 5), common_term_23);
    result_by_q.set(c_Key1(5), common_term_23);
    // q = 6
    result_by_lpq.set(c_Key3(9, 0, 6), common_term_24);
    result_by_q.set(c_Key1(6), common_term_24);
    // q = 7
    result_by_lpq.set(c_Key3(9, 0, 7), common_term_25);
    result_by_q.set(c_Key1(7), common_term_25);
    // q = 8
    result_by_lpq.set(c_Key3(9, 0, 8), common_term_26);
    result_by_q.set(c_Key1(8), common_term_26);
    // q = 9
    result_by_lpq.set(c_Key3(9, 0, 9), common_term_27);
    result_by_q.set(c_Key1(9), common_term_27);
    // q = 10
    result_by_lpq.set(c_Key3(9, 0, 10), common_term_28);
    result_by_q.set(c_Key1(10), common_term_28);
    // q = 11
    result_by_lpq.set(c_Key3(9, 0, 11), common_term_29);
    result_by_q.set(c_Key1(11), common_term_29);
    // q = 12
    result_by_lpq.set(c_Key3(9, 0, 12), common_term_30);
    result_by_q.set(c_Key1(12), common_term_30);
    // q = 13
    result_by_lpq.set(c_Key3(9, 0, 13), common_term_31);
    result_by_q.set(c_Key1(13), common_term_31);
    // q = 14
    result_by_lpq.set(c_Key3(9, 0, 14), common_term_32);
    result_by_q.set(c_Key1(14), common_term_32);
    // q = 15
    result_by_lpq.set(c_Key3(9, 0, 15), common_term_33);
    result_by_q.set(c_Key1(15), common_term_33);
    // q = 16
    result_by_lpq.set(c_Key3(9, 0, 16), common_term_34);
    result_by_q.set(c_Key1(16), common_term_34);
    // q = 17
    result_by_lpq.set(c_Key3(9, 0, 17), common_term_35);
    result_by_q.set(c_Key1(17), common_term_35);
    // q = 18
    result_by_lpq.set(c_Key3(9, 0, 18), common_term_36);
    result_by_q.set(c_Key1(18), common_term_36);
    // q = 19
    result_by_lpq.set(c_Key3(9, 0, 19), common_term_37);
    result_by_q.set(c_Key1(19), common_term_37);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 0), result_by_q);
    result_by_q.clear();

    // l , p = (9, 1).
    // q = -19
    result_by_lpq.set(c_Key3(9, 1, -19), common_term_38);
    result_by_q.set(c_Key1(-19), common_term_38);
    // q = -18
    result_by_lpq.set(c_Key3(9, 1, -18), common_term_39);
    result_by_q.set(c_Key1(-18), common_term_39);
    // q = -17
    result_by_lpq.set(c_Key3(9, 1, -17), common_term_40);
    result_by_q.set(c_Key1(-17), common_term_40);
    // q = -16
    result_by_lpq.set(c_Key3(9, 1, -16), common_term_41);
    result_by_q.set(c_Key1(-16), common_term_41);
    // q = -15
    result_by_lpq.set(c_Key3(9, 1, -15), common_term_42);
    result_by_q.set(c_Key1(-15), common_term_42);
    // q = -14
    result_by_lpq.set(c_Key3(9, 1, -14), common_term_43);
    result_by_q.set(c_Key1(-14), common_term_43);
    // q = -13
    result_by_lpq.set(c_Key3(9, 1, -13), common_term_44);
    result_by_q.set(c_Key1(-13), common_term_44);
    // q = -12
    result_by_lpq.set(c_Key3(9, 1, -12), common_term_45);
    result_by_q.set(c_Key1(-12), common_term_45);
    // q = -11
    result_by_lpq.set(c_Key3(9, 1, -11), common_term_46);
    result_by_q.set(c_Key1(-11), common_term_46);
    // q = -10
    result_by_lpq.set(c_Key3(9, 1, -10), common_term_47);
    result_by_q.set(c_Key1(-10), common_term_47);
    // q = -9
    result_by_lpq.set(c_Key3(9, 1, -9), common_term_48);
    result_by_q.set(c_Key1(-9), common_term_48);
    // q = -8
    result_by_lpq.set(c_Key3(9, 1, -8), common_term_49);
    result_by_q.set(c_Key1(-8), common_term_49);
    // q = -7
    result_by_lpq.set(c_Key3(9, 1, -7), common_term_50);
    result_by_q.set(c_Key1(-7), common_term_50);
    // q = -6
    result_by_lpq.set(c_Key3(9, 1, -6), common_term_51);
    result_by_q.set(c_Key1(-6), common_term_51);
    // q = -5
    result_by_lpq.set(c_Key3(9, 1, -5), common_term_52);
    result_by_q.set(c_Key1(-5), common_term_52);
    // q = -4
    result_by_lpq.set(c_Key3(9, 1, -4), common_term_53);
    result_by_q.set(c_Key1(-4), common_term_53);
    // q = -3
    result_by_lpq.set(c_Key3(9, 1, -3), common_term_54);
    result_by_q.set(c_Key1(-3), common_term_54);
    // q = -2
    result_by_lpq.set(c_Key3(9, 1, -2), common_term_55);
    result_by_q.set(c_Key1(-2), common_term_55);
    // q = -1
    result_by_lpq.set(c_Key3(9, 1, -1), common_term_56);
    result_by_q.set(c_Key1(-1), common_term_56);
    // q = 0
    result_by_lpq.set(c_Key3(9, 1, 0), common_term_57);
    result_by_q.set(c_Key1(0), common_term_57);
    // q = 1
    result_by_lpq.set(c_Key3(9, 1, 1), common_term_58);
    result_by_q.set(c_Key1(1), common_term_58);
    // q = 2
    result_by_lpq.set(c_Key3(9, 1, 2), common_term_59);
    result_by_q.set(c_Key1(2), common_term_59);
    // q = 3
    result_by_lpq.set(c_Key3(9, 1, 3), common_term_60);
    result_by_q.set(c_Key1(3), common_term_60);
    // q = 4
    result_by_lpq.set(c_Key3(9, 1, 4), common_term_61);
    result_by_q.set(c_Key1(4), common_term_61);
    // q = 5
    result_by_lpq.set(c_Key3(9, 1, 5), common_term_62);
    result_by_q.set(c_Key1(5), common_term_62);
    // q = 6
    result_by_lpq.set(c_Key3(9, 1, 6), common_term_63);
    result_by_q.set(c_Key1(6), common_term_63);
    // q = 7
    result_by_lpq.set(c_Key3(9, 1, 7), common_term_64);
    result_by_q.set(c_Key1(7), common_term_64);
    // q = 8
    result_by_lpq.set(c_Key3(9, 1, 8), common_term_65);
    result_by_q.set(c_Key1(8), common_term_65);
    // q = 9
    result_by_lpq.set(c_Key3(9, 1, 9), common_term_66);
    result_by_q.set(c_Key1(9), common_term_66);
    // q = 10
    result_by_lpq.set(c_Key3(9, 1, 10), common_term_67);
    result_by_q.set(c_Key1(10), common_term_67);
    // q = 11
    result_by_lpq.set(c_Key3(9, 1, 11), common_term_68);
    result_by_q.set(c_Key1(11), common_term_68);
    // q = 12
    result_by_lpq.set(c_Key3(9, 1, 12), common_term_69);
    result_by_q.set(c_Key1(12), common_term_69);
    // q = 13
    result_by_lpq.set(c_Key3(9, 1, 13), common_term_70);
    result_by_q.set(c_Key1(13), common_term_70);
    // q = 14
    result_by_lpq.set(c_Key3(9, 1, 14), common_term_71);
    result_by_q.set(c_Key1(14), common_term_71);
    // q = 15
    result_by_lpq.set(c_Key3(9, 1, 15), common_term_72);
    result_by_q.set(c_Key1(15), common_term_72);
    // q = 16
    result_by_lpq.set(c_Key3(9, 1, 16), common_term_73);
    result_by_q.set(c_Key1(16), common_term_73);
    // q = 17
    result_by_lpq.set(c_Key3(9, 1, 17), common_term_74);
    result_by_q.set(c_Key1(17), common_term_74);
    // q = 18
    result_by_lpq.set(c_Key3(9, 1, 18), common_term_75);
    result_by_q.set(c_Key1(18), common_term_75);
    // q = 19
    result_by_lpq.set(c_Key3(9, 1, 19), common_term_76);
    result_by_q.set(c_Key1(19), common_term_76);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 1), result_by_q);
    result_by_q.clear();

    // l , p = (9, 2).
    // q = -19
    result_by_lpq.set(c_Key3(9, 2, -19), common_term_77);
    result_by_q.set(c_Key1(-19), common_term_77);
    // q = -18
    result_by_lpq.set(c_Key3(9, 2, -18), common_term_78);
    result_by_q.set(c_Key1(-18), common_term_78);
    // q = -17
    result_by_lpq.set(c_Key3(9, 2, -17), common_term_79);
    result_by_q.set(c_Key1(-17), common_term_79);
    // q = -16
    result_by_lpq.set(c_Key3(9, 2, -16), common_term_80);
    result_by_q.set(c_Key1(-16), common_term_80);
    // q = -15
    result_by_lpq.set(c_Key3(9, 2, -15), common_term_81);
    result_by_q.set(c_Key1(-15), common_term_81);
    // q = -14
    result_by_lpq.set(c_Key3(9, 2, -14), common_term_82);
    result_by_q.set(c_Key1(-14), common_term_82);
    // q = -13
    result_by_lpq.set(c_Key3(9, 2, -13), common_term_83);
    result_by_q.set(c_Key1(-13), common_term_83);
    // q = -12
    result_by_lpq.set(c_Key3(9, 2, -12), common_term_84);
    result_by_q.set(c_Key1(-12), common_term_84);
    // q = -11
    result_by_lpq.set(c_Key3(9, 2, -11), common_term_85);
    result_by_q.set(c_Key1(-11), common_term_85);
    // q = -10
    result_by_lpq.set(c_Key3(9, 2, -10), common_term_86);
    result_by_q.set(c_Key1(-10), common_term_86);
    // q = -9
    result_by_lpq.set(c_Key3(9, 2, -9), common_term_87);
    result_by_q.set(c_Key1(-9), common_term_87);
    // q = -8
    result_by_lpq.set(c_Key3(9, 2, -8), common_term_88);
    result_by_q.set(c_Key1(-8), common_term_88);
    // q = -7
    result_by_lpq.set(c_Key3(9, 2, -7), common_term_89);
    result_by_q.set(c_Key1(-7), common_term_89);
    // q = -6
    result_by_lpq.set(c_Key3(9, 2, -6), common_term_90);
    result_by_q.set(c_Key1(-6), common_term_90);
    // q = -5
    result_by_lpq.set(c_Key3(9, 2, -5), common_term_91);
    result_by_q.set(c_Key1(-5), common_term_91);
    // q = -4
    result_by_lpq.set(c_Key3(9, 2, -4), common_term_92);
    result_by_q.set(c_Key1(-4), common_term_92);
    // q = -3
    result_by_lpq.set(c_Key3(9, 2, -3), common_term_93);
    result_by_q.set(c_Key1(-3), common_term_93);
    // q = -2
    result_by_lpq.set(c_Key3(9, 2, -2), common_term_94);
    result_by_q.set(c_Key1(-2), common_term_94);
    // q = -1
    result_by_lpq.set(c_Key3(9, 2, -1), common_term_95);
    result_by_q.set(c_Key1(-1), common_term_95);
    // q = 0
    result_by_lpq.set(c_Key3(9, 2, 0), common_term_96);
    result_by_q.set(c_Key1(0), common_term_96);
    // q = 1
    result_by_lpq.set(c_Key3(9, 2, 1), common_term_97);
    result_by_q.set(c_Key1(1), common_term_97);
    // q = 2
    result_by_lpq.set(c_Key3(9, 2, 2), common_term_98);
    result_by_q.set(c_Key1(2), common_term_98);
    // q = 3
    result_by_lpq.set(c_Key3(9, 2, 3), common_term_99);
    result_by_q.set(c_Key1(3), common_term_99);
    // q = 4
    result_by_lpq.set(c_Key3(9, 2, 4), common_term_100);
    result_by_q.set(c_Key1(4), common_term_100);
    // q = 5
    result_by_lpq.set(c_Key3(9, 2, 5), common_term_101);
    result_by_q.set(c_Key1(5), common_term_101);
    // q = 6
    result_by_lpq.set(c_Key3(9, 2, 6), common_term_102);
    result_by_q.set(c_Key1(6), common_term_102);
    // q = 7
    result_by_lpq.set(c_Key3(9, 2, 7), common_term_103);
    result_by_q.set(c_Key1(7), common_term_103);
    // q = 8
    result_by_lpq.set(c_Key3(9, 2, 8), common_term_104);
    result_by_q.set(c_Key1(8), common_term_104);
    // q = 9
    result_by_lpq.set(c_Key3(9, 2, 9), common_term_105);
    result_by_q.set(c_Key1(9), common_term_105);
    // q = 10
    result_by_lpq.set(c_Key3(9, 2, 10), common_term_106);
    result_by_q.set(c_Key1(10), common_term_106);
    // q = 11
    result_by_lpq.set(c_Key3(9, 2, 11), common_term_107);
    result_by_q.set(c_Key1(11), common_term_107);
    // q = 12
    result_by_lpq.set(c_Key3(9, 2, 12), common_term_108);
    result_by_q.set(c_Key1(12), common_term_108);
    // q = 13
    result_by_lpq.set(c_Key3(9, 2, 13), common_term_109);
    result_by_q.set(c_Key1(13), common_term_109);
    // q = 14
    result_by_lpq.set(c_Key3(9, 2, 14), common_term_110);
    result_by_q.set(c_Key1(14), common_term_110);
    // q = 15
    result_by_lpq.set(c_Key3(9, 2, 15), common_term_111);
    result_by_q.set(c_Key1(15), common_term_111);
    // q = 16
    result_by_lpq.set(c_Key3(9, 2, 16), common_term_112);
    result_by_q.set(c_Key1(16), common_term_112);
    // q = 17
    result_by_lpq.set(c_Key3(9, 2, 17), common_term_113);
    result_by_q.set(c_Key1(17), common_term_113);
    // q = 18
    result_by_lpq.set(c_Key3(9, 2, 18), common_term_114);
    result_by_q.set(c_Key1(18), common_term_114);
    // q = 19
    result_by_lpq.set(c_Key3(9, 2, 19), common_term_115);
    result_by_q.set(c_Key1(19), common_term_115);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 2), result_by_q);
    result_by_q.clear();

    // l , p = (9, 3).
    // q = -19
    result_by_lpq.set(c_Key3(9, 3, -19), common_term_116);
    result_by_q.set(c_Key1(-19), common_term_116);
    // q = -18
    result_by_lpq.set(c_Key3(9, 3, -18), common_term_117);
    result_by_q.set(c_Key1(-18), common_term_117);
    // q = -17
    result_by_lpq.set(c_Key3(9, 3, -17), common_term_118);
    result_by_q.set(c_Key1(-17), common_term_118);
    // q = -16
    result_by_lpq.set(c_Key3(9, 3, -16), common_term_119);
    result_by_q.set(c_Key1(-16), common_term_119);
    // q = -15
    result_by_lpq.set(c_Key3(9, 3, -15), common_term_120);
    result_by_q.set(c_Key1(-15), common_term_120);
    // q = -14
    result_by_lpq.set(c_Key3(9, 3, -14), common_term_121);
    result_by_q.set(c_Key1(-14), common_term_121);
    // q = -13
    result_by_lpq.set(c_Key3(9, 3, -13), common_term_122);
    result_by_q.set(c_Key1(-13), common_term_122);
    // q = -12
    result_by_lpq.set(c_Key3(9, 3, -12), common_term_123);
    result_by_q.set(c_Key1(-12), common_term_123);
    // q = -11
    result_by_lpq.set(c_Key3(9, 3, -11), common_term_124);
    result_by_q.set(c_Key1(-11), common_term_124);
    // q = -10
    result_by_lpq.set(c_Key3(9, 3, -10), common_term_125);
    result_by_q.set(c_Key1(-10), common_term_125);
    // q = -9
    result_by_lpq.set(c_Key3(9, 3, -9), common_term_126);
    result_by_q.set(c_Key1(-9), common_term_126);
    // q = -8
    result_by_lpq.set(c_Key3(9, 3, -8), common_term_127);
    result_by_q.set(c_Key1(-8), common_term_127);
    // q = -7
    result_by_lpq.set(c_Key3(9, 3, -7), common_term_128);
    result_by_q.set(c_Key1(-7), common_term_128);
    // q = -6
    result_by_lpq.set(c_Key3(9, 3, -6), common_term_129);
    result_by_q.set(c_Key1(-6), common_term_129);
    // q = -5
    result_by_lpq.set(c_Key3(9, 3, -5), common_term_130);
    result_by_q.set(c_Key1(-5), common_term_130);
    // q = -4
    result_by_lpq.set(c_Key3(9, 3, -4), common_term_131);
    result_by_q.set(c_Key1(-4), common_term_131);
    // q = -3
    result_by_lpq.set(c_Key3(9, 3, -3), common_term_132);
    result_by_q.set(c_Key1(-3), common_term_132);
    // q = -2
    result_by_lpq.set(c_Key3(9, 3, -2), common_term_133);
    result_by_q.set(c_Key1(-2), common_term_133);
    // q = -1
    result_by_lpq.set(c_Key3(9, 3, -1), common_term_134);
    result_by_q.set(c_Key1(-1), common_term_134);
    // q = 0
    result_by_lpq.set(c_Key3(9, 3, 0), common_term_135);
    result_by_q.set(c_Key1(0), common_term_135);
    // q = 1
    result_by_lpq.set(c_Key3(9, 3, 1), common_term_136);
    result_by_q.set(c_Key1(1), common_term_136);
    // q = 2
    result_by_lpq.set(c_Key3(9, 3, 2), common_term_137);
    result_by_q.set(c_Key1(2), common_term_137);
    // q = 3
    result_by_lpq.set(c_Key3(9, 3, 3), common_term_138);
    result_by_q.set(c_Key1(3), common_term_138);
    // q = 4
    result_by_lpq.set(c_Key3(9, 3, 4), common_term_139);
    result_by_q.set(c_Key1(4), common_term_139);
    // q = 5
    result_by_lpq.set(c_Key3(9, 3, 5), common_term_140);
    result_by_q.set(c_Key1(5), common_term_140);
    // q = 6
    result_by_lpq.set(c_Key3(9, 3, 6), common_term_141);
    result_by_q.set(c_Key1(6), common_term_141);
    // q = 7
    result_by_lpq.set(c_Key3(9, 3, 7), common_term_142);
    result_by_q.set(c_Key1(7), common_term_142);
    // q = 8
    result_by_lpq.set(c_Key3(9, 3, 8), common_term_143);
    result_by_q.set(c_Key1(8), common_term_143);
    // q = 9
    result_by_lpq.set(c_Key3(9, 3, 9), common_term_144);
    result_by_q.set(c_Key1(9), common_term_144);
    // q = 10
    result_by_lpq.set(c_Key3(9, 3, 10), common_term_145);
    result_by_q.set(c_Key1(10), common_term_145);
    // q = 11
    result_by_lpq.set(c_Key3(9, 3, 11), common_term_146);
    result_by_q.set(c_Key1(11), common_term_146);
    // q = 12
    result_by_lpq.set(c_Key3(9, 3, 12), common_term_147);
    result_by_q.set(c_Key1(12), common_term_147);
    // q = 13
    result_by_lpq.set(c_Key3(9, 3, 13), common_term_148);
    result_by_q.set(c_Key1(13), common_term_148);
    // q = 14
    result_by_lpq.set(c_Key3(9, 3, 14), common_term_149);
    result_by_q.set(c_Key1(14), common_term_149);
    // q = 15
    result_by_lpq.set(c_Key3(9, 3, 15), common_term_150);
    result_by_q.set(c_Key1(15), common_term_150);
    // q = 16
    result_by_lpq.set(c_Key3(9, 3, 16), common_term_151);
    result_by_q.set(c_Key1(16), common_term_151);
    // q = 17
    result_by_lpq.set(c_Key3(9, 3, 17), common_term_152);
    result_by_q.set(c_Key1(17), common_term_152);
    // q = 18
    result_by_lpq.set(c_Key3(9, 3, 18), common_term_153);
    result_by_q.set(c_Key1(18), common_term_153);
    // q = 19
    result_by_lpq.set(c_Key3(9, 3, 19), common_term_154);
    result_by_q.set(c_Key1(19), common_term_154);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 3), result_by_q);
    result_by_q.clear();

    // l , p = (9, 4).
    // q = -19
    result_by_lpq.set(c_Key3(9, 4, -19), common_term_155);
    result_by_q.set(c_Key1(-19), common_term_155);
    // q = -18
    result_by_lpq.set(c_Key3(9, 4, -18), common_term_156);
    result_by_q.set(c_Key1(-18), common_term_156);
    // q = -17
    result_by_lpq.set(c_Key3(9, 4, -17), common_term_157);
    result_by_q.set(c_Key1(-17), common_term_157);
    // q = -16
    result_by_lpq.set(c_Key3(9, 4, -16), common_term_158);
    result_by_q.set(c_Key1(-16), common_term_158);
    // q = -15
    result_by_lpq.set(c_Key3(9, 4, -15), common_term_159);
    result_by_q.set(c_Key1(-15), common_term_159);
    // q = -14
    result_by_lpq.set(c_Key3(9, 4, -14), common_term_160);
    result_by_q.set(c_Key1(-14), common_term_160);
    // q = -13
    result_by_lpq.set(c_Key3(9, 4, -13), common_term_161);
    result_by_q.set(c_Key1(-13), common_term_161);
    // q = -12
    result_by_lpq.set(c_Key3(9, 4, -12), common_term_162);
    result_by_q.set(c_Key1(-12), common_term_162);
    // q = -11
    result_by_lpq.set(c_Key3(9, 4, -11), common_term_163);
    result_by_q.set(c_Key1(-11), common_term_163);
    // q = -10
    result_by_lpq.set(c_Key3(9, 4, -10), common_term_164);
    result_by_q.set(c_Key1(-10), common_term_164);
    // q = -9
    result_by_lpq.set(c_Key3(9, 4, -9), common_term_165);
    result_by_q.set(c_Key1(-9), common_term_165);
    // q = -8
    result_by_lpq.set(c_Key3(9, 4, -8), common_term_166);
    result_by_q.set(c_Key1(-8), common_term_166);
    // q = -7
    result_by_lpq.set(c_Key3(9, 4, -7), common_term_167);
    result_by_q.set(c_Key1(-7), common_term_167);
    // q = -6
    result_by_lpq.set(c_Key3(9, 4, -6), common_term_168);
    result_by_q.set(c_Key1(-6), common_term_168);
    // q = -5
    result_by_lpq.set(c_Key3(9, 4, -5), common_term_169);
    result_by_q.set(c_Key1(-5), common_term_169);
    // q = -4
    result_by_lpq.set(c_Key3(9, 4, -4), common_term_170);
    result_by_q.set(c_Key1(-4), common_term_170);
    // q = -3
    result_by_lpq.set(c_Key3(9, 4, -3), common_term_171);
    result_by_q.set(c_Key1(-3), common_term_171);
    // q = -2
    result_by_lpq.set(c_Key3(9, 4, -2), common_term_172);
    result_by_q.set(c_Key1(-2), common_term_172);
    // q = -1
    result_by_lpq.set(c_Key3(9, 4, -1), common_term_173);
    result_by_q.set(c_Key1(-1), common_term_173);
    // q = 0
    result_by_lpq.set(c_Key3(9, 4, 0), common_term_174);
    result_by_q.set(c_Key1(0), common_term_174);
    // q = 1
    result_by_lpq.set(c_Key3(9, 4, 1), common_term_175);
    result_by_q.set(c_Key1(1), common_term_175);
    // q = 2
    result_by_lpq.set(c_Key3(9, 4, 2), common_term_176);
    result_by_q.set(c_Key1(2), common_term_176);
    // q = 3
    result_by_lpq.set(c_Key3(9, 4, 3), common_term_177);
    result_by_q.set(c_Key1(3), common_term_177);
    // q = 4
    result_by_lpq.set(c_Key3(9, 4, 4), common_term_178);
    result_by_q.set(c_Key1(4), common_term_178);
    // q = 5
    result_by_lpq.set(c_Key3(9, 4, 5), common_term_179);
    result_by_q.set(c_Key1(5), common_term_179);
    // q = 6
    result_by_lpq.set(c_Key3(9, 4, 6), common_term_180);
    result_by_q.set(c_Key1(6), common_term_180);
    // q = 7
    result_by_lpq.set(c_Key3(9, 4, 7), common_term_181);
    result_by_q.set(c_Key1(7), common_term_181);
    // q = 8
    result_by_lpq.set(c_Key3(9, 4, 8), common_term_182);
    result_by_q.set(c_Key1(8), common_term_182);
    // q = 9
    result_by_lpq.set(c_Key3(9, 4, 9), common_term_183);
    result_by_q.set(c_Key1(9), common_term_183);
    // q = 10
    result_by_lpq.set(c_Key3(9, 4, 10), common_term_184);
    result_by_q.set(c_Key1(10), common_term_184);
    // q = 11
    result_by_lpq.set(c_Key3(9, 4, 11), common_term_185);
    result_by_q.set(c_Key1(11), common_term_185);
    // q = 12
    result_by_lpq.set(c_Key3(9, 4, 12), common_term_186);
    result_by_q.set(c_Key1(12), common_term_186);
    // q = 13
    result_by_lpq.set(c_Key3(9, 4, 13), common_term_187);
    result_by_q.set(c_Key1(13), common_term_187);
    // q = 14
    result_by_lpq.set(c_Key3(9, 4, 14), common_term_188);
    result_by_q.set(c_Key1(14), common_term_188);
    // q = 15
    result_by_lpq.set(c_Key3(9, 4, 15), common_term_189);
    result_by_q.set(c_Key1(15), common_term_189);
    // q = 16
    result_by_lpq.set(c_Key3(9, 4, 16), common_term_190);
    result_by_q.set(c_Key1(16), common_term_190);
    // q = 17
    result_by_lpq.set(c_Key3(9, 4, 17), common_term_191);
    result_by_q.set(c_Key1(17), common_term_191);
    // q = 18
    result_by_lpq.set(c_Key3(9, 4, 18), common_term_192);
    result_by_q.set(c_Key1(18), common_term_192);
    // q = 19
    result_by_lpq.set(c_Key3(9, 4, 19), common_term_193);
    result_by_q.set(c_Key1(19), common_term_193);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 4), result_by_q);
    result_by_q.clear();

    // l , p = (9, 5).
    // q = -19
    result_by_lpq.set(c_Key3(9, 5, -19), common_term_193);
    result_by_q.set(c_Key1(-19), common_term_193);
    // q = -18
    result_by_lpq.set(c_Key3(9, 5, -18), common_term_192);
    result_by_q.set(c_Key1(-18), common_term_192);
    // q = -17
    result_by_lpq.set(c_Key3(9, 5, -17), common_term_191);
    result_by_q.set(c_Key1(-17), common_term_191);
    // q = -16
    result_by_lpq.set(c_Key3(9, 5, -16), common_term_190);
    result_by_q.set(c_Key1(-16), common_term_190);
    // q = -15
    result_by_lpq.set(c_Key3(9, 5, -15), common_term_189);
    result_by_q.set(c_Key1(-15), common_term_189);
    // q = -14
    result_by_lpq.set(c_Key3(9, 5, -14), common_term_188);
    result_by_q.set(c_Key1(-14), common_term_188);
    // q = -13
    result_by_lpq.set(c_Key3(9, 5, -13), common_term_187);
    result_by_q.set(c_Key1(-13), common_term_187);
    // q = -12
    result_by_lpq.set(c_Key3(9, 5, -12), common_term_186);
    result_by_q.set(c_Key1(-12), common_term_186);
    // q = -11
    result_by_lpq.set(c_Key3(9, 5, -11), common_term_185);
    result_by_q.set(c_Key1(-11), common_term_185);
    // q = -10
    result_by_lpq.set(c_Key3(9, 5, -10), common_term_184);
    result_by_q.set(c_Key1(-10), common_term_184);
    // q = -9
    result_by_lpq.set(c_Key3(9, 5, -9), common_term_183);
    result_by_q.set(c_Key1(-9), common_term_183);
    // q = -8
    result_by_lpq.set(c_Key3(9, 5, -8), common_term_182);
    result_by_q.set(c_Key1(-8), common_term_182);
    // q = -7
    result_by_lpq.set(c_Key3(9, 5, -7), common_term_181);
    result_by_q.set(c_Key1(-7), common_term_181);
    // q = -6
    result_by_lpq.set(c_Key3(9, 5, -6), common_term_180);
    result_by_q.set(c_Key1(-6), common_term_180);
    // q = -5
    result_by_lpq.set(c_Key3(9, 5, -5), common_term_179);
    result_by_q.set(c_Key1(-5), common_term_179);
    // q = -4
    result_by_lpq.set(c_Key3(9, 5, -4), common_term_178);
    result_by_q.set(c_Key1(-4), common_term_178);
    // q = -3
    result_by_lpq.set(c_Key3(9, 5, -3), common_term_177);
    result_by_q.set(c_Key1(-3), common_term_177);
    // q = -2
    result_by_lpq.set(c_Key3(9, 5, -2), common_term_176);
    result_by_q.set(c_Key1(-2), common_term_176);
    // q = -1
    result_by_lpq.set(c_Key3(9, 5, -1), common_term_175);
    result_by_q.set(c_Key1(-1), common_term_175);
    // q = 0
    result_by_lpq.set(c_Key3(9, 5, 0), common_term_174);
    result_by_q.set(c_Key1(0), common_term_174);
    // q = 1
    result_by_lpq.set(c_Key3(9, 5, 1), common_term_173);
    result_by_q.set(c_Key1(1), common_term_173);
    // q = 2
    result_by_lpq.set(c_Key3(9, 5, 2), common_term_172);
    result_by_q.set(c_Key1(2), common_term_172);
    // q = 3
    result_by_lpq.set(c_Key3(9, 5, 3), common_term_171);
    result_by_q.set(c_Key1(3), common_term_171);
    // q = 4
    result_by_lpq.set(c_Key3(9, 5, 4), common_term_170);
    result_by_q.set(c_Key1(4), common_term_170);
    // q = 5
    result_by_lpq.set(c_Key3(9, 5, 5), common_term_169);
    result_by_q.set(c_Key1(5), common_term_169);
    // q = 6
    result_by_lpq.set(c_Key3(9, 5, 6), common_term_168);
    result_by_q.set(c_Key1(6), common_term_168);
    // q = 7
    result_by_lpq.set(c_Key3(9, 5, 7), common_term_167);
    result_by_q.set(c_Key1(7), common_term_167);
    // q = 8
    result_by_lpq.set(c_Key3(9, 5, 8), common_term_166);
    result_by_q.set(c_Key1(8), common_term_166);
    // q = 9
    result_by_lpq.set(c_Key3(9, 5, 9), common_term_165);
    result_by_q.set(c_Key1(9), common_term_165);
    // q = 10
    result_by_lpq.set(c_Key3(9, 5, 10), common_term_164);
    result_by_q.set(c_Key1(10), common_term_164);
    // q = 11
    result_by_lpq.set(c_Key3(9, 5, 11), common_term_163);
    result_by_q.set(c_Key1(11), common_term_163);
    // q = 12
    result_by_lpq.set(c_Key3(9, 5, 12), common_term_162);
    result_by_q.set(c_Key1(12), common_term_162);
    // q = 13
    result_by_lpq.set(c_Key3(9, 5, 13), common_term_161);
    result_by_q.set(c_Key1(13), common_term_161);
    // q = 14
    result_by_lpq.set(c_Key3(9, 5, 14), common_term_160);
    result_by_q.set(c_Key1(14), common_term_160);
    // q = 15
    result_by_lpq.set(c_Key3(9, 5, 15), common_term_159);
    result_by_q.set(c_Key1(15), common_term_159);
    // q = 16
    result_by_lpq.set(c_Key3(9, 5, 16), common_term_158);
    result_by_q.set(c_Key1(16), common_term_158);
    // q = 17
    result_by_lpq.set(c_Key3(9, 5, 17), common_term_157);
    result_by_q.set(c_Key1(17), common_term_157);
    // q = 18
    result_by_lpq.set(c_Key3(9, 5, 18), common_term_156);
    result_by_q.set(c_Key1(18), common_term_156);
    // q = 19
    result_by_lpq.set(c_Key3(9, 5, 19), common_term_155);
    result_by_q.set(c_Key1(19), common_term_155);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 5), result_by_q);
    result_by_q.clear();

    // l , p = (9, 6).
    // q = -19
    result_by_lpq.set(c_Key3(9, 6, -19), common_term_154);
    result_by_q.set(c_Key1(-19), common_term_154);
    // q = -18
    result_by_lpq.set(c_Key3(9, 6, -18), common_term_153);
    result_by_q.set(c_Key1(-18), common_term_153);
    // q = -17
    result_by_lpq.set(c_Key3(9, 6, -17), common_term_152);
    result_by_q.set(c_Key1(-17), common_term_152);
    // q = -16
    result_by_lpq.set(c_Key3(9, 6, -16), common_term_151);
    result_by_q.set(c_Key1(-16), common_term_151);
    // q = -15
    result_by_lpq.set(c_Key3(9, 6, -15), common_term_150);
    result_by_q.set(c_Key1(-15), common_term_150);
    // q = -14
    result_by_lpq.set(c_Key3(9, 6, -14), common_term_149);
    result_by_q.set(c_Key1(-14), common_term_149);
    // q = -13
    result_by_lpq.set(c_Key3(9, 6, -13), common_term_148);
    result_by_q.set(c_Key1(-13), common_term_148);
    // q = -12
    result_by_lpq.set(c_Key3(9, 6, -12), common_term_147);
    result_by_q.set(c_Key1(-12), common_term_147);
    // q = -11
    result_by_lpq.set(c_Key3(9, 6, -11), common_term_146);
    result_by_q.set(c_Key1(-11), common_term_146);
    // q = -10
    result_by_lpq.set(c_Key3(9, 6, -10), common_term_145);
    result_by_q.set(c_Key1(-10), common_term_145);
    // q = -9
    result_by_lpq.set(c_Key3(9, 6, -9), common_term_144);
    result_by_q.set(c_Key1(-9), common_term_144);
    // q = -8
    result_by_lpq.set(c_Key3(9, 6, -8), common_term_143);
    result_by_q.set(c_Key1(-8), common_term_143);
    // q = -7
    result_by_lpq.set(c_Key3(9, 6, -7), common_term_142);
    result_by_q.set(c_Key1(-7), common_term_142);
    // q = -6
    result_by_lpq.set(c_Key3(9, 6, -6), common_term_141);
    result_by_q.set(c_Key1(-6), common_term_141);
    // q = -5
    result_by_lpq.set(c_Key3(9, 6, -5), common_term_140);
    result_by_q.set(c_Key1(-5), common_term_140);
    // q = -4
    result_by_lpq.set(c_Key3(9, 6, -4), common_term_139);
    result_by_q.set(c_Key1(-4), common_term_139);
    // q = -3
    result_by_lpq.set(c_Key3(9, 6, -3), common_term_138);
    result_by_q.set(c_Key1(-3), common_term_138);
    // q = -2
    result_by_lpq.set(c_Key3(9, 6, -2), common_term_137);
    result_by_q.set(c_Key1(-2), common_term_137);
    // q = -1
    result_by_lpq.set(c_Key3(9, 6, -1), common_term_136);
    result_by_q.set(c_Key1(-1), common_term_136);
    // q = 0
    result_by_lpq.set(c_Key3(9, 6, 0), common_term_135);
    result_by_q.set(c_Key1(0), common_term_135);
    // q = 1
    result_by_lpq.set(c_Key3(9, 6, 1), common_term_134);
    result_by_q.set(c_Key1(1), common_term_134);
    // q = 2
    result_by_lpq.set(c_Key3(9, 6, 2), common_term_133);
    result_by_q.set(c_Key1(2), common_term_133);
    // q = 3
    result_by_lpq.set(c_Key3(9, 6, 3), common_term_132);
    result_by_q.set(c_Key1(3), common_term_132);
    // q = 4
    result_by_lpq.set(c_Key3(9, 6, 4), common_term_131);
    result_by_q.set(c_Key1(4), common_term_131);
    // q = 5
    result_by_lpq.set(c_Key3(9, 6, 5), common_term_130);
    result_by_q.set(c_Key1(5), common_term_130);
    // q = 6
    result_by_lpq.set(c_Key3(9, 6, 6), common_term_129);
    result_by_q.set(c_Key1(6), common_term_129);
    // q = 7
    result_by_lpq.set(c_Key3(9, 6, 7), common_term_128);
    result_by_q.set(c_Key1(7), common_term_128);
    // q = 8
    result_by_lpq.set(c_Key3(9, 6, 8), common_term_127);
    result_by_q.set(c_Key1(8), common_term_127);
    // q = 9
    result_by_lpq.set(c_Key3(9, 6, 9), common_term_126);
    result_by_q.set(c_Key1(9), common_term_126);
    // q = 10
    result_by_lpq.set(c_Key3(9, 6, 10), common_term_125);
    result_by_q.set(c_Key1(10), common_term_125);
    // q = 11
    result_by_lpq.set(c_Key3(9, 6, 11), common_term_124);
    result_by_q.set(c_Key1(11), common_term_124);
    // q = 12
    result_by_lpq.set(c_Key3(9, 6, 12), common_term_123);
    result_by_q.set(c_Key1(12), common_term_123);
    // q = 13
    result_by_lpq.set(c_Key3(9, 6, 13), common_term_122);
    result_by_q.set(c_Key1(13), common_term_122);
    // q = 14
    result_by_lpq.set(c_Key3(9, 6, 14), common_term_121);
    result_by_q.set(c_Key1(14), common_term_121);
    // q = 15
    result_by_lpq.set(c_Key3(9, 6, 15), common_term_120);
    result_by_q.set(c_Key1(15), common_term_120);
    // q = 16
    result_by_lpq.set(c_Key3(9, 6, 16), common_term_119);
    result_by_q.set(c_Key1(16), common_term_119);
    // q = 17
    result_by_lpq.set(c_Key3(9, 6, 17), common_term_118);
    result_by_q.set(c_Key1(17), common_term_118);
    // q = 18
    result_by_lpq.set(c_Key3(9, 6, 18), common_term_117);
    result_by_q.set(c_Key1(18), common_term_117);
    // q = 19
    result_by_lpq.set(c_Key3(9, 6, 19), common_term_116);
    result_by_q.set(c_Key1(19), common_term_116);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 6), result_by_q);
    result_by_q.clear();

    // l , p = (9, 7).
    // q = -19
    result_by_lpq.set(c_Key3(9, 7, -19), common_term_115);
    result_by_q.set(c_Key1(-19), common_term_115);
    // q = -18
    result_by_lpq.set(c_Key3(9, 7, -18), common_term_114);
    result_by_q.set(c_Key1(-18), common_term_114);
    // q = -17
    result_by_lpq.set(c_Key3(9, 7, -17), common_term_113);
    result_by_q.set(c_Key1(-17), common_term_113);
    // q = -16
    result_by_lpq.set(c_Key3(9, 7, -16), common_term_112);
    result_by_q.set(c_Key1(-16), common_term_112);
    // q = -15
    result_by_lpq.set(c_Key3(9, 7, -15), common_term_111);
    result_by_q.set(c_Key1(-15), common_term_111);
    // q = -14
    result_by_lpq.set(c_Key3(9, 7, -14), common_term_110);
    result_by_q.set(c_Key1(-14), common_term_110);
    // q = -13
    result_by_lpq.set(c_Key3(9, 7, -13), common_term_109);
    result_by_q.set(c_Key1(-13), common_term_109);
    // q = -12
    result_by_lpq.set(c_Key3(9, 7, -12), common_term_108);
    result_by_q.set(c_Key1(-12), common_term_108);
    // q = -11
    result_by_lpq.set(c_Key3(9, 7, -11), common_term_107);
    result_by_q.set(c_Key1(-11), common_term_107);
    // q = -10
    result_by_lpq.set(c_Key3(9, 7, -10), common_term_106);
    result_by_q.set(c_Key1(-10), common_term_106);
    // q = -9
    result_by_lpq.set(c_Key3(9, 7, -9), common_term_105);
    result_by_q.set(c_Key1(-9), common_term_105);
    // q = -8
    result_by_lpq.set(c_Key3(9, 7, -8), common_term_104);
    result_by_q.set(c_Key1(-8), common_term_104);
    // q = -7
    result_by_lpq.set(c_Key3(9, 7, -7), common_term_103);
    result_by_q.set(c_Key1(-7), common_term_103);
    // q = -6
    result_by_lpq.set(c_Key3(9, 7, -6), common_term_102);
    result_by_q.set(c_Key1(-6), common_term_102);
    // q = -5
    result_by_lpq.set(c_Key3(9, 7, -5), common_term_101);
    result_by_q.set(c_Key1(-5), common_term_101);
    // q = -4
    result_by_lpq.set(c_Key3(9, 7, -4), common_term_100);
    result_by_q.set(c_Key1(-4), common_term_100);
    // q = -3
    result_by_lpq.set(c_Key3(9, 7, -3), common_term_99);
    result_by_q.set(c_Key1(-3), common_term_99);
    // q = -2
    result_by_lpq.set(c_Key3(9, 7, -2), common_term_98);
    result_by_q.set(c_Key1(-2), common_term_98);
    // q = -1
    result_by_lpq.set(c_Key3(9, 7, -1), common_term_97);
    result_by_q.set(c_Key1(-1), common_term_97);
    // q = 0
    result_by_lpq.set(c_Key3(9, 7, 0), common_term_96);
    result_by_q.set(c_Key1(0), common_term_96);
    // q = 1
    result_by_lpq.set(c_Key3(9, 7, 1), common_term_95);
    result_by_q.set(c_Key1(1), common_term_95);
    // q = 2
    result_by_lpq.set(c_Key3(9, 7, 2), common_term_94);
    result_by_q.set(c_Key1(2), common_term_94);
    // q = 3
    result_by_lpq.set(c_Key3(9, 7, 3), common_term_93);
    result_by_q.set(c_Key1(3), common_term_93);
    // q = 4
    result_by_lpq.set(c_Key3(9, 7, 4), common_term_92);
    result_by_q.set(c_Key1(4), common_term_92);
    // q = 5
    result_by_lpq.set(c_Key3(9, 7, 5), common_term_91);
    result_by_q.set(c_Key1(5), common_term_91);
    // q = 6
    result_by_lpq.set(c_Key3(9, 7, 6), common_term_90);
    result_by_q.set(c_Key1(6), common_term_90);
    // q = 7
    result_by_lpq.set(c_Key3(9, 7, 7), common_term_89);
    result_by_q.set(c_Key1(7), common_term_89);
    // q = 8
    result_by_lpq.set(c_Key3(9, 7, 8), common_term_88);
    result_by_q.set(c_Key1(8), common_term_88);
    // q = 9
    result_by_lpq.set(c_Key3(9, 7, 9), common_term_87);
    result_by_q.set(c_Key1(9), common_term_87);
    // q = 10
    result_by_lpq.set(c_Key3(9, 7, 10), common_term_86);
    result_by_q.set(c_Key1(10), common_term_86);
    // q = 11
    result_by_lpq.set(c_Key3(9, 7, 11), common_term_85);
    result_by_q.set(c_Key1(11), common_term_85);
    // q = 12
    result_by_lpq.set(c_Key3(9, 7, 12), common_term_84);
    result_by_q.set(c_Key1(12), common_term_84);
    // q = 13
    result_by_lpq.set(c_Key3(9, 7, 13), common_term_83);
    result_by_q.set(c_Key1(13), common_term_83);
    // q = 14
    result_by_lpq.set(c_Key3(9, 7, 14), common_term_82);
    result_by_q.set(c_Key1(14), common_term_82);
    // q = 15
    result_by_lpq.set(c_Key3(9, 7, 15), common_term_81);
    result_by_q.set(c_Key1(15), common_term_81);
    // q = 16
    result_by_lpq.set(c_Key3(9, 7, 16), common_term_80);
    result_by_q.set(c_Key1(16), common_term_80);
    // q = 17
    result_by_lpq.set(c_Key3(9, 7, 17), common_term_79);
    result_by_q.set(c_Key1(17), common_term_79);
    // q = 18
    result_by_lpq.set(c_Key3(9, 7, 18), common_term_78);
    result_by_q.set(c_Key1(18), common_term_78);
    // q = 19
    result_by_lpq.set(c_Key3(9, 7, 19), common_term_77);
    result_by_q.set(c_Key1(19), common_term_77);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 7), result_by_q);
    result_by_q.clear();

    // l , p = (9, 8).
    // q = -19
    result_by_lpq.set(c_Key3(9, 8, -19), common_term_76);
    result_by_q.set(c_Key1(-19), common_term_76);
    // q = -18
    result_by_lpq.set(c_Key3(9, 8, -18), common_term_75);
    result_by_q.set(c_Key1(-18), common_term_75);
    // q = -17
    result_by_lpq.set(c_Key3(9, 8, -17), common_term_74);
    result_by_q.set(c_Key1(-17), common_term_74);
    // q = -16
    result_by_lpq.set(c_Key3(9, 8, -16), common_term_73);
    result_by_q.set(c_Key1(-16), common_term_73);
    // q = -15
    result_by_lpq.set(c_Key3(9, 8, -15), common_term_72);
    result_by_q.set(c_Key1(-15), common_term_72);
    // q = -14
    result_by_lpq.set(c_Key3(9, 8, -14), common_term_71);
    result_by_q.set(c_Key1(-14), common_term_71);
    // q = -13
    result_by_lpq.set(c_Key3(9, 8, -13), common_term_70);
    result_by_q.set(c_Key1(-13), common_term_70);
    // q = -12
    result_by_lpq.set(c_Key3(9, 8, -12), common_term_69);
    result_by_q.set(c_Key1(-12), common_term_69);
    // q = -11
    result_by_lpq.set(c_Key3(9, 8, -11), common_term_68);
    result_by_q.set(c_Key1(-11), common_term_68);
    // q = -10
    result_by_lpq.set(c_Key3(9, 8, -10), common_term_67);
    result_by_q.set(c_Key1(-10), common_term_67);
    // q = -9
    result_by_lpq.set(c_Key3(9, 8, -9), common_term_66);
    result_by_q.set(c_Key1(-9), common_term_66);
    // q = -8
    result_by_lpq.set(c_Key3(9, 8, -8), common_term_65);
    result_by_q.set(c_Key1(-8), common_term_65);
    // q = -7
    result_by_lpq.set(c_Key3(9, 8, -7), common_term_64);
    result_by_q.set(c_Key1(-7), common_term_64);
    // q = -6
    result_by_lpq.set(c_Key3(9, 8, -6), common_term_63);
    result_by_q.set(c_Key1(-6), common_term_63);
    // q = -5
    result_by_lpq.set(c_Key3(9, 8, -5), common_term_62);
    result_by_q.set(c_Key1(-5), common_term_62);
    // q = -4
    result_by_lpq.set(c_Key3(9, 8, -4), common_term_61);
    result_by_q.set(c_Key1(-4), common_term_61);
    // q = -3
    result_by_lpq.set(c_Key3(9, 8, -3), common_term_60);
    result_by_q.set(c_Key1(-3), common_term_60);
    // q = -2
    result_by_lpq.set(c_Key3(9, 8, -2), common_term_59);
    result_by_q.set(c_Key1(-2), common_term_59);
    // q = -1
    result_by_lpq.set(c_Key3(9, 8, -1), common_term_58);
    result_by_q.set(c_Key1(-1), common_term_58);
    // q = 0
    result_by_lpq.set(c_Key3(9, 8, 0), common_term_57);
    result_by_q.set(c_Key1(0), common_term_57);
    // q = 1
    result_by_lpq.set(c_Key3(9, 8, 1), common_term_56);
    result_by_q.set(c_Key1(1), common_term_56);
    // q = 2
    result_by_lpq.set(c_Key3(9, 8, 2), common_term_55);
    result_by_q.set(c_Key1(2), common_term_55);
    // q = 3
    result_by_lpq.set(c_Key3(9, 8, 3), common_term_54);
    result_by_q.set(c_Key1(3), common_term_54);
    // q = 4
    result_by_lpq.set(c_Key3(9, 8, 4), common_term_53);
    result_by_q.set(c_Key1(4), common_term_53);
    // q = 5
    result_by_lpq.set(c_Key3(9, 8, 5), common_term_52);
    result_by_q.set(c_Key1(5), common_term_52);
    // q = 6
    result_by_lpq.set(c_Key3(9, 8, 6), common_term_51);
    result_by_q.set(c_Key1(6), common_term_51);
    // q = 7
    result_by_lpq.set(c_Key3(9, 8, 7), common_term_50);
    result_by_q.set(c_Key1(7), common_term_50);
    // q = 8
    result_by_lpq.set(c_Key3(9, 8, 8), common_term_49);
    result_by_q.set(c_Key1(8), common_term_49);
    // q = 9
    result_by_lpq.set(c_Key3(9, 8, 9), common_term_48);
    result_by_q.set(c_Key1(9), common_term_48);
    // q = 10
    result_by_lpq.set(c_Key3(9, 8, 10), common_term_47);
    result_by_q.set(c_Key1(10), common_term_47);
    // q = 11
    result_by_lpq.set(c_Key3(9, 8, 11), common_term_46);
    result_by_q.set(c_Key1(11), common_term_46);
    // q = 12
    result_by_lpq.set(c_Key3(9, 8, 12), common_term_45);
    result_by_q.set(c_Key1(12), common_term_45);
    // q = 13
    result_by_lpq.set(c_Key3(9, 8, 13), common_term_44);
    result_by_q.set(c_Key1(13), common_term_44);
    // q = 14
    result_by_lpq.set(c_Key3(9, 8, 14), common_term_43);
    result_by_q.set(c_Key1(14), common_term_43);
    // q = 15
    result_by_lpq.set(c_Key3(9, 8, 15), common_term_42);
    result_by_q.set(c_Key1(15), common_term_42);
    // q = 16
    result_by_lpq.set(c_Key3(9, 8, 16), common_term_41);
    result_by_q.set(c_Key1(16), common_term_41);
    // q = 17
    result_by_lpq.set(c_Key3(9, 8, 17), common_term_40);
    result_by_q.set(c_Key1(17), common_term_40);
    // q = 18
    result_by_lpq.set(c_Key3(9, 8, 18), common_term_39);
    result_by_q.set(c_Key1(18), common_term_39);
    // q = 19
    result_by_lpq.set(c_Key3(9, 8, 19), common_term_38);
    result_by_q.set(c_Key1(19), common_term_38);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 8), result_by_q);
    result_by_q.clear();

    // l , p = (9, 9).
    // q = -19
    result_by_lpq.set(c_Key3(9, 9, -19), common_term_37);
    result_by_q.set(c_Key1(-19), common_term_37);
    // q = -18
    result_by_lpq.set(c_Key3(9, 9, -18), common_term_36);
    result_by_q.set(c_Key1(-18), common_term_36);
    // q = -17
    result_by_lpq.set(c_Key3(9, 9, -17), common_term_35);
    result_by_q.set(c_Key1(-17), common_term_35);
    // q = -16
    result_by_lpq.set(c_Key3(9, 9, -16), common_term_34);
    result_by_q.set(c_Key1(-16), common_term_34);
    // q = -15
    result_by_lpq.set(c_Key3(9, 9, -15), common_term_33);
    result_by_q.set(c_Key1(-15), common_term_33);
    // q = -14
    result_by_lpq.set(c_Key3(9, 9, -14), common_term_32);
    result_by_q.set(c_Key1(-14), common_term_32);
    // q = -13
    result_by_lpq.set(c_Key3(9, 9, -13), common_term_31);
    result_by_q.set(c_Key1(-13), common_term_31);
    // q = -12
    result_by_lpq.set(c_Key3(9, 9, -12), common_term_30);
    result_by_q.set(c_Key1(-12), common_term_30);
    // q = -11
    result_by_lpq.set(c_Key3(9, 9, -11), common_term_29);
    result_by_q.set(c_Key1(-11), common_term_29);
    // q = -10
    result_by_lpq.set(c_Key3(9, 9, -10), common_term_28);
    result_by_q.set(c_Key1(-10), common_term_28);
    // q = -9
    result_by_lpq.set(c_Key3(9, 9, -9), common_term_27);
    result_by_q.set(c_Key1(-9), common_term_27);
    // q = -8
    result_by_lpq.set(c_Key3(9, 9, -8), common_term_26);
    result_by_q.set(c_Key1(-8), common_term_26);
    // q = -7
    result_by_lpq.set(c_Key3(9, 9, -7), common_term_25);
    result_by_q.set(c_Key1(-7), common_term_25);
    // q = -6
    result_by_lpq.set(c_Key3(9, 9, -6), common_term_24);
    result_by_q.set(c_Key1(-6), common_term_24);
    // q = -5
    result_by_lpq.set(c_Key3(9, 9, -5), common_term_23);
    result_by_q.set(c_Key1(-5), common_term_23);
    // q = -4
    result_by_lpq.set(c_Key3(9, 9, -4), common_term_22);
    result_by_q.set(c_Key1(-4), common_term_22);
    // q = -3
    result_by_lpq.set(c_Key3(9, 9, -3), common_term_21);
    result_by_q.set(c_Key1(-3), common_term_21);
    // q = -2
    result_by_lpq.set(c_Key3(9, 9, -2), common_term_20);
    result_by_q.set(c_Key1(-2), common_term_20);
    // q = -1
    result_by_lpq.set(c_Key3(9, 9, -1), common_term_19);
    result_by_q.set(c_Key1(-1), common_term_19);
    // q = 0
    result_by_lpq.set(c_Key3(9, 9, 0), common_term_18);
    result_by_q.set(c_Key1(0), common_term_18);
    // q = 1
    result_by_lpq.set(c_Key3(9, 9, 1), common_term_17);
    result_by_q.set(c_Key1(1), common_term_17);
    // q = 2
    result_by_lpq.set(c_Key3(9, 9, 2), common_term_16);
    result_by_q.set(c_Key1(2), common_term_16);
    // q = 3
    result_by_lpq.set(c_Key3(9, 9, 3), common_term_15);
    result_by_q.set(c_Key1(3), common_term_15);
    // q = 4
    result_by_lpq.set(c_Key3(9, 9, 4), common_term_14);
    result_by_q.set(c_Key1(4), common_term_14);
    // q = 5
    result_by_lpq.set(c_Key3(9, 9, 5), common_term_13);
    result_by_q.set(c_Key1(5), common_term_13);
    // q = 6
    result_by_lpq.set(c_Key3(9, 9, 6), common_term_12);
    result_by_q.set(c_Key1(6), common_term_12);
    // q = 7
    result_by_lpq.set(c_Key3(9, 9, 7), common_term_11);
    result_by_q.set(c_Key1(7), common_term_11);
    // q = 8
    result_by_lpq.set(c_Key3(9, 9, 8), common_term_10);
    result_by_q.set(c_Key1(8), common_term_10);
    // q = 10
    result_by_lpq.set(c_Key3(9, 9, 10), common_term_9);
    result_by_q.set(c_Key1(10), common_term_9);
    // q = 11
    result_by_lpq.set(c_Key3(9, 9, 11), common_term_8);
    result_by_q.set(c_Key1(11), common_term_8);
    // q = 12
    result_by_lpq.set(c_Key3(9, 9, 12), common_term_7);
    result_by_q.set(c_Key1(12), common_term_7);
    // q = 13
    result_by_lpq.set(c_Key3(9, 9, 13), common_term_6);
    result_by_q.set(c_Key1(13), common_term_6);
    // q = 14
    result_by_lpq.set(c_Key3(9, 9, 14), common_term_5);
    result_by_q.set(c_Key1(14), common_term_5);
    // q = 15
    result_by_lpq.set(c_Key3(9, 9, 15), common_term_4);
    result_by_q.set(c_Key1(15), common_term_4);
    // q = 16
    result_by_lpq.set(c_Key3(9, 9, 16), common_term_3);
    result_by_q.set(c_Key1(16), common_term_3);
    // q = 17
    result_by_lpq.set(c_Key3(9, 9, 17), common_term_2);
    result_by_q.set(c_Key1(17), common_term_2);
    // q = 18
    result_by_lpq.set(c_Key3(9, 9, 18), common_term_1);
    result_by_q.set(c_Key1(18), common_term_1);
    // q = 19
    result_by_lpq.set(c_Key3(9, 9, 19), common_term_0);
    result_by_q.set(c_Key1(19), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(9, 9), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}
