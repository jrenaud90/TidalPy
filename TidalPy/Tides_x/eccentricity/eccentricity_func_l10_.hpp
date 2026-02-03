#pragma once

#include "eccentricity_common_.hpp"

EccentricityFuncOutput c_eccentricity_function_l10_e1(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 10.
    // Functions are taylor expanded and truncated at eccentricity^1.

    // Functions Calculated for l = 10.

    c_IntMap<c_Key3, double> result_by_lpq(11);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(11);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;

    c_IntMap<c_Key1, double> result_by_q(1);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (10, 0).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 0, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 0), result_by_q);
    result_by_q.clear();

    // l , p = (10, 1).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 1), result_by_q);
    result_by_q.clear();

    // l , p = (10, 2).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 2), result_by_q);
    result_by_q.clear();

    // l , p = (10, 3).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 3), result_by_q);
    result_by_q.clear();

    // l , p = (10, 4).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 4, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 4), result_by_q);
    result_by_q.clear();

    // l , p = (10, 5).
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -9.5);
    result_by_lpq.set(c_Key3(10, 5, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 5), result_by_q);
    result_by_q.clear();

    // l , p = (10, 6).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 6, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 6), result_by_q);
    result_by_q.clear();

    // l , p = (10, 7).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 7, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 7), result_by_q);
    result_by_q.clear();

    // l , p = (10, 8).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 8, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 8), result_by_q);
    result_by_q.clear();

    // l , p = (10, 9).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 9, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 9), result_by_q);
    result_by_q.clear();

    // l , p = (10, 10).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 10, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 10), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l10_e2(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 10.
    // Functions are taylor expanded and truncated at eccentricity^2.

    // Functions Calculated for l = 10.

    c_IntMap<c_Key3, double> result_by_lpq(35);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(11);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = -4.5*eccentricity;
    double common_term_1 = 15.5*eccentricity;
    double common_term_2 = -2.5*eccentricity;
    double common_term_3 = 13.5*eccentricity;
    double common_term_4 = -0.5*eccentricity;
    double common_term_5 = 11.5*eccentricity;
    double common_term_6 = 1.5*eccentricity;
    double common_term_7 = 9.5*eccentricity;
    double common_term_8 = 9.0*eccentricity_2*std::pow(1.0 - eccentricity_2, -9.5);
    double common_term_9 = 3.5*eccentricity;
    double common_term_10 = 7.5*eccentricity;
    double common_term_11 = 5.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(4);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (10, 0).
    // q = -1
    result_by_lpq.set(c_Key3(10, 0, -1), common_term_0);
    result_by_q.set(c_Key1(-1), common_term_0);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 0, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(10, 0, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 0), result_by_q);
    result_by_q.clear();

    // l , p = (10, 1).
    // q = -1
    result_by_lpq.set(c_Key3(10, 1, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(10, 1, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 1), result_by_q);
    result_by_q.clear();

    // l , p = (10, 2).
    // q = -1
    result_by_lpq.set(c_Key3(10, 2, -1), common_term_4);
    result_by_q.set(c_Key1(-1), common_term_4);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(10, 2, 1), common_term_5);
    result_by_q.set(c_Key1(1), common_term_5);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 2), result_by_q);
    result_by_q.clear();

    // l , p = (10, 3).
    // q = -1
    result_by_lpq.set(c_Key3(10, 3, -1), common_term_6);
    result_by_q.set(c_Key1(-1), common_term_6);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(10, 3, 1), common_term_7);
    result_by_q.set(c_Key1(1), common_term_7);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 3), result_by_q);
    result_by_q.clear();

    // l , p = (10, 4).
    // q = -2
    result_by_lpq.set(c_Key3(10, 4, -2), common_term_8);
    result_by_q.set(c_Key1(-2), common_term_8);
    // q = -1
    result_by_lpq.set(c_Key3(10, 4, -1), common_term_9);
    result_by_q.set(c_Key1(-1), common_term_9);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 4, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(10, 4, 1), common_term_10);
    result_by_q.set(c_Key1(1), common_term_10);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 4), result_by_q);
    result_by_q.clear();

    // l , p = (10, 5).
    // q = -1
    result_by_lpq.set(c_Key3(10, 5, -1), common_term_11);
    result_by_q.set(c_Key1(-1), common_term_11);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -9.5)*(18.0*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(10, 5, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(10, 5, 1), common_term_11);
    result_by_q.set(c_Key1(1), common_term_11);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 5), result_by_q);
    result_by_q.clear();

    // l , p = (10, 6).
    // q = -1
    result_by_lpq.set(c_Key3(10, 6, -1), common_term_10);
    result_by_q.set(c_Key1(-1), common_term_10);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 6, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(10, 6, 1), common_term_9);
    result_by_q.set(c_Key1(1), common_term_9);
    // q = 2
    result_by_lpq.set(c_Key3(10, 6, 2), common_term_8);
    result_by_q.set(c_Key1(2), common_term_8);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 6), result_by_q);
    result_by_q.clear();

    // l , p = (10, 7).
    // q = -1
    result_by_lpq.set(c_Key3(10, 7, -1), common_term_7);
    result_by_q.set(c_Key1(-1), common_term_7);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 7, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(10, 7, 1), common_term_6);
    result_by_q.set(c_Key1(1), common_term_6);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 7), result_by_q);
    result_by_q.clear();

    // l , p = (10, 8).
    // q = -1
    result_by_lpq.set(c_Key3(10, 8, -1), common_term_5);
    result_by_q.set(c_Key1(-1), common_term_5);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 8, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(10, 8, 1), common_term_4);
    result_by_q.set(c_Key1(1), common_term_4);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 8), result_by_q);
    result_by_q.clear();

    // l , p = (10, 9).
    // q = -1
    result_by_lpq.set(c_Key3(10, 9, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 9, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(10, 9, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 9), result_by_q);
    result_by_q.clear();

    // l , p = (10, 10).
    // q = -1
    result_by_lpq.set(c_Key3(10, 10, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(10, 10, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(10, 10, 1), common_term_0);
    result_by_q.set(c_Key1(1), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 10), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l10_e3(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 10.
    // Functions are taylor expanded and truncated at eccentricity^3.

    // Functions Calculated for l = 10.

    c_IntMap<c_Key3, double> result_by_lpq(55);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(11);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = 8.0*eccentricity_2;
    double common_term_1 = -4.5*eccentricity;
    double common_term_2 = 1.0 - 72.5*eccentricity_2;
    double common_term_3 = 15.5*eccentricity;
    double common_term_4 = 130.5*eccentricity_2;
    double common_term_5 = 2.25*eccentricity_2;
    double common_term_6 = -2.5*eccentricity;
    double common_term_7 = 1.0 - 36.5*eccentricity_2;
    double common_term_8 = 13.5*eccentricity;
    double common_term_9 = 100.25*eccentricity_2;
    double common_term_10 = 0.5*eccentricity_2;
    double common_term_11 = -0.5*eccentricity;
    double common_term_12 = 1.0 - 8.5*eccentricity_2;
    double common_term_13 = 11.5*eccentricity;
    double common_term_14 = 74.0*eccentricity_2;
    double common_term_15 = 2.75*eccentricity_2;
    double common_term_16 = 1.5*eccentricity;
    double common_term_17 = 11.5*eccentricity_2 + 1.0;
    double common_term_18 = 9.5*eccentricity;
    double common_term_19 = 51.75*eccentricity_2;
    double common_term_20 = 9.0*eccentricity_2*std::pow(1.0 - eccentricity_2, -9.5);
    double common_term_21 = 3.5*eccentricity;
    double common_term_22 = 23.5*eccentricity_2 + 1.0;
    double common_term_23 = 7.5*eccentricity;
    double common_term_24 = 33.5*eccentricity_2;
    double common_term_25 = 19.25*eccentricity_2;
    double common_term_26 = 5.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(5);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (10, 0).
    // q = -2
    result_by_lpq.set(c_Key3(10, 0, -2), common_term_0);
    result_by_q.set(c_Key1(-2), common_term_0);
    // q = -1
    result_by_lpq.set(c_Key3(10, 0, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    result_by_lpq.set(c_Key3(10, 0, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(10, 0, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // q = 2
    result_by_lpq.set(c_Key3(10, 0, 2), common_term_4);
    result_by_q.set(c_Key1(2), common_term_4);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 0), result_by_q);
    result_by_q.clear();

    // l , p = (10, 1).
    // q = -2
    result_by_lpq.set(c_Key3(10, 1, -2), common_term_5);
    result_by_q.set(c_Key1(-2), common_term_5);
    // q = -1
    result_by_lpq.set(c_Key3(10, 1, -1), common_term_6);
    result_by_q.set(c_Key1(-1), common_term_6);
    // q = 0
    result_by_lpq.set(c_Key3(10, 1, 0), common_term_7);
    result_by_q.set(c_Key1(0), common_term_7);
    // q = 1
    result_by_lpq.set(c_Key3(10, 1, 1), common_term_8);
    result_by_q.set(c_Key1(1), common_term_8);
    // q = 2
    result_by_lpq.set(c_Key3(10, 1, 2), common_term_9);
    result_by_q.set(c_Key1(2), common_term_9);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 1), result_by_q);
    result_by_q.clear();

    // l , p = (10, 2).
    // q = -2
    result_by_lpq.set(c_Key3(10, 2, -2), common_term_10);
    result_by_q.set(c_Key1(-2), common_term_10);
    // q = -1
    result_by_lpq.set(c_Key3(10, 2, -1), common_term_11);
    result_by_q.set(c_Key1(-1), common_term_11);
    // q = 0
    result_by_lpq.set(c_Key3(10, 2, 0), common_term_12);
    result_by_q.set(c_Key1(0), common_term_12);
    // q = 1
    result_by_lpq.set(c_Key3(10, 2, 1), common_term_13);
    result_by_q.set(c_Key1(1), common_term_13);
    // q = 2
    result_by_lpq.set(c_Key3(10, 2, 2), common_term_14);
    result_by_q.set(c_Key1(2), common_term_14);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 2), result_by_q);
    result_by_q.clear();

    // l , p = (10, 3).
    // q = -2
    result_by_lpq.set(c_Key3(10, 3, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(10, 3, -1), common_term_16);
    result_by_q.set(c_Key1(-1), common_term_16);
    // q = 0
    result_by_lpq.set(c_Key3(10, 3, 0), common_term_17);
    result_by_q.set(c_Key1(0), common_term_17);
    // q = 1
    result_by_lpq.set(c_Key3(10, 3, 1), common_term_18);
    result_by_q.set(c_Key1(1), common_term_18);
    // q = 2
    result_by_lpq.set(c_Key3(10, 3, 2), common_term_19);
    result_by_q.set(c_Key1(2), common_term_19);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 3), result_by_q);
    result_by_q.clear();

    // l , p = (10, 4).
    // q = -2
    result_by_lpq.set(c_Key3(10, 4, -2), common_term_20);
    result_by_q.set(c_Key1(-2), common_term_20);
    // q = -1
    result_by_lpq.set(c_Key3(10, 4, -1), common_term_21);
    result_by_q.set(c_Key1(-1), common_term_21);
    // q = 0
    result_by_lpq.set(c_Key3(10, 4, 0), common_term_22);
    result_by_q.set(c_Key1(0), common_term_22);
    // q = 1
    result_by_lpq.set(c_Key3(10, 4, 1), common_term_23);
    result_by_q.set(c_Key1(1), common_term_23);
    // q = 2
    result_by_lpq.set(c_Key3(10, 4, 2), common_term_24);
    result_by_q.set(c_Key1(2), common_term_24);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 4), result_by_q);
    result_by_q.clear();

    // l , p = (10, 5).
    // q = -2
    result_by_lpq.set(c_Key3(10, 5, -2), common_term_25);
    result_by_q.set(c_Key1(-2), common_term_25);
    // q = -1
    result_by_lpq.set(c_Key3(10, 5, -1), common_term_26);
    result_by_q.set(c_Key1(-1), common_term_26);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -9.5)*(18.0*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(10, 5, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(10, 5, 1), common_term_26);
    result_by_q.set(c_Key1(1), common_term_26);
    // q = 2
    result_by_lpq.set(c_Key3(10, 5, 2), common_term_25);
    result_by_q.set(c_Key1(2), common_term_25);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 5), result_by_q);
    result_by_q.clear();

    // l , p = (10, 6).
    // q = -2
    result_by_lpq.set(c_Key3(10, 6, -2), common_term_24);
    result_by_q.set(c_Key1(-2), common_term_24);
    // q = -1
    result_by_lpq.set(c_Key3(10, 6, -1), common_term_23);
    result_by_q.set(c_Key1(-1), common_term_23);
    // q = 0
    result_by_lpq.set(c_Key3(10, 6, 0), common_term_22);
    result_by_q.set(c_Key1(0), common_term_22);
    // q = 1
    result_by_lpq.set(c_Key3(10, 6, 1), common_term_21);
    result_by_q.set(c_Key1(1), common_term_21);
    // q = 2
    result_by_lpq.set(c_Key3(10, 6, 2), common_term_20);
    result_by_q.set(c_Key1(2), common_term_20);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 6), result_by_q);
    result_by_q.clear();

    // l , p = (10, 7).
    // q = -2
    result_by_lpq.set(c_Key3(10, 7, -2), common_term_19);
    result_by_q.set(c_Key1(-2), common_term_19);
    // q = -1
    result_by_lpq.set(c_Key3(10, 7, -1), common_term_18);
    result_by_q.set(c_Key1(-1), common_term_18);
    // q = 0
    result_by_lpq.set(c_Key3(10, 7, 0), common_term_17);
    result_by_q.set(c_Key1(0), common_term_17);
    // q = 1
    result_by_lpq.set(c_Key3(10, 7, 1), common_term_16);
    result_by_q.set(c_Key1(1), common_term_16);
    // q = 2
    result_by_lpq.set(c_Key3(10, 7, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 7), result_by_q);
    result_by_q.clear();

    // l , p = (10, 8).
    // q = -2
    result_by_lpq.set(c_Key3(10, 8, -2), common_term_14);
    result_by_q.set(c_Key1(-2), common_term_14);
    // q = -1
    result_by_lpq.set(c_Key3(10, 8, -1), common_term_13);
    result_by_q.set(c_Key1(-1), common_term_13);
    // q = 0
    result_by_lpq.set(c_Key3(10, 8, 0), common_term_12);
    result_by_q.set(c_Key1(0), common_term_12);
    // q = 1
    result_by_lpq.set(c_Key3(10, 8, 1), common_term_11);
    result_by_q.set(c_Key1(1), common_term_11);
    // q = 2
    result_by_lpq.set(c_Key3(10, 8, 2), common_term_10);
    result_by_q.set(c_Key1(2), common_term_10);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 8), result_by_q);
    result_by_q.clear();

    // l , p = (10, 9).
    // q = -2
    result_by_lpq.set(c_Key3(10, 9, -2), common_term_9);
    result_by_q.set(c_Key1(-2), common_term_9);
    // q = -1
    result_by_lpq.set(c_Key3(10, 9, -1), common_term_8);
    result_by_q.set(c_Key1(-1), common_term_8);
    // q = 0
    result_by_lpq.set(c_Key3(10, 9, 0), common_term_7);
    result_by_q.set(c_Key1(0), common_term_7);
    // q = 1
    result_by_lpq.set(c_Key3(10, 9, 1), common_term_6);
    result_by_q.set(c_Key1(1), common_term_6);
    // q = 2
    result_by_lpq.set(c_Key3(10, 9, 2), common_term_5);
    result_by_q.set(c_Key1(2), common_term_5);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 9), result_by_q);
    result_by_q.clear();

    // l , p = (10, 10).
    // q = -2
    result_by_lpq.set(c_Key3(10, 10, -2), common_term_4);
    result_by_q.set(c_Key1(-2), common_term_4);
    // q = -1
    result_by_lpq.set(c_Key3(10, 10, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    result_by_lpq.set(c_Key3(10, 10, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(10, 10, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // q = 2
    result_by_lpq.set(c_Key3(10, 10, 2), common_term_0);
    result_by_q.set(c_Key1(2), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 10), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l10_e4(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 10.
    // Functions are taylor expanded and truncated at eccentricity^4.

    // Functions Calculated for l = 10.

    c_IntMap<c_Key3, double> result_by_lpq(79);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(11);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double common_term_0 = -7.1458333333333333*eccentricity_3;
    double common_term_1 = 8.0*eccentricity_2;
    double common_term_2 = 135.5625*eccentricity_3 - 4.5*eccentricity;
    double common_term_3 = 1.0 - 72.5*eccentricity_2;
    double common_term_4 = -633.1875*eccentricity_3 + 15.5*eccentricity;
    double common_term_5 = 130.5*eccentricity_2;
    double common_term_6 = 790.77083333333333*eccentricity_3;
    double common_term_7 = -0.85416666666666667*eccentricity_3;
    double common_term_8 = 2.25*eccentricity_2;
    double common_term_9 = 36.1875*eccentricity_3 - 2.5*eccentricity;
    double common_term_10 = 1.0 - 36.5*eccentricity_2;
    double common_term_11 = -290.8125*eccentricity_3 + 13.5*eccentricity;
    double common_term_12 = 100.25*eccentricity_2;
    double common_term_13 = 541.47916666666667*eccentricity_3;
    double common_term_14 = 0.4375*eccentricity_3;
    double common_term_15 = 0.5*eccentricity_2;
    double common_term_16 = 5.8125*eccentricity_3 - 0.5*eccentricity;
    double common_term_17 = 1.0 - 8.5*eccentricity_2;
    double common_term_18 = -71.4375*eccentricity_3 + 11.5*eccentricity;
    double common_term_19 = 74.0*eccentricity_2;
    double common_term_20 = 351.1875*eccentricity_3;
    double common_term_21 = 7.875*eccentricity_4*std::pow(1.0 - eccentricity_2, -9.5);
    double common_term_22 = 4.7291666666666667*eccentricity_3;
    double common_term_23 = 2.75*eccentricity_2;
    double common_term_24 = 20.4375*eccentricity_3 + 1.5*eccentricity;
    double common_term_25 = 11.5*eccentricity_2 + 1.0;
    double common_term_26 = 48.9375*eccentricity_3 + 9.5*eccentricity;
    double common_term_27 = 51.75*eccentricity_2;
    double common_term_28 = 211.89583333333333*eccentricity_3;
    double common_term_29 = 20.020833333333333*eccentricity_3;
    double common_term_30 = std::pow(1.0 - eccentricity_2, -9.5)*(31.5*eccentricity_4 + 9.0*eccentricity_2);
    double common_term_31 = 56.0625*eccentricity_3 + 3.5*eccentricity;
    double common_term_32 = 23.5*eccentricity_2 + 1.0;
    double common_term_33 = 94.3125*eccentricity_3 + 7.5*eccentricity;
    double common_term_34 = 33.5*eccentricity_2;
    double common_term_35 = 115.60416666666667*eccentricity_3;
    double common_term_36 = 54.3125*eccentricity_3;
    double common_term_37 = 19.25*eccentricity_2;
    double common_term_38 = 88.6875*eccentricity_3 + 5.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(8);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (10, 0).
    // q = -3
    result_by_lpq.set(c_Key3(10, 0, -3), common_term_0);
    result_by_q.set(c_Key1(-3), common_term_0);
    // q = -2
    result_by_lpq.set(c_Key3(10, 0, -2), common_term_1);
    result_by_q.set(c_Key1(-2), common_term_1);
    // q = -1
    result_by_lpq.set(c_Key3(10, 0, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    result_by_lpq.set(c_Key3(10, 0, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(10, 0, 1), common_term_4);
    result_by_q.set(c_Key1(1), common_term_4);
    // q = 2
    result_by_lpq.set(c_Key3(10, 0, 2), common_term_5);
    result_by_q.set(c_Key1(2), common_term_5);
    // q = 3
    result_by_lpq.set(c_Key3(10, 0, 3), common_term_6);
    result_by_q.set(c_Key1(3), common_term_6);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 0), result_by_q);
    result_by_q.clear();

    // l , p = (10, 1).
    // q = -3
    result_by_lpq.set(c_Key3(10, 1, -3), common_term_7);
    result_by_q.set(c_Key1(-3), common_term_7);
    // q = -2
    result_by_lpq.set(c_Key3(10, 1, -2), common_term_8);
    result_by_q.set(c_Key1(-2), common_term_8);
    // q = -1
    result_by_lpq.set(c_Key3(10, 1, -1), common_term_9);
    result_by_q.set(c_Key1(-1), common_term_9);
    // q = 0
    result_by_lpq.set(c_Key3(10, 1, 0), common_term_10);
    result_by_q.set(c_Key1(0), common_term_10);
    // q = 1
    result_by_lpq.set(c_Key3(10, 1, 1), common_term_11);
    result_by_q.set(c_Key1(1), common_term_11);
    // q = 2
    result_by_lpq.set(c_Key3(10, 1, 2), common_term_12);
    result_by_q.set(c_Key1(2), common_term_12);
    // q = 3
    result_by_lpq.set(c_Key3(10, 1, 3), common_term_13);
    result_by_q.set(c_Key1(3), common_term_13);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 1), result_by_q);
    result_by_q.clear();

    // l , p = (10, 2).
    // q = -3
    result_by_lpq.set(c_Key3(10, 2, -3), common_term_14);
    result_by_q.set(c_Key1(-3), common_term_14);
    // q = -2
    result_by_lpq.set(c_Key3(10, 2, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(10, 2, -1), common_term_16);
    result_by_q.set(c_Key1(-1), common_term_16);
    // q = 0
    result_by_lpq.set(c_Key3(10, 2, 0), common_term_17);
    result_by_q.set(c_Key1(0), common_term_17);
    // q = 1
    result_by_lpq.set(c_Key3(10, 2, 1), common_term_18);
    result_by_q.set(c_Key1(1), common_term_18);
    // q = 2
    result_by_lpq.set(c_Key3(10, 2, 2), common_term_19);
    result_by_q.set(c_Key1(2), common_term_19);
    // q = 3
    result_by_lpq.set(c_Key3(10, 2, 3), common_term_20);
    result_by_q.set(c_Key1(3), common_term_20);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 2), result_by_q);
    result_by_q.clear();

    // l , p = (10, 3).
    // q = -4
    result_by_lpq.set(c_Key3(10, 3, -4), common_term_21);
    result_by_q.set(c_Key1(-4), common_term_21);
    // q = -3
    result_by_lpq.set(c_Key3(10, 3, -3), common_term_22);
    result_by_q.set(c_Key1(-3), common_term_22);
    // q = -2
    result_by_lpq.set(c_Key3(10, 3, -2), common_term_23);
    result_by_q.set(c_Key1(-2), common_term_23);
    // q = -1
    result_by_lpq.set(c_Key3(10, 3, -1), common_term_24);
    result_by_q.set(c_Key1(-1), common_term_24);
    // q = 0
    result_by_lpq.set(c_Key3(10, 3, 0), common_term_25);
    result_by_q.set(c_Key1(0), common_term_25);
    // q = 1
    result_by_lpq.set(c_Key3(10, 3, 1), common_term_26);
    result_by_q.set(c_Key1(1), common_term_26);
    // q = 2
    result_by_lpq.set(c_Key3(10, 3, 2), common_term_27);
    result_by_q.set(c_Key1(2), common_term_27);
    // q = 3
    result_by_lpq.set(c_Key3(10, 3, 3), common_term_28);
    result_by_q.set(c_Key1(3), common_term_28);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 3), result_by_q);
    result_by_q.clear();

    // l , p = (10, 4).
    // q = -3
    result_by_lpq.set(c_Key3(10, 4, -3), common_term_29);
    result_by_q.set(c_Key1(-3), common_term_29);
    // q = -2
    result_by_lpq.set(c_Key3(10, 4, -2), common_term_30);
    result_by_q.set(c_Key1(-2), common_term_30);
    // q = -1
    result_by_lpq.set(c_Key3(10, 4, -1), common_term_31);
    result_by_q.set(c_Key1(-1), common_term_31);
    // q = 0
    result_by_lpq.set(c_Key3(10, 4, 0), common_term_32);
    result_by_q.set(c_Key1(0), common_term_32);
    // q = 1
    result_by_lpq.set(c_Key3(10, 4, 1), common_term_33);
    result_by_q.set(c_Key1(1), common_term_33);
    // q = 2
    result_by_lpq.set(c_Key3(10, 4, 2), common_term_34);
    result_by_q.set(c_Key1(2), common_term_34);
    // q = 3
    result_by_lpq.set(c_Key3(10, 4, 3), common_term_35);
    result_by_q.set(c_Key1(3), common_term_35);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 4), result_by_q);
    result_by_q.clear();

    // l , p = (10, 5).
    // q = -3
    result_by_lpq.set(c_Key3(10, 5, -3), common_term_36);
    result_by_q.set(c_Key1(-3), common_term_36);
    // q = -2
    result_by_lpq.set(c_Key3(10, 5, -2), common_term_37);
    result_by_q.set(c_Key1(-2), common_term_37);
    // q = -1
    result_by_lpq.set(c_Key3(10, 5, -1), common_term_38);
    result_by_q.set(c_Key1(-1), common_term_38);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -9.5)*(47.25*eccentricity_4 + 18.0*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(10, 5, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(10, 5, 1), common_term_38);
    result_by_q.set(c_Key1(1), common_term_38);
    // q = 2
    result_by_lpq.set(c_Key3(10, 5, 2), common_term_37);
    result_by_q.set(c_Key1(2), common_term_37);
    // q = 3
    result_by_lpq.set(c_Key3(10, 5, 3), common_term_36);
    result_by_q.set(c_Key1(3), common_term_36);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 5), result_by_q);
    result_by_q.clear();

    // l , p = (10, 6).
    // q = -3
    result_by_lpq.set(c_Key3(10, 6, -3), common_term_35);
    result_by_q.set(c_Key1(-3), common_term_35);
    // q = -2
    result_by_lpq.set(c_Key3(10, 6, -2), common_term_34);
    result_by_q.set(c_Key1(-2), common_term_34);
    // q = -1
    result_by_lpq.set(c_Key3(10, 6, -1), common_term_33);
    result_by_q.set(c_Key1(-1), common_term_33);
    // q = 0
    result_by_lpq.set(c_Key3(10, 6, 0), common_term_32);
    result_by_q.set(c_Key1(0), common_term_32);
    // q = 1
    result_by_lpq.set(c_Key3(10, 6, 1), common_term_31);
    result_by_q.set(c_Key1(1), common_term_31);
    // q = 2
    result_by_lpq.set(c_Key3(10, 6, 2), common_term_30);
    result_by_q.set(c_Key1(2), common_term_30);
    // q = 3
    result_by_lpq.set(c_Key3(10, 6, 3), common_term_29);
    result_by_q.set(c_Key1(3), common_term_29);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 6), result_by_q);
    result_by_q.clear();

    // l , p = (10, 7).
    // q = -3
    result_by_lpq.set(c_Key3(10, 7, -3), common_term_28);
    result_by_q.set(c_Key1(-3), common_term_28);
    // q = -2
    result_by_lpq.set(c_Key3(10, 7, -2), common_term_27);
    result_by_q.set(c_Key1(-2), common_term_27);
    // q = -1
    result_by_lpq.set(c_Key3(10, 7, -1), common_term_26);
    result_by_q.set(c_Key1(-1), common_term_26);
    // q = 0
    result_by_lpq.set(c_Key3(10, 7, 0), common_term_25);
    result_by_q.set(c_Key1(0), common_term_25);
    // q = 1
    result_by_lpq.set(c_Key3(10, 7, 1), common_term_24);
    result_by_q.set(c_Key1(1), common_term_24);
    // q = 2
    result_by_lpq.set(c_Key3(10, 7, 2), common_term_23);
    result_by_q.set(c_Key1(2), common_term_23);
    // q = 3
    result_by_lpq.set(c_Key3(10, 7, 3), common_term_22);
    result_by_q.set(c_Key1(3), common_term_22);
    // q = 4
    result_by_lpq.set(c_Key3(10, 7, 4), common_term_21);
    result_by_q.set(c_Key1(4), common_term_21);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 7), result_by_q);
    result_by_q.clear();

    // l , p = (10, 8).
    // q = -3
    result_by_lpq.set(c_Key3(10, 8, -3), common_term_20);
    result_by_q.set(c_Key1(-3), common_term_20);
    // q = -2
    result_by_lpq.set(c_Key3(10, 8, -2), common_term_19);
    result_by_q.set(c_Key1(-2), common_term_19);
    // q = -1
    result_by_lpq.set(c_Key3(10, 8, -1), common_term_18);
    result_by_q.set(c_Key1(-1), common_term_18);
    // q = 0
    result_by_lpq.set(c_Key3(10, 8, 0), common_term_17);
    result_by_q.set(c_Key1(0), common_term_17);
    // q = 1
    result_by_lpq.set(c_Key3(10, 8, 1), common_term_16);
    result_by_q.set(c_Key1(1), common_term_16);
    // q = 2
    result_by_lpq.set(c_Key3(10, 8, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // q = 3
    result_by_lpq.set(c_Key3(10, 8, 3), common_term_14);
    result_by_q.set(c_Key1(3), common_term_14);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 8), result_by_q);
    result_by_q.clear();

    // l , p = (10, 9).
    // q = -3
    result_by_lpq.set(c_Key3(10, 9, -3), common_term_13);
    result_by_q.set(c_Key1(-3), common_term_13);
    // q = -2
    result_by_lpq.set(c_Key3(10, 9, -2), common_term_12);
    result_by_q.set(c_Key1(-2), common_term_12);
    // q = -1
    result_by_lpq.set(c_Key3(10, 9, -1), common_term_11);
    result_by_q.set(c_Key1(-1), common_term_11);
    // q = 0
    result_by_lpq.set(c_Key3(10, 9, 0), common_term_10);
    result_by_q.set(c_Key1(0), common_term_10);
    // q = 1
    result_by_lpq.set(c_Key3(10, 9, 1), common_term_9);
    result_by_q.set(c_Key1(1), common_term_9);
    // q = 2
    result_by_lpq.set(c_Key3(10, 9, 2), common_term_8);
    result_by_q.set(c_Key1(2), common_term_8);
    // q = 3
    result_by_lpq.set(c_Key3(10, 9, 3), common_term_7);
    result_by_q.set(c_Key1(3), common_term_7);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 9), result_by_q);
    result_by_q.clear();

    // l , p = (10, 10).
    // q = -3
    result_by_lpq.set(c_Key3(10, 10, -3), common_term_6);
    result_by_q.set(c_Key1(-3), common_term_6);
    // q = -2
    result_by_lpq.set(c_Key3(10, 10, -2), common_term_5);
    result_by_q.set(c_Key1(-2), common_term_5);
    // q = -1
    result_by_lpq.set(c_Key3(10, 10, -1), common_term_4);
    result_by_q.set(c_Key1(-1), common_term_4);
    // q = 0
    result_by_lpq.set(c_Key3(10, 10, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(10, 10, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // q = 2
    result_by_lpq.set(c_Key3(10, 10, 2), common_term_1);
    result_by_q.set(c_Key1(2), common_term_1);
    // q = 3
    result_by_lpq.set(c_Key3(10, 10, 3), common_term_0);
    result_by_q.set(c_Key1(3), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 10), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l10_e5(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 10.
    // Functions are taylor expanded and truncated at eccentricity^5.

    // Functions Calculated for l = 10.

    c_IntMap<c_Key3, double> result_by_lpq(99);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(11);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double common_term_0 = 3.375*eccentricity_4;
    double common_term_1 = -7.1458333333333333*eccentricity_3;
    double common_term_2 = -129.33333333333333*eccentricity_4 + 8.0*eccentricity_2;
    double common_term_3 = 135.5625*eccentricity_3 - 4.5*eccentricity;
    double common_term_4 = 1241.5625*eccentricity_4 - 72.5*eccentricity_2 + 1.0;
    double common_term_5 = -633.1875*eccentricity_3 + 15.5*eccentricity;
    double common_term_6 = -3973.5*eccentricity_4 + 130.5*eccentricity_2;
    double common_term_7 = 790.77083333333333*eccentricity_3;
    double common_term_8 = 3858.8958333333334*eccentricity_4;
    double common_term_9 = 0.14583333333333333*eccentricity_4;
    double common_term_10 = -0.85416666666666667*eccentricity_3;
    double common_term_11 = -15.75*eccentricity_4 + 2.25*eccentricity_2;
    double common_term_12 = 36.1875*eccentricity_3 - 2.5*eccentricity;
    double common_term_13 = 315.125*eccentricity_4 - 36.5*eccentricity_2 + 1.0;
    double common_term_14 = -290.8125*eccentricity_3 + 13.5*eccentricity;
    double common_term_15 = -1675.0833333333333*eccentricity_4 + 100.25*eccentricity_2;
    double common_term_16 = 541.47916666666667*eccentricity_3;
    double common_term_17 = 2376.5625*eccentricity_4;
    double common_term_18 = 0.64583333333333333*eccentricity_4;
    double common_term_19 = 0.4375*eccentricity_3;
    double common_term_20 = 3.1666666666666667*eccentricity_4 + 0.5*eccentricity_2;
    double common_term_21 = 5.8125*eccentricity_3 - 0.5*eccentricity;
    double common_term_22 = 42.5625*eccentricity_4 - 8.5*eccentricity_2 + 1.0;
    double common_term_23 = -71.4375*eccentricity_3 + 11.5*eccentricity;
    double common_term_24 = -415.33333333333333*eccentricity_4 + 74.0*eccentricity_2;
    double common_term_25 = 351.1875*eccentricity_3;
    double common_term_26 = 1369.9583333333333*eccentricity_4;
    double common_term_27 = 7.875*eccentricity_4*std::pow(1.0 - eccentricity_2, -9.5);
    double common_term_28 = 4.7291666666666667*eccentricity_3;
    double common_term_29 = 33.416666666666667*eccentricity_4 + 2.75*eccentricity_2;
    double common_term_30 = 20.4375*eccentricity_3 + 1.5*eccentricity;
    double common_term_31 = 87.875*eccentricity_4 + 11.5*eccentricity_2 + 1.0;
    double common_term_32 = 48.9375*eccentricity_3 + 9.5*eccentricity;
    double common_term_33 = 147.75*eccentricity_4 + 51.75*eccentricity_2;
    double common_term_34 = 211.89583333333333*eccentricity_3;
    double common_term_35 = 724.08333333333333*eccentricity_4;
    double common_term_36 = 40.833333333333333*eccentricity_4;
    double common_term_37 = 20.020833333333333*eccentricity_3;
    double common_term_38 = std::pow(1.0 - eccentricity_2, -9.5)*(31.5*eccentricity_4 + 9.0*eccentricity_2);
    double common_term_39 = 56.0625*eccentricity_3 + 3.5*eccentricity;
    double common_term_40 = 211.0625*eccentricity_4 + 23.5*eccentricity_2 + 1.0;
    double common_term_41 = 94.3125*eccentricity_3 + 7.5*eccentricity;
    double common_term_42 = 292.16666666666667*eccentricity_4 + 33.5*eccentricity_2;
    double common_term_43 = 115.60416666666667*eccentricity_3;
    double common_term_44 = 339.9375*eccentricity_4;
    double common_term_45 = 134.52083333333333*eccentricity_4;
    double common_term_46 = 54.3125*eccentricity_3;
    double common_term_47 = 231.91666666666667*eccentricity_4 + 19.25*eccentricity_2;
    double common_term_48 = 88.6875*eccentricity_3 + 5.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(9);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (10, 0).
    // q = -4
    result_by_lpq.set(c_Key3(10, 0, -4), common_term_0);
    result_by_q.set(c_Key1(-4), common_term_0);
    // q = -3
    result_by_lpq.set(c_Key3(10, 0, -3), common_term_1);
    result_by_q.set(c_Key1(-3), common_term_1);
    // q = -2
    result_by_lpq.set(c_Key3(10, 0, -2), common_term_2);
    result_by_q.set(c_Key1(-2), common_term_2);
    // q = -1
    result_by_lpq.set(c_Key3(10, 0, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    result_by_lpq.set(c_Key3(10, 0, 0), common_term_4);
    result_by_q.set(c_Key1(0), common_term_4);
    // q = 1
    result_by_lpq.set(c_Key3(10, 0, 1), common_term_5);
    result_by_q.set(c_Key1(1), common_term_5);
    // q = 2
    result_by_lpq.set(c_Key3(10, 0, 2), common_term_6);
    result_by_q.set(c_Key1(2), common_term_6);
    // q = 3
    result_by_lpq.set(c_Key3(10, 0, 3), common_term_7);
    result_by_q.set(c_Key1(3), common_term_7);
    // q = 4
    result_by_lpq.set(c_Key3(10, 0, 4), common_term_8);
    result_by_q.set(c_Key1(4), common_term_8);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 0), result_by_q);
    result_by_q.clear();

    // l , p = (10, 1).
    // q = -4
    result_by_lpq.set(c_Key3(10, 1, -4), common_term_9);
    result_by_q.set(c_Key1(-4), common_term_9);
    // q = -3
    result_by_lpq.set(c_Key3(10, 1, -3), common_term_10);
    result_by_q.set(c_Key1(-3), common_term_10);
    // q = -2
    result_by_lpq.set(c_Key3(10, 1, -2), common_term_11);
    result_by_q.set(c_Key1(-2), common_term_11);
    // q = -1
    result_by_lpq.set(c_Key3(10, 1, -1), common_term_12);
    result_by_q.set(c_Key1(-1), common_term_12);
    // q = 0
    result_by_lpq.set(c_Key3(10, 1, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(10, 1, 1), common_term_14);
    result_by_q.set(c_Key1(1), common_term_14);
    // q = 2
    result_by_lpq.set(c_Key3(10, 1, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // q = 3
    result_by_lpq.set(c_Key3(10, 1, 3), common_term_16);
    result_by_q.set(c_Key1(3), common_term_16);
    // q = 4
    result_by_lpq.set(c_Key3(10, 1, 4), common_term_17);
    result_by_q.set(c_Key1(4), common_term_17);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 1), result_by_q);
    result_by_q.clear();

    // l , p = (10, 2).
    // q = -4
    result_by_lpq.set(c_Key3(10, 2, -4), common_term_18);
    result_by_q.set(c_Key1(-4), common_term_18);
    // q = -3
    result_by_lpq.set(c_Key3(10, 2, -3), common_term_19);
    result_by_q.set(c_Key1(-3), common_term_19);
    // q = -2
    result_by_lpq.set(c_Key3(10, 2, -2), common_term_20);
    result_by_q.set(c_Key1(-2), common_term_20);
    // q = -1
    result_by_lpq.set(c_Key3(10, 2, -1), common_term_21);
    result_by_q.set(c_Key1(-1), common_term_21);
    // q = 0
    result_by_lpq.set(c_Key3(10, 2, 0), common_term_22);
    result_by_q.set(c_Key1(0), common_term_22);
    // q = 1
    result_by_lpq.set(c_Key3(10, 2, 1), common_term_23);
    result_by_q.set(c_Key1(1), common_term_23);
    // q = 2
    result_by_lpq.set(c_Key3(10, 2, 2), common_term_24);
    result_by_q.set(c_Key1(2), common_term_24);
    // q = 3
    result_by_lpq.set(c_Key3(10, 2, 3), common_term_25);
    result_by_q.set(c_Key1(3), common_term_25);
    // q = 4
    result_by_lpq.set(c_Key3(10, 2, 4), common_term_26);
    result_by_q.set(c_Key1(4), common_term_26);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 2), result_by_q);
    result_by_q.clear();

    // l , p = (10, 3).
    // q = -4
    result_by_lpq.set(c_Key3(10, 3, -4), common_term_27);
    result_by_q.set(c_Key1(-4), common_term_27);
    // q = -3
    result_by_lpq.set(c_Key3(10, 3, -3), common_term_28);
    result_by_q.set(c_Key1(-3), common_term_28);
    // q = -2
    result_by_lpq.set(c_Key3(10, 3, -2), common_term_29);
    result_by_q.set(c_Key1(-2), common_term_29);
    // q = -1
    result_by_lpq.set(c_Key3(10, 3, -1), common_term_30);
    result_by_q.set(c_Key1(-1), common_term_30);
    // q = 0
    result_by_lpq.set(c_Key3(10, 3, 0), common_term_31);
    result_by_q.set(c_Key1(0), common_term_31);
    // q = 1
    result_by_lpq.set(c_Key3(10, 3, 1), common_term_32);
    result_by_q.set(c_Key1(1), common_term_32);
    // q = 2
    result_by_lpq.set(c_Key3(10, 3, 2), common_term_33);
    result_by_q.set(c_Key1(2), common_term_33);
    // q = 3
    result_by_lpq.set(c_Key3(10, 3, 3), common_term_34);
    result_by_q.set(c_Key1(3), common_term_34);
    // q = 4
    result_by_lpq.set(c_Key3(10, 3, 4), common_term_35);
    result_by_q.set(c_Key1(4), common_term_35);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 3), result_by_q);
    result_by_q.clear();

    // l , p = (10, 4).
    // q = -4
    result_by_lpq.set(c_Key3(10, 4, -4), common_term_36);
    result_by_q.set(c_Key1(-4), common_term_36);
    // q = -3
    result_by_lpq.set(c_Key3(10, 4, -3), common_term_37);
    result_by_q.set(c_Key1(-3), common_term_37);
    // q = -2
    result_by_lpq.set(c_Key3(10, 4, -2), common_term_38);
    result_by_q.set(c_Key1(-2), common_term_38);
    // q = -1
    result_by_lpq.set(c_Key3(10, 4, -1), common_term_39);
    result_by_q.set(c_Key1(-1), common_term_39);
    // q = 0
    result_by_lpq.set(c_Key3(10, 4, 0), common_term_40);
    result_by_q.set(c_Key1(0), common_term_40);
    // q = 1
    result_by_lpq.set(c_Key3(10, 4, 1), common_term_41);
    result_by_q.set(c_Key1(1), common_term_41);
    // q = 2
    result_by_lpq.set(c_Key3(10, 4, 2), common_term_42);
    result_by_q.set(c_Key1(2), common_term_42);
    // q = 3
    result_by_lpq.set(c_Key3(10, 4, 3), common_term_43);
    result_by_q.set(c_Key1(3), common_term_43);
    // q = 4
    result_by_lpq.set(c_Key3(10, 4, 4), common_term_44);
    result_by_q.set(c_Key1(4), common_term_44);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 4), result_by_q);
    result_by_q.clear();

    // l , p = (10, 5).
    // q = -4
    result_by_lpq.set(c_Key3(10, 5, -4), common_term_45);
    result_by_q.set(c_Key1(-4), common_term_45);
    // q = -3
    result_by_lpq.set(c_Key3(10, 5, -3), common_term_46);
    result_by_q.set(c_Key1(-3), common_term_46);
    // q = -2
    result_by_lpq.set(c_Key3(10, 5, -2), common_term_47);
    result_by_q.set(c_Key1(-2), common_term_47);
    // q = -1
    result_by_lpq.set(c_Key3(10, 5, -1), common_term_48);
    result_by_q.set(c_Key1(-1), common_term_48);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -9.5)*(47.25*eccentricity_4 + 18.0*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(10, 5, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(10, 5, 1), common_term_48);
    result_by_q.set(c_Key1(1), common_term_48);
    // q = 2
    result_by_lpq.set(c_Key3(10, 5, 2), common_term_47);
    result_by_q.set(c_Key1(2), common_term_47);
    // q = 3
    result_by_lpq.set(c_Key3(10, 5, 3), common_term_46);
    result_by_q.set(c_Key1(3), common_term_46);
    // q = 4
    result_by_lpq.set(c_Key3(10, 5, 4), common_term_45);
    result_by_q.set(c_Key1(4), common_term_45);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 5), result_by_q);
    result_by_q.clear();

    // l , p = (10, 6).
    // q = -4
    result_by_lpq.set(c_Key3(10, 6, -4), common_term_44);
    result_by_q.set(c_Key1(-4), common_term_44);
    // q = -3
    result_by_lpq.set(c_Key3(10, 6, -3), common_term_43);
    result_by_q.set(c_Key1(-3), common_term_43);
    // q = -2
    result_by_lpq.set(c_Key3(10, 6, -2), common_term_42);
    result_by_q.set(c_Key1(-2), common_term_42);
    // q = -1
    result_by_lpq.set(c_Key3(10, 6, -1), common_term_41);
    result_by_q.set(c_Key1(-1), common_term_41);
    // q = 0
    result_by_lpq.set(c_Key3(10, 6, 0), common_term_40);
    result_by_q.set(c_Key1(0), common_term_40);
    // q = 1
    result_by_lpq.set(c_Key3(10, 6, 1), common_term_39);
    result_by_q.set(c_Key1(1), common_term_39);
    // q = 2
    result_by_lpq.set(c_Key3(10, 6, 2), common_term_38);
    result_by_q.set(c_Key1(2), common_term_38);
    // q = 3
    result_by_lpq.set(c_Key3(10, 6, 3), common_term_37);
    result_by_q.set(c_Key1(3), common_term_37);
    // q = 4
    result_by_lpq.set(c_Key3(10, 6, 4), common_term_36);
    result_by_q.set(c_Key1(4), common_term_36);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 6), result_by_q);
    result_by_q.clear();

    // l , p = (10, 7).
    // q = -4
    result_by_lpq.set(c_Key3(10, 7, -4), common_term_35);
    result_by_q.set(c_Key1(-4), common_term_35);
    // q = -3
    result_by_lpq.set(c_Key3(10, 7, -3), common_term_34);
    result_by_q.set(c_Key1(-3), common_term_34);
    // q = -2
    result_by_lpq.set(c_Key3(10, 7, -2), common_term_33);
    result_by_q.set(c_Key1(-2), common_term_33);
    // q = -1
    result_by_lpq.set(c_Key3(10, 7, -1), common_term_32);
    result_by_q.set(c_Key1(-1), common_term_32);
    // q = 0
    result_by_lpq.set(c_Key3(10, 7, 0), common_term_31);
    result_by_q.set(c_Key1(0), common_term_31);
    // q = 1
    result_by_lpq.set(c_Key3(10, 7, 1), common_term_30);
    result_by_q.set(c_Key1(1), common_term_30);
    // q = 2
    result_by_lpq.set(c_Key3(10, 7, 2), common_term_29);
    result_by_q.set(c_Key1(2), common_term_29);
    // q = 3
    result_by_lpq.set(c_Key3(10, 7, 3), common_term_28);
    result_by_q.set(c_Key1(3), common_term_28);
    // q = 4
    result_by_lpq.set(c_Key3(10, 7, 4), common_term_27);
    result_by_q.set(c_Key1(4), common_term_27);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 7), result_by_q);
    result_by_q.clear();

    // l , p = (10, 8).
    // q = -4
    result_by_lpq.set(c_Key3(10, 8, -4), common_term_26);
    result_by_q.set(c_Key1(-4), common_term_26);
    // q = -3
    result_by_lpq.set(c_Key3(10, 8, -3), common_term_25);
    result_by_q.set(c_Key1(-3), common_term_25);
    // q = -2
    result_by_lpq.set(c_Key3(10, 8, -2), common_term_24);
    result_by_q.set(c_Key1(-2), common_term_24);
    // q = -1
    result_by_lpq.set(c_Key3(10, 8, -1), common_term_23);
    result_by_q.set(c_Key1(-1), common_term_23);
    // q = 0
    result_by_lpq.set(c_Key3(10, 8, 0), common_term_22);
    result_by_q.set(c_Key1(0), common_term_22);
    // q = 1
    result_by_lpq.set(c_Key3(10, 8, 1), common_term_21);
    result_by_q.set(c_Key1(1), common_term_21);
    // q = 2
    result_by_lpq.set(c_Key3(10, 8, 2), common_term_20);
    result_by_q.set(c_Key1(2), common_term_20);
    // q = 3
    result_by_lpq.set(c_Key3(10, 8, 3), common_term_19);
    result_by_q.set(c_Key1(3), common_term_19);
    // q = 4
    result_by_lpq.set(c_Key3(10, 8, 4), common_term_18);
    result_by_q.set(c_Key1(4), common_term_18);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 8), result_by_q);
    result_by_q.clear();

    // l , p = (10, 9).
    // q = -4
    result_by_lpq.set(c_Key3(10, 9, -4), common_term_17);
    result_by_q.set(c_Key1(-4), common_term_17);
    // q = -3
    result_by_lpq.set(c_Key3(10, 9, -3), common_term_16);
    result_by_q.set(c_Key1(-3), common_term_16);
    // q = -2
    result_by_lpq.set(c_Key3(10, 9, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(10, 9, -1), common_term_14);
    result_by_q.set(c_Key1(-1), common_term_14);
    // q = 0
    result_by_lpq.set(c_Key3(10, 9, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(10, 9, 1), common_term_12);
    result_by_q.set(c_Key1(1), common_term_12);
    // q = 2
    result_by_lpq.set(c_Key3(10, 9, 2), common_term_11);
    result_by_q.set(c_Key1(2), common_term_11);
    // q = 3
    result_by_lpq.set(c_Key3(10, 9, 3), common_term_10);
    result_by_q.set(c_Key1(3), common_term_10);
    // q = 4
    result_by_lpq.set(c_Key3(10, 9, 4), common_term_9);
    result_by_q.set(c_Key1(4), common_term_9);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 9), result_by_q);
    result_by_q.clear();

    // l , p = (10, 10).
    // q = -4
    result_by_lpq.set(c_Key3(10, 10, -4), common_term_8);
    result_by_q.set(c_Key1(-4), common_term_8);
    // q = -3
    result_by_lpq.set(c_Key3(10, 10, -3), common_term_7);
    result_by_q.set(c_Key1(-3), common_term_7);
    // q = -2
    result_by_lpq.set(c_Key3(10, 10, -2), common_term_6);
    result_by_q.set(c_Key1(-2), common_term_6);
    // q = -1
    result_by_lpq.set(c_Key3(10, 10, -1), common_term_5);
    result_by_q.set(c_Key1(-1), common_term_5);
    // q = 0
    result_by_lpq.set(c_Key3(10, 10, 0), common_term_4);
    result_by_q.set(c_Key1(0), common_term_4);
    // q = 1
    result_by_lpq.set(c_Key3(10, 10, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // q = 2
    result_by_lpq.set(c_Key3(10, 10, 2), common_term_2);
    result_by_q.set(c_Key1(2), common_term_2);
    // q = 3
    result_by_lpq.set(c_Key3(10, 10, 3), common_term_1);
    result_by_q.set(c_Key1(3), common_term_1);
    // q = 4
    result_by_lpq.set(c_Key3(10, 10, 4), common_term_0);
    result_by_q.set(c_Key1(4), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 10), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l10_e10(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 10.
    // Functions are taylor expanded and truncated at eccentricity^10.

    // Functions Calculated for l = 10.

    c_IntMap<c_Key3, double> result_by_lpq(209);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(11);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_8 = eccentricity_4 * eccentricity_4;
    double eccentricity_9 = eccentricity_8 * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double eccentricity_6 = eccentricity_3 * eccentricity_3;
    double eccentricity_7 = eccentricity_6 * eccentricity;
    double eccentricity_5 = eccentricity_4 * eccentricity;
    double common_term_0 = -5.3822889109347443e-9*eccentricity_9;
    double common_term_1 = 2.4801587301587302e-5*eccentricity_8;
    double common_term_2 = -0.0011653355189732143*eccentricity_9 - 0.0033900669642857143*eccentricity_7;
    double common_term_3 = -0.08253968253968254*eccentricity_8 + 0.088888888888888889*eccentricity_6;
    double common_term_4 = -1.6257876441592262*eccentricity_9 + 2.2040473090277778*eccentricity_7 - 0.81380208333333333*eccentricity_5;
    double common_term_5 = 30.9234375*eccentricity_8 - 17.8875*eccentricity_6 + 3.375*eccentricity_4;
    double common_term_6 = 308.33914161964699*eccentricity_9 - 212.15309244791667*eccentricity_7 + 66.545572916666667*eccentricity_5 - 7.1458333333333333*eccentricity_3;
    double common_term_7 = -1784.6055555555556*eccentricity_8 + 704.5*eccentricity_6 - 129.33333333333333*eccentricity_4 + 8.0*eccentricity_2;
    double common_term_8 = -11881.492437744141*eccentricity_9 + 5306.87548828125*eccentricity_7 - 1260.4921875*eccentricity_5 + 135.5625*eccentricity_3 - 4.5*eccentricity;
    double common_term_9 = 31844.422743055556*eccentricity_8 - 8775.9548611111111*eccentricity_6 + 1241.5625*eccentricity_4 - 72.5*eccentricity_2 + 1.0;
    double common_term_10 = 161721.53108113607*eccentricity_9 - 48884.568739149307*eccentricity_7 + 8150.3411458333333*eccentricity_5 - 633.1875*eccentricity_3 + 15.5*eccentricity;
    double common_term_11 = -231454.8*eccentricity_8 + 42937.312500000001*eccentricity_6 - 3973.5*eccentricity_4 + 130.5*eccentricity_2;
    double common_term_12 = -966863.29034333404*eccentricity_9 + 192809.99469401042*eccentricity_7 - 20053.194010416667*eccentricity_5 + 790.77083333333333*eccentricity_3;
    double common_term_13 = 765856.92907986111*eccentricity_8 - 86424.668750000001*eccentricity_6 + 3858.8958333333334*eccentricity_4;
    double common_term_14 = 2758410.2350289481*eccentricity_9 - 330019.78271484375*eccentricity_7 + 16100.61328125*eccentricity_5;
    double common_term_15 = -1144440.9515873016*eccentricity_8 + 59584.609722222222*eccentricity_6;
    double common_term_16 = -3667624.6796146938*eccentricity_9 + 200445.43370690724*eccentricity_7;
    double common_term_17 = 623706.16439732142*eccentricity_8;
    double common_term_18 = 1818308.9840970532*eccentricity_9;
    double common_term_19 = 0.047782297824005181*eccentricity_9;
    double common_term_20 = 0.03515625*eccentricity_8*std::pow(1.0 - eccentricity_2, -9.5);
    double common_term_21 = 0.25219857352120536*eccentricity_9 + 0.025866505456349206*eccentricity_7;
    double common_term_22 = 0.19037698412698413*eccentricity_8 + 0.019097222222222222*eccentricity_6;
    double common_term_23 = 0.80217459542410714*eccentricity_9 + 0.14326171875*eccentricity_7 + 0.00703125*eccentricity_5;
    double common_term_24 = 0.66111111111111111*eccentricity_8 - 0.06875*eccentricity_6 + 0.14583333333333333*eccentricity_4;
    double common_term_25 = 3.101062915943287*eccentricity_9 - 2.3777669270833333*eccentricity_7 + 3.0247395833333333*eccentricity_5 - 0.85416666666666667*eccentricity_3;
    double common_term_26 = -30.20625*eccentricity_8 + 35.015625*eccentricity_6 - 15.75*eccentricity_4 + 2.25*eccentricity_2;
    double common_term_27 = -260.62958780924479*eccentricity_9 + 288.20958116319444*eccentricity_7 - 155.44010416666667*eccentricity_5 + 36.1875*eccentricity_3 - 2.5*eccentricity;
    double common_term_28 = 1866.3055555555556*eccentricity_8 - 1089.5486111111111*eccentricity_6 + 315.125*eccentricity_4 - 36.5*eccentricity_2 + 1.0;
    double common_term_29 = 10084.092352294922*eccentricity_9 - 6078.41455078125*eccentricity_7 + 1970.2265625*eccentricity_5 - 290.8125*eccentricity_3 + 13.5*eccentricity;
    double common_term_30 = -28695.711805555556*eccentricity_8 + 9901.328125*eccentricity_6 - 1675.0833333333333*eccentricity_4 + 100.25*eccentricity_2;
    double common_term_31 = -119081.11242088035*eccentricity_9 + 42477.821712239583*eccentricity_7 - 7802.2669270833333*eccentricity_5 + 541.47916666666667*eccentricity_3;
    double common_term_32 = 161427.2625*eccentricity_8 - 31188.09375*eccentricity_6 + 2376.5625*eccentricity_4;
    double common_term_33 = 557026.68569394067*eccentricity_9 - 110947.63989800347*eccentricity_7 + 8987.4486979166667*eccentricity_5;
    double common_term_34 = -359842.25486111111*eccentricity_8 + 30349.348263888889*eccentricity_6;
    double common_term_35 = -1082391.1435743059*eccentricity_9 + 93710.231794084822*eccentricity_7;
    double common_term_36 = 269031.71821676587*eccentricity_8;
    double common_term_37 = 726995.90431837725*eccentricity_9;
    double common_term_38 = 3.7464246477399554*eccentricity_9;
    double common_term_39 = 2.6458209325396825*eccentricity_8;
    double common_term_40 = 17.51398184640067*eccentricity_9 + 1.8652793278769841*eccentricity_7;
    double common_term_41 = std::pow(1.0 - eccentricity_2, -9.5)*(0.28125*eccentricity_8 + 1.3125*eccentricity_6);
    double common_term_42 = 50.857536388578869*eccentricity_9 + 9.2587782118055556*eccentricity_7 + 0.92161458333333333*eccentricity_5;
    double common_term_43 = 37.861371527777778*eccentricity_8 + 6.70625*eccentricity_6 + 0.64583333333333333*eccentricity_4;
    double common_term_44 = 117.13199462890625*eccentricity_9 + 28.10771484375*eccentricity_7 + 4.84765625*eccentricity_5 + 0.4375*eccentricity_3;
    double common_term_45 = 88.711111111111111*eccentricity_8 + 20.9375*eccentricity_6 + 3.1666666666666667*eccentricity_4 + 0.5*eccentricity_2;
    double common_term_46 = 233.38388061523438*eccentricity_9 + 69.138726128472222*eccentricity_7 + 11.174479166666667*eccentricity_5 + 5.8125*eccentricity_3 - 0.5*eccentricity;
    double common_term_47 = 201.1015625*eccentricity_8 + 13.65625*eccentricity_6 + 42.5625*eccentricity_4 - 8.5*eccentricity_2 + 1.0;
    double common_term_48 = 582.28132527669271*eccentricity_9 - 107.68614366319444*eccentricity_7 + 245.67447916666667*eccentricity_5 - 71.4375*eccentricity_3 + 11.5*eccentricity;
    double common_term_49 = -998.73888888888889*eccentricity_8 + 1195.375*eccentricity_6 - 415.33333333333333*eccentricity_4 + 74.0*eccentricity_2;
    double common_term_50 = -5407.3172241210938*eccentricity_9 + 5068.97138671875*eccentricity_7 - 1908.99609375*eccentricity_5 + 351.1875*eccentricity_3;
    double common_term_51 = 19177.660590277778*eccentricity_8 - 7436.5625*eccentricity_6 + 1369.9583333333333*eccentricity_4;
    double common_term_52 = 65947.025514439174*eccentricity_9 - 25597.524945746528*eccentricity_7 + 4649.3778645833333*eccentricity_5;
    double common_term_53 = -80002.971428571429*eccentricity_8 + 14206.575*eccentricity_6;
    double common_term_54 = -231366.43252229236*eccentricity_9 + 39975.869204179067*eccentricity_7;
    double common_term_55 = 105238.97986111111*eccentricity_8;
    double common_term_56 = 262199.80908857073*eccentricity_9;
    double common_term_57 = 76.972772695820588*eccentricity_9;
    double common_term_58 = 50.005716765873016*eccentricity_8;
    double common_term_59 = 284.08580278669085*eccentricity_9 + 32.164327566964286*eccentricity_7;
    double common_term_60 = 190.68382936507936*eccentricity_8 + 20.443055555555556*eccentricity_6;
    double common_term_61 = 683.84254847935268*eccentricity_9 + 126.46257595486111*eccentricity_7 + 12.804947916666667*eccentricity_5;
    double common_term_62 = std::pow(1.0 - eccentricity_2, -9.5)*(0.984375*eccentricity_8 + 7.875*eccentricity_6 + 7.875*eccentricity_4);
    double common_term_63 = 1342.1533343279803*eccentricity_9 + 316.81741536458333*eccentricity_7 + 53.139322916666667*eccentricity_5 + 4.7291666666666667*eccentricity_3;
    double common_term_64 = 932.75069444444444*eccentricity_8 + 210.84375*eccentricity_6 + 33.416666666666667*eccentricity_4 + 2.75*eccentricity_2;
    double common_term_65 = 2326.2882019042969*eccentricity_9 + 638.93115234375*eccentricity_7 + 137.6015625*eccentricity_5 + 20.4375*eccentricity_3 + 1.5*eccentricity;
    double common_term_66 = 1632.9787326388889*eccentricity_8 + 429.99305555555556*eccentricity_6 + 87.875*eccentricity_4 + 11.5*eccentricity_2 + 1.0;
    double common_term_67 = 3702.8700052897135*eccentricity_9 + 1127.3065863715278*eccentricity_7 + 286.43489583333333*eccentricity_5 + 48.9375*eccentricity_3 + 9.5*eccentricity;
    double common_term_68 = 2607.13125*eccentricity_8 + 789.375*eccentricity_6 + 147.75*eccentricity_4 + 51.75*eccentricity_2;
    double common_term_69 = 5465.05249520761*eccentricity_9 + 1960.3046549479167*eccentricity_7 + 333.49348958333333*eccentricity_5 + 211.89583333333333*eccentricity_3;
    double common_term_70 = 4597.7611111111111*eccentricity_8 + 506.675*eccentricity_6 + 724.08333333333333*eccentricity_4;
    double common_term_71 = 10579.168385532924*eccentricity_9 + 65.10849609375*eccentricity_7 + 2179.52578125*eccentricity_5;
    double common_term_72 = -3120.3387896825397*eccentricity_8 + 5967.9826388888889*eccentricity_6;
    double common_term_73 = -15211.783884248279*eccentricity_9 + 15180.904791356647*eccentricity_7;
    double common_term_74 = 36398.464955357143*eccentricity_8;
    double common_term_75 = 83129.35033802389*eccentricity_9;
    double common_term_76 = 763.40796896852093*eccentricity_9;
    double common_term_77 = 448.60078125*eccentricity_8;
    double common_term_78 = 1979.6183682396298*eccentricity_9 + 258.11415705605159*eccentricity_7;
    double common_term_79 = 1204.2793650793651*eccentricity_8 + 144.71805555555556*eccentricity_6;
    double common_term_80 = 3636.8901907784598*eccentricity_9 + 713.45830078125*eccentricity_7 + 78.53203125*eccentricity_5;
    double common_term_81 = 2240.1704861111111*eccentricity_8 + 409.125*eccentricity_6 + 40.833333333333333*eccentricity_4;
    double common_term_82 = 5642.6422087492766*eccentricity_9 + 1337.6403971354167*eccentricity_7 + 225.02473958333333*eccentricity_5 + 20.020833333333333*eccentricity_3;
    double common_term_83 = std::pow(1.0 - eccentricity_2, -9.5)*(1.96875*eccentricity_8 + 19.6875*eccentricity_6 + 31.5*eccentricity_4 + 9.0*eccentricity_2);
    double common_term_84 = 7862.1723083496094*eccentricity_9 + 2069.1467556423611*eccentricity_7 + 418.09114583333333*eccentricity_5 + 56.0625*eccentricity_3 + 3.5*eccentricity;
    double common_term_85 = 4815.6284722222222*eccentricity_8 + 1171.1284722222222*eccentricity_6 + 211.0625*eccentricity_4 + 23.5*eccentricity_2 + 1.0;
    double common_term_86 = 10120.135345458984*eccentricity_9 + 2818.54541015625*eccentricity_7 + 618.2578125*eccentricity_5 + 94.3125*eccentricity_3 + 7.5*eccentricity;
    double common_term_87 = 6101.6444444444444*eccentricity_8 + 1547.75*eccentricity_6 + 292.16666666666667*eccentricity_4 + 33.5*eccentricity_2;
    double common_term_88 = 12203.147214536314*eccentricity_9 + 3471.0574544270833*eccentricity_7 + 768.07682291666667*eccentricity_5 + 115.60416666666667*eccentricity_3;
    double common_term_89 = 7169.06953125*eccentricity_8 + 1795.36875*eccentricity_6 + 339.9375*eccentricity_4;
    double common_term_90 = 13885.421431477865*eccentricity_9 + 3829.5984483506944*eccentricity_7 + 895.01744791666667*eccentricity_5;
    double common_term_91 = 7566.1055555555556*eccentricity_8 + 2172.0972222222222*eccentricity_6;
    double common_term_92 = 13952.05462777274*eccentricity_9 + 4950.4390764508929*eccentricity_7;
    double common_term_93 = 10731.751364087302*eccentricity_8;
    double common_term_94 = 22332.510335668601*eccentricity_9;
    double common_term_95 = 4813.0695709228516*eccentricity_9;
    double common_term_96 = 2552.854191468254*eccentricity_8;
    double common_term_97 = 7470.4682444254557*eccentricity_9 + 1311.7809136284722*eccentricity_7;
    double common_term_98 = 4157.7741071428571*eccentricity_8 + 648.278125*eccentricity_6;
    double common_term_99 = 9865.7245890299479*eccentricity_9 + 2218.4079318576389*eccentricity_7 + 304.97786458333333*eccentricity_5;
    double common_term_100 = 5493.2090277777778*eccentricity_8 + 1124.26875*eccentricity_6 + 134.52083333333333*eccentricity_4;
    double common_term_101 = 11575.071057128906*eccentricity_9 + 2911.84072265625*eccentricity_7 + 533.62890625*eccentricity_5 + 54.3125*eccentricity_3;
    double common_term_102 = 6351.5298611111111*eccentricity_8 + 1449.421875*eccentricity_6 + 231.91666666666667*eccentricity_4 + 19.25*eccentricity_2;
    double common_term_103 = 12464.666538492839*eccentricity_9 + 3281.9704318576389*eccentricity_7 + 662.89322916666667*eccentricity_5 + 88.6875*eccentricity_3 + 5.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(19);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (10, 0).
    // q = -9
    result_by_lpq.set(c_Key3(10, 0, -9), common_term_0);
    result_by_q.set(c_Key1(-9), common_term_0);
    // q = -8
    result_by_lpq.set(c_Key3(10, 0, -8), common_term_1);
    result_by_q.set(c_Key1(-8), common_term_1);
    // q = -7
    result_by_lpq.set(c_Key3(10, 0, -7), common_term_2);
    result_by_q.set(c_Key1(-7), common_term_2);
    // q = -6
    result_by_lpq.set(c_Key3(10, 0, -6), common_term_3);
    result_by_q.set(c_Key1(-6), common_term_3);
    // q = -5
    result_by_lpq.set(c_Key3(10, 0, -5), common_term_4);
    result_by_q.set(c_Key1(-5), common_term_4);
    // q = -4
    result_by_lpq.set(c_Key3(10, 0, -4), common_term_5);
    result_by_q.set(c_Key1(-4), common_term_5);
    // q = -3
    result_by_lpq.set(c_Key3(10, 0, -3), common_term_6);
    result_by_q.set(c_Key1(-3), common_term_6);
    // q = -2
    result_by_lpq.set(c_Key3(10, 0, -2), common_term_7);
    result_by_q.set(c_Key1(-2), common_term_7);
    // q = -1
    result_by_lpq.set(c_Key3(10, 0, -1), common_term_8);
    result_by_q.set(c_Key1(-1), common_term_8);
    // q = 0
    result_by_lpq.set(c_Key3(10, 0, 0), common_term_9);
    result_by_q.set(c_Key1(0), common_term_9);
    // q = 1
    result_by_lpq.set(c_Key3(10, 0, 1), common_term_10);
    result_by_q.set(c_Key1(1), common_term_10);
    // q = 2
    result_by_lpq.set(c_Key3(10, 0, 2), common_term_11);
    result_by_q.set(c_Key1(2), common_term_11);
    // q = 3
    result_by_lpq.set(c_Key3(10, 0, 3), common_term_12);
    result_by_q.set(c_Key1(3), common_term_12);
    // q = 4
    result_by_lpq.set(c_Key3(10, 0, 4), common_term_13);
    result_by_q.set(c_Key1(4), common_term_13);
    // q = 5
    result_by_lpq.set(c_Key3(10, 0, 5), common_term_14);
    result_by_q.set(c_Key1(5), common_term_14);
    // q = 6
    result_by_lpq.set(c_Key3(10, 0, 6), common_term_15);
    result_by_q.set(c_Key1(6), common_term_15);
    // q = 7
    result_by_lpq.set(c_Key3(10, 0, 7), common_term_16);
    result_by_q.set(c_Key1(7), common_term_16);
    // q = 8
    result_by_lpq.set(c_Key3(10, 0, 8), common_term_17);
    result_by_q.set(c_Key1(8), common_term_17);
    // q = 9
    result_by_lpq.set(c_Key3(10, 0, 9), common_term_18);
    result_by_q.set(c_Key1(9), common_term_18);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 0), result_by_q);
    result_by_q.clear();

    // l , p = (10, 1).
    // q = -9
    result_by_lpq.set(c_Key3(10, 1, -9), common_term_19);
    result_by_q.set(c_Key1(-9), common_term_19);
    // q = -8
    result_by_lpq.set(c_Key3(10, 1, -8), common_term_20);
    result_by_q.set(c_Key1(-8), common_term_20);
    // q = -7
    result_by_lpq.set(c_Key3(10, 1, -7), common_term_21);
    result_by_q.set(c_Key1(-7), common_term_21);
    // q = -6
    result_by_lpq.set(c_Key3(10, 1, -6), common_term_22);
    result_by_q.set(c_Key1(-6), common_term_22);
    // q = -5
    result_by_lpq.set(c_Key3(10, 1, -5), common_term_23);
    result_by_q.set(c_Key1(-5), common_term_23);
    // q = -4
    result_by_lpq.set(c_Key3(10, 1, -4), common_term_24);
    result_by_q.set(c_Key1(-4), common_term_24);
    // q = -3
    result_by_lpq.set(c_Key3(10, 1, -3), common_term_25);
    result_by_q.set(c_Key1(-3), common_term_25);
    // q = -2
    result_by_lpq.set(c_Key3(10, 1, -2), common_term_26);
    result_by_q.set(c_Key1(-2), common_term_26);
    // q = -1
    result_by_lpq.set(c_Key3(10, 1, -1), common_term_27);
    result_by_q.set(c_Key1(-1), common_term_27);
    // q = 0
    result_by_lpq.set(c_Key3(10, 1, 0), common_term_28);
    result_by_q.set(c_Key1(0), common_term_28);
    // q = 1
    result_by_lpq.set(c_Key3(10, 1, 1), common_term_29);
    result_by_q.set(c_Key1(1), common_term_29);
    // q = 2
    result_by_lpq.set(c_Key3(10, 1, 2), common_term_30);
    result_by_q.set(c_Key1(2), common_term_30);
    // q = 3
    result_by_lpq.set(c_Key3(10, 1, 3), common_term_31);
    result_by_q.set(c_Key1(3), common_term_31);
    // q = 4
    result_by_lpq.set(c_Key3(10, 1, 4), common_term_32);
    result_by_q.set(c_Key1(4), common_term_32);
    // q = 5
    result_by_lpq.set(c_Key3(10, 1, 5), common_term_33);
    result_by_q.set(c_Key1(5), common_term_33);
    // q = 6
    result_by_lpq.set(c_Key3(10, 1, 6), common_term_34);
    result_by_q.set(c_Key1(6), common_term_34);
    // q = 7
    result_by_lpq.set(c_Key3(10, 1, 7), common_term_35);
    result_by_q.set(c_Key1(7), common_term_35);
    // q = 8
    result_by_lpq.set(c_Key3(10, 1, 8), common_term_36);
    result_by_q.set(c_Key1(8), common_term_36);
    // q = 9
    result_by_lpq.set(c_Key3(10, 1, 9), common_term_37);
    result_by_q.set(c_Key1(9), common_term_37);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 1), result_by_q);
    result_by_q.clear();

    // l , p = (10, 2).
    // q = -9
    result_by_lpq.set(c_Key3(10, 2, -9), common_term_38);
    result_by_q.set(c_Key1(-9), common_term_38);
    // q = -8
    result_by_lpq.set(c_Key3(10, 2, -8), common_term_39);
    result_by_q.set(c_Key1(-8), common_term_39);
    // q = -7
    result_by_lpq.set(c_Key3(10, 2, -7), common_term_40);
    result_by_q.set(c_Key1(-7), common_term_40);
    // q = -6
    result_by_lpq.set(c_Key3(10, 2, -6), common_term_41);
    result_by_q.set(c_Key1(-6), common_term_41);
    // q = -5
    result_by_lpq.set(c_Key3(10, 2, -5), common_term_42);
    result_by_q.set(c_Key1(-5), common_term_42);
    // q = -4
    result_by_lpq.set(c_Key3(10, 2, -4), common_term_43);
    result_by_q.set(c_Key1(-4), common_term_43);
    // q = -3
    result_by_lpq.set(c_Key3(10, 2, -3), common_term_44);
    result_by_q.set(c_Key1(-3), common_term_44);
    // q = -2
    result_by_lpq.set(c_Key3(10, 2, -2), common_term_45);
    result_by_q.set(c_Key1(-2), common_term_45);
    // q = -1
    result_by_lpq.set(c_Key3(10, 2, -1), common_term_46);
    result_by_q.set(c_Key1(-1), common_term_46);
    // q = 0
    result_by_lpq.set(c_Key3(10, 2, 0), common_term_47);
    result_by_q.set(c_Key1(0), common_term_47);
    // q = 1
    result_by_lpq.set(c_Key3(10, 2, 1), common_term_48);
    result_by_q.set(c_Key1(1), common_term_48);
    // q = 2
    result_by_lpq.set(c_Key3(10, 2, 2), common_term_49);
    result_by_q.set(c_Key1(2), common_term_49);
    // q = 3
    result_by_lpq.set(c_Key3(10, 2, 3), common_term_50);
    result_by_q.set(c_Key1(3), common_term_50);
    // q = 4
    result_by_lpq.set(c_Key3(10, 2, 4), common_term_51);
    result_by_q.set(c_Key1(4), common_term_51);
    // q = 5
    result_by_lpq.set(c_Key3(10, 2, 5), common_term_52);
    result_by_q.set(c_Key1(5), common_term_52);
    // q = 6
    result_by_lpq.set(c_Key3(10, 2, 6), common_term_53);
    result_by_q.set(c_Key1(6), common_term_53);
    // q = 7
    result_by_lpq.set(c_Key3(10, 2, 7), common_term_54);
    result_by_q.set(c_Key1(7), common_term_54);
    // q = 8
    result_by_lpq.set(c_Key3(10, 2, 8), common_term_55);
    result_by_q.set(c_Key1(8), common_term_55);
    // q = 9
    result_by_lpq.set(c_Key3(10, 2, 9), common_term_56);
    result_by_q.set(c_Key1(9), common_term_56);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 2), result_by_q);
    result_by_q.clear();

    // l , p = (10, 3).
    // q = -9
    result_by_lpq.set(c_Key3(10, 3, -9), common_term_57);
    result_by_q.set(c_Key1(-9), common_term_57);
    // q = -8
    result_by_lpq.set(c_Key3(10, 3, -8), common_term_58);
    result_by_q.set(c_Key1(-8), common_term_58);
    // q = -7
    result_by_lpq.set(c_Key3(10, 3, -7), common_term_59);
    result_by_q.set(c_Key1(-7), common_term_59);
    // q = -6
    result_by_lpq.set(c_Key3(10, 3, -6), common_term_60);
    result_by_q.set(c_Key1(-6), common_term_60);
    // q = -5
    result_by_lpq.set(c_Key3(10, 3, -5), common_term_61);
    result_by_q.set(c_Key1(-5), common_term_61);
    // q = -4
    result_by_lpq.set(c_Key3(10, 3, -4), common_term_62);
    result_by_q.set(c_Key1(-4), common_term_62);
    // q = -3
    result_by_lpq.set(c_Key3(10, 3, -3), common_term_63);
    result_by_q.set(c_Key1(-3), common_term_63);
    // q = -2
    result_by_lpq.set(c_Key3(10, 3, -2), common_term_64);
    result_by_q.set(c_Key1(-2), common_term_64);
    // q = -1
    result_by_lpq.set(c_Key3(10, 3, -1), common_term_65);
    result_by_q.set(c_Key1(-1), common_term_65);
    // q = 0
    result_by_lpq.set(c_Key3(10, 3, 0), common_term_66);
    result_by_q.set(c_Key1(0), common_term_66);
    // q = 1
    result_by_lpq.set(c_Key3(10, 3, 1), common_term_67);
    result_by_q.set(c_Key1(1), common_term_67);
    // q = 2
    result_by_lpq.set(c_Key3(10, 3, 2), common_term_68);
    result_by_q.set(c_Key1(2), common_term_68);
    // q = 3
    result_by_lpq.set(c_Key3(10, 3, 3), common_term_69);
    result_by_q.set(c_Key1(3), common_term_69);
    // q = 4
    result_by_lpq.set(c_Key3(10, 3, 4), common_term_70);
    result_by_q.set(c_Key1(4), common_term_70);
    // q = 5
    result_by_lpq.set(c_Key3(10, 3, 5), common_term_71);
    result_by_q.set(c_Key1(5), common_term_71);
    // q = 6
    result_by_lpq.set(c_Key3(10, 3, 6), common_term_72);
    result_by_q.set(c_Key1(6), common_term_72);
    // q = 7
    result_by_lpq.set(c_Key3(10, 3, 7), common_term_73);
    result_by_q.set(c_Key1(7), common_term_73);
    // q = 8
    result_by_lpq.set(c_Key3(10, 3, 8), common_term_74);
    result_by_q.set(c_Key1(8), common_term_74);
    // q = 9
    result_by_lpq.set(c_Key3(10, 3, 9), common_term_75);
    result_by_q.set(c_Key1(9), common_term_75);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 3), result_by_q);
    result_by_q.clear();

    // l , p = (10, 4).
    // q = -9
    result_by_lpq.set(c_Key3(10, 4, -9), common_term_76);
    result_by_q.set(c_Key1(-9), common_term_76);
    // q = -8
    result_by_lpq.set(c_Key3(10, 4, -8), common_term_77);
    result_by_q.set(c_Key1(-8), common_term_77);
    // q = -7
    result_by_lpq.set(c_Key3(10, 4, -7), common_term_78);
    result_by_q.set(c_Key1(-7), common_term_78);
    // q = -6
    result_by_lpq.set(c_Key3(10, 4, -6), common_term_79);
    result_by_q.set(c_Key1(-6), common_term_79);
    // q = -5
    result_by_lpq.set(c_Key3(10, 4, -5), common_term_80);
    result_by_q.set(c_Key1(-5), common_term_80);
    // q = -4
    result_by_lpq.set(c_Key3(10, 4, -4), common_term_81);
    result_by_q.set(c_Key1(-4), common_term_81);
    // q = -3
    result_by_lpq.set(c_Key3(10, 4, -3), common_term_82);
    result_by_q.set(c_Key1(-3), common_term_82);
    // q = -2
    result_by_lpq.set(c_Key3(10, 4, -2), common_term_83);
    result_by_q.set(c_Key1(-2), common_term_83);
    // q = -1
    result_by_lpq.set(c_Key3(10, 4, -1), common_term_84);
    result_by_q.set(c_Key1(-1), common_term_84);
    // q = 0
    result_by_lpq.set(c_Key3(10, 4, 0), common_term_85);
    result_by_q.set(c_Key1(0), common_term_85);
    // q = 1
    result_by_lpq.set(c_Key3(10, 4, 1), common_term_86);
    result_by_q.set(c_Key1(1), common_term_86);
    // q = 2
    result_by_lpq.set(c_Key3(10, 4, 2), common_term_87);
    result_by_q.set(c_Key1(2), common_term_87);
    // q = 3
    result_by_lpq.set(c_Key3(10, 4, 3), common_term_88);
    result_by_q.set(c_Key1(3), common_term_88);
    // q = 4
    result_by_lpq.set(c_Key3(10, 4, 4), common_term_89);
    result_by_q.set(c_Key1(4), common_term_89);
    // q = 5
    result_by_lpq.set(c_Key3(10, 4, 5), common_term_90);
    result_by_q.set(c_Key1(5), common_term_90);
    // q = 6
    result_by_lpq.set(c_Key3(10, 4, 6), common_term_91);
    result_by_q.set(c_Key1(6), common_term_91);
    // q = 7
    result_by_lpq.set(c_Key3(10, 4, 7), common_term_92);
    result_by_q.set(c_Key1(7), common_term_92);
    // q = 8
    result_by_lpq.set(c_Key3(10, 4, 8), common_term_93);
    result_by_q.set(c_Key1(8), common_term_93);
    // q = 9
    result_by_lpq.set(c_Key3(10, 4, 9), common_term_94);
    result_by_q.set(c_Key1(9), common_term_94);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 4), result_by_q);
    result_by_q.clear();

    // l , p = (10, 5).
    // q = -9
    result_by_lpq.set(c_Key3(10, 5, -9), common_term_95);
    result_by_q.set(c_Key1(-9), common_term_95);
    // q = -8
    result_by_lpq.set(c_Key3(10, 5, -8), common_term_96);
    result_by_q.set(c_Key1(-8), common_term_96);
    // q = -7
    result_by_lpq.set(c_Key3(10, 5, -7), common_term_97);
    result_by_q.set(c_Key1(-7), common_term_97);
    // q = -6
    result_by_lpq.set(c_Key3(10, 5, -6), common_term_98);
    result_by_q.set(c_Key1(-6), common_term_98);
    // q = -5
    result_by_lpq.set(c_Key3(10, 5, -5), common_term_99);
    result_by_q.set(c_Key1(-5), common_term_99);
    // q = -4
    result_by_lpq.set(c_Key3(10, 5, -4), common_term_100);
    result_by_q.set(c_Key1(-4), common_term_100);
    // q = -3
    result_by_lpq.set(c_Key3(10, 5, -3), common_term_101);
    result_by_q.set(c_Key1(-3), common_term_101);
    // q = -2
    result_by_lpq.set(c_Key3(10, 5, -2), common_term_102);
    result_by_q.set(c_Key1(-2), common_term_102);
    // q = -1
    result_by_lpq.set(c_Key3(10, 5, -1), common_term_103);
    result_by_q.set(c_Key1(-1), common_term_103);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -9.5)*(2.4609375*eccentricity_8 + 26.25*eccentricity_6 + 47.25*eccentricity_4 + 18.0*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(10, 5, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(10, 5, 1), common_term_103);
    result_by_q.set(c_Key1(1), common_term_103);
    // q = 2
    result_by_lpq.set(c_Key3(10, 5, 2), common_term_102);
    result_by_q.set(c_Key1(2), common_term_102);
    // q = 3
    result_by_lpq.set(c_Key3(10, 5, 3), common_term_101);
    result_by_q.set(c_Key1(3), common_term_101);
    // q = 4
    result_by_lpq.set(c_Key3(10, 5, 4), common_term_100);
    result_by_q.set(c_Key1(4), common_term_100);
    // q = 5
    result_by_lpq.set(c_Key3(10, 5, 5), common_term_99);
    result_by_q.set(c_Key1(5), common_term_99);
    // q = 6
    result_by_lpq.set(c_Key3(10, 5, 6), common_term_98);
    result_by_q.set(c_Key1(6), common_term_98);
    // q = 7
    result_by_lpq.set(c_Key3(10, 5, 7), common_term_97);
    result_by_q.set(c_Key1(7), common_term_97);
    // q = 8
    result_by_lpq.set(c_Key3(10, 5, 8), common_term_96);
    result_by_q.set(c_Key1(8), common_term_96);
    // q = 9
    result_by_lpq.set(c_Key3(10, 5, 9), common_term_95);
    result_by_q.set(c_Key1(9), common_term_95);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 5), result_by_q);
    result_by_q.clear();

    // l , p = (10, 6).
    // q = -9
    result_by_lpq.set(c_Key3(10, 6, -9), common_term_94);
    result_by_q.set(c_Key1(-9), common_term_94);
    // q = -8
    result_by_lpq.set(c_Key3(10, 6, -8), common_term_93);
    result_by_q.set(c_Key1(-8), common_term_93);
    // q = -7
    result_by_lpq.set(c_Key3(10, 6, -7), common_term_92);
    result_by_q.set(c_Key1(-7), common_term_92);
    // q = -6
    result_by_lpq.set(c_Key3(10, 6, -6), common_term_91);
    result_by_q.set(c_Key1(-6), common_term_91);
    // q = -5
    result_by_lpq.set(c_Key3(10, 6, -5), common_term_90);
    result_by_q.set(c_Key1(-5), common_term_90);
    // q = -4
    result_by_lpq.set(c_Key3(10, 6, -4), common_term_89);
    result_by_q.set(c_Key1(-4), common_term_89);
    // q = -3
    result_by_lpq.set(c_Key3(10, 6, -3), common_term_88);
    result_by_q.set(c_Key1(-3), common_term_88);
    // q = -2
    result_by_lpq.set(c_Key3(10, 6, -2), common_term_87);
    result_by_q.set(c_Key1(-2), common_term_87);
    // q = -1
    result_by_lpq.set(c_Key3(10, 6, -1), common_term_86);
    result_by_q.set(c_Key1(-1), common_term_86);
    // q = 0
    result_by_lpq.set(c_Key3(10, 6, 0), common_term_85);
    result_by_q.set(c_Key1(0), common_term_85);
    // q = 1
    result_by_lpq.set(c_Key3(10, 6, 1), common_term_84);
    result_by_q.set(c_Key1(1), common_term_84);
    // q = 2
    result_by_lpq.set(c_Key3(10, 6, 2), common_term_83);
    result_by_q.set(c_Key1(2), common_term_83);
    // q = 3
    result_by_lpq.set(c_Key3(10, 6, 3), common_term_82);
    result_by_q.set(c_Key1(3), common_term_82);
    // q = 4
    result_by_lpq.set(c_Key3(10, 6, 4), common_term_81);
    result_by_q.set(c_Key1(4), common_term_81);
    // q = 5
    result_by_lpq.set(c_Key3(10, 6, 5), common_term_80);
    result_by_q.set(c_Key1(5), common_term_80);
    // q = 6
    result_by_lpq.set(c_Key3(10, 6, 6), common_term_79);
    result_by_q.set(c_Key1(6), common_term_79);
    // q = 7
    result_by_lpq.set(c_Key3(10, 6, 7), common_term_78);
    result_by_q.set(c_Key1(7), common_term_78);
    // q = 8
    result_by_lpq.set(c_Key3(10, 6, 8), common_term_77);
    result_by_q.set(c_Key1(8), common_term_77);
    // q = 9
    result_by_lpq.set(c_Key3(10, 6, 9), common_term_76);
    result_by_q.set(c_Key1(9), common_term_76);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 6), result_by_q);
    result_by_q.clear();

    // l , p = (10, 7).
    // q = -9
    result_by_lpq.set(c_Key3(10, 7, -9), common_term_75);
    result_by_q.set(c_Key1(-9), common_term_75);
    // q = -8
    result_by_lpq.set(c_Key3(10, 7, -8), common_term_74);
    result_by_q.set(c_Key1(-8), common_term_74);
    // q = -7
    result_by_lpq.set(c_Key3(10, 7, -7), common_term_73);
    result_by_q.set(c_Key1(-7), common_term_73);
    // q = -6
    result_by_lpq.set(c_Key3(10, 7, -6), common_term_72);
    result_by_q.set(c_Key1(-6), common_term_72);
    // q = -5
    result_by_lpq.set(c_Key3(10, 7, -5), common_term_71);
    result_by_q.set(c_Key1(-5), common_term_71);
    // q = -4
    result_by_lpq.set(c_Key3(10, 7, -4), common_term_70);
    result_by_q.set(c_Key1(-4), common_term_70);
    // q = -3
    result_by_lpq.set(c_Key3(10, 7, -3), common_term_69);
    result_by_q.set(c_Key1(-3), common_term_69);
    // q = -2
    result_by_lpq.set(c_Key3(10, 7, -2), common_term_68);
    result_by_q.set(c_Key1(-2), common_term_68);
    // q = -1
    result_by_lpq.set(c_Key3(10, 7, -1), common_term_67);
    result_by_q.set(c_Key1(-1), common_term_67);
    // q = 0
    result_by_lpq.set(c_Key3(10, 7, 0), common_term_66);
    result_by_q.set(c_Key1(0), common_term_66);
    // q = 1
    result_by_lpq.set(c_Key3(10, 7, 1), common_term_65);
    result_by_q.set(c_Key1(1), common_term_65);
    // q = 2
    result_by_lpq.set(c_Key3(10, 7, 2), common_term_64);
    result_by_q.set(c_Key1(2), common_term_64);
    // q = 3
    result_by_lpq.set(c_Key3(10, 7, 3), common_term_63);
    result_by_q.set(c_Key1(3), common_term_63);
    // q = 4
    result_by_lpq.set(c_Key3(10, 7, 4), common_term_62);
    result_by_q.set(c_Key1(4), common_term_62);
    // q = 5
    result_by_lpq.set(c_Key3(10, 7, 5), common_term_61);
    result_by_q.set(c_Key1(5), common_term_61);
    // q = 6
    result_by_lpq.set(c_Key3(10, 7, 6), common_term_60);
    result_by_q.set(c_Key1(6), common_term_60);
    // q = 7
    result_by_lpq.set(c_Key3(10, 7, 7), common_term_59);
    result_by_q.set(c_Key1(7), common_term_59);
    // q = 8
    result_by_lpq.set(c_Key3(10, 7, 8), common_term_58);
    result_by_q.set(c_Key1(8), common_term_58);
    // q = 9
    result_by_lpq.set(c_Key3(10, 7, 9), common_term_57);
    result_by_q.set(c_Key1(9), common_term_57);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 7), result_by_q);
    result_by_q.clear();

    // l , p = (10, 8).
    // q = -9
    result_by_lpq.set(c_Key3(10, 8, -9), common_term_56);
    result_by_q.set(c_Key1(-9), common_term_56);
    // q = -8
    result_by_lpq.set(c_Key3(10, 8, -8), common_term_55);
    result_by_q.set(c_Key1(-8), common_term_55);
    // q = -7
    result_by_lpq.set(c_Key3(10, 8, -7), common_term_54);
    result_by_q.set(c_Key1(-7), common_term_54);
    // q = -6
    result_by_lpq.set(c_Key3(10, 8, -6), common_term_53);
    result_by_q.set(c_Key1(-6), common_term_53);
    // q = -5
    result_by_lpq.set(c_Key3(10, 8, -5), common_term_52);
    result_by_q.set(c_Key1(-5), common_term_52);
    // q = -4
    result_by_lpq.set(c_Key3(10, 8, -4), common_term_51);
    result_by_q.set(c_Key1(-4), common_term_51);
    // q = -3
    result_by_lpq.set(c_Key3(10, 8, -3), common_term_50);
    result_by_q.set(c_Key1(-3), common_term_50);
    // q = -2
    result_by_lpq.set(c_Key3(10, 8, -2), common_term_49);
    result_by_q.set(c_Key1(-2), common_term_49);
    // q = -1
    result_by_lpq.set(c_Key3(10, 8, -1), common_term_48);
    result_by_q.set(c_Key1(-1), common_term_48);
    // q = 0
    result_by_lpq.set(c_Key3(10, 8, 0), common_term_47);
    result_by_q.set(c_Key1(0), common_term_47);
    // q = 1
    result_by_lpq.set(c_Key3(10, 8, 1), common_term_46);
    result_by_q.set(c_Key1(1), common_term_46);
    // q = 2
    result_by_lpq.set(c_Key3(10, 8, 2), common_term_45);
    result_by_q.set(c_Key1(2), common_term_45);
    // q = 3
    result_by_lpq.set(c_Key3(10, 8, 3), common_term_44);
    result_by_q.set(c_Key1(3), common_term_44);
    // q = 4
    result_by_lpq.set(c_Key3(10, 8, 4), common_term_43);
    result_by_q.set(c_Key1(4), common_term_43);
    // q = 5
    result_by_lpq.set(c_Key3(10, 8, 5), common_term_42);
    result_by_q.set(c_Key1(5), common_term_42);
    // q = 6
    result_by_lpq.set(c_Key3(10, 8, 6), common_term_41);
    result_by_q.set(c_Key1(6), common_term_41);
    // q = 7
    result_by_lpq.set(c_Key3(10, 8, 7), common_term_40);
    result_by_q.set(c_Key1(7), common_term_40);
    // q = 8
    result_by_lpq.set(c_Key3(10, 8, 8), common_term_39);
    result_by_q.set(c_Key1(8), common_term_39);
    // q = 9
    result_by_lpq.set(c_Key3(10, 8, 9), common_term_38);
    result_by_q.set(c_Key1(9), common_term_38);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 8), result_by_q);
    result_by_q.clear();

    // l , p = (10, 9).
    // q = -9
    result_by_lpq.set(c_Key3(10, 9, -9), common_term_37);
    result_by_q.set(c_Key1(-9), common_term_37);
    // q = -8
    result_by_lpq.set(c_Key3(10, 9, -8), common_term_36);
    result_by_q.set(c_Key1(-8), common_term_36);
    // q = -7
    result_by_lpq.set(c_Key3(10, 9, -7), common_term_35);
    result_by_q.set(c_Key1(-7), common_term_35);
    // q = -6
    result_by_lpq.set(c_Key3(10, 9, -6), common_term_34);
    result_by_q.set(c_Key1(-6), common_term_34);
    // q = -5
    result_by_lpq.set(c_Key3(10, 9, -5), common_term_33);
    result_by_q.set(c_Key1(-5), common_term_33);
    // q = -4
    result_by_lpq.set(c_Key3(10, 9, -4), common_term_32);
    result_by_q.set(c_Key1(-4), common_term_32);
    // q = -3
    result_by_lpq.set(c_Key3(10, 9, -3), common_term_31);
    result_by_q.set(c_Key1(-3), common_term_31);
    // q = -2
    result_by_lpq.set(c_Key3(10, 9, -2), common_term_30);
    result_by_q.set(c_Key1(-2), common_term_30);
    // q = -1
    result_by_lpq.set(c_Key3(10, 9, -1), common_term_29);
    result_by_q.set(c_Key1(-1), common_term_29);
    // q = 0
    result_by_lpq.set(c_Key3(10, 9, 0), common_term_28);
    result_by_q.set(c_Key1(0), common_term_28);
    // q = 1
    result_by_lpq.set(c_Key3(10, 9, 1), common_term_27);
    result_by_q.set(c_Key1(1), common_term_27);
    // q = 2
    result_by_lpq.set(c_Key3(10, 9, 2), common_term_26);
    result_by_q.set(c_Key1(2), common_term_26);
    // q = 3
    result_by_lpq.set(c_Key3(10, 9, 3), common_term_25);
    result_by_q.set(c_Key1(3), common_term_25);
    // q = 4
    result_by_lpq.set(c_Key3(10, 9, 4), common_term_24);
    result_by_q.set(c_Key1(4), common_term_24);
    // q = 5
    result_by_lpq.set(c_Key3(10, 9, 5), common_term_23);
    result_by_q.set(c_Key1(5), common_term_23);
    // q = 6
    result_by_lpq.set(c_Key3(10, 9, 6), common_term_22);
    result_by_q.set(c_Key1(6), common_term_22);
    // q = 7
    result_by_lpq.set(c_Key3(10, 9, 7), common_term_21);
    result_by_q.set(c_Key1(7), common_term_21);
    // q = 8
    result_by_lpq.set(c_Key3(10, 9, 8), common_term_20);
    result_by_q.set(c_Key1(8), common_term_20);
    // q = 9
    result_by_lpq.set(c_Key3(10, 9, 9), common_term_19);
    result_by_q.set(c_Key1(9), common_term_19);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 9), result_by_q);
    result_by_q.clear();

    // l , p = (10, 10).
    // q = -9
    result_by_lpq.set(c_Key3(10, 10, -9), common_term_18);
    result_by_q.set(c_Key1(-9), common_term_18);
    // q = -8
    result_by_lpq.set(c_Key3(10, 10, -8), common_term_17);
    result_by_q.set(c_Key1(-8), common_term_17);
    // q = -7
    result_by_lpq.set(c_Key3(10, 10, -7), common_term_16);
    result_by_q.set(c_Key1(-7), common_term_16);
    // q = -6
    result_by_lpq.set(c_Key3(10, 10, -6), common_term_15);
    result_by_q.set(c_Key1(-6), common_term_15);
    // q = -5
    result_by_lpq.set(c_Key3(10, 10, -5), common_term_14);
    result_by_q.set(c_Key1(-5), common_term_14);
    // q = -4
    result_by_lpq.set(c_Key3(10, 10, -4), common_term_13);
    result_by_q.set(c_Key1(-4), common_term_13);
    // q = -3
    result_by_lpq.set(c_Key3(10, 10, -3), common_term_12);
    result_by_q.set(c_Key1(-3), common_term_12);
    // q = -2
    result_by_lpq.set(c_Key3(10, 10, -2), common_term_11);
    result_by_q.set(c_Key1(-2), common_term_11);
    // q = -1
    result_by_lpq.set(c_Key3(10, 10, -1), common_term_10);
    result_by_q.set(c_Key1(-1), common_term_10);
    // q = 0
    result_by_lpq.set(c_Key3(10, 10, 0), common_term_9);
    result_by_q.set(c_Key1(0), common_term_9);
    // q = 1
    result_by_lpq.set(c_Key3(10, 10, 1), common_term_8);
    result_by_q.set(c_Key1(1), common_term_8);
    // q = 2
    result_by_lpq.set(c_Key3(10, 10, 2), common_term_7);
    result_by_q.set(c_Key1(2), common_term_7);
    // q = 3
    result_by_lpq.set(c_Key3(10, 10, 3), common_term_6);
    result_by_q.set(c_Key1(3), common_term_6);
    // q = 4
    result_by_lpq.set(c_Key3(10, 10, 4), common_term_5);
    result_by_q.set(c_Key1(4), common_term_5);
    // q = 5
    result_by_lpq.set(c_Key3(10, 10, 5), common_term_4);
    result_by_q.set(c_Key1(5), common_term_4);
    // q = 6
    result_by_lpq.set(c_Key3(10, 10, 6), common_term_3);
    result_by_q.set(c_Key1(6), common_term_3);
    // q = 7
    result_by_lpq.set(c_Key3(10, 10, 7), common_term_2);
    result_by_q.set(c_Key1(7), common_term_2);
    // q = 8
    result_by_lpq.set(c_Key3(10, 10, 8), common_term_1);
    result_by_q.set(c_Key1(8), common_term_1);
    // q = 9
    result_by_lpq.set(c_Key3(10, 10, 9), common_term_0);
    result_by_q.set(c_Key1(9), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 10), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l10_e15(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 10.
    // Functions are taylor expanded and truncated at eccentricity^15.

    // Functions Calculated for l = 10.

    c_IntMap<c_Key3, double> result_by_lpq(317);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(11);
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
    double common_term_0 = 1.8793669587320381e-7*eccentricity_14;
    double common_term_1 = 3.1254031917789242e-8*eccentricity_13;
    double common_term_2 = 6.6645031922809701e-9*eccentricity_14 + 2.0876756987868099e-9*eccentricity_12;
    double common_term_3 = 3.5423208267989084e-11*eccentricity_13 + 1.2232474797578964e-11*eccentricity_11;
    double common_term_4 = -1.4833404751414192e-8*eccentricity_13 - 1.063002059909612e-8*eccentricity_11 - 5.3822889109347443e-9*eccentricity_9;
    double common_term_5 = 3.6877225462294907e-5*eccentricity_14 + 3.5445601851851852e-5*eccentricity_12 + 3.1690917107583774e-5*eccentricity_10 + 2.4801587301587302e-5*eccentricity_8;
    double common_term_6 = -0.00093127523149762835*eccentricity_13 - 0.0011785779680524554*eccentricity_11 - 0.0011653355189732143*eccentricity_9 - 0.0033900669642857143*eccentricity_7;
    double common_term_7 = -0.0069032921810699493*eccentricity_14 - 0.011360964138741917*eccentricity_12 + 0.0087301587301587302*eccentricity_10 - 0.08253968253968254*eccentricity_8 + 0.088888888888888889*eccentricity_6;
    double common_term_8 = -0.14282798045246444*eccentricity_13 + 0.49970137379156849*eccentricity_11 - 1.6257876441592262*eccentricity_9 + 2.2040473090277778*eccentricity_7 - 0.81380208333333333*eccentricity_5;
    double common_term_9 = -2.40773856026784*eccentricity_14 + 8.9207728794642912*eccentricity_12 - 22.639620535714286*eccentricity_10 + 30.9234375*eccentricity_8 - 17.8875*eccentricity_6 + 3.375*eccentricity_4;
    double common_term_10 = 105.30185376997335*eccentricity_13 - 232.47241801509155*eccentricity_11 + 308.33914161964699*eccentricity_9 - 212.15309244791667*eccentricity_7 + 66.545572916666667*eccentricity_5 - 7.1458333333333333*eccentricity_3;
    double common_term_11 = 946.38402626212531*eccentricity_14 - 1881.2400396825392*eccentricity_12 + 2402.6212962962965*eccentricity_10 - 1784.6055555555556*eccentricity_8 + 704.5*eccentricity_6 - 129.33333333333333*eccentricity_4 + 8.0*eccentricity_2;
    double common_term_12 = -12623.381456181661*eccentricity_13 + 15486.722428436279*eccentricity_11 - 11881.492437744141*eccentricity_9 + 5306.87548828125*eccentricity_7 - 1260.4921875*eccentricity_5 + 135.5625*eccentricity_3 - 4.5*eccentricity;
    double common_term_13 = -72857.953256235602*eccentricity_14 + 85822.01504931037*eccentricity_12 - 66460.387586805555*eccentricity_10 + 31844.422743055556*eccentricity_8 - 8775.9548611111111*eccentricity_6 + 1241.5625*eccentricity_4 - 72.5*eccentricity_2 + 1.0;
    double common_term_14 = 420362.94863797451*eccentricity_13 - 324311.01459763702*eccentricity_11 + 161721.53108113607*eccentricity_9 - 48884.568739149307*eccentricity_7 + 8150.3411458333333*eccentricity_5 - 633.1875*eccentricity_3 + 15.5*eccentricity;
    double common_term_15 = 1857749.3601102121*eccentricity_14 - 1416734.9751897321*eccentricity_12 + 721640.82304687498*eccentricity_10 - 231454.8*eccentricity_8 + 42937.312500000001*eccentricity_6 - 3973.5*eccentricity_4 + 130.5*eccentricity_2;
    double common_term_16 = -5645298.2614953222*eccentricity_13 + 2901880.221849056*eccentricity_11 - 966863.29034333404*eccentricity_9 + 192809.99469401042*eccentricity_7 - 20053.194010416667*eccentricity_5 + 790.77083333333333*eccentricity_3;
    double common_term_17 = -20811417.030628115*eccentricity_14 + 10707868.584195964*eccentricity_12 - 3653666.6625723379*eccentricity_10 + 765856.92907986111*eccentricity_8 - 86424.668750000001*eccentricity_6 + 3858.8958333333334*eccentricity_4;
    double common_term_18 = 36749951.747007455*eccentricity_13 - 12713999.931560244*eccentricity_11 + 2758410.2350289481*eccentricity_9 - 330019.78271484375*eccentricity_7 + 16100.61328125*eccentricity_5;
    double common_term_19 = 118538338.34825096*eccentricity_14 - 41281660.227906379*eccentricity_12 + 9168641.8137896825*eccentricity_10 - 1144440.9515873016*eccentricity_8 + 59584.609722222222*eccentricity_6;
    double common_term_20 = -126341976.86238201*eccentricity_13 + 28493340.940017857*eccentricity_11 - 3667624.6796146938*eccentricity_9 + 200445.43370690724*eccentricity_7;
    double common_term_21 = -367374606.24462491*eccentricity_14 + 83618046.254594029*eccentricity_12 - 11003358.154520089*eccentricity_10 + 623706.16439732142*eccentricity_8;
    double common_term_22 = 233544169.4545519*eccentricity_13 - 31210475.019879966*eccentricity_11 + 1818308.9840970532*eccentricity_9;
    double common_term_23 = 624697168.85578806*eccentricity_14 - 84346904.447915727*eccentricity_12 + 5015578.577405065*eccentricity_10;
    double common_term_24 = -218534200.51042232*eccentricity_13 + 13191026.055637727*eccentricity_11;
    double common_term_25 = -545556461.6433747*eccentricity_14 + 33282081.127764854*eccentricity_12;
    double common_term_26 = 80963575.707463975*eccentricity_13;
    double common_term_27 = 190682398.33831324*eccentricity_14;
    double common_term_28 = 0.22162890785274993*eccentricity_14;
    double common_term_29 = 0.16305740723786139*eccentricity_13;
    double common_term_30 = 1.0197201563667492*eccentricity_14 + 0.11996800584799717*eccentricity_12;
    double common_term_31 = 0.77233257646684523*eccentricity_13 + 0.088266719471324574*eccentricity_11;
    double common_term_32 = 2.9346062504018776*eccentricity_14 + 0.58448591194684945*eccentricity_12 + 0.064942887455908289*eccentricity_10;
    double common_term_33 = 2.2696591440821298*eccentricity_13 + 0.44198625446837625*eccentricity_11 + 0.047782297824005181*eccentricity_9;
    double common_term_34 = 0.03515625*eccentricity_8*std::pow(1.0 - eccentricity_2, -9.5);
    double common_term_35 = 5.2837762153744347*eccentricity_13 + 1.3531424621635854*eccentricity_11 + 0.25219857352120536*eccentricity_9 + 0.025866505456349206*eccentricity_7;
    double common_term_36 = 13.381717950550687*eccentricity_14 + 4.1485331468621399*eccentricity_12 + 1.043226066468254*eccentricity_10 + 0.19037698412698413*eccentricity_8 + 0.019097222222222222*eccentricity_6;
    double common_term_37 = 10.65913347857339*eccentricity_13 + 3.2524068014962333*eccentricity_11 + 0.80217459542410714*eccentricity_9 + 0.14326171875*eccentricity_7 + 0.00703125*eccentricity_5;
    double common_term_38 = 24.189023784837044*eccentricity_14 + 8.4745465959821429*eccentricity_12 + 2.5287037037037037*eccentricity_10 + 0.66111111111111111*eccentricity_8 - 0.06875*eccentricity_6 + 0.14583333333333333*eccentricity_4;
    double common_term_39 = 19.503868932446475*eccentricity_13 + 6.4177327433591166*eccentricity_11 + 3.101062915943287*eccentricity_9 - 2.3777669270833333*eccentricity_7 + 3.0247395833333333*eccentricity_5 - 0.85416666666666667*eccentricity_3;
    double common_term_40 = 41.58755545479911*eccentricity_14 + 11.053727678571429*eccentricity_12 + 20.2921875*eccentricity_10 - 30.20625*eccentricity_8 + 35.015625*eccentricity_6 - 15.75*eccentricity_4 + 2.25*eccentricity_2;
    double common_term_41 = -16.790848970530973*eccentricity_13 + 153.76848819873951*eccentricity_11 - 260.62958780924479*eccentricity_9 + 288.20958116319444*eccentricity_7 - 155.44010416666667*eccentricity_5 + 36.1875*eccentricity_3 - 2.5*eccentricity;
    double common_term_42 = -352.56799445704168*eccentricity_14 + 1067.3436419753088*eccentricity_12 - 1758.5455555555556*eccentricity_10 + 1866.3055555555556*eccentricity_8 - 1089.5486111111111*eccentricity_6 + 315.125*eccentricity_4 - 36.5*eccentricity_2 + 1.0;
    double common_term_43 = 6407.1587927218836*eccentricity_13 - 9903.7892381286627*eccentricity_11 + 10084.092352294922*eccentricity_9 - 6078.41455078125*eccentricity_7 + 1970.2265625*eccentricity_5 - 290.8125*eccentricity_3 + 13.5*eccentricity;
    double common_term_44 = 33541.5456548351*eccentricity_14 - 48417.402591765872*eccentricity_12 + 47262.871817129631*eccentricity_10 - 28695.711805555556*eccentricity_8 + 9901.328125*eccentricity_6 - 1675.0833333333333*eccentricity_4 + 100.25*eccentricity_2;
    double common_term_45 = -211038.32953387519*eccentricity_13 + 197461.10533111007*eccentricity_11 - 119081.11242088035*eccentricity_9 + 42477.821712239583*eccentricity_7 - 7802.2669270833333*eccentricity_5 + 541.47916666666667*eccentricity_3;
    double common_term_46 = -836244.99426269532*eccentricity_14 + 750256.73368443079*eccentricity_12 - 445681.59910714286*eccentricity_10 + 161427.2625*eccentricity_8 - 31188.09375*eccentricity_6 + 2376.5625*eccentricity_4;
    double common_term_47 = 2631957.8339697471*eccentricity_13 - 1532288.5636892672*eccentricity_11 + 557026.68569394067*eccentricity_9 - 110947.63989800347*eccentricity_7 + 8987.4486979166667*eccentricity_5;
    double common_term_48 = 8625578.5376998437*eccentricity_14 - 4906074.2055129565*eccentricity_12 + 1776120.1692491319*eccentricity_10 - 359842.25486111111*eccentricity_8 + 30349.348263888889*eccentricity_6;
    double common_term_49 = -14783356.235550846*eccentricity_13 + 5301464.0178067343*eccentricity_11 - 1082391.1435743059*eccentricity_9 + 93710.231794084822*eccentricity_7;
    double common_term_50 = -42272175.433449653*eccentricity_14 + 14960518.251003327*eccentricity_12 - 3057694.4276902833*eccentricity_10 + 269031.71821676587*eccentricity_8;
    double common_term_51 = 40224855.413677096*eccentricity_13 - 8190235.4274458904*eccentricity_11 + 726995.90431837725*eccentricity_9;
    double common_term_52 = 103689654.78311793*eccentricity_14 - 20957317.37134791*eccentricity_12 + 1866483.2682421875*eccentricity_10;
    double common_term_53 = -51534719.673596271*eccentricity_13 + 4586053.3352122621*eccentricity_11;
    double common_term_54 = -122374212.7556137*eccentricity_14 + 10846586.516443591*eccentricity_12;
    double common_term_55 = 24809937.971262643*eccentricity_13;
    double common_term_56 = 55095376.82915985*eccentricity_14;
    double common_term_57 = 20.864082391798476*eccentricity_14;
    double common_term_58 = 14.837341078135081*eccentricity_13;
    double common_term_59 = 82.685478052221216*eccentricity_14 + 10.538954317167208*eccentricity_12;
    double common_term_60 = 60.899847742333952*eccentricity_13 + 7.4762625074241336*eccentricity_11;
    double common_term_61 = 213.86499566859484*eccentricity_14 + 44.752017909752285*eccentricity_12 + 5.2963146219135802*eccentricity_10;
    double common_term_62 = 161.28867905393823*eccentricity_13 + 32.810629239763532*eccentricity_11 + 3.7464246477399554*eccentricity_9;
    double common_term_63 = 450.6495398614013*eccentricity_14 + 121.32952799479167*eccentricity_12 + 24.000000688932981*eccentricity_10 + 2.6458209325396825*eccentricity_8;
    double common_term_64 = 345.99517598736363*eccentricity_13 + 91.03766524770693*eccentricity_11 + 17.51398184640067*eccentricity_9 + 1.8652793278769841*eccentricity_7;
    double common_term_65 = std::pow(1.0 - eccentricity_2, -9.5)*(0.28125*eccentricity_8 + 1.3125*eccentricity_6);
    double common_term_66 = 652.62349498479385*eccentricity_13 + 202.37828039946379*eccentricity_11 + 50.857536388578869*eccentricity_9 + 9.2587782118055556*eccentricity_7 + 0.92161458333333333*eccentricity_5;
    double common_term_67 = 1430.0388710656186*eccentricity_14 + 507.00156637524802*eccentricity_12 + 154.17078579695767*eccentricity_10 + 37.861371527777778*eccentricity_8 + 6.70625*eccentricity_6 + 0.64583333333333333*eccentricity_4;
    double common_term_68 = 1127.1018888698305*eccentricity_13 + 392.86032180786133*eccentricity_11 + 117.13199462890625*eccentricity_9 + 28.10771484375*eccentricity_7 + 4.84765625*eccentricity_5 + 0.4375*eccentricity_3;
    double common_term_69 = 2291.4185359657463*eccentricity_14 + 886.1375037202381*eccentricity_12 + 303.62125289351852*eccentricity_10 + 88.711111111111111*eccentricity_8 + 20.9375*eccentricity_6 + 3.1666666666666667*eccentricity_4 + 0.5*eccentricity_2;
    double common_term_70 = 1824.3681644148717*eccentricity_13 + 695.03978884661639*eccentricity_11 + 233.38388061523438*eccentricity_9 + 69.138726128472222*eccentricity_7 + 11.174479166666667*eccentricity_5 + 5.8125*eccentricity_3 - 0.5*eccentricity;
    double common_term_71 = 3497.1331227080676*eccentricity_14 + 1450.773330078125*eccentricity_12 + 536.04703125*eccentricity_10 + 201.1015625*eccentricity_8 + 13.65625*eccentricity_6 + 42.5625*eccentricity_4 - 8.5*eccentricity_2 + 1.0;
    double common_term_72 = 2827.4208128707792*eccentricity_13 + 1083.1575818831832*eccentricity_11 + 582.28132527669271*eccentricity_9 - 107.68614366319444*eccentricity_7 + 245.67447916666667*eccentricity_5 - 71.4375*eccentricity_3 + 11.5*eccentricity;
    double common_term_73 = 5279.929886153825*eccentricity_14 + 1811.311626984127*eccentricity_12 + 1850.3712962962963*eccentricity_10 - 998.73888888888889*eccentricity_8 + 1195.375*eccentricity_6 - 415.33333333333333*eccentricity_4 + 74.0*eccentricity_2;
    double common_term_74 = 1684.3778098848886*eccentricity_13 + 6538.4422776358468*eccentricity_11 - 5407.3172241210938*eccentricity_9 + 5068.97138671875*eccentricity_7 - 1908.99609375*eccentricity_5 + 351.1875*eccentricity_3;
    double common_term_75 = -4689.4129891998266*eccentricity_14 + 24165.240187872024*eccentricity_12 - 23403.038917824074*eccentricity_10 + 19177.660590277778*eccentricity_8 - 7436.5625*eccentricity_6 + 1369.9583333333333*eccentricity_4;
    double common_term_76 = 87843.350422844994*eccentricity_13 - 88156.786500256655*eccentricity_11 + 65947.025514439174*eccentricity_9 - 25597.524945746528*eccentricity_7 + 4649.3778645833333*eccentricity_5;
    double common_term_77 = 304441.09340401785*eccentricity_14 - 300270.63828125*eccentricity_12 + 209230.68113839286*eccentricity_10 - 80002.971428571429*eccentricity_8 + 14206.575*eccentricity_6;
    double common_term_78 = -945071.06842828232*eccentricity_13 + 619928.1896216586*eccentricity_11 - 231366.43252229236*eccentricity_9 + 39975.869204179067*eccentricity_7;
    double common_term_79 = -2787795.5507886884*eccentricity_14 + 1732188.5863252315*eccentricity_12 - 627670.41844135803*eccentricity_10 + 105238.97986111111*eccentricity_8;
    double common_term_80 = 4600908.8919337972*eccentricity_13 - 1613867.4114398956*eccentricity_11 + 262199.80908857073*eccentricity_9;
    double common_term_81 = 11692765.550166911*eccentricity_14 - 3964207.7674797078*eccentricity_12 + 623657.13672674163*eccentricity_10;
    double common_term_82 = -9360904.4284804557*eccentricity_13 + 1425768.8055109258*eccentricity_11;
    double common_term_83 = -21356956.851246468*eccentricity_14 + 3149673.7737829748*eccentricity_12;
    double common_term_84 = 6752590.3413372089*eccentricity_13;
    double common_term_85 = 14099386.647976717*eccentricity_14;
    double common_term_86 = 597.8781373447734*eccentricity_14;
    double common_term_87 = 401.30327557254624*eccentricity_13;
    double common_term_88 = 1820.8968782889725*eccentricity_14 + 268.0270107093588*eccentricity_12;
    double common_term_89 = 1275.1064597254856*eccentricity_13 + 178.01548490179642*eccentricity_11;
    double common_term_90 = 3907.4299074168019*eccentricity_14 + 886.67794693587662*eccentricity_12 + 117.48426897321429*eccentricity_10;
    double common_term_91 = 2802.9247542444055*eccentricity_13 + 611.94196460528769*eccentricity_11 + 76.972772695820588*eccentricity_9;
    double common_term_92 = 7079.0373773140841*eccentricity_14 + 1995.7419670758929*eccentricity_12 + 418.86314966380071*eccentricity_10 + 50.005716765873016*eccentricity_8;
    double common_term_93 = 5163.7462441757747*eccentricity_13 + 1409.5655752454485*eccentricity_11 + 284.08580278669085*eccentricity_9 + 32.164327566964286*eccentricity_7;
    double common_term_94 = 11564.737798184547*eccentricity_14 + 3737.825384998714*eccentricity_12 + 986.71463603670635*eccentricity_10 + 190.68382936507936*eccentricity_8 + 20.443055555555556*eccentricity_6;
    double common_term_95 = 8542.0413697642904*eccentricity_13 + 2682.9199675282473*eccentricity_11 + 683.84254847935268*eccentricity_9 + 126.46257595486111*eccentricity_7 + 12.804947916666667*eccentricity_5;
    double common_term_96 = std::pow(1.0 - eccentricity_2, -9.5)*(0.984375*eccentricity_8 + 7.875*eccentricity_6 + 7.875*eccentricity_4);
    double common_term_97 = 13124.511923522621*eccentricity_13 + 4547.1542784473883*eccentricity_11 + 1342.1533343279803*eccentricity_9 + 316.81741536458333*eccentricity_7 + 53.139322916666667*eccentricity_5 + 4.7291666666666667*eccentricity_3;
    double common_term_98 = 25407.207472106826*eccentricity_14 + 9709.4953354414683*eccentricity_12 + 3270.7477068865741*eccentricity_10 + 932.75069444444444*eccentricity_8 + 210.84375*eccentricity_6 + 33.416666666666667*eccentricity_4 + 2.75*eccentricity_2;
    double common_term_99 = 19097.12023960386*eccentricity_13 + 7117.9510542297363*eccentricity_11 + 2326.2882019042969*eccentricity_9 + 638.93115234375*eccentricity_7 + 137.6015625*eccentricity_5 + 20.4375*eccentricity_3 + 1.5*eccentricity;
    double common_term_100 = 35225.259770100604*eccentricity_14 + 14233.868095823688*eccentricity_12 + 5164.4587934027778*eccentricity_10 + 1632.9787326388889*eccentricity_8 + 429.99305555555556*eccentricity_6 + 87.875*eccentricity_4 + 11.5*eccentricity_2 + 1.0;
    double common_term_101 = 26642.034517889511*eccentricity_13 + 10508.832178045202*eccentricity_11 + 3702.8700052897135*eccentricity_9 + 1127.3065863715278*eccentricity_7 + 286.43489583333333*eccentricity_5 + 48.9375*eccentricity_3 + 9.5*eccentricity;
    double common_term_102 = 47273.520184151786*eccentricity_14 + 19974.563219866071*eccentricity_12 + 7677.5900390625*eccentricity_10 + 2607.13125*eccentricity_8 + 789.375*eccentricity_6 + 147.75*eccentricity_4 + 51.75*eccentricity_2;
    double common_term_103 = 35929.081201025822*eccentricity_13 + 14851.215758062292*eccentricity_11 + 5465.05249520761*eccentricity_9 + 1960.3046549479167*eccentricity_7 + 333.49348958333333*eccentricity_5 + 211.89583333333333*eccentricity_3;
    double common_term_104 = 61724.547504213973*eccentricity_14 + 27210.794546750992*eccentricity_12 + 10508.673544973545*eccentricity_10 + 4597.7611111111111*eccentricity_8 + 506.675*eccentricity_6 + 724.08333333333333*eccentricity_4;
    double common_term_105 = 47884.189076154573*eccentricity_13 + 18473.2914163317*eccentricity_11 + 10579.168385532924*eccentricity_9 + 65.10849609375*eccentricity_7 + 2179.52578125*eccentricity_5;
    double common_term_106 = 82251.798673904482*eccentricity_14 + 28825.261025224133*eccentricity_12 + 24545.98294890873*eccentricity_10 - 3120.3387896825397*eccentricity_8 + 5967.9826388888889*eccentricity_6;
    double common_term_107 = 35951.121408174835*eccentricity_13 + 58159.198572374105*eccentricity_11 - 15211.783884248279*eccentricity_9 + 15180.904791356647*eccentricity_7;
    double common_term_108 = 18466.207672737418*eccentricity_14 + 140388.0997265625*eccentricity_12 - 52121.822209821429*eccentricity_10 + 36398.464955357143*eccentricity_8;
    double common_term_109 = 341351.77715522546*eccentricity_13 - 151901.36367563628*eccentricity_11 + 83129.35033802389*eccentricity_9;
    double common_term_110 = 825479.23414268188*eccentricity_14 - 400520.69452809344*eccentricity_12 + 182279.58837191358*eccentricity_10;
    double common_term_111 = -983804.59866731663*eccentricity_13 + 386077.27742586507*eccentricity_11;
    double common_term_112 = -2289259.6301116501*eccentricity_14 + 793676.63812026724*eccentricity_12;
    double common_term_113 = 1589705.6410422374*eccentricity_13;
    double common_term_114 = 3112121.1563314852*eccentricity_14;
    double common_term_115 = 8770.7349815502578*eccentricity_14;
    double common_term_116 = 5504.6219770765671*eccentricity_13;
    double common_term_117 = 17901.450702323614*eccentricity_14 + 3421.4312731157892*eccentricity_12;
    double common_term_118 = 11887.858538566255*eccentricity_13 + 2103.2376849850741*eccentricity_11;
    double common_term_119 = 29952.798654308253*eccentricity_14 + 7783.9591872594998*eccentricity_12 + 1276.5633008156966*eccentricity_10;
    double common_term_120 = 20304.147889377378*eccentricity_13 + 5019.5518635455473*eccentricity_11 + 763.40796896852093*eccentricity_9;
    double common_term_121 = 44359.788080166903*eccentricity_14 + 13566.126449497768*eccentricity_12 + 3182.6087611607143*eccentricity_10 + 448.60078125*eccentricity_8;
    double common_term_122 = 30434.914156567132*eccentricity_13 + 8918.7079399182893*eccentricity_11 + 1979.6183682396298*eccentricity_9 + 258.11415705605159*eccentricity_7;
    double common_term_123 = 60807.581296233144*eccentricity_14 + 20563.361818185994*eccentricity_12 + 5756.4458023313492*eccentricity_10 + 1204.2793650793651*eccentricity_8 + 144.71805555555556*eccentricity_6;
    double common_term_124 = 42006.00026381697*eccentricity_13 + 13653.1116166796*eccentricity_11 + 3636.8901907784598*eccentricity_9 + 713.45830078125*eccentricity_7 + 78.53203125*eccentricity_5;
    double common_term_125 = 78872.207621392862*eccentricity_14 + 28542.739133029514*eccentricity_12 + 8883.2707671957672*eccentricity_10 + 2240.1704861111111*eccentricity_8 + 409.125*eccentricity_6 + 40.833333333333333*eccentricity_4;
    double common_term_126 = 54657.662972711104*eccentricity_13 + 19027.201894096092*eccentricity_11 + 5642.6422087492766*eccentricity_9 + 1337.6403971354167*eccentricity_7 + 225.02473958333333*eccentricity_5 + 20.020833333333333*eccentricity_3;
    double common_term_127 = std::pow(1.0 - eccentricity_2, -9.5)*(1.96875*eccentricity_8 + 19.6875*eccentricity_6 + 31.5*eccentricity_4 + 9.0*eccentricity_2);
    double common_term_128 = 67951.780440018097*eccentricity_13 + 24786.066802034731*eccentricity_11 + 7862.1723083496094*eccentricity_9 + 2069.1467556423611*eccentricity_7 + 418.09114583333333*eccentricity_5 + 56.0625*eccentricity_3 + 3.5*eccentricity;
    double common_term_129 = 117705.31147752381*eccentricity_14 + 46166.126637128665*eccentricity_12 + 16094.778923611111*eccentricity_10 + 4815.6284722222222*eccentricity_8 + 1171.1284722222222*eccentricity_6 + 211.0625*eccentricity_4 + 23.5*eccentricity_2 + 1.0;
    double common_term_130 = 81375.283556417738*eccentricity_13 + 30618.552248382568*eccentricity_11 + 10120.135345458984*eccentricity_9 + 2818.54541015625*eccentricity_7 + 618.2578125*eccentricity_5 + 94.3125*eccentricity_3 + 7.5*eccentricity;
    double common_term_131 = 137185.77507562762*eccentricity_14 + 55005.507094494048*eccentricity_12 + 19708.091304976852*eccentricity_10 + 6101.6444444444444*eccentricity_8 + 1547.75*eccentricity_6 + 292.16666666666667*eccentricity_4 + 33.5*eccentricity_2;
    double common_term_132 = 94343.864676982007*eccentricity_13 + 36160.464863953767*eccentricity_11 + 12203.147214536314*eccentricity_9 + 3471.0574544270833*eccentricity_7 + 768.07682291666667*eccentricity_5 + 115.60416666666667*eccentricity_3;
    double common_term_133 = 155714.95866420201*eccentricity_14 + 63224.796323939732*eccentricity_12 + 22936.734877232143*eccentricity_10 + 7169.06953125*eccentricity_8 + 1795.36875*eccentricity_6 + 339.9375*eccentricity_4;
    double common_term_134 = 106207.14641562381*eccentricity_13 + 40991.901643187911*eccentricity_11 + 13885.421431477865*eccentricity_9 + 3829.5984483506944*eccentricity_7 + 895.01744791666667*eccentricity_5;
    double common_term_135 = 172461.03920998907*eccentricity_14 + 70237.92071851117*eccentricity_12 + 25550.523834325397*eccentricity_10 + 7566.1055555555556*eccentricity_8 + 2172.0972222222222*eccentricity_6;
    double common_term_136 = 116073.89904596737*eccentricity_13 + 45131.747396033151*eccentricity_11 + 13952.05462777274*eccentricity_9 + 4950.4390764508929*eccentricity_7;
    double common_term_137 = 185762.48939678646*eccentricity_14 + 77236.197046027612*eccentricity_12 + 24050.714995315256*eccentricity_10 + 10731.751364087302*eccentricity_8;
    double common_term_138 = 129249.85086463271*eccentricity_13 + 38539.346853018328*eccentricity_11 + 22332.510335668601*eccentricity_9;
    double common_term_139 = 213629.86151481331*eccentricity_14 + 56401.696065340909*eccentricity_12 + 44916.095558035714*eccentricity_10;
    double common_term_140 = 71970.9917549521*eccentricity_13 + 87766.10919416494*eccentricity_11;
    double common_term_141 = 68751.752767771307*eccentricity_14 + 167295.6226284693*eccentricity_12;
    double common_term_142 = 312097.97744729076*eccentricity_13;
    double common_term_143 = 571341.39167059263*eccentricity_14;
    double common_term_144 = 82829.353837736057*eccentricity_14;
    double common_term_145 = 48483.878171030515*eccentricity_13;
    double common_term_146 = 87060.351743539664*eccentricity_14 + 27972.480186941964*eccentricity_12;
    double common_term_147 = 56324.621720736673*eccentricity_13 + 15874.286491991612*eccentricity_11;
    double common_term_148 = 112955.73072896573*eccentricity_14 + 35477.132213128307*eccentricity_12 + 8838.6522851906966*eccentricity_10;
    double common_term_149 = 72943.453404656819*eccentricity_13 + 21744.861643011911*eccentricity_11 + 4813.0695709228516*eccentricity_9;
    double common_term_150 = 132926.16607212899*eccentricity_14 + 46094.235534267526*eccentricity_12 + 12948.758050870811*eccentricity_10 + 2552.854191468254*eccentricity_8;
    double common_term_151 = 86225.641583319809*eccentricity_13 + 28415.779286085529*eccentricity_11 + 7470.4682444254557*eccentricity_9 + 1311.7809136284722*eccentricity_7;
    double common_term_152 = 149262.66018554688*eccentricity_14 + 54632.448995535714*eccentricity_12 + 17025.143805803571*eccentricity_10 + 4157.7741071428571*eccentricity_8 + 648.278125*eccentricity_6;
    double common_term_153 = 96674.068711462405*eccentricity_13 + 33690.52670468729*eccentricity_11 + 9865.7245890299479*eccentricity_9 + 2218.4079318576389*eccentricity_7 + 304.97786458333333*eccentricity_5;
    double common_term_154 = 161311.29156125418*eccentricity_14 + 61022.377556113591*eccentricity_12 + 20127.049388227513*eccentricity_10 + 5493.2090277777778*eccentricity_8 + 1124.26875*eccentricity_6 + 134.52083333333333*eccentricity_4;
    double common_term_155 = 103860.76791643483*eccentricity_13 + 37371.455359431676*eccentricity_11 + 11575.071057128906*eccentricity_9 + 2911.84072265625*eccentricity_7 + 533.62890625*eccentricity_5 + 54.3125*eccentricity_3;
    double common_term_156 = 168692.96589035976*eccentricity_14 + 64971.088209325397*eccentricity_12 + 22071.572251157407*eccentricity_10 + 6351.5298611111111*eccentricity_8 + 1449.421875*eccentricity_6 + 231.91666666666667*eccentricity_4 + 19.25*eccentricity_2;
    double common_term_157 = 107520.56849235895*eccentricity_13 + 39261.090514791983*eccentricity_11 + 12464.666538492839*eccentricity_9 + 3281.9704318576389*eccentricity_7 + 662.89322916666667*eccentricity_5 + 88.6875*eccentricity_3 + 5.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(29);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (10, 0).
    // q = -14
    result_by_lpq.set(c_Key3(10, 0, -14), common_term_0);
    result_by_q.set(c_Key1(-14), common_term_0);
    // q = -13
    result_by_lpq.set(c_Key3(10, 0, -13), common_term_1);
    result_by_q.set(c_Key1(-13), common_term_1);
    // q = -12
    result_by_lpq.set(c_Key3(10, 0, -12), common_term_2);
    result_by_q.set(c_Key1(-12), common_term_2);
    // q = -11
    result_by_lpq.set(c_Key3(10, 0, -11), common_term_3);
    result_by_q.set(c_Key1(-11), common_term_3);
    // q = -9
    result_by_lpq.set(c_Key3(10, 0, -9), common_term_4);
    result_by_q.set(c_Key1(-9), common_term_4);
    // q = -8
    result_by_lpq.set(c_Key3(10, 0, -8), common_term_5);
    result_by_q.set(c_Key1(-8), common_term_5);
    // q = -7
    result_by_lpq.set(c_Key3(10, 0, -7), common_term_6);
    result_by_q.set(c_Key1(-7), common_term_6);
    // q = -6
    result_by_lpq.set(c_Key3(10, 0, -6), common_term_7);
    result_by_q.set(c_Key1(-6), common_term_7);
    // q = -5
    result_by_lpq.set(c_Key3(10, 0, -5), common_term_8);
    result_by_q.set(c_Key1(-5), common_term_8);
    // q = -4
    result_by_lpq.set(c_Key3(10, 0, -4), common_term_9);
    result_by_q.set(c_Key1(-4), common_term_9);
    // q = -3
    result_by_lpq.set(c_Key3(10, 0, -3), common_term_10);
    result_by_q.set(c_Key1(-3), common_term_10);
    // q = -2
    result_by_lpq.set(c_Key3(10, 0, -2), common_term_11);
    result_by_q.set(c_Key1(-2), common_term_11);
    // q = -1
    result_by_lpq.set(c_Key3(10, 0, -1), common_term_12);
    result_by_q.set(c_Key1(-1), common_term_12);
    // q = 0
    result_by_lpq.set(c_Key3(10, 0, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(10, 0, 1), common_term_14);
    result_by_q.set(c_Key1(1), common_term_14);
    // q = 2
    result_by_lpq.set(c_Key3(10, 0, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // q = 3
    result_by_lpq.set(c_Key3(10, 0, 3), common_term_16);
    result_by_q.set(c_Key1(3), common_term_16);
    // q = 4
    result_by_lpq.set(c_Key3(10, 0, 4), common_term_17);
    result_by_q.set(c_Key1(4), common_term_17);
    // q = 5
    result_by_lpq.set(c_Key3(10, 0, 5), common_term_18);
    result_by_q.set(c_Key1(5), common_term_18);
    // q = 6
    result_by_lpq.set(c_Key3(10, 0, 6), common_term_19);
    result_by_q.set(c_Key1(6), common_term_19);
    // q = 7
    result_by_lpq.set(c_Key3(10, 0, 7), common_term_20);
    result_by_q.set(c_Key1(7), common_term_20);
    // q = 8
    result_by_lpq.set(c_Key3(10, 0, 8), common_term_21);
    result_by_q.set(c_Key1(8), common_term_21);
    // q = 9
    result_by_lpq.set(c_Key3(10, 0, 9), common_term_22);
    result_by_q.set(c_Key1(9), common_term_22);
    // q = 10
    result_by_lpq.set(c_Key3(10, 0, 10), common_term_23);
    result_by_q.set(c_Key1(10), common_term_23);
    // q = 11
    result_by_lpq.set(c_Key3(10, 0, 11), common_term_24);
    result_by_q.set(c_Key1(11), common_term_24);
    // q = 12
    result_by_lpq.set(c_Key3(10, 0, 12), common_term_25);
    result_by_q.set(c_Key1(12), common_term_25);
    // q = 13
    result_by_lpq.set(c_Key3(10, 0, 13), common_term_26);
    result_by_q.set(c_Key1(13), common_term_26);
    // q = 14
    result_by_lpq.set(c_Key3(10, 0, 14), common_term_27);
    result_by_q.set(c_Key1(14), common_term_27);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 0), result_by_q);
    result_by_q.clear();

    // l , p = (10, 1).
    // q = -14
    result_by_lpq.set(c_Key3(10, 1, -14), common_term_28);
    result_by_q.set(c_Key1(-14), common_term_28);
    // q = -13
    result_by_lpq.set(c_Key3(10, 1, -13), common_term_29);
    result_by_q.set(c_Key1(-13), common_term_29);
    // q = -12
    result_by_lpq.set(c_Key3(10, 1, -12), common_term_30);
    result_by_q.set(c_Key1(-12), common_term_30);
    // q = -11
    result_by_lpq.set(c_Key3(10, 1, -11), common_term_31);
    result_by_q.set(c_Key1(-11), common_term_31);
    // q = -10
    result_by_lpq.set(c_Key3(10, 1, -10), common_term_32);
    result_by_q.set(c_Key1(-10), common_term_32);
    // q = -9
    result_by_lpq.set(c_Key3(10, 1, -9), common_term_33);
    result_by_q.set(c_Key1(-9), common_term_33);
    // q = -8
    result_by_lpq.set(c_Key3(10, 1, -8), common_term_34);
    result_by_q.set(c_Key1(-8), common_term_34);
    // q = -7
    result_by_lpq.set(c_Key3(10, 1, -7), common_term_35);
    result_by_q.set(c_Key1(-7), common_term_35);
    // q = -6
    result_by_lpq.set(c_Key3(10, 1, -6), common_term_36);
    result_by_q.set(c_Key1(-6), common_term_36);
    // q = -5
    result_by_lpq.set(c_Key3(10, 1, -5), common_term_37);
    result_by_q.set(c_Key1(-5), common_term_37);
    // q = -4
    result_by_lpq.set(c_Key3(10, 1, -4), common_term_38);
    result_by_q.set(c_Key1(-4), common_term_38);
    // q = -3
    result_by_lpq.set(c_Key3(10, 1, -3), common_term_39);
    result_by_q.set(c_Key1(-3), common_term_39);
    // q = -2
    result_by_lpq.set(c_Key3(10, 1, -2), common_term_40);
    result_by_q.set(c_Key1(-2), common_term_40);
    // q = -1
    result_by_lpq.set(c_Key3(10, 1, -1), common_term_41);
    result_by_q.set(c_Key1(-1), common_term_41);
    // q = 0
    result_by_lpq.set(c_Key3(10, 1, 0), common_term_42);
    result_by_q.set(c_Key1(0), common_term_42);
    // q = 1
    result_by_lpq.set(c_Key3(10, 1, 1), common_term_43);
    result_by_q.set(c_Key1(1), common_term_43);
    // q = 2
    result_by_lpq.set(c_Key3(10, 1, 2), common_term_44);
    result_by_q.set(c_Key1(2), common_term_44);
    // q = 3
    result_by_lpq.set(c_Key3(10, 1, 3), common_term_45);
    result_by_q.set(c_Key1(3), common_term_45);
    // q = 4
    result_by_lpq.set(c_Key3(10, 1, 4), common_term_46);
    result_by_q.set(c_Key1(4), common_term_46);
    // q = 5
    result_by_lpq.set(c_Key3(10, 1, 5), common_term_47);
    result_by_q.set(c_Key1(5), common_term_47);
    // q = 6
    result_by_lpq.set(c_Key3(10, 1, 6), common_term_48);
    result_by_q.set(c_Key1(6), common_term_48);
    // q = 7
    result_by_lpq.set(c_Key3(10, 1, 7), common_term_49);
    result_by_q.set(c_Key1(7), common_term_49);
    // q = 8
    result_by_lpq.set(c_Key3(10, 1, 8), common_term_50);
    result_by_q.set(c_Key1(8), common_term_50);
    // q = 9
    result_by_lpq.set(c_Key3(10, 1, 9), common_term_51);
    result_by_q.set(c_Key1(9), common_term_51);
    // q = 10
    result_by_lpq.set(c_Key3(10, 1, 10), common_term_52);
    result_by_q.set(c_Key1(10), common_term_52);
    // q = 11
    result_by_lpq.set(c_Key3(10, 1, 11), common_term_53);
    result_by_q.set(c_Key1(11), common_term_53);
    // q = 12
    result_by_lpq.set(c_Key3(10, 1, 12), common_term_54);
    result_by_q.set(c_Key1(12), common_term_54);
    // q = 13
    result_by_lpq.set(c_Key3(10, 1, 13), common_term_55);
    result_by_q.set(c_Key1(13), common_term_55);
    // q = 14
    result_by_lpq.set(c_Key3(10, 1, 14), common_term_56);
    result_by_q.set(c_Key1(14), common_term_56);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 1), result_by_q);
    result_by_q.clear();

    // l , p = (10, 2).
    // q = -14
    result_by_lpq.set(c_Key3(10, 2, -14), common_term_57);
    result_by_q.set(c_Key1(-14), common_term_57);
    // q = -13
    result_by_lpq.set(c_Key3(10, 2, -13), common_term_58);
    result_by_q.set(c_Key1(-13), common_term_58);
    // q = -12
    result_by_lpq.set(c_Key3(10, 2, -12), common_term_59);
    result_by_q.set(c_Key1(-12), common_term_59);
    // q = -11
    result_by_lpq.set(c_Key3(10, 2, -11), common_term_60);
    result_by_q.set(c_Key1(-11), common_term_60);
    // q = -10
    result_by_lpq.set(c_Key3(10, 2, -10), common_term_61);
    result_by_q.set(c_Key1(-10), common_term_61);
    // q = -9
    result_by_lpq.set(c_Key3(10, 2, -9), common_term_62);
    result_by_q.set(c_Key1(-9), common_term_62);
    // q = -8
    result_by_lpq.set(c_Key3(10, 2, -8), common_term_63);
    result_by_q.set(c_Key1(-8), common_term_63);
    // q = -7
    result_by_lpq.set(c_Key3(10, 2, -7), common_term_64);
    result_by_q.set(c_Key1(-7), common_term_64);
    // q = -6
    result_by_lpq.set(c_Key3(10, 2, -6), common_term_65);
    result_by_q.set(c_Key1(-6), common_term_65);
    // q = -5
    result_by_lpq.set(c_Key3(10, 2, -5), common_term_66);
    result_by_q.set(c_Key1(-5), common_term_66);
    // q = -4
    result_by_lpq.set(c_Key3(10, 2, -4), common_term_67);
    result_by_q.set(c_Key1(-4), common_term_67);
    // q = -3
    result_by_lpq.set(c_Key3(10, 2, -3), common_term_68);
    result_by_q.set(c_Key1(-3), common_term_68);
    // q = -2
    result_by_lpq.set(c_Key3(10, 2, -2), common_term_69);
    result_by_q.set(c_Key1(-2), common_term_69);
    // q = -1
    result_by_lpq.set(c_Key3(10, 2, -1), common_term_70);
    result_by_q.set(c_Key1(-1), common_term_70);
    // q = 0
    result_by_lpq.set(c_Key3(10, 2, 0), common_term_71);
    result_by_q.set(c_Key1(0), common_term_71);
    // q = 1
    result_by_lpq.set(c_Key3(10, 2, 1), common_term_72);
    result_by_q.set(c_Key1(1), common_term_72);
    // q = 2
    result_by_lpq.set(c_Key3(10, 2, 2), common_term_73);
    result_by_q.set(c_Key1(2), common_term_73);
    // q = 3
    result_by_lpq.set(c_Key3(10, 2, 3), common_term_74);
    result_by_q.set(c_Key1(3), common_term_74);
    // q = 4
    result_by_lpq.set(c_Key3(10, 2, 4), common_term_75);
    result_by_q.set(c_Key1(4), common_term_75);
    // q = 5
    result_by_lpq.set(c_Key3(10, 2, 5), common_term_76);
    result_by_q.set(c_Key1(5), common_term_76);
    // q = 6
    result_by_lpq.set(c_Key3(10, 2, 6), common_term_77);
    result_by_q.set(c_Key1(6), common_term_77);
    // q = 7
    result_by_lpq.set(c_Key3(10, 2, 7), common_term_78);
    result_by_q.set(c_Key1(7), common_term_78);
    // q = 8
    result_by_lpq.set(c_Key3(10, 2, 8), common_term_79);
    result_by_q.set(c_Key1(8), common_term_79);
    // q = 9
    result_by_lpq.set(c_Key3(10, 2, 9), common_term_80);
    result_by_q.set(c_Key1(9), common_term_80);
    // q = 10
    result_by_lpq.set(c_Key3(10, 2, 10), common_term_81);
    result_by_q.set(c_Key1(10), common_term_81);
    // q = 11
    result_by_lpq.set(c_Key3(10, 2, 11), common_term_82);
    result_by_q.set(c_Key1(11), common_term_82);
    // q = 12
    result_by_lpq.set(c_Key3(10, 2, 12), common_term_83);
    result_by_q.set(c_Key1(12), common_term_83);
    // q = 13
    result_by_lpq.set(c_Key3(10, 2, 13), common_term_84);
    result_by_q.set(c_Key1(13), common_term_84);
    // q = 14
    result_by_lpq.set(c_Key3(10, 2, 14), common_term_85);
    result_by_q.set(c_Key1(14), common_term_85);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 2), result_by_q);
    result_by_q.clear();

    // l , p = (10, 3).
    // q = -14
    result_by_lpq.set(c_Key3(10, 3, -14), common_term_86);
    result_by_q.set(c_Key1(-14), common_term_86);
    // q = -13
    result_by_lpq.set(c_Key3(10, 3, -13), common_term_87);
    result_by_q.set(c_Key1(-13), common_term_87);
    // q = -12
    result_by_lpq.set(c_Key3(10, 3, -12), common_term_88);
    result_by_q.set(c_Key1(-12), common_term_88);
    // q = -11
    result_by_lpq.set(c_Key3(10, 3, -11), common_term_89);
    result_by_q.set(c_Key1(-11), common_term_89);
    // q = -10
    result_by_lpq.set(c_Key3(10, 3, -10), common_term_90);
    result_by_q.set(c_Key1(-10), common_term_90);
    // q = -9
    result_by_lpq.set(c_Key3(10, 3, -9), common_term_91);
    result_by_q.set(c_Key1(-9), common_term_91);
    // q = -8
    result_by_lpq.set(c_Key3(10, 3, -8), common_term_92);
    result_by_q.set(c_Key1(-8), common_term_92);
    // q = -7
    result_by_lpq.set(c_Key3(10, 3, -7), common_term_93);
    result_by_q.set(c_Key1(-7), common_term_93);
    // q = -6
    result_by_lpq.set(c_Key3(10, 3, -6), common_term_94);
    result_by_q.set(c_Key1(-6), common_term_94);
    // q = -5
    result_by_lpq.set(c_Key3(10, 3, -5), common_term_95);
    result_by_q.set(c_Key1(-5), common_term_95);
    // q = -4
    result_by_lpq.set(c_Key3(10, 3, -4), common_term_96);
    result_by_q.set(c_Key1(-4), common_term_96);
    // q = -3
    result_by_lpq.set(c_Key3(10, 3, -3), common_term_97);
    result_by_q.set(c_Key1(-3), common_term_97);
    // q = -2
    result_by_lpq.set(c_Key3(10, 3, -2), common_term_98);
    result_by_q.set(c_Key1(-2), common_term_98);
    // q = -1
    result_by_lpq.set(c_Key3(10, 3, -1), common_term_99);
    result_by_q.set(c_Key1(-1), common_term_99);
    // q = 0
    result_by_lpq.set(c_Key3(10, 3, 0), common_term_100);
    result_by_q.set(c_Key1(0), common_term_100);
    // q = 1
    result_by_lpq.set(c_Key3(10, 3, 1), common_term_101);
    result_by_q.set(c_Key1(1), common_term_101);
    // q = 2
    result_by_lpq.set(c_Key3(10, 3, 2), common_term_102);
    result_by_q.set(c_Key1(2), common_term_102);
    // q = 3
    result_by_lpq.set(c_Key3(10, 3, 3), common_term_103);
    result_by_q.set(c_Key1(3), common_term_103);
    // q = 4
    result_by_lpq.set(c_Key3(10, 3, 4), common_term_104);
    result_by_q.set(c_Key1(4), common_term_104);
    // q = 5
    result_by_lpq.set(c_Key3(10, 3, 5), common_term_105);
    result_by_q.set(c_Key1(5), common_term_105);
    // q = 6
    result_by_lpq.set(c_Key3(10, 3, 6), common_term_106);
    result_by_q.set(c_Key1(6), common_term_106);
    // q = 7
    result_by_lpq.set(c_Key3(10, 3, 7), common_term_107);
    result_by_q.set(c_Key1(7), common_term_107);
    // q = 8
    result_by_lpq.set(c_Key3(10, 3, 8), common_term_108);
    result_by_q.set(c_Key1(8), common_term_108);
    // q = 9
    result_by_lpq.set(c_Key3(10, 3, 9), common_term_109);
    result_by_q.set(c_Key1(9), common_term_109);
    // q = 10
    result_by_lpq.set(c_Key3(10, 3, 10), common_term_110);
    result_by_q.set(c_Key1(10), common_term_110);
    // q = 11
    result_by_lpq.set(c_Key3(10, 3, 11), common_term_111);
    result_by_q.set(c_Key1(11), common_term_111);
    // q = 12
    result_by_lpq.set(c_Key3(10, 3, 12), common_term_112);
    result_by_q.set(c_Key1(12), common_term_112);
    // q = 13
    result_by_lpq.set(c_Key3(10, 3, 13), common_term_113);
    result_by_q.set(c_Key1(13), common_term_113);
    // q = 14
    result_by_lpq.set(c_Key3(10, 3, 14), common_term_114);
    result_by_q.set(c_Key1(14), common_term_114);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 3), result_by_q);
    result_by_q.clear();

    // l , p = (10, 4).
    // q = -14
    result_by_lpq.set(c_Key3(10, 4, -14), common_term_115);
    result_by_q.set(c_Key1(-14), common_term_115);
    // q = -13
    result_by_lpq.set(c_Key3(10, 4, -13), common_term_116);
    result_by_q.set(c_Key1(-13), common_term_116);
    // q = -12
    result_by_lpq.set(c_Key3(10, 4, -12), common_term_117);
    result_by_q.set(c_Key1(-12), common_term_117);
    // q = -11
    result_by_lpq.set(c_Key3(10, 4, -11), common_term_118);
    result_by_q.set(c_Key1(-11), common_term_118);
    // q = -10
    result_by_lpq.set(c_Key3(10, 4, -10), common_term_119);
    result_by_q.set(c_Key1(-10), common_term_119);
    // q = -9
    result_by_lpq.set(c_Key3(10, 4, -9), common_term_120);
    result_by_q.set(c_Key1(-9), common_term_120);
    // q = -8
    result_by_lpq.set(c_Key3(10, 4, -8), common_term_121);
    result_by_q.set(c_Key1(-8), common_term_121);
    // q = -7
    result_by_lpq.set(c_Key3(10, 4, -7), common_term_122);
    result_by_q.set(c_Key1(-7), common_term_122);
    // q = -6
    result_by_lpq.set(c_Key3(10, 4, -6), common_term_123);
    result_by_q.set(c_Key1(-6), common_term_123);
    // q = -5
    result_by_lpq.set(c_Key3(10, 4, -5), common_term_124);
    result_by_q.set(c_Key1(-5), common_term_124);
    // q = -4
    result_by_lpq.set(c_Key3(10, 4, -4), common_term_125);
    result_by_q.set(c_Key1(-4), common_term_125);
    // q = -3
    result_by_lpq.set(c_Key3(10, 4, -3), common_term_126);
    result_by_q.set(c_Key1(-3), common_term_126);
    // q = -2
    result_by_lpq.set(c_Key3(10, 4, -2), common_term_127);
    result_by_q.set(c_Key1(-2), common_term_127);
    // q = -1
    result_by_lpq.set(c_Key3(10, 4, -1), common_term_128);
    result_by_q.set(c_Key1(-1), common_term_128);
    // q = 0
    result_by_lpq.set(c_Key3(10, 4, 0), common_term_129);
    result_by_q.set(c_Key1(0), common_term_129);
    // q = 1
    result_by_lpq.set(c_Key3(10, 4, 1), common_term_130);
    result_by_q.set(c_Key1(1), common_term_130);
    // q = 2
    result_by_lpq.set(c_Key3(10, 4, 2), common_term_131);
    result_by_q.set(c_Key1(2), common_term_131);
    // q = 3
    result_by_lpq.set(c_Key3(10, 4, 3), common_term_132);
    result_by_q.set(c_Key1(3), common_term_132);
    // q = 4
    result_by_lpq.set(c_Key3(10, 4, 4), common_term_133);
    result_by_q.set(c_Key1(4), common_term_133);
    // q = 5
    result_by_lpq.set(c_Key3(10, 4, 5), common_term_134);
    result_by_q.set(c_Key1(5), common_term_134);
    // q = 6
    result_by_lpq.set(c_Key3(10, 4, 6), common_term_135);
    result_by_q.set(c_Key1(6), common_term_135);
    // q = 7
    result_by_lpq.set(c_Key3(10, 4, 7), common_term_136);
    result_by_q.set(c_Key1(7), common_term_136);
    // q = 8
    result_by_lpq.set(c_Key3(10, 4, 8), common_term_137);
    result_by_q.set(c_Key1(8), common_term_137);
    // q = 9
    result_by_lpq.set(c_Key3(10, 4, 9), common_term_138);
    result_by_q.set(c_Key1(9), common_term_138);
    // q = 10
    result_by_lpq.set(c_Key3(10, 4, 10), common_term_139);
    result_by_q.set(c_Key1(10), common_term_139);
    // q = 11
    result_by_lpq.set(c_Key3(10, 4, 11), common_term_140);
    result_by_q.set(c_Key1(11), common_term_140);
    // q = 12
    result_by_lpq.set(c_Key3(10, 4, 12), common_term_141);
    result_by_q.set(c_Key1(12), common_term_141);
    // q = 13
    result_by_lpq.set(c_Key3(10, 4, 13), common_term_142);
    result_by_q.set(c_Key1(13), common_term_142);
    // q = 14
    result_by_lpq.set(c_Key3(10, 4, 14), common_term_143);
    result_by_q.set(c_Key1(14), common_term_143);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 4), result_by_q);
    result_by_q.clear();

    // l , p = (10, 5).
    // q = -14
    result_by_lpq.set(c_Key3(10, 5, -14), common_term_144);
    result_by_q.set(c_Key1(-14), common_term_144);
    // q = -13
    result_by_lpq.set(c_Key3(10, 5, -13), common_term_145);
    result_by_q.set(c_Key1(-13), common_term_145);
    // q = -12
    result_by_lpq.set(c_Key3(10, 5, -12), common_term_146);
    result_by_q.set(c_Key1(-12), common_term_146);
    // q = -11
    result_by_lpq.set(c_Key3(10, 5, -11), common_term_147);
    result_by_q.set(c_Key1(-11), common_term_147);
    // q = -10
    result_by_lpq.set(c_Key3(10, 5, -10), common_term_148);
    result_by_q.set(c_Key1(-10), common_term_148);
    // q = -9
    result_by_lpq.set(c_Key3(10, 5, -9), common_term_149);
    result_by_q.set(c_Key1(-9), common_term_149);
    // q = -8
    result_by_lpq.set(c_Key3(10, 5, -8), common_term_150);
    result_by_q.set(c_Key1(-8), common_term_150);
    // q = -7
    result_by_lpq.set(c_Key3(10, 5, -7), common_term_151);
    result_by_q.set(c_Key1(-7), common_term_151);
    // q = -6
    result_by_lpq.set(c_Key3(10, 5, -6), common_term_152);
    result_by_q.set(c_Key1(-6), common_term_152);
    // q = -5
    result_by_lpq.set(c_Key3(10, 5, -5), common_term_153);
    result_by_q.set(c_Key1(-5), common_term_153);
    // q = -4
    result_by_lpq.set(c_Key3(10, 5, -4), common_term_154);
    result_by_q.set(c_Key1(-4), common_term_154);
    // q = -3
    result_by_lpq.set(c_Key3(10, 5, -3), common_term_155);
    result_by_q.set(c_Key1(-3), common_term_155);
    // q = -2
    result_by_lpq.set(c_Key3(10, 5, -2), common_term_156);
    result_by_q.set(c_Key1(-2), common_term_156);
    // q = -1
    result_by_lpq.set(c_Key3(10, 5, -1), common_term_157);
    result_by_q.set(c_Key1(-1), common_term_157);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -9.5)*(2.4609375*eccentricity_8 + 26.25*eccentricity_6 + 47.25*eccentricity_4 + 18.0*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(10, 5, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(10, 5, 1), common_term_157);
    result_by_q.set(c_Key1(1), common_term_157);
    // q = 2
    result_by_lpq.set(c_Key3(10, 5, 2), common_term_156);
    result_by_q.set(c_Key1(2), common_term_156);
    // q = 3
    result_by_lpq.set(c_Key3(10, 5, 3), common_term_155);
    result_by_q.set(c_Key1(3), common_term_155);
    // q = 4
    result_by_lpq.set(c_Key3(10, 5, 4), common_term_154);
    result_by_q.set(c_Key1(4), common_term_154);
    // q = 5
    result_by_lpq.set(c_Key3(10, 5, 5), common_term_153);
    result_by_q.set(c_Key1(5), common_term_153);
    // q = 6
    result_by_lpq.set(c_Key3(10, 5, 6), common_term_152);
    result_by_q.set(c_Key1(6), common_term_152);
    // q = 7
    result_by_lpq.set(c_Key3(10, 5, 7), common_term_151);
    result_by_q.set(c_Key1(7), common_term_151);
    // q = 8
    result_by_lpq.set(c_Key3(10, 5, 8), common_term_150);
    result_by_q.set(c_Key1(8), common_term_150);
    // q = 9
    result_by_lpq.set(c_Key3(10, 5, 9), common_term_149);
    result_by_q.set(c_Key1(9), common_term_149);
    // q = 10
    result_by_lpq.set(c_Key3(10, 5, 10), common_term_148);
    result_by_q.set(c_Key1(10), common_term_148);
    // q = 11
    result_by_lpq.set(c_Key3(10, 5, 11), common_term_147);
    result_by_q.set(c_Key1(11), common_term_147);
    // q = 12
    result_by_lpq.set(c_Key3(10, 5, 12), common_term_146);
    result_by_q.set(c_Key1(12), common_term_146);
    // q = 13
    result_by_lpq.set(c_Key3(10, 5, 13), common_term_145);
    result_by_q.set(c_Key1(13), common_term_145);
    // q = 14
    result_by_lpq.set(c_Key3(10, 5, 14), common_term_144);
    result_by_q.set(c_Key1(14), common_term_144);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 5), result_by_q);
    result_by_q.clear();

    // l , p = (10, 6).
    // q = -14
    result_by_lpq.set(c_Key3(10, 6, -14), common_term_143);
    result_by_q.set(c_Key1(-14), common_term_143);
    // q = -13
    result_by_lpq.set(c_Key3(10, 6, -13), common_term_142);
    result_by_q.set(c_Key1(-13), common_term_142);
    // q = -12
    result_by_lpq.set(c_Key3(10, 6, -12), common_term_141);
    result_by_q.set(c_Key1(-12), common_term_141);
    // q = -11
    result_by_lpq.set(c_Key3(10, 6, -11), common_term_140);
    result_by_q.set(c_Key1(-11), common_term_140);
    // q = -10
    result_by_lpq.set(c_Key3(10, 6, -10), common_term_139);
    result_by_q.set(c_Key1(-10), common_term_139);
    // q = -9
    result_by_lpq.set(c_Key3(10, 6, -9), common_term_138);
    result_by_q.set(c_Key1(-9), common_term_138);
    // q = -8
    result_by_lpq.set(c_Key3(10, 6, -8), common_term_137);
    result_by_q.set(c_Key1(-8), common_term_137);
    // q = -7
    result_by_lpq.set(c_Key3(10, 6, -7), common_term_136);
    result_by_q.set(c_Key1(-7), common_term_136);
    // q = -6
    result_by_lpq.set(c_Key3(10, 6, -6), common_term_135);
    result_by_q.set(c_Key1(-6), common_term_135);
    // q = -5
    result_by_lpq.set(c_Key3(10, 6, -5), common_term_134);
    result_by_q.set(c_Key1(-5), common_term_134);
    // q = -4
    result_by_lpq.set(c_Key3(10, 6, -4), common_term_133);
    result_by_q.set(c_Key1(-4), common_term_133);
    // q = -3
    result_by_lpq.set(c_Key3(10, 6, -3), common_term_132);
    result_by_q.set(c_Key1(-3), common_term_132);
    // q = -2
    result_by_lpq.set(c_Key3(10, 6, -2), common_term_131);
    result_by_q.set(c_Key1(-2), common_term_131);
    // q = -1
    result_by_lpq.set(c_Key3(10, 6, -1), common_term_130);
    result_by_q.set(c_Key1(-1), common_term_130);
    // q = 0
    result_by_lpq.set(c_Key3(10, 6, 0), common_term_129);
    result_by_q.set(c_Key1(0), common_term_129);
    // q = 1
    result_by_lpq.set(c_Key3(10, 6, 1), common_term_128);
    result_by_q.set(c_Key1(1), common_term_128);
    // q = 2
    result_by_lpq.set(c_Key3(10, 6, 2), common_term_127);
    result_by_q.set(c_Key1(2), common_term_127);
    // q = 3
    result_by_lpq.set(c_Key3(10, 6, 3), common_term_126);
    result_by_q.set(c_Key1(3), common_term_126);
    // q = 4
    result_by_lpq.set(c_Key3(10, 6, 4), common_term_125);
    result_by_q.set(c_Key1(4), common_term_125);
    // q = 5
    result_by_lpq.set(c_Key3(10, 6, 5), common_term_124);
    result_by_q.set(c_Key1(5), common_term_124);
    // q = 6
    result_by_lpq.set(c_Key3(10, 6, 6), common_term_123);
    result_by_q.set(c_Key1(6), common_term_123);
    // q = 7
    result_by_lpq.set(c_Key3(10, 6, 7), common_term_122);
    result_by_q.set(c_Key1(7), common_term_122);
    // q = 8
    result_by_lpq.set(c_Key3(10, 6, 8), common_term_121);
    result_by_q.set(c_Key1(8), common_term_121);
    // q = 9
    result_by_lpq.set(c_Key3(10, 6, 9), common_term_120);
    result_by_q.set(c_Key1(9), common_term_120);
    // q = 10
    result_by_lpq.set(c_Key3(10, 6, 10), common_term_119);
    result_by_q.set(c_Key1(10), common_term_119);
    // q = 11
    result_by_lpq.set(c_Key3(10, 6, 11), common_term_118);
    result_by_q.set(c_Key1(11), common_term_118);
    // q = 12
    result_by_lpq.set(c_Key3(10, 6, 12), common_term_117);
    result_by_q.set(c_Key1(12), common_term_117);
    // q = 13
    result_by_lpq.set(c_Key3(10, 6, 13), common_term_116);
    result_by_q.set(c_Key1(13), common_term_116);
    // q = 14
    result_by_lpq.set(c_Key3(10, 6, 14), common_term_115);
    result_by_q.set(c_Key1(14), common_term_115);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 6), result_by_q);
    result_by_q.clear();

    // l , p = (10, 7).
    // q = -14
    result_by_lpq.set(c_Key3(10, 7, -14), common_term_114);
    result_by_q.set(c_Key1(-14), common_term_114);
    // q = -13
    result_by_lpq.set(c_Key3(10, 7, -13), common_term_113);
    result_by_q.set(c_Key1(-13), common_term_113);
    // q = -12
    result_by_lpq.set(c_Key3(10, 7, -12), common_term_112);
    result_by_q.set(c_Key1(-12), common_term_112);
    // q = -11
    result_by_lpq.set(c_Key3(10, 7, -11), common_term_111);
    result_by_q.set(c_Key1(-11), common_term_111);
    // q = -10
    result_by_lpq.set(c_Key3(10, 7, -10), common_term_110);
    result_by_q.set(c_Key1(-10), common_term_110);
    // q = -9
    result_by_lpq.set(c_Key3(10, 7, -9), common_term_109);
    result_by_q.set(c_Key1(-9), common_term_109);
    // q = -8
    result_by_lpq.set(c_Key3(10, 7, -8), common_term_108);
    result_by_q.set(c_Key1(-8), common_term_108);
    // q = -7
    result_by_lpq.set(c_Key3(10, 7, -7), common_term_107);
    result_by_q.set(c_Key1(-7), common_term_107);
    // q = -6
    result_by_lpq.set(c_Key3(10, 7, -6), common_term_106);
    result_by_q.set(c_Key1(-6), common_term_106);
    // q = -5
    result_by_lpq.set(c_Key3(10, 7, -5), common_term_105);
    result_by_q.set(c_Key1(-5), common_term_105);
    // q = -4
    result_by_lpq.set(c_Key3(10, 7, -4), common_term_104);
    result_by_q.set(c_Key1(-4), common_term_104);
    // q = -3
    result_by_lpq.set(c_Key3(10, 7, -3), common_term_103);
    result_by_q.set(c_Key1(-3), common_term_103);
    // q = -2
    result_by_lpq.set(c_Key3(10, 7, -2), common_term_102);
    result_by_q.set(c_Key1(-2), common_term_102);
    // q = -1
    result_by_lpq.set(c_Key3(10, 7, -1), common_term_101);
    result_by_q.set(c_Key1(-1), common_term_101);
    // q = 0
    result_by_lpq.set(c_Key3(10, 7, 0), common_term_100);
    result_by_q.set(c_Key1(0), common_term_100);
    // q = 1
    result_by_lpq.set(c_Key3(10, 7, 1), common_term_99);
    result_by_q.set(c_Key1(1), common_term_99);
    // q = 2
    result_by_lpq.set(c_Key3(10, 7, 2), common_term_98);
    result_by_q.set(c_Key1(2), common_term_98);
    // q = 3
    result_by_lpq.set(c_Key3(10, 7, 3), common_term_97);
    result_by_q.set(c_Key1(3), common_term_97);
    // q = 4
    result_by_lpq.set(c_Key3(10, 7, 4), common_term_96);
    result_by_q.set(c_Key1(4), common_term_96);
    // q = 5
    result_by_lpq.set(c_Key3(10, 7, 5), common_term_95);
    result_by_q.set(c_Key1(5), common_term_95);
    // q = 6
    result_by_lpq.set(c_Key3(10, 7, 6), common_term_94);
    result_by_q.set(c_Key1(6), common_term_94);
    // q = 7
    result_by_lpq.set(c_Key3(10, 7, 7), common_term_93);
    result_by_q.set(c_Key1(7), common_term_93);
    // q = 8
    result_by_lpq.set(c_Key3(10, 7, 8), common_term_92);
    result_by_q.set(c_Key1(8), common_term_92);
    // q = 9
    result_by_lpq.set(c_Key3(10, 7, 9), common_term_91);
    result_by_q.set(c_Key1(9), common_term_91);
    // q = 10
    result_by_lpq.set(c_Key3(10, 7, 10), common_term_90);
    result_by_q.set(c_Key1(10), common_term_90);
    // q = 11
    result_by_lpq.set(c_Key3(10, 7, 11), common_term_89);
    result_by_q.set(c_Key1(11), common_term_89);
    // q = 12
    result_by_lpq.set(c_Key3(10, 7, 12), common_term_88);
    result_by_q.set(c_Key1(12), common_term_88);
    // q = 13
    result_by_lpq.set(c_Key3(10, 7, 13), common_term_87);
    result_by_q.set(c_Key1(13), common_term_87);
    // q = 14
    result_by_lpq.set(c_Key3(10, 7, 14), common_term_86);
    result_by_q.set(c_Key1(14), common_term_86);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 7), result_by_q);
    result_by_q.clear();

    // l , p = (10, 8).
    // q = -14
    result_by_lpq.set(c_Key3(10, 8, -14), common_term_85);
    result_by_q.set(c_Key1(-14), common_term_85);
    // q = -13
    result_by_lpq.set(c_Key3(10, 8, -13), common_term_84);
    result_by_q.set(c_Key1(-13), common_term_84);
    // q = -12
    result_by_lpq.set(c_Key3(10, 8, -12), common_term_83);
    result_by_q.set(c_Key1(-12), common_term_83);
    // q = -11
    result_by_lpq.set(c_Key3(10, 8, -11), common_term_82);
    result_by_q.set(c_Key1(-11), common_term_82);
    // q = -10
    result_by_lpq.set(c_Key3(10, 8, -10), common_term_81);
    result_by_q.set(c_Key1(-10), common_term_81);
    // q = -9
    result_by_lpq.set(c_Key3(10, 8, -9), common_term_80);
    result_by_q.set(c_Key1(-9), common_term_80);
    // q = -8
    result_by_lpq.set(c_Key3(10, 8, -8), common_term_79);
    result_by_q.set(c_Key1(-8), common_term_79);
    // q = -7
    result_by_lpq.set(c_Key3(10, 8, -7), common_term_78);
    result_by_q.set(c_Key1(-7), common_term_78);
    // q = -6
    result_by_lpq.set(c_Key3(10, 8, -6), common_term_77);
    result_by_q.set(c_Key1(-6), common_term_77);
    // q = -5
    result_by_lpq.set(c_Key3(10, 8, -5), common_term_76);
    result_by_q.set(c_Key1(-5), common_term_76);
    // q = -4
    result_by_lpq.set(c_Key3(10, 8, -4), common_term_75);
    result_by_q.set(c_Key1(-4), common_term_75);
    // q = -3
    result_by_lpq.set(c_Key3(10, 8, -3), common_term_74);
    result_by_q.set(c_Key1(-3), common_term_74);
    // q = -2
    result_by_lpq.set(c_Key3(10, 8, -2), common_term_73);
    result_by_q.set(c_Key1(-2), common_term_73);
    // q = -1
    result_by_lpq.set(c_Key3(10, 8, -1), common_term_72);
    result_by_q.set(c_Key1(-1), common_term_72);
    // q = 0
    result_by_lpq.set(c_Key3(10, 8, 0), common_term_71);
    result_by_q.set(c_Key1(0), common_term_71);
    // q = 1
    result_by_lpq.set(c_Key3(10, 8, 1), common_term_70);
    result_by_q.set(c_Key1(1), common_term_70);
    // q = 2
    result_by_lpq.set(c_Key3(10, 8, 2), common_term_69);
    result_by_q.set(c_Key1(2), common_term_69);
    // q = 3
    result_by_lpq.set(c_Key3(10, 8, 3), common_term_68);
    result_by_q.set(c_Key1(3), common_term_68);
    // q = 4
    result_by_lpq.set(c_Key3(10, 8, 4), common_term_67);
    result_by_q.set(c_Key1(4), common_term_67);
    // q = 5
    result_by_lpq.set(c_Key3(10, 8, 5), common_term_66);
    result_by_q.set(c_Key1(5), common_term_66);
    // q = 6
    result_by_lpq.set(c_Key3(10, 8, 6), common_term_65);
    result_by_q.set(c_Key1(6), common_term_65);
    // q = 7
    result_by_lpq.set(c_Key3(10, 8, 7), common_term_64);
    result_by_q.set(c_Key1(7), common_term_64);
    // q = 8
    result_by_lpq.set(c_Key3(10, 8, 8), common_term_63);
    result_by_q.set(c_Key1(8), common_term_63);
    // q = 9
    result_by_lpq.set(c_Key3(10, 8, 9), common_term_62);
    result_by_q.set(c_Key1(9), common_term_62);
    // q = 10
    result_by_lpq.set(c_Key3(10, 8, 10), common_term_61);
    result_by_q.set(c_Key1(10), common_term_61);
    // q = 11
    result_by_lpq.set(c_Key3(10, 8, 11), common_term_60);
    result_by_q.set(c_Key1(11), common_term_60);
    // q = 12
    result_by_lpq.set(c_Key3(10, 8, 12), common_term_59);
    result_by_q.set(c_Key1(12), common_term_59);
    // q = 13
    result_by_lpq.set(c_Key3(10, 8, 13), common_term_58);
    result_by_q.set(c_Key1(13), common_term_58);
    // q = 14
    result_by_lpq.set(c_Key3(10, 8, 14), common_term_57);
    result_by_q.set(c_Key1(14), common_term_57);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 8), result_by_q);
    result_by_q.clear();

    // l , p = (10, 9).
    // q = -14
    result_by_lpq.set(c_Key3(10, 9, -14), common_term_56);
    result_by_q.set(c_Key1(-14), common_term_56);
    // q = -13
    result_by_lpq.set(c_Key3(10, 9, -13), common_term_55);
    result_by_q.set(c_Key1(-13), common_term_55);
    // q = -12
    result_by_lpq.set(c_Key3(10, 9, -12), common_term_54);
    result_by_q.set(c_Key1(-12), common_term_54);
    // q = -11
    result_by_lpq.set(c_Key3(10, 9, -11), common_term_53);
    result_by_q.set(c_Key1(-11), common_term_53);
    // q = -10
    result_by_lpq.set(c_Key3(10, 9, -10), common_term_52);
    result_by_q.set(c_Key1(-10), common_term_52);
    // q = -9
    result_by_lpq.set(c_Key3(10, 9, -9), common_term_51);
    result_by_q.set(c_Key1(-9), common_term_51);
    // q = -8
    result_by_lpq.set(c_Key3(10, 9, -8), common_term_50);
    result_by_q.set(c_Key1(-8), common_term_50);
    // q = -7
    result_by_lpq.set(c_Key3(10, 9, -7), common_term_49);
    result_by_q.set(c_Key1(-7), common_term_49);
    // q = -6
    result_by_lpq.set(c_Key3(10, 9, -6), common_term_48);
    result_by_q.set(c_Key1(-6), common_term_48);
    // q = -5
    result_by_lpq.set(c_Key3(10, 9, -5), common_term_47);
    result_by_q.set(c_Key1(-5), common_term_47);
    // q = -4
    result_by_lpq.set(c_Key3(10, 9, -4), common_term_46);
    result_by_q.set(c_Key1(-4), common_term_46);
    // q = -3
    result_by_lpq.set(c_Key3(10, 9, -3), common_term_45);
    result_by_q.set(c_Key1(-3), common_term_45);
    // q = -2
    result_by_lpq.set(c_Key3(10, 9, -2), common_term_44);
    result_by_q.set(c_Key1(-2), common_term_44);
    // q = -1
    result_by_lpq.set(c_Key3(10, 9, -1), common_term_43);
    result_by_q.set(c_Key1(-1), common_term_43);
    // q = 0
    result_by_lpq.set(c_Key3(10, 9, 0), common_term_42);
    result_by_q.set(c_Key1(0), common_term_42);
    // q = 1
    result_by_lpq.set(c_Key3(10, 9, 1), common_term_41);
    result_by_q.set(c_Key1(1), common_term_41);
    // q = 2
    result_by_lpq.set(c_Key3(10, 9, 2), common_term_40);
    result_by_q.set(c_Key1(2), common_term_40);
    // q = 3
    result_by_lpq.set(c_Key3(10, 9, 3), common_term_39);
    result_by_q.set(c_Key1(3), common_term_39);
    // q = 4
    result_by_lpq.set(c_Key3(10, 9, 4), common_term_38);
    result_by_q.set(c_Key1(4), common_term_38);
    // q = 5
    result_by_lpq.set(c_Key3(10, 9, 5), common_term_37);
    result_by_q.set(c_Key1(5), common_term_37);
    // q = 6
    result_by_lpq.set(c_Key3(10, 9, 6), common_term_36);
    result_by_q.set(c_Key1(6), common_term_36);
    // q = 7
    result_by_lpq.set(c_Key3(10, 9, 7), common_term_35);
    result_by_q.set(c_Key1(7), common_term_35);
    // q = 8
    result_by_lpq.set(c_Key3(10, 9, 8), common_term_34);
    result_by_q.set(c_Key1(8), common_term_34);
    // q = 9
    result_by_lpq.set(c_Key3(10, 9, 9), common_term_33);
    result_by_q.set(c_Key1(9), common_term_33);
    // q = 10
    result_by_lpq.set(c_Key3(10, 9, 10), common_term_32);
    result_by_q.set(c_Key1(10), common_term_32);
    // q = 11
    result_by_lpq.set(c_Key3(10, 9, 11), common_term_31);
    result_by_q.set(c_Key1(11), common_term_31);
    // q = 12
    result_by_lpq.set(c_Key3(10, 9, 12), common_term_30);
    result_by_q.set(c_Key1(12), common_term_30);
    // q = 13
    result_by_lpq.set(c_Key3(10, 9, 13), common_term_29);
    result_by_q.set(c_Key1(13), common_term_29);
    // q = 14
    result_by_lpq.set(c_Key3(10, 9, 14), common_term_28);
    result_by_q.set(c_Key1(14), common_term_28);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 9), result_by_q);
    result_by_q.clear();

    // l , p = (10, 10).
    // q = -14
    result_by_lpq.set(c_Key3(10, 10, -14), common_term_27);
    result_by_q.set(c_Key1(-14), common_term_27);
    // q = -13
    result_by_lpq.set(c_Key3(10, 10, -13), common_term_26);
    result_by_q.set(c_Key1(-13), common_term_26);
    // q = -12
    result_by_lpq.set(c_Key3(10, 10, -12), common_term_25);
    result_by_q.set(c_Key1(-12), common_term_25);
    // q = -11
    result_by_lpq.set(c_Key3(10, 10, -11), common_term_24);
    result_by_q.set(c_Key1(-11), common_term_24);
    // q = -10
    result_by_lpq.set(c_Key3(10, 10, -10), common_term_23);
    result_by_q.set(c_Key1(-10), common_term_23);
    // q = -9
    result_by_lpq.set(c_Key3(10, 10, -9), common_term_22);
    result_by_q.set(c_Key1(-9), common_term_22);
    // q = -8
    result_by_lpq.set(c_Key3(10, 10, -8), common_term_21);
    result_by_q.set(c_Key1(-8), common_term_21);
    // q = -7
    result_by_lpq.set(c_Key3(10, 10, -7), common_term_20);
    result_by_q.set(c_Key1(-7), common_term_20);
    // q = -6
    result_by_lpq.set(c_Key3(10, 10, -6), common_term_19);
    result_by_q.set(c_Key1(-6), common_term_19);
    // q = -5
    result_by_lpq.set(c_Key3(10, 10, -5), common_term_18);
    result_by_q.set(c_Key1(-5), common_term_18);
    // q = -4
    result_by_lpq.set(c_Key3(10, 10, -4), common_term_17);
    result_by_q.set(c_Key1(-4), common_term_17);
    // q = -3
    result_by_lpq.set(c_Key3(10, 10, -3), common_term_16);
    result_by_q.set(c_Key1(-3), common_term_16);
    // q = -2
    result_by_lpq.set(c_Key3(10, 10, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(10, 10, -1), common_term_14);
    result_by_q.set(c_Key1(-1), common_term_14);
    // q = 0
    result_by_lpq.set(c_Key3(10, 10, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(10, 10, 1), common_term_12);
    result_by_q.set(c_Key1(1), common_term_12);
    // q = 2
    result_by_lpq.set(c_Key3(10, 10, 2), common_term_11);
    result_by_q.set(c_Key1(2), common_term_11);
    // q = 3
    result_by_lpq.set(c_Key3(10, 10, 3), common_term_10);
    result_by_q.set(c_Key1(3), common_term_10);
    // q = 4
    result_by_lpq.set(c_Key3(10, 10, 4), common_term_9);
    result_by_q.set(c_Key1(4), common_term_9);
    // q = 5
    result_by_lpq.set(c_Key3(10, 10, 5), common_term_8);
    result_by_q.set(c_Key1(5), common_term_8);
    // q = 6
    result_by_lpq.set(c_Key3(10, 10, 6), common_term_7);
    result_by_q.set(c_Key1(6), common_term_7);
    // q = 7
    result_by_lpq.set(c_Key3(10, 10, 7), common_term_6);
    result_by_q.set(c_Key1(7), common_term_6);
    // q = 8
    result_by_lpq.set(c_Key3(10, 10, 8), common_term_5);
    result_by_q.set(c_Key1(8), common_term_5);
    // q = 9
    result_by_lpq.set(c_Key3(10, 10, 9), common_term_4);
    result_by_q.set(c_Key1(9), common_term_4);
    // q = 11
    result_by_lpq.set(c_Key3(10, 10, 11), common_term_3);
    result_by_q.set(c_Key1(11), common_term_3);
    // q = 12
    result_by_lpq.set(c_Key3(10, 10, 12), common_term_2);
    result_by_q.set(c_Key1(12), common_term_2);
    // q = 13
    result_by_lpq.set(c_Key3(10, 10, 13), common_term_1);
    result_by_q.set(c_Key1(13), common_term_1);
    // q = 14
    result_by_lpq.set(c_Key3(10, 10, 14), common_term_0);
    result_by_q.set(c_Key1(14), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 10), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l10_e20(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 10.
    // Functions are taylor expanded and truncated at eccentricity^20.

    // Functions Calculated for l = 10.

    c_IntMap<c_Key3, double> result_by_lpq(427);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(11);
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
    double common_term_0 = 2.1180838095217927e-5*eccentricity_19;
    double common_term_1 = 1.0733437299125303e-5*eccentricity_18;
    double common_term_2 = 1.8781273333403877e-5*eccentricity_19 + 4.98985859780472e-6*eccentricity_17;
    double common_term_3 = 7.685025017983592e-6*eccentricity_18 + 2.0574082725310404e-6*eccentricity_16;
    double common_term_4 = 5.7591945761555538e-6*eccentricity_19 + 2.6151006415277616e-6*eccentricity_17 + 7.1219762152245422e-7*eccentricity_15;
    double common_term_5 = 1.4353665147315941e-6*eccentricity_18 + 6.7030754861442692e-7*eccentricity_16 + 1.8793669587320381e-7*eccentricity_14;
    double common_term_6 = 3.554446315899251e-7*eccentricity_19 + 2.1943399998708545e-7*eccentricity_17 + 1.0659857314817402e-7*eccentricity_15 + 3.1254031917789242e-8*eccentricity_13;
    double common_term_7 = 2.0109284449424441e-8*eccentricity_18 + 1.2989185546229197e-8*eccentricity_16 + 6.6645031922809701e-9*eccentricity_14 + 2.0876756987868099e-9*eccentricity_12;
    double common_term_8 = 1.2060196971832842e-10*eccentricity_19 + 9.3331065661345743e-11*eccentricity_17 + 6.4007306527956756e-11*eccentricity_15 + 3.5423208267989084e-11*eccentricity_13 + 1.2232474797578964e-11*eccentricity_11;
    double common_term_9 = -2.1477344726539044e-8*eccentricity_19 - 2.0069322628497502e-8*eccentricity_17 - 1.7923113392234972e-8*eccentricity_15 - 1.4833404751414192e-8*eccentricity_13 - 1.063002059909612e-8*eccentricity_11 - 5.3822889109347443e-9*eccentricity_9;
    double common_term_10 = 3.6413216453364023e-5*eccentricity_18 + 3.7014794592196213e-5*eccentricity_16 + 3.6877225462294907e-5*eccentricity_14 + 3.5445601851851852e-5*eccentricity_12 + 3.1690917107583774e-5*eccentricity_10 + 2.4801587301587302e-5*eccentricity_8;
    double common_term_11 = -0.00054751600500858712*eccentricity_19 - 0.00064538436982925862*eccentricity_17 - 0.00077047540769948595*eccentricity_15 - 0.00093127523149762835*eccentricity_13 - 0.0011785779680524554*eccentricity_11 - 0.0011653355189732143*eccentricity_9 - 0.0033900669642857143*eccentricity_7;
    double common_term_12 = -0.0054140144370660676*eccentricity_18 - 0.0062333520389075614*eccentricity_16 - 0.0069032921810699493*eccentricity_14 - 0.011360964138741917*eccentricity_12 + 0.0087301587301587302*eccentricity_10 - 0.08253968253968254*eccentricity_8 + 0.088888888888888889*eccentricity_6;
    double common_term_13 = -0.021692770731941891*eccentricity_19 - 0.026166704323142935*eccentricity_17 - 0.012011810926954924*eccentricity_15 - 0.14282798045246444*eccentricity_13 + 0.49970137379156849*eccentricity_11 - 1.6257876441592262*eccentricity_9 + 2.2040473090277778*eccentricity_7 - 0.81380208333333333*eccentricity_5;
    double common_term_14 = -0.12578388526927644*eccentricity_18 + 0.34545347377234916*eccentricity_16 - 2.40773856026784*eccentricity_14 + 8.9207728794642912*eccentricity_12 - 22.639620535714286*eccentricity_10 + 30.9234375*eccentricity_8 - 17.8875*eccentricity_6 + 3.375*eccentricity_4;
    double common_term_15 = -1.3459285472799216*eccentricity_19 + 6.9418129501327581*eccentricity_17 - 32.406303399228239*eccentricity_15 + 105.30185376997335*eccentricity_13 - 232.47241801509155*eccentricity_11 + 308.33914161964699*eccentricity_9 - 212.15309244791667*eccentricity_7 + 66.545572916666667*eccentricity_5 - 7.1458333333333333*eccentricity_3;
    double common_term_16 = 85.636322955785066*eccentricity_18 - 332.30942524071214*eccentricity_16 + 946.38402626212531*eccentricity_14 - 1881.2400396825392*eccentricity_12 + 2402.6212962962965*eccentricity_10 - 1784.6055555555556*eccentricity_8 + 704.5*eccentricity_6 - 129.33333333333333*eccentricity_4 + 8.0*eccentricity_2;
    double common_term_17 = 807.34050502441476*eccentricity_19 - 2728.1220725185879*eccentricity_17 + 6929.1780014584211*eccentricity_15 - 12623.381456181661*eccentricity_13 + 15486.722428436279*eccentricity_11 - 11881.492437744141*eccentricity_9 + 5306.87548828125*eccentricity_7 - 1260.4921875*eccentricity_5 + 135.5625*eccentricity_3 - 4.5*eccentricity;
    double common_term_18 = -18752.42272595497*eccentricity_18 + 43140.949632659524*eccentricity_16 - 72857.953256235602*eccentricity_14 + 85822.01504931037*eccentricity_12 - 66460.387586805555*eccentricity_10 + 31844.422743055556*eccentricity_8 - 8775.9548611111111*eccentricity_6 + 1241.5625*eccentricity_4 - 72.5*eccentricity_2 + 1.0;
    double common_term_19 = -111586.78942599429*eccentricity_19 + 235358.66485664334*eccentricity_17 - 371537.34974909366*eccentricity_15 + 420362.94863797451*eccentricity_13 - 324311.01459763702*eccentricity_11 + 161721.53108113607*eccentricity_9 - 48884.568739149307*eccentricity_7 + 8150.3411458333333*eccentricity_5 - 633.1875*eccentricity_3 + 15.5*eccentricity;
    double common_term_20 = 1150299.5211249165*eccentricity_18 - 1708119.2174701051*eccentricity_16 + 1857749.3601102121*eccentricity_14 - 1416734.9751897321*eccentricity_12 + 721640.82304687498*eccentricity_10 - 231454.8*eccentricity_8 + 42937.312500000001*eccentricity_6 - 3973.5*eccentricity_4 + 130.5*eccentricity_2;
    double common_term_21 = 5122131.8797949193*eccentricity_19 - 7190885.2407993091*eccentricity_17 + 7525559.9876692543*eccentricity_15 - 5645298.2614953222*eccentricity_13 + 2901880.221849056*eccentricity_11 - 966863.29034333404*eccentricity_9 + 192809.99469401042*eccentricity_7 - 20053.194010416667*eccentricity_5 + 790.77083333333333*eccentricity_3;
    double common_term_22 = -28061944.347515996*eccentricity_18 + 28290833.563823307*eccentricity_16 - 20811417.030628115*eccentricity_14 + 10707868.584195964*eccentricity_12 - 3653666.6625723379*eccentricity_10 + 765856.92907986111*eccentricity_8 - 86424.668750000001*eccentricity_6 + 3858.8958333333334*eccentricity_4;
    double common_term_23 = -102516456.18715353*eccentricity_19 + 99676225.161767525*eccentricity_17 - 71764218.113639083*eccentricity_15 + 36749951.747007455*eccentricity_13 - 12713999.931560244*eccentricity_11 + 2758410.2350289481*eccentricity_9 - 330019.78271484375*eccentricity_7 + 16100.61328125*eccentricity_5;
    double common_term_24 = 331781573.83610157*eccentricity_18 - 233508392.56055399*eccentricity_16 + 118538338.34825096*eccentricity_14 - 41281660.227906379*eccentricity_12 + 9168641.8137896825*eccentricity_10 - 1144440.9515873016*eccentricity_8 + 59584.609722222222*eccentricity_6;
    double common_term_25 = 1050239719.9822683*eccentricity_19 - 722039341.90925662*eccentricity_17 + 362308093.33481603*eccentricity_15 - 126341976.86238201*eccentricity_13 + 28493340.940017857*eccentricity_11 - 3667624.6796146938*eccentricity_9 + 200445.43370690724*eccentricity_7;
    double common_term_26 = -2134093960.9808699*eccentricity_18 + 1056309105.3580081*eccentricity_16 - 367374606.24462491*eccentricity_14 + 83618046.254594029*eccentricity_12 - 11003358.154520089*eccentricity_10 + 623706.16439732142*eccentricity_8;
    double common_term_27 = -6058570430.2100537*eccentricity_19 + 2953627505.2667947*eccentricity_17 - 1021447407.8414878*eccentricity_15 + 233544169.4545519*eccentricity_13 - 31210475.019879966*eccentricity_11 + 1818308.9840970532*eccentricity_9;
    double common_term_28 = 7956689158.8754096*eccentricity_18 - 2729824988.5338813*eccentricity_16 + 624697168.85578806*eccentricity_14 - 84346904.447915727*eccentricity_12 + 5015578.577405065*eccentricity_10;
    double common_term_29 = 20728626830.326469*eccentricity_19 - 7042782235.7584779*eccentricity_17 + 1608488960.9014967*eccentricity_15 - 218534200.51042232*eccentricity_13 + 13191026.055637727*eccentricity_11;
    double common_term_30 = -17604401987.420231*eccentricity_18 + 4003577809.9578803*eccentricity_16 - 545556461.6433747*eccentricity_14 + 33282081.127764854*eccentricity_12;
    double common_term_31 = -42766511936.392622*eccentricity_19 + 9667099187.8532077*eccentricity_17 - 1317766134.1665586*eccentricity_15 + 80963575.707463975*eccentricity_13;
    double common_term_32 = 22712386893.961056*eccentricity_18 - 3090497880.0471204*eccentricity_16 + 190682398.33831324*eccentricity_14;
    double common_term_33 = 52054779222.892677*eccentricity_19 - 7058137015.9968525*eccentricity_17 + 436288546.48853598*eccentricity_15;
    double common_term_34 = -15736806477.687374*eccentricity_18 + 972630557.49984691*eccentricity_16;
    double common_term_35 = -34328109442.710094*eccentricity_19 + 2117950208.1658236*eccentricity_17;
    double common_term_36 = 4514498221.1598429*eccentricity_18;
    double common_term_37 = 9437081151.2357576*eccentricity_19;
    double common_term_38 = 1.0288084466013607*eccentricity_19;
    double common_term_39 = 0.75674237506694758*eccentricity_18;
    double common_term_40 = 4.0345365441784855*eccentricity_19 + 0.55665432440688455*eccentricity_17;
    double common_term_41 = 3.0706137184289076*eccentricity_18 + 0.40949244320574138*eccentricity_16;
    double common_term_42 = 10.410449274432229*eccentricity_19 + 2.3344330100511606*eccentricity_17 + 0.30124999907608369*eccentricity_15;
    double common_term_43 = 8.1027108093192034*eccentricity_18 + 1.7729325072249179*eccentricity_16 + 0.22162890785274993*eccentricity_14;
    double common_term_44 = 21.946612279616421*eccentricity_19 + 6.2978953428953424*eccentricity_17 + 1.3451915608193879*eccentricity_15 + 0.16305740723786139*eccentricity_13;
    double common_term_45 = 17.37017869758595*eccentricity_18 + 4.8886469988915658*eccentricity_16 + 1.0197201563667492*eccentricity_14 + 0.11996800584799717*eccentricity_12;
    double common_term_46 = 40.878433451408526*eccentricity_19 + 13.728203373475821*eccentricity_17 + 3.7899445827149011*eccentricity_15 + 0.77233257646684523*eccentricity_13 + 0.088266719471324574*eccentricity_11;
    double common_term_47 = 32.788541895858795*eccentricity_18 + 10.834636520190629*eccentricity_16 + 2.9346062504018776*eccentricity_14 + 0.58448591194684945*eccentricity_12 + 0.064942887455908289*eccentricity_10;
    double common_term_48 = 69.929069254715182*eccentricity_19 + 26.261723244532788*eccentricity_17 + 8.5393436629686325*eccentricity_15 + 2.2696591440821298*eccentricity_13 + 0.44198625446837625*eccentricity_11 + 0.047782297824005181*eccentricity_9;
    double common_term_49 = 0.03515625*eccentricity_8*std::pow(1.0 - eccentricity_2, -9.5);
    double common_term_50 = 112.3524430982397*eccentricity_19 + 45.928325759795107*eccentricity_17 + 16.776562126780791*eccentricity_15 + 5.2837762153744347*eccentricity_13 + 1.3531424621635854*eccentricity_11 + 0.25219857352120536*eccentricity_9 + 0.025866505456349206*eccentricity_7;
    double common_term_51 = 91.974881293084915*eccentricity_18 + 37.143598064181832*eccentricity_16 + 13.381717950550687*eccentricity_14 + 4.1485331468621399*eccentricity_12 + 1.043226066468254*eccentricity_10 + 0.19037698412698413*eccentricity_8 + 0.019097222222222222*eccentricity_6;
    double common_term_52 = 171.97625909506005*eccentricity_19 + 75.189718712027777*eccentricity_17 + 29.997713000587055*eccentricity_15 + 10.65913347857339*eccentricity_13 + 3.2524068014962333*eccentricity_11 + 0.80217459542410714*eccentricity_9 + 0.14326171875*eccentricity_7 + 0.00703125*eccentricity_5;
    double common_term_53 = 141.92732060695361*eccentricity_18 + 61.38062436436364*eccentricity_16 + 24.189023784837044*eccentricity_14 + 8.4745465959821429*eccentricity_12 + 2.5287037037037037*eccentricity_10 + 0.66111111111111111*eccentricity_8 - 0.06875*eccentricity_6 + 0.14583333333333333*eccentricity_4;
    double common_term_54 = 253.22851290423462*eccentricity_19 + 116.95918901036959*eccentricity_17 + 50.017628494379288*eccentricity_15 + 19.503868932446475*eccentricity_13 + 6.4177327433591166*eccentricity_11 + 3.101062915943287*eccentricity_9 - 2.3777669270833333*eccentricity_7 + 3.0247395833333333*eccentricity_5 - 0.85416666666666667*eccentricity_3;
    double common_term_55 = 210.44835771397184*eccentricity_18 + 96.069572305484702*eccentricity_16 + 41.58755545479911*eccentricity_14 + 11.053727678571429*eccentricity_12 + 20.2921875*eccentricity_10 - 30.20625*eccentricity_8 + 35.015625*eccentricity_6 - 15.75*eccentricity_4 + 2.25*eccentricity_2;
    double common_term_56 = 361.43738721090084*eccentricity_19 + 172.24542382768535*eccentricity_17 + 91.274410893585868*eccentricity_15 - 16.790848970530973*eccentricity_13 + 153.76848819873951*eccentricity_11 - 260.62958780924479*eccentricity_9 + 288.20958116319444*eccentricity_7 - 155.44010416666667*eccentricity_5 + 36.1875*eccentricity_3 - 2.5*eccentricity;
    double common_term_57 = 275.00062984554779*eccentricity_18 + 265.52491119300234*eccentricity_16 - 352.56799445704168*eccentricity_14 + 1067.3436419753088*eccentricity_12 - 1758.5455555555556*eccentricity_10 + 1866.3055555555556*eccentricity_8 - 1089.5486111111111*eccentricity_6 + 315.125*eccentricity_4 - 36.5*eccentricity_2 + 1.0;
    double common_term_58 = 259.81781863683375*eccentricity_19 + 1195.579425712389*eccentricity_17 - 2729.5336153178436*eccentricity_15 + 6407.1587927218836*eccentricity_13 - 9903.7892381286627*eccentricity_11 + 10084.092352294922*eccentricity_9 - 6078.41455078125*eccentricity_7 + 1970.2265625*eccentricity_5 - 290.8125*eccentricity_3 + 13.5*eccentricity;
    double common_term_59 = 6562.2249977234656*eccentricity_18 - 16323.953564591196*eccentricity_16 + 33541.5456548351*eccentricity_14 - 48417.402591765872*eccentricity_12 + 47262.871817129631*eccentricity_10 - 28695.711805555556*eccentricity_8 + 9901.328125*eccentricity_6 - 1675.0833333333333*eccentricity_4 + 100.25*eccentricity_2;
    double common_term_60 = 35220.642963427729*eccentricity_19 - 83740.935562125634*eccentricity_17 + 156136.48593028098*eccentricity_15 - 211038.32953387519*eccentricity_13 + 197461.10533111007*eccentricity_11 - 119081.11242088035*eccentricity_9 + 42477.821712239583*eccentricity_7 - 7802.2669270833333*eccentricity_5 + 541.47916666666667*eccentricity_3;
    double common_term_61 = -382886.57819176992*eccentricity_18 + 658062.75218994142*eccentricity_16 - 836244.99426269532*eccentricity_14 + 750256.73368443079*eccentricity_12 - 445681.59910714286*eccentricity_10 + 161427.2625*eccentricity_8 - 31188.09375*eccentricity_6 + 2376.5625*eccentricity_4;
    double common_term_62 = -1593926.3228540845*eccentricity_19 + 2548625.8142908708*eccentricity_17 - 3057232.2783962995*eccentricity_15 + 2631957.8339697471*eccentricity_13 - 1532288.5636892672*eccentricity_11 + 557026.68569394067*eccentricity_9 - 110947.63989800347*eccentricity_7 + 8987.4486979166667*eccentricity_5;
    double common_term_63 = 9178882.1462976283*eccentricity_18 - 10431492.289385038*eccentricity_16 + 8625578.5376998437*eccentricity_14 - 4906074.2055129565*eccentricity_12 + 1776120.1692491319*eccentricity_10 - 359842.25486111111*eccentricity_8 + 30349.348263888889*eccentricity_6;
    double common_term_64 = 31037421.969894867*eccentricity_19 - 33525279.842521398*eccentricity_17 + 26654867.601642882*eccentricity_15 - 14783356.235550846*eccentricity_13 + 5301464.0178067343*eccentricity_11 - 1082391.1435743059*eccentricity_9 + 93710.231794084822*eccentricity_7;
    double common_term_65 = -102243859.86769951*eccentricity_18 + 78252604.649985213*eccentricity_16 - 42272175.433449653*eccentricity_14 + 14960518.251003327*eccentricity_12 - 3057694.4276902833*eccentricity_10 + 269031.71821676587*eccentricity_8;
    double common_term_66 = -297713606.94825665*eccentricity_19 + 219592801.34992108*eccentricity_17 - 115469875.23999036*eccentricity_15 + 40224855.413677096*eccentricity_13 - 8190235.4274458904*eccentricity_11 + 726995.90431837725*eccentricity_9;
    double common_term_67 = 592027359.1881786*eccentricity_18 - 302953343.28624617*eccentricity_16 + 103689654.78311793*eccentricity_14 - 20957317.37134791*eccentricity_12 + 1866483.2682421875*eccentricity_10;
    double common_term_68 = 1539999274.1544579*eccentricity_19 - 766886207.75928423*eccentricity_17 + 257551566.72457339*eccentricity_15 - 51534719.673596271*eccentricity_13 + 4586053.3352122621*eccentricity_11;
    double common_term_69 = -1880067356.4129967*eccentricity_18 + 619004518.39645683*eccentricity_16 - 122374212.7556137*eccentricity_14 + 10846586.516443591*eccentricity_12;
    double common_term_70 = -4478081538.4769334*eccentricity_19 + 1444577850.2083849*eccentricity_17 - 281732359.2905232*eccentricity_15 + 24809937.971262643*eccentricity_13;
    double common_term_71 = 3283131451.112348*eccentricity_18 - 630935770.66319587*eccentricity_16 + 55095376.82915985*eccentricity_14;
    double common_term_72 = 7285024950.3031*eccentricity_19 - 1378334101.8391283*eccentricity_17 + 119168183.4883343*eccentricity_15;
    double common_term_73 = -2944313061.3853878*eccentricity_18 + 251732980.36485921*eccentricity_16;
    double common_term_74 = -6162642757.3730721*eccentricity_19 + 520544682.92290683*eccentricity_17;
    double common_term_75 = 1055785181.4203645*eccentricity_18;
    double common_term_76 = 2103963665.0103699*eccentricity_19;
    double common_term_77 = 112.98191736461063*eccentricity_19;
    double common_term_78 = 80.737549116004609*eccentricity_18;
    double common_term_79 = 368.37186307628215*eccentricity_19 + 57.64650582506177*eccentricity_17;
    double common_term_80 = 274.54455524536473*eccentricity_18 + 41.122295597671707*eccentricity_16;
    double common_term_81 = 843.77873716841221*eccentricity_19 + 204.10867910842715*eccentricity_17 + 29.306396127461009*eccentricity_15;
    double common_term_82 = 644.53680197606078*eccentricity_18 + 151.37887298870417*eccentricity_16 + 20.864082391798476*eccentricity_14;
    double common_term_83 = 1625.7041565383165*eccentricity_19 + 491.06340377946205*eccentricity_17 + 112.0079527130592*eccentricity_15 + 14.837341078135081*eccentricity_13;
    double common_term_84 = 1264.2089183228271*eccentricity_18 + 373.17328427841341*eccentricity_16 + 82.685478052221216*eccentricity_14 + 10.538954317167208*eccentricity_12;
    double common_term_85 = 2816.4034988855777*eccentricity_19 + 980.58515516862394*eccentricity_17 + 282.86288363916402*eccentricity_15 + 60.899847742333952*eccentricity_13 + 7.4762625074241336*eccentricity_11;
    double common_term_86 = 2221.0071568959703*eccentricity_18 + 758.65423937690503*eccentricity_16 + 213.86499566859484*eccentricity_14 + 44.752017909752285*eccentricity_12 + 5.2963146219135802*eccentricity_10;
    double common_term_87 = 4534.4679592122625*eccentricity_19 + 1747.1597537561566*eccentricity_17 + 585.4577306037522*eccentricity_15 + 161.28867905393823*eccentricity_13 + 32.810629239763532*eccentricity_11 + 3.7464246477399554*eccentricity_9;
    double common_term_88 = 3617.1157281978762*eccentricity_18 + 1371.0089774932163*eccentricity_16 + 450.6495398614013*eccentricity_14 + 121.32952799479167*eccentricity_12 + 24.000000688932981*eccentricity_10 + 2.6458209325396825*eccentricity_8;
    double common_term_89 = 6915.8077877422675*eccentricity_19 + 2878.5201308646657*eccentricity_17 + 1073.1708764394303*eccentricity_15 + 345.99517598736363*eccentricity_13 + 91.03766524770693*eccentricity_11 + 17.51398184640067*eccentricity_9 + 1.8652793278769841*eccentricity_7;
    double common_term_90 = std::pow(1.0 - eccentricity_2, -9.5)*(0.28125*eccentricity_8 + 1.3125*eccentricity_6);
    double common_term_91 = 10114.631608995727*eccentricity_19 + 4476.5345796185755*eccentricity_17 + 1809.964754093561*eccentricity_15 + 652.62349498479385*eccentricity_13 + 202.37828039946379*eccentricity_11 + 50.857536388578869*eccentricity_9 + 9.2587782118055556*eccentricity_7 + 0.92161458333333333*eccentricity_5;
    double common_term_92 = 8215.2922533077255*eccentricity_18 + 3589.2510670154294*eccentricity_16 + 1430.0388710656186*eccentricity_14 + 507.00156637524802*eccentricity_12 + 154.17078579695767*eccentricity_10 + 37.861371527777778*eccentricity_8 + 6.70625*eccentricity_6 + 0.64583333333333333*eccentricity_4;
    double common_term_93 = 14304.43332593057*eccentricity_19 + 6658.1018719744533*eccentricity_17 + 2871.1891703432798*eccentricity_15 + 1127.1018888698305*eccentricity_13 + 392.86032180786133*eccentricity_11 + 117.13199462890625*eccentricity_9 + 28.10771484375*eccentricity_7 + 4.84765625*eccentricity_5 + 0.4375*eccentricity_3;
    double common_term_94 = 11703.836967576297*eccentricity_18 + 5384.2165513332712*eccentricity_16 + 2291.4185359657463*eccentricity_14 + 886.1375037202381*eccentricity_12 + 303.62125289351852*eccentricity_10 + 88.711111111111111*eccentricity_8 + 20.9375*eccentricity_6 + 3.1666666666666667*eccentricity_4 + 0.5*eccentricity_2;
    double common_term_95 = 19678.975043607169*eccentricity_19 + 9556.0406493825796*eccentricity_17 + 4344.3784394099216*eccentricity_15 + 1824.3681644148717*eccentricity_13 + 695.03978884661639*eccentricity_11 + 233.38388061523438*eccentricity_9 + 69.138726128472222*eccentricity_7 + 11.174479166666667*eccentricity_5 + 5.8125*eccentricity_3 - 0.5*eccentricity;
    double common_term_96 = 16206.502870589276*eccentricity_18 + 7785.923318780393*eccentricity_16 + 3497.1331227080676*eccentricity_14 + 1450.773330078125*eccentricity_12 + 536.04703125*eccentricity_10 + 201.1015625*eccentricity_8 + 13.65625*eccentricity_6 + 42.5625*eccentricity_4 - 8.5*eccentricity_2 + 1.0;
    double common_term_97 = 26453.149897806083*eccentricity_19 + 13320.544781455107*eccentricity_17 + 6326.1527697419016*eccentricity_15 + 2827.4208128707792*eccentricity_13 + 1083.1575818831832*eccentricity_11 + 582.28132527669271*eccentricity_9 - 107.68614366319444*eccentricity_7 + 245.67447916666667*eccentricity_5 - 71.4375*eccentricity_3 + 11.5*eccentricity;
    double common_term_98 = 21919.980312267055*eccentricity_18 + 10889.303889415977*eccentricity_16 + 5279.929886153825*eccentricity_14 + 1811.311626984127*eccentricity_12 + 1850.3712962962963*eccentricity_10 - 998.73888888888889*eccentricity_8 + 1195.375*eccentricity_6 - 415.33333333333333*eccentricity_4 + 74.0*eccentricity_2;
    double common_term_99 = 34923.580993740338*eccentricity_19 + 17851.89538786889*eccentricity_17 + 9871.6870311158891*eccentricity_15 + 1684.3778098848886*eccentricity_13 + 6538.4422776358468*eccentricity_11 - 5407.3172241210938*eccentricity_9 + 5068.97138671875*eccentricity_7 - 1908.99609375*eccentricity_5 + 351.1875*eccentricity_3;
    double common_term_100 = 27398.194515506151*eccentricity_18 + 20008.281768129032*eccentricity_16 - 4689.4129891998266*eccentricity_14 + 24165.240187872024*eccentricity_12 - 23403.038917824074*eccentricity_10 + 19177.660590277778*eccentricity_8 - 7436.5625*eccentricity_6 + 1369.9583333333333*eccentricity_4;
    double common_term_101 = 36455.424688022926*eccentricity_19 + 48257.145827312007*eccentricity_17 - 39318.263637027362*eccentricity_15 + 87843.350422844994*eccentricity_13 - 88156.786500256655*eccentricity_11 + 65947.025514439174*eccentricity_9 - 25597.524945746528*eccentricity_7 + 4649.3778645833333*eccentricity_5;
    double common_term_102 = 140939.50338318188*eccentricity_18 - 181306.44723457793*eccentricity_16 + 304441.09340401785*eccentricity_14 - 300270.63828125*eccentricity_12 + 209230.68113839286*eccentricity_10 - 80002.971428571429*eccentricity_8 + 14206.575*eccentricity_6;
    double common_term_103 = 460078.34762914776*eccentricity_19 - 690886.21944926894*eccentricity_17 + 996598.25324372381*eccentricity_15 - 945071.06842828232*eccentricity_13 + 619928.1896216586*eccentricity_11 - 231366.43252229236*eccentricity_9 + 39975.869204179067*eccentricity_7;
    double common_term_104 = -2366625.8994805733*eccentricity_18 + 3083403.2960288772*eccentricity_16 - 2787795.5507886884*eccentricity_14 + 1732188.5863252315*eccentricity_12 - 627670.41844135803*eccentricity_10 + 105238.97986111111*eccentricity_8;
    double common_term_105 = -7521529.0678065039*eccentricity_19 + 9054229.1342142199*eccentricity_17 - 7785523.9640435582*eccentricity_15 + 4600908.8919337972*eccentricity_13 - 1613867.4114398956*eccentricity_11 + 262199.80908857073*eccentricity_9;
    double common_term_106 = 25361397.390881119*eccentricity_18 - 20742423.997484287*eccentricity_16 + 11692765.550166911*eccentricity_14 - 3964207.7674797078*eccentricity_12 + 623657.13672674163*eccentricity_10;
    double common_term_107 = 68095320.0149021*eccentricity_19 - 53037294.209015826*eccentricity_17 + 28585948.363469411*eccentricity_15 - 9360904.4284804557*eccentricity_13 + 1425768.8055109258*eccentricity_11;
    double common_term_108 = -130785441.92147947*eccentricity_18 + 67529760.1599794*eccentricity_16 - 21356956.851246468*eccentricity_14 + 3149673.7737829748*eccentricity_12;
    double common_term_109 = -312270550.99270186*eccentricity_19 + 154731480.30733839*eccentricity_17 - 47272418.196230349*eccentricity_15 + 6752590.3413372089*eccentricity_13;
    double common_term_110 = 344973743.08600822*eccentricity_18 - 101859771.81369891*eccentricity_16 + 14099386.647976717*eccentricity_14;
    double common_term_111 = 750407602.20340997*eccentricity_19 - 214271185.23442369*eccentricity_17 + 28756167.629193573*eccentricity_15;
    double common_term_112 = -441103737.0843552*eccentricity_18 + 57429436.20005915*eccentricity_16;
    double common_term_113 = -890493160.4594087*eccentricity_19 + 112543461.87461006*eccentricity_17;
    double common_term_114 = 216803718.57770425*eccentricity_18;
    double common_term_115 = 411193102.24249776*eccentricity_19;
    double common_term_116 = 4138.6405281912606*eccentricity_19;
    double common_term_117 = 2829.3642413238161*eccentricity_18;
    double common_term_118 = 9835.7341477788659*eccentricity_19 + 1928.4650025296497*eccentricity_17;
    double common_term_119 = 7105.4019251714062*eccentricity_18 + 1310.0997631262258*eccentricity_16;
    double common_term_120 = 18680.067123371184*eccentricity_19 + 5102.0903109296229*eccentricity_17 + 886.8004370999841*eccentricity_15;
    double common_term_121 = 13823.700525493026*eccentricity_18 + 3641.5318466166829*eccentricity_16 + 597.8781373447734*eccentricity_14;
    double common_term_122 = 31091.857706138587*eccentricity_19 + 10171.974739526181*eccentricity_17 + 2583.1913630254973*eccentricity_15 + 401.30327557254624*eccentricity_13;
    double common_term_123 = 23400.328077968664*eccentricity_18 + 7440.906681033873*eccentricity_16 + 1820.8968782889725*eccentricity_14 + 268.0270107093588*eccentricity_12;
    double common_term_124 = 47665.459288129598*eccentricity_19 + 17513.57750786246*eccentricity_17 + 5409.6519473549704*eccentricity_15 + 1275.1064597254856*eccentricity_13 + 178.01548490179642*eccentricity_11;
    double common_term_125 = 36331.264339506515*eccentricity_18 + 13031.2340158328*eccentricity_16 + 3907.4299074168019*eccentricity_14 + 886.67794693587662*eccentricity_12 + 117.48426897321429*eccentricity_10;
    double common_term_126 = 68994.15701440329*eccentricity_19 + 27540.031502165887*eccentricity_17 + 9636.2947931303867*eccentricity_15 + 2802.9247542444055*eccentricity_13 + 611.94196460528769*eccentricity_11 + 76.972772695820588*eccentricity_9;
    double common_term_127 = 53116.720900924527*eccentricity_18 + 20754.732550441192*eccentricity_16 + 7079.0373773140841*eccentricity_14 + 1995.7419670758929*eccentricity_12 + 418.86314966380071*eccentricity_10 + 50.005716765873016*eccentricity_8;
    double common_term_128 = 95676.24463123057*eccentricity_19 + 40669.822642261914*eccentricity_17 + 15544.392444987491*eccentricity_15 + 5163.7462441757747*eccentricity_13 + 1409.5655752454485*eccentricity_11 + 284.08580278669085*eccentricity_9 + 32.164327566964286*eccentricity_7;
    double common_term_129 = 74260.695389085035*eccentricity_18 + 30958.595060792281*eccentricity_16 + 11564.737798184547*eccentricity_14 + 3737.825384998714*eccentricity_12 + 986.71463603670635*eccentricity_10 + 190.68382936507936*eccentricity_8 + 20.443055555555556*eccentricity_6;
    double common_term_130 = 128310.71458485333*eccentricity_19 + 57324.178242893404*eccentricity_17 + 23419.340772264431*eccentricity_15 + 8542.0413697642904*eccentricity_13 + 2682.9199675282473*eccentricity_11 + 683.84254847935268*eccentricity_9 + 126.46257595486111*eccentricity_7 + 12.804947916666667*eccentricity_5;
    double common_term_131 = std::pow(1.0 - eccentricity_2, -9.5)*(0.984375*eccentricity_8 + 7.875*eccentricity_6 + 7.875*eccentricity_4);
    double common_term_132 = 167493.12306250129*eccentricity_19 + 77923.329423605725*eccentricity_17 + 33547.418077676931*eccentricity_15 + 13124.511923522621*eccentricity_13 + 4547.1542784473883*eccentricity_11 + 1342.1533343279803*eccentricity_9 + 316.81741536458333*eccentricity_7 + 53.139322916666667*eccentricity_5 + 4.7291666666666667*eccentricity_3;
    double common_term_133 = 131635.37311374896*eccentricity_18 + 60201.534750570721*eccentricity_16 + 25407.207472106826*eccentricity_14 + 9709.4953354414683*eccentricity_12 + 3270.7477068865741*eccentricity_10 + 932.75069444444444*eccentricity_8 + 210.84375*eccentricity_6 + 33.416666666666667*eccentricity_4 + 2.75*eccentricity_2;
    double common_term_134 = 213811.39227449093*eccentricity_19 + 102882.72085181049*eccentricity_17 + 46212.40254205902*eccentricity_15 + 19097.12023960386*eccentricity_13 + 7117.9510542297363*eccentricity_11 + 2326.2882019042969*eccentricity_9 + 638.93115234375*eccentricity_7 + 137.6015625*eccentricity_5 + 20.4375*eccentricity_3 + 1.5*eccentricity;
    double common_term_135 = 168857.30383701822*eccentricity_18 + 79929.617104228148*eccentricity_16 + 35225.259770100604*eccentricity_14 + 14233.868095823688*eccentricity_12 + 5164.4587934027778*eccentricity_10 + 1632.9787326388889*eccentricity_8 + 429.99305555555556*eccentricity_6 + 87.875*eccentricity_4 + 11.5*eccentricity_2 + 1.0;
    double common_term_136 = 267841.54221468846*eccentricity_19 + 132609.14588536614*eccentricity_17 + 61692.110116913397*eccentricity_15 + 26642.034517889511*eccentricity_13 + 10508.832178045202*eccentricity_11 + 3702.8700052897135*eccentricity_9 + 1127.3065863715278*eccentricity_7 + 286.43489583333333*eccentricity_5 + 48.9375*eccentricity_3 + 9.5*eccentricity;
    double common_term_137 = 212412.06817026504*eccentricity_18 + 103508.73039251834*eccentricity_16 + 47273.520184151786*eccentricity_14 + 19974.563219866071*eccentricity_12 + 7677.5900390625*eccentricity_10 + 2607.13125*eccentricity_8 + 789.375*eccentricity_6 + 147.75*eccentricity_4 + 51.75*eccentricity_2;
    double common_term_138 = 330143.34347561482*eccentricity_19 + 167496.64589778274*eccentricity_17 + 80255.781040617116*eccentricity_15 + 35929.081201025822*eccentricity_13 + 14851.215758062292*eccentricity_11 + 5465.05249520761*eccentricity_9 + 1960.3046549479167*eccentricity_7 + 333.49348958333333*eccentricity_5 + 211.89583333333333*eccentricity_3;
    double common_term_139 = 262760.98954544106*eccentricity_18 + 131267.17549175998*eccentricity_16 + 61724.547504213973*eccentricity_14 + 27210.794546750992*eccentricity_12 + 10508.673544973545*eccentricity_10 + 4597.7611111111111*eccentricity_8 + 506.675*eccentricity_6 + 724.08333333333333*eccentricity_4;
    double common_term_140 = 401243.84246810321*eccentricity_19 + 207983.05790368313*eccentricity_17 + 101915.08377083501*eccentricity_15 + 47884.189076154573*eccentricity_13 + 18473.2914163317*eccentricity_11 + 10579.168385532924*eccentricity_9 + 65.10849609375*eccentricity_7 + 2179.52578125*eccentricity_5;
    double common_term_141 = 320701.76416126877*eccentricity_18 + 162251.0968813113*eccentricity_16 + 82251.798673904482*eccentricity_14 + 28825.261025224133*eccentricity_12 + 24545.98294890873*eccentricity_10 - 3120.3387896825397*eccentricity_8 + 5967.9826388888889*eccentricity_6;
    double common_term_142 = 483470.11189502292*eccentricity_19 + 248725.39500089955*eccentricity_17 + 141041.9657555091*eccentricity_15 + 35951.121408174835*eccentricity_13 + 58159.198572374105*eccentricity_11 - 15211.783884248279*eccentricity_9 + 15180.904791356647*eccentricity_7;
    double common_term_143 = 363403.99318024883*eccentricity_18 + 249061.24489227691*eccentricity_16 + 18466.207672737418*eccentricity_14 + 140388.0997265625*eccentricity_12 - 52121.822209821429*eccentricity_10 + 36398.464955357143*eccentricity_8;
    double common_term_144 = 490502.0887191879*eccentricity_19 + 469059.31964888346*eccentricity_17 - 88489.882787471435*eccentricity_15 + 341351.77715522546*eccentricity_13 - 151901.36367563628*eccentricity_11 + 83129.35033802389*eccentricity_9;
    double common_term_145 = 964202.03209421113*eccentricity_18 - 460999.53504451267*eccentricity_16 + 825479.23414268188*eccentricity_14 - 400520.69452809344*eccentricity_12 + 182279.58837191358*eccentricity_10;
    double common_term_146 = 2154675.1202137102*eccentricity_19 - 1544600.2023080229*eccentricity_17 + 1966077.6769274138*eccentricity_15 - 983804.59866731663*eccentricity_13 + 386077.27742586507*eccentricity_11;
    double common_term_147 = -4410097.7830178756*eccentricity_18 + 4584985.0398818621*eccentricity_16 - 2289259.6301116501*eccentricity_14 + 793676.63812026724*eccentricity_12;
    double common_term_148 = -11528251.034230587*eccentricity_19 + 10440907.745321757*eccentricity_17 - 5101189.050751141*eccentricity_15 + 1589705.6410422374*eccentricity_13;
    double common_term_149 = 23201780.027656546*eccentricity_18 - 10967225.474120779*eccentricity_16 + 3112121.1563314852*eccentricity_14;
    double common_term_150 = 50341647.989848946*eccentricity_19 - 22875061.319972804*eccentricity_17 + 5970216.8235497252*eccentricity_15;
    double common_term_151 = -46482933.247512121*eccentricity_18 + 11247626.674198625*eccentricity_16;
    double common_term_152 = -92326682.164787853*eccentricity_19 + 20848109.022398869*eccentricity_17;
    double common_term_153 = 38079252.988601268*eccentricity_18;
    double common_term_154 = 68629721.769822079*eccentricity_19;
    double common_term_155 = 80509.808802793336*eccentricity_19;
    double common_term_156 = 52326.485955429111*eccentricity_18;
    double common_term_157 = 114173.64627719518*eccentricity_19 + 33817.842897087828*eccentricity_17;
    double common_term_158 = 80877.105913393353*eccentricity_18 + 21720.709475718762*eccentricity_16;
    double common_term_159 = 175380.05996729813*eccentricity_19 + 56540.645051780581*eccentricity_17 + 13855.339375901889*eccentricity_15;
    double common_term_160 = 125725.37115503593*eccentricity_18 + 39027.468688471127*eccentricity_16 + 8770.7349815502578*eccentricity_14;
    double common_term_161 = 243219.50039692118*eccentricity_19 + 89274.024661040498*eccentricity_17 + 26601.523448693696*eccentricity_15 + 5504.6219770765671*eccentricity_13;
    double common_term_162 = 176723.74687709893*eccentricity_18 + 62743.436603254658*eccentricity_16 + 17901.450702323614*eccentricity_14 + 3421.4312731157892*eccentricity_12;
    double common_term_163 = 319602.64159973512*eccentricity_19 + 127173.22323520637*eccentricity_17 + 43611.956193055484*eccentricity_15 + 11887.858538566255*eccentricity_13 + 2103.2376849850741*eccentricity_11;
    double common_term_164 = 234307.54926008655*eccentricity_18 + 90569.30866500826*eccentricity_16 + 29952.798654308253*eccentricity_14 + 7783.9591872594998*eccentricity_12 + 1276.5633008156966*eccentricity_10;
    double common_term_165 = 403469.38357003145*eccentricity_19 + 170106.88567987157*eccentricity_17 + 63777.513981546737*eccentricity_15 + 20304.147889377378*eccentricity_13 + 5019.5518635455473*eccentricity_11 + 763.40796896852093*eccentricity_9;
    double common_term_166 = 297643.81710957402*eccentricity_18 + 122192.66030397093*eccentricity_16 + 44359.788080166903*eccentricity_14 + 13566.126449497768*eccentricity_12 + 3182.6087611607143*eccentricity_10 + 448.60078125*eccentricity_8;
    double common_term_167 = 493746.88522566438*eccentricity_19 + 217384.87587431255*eccentricity_17 + 86757.508682937579*eccentricity_15 + 30434.914156567132*eccentricity_13 + 8918.7079399182893*eccentricity_11 + 1979.6183682396298*eccentricity_9 + 258.11415705605159*eccentricity_7;
    double common_term_168 = 365797.48265010547*eccentricity_18 + 157029.09361794554*eccentricity_16 + 60807.581296233144*eccentricity_14 + 20563.361818185994*eccentricity_12 + 5756.4458023313492*eccentricity_10 + 1204.2793650793651*eccentricity_8 + 144.71805555555556*eccentricity_6;
    double common_term_169 = 589192.84976935149*eccentricity_19 + 268196.44295711927*eccentricity_17 + 112054.0201218079*eccentricity_15 + 42006.00026381697*eccentricity_13 + 13653.1116166796*eccentricity_11 + 3636.8901907784598*eccentricity_9 + 713.45830078125*eccentricity_7 + 78.53203125*eccentricity_5;
    double common_term_170 = 437684.67519404149*eccentricity_18 + 194378.6562848673*eccentricity_16 + 78872.207621392862*eccentricity_14 + 28542.739133029514*eccentricity_12 + 8883.2707671957672*eccentricity_10 + 2240.1704861111111*eccentricity_8 + 409.125*eccentricity_6 + 40.833333333333333*eccentricity_4;
    double common_term_171 = 688403.66054336277*eccentricity_19 + 321599.47985147387*eccentricity_17 + 139065.34923841712*eccentricity_15 + 54657.662972711104*eccentricity_13 + 19027.201894096092*eccentricity_11 + 5642.6422087492766*eccentricity_9 + 1337.6403971354167*eccentricity_7 + 225.02473958333333*eccentricity_5 + 20.020833333333333*eccentricity_3;
    double common_term_172 = std::pow(1.0 - eccentricity_2, -9.5)*(1.96875*eccentricity_8 + 19.6875*eccentricity_6 + 31.5*eccentricity_4 + 9.0*eccentricity_2);
    double common_term_173 = 789819.1786291366*eccentricity_19 + 376524.96773389908*eccentricity_17 + 167088.86412091259*eccentricity_15 + 67951.780440018097*eccentricity_13 + 24786.066802034731*eccentricity_11 + 7862.1723083496094*eccentricity_9 + 2069.1467556423611*eccentricity_7 + 418.09114583333333*eccentricity_5 + 56.0625*eccentricity_3 + 3.5*eccentricity;
    double common_term_174 = 587611.04435552879*eccentricity_18 + 273243.81438047871*eccentricity_16 + 117705.31147752381*eccentricity_14 + 46166.126637128665*eccentricity_12 + 16094.778923611111*eccentricity_10 + 4815.6284722222222*eccentricity_8 + 1171.1284722222222*eccentricity_6 + 211.0625*eccentricity_4 + 23.5*eccentricity_2 + 1.0;
    double common_term_175 = 891727.78670617842*eccentricity_19 + 431781.4933825241*eccentricity_17 + 195324.99331182875*eccentricity_15 + 81375.283556417738*eccentricity_13 + 30618.552248382568*eccentricity_11 + 10120.135345458984*eccentricity_9 + 2818.54541015625*eccentricity_7 + 618.2578125*eccentricity_5 + 94.3125*eccentricity_3 + 7.5*eccentricity;
    double common_term_176 = 662783.58640818212*eccentricity_18 + 312798.65258819381*eccentricity_16 + 137185.77507562762*eccentricity_14 + 55005.507094494048*eccentricity_12 + 19708.091304976852*eccentricity_10 + 6101.6444444444444*eccentricity_8 + 1547.75*eccentricity_6 + 292.16666666666667*eccentricity_4 + 33.5*eccentricity_2;
    double common_term_177 = 992271.64644624576*eccentricity_19 + 486059.99073301221*eccentricity_17 + 222881.45152326412*eccentricity_15 + 94343.864676982007*eccentricity_13 + 36160.464863953767*eccentricity_11 + 12203.147214536314*eccentricity_9 + 3471.0574544270833*eccentricity_7 + 768.07682291666667*eccentricity_5 + 115.60416666666667*eccentricity_3;
    double common_term_178 = 735966.3817457862*eccentricity_18 + 350954.07066671317*eccentricity_16 + 155714.95866420201*eccentricity_14 + 63224.796323939732*eccentricity_12 + 22936.734877232143*eccentricity_10 + 7169.06953125*eccentricity_8 + 1795.36875*eccentricity_6 + 339.9375*eccentricity_4;
    double common_term_179 = 1089452.1553185649*eccentricity_19 + 537938.7168549791*eccentricity_17 + 248777.49520401884*eccentricity_15 + 106207.14641562381*eccentricity_13 + 40991.901643187911*eccentricity_11 + 13885.421431477865*eccentricity_9 + 3829.5984483506944*eccentricity_7 + 895.01744791666667*eccentricity_5;
    double common_term_180 = 805406.87271481334*eccentricity_18 + 386474.53720636006*eccentricity_16 + 172461.03920998907*eccentricity_14 + 70237.92071851117*eccentricity_12 + 25550.523834325397*eccentricity_10 + 7566.1055555555556*eccentricity_8 + 2172.0972222222222*eccentricity_6;
    double common_term_181 = 1181137.5442669622*eccentricity_19 + 585877.49210235729*eccentricity_17 + 271998.33513745087*eccentricity_15 + 116073.89904596737*eccentricity_13 + 45131.747396033151*eccentricity_11 + 13952.05462777274*eccentricity_9 + 4950.4390764508929*eccentricity_7;
    double common_term_182 = 869173.34939077183*eccentricity_18 + 418273.58801654392*eccentricity_16 + 185762.48939678646*eccentricity_14 + 77236.197046027612*eccentricity_12 + 24050.714995315256*eccentricity_10 + 10731.751364087302*eccentricity_8;
    double common_term_183 = 1264770.5665176078*eccentricity_19 + 629270.03590815109*eccentricity_17 + 288555.94396704217*eccentricity_15 + 129249.85086463271*eccentricity_13 + 38539.346853018328*eccentricity_11 + 22332.510335668601*eccentricity_9;
    double common_term_184 = 929221.76420179613*eccentricity_18 + 435112.20387268981*eccentricity_16 + 213629.86151481331*eccentricity_14 + 56401.696065340909*eccentricity_12 + 44916.095558035714*eccentricity_10;
    double common_term_185 = 1351829.2749774022*eccentricity_19 + 634990.42284656349*eccentricity_17 + 352644.84690263587*eccentricity_15 + 71970.9917549521*eccentricity_13 + 87766.10919416494*eccentricity_11;
    double common_term_186 = 889473.81433039113*eccentricity_18 + 588194.66599826298*eccentricity_16 + 68751.752767771307*eccentricity_14 + 167295.6226284693*eccentricity_12;
    double common_term_187 = 1173734.9388192032*eccentricity_19 + 1001906.0421649719*eccentricity_17 + 7166.5252516924226*eccentricity_15 + 312097.97744729076*eccentricity_13;
    double common_term_188 = 1755510.4342203288*eccentricity_18 - 198850.21150730427*eccentricity_16 + 571341.39167059263*eccentricity_14;
    double common_term_189 = 3170355.2402208336*eccentricity_19 - 725293.97401762842*eccentricity_17 + 1028598.333538794*eccentricity_15;
    double common_term_190 = -1916885.9304023117*eccentricity_18 + 1824446.4906466377*eccentricity_16;
    double common_term_191 = -4427310.0109400648*eccentricity_19 + 3193136.3920832776*eccentricity_17;
    double common_term_192 = 5521684.4806277122*eccentricity_18;
    double common_term_193 = 9444523.2428855771*eccentricity_19;
    double common_term_194 = 1020324.9362591622*eccentricity_19;
    double common_term_195 = 629038.18850110831*eccentricity_18;
    double common_term_196 = 467115.15602483617*eccentricity_19 + 384588.48078326111*eccentricity_17;
    double common_term_197 = 365104.07238555786*eccentricity_18 + 232982.5585813749*eccentricity_16;
    double common_term_198 = 796277.62880306136*eccentricity_19 + 269740.82610943422*eccentricity_17 + 139708.68784218533*eccentricity_15;
    double common_term_199 = 550729.40311594368*eccentricity_18 + 191137.75037063149*eccentricity_16 + 82829.353837736057*eccentricity_14;
    double common_term_200 = 880585.13279601822*eccentricity_19 + 377931.40039411206*eccentricity_17 + 130923.06267031424*eccentricity_15 + 48483.878171030515*eccentricity_13;
    double common_term_201 = 622416.57077588657*eccentricity_18 + 256466.68845654833*eccentricity_16 + 87060.351743539664*eccentricity_14 + 27972.480186941964*eccentricity_12;
    double common_term_202 = 981880.19771479974*eccentricity_19 + 433645.93144465831*eccentricity_17 + 171618.37500923486*eccentricity_15 + 56324.621720736673*eccentricity_13 + 15874.286491991612*eccentricity_11;
    double common_term_203 = 695354.81434512977*eccentricity_18 + 297508.91555410983*eccentricity_16 + 112955.73072896573*eccentricity_14 + 35477.132213128307*eccentricity_12 + 8838.6522851906966*eccentricity_10;
    double common_term_204 = 1067424.5495830361*eccentricity_19 + 485487.18598420322*eccentricity_17 + 200714.73108897116*eccentricity_15 + 72943.453404656819*eccentricity_13 + 21744.861643011911*eccentricity_11 + 4813.0695709228516*eccentricity_9;
    double common_term_205 = 756471.96248446309*eccentricity_18 + 333681.45621316457*eccentricity_16 + 132926.16607212899*eccentricity_14 + 46094.235534267526*eccentricity_12 + 12948.758050870811*eccentricity_10 + 2552.854191468254*eccentricity_8;
    double common_term_206 = 1137544.589899778*eccentricity_19 + 528173.48305977826*eccentricity_17 + 225373.22710565397*eccentricity_15 + 86225.641583319809*eccentricity_13 + 28415.779286085529*eccentricity_11 + 7470.4682444254557*eccentricity_9 + 1311.7809136284722*eccentricity_7;
    double common_term_207 = 805169.44942701939*eccentricity_18 + 362704.71270117187*eccentricity_16 + 149262.66018554688*eccentricity_14 + 54632.448995535714*eccentricity_12 + 17025.143805803571*eccentricity_10 + 4157.7741071428571*eccentricity_8 + 648.278125*eccentricity_6;
    double common_term_208 = 1191100.8810914528*eccentricity_19 + 560963.72188863101*eccentricity_17 + 244468.53027659426*eccentricity_15 + 96674.068711462405*eccentricity_13 + 33690.52670468729*eccentricity_11 + 9865.7245890299479*eccentricity_9 + 2218.4079318576389*eccentricity_7 + 304.97786458333333*eccentricity_5;
    double common_term_209 = 840573.09465315517*eccentricity_18 + 383927.34040087422*eccentricity_16 + 161311.29156125418*eccentricity_14 + 61022.377556113591*eccentricity_12 + 20127.049388227513*eccentricity_10 + 5493.2090277777778*eccentricity_8 + 1124.26875*eccentricity_6 + 134.52083333333333*eccentricity_4;
    double common_term_210 = 1227265.6495161146*eccentricity_19 + 583193.0069256193*eccentricity_17 + 257488.28485027066*eccentricity_15 + 103860.76791643483*eccentricity_13 + 37371.455359431676*eccentricity_11 + 11575.071057128906*eccentricity_9 + 2911.84072265625*eccentricity_7 + 533.62890625*eccentricity_5 + 54.3125*eccentricity_3;
    double common_term_211 = 862063.36628591829*eccentricity_18 + 396857.73401790116*eccentricity_16 + 168692.96589035976*eccentricity_14 + 64971.088209325397*eccentricity_12 + 22071.572251157407*eccentricity_10 + 6351.5298611111111*eccentricity_8 + 1449.421875*eccentricity_6 + 231.91666666666667*eccentricity_4 + 19.25*eccentricity_2;
    double common_term_212 = 1245486.2783117408*eccentricity_19 + 594418.46341030252*eccentricity_17 + 264085.10115677552*eccentricity_15 + 107520.56849235895*eccentricity_13 + 39261.090514791983*eccentricity_11 + 12464.666538492839*eccentricity_9 + 3281.9704318576389*eccentricity_7 + 662.89322916666667*eccentricity_5 + 88.6875*eccentricity_3 + 5.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(39);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (10, 0).
    // q = -19
    result_by_lpq.set(c_Key3(10, 0, -19), common_term_0);
    result_by_q.set(c_Key1(-19), common_term_0);
    // q = -18
    result_by_lpq.set(c_Key3(10, 0, -18), common_term_1);
    result_by_q.set(c_Key1(-18), common_term_1);
    // q = -17
    result_by_lpq.set(c_Key3(10, 0, -17), common_term_2);
    result_by_q.set(c_Key1(-17), common_term_2);
    // q = -16
    result_by_lpq.set(c_Key3(10, 0, -16), common_term_3);
    result_by_q.set(c_Key1(-16), common_term_3);
    // q = -15
    result_by_lpq.set(c_Key3(10, 0, -15), common_term_4);
    result_by_q.set(c_Key1(-15), common_term_4);
    // q = -14
    result_by_lpq.set(c_Key3(10, 0, -14), common_term_5);
    result_by_q.set(c_Key1(-14), common_term_5);
    // q = -13
    result_by_lpq.set(c_Key3(10, 0, -13), common_term_6);
    result_by_q.set(c_Key1(-13), common_term_6);
    // q = -12
    result_by_lpq.set(c_Key3(10, 0, -12), common_term_7);
    result_by_q.set(c_Key1(-12), common_term_7);
    // q = -11
    result_by_lpq.set(c_Key3(10, 0, -11), common_term_8);
    result_by_q.set(c_Key1(-11), common_term_8);
    // q = -9
    result_by_lpq.set(c_Key3(10, 0, -9), common_term_9);
    result_by_q.set(c_Key1(-9), common_term_9);
    // q = -8
    result_by_lpq.set(c_Key3(10, 0, -8), common_term_10);
    result_by_q.set(c_Key1(-8), common_term_10);
    // q = -7
    result_by_lpq.set(c_Key3(10, 0, -7), common_term_11);
    result_by_q.set(c_Key1(-7), common_term_11);
    // q = -6
    result_by_lpq.set(c_Key3(10, 0, -6), common_term_12);
    result_by_q.set(c_Key1(-6), common_term_12);
    // q = -5
    result_by_lpq.set(c_Key3(10, 0, -5), common_term_13);
    result_by_q.set(c_Key1(-5), common_term_13);
    // q = -4
    result_by_lpq.set(c_Key3(10, 0, -4), common_term_14);
    result_by_q.set(c_Key1(-4), common_term_14);
    // q = -3
    result_by_lpq.set(c_Key3(10, 0, -3), common_term_15);
    result_by_q.set(c_Key1(-3), common_term_15);
    // q = -2
    result_by_lpq.set(c_Key3(10, 0, -2), common_term_16);
    result_by_q.set(c_Key1(-2), common_term_16);
    // q = -1
    result_by_lpq.set(c_Key3(10, 0, -1), common_term_17);
    result_by_q.set(c_Key1(-1), common_term_17);
    // q = 0
    result_by_lpq.set(c_Key3(10, 0, 0), common_term_18);
    result_by_q.set(c_Key1(0), common_term_18);
    // q = 1
    result_by_lpq.set(c_Key3(10, 0, 1), common_term_19);
    result_by_q.set(c_Key1(1), common_term_19);
    // q = 2
    result_by_lpq.set(c_Key3(10, 0, 2), common_term_20);
    result_by_q.set(c_Key1(2), common_term_20);
    // q = 3
    result_by_lpq.set(c_Key3(10, 0, 3), common_term_21);
    result_by_q.set(c_Key1(3), common_term_21);
    // q = 4
    result_by_lpq.set(c_Key3(10, 0, 4), common_term_22);
    result_by_q.set(c_Key1(4), common_term_22);
    // q = 5
    result_by_lpq.set(c_Key3(10, 0, 5), common_term_23);
    result_by_q.set(c_Key1(5), common_term_23);
    // q = 6
    result_by_lpq.set(c_Key3(10, 0, 6), common_term_24);
    result_by_q.set(c_Key1(6), common_term_24);
    // q = 7
    result_by_lpq.set(c_Key3(10, 0, 7), common_term_25);
    result_by_q.set(c_Key1(7), common_term_25);
    // q = 8
    result_by_lpq.set(c_Key3(10, 0, 8), common_term_26);
    result_by_q.set(c_Key1(8), common_term_26);
    // q = 9
    result_by_lpq.set(c_Key3(10, 0, 9), common_term_27);
    result_by_q.set(c_Key1(9), common_term_27);
    // q = 10
    result_by_lpq.set(c_Key3(10, 0, 10), common_term_28);
    result_by_q.set(c_Key1(10), common_term_28);
    // q = 11
    result_by_lpq.set(c_Key3(10, 0, 11), common_term_29);
    result_by_q.set(c_Key1(11), common_term_29);
    // q = 12
    result_by_lpq.set(c_Key3(10, 0, 12), common_term_30);
    result_by_q.set(c_Key1(12), common_term_30);
    // q = 13
    result_by_lpq.set(c_Key3(10, 0, 13), common_term_31);
    result_by_q.set(c_Key1(13), common_term_31);
    // q = 14
    result_by_lpq.set(c_Key3(10, 0, 14), common_term_32);
    result_by_q.set(c_Key1(14), common_term_32);
    // q = 15
    result_by_lpq.set(c_Key3(10, 0, 15), common_term_33);
    result_by_q.set(c_Key1(15), common_term_33);
    // q = 16
    result_by_lpq.set(c_Key3(10, 0, 16), common_term_34);
    result_by_q.set(c_Key1(16), common_term_34);
    // q = 17
    result_by_lpq.set(c_Key3(10, 0, 17), common_term_35);
    result_by_q.set(c_Key1(17), common_term_35);
    // q = 18
    result_by_lpq.set(c_Key3(10, 0, 18), common_term_36);
    result_by_q.set(c_Key1(18), common_term_36);
    // q = 19
    result_by_lpq.set(c_Key3(10, 0, 19), common_term_37);
    result_by_q.set(c_Key1(19), common_term_37);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 0), result_by_q);
    result_by_q.clear();

    // l , p = (10, 1).
    // q = -19
    result_by_lpq.set(c_Key3(10, 1, -19), common_term_38);
    result_by_q.set(c_Key1(-19), common_term_38);
    // q = -18
    result_by_lpq.set(c_Key3(10, 1, -18), common_term_39);
    result_by_q.set(c_Key1(-18), common_term_39);
    // q = -17
    result_by_lpq.set(c_Key3(10, 1, -17), common_term_40);
    result_by_q.set(c_Key1(-17), common_term_40);
    // q = -16
    result_by_lpq.set(c_Key3(10, 1, -16), common_term_41);
    result_by_q.set(c_Key1(-16), common_term_41);
    // q = -15
    result_by_lpq.set(c_Key3(10, 1, -15), common_term_42);
    result_by_q.set(c_Key1(-15), common_term_42);
    // q = -14
    result_by_lpq.set(c_Key3(10, 1, -14), common_term_43);
    result_by_q.set(c_Key1(-14), common_term_43);
    // q = -13
    result_by_lpq.set(c_Key3(10, 1, -13), common_term_44);
    result_by_q.set(c_Key1(-13), common_term_44);
    // q = -12
    result_by_lpq.set(c_Key3(10, 1, -12), common_term_45);
    result_by_q.set(c_Key1(-12), common_term_45);
    // q = -11
    result_by_lpq.set(c_Key3(10, 1, -11), common_term_46);
    result_by_q.set(c_Key1(-11), common_term_46);
    // q = -10
    result_by_lpq.set(c_Key3(10, 1, -10), common_term_47);
    result_by_q.set(c_Key1(-10), common_term_47);
    // q = -9
    result_by_lpq.set(c_Key3(10, 1, -9), common_term_48);
    result_by_q.set(c_Key1(-9), common_term_48);
    // q = -8
    result_by_lpq.set(c_Key3(10, 1, -8), common_term_49);
    result_by_q.set(c_Key1(-8), common_term_49);
    // q = -7
    result_by_lpq.set(c_Key3(10, 1, -7), common_term_50);
    result_by_q.set(c_Key1(-7), common_term_50);
    // q = -6
    result_by_lpq.set(c_Key3(10, 1, -6), common_term_51);
    result_by_q.set(c_Key1(-6), common_term_51);
    // q = -5
    result_by_lpq.set(c_Key3(10, 1, -5), common_term_52);
    result_by_q.set(c_Key1(-5), common_term_52);
    // q = -4
    result_by_lpq.set(c_Key3(10, 1, -4), common_term_53);
    result_by_q.set(c_Key1(-4), common_term_53);
    // q = -3
    result_by_lpq.set(c_Key3(10, 1, -3), common_term_54);
    result_by_q.set(c_Key1(-3), common_term_54);
    // q = -2
    result_by_lpq.set(c_Key3(10, 1, -2), common_term_55);
    result_by_q.set(c_Key1(-2), common_term_55);
    // q = -1
    result_by_lpq.set(c_Key3(10, 1, -1), common_term_56);
    result_by_q.set(c_Key1(-1), common_term_56);
    // q = 0
    result_by_lpq.set(c_Key3(10, 1, 0), common_term_57);
    result_by_q.set(c_Key1(0), common_term_57);
    // q = 1
    result_by_lpq.set(c_Key3(10, 1, 1), common_term_58);
    result_by_q.set(c_Key1(1), common_term_58);
    // q = 2
    result_by_lpq.set(c_Key3(10, 1, 2), common_term_59);
    result_by_q.set(c_Key1(2), common_term_59);
    // q = 3
    result_by_lpq.set(c_Key3(10, 1, 3), common_term_60);
    result_by_q.set(c_Key1(3), common_term_60);
    // q = 4
    result_by_lpq.set(c_Key3(10, 1, 4), common_term_61);
    result_by_q.set(c_Key1(4), common_term_61);
    // q = 5
    result_by_lpq.set(c_Key3(10, 1, 5), common_term_62);
    result_by_q.set(c_Key1(5), common_term_62);
    // q = 6
    result_by_lpq.set(c_Key3(10, 1, 6), common_term_63);
    result_by_q.set(c_Key1(6), common_term_63);
    // q = 7
    result_by_lpq.set(c_Key3(10, 1, 7), common_term_64);
    result_by_q.set(c_Key1(7), common_term_64);
    // q = 8
    result_by_lpq.set(c_Key3(10, 1, 8), common_term_65);
    result_by_q.set(c_Key1(8), common_term_65);
    // q = 9
    result_by_lpq.set(c_Key3(10, 1, 9), common_term_66);
    result_by_q.set(c_Key1(9), common_term_66);
    // q = 10
    result_by_lpq.set(c_Key3(10, 1, 10), common_term_67);
    result_by_q.set(c_Key1(10), common_term_67);
    // q = 11
    result_by_lpq.set(c_Key3(10, 1, 11), common_term_68);
    result_by_q.set(c_Key1(11), common_term_68);
    // q = 12
    result_by_lpq.set(c_Key3(10, 1, 12), common_term_69);
    result_by_q.set(c_Key1(12), common_term_69);
    // q = 13
    result_by_lpq.set(c_Key3(10, 1, 13), common_term_70);
    result_by_q.set(c_Key1(13), common_term_70);
    // q = 14
    result_by_lpq.set(c_Key3(10, 1, 14), common_term_71);
    result_by_q.set(c_Key1(14), common_term_71);
    // q = 15
    result_by_lpq.set(c_Key3(10, 1, 15), common_term_72);
    result_by_q.set(c_Key1(15), common_term_72);
    // q = 16
    result_by_lpq.set(c_Key3(10, 1, 16), common_term_73);
    result_by_q.set(c_Key1(16), common_term_73);
    // q = 17
    result_by_lpq.set(c_Key3(10, 1, 17), common_term_74);
    result_by_q.set(c_Key1(17), common_term_74);
    // q = 18
    result_by_lpq.set(c_Key3(10, 1, 18), common_term_75);
    result_by_q.set(c_Key1(18), common_term_75);
    // q = 19
    result_by_lpq.set(c_Key3(10, 1, 19), common_term_76);
    result_by_q.set(c_Key1(19), common_term_76);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 1), result_by_q);
    result_by_q.clear();

    // l , p = (10, 2).
    // q = -19
    result_by_lpq.set(c_Key3(10, 2, -19), common_term_77);
    result_by_q.set(c_Key1(-19), common_term_77);
    // q = -18
    result_by_lpq.set(c_Key3(10, 2, -18), common_term_78);
    result_by_q.set(c_Key1(-18), common_term_78);
    // q = -17
    result_by_lpq.set(c_Key3(10, 2, -17), common_term_79);
    result_by_q.set(c_Key1(-17), common_term_79);
    // q = -16
    result_by_lpq.set(c_Key3(10, 2, -16), common_term_80);
    result_by_q.set(c_Key1(-16), common_term_80);
    // q = -15
    result_by_lpq.set(c_Key3(10, 2, -15), common_term_81);
    result_by_q.set(c_Key1(-15), common_term_81);
    // q = -14
    result_by_lpq.set(c_Key3(10, 2, -14), common_term_82);
    result_by_q.set(c_Key1(-14), common_term_82);
    // q = -13
    result_by_lpq.set(c_Key3(10, 2, -13), common_term_83);
    result_by_q.set(c_Key1(-13), common_term_83);
    // q = -12
    result_by_lpq.set(c_Key3(10, 2, -12), common_term_84);
    result_by_q.set(c_Key1(-12), common_term_84);
    // q = -11
    result_by_lpq.set(c_Key3(10, 2, -11), common_term_85);
    result_by_q.set(c_Key1(-11), common_term_85);
    // q = -10
    result_by_lpq.set(c_Key3(10, 2, -10), common_term_86);
    result_by_q.set(c_Key1(-10), common_term_86);
    // q = -9
    result_by_lpq.set(c_Key3(10, 2, -9), common_term_87);
    result_by_q.set(c_Key1(-9), common_term_87);
    // q = -8
    result_by_lpq.set(c_Key3(10, 2, -8), common_term_88);
    result_by_q.set(c_Key1(-8), common_term_88);
    // q = -7
    result_by_lpq.set(c_Key3(10, 2, -7), common_term_89);
    result_by_q.set(c_Key1(-7), common_term_89);
    // q = -6
    result_by_lpq.set(c_Key3(10, 2, -6), common_term_90);
    result_by_q.set(c_Key1(-6), common_term_90);
    // q = -5
    result_by_lpq.set(c_Key3(10, 2, -5), common_term_91);
    result_by_q.set(c_Key1(-5), common_term_91);
    // q = -4
    result_by_lpq.set(c_Key3(10, 2, -4), common_term_92);
    result_by_q.set(c_Key1(-4), common_term_92);
    // q = -3
    result_by_lpq.set(c_Key3(10, 2, -3), common_term_93);
    result_by_q.set(c_Key1(-3), common_term_93);
    // q = -2
    result_by_lpq.set(c_Key3(10, 2, -2), common_term_94);
    result_by_q.set(c_Key1(-2), common_term_94);
    // q = -1
    result_by_lpq.set(c_Key3(10, 2, -1), common_term_95);
    result_by_q.set(c_Key1(-1), common_term_95);
    // q = 0
    result_by_lpq.set(c_Key3(10, 2, 0), common_term_96);
    result_by_q.set(c_Key1(0), common_term_96);
    // q = 1
    result_by_lpq.set(c_Key3(10, 2, 1), common_term_97);
    result_by_q.set(c_Key1(1), common_term_97);
    // q = 2
    result_by_lpq.set(c_Key3(10, 2, 2), common_term_98);
    result_by_q.set(c_Key1(2), common_term_98);
    // q = 3
    result_by_lpq.set(c_Key3(10, 2, 3), common_term_99);
    result_by_q.set(c_Key1(3), common_term_99);
    // q = 4
    result_by_lpq.set(c_Key3(10, 2, 4), common_term_100);
    result_by_q.set(c_Key1(4), common_term_100);
    // q = 5
    result_by_lpq.set(c_Key3(10, 2, 5), common_term_101);
    result_by_q.set(c_Key1(5), common_term_101);
    // q = 6
    result_by_lpq.set(c_Key3(10, 2, 6), common_term_102);
    result_by_q.set(c_Key1(6), common_term_102);
    // q = 7
    result_by_lpq.set(c_Key3(10, 2, 7), common_term_103);
    result_by_q.set(c_Key1(7), common_term_103);
    // q = 8
    result_by_lpq.set(c_Key3(10, 2, 8), common_term_104);
    result_by_q.set(c_Key1(8), common_term_104);
    // q = 9
    result_by_lpq.set(c_Key3(10, 2, 9), common_term_105);
    result_by_q.set(c_Key1(9), common_term_105);
    // q = 10
    result_by_lpq.set(c_Key3(10, 2, 10), common_term_106);
    result_by_q.set(c_Key1(10), common_term_106);
    // q = 11
    result_by_lpq.set(c_Key3(10, 2, 11), common_term_107);
    result_by_q.set(c_Key1(11), common_term_107);
    // q = 12
    result_by_lpq.set(c_Key3(10, 2, 12), common_term_108);
    result_by_q.set(c_Key1(12), common_term_108);
    // q = 13
    result_by_lpq.set(c_Key3(10, 2, 13), common_term_109);
    result_by_q.set(c_Key1(13), common_term_109);
    // q = 14
    result_by_lpq.set(c_Key3(10, 2, 14), common_term_110);
    result_by_q.set(c_Key1(14), common_term_110);
    // q = 15
    result_by_lpq.set(c_Key3(10, 2, 15), common_term_111);
    result_by_q.set(c_Key1(15), common_term_111);
    // q = 16
    result_by_lpq.set(c_Key3(10, 2, 16), common_term_112);
    result_by_q.set(c_Key1(16), common_term_112);
    // q = 17
    result_by_lpq.set(c_Key3(10, 2, 17), common_term_113);
    result_by_q.set(c_Key1(17), common_term_113);
    // q = 18
    result_by_lpq.set(c_Key3(10, 2, 18), common_term_114);
    result_by_q.set(c_Key1(18), common_term_114);
    // q = 19
    result_by_lpq.set(c_Key3(10, 2, 19), common_term_115);
    result_by_q.set(c_Key1(19), common_term_115);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 2), result_by_q);
    result_by_q.clear();

    // l , p = (10, 3).
    // q = -19
    result_by_lpq.set(c_Key3(10, 3, -19), common_term_116);
    result_by_q.set(c_Key1(-19), common_term_116);
    // q = -18
    result_by_lpq.set(c_Key3(10, 3, -18), common_term_117);
    result_by_q.set(c_Key1(-18), common_term_117);
    // q = -17
    result_by_lpq.set(c_Key3(10, 3, -17), common_term_118);
    result_by_q.set(c_Key1(-17), common_term_118);
    // q = -16
    result_by_lpq.set(c_Key3(10, 3, -16), common_term_119);
    result_by_q.set(c_Key1(-16), common_term_119);
    // q = -15
    result_by_lpq.set(c_Key3(10, 3, -15), common_term_120);
    result_by_q.set(c_Key1(-15), common_term_120);
    // q = -14
    result_by_lpq.set(c_Key3(10, 3, -14), common_term_121);
    result_by_q.set(c_Key1(-14), common_term_121);
    // q = -13
    result_by_lpq.set(c_Key3(10, 3, -13), common_term_122);
    result_by_q.set(c_Key1(-13), common_term_122);
    // q = -12
    result_by_lpq.set(c_Key3(10, 3, -12), common_term_123);
    result_by_q.set(c_Key1(-12), common_term_123);
    // q = -11
    result_by_lpq.set(c_Key3(10, 3, -11), common_term_124);
    result_by_q.set(c_Key1(-11), common_term_124);
    // q = -10
    result_by_lpq.set(c_Key3(10, 3, -10), common_term_125);
    result_by_q.set(c_Key1(-10), common_term_125);
    // q = -9
    result_by_lpq.set(c_Key3(10, 3, -9), common_term_126);
    result_by_q.set(c_Key1(-9), common_term_126);
    // q = -8
    result_by_lpq.set(c_Key3(10, 3, -8), common_term_127);
    result_by_q.set(c_Key1(-8), common_term_127);
    // q = -7
    result_by_lpq.set(c_Key3(10, 3, -7), common_term_128);
    result_by_q.set(c_Key1(-7), common_term_128);
    // q = -6
    result_by_lpq.set(c_Key3(10, 3, -6), common_term_129);
    result_by_q.set(c_Key1(-6), common_term_129);
    // q = -5
    result_by_lpq.set(c_Key3(10, 3, -5), common_term_130);
    result_by_q.set(c_Key1(-5), common_term_130);
    // q = -4
    result_by_lpq.set(c_Key3(10, 3, -4), common_term_131);
    result_by_q.set(c_Key1(-4), common_term_131);
    // q = -3
    result_by_lpq.set(c_Key3(10, 3, -3), common_term_132);
    result_by_q.set(c_Key1(-3), common_term_132);
    // q = -2
    result_by_lpq.set(c_Key3(10, 3, -2), common_term_133);
    result_by_q.set(c_Key1(-2), common_term_133);
    // q = -1
    result_by_lpq.set(c_Key3(10, 3, -1), common_term_134);
    result_by_q.set(c_Key1(-1), common_term_134);
    // q = 0
    result_by_lpq.set(c_Key3(10, 3, 0), common_term_135);
    result_by_q.set(c_Key1(0), common_term_135);
    // q = 1
    result_by_lpq.set(c_Key3(10, 3, 1), common_term_136);
    result_by_q.set(c_Key1(1), common_term_136);
    // q = 2
    result_by_lpq.set(c_Key3(10, 3, 2), common_term_137);
    result_by_q.set(c_Key1(2), common_term_137);
    // q = 3
    result_by_lpq.set(c_Key3(10, 3, 3), common_term_138);
    result_by_q.set(c_Key1(3), common_term_138);
    // q = 4
    result_by_lpq.set(c_Key3(10, 3, 4), common_term_139);
    result_by_q.set(c_Key1(4), common_term_139);
    // q = 5
    result_by_lpq.set(c_Key3(10, 3, 5), common_term_140);
    result_by_q.set(c_Key1(5), common_term_140);
    // q = 6
    result_by_lpq.set(c_Key3(10, 3, 6), common_term_141);
    result_by_q.set(c_Key1(6), common_term_141);
    // q = 7
    result_by_lpq.set(c_Key3(10, 3, 7), common_term_142);
    result_by_q.set(c_Key1(7), common_term_142);
    // q = 8
    result_by_lpq.set(c_Key3(10, 3, 8), common_term_143);
    result_by_q.set(c_Key1(8), common_term_143);
    // q = 9
    result_by_lpq.set(c_Key3(10, 3, 9), common_term_144);
    result_by_q.set(c_Key1(9), common_term_144);
    // q = 10
    result_by_lpq.set(c_Key3(10, 3, 10), common_term_145);
    result_by_q.set(c_Key1(10), common_term_145);
    // q = 11
    result_by_lpq.set(c_Key3(10, 3, 11), common_term_146);
    result_by_q.set(c_Key1(11), common_term_146);
    // q = 12
    result_by_lpq.set(c_Key3(10, 3, 12), common_term_147);
    result_by_q.set(c_Key1(12), common_term_147);
    // q = 13
    result_by_lpq.set(c_Key3(10, 3, 13), common_term_148);
    result_by_q.set(c_Key1(13), common_term_148);
    // q = 14
    result_by_lpq.set(c_Key3(10, 3, 14), common_term_149);
    result_by_q.set(c_Key1(14), common_term_149);
    // q = 15
    result_by_lpq.set(c_Key3(10, 3, 15), common_term_150);
    result_by_q.set(c_Key1(15), common_term_150);
    // q = 16
    result_by_lpq.set(c_Key3(10, 3, 16), common_term_151);
    result_by_q.set(c_Key1(16), common_term_151);
    // q = 17
    result_by_lpq.set(c_Key3(10, 3, 17), common_term_152);
    result_by_q.set(c_Key1(17), common_term_152);
    // q = 18
    result_by_lpq.set(c_Key3(10, 3, 18), common_term_153);
    result_by_q.set(c_Key1(18), common_term_153);
    // q = 19
    result_by_lpq.set(c_Key3(10, 3, 19), common_term_154);
    result_by_q.set(c_Key1(19), common_term_154);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 3), result_by_q);
    result_by_q.clear();

    // l , p = (10, 4).
    // q = -19
    result_by_lpq.set(c_Key3(10, 4, -19), common_term_155);
    result_by_q.set(c_Key1(-19), common_term_155);
    // q = -18
    result_by_lpq.set(c_Key3(10, 4, -18), common_term_156);
    result_by_q.set(c_Key1(-18), common_term_156);
    // q = -17
    result_by_lpq.set(c_Key3(10, 4, -17), common_term_157);
    result_by_q.set(c_Key1(-17), common_term_157);
    // q = -16
    result_by_lpq.set(c_Key3(10, 4, -16), common_term_158);
    result_by_q.set(c_Key1(-16), common_term_158);
    // q = -15
    result_by_lpq.set(c_Key3(10, 4, -15), common_term_159);
    result_by_q.set(c_Key1(-15), common_term_159);
    // q = -14
    result_by_lpq.set(c_Key3(10, 4, -14), common_term_160);
    result_by_q.set(c_Key1(-14), common_term_160);
    // q = -13
    result_by_lpq.set(c_Key3(10, 4, -13), common_term_161);
    result_by_q.set(c_Key1(-13), common_term_161);
    // q = -12
    result_by_lpq.set(c_Key3(10, 4, -12), common_term_162);
    result_by_q.set(c_Key1(-12), common_term_162);
    // q = -11
    result_by_lpq.set(c_Key3(10, 4, -11), common_term_163);
    result_by_q.set(c_Key1(-11), common_term_163);
    // q = -10
    result_by_lpq.set(c_Key3(10, 4, -10), common_term_164);
    result_by_q.set(c_Key1(-10), common_term_164);
    // q = -9
    result_by_lpq.set(c_Key3(10, 4, -9), common_term_165);
    result_by_q.set(c_Key1(-9), common_term_165);
    // q = -8
    result_by_lpq.set(c_Key3(10, 4, -8), common_term_166);
    result_by_q.set(c_Key1(-8), common_term_166);
    // q = -7
    result_by_lpq.set(c_Key3(10, 4, -7), common_term_167);
    result_by_q.set(c_Key1(-7), common_term_167);
    // q = -6
    result_by_lpq.set(c_Key3(10, 4, -6), common_term_168);
    result_by_q.set(c_Key1(-6), common_term_168);
    // q = -5
    result_by_lpq.set(c_Key3(10, 4, -5), common_term_169);
    result_by_q.set(c_Key1(-5), common_term_169);
    // q = -4
    result_by_lpq.set(c_Key3(10, 4, -4), common_term_170);
    result_by_q.set(c_Key1(-4), common_term_170);
    // q = -3
    result_by_lpq.set(c_Key3(10, 4, -3), common_term_171);
    result_by_q.set(c_Key1(-3), common_term_171);
    // q = -2
    result_by_lpq.set(c_Key3(10, 4, -2), common_term_172);
    result_by_q.set(c_Key1(-2), common_term_172);
    // q = -1
    result_by_lpq.set(c_Key3(10, 4, -1), common_term_173);
    result_by_q.set(c_Key1(-1), common_term_173);
    // q = 0
    result_by_lpq.set(c_Key3(10, 4, 0), common_term_174);
    result_by_q.set(c_Key1(0), common_term_174);
    // q = 1
    result_by_lpq.set(c_Key3(10, 4, 1), common_term_175);
    result_by_q.set(c_Key1(1), common_term_175);
    // q = 2
    result_by_lpq.set(c_Key3(10, 4, 2), common_term_176);
    result_by_q.set(c_Key1(2), common_term_176);
    // q = 3
    result_by_lpq.set(c_Key3(10, 4, 3), common_term_177);
    result_by_q.set(c_Key1(3), common_term_177);
    // q = 4
    result_by_lpq.set(c_Key3(10, 4, 4), common_term_178);
    result_by_q.set(c_Key1(4), common_term_178);
    // q = 5
    result_by_lpq.set(c_Key3(10, 4, 5), common_term_179);
    result_by_q.set(c_Key1(5), common_term_179);
    // q = 6
    result_by_lpq.set(c_Key3(10, 4, 6), common_term_180);
    result_by_q.set(c_Key1(6), common_term_180);
    // q = 7
    result_by_lpq.set(c_Key3(10, 4, 7), common_term_181);
    result_by_q.set(c_Key1(7), common_term_181);
    // q = 8
    result_by_lpq.set(c_Key3(10, 4, 8), common_term_182);
    result_by_q.set(c_Key1(8), common_term_182);
    // q = 9
    result_by_lpq.set(c_Key3(10, 4, 9), common_term_183);
    result_by_q.set(c_Key1(9), common_term_183);
    // q = 10
    result_by_lpq.set(c_Key3(10, 4, 10), common_term_184);
    result_by_q.set(c_Key1(10), common_term_184);
    // q = 11
    result_by_lpq.set(c_Key3(10, 4, 11), common_term_185);
    result_by_q.set(c_Key1(11), common_term_185);
    // q = 12
    result_by_lpq.set(c_Key3(10, 4, 12), common_term_186);
    result_by_q.set(c_Key1(12), common_term_186);
    // q = 13
    result_by_lpq.set(c_Key3(10, 4, 13), common_term_187);
    result_by_q.set(c_Key1(13), common_term_187);
    // q = 14
    result_by_lpq.set(c_Key3(10, 4, 14), common_term_188);
    result_by_q.set(c_Key1(14), common_term_188);
    // q = 15
    result_by_lpq.set(c_Key3(10, 4, 15), common_term_189);
    result_by_q.set(c_Key1(15), common_term_189);
    // q = 16
    result_by_lpq.set(c_Key3(10, 4, 16), common_term_190);
    result_by_q.set(c_Key1(16), common_term_190);
    // q = 17
    result_by_lpq.set(c_Key3(10, 4, 17), common_term_191);
    result_by_q.set(c_Key1(17), common_term_191);
    // q = 18
    result_by_lpq.set(c_Key3(10, 4, 18), common_term_192);
    result_by_q.set(c_Key1(18), common_term_192);
    // q = 19
    result_by_lpq.set(c_Key3(10, 4, 19), common_term_193);
    result_by_q.set(c_Key1(19), common_term_193);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 4), result_by_q);
    result_by_q.clear();

    // l , p = (10, 5).
    // q = -19
    result_by_lpq.set(c_Key3(10, 5, -19), common_term_194);
    result_by_q.set(c_Key1(-19), common_term_194);
    // q = -18
    result_by_lpq.set(c_Key3(10, 5, -18), common_term_195);
    result_by_q.set(c_Key1(-18), common_term_195);
    // q = -17
    result_by_lpq.set(c_Key3(10, 5, -17), common_term_196);
    result_by_q.set(c_Key1(-17), common_term_196);
    // q = -16
    result_by_lpq.set(c_Key3(10, 5, -16), common_term_197);
    result_by_q.set(c_Key1(-16), common_term_197);
    // q = -15
    result_by_lpq.set(c_Key3(10, 5, -15), common_term_198);
    result_by_q.set(c_Key1(-15), common_term_198);
    // q = -14
    result_by_lpq.set(c_Key3(10, 5, -14), common_term_199);
    result_by_q.set(c_Key1(-14), common_term_199);
    // q = -13
    result_by_lpq.set(c_Key3(10, 5, -13), common_term_200);
    result_by_q.set(c_Key1(-13), common_term_200);
    // q = -12
    result_by_lpq.set(c_Key3(10, 5, -12), common_term_201);
    result_by_q.set(c_Key1(-12), common_term_201);
    // q = -11
    result_by_lpq.set(c_Key3(10, 5, -11), common_term_202);
    result_by_q.set(c_Key1(-11), common_term_202);
    // q = -10
    result_by_lpq.set(c_Key3(10, 5, -10), common_term_203);
    result_by_q.set(c_Key1(-10), common_term_203);
    // q = -9
    result_by_lpq.set(c_Key3(10, 5, -9), common_term_204);
    result_by_q.set(c_Key1(-9), common_term_204);
    // q = -8
    result_by_lpq.set(c_Key3(10, 5, -8), common_term_205);
    result_by_q.set(c_Key1(-8), common_term_205);
    // q = -7
    result_by_lpq.set(c_Key3(10, 5, -7), common_term_206);
    result_by_q.set(c_Key1(-7), common_term_206);
    // q = -6
    result_by_lpq.set(c_Key3(10, 5, -6), common_term_207);
    result_by_q.set(c_Key1(-6), common_term_207);
    // q = -5
    result_by_lpq.set(c_Key3(10, 5, -5), common_term_208);
    result_by_q.set(c_Key1(-5), common_term_208);
    // q = -4
    result_by_lpq.set(c_Key3(10, 5, -4), common_term_209);
    result_by_q.set(c_Key1(-4), common_term_209);
    // q = -3
    result_by_lpq.set(c_Key3(10, 5, -3), common_term_210);
    result_by_q.set(c_Key1(-3), common_term_210);
    // q = -2
    result_by_lpq.set(c_Key3(10, 5, -2), common_term_211);
    result_by_q.set(c_Key1(-2), common_term_211);
    // q = -1
    result_by_lpq.set(c_Key3(10, 5, -1), common_term_212);
    result_by_q.set(c_Key1(-1), common_term_212);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -9.5)*(2.4609375*eccentricity_8 + 26.25*eccentricity_6 + 47.25*eccentricity_4 + 18.0*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(10, 5, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(10, 5, 1), common_term_212);
    result_by_q.set(c_Key1(1), common_term_212);
    // q = 2
    result_by_lpq.set(c_Key3(10, 5, 2), common_term_211);
    result_by_q.set(c_Key1(2), common_term_211);
    // q = 3
    result_by_lpq.set(c_Key3(10, 5, 3), common_term_210);
    result_by_q.set(c_Key1(3), common_term_210);
    // q = 4
    result_by_lpq.set(c_Key3(10, 5, 4), common_term_209);
    result_by_q.set(c_Key1(4), common_term_209);
    // q = 5
    result_by_lpq.set(c_Key3(10, 5, 5), common_term_208);
    result_by_q.set(c_Key1(5), common_term_208);
    // q = 6
    result_by_lpq.set(c_Key3(10, 5, 6), common_term_207);
    result_by_q.set(c_Key1(6), common_term_207);
    // q = 7
    result_by_lpq.set(c_Key3(10, 5, 7), common_term_206);
    result_by_q.set(c_Key1(7), common_term_206);
    // q = 8
    result_by_lpq.set(c_Key3(10, 5, 8), common_term_205);
    result_by_q.set(c_Key1(8), common_term_205);
    // q = 9
    result_by_lpq.set(c_Key3(10, 5, 9), common_term_204);
    result_by_q.set(c_Key1(9), common_term_204);
    // q = 10
    result_by_lpq.set(c_Key3(10, 5, 10), common_term_203);
    result_by_q.set(c_Key1(10), common_term_203);
    // q = 11
    result_by_lpq.set(c_Key3(10, 5, 11), common_term_202);
    result_by_q.set(c_Key1(11), common_term_202);
    // q = 12
    result_by_lpq.set(c_Key3(10, 5, 12), common_term_201);
    result_by_q.set(c_Key1(12), common_term_201);
    // q = 13
    result_by_lpq.set(c_Key3(10, 5, 13), common_term_200);
    result_by_q.set(c_Key1(13), common_term_200);
    // q = 14
    result_by_lpq.set(c_Key3(10, 5, 14), common_term_199);
    result_by_q.set(c_Key1(14), common_term_199);
    // q = 15
    result_by_lpq.set(c_Key3(10, 5, 15), common_term_198);
    result_by_q.set(c_Key1(15), common_term_198);
    // q = 16
    result_by_lpq.set(c_Key3(10, 5, 16), common_term_197);
    result_by_q.set(c_Key1(16), common_term_197);
    // q = 17
    result_by_lpq.set(c_Key3(10, 5, 17), common_term_196);
    result_by_q.set(c_Key1(17), common_term_196);
    // q = 18
    result_by_lpq.set(c_Key3(10, 5, 18), common_term_195);
    result_by_q.set(c_Key1(18), common_term_195);
    // q = 19
    result_by_lpq.set(c_Key3(10, 5, 19), common_term_194);
    result_by_q.set(c_Key1(19), common_term_194);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 5), result_by_q);
    result_by_q.clear();

    // l , p = (10, 6).
    // q = -19
    result_by_lpq.set(c_Key3(10, 6, -19), common_term_193);
    result_by_q.set(c_Key1(-19), common_term_193);
    // q = -18
    result_by_lpq.set(c_Key3(10, 6, -18), common_term_192);
    result_by_q.set(c_Key1(-18), common_term_192);
    // q = -17
    result_by_lpq.set(c_Key3(10, 6, -17), common_term_191);
    result_by_q.set(c_Key1(-17), common_term_191);
    // q = -16
    result_by_lpq.set(c_Key3(10, 6, -16), common_term_190);
    result_by_q.set(c_Key1(-16), common_term_190);
    // q = -15
    result_by_lpq.set(c_Key3(10, 6, -15), common_term_189);
    result_by_q.set(c_Key1(-15), common_term_189);
    // q = -14
    result_by_lpq.set(c_Key3(10, 6, -14), common_term_188);
    result_by_q.set(c_Key1(-14), common_term_188);
    // q = -13
    result_by_lpq.set(c_Key3(10, 6, -13), common_term_187);
    result_by_q.set(c_Key1(-13), common_term_187);
    // q = -12
    result_by_lpq.set(c_Key3(10, 6, -12), common_term_186);
    result_by_q.set(c_Key1(-12), common_term_186);
    // q = -11
    result_by_lpq.set(c_Key3(10, 6, -11), common_term_185);
    result_by_q.set(c_Key1(-11), common_term_185);
    // q = -10
    result_by_lpq.set(c_Key3(10, 6, -10), common_term_184);
    result_by_q.set(c_Key1(-10), common_term_184);
    // q = -9
    result_by_lpq.set(c_Key3(10, 6, -9), common_term_183);
    result_by_q.set(c_Key1(-9), common_term_183);
    // q = -8
    result_by_lpq.set(c_Key3(10, 6, -8), common_term_182);
    result_by_q.set(c_Key1(-8), common_term_182);
    // q = -7
    result_by_lpq.set(c_Key3(10, 6, -7), common_term_181);
    result_by_q.set(c_Key1(-7), common_term_181);
    // q = -6
    result_by_lpq.set(c_Key3(10, 6, -6), common_term_180);
    result_by_q.set(c_Key1(-6), common_term_180);
    // q = -5
    result_by_lpq.set(c_Key3(10, 6, -5), common_term_179);
    result_by_q.set(c_Key1(-5), common_term_179);
    // q = -4
    result_by_lpq.set(c_Key3(10, 6, -4), common_term_178);
    result_by_q.set(c_Key1(-4), common_term_178);
    // q = -3
    result_by_lpq.set(c_Key3(10, 6, -3), common_term_177);
    result_by_q.set(c_Key1(-3), common_term_177);
    // q = -2
    result_by_lpq.set(c_Key3(10, 6, -2), common_term_176);
    result_by_q.set(c_Key1(-2), common_term_176);
    // q = -1
    result_by_lpq.set(c_Key3(10, 6, -1), common_term_175);
    result_by_q.set(c_Key1(-1), common_term_175);
    // q = 0
    result_by_lpq.set(c_Key3(10, 6, 0), common_term_174);
    result_by_q.set(c_Key1(0), common_term_174);
    // q = 1
    result_by_lpq.set(c_Key3(10, 6, 1), common_term_173);
    result_by_q.set(c_Key1(1), common_term_173);
    // q = 2
    result_by_lpq.set(c_Key3(10, 6, 2), common_term_172);
    result_by_q.set(c_Key1(2), common_term_172);
    // q = 3
    result_by_lpq.set(c_Key3(10, 6, 3), common_term_171);
    result_by_q.set(c_Key1(3), common_term_171);
    // q = 4
    result_by_lpq.set(c_Key3(10, 6, 4), common_term_170);
    result_by_q.set(c_Key1(4), common_term_170);
    // q = 5
    result_by_lpq.set(c_Key3(10, 6, 5), common_term_169);
    result_by_q.set(c_Key1(5), common_term_169);
    // q = 6
    result_by_lpq.set(c_Key3(10, 6, 6), common_term_168);
    result_by_q.set(c_Key1(6), common_term_168);
    // q = 7
    result_by_lpq.set(c_Key3(10, 6, 7), common_term_167);
    result_by_q.set(c_Key1(7), common_term_167);
    // q = 8
    result_by_lpq.set(c_Key3(10, 6, 8), common_term_166);
    result_by_q.set(c_Key1(8), common_term_166);
    // q = 9
    result_by_lpq.set(c_Key3(10, 6, 9), common_term_165);
    result_by_q.set(c_Key1(9), common_term_165);
    // q = 10
    result_by_lpq.set(c_Key3(10, 6, 10), common_term_164);
    result_by_q.set(c_Key1(10), common_term_164);
    // q = 11
    result_by_lpq.set(c_Key3(10, 6, 11), common_term_163);
    result_by_q.set(c_Key1(11), common_term_163);
    // q = 12
    result_by_lpq.set(c_Key3(10, 6, 12), common_term_162);
    result_by_q.set(c_Key1(12), common_term_162);
    // q = 13
    result_by_lpq.set(c_Key3(10, 6, 13), common_term_161);
    result_by_q.set(c_Key1(13), common_term_161);
    // q = 14
    result_by_lpq.set(c_Key3(10, 6, 14), common_term_160);
    result_by_q.set(c_Key1(14), common_term_160);
    // q = 15
    result_by_lpq.set(c_Key3(10, 6, 15), common_term_159);
    result_by_q.set(c_Key1(15), common_term_159);
    // q = 16
    result_by_lpq.set(c_Key3(10, 6, 16), common_term_158);
    result_by_q.set(c_Key1(16), common_term_158);
    // q = 17
    result_by_lpq.set(c_Key3(10, 6, 17), common_term_157);
    result_by_q.set(c_Key1(17), common_term_157);
    // q = 18
    result_by_lpq.set(c_Key3(10, 6, 18), common_term_156);
    result_by_q.set(c_Key1(18), common_term_156);
    // q = 19
    result_by_lpq.set(c_Key3(10, 6, 19), common_term_155);
    result_by_q.set(c_Key1(19), common_term_155);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 6), result_by_q);
    result_by_q.clear();

    // l , p = (10, 7).
    // q = -19
    result_by_lpq.set(c_Key3(10, 7, -19), common_term_154);
    result_by_q.set(c_Key1(-19), common_term_154);
    // q = -18
    result_by_lpq.set(c_Key3(10, 7, -18), common_term_153);
    result_by_q.set(c_Key1(-18), common_term_153);
    // q = -17
    result_by_lpq.set(c_Key3(10, 7, -17), common_term_152);
    result_by_q.set(c_Key1(-17), common_term_152);
    // q = -16
    result_by_lpq.set(c_Key3(10, 7, -16), common_term_151);
    result_by_q.set(c_Key1(-16), common_term_151);
    // q = -15
    result_by_lpq.set(c_Key3(10, 7, -15), common_term_150);
    result_by_q.set(c_Key1(-15), common_term_150);
    // q = -14
    result_by_lpq.set(c_Key3(10, 7, -14), common_term_149);
    result_by_q.set(c_Key1(-14), common_term_149);
    // q = -13
    result_by_lpq.set(c_Key3(10, 7, -13), common_term_148);
    result_by_q.set(c_Key1(-13), common_term_148);
    // q = -12
    result_by_lpq.set(c_Key3(10, 7, -12), common_term_147);
    result_by_q.set(c_Key1(-12), common_term_147);
    // q = -11
    result_by_lpq.set(c_Key3(10, 7, -11), common_term_146);
    result_by_q.set(c_Key1(-11), common_term_146);
    // q = -10
    result_by_lpq.set(c_Key3(10, 7, -10), common_term_145);
    result_by_q.set(c_Key1(-10), common_term_145);
    // q = -9
    result_by_lpq.set(c_Key3(10, 7, -9), common_term_144);
    result_by_q.set(c_Key1(-9), common_term_144);
    // q = -8
    result_by_lpq.set(c_Key3(10, 7, -8), common_term_143);
    result_by_q.set(c_Key1(-8), common_term_143);
    // q = -7
    result_by_lpq.set(c_Key3(10, 7, -7), common_term_142);
    result_by_q.set(c_Key1(-7), common_term_142);
    // q = -6
    result_by_lpq.set(c_Key3(10, 7, -6), common_term_141);
    result_by_q.set(c_Key1(-6), common_term_141);
    // q = -5
    result_by_lpq.set(c_Key3(10, 7, -5), common_term_140);
    result_by_q.set(c_Key1(-5), common_term_140);
    // q = -4
    result_by_lpq.set(c_Key3(10, 7, -4), common_term_139);
    result_by_q.set(c_Key1(-4), common_term_139);
    // q = -3
    result_by_lpq.set(c_Key3(10, 7, -3), common_term_138);
    result_by_q.set(c_Key1(-3), common_term_138);
    // q = -2
    result_by_lpq.set(c_Key3(10, 7, -2), common_term_137);
    result_by_q.set(c_Key1(-2), common_term_137);
    // q = -1
    result_by_lpq.set(c_Key3(10, 7, -1), common_term_136);
    result_by_q.set(c_Key1(-1), common_term_136);
    // q = 0
    result_by_lpq.set(c_Key3(10, 7, 0), common_term_135);
    result_by_q.set(c_Key1(0), common_term_135);
    // q = 1
    result_by_lpq.set(c_Key3(10, 7, 1), common_term_134);
    result_by_q.set(c_Key1(1), common_term_134);
    // q = 2
    result_by_lpq.set(c_Key3(10, 7, 2), common_term_133);
    result_by_q.set(c_Key1(2), common_term_133);
    // q = 3
    result_by_lpq.set(c_Key3(10, 7, 3), common_term_132);
    result_by_q.set(c_Key1(3), common_term_132);
    // q = 4
    result_by_lpq.set(c_Key3(10, 7, 4), common_term_131);
    result_by_q.set(c_Key1(4), common_term_131);
    // q = 5
    result_by_lpq.set(c_Key3(10, 7, 5), common_term_130);
    result_by_q.set(c_Key1(5), common_term_130);
    // q = 6
    result_by_lpq.set(c_Key3(10, 7, 6), common_term_129);
    result_by_q.set(c_Key1(6), common_term_129);
    // q = 7
    result_by_lpq.set(c_Key3(10, 7, 7), common_term_128);
    result_by_q.set(c_Key1(7), common_term_128);
    // q = 8
    result_by_lpq.set(c_Key3(10, 7, 8), common_term_127);
    result_by_q.set(c_Key1(8), common_term_127);
    // q = 9
    result_by_lpq.set(c_Key3(10, 7, 9), common_term_126);
    result_by_q.set(c_Key1(9), common_term_126);
    // q = 10
    result_by_lpq.set(c_Key3(10, 7, 10), common_term_125);
    result_by_q.set(c_Key1(10), common_term_125);
    // q = 11
    result_by_lpq.set(c_Key3(10, 7, 11), common_term_124);
    result_by_q.set(c_Key1(11), common_term_124);
    // q = 12
    result_by_lpq.set(c_Key3(10, 7, 12), common_term_123);
    result_by_q.set(c_Key1(12), common_term_123);
    // q = 13
    result_by_lpq.set(c_Key3(10, 7, 13), common_term_122);
    result_by_q.set(c_Key1(13), common_term_122);
    // q = 14
    result_by_lpq.set(c_Key3(10, 7, 14), common_term_121);
    result_by_q.set(c_Key1(14), common_term_121);
    // q = 15
    result_by_lpq.set(c_Key3(10, 7, 15), common_term_120);
    result_by_q.set(c_Key1(15), common_term_120);
    // q = 16
    result_by_lpq.set(c_Key3(10, 7, 16), common_term_119);
    result_by_q.set(c_Key1(16), common_term_119);
    // q = 17
    result_by_lpq.set(c_Key3(10, 7, 17), common_term_118);
    result_by_q.set(c_Key1(17), common_term_118);
    // q = 18
    result_by_lpq.set(c_Key3(10, 7, 18), common_term_117);
    result_by_q.set(c_Key1(18), common_term_117);
    // q = 19
    result_by_lpq.set(c_Key3(10, 7, 19), common_term_116);
    result_by_q.set(c_Key1(19), common_term_116);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 7), result_by_q);
    result_by_q.clear();

    // l , p = (10, 8).
    // q = -19
    result_by_lpq.set(c_Key3(10, 8, -19), common_term_115);
    result_by_q.set(c_Key1(-19), common_term_115);
    // q = -18
    result_by_lpq.set(c_Key3(10, 8, -18), common_term_114);
    result_by_q.set(c_Key1(-18), common_term_114);
    // q = -17
    result_by_lpq.set(c_Key3(10, 8, -17), common_term_113);
    result_by_q.set(c_Key1(-17), common_term_113);
    // q = -16
    result_by_lpq.set(c_Key3(10, 8, -16), common_term_112);
    result_by_q.set(c_Key1(-16), common_term_112);
    // q = -15
    result_by_lpq.set(c_Key3(10, 8, -15), common_term_111);
    result_by_q.set(c_Key1(-15), common_term_111);
    // q = -14
    result_by_lpq.set(c_Key3(10, 8, -14), common_term_110);
    result_by_q.set(c_Key1(-14), common_term_110);
    // q = -13
    result_by_lpq.set(c_Key3(10, 8, -13), common_term_109);
    result_by_q.set(c_Key1(-13), common_term_109);
    // q = -12
    result_by_lpq.set(c_Key3(10, 8, -12), common_term_108);
    result_by_q.set(c_Key1(-12), common_term_108);
    // q = -11
    result_by_lpq.set(c_Key3(10, 8, -11), common_term_107);
    result_by_q.set(c_Key1(-11), common_term_107);
    // q = -10
    result_by_lpq.set(c_Key3(10, 8, -10), common_term_106);
    result_by_q.set(c_Key1(-10), common_term_106);
    // q = -9
    result_by_lpq.set(c_Key3(10, 8, -9), common_term_105);
    result_by_q.set(c_Key1(-9), common_term_105);
    // q = -8
    result_by_lpq.set(c_Key3(10, 8, -8), common_term_104);
    result_by_q.set(c_Key1(-8), common_term_104);
    // q = -7
    result_by_lpq.set(c_Key3(10, 8, -7), common_term_103);
    result_by_q.set(c_Key1(-7), common_term_103);
    // q = -6
    result_by_lpq.set(c_Key3(10, 8, -6), common_term_102);
    result_by_q.set(c_Key1(-6), common_term_102);
    // q = -5
    result_by_lpq.set(c_Key3(10, 8, -5), common_term_101);
    result_by_q.set(c_Key1(-5), common_term_101);
    // q = -4
    result_by_lpq.set(c_Key3(10, 8, -4), common_term_100);
    result_by_q.set(c_Key1(-4), common_term_100);
    // q = -3
    result_by_lpq.set(c_Key3(10, 8, -3), common_term_99);
    result_by_q.set(c_Key1(-3), common_term_99);
    // q = -2
    result_by_lpq.set(c_Key3(10, 8, -2), common_term_98);
    result_by_q.set(c_Key1(-2), common_term_98);
    // q = -1
    result_by_lpq.set(c_Key3(10, 8, -1), common_term_97);
    result_by_q.set(c_Key1(-1), common_term_97);
    // q = 0
    result_by_lpq.set(c_Key3(10, 8, 0), common_term_96);
    result_by_q.set(c_Key1(0), common_term_96);
    // q = 1
    result_by_lpq.set(c_Key3(10, 8, 1), common_term_95);
    result_by_q.set(c_Key1(1), common_term_95);
    // q = 2
    result_by_lpq.set(c_Key3(10, 8, 2), common_term_94);
    result_by_q.set(c_Key1(2), common_term_94);
    // q = 3
    result_by_lpq.set(c_Key3(10, 8, 3), common_term_93);
    result_by_q.set(c_Key1(3), common_term_93);
    // q = 4
    result_by_lpq.set(c_Key3(10, 8, 4), common_term_92);
    result_by_q.set(c_Key1(4), common_term_92);
    // q = 5
    result_by_lpq.set(c_Key3(10, 8, 5), common_term_91);
    result_by_q.set(c_Key1(5), common_term_91);
    // q = 6
    result_by_lpq.set(c_Key3(10, 8, 6), common_term_90);
    result_by_q.set(c_Key1(6), common_term_90);
    // q = 7
    result_by_lpq.set(c_Key3(10, 8, 7), common_term_89);
    result_by_q.set(c_Key1(7), common_term_89);
    // q = 8
    result_by_lpq.set(c_Key3(10, 8, 8), common_term_88);
    result_by_q.set(c_Key1(8), common_term_88);
    // q = 9
    result_by_lpq.set(c_Key3(10, 8, 9), common_term_87);
    result_by_q.set(c_Key1(9), common_term_87);
    // q = 10
    result_by_lpq.set(c_Key3(10, 8, 10), common_term_86);
    result_by_q.set(c_Key1(10), common_term_86);
    // q = 11
    result_by_lpq.set(c_Key3(10, 8, 11), common_term_85);
    result_by_q.set(c_Key1(11), common_term_85);
    // q = 12
    result_by_lpq.set(c_Key3(10, 8, 12), common_term_84);
    result_by_q.set(c_Key1(12), common_term_84);
    // q = 13
    result_by_lpq.set(c_Key3(10, 8, 13), common_term_83);
    result_by_q.set(c_Key1(13), common_term_83);
    // q = 14
    result_by_lpq.set(c_Key3(10, 8, 14), common_term_82);
    result_by_q.set(c_Key1(14), common_term_82);
    // q = 15
    result_by_lpq.set(c_Key3(10, 8, 15), common_term_81);
    result_by_q.set(c_Key1(15), common_term_81);
    // q = 16
    result_by_lpq.set(c_Key3(10, 8, 16), common_term_80);
    result_by_q.set(c_Key1(16), common_term_80);
    // q = 17
    result_by_lpq.set(c_Key3(10, 8, 17), common_term_79);
    result_by_q.set(c_Key1(17), common_term_79);
    // q = 18
    result_by_lpq.set(c_Key3(10, 8, 18), common_term_78);
    result_by_q.set(c_Key1(18), common_term_78);
    // q = 19
    result_by_lpq.set(c_Key3(10, 8, 19), common_term_77);
    result_by_q.set(c_Key1(19), common_term_77);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 8), result_by_q);
    result_by_q.clear();

    // l , p = (10, 9).
    // q = -19
    result_by_lpq.set(c_Key3(10, 9, -19), common_term_76);
    result_by_q.set(c_Key1(-19), common_term_76);
    // q = -18
    result_by_lpq.set(c_Key3(10, 9, -18), common_term_75);
    result_by_q.set(c_Key1(-18), common_term_75);
    // q = -17
    result_by_lpq.set(c_Key3(10, 9, -17), common_term_74);
    result_by_q.set(c_Key1(-17), common_term_74);
    // q = -16
    result_by_lpq.set(c_Key3(10, 9, -16), common_term_73);
    result_by_q.set(c_Key1(-16), common_term_73);
    // q = -15
    result_by_lpq.set(c_Key3(10, 9, -15), common_term_72);
    result_by_q.set(c_Key1(-15), common_term_72);
    // q = -14
    result_by_lpq.set(c_Key3(10, 9, -14), common_term_71);
    result_by_q.set(c_Key1(-14), common_term_71);
    // q = -13
    result_by_lpq.set(c_Key3(10, 9, -13), common_term_70);
    result_by_q.set(c_Key1(-13), common_term_70);
    // q = -12
    result_by_lpq.set(c_Key3(10, 9, -12), common_term_69);
    result_by_q.set(c_Key1(-12), common_term_69);
    // q = -11
    result_by_lpq.set(c_Key3(10, 9, -11), common_term_68);
    result_by_q.set(c_Key1(-11), common_term_68);
    // q = -10
    result_by_lpq.set(c_Key3(10, 9, -10), common_term_67);
    result_by_q.set(c_Key1(-10), common_term_67);
    // q = -9
    result_by_lpq.set(c_Key3(10, 9, -9), common_term_66);
    result_by_q.set(c_Key1(-9), common_term_66);
    // q = -8
    result_by_lpq.set(c_Key3(10, 9, -8), common_term_65);
    result_by_q.set(c_Key1(-8), common_term_65);
    // q = -7
    result_by_lpq.set(c_Key3(10, 9, -7), common_term_64);
    result_by_q.set(c_Key1(-7), common_term_64);
    // q = -6
    result_by_lpq.set(c_Key3(10, 9, -6), common_term_63);
    result_by_q.set(c_Key1(-6), common_term_63);
    // q = -5
    result_by_lpq.set(c_Key3(10, 9, -5), common_term_62);
    result_by_q.set(c_Key1(-5), common_term_62);
    // q = -4
    result_by_lpq.set(c_Key3(10, 9, -4), common_term_61);
    result_by_q.set(c_Key1(-4), common_term_61);
    // q = -3
    result_by_lpq.set(c_Key3(10, 9, -3), common_term_60);
    result_by_q.set(c_Key1(-3), common_term_60);
    // q = -2
    result_by_lpq.set(c_Key3(10, 9, -2), common_term_59);
    result_by_q.set(c_Key1(-2), common_term_59);
    // q = -1
    result_by_lpq.set(c_Key3(10, 9, -1), common_term_58);
    result_by_q.set(c_Key1(-1), common_term_58);
    // q = 0
    result_by_lpq.set(c_Key3(10, 9, 0), common_term_57);
    result_by_q.set(c_Key1(0), common_term_57);
    // q = 1
    result_by_lpq.set(c_Key3(10, 9, 1), common_term_56);
    result_by_q.set(c_Key1(1), common_term_56);
    // q = 2
    result_by_lpq.set(c_Key3(10, 9, 2), common_term_55);
    result_by_q.set(c_Key1(2), common_term_55);
    // q = 3
    result_by_lpq.set(c_Key3(10, 9, 3), common_term_54);
    result_by_q.set(c_Key1(3), common_term_54);
    // q = 4
    result_by_lpq.set(c_Key3(10, 9, 4), common_term_53);
    result_by_q.set(c_Key1(4), common_term_53);
    // q = 5
    result_by_lpq.set(c_Key3(10, 9, 5), common_term_52);
    result_by_q.set(c_Key1(5), common_term_52);
    // q = 6
    result_by_lpq.set(c_Key3(10, 9, 6), common_term_51);
    result_by_q.set(c_Key1(6), common_term_51);
    // q = 7
    result_by_lpq.set(c_Key3(10, 9, 7), common_term_50);
    result_by_q.set(c_Key1(7), common_term_50);
    // q = 8
    result_by_lpq.set(c_Key3(10, 9, 8), common_term_49);
    result_by_q.set(c_Key1(8), common_term_49);
    // q = 9
    result_by_lpq.set(c_Key3(10, 9, 9), common_term_48);
    result_by_q.set(c_Key1(9), common_term_48);
    // q = 10
    result_by_lpq.set(c_Key3(10, 9, 10), common_term_47);
    result_by_q.set(c_Key1(10), common_term_47);
    // q = 11
    result_by_lpq.set(c_Key3(10, 9, 11), common_term_46);
    result_by_q.set(c_Key1(11), common_term_46);
    // q = 12
    result_by_lpq.set(c_Key3(10, 9, 12), common_term_45);
    result_by_q.set(c_Key1(12), common_term_45);
    // q = 13
    result_by_lpq.set(c_Key3(10, 9, 13), common_term_44);
    result_by_q.set(c_Key1(13), common_term_44);
    // q = 14
    result_by_lpq.set(c_Key3(10, 9, 14), common_term_43);
    result_by_q.set(c_Key1(14), common_term_43);
    // q = 15
    result_by_lpq.set(c_Key3(10, 9, 15), common_term_42);
    result_by_q.set(c_Key1(15), common_term_42);
    // q = 16
    result_by_lpq.set(c_Key3(10, 9, 16), common_term_41);
    result_by_q.set(c_Key1(16), common_term_41);
    // q = 17
    result_by_lpq.set(c_Key3(10, 9, 17), common_term_40);
    result_by_q.set(c_Key1(17), common_term_40);
    // q = 18
    result_by_lpq.set(c_Key3(10, 9, 18), common_term_39);
    result_by_q.set(c_Key1(18), common_term_39);
    // q = 19
    result_by_lpq.set(c_Key3(10, 9, 19), common_term_38);
    result_by_q.set(c_Key1(19), common_term_38);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 9), result_by_q);
    result_by_q.clear();

    // l , p = (10, 10).
    // q = -19
    result_by_lpq.set(c_Key3(10, 10, -19), common_term_37);
    result_by_q.set(c_Key1(-19), common_term_37);
    // q = -18
    result_by_lpq.set(c_Key3(10, 10, -18), common_term_36);
    result_by_q.set(c_Key1(-18), common_term_36);
    // q = -17
    result_by_lpq.set(c_Key3(10, 10, -17), common_term_35);
    result_by_q.set(c_Key1(-17), common_term_35);
    // q = -16
    result_by_lpq.set(c_Key3(10, 10, -16), common_term_34);
    result_by_q.set(c_Key1(-16), common_term_34);
    // q = -15
    result_by_lpq.set(c_Key3(10, 10, -15), common_term_33);
    result_by_q.set(c_Key1(-15), common_term_33);
    // q = -14
    result_by_lpq.set(c_Key3(10, 10, -14), common_term_32);
    result_by_q.set(c_Key1(-14), common_term_32);
    // q = -13
    result_by_lpq.set(c_Key3(10, 10, -13), common_term_31);
    result_by_q.set(c_Key1(-13), common_term_31);
    // q = -12
    result_by_lpq.set(c_Key3(10, 10, -12), common_term_30);
    result_by_q.set(c_Key1(-12), common_term_30);
    // q = -11
    result_by_lpq.set(c_Key3(10, 10, -11), common_term_29);
    result_by_q.set(c_Key1(-11), common_term_29);
    // q = -10
    result_by_lpq.set(c_Key3(10, 10, -10), common_term_28);
    result_by_q.set(c_Key1(-10), common_term_28);
    // q = -9
    result_by_lpq.set(c_Key3(10, 10, -9), common_term_27);
    result_by_q.set(c_Key1(-9), common_term_27);
    // q = -8
    result_by_lpq.set(c_Key3(10, 10, -8), common_term_26);
    result_by_q.set(c_Key1(-8), common_term_26);
    // q = -7
    result_by_lpq.set(c_Key3(10, 10, -7), common_term_25);
    result_by_q.set(c_Key1(-7), common_term_25);
    // q = -6
    result_by_lpq.set(c_Key3(10, 10, -6), common_term_24);
    result_by_q.set(c_Key1(-6), common_term_24);
    // q = -5
    result_by_lpq.set(c_Key3(10, 10, -5), common_term_23);
    result_by_q.set(c_Key1(-5), common_term_23);
    // q = -4
    result_by_lpq.set(c_Key3(10, 10, -4), common_term_22);
    result_by_q.set(c_Key1(-4), common_term_22);
    // q = -3
    result_by_lpq.set(c_Key3(10, 10, -3), common_term_21);
    result_by_q.set(c_Key1(-3), common_term_21);
    // q = -2
    result_by_lpq.set(c_Key3(10, 10, -2), common_term_20);
    result_by_q.set(c_Key1(-2), common_term_20);
    // q = -1
    result_by_lpq.set(c_Key3(10, 10, -1), common_term_19);
    result_by_q.set(c_Key1(-1), common_term_19);
    // q = 0
    result_by_lpq.set(c_Key3(10, 10, 0), common_term_18);
    result_by_q.set(c_Key1(0), common_term_18);
    // q = 1
    result_by_lpq.set(c_Key3(10, 10, 1), common_term_17);
    result_by_q.set(c_Key1(1), common_term_17);
    // q = 2
    result_by_lpq.set(c_Key3(10, 10, 2), common_term_16);
    result_by_q.set(c_Key1(2), common_term_16);
    // q = 3
    result_by_lpq.set(c_Key3(10, 10, 3), common_term_15);
    result_by_q.set(c_Key1(3), common_term_15);
    // q = 4
    result_by_lpq.set(c_Key3(10, 10, 4), common_term_14);
    result_by_q.set(c_Key1(4), common_term_14);
    // q = 5
    result_by_lpq.set(c_Key3(10, 10, 5), common_term_13);
    result_by_q.set(c_Key1(5), common_term_13);
    // q = 6
    result_by_lpq.set(c_Key3(10, 10, 6), common_term_12);
    result_by_q.set(c_Key1(6), common_term_12);
    // q = 7
    result_by_lpq.set(c_Key3(10, 10, 7), common_term_11);
    result_by_q.set(c_Key1(7), common_term_11);
    // q = 8
    result_by_lpq.set(c_Key3(10, 10, 8), common_term_10);
    result_by_q.set(c_Key1(8), common_term_10);
    // q = 9
    result_by_lpq.set(c_Key3(10, 10, 9), common_term_9);
    result_by_q.set(c_Key1(9), common_term_9);
    // q = 11
    result_by_lpq.set(c_Key3(10, 10, 11), common_term_8);
    result_by_q.set(c_Key1(11), common_term_8);
    // q = 12
    result_by_lpq.set(c_Key3(10, 10, 12), common_term_7);
    result_by_q.set(c_Key1(12), common_term_7);
    // q = 13
    result_by_lpq.set(c_Key3(10, 10, 13), common_term_6);
    result_by_q.set(c_Key1(13), common_term_6);
    // q = 14
    result_by_lpq.set(c_Key3(10, 10, 14), common_term_5);
    result_by_q.set(c_Key1(14), common_term_5);
    // q = 15
    result_by_lpq.set(c_Key3(10, 10, 15), common_term_4);
    result_by_q.set(c_Key1(15), common_term_4);
    // q = 16
    result_by_lpq.set(c_Key3(10, 10, 16), common_term_3);
    result_by_q.set(c_Key1(16), common_term_3);
    // q = 17
    result_by_lpq.set(c_Key3(10, 10, 17), common_term_2);
    result_by_q.set(c_Key1(17), common_term_2);
    // q = 18
    result_by_lpq.set(c_Key3(10, 10, 18), common_term_1);
    result_by_q.set(c_Key1(18), common_term_1);
    // q = 19
    result_by_lpq.set(c_Key3(10, 10, 19), common_term_0);
    result_by_q.set(c_Key1(19), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(10, 10), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}
