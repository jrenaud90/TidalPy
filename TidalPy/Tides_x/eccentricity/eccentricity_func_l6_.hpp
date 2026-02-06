#pragma once

#include "eccentricity_common_.hpp"

EccentricityFuncOutput c_eccentricity_function_l6_e1(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 6.
    // Functions are taylor expanded and truncated at eccentricity^1.

    // Functions Calculated for l = 6.

    c_IntMap<c_Key3, double> result_by_lpq(7);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(7);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;

    c_IntMap<c_Key1, double> result_by_q(1);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (6, 0).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(6, 0, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 0), result_by_q);
    result_by_q.clear();

    // l , p = (6, 1).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(6, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 1), result_by_q);
    result_by_q.clear();

    // l , p = (6, 2).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(6, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 2), result_by_q);
    result_by_q.clear();

    // l , p = (6, 3).
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -5.5);
    result_by_lpq.set(c_Key3(6, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 3), result_by_q);
    result_by_q.clear();

    // l , p = (6, 4).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(6, 4, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 4), result_by_q);
    result_by_q.clear();

    // l , p = (6, 5).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(6, 5, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 5), result_by_q);
    result_by_q.clear();

    // l , p = (6, 6).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(6, 6, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 6), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l6_e2(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 6.
    // Functions are taylor expanded and truncated at eccentricity^2.

    // Functions Calculated for l = 6.

    c_IntMap<c_Key3, double> result_by_lpq(23);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(7);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = -2.5*eccentricity;
    double common_term_1 = 9.5*eccentricity;
    double common_term_2 = -0.5*eccentricity;
    double common_term_3 = 7.5*eccentricity;
    double common_term_4 = 2.5*eccentricity_2*std::pow(1.0 - eccentricity_2, -5.5);
    double common_term_5 = 1.5*eccentricity;
    double common_term_6 = 5.5*eccentricity;
    double common_term_7 = 3.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(4);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (6, 0).
    // q = -1
    result_by_lpq.set(c_Key3(6, 0, -1), common_term_0);
    result_by_q.set(c_Key1(-1), common_term_0);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(6, 0, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(6, 0, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 0), result_by_q);
    result_by_q.clear();

    // l , p = (6, 1).
    // q = -1
    result_by_lpq.set(c_Key3(6, 1, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(6, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(6, 1, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 1), result_by_q);
    result_by_q.clear();

    // l , p = (6, 2).
    // q = -2
    result_by_lpq.set(c_Key3(6, 2, -2), common_term_4);
    result_by_q.set(c_Key1(-2), common_term_4);
    // q = -1
    result_by_lpq.set(c_Key3(6, 2, -1), common_term_5);
    result_by_q.set(c_Key1(-1), common_term_5);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(6, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(6, 2, 1), common_term_6);
    result_by_q.set(c_Key1(1), common_term_6);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 2), result_by_q);
    result_by_q.clear();

    // l , p = (6, 3).
    // q = -1
    result_by_lpq.set(c_Key3(6, 3, -1), common_term_7);
    result_by_q.set(c_Key1(-1), common_term_7);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -5.5)*(5.0*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(6, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(6, 3, 1), common_term_7);
    result_by_q.set(c_Key1(1), common_term_7);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 3), result_by_q);
    result_by_q.clear();

    // l , p = (6, 4).
    // q = -1
    result_by_lpq.set(c_Key3(6, 4, -1), common_term_6);
    result_by_q.set(c_Key1(-1), common_term_6);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(6, 4, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(6, 4, 1), common_term_5);
    result_by_q.set(c_Key1(1), common_term_5);
    // q = 2
    result_by_lpq.set(c_Key3(6, 4, 2), common_term_4);
    result_by_q.set(c_Key1(2), common_term_4);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 4), result_by_q);
    result_by_q.clear();

    // l , p = (6, 5).
    // q = -1
    result_by_lpq.set(c_Key3(6, 5, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(6, 5, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(6, 5, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 5), result_by_q);
    result_by_q.clear();

    // l , p = (6, 6).
    // q = -1
    result_by_lpq.set(c_Key3(6, 6, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(6, 6, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(6, 6, 1), common_term_0);
    result_by_q.set(c_Key1(1), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 6), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l6_e3(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 6.
    // Functions are taylor expanded and truncated at eccentricity^3.

    // Functions Calculated for l = 6.

    c_IntMap<c_Key3, double> result_by_lpq(35);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(7);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = 2.0*eccentricity_2;
    double common_term_1 = -2.5*eccentricity;
    double common_term_2 = 1.0 - 25.5*eccentricity_2;
    double common_term_3 = 9.5*eccentricity;
    double common_term_4 = 51.5*eccentricity_2;
    double common_term_5 = 0.25*eccentricity_2;
    double common_term_6 = -0.5*eccentricity;
    double common_term_7 = 1.0 - 5.5*eccentricity_2;
    double common_term_8 = 7.5*eccentricity;
    double common_term_9 = 33.25*eccentricity_2;
    double common_term_10 = 2.5*eccentricity_2*std::pow(1.0 - eccentricity_2, -5.5);
    double common_term_11 = 1.5*eccentricity;
    double common_term_12 = 6.5*eccentricity_2 + 1.0;
    double common_term_13 = 5.5*eccentricity;
    double common_term_14 = 19.0*eccentricity_2;
    double common_term_15 = 8.75*eccentricity_2;
    double common_term_16 = 3.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(5);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (6, 0).
    // q = -2
    result_by_lpq.set(c_Key3(6, 0, -2), common_term_0);
    result_by_q.set(c_Key1(-2), common_term_0);
    // q = -1
    result_by_lpq.set(c_Key3(6, 0, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    result_by_lpq.set(c_Key3(6, 0, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(6, 0, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // q = 2
    result_by_lpq.set(c_Key3(6, 0, 2), common_term_4);
    result_by_q.set(c_Key1(2), common_term_4);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 0), result_by_q);
    result_by_q.clear();

    // l , p = (6, 1).
    // q = -2
    result_by_lpq.set(c_Key3(6, 1, -2), common_term_5);
    result_by_q.set(c_Key1(-2), common_term_5);
    // q = -1
    result_by_lpq.set(c_Key3(6, 1, -1), common_term_6);
    result_by_q.set(c_Key1(-1), common_term_6);
    // q = 0
    result_by_lpq.set(c_Key3(6, 1, 0), common_term_7);
    result_by_q.set(c_Key1(0), common_term_7);
    // q = 1
    result_by_lpq.set(c_Key3(6, 1, 1), common_term_8);
    result_by_q.set(c_Key1(1), common_term_8);
    // q = 2
    result_by_lpq.set(c_Key3(6, 1, 2), common_term_9);
    result_by_q.set(c_Key1(2), common_term_9);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 1), result_by_q);
    result_by_q.clear();

    // l , p = (6, 2).
    // q = -2
    result_by_lpq.set(c_Key3(6, 2, -2), common_term_10);
    result_by_q.set(c_Key1(-2), common_term_10);
    // q = -1
    result_by_lpq.set(c_Key3(6, 2, -1), common_term_11);
    result_by_q.set(c_Key1(-1), common_term_11);
    // q = 0
    result_by_lpq.set(c_Key3(6, 2, 0), common_term_12);
    result_by_q.set(c_Key1(0), common_term_12);
    // q = 1
    result_by_lpq.set(c_Key3(6, 2, 1), common_term_13);
    result_by_q.set(c_Key1(1), common_term_13);
    // q = 2
    result_by_lpq.set(c_Key3(6, 2, 2), common_term_14);
    result_by_q.set(c_Key1(2), common_term_14);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 2), result_by_q);
    result_by_q.clear();

    // l , p = (6, 3).
    // q = -2
    result_by_lpq.set(c_Key3(6, 3, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(6, 3, -1), common_term_16);
    result_by_q.set(c_Key1(-1), common_term_16);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -5.5)*(5.0*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(6, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(6, 3, 1), common_term_16);
    result_by_q.set(c_Key1(1), common_term_16);
    // q = 2
    result_by_lpq.set(c_Key3(6, 3, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 3), result_by_q);
    result_by_q.clear();

    // l , p = (6, 4).
    // q = -2
    result_by_lpq.set(c_Key3(6, 4, -2), common_term_14);
    result_by_q.set(c_Key1(-2), common_term_14);
    // q = -1
    result_by_lpq.set(c_Key3(6, 4, -1), common_term_13);
    result_by_q.set(c_Key1(-1), common_term_13);
    // q = 0
    result_by_lpq.set(c_Key3(6, 4, 0), common_term_12);
    result_by_q.set(c_Key1(0), common_term_12);
    // q = 1
    result_by_lpq.set(c_Key3(6, 4, 1), common_term_11);
    result_by_q.set(c_Key1(1), common_term_11);
    // q = 2
    result_by_lpq.set(c_Key3(6, 4, 2), common_term_10);
    result_by_q.set(c_Key1(2), common_term_10);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 4), result_by_q);
    result_by_q.clear();

    // l , p = (6, 5).
    // q = -2
    result_by_lpq.set(c_Key3(6, 5, -2), common_term_9);
    result_by_q.set(c_Key1(-2), common_term_9);
    // q = -1
    result_by_lpq.set(c_Key3(6, 5, -1), common_term_8);
    result_by_q.set(c_Key1(-1), common_term_8);
    // q = 0
    result_by_lpq.set(c_Key3(6, 5, 0), common_term_7);
    result_by_q.set(c_Key1(0), common_term_7);
    // q = 1
    result_by_lpq.set(c_Key3(6, 5, 1), common_term_6);
    result_by_q.set(c_Key1(1), common_term_6);
    // q = 2
    result_by_lpq.set(c_Key3(6, 5, 2), common_term_5);
    result_by_q.set(c_Key1(2), common_term_5);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 5), result_by_q);
    result_by_q.clear();

    // l , p = (6, 6).
    // q = -2
    result_by_lpq.set(c_Key3(6, 6, -2), common_term_4);
    result_by_q.set(c_Key1(-2), common_term_4);
    // q = -1
    result_by_lpq.set(c_Key3(6, 6, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    result_by_lpq.set(c_Key3(6, 6, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(6, 6, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // q = 2
    result_by_lpq.set(c_Key3(6, 6, 2), common_term_0);
    result_by_q.set(c_Key1(2), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 6), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l6_e4(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 6.
    // Functions are taylor expanded and truncated at eccentricity^4.

    // Functions Calculated for l = 6.

    c_IntMap<c_Key3, double> result_by_lpq(51);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(7);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double common_term_0 = -0.5625*eccentricity_3;
    double common_term_1 = 2.0*eccentricity_2;
    double common_term_2 = 22.8125*eccentricity_3 - 2.5*eccentricity;
    double common_term_3 = 1.0 - 25.5*eccentricity_2;
    double common_term_4 = -147.4375*eccentricity_3 + 9.5*eccentricity;
    double common_term_5 = 51.5*eccentricity_2;
    double common_term_6 = 209.1875*eccentricity_3;
    double common_term_7 = 0.3125*eccentricity_4*std::pow(1.0 - eccentricity_2, -5.5);
    double common_term_8 = 0.22916666666666667*eccentricity_3;
    double common_term_9 = 0.25*eccentricity_2;
    double common_term_10 = 1.9375*eccentricity_3 - 0.5*eccentricity;
    double common_term_11 = 1.0 - 5.5*eccentricity_2;
    double common_term_12 = -31.5625*eccentricity_3 + 7.5*eccentricity;
    double common_term_13 = 33.25*eccentricity_2;
    double common_term_14 = 113.39583333333333*eccentricity_3;
    double common_term_15 = 4.0208333333333333*eccentricity_3;
    double common_term_16 = std::pow(1.0 - eccentricity_2, -5.5)*(1.25*eccentricity_4 + 2.5*eccentricity_2);
    double common_term_17 = 10.0625*eccentricity_3 + 1.5*eccentricity;
    double common_term_18 = 6.5*eccentricity_2 + 1.0;
    double common_term_19 = 17.3125*eccentricity_3 + 5.5*eccentricity;
    double common_term_20 = 19.0*eccentricity_2;
    double common_term_21 = 52.604166666666667*eccentricity_3;
    double common_term_22 = 18.8125*eccentricity_3;
    double common_term_23 = 8.75*eccentricity_2;
    double common_term_24 = 23.1875*eccentricity_3 + 3.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(8);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (6, 0).
    // q = -3
    result_by_lpq.set(c_Key3(6, 0, -3), common_term_0);
    result_by_q.set(c_Key1(-3), common_term_0);
    // q = -2
    result_by_lpq.set(c_Key3(6, 0, -2), common_term_1);
    result_by_q.set(c_Key1(-2), common_term_1);
    // q = -1
    result_by_lpq.set(c_Key3(6, 0, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    result_by_lpq.set(c_Key3(6, 0, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(6, 0, 1), common_term_4);
    result_by_q.set(c_Key1(1), common_term_4);
    // q = 2
    result_by_lpq.set(c_Key3(6, 0, 2), common_term_5);
    result_by_q.set(c_Key1(2), common_term_5);
    // q = 3
    result_by_lpq.set(c_Key3(6, 0, 3), common_term_6);
    result_by_q.set(c_Key1(3), common_term_6);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 0), result_by_q);
    result_by_q.clear();

    // l , p = (6, 1).
    // q = -4
    result_by_lpq.set(c_Key3(6, 1, -4), common_term_7);
    result_by_q.set(c_Key1(-4), common_term_7);
    // q = -3
    result_by_lpq.set(c_Key3(6, 1, -3), common_term_8);
    result_by_q.set(c_Key1(-3), common_term_8);
    // q = -2
    result_by_lpq.set(c_Key3(6, 1, -2), common_term_9);
    result_by_q.set(c_Key1(-2), common_term_9);
    // q = -1
    result_by_lpq.set(c_Key3(6, 1, -1), common_term_10);
    result_by_q.set(c_Key1(-1), common_term_10);
    // q = 0
    result_by_lpq.set(c_Key3(6, 1, 0), common_term_11);
    result_by_q.set(c_Key1(0), common_term_11);
    // q = 1
    result_by_lpq.set(c_Key3(6, 1, 1), common_term_12);
    result_by_q.set(c_Key1(1), common_term_12);
    // q = 2
    result_by_lpq.set(c_Key3(6, 1, 2), common_term_13);
    result_by_q.set(c_Key1(2), common_term_13);
    // q = 3
    result_by_lpq.set(c_Key3(6, 1, 3), common_term_14);
    result_by_q.set(c_Key1(3), common_term_14);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 1), result_by_q);
    result_by_q.clear();

    // l , p = (6, 2).
    // q = -3
    result_by_lpq.set(c_Key3(6, 2, -3), common_term_15);
    result_by_q.set(c_Key1(-3), common_term_15);
    // q = -2
    result_by_lpq.set(c_Key3(6, 2, -2), common_term_16);
    result_by_q.set(c_Key1(-2), common_term_16);
    // q = -1
    result_by_lpq.set(c_Key3(6, 2, -1), common_term_17);
    result_by_q.set(c_Key1(-1), common_term_17);
    // q = 0
    result_by_lpq.set(c_Key3(6, 2, 0), common_term_18);
    result_by_q.set(c_Key1(0), common_term_18);
    // q = 1
    result_by_lpq.set(c_Key3(6, 2, 1), common_term_19);
    result_by_q.set(c_Key1(1), common_term_19);
    // q = 2
    result_by_lpq.set(c_Key3(6, 2, 2), common_term_20);
    result_by_q.set(c_Key1(2), common_term_20);
    // q = 3
    result_by_lpq.set(c_Key3(6, 2, 3), common_term_21);
    result_by_q.set(c_Key1(3), common_term_21);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 2), result_by_q);
    result_by_q.clear();

    // l , p = (6, 3).
    // q = -3
    result_by_lpq.set(c_Key3(6, 3, -3), common_term_22);
    result_by_q.set(c_Key1(-3), common_term_22);
    // q = -2
    result_by_lpq.set(c_Key3(6, 3, -2), common_term_23);
    result_by_q.set(c_Key1(-2), common_term_23);
    // q = -1
    result_by_lpq.set(c_Key3(6, 3, -1), common_term_24);
    result_by_q.set(c_Key1(-1), common_term_24);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -5.5)*(1.875*eccentricity_4 + 5.0*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(6, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(6, 3, 1), common_term_24);
    result_by_q.set(c_Key1(1), common_term_24);
    // q = 2
    result_by_lpq.set(c_Key3(6, 3, 2), common_term_23);
    result_by_q.set(c_Key1(2), common_term_23);
    // q = 3
    result_by_lpq.set(c_Key3(6, 3, 3), common_term_22);
    result_by_q.set(c_Key1(3), common_term_22);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 3), result_by_q);
    result_by_q.clear();

    // l , p = (6, 4).
    // q = -3
    result_by_lpq.set(c_Key3(6, 4, -3), common_term_21);
    result_by_q.set(c_Key1(-3), common_term_21);
    // q = -2
    result_by_lpq.set(c_Key3(6, 4, -2), common_term_20);
    result_by_q.set(c_Key1(-2), common_term_20);
    // q = -1
    result_by_lpq.set(c_Key3(6, 4, -1), common_term_19);
    result_by_q.set(c_Key1(-1), common_term_19);
    // q = 0
    result_by_lpq.set(c_Key3(6, 4, 0), common_term_18);
    result_by_q.set(c_Key1(0), common_term_18);
    // q = 1
    result_by_lpq.set(c_Key3(6, 4, 1), common_term_17);
    result_by_q.set(c_Key1(1), common_term_17);
    // q = 2
    result_by_lpq.set(c_Key3(6, 4, 2), common_term_16);
    result_by_q.set(c_Key1(2), common_term_16);
    // q = 3
    result_by_lpq.set(c_Key3(6, 4, 3), common_term_15);
    result_by_q.set(c_Key1(3), common_term_15);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 4), result_by_q);
    result_by_q.clear();

    // l , p = (6, 5).
    // q = -3
    result_by_lpq.set(c_Key3(6, 5, -3), common_term_14);
    result_by_q.set(c_Key1(-3), common_term_14);
    // q = -2
    result_by_lpq.set(c_Key3(6, 5, -2), common_term_13);
    result_by_q.set(c_Key1(-2), common_term_13);
    // q = -1
    result_by_lpq.set(c_Key3(6, 5, -1), common_term_12);
    result_by_q.set(c_Key1(-1), common_term_12);
    // q = 0
    result_by_lpq.set(c_Key3(6, 5, 0), common_term_11);
    result_by_q.set(c_Key1(0), common_term_11);
    // q = 1
    result_by_lpq.set(c_Key3(6, 5, 1), common_term_10);
    result_by_q.set(c_Key1(1), common_term_10);
    // q = 2
    result_by_lpq.set(c_Key3(6, 5, 2), common_term_9);
    result_by_q.set(c_Key1(2), common_term_9);
    // q = 3
    result_by_lpq.set(c_Key3(6, 5, 3), common_term_8);
    result_by_q.set(c_Key1(3), common_term_8);
    // q = 4
    result_by_lpq.set(c_Key3(6, 5, 4), common_term_7);
    result_by_q.set(c_Key1(4), common_term_7);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 5), result_by_q);
    result_by_q.clear();

    // l , p = (6, 6).
    // q = -3
    result_by_lpq.set(c_Key3(6, 6, -3), common_term_6);
    result_by_q.set(c_Key1(-3), common_term_6);
    // q = -2
    result_by_lpq.set(c_Key3(6, 6, -2), common_term_5);
    result_by_q.set(c_Key1(-2), common_term_5);
    // q = -1
    result_by_lpq.set(c_Key3(6, 6, -1), common_term_4);
    result_by_q.set(c_Key1(-1), common_term_4);
    // q = 0
    result_by_lpq.set(c_Key3(6, 6, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(6, 6, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // q = 2
    result_by_lpq.set(c_Key3(6, 6, 2), common_term_1);
    result_by_q.set(c_Key1(2), common_term_1);
    // q = 3
    result_by_lpq.set(c_Key3(6, 6, 3), common_term_0);
    result_by_q.set(c_Key1(3), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 6), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l6_e5(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 6.
    // Functions are taylor expanded and truncated at eccentricity^5.

    // Functions Calculated for l = 6.

    c_IntMap<c_Key3, double> result_by_lpq(63);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(7);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double common_term_0 = 0.041666666666666667*eccentricity_4;
    double common_term_1 = -0.5625*eccentricity_3;
    double common_term_2 = -7.6666666666666667*eccentricity_4 + 2.0*eccentricity_2;
    double common_term_3 = 22.8125*eccentricity_3 - 2.5*eccentricity;
    double common_term_4 = 145.6875*eccentricity_4 - 25.5*eccentricity_2 + 1.0;
    double common_term_5 = -147.4375*eccentricity_3 + 9.5*eccentricity;
    double common_term_6 = -635.66666666666667*eccentricity_4 + 51.5*eccentricity_2;
    double common_term_7 = 209.1875*eccentricity_3;
    double common_term_8 = 707.60416666666667*eccentricity_4;
    double common_term_9 = 0.3125*eccentricity_4*std::pow(1.0 - eccentricity_2, -5.5);
    double common_term_10 = 0.22916666666666667*eccentricity_3;
    double common_term_11 = eccentricity_4 + 0.25*eccentricity_2;
    double common_term_12 = 1.9375*eccentricity_3 - 0.5*eccentricity;
    double common_term_13 = 11.0*eccentricity_4 - 5.5*eccentricity_2 + 1.0;
    double common_term_14 = -31.5625*eccentricity_3 + 7.5*eccentricity;
    double common_term_15 = -131.0*eccentricity_4 + 33.25*eccentricity_2;
    double common_term_16 = 113.39583333333333*eccentricity_3;
    double common_term_17 = 328.6875*eccentricity_4;
    double common_term_18 = 6.3125*eccentricity_4;
    double common_term_19 = 4.0208333333333333*eccentricity_3;
    double common_term_20 = std::pow(1.0 - eccentricity_2, -5.5)*(1.25*eccentricity_4 + 2.5*eccentricity_2);
    double common_term_21 = 10.0625*eccentricity_3 + 1.5*eccentricity;
    double common_term_22 = 26.1875*eccentricity_4 + 6.5*eccentricity_2 + 1.0;
    double common_term_23 = 17.3125*eccentricity_3 + 5.5*eccentricity;
    double common_term_24 = 35.0*eccentricity_4 + 19.0*eccentricity_2;
    double common_term_25 = 52.604166666666667*eccentricity_3;
    double common_term_26 = 127.5*eccentricity_4;
    double common_term_27 = 37.041666666666667*eccentricity_4;
    double common_term_28 = 18.8125*eccentricity_3;
    double common_term_29 = 44.333333333333333*eccentricity_4 + 8.75*eccentricity_2;
    double common_term_30 = 23.1875*eccentricity_3 + 3.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(9);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (6, 0).
    // q = -4
    result_by_lpq.set(c_Key3(6, 0, -4), common_term_0);
    result_by_q.set(c_Key1(-4), common_term_0);
    // q = -3
    result_by_lpq.set(c_Key3(6, 0, -3), common_term_1);
    result_by_q.set(c_Key1(-3), common_term_1);
    // q = -2
    result_by_lpq.set(c_Key3(6, 0, -2), common_term_2);
    result_by_q.set(c_Key1(-2), common_term_2);
    // q = -1
    result_by_lpq.set(c_Key3(6, 0, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    result_by_lpq.set(c_Key3(6, 0, 0), common_term_4);
    result_by_q.set(c_Key1(0), common_term_4);
    // q = 1
    result_by_lpq.set(c_Key3(6, 0, 1), common_term_5);
    result_by_q.set(c_Key1(1), common_term_5);
    // q = 2
    result_by_lpq.set(c_Key3(6, 0, 2), common_term_6);
    result_by_q.set(c_Key1(2), common_term_6);
    // q = 3
    result_by_lpq.set(c_Key3(6, 0, 3), common_term_7);
    result_by_q.set(c_Key1(3), common_term_7);
    // q = 4
    result_by_lpq.set(c_Key3(6, 0, 4), common_term_8);
    result_by_q.set(c_Key1(4), common_term_8);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 0), result_by_q);
    result_by_q.clear();

    // l , p = (6, 1).
    // q = -4
    result_by_lpq.set(c_Key3(6, 1, -4), common_term_9);
    result_by_q.set(c_Key1(-4), common_term_9);
    // q = -3
    result_by_lpq.set(c_Key3(6, 1, -3), common_term_10);
    result_by_q.set(c_Key1(-3), common_term_10);
    // q = -2
    result_by_lpq.set(c_Key3(6, 1, -2), common_term_11);
    result_by_q.set(c_Key1(-2), common_term_11);
    // q = -1
    result_by_lpq.set(c_Key3(6, 1, -1), common_term_12);
    result_by_q.set(c_Key1(-1), common_term_12);
    // q = 0
    result_by_lpq.set(c_Key3(6, 1, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(6, 1, 1), common_term_14);
    result_by_q.set(c_Key1(1), common_term_14);
    // q = 2
    result_by_lpq.set(c_Key3(6, 1, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // q = 3
    result_by_lpq.set(c_Key3(6, 1, 3), common_term_16);
    result_by_q.set(c_Key1(3), common_term_16);
    // q = 4
    result_by_lpq.set(c_Key3(6, 1, 4), common_term_17);
    result_by_q.set(c_Key1(4), common_term_17);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 1), result_by_q);
    result_by_q.clear();

    // l , p = (6, 2).
    // q = -4
    result_by_lpq.set(c_Key3(6, 2, -4), common_term_18);
    result_by_q.set(c_Key1(-4), common_term_18);
    // q = -3
    result_by_lpq.set(c_Key3(6, 2, -3), common_term_19);
    result_by_q.set(c_Key1(-3), common_term_19);
    // q = -2
    result_by_lpq.set(c_Key3(6, 2, -2), common_term_20);
    result_by_q.set(c_Key1(-2), common_term_20);
    // q = -1
    result_by_lpq.set(c_Key3(6, 2, -1), common_term_21);
    result_by_q.set(c_Key1(-1), common_term_21);
    // q = 0
    result_by_lpq.set(c_Key3(6, 2, 0), common_term_22);
    result_by_q.set(c_Key1(0), common_term_22);
    // q = 1
    result_by_lpq.set(c_Key3(6, 2, 1), common_term_23);
    result_by_q.set(c_Key1(1), common_term_23);
    // q = 2
    result_by_lpq.set(c_Key3(6, 2, 2), common_term_24);
    result_by_q.set(c_Key1(2), common_term_24);
    // q = 3
    result_by_lpq.set(c_Key3(6, 2, 3), common_term_25);
    result_by_q.set(c_Key1(3), common_term_25);
    // q = 4
    result_by_lpq.set(c_Key3(6, 2, 4), common_term_26);
    result_by_q.set(c_Key1(4), common_term_26);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 2), result_by_q);
    result_by_q.clear();

    // l , p = (6, 3).
    // q = -4
    result_by_lpq.set(c_Key3(6, 3, -4), common_term_27);
    result_by_q.set(c_Key1(-4), common_term_27);
    // q = -3
    result_by_lpq.set(c_Key3(6, 3, -3), common_term_28);
    result_by_q.set(c_Key1(-3), common_term_28);
    // q = -2
    result_by_lpq.set(c_Key3(6, 3, -2), common_term_29);
    result_by_q.set(c_Key1(-2), common_term_29);
    // q = -1
    result_by_lpq.set(c_Key3(6, 3, -1), common_term_30);
    result_by_q.set(c_Key1(-1), common_term_30);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -5.5)*(1.875*eccentricity_4 + 5.0*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(6, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(6, 3, 1), common_term_30);
    result_by_q.set(c_Key1(1), common_term_30);
    // q = 2
    result_by_lpq.set(c_Key3(6, 3, 2), common_term_29);
    result_by_q.set(c_Key1(2), common_term_29);
    // q = 3
    result_by_lpq.set(c_Key3(6, 3, 3), common_term_28);
    result_by_q.set(c_Key1(3), common_term_28);
    // q = 4
    result_by_lpq.set(c_Key3(6, 3, 4), common_term_27);
    result_by_q.set(c_Key1(4), common_term_27);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 3), result_by_q);
    result_by_q.clear();

    // l , p = (6, 4).
    // q = -4
    result_by_lpq.set(c_Key3(6, 4, -4), common_term_26);
    result_by_q.set(c_Key1(-4), common_term_26);
    // q = -3
    result_by_lpq.set(c_Key3(6, 4, -3), common_term_25);
    result_by_q.set(c_Key1(-3), common_term_25);
    // q = -2
    result_by_lpq.set(c_Key3(6, 4, -2), common_term_24);
    result_by_q.set(c_Key1(-2), common_term_24);
    // q = -1
    result_by_lpq.set(c_Key3(6, 4, -1), common_term_23);
    result_by_q.set(c_Key1(-1), common_term_23);
    // q = 0
    result_by_lpq.set(c_Key3(6, 4, 0), common_term_22);
    result_by_q.set(c_Key1(0), common_term_22);
    // q = 1
    result_by_lpq.set(c_Key3(6, 4, 1), common_term_21);
    result_by_q.set(c_Key1(1), common_term_21);
    // q = 2
    result_by_lpq.set(c_Key3(6, 4, 2), common_term_20);
    result_by_q.set(c_Key1(2), common_term_20);
    // q = 3
    result_by_lpq.set(c_Key3(6, 4, 3), common_term_19);
    result_by_q.set(c_Key1(3), common_term_19);
    // q = 4
    result_by_lpq.set(c_Key3(6, 4, 4), common_term_18);
    result_by_q.set(c_Key1(4), common_term_18);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 4), result_by_q);
    result_by_q.clear();

    // l , p = (6, 5).
    // q = -4
    result_by_lpq.set(c_Key3(6, 5, -4), common_term_17);
    result_by_q.set(c_Key1(-4), common_term_17);
    // q = -3
    result_by_lpq.set(c_Key3(6, 5, -3), common_term_16);
    result_by_q.set(c_Key1(-3), common_term_16);
    // q = -2
    result_by_lpq.set(c_Key3(6, 5, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(6, 5, -1), common_term_14);
    result_by_q.set(c_Key1(-1), common_term_14);
    // q = 0
    result_by_lpq.set(c_Key3(6, 5, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(6, 5, 1), common_term_12);
    result_by_q.set(c_Key1(1), common_term_12);
    // q = 2
    result_by_lpq.set(c_Key3(6, 5, 2), common_term_11);
    result_by_q.set(c_Key1(2), common_term_11);
    // q = 3
    result_by_lpq.set(c_Key3(6, 5, 3), common_term_10);
    result_by_q.set(c_Key1(3), common_term_10);
    // q = 4
    result_by_lpq.set(c_Key3(6, 5, 4), common_term_9);
    result_by_q.set(c_Key1(4), common_term_9);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 5), result_by_q);
    result_by_q.clear();

    // l , p = (6, 6).
    // q = -4
    result_by_lpq.set(c_Key3(6, 6, -4), common_term_8);
    result_by_q.set(c_Key1(-4), common_term_8);
    // q = -3
    result_by_lpq.set(c_Key3(6, 6, -3), common_term_7);
    result_by_q.set(c_Key1(-3), common_term_7);
    // q = -2
    result_by_lpq.set(c_Key3(6, 6, -2), common_term_6);
    result_by_q.set(c_Key1(-2), common_term_6);
    // q = -1
    result_by_lpq.set(c_Key3(6, 6, -1), common_term_5);
    result_by_q.set(c_Key1(-1), common_term_5);
    // q = 0
    result_by_lpq.set(c_Key3(6, 6, 0), common_term_4);
    result_by_q.set(c_Key1(0), common_term_4);
    // q = 1
    result_by_lpq.set(c_Key3(6, 6, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // q = 2
    result_by_lpq.set(c_Key3(6, 6, 2), common_term_2);
    result_by_q.set(c_Key1(2), common_term_2);
    // q = 3
    result_by_lpq.set(c_Key3(6, 6, 3), common_term_1);
    result_by_q.set(c_Key1(3), common_term_1);
    // q = 4
    result_by_lpq.set(c_Key3(6, 6, 4), common_term_0);
    result_by_q.set(c_Key1(4), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 6), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l6_e10(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 6.
    // Functions are taylor expanded and truncated at eccentricity^10.

    // Functions Calculated for l = 6.

    c_IntMap<c_Key3, double> result_by_lpq(131);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(7);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_8 = eccentricity_4 * eccentricity_4;
    double eccentricity_9 = eccentricity_8 * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double eccentricity_6 = eccentricity_3 * eccentricity_3;
    double eccentricity_7 = eccentricity_6 * eccentricity;
    double eccentricity_5 = eccentricity_4 * eccentricity;
    double common_term_0 = 0.00010593959263392857*eccentricity_9;
    double common_term_1 = 2.4801587301587302e-5*eccentricity_8;
    double common_term_2 = 2.8579954117063492e-6*eccentricity_9 + 1.5500992063492063e-6*eccentricity_7;
    double common_term_3 = -0.00023038349454365079*eccentricity_9 - 0.00024956597222222222*eccentricity_7 - 0.00026041666666666667*eccentricity_5;
    double common_term_4 = 0.0064236111111111111*eccentricity_8 + 0.0041666666666666667*eccentricity_6 + 0.041666666666666667*eccentricity_4;
    double common_term_5 = 0.06251220703125*eccentricity_9 - 0.13447265625*eccentricity_7 + 0.73828125*eccentricity_5 - 0.5625*eccentricity_3;
    double common_term_6 = -2.6513888888888889*eccentricity_8 + 7.7083333333333333*eccentricity_6 - 7.6666666666666667*eccentricity_4 + 2.0*eccentricity_2;
    double common_term_7 = -25.022871229383681*eccentricity_9 + 55.931260850694444*eccentricity_7 - 57.669270833333333*eccentricity_5 + 22.8125*eccentricity_3 - 2.5*eccentricity;
    double common_term_8 = 312.66796875*eccentricity_8 - 314.96875*eccentricity_6 + 145.6875*eccentricity_4 - 25.5*eccentricity_2 + 1.0;
    double common_term_9 = 1443.2193650987413*eccentricity_9 - 1393.6566297743056*eccentricity_7 + 687.33072916666667*eccentricity_5 - 147.4375*eccentricity_3 + 9.5*eccentricity;
    double common_term_10 = -5298.1472222222222*eccentricity_8 + 2669.9583333333333*eccentricity_6 - 635.66666666666667*eccentricity_4 + 51.5*eccentricity_2;
    double common_term_11 = -17943.421472167969*eccentricity_9 + 9036.48076171875*eccentricity_7 - 2273.73046875*eccentricity_5 + 209.1875*eccentricity_3;
    double common_term_12 = 27569.268663194444*eccentricity_8 - 7131.3020833333333*eccentricity_6 + 707.60416666666667*eccentricity_4;
    double common_term_13 = 77533.185511804006*eccentricity_9 - 20271.980327690972*eccentricity_7 + 2105.3309895833333*eccentricity_5;
    double common_term_14 = -53371.869642857143*eccentricity_8 + 5692.6*eccentricity_6;
    double common_term_15 = -132109.62696073017*eccentricity_9 + 14289.359655567956*eccentricity_7;
    double common_term_16 = 33793.725173611111*eccentricity_8;
    double common_term_17 = 76106.604036603655*eccentricity_9;
    double common_term_18 = 1.4557787644589809*eccentricity_9;
    double common_term_19 = 1.0689856150793651*eccentricity_8;
    double common_term_20 = 3.7278908865792411*eccentricity_9 + 0.78542131696428571*eccentricity_7;
    double common_term_21 = 2.8865575396825397*eccentricity_8 + 0.57743055555555556*eccentricity_6;
    double common_term_22 = 7.0080796983506944*eccentricity_9 + 2.2298502604166667*eccentricity_7 + 0.42473958333333333*eccentricity_5;
    double common_term_23 = 0.3125*eccentricity_4*std::pow(1.0 - eccentricity_2, -5.5);
    double common_term_24 = 11.385783329716435*eccentricity_9 + 4.4400716145833333*eccentricity_7 + 1.3216145833333333*eccentricity_5 + 0.22916666666666667*eccentricity_3;
    double common_term_25 = 9.2607638888888889*eccentricity_8 + 3.5260416666666667*eccentricity_6 + eccentricity_4 + 0.25*eccentricity_2;
    double common_term_26 = 16.961810302734375*eccentricity_9 + 7.60693359375*eccentricity_7 + 2.3828125*eccentricity_5 + 1.9375*eccentricity_3 - 0.5*eccentricity;
    double common_term_27 = 15.146267361111111*eccentricity_8 + 1.8888888888888889*eccentricity_6 + 11.0*eccentricity_4 - 5.5*eccentricity_2 + 1.0;
    double common_term_28 = 32.950558132595486*eccentricity_9 - 15.68603515625*eccentricity_7 + 51.549479166666667*eccentricity_5 - 31.5625*eccentricity_3 + 7.5*eccentricity;
    double common_term_29 = -113.184375*eccentricity_8 + 204.765625*eccentricity_6 - 131.0*eccentricity_4 + 33.25*eccentricity_2;
    double common_term_30 = -504.1496912073206*eccentricity_9 + 708.82418619791667*eccentricity_7 - 445.40755208333333*eccentricity_5 + 113.39583333333333*eccentricity_3;
    double common_term_31 = 2195.2217013888889*eccentricity_8 - 1320.9854166666667*eccentricity_6 + 328.6875*eccentricity_4;
    double common_term_32 = 6214.475967843192*eccentricity_9 - 3544.72041015625*eccentricity_7 + 851.72890625*eccentricity_5;
    double common_term_33 = -8809.5994543650794*eccentricity_8 + 2032.5086805555556*eccentricity_6;
    double common_term_34 = -20602.832414947994*eccentricity_9 + 4552.519982328869*eccentricity_7;
    double common_term_35 = 9696.6136160714286*eccentricity_8;
    double common_term_36 = 19824.016237197688*eccentricity_9;
    double common_term_37 = 49.231494630413291*eccentricity_9;
    double common_term_38 = 33.245424107142857*eccentricity_8;
    double common_term_39 = 84.318628462534102*eccentricity_9 + 22.286900111607143*eccentricity_7;
    double common_term_40 = 61.431349206349206*eccentricity_8 + 14.809722222222222*eccentricity_6;
    double common_term_41 = 127.20421840122768*eccentricity_9 + 44.17548828125*eccentricity_7 + 9.73515625*eccentricity_5;
    double common_term_42 = 95.336545138888889*eccentricity_8 + 31.327083333333333*eccentricity_6 + 6.3125*eccentricity_4;
    double common_term_43 = 175.57016556351273*eccentricity_9 + 70.590397135416667*eccentricity_7 + 21.873697916666667*eccentricity_5 + 4.0208333333333333*eccentricity_3;
    double common_term_44 = std::pow(1.0 - eccentricity_2, -5.5)*(1.25*eccentricity_4 + 2.5*eccentricity_2);
    double common_term_45 = 229.31607733832465*eccentricity_9 + 101.17106119791667*eccentricity_7 + 37.080729166666667*eccentricity_5 + 10.0625*eccentricity_3 + 1.5*eccentricity;
    double common_term_46 = 177.49001736111111*eccentricity_8 + 75.336805555555556*eccentricity_6 + 26.1875*eccentricity_4 + 6.5*eccentricity_2 + 1.0;
    double common_term_47 = 288.33571166992188*eccentricity_9 + 135.78759765625*eccentricity_7 + 55.4140625*eccentricity_5 + 17.3125*eccentricity_3 + 5.5*eccentricity;
    double common_term_48 = 225.23888888888889*eccentricity_8 + 104.10416666666667*eccentricity_6 + 35.0*eccentricity_4 + 19.0*eccentricity_2;
    double common_term_49 = 350.33441614221644*eccentricity_9 + 183.17464192708333*eccentricity_7 + 55.384114583333333*eccentricity_5 + 52.604166666666667*eccentricity_3;
    double common_term_50 = 315.759375*eccentricity_8 + 60.225*eccentricity_6 + 127.5*eccentricity_4;
    double common_term_51 = 556.04554714626736*eccentricity_9 - 0.72203776041666667*eccentricity_7 + 282.51223958333333*eccentricity_5;
    double common_term_52 = -244.28353174603175*eccentricity_8 + 586.63472222222222*eccentricity_6;
    double common_term_53 = -917.96151515415737*eccentricity_9 + 1159.5359514508929*eccentricity_7;
    double common_term_54 = 2204.7194320436508*eccentricity_8;
    double common_term_55 = 4062.6052557997392*eccentricity_9;
    double common_term_56 = 588.71616516113281*eccentricity_9;
    double common_term_57 = 356.50512152777778*eccentricity_8;
    double common_term_58 = 416.56784430609809*eccentricity_9 + 211.64350043402778*eccentricity_7;
    double common_term_59 = 293.68125*eccentricity_8 + 122.5875*eccentricity_6;
    double common_term_60 = 477.65137396918403*eccentricity_9 + 197.73822699652778*eccentricity_7 + 68.805989583333333*eccentricity_5;
    double common_term_61 = 329.14340277777778*eccentricity_8 + 127.07916666666667*eccentricity_6 + 37.041666666666667*eccentricity_4;
    double common_term_62 = 503.28966064453125*eccentricity_9 + 219.36181640625*eccentricity_7 + 77.51953125*eccentricity_5 + 18.8125*eccentricity_3;
    double common_term_63 = 346.21319444444444*eccentricity_8 + 140.14583333333333*eccentricity_6 + 44.333333333333333*eccentricity_4 + 8.75*eccentricity_2;
    double common_term_64 = 516.27495591905382*eccentricity_9 + 229.25645616319444*eccentricity_7 + 84.674479166666667*eccentricity_5 + 23.1875*eccentricity_3 + 3.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(19);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (6, 0).
    // q = -9
    result_by_lpq.set(c_Key3(6, 0, -9), common_term_0);
    result_by_q.set(c_Key1(-9), common_term_0);
    // q = -8
    result_by_lpq.set(c_Key3(6, 0, -8), common_term_1);
    result_by_q.set(c_Key1(-8), common_term_1);
    // q = -7
    result_by_lpq.set(c_Key3(6, 0, -7), common_term_2);
    result_by_q.set(c_Key1(-7), common_term_2);
    // q = -5
    result_by_lpq.set(c_Key3(6, 0, -5), common_term_3);
    result_by_q.set(c_Key1(-5), common_term_3);
    // q = -4
    result_by_lpq.set(c_Key3(6, 0, -4), common_term_4);
    result_by_q.set(c_Key1(-4), common_term_4);
    // q = -3
    result_by_lpq.set(c_Key3(6, 0, -3), common_term_5);
    result_by_q.set(c_Key1(-3), common_term_5);
    // q = -2
    result_by_lpq.set(c_Key3(6, 0, -2), common_term_6);
    result_by_q.set(c_Key1(-2), common_term_6);
    // q = -1
    result_by_lpq.set(c_Key3(6, 0, -1), common_term_7);
    result_by_q.set(c_Key1(-1), common_term_7);
    // q = 0
    result_by_lpq.set(c_Key3(6, 0, 0), common_term_8);
    result_by_q.set(c_Key1(0), common_term_8);
    // q = 1
    result_by_lpq.set(c_Key3(6, 0, 1), common_term_9);
    result_by_q.set(c_Key1(1), common_term_9);
    // q = 2
    result_by_lpq.set(c_Key3(6, 0, 2), common_term_10);
    result_by_q.set(c_Key1(2), common_term_10);
    // q = 3
    result_by_lpq.set(c_Key3(6, 0, 3), common_term_11);
    result_by_q.set(c_Key1(3), common_term_11);
    // q = 4
    result_by_lpq.set(c_Key3(6, 0, 4), common_term_12);
    result_by_q.set(c_Key1(4), common_term_12);
    // q = 5
    result_by_lpq.set(c_Key3(6, 0, 5), common_term_13);
    result_by_q.set(c_Key1(5), common_term_13);
    // q = 6
    result_by_lpq.set(c_Key3(6, 0, 6), common_term_14);
    result_by_q.set(c_Key1(6), common_term_14);
    // q = 7
    result_by_lpq.set(c_Key3(6, 0, 7), common_term_15);
    result_by_q.set(c_Key1(7), common_term_15);
    // q = 8
    result_by_lpq.set(c_Key3(6, 0, 8), common_term_16);
    result_by_q.set(c_Key1(8), common_term_16);
    // q = 9
    result_by_lpq.set(c_Key3(6, 0, 9), common_term_17);
    result_by_q.set(c_Key1(9), common_term_17);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 0), result_by_q);
    result_by_q.clear();

    // l , p = (6, 1).
    // q = -9
    result_by_lpq.set(c_Key3(6, 1, -9), common_term_18);
    result_by_q.set(c_Key1(-9), common_term_18);
    // q = -8
    result_by_lpq.set(c_Key3(6, 1, -8), common_term_19);
    result_by_q.set(c_Key1(-8), common_term_19);
    // q = -7
    result_by_lpq.set(c_Key3(6, 1, -7), common_term_20);
    result_by_q.set(c_Key1(-7), common_term_20);
    // q = -6
    result_by_lpq.set(c_Key3(6, 1, -6), common_term_21);
    result_by_q.set(c_Key1(-6), common_term_21);
    // q = -5
    result_by_lpq.set(c_Key3(6, 1, -5), common_term_22);
    result_by_q.set(c_Key1(-5), common_term_22);
    // q = -4
    result_by_lpq.set(c_Key3(6, 1, -4), common_term_23);
    result_by_q.set(c_Key1(-4), common_term_23);
    // q = -3
    result_by_lpq.set(c_Key3(6, 1, -3), common_term_24);
    result_by_q.set(c_Key1(-3), common_term_24);
    // q = -2
    result_by_lpq.set(c_Key3(6, 1, -2), common_term_25);
    result_by_q.set(c_Key1(-2), common_term_25);
    // q = -1
    result_by_lpq.set(c_Key3(6, 1, -1), common_term_26);
    result_by_q.set(c_Key1(-1), common_term_26);
    // q = 0
    result_by_lpq.set(c_Key3(6, 1, 0), common_term_27);
    result_by_q.set(c_Key1(0), common_term_27);
    // q = 1
    result_by_lpq.set(c_Key3(6, 1, 1), common_term_28);
    result_by_q.set(c_Key1(1), common_term_28);
    // q = 2
    result_by_lpq.set(c_Key3(6, 1, 2), common_term_29);
    result_by_q.set(c_Key1(2), common_term_29);
    // q = 3
    result_by_lpq.set(c_Key3(6, 1, 3), common_term_30);
    result_by_q.set(c_Key1(3), common_term_30);
    // q = 4
    result_by_lpq.set(c_Key3(6, 1, 4), common_term_31);
    result_by_q.set(c_Key1(4), common_term_31);
    // q = 5
    result_by_lpq.set(c_Key3(6, 1, 5), common_term_32);
    result_by_q.set(c_Key1(5), common_term_32);
    // q = 6
    result_by_lpq.set(c_Key3(6, 1, 6), common_term_33);
    result_by_q.set(c_Key1(6), common_term_33);
    // q = 7
    result_by_lpq.set(c_Key3(6, 1, 7), common_term_34);
    result_by_q.set(c_Key1(7), common_term_34);
    // q = 8
    result_by_lpq.set(c_Key3(6, 1, 8), common_term_35);
    result_by_q.set(c_Key1(8), common_term_35);
    // q = 9
    result_by_lpq.set(c_Key3(6, 1, 9), common_term_36);
    result_by_q.set(c_Key1(9), common_term_36);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 1), result_by_q);
    result_by_q.clear();

    // l , p = (6, 2).
    // q = -9
    result_by_lpq.set(c_Key3(6, 2, -9), common_term_37);
    result_by_q.set(c_Key1(-9), common_term_37);
    // q = -8
    result_by_lpq.set(c_Key3(6, 2, -8), common_term_38);
    result_by_q.set(c_Key1(-8), common_term_38);
    // q = -7
    result_by_lpq.set(c_Key3(6, 2, -7), common_term_39);
    result_by_q.set(c_Key1(-7), common_term_39);
    // q = -6
    result_by_lpq.set(c_Key3(6, 2, -6), common_term_40);
    result_by_q.set(c_Key1(-6), common_term_40);
    // q = -5
    result_by_lpq.set(c_Key3(6, 2, -5), common_term_41);
    result_by_q.set(c_Key1(-5), common_term_41);
    // q = -4
    result_by_lpq.set(c_Key3(6, 2, -4), common_term_42);
    result_by_q.set(c_Key1(-4), common_term_42);
    // q = -3
    result_by_lpq.set(c_Key3(6, 2, -3), common_term_43);
    result_by_q.set(c_Key1(-3), common_term_43);
    // q = -2
    result_by_lpq.set(c_Key3(6, 2, -2), common_term_44);
    result_by_q.set(c_Key1(-2), common_term_44);
    // q = -1
    result_by_lpq.set(c_Key3(6, 2, -1), common_term_45);
    result_by_q.set(c_Key1(-1), common_term_45);
    // q = 0
    result_by_lpq.set(c_Key3(6, 2, 0), common_term_46);
    result_by_q.set(c_Key1(0), common_term_46);
    // q = 1
    result_by_lpq.set(c_Key3(6, 2, 1), common_term_47);
    result_by_q.set(c_Key1(1), common_term_47);
    // q = 2
    result_by_lpq.set(c_Key3(6, 2, 2), common_term_48);
    result_by_q.set(c_Key1(2), common_term_48);
    // q = 3
    result_by_lpq.set(c_Key3(6, 2, 3), common_term_49);
    result_by_q.set(c_Key1(3), common_term_49);
    // q = 4
    result_by_lpq.set(c_Key3(6, 2, 4), common_term_50);
    result_by_q.set(c_Key1(4), common_term_50);
    // q = 5
    result_by_lpq.set(c_Key3(6, 2, 5), common_term_51);
    result_by_q.set(c_Key1(5), common_term_51);
    // q = 6
    result_by_lpq.set(c_Key3(6, 2, 6), common_term_52);
    result_by_q.set(c_Key1(6), common_term_52);
    // q = 7
    result_by_lpq.set(c_Key3(6, 2, 7), common_term_53);
    result_by_q.set(c_Key1(7), common_term_53);
    // q = 8
    result_by_lpq.set(c_Key3(6, 2, 8), common_term_54);
    result_by_q.set(c_Key1(8), common_term_54);
    // q = 9
    result_by_lpq.set(c_Key3(6, 2, 9), common_term_55);
    result_by_q.set(c_Key1(9), common_term_55);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 2), result_by_q);
    result_by_q.clear();

    // l , p = (6, 3).
    // q = -9
    result_by_lpq.set(c_Key3(6, 3, -9), common_term_56);
    result_by_q.set(c_Key1(-9), common_term_56);
    // q = -8
    result_by_lpq.set(c_Key3(6, 3, -8), common_term_57);
    result_by_q.set(c_Key1(-8), common_term_57);
    // q = -7
    result_by_lpq.set(c_Key3(6, 3, -7), common_term_58);
    result_by_q.set(c_Key1(-7), common_term_58);
    // q = -6
    result_by_lpq.set(c_Key3(6, 3, -6), common_term_59);
    result_by_q.set(c_Key1(-6), common_term_59);
    // q = -5
    result_by_lpq.set(c_Key3(6, 3, -5), common_term_60);
    result_by_q.set(c_Key1(-5), common_term_60);
    // q = -4
    result_by_lpq.set(c_Key3(6, 3, -4), common_term_61);
    result_by_q.set(c_Key1(-4), common_term_61);
    // q = -3
    result_by_lpq.set(c_Key3(6, 3, -3), common_term_62);
    result_by_q.set(c_Key1(-3), common_term_62);
    // q = -2
    result_by_lpq.set(c_Key3(6, 3, -2), common_term_63);
    result_by_q.set(c_Key1(-2), common_term_63);
    // q = -1
    result_by_lpq.set(c_Key3(6, 3, -1), common_term_64);
    result_by_q.set(c_Key1(-1), common_term_64);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -5.5)*(1.875*eccentricity_4 + 5.0*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(6, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(6, 3, 1), common_term_64);
    result_by_q.set(c_Key1(1), common_term_64);
    // q = 2
    result_by_lpq.set(c_Key3(6, 3, 2), common_term_63);
    result_by_q.set(c_Key1(2), common_term_63);
    // q = 3
    result_by_lpq.set(c_Key3(6, 3, 3), common_term_62);
    result_by_q.set(c_Key1(3), common_term_62);
    // q = 4
    result_by_lpq.set(c_Key3(6, 3, 4), common_term_61);
    result_by_q.set(c_Key1(4), common_term_61);
    // q = 5
    result_by_lpq.set(c_Key3(6, 3, 5), common_term_60);
    result_by_q.set(c_Key1(5), common_term_60);
    // q = 6
    result_by_lpq.set(c_Key3(6, 3, 6), common_term_59);
    result_by_q.set(c_Key1(6), common_term_59);
    // q = 7
    result_by_lpq.set(c_Key3(6, 3, 7), common_term_58);
    result_by_q.set(c_Key1(7), common_term_58);
    // q = 8
    result_by_lpq.set(c_Key3(6, 3, 8), common_term_57);
    result_by_q.set(c_Key1(8), common_term_57);
    // q = 9
    result_by_lpq.set(c_Key3(6, 3, 9), common_term_56);
    result_by_q.set(c_Key1(9), common_term_56);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 3), result_by_q);
    result_by_q.clear();

    // l , p = (6, 4).
    // q = -9
    result_by_lpq.set(c_Key3(6, 4, -9), common_term_55);
    result_by_q.set(c_Key1(-9), common_term_55);
    // q = -8
    result_by_lpq.set(c_Key3(6, 4, -8), common_term_54);
    result_by_q.set(c_Key1(-8), common_term_54);
    // q = -7
    result_by_lpq.set(c_Key3(6, 4, -7), common_term_53);
    result_by_q.set(c_Key1(-7), common_term_53);
    // q = -6
    result_by_lpq.set(c_Key3(6, 4, -6), common_term_52);
    result_by_q.set(c_Key1(-6), common_term_52);
    // q = -5
    result_by_lpq.set(c_Key3(6, 4, -5), common_term_51);
    result_by_q.set(c_Key1(-5), common_term_51);
    // q = -4
    result_by_lpq.set(c_Key3(6, 4, -4), common_term_50);
    result_by_q.set(c_Key1(-4), common_term_50);
    // q = -3
    result_by_lpq.set(c_Key3(6, 4, -3), common_term_49);
    result_by_q.set(c_Key1(-3), common_term_49);
    // q = -2
    result_by_lpq.set(c_Key3(6, 4, -2), common_term_48);
    result_by_q.set(c_Key1(-2), common_term_48);
    // q = -1
    result_by_lpq.set(c_Key3(6, 4, -1), common_term_47);
    result_by_q.set(c_Key1(-1), common_term_47);
    // q = 0
    result_by_lpq.set(c_Key3(6, 4, 0), common_term_46);
    result_by_q.set(c_Key1(0), common_term_46);
    // q = 1
    result_by_lpq.set(c_Key3(6, 4, 1), common_term_45);
    result_by_q.set(c_Key1(1), common_term_45);
    // q = 2
    result_by_lpq.set(c_Key3(6, 4, 2), common_term_44);
    result_by_q.set(c_Key1(2), common_term_44);
    // q = 3
    result_by_lpq.set(c_Key3(6, 4, 3), common_term_43);
    result_by_q.set(c_Key1(3), common_term_43);
    // q = 4
    result_by_lpq.set(c_Key3(6, 4, 4), common_term_42);
    result_by_q.set(c_Key1(4), common_term_42);
    // q = 5
    result_by_lpq.set(c_Key3(6, 4, 5), common_term_41);
    result_by_q.set(c_Key1(5), common_term_41);
    // q = 6
    result_by_lpq.set(c_Key3(6, 4, 6), common_term_40);
    result_by_q.set(c_Key1(6), common_term_40);
    // q = 7
    result_by_lpq.set(c_Key3(6, 4, 7), common_term_39);
    result_by_q.set(c_Key1(7), common_term_39);
    // q = 8
    result_by_lpq.set(c_Key3(6, 4, 8), common_term_38);
    result_by_q.set(c_Key1(8), common_term_38);
    // q = 9
    result_by_lpq.set(c_Key3(6, 4, 9), common_term_37);
    result_by_q.set(c_Key1(9), common_term_37);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 4), result_by_q);
    result_by_q.clear();

    // l , p = (6, 5).
    // q = -9
    result_by_lpq.set(c_Key3(6, 5, -9), common_term_36);
    result_by_q.set(c_Key1(-9), common_term_36);
    // q = -8
    result_by_lpq.set(c_Key3(6, 5, -8), common_term_35);
    result_by_q.set(c_Key1(-8), common_term_35);
    // q = -7
    result_by_lpq.set(c_Key3(6, 5, -7), common_term_34);
    result_by_q.set(c_Key1(-7), common_term_34);
    // q = -6
    result_by_lpq.set(c_Key3(6, 5, -6), common_term_33);
    result_by_q.set(c_Key1(-6), common_term_33);
    // q = -5
    result_by_lpq.set(c_Key3(6, 5, -5), common_term_32);
    result_by_q.set(c_Key1(-5), common_term_32);
    // q = -4
    result_by_lpq.set(c_Key3(6, 5, -4), common_term_31);
    result_by_q.set(c_Key1(-4), common_term_31);
    // q = -3
    result_by_lpq.set(c_Key3(6, 5, -3), common_term_30);
    result_by_q.set(c_Key1(-3), common_term_30);
    // q = -2
    result_by_lpq.set(c_Key3(6, 5, -2), common_term_29);
    result_by_q.set(c_Key1(-2), common_term_29);
    // q = -1
    result_by_lpq.set(c_Key3(6, 5, -1), common_term_28);
    result_by_q.set(c_Key1(-1), common_term_28);
    // q = 0
    result_by_lpq.set(c_Key3(6, 5, 0), common_term_27);
    result_by_q.set(c_Key1(0), common_term_27);
    // q = 1
    result_by_lpq.set(c_Key3(6, 5, 1), common_term_26);
    result_by_q.set(c_Key1(1), common_term_26);
    // q = 2
    result_by_lpq.set(c_Key3(6, 5, 2), common_term_25);
    result_by_q.set(c_Key1(2), common_term_25);
    // q = 3
    result_by_lpq.set(c_Key3(6, 5, 3), common_term_24);
    result_by_q.set(c_Key1(3), common_term_24);
    // q = 4
    result_by_lpq.set(c_Key3(6, 5, 4), common_term_23);
    result_by_q.set(c_Key1(4), common_term_23);
    // q = 5
    result_by_lpq.set(c_Key3(6, 5, 5), common_term_22);
    result_by_q.set(c_Key1(5), common_term_22);
    // q = 6
    result_by_lpq.set(c_Key3(6, 5, 6), common_term_21);
    result_by_q.set(c_Key1(6), common_term_21);
    // q = 7
    result_by_lpq.set(c_Key3(6, 5, 7), common_term_20);
    result_by_q.set(c_Key1(7), common_term_20);
    // q = 8
    result_by_lpq.set(c_Key3(6, 5, 8), common_term_19);
    result_by_q.set(c_Key1(8), common_term_19);
    // q = 9
    result_by_lpq.set(c_Key3(6, 5, 9), common_term_18);
    result_by_q.set(c_Key1(9), common_term_18);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 5), result_by_q);
    result_by_q.clear();

    // l , p = (6, 6).
    // q = -9
    result_by_lpq.set(c_Key3(6, 6, -9), common_term_17);
    result_by_q.set(c_Key1(-9), common_term_17);
    // q = -8
    result_by_lpq.set(c_Key3(6, 6, -8), common_term_16);
    result_by_q.set(c_Key1(-8), common_term_16);
    // q = -7
    result_by_lpq.set(c_Key3(6, 6, -7), common_term_15);
    result_by_q.set(c_Key1(-7), common_term_15);
    // q = -6
    result_by_lpq.set(c_Key3(6, 6, -6), common_term_14);
    result_by_q.set(c_Key1(-6), common_term_14);
    // q = -5
    result_by_lpq.set(c_Key3(6, 6, -5), common_term_13);
    result_by_q.set(c_Key1(-5), common_term_13);
    // q = -4
    result_by_lpq.set(c_Key3(6, 6, -4), common_term_12);
    result_by_q.set(c_Key1(-4), common_term_12);
    // q = -3
    result_by_lpq.set(c_Key3(6, 6, -3), common_term_11);
    result_by_q.set(c_Key1(-3), common_term_11);
    // q = -2
    result_by_lpq.set(c_Key3(6, 6, -2), common_term_10);
    result_by_q.set(c_Key1(-2), common_term_10);
    // q = -1
    result_by_lpq.set(c_Key3(6, 6, -1), common_term_9);
    result_by_q.set(c_Key1(-1), common_term_9);
    // q = 0
    result_by_lpq.set(c_Key3(6, 6, 0), common_term_8);
    result_by_q.set(c_Key1(0), common_term_8);
    // q = 1
    result_by_lpq.set(c_Key3(6, 6, 1), common_term_7);
    result_by_q.set(c_Key1(1), common_term_7);
    // q = 2
    result_by_lpq.set(c_Key3(6, 6, 2), common_term_6);
    result_by_q.set(c_Key1(2), common_term_6);
    // q = 3
    result_by_lpq.set(c_Key3(6, 6, 3), common_term_5);
    result_by_q.set(c_Key1(3), common_term_5);
    // q = 4
    result_by_lpq.set(c_Key3(6, 6, 4), common_term_4);
    result_by_q.set(c_Key1(4), common_term_4);
    // q = 5
    result_by_lpq.set(c_Key3(6, 6, 5), common_term_3);
    result_by_q.set(c_Key1(5), common_term_3);
    // q = 7
    result_by_lpq.set(c_Key3(6, 6, 7), common_term_2);
    result_by_q.set(c_Key1(7), common_term_2);
    // q = 8
    result_by_lpq.set(c_Key3(6, 6, 8), common_term_1);
    result_by_q.set(c_Key1(8), common_term_1);
    // q = 9
    result_by_lpq.set(c_Key3(6, 6, 9), common_term_0);
    result_by_q.set(c_Key1(9), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 6), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l6_e15(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 6.
    // Functions are taylor expanded and truncated at eccentricity^15.

    // Functions Calculated for l = 6.

    c_IntMap<c_Key3, double> result_by_lpq(201);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(7);
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
    double common_term_0 = 0.0030791548251865712*eccentricity_14;
    double common_term_1 = 0.001899346759560887*eccentricity_13;
    double common_term_2 = 0.00243231377997003*eccentricity_14 + 0.001109476461038961*eccentricity_12;
    double common_term_3 = 0.0013314563022267142*eccentricity_13 + 0.00059728880847553536*eccentricity_11;
    double common_term_4 = 0.00094222649778205334*eccentricity_14 + 0.00062850729517396184*eccentricity_12 + 0.00028218694885361552*eccentricity_10;
    double common_term_5 = 0.0003350640581799792*eccentricity_13 + 0.00023041861397879464*eccentricity_11 + 0.00010593959263392857*eccentricity_9;
    double common_term_6 = 8.2929263701833146e-5*eccentricity_14 + 7.0512290564373898e-5*eccentricity_12 + 5.0981040564373898e-5*eccentricity_10 + 2.4801587301587302e-5*eccentricity_8;
    double common_term_7 = 4.0465337539125652e-6*eccentricity_13 + 3.6498646677276235e-6*eccentricity_11 + 2.8579954117063492e-6*eccentricity_9 + 1.5500992063492063e-6*eccentricity_7;
    double common_term_8 = -0.00018671935337701312*eccentricity_13 - 0.0002079050376932457*eccentricity_11 - 0.00023038349454365079*eccentricity_9 - 0.00024956597222222222*eccentricity_7 - 0.00026041666666666667*eccentricity_5;
    double common_term_9 = 0.0027452027300117578*eccentricity_14 + 0.0034391534391534392*eccentricity_12 + 0.0044436177248677249*eccentricity_10 + 0.0064236111111111111*eccentricity_8 + 0.0041666666666666667*eccentricity_6 + 0.041666666666666667*eccentricity_4;
    double common_term_10 = 0.022001722199576242*eccentricity_13 + 0.025369589669363839*eccentricity_11 + 0.06251220703125*eccentricity_9 - 0.13447265625*eccentricity_7 + 0.73828125*eccentricity_5 - 0.5625*eccentricity_3;
    double common_term_11 = 0.077645640432098765*eccentricity_14 + 0.011097056878306878*eccentricity_12 + 0.69415509259259259*eccentricity_10 - 2.6513888888888889*eccentricity_8 + 7.7083333333333333*eccentricity_6 - 7.6666666666666667*eccentricity_4 + 2.0*eccentricity_2;
    double common_term_12 = -1.0275152036542405*eccentricity_13 + 7.129992025869864*eccentricity_11 - 25.022871229383681*eccentricity_9 + 55.931260850694444*eccentricity_7 - 57.669270833333333*eccentricity_5 + 22.8125*eccentricity_3 - 2.5*eccentricity;
    double common_term_13 = -12.093317721619898*eccentricity_14 + 55.52111328125*eccentricity_12 - 165.161953125*eccentricity_10 + 312.66796875*eccentricity_8 - 314.96875*eccentricity_6 + 145.6875*eccentricity_4 - 25.5*eccentricity_2 + 1.0;
    double common_term_14 = 338.7226141924917*eccentricity_13 - 864.59438154997649*eccentricity_11 + 1443.2193650987413*eccentricity_9 - 1393.6566297743056*eccentricity_7 + 687.33072916666667*eccentricity_5 - 147.4375*eccentricity_3 + 9.5*eccentricity;
    double common_term_15 = 1711.4160550319665*eccentricity_14 - 3828.2373974867725*eccentricity_12 + 5759.9242766203704*eccentricity_10 - 5298.1472222222222*eccentricity_8 + 2669.9583333333333*eccentricity_6 - 635.66666666666667*eccentricity_4 + 51.5*eccentricity_2;
    double common_term_16 = -14912.162450819697*eccentricity_13 + 20512.371441977365*eccentricity_11 - 17943.421472167969*eccentricity_9 + 9036.48076171875*eccentricity_7 - 2273.73046875*eccentricity_5 + 209.1875*eccentricity_3;
    double common_term_17 = -52466.461716998962*eccentricity_14 + 66651.739275896991*eccentricity_12 - 55466.796926669974*eccentricity_10 + 27569.268663194444*eccentricity_8 - 7131.3020833333333*eccentricity_6 + 707.60416666666667*eccentricity_4;
    double common_term_18 = 200856.68874116792*eccentricity_13 - 159206.19147035912*eccentricity_11 + 77533.185511804006*eccentricity_9 - 20271.980327690972*eccentricity_7 + 2105.3309895833333*eccentricity_5;
    double common_term_19 = 568293.86935546875*eccentricity_14 - 429709.44654017857*eccentricity_12 + 204154.14174107143*eccentricity_10 - 53371.869642857143*eccentricity_8 + 5692.6*eccentricity_6;
    double common_term_20 = -1101154.6901831306*eccentricity_13 + 509067.61315316565*eccentricity_11 - 132109.62696073017*eccentricity_9 + 14289.359655567956*eccentricity_7;
    double common_term_21 = -2699160.2059729639*eccentricity_14 + 1212454.643498505*eccentricity_12 - 310806.6863908179*eccentricity_10 + 33793.725173611111*eccentricity_8;
    double common_term_22 = 2776620.1823720529*eccentricity_13 - 700696.61881991795*eccentricity_11 + 76106.604036603655*eccentricity_9;
    double common_term_23 = 6146345.1867708125*eccentricity_14 - 1523325.8970551748*eccentricity_12 + 164531.88032104277*eccentricity_10;
    double common_term_24 = -3209522.8770327859*eccentricity_13 + 343556.26409161624*eccentricity_11;
    double common_term_25 = -6579774.7889591073*eccentricity_14 + 696265.66076045049*eccentricity_12;
    double common_term_26 = 1374907.1380696101*eccentricity_13;
    double common_term_27 = 2653822.990858011*eccentricity_14;
    double common_term_28 = 6.8610594833284023*eccentricity_14;
    double common_term_29 = 5.0289451730453766*eccentricity_13;
    double common_term_30 = 12.774890073295564*eccentricity_14 + 3.6869184502723999*eccentricity_12;
    double common_term_31 = 10.065197431835127*eccentricity_13 + 2.7038058547888005*eccentricity_11;
    double common_term_32 = 20.823311339539367*eccentricity_14 + 7.8942329545454545*eccentricity_12 + 1.9835435267857143*eccentricity_10;
    double common_term_33 = 16.857054417993809*eccentricity_13 + 6.1673492182709763*eccentricity_11 + 1.4557787644589809*eccentricity_9;
    double common_term_34 = 30.706129331300772*eccentricity_14 + 13.599795869157848*eccentricity_12 + 4.8019696593915344*eccentricity_10 + 1.0689856150793651*eccentricity_8;
    double common_term_35 = 25.328705772672381*eccentricity_13 + 10.936387797764369*eccentricity_11 + 3.7278908865792411*eccentricity_9 + 0.78542131696428571*eccentricity_7;
    double common_term_36 = 42.530252436526308*eccentricity_14 + 20.830027373603762*eccentricity_12 + 8.7675006200396825*eccentricity_10 + 2.8865575396825397*eccentricity_8 + 0.57743055555555556*eccentricity_6;
    double common_term_37 = 35.574591930202706*eccentricity_13 + 17.079967815661557*eccentricity_11 + 7.0080796983506944*eccentricity_9 + 2.2298502604166667*eccentricity_7 + 0.42473958333333333*eccentricity_5;
    double common_term_38 = 0.3125*eccentricity_4*std::pow(1.0 - eccentricity_2, -5.5);
    double common_term_39 = 47.683038475374696*eccentricity_13 + 24.687632868529628*eccentricity_11 + 11.385783329716435*eccentricity_9 + 4.4400716145833333*eccentricity_7 + 1.3216145833333333*eccentricity_5 + 0.22916666666666667*eccentricity_3;
    double common_term_40 = 72.355813095927028*eccentricity_14 + 40.227262318121693*eccentricity_12 + 20.4865234375*eccentricity_10 + 9.2607638888888889*eccentricity_8 + 3.5260416666666667*eccentricity_6 + eccentricity_4 + 0.25*eccentricity_2;
    double common_term_41 = 61.756888719286237*eccentricity_13 + 33.865264434814453*eccentricity_11 + 16.961810302734375*eccentricity_9 + 7.60693359375*eccentricity_7 + 2.3828125*eccentricity_5 + 1.9375*eccentricity_3 - 0.5*eccentricity;
    double common_term_42 = 90.570352541769022*eccentricity_14 + 52.633672357253086*eccentricity_12 + 28.299939236111111*eccentricity_10 + 15.146267361111111*eccentricity_8 + 1.8888888888888889*eccentricity_6 + 11.0*eccentricity_4 - 5.5*eccentricity_2 + 1.0;
    double common_term_43 = 78.23271772848866*eccentricity_13 + 42.88214436283818*eccentricity_11 + 32.950558132595486*eccentricity_9 - 15.68603515625*eccentricity_7 + 51.549479166666667*eccentricity_5 - 31.5625*eccentricity_3 + 7.5*eccentricity;
    double common_term_44 = 113.89916190011161*eccentricity_14 + 52.828828125*eccentricity_12 + 91.6787109375*eccentricity_10 - 113.184375*eccentricity_8 + 204.765625*eccentricity_6 - 131.0*eccentricity_4 + 33.25*eccentricity_2;
    double common_term_45 = 15.369647119486773*eccentricity_13 + 311.31128373322663*eccentricity_11 - 504.1496912073206*eccentricity_9 + 708.82418619791667*eccentricity_7 - 445.40755208333333*eccentricity_5 + 113.39583333333333*eccentricity_3;
    double common_term_46 = -246.34788596034961*eccentricity_14 + 1101.6653302228009*eccentricity_12 - 1826.6253100198413*eccentricity_10 + 2195.2217013888889*eccentricity_8 - 1320.9854166666667*eccentricity_6 + 328.6875*eccentricity_4;
    double common_term_47 = 3718.3121417897088*eccentricity_13 - 5816.4672563280378*eccentricity_11 + 6214.475967843192*eccentricity_9 - 3544.72041015625*eccentricity_7 + 851.72890625*eccentricity_5;
    double common_term_48 = 11691.499287349066*eccentricity_14 - 16867.581598370444*eccentricity_12 + 16351.393849206349*eccentricity_10 - 8809.5994543650794*eccentricity_8 + 2032.5086805555556*eccentricity_6;
    double common_term_49 = -45480.198617447125*eccentricity_13 + 40502.908957301414*eccentricity_11 - 20602.832414947994*eccentricity_9 + 4552.519982328869*eccentricity_7;
    double common_term_50 = -115592.32271763393*eccentricity_14 + 95387.253716517857*eccentricity_12 - 45858.797209821429*eccentricity_10 + 9696.6136160714286*eccentricity_8;
    double common_term_51 = 215239.10583281955*eccentricity_13 - 97969.551795844964*eccentricity_11 + 19824.016237197688*eccentricity_9;
    double common_term_52 = 468205.26353690259*eccentricity_14 - 202167.86227439324*eccentricity_12 + 39172.045639467593*eccentricity_10;
    double common_term_53 = -405000.10163826602*eccentricity_13 + 75209.476258785694*eccentricity_11;
    double common_term_54 = -790764.76856383999*eccentricity_14 + 140889.14879560632*eccentricity_12;
    double common_term_55 = 258358.9854943355*eccentricity_13;
    double common_term_56 = 465019.77231634019*eccentricity_14;
    double common_term_57 = 324.51430162471457*eccentricity_14;
    double common_term_58 = 224.37256107712711*eccentricity_13;
    double common_term_59 = 334.65154916808291*eccentricity_14 + 154.58385662229521*eccentricity_12;
    double common_term_60 = 262.05906939797587*eccentricity_13 + 106.07548423172592*eccentricity_11;
    double common_term_61 = 462.35437831366325*eccentricity_14 + 201.66168540564374*eccentricity_12 + 72.456307457010582*eccentricity_10;
    double common_term_62 = 363.61292241894971*eccentricity_13 + 152.84656421708472*eccentricity_11 + 49.231494630413291*eccentricity_9;
    double common_term_63 = 586.50778618354302*eccentricity_14 + 283.60770228794643*eccentricity_12 + 114.264453125*eccentricity_10 + 33.245424107142857*eccentricity_8;
    double common_term_64 = 468.57150143626716*eccentricity_13 + 219.23386247490239*eccentricity_11 + 84.318628462534102*eccentricity_9 + 22.286900111607143*eccentricity_7;
    double common_term_65 = 720.3200735550779*eccentricity_14 + 371.4345467739565*eccentricity_12 + 167.85452628968254*eccentricity_10 + 61.431349206349206*eccentricity_8 + 14.809722222222222*eccentricity_6;
    double common_term_66 = 582.03084909234728*eccentricity_13 + 291.99724464416504*eccentricity_11 + 127.20421840122768*eccentricity_9 + 44.17548828125*eccentricity_7 + 9.73515625*eccentricity_5;
    double common_term_67 = 863.15967910075875*eccentricity_14 + 466.87657231729497*eccentricity_12 + 227.51363467261905*eccentricity_10 + 95.336545138888889*eccentricity_8 + 31.327083333333333*eccentricity_6 + 6.3125*eccentricity_4;
    double common_term_68 = 703.71636338630051*eccentricity_13 + 371.5901080459514*eccentricity_11 + 175.57016556351273*eccentricity_9 + 70.590397135416667*eccentricity_7 + 21.873697916666667*eccentricity_5 + 4.0208333333333333*eccentricity_3;
    double common_term_69 = std::pow(1.0 - eccentricity_2, -5.5)*(1.25*eccentricity_4 + 2.5*eccentricity_2);
    double common_term_70 = 833.50103988374983*eccentricity_13 + 457.87624250058775*eccentricity_11 + 229.31607733832465*eccentricity_9 + 101.17106119791667*eccentricity_7 + 37.080729166666667*eccentricity_5 + 10.0625*eccentricity_3 + 1.5*eccentricity;
    double common_term_71 = 1175.4085234394684*eccentricity_14 + 679.98254677854938*eccentricity_12 + 364.97671006944444*eccentricity_10 + 177.49001736111111*eccentricity_8 + 75.336805555555556*eccentricity_6 + 26.1875*eccentricity_4 + 6.5*eccentricity_2 + 1.0;
    double common_term_72 = 971.27648251942226*eccentricity_13 + 550.74737289428711*eccentricity_11 + 288.33571166992188*eccentricity_9 + 135.78759765625*eccentricity_7 + 55.4140625*eccentricity_5 + 17.3125*eccentricity_3 + 5.5*eccentricity;
    double common_term_73 = 1344.6030678185626*eccentricity_14 + 797.42890294312169*eccentricity_12 + 442.60251736111111*eccentricity_10 + 225.23888888888889*eccentricity_8 + 104.10416666666667*eccentricity_6 + 35.0*eccentricity_4 + 19.0*eccentricity_2;
    double common_term_74 = 1116.9114152054307*eccentricity_13 + 650.50363742485248*eccentricity_11 + 350.33441614221644*eccentricity_9 + 183.17464192708333*eccentricity_7 + 55.384114583333333*eccentricity_5 + 52.604166666666667*eccentricity_3;
    double common_term_75 = 1521.9739501953125*eccentricity_14 + 924.69172014508929*eccentricity_12 + 513.92678571428571*eccentricity_10 + 315.759375*eccentricity_8 + 60.225*eccentricity_6 + 127.5*eccentricity_4;
    double common_term_76 = 1284.5458838412791*eccentricity_13 + 703.90580623767994*eccentricity_11 + 556.04554714626736*eccentricity_9 - 0.72203776041666667*eccentricity_7 + 282.51223958333333*eccentricity_5;
    double common_term_77 = 1770.0839864016112*eccentricity_14 + 863.15170120517343*eccentricity_12 + 1031.6815538194444*eccentricity_10 - 244.28353174603175*eccentricity_8 + 586.63472222222222*eccentricity_6;
    double common_term_78 = 815.92556585720607*eccentricity_13 + 2033.4356371743338*eccentricity_11 - 917.96151515415737*eccentricity_9 + 1159.5359514508929*eccentricity_7;
    double common_term_79 = 99.537103549080003*eccentricity_14 + 4195.9863647245646*eccentricity_12 - 2515.2845465443122*eccentricity_10 + 2204.7194320436508*eccentricity_8;
    double common_term_80 = 8849.209973816325*eccentricity_13 - 5978.6031004550133*eccentricity_11 + 4062.6052557997392*eccentricity_9;
    double common_term_81 = 18680.647027698864*eccentricity_14 - 13048.854811282468*eccentricity_12 + 7294.7463392857143*eccentricity_10;
    double common_term_82 = -26859.214041007832*eccentricity_13 + 12816.391140609443*eccentricity_11;
    double common_term_83 = -52928.991878443501*eccentricity_14 + 22103.683587078832*eccentricity_12;
    double common_term_84 = 37515.480413075451*eccentricity_13;
    double common_term_85 = 62790.335431863512*eccentricity_14;
    double common_term_86 = 5914.2904116453048*eccentricity_14;
    double common_term_87 = 3805.5325424087962*eccentricity_13;
    double common_term_88 = 487.45542326813811*eccentricity_14 + 2427.2028231534091*eccentricity_12;
    double common_term_89 = 803.49230794009848*eccentricity_13 + 1532.5469050381507*eccentricity_11;
    double common_term_90 = 2351.9291003477337*eccentricity_14 + 821.53377744458474*eccentricity_12 + 956.40387490354938*eccentricity_10;
    double common_term_91 = 1728.3684608225389*eccentricity_13 + 711.15719764709473*eccentricity_11 + 588.71616516113281*eccentricity_9;
    double common_term_92 = 2254.0798681403152*eccentricity_14 + 1274.7245322145062*eccentricity_12 + 561.71428915895062*eccentricity_10 + 356.50512152777778*eccentricity_8;
    double common_term_93 = 1736.0021989555908*eccentricity_13 + 934.16348894378285*eccentricity_11 + 416.56784430609809*eccentricity_9 + 211.64350043402778*eccentricity_7;
    double common_term_94 = 2337.4749462890625*eccentricity_14 + 1314.013671875*eccentricity_12 + 674.9208984375*eccentricity_10 + 293.68125*eccentricity_8 + 122.5875*eccentricity_6;
    double common_term_95 = 1797.6122005314493*eccentricity_13 + 976.08916423938892*eccentricity_11 + 477.65137396918403*eccentricity_9 + 197.73822699652778*eccentricity_7 + 68.805989583333333*eccentricity_5;
    double common_term_96 = 2392.0661497717335*eccentricity_14 + 1359.6920211226852*eccentricity_12 + 709.74508101851852*eccentricity_10 + 329.14340277777778*eccentricity_8 + 127.07916666666667*eccentricity_6 + 37.041666666666667*eccentricity_4;
    double common_term_97 = 1837.5225867605209*eccentricity_13 + 1008.8510215759277*eccentricity_11 + 503.28966064453125*eccentricity_9 + 219.36181640625*eccentricity_7 + 77.51953125*eccentricity_5 + 18.8125*eccentricity_3;
    double common_term_98 = 2424.8300771002122*eccentricity_14 + 1386.8702835648148*eccentricity_12 + 731.70609085648148*eccentricity_10 + 346.21319444444444*eccentricity_8 + 140.14583333333333*eccentricity_6 + 44.333333333333333*eccentricity_4 + 8.75*eccentricity_2;
    double common_term_99 = 1857.4697025810053*eccentricity_13 + 1025.192222555655*eccentricity_11 + 516.27495591905382*eccentricity_9 + 229.25645616319444*eccentricity_7 + 84.674479166666667*eccentricity_5 + 23.1875*eccentricity_3 + 3.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(29);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (6, 0).
    // q = -14
    result_by_lpq.set(c_Key3(6, 0, -14), common_term_0);
    result_by_q.set(c_Key1(-14), common_term_0);
    // q = -13
    result_by_lpq.set(c_Key3(6, 0, -13), common_term_1);
    result_by_q.set(c_Key1(-13), common_term_1);
    // q = -12
    result_by_lpq.set(c_Key3(6, 0, -12), common_term_2);
    result_by_q.set(c_Key1(-12), common_term_2);
    // q = -11
    result_by_lpq.set(c_Key3(6, 0, -11), common_term_3);
    result_by_q.set(c_Key1(-11), common_term_3);
    // q = -10
    result_by_lpq.set(c_Key3(6, 0, -10), common_term_4);
    result_by_q.set(c_Key1(-10), common_term_4);
    // q = -9
    result_by_lpq.set(c_Key3(6, 0, -9), common_term_5);
    result_by_q.set(c_Key1(-9), common_term_5);
    // q = -8
    result_by_lpq.set(c_Key3(6, 0, -8), common_term_6);
    result_by_q.set(c_Key1(-8), common_term_6);
    // q = -7
    result_by_lpq.set(c_Key3(6, 0, -7), common_term_7);
    result_by_q.set(c_Key1(-7), common_term_7);
    // q = -5
    result_by_lpq.set(c_Key3(6, 0, -5), common_term_8);
    result_by_q.set(c_Key1(-5), common_term_8);
    // q = -4
    result_by_lpq.set(c_Key3(6, 0, -4), common_term_9);
    result_by_q.set(c_Key1(-4), common_term_9);
    // q = -3
    result_by_lpq.set(c_Key3(6, 0, -3), common_term_10);
    result_by_q.set(c_Key1(-3), common_term_10);
    // q = -2
    result_by_lpq.set(c_Key3(6, 0, -2), common_term_11);
    result_by_q.set(c_Key1(-2), common_term_11);
    // q = -1
    result_by_lpq.set(c_Key3(6, 0, -1), common_term_12);
    result_by_q.set(c_Key1(-1), common_term_12);
    // q = 0
    result_by_lpq.set(c_Key3(6, 0, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(6, 0, 1), common_term_14);
    result_by_q.set(c_Key1(1), common_term_14);
    // q = 2
    result_by_lpq.set(c_Key3(6, 0, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // q = 3
    result_by_lpq.set(c_Key3(6, 0, 3), common_term_16);
    result_by_q.set(c_Key1(3), common_term_16);
    // q = 4
    result_by_lpq.set(c_Key3(6, 0, 4), common_term_17);
    result_by_q.set(c_Key1(4), common_term_17);
    // q = 5
    result_by_lpq.set(c_Key3(6, 0, 5), common_term_18);
    result_by_q.set(c_Key1(5), common_term_18);
    // q = 6
    result_by_lpq.set(c_Key3(6, 0, 6), common_term_19);
    result_by_q.set(c_Key1(6), common_term_19);
    // q = 7
    result_by_lpq.set(c_Key3(6, 0, 7), common_term_20);
    result_by_q.set(c_Key1(7), common_term_20);
    // q = 8
    result_by_lpq.set(c_Key3(6, 0, 8), common_term_21);
    result_by_q.set(c_Key1(8), common_term_21);
    // q = 9
    result_by_lpq.set(c_Key3(6, 0, 9), common_term_22);
    result_by_q.set(c_Key1(9), common_term_22);
    // q = 10
    result_by_lpq.set(c_Key3(6, 0, 10), common_term_23);
    result_by_q.set(c_Key1(10), common_term_23);
    // q = 11
    result_by_lpq.set(c_Key3(6, 0, 11), common_term_24);
    result_by_q.set(c_Key1(11), common_term_24);
    // q = 12
    result_by_lpq.set(c_Key3(6, 0, 12), common_term_25);
    result_by_q.set(c_Key1(12), common_term_25);
    // q = 13
    result_by_lpq.set(c_Key3(6, 0, 13), common_term_26);
    result_by_q.set(c_Key1(13), common_term_26);
    // q = 14
    result_by_lpq.set(c_Key3(6, 0, 14), common_term_27);
    result_by_q.set(c_Key1(14), common_term_27);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 0), result_by_q);
    result_by_q.clear();

    // l , p = (6, 1).
    // q = -14
    result_by_lpq.set(c_Key3(6, 1, -14), common_term_28);
    result_by_q.set(c_Key1(-14), common_term_28);
    // q = -13
    result_by_lpq.set(c_Key3(6, 1, -13), common_term_29);
    result_by_q.set(c_Key1(-13), common_term_29);
    // q = -12
    result_by_lpq.set(c_Key3(6, 1, -12), common_term_30);
    result_by_q.set(c_Key1(-12), common_term_30);
    // q = -11
    result_by_lpq.set(c_Key3(6, 1, -11), common_term_31);
    result_by_q.set(c_Key1(-11), common_term_31);
    // q = -10
    result_by_lpq.set(c_Key3(6, 1, -10), common_term_32);
    result_by_q.set(c_Key1(-10), common_term_32);
    // q = -9
    result_by_lpq.set(c_Key3(6, 1, -9), common_term_33);
    result_by_q.set(c_Key1(-9), common_term_33);
    // q = -8
    result_by_lpq.set(c_Key3(6, 1, -8), common_term_34);
    result_by_q.set(c_Key1(-8), common_term_34);
    // q = -7
    result_by_lpq.set(c_Key3(6, 1, -7), common_term_35);
    result_by_q.set(c_Key1(-7), common_term_35);
    // q = -6
    result_by_lpq.set(c_Key3(6, 1, -6), common_term_36);
    result_by_q.set(c_Key1(-6), common_term_36);
    // q = -5
    result_by_lpq.set(c_Key3(6, 1, -5), common_term_37);
    result_by_q.set(c_Key1(-5), common_term_37);
    // q = -4
    result_by_lpq.set(c_Key3(6, 1, -4), common_term_38);
    result_by_q.set(c_Key1(-4), common_term_38);
    // q = -3
    result_by_lpq.set(c_Key3(6, 1, -3), common_term_39);
    result_by_q.set(c_Key1(-3), common_term_39);
    // q = -2
    result_by_lpq.set(c_Key3(6, 1, -2), common_term_40);
    result_by_q.set(c_Key1(-2), common_term_40);
    // q = -1
    result_by_lpq.set(c_Key3(6, 1, -1), common_term_41);
    result_by_q.set(c_Key1(-1), common_term_41);
    // q = 0
    result_by_lpq.set(c_Key3(6, 1, 0), common_term_42);
    result_by_q.set(c_Key1(0), common_term_42);
    // q = 1
    result_by_lpq.set(c_Key3(6, 1, 1), common_term_43);
    result_by_q.set(c_Key1(1), common_term_43);
    // q = 2
    result_by_lpq.set(c_Key3(6, 1, 2), common_term_44);
    result_by_q.set(c_Key1(2), common_term_44);
    // q = 3
    result_by_lpq.set(c_Key3(6, 1, 3), common_term_45);
    result_by_q.set(c_Key1(3), common_term_45);
    // q = 4
    result_by_lpq.set(c_Key3(6, 1, 4), common_term_46);
    result_by_q.set(c_Key1(4), common_term_46);
    // q = 5
    result_by_lpq.set(c_Key3(6, 1, 5), common_term_47);
    result_by_q.set(c_Key1(5), common_term_47);
    // q = 6
    result_by_lpq.set(c_Key3(6, 1, 6), common_term_48);
    result_by_q.set(c_Key1(6), common_term_48);
    // q = 7
    result_by_lpq.set(c_Key3(6, 1, 7), common_term_49);
    result_by_q.set(c_Key1(7), common_term_49);
    // q = 8
    result_by_lpq.set(c_Key3(6, 1, 8), common_term_50);
    result_by_q.set(c_Key1(8), common_term_50);
    // q = 9
    result_by_lpq.set(c_Key3(6, 1, 9), common_term_51);
    result_by_q.set(c_Key1(9), common_term_51);
    // q = 10
    result_by_lpq.set(c_Key3(6, 1, 10), common_term_52);
    result_by_q.set(c_Key1(10), common_term_52);
    // q = 11
    result_by_lpq.set(c_Key3(6, 1, 11), common_term_53);
    result_by_q.set(c_Key1(11), common_term_53);
    // q = 12
    result_by_lpq.set(c_Key3(6, 1, 12), common_term_54);
    result_by_q.set(c_Key1(12), common_term_54);
    // q = 13
    result_by_lpq.set(c_Key3(6, 1, 13), common_term_55);
    result_by_q.set(c_Key1(13), common_term_55);
    // q = 14
    result_by_lpq.set(c_Key3(6, 1, 14), common_term_56);
    result_by_q.set(c_Key1(14), common_term_56);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 1), result_by_q);
    result_by_q.clear();

    // l , p = (6, 2).
    // q = -14
    result_by_lpq.set(c_Key3(6, 2, -14), common_term_57);
    result_by_q.set(c_Key1(-14), common_term_57);
    // q = -13
    result_by_lpq.set(c_Key3(6, 2, -13), common_term_58);
    result_by_q.set(c_Key1(-13), common_term_58);
    // q = -12
    result_by_lpq.set(c_Key3(6, 2, -12), common_term_59);
    result_by_q.set(c_Key1(-12), common_term_59);
    // q = -11
    result_by_lpq.set(c_Key3(6, 2, -11), common_term_60);
    result_by_q.set(c_Key1(-11), common_term_60);
    // q = -10
    result_by_lpq.set(c_Key3(6, 2, -10), common_term_61);
    result_by_q.set(c_Key1(-10), common_term_61);
    // q = -9
    result_by_lpq.set(c_Key3(6, 2, -9), common_term_62);
    result_by_q.set(c_Key1(-9), common_term_62);
    // q = -8
    result_by_lpq.set(c_Key3(6, 2, -8), common_term_63);
    result_by_q.set(c_Key1(-8), common_term_63);
    // q = -7
    result_by_lpq.set(c_Key3(6, 2, -7), common_term_64);
    result_by_q.set(c_Key1(-7), common_term_64);
    // q = -6
    result_by_lpq.set(c_Key3(6, 2, -6), common_term_65);
    result_by_q.set(c_Key1(-6), common_term_65);
    // q = -5
    result_by_lpq.set(c_Key3(6, 2, -5), common_term_66);
    result_by_q.set(c_Key1(-5), common_term_66);
    // q = -4
    result_by_lpq.set(c_Key3(6, 2, -4), common_term_67);
    result_by_q.set(c_Key1(-4), common_term_67);
    // q = -3
    result_by_lpq.set(c_Key3(6, 2, -3), common_term_68);
    result_by_q.set(c_Key1(-3), common_term_68);
    // q = -2
    result_by_lpq.set(c_Key3(6, 2, -2), common_term_69);
    result_by_q.set(c_Key1(-2), common_term_69);
    // q = -1
    result_by_lpq.set(c_Key3(6, 2, -1), common_term_70);
    result_by_q.set(c_Key1(-1), common_term_70);
    // q = 0
    result_by_lpq.set(c_Key3(6, 2, 0), common_term_71);
    result_by_q.set(c_Key1(0), common_term_71);
    // q = 1
    result_by_lpq.set(c_Key3(6, 2, 1), common_term_72);
    result_by_q.set(c_Key1(1), common_term_72);
    // q = 2
    result_by_lpq.set(c_Key3(6, 2, 2), common_term_73);
    result_by_q.set(c_Key1(2), common_term_73);
    // q = 3
    result_by_lpq.set(c_Key3(6, 2, 3), common_term_74);
    result_by_q.set(c_Key1(3), common_term_74);
    // q = 4
    result_by_lpq.set(c_Key3(6, 2, 4), common_term_75);
    result_by_q.set(c_Key1(4), common_term_75);
    // q = 5
    result_by_lpq.set(c_Key3(6, 2, 5), common_term_76);
    result_by_q.set(c_Key1(5), common_term_76);
    // q = 6
    result_by_lpq.set(c_Key3(6, 2, 6), common_term_77);
    result_by_q.set(c_Key1(6), common_term_77);
    // q = 7
    result_by_lpq.set(c_Key3(6, 2, 7), common_term_78);
    result_by_q.set(c_Key1(7), common_term_78);
    // q = 8
    result_by_lpq.set(c_Key3(6, 2, 8), common_term_79);
    result_by_q.set(c_Key1(8), common_term_79);
    // q = 9
    result_by_lpq.set(c_Key3(6, 2, 9), common_term_80);
    result_by_q.set(c_Key1(9), common_term_80);
    // q = 10
    result_by_lpq.set(c_Key3(6, 2, 10), common_term_81);
    result_by_q.set(c_Key1(10), common_term_81);
    // q = 11
    result_by_lpq.set(c_Key3(6, 2, 11), common_term_82);
    result_by_q.set(c_Key1(11), common_term_82);
    // q = 12
    result_by_lpq.set(c_Key3(6, 2, 12), common_term_83);
    result_by_q.set(c_Key1(12), common_term_83);
    // q = 13
    result_by_lpq.set(c_Key3(6, 2, 13), common_term_84);
    result_by_q.set(c_Key1(13), common_term_84);
    // q = 14
    result_by_lpq.set(c_Key3(6, 2, 14), common_term_85);
    result_by_q.set(c_Key1(14), common_term_85);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 2), result_by_q);
    result_by_q.clear();

    // l , p = (6, 3).
    // q = -14
    result_by_lpq.set(c_Key3(6, 3, -14), common_term_86);
    result_by_q.set(c_Key1(-14), common_term_86);
    // q = -13
    result_by_lpq.set(c_Key3(6, 3, -13), common_term_87);
    result_by_q.set(c_Key1(-13), common_term_87);
    // q = -12
    result_by_lpq.set(c_Key3(6, 3, -12), common_term_88);
    result_by_q.set(c_Key1(-12), common_term_88);
    // q = -11
    result_by_lpq.set(c_Key3(6, 3, -11), common_term_89);
    result_by_q.set(c_Key1(-11), common_term_89);
    // q = -10
    result_by_lpq.set(c_Key3(6, 3, -10), common_term_90);
    result_by_q.set(c_Key1(-10), common_term_90);
    // q = -9
    result_by_lpq.set(c_Key3(6, 3, -9), common_term_91);
    result_by_q.set(c_Key1(-9), common_term_91);
    // q = -8
    result_by_lpq.set(c_Key3(6, 3, -8), common_term_92);
    result_by_q.set(c_Key1(-8), common_term_92);
    // q = -7
    result_by_lpq.set(c_Key3(6, 3, -7), common_term_93);
    result_by_q.set(c_Key1(-7), common_term_93);
    // q = -6
    result_by_lpq.set(c_Key3(6, 3, -6), common_term_94);
    result_by_q.set(c_Key1(-6), common_term_94);
    // q = -5
    result_by_lpq.set(c_Key3(6, 3, -5), common_term_95);
    result_by_q.set(c_Key1(-5), common_term_95);
    // q = -4
    result_by_lpq.set(c_Key3(6, 3, -4), common_term_96);
    result_by_q.set(c_Key1(-4), common_term_96);
    // q = -3
    result_by_lpq.set(c_Key3(6, 3, -3), common_term_97);
    result_by_q.set(c_Key1(-3), common_term_97);
    // q = -2
    result_by_lpq.set(c_Key3(6, 3, -2), common_term_98);
    result_by_q.set(c_Key1(-2), common_term_98);
    // q = -1
    result_by_lpq.set(c_Key3(6, 3, -1), common_term_99);
    result_by_q.set(c_Key1(-1), common_term_99);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -5.5)*(1.875*eccentricity_4 + 5.0*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(6, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(6, 3, 1), common_term_99);
    result_by_q.set(c_Key1(1), common_term_99);
    // q = 2
    result_by_lpq.set(c_Key3(6, 3, 2), common_term_98);
    result_by_q.set(c_Key1(2), common_term_98);
    // q = 3
    result_by_lpq.set(c_Key3(6, 3, 3), common_term_97);
    result_by_q.set(c_Key1(3), common_term_97);
    // q = 4
    result_by_lpq.set(c_Key3(6, 3, 4), common_term_96);
    result_by_q.set(c_Key1(4), common_term_96);
    // q = 5
    result_by_lpq.set(c_Key3(6, 3, 5), common_term_95);
    result_by_q.set(c_Key1(5), common_term_95);
    // q = 6
    result_by_lpq.set(c_Key3(6, 3, 6), common_term_94);
    result_by_q.set(c_Key1(6), common_term_94);
    // q = 7
    result_by_lpq.set(c_Key3(6, 3, 7), common_term_93);
    result_by_q.set(c_Key1(7), common_term_93);
    // q = 8
    result_by_lpq.set(c_Key3(6, 3, 8), common_term_92);
    result_by_q.set(c_Key1(8), common_term_92);
    // q = 9
    result_by_lpq.set(c_Key3(6, 3, 9), common_term_91);
    result_by_q.set(c_Key1(9), common_term_91);
    // q = 10
    result_by_lpq.set(c_Key3(6, 3, 10), common_term_90);
    result_by_q.set(c_Key1(10), common_term_90);
    // q = 11
    result_by_lpq.set(c_Key3(6, 3, 11), common_term_89);
    result_by_q.set(c_Key1(11), common_term_89);
    // q = 12
    result_by_lpq.set(c_Key3(6, 3, 12), common_term_88);
    result_by_q.set(c_Key1(12), common_term_88);
    // q = 13
    result_by_lpq.set(c_Key3(6, 3, 13), common_term_87);
    result_by_q.set(c_Key1(13), common_term_87);
    // q = 14
    result_by_lpq.set(c_Key3(6, 3, 14), common_term_86);
    result_by_q.set(c_Key1(14), common_term_86);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 3), result_by_q);
    result_by_q.clear();

    // l , p = (6, 4).
    // q = -14
    result_by_lpq.set(c_Key3(6, 4, -14), common_term_85);
    result_by_q.set(c_Key1(-14), common_term_85);
    // q = -13
    result_by_lpq.set(c_Key3(6, 4, -13), common_term_84);
    result_by_q.set(c_Key1(-13), common_term_84);
    // q = -12
    result_by_lpq.set(c_Key3(6, 4, -12), common_term_83);
    result_by_q.set(c_Key1(-12), common_term_83);
    // q = -11
    result_by_lpq.set(c_Key3(6, 4, -11), common_term_82);
    result_by_q.set(c_Key1(-11), common_term_82);
    // q = -10
    result_by_lpq.set(c_Key3(6, 4, -10), common_term_81);
    result_by_q.set(c_Key1(-10), common_term_81);
    // q = -9
    result_by_lpq.set(c_Key3(6, 4, -9), common_term_80);
    result_by_q.set(c_Key1(-9), common_term_80);
    // q = -8
    result_by_lpq.set(c_Key3(6, 4, -8), common_term_79);
    result_by_q.set(c_Key1(-8), common_term_79);
    // q = -7
    result_by_lpq.set(c_Key3(6, 4, -7), common_term_78);
    result_by_q.set(c_Key1(-7), common_term_78);
    // q = -6
    result_by_lpq.set(c_Key3(6, 4, -6), common_term_77);
    result_by_q.set(c_Key1(-6), common_term_77);
    // q = -5
    result_by_lpq.set(c_Key3(6, 4, -5), common_term_76);
    result_by_q.set(c_Key1(-5), common_term_76);
    // q = -4
    result_by_lpq.set(c_Key3(6, 4, -4), common_term_75);
    result_by_q.set(c_Key1(-4), common_term_75);
    // q = -3
    result_by_lpq.set(c_Key3(6, 4, -3), common_term_74);
    result_by_q.set(c_Key1(-3), common_term_74);
    // q = -2
    result_by_lpq.set(c_Key3(6, 4, -2), common_term_73);
    result_by_q.set(c_Key1(-2), common_term_73);
    // q = -1
    result_by_lpq.set(c_Key3(6, 4, -1), common_term_72);
    result_by_q.set(c_Key1(-1), common_term_72);
    // q = 0
    result_by_lpq.set(c_Key3(6, 4, 0), common_term_71);
    result_by_q.set(c_Key1(0), common_term_71);
    // q = 1
    result_by_lpq.set(c_Key3(6, 4, 1), common_term_70);
    result_by_q.set(c_Key1(1), common_term_70);
    // q = 2
    result_by_lpq.set(c_Key3(6, 4, 2), common_term_69);
    result_by_q.set(c_Key1(2), common_term_69);
    // q = 3
    result_by_lpq.set(c_Key3(6, 4, 3), common_term_68);
    result_by_q.set(c_Key1(3), common_term_68);
    // q = 4
    result_by_lpq.set(c_Key3(6, 4, 4), common_term_67);
    result_by_q.set(c_Key1(4), common_term_67);
    // q = 5
    result_by_lpq.set(c_Key3(6, 4, 5), common_term_66);
    result_by_q.set(c_Key1(5), common_term_66);
    // q = 6
    result_by_lpq.set(c_Key3(6, 4, 6), common_term_65);
    result_by_q.set(c_Key1(6), common_term_65);
    // q = 7
    result_by_lpq.set(c_Key3(6, 4, 7), common_term_64);
    result_by_q.set(c_Key1(7), common_term_64);
    // q = 8
    result_by_lpq.set(c_Key3(6, 4, 8), common_term_63);
    result_by_q.set(c_Key1(8), common_term_63);
    // q = 9
    result_by_lpq.set(c_Key3(6, 4, 9), common_term_62);
    result_by_q.set(c_Key1(9), common_term_62);
    // q = 10
    result_by_lpq.set(c_Key3(6, 4, 10), common_term_61);
    result_by_q.set(c_Key1(10), common_term_61);
    // q = 11
    result_by_lpq.set(c_Key3(6, 4, 11), common_term_60);
    result_by_q.set(c_Key1(11), common_term_60);
    // q = 12
    result_by_lpq.set(c_Key3(6, 4, 12), common_term_59);
    result_by_q.set(c_Key1(12), common_term_59);
    // q = 13
    result_by_lpq.set(c_Key3(6, 4, 13), common_term_58);
    result_by_q.set(c_Key1(13), common_term_58);
    // q = 14
    result_by_lpq.set(c_Key3(6, 4, 14), common_term_57);
    result_by_q.set(c_Key1(14), common_term_57);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 4), result_by_q);
    result_by_q.clear();

    // l , p = (6, 5).
    // q = -14
    result_by_lpq.set(c_Key3(6, 5, -14), common_term_56);
    result_by_q.set(c_Key1(-14), common_term_56);
    // q = -13
    result_by_lpq.set(c_Key3(6, 5, -13), common_term_55);
    result_by_q.set(c_Key1(-13), common_term_55);
    // q = -12
    result_by_lpq.set(c_Key3(6, 5, -12), common_term_54);
    result_by_q.set(c_Key1(-12), common_term_54);
    // q = -11
    result_by_lpq.set(c_Key3(6, 5, -11), common_term_53);
    result_by_q.set(c_Key1(-11), common_term_53);
    // q = -10
    result_by_lpq.set(c_Key3(6, 5, -10), common_term_52);
    result_by_q.set(c_Key1(-10), common_term_52);
    // q = -9
    result_by_lpq.set(c_Key3(6, 5, -9), common_term_51);
    result_by_q.set(c_Key1(-9), common_term_51);
    // q = -8
    result_by_lpq.set(c_Key3(6, 5, -8), common_term_50);
    result_by_q.set(c_Key1(-8), common_term_50);
    // q = -7
    result_by_lpq.set(c_Key3(6, 5, -7), common_term_49);
    result_by_q.set(c_Key1(-7), common_term_49);
    // q = -6
    result_by_lpq.set(c_Key3(6, 5, -6), common_term_48);
    result_by_q.set(c_Key1(-6), common_term_48);
    // q = -5
    result_by_lpq.set(c_Key3(6, 5, -5), common_term_47);
    result_by_q.set(c_Key1(-5), common_term_47);
    // q = -4
    result_by_lpq.set(c_Key3(6, 5, -4), common_term_46);
    result_by_q.set(c_Key1(-4), common_term_46);
    // q = -3
    result_by_lpq.set(c_Key3(6, 5, -3), common_term_45);
    result_by_q.set(c_Key1(-3), common_term_45);
    // q = -2
    result_by_lpq.set(c_Key3(6, 5, -2), common_term_44);
    result_by_q.set(c_Key1(-2), common_term_44);
    // q = -1
    result_by_lpq.set(c_Key3(6, 5, -1), common_term_43);
    result_by_q.set(c_Key1(-1), common_term_43);
    // q = 0
    result_by_lpq.set(c_Key3(6, 5, 0), common_term_42);
    result_by_q.set(c_Key1(0), common_term_42);
    // q = 1
    result_by_lpq.set(c_Key3(6, 5, 1), common_term_41);
    result_by_q.set(c_Key1(1), common_term_41);
    // q = 2
    result_by_lpq.set(c_Key3(6, 5, 2), common_term_40);
    result_by_q.set(c_Key1(2), common_term_40);
    // q = 3
    result_by_lpq.set(c_Key3(6, 5, 3), common_term_39);
    result_by_q.set(c_Key1(3), common_term_39);
    // q = 4
    result_by_lpq.set(c_Key3(6, 5, 4), common_term_38);
    result_by_q.set(c_Key1(4), common_term_38);
    // q = 5
    result_by_lpq.set(c_Key3(6, 5, 5), common_term_37);
    result_by_q.set(c_Key1(5), common_term_37);
    // q = 6
    result_by_lpq.set(c_Key3(6, 5, 6), common_term_36);
    result_by_q.set(c_Key1(6), common_term_36);
    // q = 7
    result_by_lpq.set(c_Key3(6, 5, 7), common_term_35);
    result_by_q.set(c_Key1(7), common_term_35);
    // q = 8
    result_by_lpq.set(c_Key3(6, 5, 8), common_term_34);
    result_by_q.set(c_Key1(8), common_term_34);
    // q = 9
    result_by_lpq.set(c_Key3(6, 5, 9), common_term_33);
    result_by_q.set(c_Key1(9), common_term_33);
    // q = 10
    result_by_lpq.set(c_Key3(6, 5, 10), common_term_32);
    result_by_q.set(c_Key1(10), common_term_32);
    // q = 11
    result_by_lpq.set(c_Key3(6, 5, 11), common_term_31);
    result_by_q.set(c_Key1(11), common_term_31);
    // q = 12
    result_by_lpq.set(c_Key3(6, 5, 12), common_term_30);
    result_by_q.set(c_Key1(12), common_term_30);
    // q = 13
    result_by_lpq.set(c_Key3(6, 5, 13), common_term_29);
    result_by_q.set(c_Key1(13), common_term_29);
    // q = 14
    result_by_lpq.set(c_Key3(6, 5, 14), common_term_28);
    result_by_q.set(c_Key1(14), common_term_28);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 5), result_by_q);
    result_by_q.clear();

    // l , p = (6, 6).
    // q = -14
    result_by_lpq.set(c_Key3(6, 6, -14), common_term_27);
    result_by_q.set(c_Key1(-14), common_term_27);
    // q = -13
    result_by_lpq.set(c_Key3(6, 6, -13), common_term_26);
    result_by_q.set(c_Key1(-13), common_term_26);
    // q = -12
    result_by_lpq.set(c_Key3(6, 6, -12), common_term_25);
    result_by_q.set(c_Key1(-12), common_term_25);
    // q = -11
    result_by_lpq.set(c_Key3(6, 6, -11), common_term_24);
    result_by_q.set(c_Key1(-11), common_term_24);
    // q = -10
    result_by_lpq.set(c_Key3(6, 6, -10), common_term_23);
    result_by_q.set(c_Key1(-10), common_term_23);
    // q = -9
    result_by_lpq.set(c_Key3(6, 6, -9), common_term_22);
    result_by_q.set(c_Key1(-9), common_term_22);
    // q = -8
    result_by_lpq.set(c_Key3(6, 6, -8), common_term_21);
    result_by_q.set(c_Key1(-8), common_term_21);
    // q = -7
    result_by_lpq.set(c_Key3(6, 6, -7), common_term_20);
    result_by_q.set(c_Key1(-7), common_term_20);
    // q = -6
    result_by_lpq.set(c_Key3(6, 6, -6), common_term_19);
    result_by_q.set(c_Key1(-6), common_term_19);
    // q = -5
    result_by_lpq.set(c_Key3(6, 6, -5), common_term_18);
    result_by_q.set(c_Key1(-5), common_term_18);
    // q = -4
    result_by_lpq.set(c_Key3(6, 6, -4), common_term_17);
    result_by_q.set(c_Key1(-4), common_term_17);
    // q = -3
    result_by_lpq.set(c_Key3(6, 6, -3), common_term_16);
    result_by_q.set(c_Key1(-3), common_term_16);
    // q = -2
    result_by_lpq.set(c_Key3(6, 6, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(6, 6, -1), common_term_14);
    result_by_q.set(c_Key1(-1), common_term_14);
    // q = 0
    result_by_lpq.set(c_Key3(6, 6, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(6, 6, 1), common_term_12);
    result_by_q.set(c_Key1(1), common_term_12);
    // q = 2
    result_by_lpq.set(c_Key3(6, 6, 2), common_term_11);
    result_by_q.set(c_Key1(2), common_term_11);
    // q = 3
    result_by_lpq.set(c_Key3(6, 6, 3), common_term_10);
    result_by_q.set(c_Key1(3), common_term_10);
    // q = 4
    result_by_lpq.set(c_Key3(6, 6, 4), common_term_9);
    result_by_q.set(c_Key1(4), common_term_9);
    // q = 5
    result_by_lpq.set(c_Key3(6, 6, 5), common_term_8);
    result_by_q.set(c_Key1(5), common_term_8);
    // q = 7
    result_by_lpq.set(c_Key3(6, 6, 7), common_term_7);
    result_by_q.set(c_Key1(7), common_term_7);
    // q = 8
    result_by_lpq.set(c_Key3(6, 6, 8), common_term_6);
    result_by_q.set(c_Key1(8), common_term_6);
    // q = 9
    result_by_lpq.set(c_Key3(6, 6, 9), common_term_5);
    result_by_q.set(c_Key1(9), common_term_5);
    // q = 10
    result_by_lpq.set(c_Key3(6, 6, 10), common_term_4);
    result_by_q.set(c_Key1(10), common_term_4);
    // q = 11
    result_by_lpq.set(c_Key3(6, 6, 11), common_term_3);
    result_by_q.set(c_Key1(11), common_term_3);
    // q = 12
    result_by_lpq.set(c_Key3(6, 6, 12), common_term_2);
    result_by_q.set(c_Key1(12), common_term_2);
    // q = 13
    result_by_lpq.set(c_Key3(6, 6, 13), common_term_1);
    result_by_q.set(c_Key1(13), common_term_1);
    // q = 14
    result_by_lpq.set(c_Key3(6, 6, 14), common_term_0);
    result_by_q.set(c_Key1(14), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 6), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l6_e20(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 6.
    // Functions are taylor expanded and truncated at eccentricity^20.

    // Functions Calculated for l = 6.

    c_IntMap<c_Key3, double> result_by_lpq(271);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(7);
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
    double common_term_0 = 0.022922350820189795*eccentricity_19;
    double common_term_1 = 0.015862859829246384*eccentricity_18;
    double common_term_2 = 0.017918913425916616*eccentricity_19 + 0.01084169551820165*eccentricity_17;
    double common_term_3 = 0.013084327126699582*eccentricity_18 + 0.0072929036443899312*eccentricity_16;
    double common_term_4 = 0.014115067401750228*eccentricity_19 + 0.0092344967537439017*eccentricity_17 + 0.0048049414003220301*eccentricity_15;
    double common_term_5 = 0.0096095290169364243*eccentricity_18 + 0.0062609481445460281*eccentricity_16 + 0.0030791548251865712*eccentricity_14;
    double common_term_6 = 0.0081255707234367575*eccentricity_19 + 0.0062104682065225253*eccentricity_17 + 0.0040361118640668849*eccentricity_15 + 0.001899346759560887*eccentricity_13;
    double common_term_7 = 0.0048572361302202262*eccentricity_18 + 0.00373762502843139*eccentricity_16 + 0.00243231377997003*eccentricity_14 + 0.001109476461038961*eccentricity_12;
    double common_term_8 = 0.0030445391448489479*eccentricity_19 + 0.0026072237358337232*eccentricity_17 + 0.0020303272818552403*eccentricity_15 + 0.0013314563022267142*eccentricity_13 + 0.00059728880847553536*eccentricity_11;
    double common_term_9 = 0.0013621491432337993*eccentricity_18 + 0.0011874065577769281*eccentricity_16 + 0.00094222649778205334*eccentricity_14 + 0.00062850729517396184*eccentricity_12 + 0.00028218694885361552*eccentricity_10;
    double common_term_10 = 0.00048813635097630242*eccentricity_19 + 0.00045961818506161531*eccentricity_17 + 0.00041054413690195455*eccentricity_15 + 0.0003350640581799792*eccentricity_13 + 0.00023041861397879464*eccentricity_11 + 0.00010593959263392857*eccentricity_9;
    double common_term_11 = 9.2860865109309843e-5*eccentricity_18 + 8.9819137173376178e-5*eccentricity_16 + 8.2929263701833146e-5*eccentricity_14 + 7.0512290564373898e-5*eccentricity_12 + 5.0981040564373898e-5*eccentricity_10 + 2.4801587301587302e-5*eccentricity_8;
    double common_term_12 = 4.0924552647722061e-6*eccentricity_19 + 4.1848215867039857e-6*eccentricity_17 + 4.1909728094319812e-6*eccentricity_15 + 4.0465337539125652e-6*eccentricity_13 + 3.6498646677276235e-6*eccentricity_11 + 2.8579954117063492e-6*eccentricity_9 + 1.5500992063492063e-6*eccentricity_7;
    double common_term_13 = -0.00013739526612460433*eccentricity_19 - 0.0001515485574898391*eccentricity_17 - 0.00016791066998667327*eccentricity_15 - 0.00018671935337701312*eccentricity_13 - 0.0002079050376932457*eccentricity_11 - 0.00023038349454365079*eccentricity_9 - 0.00024956597222222222*eccentricity_7 - 0.00026041666666666667*eccentricity_5;
    double common_term_14 = 0.0018851574702340193*eccentricity_18 + 0.0022509684771059916*eccentricity_16 + 0.0027452027300117578*eccentricity_14 + 0.0034391534391534392*eccentricity_12 + 0.0044436177248677249*eccentricity_10 + 0.0064236111111111111*eccentricity_8 + 0.0041666666666666667*eccentricity_6 + 0.041666666666666667*eccentricity_4;
    double common_term_15 = 0.012662847740569624*eccentricity_19 + 0.014886944759561091*eccentricity_17 + 0.017825560314314706*eccentricity_15 + 0.022001722199576242*eccentricity_13 + 0.025369589669363839*eccentricity_11 + 0.06251220703125*eccentricity_9 - 0.13447265625*eccentricity_7 + 0.73828125*eccentricity_5 - 0.5625*eccentricity_3;
    double common_term_16 = 0.048817643281110026*eccentricity_18 + 0.057518156117899835*eccentricity_16 + 0.077645640432098765*eccentricity_14 + 0.011097056878306878*eccentricity_12 + 0.69415509259259259*eccentricity_10 - 2.6513888888888889*eccentricity_8 + 7.7083333333333333*eccentricity_6 - 7.6666666666666667*eccentricity_4 + 2.0*eccentricity_2;
    double common_term_17 = 0.11870338844372974*eccentricity_19 + 0.12349315006894169*eccentricity_17 + 0.32680401728251744*eccentricity_15 - 1.0275152036542405*eccentricity_13 + 7.129992025869864*eccentricity_11 - 25.022871229383681*eccentricity_9 + 55.931260850694444*eccentricity_7 - 57.669270833333333*eccentricity_5 + 22.8125*eccentricity_3 - 2.5*eccentricity;
    double common_term_18 = 0.014733101981026786*eccentricity_18 + 2.4089809121890944*eccentricity_16 - 12.093317721619898*eccentricity_14 + 55.52111328125*eccentricity_12 - 165.161953125*eccentricity_10 + 312.66796875*eccentricity_8 - 314.96875*eccentricity_6 + 145.6875*eccentricity_4 - 25.5*eccentricity_2 + 1.0;
    double common_term_19 = -2.4799441857411133*eccentricity_19 + 19.458961742712141*eccentricity_17 - 91.923493899644157*eccentricity_15 + 338.7226141924917*eccentricity_13 - 864.59438154997649*eccentricity_11 + 1443.2193650987413*eccentricity_9 - 1393.6566297743056*eccentricity_7 + 687.33072916666667*eccentricity_5 - 147.4375*eccentricity_3 + 9.5*eccentricity;
    double common_term_20 = 133.27556538528789*eccentricity_18 - 546.74611371155973*eccentricity_16 + 1711.4160550319665*eccentricity_14 - 3828.2373974867725*eccentricity_12 + 5759.9242766203704*eccentricity_10 - 5298.1472222222222*eccentricity_8 + 2669.9583333333333*eccentricity_6 - 635.66666666666667*eccentricity_4 + 51.5*eccentricity_2;
    double common_term_21 = 764.31647945043453*eccentricity_19 - 2730.4086346368464*eccentricity_17 + 7471.4286582968065*eccentricity_15 - 14912.162450819697*eccentricity_13 + 20512.371441977365*eccentricity_11 - 17943.421472167969*eccentricity_9 + 9036.48076171875*eccentricity_7 - 2273.73046875*eccentricity_5 + 209.1875*eccentricity_3;
    double common_term_22 = -11915.753622214634*eccentricity_18 + 29046.852838978589*eccentricity_16 - 52466.461716998962*eccentricity_14 + 66651.739275896991*eccentricity_12 - 55466.796926669974*eccentricity_10 + 27569.268663194444*eccentricity_8 - 7131.3020833333333*eccentricity_6 + 707.60416666666667*eccentricity_4;
    double common_term_23 = -46673.034149869906*eccentricity_19 + 102776.8674300622*eccentricity_17 - 169884.58606799649*eccentricity_15 + 200856.68874116792*eccentricity_13 - 159206.19147035912*eccentricity_11 + 77533.185511804006*eccentricity_9 - 20271.980327690972*eccentricity_7 + 2105.3309895833333*eccentricity_5;
    double common_term_24 = 336348.99577181095*eccentricity_18 - 513305.51556331169*eccentricity_16 + 568293.86935546875*eccentricity_14 - 429709.44654017857*eccentricity_12 + 204154.14174107143*eccentricity_10 - 53371.869642857143*eccentricity_8 + 5692.6*eccentricity_6;
    double common_term_25 = 1030625.0000255625*eccentricity_19 - 1462655.5248712089*eccentricity_17 + 1524005.4527814355*eccentricity_15 - 1101154.6901831306*eccentricity_13 + 509067.61315316565*eccentricity_11 - 132109.62696073017*eccentricity_9 + 14289.359655567956*eccentricity_7;
    double common_term_26 = -3963205.0956288945*eccentricity_18 + 3902740.1575198022*eccentricity_16 - 2699160.2059729639*eccentricity_14 + 1212454.643498505*eccentricity_12 - 310806.6863908179*eccentricity_10 + 33793.725173611111*eccentricity_8;
    double common_term_27 = -10279174.999861816*eccentricity_19 + 9601164.9413476062*eccentricity_17 - 6366468.7414403115*eccentricity_15 + 2776620.1823720529*eccentricity_13 - 700696.61881991795*eccentricity_11 + 76106.604036603655*eccentricity_9;
    double common_term_28 = 22801809.408598651*eccentricity_18 - 14519357.447903737*eccentricity_16 + 6146345.1867708125*eccentricity_14 - 1523325.8970551748*eccentricity_12 + 164531.88032104277*eccentricity_10;
    double common_term_29 = 52487593.808247077*eccentricity_19 - 32143278.689166799*eccentricity_17 + 13207295.111706393*eccentricity_15 - 3209522.8770327859*eccentricity_13 + 343556.26409161624*eccentricity_11;
    double common_term_30 = -69303174.219026244*eccentricity_18 + 27645292.495112135*eccentricity_16 - 6579774.7889591073*eccentricity_14 + 696265.66076045049*eccentricity_12;
    double common_term_31 = -145926951.25053421*eccentricity_19 + 56532270.950514135*eccentricity_17 - 13168134.28329271*eccentricity_15 + 1374907.1380696101*eccentricity_13;
    double common_term_32 = 113213027.77607351*eccentricity_18 - 25796204.09783613*eccentricity_16 + 2653822.990858011*eccentricity_14;
    double common_term_33 = 222493980.50529293*eccentricity_19 - 49578212.632336916*eccentricity_17 + 5020071.1854010095*eccentricity_15;
    double common_term_34 = -93662035.610407658*eccentricity_18 + 9326944.8146272347*eccentricity_16;
    double common_term_35 = -174215515.00922377*eccentricity_19 + 17051591.869169781*eccentricity_17;
    double common_term_36 = 30723607.112370791*eccentricity_18;
    double common_term_37 = 54632452.945821325*eccentricity_19;
    double common_term_38 = 32.499502814864729*eccentricity_19;
    double common_term_39 = 23.806769508069648*eccentricity_18;
    double common_term_40 = 37.934044841615997*eccentricity_19 + 17.44027503927414*eccentricity_17;
    double common_term_41 = 31.087069665787842*eccentricity_18 + 12.777465638626859*eccentricity_16;
    double common_term_42 = 56.700521426502138*eccentricity_19 + 25.193606264517605*eccentricity_17 + 9.3623755403439596*eccentricity_15;
    double common_term_43 = 46.747153665108266*eccentricity_18 + 20.23311907505449*eccentricity_16 + 6.8610594833284023*eccentricity_14;
    double common_term_44 = 76.654507201453582*eccentricity_19 + 38.404338917488468*eccentricity_17 + 16.127849449382982*eccentricity_15 + 5.0289451730453766*eccentricity_13;
    double common_term_45 = 64.258903808424788*eccentricity_18 + 31.433252209402247*eccentricity_16 + 12.774890073295564*eccentricity_14 + 3.6869184502723999*eccentricity_12;
    double common_term_46 = 99.581952046703531*eccentricity_19 + 53.687773503069497*eccentricity_17 + 25.631257987134666*eccentricity_15 + 10.065197431835127*eccentricity_13 + 2.7038058547888005*eccentricity_11;
    double common_term_47 = 84.482214257104328*eccentricity_18 + 44.708165903920298*eccentricity_16 + 20.823311339539367*eccentricity_14 + 7.8942329545454545*eccentricity_12 + 1.9835435267857143*eccentricity_10;
    double common_term_48 = 125.47682746354936*eccentricity_19 + 71.465923451466082*eccentricity_17 + 37.110228392040713*eccentricity_15 + 16.857054417993809*eccentricity_13 + 6.1673492182709763*eccentricity_11 + 1.4557787644589809*eccentricity_9;
    double common_term_49 = 107.47290128583297*eccentricity_18 + 60.282476891289351*eccentricity_16 + 30.706129331300772*eccentricity_14 + 13.599795869157848*eccentricity_12 + 4.8019696593915344*eccentricity_10 + 1.0689856150793651*eccentricity_8;
    double common_term_50 = 154.42704547593751*eccentricity_19 + 91.816699185913163*eccentricity_17 + 50.705265013281401*eccentricity_15 + 25.328705772672381*eccentricity_13 + 10.936387797764369*eccentricity_11 + 3.7278908865792411*eccentricity_9 + 0.78542131696428571*eccentricity_7;
    double common_term_51 = 133.31718833349601*eccentricity_18 + 78.241000453965081*eccentricity_16 + 42.530252436526308*eccentricity_14 + 20.830027373603762*eccentricity_12 + 8.7675006200396825*eccentricity_10 + 2.8865575396825397*eccentricity_8 + 0.57743055555555556*eccentricity_6;
    double common_term_52 = 186.51772998109614*eccentricity_19 + 114.82607533127794*eccentricity_17 + 66.503010614802668*eccentricity_15 + 35.574591930202706*eccentricity_13 + 17.079967815661557*eccentricity_11 + 7.0080796983506944*eccentricity_9 + 2.2298502604166667*eccentricity_7 + 0.42473958333333333*eccentricity_5;
    double common_term_53 = 0.3125*eccentricity_4*std::pow(1.0 - eccentricity_2, -5.5);
    double common_term_54 = 221.83681072248169*eccentricity_19 + 140.58189417757143*eccentricity_17 + 84.591396417704583*eccentricity_15 + 47.683038475374696*eccentricity_13 + 24.687632868529628*eccentricity_11 + 11.385783329716435*eccentricity_9 + 4.4400716145833333*eccentricity_7 + 1.3216145833333333*eccentricity_5 + 0.22916666666666667*eccentricity_3;
    double common_term_55 = 193.91624807174203*eccentricity_18 + 121.66162162931338*eccentricity_16 + 72.355813095927028*eccentricity_14 + 40.227262318121693*eccentricity_12 + 20.4865234375*eccentricity_10 + 9.2607638888888889*eccentricity_8 + 3.5260416666666667*eccentricity_6 + eccentricity_4 + 0.25*eccentricity_2;
    double common_term_56 = 260.48306421527342*eccentricity_19 + 169.18394301255405*eccentricity_17 + 105.07152791990309*eccentricity_15 + 61.756888719286237*eccentricity_13 + 33.865264434814453*eccentricity_11 + 16.961810302734375*eccentricity_9 + 7.60693359375*eccentricity_7 + 2.3828125*eccentricity_5 + 1.9375*eccentricity_3 - 0.5*eccentricity;
    double common_term_57 = 228.87702869884727*eccentricity_18 + 147.33348536100457*eccentricity_16 + 90.570352541769022*eccentricity_14 + 52.633672357253086*eccentricity_12 + 28.299939236111111*eccentricity_10 + 15.146267361111111*eccentricity_8 + 1.8888888888888889*eccentricity_6 + 11.0*eccentricity_4 - 5.5*eccentricity_2 + 1.0;
    double common_term_58 = 302.58869893316846*eccentricity_19 + 200.77331408659388*eccentricity_17 + 128.05906230778546*eccentricity_15 + 78.23271772848866*eccentricity_13 + 42.88214436283818*eccentricity_11 + 32.950558132595486*eccentricity_9 - 15.68603515625*eccentricity_7 + 51.549479166666667*eccentricity_5 - 31.5625*eccentricity_3 + 7.5*eccentricity;
    double common_term_59 = 267.18237103956573*eccentricity_18 + 175.46965506417411*eccentricity_16 + 113.89916190011161*eccentricity_14 + 52.828828125*eccentricity_12 + 91.6787109375*eccentricity_10 - 113.184375*eccentricity_8 + 204.765625*eccentricity_6 - 131.0*eccentricity_4 + 33.25*eccentricity_2;
    double common_term_60 = 348.78616394559121*eccentricity_19 + 232.20966169124055*eccentricity_17 + 172.6497450068172*eccentricity_15 + 15.369647119486773*eccentricity_13 + 311.31128373322663*eccentricity_11 - 504.1496912073206*eccentricity_9 + 708.82418619791667*eccentricity_7 - 445.40755208333333*eccentricity_5 + 113.39583333333333*eccentricity_3;
    double common_term_61 = 286.79337526104925*eccentricity_18 + 312.12277597607451*eccentricity_16 - 246.34788596034961*eccentricity_14 + 1101.6653302228009*eccentricity_12 - 1826.6253100198413*eccentricity_10 + 2195.2217013888889*eccentricity_8 - 1320.9854166666667*eccentricity_6 + 328.6875*eccentricity_4;
    double common_term_62 = 276.91985403468286*eccentricity_19 + 764.30016386689826*eccentricity_17 - 1359.3114061362914*eccentricity_15 + 3718.3121417897088*eccentricity_13 - 5816.4672563280378*eccentricity_11 + 6214.475967843192*eccentricity_9 - 3544.72041015625*eccentricity_7 + 851.72890625*eccentricity_5;
    double common_term_63 = 2364.3348192932344*eccentricity_18 - 5315.6236015997911*eccentricity_16 + 11691.499287349066*eccentricity_14 - 16867.581598370444*eccentricity_12 + 16351.393849206349*eccentricity_10 - 8809.5994543650794*eccentricity_8 + 2032.5086805555556*eccentricity_6;
    double common_term_64 = 7843.1667342389274*eccentricity_19 - 17930.403401144712*eccentricity_17 + 34286.309547444649*eccentricity_15 - 45480.198617447125*eccentricity_13 + 40502.908957301414*eccentricity_11 - 20602.832414947994*eccentricity_9 + 4552.519982328869*eccentricity_7;
    double common_term_65 = -55058.978267031674*eccentricity_18 + 94485.756624026608*eccentricity_16 - 115592.32271763393*eccentricity_14 + 95387.253716517857*eccentricity_12 - 45858.797209821429*eccentricity_10 + 9696.6136160714286*eccentricity_8;
    double common_term_66 = -157563.67113880886*eccentricity_19 + 246684.59148802252*eccentricity_17 - 279670.83386660374*eccentricity_15 + 215239.10583281955*eccentricity_13 - 97969.551795844964*eccentricity_11 + 19824.016237197688*eccentricity_9;
    double common_term_67 = 614599.52279763078*eccentricity_18 - 648932.12729824959*eccentricity_16 + 468205.26353690259*eccentricity_14 - 202167.86227439324*eccentricity_12 + 39172.045639467593*eccentricity_10;
    double common_term_68 = 1470169.4314044627*eccentricity_19 - 1452441.8929208383*eccentricity_17 + 986666.93010043558*eccentricity_15 - 405000.10163826602*eccentricity_13 + 75209.476258785694*eccentricity_11;
    double common_term_69 = -3150324.4340984275*eccentricity_18 + 2022374.0993700757*eccentricity_16 - 790764.76856383999*eccentricity_14 + 140889.14879560632*eccentricity_12;
    double common_term_70 = -6646703.9632778247*eccentricity_19 + 4045201.946674794*eccentricity_17 - 1509687.3642940114*eccentricity_15 + 258358.9854943355*eccentricity_13;
    double common_term_71 = 7917691.3983317176*eccentricity_18 - 2825679.5997690245*eccentricity_16 + 465019.77231634019*eccentricity_14;
    double common_term_72 = 15199844.212386322*eccentricity_19 - 5196494.4906410846*eccentricity_17 + 823332.45820022203*eccentricity_15;
    double common_term_73 = -9407045.9366181644*eccentricity_18 + 1436577.5613074477*eccentricity_16;
    double common_term_74 = -16789344.846203008*eccentricity_19 + 2474007.3651317233*eccentricity_17;
    double common_term_75 = 4210740.050523545*eccentricity_18;
    double common_term_76 = 7090675.2633715155*eccentricity_19;
    double common_term_77 = 1970.0581517951859*eccentricity_19;
    double common_term_78 = 1379.9519990815706*eccentricity_18;
    double common_term_79 = 682.92324856619946*eccentricity_19 + 964.55733379005519*eccentricity_17;
    double common_term_80 = 667.37643511432351*eccentricity_18 + 672.63959971401422*eccentricity_16;
    double common_term_81 = 1425.3455366269372*eccentricity_19 + 598.58843041142969*eccentricity_17 + 467.86907427341745*eccentricity_15;
    double common_term_82 = 1140.51596159684*eccentricity_18 + 509.53761141174674*eccentricity_16 + 324.51430162471457*eccentricity_14;
    double common_term_83 = 1620.8799199100215*eccentricity_19 + 914.54633828437*eccentricity_17 + 418.46794526810884*eccentricity_15 + 224.37256107712711*eccentricity_13;
    double common_term_84 = 1340.8271825096227*eccentricity_18 + 732.13784381473508*eccentricity_16 + 334.65154916808291*eccentricity_14 + 154.58385662229521*eccentricity_12;
    double common_term_85 = 1900.9690804746948*eccentricity_19 + 1101.6096257416986*eccentricity_17 + 583.59585452682548*eccentricity_15 + 262.05906939797587*eccentricity_13 + 106.07548423172592*eccentricity_11;
    double common_term_86 = 1583.9170905567324*eccentricity_18 + 899.024195717758*eccentricity_16 + 462.35437831366325*eccentricity_14 + 201.66168540564374*eccentricity_12 + 72.456307457010582*eccentricity_10;
    double common_term_87 = 2190.7161050994606*eccentricity_19 + 1312.4932872677859*eccentricity_17 + 728.71253177272477*eccentricity_15 + 363.61292241894971*eccentricity_13 + 152.84656421708472*eccentricity_11 + 49.231494630413291*eccentricity_9;
    double common_term_88 = 1838.2648011325223*eccentricity_18 + 1081.2914443058163*eccentricity_16 + 586.50778618354302*eccentricity_14 + 283.60770228794643*eccentricity_12 + 114.264453125*eccentricity_10 + 33.245424107142857*eccentricity_8;
    double common_term_89 = 2493.7485403906004*eccentricity_19 + 1534.5999084893122*eccentricity_17 + 885.39092488874398*eccentricity_15 + 468.57150143626716*eccentricity_13 + 219.23386247490239*eccentricity_11 + 84.318628462534102*eccentricity_9 + 22.286900111607143*eccentricity_7;
    double common_term_90 = 2104.9393267147584*eccentricity_18 + 1274.1749809270992*eccentricity_16 + 720.3200735550779*eccentricity_14 + 371.4345467739565*eccentricity_12 + 167.85452628968254*eccentricity_10 + 61.431349206349206*eccentricity_8 + 14.809722222222222*eccentricity_6;
    double common_term_91 = 2809.8535036986827*eccentricity_19 + 1768.1347850554227*eccentricity_17 + 1051.9076581006561*eccentricity_15 + 582.03084909234728*eccentricity_13 + 291.99724464416504*eccentricity_11 + 127.20421840122768*eccentricity_9 + 44.17548828125*eccentricity_7 + 9.73515625*eccentricity_5;
    double common_term_92 = 2383.7840763945609*eccentricity_18 + 1477.6188254402425*eccentricity_16 + 863.15967910075875*eccentricity_14 + 466.87657231729497*eccentricity_12 + 227.51363467261905*eccentricity_10 + 95.336545138888889*eccentricity_8 + 31.327083333333333*eccentricity_6 + 6.3125*eccentricity_4;
    double common_term_93 = 3138.8980659529389*eccentricity_19 + 2012.9583797879042*eccentricity_17 + 1228.1407056809174*eccentricity_15 + 703.71636338630051*eccentricity_13 + 371.5901080459514*eccentricity_11 + 175.57016556351273*eccentricity_9 + 70.590397135416667*eccentricity_7 + 21.873697916666667*eccentricity_5 + 4.0208333333333333*eccentricity_3;
    double common_term_94 = std::pow(1.0 - eccentricity_2, -5.5)*(1.25*eccentricity_4 + 2.5*eccentricity_2);
    double common_term_95 = 3480.7616618177654*eccentricity_19 + 2268.94845127767*eccentricity_17 + 1413.9655824127353*eccentricity_15 + 833.50103988374983*eccentricity_13 + 457.87624250058775*eccentricity_11 + 229.31607733832465*eccentricity_9 + 101.17106119791667*eccentricity_7 + 37.080729166666667*eccentricity_5 + 10.0625*eccentricity_3 + 1.5*eccentricity;
    double common_term_96 = 2977.4852960165959*eccentricity_18 + 1915.6771936070526*eccentricity_16 + 1175.4085234394684*eccentricity_14 + 679.98254677854938*eccentricity_12 + 364.97671006944444*eccentricity_10 + 177.49001736111111*eccentricity_8 + 75.336805555555556*eccentricity_6 + 26.1875*eccentricity_4 + 6.5*eccentricity_2 + 1.0;
    double common_term_97 = 3835.3352138753201*eccentricity_19 + 2535.9961280482471*eccentricity_17 + 1609.2736452334876*eccentricity_15 + 971.27648251942226*eccentricity_13 + 550.74737289428711*eccentricity_11 + 288.33571166992188*eccentricity_9 + 135.78759765625*eccentricity_7 + 55.4140625*eccentricity_5 + 17.3125*eccentricity_3 + 5.5*eccentricity;
    double common_term_98 = 3292.1257042683656*eccentricity_18 + 2150.0762100640314*eccentricity_16 + 1344.6030678185626*eccentricity_14 + 797.42890294312169*eccentricity_12 + 442.60251736111111*eccentricity_10 + 225.23888888888889*eccentricity_8 + 104.10416666666667*eccentricity_6 + 35.0*eccentricity_4 + 19.0*eccentricity_2;
    double common_term_99 = 4202.5261790501734*eccentricity_19 + 2814.0108402993383*eccentricity_17 + 1813.9819219717349*eccentricity_15 + 1116.9114152054307*eccentricity_13 + 650.50363742485248*eccentricity_11 + 350.33441614221644*eccentricity_9 + 183.17464192708333*eccentricity_7 + 55.384114583333333*eccentricity_5 + 52.604166666666667*eccentricity_3;
    double common_term_100 = 3618.5038254789731*eccentricity_18 + 2394.6618194231306*eccentricity_16 + 1521.9739501953125*eccentricity_14 + 924.69172014508929*eccentricity_12 + 513.92678571428571*eccentricity_10 + 315.759375*eccentricity_8 + 60.225*eccentricity_6 + 127.5*eccentricity_4;
    double common_term_101 = 4582.2019599417241*eccentricity_19 + 3103.3563236145644*eccentricity_17 + 2025.1875268339969*eccentricity_15 + 1284.5458838412791*eccentricity_13 + 703.90580623767994*eccentricity_11 + 556.04554714626736*eccentricity_9 - 0.72203776041666667*eccentricity_7 + 282.51223958333333*eccentricity_5;
    double common_term_102 = 3959.3322305698486*eccentricity_18 + 2634.4506103112182*eccentricity_16 + 1770.0839864016112*eccentricity_14 + 863.15170120517343*eccentricity_12 + 1031.6815538194444*eccentricity_10 - 244.28353174603175*eccentricity_8 + 586.63472222222222*eccentricity_6;
    double common_term_103 = 4988.8131773200446*eccentricity_19 + 3337.6673528652031*eccentricity_17 + 2481.9308368206256*eccentricity_15 + 815.92556585720607*eccentricity_13 + 2033.4356371743338*eccentricity_11 - 917.96151515415737*eccentricity_9 + 1159.5359514508929*eccentricity_7;
    double common_term_104 = 4056.1140025056224*eccentricity_18 + 3686.1121246735854*eccentricity_16 + 99.537103549080003*eccentricity_14 + 4195.9863647245646*eccentricity_12 - 2515.2845465443122*eccentricity_10 + 2204.7194320436508*eccentricity_8;
    double common_term_105 = 4514.8317519214222*eccentricity_19 + 6070.0886626209522*eccentricity_17 - 2400.6970032988403*eccentricity_15 + 8849.209973816325*eccentricity_13 - 5978.6031004550133*eccentricity_11 + 4062.6052557997392*eccentricity_9;
    double common_term_106 = 11334.191663178479*eccentricity_18 - 9210.2769493494786*eccentricity_16 + 18680.647027698864*eccentricity_14 - 13048.854811282468*eccentricity_12 + 7294.7463392857143*eccentricity_10;
    double common_term_107 = 23496.906469001575*eccentricity_19 - 25789.274531413559*eccentricity_17 + 38954.317796224522*eccentricity_15 - 26859.214041007832*eccentricity_13 + 12816.391140609443*eccentricity_11;
    double common_term_108 = -63505.036264258579*eccentricity_18 + 79706.951875949787*eccentricity_16 - 52928.991878443501*eccentricity_14 + 22103.683587078832*eccentricity_12;
    double common_term_109 = -145291.30699695886*eccentricity_19 + 159635.651941569*eccentricity_17 - 100808.05181331518*eccentricity_15 + 37515.480413075451*eccentricity_13;
    double common_term_110 = 312874.30504545611*eccentricity_18 - 186769.41197038686*eccentricity_16 + 62790.335431863512*eccentricity_14;
    double common_term_111 = 600646.58476021791*eccentricity_19 - 338177.82659816682*eccentricity_17 + 103810.00681260552*eccentricity_15;
    double common_term_112 = -600519.57383984692*eccentricity_18 + 169766.84288515738*eccentricity_16;
    double common_term_113 = -1048631.2712161361*eccentricity_19 + 274940.89334045742*eccentricity_17;
    double common_term_114 = 441394.30421414092*eccentricity_18;
    double common_term_115 = 703041.22488743427*eccentricity_19;
    double common_term_116 = 48485.244416646629*eccentricity_19;
    double common_term_117 = 32188.934890769995*eccentricity_18;
    double common_term_118 = -27837.987495213723*eccentricity_19 + 21263.489026116978*eccentricity_17;
    double common_term_119 = -14214.933231771154*eccentricity_18 + 13968.951714812345*eccentricity_16;
    double common_term_120 = 17270.193626971695*eccentricity_19 - 6582.4139006824708*eccentricity_17 + 9120.6484354619648*eccentricity_15;
    double common_term_121 = 10650.224142819642*eccentricity_18 - 2488.1435316882245*eccentricity_16 + 5914.2904116453048*eccentricity_14;
    double common_term_122 = 5640.5869828272548*eccentricity_19 + 6861.454586724783*eccentricity_17 - 432.35112940257249*eccentricity_15 + 3805.5325424087962*eccentricity_13;
    double common_term_123 = 5151.5794015705974*eccentricity_18 + 4628.8066493028456*eccentricity_16 + 487.45542326813811*eccentricity_14 + 2427.2028231534091*eccentricity_12;
    double common_term_124 = 7262.8269905151646*eccentricity_19 + 4391.2025251699104*eccentricity_17 + 3251.694003623174*eccentricity_15 + 803.49230794009848*eccentricity_13 + 1532.5469050381507*eccentricity_11;
    double common_term_125 = 5887.8593412949735*eccentricity_18 + 3599.7910864969061*eccentricity_16 + 2351.9291003477337*eccentricity_14 + 821.53377744458474*eccentricity_12 + 956.40387490354938*eccentricity_10;
    double common_term_126 = 7352.0443589849854*eccentricity_19 + 4745.5366468403684*eccentricity_17 + 2876.4945120911707*eccentricity_15 + 1728.3684608225389*eccentricity_13 + 711.15719764709473*eccentricity_11 + 588.71616516113281*eccentricity_9;
    double common_term_127 = 6000.8984540274657*eccentricity_18 + 3791.1443781398208*eccentricity_16 + 2254.0798681403152*eccentricity_14 + 1274.7245322145062*eccentricity_12 + 561.71428915895062*eccentricity_10 + 356.50512152777778*eccentricity_8;
    double common_term_128 = 7481.8360102649544*eccentricity_19 + 4851.0151231210192*eccentricity_17 + 2995.7861820289237*eccentricity_15 + 1736.0021989555908*eccentricity_13 + 934.16348894378285*eccentricity_11 + 416.56784430609809*eccentricity_9 + 211.64350043402778*eccentricity_7;
    double common_term_129 = 6106.0745450217507*eccentricity_18 + 3880.0526324573864*eccentricity_16 + 2337.4749462890625*eccentricity_14 + 1314.013671875*eccentricity_12 + 674.9208984375*eccentricity_10 + 293.68125*eccentricity_8 + 122.5875*eccentricity_6;
    double common_term_130 = 7578.2218510225676*eccentricity_19 + 4934.6530023668754*eccentricity_17 + 3066.9197633159136*eccentricity_15 + 1797.6122005314493*eccentricity_13 + 976.08916423938892*eccentricity_11 + 477.65137396918403*eccentricity_9 + 197.73822699652778*eccentricity_7 + 68.805989583333333*eccentricity_5;
    double common_term_131 = 6180.9889529862131*eccentricity_18 + 3944.566841425808*eccentricity_16 + 2392.0661497717335*eccentricity_14 + 1359.6920211226852*eccentricity_12 + 709.74508101851852*eccentricity_10 + 329.14340277777778*eccentricity_8 + 127.07916666666667*eccentricity_6 + 37.041666666666667*eccentricity_4;
    double common_term_132 = 7642.4603217811706*eccentricity_19 + 4990.3588658984325*eccentricity_17 + 3114.5109671384096*eccentricity_15 + 1837.5225867605209*eccentricity_13 + 1008.8510215759277*eccentricity_11 + 503.28966064453125*eccentricity_9 + 219.36181640625*eccentricity_7 + 77.51953125*eccentricity_5 + 18.8125*eccentricity_3;
    double common_term_133 = 6225.9237123259272*eccentricity_18 + 3983.2556593509639*eccentricity_16 + 2424.8300771002122*eccentricity_14 + 1386.8702835648148*eccentricity_12 + 731.70609085648148*eccentricity_10 + 346.21319444444444*eccentricity_8 + 140.14583333333333*eccentricity_6 + 44.333333333333333*eccentricity_4 + 8.75*eccentricity_2;
    double common_term_134 = 7674.572936199302*eccentricity_19 + 5018.2041244441034*eccentricity_17 + 3138.2972657430161*eccentricity_15 + 1857.4697025810053*eccentricity_13 + 1025.192222555655*eccentricity_11 + 516.27495591905382*eccentricity_9 + 229.25645616319444*eccentricity_7 + 84.674479166666667*eccentricity_5 + 23.1875*eccentricity_3 + 3.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(39);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (6, 0).
    // q = -19
    result_by_lpq.set(c_Key3(6, 0, -19), common_term_0);
    result_by_q.set(c_Key1(-19), common_term_0);
    // q = -18
    result_by_lpq.set(c_Key3(6, 0, -18), common_term_1);
    result_by_q.set(c_Key1(-18), common_term_1);
    // q = -17
    result_by_lpq.set(c_Key3(6, 0, -17), common_term_2);
    result_by_q.set(c_Key1(-17), common_term_2);
    // q = -16
    result_by_lpq.set(c_Key3(6, 0, -16), common_term_3);
    result_by_q.set(c_Key1(-16), common_term_3);
    // q = -15
    result_by_lpq.set(c_Key3(6, 0, -15), common_term_4);
    result_by_q.set(c_Key1(-15), common_term_4);
    // q = -14
    result_by_lpq.set(c_Key3(6, 0, -14), common_term_5);
    result_by_q.set(c_Key1(-14), common_term_5);
    // q = -13
    result_by_lpq.set(c_Key3(6, 0, -13), common_term_6);
    result_by_q.set(c_Key1(-13), common_term_6);
    // q = -12
    result_by_lpq.set(c_Key3(6, 0, -12), common_term_7);
    result_by_q.set(c_Key1(-12), common_term_7);
    // q = -11
    result_by_lpq.set(c_Key3(6, 0, -11), common_term_8);
    result_by_q.set(c_Key1(-11), common_term_8);
    // q = -10
    result_by_lpq.set(c_Key3(6, 0, -10), common_term_9);
    result_by_q.set(c_Key1(-10), common_term_9);
    // q = -9
    result_by_lpq.set(c_Key3(6, 0, -9), common_term_10);
    result_by_q.set(c_Key1(-9), common_term_10);
    // q = -8
    result_by_lpq.set(c_Key3(6, 0, -8), common_term_11);
    result_by_q.set(c_Key1(-8), common_term_11);
    // q = -7
    result_by_lpq.set(c_Key3(6, 0, -7), common_term_12);
    result_by_q.set(c_Key1(-7), common_term_12);
    // q = -5
    result_by_lpq.set(c_Key3(6, 0, -5), common_term_13);
    result_by_q.set(c_Key1(-5), common_term_13);
    // q = -4
    result_by_lpq.set(c_Key3(6, 0, -4), common_term_14);
    result_by_q.set(c_Key1(-4), common_term_14);
    // q = -3
    result_by_lpq.set(c_Key3(6, 0, -3), common_term_15);
    result_by_q.set(c_Key1(-3), common_term_15);
    // q = -2
    result_by_lpq.set(c_Key3(6, 0, -2), common_term_16);
    result_by_q.set(c_Key1(-2), common_term_16);
    // q = -1
    result_by_lpq.set(c_Key3(6, 0, -1), common_term_17);
    result_by_q.set(c_Key1(-1), common_term_17);
    // q = 0
    result_by_lpq.set(c_Key3(6, 0, 0), common_term_18);
    result_by_q.set(c_Key1(0), common_term_18);
    // q = 1
    result_by_lpq.set(c_Key3(6, 0, 1), common_term_19);
    result_by_q.set(c_Key1(1), common_term_19);
    // q = 2
    result_by_lpq.set(c_Key3(6, 0, 2), common_term_20);
    result_by_q.set(c_Key1(2), common_term_20);
    // q = 3
    result_by_lpq.set(c_Key3(6, 0, 3), common_term_21);
    result_by_q.set(c_Key1(3), common_term_21);
    // q = 4
    result_by_lpq.set(c_Key3(6, 0, 4), common_term_22);
    result_by_q.set(c_Key1(4), common_term_22);
    // q = 5
    result_by_lpq.set(c_Key3(6, 0, 5), common_term_23);
    result_by_q.set(c_Key1(5), common_term_23);
    // q = 6
    result_by_lpq.set(c_Key3(6, 0, 6), common_term_24);
    result_by_q.set(c_Key1(6), common_term_24);
    // q = 7
    result_by_lpq.set(c_Key3(6, 0, 7), common_term_25);
    result_by_q.set(c_Key1(7), common_term_25);
    // q = 8
    result_by_lpq.set(c_Key3(6, 0, 8), common_term_26);
    result_by_q.set(c_Key1(8), common_term_26);
    // q = 9
    result_by_lpq.set(c_Key3(6, 0, 9), common_term_27);
    result_by_q.set(c_Key1(9), common_term_27);
    // q = 10
    result_by_lpq.set(c_Key3(6, 0, 10), common_term_28);
    result_by_q.set(c_Key1(10), common_term_28);
    // q = 11
    result_by_lpq.set(c_Key3(6, 0, 11), common_term_29);
    result_by_q.set(c_Key1(11), common_term_29);
    // q = 12
    result_by_lpq.set(c_Key3(6, 0, 12), common_term_30);
    result_by_q.set(c_Key1(12), common_term_30);
    // q = 13
    result_by_lpq.set(c_Key3(6, 0, 13), common_term_31);
    result_by_q.set(c_Key1(13), common_term_31);
    // q = 14
    result_by_lpq.set(c_Key3(6, 0, 14), common_term_32);
    result_by_q.set(c_Key1(14), common_term_32);
    // q = 15
    result_by_lpq.set(c_Key3(6, 0, 15), common_term_33);
    result_by_q.set(c_Key1(15), common_term_33);
    // q = 16
    result_by_lpq.set(c_Key3(6, 0, 16), common_term_34);
    result_by_q.set(c_Key1(16), common_term_34);
    // q = 17
    result_by_lpq.set(c_Key3(6, 0, 17), common_term_35);
    result_by_q.set(c_Key1(17), common_term_35);
    // q = 18
    result_by_lpq.set(c_Key3(6, 0, 18), common_term_36);
    result_by_q.set(c_Key1(18), common_term_36);
    // q = 19
    result_by_lpq.set(c_Key3(6, 0, 19), common_term_37);
    result_by_q.set(c_Key1(19), common_term_37);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 0), result_by_q);
    result_by_q.clear();

    // l , p = (6, 1).
    // q = -19
    result_by_lpq.set(c_Key3(6, 1, -19), common_term_38);
    result_by_q.set(c_Key1(-19), common_term_38);
    // q = -18
    result_by_lpq.set(c_Key3(6, 1, -18), common_term_39);
    result_by_q.set(c_Key1(-18), common_term_39);
    // q = -17
    result_by_lpq.set(c_Key3(6, 1, -17), common_term_40);
    result_by_q.set(c_Key1(-17), common_term_40);
    // q = -16
    result_by_lpq.set(c_Key3(6, 1, -16), common_term_41);
    result_by_q.set(c_Key1(-16), common_term_41);
    // q = -15
    result_by_lpq.set(c_Key3(6, 1, -15), common_term_42);
    result_by_q.set(c_Key1(-15), common_term_42);
    // q = -14
    result_by_lpq.set(c_Key3(6, 1, -14), common_term_43);
    result_by_q.set(c_Key1(-14), common_term_43);
    // q = -13
    result_by_lpq.set(c_Key3(6, 1, -13), common_term_44);
    result_by_q.set(c_Key1(-13), common_term_44);
    // q = -12
    result_by_lpq.set(c_Key3(6, 1, -12), common_term_45);
    result_by_q.set(c_Key1(-12), common_term_45);
    // q = -11
    result_by_lpq.set(c_Key3(6, 1, -11), common_term_46);
    result_by_q.set(c_Key1(-11), common_term_46);
    // q = -10
    result_by_lpq.set(c_Key3(6, 1, -10), common_term_47);
    result_by_q.set(c_Key1(-10), common_term_47);
    // q = -9
    result_by_lpq.set(c_Key3(6, 1, -9), common_term_48);
    result_by_q.set(c_Key1(-9), common_term_48);
    // q = -8
    result_by_lpq.set(c_Key3(6, 1, -8), common_term_49);
    result_by_q.set(c_Key1(-8), common_term_49);
    // q = -7
    result_by_lpq.set(c_Key3(6, 1, -7), common_term_50);
    result_by_q.set(c_Key1(-7), common_term_50);
    // q = -6
    result_by_lpq.set(c_Key3(6, 1, -6), common_term_51);
    result_by_q.set(c_Key1(-6), common_term_51);
    // q = -5
    result_by_lpq.set(c_Key3(6, 1, -5), common_term_52);
    result_by_q.set(c_Key1(-5), common_term_52);
    // q = -4
    result_by_lpq.set(c_Key3(6, 1, -4), common_term_53);
    result_by_q.set(c_Key1(-4), common_term_53);
    // q = -3
    result_by_lpq.set(c_Key3(6, 1, -3), common_term_54);
    result_by_q.set(c_Key1(-3), common_term_54);
    // q = -2
    result_by_lpq.set(c_Key3(6, 1, -2), common_term_55);
    result_by_q.set(c_Key1(-2), common_term_55);
    // q = -1
    result_by_lpq.set(c_Key3(6, 1, -1), common_term_56);
    result_by_q.set(c_Key1(-1), common_term_56);
    // q = 0
    result_by_lpq.set(c_Key3(6, 1, 0), common_term_57);
    result_by_q.set(c_Key1(0), common_term_57);
    // q = 1
    result_by_lpq.set(c_Key3(6, 1, 1), common_term_58);
    result_by_q.set(c_Key1(1), common_term_58);
    // q = 2
    result_by_lpq.set(c_Key3(6, 1, 2), common_term_59);
    result_by_q.set(c_Key1(2), common_term_59);
    // q = 3
    result_by_lpq.set(c_Key3(6, 1, 3), common_term_60);
    result_by_q.set(c_Key1(3), common_term_60);
    // q = 4
    result_by_lpq.set(c_Key3(6, 1, 4), common_term_61);
    result_by_q.set(c_Key1(4), common_term_61);
    // q = 5
    result_by_lpq.set(c_Key3(6, 1, 5), common_term_62);
    result_by_q.set(c_Key1(5), common_term_62);
    // q = 6
    result_by_lpq.set(c_Key3(6, 1, 6), common_term_63);
    result_by_q.set(c_Key1(6), common_term_63);
    // q = 7
    result_by_lpq.set(c_Key3(6, 1, 7), common_term_64);
    result_by_q.set(c_Key1(7), common_term_64);
    // q = 8
    result_by_lpq.set(c_Key3(6, 1, 8), common_term_65);
    result_by_q.set(c_Key1(8), common_term_65);
    // q = 9
    result_by_lpq.set(c_Key3(6, 1, 9), common_term_66);
    result_by_q.set(c_Key1(9), common_term_66);
    // q = 10
    result_by_lpq.set(c_Key3(6, 1, 10), common_term_67);
    result_by_q.set(c_Key1(10), common_term_67);
    // q = 11
    result_by_lpq.set(c_Key3(6, 1, 11), common_term_68);
    result_by_q.set(c_Key1(11), common_term_68);
    // q = 12
    result_by_lpq.set(c_Key3(6, 1, 12), common_term_69);
    result_by_q.set(c_Key1(12), common_term_69);
    // q = 13
    result_by_lpq.set(c_Key3(6, 1, 13), common_term_70);
    result_by_q.set(c_Key1(13), common_term_70);
    // q = 14
    result_by_lpq.set(c_Key3(6, 1, 14), common_term_71);
    result_by_q.set(c_Key1(14), common_term_71);
    // q = 15
    result_by_lpq.set(c_Key3(6, 1, 15), common_term_72);
    result_by_q.set(c_Key1(15), common_term_72);
    // q = 16
    result_by_lpq.set(c_Key3(6, 1, 16), common_term_73);
    result_by_q.set(c_Key1(16), common_term_73);
    // q = 17
    result_by_lpq.set(c_Key3(6, 1, 17), common_term_74);
    result_by_q.set(c_Key1(17), common_term_74);
    // q = 18
    result_by_lpq.set(c_Key3(6, 1, 18), common_term_75);
    result_by_q.set(c_Key1(18), common_term_75);
    // q = 19
    result_by_lpq.set(c_Key3(6, 1, 19), common_term_76);
    result_by_q.set(c_Key1(19), common_term_76);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 1), result_by_q);
    result_by_q.clear();

    // l , p = (6, 2).
    // q = -19
    result_by_lpq.set(c_Key3(6, 2, -19), common_term_77);
    result_by_q.set(c_Key1(-19), common_term_77);
    // q = -18
    result_by_lpq.set(c_Key3(6, 2, -18), common_term_78);
    result_by_q.set(c_Key1(-18), common_term_78);
    // q = -17
    result_by_lpq.set(c_Key3(6, 2, -17), common_term_79);
    result_by_q.set(c_Key1(-17), common_term_79);
    // q = -16
    result_by_lpq.set(c_Key3(6, 2, -16), common_term_80);
    result_by_q.set(c_Key1(-16), common_term_80);
    // q = -15
    result_by_lpq.set(c_Key3(6, 2, -15), common_term_81);
    result_by_q.set(c_Key1(-15), common_term_81);
    // q = -14
    result_by_lpq.set(c_Key3(6, 2, -14), common_term_82);
    result_by_q.set(c_Key1(-14), common_term_82);
    // q = -13
    result_by_lpq.set(c_Key3(6, 2, -13), common_term_83);
    result_by_q.set(c_Key1(-13), common_term_83);
    // q = -12
    result_by_lpq.set(c_Key3(6, 2, -12), common_term_84);
    result_by_q.set(c_Key1(-12), common_term_84);
    // q = -11
    result_by_lpq.set(c_Key3(6, 2, -11), common_term_85);
    result_by_q.set(c_Key1(-11), common_term_85);
    // q = -10
    result_by_lpq.set(c_Key3(6, 2, -10), common_term_86);
    result_by_q.set(c_Key1(-10), common_term_86);
    // q = -9
    result_by_lpq.set(c_Key3(6, 2, -9), common_term_87);
    result_by_q.set(c_Key1(-9), common_term_87);
    // q = -8
    result_by_lpq.set(c_Key3(6, 2, -8), common_term_88);
    result_by_q.set(c_Key1(-8), common_term_88);
    // q = -7
    result_by_lpq.set(c_Key3(6, 2, -7), common_term_89);
    result_by_q.set(c_Key1(-7), common_term_89);
    // q = -6
    result_by_lpq.set(c_Key3(6, 2, -6), common_term_90);
    result_by_q.set(c_Key1(-6), common_term_90);
    // q = -5
    result_by_lpq.set(c_Key3(6, 2, -5), common_term_91);
    result_by_q.set(c_Key1(-5), common_term_91);
    // q = -4
    result_by_lpq.set(c_Key3(6, 2, -4), common_term_92);
    result_by_q.set(c_Key1(-4), common_term_92);
    // q = -3
    result_by_lpq.set(c_Key3(6, 2, -3), common_term_93);
    result_by_q.set(c_Key1(-3), common_term_93);
    // q = -2
    result_by_lpq.set(c_Key3(6, 2, -2), common_term_94);
    result_by_q.set(c_Key1(-2), common_term_94);
    // q = -1
    result_by_lpq.set(c_Key3(6, 2, -1), common_term_95);
    result_by_q.set(c_Key1(-1), common_term_95);
    // q = 0
    result_by_lpq.set(c_Key3(6, 2, 0), common_term_96);
    result_by_q.set(c_Key1(0), common_term_96);
    // q = 1
    result_by_lpq.set(c_Key3(6, 2, 1), common_term_97);
    result_by_q.set(c_Key1(1), common_term_97);
    // q = 2
    result_by_lpq.set(c_Key3(6, 2, 2), common_term_98);
    result_by_q.set(c_Key1(2), common_term_98);
    // q = 3
    result_by_lpq.set(c_Key3(6, 2, 3), common_term_99);
    result_by_q.set(c_Key1(3), common_term_99);
    // q = 4
    result_by_lpq.set(c_Key3(6, 2, 4), common_term_100);
    result_by_q.set(c_Key1(4), common_term_100);
    // q = 5
    result_by_lpq.set(c_Key3(6, 2, 5), common_term_101);
    result_by_q.set(c_Key1(5), common_term_101);
    // q = 6
    result_by_lpq.set(c_Key3(6, 2, 6), common_term_102);
    result_by_q.set(c_Key1(6), common_term_102);
    // q = 7
    result_by_lpq.set(c_Key3(6, 2, 7), common_term_103);
    result_by_q.set(c_Key1(7), common_term_103);
    // q = 8
    result_by_lpq.set(c_Key3(6, 2, 8), common_term_104);
    result_by_q.set(c_Key1(8), common_term_104);
    // q = 9
    result_by_lpq.set(c_Key3(6, 2, 9), common_term_105);
    result_by_q.set(c_Key1(9), common_term_105);
    // q = 10
    result_by_lpq.set(c_Key3(6, 2, 10), common_term_106);
    result_by_q.set(c_Key1(10), common_term_106);
    // q = 11
    result_by_lpq.set(c_Key3(6, 2, 11), common_term_107);
    result_by_q.set(c_Key1(11), common_term_107);
    // q = 12
    result_by_lpq.set(c_Key3(6, 2, 12), common_term_108);
    result_by_q.set(c_Key1(12), common_term_108);
    // q = 13
    result_by_lpq.set(c_Key3(6, 2, 13), common_term_109);
    result_by_q.set(c_Key1(13), common_term_109);
    // q = 14
    result_by_lpq.set(c_Key3(6, 2, 14), common_term_110);
    result_by_q.set(c_Key1(14), common_term_110);
    // q = 15
    result_by_lpq.set(c_Key3(6, 2, 15), common_term_111);
    result_by_q.set(c_Key1(15), common_term_111);
    // q = 16
    result_by_lpq.set(c_Key3(6, 2, 16), common_term_112);
    result_by_q.set(c_Key1(16), common_term_112);
    // q = 17
    result_by_lpq.set(c_Key3(6, 2, 17), common_term_113);
    result_by_q.set(c_Key1(17), common_term_113);
    // q = 18
    result_by_lpq.set(c_Key3(6, 2, 18), common_term_114);
    result_by_q.set(c_Key1(18), common_term_114);
    // q = 19
    result_by_lpq.set(c_Key3(6, 2, 19), common_term_115);
    result_by_q.set(c_Key1(19), common_term_115);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 2), result_by_q);
    result_by_q.clear();

    // l , p = (6, 3).
    // q = -19
    result_by_lpq.set(c_Key3(6, 3, -19), common_term_116);
    result_by_q.set(c_Key1(-19), common_term_116);
    // q = -18
    result_by_lpq.set(c_Key3(6, 3, -18), common_term_117);
    result_by_q.set(c_Key1(-18), common_term_117);
    // q = -17
    result_by_lpq.set(c_Key3(6, 3, -17), common_term_118);
    result_by_q.set(c_Key1(-17), common_term_118);
    // q = -16
    result_by_lpq.set(c_Key3(6, 3, -16), common_term_119);
    result_by_q.set(c_Key1(-16), common_term_119);
    // q = -15
    result_by_lpq.set(c_Key3(6, 3, -15), common_term_120);
    result_by_q.set(c_Key1(-15), common_term_120);
    // q = -14
    result_by_lpq.set(c_Key3(6, 3, -14), common_term_121);
    result_by_q.set(c_Key1(-14), common_term_121);
    // q = -13
    result_by_lpq.set(c_Key3(6, 3, -13), common_term_122);
    result_by_q.set(c_Key1(-13), common_term_122);
    // q = -12
    result_by_lpq.set(c_Key3(6, 3, -12), common_term_123);
    result_by_q.set(c_Key1(-12), common_term_123);
    // q = -11
    result_by_lpq.set(c_Key3(6, 3, -11), common_term_124);
    result_by_q.set(c_Key1(-11), common_term_124);
    // q = -10
    result_by_lpq.set(c_Key3(6, 3, -10), common_term_125);
    result_by_q.set(c_Key1(-10), common_term_125);
    // q = -9
    result_by_lpq.set(c_Key3(6, 3, -9), common_term_126);
    result_by_q.set(c_Key1(-9), common_term_126);
    // q = -8
    result_by_lpq.set(c_Key3(6, 3, -8), common_term_127);
    result_by_q.set(c_Key1(-8), common_term_127);
    // q = -7
    result_by_lpq.set(c_Key3(6, 3, -7), common_term_128);
    result_by_q.set(c_Key1(-7), common_term_128);
    // q = -6
    result_by_lpq.set(c_Key3(6, 3, -6), common_term_129);
    result_by_q.set(c_Key1(-6), common_term_129);
    // q = -5
    result_by_lpq.set(c_Key3(6, 3, -5), common_term_130);
    result_by_q.set(c_Key1(-5), common_term_130);
    // q = -4
    result_by_lpq.set(c_Key3(6, 3, -4), common_term_131);
    result_by_q.set(c_Key1(-4), common_term_131);
    // q = -3
    result_by_lpq.set(c_Key3(6, 3, -3), common_term_132);
    result_by_q.set(c_Key1(-3), common_term_132);
    // q = -2
    result_by_lpq.set(c_Key3(6, 3, -2), common_term_133);
    result_by_q.set(c_Key1(-2), common_term_133);
    // q = -1
    result_by_lpq.set(c_Key3(6, 3, -1), common_term_134);
    result_by_q.set(c_Key1(-1), common_term_134);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -5.5)*(1.875*eccentricity_4 + 5.0*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(6, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(6, 3, 1), common_term_134);
    result_by_q.set(c_Key1(1), common_term_134);
    // q = 2
    result_by_lpq.set(c_Key3(6, 3, 2), common_term_133);
    result_by_q.set(c_Key1(2), common_term_133);
    // q = 3
    result_by_lpq.set(c_Key3(6, 3, 3), common_term_132);
    result_by_q.set(c_Key1(3), common_term_132);
    // q = 4
    result_by_lpq.set(c_Key3(6, 3, 4), common_term_131);
    result_by_q.set(c_Key1(4), common_term_131);
    // q = 5
    result_by_lpq.set(c_Key3(6, 3, 5), common_term_130);
    result_by_q.set(c_Key1(5), common_term_130);
    // q = 6
    result_by_lpq.set(c_Key3(6, 3, 6), common_term_129);
    result_by_q.set(c_Key1(6), common_term_129);
    // q = 7
    result_by_lpq.set(c_Key3(6, 3, 7), common_term_128);
    result_by_q.set(c_Key1(7), common_term_128);
    // q = 8
    result_by_lpq.set(c_Key3(6, 3, 8), common_term_127);
    result_by_q.set(c_Key1(8), common_term_127);
    // q = 9
    result_by_lpq.set(c_Key3(6, 3, 9), common_term_126);
    result_by_q.set(c_Key1(9), common_term_126);
    // q = 10
    result_by_lpq.set(c_Key3(6, 3, 10), common_term_125);
    result_by_q.set(c_Key1(10), common_term_125);
    // q = 11
    result_by_lpq.set(c_Key3(6, 3, 11), common_term_124);
    result_by_q.set(c_Key1(11), common_term_124);
    // q = 12
    result_by_lpq.set(c_Key3(6, 3, 12), common_term_123);
    result_by_q.set(c_Key1(12), common_term_123);
    // q = 13
    result_by_lpq.set(c_Key3(6, 3, 13), common_term_122);
    result_by_q.set(c_Key1(13), common_term_122);
    // q = 14
    result_by_lpq.set(c_Key3(6, 3, 14), common_term_121);
    result_by_q.set(c_Key1(14), common_term_121);
    // q = 15
    result_by_lpq.set(c_Key3(6, 3, 15), common_term_120);
    result_by_q.set(c_Key1(15), common_term_120);
    // q = 16
    result_by_lpq.set(c_Key3(6, 3, 16), common_term_119);
    result_by_q.set(c_Key1(16), common_term_119);
    // q = 17
    result_by_lpq.set(c_Key3(6, 3, 17), common_term_118);
    result_by_q.set(c_Key1(17), common_term_118);
    // q = 18
    result_by_lpq.set(c_Key3(6, 3, 18), common_term_117);
    result_by_q.set(c_Key1(18), common_term_117);
    // q = 19
    result_by_lpq.set(c_Key3(6, 3, 19), common_term_116);
    result_by_q.set(c_Key1(19), common_term_116);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 3), result_by_q);
    result_by_q.clear();

    // l , p = (6, 4).
    // q = -19
    result_by_lpq.set(c_Key3(6, 4, -19), common_term_115);
    result_by_q.set(c_Key1(-19), common_term_115);
    // q = -18
    result_by_lpq.set(c_Key3(6, 4, -18), common_term_114);
    result_by_q.set(c_Key1(-18), common_term_114);
    // q = -17
    result_by_lpq.set(c_Key3(6, 4, -17), common_term_113);
    result_by_q.set(c_Key1(-17), common_term_113);
    // q = -16
    result_by_lpq.set(c_Key3(6, 4, -16), common_term_112);
    result_by_q.set(c_Key1(-16), common_term_112);
    // q = -15
    result_by_lpq.set(c_Key3(6, 4, -15), common_term_111);
    result_by_q.set(c_Key1(-15), common_term_111);
    // q = -14
    result_by_lpq.set(c_Key3(6, 4, -14), common_term_110);
    result_by_q.set(c_Key1(-14), common_term_110);
    // q = -13
    result_by_lpq.set(c_Key3(6, 4, -13), common_term_109);
    result_by_q.set(c_Key1(-13), common_term_109);
    // q = -12
    result_by_lpq.set(c_Key3(6, 4, -12), common_term_108);
    result_by_q.set(c_Key1(-12), common_term_108);
    // q = -11
    result_by_lpq.set(c_Key3(6, 4, -11), common_term_107);
    result_by_q.set(c_Key1(-11), common_term_107);
    // q = -10
    result_by_lpq.set(c_Key3(6, 4, -10), common_term_106);
    result_by_q.set(c_Key1(-10), common_term_106);
    // q = -9
    result_by_lpq.set(c_Key3(6, 4, -9), common_term_105);
    result_by_q.set(c_Key1(-9), common_term_105);
    // q = -8
    result_by_lpq.set(c_Key3(6, 4, -8), common_term_104);
    result_by_q.set(c_Key1(-8), common_term_104);
    // q = -7
    result_by_lpq.set(c_Key3(6, 4, -7), common_term_103);
    result_by_q.set(c_Key1(-7), common_term_103);
    // q = -6
    result_by_lpq.set(c_Key3(6, 4, -6), common_term_102);
    result_by_q.set(c_Key1(-6), common_term_102);
    // q = -5
    result_by_lpq.set(c_Key3(6, 4, -5), common_term_101);
    result_by_q.set(c_Key1(-5), common_term_101);
    // q = -4
    result_by_lpq.set(c_Key3(6, 4, -4), common_term_100);
    result_by_q.set(c_Key1(-4), common_term_100);
    // q = -3
    result_by_lpq.set(c_Key3(6, 4, -3), common_term_99);
    result_by_q.set(c_Key1(-3), common_term_99);
    // q = -2
    result_by_lpq.set(c_Key3(6, 4, -2), common_term_98);
    result_by_q.set(c_Key1(-2), common_term_98);
    // q = -1
    result_by_lpq.set(c_Key3(6, 4, -1), common_term_97);
    result_by_q.set(c_Key1(-1), common_term_97);
    // q = 0
    result_by_lpq.set(c_Key3(6, 4, 0), common_term_96);
    result_by_q.set(c_Key1(0), common_term_96);
    // q = 1
    result_by_lpq.set(c_Key3(6, 4, 1), common_term_95);
    result_by_q.set(c_Key1(1), common_term_95);
    // q = 2
    result_by_lpq.set(c_Key3(6, 4, 2), common_term_94);
    result_by_q.set(c_Key1(2), common_term_94);
    // q = 3
    result_by_lpq.set(c_Key3(6, 4, 3), common_term_93);
    result_by_q.set(c_Key1(3), common_term_93);
    // q = 4
    result_by_lpq.set(c_Key3(6, 4, 4), common_term_92);
    result_by_q.set(c_Key1(4), common_term_92);
    // q = 5
    result_by_lpq.set(c_Key3(6, 4, 5), common_term_91);
    result_by_q.set(c_Key1(5), common_term_91);
    // q = 6
    result_by_lpq.set(c_Key3(6, 4, 6), common_term_90);
    result_by_q.set(c_Key1(6), common_term_90);
    // q = 7
    result_by_lpq.set(c_Key3(6, 4, 7), common_term_89);
    result_by_q.set(c_Key1(7), common_term_89);
    // q = 8
    result_by_lpq.set(c_Key3(6, 4, 8), common_term_88);
    result_by_q.set(c_Key1(8), common_term_88);
    // q = 9
    result_by_lpq.set(c_Key3(6, 4, 9), common_term_87);
    result_by_q.set(c_Key1(9), common_term_87);
    // q = 10
    result_by_lpq.set(c_Key3(6, 4, 10), common_term_86);
    result_by_q.set(c_Key1(10), common_term_86);
    // q = 11
    result_by_lpq.set(c_Key3(6, 4, 11), common_term_85);
    result_by_q.set(c_Key1(11), common_term_85);
    // q = 12
    result_by_lpq.set(c_Key3(6, 4, 12), common_term_84);
    result_by_q.set(c_Key1(12), common_term_84);
    // q = 13
    result_by_lpq.set(c_Key3(6, 4, 13), common_term_83);
    result_by_q.set(c_Key1(13), common_term_83);
    // q = 14
    result_by_lpq.set(c_Key3(6, 4, 14), common_term_82);
    result_by_q.set(c_Key1(14), common_term_82);
    // q = 15
    result_by_lpq.set(c_Key3(6, 4, 15), common_term_81);
    result_by_q.set(c_Key1(15), common_term_81);
    // q = 16
    result_by_lpq.set(c_Key3(6, 4, 16), common_term_80);
    result_by_q.set(c_Key1(16), common_term_80);
    // q = 17
    result_by_lpq.set(c_Key3(6, 4, 17), common_term_79);
    result_by_q.set(c_Key1(17), common_term_79);
    // q = 18
    result_by_lpq.set(c_Key3(6, 4, 18), common_term_78);
    result_by_q.set(c_Key1(18), common_term_78);
    // q = 19
    result_by_lpq.set(c_Key3(6, 4, 19), common_term_77);
    result_by_q.set(c_Key1(19), common_term_77);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 4), result_by_q);
    result_by_q.clear();

    // l , p = (6, 5).
    // q = -19
    result_by_lpq.set(c_Key3(6, 5, -19), common_term_76);
    result_by_q.set(c_Key1(-19), common_term_76);
    // q = -18
    result_by_lpq.set(c_Key3(6, 5, -18), common_term_75);
    result_by_q.set(c_Key1(-18), common_term_75);
    // q = -17
    result_by_lpq.set(c_Key3(6, 5, -17), common_term_74);
    result_by_q.set(c_Key1(-17), common_term_74);
    // q = -16
    result_by_lpq.set(c_Key3(6, 5, -16), common_term_73);
    result_by_q.set(c_Key1(-16), common_term_73);
    // q = -15
    result_by_lpq.set(c_Key3(6, 5, -15), common_term_72);
    result_by_q.set(c_Key1(-15), common_term_72);
    // q = -14
    result_by_lpq.set(c_Key3(6, 5, -14), common_term_71);
    result_by_q.set(c_Key1(-14), common_term_71);
    // q = -13
    result_by_lpq.set(c_Key3(6, 5, -13), common_term_70);
    result_by_q.set(c_Key1(-13), common_term_70);
    // q = -12
    result_by_lpq.set(c_Key3(6, 5, -12), common_term_69);
    result_by_q.set(c_Key1(-12), common_term_69);
    // q = -11
    result_by_lpq.set(c_Key3(6, 5, -11), common_term_68);
    result_by_q.set(c_Key1(-11), common_term_68);
    // q = -10
    result_by_lpq.set(c_Key3(6, 5, -10), common_term_67);
    result_by_q.set(c_Key1(-10), common_term_67);
    // q = -9
    result_by_lpq.set(c_Key3(6, 5, -9), common_term_66);
    result_by_q.set(c_Key1(-9), common_term_66);
    // q = -8
    result_by_lpq.set(c_Key3(6, 5, -8), common_term_65);
    result_by_q.set(c_Key1(-8), common_term_65);
    // q = -7
    result_by_lpq.set(c_Key3(6, 5, -7), common_term_64);
    result_by_q.set(c_Key1(-7), common_term_64);
    // q = -6
    result_by_lpq.set(c_Key3(6, 5, -6), common_term_63);
    result_by_q.set(c_Key1(-6), common_term_63);
    // q = -5
    result_by_lpq.set(c_Key3(6, 5, -5), common_term_62);
    result_by_q.set(c_Key1(-5), common_term_62);
    // q = -4
    result_by_lpq.set(c_Key3(6, 5, -4), common_term_61);
    result_by_q.set(c_Key1(-4), common_term_61);
    // q = -3
    result_by_lpq.set(c_Key3(6, 5, -3), common_term_60);
    result_by_q.set(c_Key1(-3), common_term_60);
    // q = -2
    result_by_lpq.set(c_Key3(6, 5, -2), common_term_59);
    result_by_q.set(c_Key1(-2), common_term_59);
    // q = -1
    result_by_lpq.set(c_Key3(6, 5, -1), common_term_58);
    result_by_q.set(c_Key1(-1), common_term_58);
    // q = 0
    result_by_lpq.set(c_Key3(6, 5, 0), common_term_57);
    result_by_q.set(c_Key1(0), common_term_57);
    // q = 1
    result_by_lpq.set(c_Key3(6, 5, 1), common_term_56);
    result_by_q.set(c_Key1(1), common_term_56);
    // q = 2
    result_by_lpq.set(c_Key3(6, 5, 2), common_term_55);
    result_by_q.set(c_Key1(2), common_term_55);
    // q = 3
    result_by_lpq.set(c_Key3(6, 5, 3), common_term_54);
    result_by_q.set(c_Key1(3), common_term_54);
    // q = 4
    result_by_lpq.set(c_Key3(6, 5, 4), common_term_53);
    result_by_q.set(c_Key1(4), common_term_53);
    // q = 5
    result_by_lpq.set(c_Key3(6, 5, 5), common_term_52);
    result_by_q.set(c_Key1(5), common_term_52);
    // q = 6
    result_by_lpq.set(c_Key3(6, 5, 6), common_term_51);
    result_by_q.set(c_Key1(6), common_term_51);
    // q = 7
    result_by_lpq.set(c_Key3(6, 5, 7), common_term_50);
    result_by_q.set(c_Key1(7), common_term_50);
    // q = 8
    result_by_lpq.set(c_Key3(6, 5, 8), common_term_49);
    result_by_q.set(c_Key1(8), common_term_49);
    // q = 9
    result_by_lpq.set(c_Key3(6, 5, 9), common_term_48);
    result_by_q.set(c_Key1(9), common_term_48);
    // q = 10
    result_by_lpq.set(c_Key3(6, 5, 10), common_term_47);
    result_by_q.set(c_Key1(10), common_term_47);
    // q = 11
    result_by_lpq.set(c_Key3(6, 5, 11), common_term_46);
    result_by_q.set(c_Key1(11), common_term_46);
    // q = 12
    result_by_lpq.set(c_Key3(6, 5, 12), common_term_45);
    result_by_q.set(c_Key1(12), common_term_45);
    // q = 13
    result_by_lpq.set(c_Key3(6, 5, 13), common_term_44);
    result_by_q.set(c_Key1(13), common_term_44);
    // q = 14
    result_by_lpq.set(c_Key3(6, 5, 14), common_term_43);
    result_by_q.set(c_Key1(14), common_term_43);
    // q = 15
    result_by_lpq.set(c_Key3(6, 5, 15), common_term_42);
    result_by_q.set(c_Key1(15), common_term_42);
    // q = 16
    result_by_lpq.set(c_Key3(6, 5, 16), common_term_41);
    result_by_q.set(c_Key1(16), common_term_41);
    // q = 17
    result_by_lpq.set(c_Key3(6, 5, 17), common_term_40);
    result_by_q.set(c_Key1(17), common_term_40);
    // q = 18
    result_by_lpq.set(c_Key3(6, 5, 18), common_term_39);
    result_by_q.set(c_Key1(18), common_term_39);
    // q = 19
    result_by_lpq.set(c_Key3(6, 5, 19), common_term_38);
    result_by_q.set(c_Key1(19), common_term_38);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 5), result_by_q);
    result_by_q.clear();

    // l , p = (6, 6).
    // q = -19
    result_by_lpq.set(c_Key3(6, 6, -19), common_term_37);
    result_by_q.set(c_Key1(-19), common_term_37);
    // q = -18
    result_by_lpq.set(c_Key3(6, 6, -18), common_term_36);
    result_by_q.set(c_Key1(-18), common_term_36);
    // q = -17
    result_by_lpq.set(c_Key3(6, 6, -17), common_term_35);
    result_by_q.set(c_Key1(-17), common_term_35);
    // q = -16
    result_by_lpq.set(c_Key3(6, 6, -16), common_term_34);
    result_by_q.set(c_Key1(-16), common_term_34);
    // q = -15
    result_by_lpq.set(c_Key3(6, 6, -15), common_term_33);
    result_by_q.set(c_Key1(-15), common_term_33);
    // q = -14
    result_by_lpq.set(c_Key3(6, 6, -14), common_term_32);
    result_by_q.set(c_Key1(-14), common_term_32);
    // q = -13
    result_by_lpq.set(c_Key3(6, 6, -13), common_term_31);
    result_by_q.set(c_Key1(-13), common_term_31);
    // q = -12
    result_by_lpq.set(c_Key3(6, 6, -12), common_term_30);
    result_by_q.set(c_Key1(-12), common_term_30);
    // q = -11
    result_by_lpq.set(c_Key3(6, 6, -11), common_term_29);
    result_by_q.set(c_Key1(-11), common_term_29);
    // q = -10
    result_by_lpq.set(c_Key3(6, 6, -10), common_term_28);
    result_by_q.set(c_Key1(-10), common_term_28);
    // q = -9
    result_by_lpq.set(c_Key3(6, 6, -9), common_term_27);
    result_by_q.set(c_Key1(-9), common_term_27);
    // q = -8
    result_by_lpq.set(c_Key3(6, 6, -8), common_term_26);
    result_by_q.set(c_Key1(-8), common_term_26);
    // q = -7
    result_by_lpq.set(c_Key3(6, 6, -7), common_term_25);
    result_by_q.set(c_Key1(-7), common_term_25);
    // q = -6
    result_by_lpq.set(c_Key3(6, 6, -6), common_term_24);
    result_by_q.set(c_Key1(-6), common_term_24);
    // q = -5
    result_by_lpq.set(c_Key3(6, 6, -5), common_term_23);
    result_by_q.set(c_Key1(-5), common_term_23);
    // q = -4
    result_by_lpq.set(c_Key3(6, 6, -4), common_term_22);
    result_by_q.set(c_Key1(-4), common_term_22);
    // q = -3
    result_by_lpq.set(c_Key3(6, 6, -3), common_term_21);
    result_by_q.set(c_Key1(-3), common_term_21);
    // q = -2
    result_by_lpq.set(c_Key3(6, 6, -2), common_term_20);
    result_by_q.set(c_Key1(-2), common_term_20);
    // q = -1
    result_by_lpq.set(c_Key3(6, 6, -1), common_term_19);
    result_by_q.set(c_Key1(-1), common_term_19);
    // q = 0
    result_by_lpq.set(c_Key3(6, 6, 0), common_term_18);
    result_by_q.set(c_Key1(0), common_term_18);
    // q = 1
    result_by_lpq.set(c_Key3(6, 6, 1), common_term_17);
    result_by_q.set(c_Key1(1), common_term_17);
    // q = 2
    result_by_lpq.set(c_Key3(6, 6, 2), common_term_16);
    result_by_q.set(c_Key1(2), common_term_16);
    // q = 3
    result_by_lpq.set(c_Key3(6, 6, 3), common_term_15);
    result_by_q.set(c_Key1(3), common_term_15);
    // q = 4
    result_by_lpq.set(c_Key3(6, 6, 4), common_term_14);
    result_by_q.set(c_Key1(4), common_term_14);
    // q = 5
    result_by_lpq.set(c_Key3(6, 6, 5), common_term_13);
    result_by_q.set(c_Key1(5), common_term_13);
    // q = 7
    result_by_lpq.set(c_Key3(6, 6, 7), common_term_12);
    result_by_q.set(c_Key1(7), common_term_12);
    // q = 8
    result_by_lpq.set(c_Key3(6, 6, 8), common_term_11);
    result_by_q.set(c_Key1(8), common_term_11);
    // q = 9
    result_by_lpq.set(c_Key3(6, 6, 9), common_term_10);
    result_by_q.set(c_Key1(9), common_term_10);
    // q = 10
    result_by_lpq.set(c_Key3(6, 6, 10), common_term_9);
    result_by_q.set(c_Key1(10), common_term_9);
    // q = 11
    result_by_lpq.set(c_Key3(6, 6, 11), common_term_8);
    result_by_q.set(c_Key1(11), common_term_8);
    // q = 12
    result_by_lpq.set(c_Key3(6, 6, 12), common_term_7);
    result_by_q.set(c_Key1(12), common_term_7);
    // q = 13
    result_by_lpq.set(c_Key3(6, 6, 13), common_term_6);
    result_by_q.set(c_Key1(13), common_term_6);
    // q = 14
    result_by_lpq.set(c_Key3(6, 6, 14), common_term_5);
    result_by_q.set(c_Key1(14), common_term_5);
    // q = 15
    result_by_lpq.set(c_Key3(6, 6, 15), common_term_4);
    result_by_q.set(c_Key1(15), common_term_4);
    // q = 16
    result_by_lpq.set(c_Key3(6, 6, 16), common_term_3);
    result_by_q.set(c_Key1(16), common_term_3);
    // q = 17
    result_by_lpq.set(c_Key3(6, 6, 17), common_term_2);
    result_by_q.set(c_Key1(17), common_term_2);
    // q = 18
    result_by_lpq.set(c_Key3(6, 6, 18), common_term_1);
    result_by_q.set(c_Key1(18), common_term_1);
    // q = 19
    result_by_lpq.set(c_Key3(6, 6, 19), common_term_0);
    result_by_q.set(c_Key1(19), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(6, 6), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}
