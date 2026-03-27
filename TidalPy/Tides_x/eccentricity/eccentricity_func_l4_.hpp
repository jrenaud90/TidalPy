#pragma once

#include "eccentricity_common_.hpp"

EccentricityFuncOutput c_eccentricity_function_l4_e1(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 4.
    // Functions are taylor expanded and truncated at eccentricity^1.

    // Functions Calculated for l = 4.

    c_IntMap<c_Key3, double> result_by_lpq(5);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(5);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;

    c_IntMap<c_Key1, double> result_by_q(1);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (4, 0).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(4, 0, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 0), result_by_q);
    result_by_q.clear();

    // l , p = (4, 1).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(4, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 1), result_by_q);
    result_by_q.clear();

    // l , p = (4, 2).
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -3.5);
    result_by_lpq.set(c_Key3(4, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 2), result_by_q);
    result_by_q.clear();

    // l , p = (4, 3).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(4, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 3), result_by_q);
    result_by_q.clear();

    // l , p = (4, 4).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(4, 4, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 4), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l4_e2(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 4.
    // Functions are taylor expanded and truncated at eccentricity^2.

    // Functions Calculated for l = 4.

    c_IntMap<c_Key3, double> result_by_lpq(17);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(5);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = -1.5*eccentricity;
    double common_term_1 = 6.5*eccentricity;
    double common_term_2 = 0.75*eccentricity_2*std::pow(1.0 - eccentricity_2, -3.5);
    double common_term_3 = 0.5*eccentricity;
    double common_term_4 = 4.5*eccentricity;
    double common_term_5 = 2.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(4);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (4, 0).
    // q = -1
    result_by_lpq.set(c_Key3(4, 0, -1), common_term_0);
    result_by_q.set(c_Key1(-1), common_term_0);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(4, 0, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(4, 0, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 0), result_by_q);
    result_by_q.clear();

    // l , p = (4, 1).
    // q = -2
    result_by_lpq.set(c_Key3(4, 1, -2), common_term_2);
    result_by_q.set(c_Key1(-2), common_term_2);
    // q = -1
    result_by_lpq.set(c_Key3(4, 1, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(4, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(4, 1, 1), common_term_4);
    result_by_q.set(c_Key1(1), common_term_4);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 1), result_by_q);
    result_by_q.clear();

    // l , p = (4, 2).
    // q = -1
    result_by_lpq.set(c_Key3(4, 2, -1), common_term_5);
    result_by_q.set(c_Key1(-1), common_term_5);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -3.5)*(1.5*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(4, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(4, 2, 1), common_term_5);
    result_by_q.set(c_Key1(1), common_term_5);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 2), result_by_q);
    result_by_q.clear();

    // l , p = (4, 3).
    // q = -1
    result_by_lpq.set(c_Key3(4, 3, -1), common_term_4);
    result_by_q.set(c_Key1(-1), common_term_4);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(4, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(4, 3, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // q = 2
    result_by_lpq.set(c_Key3(4, 3, 2), common_term_2);
    result_by_q.set(c_Key1(2), common_term_2);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 3), result_by_q);
    result_by_q.clear();

    // l , p = (4, 4).
    // q = -1
    result_by_lpq.set(c_Key3(4, 4, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(4, 4, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(4, 4, 1), common_term_0);
    result_by_q.set(c_Key1(1), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 4), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l4_e3(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 4.
    // Functions are taylor expanded and truncated at eccentricity^3.

    // Functions Calculated for l = 4.

    c_IntMap<c_Key3, double> result_by_lpq(25);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(5);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = 0.5*eccentricity_2;
    double common_term_1 = -1.5*eccentricity;
    double common_term_2 = 1.0 - 11.0*eccentricity_2;
    double common_term_3 = 6.5*eccentricity;
    double common_term_4 = 25.5*eccentricity_2;
    double common_term_5 = 0.75*eccentricity_2*std::pow(1.0 - eccentricity_2, -3.5);
    double common_term_6 = 0.5*eccentricity;
    double common_term_7 = eccentricity_2 + 1.0;
    double common_term_8 = 4.5*eccentricity;
    double common_term_9 = 13.25*eccentricity_2;
    double common_term_10 = 5.0*eccentricity_2;
    double common_term_11 = 2.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(5);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (4, 0).
    // q = -2
    result_by_lpq.set(c_Key3(4, 0, -2), common_term_0);
    result_by_q.set(c_Key1(-2), common_term_0);
    // q = -1
    result_by_lpq.set(c_Key3(4, 0, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    result_by_lpq.set(c_Key3(4, 0, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(4, 0, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // q = 2
    result_by_lpq.set(c_Key3(4, 0, 2), common_term_4);
    result_by_q.set(c_Key1(2), common_term_4);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 0), result_by_q);
    result_by_q.clear();

    // l , p = (4, 1).
    // q = -2
    result_by_lpq.set(c_Key3(4, 1, -2), common_term_5);
    result_by_q.set(c_Key1(-2), common_term_5);
    // q = -1
    result_by_lpq.set(c_Key3(4, 1, -1), common_term_6);
    result_by_q.set(c_Key1(-1), common_term_6);
    // q = 0
    result_by_lpq.set(c_Key3(4, 1, 0), common_term_7);
    result_by_q.set(c_Key1(0), common_term_7);
    // q = 1
    result_by_lpq.set(c_Key3(4, 1, 1), common_term_8);
    result_by_q.set(c_Key1(1), common_term_8);
    // q = 2
    result_by_lpq.set(c_Key3(4, 1, 2), common_term_9);
    result_by_q.set(c_Key1(2), common_term_9);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 1), result_by_q);
    result_by_q.clear();

    // l , p = (4, 2).
    // q = -2
    result_by_lpq.set(c_Key3(4, 2, -2), common_term_10);
    result_by_q.set(c_Key1(-2), common_term_10);
    // q = -1
    result_by_lpq.set(c_Key3(4, 2, -1), common_term_11);
    result_by_q.set(c_Key1(-1), common_term_11);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -3.5)*(1.5*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(4, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(4, 2, 1), common_term_11);
    result_by_q.set(c_Key1(1), common_term_11);
    // q = 2
    result_by_lpq.set(c_Key3(4, 2, 2), common_term_10);
    result_by_q.set(c_Key1(2), common_term_10);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 2), result_by_q);
    result_by_q.clear();

    // l , p = (4, 3).
    // q = -2
    result_by_lpq.set(c_Key3(4, 3, -2), common_term_9);
    result_by_q.set(c_Key1(-2), common_term_9);
    // q = -1
    result_by_lpq.set(c_Key3(4, 3, -1), common_term_8);
    result_by_q.set(c_Key1(-1), common_term_8);
    // q = 0
    result_by_lpq.set(c_Key3(4, 3, 0), common_term_7);
    result_by_q.set(c_Key1(0), common_term_7);
    // q = 1
    result_by_lpq.set(c_Key3(4, 3, 1), common_term_6);
    result_by_q.set(c_Key1(1), common_term_6);
    // q = 2
    result_by_lpq.set(c_Key3(4, 3, 2), common_term_5);
    result_by_q.set(c_Key1(2), common_term_5);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 3), result_by_q);
    result_by_q.clear();

    // l , p = (4, 4).
    // q = -2
    result_by_lpq.set(c_Key3(4, 4, -2), common_term_4);
    result_by_q.set(c_Key1(-2), common_term_4);
    // q = -1
    result_by_lpq.set(c_Key3(4, 4, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    result_by_lpq.set(c_Key3(4, 4, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(4, 4, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // q = 2
    result_by_lpq.set(c_Key3(4, 4, 2), common_term_0);
    result_by_q.set(c_Key1(2), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 4), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l4_e4(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 4.
    // Functions are taylor expanded and truncated at eccentricity^4.

    // Functions Calculated for l = 4.

    c_IntMap<c_Key3, double> result_by_lpq(35);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(5);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double common_term_0 = -0.020833333333333333*eccentricity_3;
    double common_term_1 = 0.5*eccentricity_2;
    double common_term_2 = 4.6875*eccentricity_3 - 1.5*eccentricity;
    double common_term_3 = 1.0 - 11.0*eccentricity_2;
    double common_term_4 = -47.8125*eccentricity_3 + 6.5*eccentricity;
    double common_term_5 = 25.5*eccentricity_2;
    double common_term_6 = 78.145833333333333*eccentricity_3;
    double common_term_7 = 1.0208333333333333*eccentricity_3;
    double common_term_8 = 0.75*eccentricity_2*std::pow(1.0 - eccentricity_2, -3.5);
    double common_term_9 = 2.0625*eccentricity_3 + 0.5*eccentricity;
    double common_term_10 = eccentricity_2 + 1.0;
    double common_term_11 = -0.1875*eccentricity_3 + 4.5*eccentricity;
    double common_term_12 = 13.25*eccentricity_2;
    double common_term_13 = 32.104166666666667*eccentricity_3;
    double common_term_14 = 9.0625*eccentricity_3;
    double common_term_15 = 5.0*eccentricity_2;
    double common_term_16 = 8.4375*eccentricity_3 + 2.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(7);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (4, 0).
    // q = -3
    result_by_lpq.set(c_Key3(4, 0, -3), common_term_0);
    result_by_q.set(c_Key1(-3), common_term_0);
    // q = -2
    result_by_lpq.set(c_Key3(4, 0, -2), common_term_1);
    result_by_q.set(c_Key1(-2), common_term_1);
    // q = -1
    result_by_lpq.set(c_Key3(4, 0, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    result_by_lpq.set(c_Key3(4, 0, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(4, 0, 1), common_term_4);
    result_by_q.set(c_Key1(1), common_term_4);
    // q = 2
    result_by_lpq.set(c_Key3(4, 0, 2), common_term_5);
    result_by_q.set(c_Key1(2), common_term_5);
    // q = 3
    result_by_lpq.set(c_Key3(4, 0, 3), common_term_6);
    result_by_q.set(c_Key1(3), common_term_6);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 0), result_by_q);
    result_by_q.clear();

    // l , p = (4, 1).
    // q = -3
    result_by_lpq.set(c_Key3(4, 1, -3), common_term_7);
    result_by_q.set(c_Key1(-3), common_term_7);
    // q = -2
    result_by_lpq.set(c_Key3(4, 1, -2), common_term_8);
    result_by_q.set(c_Key1(-2), common_term_8);
    // q = -1
    result_by_lpq.set(c_Key3(4, 1, -1), common_term_9);
    result_by_q.set(c_Key1(-1), common_term_9);
    // q = 0
    result_by_lpq.set(c_Key3(4, 1, 0), common_term_10);
    result_by_q.set(c_Key1(0), common_term_10);
    // q = 1
    result_by_lpq.set(c_Key3(4, 1, 1), common_term_11);
    result_by_q.set(c_Key1(1), common_term_11);
    // q = 2
    result_by_lpq.set(c_Key3(4, 1, 2), common_term_12);
    result_by_q.set(c_Key1(2), common_term_12);
    // q = 3
    result_by_lpq.set(c_Key3(4, 1, 3), common_term_13);
    result_by_q.set(c_Key1(3), common_term_13);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 1), result_by_q);
    result_by_q.clear();

    // l , p = (4, 2).
    // q = -3
    result_by_lpq.set(c_Key3(4, 2, -3), common_term_14);
    result_by_q.set(c_Key1(-3), common_term_14);
    // q = -2
    result_by_lpq.set(c_Key3(4, 2, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(4, 2, -1), common_term_16);
    result_by_q.set(c_Key1(-1), common_term_16);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -3.5)*(1.5*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(4, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(4, 2, 1), common_term_16);
    result_by_q.set(c_Key1(1), common_term_16);
    // q = 2
    result_by_lpq.set(c_Key3(4, 2, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // q = 3
    result_by_lpq.set(c_Key3(4, 2, 3), common_term_14);
    result_by_q.set(c_Key1(3), common_term_14);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 2), result_by_q);
    result_by_q.clear();

    // l , p = (4, 3).
    // q = -3
    result_by_lpq.set(c_Key3(4, 3, -3), common_term_13);
    result_by_q.set(c_Key1(-3), common_term_13);
    // q = -2
    result_by_lpq.set(c_Key3(4, 3, -2), common_term_12);
    result_by_q.set(c_Key1(-2), common_term_12);
    // q = -1
    result_by_lpq.set(c_Key3(4, 3, -1), common_term_11);
    result_by_q.set(c_Key1(-1), common_term_11);
    // q = 0
    result_by_lpq.set(c_Key3(4, 3, 0), common_term_10);
    result_by_q.set(c_Key1(0), common_term_10);
    // q = 1
    result_by_lpq.set(c_Key3(4, 3, 1), common_term_9);
    result_by_q.set(c_Key1(1), common_term_9);
    // q = 2
    result_by_lpq.set(c_Key3(4, 3, 2), common_term_8);
    result_by_q.set(c_Key1(2), common_term_8);
    // q = 3
    result_by_lpq.set(c_Key3(4, 3, 3), common_term_7);
    result_by_q.set(c_Key1(3), common_term_7);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 3), result_by_q);
    result_by_q.clear();

    // l , p = (4, 4).
    // q = -3
    result_by_lpq.set(c_Key3(4, 4, -3), common_term_6);
    result_by_q.set(c_Key1(-3), common_term_6);
    // q = -2
    result_by_lpq.set(c_Key3(4, 4, -2), common_term_5);
    result_by_q.set(c_Key1(-2), common_term_5);
    // q = -1
    result_by_lpq.set(c_Key3(4, 4, -1), common_term_4);
    result_by_q.set(c_Key1(-1), common_term_4);
    // q = 0
    result_by_lpq.set(c_Key3(4, 4, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(4, 4, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // q = 2
    result_by_lpq.set(c_Key3(4, 4, 2), common_term_1);
    result_by_q.set(c_Key1(2), common_term_1);
    // q = 3
    result_by_lpq.set(c_Key3(4, 4, 3), common_term_0);
    result_by_q.set(c_Key1(3), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 4), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l4_e5(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 4.
    // Functions are taylor expanded and truncated at eccentricity^5.

    // Functions Calculated for l = 4.

    c_IntMap<c_Key3, double> result_by_lpq(43);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(5);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double common_term_0 = -0.020833333333333333*eccentricity_3;
    double common_term_1 = -0.33333333333333333*eccentricity_4 + 0.5*eccentricity_2;
    double common_term_2 = 4.6875*eccentricity_3 - 1.5*eccentricity;
    double common_term_3 = 24.875*eccentricity_4 - 11.0*eccentricity_2 + 1.0;
    double common_term_4 = -47.8125*eccentricity_3 + 6.5*eccentricity;
    double common_term_5 = -160.5*eccentricity_4 + 25.5*eccentricity_2;
    double common_term_6 = 78.145833333333333*eccentricity_3;
    double common_term_7 = 205.95833333333333*eccentricity_4;
    double common_term_8 = 1.3958333333333333*eccentricity_4;
    double common_term_9 = 1.0208333333333333*eccentricity_3;
    double common_term_10 = 0.75*eccentricity_2*std::pow(1.0 - eccentricity_2, -3.5);
    double common_term_11 = 2.0625*eccentricity_3 + 0.5*eccentricity;
    double common_term_12 = 4.0625*eccentricity_4 + eccentricity_2 + 1.0;
    double common_term_13 = -0.1875*eccentricity_3 + 4.5*eccentricity;
    double common_term_14 = -7.4583333333333333*eccentricity_4 + 13.25*eccentricity_2;
    double common_term_15 = 32.104166666666667*eccentricity_3;
    double common_term_16 = 69.375*eccentricity_4;
    double common_term_17 = 15.520833333333333*eccentricity_4;
    double common_term_18 = 9.0625*eccentricity_3;
    double common_term_19 = 12.916666666666667*eccentricity_4 + 5.0*eccentricity_2;
    double common_term_20 = 8.4375*eccentricity_3 + 2.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(9);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (4, 0).
    // q = -3
    result_by_lpq.set(c_Key3(4, 0, -3), common_term_0);
    result_by_q.set(c_Key1(-3), common_term_0);
    // q = -2
    result_by_lpq.set(c_Key3(4, 0, -2), common_term_1);
    result_by_q.set(c_Key1(-2), common_term_1);
    // q = -1
    result_by_lpq.set(c_Key3(4, 0, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    result_by_lpq.set(c_Key3(4, 0, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(4, 0, 1), common_term_4);
    result_by_q.set(c_Key1(1), common_term_4);
    // q = 2
    result_by_lpq.set(c_Key3(4, 0, 2), common_term_5);
    result_by_q.set(c_Key1(2), common_term_5);
    // q = 3
    result_by_lpq.set(c_Key3(4, 0, 3), common_term_6);
    result_by_q.set(c_Key1(3), common_term_6);
    // q = 4
    result_by_lpq.set(c_Key3(4, 0, 4), common_term_7);
    result_by_q.set(c_Key1(4), common_term_7);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 0), result_by_q);
    result_by_q.clear();

    // l , p = (4, 1).
    // q = -4
    result_by_lpq.set(c_Key3(4, 1, -4), common_term_8);
    result_by_q.set(c_Key1(-4), common_term_8);
    // q = -3
    result_by_lpq.set(c_Key3(4, 1, -3), common_term_9);
    result_by_q.set(c_Key1(-3), common_term_9);
    // q = -2
    result_by_lpq.set(c_Key3(4, 1, -2), common_term_10);
    result_by_q.set(c_Key1(-2), common_term_10);
    // q = -1
    result_by_lpq.set(c_Key3(4, 1, -1), common_term_11);
    result_by_q.set(c_Key1(-1), common_term_11);
    // q = 0
    result_by_lpq.set(c_Key3(4, 1, 0), common_term_12);
    result_by_q.set(c_Key1(0), common_term_12);
    // q = 1
    result_by_lpq.set(c_Key3(4, 1, 1), common_term_13);
    result_by_q.set(c_Key1(1), common_term_13);
    // q = 2
    result_by_lpq.set(c_Key3(4, 1, 2), common_term_14);
    result_by_q.set(c_Key1(2), common_term_14);
    // q = 3
    result_by_lpq.set(c_Key3(4, 1, 3), common_term_15);
    result_by_q.set(c_Key1(3), common_term_15);
    // q = 4
    result_by_lpq.set(c_Key3(4, 1, 4), common_term_16);
    result_by_q.set(c_Key1(4), common_term_16);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 1), result_by_q);
    result_by_q.clear();

    // l , p = (4, 2).
    // q = -4
    result_by_lpq.set(c_Key3(4, 2, -4), common_term_17);
    result_by_q.set(c_Key1(-4), common_term_17);
    // q = -3
    result_by_lpq.set(c_Key3(4, 2, -3), common_term_18);
    result_by_q.set(c_Key1(-3), common_term_18);
    // q = -2
    result_by_lpq.set(c_Key3(4, 2, -2), common_term_19);
    result_by_q.set(c_Key1(-2), common_term_19);
    // q = -1
    result_by_lpq.set(c_Key3(4, 2, -1), common_term_20);
    result_by_q.set(c_Key1(-1), common_term_20);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -3.5)*(1.5*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(4, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(4, 2, 1), common_term_20);
    result_by_q.set(c_Key1(1), common_term_20);
    // q = 2
    result_by_lpq.set(c_Key3(4, 2, 2), common_term_19);
    result_by_q.set(c_Key1(2), common_term_19);
    // q = 3
    result_by_lpq.set(c_Key3(4, 2, 3), common_term_18);
    result_by_q.set(c_Key1(3), common_term_18);
    // q = 4
    result_by_lpq.set(c_Key3(4, 2, 4), common_term_17);
    result_by_q.set(c_Key1(4), common_term_17);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 2), result_by_q);
    result_by_q.clear();

    // l , p = (4, 3).
    // q = -4
    result_by_lpq.set(c_Key3(4, 3, -4), common_term_16);
    result_by_q.set(c_Key1(-4), common_term_16);
    // q = -3
    result_by_lpq.set(c_Key3(4, 3, -3), common_term_15);
    result_by_q.set(c_Key1(-3), common_term_15);
    // q = -2
    result_by_lpq.set(c_Key3(4, 3, -2), common_term_14);
    result_by_q.set(c_Key1(-2), common_term_14);
    // q = -1
    result_by_lpq.set(c_Key3(4, 3, -1), common_term_13);
    result_by_q.set(c_Key1(-1), common_term_13);
    // q = 0
    result_by_lpq.set(c_Key3(4, 3, 0), common_term_12);
    result_by_q.set(c_Key1(0), common_term_12);
    // q = 1
    result_by_lpq.set(c_Key3(4, 3, 1), common_term_11);
    result_by_q.set(c_Key1(1), common_term_11);
    // q = 2
    result_by_lpq.set(c_Key3(4, 3, 2), common_term_10);
    result_by_q.set(c_Key1(2), common_term_10);
    // q = 3
    result_by_lpq.set(c_Key3(4, 3, 3), common_term_9);
    result_by_q.set(c_Key1(3), common_term_9);
    // q = 4
    result_by_lpq.set(c_Key3(4, 3, 4), common_term_8);
    result_by_q.set(c_Key1(4), common_term_8);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 3), result_by_q);
    result_by_q.clear();

    // l , p = (4, 4).
    // q = -4
    result_by_lpq.set(c_Key3(4, 4, -4), common_term_7);
    result_by_q.set(c_Key1(-4), common_term_7);
    // q = -3
    result_by_lpq.set(c_Key3(4, 4, -3), common_term_6);
    result_by_q.set(c_Key1(-3), common_term_6);
    // q = -2
    result_by_lpq.set(c_Key3(4, 4, -2), common_term_5);
    result_by_q.set(c_Key1(-2), common_term_5);
    // q = -1
    result_by_lpq.set(c_Key3(4, 4, -1), common_term_4);
    result_by_q.set(c_Key1(-1), common_term_4);
    // q = 0
    result_by_lpq.set(c_Key3(4, 4, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(4, 4, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // q = 2
    result_by_lpq.set(c_Key3(4, 4, 2), common_term_1);
    result_by_q.set(c_Key1(2), common_term_1);
    // q = 3
    result_by_lpq.set(c_Key3(4, 4, 3), common_term_0);
    result_by_q.set(c_Key1(3), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 4), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l4_e10(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 4.
    // Functions are taylor expanded and truncated at eccentricity^10.

    // Functions Calculated for l = 4.

    c_IntMap<c_Key3, double> result_by_lpq(93);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(5);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_8 = eccentricity_4 * eccentricity_4;
    double eccentricity_9 = eccentricity_8 * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double eccentricity_6 = eccentricity_3 * eccentricity_3;
    double eccentricity_7 = eccentricity_6 * eccentricity;
    double eccentricity_5 = eccentricity_4 * eccentricity;
    double common_term_0 = 0.010512283029169422*eccentricity_9;
    double common_term_1 = 0.0063492063492063492*eccentricity_8;
    double common_term_2 = 0.0049791608537946429*eccentricity_9 + 0.0033900669642857143*eccentricity_7;
    double common_term_3 = 0.0019841269841269841*eccentricity_8 + 0.0013888888888888889*eccentricity_6;
    double common_term_4 = 0.00033656529017857143*eccentricity_9 + 0.00033637152777777778*eccentricity_7 + 0.00026041666666666667*eccentricity_5;
    double common_term_5 = -0.0052295826099537037*eccentricity_9 - 0.0068033854166666667*eccentricity_7 - 0.0091145833333333333*eccentricity_5 - 0.020833333333333333*eccentricity_3;
    double common_term_6 = -0.030555555555555556*eccentricity_8 - 0.33333333333333333*eccentricity_4 + 0.5*eccentricity_2;
    double common_term_7 = -0.154852294921875*eccentricity_9 + 0.61083984375*eccentricity_7 - 3.0703125*eccentricity_5 + 4.6875*eccentricity_3 - 1.5*eccentricity;
    double common_term_8 = 5.7834201388888889*eccentricity_8 - 18.194444444444444*eccentricity_6 + 24.875*eccentricity_4 - 11.0*eccentricity_2 + 1.0;
    double common_term_9 = 33.458201090494792*eccentricity_9 - 81.000705295138889*eccentricity_7 + 98.763020833333333*eccentricity_5 - 47.8125*eccentricity_3 + 6.5*eccentricity;
    double common_term_10 = -297.39375*eccentricity_8 + 326.625*eccentricity_6 - 160.5*eccentricity_4 + 25.5*eccentricity_2;
    double common_term_11 = -951.35203224464699*eccentricity_9 + 950.95387369791667*eccentricity_7 - 459.52994791666667*eccentricity_5 + 78.145833333333333*eccentricity_3;
    double common_term_12 = 2518.9965277777778*eccentricity_8 - 1179.55*eccentricity_6 + 205.95833333333333*eccentricity_4;
    double common_term_13 = 6202.1901803152902*eccentricity_9 - 2794.46572265625*eccentricity_7 + 489.84609375*eccentricity_5;
    double common_term_14 = -6224.7792658730159*eccentricity_8 + 1081.1180555555556*eccentricity_6;
    double common_term_15 = -13203.441354079474*eccentricity_9 + 2253.9047913566468*eccentricity_7;
    double common_term_16 = 4492.4152901785714*eccentricity_8;
    double common_term_17 = 8634.2892587866943*eccentricity_9;
    double common_term_18 = 6.8567855646580826*eccentricity_9;
    double common_term_19 = 4.9818080357142857*eccentricity_8;
    double common_term_20 = 7.861205328078497*eccentricity_9 + 3.6200164310515873*eccentricity_7;
    double common_term_21 = 6.4266121031746032*eccentricity_8 + 2.6315972222222222*eccentricity_6;
    double common_term_22 = 10.252205984933036*eccentricity_9 + 5.19462890625*eccentricity_7 + 1.91484375*eccentricity_5;
    double common_term_23 = 8.5769965277777778*eccentricity_8 + 4.1625*eccentricity_6 + 1.3958333333333333*eccentricity_4;
    double common_term_24 = 12.705226869936343*eccentricity_9 + 7.1361002604166667*eccentricity_7 + 3.3138020833333333*eccentricity_5 + 1.0208333333333333*eccentricity_3;
    double common_term_25 = 0.75*eccentricity_2*std::pow(1.0 - eccentricity_2, -3.5);
    double common_term_26 = 15.251617431640625*eccentricity_9 + 9.1812608506944444*eccentricity_7 + 4.8567708333333333*eccentricity_5 + 2.0625*eccentricity_3 + 0.5*eccentricity;
    double common_term_27 = 13.155815972222222*eccentricity_8 + 7.7222222222222222*eccentricity_6 + 4.0625*eccentricity_4 + eccentricity_2 + 1.0;
    double common_term_28 = 17.867779541015625*eccentricity_9 + 11.10205078125*eccentricity_7 + 7.5234375*eccentricity_5 - 0.1875*eccentricity_3 + 4.5*eccentricity;
    double common_term_29 = 14.158506944444444*eccentricity_8 + 15.265625*eccentricity_6 - 7.4583333333333333*eccentricity_4 + 13.25*eccentricity_2;
    double common_term_30 = 13.520844636140046*eccentricity_9 + 35.305501302083333*eccentricity_7 - 30.571614583333333*eccentricity_5 + 32.104166666666667*eccentricity_3;
    double common_term_31 = 87.3796875*eccentricity_8 - 89.3625*eccentricity_6 + 69.375*eccentricity_4;
    double common_term_32 = 215.74904378255208*eccentricity_9 - 221.89308810763889*eccentricity_7 + 138.96276041666667*eccentricity_5;
    double common_term_33 = -498.2624503968254*eccentricity_8 + 263.67430555555556*eccentricity_6;
    double common_term_34 = -1043.61390511649*eccentricity_9 + 480.35352957589286*eccentricity_7;
    double common_term_35 = 847.75489831349206*eccentricity_8;
    double common_term_36 = 1458.5491959775356*eccentricity_9;
    double common_term_37 = 154.08794512067522*eccentricity_9;
    double common_term_38 = 100.57707093253968*eccentricity_8;
    double common_term_39 = 32.994084676106771*eccentricity_9 + 64.837944878472222*eccentricity_7;
    double common_term_40 = 34.647321428571429*eccentricity_8 + 41.15625*eccentricity_6;
    double common_term_41 = 56.933056059337798*eccentricity_9 + 30.660536024305556*eccentricity_7 + 25.610677083333333*eccentricity_5;
    double common_term_42 = 44.741319444444444*eccentricity_8 + 24.625*eccentricity_6 + 15.520833333333333*eccentricity_4;
    double common_term_43 = 56.60687255859375*eccentricity_9 + 34.61572265625*eccentricity_7 + 18.41796875*eccentricity_5 + 9.0625*eccentricity_3;
    double common_term_44 = 44.975694444444444*eccentricity_8 + 26.09375*eccentricity_6 + 12.916666666666667*eccentricity_4 + 5.0*eccentricity_2;
    double common_term_45 = 56.878041585286458*eccentricity_9 + 34.885796440972222*eccentricity_7 + 18.971354166666667*eccentricity_5 + 8.4375*eccentricity_3 + 2.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(19);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (4, 0).
    // q = -9
    result_by_lpq.set(c_Key3(4, 0, -9), common_term_0);
    result_by_q.set(c_Key1(-9), common_term_0);
    // q = -8
    result_by_lpq.set(c_Key3(4, 0, -8), common_term_1);
    result_by_q.set(c_Key1(-8), common_term_1);
    // q = -7
    result_by_lpq.set(c_Key3(4, 0, -7), common_term_2);
    result_by_q.set(c_Key1(-7), common_term_2);
    // q = -6
    result_by_lpq.set(c_Key3(4, 0, -6), common_term_3);
    result_by_q.set(c_Key1(-6), common_term_3);
    // q = -5
    result_by_lpq.set(c_Key3(4, 0, -5), common_term_4);
    result_by_q.set(c_Key1(-5), common_term_4);
    // q = -3
    result_by_lpq.set(c_Key3(4, 0, -3), common_term_5);
    result_by_q.set(c_Key1(-3), common_term_5);
    // q = -2
    result_by_lpq.set(c_Key3(4, 0, -2), common_term_6);
    result_by_q.set(c_Key1(-2), common_term_6);
    // q = -1
    result_by_lpq.set(c_Key3(4, 0, -1), common_term_7);
    result_by_q.set(c_Key1(-1), common_term_7);
    // q = 0
    result_by_lpq.set(c_Key3(4, 0, 0), common_term_8);
    result_by_q.set(c_Key1(0), common_term_8);
    // q = 1
    result_by_lpq.set(c_Key3(4, 0, 1), common_term_9);
    result_by_q.set(c_Key1(1), common_term_9);
    // q = 2
    result_by_lpq.set(c_Key3(4, 0, 2), common_term_10);
    result_by_q.set(c_Key1(2), common_term_10);
    // q = 3
    result_by_lpq.set(c_Key3(4, 0, 3), common_term_11);
    result_by_q.set(c_Key1(3), common_term_11);
    // q = 4
    result_by_lpq.set(c_Key3(4, 0, 4), common_term_12);
    result_by_q.set(c_Key1(4), common_term_12);
    // q = 5
    result_by_lpq.set(c_Key3(4, 0, 5), common_term_13);
    result_by_q.set(c_Key1(5), common_term_13);
    // q = 6
    result_by_lpq.set(c_Key3(4, 0, 6), common_term_14);
    result_by_q.set(c_Key1(6), common_term_14);
    // q = 7
    result_by_lpq.set(c_Key3(4, 0, 7), common_term_15);
    result_by_q.set(c_Key1(7), common_term_15);
    // q = 8
    result_by_lpq.set(c_Key3(4, 0, 8), common_term_16);
    result_by_q.set(c_Key1(8), common_term_16);
    // q = 9
    result_by_lpq.set(c_Key3(4, 0, 9), common_term_17);
    result_by_q.set(c_Key1(9), common_term_17);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 0), result_by_q);
    result_by_q.clear();

    // l , p = (4, 1).
    // q = -9
    result_by_lpq.set(c_Key3(4, 1, -9), common_term_18);
    result_by_q.set(c_Key1(-9), common_term_18);
    // q = -8
    result_by_lpq.set(c_Key3(4, 1, -8), common_term_19);
    result_by_q.set(c_Key1(-8), common_term_19);
    // q = -7
    result_by_lpq.set(c_Key3(4, 1, -7), common_term_20);
    result_by_q.set(c_Key1(-7), common_term_20);
    // q = -6
    result_by_lpq.set(c_Key3(4, 1, -6), common_term_21);
    result_by_q.set(c_Key1(-6), common_term_21);
    // q = -5
    result_by_lpq.set(c_Key3(4, 1, -5), common_term_22);
    result_by_q.set(c_Key1(-5), common_term_22);
    // q = -4
    result_by_lpq.set(c_Key3(4, 1, -4), common_term_23);
    result_by_q.set(c_Key1(-4), common_term_23);
    // q = -3
    result_by_lpq.set(c_Key3(4, 1, -3), common_term_24);
    result_by_q.set(c_Key1(-3), common_term_24);
    // q = -2
    result_by_lpq.set(c_Key3(4, 1, -2), common_term_25);
    result_by_q.set(c_Key1(-2), common_term_25);
    // q = -1
    result_by_lpq.set(c_Key3(4, 1, -1), common_term_26);
    result_by_q.set(c_Key1(-1), common_term_26);
    // q = 0
    result_by_lpq.set(c_Key3(4, 1, 0), common_term_27);
    result_by_q.set(c_Key1(0), common_term_27);
    // q = 1
    result_by_lpq.set(c_Key3(4, 1, 1), common_term_28);
    result_by_q.set(c_Key1(1), common_term_28);
    // q = 2
    result_by_lpq.set(c_Key3(4, 1, 2), common_term_29);
    result_by_q.set(c_Key1(2), common_term_29);
    // q = 3
    result_by_lpq.set(c_Key3(4, 1, 3), common_term_30);
    result_by_q.set(c_Key1(3), common_term_30);
    // q = 4
    result_by_lpq.set(c_Key3(4, 1, 4), common_term_31);
    result_by_q.set(c_Key1(4), common_term_31);
    // q = 5
    result_by_lpq.set(c_Key3(4, 1, 5), common_term_32);
    result_by_q.set(c_Key1(5), common_term_32);
    // q = 6
    result_by_lpq.set(c_Key3(4, 1, 6), common_term_33);
    result_by_q.set(c_Key1(6), common_term_33);
    // q = 7
    result_by_lpq.set(c_Key3(4, 1, 7), common_term_34);
    result_by_q.set(c_Key1(7), common_term_34);
    // q = 8
    result_by_lpq.set(c_Key3(4, 1, 8), common_term_35);
    result_by_q.set(c_Key1(8), common_term_35);
    // q = 9
    result_by_lpq.set(c_Key3(4, 1, 9), common_term_36);
    result_by_q.set(c_Key1(9), common_term_36);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 1), result_by_q);
    result_by_q.clear();

    // l , p = (4, 2).
    // q = -9
    result_by_lpq.set(c_Key3(4, 2, -9), common_term_37);
    result_by_q.set(c_Key1(-9), common_term_37);
    // q = -8
    result_by_lpq.set(c_Key3(4, 2, -8), common_term_38);
    result_by_q.set(c_Key1(-8), common_term_38);
    // q = -7
    result_by_lpq.set(c_Key3(4, 2, -7), common_term_39);
    result_by_q.set(c_Key1(-7), common_term_39);
    // q = -6
    result_by_lpq.set(c_Key3(4, 2, -6), common_term_40);
    result_by_q.set(c_Key1(-6), common_term_40);
    // q = -5
    result_by_lpq.set(c_Key3(4, 2, -5), common_term_41);
    result_by_q.set(c_Key1(-5), common_term_41);
    // q = -4
    result_by_lpq.set(c_Key3(4, 2, -4), common_term_42);
    result_by_q.set(c_Key1(-4), common_term_42);
    // q = -3
    result_by_lpq.set(c_Key3(4, 2, -3), common_term_43);
    result_by_q.set(c_Key1(-3), common_term_43);
    // q = -2
    result_by_lpq.set(c_Key3(4, 2, -2), common_term_44);
    result_by_q.set(c_Key1(-2), common_term_44);
    // q = -1
    result_by_lpq.set(c_Key3(4, 2, -1), common_term_45);
    result_by_q.set(c_Key1(-1), common_term_45);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -3.5)*(1.5*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(4, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(4, 2, 1), common_term_45);
    result_by_q.set(c_Key1(1), common_term_45);
    // q = 2
    result_by_lpq.set(c_Key3(4, 2, 2), common_term_44);
    result_by_q.set(c_Key1(2), common_term_44);
    // q = 3
    result_by_lpq.set(c_Key3(4, 2, 3), common_term_43);
    result_by_q.set(c_Key1(3), common_term_43);
    // q = 4
    result_by_lpq.set(c_Key3(4, 2, 4), common_term_42);
    result_by_q.set(c_Key1(4), common_term_42);
    // q = 5
    result_by_lpq.set(c_Key3(4, 2, 5), common_term_41);
    result_by_q.set(c_Key1(5), common_term_41);
    // q = 6
    result_by_lpq.set(c_Key3(4, 2, 6), common_term_40);
    result_by_q.set(c_Key1(6), common_term_40);
    // q = 7
    result_by_lpq.set(c_Key3(4, 2, 7), common_term_39);
    result_by_q.set(c_Key1(7), common_term_39);
    // q = 8
    result_by_lpq.set(c_Key3(4, 2, 8), common_term_38);
    result_by_q.set(c_Key1(8), common_term_38);
    // q = 9
    result_by_lpq.set(c_Key3(4, 2, 9), common_term_37);
    result_by_q.set(c_Key1(9), common_term_37);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 2), result_by_q);
    result_by_q.clear();

    // l , p = (4, 3).
    // q = -9
    result_by_lpq.set(c_Key3(4, 3, -9), common_term_36);
    result_by_q.set(c_Key1(-9), common_term_36);
    // q = -8
    result_by_lpq.set(c_Key3(4, 3, -8), common_term_35);
    result_by_q.set(c_Key1(-8), common_term_35);
    // q = -7
    result_by_lpq.set(c_Key3(4, 3, -7), common_term_34);
    result_by_q.set(c_Key1(-7), common_term_34);
    // q = -6
    result_by_lpq.set(c_Key3(4, 3, -6), common_term_33);
    result_by_q.set(c_Key1(-6), common_term_33);
    // q = -5
    result_by_lpq.set(c_Key3(4, 3, -5), common_term_32);
    result_by_q.set(c_Key1(-5), common_term_32);
    // q = -4
    result_by_lpq.set(c_Key3(4, 3, -4), common_term_31);
    result_by_q.set(c_Key1(-4), common_term_31);
    // q = -3
    result_by_lpq.set(c_Key3(4, 3, -3), common_term_30);
    result_by_q.set(c_Key1(-3), common_term_30);
    // q = -2
    result_by_lpq.set(c_Key3(4, 3, -2), common_term_29);
    result_by_q.set(c_Key1(-2), common_term_29);
    // q = -1
    result_by_lpq.set(c_Key3(4, 3, -1), common_term_28);
    result_by_q.set(c_Key1(-1), common_term_28);
    // q = 0
    result_by_lpq.set(c_Key3(4, 3, 0), common_term_27);
    result_by_q.set(c_Key1(0), common_term_27);
    // q = 1
    result_by_lpq.set(c_Key3(4, 3, 1), common_term_26);
    result_by_q.set(c_Key1(1), common_term_26);
    // q = 2
    result_by_lpq.set(c_Key3(4, 3, 2), common_term_25);
    result_by_q.set(c_Key1(2), common_term_25);
    // q = 3
    result_by_lpq.set(c_Key3(4, 3, 3), common_term_24);
    result_by_q.set(c_Key1(3), common_term_24);
    // q = 4
    result_by_lpq.set(c_Key3(4, 3, 4), common_term_23);
    result_by_q.set(c_Key1(4), common_term_23);
    // q = 5
    result_by_lpq.set(c_Key3(4, 3, 5), common_term_22);
    result_by_q.set(c_Key1(5), common_term_22);
    // q = 6
    result_by_lpq.set(c_Key3(4, 3, 6), common_term_21);
    result_by_q.set(c_Key1(6), common_term_21);
    // q = 7
    result_by_lpq.set(c_Key3(4, 3, 7), common_term_20);
    result_by_q.set(c_Key1(7), common_term_20);
    // q = 8
    result_by_lpq.set(c_Key3(4, 3, 8), common_term_19);
    result_by_q.set(c_Key1(8), common_term_19);
    // q = 9
    result_by_lpq.set(c_Key3(4, 3, 9), common_term_18);
    result_by_q.set(c_Key1(9), common_term_18);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 3), result_by_q);
    result_by_q.clear();

    // l , p = (4, 4).
    // q = -9
    result_by_lpq.set(c_Key3(4, 4, -9), common_term_17);
    result_by_q.set(c_Key1(-9), common_term_17);
    // q = -8
    result_by_lpq.set(c_Key3(4, 4, -8), common_term_16);
    result_by_q.set(c_Key1(-8), common_term_16);
    // q = -7
    result_by_lpq.set(c_Key3(4, 4, -7), common_term_15);
    result_by_q.set(c_Key1(-7), common_term_15);
    // q = -6
    result_by_lpq.set(c_Key3(4, 4, -6), common_term_14);
    result_by_q.set(c_Key1(-6), common_term_14);
    // q = -5
    result_by_lpq.set(c_Key3(4, 4, -5), common_term_13);
    result_by_q.set(c_Key1(-5), common_term_13);
    // q = -4
    result_by_lpq.set(c_Key3(4, 4, -4), common_term_12);
    result_by_q.set(c_Key1(-4), common_term_12);
    // q = -3
    result_by_lpq.set(c_Key3(4, 4, -3), common_term_11);
    result_by_q.set(c_Key1(-3), common_term_11);
    // q = -2
    result_by_lpq.set(c_Key3(4, 4, -2), common_term_10);
    result_by_q.set(c_Key1(-2), common_term_10);
    // q = -1
    result_by_lpq.set(c_Key3(4, 4, -1), common_term_9);
    result_by_q.set(c_Key1(-1), common_term_9);
    // q = 0
    result_by_lpq.set(c_Key3(4, 4, 0), common_term_8);
    result_by_q.set(c_Key1(0), common_term_8);
    // q = 1
    result_by_lpq.set(c_Key3(4, 4, 1), common_term_7);
    result_by_q.set(c_Key1(1), common_term_7);
    // q = 2
    result_by_lpq.set(c_Key3(4, 4, 2), common_term_6);
    result_by_q.set(c_Key1(2), common_term_6);
    // q = 3
    result_by_lpq.set(c_Key3(4, 4, 3), common_term_5);
    result_by_q.set(c_Key1(3), common_term_5);
    // q = 5
    result_by_lpq.set(c_Key3(4, 4, 5), common_term_4);
    result_by_q.set(c_Key1(5), common_term_4);
    // q = 6
    result_by_lpq.set(c_Key3(4, 4, 6), common_term_3);
    result_by_q.set(c_Key1(6), common_term_3);
    // q = 7
    result_by_lpq.set(c_Key3(4, 4, 7), common_term_2);
    result_by_q.set(c_Key1(7), common_term_2);
    // q = 8
    result_by_lpq.set(c_Key3(4, 4, 8), common_term_1);
    result_by_q.set(c_Key1(8), common_term_1);
    // q = 9
    result_by_lpq.set(c_Key3(4, 4, 9), common_term_0);
    result_by_q.set(c_Key1(9), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 4), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l4_e15(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 4.
    // Functions are taylor expanded and truncated at eccentricity^15.

    // Functions Calculated for l = 4.

    c_IntMap<c_Key3, double> result_by_lpq(143);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(5);
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
    double common_term_0 = 0.070011874986143339*eccentricity_14;
    double common_term_1 = 0.049829021929265498*eccentricity_13;
    double common_term_2 = 0.035025386136497248*eccentricity_14 + 0.035025386136497248*eccentricity_12;
    double common_term_3 = 0.027714957818082331*eccentricity_13 + 0.024187599550326398*eccentricity_11;
    double common_term_4 = 0.02477830762987013*eccentricity_14 + 0.020710227272727273*eccentricity_12 + 0.016272321428571429*eccentricity_10;
    double common_term_5 = 0.017037663261764646*eccentricity_13 + 0.014454389165107956*eccentricity_11 + 0.010512283029169422*eccentricity_9;
    double common_term_6 = 0.011087595532039976*eccentricity_14 + 0.010634920634920635*eccentricity_12 + 0.0091710758377425044*eccentricity_10 + 0.0063492063492063492*eccentricity_8;
    double common_term_7 = 0.0057276902879987444*eccentricity_13 + 0.0056280408586774554*eccentricity_11 + 0.0049791608537946429*eccentricity_9 + 0.0033900669642857143*eccentricity_7;
    double common_term_8 = 0.0019841327252351558*eccentricity_14 + 0.0021049199000587889*eccentricity_12 + 0.0021453373015873016*eccentricity_10 + 0.0019841269841269841*eccentricity_8 + 0.0013888888888888889*eccentricity_6;
    double common_term_9 = 0.00028329432396101012*eccentricity_13 + 0.00031263361532221395*eccentricity_11 + 0.00033656529017857143*eccentricity_9 + 0.00033637152777777778*eccentricity_7 + 0.00026041666666666667*eccentricity_5;
    double common_term_10 = -0.0034640831669802388*eccentricity_13 - 0.0041940557893621858*eccentricity_11 - 0.0052295826099537037*eccentricity_9 - 0.0068033854166666667*eccentricity_7 - 0.0091145833333333333*eccentricity_5 - 0.020833333333333333*eccentricity_3;
    double common_term_11 = -0.014401730599647266*eccentricity_14 - 0.017333829365079365*eccentricity_12 - 0.021368634259259259*eccentricity_10 - 0.030555555555555556*eccentricity_8 - 0.33333333333333333*eccentricity_4 + 0.5*eccentricity_2;
    double common_term_12 = -0.040518174852643694*eccentricity_13 - 0.042028656005859375*eccentricity_11 - 0.154852294921875*eccentricity_9 + 0.61083984375*eccentricity_7 - 3.0703125*eccentricity_5 + 4.6875*eccentricity_3 - 1.5*eccentricity;
    double common_term_13 = -0.085807311350466112*eccentricity_14 + 0.055096691743827161*eccentricity_12 - 1.2580034722222222*eccentricity_10 + 5.7834201388888889*eccentricity_8 - 18.194444444444444*eccentricity_6 + 24.875*eccentricity_4 - 11.0*eccentricity_2 + 1.0;
    double common_term_14 = 1.3447142614468906*eccentricity_13 - 8.7808166786476418*eccentricity_11 + 33.458201090494792*eccentricity_9 - 81.000705295138889*eccentricity_7 + 98.763020833333333*eccentricity_5 - 47.8125*eccentricity_3 + 6.5*eccentricity;
    double common_term_15 = 10.165484095982143*eccentricity_14 - 47.400033482142857*eccentricity_12 + 148.210546875*eccentricity_10 - 297.39375*eccentricity_8 + 326.625*eccentricity_6 - 160.5*eccentricity_4 + 25.5*eccentricity_2;
    double common_term_16 = -208.7227470004117*eccentricity_13 + 550.15952060840748*eccentricity_11 - 951.35203224464699*eccentricity_9 + 950.95387369791667*eccentricity_7 - 459.52994791666667*eccentricity_5 + 78.145833333333333*eccentricity_3;
    double common_term_17 = -788.61598664388595*eccentricity_14 + 1797.0557074652778*eccentricity_12 - 2744.5086970899471*eccentricity_10 + 2518.9965277777778*eccentricity_8 - 1179.55*eccentricity_6 + 205.95833333333333*eccentricity_4;
    double common_term_18 = 5324.4495974779129*eccentricity_13 - 7305.2158461979457*eccentricity_11 + 6202.1901803152902*eccentricity_9 - 2794.46572265625*eccentricity_7 + 489.84609375*eccentricity_5;
    double common_term_19 = 14604.241807008056*eccentricity_14 - 18229.565656461273*eccentricity_12 + 14405.240110367063*eccentricity_10 - 6224.7792658730159*eccentricity_8 + 1081.1180555555556*eccentricity_6;
    double common_term_20 = -43144.437509774074*eccentricity_13 + 31898.767145541694*eccentricity_11 - 13203.441354079474*eccentricity_9 + 2253.9047913566468*eccentricity_7;
    double common_term_21 = -97689.773614042208*eccentricity_14 + 67881.609400111607*eccentricity_12 - 26910.148325892857*eccentricity_10 + 4492.4152901785714*eccentricity_8;
    double common_term_22 = 139668.49558447324*eccentricity_13 - 53054.607669626737*eccentricity_11 + 8634.2892587866943*eccentricity_9;
    double common_term_23 = 279181.24935633096*eccentricity_14 - 101702.49863162879*eccentricity_12 + 16103.458149112654*eccentricity_10;
    double common_term_24 = -190318.12322846481*eccentricity_13 + 29284.679306401835*eccentricity_11;
    double common_term_25 = -348782.98930667648*eccentricity_14 + 52120.692975678996*eccentricity_12;
    double common_term_26 = 91056.936978657089*eccentricity_13;
    double common_term_27 = 156525.67958167251*eccentricity_14;
    double common_term_28 = 33.794073491031848*eccentricity_14;
    double common_term_29 = 24.575391528512721*eccentricity_13;
    double common_term_30 = 14.926281510293012*eccentricity_14 + 17.866780259711032*eccentricity_12;
    double common_term_31 = 14.289571168267882*eccentricity_13 + 12.986322418807389*eccentricity_11;
    double common_term_32 = 23.187307042659148*eccentricity_14 + 12.894632250706469*eccentricity_12 + 9.4370307677469136*eccentricity_10;
    double common_term_33 = 19.822783763630757*eccentricity_13 + 11.201183952990873*eccentricity_11 + 6.8567855646580826*eccentricity_9;
    double common_term_34 = 26.351643795657468*eccentricity_14 + 16.925938895089286*eccentricity_12 + 9.4754464285714286*eccentricity_10 + 4.9818080357142857*eccentricity_8;
    double common_term_35 = 22.990468170557353*eccentricity_13 + 14.397755378978803*eccentricity_11 + 7.861205328078497*eccentricity_9 + 3.6200164310515873*eccentricity_7;
    double common_term_36 = 29.945886675404633*eccentricity_14 + 19.962523653365667*eccentricity_12 + 12.18421378968254*eccentricity_10 + 6.4266121031746032*eccentricity_8 + 2.6315972222222222*eccentricity_6;
    double common_term_37 = 26.370534416607448*eccentricity_13 + 17.25169061933245*eccentricity_11 + 10.252205984933036*eccentricity_9 + 5.19462890625*eccentricity_7 + 1.91484375*eccentricity_5;
    double common_term_38 = 33.627440177652851*eccentricity_14 + 23.133182973710317*eccentricity_12 + 14.83926917989418*eccentricity_10 + 8.5769965277777778*eccentricity_8 + 4.1625*eccentricity_6 + 1.3958333333333333*eccentricity_4;
    double common_term_39 = 29.841378941636868*eccentricity_13 + 20.214907857097646*eccentricity_11 + 12.705226869936343*eccentricity_9 + 7.1361002604166667*eccentricity_7 + 3.3138020833333333*eccentricity_5 + 1.0208333333333333*eccentricity_3;
    double common_term_40 = 0.75*eccentricity_2*std::pow(1.0 - eccentricity_2, -3.5);
    double common_term_41 = 33.38909522085173*eccentricity_13 + 23.262085916024667*eccentricity_11 + 15.251617431640625*eccentricity_9 + 9.1812608506944444*eccentricity_7 + 4.8567708333333333*eccentricity_5 + 2.0625*eccentricity_3 + 0.5*eccentricity;
    double common_term_42 = 41.216691825022833*eccentricity_14 + 29.720034481095679*eccentricity_12 + 20.422673611111111*eccentricity_10 + 13.155815972222222*eccentricity_8 + 7.7222222222222222*eccentricity_6 + 4.0625*eccentricity_4 + eccentricity_2 + 1.0;
    double common_term_43 = 36.988072577885219*eccentricity_13 + 26.360151214599609*eccentricity_11 + 17.867779541015625*eccentricity_9 + 11.10205078125*eccentricity_7 + 7.5234375*eccentricity_5 - 0.1875*eccentricity_3 + 4.5*eccentricity;
    double common_term_44 = 45.08621452425733*eccentricity_14 + 33.066267671130952*eccentricity_12 + 23.493453414351852*eccentricity_10 + 14.158506944444444*eccentricity_8 + 15.265625*eccentricity_6 - 7.4583333333333333*eccentricity_4 + 13.25*eccentricity_2;
    double common_term_45 = 40.415334258129988*eccentricity_13 + 30.918859915758567*eccentricity_11 + 13.520844636140046*eccentricity_9 + 35.305501302083333*eccentricity_7 - 30.571614583333333*eccentricity_5 + 32.104166666666667*eccentricity_3;
    double common_term_46 = 47.635421316964286*eccentricity_14 + 43.692543247767857*eccentricity_12 - 1.39453125*eccentricity_10 + 87.3796875*eccentricity_8 - 89.3625*eccentricity_6 + 69.375*eccentricity_4;
    double common_term_47 = 73.611281795168119*eccentricity_13 - 59.243043023568613*eccentricity_11 + 215.74904378255208*eccentricity_9 - 221.89308810763889*eccentricity_7 + 138.96276041666667*eccentricity_5;
    double common_term_48 = 154.96848346618258*eccentricity_14 - 230.21997073183054*eccentricity_12 + 513.34901103670635*eccentricity_10 - 498.2624503968254*eccentricity_8 + 263.67430555555556*eccentricity_6;
    double common_term_49 = -674.21670128549848*eccentricity_13 + 1165.9634992871966*eccentricity_11 - 1043.61390511649*eccentricity_9 + 480.35352957589286*eccentricity_7;
    double common_term_50 = -1735.7767232666864*eccentricity_14 + 2531.0599901310351*eccentricity_12 - 2075.9276551477072*eccentricity_10 + 847.75489831349206*eccentricity_8;
    double common_term_51 = 5275.6588546423769*eccentricity_13 - 3967.35798898781*eccentricity_11 + 1458.5491959775356*eccentricity_9;
    double common_term_52 = 10613.364179338728*eccentricity_14 - 7342.7522947950487*eccentricity_12 + 2457.5253850446429*eccentricity_10;
    double common_term_53 = -13236.455289503308*eccentricity_13 + 4069.044575450592*eccentricity_11;
    double common_term_54 = -23339.876825966721*eccentricity_14 + 6638.2568172757669*eccentricity_12;
    double common_term_55 = 10692.719224305181*eccentricity_13;
    double common_term_56 = 17034.159056928928*eccentricity_14;
    double common_term_57 = 1147.5537296915525*eccentricity_14;
    double common_term_58 = 777.86994939862003*eccentricity_13;
    double common_term_59 = -518.67570593469031*eccentricity_14 + 524.41106369216721*eccentricity_12;
    double common_term_60 = -247.08419345608184*eccentricity_13 + 351.33040815417308*eccentricity_11;
    double common_term_61 = 271.50264519601813*eccentricity_14 - 96.335783899260462*eccentricity_12 + 233.66851989638448*eccentricity_10;
    double common_term_62 = 177.89107343562238*eccentricity_13 - 17.670006288800921*eccentricity_11 + 154.08794512067522*eccentricity_9;
    double common_term_63 = 132.52199727516568*eccentricity_14 + 125.02479693700397*eccentricity_12 + 19.254719190917108*eccentricity_10 + 100.57707093253968*eccentricity_8;
    double common_term_64 = 117.61316365351893*eccentricity_13 + 93.258982458232362*eccentricity_11 + 32.994084676106771*eccentricity_9 + 64.837944878472222*eccentricity_7;
    double common_term_65 = 141.55147530691964*eccentricity_14 + 100.88992745535714*eccentricity_12 + 72.286272321428571*eccentricity_10 + 34.647321428571429*eccentricity_8 + 41.15625*eccentricity_6;
    double common_term_66 = 120.86578909086522*eccentricity_13 + 84.684386228127454*eccentricity_11 + 56.933056059337798*eccentricity_9 + 30.660536024305556*eccentricity_7 + 25.610677083333333*eccentricity_5;
    double common_term_67 = 141.99537978578777*eccentricity_14 + 102.16700381324405*eccentricity_12 + 69.845320767195767*eccentricity_10 + 44.741319444444444*eccentricity_8 + 24.625*eccentricity_6 + 15.520833333333333*eccentricity_4;
    double common_term_68 = 121.30688092027392*eccentricity_13 + 85.32438714163644*eccentricity_11 + 56.60687255859375*eccentricity_9 + 34.61572265625*eccentricity_7 + 18.41796875*eccentricity_5 + 9.0625*eccentricity_3;
    double common_term_69 = 142.32977555424658*eccentricity_14 + 102.5221912202381*eccentricity_12 + 70.253761574074074*eccentricity_10 + 44.975694444444444*eccentricity_8 + 26.09375*eccentricity_6 + 12.916666666666667*eccentricity_4 + 5.0*eccentricity_2;
    double common_term_70 = 121.53733695718136*eccentricity_13 + 85.571842306631583*eccentricity_11 + 56.878041585286458*eccentricity_9 + 34.885796440972222*eccentricity_7 + 18.971354166666667*eccentricity_5 + 8.4375*eccentricity_3 + 2.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(29);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (4, 0).
    // q = -14
    result_by_lpq.set(c_Key3(4, 0, -14), common_term_0);
    result_by_q.set(c_Key1(-14), common_term_0);
    // q = -13
    result_by_lpq.set(c_Key3(4, 0, -13), common_term_1);
    result_by_q.set(c_Key1(-13), common_term_1);
    // q = -12
    result_by_lpq.set(c_Key3(4, 0, -12), common_term_2);
    result_by_q.set(c_Key1(-12), common_term_2);
    // q = -11
    result_by_lpq.set(c_Key3(4, 0, -11), common_term_3);
    result_by_q.set(c_Key1(-11), common_term_3);
    // q = -10
    result_by_lpq.set(c_Key3(4, 0, -10), common_term_4);
    result_by_q.set(c_Key1(-10), common_term_4);
    // q = -9
    result_by_lpq.set(c_Key3(4, 0, -9), common_term_5);
    result_by_q.set(c_Key1(-9), common_term_5);
    // q = -8
    result_by_lpq.set(c_Key3(4, 0, -8), common_term_6);
    result_by_q.set(c_Key1(-8), common_term_6);
    // q = -7
    result_by_lpq.set(c_Key3(4, 0, -7), common_term_7);
    result_by_q.set(c_Key1(-7), common_term_7);
    // q = -6
    result_by_lpq.set(c_Key3(4, 0, -6), common_term_8);
    result_by_q.set(c_Key1(-6), common_term_8);
    // q = -5
    result_by_lpq.set(c_Key3(4, 0, -5), common_term_9);
    result_by_q.set(c_Key1(-5), common_term_9);
    // q = -3
    result_by_lpq.set(c_Key3(4, 0, -3), common_term_10);
    result_by_q.set(c_Key1(-3), common_term_10);
    // q = -2
    result_by_lpq.set(c_Key3(4, 0, -2), common_term_11);
    result_by_q.set(c_Key1(-2), common_term_11);
    // q = -1
    result_by_lpq.set(c_Key3(4, 0, -1), common_term_12);
    result_by_q.set(c_Key1(-1), common_term_12);
    // q = 0
    result_by_lpq.set(c_Key3(4, 0, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(4, 0, 1), common_term_14);
    result_by_q.set(c_Key1(1), common_term_14);
    // q = 2
    result_by_lpq.set(c_Key3(4, 0, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // q = 3
    result_by_lpq.set(c_Key3(4, 0, 3), common_term_16);
    result_by_q.set(c_Key1(3), common_term_16);
    // q = 4
    result_by_lpq.set(c_Key3(4, 0, 4), common_term_17);
    result_by_q.set(c_Key1(4), common_term_17);
    // q = 5
    result_by_lpq.set(c_Key3(4, 0, 5), common_term_18);
    result_by_q.set(c_Key1(5), common_term_18);
    // q = 6
    result_by_lpq.set(c_Key3(4, 0, 6), common_term_19);
    result_by_q.set(c_Key1(6), common_term_19);
    // q = 7
    result_by_lpq.set(c_Key3(4, 0, 7), common_term_20);
    result_by_q.set(c_Key1(7), common_term_20);
    // q = 8
    result_by_lpq.set(c_Key3(4, 0, 8), common_term_21);
    result_by_q.set(c_Key1(8), common_term_21);
    // q = 9
    result_by_lpq.set(c_Key3(4, 0, 9), common_term_22);
    result_by_q.set(c_Key1(9), common_term_22);
    // q = 10
    result_by_lpq.set(c_Key3(4, 0, 10), common_term_23);
    result_by_q.set(c_Key1(10), common_term_23);
    // q = 11
    result_by_lpq.set(c_Key3(4, 0, 11), common_term_24);
    result_by_q.set(c_Key1(11), common_term_24);
    // q = 12
    result_by_lpq.set(c_Key3(4, 0, 12), common_term_25);
    result_by_q.set(c_Key1(12), common_term_25);
    // q = 13
    result_by_lpq.set(c_Key3(4, 0, 13), common_term_26);
    result_by_q.set(c_Key1(13), common_term_26);
    // q = 14
    result_by_lpq.set(c_Key3(4, 0, 14), common_term_27);
    result_by_q.set(c_Key1(14), common_term_27);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 0), result_by_q);
    result_by_q.clear();

    // l , p = (4, 1).
    // q = -14
    result_by_lpq.set(c_Key3(4, 1, -14), common_term_28);
    result_by_q.set(c_Key1(-14), common_term_28);
    // q = -13
    result_by_lpq.set(c_Key3(4, 1, -13), common_term_29);
    result_by_q.set(c_Key1(-13), common_term_29);
    // q = -12
    result_by_lpq.set(c_Key3(4, 1, -12), common_term_30);
    result_by_q.set(c_Key1(-12), common_term_30);
    // q = -11
    result_by_lpq.set(c_Key3(4, 1, -11), common_term_31);
    result_by_q.set(c_Key1(-11), common_term_31);
    // q = -10
    result_by_lpq.set(c_Key3(4, 1, -10), common_term_32);
    result_by_q.set(c_Key1(-10), common_term_32);
    // q = -9
    result_by_lpq.set(c_Key3(4, 1, -9), common_term_33);
    result_by_q.set(c_Key1(-9), common_term_33);
    // q = -8
    result_by_lpq.set(c_Key3(4, 1, -8), common_term_34);
    result_by_q.set(c_Key1(-8), common_term_34);
    // q = -7
    result_by_lpq.set(c_Key3(4, 1, -7), common_term_35);
    result_by_q.set(c_Key1(-7), common_term_35);
    // q = -6
    result_by_lpq.set(c_Key3(4, 1, -6), common_term_36);
    result_by_q.set(c_Key1(-6), common_term_36);
    // q = -5
    result_by_lpq.set(c_Key3(4, 1, -5), common_term_37);
    result_by_q.set(c_Key1(-5), common_term_37);
    // q = -4
    result_by_lpq.set(c_Key3(4, 1, -4), common_term_38);
    result_by_q.set(c_Key1(-4), common_term_38);
    // q = -3
    result_by_lpq.set(c_Key3(4, 1, -3), common_term_39);
    result_by_q.set(c_Key1(-3), common_term_39);
    // q = -2
    result_by_lpq.set(c_Key3(4, 1, -2), common_term_40);
    result_by_q.set(c_Key1(-2), common_term_40);
    // q = -1
    result_by_lpq.set(c_Key3(4, 1, -1), common_term_41);
    result_by_q.set(c_Key1(-1), common_term_41);
    // q = 0
    result_by_lpq.set(c_Key3(4, 1, 0), common_term_42);
    result_by_q.set(c_Key1(0), common_term_42);
    // q = 1
    result_by_lpq.set(c_Key3(4, 1, 1), common_term_43);
    result_by_q.set(c_Key1(1), common_term_43);
    // q = 2
    result_by_lpq.set(c_Key3(4, 1, 2), common_term_44);
    result_by_q.set(c_Key1(2), common_term_44);
    // q = 3
    result_by_lpq.set(c_Key3(4, 1, 3), common_term_45);
    result_by_q.set(c_Key1(3), common_term_45);
    // q = 4
    result_by_lpq.set(c_Key3(4, 1, 4), common_term_46);
    result_by_q.set(c_Key1(4), common_term_46);
    // q = 5
    result_by_lpq.set(c_Key3(4, 1, 5), common_term_47);
    result_by_q.set(c_Key1(5), common_term_47);
    // q = 6
    result_by_lpq.set(c_Key3(4, 1, 6), common_term_48);
    result_by_q.set(c_Key1(6), common_term_48);
    // q = 7
    result_by_lpq.set(c_Key3(4, 1, 7), common_term_49);
    result_by_q.set(c_Key1(7), common_term_49);
    // q = 8
    result_by_lpq.set(c_Key3(4, 1, 8), common_term_50);
    result_by_q.set(c_Key1(8), common_term_50);
    // q = 9
    result_by_lpq.set(c_Key3(4, 1, 9), common_term_51);
    result_by_q.set(c_Key1(9), common_term_51);
    // q = 10
    result_by_lpq.set(c_Key3(4, 1, 10), common_term_52);
    result_by_q.set(c_Key1(10), common_term_52);
    // q = 11
    result_by_lpq.set(c_Key3(4, 1, 11), common_term_53);
    result_by_q.set(c_Key1(11), common_term_53);
    // q = 12
    result_by_lpq.set(c_Key3(4, 1, 12), common_term_54);
    result_by_q.set(c_Key1(12), common_term_54);
    // q = 13
    result_by_lpq.set(c_Key3(4, 1, 13), common_term_55);
    result_by_q.set(c_Key1(13), common_term_55);
    // q = 14
    result_by_lpq.set(c_Key3(4, 1, 14), common_term_56);
    result_by_q.set(c_Key1(14), common_term_56);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 1), result_by_q);
    result_by_q.clear();

    // l , p = (4, 2).
    // q = -14
    result_by_lpq.set(c_Key3(4, 2, -14), common_term_57);
    result_by_q.set(c_Key1(-14), common_term_57);
    // q = -13
    result_by_lpq.set(c_Key3(4, 2, -13), common_term_58);
    result_by_q.set(c_Key1(-13), common_term_58);
    // q = -12
    result_by_lpq.set(c_Key3(4, 2, -12), common_term_59);
    result_by_q.set(c_Key1(-12), common_term_59);
    // q = -11
    result_by_lpq.set(c_Key3(4, 2, -11), common_term_60);
    result_by_q.set(c_Key1(-11), common_term_60);
    // q = -10
    result_by_lpq.set(c_Key3(4, 2, -10), common_term_61);
    result_by_q.set(c_Key1(-10), common_term_61);
    // q = -9
    result_by_lpq.set(c_Key3(4, 2, -9), common_term_62);
    result_by_q.set(c_Key1(-9), common_term_62);
    // q = -8
    result_by_lpq.set(c_Key3(4, 2, -8), common_term_63);
    result_by_q.set(c_Key1(-8), common_term_63);
    // q = -7
    result_by_lpq.set(c_Key3(4, 2, -7), common_term_64);
    result_by_q.set(c_Key1(-7), common_term_64);
    // q = -6
    result_by_lpq.set(c_Key3(4, 2, -6), common_term_65);
    result_by_q.set(c_Key1(-6), common_term_65);
    // q = -5
    result_by_lpq.set(c_Key3(4, 2, -5), common_term_66);
    result_by_q.set(c_Key1(-5), common_term_66);
    // q = -4
    result_by_lpq.set(c_Key3(4, 2, -4), common_term_67);
    result_by_q.set(c_Key1(-4), common_term_67);
    // q = -3
    result_by_lpq.set(c_Key3(4, 2, -3), common_term_68);
    result_by_q.set(c_Key1(-3), common_term_68);
    // q = -2
    result_by_lpq.set(c_Key3(4, 2, -2), common_term_69);
    result_by_q.set(c_Key1(-2), common_term_69);
    // q = -1
    result_by_lpq.set(c_Key3(4, 2, -1), common_term_70);
    result_by_q.set(c_Key1(-1), common_term_70);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -3.5)*(1.5*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(4, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(4, 2, 1), common_term_70);
    result_by_q.set(c_Key1(1), common_term_70);
    // q = 2
    result_by_lpq.set(c_Key3(4, 2, 2), common_term_69);
    result_by_q.set(c_Key1(2), common_term_69);
    // q = 3
    result_by_lpq.set(c_Key3(4, 2, 3), common_term_68);
    result_by_q.set(c_Key1(3), common_term_68);
    // q = 4
    result_by_lpq.set(c_Key3(4, 2, 4), common_term_67);
    result_by_q.set(c_Key1(4), common_term_67);
    // q = 5
    result_by_lpq.set(c_Key3(4, 2, 5), common_term_66);
    result_by_q.set(c_Key1(5), common_term_66);
    // q = 6
    result_by_lpq.set(c_Key3(4, 2, 6), common_term_65);
    result_by_q.set(c_Key1(6), common_term_65);
    // q = 7
    result_by_lpq.set(c_Key3(4, 2, 7), common_term_64);
    result_by_q.set(c_Key1(7), common_term_64);
    // q = 8
    result_by_lpq.set(c_Key3(4, 2, 8), common_term_63);
    result_by_q.set(c_Key1(8), common_term_63);
    // q = 9
    result_by_lpq.set(c_Key3(4, 2, 9), common_term_62);
    result_by_q.set(c_Key1(9), common_term_62);
    // q = 10
    result_by_lpq.set(c_Key3(4, 2, 10), common_term_61);
    result_by_q.set(c_Key1(10), common_term_61);
    // q = 11
    result_by_lpq.set(c_Key3(4, 2, 11), common_term_60);
    result_by_q.set(c_Key1(11), common_term_60);
    // q = 12
    result_by_lpq.set(c_Key3(4, 2, 12), common_term_59);
    result_by_q.set(c_Key1(12), common_term_59);
    // q = 13
    result_by_lpq.set(c_Key3(4, 2, 13), common_term_58);
    result_by_q.set(c_Key1(13), common_term_58);
    // q = 14
    result_by_lpq.set(c_Key3(4, 2, 14), common_term_57);
    result_by_q.set(c_Key1(14), common_term_57);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 2), result_by_q);
    result_by_q.clear();

    // l , p = (4, 3).
    // q = -14
    result_by_lpq.set(c_Key3(4, 3, -14), common_term_56);
    result_by_q.set(c_Key1(-14), common_term_56);
    // q = -13
    result_by_lpq.set(c_Key3(4, 3, -13), common_term_55);
    result_by_q.set(c_Key1(-13), common_term_55);
    // q = -12
    result_by_lpq.set(c_Key3(4, 3, -12), common_term_54);
    result_by_q.set(c_Key1(-12), common_term_54);
    // q = -11
    result_by_lpq.set(c_Key3(4, 3, -11), common_term_53);
    result_by_q.set(c_Key1(-11), common_term_53);
    // q = -10
    result_by_lpq.set(c_Key3(4, 3, -10), common_term_52);
    result_by_q.set(c_Key1(-10), common_term_52);
    // q = -9
    result_by_lpq.set(c_Key3(4, 3, -9), common_term_51);
    result_by_q.set(c_Key1(-9), common_term_51);
    // q = -8
    result_by_lpq.set(c_Key3(4, 3, -8), common_term_50);
    result_by_q.set(c_Key1(-8), common_term_50);
    // q = -7
    result_by_lpq.set(c_Key3(4, 3, -7), common_term_49);
    result_by_q.set(c_Key1(-7), common_term_49);
    // q = -6
    result_by_lpq.set(c_Key3(4, 3, -6), common_term_48);
    result_by_q.set(c_Key1(-6), common_term_48);
    // q = -5
    result_by_lpq.set(c_Key3(4, 3, -5), common_term_47);
    result_by_q.set(c_Key1(-5), common_term_47);
    // q = -4
    result_by_lpq.set(c_Key3(4, 3, -4), common_term_46);
    result_by_q.set(c_Key1(-4), common_term_46);
    // q = -3
    result_by_lpq.set(c_Key3(4, 3, -3), common_term_45);
    result_by_q.set(c_Key1(-3), common_term_45);
    // q = -2
    result_by_lpq.set(c_Key3(4, 3, -2), common_term_44);
    result_by_q.set(c_Key1(-2), common_term_44);
    // q = -1
    result_by_lpq.set(c_Key3(4, 3, -1), common_term_43);
    result_by_q.set(c_Key1(-1), common_term_43);
    // q = 0
    result_by_lpq.set(c_Key3(4, 3, 0), common_term_42);
    result_by_q.set(c_Key1(0), common_term_42);
    // q = 1
    result_by_lpq.set(c_Key3(4, 3, 1), common_term_41);
    result_by_q.set(c_Key1(1), common_term_41);
    // q = 2
    result_by_lpq.set(c_Key3(4, 3, 2), common_term_40);
    result_by_q.set(c_Key1(2), common_term_40);
    // q = 3
    result_by_lpq.set(c_Key3(4, 3, 3), common_term_39);
    result_by_q.set(c_Key1(3), common_term_39);
    // q = 4
    result_by_lpq.set(c_Key3(4, 3, 4), common_term_38);
    result_by_q.set(c_Key1(4), common_term_38);
    // q = 5
    result_by_lpq.set(c_Key3(4, 3, 5), common_term_37);
    result_by_q.set(c_Key1(5), common_term_37);
    // q = 6
    result_by_lpq.set(c_Key3(4, 3, 6), common_term_36);
    result_by_q.set(c_Key1(6), common_term_36);
    // q = 7
    result_by_lpq.set(c_Key3(4, 3, 7), common_term_35);
    result_by_q.set(c_Key1(7), common_term_35);
    // q = 8
    result_by_lpq.set(c_Key3(4, 3, 8), common_term_34);
    result_by_q.set(c_Key1(8), common_term_34);
    // q = 9
    result_by_lpq.set(c_Key3(4, 3, 9), common_term_33);
    result_by_q.set(c_Key1(9), common_term_33);
    // q = 10
    result_by_lpq.set(c_Key3(4, 3, 10), common_term_32);
    result_by_q.set(c_Key1(10), common_term_32);
    // q = 11
    result_by_lpq.set(c_Key3(4, 3, 11), common_term_31);
    result_by_q.set(c_Key1(11), common_term_31);
    // q = 12
    result_by_lpq.set(c_Key3(4, 3, 12), common_term_30);
    result_by_q.set(c_Key1(12), common_term_30);
    // q = 13
    result_by_lpq.set(c_Key3(4, 3, 13), common_term_29);
    result_by_q.set(c_Key1(13), common_term_29);
    // q = 14
    result_by_lpq.set(c_Key3(4, 3, 14), common_term_28);
    result_by_q.set(c_Key1(14), common_term_28);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 3), result_by_q);
    result_by_q.clear();

    // l , p = (4, 4).
    // q = -14
    result_by_lpq.set(c_Key3(4, 4, -14), common_term_27);
    result_by_q.set(c_Key1(-14), common_term_27);
    // q = -13
    result_by_lpq.set(c_Key3(4, 4, -13), common_term_26);
    result_by_q.set(c_Key1(-13), common_term_26);
    // q = -12
    result_by_lpq.set(c_Key3(4, 4, -12), common_term_25);
    result_by_q.set(c_Key1(-12), common_term_25);
    // q = -11
    result_by_lpq.set(c_Key3(4, 4, -11), common_term_24);
    result_by_q.set(c_Key1(-11), common_term_24);
    // q = -10
    result_by_lpq.set(c_Key3(4, 4, -10), common_term_23);
    result_by_q.set(c_Key1(-10), common_term_23);
    // q = -9
    result_by_lpq.set(c_Key3(4, 4, -9), common_term_22);
    result_by_q.set(c_Key1(-9), common_term_22);
    // q = -8
    result_by_lpq.set(c_Key3(4, 4, -8), common_term_21);
    result_by_q.set(c_Key1(-8), common_term_21);
    // q = -7
    result_by_lpq.set(c_Key3(4, 4, -7), common_term_20);
    result_by_q.set(c_Key1(-7), common_term_20);
    // q = -6
    result_by_lpq.set(c_Key3(4, 4, -6), common_term_19);
    result_by_q.set(c_Key1(-6), common_term_19);
    // q = -5
    result_by_lpq.set(c_Key3(4, 4, -5), common_term_18);
    result_by_q.set(c_Key1(-5), common_term_18);
    // q = -4
    result_by_lpq.set(c_Key3(4, 4, -4), common_term_17);
    result_by_q.set(c_Key1(-4), common_term_17);
    // q = -3
    result_by_lpq.set(c_Key3(4, 4, -3), common_term_16);
    result_by_q.set(c_Key1(-3), common_term_16);
    // q = -2
    result_by_lpq.set(c_Key3(4, 4, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(4, 4, -1), common_term_14);
    result_by_q.set(c_Key1(-1), common_term_14);
    // q = 0
    result_by_lpq.set(c_Key3(4, 4, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(4, 4, 1), common_term_12);
    result_by_q.set(c_Key1(1), common_term_12);
    // q = 2
    result_by_lpq.set(c_Key3(4, 4, 2), common_term_11);
    result_by_q.set(c_Key1(2), common_term_11);
    // q = 3
    result_by_lpq.set(c_Key3(4, 4, 3), common_term_10);
    result_by_q.set(c_Key1(3), common_term_10);
    // q = 5
    result_by_lpq.set(c_Key3(4, 4, 5), common_term_9);
    result_by_q.set(c_Key1(5), common_term_9);
    // q = 6
    result_by_lpq.set(c_Key3(4, 4, 6), common_term_8);
    result_by_q.set(c_Key1(6), common_term_8);
    // q = 7
    result_by_lpq.set(c_Key3(4, 4, 7), common_term_7);
    result_by_q.set(c_Key1(7), common_term_7);
    // q = 8
    result_by_lpq.set(c_Key3(4, 4, 8), common_term_6);
    result_by_q.set(c_Key1(8), common_term_6);
    // q = 9
    result_by_lpq.set(c_Key3(4, 4, 9), common_term_5);
    result_by_q.set(c_Key1(9), common_term_5);
    // q = 10
    result_by_lpq.set(c_Key3(4, 4, 10), common_term_4);
    result_by_q.set(c_Key1(10), common_term_4);
    // q = 11
    result_by_lpq.set(c_Key3(4, 4, 11), common_term_3);
    result_by_q.set(c_Key1(11), common_term_3);
    // q = 12
    result_by_lpq.set(c_Key3(4, 4, 12), common_term_2);
    result_by_q.set(c_Key1(12), common_term_2);
    // q = 13
    result_by_lpq.set(c_Key3(4, 4, 13), common_term_1);
    result_by_q.set(c_Key1(13), common_term_1);
    // q = 14
    result_by_lpq.set(c_Key3(4, 4, 14), common_term_0);
    result_by_q.set(c_Key1(14), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 4), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l4_e20(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 4.
    // Functions are taylor expanded and truncated at eccentricity^20.

    // Functions Calculated for l = 4.

    c_IntMap<c_Key3, double> result_by_lpq(193);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(5);
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
    double common_term_0 = 0.34759168852955555*eccentricity_19;
    double common_term_1 = 0.25434529016223454*eccentricity_18;
    double common_term_2 = 0.018039483189853507*eccentricity_19 + 0.18554896995277893*eccentricity_17;
    double common_term_3 = 0.03965714957311596*eccentricity_18 + 0.13483430854859426*eccentricity_16;
    double common_term_4 = 0.091583206872246788*eccentricity_19 + 0.047219616017291483*eccentricity_17 + 0.097485658874408223*eccentricity_15;
    double common_term_5 = 0.072199746079460319*eccentricity_18 + 0.04667458332409556*eccentricity_16 + 0.070011874986143339*eccentricity_14;
    double common_term_6 = 0.05954144073346099*eccentricity_19 + 0.057014188930672978*eccentricity_17 + 0.041820786262062114*eccentricity_15 + 0.049829021929265498*eccentricity_13;
    double common_term_7 = 0.047327037323509987*eccentricity_18 + 0.044503409582774662*eccentricity_16 + 0.035025386136497248*eccentricity_14 + 0.035025386136497248*eccentricity_12;
    double common_term_8 = 0.037206138900130733*eccentricity_19 + 0.036250914486980488*eccentricity_17 + 0.033882989514309395*eccentricity_15 + 0.027714957818082331*eccentricity_13 + 0.024187599550326398*eccentricity_11;
    double common_term_9 = 0.02692442575909024*eccentricity_18 + 0.026456746378621379*eccentricity_16 + 0.02477830762987013*eccentricity_14 + 0.020710227272727273*eccentricity_12 + 0.016272321428571429*eccentricity_10;
    double common_term_10 = 0.017752301290754539*eccentricity_19 + 0.018151244851032595*eccentricity_17 + 0.0180376035360927*eccentricity_15 + 0.017037663261764646*eccentricity_13 + 0.014454389165107956*eccentricity_11 + 0.010512283029169422*eccentricity_9;
    double common_term_11 = 0.010600084843140399*eccentricity_18 + 0.010987314725740652*eccentricity_16 + 0.011087595532039976*eccentricity_14 + 0.010634920634920635*eccentricity_12 + 0.0091710758377425044*eccentricity_10 + 0.0063492063492063492*eccentricity_8;
    double common_term_12 = 0.0049358435204554048*eccentricity_19 + 0.0052676870850967122*eccentricity_17 + 0.0055574100745188725*eccentricity_15 + 0.0057276902879987444*eccentricity_13 + 0.0056280408586774554*eccentricity_11 + 0.0049791608537946429*eccentricity_9 + 0.0033900669642857143*eccentricity_7;
    double common_term_13 = 0.0016913935072044761*eccentricity_18 + 0.0018382791414475442*eccentricity_16 + 0.0019841327252351558*eccentricity_14 + 0.0021049199000587889*eccentricity_12 + 0.0021453373015873016*eccentricity_10 + 0.0019841269841269841*eccentricity_8 + 0.0013888888888888889*eccentricity_6;
    double common_term_14 = 0.00020686434965186921*eccentricity_19 + 0.00022929149984314218*eccentricity_17 + 0.00025487830252033483*eccentricity_15 + 0.00028329432396101012*eccentricity_13 + 0.00031263361532221395*eccentricity_11 + 0.00033656529017857143*eccentricity_9 + 0.00033637152777777778*eccentricity_7 + 0.00026041666666666667*eccentricity_5;
    double common_term_15 = -0.002191463474462603*eccentricity_19 - 0.002514204690018195*eccentricity_17 - 0.0029255062409435887*eccentricity_15 - 0.0034640831669802388*eccentricity_13 - 0.0041940557893621858*eccentricity_11 - 0.0052295826099537037*eccentricity_9 - 0.0068033854166666667*eccentricity_7 - 0.0091145833333333333*eccentricity_5 - 0.020833333333333333*eccentricity_3;
    double common_term_16 = -0.010573057549360408*eccentricity_18 - 0.012235103902028919*eccentricity_16 - 0.014401730599647266*eccentricity_14 - 0.017333829365079365*eccentricity_12 - 0.021368634259259259*eccentricity_10 - 0.030555555555555556*eccentricity_8 - 0.33333333333333333*eccentricity_4 + 0.5*eccentricity_2;
    double common_term_17 = -0.024711372190791809*eccentricity_19 - 0.028496080233658455*eccentricity_17 - 0.033405401134977535*eccentricity_15 - 0.040518174852643694*eccentricity_13 - 0.042028656005859375*eccentricity_11 - 0.154852294921875*eccentricity_9 + 0.61083984375*eccentricity_7 - 3.0703125*eccentricity_5 + 4.6875*eccentricity_3 - 1.5*eccentricity;
    double common_term_18 = -0.052068951078972327*eccentricity_18 - 0.05997888512682272*eccentricity_16 - 0.085807311350466112*eccentricity_14 + 0.055096691743827161*eccentricity_12 - 1.2580034722222222*eccentricity_10 + 5.7834201388888889*eccentricity_8 - 18.194444444444444*eccentricity_6 + 24.875*eccentricity_4 - 11.0*eccentricity_2 + 1.0;
    double common_term_19 = -0.084272027012395294*eccentricity_19 - 0.079110496370550202*eccentricity_17 - 0.30139996413245673*eccentricity_15 + 1.3447142614468906*eccentricity_13 - 8.7808166786476418*eccentricity_11 + 33.458201090494792*eccentricity_9 - 81.000705295138889*eccentricity_7 + 98.763020833333333*eccentricity_5 - 47.8125*eccentricity_3 + 6.5*eccentricity;
    double common_term_20 = 0.067820352359693878*eccentricity_18 - 1.8467147640306122*eccentricity_16 + 10.165484095982143*eccentricity_14 - 47.400033482142857*eccentricity_12 + 148.210546875*eccentricity_10 - 297.39375*eccentricity_8 + 326.625*eccentricity_6 - 160.5*eccentricity_4 + 25.5*eccentricity_2;
    double common_term_21 = 1.5313601147591457*eccentricity_19 - 11.356243537181669*eccentricity_17 + 55.460092323955522*eccentricity_15 - 208.7227470004117*eccentricity_13 + 550.15952060840748*eccentricity_11 - 951.35203224464699*eccentricity_9 + 950.95387369791667*eccentricity_7 - 459.52994791666667*eccentricity_5 + 78.145833333333333*eccentricity_3;
    double common_term_22 = -59.362379564225732*eccentricity_18 + 247.90848736573462*eccentricity_16 - 788.61598664388595*eccentricity_14 + 1797.0557074652778*eccentricity_12 - 2744.5086970899471*eccentricity_10 + 2518.9965277777778*eccentricity_8 - 1179.55*eccentricity_6 + 205.95833333333333*eccentricity_4;
    double common_term_23 = -265.69890860753507*eccentricity_19 + 959.25634577632899*eccentricity_17 - 2649.4054827976738*eccentricity_15 + 5324.4495974779129*eccentricity_13 - 7305.2158461979457*eccentricity_11 + 6202.1901803152902*eccentricity_9 - 2794.46572265625*eccentricity_7 + 489.84609375*eccentricity_5;
    double common_term_24 = 3320.3409847825183*eccentricity_18 - 8113.2781843526201*eccentricity_16 + 14604.241807008056*eccentricity_14 - 18229.565656461273*eccentricity_12 + 14405.240110367063*eccentricity_10 - 6224.7792658730159*eccentricity_8 + 1081.1180555555556*eccentricity_6;
    double common_term_25 = 10512.409199683637*eccentricity_19 - 23052.058798035141*eccentricity_17 + 37622.37391105637*eccentricity_15 - 43144.437509774074*eccentricity_13 + 31898.767145541694*eccentricity_11 - 13203.441354079474*eccentricity_9 + 2253.9047913566468*eccentricity_7;
    double common_term_26 = -61568.428132416778*eccentricity_18 + 92002.489285119787*eccentricity_16 - 97689.773614042208*eccentricity_14 + 67881.609400111607*eccentricity_12 - 26910.148325892857*eccentricity_10 + 4492.4152901785714*eccentricity_8;
    double common_term_27 = -156110.50470598864*eccentricity_19 + 215304.84190499105*eccentricity_17 - 213034.97742998696*eccentricity_15 + 139668.49558447324*eccentricity_13 - 53054.607669626737*eccentricity_11 + 8634.2892587866943*eccentricity_9;
    double common_term_28 = 485230.18555370442*eccentricity_18 - 449793.21605362163*eccentricity_16 + 279181.24935633096*eccentricity_14 - 101702.49863162879*eccentricity_12 + 16103.458149112654*eccentricity_10;
    double common_term_29 = 1058425.6598960504*eccentricity_19 - 923345.78202836643*eccentricity_17 + 544218.13637855594*eccentricity_15 - 190318.12322846481*eccentricity_13 + 29284.679306401835*eccentricity_11;
    double common_term_30 = -1849246.4050413751*eccentricity_18 + 1037780.2303366452*eccentricity_16 - 348782.98930667648*eccentricity_14 + 52120.692975678996*eccentricity_12;
    double common_term_31 = -3623545.0797993803*eccentricity_19 + 1940853.8993785101*eccentricity_17 - 627598.31024167091*eccentricity_15 + 91056.936978657089*eccentricity_13;
    double common_term_32 = 3567472.1113547689*eccentricity_18 - 1111181.6067928606*eccentricity_16 + 156525.67958167251*eccentricity_14;
    double common_term_33 = 6456371.4660353803*eccentricity_19 - 1939261.0315399098*eccentricity_17 + 265262.86346406422*eccentricity_15;
    double common_term_34 = -3341065.9665771337*eccentricity_18 + 443904.12655581547*eccentricity_16;
    double common_term_35 = -5689628.6949553253*eccentricity_19 + 734536.99141943954*eccentricity_17;
    double common_term_36 = 1203232.066313763*eccentricity_18;
    double common_term_37 = 1953093.446195968*eccentricity_19;
    double common_term_38 = 165.49685092188331*eccentricity_19;
    double common_term_39 = 120.51135552231316*eccentricity_18;
    double common_term_40 = -41.745569754810392*eccentricity_19 + 87.731374034502904*eccentricity_17;
    double common_term_41 = -13.734645668760665*eccentricity_18 + 63.850962479650147*eccentricity_16;
    double common_term_42 = 62.790027306521848*eccentricity_19 + 2.1515711696799881*eccentricity_17 + 46.458274549086854*eccentricity_15;
    double common_term_43 = 48.432635568358818*eccentricity_18 + 10.424989935183399*eccentricity_16 + 33.794073491031848*eccentricity_14;
    double common_term_44 = 47.224433860567189*eccentricity_19 + 38.924160366271301*eccentricity_17 + 14.04468899703617*eccentricity_15 + 24.575391528512721*eccentricity_13;
    double common_term_45 = 42.902415216070794*eccentricity_18 + 32.236483003661252*eccentricity_16 + 14.926281510293012*eccentricity_14 + 17.866780259711032*eccentricity_12;
    double common_term_46 = 53.557910905421605*eccentricity_19 + 38.421509298650324*eccentricity_17 + 27.207378681922531*eccentricity_15 + 14.289571168267882*eccentricity_13 + 12.986322418807389*eccentricity_11;
    double common_term_47 = 47.995376808404123*eccentricity_18 + 34.094525600625771*eccentricity_16 + 23.187307042659148*eccentricity_14 + 12.894632250706469*eccentricity_12 + 9.4370307677469136*eccentricity_10;
    double common_term_48 = 58.261648537898592*eccentricity_19 + 42.881729557807514*eccentricity_17 + 30.055493095843315*eccentricity_15 + 19.822783763630757*eccentricity_13 + 11.201183952990873*eccentricity_11 + 6.8567855646580826*eccentricity_9;
    double common_term_49 = 52.517695581761988*eccentricity_18 + 38.183190973455256*eccentricity_16 + 26.351643795657468*eccentricity_14 + 16.925938895089286*eccentricity_12 + 9.4754464285714286*eccentricity_10 + 4.9818080357142857*eccentricity_8;
    double common_term_50 = 63.09722029289597*eccentricity_19 + 47.198688357392816*eccentricity_17 + 33.877461301861368*eccentricity_15 + 22.990468170557353*eccentricity_13 + 14.397755378978803*eccentricity_11 + 7.861205328078497*eccentricity_9 + 3.6200164310515873*eccentricity_7;
    double common_term_51 = 57.13707001496127*eccentricity_18 + 42.288530259750176*eccentricity_16 + 29.945886675404633*eccentricity_14 + 19.962523653365667*eccentricity_12 + 12.18421378968254*eccentricity_10 + 6.4266121031746032*eccentricity_8 + 2.6315972222222222*eccentricity_6;
    double common_term_52 = 68.010924860559679*eccentricity_19 + 51.603097574469525*eccentricity_17 + 37.770526388219425*eccentricity_15 + 26.370534416607448*eccentricity_13 + 17.25169061933245*eccentricity_11 + 10.252205984933036*eccentricity_9 + 5.19462890625*eccentricity_7 + 1.91484375*eccentricity_5;
    double common_term_53 = 61.834575263639143*eccentricity_18 + 46.478866287413845*eccentricity_16 + 33.627440177652851*eccentricity_14 + 23.133182973710317*eccentricity_12 + 14.83926917989418*eccentricity_10 + 8.5769965277777778*eccentricity_8 + 4.1625*eccentricity_6 + 1.3958333333333333*eccentricity_4;
    double common_term_54 = 72.99665718578988*eccentricity_19 + 56.084803159772704*eccentricity_17 + 41.747373619585698*eccentricity_15 + 29.841378941636868*eccentricity_13 + 20.214907857097646*eccentricity_11 + 12.705226869936343*eccentricity_9 + 7.1361002604166667*eccentricity_7 + 3.3138020833333333*eccentricity_5 + 1.0208333333333333*eccentricity_3;
    double common_term_55 = 0.75*eccentricity_2*std::pow(1.0 - eccentricity_2, -3.5);
    double common_term_56 = 78.045271160815269*eccentricity_19 + 60.633194854714707*eccentricity_17 + 45.79547658401597*eccentricity_15 + 33.38909522085173*eccentricity_13 + 23.262085916024667*eccentricity_11 + 15.251617431640625*eccentricity_9 + 9.1812608506944444*eccentricity_7 + 4.8567708333333333*eccentricity_5 + 2.0625*eccentricity_3 + 0.5*eccentricity;
    double common_term_57 = 71.426906330552456*eccentricity_18 + 55.069871777116156*eccentricity_16 + 41.216691825022833*eccentricity_14 + 29.720034481095679*eccentricity_12 + 20.422673611111111*eccentricity_10 + 13.155815972222222*eccentricity_8 + 7.7222222222222222*eccentricity_6 + 4.0625*eccentricity_4 + eccentricity_2 + 1.0;
    double common_term_58 = 83.141477710302206*eccentricity_19 + 65.230430075775984*eccentricity_17 + 49.893653763179876*eccentricity_15 + 36.988072577885219*eccentricity_13 + 26.360151214599609*eccentricity_11 + 17.867779541015625*eccentricity_9 + 11.10205078125*eccentricity_7 + 7.5234375*eccentricity_5 - 0.1875*eccentricity_3 + 4.5*eccentricity;
    double common_term_59 = 76.293541883408878*eccentricity_18 + 59.437244162516392*eccentricity_16 + 45.08621452425733*eccentricity_14 + 33.066267671130952*eccentricity_12 + 23.493453414351852*eccentricity_10 + 14.158506944444444*eccentricity_8 + 15.265625*eccentricity_6 - 7.4583333333333333*eccentricity_4 + 13.25*eccentricity_2;
    double common_term_60 = 88.275575467188025*eccentricity_19 + 69.863104645437882*eccentricity_17 + 54.050241583708966*eccentricity_15 + 40.415334258129988*eccentricity_13 + 30.918859915758567*eccentricity_11 + 13.520844636140046*eccentricity_9 + 35.305501302083333*eccentricity_7 - 30.571614583333333*eccentricity_5 + 32.104166666666667*eccentricity_3;
    double common_term_61 = 81.17263885033772*eccentricity_18 + 64.025873866489955*eccentricity_16 + 47.635421316964286*eccentricity_14 + 43.692543247767857*eccentricity_12 - 1.39453125*eccentricity_10 + 87.3796875*eccentricity_8 - 89.3625*eccentricity_6 + 69.375*eccentricity_4;
    double common_term_62 = 93.275950507078526*eccentricity_19 + 75.712667099159543*eccentricity_17 + 51.424337092412843*eccentricity_15 + 73.611281795168119*eccentricity_13 - 59.243043023568613*eccentricity_11 + 215.74904378255208*eccentricity_9 - 221.89308810763889*eccentricity_7 + 138.96276041666667*eccentricity_5;
    double common_term_63 = 92.052192126398807*eccentricity_18 + 40.133964053329884*eccentricity_16 + 154.96848346618258*eccentricity_14 - 230.21997073183054*eccentricity_12 + 513.34901103670635*eccentricity_10 - 498.2624503968254*eccentricity_8 + 263.67430555555556*eccentricity_6;
    double common_term_64 = 123.7530766704035*eccentricity_19 - 22.245418892679954*eccentricity_17 + 378.09074706857081*eccentricity_15 - 674.21670128549848*eccentricity_13 + 1165.9634992871966*eccentricity_11 - 1043.61390511649*eccentricity_9 + 480.35352957589286*eccentricity_7;
    double common_term_65 = -236.63325025247226*eccentricity_18 + 963.67508489788752*eccentricity_16 - 1735.7767232666864*eccentricity_14 + 2531.0599901310351*eccentricity_12 - 2075.9276551477072*eccentricity_10 + 847.75489831349206*eccentricity_8;
    double common_term_66 = -865.09322459488105*eccentricity_19 + 2420.0394311237091*eccentricity_17 - 4125.4724093224525*eccentricity_15 + 5275.6588546423769*eccentricity_13 - 3967.35798898781*eccentricity_11 + 1458.5491959775356*eccentricity_9;
    double common_term_67 = 5863.4897828675187*eccentricity_18 - 9257.4324250737056*eccentricity_16 + 10613.364179338728*eccentricity_14 - 7342.7522947950487*eccentricity_12 + 2457.5253850446429*eccentricity_10;
    double common_term_68 = 13651.079613011618*eccentricity_19 - 19863.860470018095*eccentricity_17 + 20705.968653969583*eccentricity_15 - 13236.455289503308*eccentricity_13 + 4069.044575450592*eccentricity_11;
    double common_term_69 = -41093.171695668956*eccentricity_18 + 39335.328647421242*eccentricity_16 - 23339.876825966721*eccentricity_14 + 6638.2568172757669*eccentricity_12;
    double common_term_70 = -82440.918105982611*eccentricity_19 + 73014.957422817603*eccentricity_17 - 40389.55888549491*eccentricity_15 + 10692.719224305181*eccentricity_13;
    double common_term_71 = 132811.0673084307*eccentricity_18 - 68771.557278755496*eccentricity_16 + 17034.159056928928*eccentricity_14;
    double common_term_72 = 237299.66824546371*eccentricity_19 - 115457.22542067062*eccentricity_17 + 26874.646988122335*eccentricity_15;
    double common_term_73 = -191444.42741346514*eccentricity_18 + 42038.016718134318*eccentricity_16;
    double common_term_74 = -313966.39166993903*eccentricity_19 + 65256.898125446444*eccentricity_17;
    double common_term_75 = 100609.53102493417*eccentricity_18;
    double common_term_76 = 154160.47368640868*eccentricity_19;
    double common_term_77 = 7534.8219706986326*eccentricity_19;
    double common_term_78 = 5207.1267158906798*eccentricity_18;
    double common_term_79 = -8493.1219468757066*eccentricity_19 + 3587.4552273266105*eccentricity_17;
    double common_term_80 = -5164.612410260701*eccentricity_18 + 2463.1745241952132*eccentricity_16;
    double common_term_81 = 4116.0665122883993*eccentricity_19 - 3073.0369941549854*eccentricity_17 + 1684.8326059340372*eccentricity_15;
    double common_term_82 = 2332.6345609577105*eccentricity_18 - 1777.2453797773032*eccentricity_16 + 1147.5537296915525*eccentricity_14;
    double common_term_83 = -447.06291502582577*eccentricity_19 + 1319.9122569198474*eccentricity_17 - 988.3816178605511*eccentricity_15 + 777.86994939862003*eccentricity_13;
    double common_term_84 = -88.81965093430565*eccentricity_18 + 753.91437607105506*eccentricity_16 - 518.67570593469031*eccentricity_14 + 524.41106369216721*eccentricity_12;
    double common_term_85 = 331.26306900988307*eccentricity_19 + 68.767883810134465*eccentricity_17 + 441.99070737488338*eccentricity_15 - 247.08419345608184*eccentricity_13 + 351.33040815417308*eccentricity_11;
    double common_term_86 = 264.49509046288521*eccentricity_18 + 127.48921466251834*eccentricity_16 + 271.50264519601813*eccentricity_14 - 96.335783899260462*eccentricity_12 + 233.66851989638448*eccentricity_10;
    double common_term_87 = 274.72474204132075*eccentricity_19 + 222.07166870663651*eccentricity_17 + 139.98744630670393*eccentricity_15 + 177.89107343562238*eccentricity_13 - 17.670006288800921*eccentricity_11 + 154.08794512067522*eccentricity_9;
    double common_term_88 = 244.70317070665839*eccentricity_18 + 190.62598717289427*eccentricity_16 + 132.52199727516568*eccentricity_14 + 125.02479693700397*eccentricity_12 + 19.254719190917108*eccentricity_10 + 100.57707093253968*eccentricity_8;
    double common_term_89 = 277.03037965589359*eccentricity_19 + 216.0576328526854*eccentricity_17 + 164.52840318232607*eccentricity_15 + 117.61316365351893*eccentricity_13 + 93.258982458232362*eccentricity_11 + 32.994084676106771*eccentricity_9 + 64.837944878472222*eccentricity_7;
    double common_term_90 = 245.7812552632914*eccentricity_18 + 189.33475370332792*eccentricity_16 + 141.55147530691964*eccentricity_14 + 100.88992745535714*eccentricity_12 + 72.286272321428571*eccentricity_10 + 34.647321428571429*eccentricity_8 + 41.15625*eccentricity_6;
    double common_term_91 = 277.60209340728307*eccentricity_19 + 216.75026649961592*eccentricity_17 + 164.65299641678678*eccentricity_15 + 120.86578909086522*eccentricity_13 + 84.684386228127454*eccentricity_11 + 56.933056059337798*eccentricity_9 + 30.660536024305556*eccentricity_7 + 25.610677083333333*eccentricity_5;
    double common_term_92 = 246.27719816209531*eccentricity_18 + 189.87248716192528*eccentricity_16 + 141.99537978578777*eccentricity_14 + 102.16700381324405*eccentricity_12 + 69.845320767195767*eccentricity_10 + 44.741319444444444*eccentricity_8 + 24.625*eccentricity_6 + 15.520833333333333*eccentricity_4;
    double common_term_93 = 277.99076383546983*eccentricity_19 + 217.1587164314745*eccentricity_17 + 165.08652688392571*eccentricity_15 + 121.30688092027392*eccentricity_13 + 85.32438714163644*eccentricity_11 + 56.60687255859375*eccentricity_9 + 34.61572265625*eccentricity_7 + 18.41796875*eccentricity_5 + 9.0625*eccentricity_3;
    double common_term_94 = 246.57613335031881*eccentricity_18 + 190.18762872931451*eccentricity_16 + 142.32977555424658*eccentricity_14 + 102.5221912202381*eccentricity_12 + 70.253761574074074*eccentricity_10 + 44.975694444444444*eccentricity_8 + 26.09375*eccentricity_6 + 12.916666666666667*eccentricity_4 + 5.0*eccentricity_2;
    double common_term_95 = 278.18533096465024*eccentricity_19 + 217.36334160798068*eccentricity_17 + 165.30293176188806*eccentricity_15 + 121.53733695718136*eccentricity_13 + 85.571842306631583*eccentricity_11 + 56.878041585286458*eccentricity_9 + 34.885796440972222*eccentricity_7 + 18.971354166666667*eccentricity_5 + 8.4375*eccentricity_3 + 2.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(39);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (4, 0).
    // q = -19
    result_by_lpq.set(c_Key3(4, 0, -19), common_term_0);
    result_by_q.set(c_Key1(-19), common_term_0);
    // q = -18
    result_by_lpq.set(c_Key3(4, 0, -18), common_term_1);
    result_by_q.set(c_Key1(-18), common_term_1);
    // q = -17
    result_by_lpq.set(c_Key3(4, 0, -17), common_term_2);
    result_by_q.set(c_Key1(-17), common_term_2);
    // q = -16
    result_by_lpq.set(c_Key3(4, 0, -16), common_term_3);
    result_by_q.set(c_Key1(-16), common_term_3);
    // q = -15
    result_by_lpq.set(c_Key3(4, 0, -15), common_term_4);
    result_by_q.set(c_Key1(-15), common_term_4);
    // q = -14
    result_by_lpq.set(c_Key3(4, 0, -14), common_term_5);
    result_by_q.set(c_Key1(-14), common_term_5);
    // q = -13
    result_by_lpq.set(c_Key3(4, 0, -13), common_term_6);
    result_by_q.set(c_Key1(-13), common_term_6);
    // q = -12
    result_by_lpq.set(c_Key3(4, 0, -12), common_term_7);
    result_by_q.set(c_Key1(-12), common_term_7);
    // q = -11
    result_by_lpq.set(c_Key3(4, 0, -11), common_term_8);
    result_by_q.set(c_Key1(-11), common_term_8);
    // q = -10
    result_by_lpq.set(c_Key3(4, 0, -10), common_term_9);
    result_by_q.set(c_Key1(-10), common_term_9);
    // q = -9
    result_by_lpq.set(c_Key3(4, 0, -9), common_term_10);
    result_by_q.set(c_Key1(-9), common_term_10);
    // q = -8
    result_by_lpq.set(c_Key3(4, 0, -8), common_term_11);
    result_by_q.set(c_Key1(-8), common_term_11);
    // q = -7
    result_by_lpq.set(c_Key3(4, 0, -7), common_term_12);
    result_by_q.set(c_Key1(-7), common_term_12);
    // q = -6
    result_by_lpq.set(c_Key3(4, 0, -6), common_term_13);
    result_by_q.set(c_Key1(-6), common_term_13);
    // q = -5
    result_by_lpq.set(c_Key3(4, 0, -5), common_term_14);
    result_by_q.set(c_Key1(-5), common_term_14);
    // q = -3
    result_by_lpq.set(c_Key3(4, 0, -3), common_term_15);
    result_by_q.set(c_Key1(-3), common_term_15);
    // q = -2
    result_by_lpq.set(c_Key3(4, 0, -2), common_term_16);
    result_by_q.set(c_Key1(-2), common_term_16);
    // q = -1
    result_by_lpq.set(c_Key3(4, 0, -1), common_term_17);
    result_by_q.set(c_Key1(-1), common_term_17);
    // q = 0
    result_by_lpq.set(c_Key3(4, 0, 0), common_term_18);
    result_by_q.set(c_Key1(0), common_term_18);
    // q = 1
    result_by_lpq.set(c_Key3(4, 0, 1), common_term_19);
    result_by_q.set(c_Key1(1), common_term_19);
    // q = 2
    result_by_lpq.set(c_Key3(4, 0, 2), common_term_20);
    result_by_q.set(c_Key1(2), common_term_20);
    // q = 3
    result_by_lpq.set(c_Key3(4, 0, 3), common_term_21);
    result_by_q.set(c_Key1(3), common_term_21);
    // q = 4
    result_by_lpq.set(c_Key3(4, 0, 4), common_term_22);
    result_by_q.set(c_Key1(4), common_term_22);
    // q = 5
    result_by_lpq.set(c_Key3(4, 0, 5), common_term_23);
    result_by_q.set(c_Key1(5), common_term_23);
    // q = 6
    result_by_lpq.set(c_Key3(4, 0, 6), common_term_24);
    result_by_q.set(c_Key1(6), common_term_24);
    // q = 7
    result_by_lpq.set(c_Key3(4, 0, 7), common_term_25);
    result_by_q.set(c_Key1(7), common_term_25);
    // q = 8
    result_by_lpq.set(c_Key3(4, 0, 8), common_term_26);
    result_by_q.set(c_Key1(8), common_term_26);
    // q = 9
    result_by_lpq.set(c_Key3(4, 0, 9), common_term_27);
    result_by_q.set(c_Key1(9), common_term_27);
    // q = 10
    result_by_lpq.set(c_Key3(4, 0, 10), common_term_28);
    result_by_q.set(c_Key1(10), common_term_28);
    // q = 11
    result_by_lpq.set(c_Key3(4, 0, 11), common_term_29);
    result_by_q.set(c_Key1(11), common_term_29);
    // q = 12
    result_by_lpq.set(c_Key3(4, 0, 12), common_term_30);
    result_by_q.set(c_Key1(12), common_term_30);
    // q = 13
    result_by_lpq.set(c_Key3(4, 0, 13), common_term_31);
    result_by_q.set(c_Key1(13), common_term_31);
    // q = 14
    result_by_lpq.set(c_Key3(4, 0, 14), common_term_32);
    result_by_q.set(c_Key1(14), common_term_32);
    // q = 15
    result_by_lpq.set(c_Key3(4, 0, 15), common_term_33);
    result_by_q.set(c_Key1(15), common_term_33);
    // q = 16
    result_by_lpq.set(c_Key3(4, 0, 16), common_term_34);
    result_by_q.set(c_Key1(16), common_term_34);
    // q = 17
    result_by_lpq.set(c_Key3(4, 0, 17), common_term_35);
    result_by_q.set(c_Key1(17), common_term_35);
    // q = 18
    result_by_lpq.set(c_Key3(4, 0, 18), common_term_36);
    result_by_q.set(c_Key1(18), common_term_36);
    // q = 19
    result_by_lpq.set(c_Key3(4, 0, 19), common_term_37);
    result_by_q.set(c_Key1(19), common_term_37);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 0), result_by_q);
    result_by_q.clear();

    // l , p = (4, 1).
    // q = -19
    result_by_lpq.set(c_Key3(4, 1, -19), common_term_38);
    result_by_q.set(c_Key1(-19), common_term_38);
    // q = -18
    result_by_lpq.set(c_Key3(4, 1, -18), common_term_39);
    result_by_q.set(c_Key1(-18), common_term_39);
    // q = -17
    result_by_lpq.set(c_Key3(4, 1, -17), common_term_40);
    result_by_q.set(c_Key1(-17), common_term_40);
    // q = -16
    result_by_lpq.set(c_Key3(4, 1, -16), common_term_41);
    result_by_q.set(c_Key1(-16), common_term_41);
    // q = -15
    result_by_lpq.set(c_Key3(4, 1, -15), common_term_42);
    result_by_q.set(c_Key1(-15), common_term_42);
    // q = -14
    result_by_lpq.set(c_Key3(4, 1, -14), common_term_43);
    result_by_q.set(c_Key1(-14), common_term_43);
    // q = -13
    result_by_lpq.set(c_Key3(4, 1, -13), common_term_44);
    result_by_q.set(c_Key1(-13), common_term_44);
    // q = -12
    result_by_lpq.set(c_Key3(4, 1, -12), common_term_45);
    result_by_q.set(c_Key1(-12), common_term_45);
    // q = -11
    result_by_lpq.set(c_Key3(4, 1, -11), common_term_46);
    result_by_q.set(c_Key1(-11), common_term_46);
    // q = -10
    result_by_lpq.set(c_Key3(4, 1, -10), common_term_47);
    result_by_q.set(c_Key1(-10), common_term_47);
    // q = -9
    result_by_lpq.set(c_Key3(4, 1, -9), common_term_48);
    result_by_q.set(c_Key1(-9), common_term_48);
    // q = -8
    result_by_lpq.set(c_Key3(4, 1, -8), common_term_49);
    result_by_q.set(c_Key1(-8), common_term_49);
    // q = -7
    result_by_lpq.set(c_Key3(4, 1, -7), common_term_50);
    result_by_q.set(c_Key1(-7), common_term_50);
    // q = -6
    result_by_lpq.set(c_Key3(4, 1, -6), common_term_51);
    result_by_q.set(c_Key1(-6), common_term_51);
    // q = -5
    result_by_lpq.set(c_Key3(4, 1, -5), common_term_52);
    result_by_q.set(c_Key1(-5), common_term_52);
    // q = -4
    result_by_lpq.set(c_Key3(4, 1, -4), common_term_53);
    result_by_q.set(c_Key1(-4), common_term_53);
    // q = -3
    result_by_lpq.set(c_Key3(4, 1, -3), common_term_54);
    result_by_q.set(c_Key1(-3), common_term_54);
    // q = -2
    result_by_lpq.set(c_Key3(4, 1, -2), common_term_55);
    result_by_q.set(c_Key1(-2), common_term_55);
    // q = -1
    result_by_lpq.set(c_Key3(4, 1, -1), common_term_56);
    result_by_q.set(c_Key1(-1), common_term_56);
    // q = 0
    result_by_lpq.set(c_Key3(4, 1, 0), common_term_57);
    result_by_q.set(c_Key1(0), common_term_57);
    // q = 1
    result_by_lpq.set(c_Key3(4, 1, 1), common_term_58);
    result_by_q.set(c_Key1(1), common_term_58);
    // q = 2
    result_by_lpq.set(c_Key3(4, 1, 2), common_term_59);
    result_by_q.set(c_Key1(2), common_term_59);
    // q = 3
    result_by_lpq.set(c_Key3(4, 1, 3), common_term_60);
    result_by_q.set(c_Key1(3), common_term_60);
    // q = 4
    result_by_lpq.set(c_Key3(4, 1, 4), common_term_61);
    result_by_q.set(c_Key1(4), common_term_61);
    // q = 5
    result_by_lpq.set(c_Key3(4, 1, 5), common_term_62);
    result_by_q.set(c_Key1(5), common_term_62);
    // q = 6
    result_by_lpq.set(c_Key3(4, 1, 6), common_term_63);
    result_by_q.set(c_Key1(6), common_term_63);
    // q = 7
    result_by_lpq.set(c_Key3(4, 1, 7), common_term_64);
    result_by_q.set(c_Key1(7), common_term_64);
    // q = 8
    result_by_lpq.set(c_Key3(4, 1, 8), common_term_65);
    result_by_q.set(c_Key1(8), common_term_65);
    // q = 9
    result_by_lpq.set(c_Key3(4, 1, 9), common_term_66);
    result_by_q.set(c_Key1(9), common_term_66);
    // q = 10
    result_by_lpq.set(c_Key3(4, 1, 10), common_term_67);
    result_by_q.set(c_Key1(10), common_term_67);
    // q = 11
    result_by_lpq.set(c_Key3(4, 1, 11), common_term_68);
    result_by_q.set(c_Key1(11), common_term_68);
    // q = 12
    result_by_lpq.set(c_Key3(4, 1, 12), common_term_69);
    result_by_q.set(c_Key1(12), common_term_69);
    // q = 13
    result_by_lpq.set(c_Key3(4, 1, 13), common_term_70);
    result_by_q.set(c_Key1(13), common_term_70);
    // q = 14
    result_by_lpq.set(c_Key3(4, 1, 14), common_term_71);
    result_by_q.set(c_Key1(14), common_term_71);
    // q = 15
    result_by_lpq.set(c_Key3(4, 1, 15), common_term_72);
    result_by_q.set(c_Key1(15), common_term_72);
    // q = 16
    result_by_lpq.set(c_Key3(4, 1, 16), common_term_73);
    result_by_q.set(c_Key1(16), common_term_73);
    // q = 17
    result_by_lpq.set(c_Key3(4, 1, 17), common_term_74);
    result_by_q.set(c_Key1(17), common_term_74);
    // q = 18
    result_by_lpq.set(c_Key3(4, 1, 18), common_term_75);
    result_by_q.set(c_Key1(18), common_term_75);
    // q = 19
    result_by_lpq.set(c_Key3(4, 1, 19), common_term_76);
    result_by_q.set(c_Key1(19), common_term_76);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 1), result_by_q);
    result_by_q.clear();

    // l , p = (4, 2).
    // q = -19
    result_by_lpq.set(c_Key3(4, 2, -19), common_term_77);
    result_by_q.set(c_Key1(-19), common_term_77);
    // q = -18
    result_by_lpq.set(c_Key3(4, 2, -18), common_term_78);
    result_by_q.set(c_Key1(-18), common_term_78);
    // q = -17
    result_by_lpq.set(c_Key3(4, 2, -17), common_term_79);
    result_by_q.set(c_Key1(-17), common_term_79);
    // q = -16
    result_by_lpq.set(c_Key3(4, 2, -16), common_term_80);
    result_by_q.set(c_Key1(-16), common_term_80);
    // q = -15
    result_by_lpq.set(c_Key3(4, 2, -15), common_term_81);
    result_by_q.set(c_Key1(-15), common_term_81);
    // q = -14
    result_by_lpq.set(c_Key3(4, 2, -14), common_term_82);
    result_by_q.set(c_Key1(-14), common_term_82);
    // q = -13
    result_by_lpq.set(c_Key3(4, 2, -13), common_term_83);
    result_by_q.set(c_Key1(-13), common_term_83);
    // q = -12
    result_by_lpq.set(c_Key3(4, 2, -12), common_term_84);
    result_by_q.set(c_Key1(-12), common_term_84);
    // q = -11
    result_by_lpq.set(c_Key3(4, 2, -11), common_term_85);
    result_by_q.set(c_Key1(-11), common_term_85);
    // q = -10
    result_by_lpq.set(c_Key3(4, 2, -10), common_term_86);
    result_by_q.set(c_Key1(-10), common_term_86);
    // q = -9
    result_by_lpq.set(c_Key3(4, 2, -9), common_term_87);
    result_by_q.set(c_Key1(-9), common_term_87);
    // q = -8
    result_by_lpq.set(c_Key3(4, 2, -8), common_term_88);
    result_by_q.set(c_Key1(-8), common_term_88);
    // q = -7
    result_by_lpq.set(c_Key3(4, 2, -7), common_term_89);
    result_by_q.set(c_Key1(-7), common_term_89);
    // q = -6
    result_by_lpq.set(c_Key3(4, 2, -6), common_term_90);
    result_by_q.set(c_Key1(-6), common_term_90);
    // q = -5
    result_by_lpq.set(c_Key3(4, 2, -5), common_term_91);
    result_by_q.set(c_Key1(-5), common_term_91);
    // q = -4
    result_by_lpq.set(c_Key3(4, 2, -4), common_term_92);
    result_by_q.set(c_Key1(-4), common_term_92);
    // q = -3
    result_by_lpq.set(c_Key3(4, 2, -3), common_term_93);
    result_by_q.set(c_Key1(-3), common_term_93);
    // q = -2
    result_by_lpq.set(c_Key3(4, 2, -2), common_term_94);
    result_by_q.set(c_Key1(-2), common_term_94);
    // q = -1
    result_by_lpq.set(c_Key3(4, 2, -1), common_term_95);
    result_by_q.set(c_Key1(-1), common_term_95);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -3.5)*(1.5*eccentricity_2 + 1.0);
    result_by_lpq.set(c_Key3(4, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(4, 2, 1), common_term_95);
    result_by_q.set(c_Key1(1), common_term_95);
    // q = 2
    result_by_lpq.set(c_Key3(4, 2, 2), common_term_94);
    result_by_q.set(c_Key1(2), common_term_94);
    // q = 3
    result_by_lpq.set(c_Key3(4, 2, 3), common_term_93);
    result_by_q.set(c_Key1(3), common_term_93);
    // q = 4
    result_by_lpq.set(c_Key3(4, 2, 4), common_term_92);
    result_by_q.set(c_Key1(4), common_term_92);
    // q = 5
    result_by_lpq.set(c_Key3(4, 2, 5), common_term_91);
    result_by_q.set(c_Key1(5), common_term_91);
    // q = 6
    result_by_lpq.set(c_Key3(4, 2, 6), common_term_90);
    result_by_q.set(c_Key1(6), common_term_90);
    // q = 7
    result_by_lpq.set(c_Key3(4, 2, 7), common_term_89);
    result_by_q.set(c_Key1(7), common_term_89);
    // q = 8
    result_by_lpq.set(c_Key3(4, 2, 8), common_term_88);
    result_by_q.set(c_Key1(8), common_term_88);
    // q = 9
    result_by_lpq.set(c_Key3(4, 2, 9), common_term_87);
    result_by_q.set(c_Key1(9), common_term_87);
    // q = 10
    result_by_lpq.set(c_Key3(4, 2, 10), common_term_86);
    result_by_q.set(c_Key1(10), common_term_86);
    // q = 11
    result_by_lpq.set(c_Key3(4, 2, 11), common_term_85);
    result_by_q.set(c_Key1(11), common_term_85);
    // q = 12
    result_by_lpq.set(c_Key3(4, 2, 12), common_term_84);
    result_by_q.set(c_Key1(12), common_term_84);
    // q = 13
    result_by_lpq.set(c_Key3(4, 2, 13), common_term_83);
    result_by_q.set(c_Key1(13), common_term_83);
    // q = 14
    result_by_lpq.set(c_Key3(4, 2, 14), common_term_82);
    result_by_q.set(c_Key1(14), common_term_82);
    // q = 15
    result_by_lpq.set(c_Key3(4, 2, 15), common_term_81);
    result_by_q.set(c_Key1(15), common_term_81);
    // q = 16
    result_by_lpq.set(c_Key3(4, 2, 16), common_term_80);
    result_by_q.set(c_Key1(16), common_term_80);
    // q = 17
    result_by_lpq.set(c_Key3(4, 2, 17), common_term_79);
    result_by_q.set(c_Key1(17), common_term_79);
    // q = 18
    result_by_lpq.set(c_Key3(4, 2, 18), common_term_78);
    result_by_q.set(c_Key1(18), common_term_78);
    // q = 19
    result_by_lpq.set(c_Key3(4, 2, 19), common_term_77);
    result_by_q.set(c_Key1(19), common_term_77);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 2), result_by_q);
    result_by_q.clear();

    // l , p = (4, 3).
    // q = -19
    result_by_lpq.set(c_Key3(4, 3, -19), common_term_76);
    result_by_q.set(c_Key1(-19), common_term_76);
    // q = -18
    result_by_lpq.set(c_Key3(4, 3, -18), common_term_75);
    result_by_q.set(c_Key1(-18), common_term_75);
    // q = -17
    result_by_lpq.set(c_Key3(4, 3, -17), common_term_74);
    result_by_q.set(c_Key1(-17), common_term_74);
    // q = -16
    result_by_lpq.set(c_Key3(4, 3, -16), common_term_73);
    result_by_q.set(c_Key1(-16), common_term_73);
    // q = -15
    result_by_lpq.set(c_Key3(4, 3, -15), common_term_72);
    result_by_q.set(c_Key1(-15), common_term_72);
    // q = -14
    result_by_lpq.set(c_Key3(4, 3, -14), common_term_71);
    result_by_q.set(c_Key1(-14), common_term_71);
    // q = -13
    result_by_lpq.set(c_Key3(4, 3, -13), common_term_70);
    result_by_q.set(c_Key1(-13), common_term_70);
    // q = -12
    result_by_lpq.set(c_Key3(4, 3, -12), common_term_69);
    result_by_q.set(c_Key1(-12), common_term_69);
    // q = -11
    result_by_lpq.set(c_Key3(4, 3, -11), common_term_68);
    result_by_q.set(c_Key1(-11), common_term_68);
    // q = -10
    result_by_lpq.set(c_Key3(4, 3, -10), common_term_67);
    result_by_q.set(c_Key1(-10), common_term_67);
    // q = -9
    result_by_lpq.set(c_Key3(4, 3, -9), common_term_66);
    result_by_q.set(c_Key1(-9), common_term_66);
    // q = -8
    result_by_lpq.set(c_Key3(4, 3, -8), common_term_65);
    result_by_q.set(c_Key1(-8), common_term_65);
    // q = -7
    result_by_lpq.set(c_Key3(4, 3, -7), common_term_64);
    result_by_q.set(c_Key1(-7), common_term_64);
    // q = -6
    result_by_lpq.set(c_Key3(4, 3, -6), common_term_63);
    result_by_q.set(c_Key1(-6), common_term_63);
    // q = -5
    result_by_lpq.set(c_Key3(4, 3, -5), common_term_62);
    result_by_q.set(c_Key1(-5), common_term_62);
    // q = -4
    result_by_lpq.set(c_Key3(4, 3, -4), common_term_61);
    result_by_q.set(c_Key1(-4), common_term_61);
    // q = -3
    result_by_lpq.set(c_Key3(4, 3, -3), common_term_60);
    result_by_q.set(c_Key1(-3), common_term_60);
    // q = -2
    result_by_lpq.set(c_Key3(4, 3, -2), common_term_59);
    result_by_q.set(c_Key1(-2), common_term_59);
    // q = -1
    result_by_lpq.set(c_Key3(4, 3, -1), common_term_58);
    result_by_q.set(c_Key1(-1), common_term_58);
    // q = 0
    result_by_lpq.set(c_Key3(4, 3, 0), common_term_57);
    result_by_q.set(c_Key1(0), common_term_57);
    // q = 1
    result_by_lpq.set(c_Key3(4, 3, 1), common_term_56);
    result_by_q.set(c_Key1(1), common_term_56);
    // q = 2
    result_by_lpq.set(c_Key3(4, 3, 2), common_term_55);
    result_by_q.set(c_Key1(2), common_term_55);
    // q = 3
    result_by_lpq.set(c_Key3(4, 3, 3), common_term_54);
    result_by_q.set(c_Key1(3), common_term_54);
    // q = 4
    result_by_lpq.set(c_Key3(4, 3, 4), common_term_53);
    result_by_q.set(c_Key1(4), common_term_53);
    // q = 5
    result_by_lpq.set(c_Key3(4, 3, 5), common_term_52);
    result_by_q.set(c_Key1(5), common_term_52);
    // q = 6
    result_by_lpq.set(c_Key3(4, 3, 6), common_term_51);
    result_by_q.set(c_Key1(6), common_term_51);
    // q = 7
    result_by_lpq.set(c_Key3(4, 3, 7), common_term_50);
    result_by_q.set(c_Key1(7), common_term_50);
    // q = 8
    result_by_lpq.set(c_Key3(4, 3, 8), common_term_49);
    result_by_q.set(c_Key1(8), common_term_49);
    // q = 9
    result_by_lpq.set(c_Key3(4, 3, 9), common_term_48);
    result_by_q.set(c_Key1(9), common_term_48);
    // q = 10
    result_by_lpq.set(c_Key3(4, 3, 10), common_term_47);
    result_by_q.set(c_Key1(10), common_term_47);
    // q = 11
    result_by_lpq.set(c_Key3(4, 3, 11), common_term_46);
    result_by_q.set(c_Key1(11), common_term_46);
    // q = 12
    result_by_lpq.set(c_Key3(4, 3, 12), common_term_45);
    result_by_q.set(c_Key1(12), common_term_45);
    // q = 13
    result_by_lpq.set(c_Key3(4, 3, 13), common_term_44);
    result_by_q.set(c_Key1(13), common_term_44);
    // q = 14
    result_by_lpq.set(c_Key3(4, 3, 14), common_term_43);
    result_by_q.set(c_Key1(14), common_term_43);
    // q = 15
    result_by_lpq.set(c_Key3(4, 3, 15), common_term_42);
    result_by_q.set(c_Key1(15), common_term_42);
    // q = 16
    result_by_lpq.set(c_Key3(4, 3, 16), common_term_41);
    result_by_q.set(c_Key1(16), common_term_41);
    // q = 17
    result_by_lpq.set(c_Key3(4, 3, 17), common_term_40);
    result_by_q.set(c_Key1(17), common_term_40);
    // q = 18
    result_by_lpq.set(c_Key3(4, 3, 18), common_term_39);
    result_by_q.set(c_Key1(18), common_term_39);
    // q = 19
    result_by_lpq.set(c_Key3(4, 3, 19), common_term_38);
    result_by_q.set(c_Key1(19), common_term_38);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 3), result_by_q);
    result_by_q.clear();

    // l , p = (4, 4).
    // q = -19
    result_by_lpq.set(c_Key3(4, 4, -19), common_term_37);
    result_by_q.set(c_Key1(-19), common_term_37);
    // q = -18
    result_by_lpq.set(c_Key3(4, 4, -18), common_term_36);
    result_by_q.set(c_Key1(-18), common_term_36);
    // q = -17
    result_by_lpq.set(c_Key3(4, 4, -17), common_term_35);
    result_by_q.set(c_Key1(-17), common_term_35);
    // q = -16
    result_by_lpq.set(c_Key3(4, 4, -16), common_term_34);
    result_by_q.set(c_Key1(-16), common_term_34);
    // q = -15
    result_by_lpq.set(c_Key3(4, 4, -15), common_term_33);
    result_by_q.set(c_Key1(-15), common_term_33);
    // q = -14
    result_by_lpq.set(c_Key3(4, 4, -14), common_term_32);
    result_by_q.set(c_Key1(-14), common_term_32);
    // q = -13
    result_by_lpq.set(c_Key3(4, 4, -13), common_term_31);
    result_by_q.set(c_Key1(-13), common_term_31);
    // q = -12
    result_by_lpq.set(c_Key3(4, 4, -12), common_term_30);
    result_by_q.set(c_Key1(-12), common_term_30);
    // q = -11
    result_by_lpq.set(c_Key3(4, 4, -11), common_term_29);
    result_by_q.set(c_Key1(-11), common_term_29);
    // q = -10
    result_by_lpq.set(c_Key3(4, 4, -10), common_term_28);
    result_by_q.set(c_Key1(-10), common_term_28);
    // q = -9
    result_by_lpq.set(c_Key3(4, 4, -9), common_term_27);
    result_by_q.set(c_Key1(-9), common_term_27);
    // q = -8
    result_by_lpq.set(c_Key3(4, 4, -8), common_term_26);
    result_by_q.set(c_Key1(-8), common_term_26);
    // q = -7
    result_by_lpq.set(c_Key3(4, 4, -7), common_term_25);
    result_by_q.set(c_Key1(-7), common_term_25);
    // q = -6
    result_by_lpq.set(c_Key3(4, 4, -6), common_term_24);
    result_by_q.set(c_Key1(-6), common_term_24);
    // q = -5
    result_by_lpq.set(c_Key3(4, 4, -5), common_term_23);
    result_by_q.set(c_Key1(-5), common_term_23);
    // q = -4
    result_by_lpq.set(c_Key3(4, 4, -4), common_term_22);
    result_by_q.set(c_Key1(-4), common_term_22);
    // q = -3
    result_by_lpq.set(c_Key3(4, 4, -3), common_term_21);
    result_by_q.set(c_Key1(-3), common_term_21);
    // q = -2
    result_by_lpq.set(c_Key3(4, 4, -2), common_term_20);
    result_by_q.set(c_Key1(-2), common_term_20);
    // q = -1
    result_by_lpq.set(c_Key3(4, 4, -1), common_term_19);
    result_by_q.set(c_Key1(-1), common_term_19);
    // q = 0
    result_by_lpq.set(c_Key3(4, 4, 0), common_term_18);
    result_by_q.set(c_Key1(0), common_term_18);
    // q = 1
    result_by_lpq.set(c_Key3(4, 4, 1), common_term_17);
    result_by_q.set(c_Key1(1), common_term_17);
    // q = 2
    result_by_lpq.set(c_Key3(4, 4, 2), common_term_16);
    result_by_q.set(c_Key1(2), common_term_16);
    // q = 3
    result_by_lpq.set(c_Key3(4, 4, 3), common_term_15);
    result_by_q.set(c_Key1(3), common_term_15);
    // q = 5
    result_by_lpq.set(c_Key3(4, 4, 5), common_term_14);
    result_by_q.set(c_Key1(5), common_term_14);
    // q = 6
    result_by_lpq.set(c_Key3(4, 4, 6), common_term_13);
    result_by_q.set(c_Key1(6), common_term_13);
    // q = 7
    result_by_lpq.set(c_Key3(4, 4, 7), common_term_12);
    result_by_q.set(c_Key1(7), common_term_12);
    // q = 8
    result_by_lpq.set(c_Key3(4, 4, 8), common_term_11);
    result_by_q.set(c_Key1(8), common_term_11);
    // q = 9
    result_by_lpq.set(c_Key3(4, 4, 9), common_term_10);
    result_by_q.set(c_Key1(9), common_term_10);
    // q = 10
    result_by_lpq.set(c_Key3(4, 4, 10), common_term_9);
    result_by_q.set(c_Key1(10), common_term_9);
    // q = 11
    result_by_lpq.set(c_Key3(4, 4, 11), common_term_8);
    result_by_q.set(c_Key1(11), common_term_8);
    // q = 12
    result_by_lpq.set(c_Key3(4, 4, 12), common_term_7);
    result_by_q.set(c_Key1(12), common_term_7);
    // q = 13
    result_by_lpq.set(c_Key3(4, 4, 13), common_term_6);
    result_by_q.set(c_Key1(13), common_term_6);
    // q = 14
    result_by_lpq.set(c_Key3(4, 4, 14), common_term_5);
    result_by_q.set(c_Key1(14), common_term_5);
    // q = 15
    result_by_lpq.set(c_Key3(4, 4, 15), common_term_4);
    result_by_q.set(c_Key1(15), common_term_4);
    // q = 16
    result_by_lpq.set(c_Key3(4, 4, 16), common_term_3);
    result_by_q.set(c_Key1(16), common_term_3);
    // q = 17
    result_by_lpq.set(c_Key3(4, 4, 17), common_term_2);
    result_by_q.set(c_Key1(17), common_term_2);
    // q = 18
    result_by_lpq.set(c_Key3(4, 4, 18), common_term_1);
    result_by_q.set(c_Key1(18), common_term_1);
    // q = 19
    result_by_lpq.set(c_Key3(4, 4, 19), common_term_0);
    result_by_q.set(c_Key1(19), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(4, 4), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}
