#pragma once

#include "eccentricity_common_.hpp"

EccentricityFuncOutput c_eccentricity_function_l5_e1(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 5.
    // Functions are taylor expanded and truncated at eccentricity^1.

    // Functions Calculated for l = 5.

    c_IntMap<c_Key3, double> result_by_lpq(16);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(6);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = -2.0*eccentricity;
    double common_term_1 = 8.0*eccentricity;
    double common_term_2 = 6.0*eccentricity;
    double common_term_3 = 2.0*eccentricity*std::pow(1.0 - eccentricity_2, -4.5);
    double common_term_4 = 4.0*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(3);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (5, 0).
    // q = -1
    result_by_lpq.set(c_Key3(5, 0, -1), common_term_0);
    result_by_q.set(c_Key1(-1), common_term_0);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(5, 0, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(5, 0, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 0), result_by_q);
    result_by_q.clear();

    // l , p = (5, 1).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(5, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(5, 1, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 1), result_by_q);
    result_by_q.clear();

    // l , p = (5, 2).
    // q = -1
    result_by_lpq.set(c_Key3(5, 2, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(5, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(5, 2, 1), common_term_4);
    result_by_q.set(c_Key1(1), common_term_4);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 2), result_by_q);
    result_by_q.clear();

    // l , p = (5, 3).
    // q = -1
    result_by_lpq.set(c_Key3(5, 3, -1), common_term_4);
    result_by_q.set(c_Key1(-1), common_term_4);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(5, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(5, 3, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 3), result_by_q);
    result_by_q.clear();

    // l , p = (5, 4).
    // q = -1
    result_by_lpq.set(c_Key3(5, 4, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(5, 4, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 4), result_by_q);
    result_by_q.clear();

    // l , p = (5, 5).
    // q = -1
    result_by_lpq.set(c_Key3(5, 5, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(5, 5, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(5, 5, 1), common_term_0);
    result_by_q.set(c_Key1(1), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 5), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l5_e2(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 5.
    // Functions are taylor expanded and truncated at eccentricity^2.

    // Functions Calculated for l = 5.

    c_IntMap<c_Key3, double> result_by_lpq(28);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(6);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = 1.125*eccentricity_2;
    double common_term_1 = -2.0*eccentricity;
    double common_term_2 = 1.0 - 17.5*eccentricity_2;
    double common_term_3 = 8.0*eccentricity;
    double common_term_4 = 37.375*eccentricity_2;
    double common_term_5 = 0.375*eccentricity_2;
    double common_term_6 = 1.0 - 1.5*eccentricity_2;
    double common_term_7 = 6.0*eccentricity;
    double common_term_8 = 22.125*eccentricity_2;
    double common_term_9 = 3.625*eccentricity_2;
    double common_term_10 = 2.0*eccentricity*std::pow(1.0 - eccentricity_2, -4.5);
    double common_term_11 = 6.5*eccentricity_2 + 1.0;
    double common_term_12 = 4.0*eccentricity;
    double common_term_13 = 10.875*eccentricity_2;

    c_IntMap<c_Key1, double> result_by_q(5);

    // eccentricity function by mode:
    // l , p = (5, 0).
    // q = -2
    result_by_lpq.set(c_Key3(5, 0, -2), common_term_0);
    result_by_q.set(c_Key1(-2), common_term_0);
    // q = -1
    result_by_lpq.set(c_Key3(5, 0, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    result_by_lpq.set(c_Key3(5, 0, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(5, 0, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // q = 2
    result_by_lpq.set(c_Key3(5, 0, 2), common_term_4);
    result_by_q.set(c_Key1(2), common_term_4);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 0), result_by_q);
    result_by_q.clear();

    // l , p = (5, 1).
    // q = -2
    result_by_lpq.set(c_Key3(5, 1, -2), common_term_5);
    result_by_q.set(c_Key1(-2), common_term_5);
    // q = 0
    result_by_lpq.set(c_Key3(5, 1, 0), common_term_6);
    result_by_q.set(c_Key1(0), common_term_6);
    // q = 1
    result_by_lpq.set(c_Key3(5, 1, 1), common_term_7);
    result_by_q.set(c_Key1(1), common_term_7);
    // q = 2
    result_by_lpq.set(c_Key3(5, 1, 2), common_term_8);
    result_by_q.set(c_Key1(2), common_term_8);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 1), result_by_q);
    result_by_q.clear();

    // l , p = (5, 2).
    // q = -2
    result_by_lpq.set(c_Key3(5, 2, -2), common_term_9);
    result_by_q.set(c_Key1(-2), common_term_9);
    // q = -1
    result_by_lpq.set(c_Key3(5, 2, -1), common_term_10);
    result_by_q.set(c_Key1(-1), common_term_10);
    // q = 0
    result_by_lpq.set(c_Key3(5, 2, 0), common_term_11);
    result_by_q.set(c_Key1(0), common_term_11);
    // q = 1
    result_by_lpq.set(c_Key3(5, 2, 1), common_term_12);
    result_by_q.set(c_Key1(1), common_term_12);
    // q = 2
    result_by_lpq.set(c_Key3(5, 2, 2), common_term_13);
    result_by_q.set(c_Key1(2), common_term_13);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 2), result_by_q);
    result_by_q.clear();

    // l , p = (5, 3).
    // q = -2
    result_by_lpq.set(c_Key3(5, 3, -2), common_term_13);
    result_by_q.set(c_Key1(-2), common_term_13);
    // q = -1
    result_by_lpq.set(c_Key3(5, 3, -1), common_term_12);
    result_by_q.set(c_Key1(-1), common_term_12);
    // q = 0
    result_by_lpq.set(c_Key3(5, 3, 0), common_term_11);
    result_by_q.set(c_Key1(0), common_term_11);
    // q = 1
    result_by_lpq.set(c_Key3(5, 3, 1), common_term_10);
    result_by_q.set(c_Key1(1), common_term_10);
    // q = 2
    result_by_lpq.set(c_Key3(5, 3, 2), common_term_9);
    result_by_q.set(c_Key1(2), common_term_9);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 3), result_by_q);
    result_by_q.clear();

    // l , p = (5, 4).
    // q = -2
    result_by_lpq.set(c_Key3(5, 4, -2), common_term_8);
    result_by_q.set(c_Key1(-2), common_term_8);
    // q = -1
    result_by_lpq.set(c_Key3(5, 4, -1), common_term_7);
    result_by_q.set(c_Key1(-1), common_term_7);
    // q = 0
    result_by_lpq.set(c_Key3(5, 4, 0), common_term_6);
    result_by_q.set(c_Key1(0), common_term_6);
    // q = 2
    result_by_lpq.set(c_Key3(5, 4, 2), common_term_5);
    result_by_q.set(c_Key1(2), common_term_5);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 4), result_by_q);
    result_by_q.clear();

    // l , p = (5, 5).
    // q = -2
    result_by_lpq.set(c_Key3(5, 5, -2), common_term_4);
    result_by_q.set(c_Key1(-2), common_term_4);
    // q = -1
    result_by_lpq.set(c_Key3(5, 5, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    result_by_lpq.set(c_Key3(5, 5, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(5, 5, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // q = 2
    result_by_lpq.set(c_Key3(5, 5, 2), common_term_0);
    result_by_q.set(c_Key1(2), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 5), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l5_e3(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 5.
    // Functions are taylor expanded and truncated at eccentricity^3.

    // Functions Calculated for l = 5.

    c_IntMap<c_Key3, double> result_by_lpq(42);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(6);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double common_term_0 = -0.16666666666666667*eccentricity_3;
    double common_term_1 = 1.125*eccentricity_2;
    double common_term_2 = 11.5*eccentricity_3 - 2.0*eccentricity;
    double common_term_3 = 1.0 - 17.5*eccentricity_2;
    double common_term_4 = -88.5*eccentricity_3 + 8.0*eccentricity;
    double common_term_5 = 37.375*eccentricity_2;
    double common_term_6 = 133.16666666666667*eccentricity_3;
    double common_term_7 = 0.5*eccentricity_3*std::pow(1.0 - eccentricity_2, -4.5);
    double common_term_8 = 0.375*eccentricity_2;
    double common_term_9 = 1.5*eccentricity_3;
    double common_term_10 = 1.0 - 1.5*eccentricity_2;
    double common_term_11 = -10.5*eccentricity_3 + 6.0*eccentricity;
    double common_term_12 = 22.125*eccentricity_2;
    double common_term_13 = 64.5*eccentricity_3;
    double common_term_14 = 6.1666666666666667*eccentricity_3;
    double common_term_15 = 3.625*eccentricity_2;
    double common_term_16 = std::pow(1.0 - eccentricity_2, -4.5)*(1.5*eccentricity_3 + 2.0*eccentricity);
    double common_term_17 = 6.5*eccentricity_2 + 1.0;
    double common_term_18 = 14.5*eccentricity_3 + 4.0*eccentricity;
    double common_term_19 = 10.875*eccentricity_2;
    double common_term_20 = 24.833333333333333*eccentricity_3;

    c_IntMap<c_Key1, double> result_by_q(7);

    // eccentricity function by mode:
    // l , p = (5, 0).
    // q = -3
    result_by_lpq.set(c_Key3(5, 0, -3), common_term_0);
    result_by_q.set(c_Key1(-3), common_term_0);
    // q = -2
    result_by_lpq.set(c_Key3(5, 0, -2), common_term_1);
    result_by_q.set(c_Key1(-2), common_term_1);
    // q = -1
    result_by_lpq.set(c_Key3(5, 0, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    result_by_lpq.set(c_Key3(5, 0, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(5, 0, 1), common_term_4);
    result_by_q.set(c_Key1(1), common_term_4);
    // q = 2
    result_by_lpq.set(c_Key3(5, 0, 2), common_term_5);
    result_by_q.set(c_Key1(2), common_term_5);
    // q = 3
    result_by_lpq.set(c_Key3(5, 0, 3), common_term_6);
    result_by_q.set(c_Key1(3), common_term_6);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 0), result_by_q);
    result_by_q.clear();

    // l , p = (5, 1).
    // q = -3
    result_by_lpq.set(c_Key3(5, 1, -3), common_term_7);
    result_by_q.set(c_Key1(-3), common_term_7);
    // q = -2
    result_by_lpq.set(c_Key3(5, 1, -2), common_term_8);
    result_by_q.set(c_Key1(-2), common_term_8);
    // q = -1
    result_by_lpq.set(c_Key3(5, 1, -1), common_term_9);
    result_by_q.set(c_Key1(-1), common_term_9);
    // q = 0
    result_by_lpq.set(c_Key3(5, 1, 0), common_term_10);
    result_by_q.set(c_Key1(0), common_term_10);
    // q = 1
    result_by_lpq.set(c_Key3(5, 1, 1), common_term_11);
    result_by_q.set(c_Key1(1), common_term_11);
    // q = 2
    result_by_lpq.set(c_Key3(5, 1, 2), common_term_12);
    result_by_q.set(c_Key1(2), common_term_12);
    // q = 3
    result_by_lpq.set(c_Key3(5, 1, 3), common_term_13);
    result_by_q.set(c_Key1(3), common_term_13);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 1), result_by_q);
    result_by_q.clear();

    // l , p = (5, 2).
    // q = -3
    result_by_lpq.set(c_Key3(5, 2, -3), common_term_14);
    result_by_q.set(c_Key1(-3), common_term_14);
    // q = -2
    result_by_lpq.set(c_Key3(5, 2, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(5, 2, -1), common_term_16);
    result_by_q.set(c_Key1(-1), common_term_16);
    // q = 0
    result_by_lpq.set(c_Key3(5, 2, 0), common_term_17);
    result_by_q.set(c_Key1(0), common_term_17);
    // q = 1
    result_by_lpq.set(c_Key3(5, 2, 1), common_term_18);
    result_by_q.set(c_Key1(1), common_term_18);
    // q = 2
    result_by_lpq.set(c_Key3(5, 2, 2), common_term_19);
    result_by_q.set(c_Key1(2), common_term_19);
    // q = 3
    result_by_lpq.set(c_Key3(5, 2, 3), common_term_20);
    result_by_q.set(c_Key1(3), common_term_20);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 2), result_by_q);
    result_by_q.clear();

    // l , p = (5, 3).
    // q = -3
    result_by_lpq.set(c_Key3(5, 3, -3), common_term_20);
    result_by_q.set(c_Key1(-3), common_term_20);
    // q = -2
    result_by_lpq.set(c_Key3(5, 3, -2), common_term_19);
    result_by_q.set(c_Key1(-2), common_term_19);
    // q = -1
    result_by_lpq.set(c_Key3(5, 3, -1), common_term_18);
    result_by_q.set(c_Key1(-1), common_term_18);
    // q = 0
    result_by_lpq.set(c_Key3(5, 3, 0), common_term_17);
    result_by_q.set(c_Key1(0), common_term_17);
    // q = 1
    result_by_lpq.set(c_Key3(5, 3, 1), common_term_16);
    result_by_q.set(c_Key1(1), common_term_16);
    // q = 2
    result_by_lpq.set(c_Key3(5, 3, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // q = 3
    result_by_lpq.set(c_Key3(5, 3, 3), common_term_14);
    result_by_q.set(c_Key1(3), common_term_14);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 3), result_by_q);
    result_by_q.clear();

    // l , p = (5, 4).
    // q = -3
    result_by_lpq.set(c_Key3(5, 4, -3), common_term_13);
    result_by_q.set(c_Key1(-3), common_term_13);
    // q = -2
    result_by_lpq.set(c_Key3(5, 4, -2), common_term_12);
    result_by_q.set(c_Key1(-2), common_term_12);
    // q = -1
    result_by_lpq.set(c_Key3(5, 4, -1), common_term_11);
    result_by_q.set(c_Key1(-1), common_term_11);
    // q = 0
    result_by_lpq.set(c_Key3(5, 4, 0), common_term_10);
    result_by_q.set(c_Key1(0), common_term_10);
    // q = 1
    result_by_lpq.set(c_Key3(5, 4, 1), common_term_9);
    result_by_q.set(c_Key1(1), common_term_9);
    // q = 2
    result_by_lpq.set(c_Key3(5, 4, 2), common_term_8);
    result_by_q.set(c_Key1(2), common_term_8);
    // q = 3
    result_by_lpq.set(c_Key3(5, 4, 3), common_term_7);
    result_by_q.set(c_Key1(3), common_term_7);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 4), result_by_q);
    result_by_q.clear();

    // l , p = (5, 5).
    // q = -3
    result_by_lpq.set(c_Key3(5, 5, -3), common_term_6);
    result_by_q.set(c_Key1(-3), common_term_6);
    // q = -2
    result_by_lpq.set(c_Key3(5, 5, -2), common_term_5);
    result_by_q.set(c_Key1(-2), common_term_5);
    // q = -1
    result_by_lpq.set(c_Key3(5, 5, -1), common_term_4);
    result_by_q.set(c_Key1(-1), common_term_4);
    // q = 0
    result_by_lpq.set(c_Key3(5, 5, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(5, 5, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // q = 2
    result_by_lpq.set(c_Key3(5, 5, 2), common_term_1);
    result_by_q.set(c_Key1(2), common_term_1);
    // q = 3
    result_by_lpq.set(c_Key3(5, 5, 3), common_term_0);
    result_by_q.set(c_Key1(3), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 5), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l5_e4(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 5.
    // Functions are taylor expanded and truncated at eccentricity^4.

    // Functions Calculated for l = 5.

    c_IntMap<c_Key3, double> result_by_lpq(54);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(6);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double common_term_0 = 0.0026041666666666667*eccentricity_4;
    double common_term_1 = -0.16666666666666667*eccentricity_3;
    double common_term_2 = -2.25*eccentricity_4 + 1.125*eccentricity_2;
    double common_term_3 = 11.5*eccentricity_3 - 2.0*eccentricity;
    double common_term_4 = 66.484375*eccentricity_4 - 17.5*eccentricity_2 + 1.0;
    double common_term_5 = -88.5*eccentricity_3 + 8.0*eccentricity;
    double common_term_6 = -338.91666666666667*eccentricity_4 + 37.375*eccentricity_2;
    double common_term_7 = 133.16666666666667*eccentricity_3;
    double common_term_8 = 400.6796875*eccentricity_4;
    double common_term_9 = 0.6796875*eccentricity_4;
    double common_term_10 = 0.5*eccentricity_3*std::pow(1.0 - eccentricity_2, -4.5);
    double common_term_11 = 1.75*eccentricity_4 + 0.375*eccentricity_2;
    double common_term_12 = 1.5*eccentricity_3;
    double common_term_13 = 4.734375*eccentricity_4 - 1.5*eccentricity_2 + 1.0;
    double common_term_14 = -10.5*eccentricity_3 + 6.0*eccentricity;
    double common_term_15 = -44.25*eccentricity_4 + 22.125*eccentricity_2;
    double common_term_16 = 64.5*eccentricity_3;
    double common_term_17 = 163.0859375*eccentricity_4;
    double common_term_18 = 10.0859375*eccentricity_4;
    double common_term_19 = 6.1666666666666667*eccentricity_3;
    double common_term_20 = 16.083333333333333*eccentricity_4 + 3.625*eccentricity_2;
    double common_term_21 = std::pow(1.0 - eccentricity_2, -4.5)*(1.5*eccentricity_3 + 2.0*eccentricity);
    double common_term_22 = 21.859375*eccentricity_4 + 6.5*eccentricity_2 + 1.0;
    double common_term_23 = 14.5*eccentricity_3 + 4.0*eccentricity;
    double common_term_24 = 26.75*eccentricity_4 + 10.875*eccentricity_2;
    double common_term_25 = 24.833333333333333*eccentricity_3;
    double common_term_26 = 51.221354166666667*eccentricity_4;

    c_IntMap<c_Key1, double> result_by_q(9);

    // eccentricity function by mode:
    // l , p = (5, 0).
    // q = -4
    result_by_lpq.set(c_Key3(5, 0, -4), common_term_0);
    result_by_q.set(c_Key1(-4), common_term_0);
    // q = -3
    result_by_lpq.set(c_Key3(5, 0, -3), common_term_1);
    result_by_q.set(c_Key1(-3), common_term_1);
    // q = -2
    result_by_lpq.set(c_Key3(5, 0, -2), common_term_2);
    result_by_q.set(c_Key1(-2), common_term_2);
    // q = -1
    result_by_lpq.set(c_Key3(5, 0, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    result_by_lpq.set(c_Key3(5, 0, 0), common_term_4);
    result_by_q.set(c_Key1(0), common_term_4);
    // q = 1
    result_by_lpq.set(c_Key3(5, 0, 1), common_term_5);
    result_by_q.set(c_Key1(1), common_term_5);
    // q = 2
    result_by_lpq.set(c_Key3(5, 0, 2), common_term_6);
    result_by_q.set(c_Key1(2), common_term_6);
    // q = 3
    result_by_lpq.set(c_Key3(5, 0, 3), common_term_7);
    result_by_q.set(c_Key1(3), common_term_7);
    // q = 4
    result_by_lpq.set(c_Key3(5, 0, 4), common_term_8);
    result_by_q.set(c_Key1(4), common_term_8);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 0), result_by_q);
    result_by_q.clear();

    // l , p = (5, 1).
    // q = -4
    result_by_lpq.set(c_Key3(5, 1, -4), common_term_9);
    result_by_q.set(c_Key1(-4), common_term_9);
    // q = -3
    result_by_lpq.set(c_Key3(5, 1, -3), common_term_10);
    result_by_q.set(c_Key1(-3), common_term_10);
    // q = -2
    result_by_lpq.set(c_Key3(5, 1, -2), common_term_11);
    result_by_q.set(c_Key1(-2), common_term_11);
    // q = -1
    result_by_lpq.set(c_Key3(5, 1, -1), common_term_12);
    result_by_q.set(c_Key1(-1), common_term_12);
    // q = 0
    result_by_lpq.set(c_Key3(5, 1, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(5, 1, 1), common_term_14);
    result_by_q.set(c_Key1(1), common_term_14);
    // q = 2
    result_by_lpq.set(c_Key3(5, 1, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // q = 3
    result_by_lpq.set(c_Key3(5, 1, 3), common_term_16);
    result_by_q.set(c_Key1(3), common_term_16);
    // q = 4
    result_by_lpq.set(c_Key3(5, 1, 4), common_term_17);
    result_by_q.set(c_Key1(4), common_term_17);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 1), result_by_q);
    result_by_q.clear();

    // l , p = (5, 2).
    // q = -4
    result_by_lpq.set(c_Key3(5, 2, -4), common_term_18);
    result_by_q.set(c_Key1(-4), common_term_18);
    // q = -3
    result_by_lpq.set(c_Key3(5, 2, -3), common_term_19);
    result_by_q.set(c_Key1(-3), common_term_19);
    // q = -2
    result_by_lpq.set(c_Key3(5, 2, -2), common_term_20);
    result_by_q.set(c_Key1(-2), common_term_20);
    // q = -1
    result_by_lpq.set(c_Key3(5, 2, -1), common_term_21);
    result_by_q.set(c_Key1(-1), common_term_21);
    // q = 0
    result_by_lpq.set(c_Key3(5, 2, 0), common_term_22);
    result_by_q.set(c_Key1(0), common_term_22);
    // q = 1
    result_by_lpq.set(c_Key3(5, 2, 1), common_term_23);
    result_by_q.set(c_Key1(1), common_term_23);
    // q = 2
    result_by_lpq.set(c_Key3(5, 2, 2), common_term_24);
    result_by_q.set(c_Key1(2), common_term_24);
    // q = 3
    result_by_lpq.set(c_Key3(5, 2, 3), common_term_25);
    result_by_q.set(c_Key1(3), common_term_25);
    // q = 4
    result_by_lpq.set(c_Key3(5, 2, 4), common_term_26);
    result_by_q.set(c_Key1(4), common_term_26);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 2), result_by_q);
    result_by_q.clear();

    // l , p = (5, 3).
    // q = -4
    result_by_lpq.set(c_Key3(5, 3, -4), common_term_26);
    result_by_q.set(c_Key1(-4), common_term_26);
    // q = -3
    result_by_lpq.set(c_Key3(5, 3, -3), common_term_25);
    result_by_q.set(c_Key1(-3), common_term_25);
    // q = -2
    result_by_lpq.set(c_Key3(5, 3, -2), common_term_24);
    result_by_q.set(c_Key1(-2), common_term_24);
    // q = -1
    result_by_lpq.set(c_Key3(5, 3, -1), common_term_23);
    result_by_q.set(c_Key1(-1), common_term_23);
    // q = 0
    result_by_lpq.set(c_Key3(5, 3, 0), common_term_22);
    result_by_q.set(c_Key1(0), common_term_22);
    // q = 1
    result_by_lpq.set(c_Key3(5, 3, 1), common_term_21);
    result_by_q.set(c_Key1(1), common_term_21);
    // q = 2
    result_by_lpq.set(c_Key3(5, 3, 2), common_term_20);
    result_by_q.set(c_Key1(2), common_term_20);
    // q = 3
    result_by_lpq.set(c_Key3(5, 3, 3), common_term_19);
    result_by_q.set(c_Key1(3), common_term_19);
    // q = 4
    result_by_lpq.set(c_Key3(5, 3, 4), common_term_18);
    result_by_q.set(c_Key1(4), common_term_18);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 3), result_by_q);
    result_by_q.clear();

    // l , p = (5, 4).
    // q = -4
    result_by_lpq.set(c_Key3(5, 4, -4), common_term_17);
    result_by_q.set(c_Key1(-4), common_term_17);
    // q = -3
    result_by_lpq.set(c_Key3(5, 4, -3), common_term_16);
    result_by_q.set(c_Key1(-3), common_term_16);
    // q = -2
    result_by_lpq.set(c_Key3(5, 4, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(5, 4, -1), common_term_14);
    result_by_q.set(c_Key1(-1), common_term_14);
    // q = 0
    result_by_lpq.set(c_Key3(5, 4, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(5, 4, 1), common_term_12);
    result_by_q.set(c_Key1(1), common_term_12);
    // q = 2
    result_by_lpq.set(c_Key3(5, 4, 2), common_term_11);
    result_by_q.set(c_Key1(2), common_term_11);
    // q = 3
    result_by_lpq.set(c_Key3(5, 4, 3), common_term_10);
    result_by_q.set(c_Key1(3), common_term_10);
    // q = 4
    result_by_lpq.set(c_Key3(5, 4, 4), common_term_9);
    result_by_q.set(c_Key1(4), common_term_9);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 4), result_by_q);
    result_by_q.clear();

    // l , p = (5, 5).
    // q = -4
    result_by_lpq.set(c_Key3(5, 5, -4), common_term_8);
    result_by_q.set(c_Key1(-4), common_term_8);
    // q = -3
    result_by_lpq.set(c_Key3(5, 5, -3), common_term_7);
    result_by_q.set(c_Key1(-3), common_term_7);
    // q = -2
    result_by_lpq.set(c_Key3(5, 5, -2), common_term_6);
    result_by_q.set(c_Key1(-2), common_term_6);
    // q = -1
    result_by_lpq.set(c_Key3(5, 5, -1), common_term_5);
    result_by_q.set(c_Key1(-1), common_term_5);
    // q = 0
    result_by_lpq.set(c_Key3(5, 5, 0), common_term_4);
    result_by_q.set(c_Key1(0), common_term_4);
    // q = 1
    result_by_lpq.set(c_Key3(5, 5, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // q = 2
    result_by_lpq.set(c_Key3(5, 5, 2), common_term_2);
    result_by_q.set(c_Key1(2), common_term_2);
    // q = 3
    result_by_lpq.set(c_Key3(5, 5, 3), common_term_1);
    result_by_q.set(c_Key1(3), common_term_1);
    // q = 4
    result_by_lpq.set(c_Key3(5, 5, 4), common_term_0);
    result_by_q.set(c_Key1(4), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 5), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l5_e5(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 5.
    // Functions are taylor expanded and truncated at eccentricity^5.

    // Functions Calculated for l = 5.

    c_IntMap<c_Key3, double> result_by_lpq(64);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(6);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double eccentricity_5 = eccentricity_4 * eccentricity;
    double common_term_0 = 0.0026041666666666667*eccentricity_4;
    double common_term_1 = 0.041666666666666667*eccentricity_5 - 0.16666666666666667*eccentricity_3;
    double common_term_2 = -2.25*eccentricity_4 + 1.125*eccentricity_2;
    double common_term_3 = -16.833333333333333*eccentricity_5 + 11.5*eccentricity_3 - 2.0*eccentricity;
    double common_term_4 = 66.484375*eccentricity_4 - 17.5*eccentricity_2 + 1.0;
    double common_term_5 = 286.5*eccentricity_5 - 88.5*eccentricity_3 + 8.0*eccentricity;
    double common_term_6 = -338.91666666666667*eccentricity_4 + 37.375*eccentricity_2;
    double common_term_7 = -1090.1666666666667*eccentricity_5 + 133.16666666666667*eccentricity_3;
    double common_term_8 = 400.6796875*eccentricity_4;
    double common_term_9 = 1072.4583333333333*eccentricity_5;
    double common_term_10 = 0.925*eccentricity_5;
    double common_term_11 = 0.6796875*eccentricity_4;
    double common_term_12 = 0.5*eccentricity_3*std::pow(1.0 - eccentricity_2, -4.5);
    double common_term_13 = 1.75*eccentricity_4 + 0.375*eccentricity_2;
    double common_term_14 = 4.0*eccentricity_5 + 1.5*eccentricity_3;
    double common_term_15 = 4.734375*eccentricity_4 - 1.5*eccentricity_2 + 1.0;
    double common_term_16 = 15.5*eccentricity_5 - 10.5*eccentricity_3 + 6.0*eccentricity;
    double common_term_17 = -44.25*eccentricity_4 + 22.125*eccentricity_2;
    double common_term_18 = -145.125*eccentricity_5 + 64.5*eccentricity_3;
    double common_term_19 = 163.0859375*eccentricity_4;
    double common_term_20 = 374.45*eccentricity_5;
    double common_term_21 = 16.058333333333333*eccentricity_5;
    double common_term_22 = 10.0859375*eccentricity_4;
    double common_term_23 = 23.708333333333333*eccentricity_5 + 6.1666666666666667*eccentricity_3;
    double common_term_24 = 16.083333333333333*eccentricity_4 + 3.625*eccentricity_2;
    double common_term_25 = std::pow(1.0 - eccentricity_2, -4.5)*(1.5*eccentricity_3 + 2.0*eccentricity);
    double common_term_26 = 21.859375*eccentricity_4 + 6.5*eccentricity_2 + 1.0;
    double common_term_27 = 38.916666666666667*eccentricity_5 + 14.5*eccentricity_3 + 4.0*eccentricity;
    double common_term_28 = 26.75*eccentricity_4 + 10.875*eccentricity_2;
    double common_term_29 = 43.041666666666667*eccentricity_5 + 24.833333333333333*eccentricity_3;
    double common_term_30 = 51.221354166666667*eccentricity_4;
    double common_term_31 = 98.775*eccentricity_5;

    c_IntMap<c_Key1, double> result_by_q(11);

    // eccentricity function by mode:
    // l , p = (5, 0).
    // q = -4
    result_by_lpq.set(c_Key3(5, 0, -4), common_term_0);
    result_by_q.set(c_Key1(-4), common_term_0);
    // q = -3
    result_by_lpq.set(c_Key3(5, 0, -3), common_term_1);
    result_by_q.set(c_Key1(-3), common_term_1);
    // q = -2
    result_by_lpq.set(c_Key3(5, 0, -2), common_term_2);
    result_by_q.set(c_Key1(-2), common_term_2);
    // q = -1
    result_by_lpq.set(c_Key3(5, 0, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    result_by_lpq.set(c_Key3(5, 0, 0), common_term_4);
    result_by_q.set(c_Key1(0), common_term_4);
    // q = 1
    result_by_lpq.set(c_Key3(5, 0, 1), common_term_5);
    result_by_q.set(c_Key1(1), common_term_5);
    // q = 2
    result_by_lpq.set(c_Key3(5, 0, 2), common_term_6);
    result_by_q.set(c_Key1(2), common_term_6);
    // q = 3
    result_by_lpq.set(c_Key3(5, 0, 3), common_term_7);
    result_by_q.set(c_Key1(3), common_term_7);
    // q = 4
    result_by_lpq.set(c_Key3(5, 0, 4), common_term_8);
    result_by_q.set(c_Key1(4), common_term_8);
    // q = 5
    result_by_lpq.set(c_Key3(5, 0, 5), common_term_9);
    result_by_q.set(c_Key1(5), common_term_9);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 0), result_by_q);
    result_by_q.clear();

    // l , p = (5, 1).
    // q = -5
    result_by_lpq.set(c_Key3(5, 1, -5), common_term_10);
    result_by_q.set(c_Key1(-5), common_term_10);
    // q = -4
    result_by_lpq.set(c_Key3(5, 1, -4), common_term_11);
    result_by_q.set(c_Key1(-4), common_term_11);
    // q = -3
    result_by_lpq.set(c_Key3(5, 1, -3), common_term_12);
    result_by_q.set(c_Key1(-3), common_term_12);
    // q = -2
    result_by_lpq.set(c_Key3(5, 1, -2), common_term_13);
    result_by_q.set(c_Key1(-2), common_term_13);
    // q = -1
    result_by_lpq.set(c_Key3(5, 1, -1), common_term_14);
    result_by_q.set(c_Key1(-1), common_term_14);
    // q = 0
    result_by_lpq.set(c_Key3(5, 1, 0), common_term_15);
    result_by_q.set(c_Key1(0), common_term_15);
    // q = 1
    result_by_lpq.set(c_Key3(5, 1, 1), common_term_16);
    result_by_q.set(c_Key1(1), common_term_16);
    // q = 2
    result_by_lpq.set(c_Key3(5, 1, 2), common_term_17);
    result_by_q.set(c_Key1(2), common_term_17);
    // q = 3
    result_by_lpq.set(c_Key3(5, 1, 3), common_term_18);
    result_by_q.set(c_Key1(3), common_term_18);
    // q = 4
    result_by_lpq.set(c_Key3(5, 1, 4), common_term_19);
    result_by_q.set(c_Key1(4), common_term_19);
    // q = 5
    result_by_lpq.set(c_Key3(5, 1, 5), common_term_20);
    result_by_q.set(c_Key1(5), common_term_20);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 1), result_by_q);
    result_by_q.clear();

    // l , p = (5, 2).
    // q = -5
    result_by_lpq.set(c_Key3(5, 2, -5), common_term_21);
    result_by_q.set(c_Key1(-5), common_term_21);
    // q = -4
    result_by_lpq.set(c_Key3(5, 2, -4), common_term_22);
    result_by_q.set(c_Key1(-4), common_term_22);
    // q = -3
    result_by_lpq.set(c_Key3(5, 2, -3), common_term_23);
    result_by_q.set(c_Key1(-3), common_term_23);
    // q = -2
    result_by_lpq.set(c_Key3(5, 2, -2), common_term_24);
    result_by_q.set(c_Key1(-2), common_term_24);
    // q = -1
    result_by_lpq.set(c_Key3(5, 2, -1), common_term_25);
    result_by_q.set(c_Key1(-1), common_term_25);
    // q = 0
    result_by_lpq.set(c_Key3(5, 2, 0), common_term_26);
    result_by_q.set(c_Key1(0), common_term_26);
    // q = 1
    result_by_lpq.set(c_Key3(5, 2, 1), common_term_27);
    result_by_q.set(c_Key1(1), common_term_27);
    // q = 2
    result_by_lpq.set(c_Key3(5, 2, 2), common_term_28);
    result_by_q.set(c_Key1(2), common_term_28);
    // q = 3
    result_by_lpq.set(c_Key3(5, 2, 3), common_term_29);
    result_by_q.set(c_Key1(3), common_term_29);
    // q = 4
    result_by_lpq.set(c_Key3(5, 2, 4), common_term_30);
    result_by_q.set(c_Key1(4), common_term_30);
    // q = 5
    result_by_lpq.set(c_Key3(5, 2, 5), common_term_31);
    result_by_q.set(c_Key1(5), common_term_31);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 2), result_by_q);
    result_by_q.clear();

    // l , p = (5, 3).
    // q = -5
    result_by_lpq.set(c_Key3(5, 3, -5), common_term_31);
    result_by_q.set(c_Key1(-5), common_term_31);
    // q = -4
    result_by_lpq.set(c_Key3(5, 3, -4), common_term_30);
    result_by_q.set(c_Key1(-4), common_term_30);
    // q = -3
    result_by_lpq.set(c_Key3(5, 3, -3), common_term_29);
    result_by_q.set(c_Key1(-3), common_term_29);
    // q = -2
    result_by_lpq.set(c_Key3(5, 3, -2), common_term_28);
    result_by_q.set(c_Key1(-2), common_term_28);
    // q = -1
    result_by_lpq.set(c_Key3(5, 3, -1), common_term_27);
    result_by_q.set(c_Key1(-1), common_term_27);
    // q = 0
    result_by_lpq.set(c_Key3(5, 3, 0), common_term_26);
    result_by_q.set(c_Key1(0), common_term_26);
    // q = 1
    result_by_lpq.set(c_Key3(5, 3, 1), common_term_25);
    result_by_q.set(c_Key1(1), common_term_25);
    // q = 2
    result_by_lpq.set(c_Key3(5, 3, 2), common_term_24);
    result_by_q.set(c_Key1(2), common_term_24);
    // q = 3
    result_by_lpq.set(c_Key3(5, 3, 3), common_term_23);
    result_by_q.set(c_Key1(3), common_term_23);
    // q = 4
    result_by_lpq.set(c_Key3(5, 3, 4), common_term_22);
    result_by_q.set(c_Key1(4), common_term_22);
    // q = 5
    result_by_lpq.set(c_Key3(5, 3, 5), common_term_21);
    result_by_q.set(c_Key1(5), common_term_21);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 3), result_by_q);
    result_by_q.clear();

    // l , p = (5, 4).
    // q = -5
    result_by_lpq.set(c_Key3(5, 4, -5), common_term_20);
    result_by_q.set(c_Key1(-5), common_term_20);
    // q = -4
    result_by_lpq.set(c_Key3(5, 4, -4), common_term_19);
    result_by_q.set(c_Key1(-4), common_term_19);
    // q = -3
    result_by_lpq.set(c_Key3(5, 4, -3), common_term_18);
    result_by_q.set(c_Key1(-3), common_term_18);
    // q = -2
    result_by_lpq.set(c_Key3(5, 4, -2), common_term_17);
    result_by_q.set(c_Key1(-2), common_term_17);
    // q = -1
    result_by_lpq.set(c_Key3(5, 4, -1), common_term_16);
    result_by_q.set(c_Key1(-1), common_term_16);
    // q = 0
    result_by_lpq.set(c_Key3(5, 4, 0), common_term_15);
    result_by_q.set(c_Key1(0), common_term_15);
    // q = 1
    result_by_lpq.set(c_Key3(5, 4, 1), common_term_14);
    result_by_q.set(c_Key1(1), common_term_14);
    // q = 2
    result_by_lpq.set(c_Key3(5, 4, 2), common_term_13);
    result_by_q.set(c_Key1(2), common_term_13);
    // q = 3
    result_by_lpq.set(c_Key3(5, 4, 3), common_term_12);
    result_by_q.set(c_Key1(3), common_term_12);
    // q = 4
    result_by_lpq.set(c_Key3(5, 4, 4), common_term_11);
    result_by_q.set(c_Key1(4), common_term_11);
    // q = 5
    result_by_lpq.set(c_Key3(5, 4, 5), common_term_10);
    result_by_q.set(c_Key1(5), common_term_10);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 4), result_by_q);
    result_by_q.clear();

    // l , p = (5, 5).
    // q = -5
    result_by_lpq.set(c_Key3(5, 5, -5), common_term_9);
    result_by_q.set(c_Key1(-5), common_term_9);
    // q = -4
    result_by_lpq.set(c_Key3(5, 5, -4), common_term_8);
    result_by_q.set(c_Key1(-4), common_term_8);
    // q = -3
    result_by_lpq.set(c_Key3(5, 5, -3), common_term_7);
    result_by_q.set(c_Key1(-3), common_term_7);
    // q = -2
    result_by_lpq.set(c_Key3(5, 5, -2), common_term_6);
    result_by_q.set(c_Key1(-2), common_term_6);
    // q = -1
    result_by_lpq.set(c_Key3(5, 5, -1), common_term_5);
    result_by_q.set(c_Key1(-1), common_term_5);
    // q = 0
    result_by_lpq.set(c_Key3(5, 5, 0), common_term_4);
    result_by_q.set(c_Key1(0), common_term_4);
    // q = 1
    result_by_lpq.set(c_Key3(5, 5, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // q = 2
    result_by_lpq.set(c_Key3(5, 5, 2), common_term_2);
    result_by_q.set(c_Key1(2), common_term_2);
    // q = 3
    result_by_lpq.set(c_Key3(5, 5, 3), common_term_1);
    result_by_q.set(c_Key1(3), common_term_1);
    // q = 4
    result_by_lpq.set(c_Key3(5, 5, 4), common_term_0);
    result_by_q.set(c_Key1(4), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 5), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l5_e10(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 5.
    // Functions are taylor expanded and truncated at eccentricity^10.

    // Functions Calculated for l = 5.

    c_IntMap<c_Key3, double> result_by_lpq(124);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(6);
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
    double common_term_0 = 0.0026280707572923556*eccentricity_10;
    double common_term_1 = 0.0014109347442680776*eccentricity_9;
    double common_term_2 = 0.0011653355189732143*eccentricity_10 + 0.00063563755580357143*eccentricity_8;
    double common_term_3 = 0.00034722222222222222*eccentricity_9 + 0.00019841269841269841*eccentricity_7;
    double common_term_4 = 3.8885691809275794e-5*eccentricity_10 + 3.410218253968254e-5*eccentricity_8 + 2.1701388888888889e-5*eccentricity_6;
    double common_term_5 = 0.0012444067253637566*eccentricity_10 + 0.0014946831597222222*eccentricity_8 + 0.0018229166666666667*eccentricity_6 + 0.0026041666666666667*eccentricity_4;
    double common_term_6 = 0.001099537037037037*eccentricity_9 - 0.0041666666666666667*eccentricity_7 + 0.041666666666666667*eccentricity_5 - 0.16666666666666667*eccentricity_3;
    double common_term_7 = 0.00580902099609375*eccentricity_10 - 0.1740234375*eccentricity_8 + 0.8525390625*eccentricity_6 - 2.25*eccentricity_4 + 1.125*eccentricity_2;
    double common_term_8 = -2.2220486111111111*eccentricity_9 + 8.3368055555555556*eccentricity_7 - 16.833333333333333*eccentricity_5 + 11.5*eccentricity_3 - 2.0*eccentricity;
    double common_term_9 = -17.866980658637153*eccentricity_10 + 53.395758734809028*eccentricity_8 - 90.342881944444444*eccentricity_6 + 66.484375*eccentricity_4 - 17.5*eccentricity_2 + 1.0;
    double common_term_10 = 262.5140625*eccentricity_9 - 388.625*eccentricity_7 + 286.5*eccentricity_5 - 88.5*eccentricity_3 + 8.0*eccentricity;
    double common_term_11 = 1075.0551734641746*eccentricity_10 - 1426.3118272569444*eccentricity_8 + 1024.1761067708333*eccentricity_6 - 338.91666666666667*eccentricity_4 + 37.375*eccentricity_2;
    double common_term_12 = -4644.6884259259259*eccentricity_9 + 3210.0166666666667*eccentricity_7 - 1090.1666666666667*eccentricity_5 + 133.16666666666667*eccentricity_3;
    double common_term_13 = -13776.578367396763*eccentricity_10 + 9117.8461669921875*eccentricity_8 - 3105.83671875*eccentricity_6 + 400.6796875*eccentricity_4;
    double common_term_14 = 23983.827504960317*eccentricity_9 - 8086.7361111111111*eccentricity_7 + 1072.4583333333333*eccentricity_5;
    double common_term_15 = 59306.949686698308*eccentricity_10 - 19637.63640562996*eccentricity_8 + 2633.1491102430556*eccentricity_6;
    double common_term_16 = -45101.5265625*eccentricity_9 + 6049.0116071428571*eccentricity_7;
    double common_term_17 = -98957.484222180671*eccentricity_10 + 13179.87816617451*eccentricity_8;
    double common_term_18 = 27504.092611882716*eccentricity_9;
    double common_term_19 = 55372.913824353899*eccentricity_10;
    double common_term_20 = 4.3938146464029948*eccentricity_10;
    double common_term_21 = 3.2127232142857143*eccentricity_9;
    double common_term_22 = 7.5595372921575314*eccentricity_10 + 2.3502705891927083*eccentricity_8;
    double common_term_23 = 5.9837797619047619*eccentricity_9 + 1.7205357142857143*eccentricity_7;
    double common_term_24 = 11.538018035888672*eccentricity_10 + 4.7125279017857143*eccentricity_8 + 1.2607421875*eccentricity_6;
    double common_term_25 = 9.410639880952381*eccentricity_9 + 3.6958333333333333*eccentricity_7 + 0.925*eccentricity_5;
    double common_term_26 = 16.115094284784226*eccentricity_10 + 7.6452555338541667*eccentricity_8 + 2.88828125*eccentricity_6 + 0.6796875*eccentricity_4;
    double common_term_27 = 0.5*eccentricity_3*std::pow(1.0 - eccentricity_2, -4.5);
    double common_term_28 = 21.283941819932726*eccentricity_10 + 11.114388020833333*eccentricity_8 + 4.9912109375*eccentricity_6 + 1.75*eccentricity_4 + 0.375*eccentricity_2;
    double common_term_29 = 17.966145833333333*eccentricity_9 + 9.1875*eccentricity_7 + 4.0*eccentricity_5 + 1.5*eccentricity_3;
    double common_term_30 = 27.040718994140625*eccentricity_10 + 15.16680908203125*eccentricity_8 + 7.2265625*eccentricity_6 + 4.734375*eccentricity_4 - 1.5*eccentricity_2 + 1.0;
    double common_term_31 = 23.668229166666667*eccentricity_9 + 9.6770833333333333*eccentricity_7 + 15.5*eccentricity_5 - 10.5*eccentricity_3 + 6.0*eccentricity;
    double common_term_32 = 37.758710225423177*eccentricity_10 + 2.5192057291666667*eccentricity_8 + 51.7041015625*eccentricity_6 - 44.25*eccentricity_4 + 22.125*eccentricity_2;
    double common_term_33 = -45.9796875*eccentricity_9 + 164.025*eccentricity_7 - 145.125*eccentricity_5 + 64.5*eccentricity_3;
    double common_term_34 = -231.58345404730903*eccentricity_10 + 480.89112141927083*eccentricity_8 - 407.71484375*eccentricity_6 + 163.0859375*eccentricity_4;
    double common_term_35 = 1302.5502976190476*eccentricity_9 - 1029.7375*eccentricity_7 + 374.45*eccentricity_5;
    double common_term_36 = 3289.3376407078334*eccentricity_10 - 2404.9318359375*eccentricity_8 + 801.6439453125*eccentricity_6;
    double common_term_37 = -5288.5431547619048*eccentricity_9 + 1627.2440476190476*eccentricity_7;
    double common_term_38 = -11086.362901814779*eccentricity_10 + 3167.5322576613653*eccentricity_8;
    double common_term_39 = 5960.3247767857143*eccentricity_9;
    double common_term_40 = 10905.918879855524*eccentricity_10;
    double common_term_41 = 131.61114809308733*eccentricity_10;
    double common_term_42 = 88.163544422398589*eccentricity_9;
    double common_term_43 = 109.04968559359327*eccentricity_10 + 58.584195285373264*eccentricity_8;
    double common_term_44 = 84.790848214285714*eccentricity_9 + 38.552678571428571*eccentricity_7;
    double common_term_45 = 135.52102582416837*eccentricity_10 + 64.078047495039683*eccentricity_8 + 25.069552951388889*eccentricity_6;
    double common_term_46 = 104.94945436507937*eccentricity_9 + 47.195138888888889*eccentricity_7 + 16.058333333333333*eccentricity_5;
    double common_term_47 = 157.00033656529018*eccentricity_10 + 80.0341552734375*eccentricity_8 + 33.89765625*eccentricity_6 + 10.0859375*eccentricity_4;
    double common_term_48 = 123.10688657407407*eccentricity_9 + 59.960416666666667*eccentricity_7 + 23.708333333333333*eccentricity_5 + 6.1666666666666667*eccentricity_3;
    double common_term_49 = 178.02036545364945*eccentricity_10 + 95.092730034722222*eccentricity_8 + 44.002278645833333*eccentricity_6 + 16.083333333333333*eccentricity_4 + 3.625*eccentricity_2;
    double common_term_50 = std::pow(1.0 - eccentricity_2, -4.5)*(1.5*eccentricity_3 + 2.0*eccentricity);
    double common_term_51 = 198.49574069552951*eccentricity_10 + 109.66821967230903*eccentricity_8 + 53.677951388888889*eccentricity_6 + 21.859375*eccentricity_4 + 6.5*eccentricity_2 + 1.0;
    double common_term_52 = 157.90815972222222*eccentricity_9 + 83.961805555555556*eccentricity_7 + 38.916666666666667*eccentricity_5 + 14.5*eccentricity_3 + 4.0*eccentricity;
    double common_term_53 = 218.42696685791016*eccentricity_10 + 123.7505859375*eccentricity_8 + 63.0263671875*eccentricity_6 + 26.75*eccentricity_4 + 10.875*eccentricity_2;
    double common_term_54 = 174.4318287037037*eccentricity_9 + 95.977083333333333*eccentricity_7 + 43.041666666666667*eccentricity_5 + 24.833333333333333*eccentricity_3;
    double common_term_55 = 237.09027357959243*eccentricity_10 + 140.88095431857639*eccentricity_8 + 60.779947916666667*eccentricity_6 + 51.221354166666667*eccentricity_4;
    double common_term_56 = 204.35424107142857*eccentricity_9 + 72.31875*eccentricity_7 + 98.775*eccentricity_5;
    double common_term_57 = 301.28983722262912*eccentricity_10 + 60.194509548611111*eccentricity_8 + 181.54624565972222*eccentricity_6;
    double common_term_58 = -11.144667658730159*eccentricity_9 + 321.84662698412698*eccentricity_7;
    double common_term_59 = -209.35917576381138*eccentricity_10 + 554.7166024344308*eccentricity_8;
    double common_term_60 = 934.67407958553792*eccentricity_9;
    double common_term_61 = 1545.8537958783616*eccentricity_10;

    c_IntMap<c_Key1, double> result_by_q(21);

    // eccentricity function by mode:
    // l , p = (5, 0).
    // q = -10
    result_by_lpq.set(c_Key3(5, 0, -10), common_term_0);
    result_by_q.set(c_Key1(-10), common_term_0);
    // q = -9
    result_by_lpq.set(c_Key3(5, 0, -9), common_term_1);
    result_by_q.set(c_Key1(-9), common_term_1);
    // q = -8
    result_by_lpq.set(c_Key3(5, 0, -8), common_term_2);
    result_by_q.set(c_Key1(-8), common_term_2);
    // q = -7
    result_by_lpq.set(c_Key3(5, 0, -7), common_term_3);
    result_by_q.set(c_Key1(-7), common_term_3);
    // q = -6
    result_by_lpq.set(c_Key3(5, 0, -6), common_term_4);
    result_by_q.set(c_Key1(-6), common_term_4);
    // q = -4
    result_by_lpq.set(c_Key3(5, 0, -4), common_term_5);
    result_by_q.set(c_Key1(-4), common_term_5);
    // q = -3
    result_by_lpq.set(c_Key3(5, 0, -3), common_term_6);
    result_by_q.set(c_Key1(-3), common_term_6);
    // q = -2
    result_by_lpq.set(c_Key3(5, 0, -2), common_term_7);
    result_by_q.set(c_Key1(-2), common_term_7);
    // q = -1
    result_by_lpq.set(c_Key3(5, 0, -1), common_term_8);
    result_by_q.set(c_Key1(-1), common_term_8);
    // q = 0
    result_by_lpq.set(c_Key3(5, 0, 0), common_term_9);
    result_by_q.set(c_Key1(0), common_term_9);
    // q = 1
    result_by_lpq.set(c_Key3(5, 0, 1), common_term_10);
    result_by_q.set(c_Key1(1), common_term_10);
    // q = 2
    result_by_lpq.set(c_Key3(5, 0, 2), common_term_11);
    result_by_q.set(c_Key1(2), common_term_11);
    // q = 3
    result_by_lpq.set(c_Key3(5, 0, 3), common_term_12);
    result_by_q.set(c_Key1(3), common_term_12);
    // q = 4
    result_by_lpq.set(c_Key3(5, 0, 4), common_term_13);
    result_by_q.set(c_Key1(4), common_term_13);
    // q = 5
    result_by_lpq.set(c_Key3(5, 0, 5), common_term_14);
    result_by_q.set(c_Key1(5), common_term_14);
    // q = 6
    result_by_lpq.set(c_Key3(5, 0, 6), common_term_15);
    result_by_q.set(c_Key1(6), common_term_15);
    // q = 7
    result_by_lpq.set(c_Key3(5, 0, 7), common_term_16);
    result_by_q.set(c_Key1(7), common_term_16);
    // q = 8
    result_by_lpq.set(c_Key3(5, 0, 8), common_term_17);
    result_by_q.set(c_Key1(8), common_term_17);
    // q = 9
    result_by_lpq.set(c_Key3(5, 0, 9), common_term_18);
    result_by_q.set(c_Key1(9), common_term_18);
    // q = 10
    result_by_lpq.set(c_Key3(5, 0, 10), common_term_19);
    result_by_q.set(c_Key1(10), common_term_19);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 0), result_by_q);
    result_by_q.clear();

    // l , p = (5, 1).
    // q = -10
    result_by_lpq.set(c_Key3(5, 1, -10), common_term_20);
    result_by_q.set(c_Key1(-10), common_term_20);
    // q = -9
    result_by_lpq.set(c_Key3(5, 1, -9), common_term_21);
    result_by_q.set(c_Key1(-9), common_term_21);
    // q = -8
    result_by_lpq.set(c_Key3(5, 1, -8), common_term_22);
    result_by_q.set(c_Key1(-8), common_term_22);
    // q = -7
    result_by_lpq.set(c_Key3(5, 1, -7), common_term_23);
    result_by_q.set(c_Key1(-7), common_term_23);
    // q = -6
    result_by_lpq.set(c_Key3(5, 1, -6), common_term_24);
    result_by_q.set(c_Key1(-6), common_term_24);
    // q = -5
    result_by_lpq.set(c_Key3(5, 1, -5), common_term_25);
    result_by_q.set(c_Key1(-5), common_term_25);
    // q = -4
    result_by_lpq.set(c_Key3(5, 1, -4), common_term_26);
    result_by_q.set(c_Key1(-4), common_term_26);
    // q = -3
    result_by_lpq.set(c_Key3(5, 1, -3), common_term_27);
    result_by_q.set(c_Key1(-3), common_term_27);
    // q = -2
    result_by_lpq.set(c_Key3(5, 1, -2), common_term_28);
    result_by_q.set(c_Key1(-2), common_term_28);
    // q = -1
    result_by_lpq.set(c_Key3(5, 1, -1), common_term_29);
    result_by_q.set(c_Key1(-1), common_term_29);
    // q = 0
    result_by_lpq.set(c_Key3(5, 1, 0), common_term_30);
    result_by_q.set(c_Key1(0), common_term_30);
    // q = 1
    result_by_lpq.set(c_Key3(5, 1, 1), common_term_31);
    result_by_q.set(c_Key1(1), common_term_31);
    // q = 2
    result_by_lpq.set(c_Key3(5, 1, 2), common_term_32);
    result_by_q.set(c_Key1(2), common_term_32);
    // q = 3
    result_by_lpq.set(c_Key3(5, 1, 3), common_term_33);
    result_by_q.set(c_Key1(3), common_term_33);
    // q = 4
    result_by_lpq.set(c_Key3(5, 1, 4), common_term_34);
    result_by_q.set(c_Key1(4), common_term_34);
    // q = 5
    result_by_lpq.set(c_Key3(5, 1, 5), common_term_35);
    result_by_q.set(c_Key1(5), common_term_35);
    // q = 6
    result_by_lpq.set(c_Key3(5, 1, 6), common_term_36);
    result_by_q.set(c_Key1(6), common_term_36);
    // q = 7
    result_by_lpq.set(c_Key3(5, 1, 7), common_term_37);
    result_by_q.set(c_Key1(7), common_term_37);
    // q = 8
    result_by_lpq.set(c_Key3(5, 1, 8), common_term_38);
    result_by_q.set(c_Key1(8), common_term_38);
    // q = 9
    result_by_lpq.set(c_Key3(5, 1, 9), common_term_39);
    result_by_q.set(c_Key1(9), common_term_39);
    // q = 10
    result_by_lpq.set(c_Key3(5, 1, 10), common_term_40);
    result_by_q.set(c_Key1(10), common_term_40);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 1), result_by_q);
    result_by_q.clear();

    // l , p = (5, 2).
    // q = -10
    result_by_lpq.set(c_Key3(5, 2, -10), common_term_41);
    result_by_q.set(c_Key1(-10), common_term_41);
    // q = -9
    result_by_lpq.set(c_Key3(5, 2, -9), common_term_42);
    result_by_q.set(c_Key1(-9), common_term_42);
    // q = -8
    result_by_lpq.set(c_Key3(5, 2, -8), common_term_43);
    result_by_q.set(c_Key1(-8), common_term_43);
    // q = -7
    result_by_lpq.set(c_Key3(5, 2, -7), common_term_44);
    result_by_q.set(c_Key1(-7), common_term_44);
    // q = -6
    result_by_lpq.set(c_Key3(5, 2, -6), common_term_45);
    result_by_q.set(c_Key1(-6), common_term_45);
    // q = -5
    result_by_lpq.set(c_Key3(5, 2, -5), common_term_46);
    result_by_q.set(c_Key1(-5), common_term_46);
    // q = -4
    result_by_lpq.set(c_Key3(5, 2, -4), common_term_47);
    result_by_q.set(c_Key1(-4), common_term_47);
    // q = -3
    result_by_lpq.set(c_Key3(5, 2, -3), common_term_48);
    result_by_q.set(c_Key1(-3), common_term_48);
    // q = -2
    result_by_lpq.set(c_Key3(5, 2, -2), common_term_49);
    result_by_q.set(c_Key1(-2), common_term_49);
    // q = -1
    result_by_lpq.set(c_Key3(5, 2, -1), common_term_50);
    result_by_q.set(c_Key1(-1), common_term_50);
    // q = 0
    result_by_lpq.set(c_Key3(5, 2, 0), common_term_51);
    result_by_q.set(c_Key1(0), common_term_51);
    // q = 1
    result_by_lpq.set(c_Key3(5, 2, 1), common_term_52);
    result_by_q.set(c_Key1(1), common_term_52);
    // q = 2
    result_by_lpq.set(c_Key3(5, 2, 2), common_term_53);
    result_by_q.set(c_Key1(2), common_term_53);
    // q = 3
    result_by_lpq.set(c_Key3(5, 2, 3), common_term_54);
    result_by_q.set(c_Key1(3), common_term_54);
    // q = 4
    result_by_lpq.set(c_Key3(5, 2, 4), common_term_55);
    result_by_q.set(c_Key1(4), common_term_55);
    // q = 5
    result_by_lpq.set(c_Key3(5, 2, 5), common_term_56);
    result_by_q.set(c_Key1(5), common_term_56);
    // q = 6
    result_by_lpq.set(c_Key3(5, 2, 6), common_term_57);
    result_by_q.set(c_Key1(6), common_term_57);
    // q = 7
    result_by_lpq.set(c_Key3(5, 2, 7), common_term_58);
    result_by_q.set(c_Key1(7), common_term_58);
    // q = 8
    result_by_lpq.set(c_Key3(5, 2, 8), common_term_59);
    result_by_q.set(c_Key1(8), common_term_59);
    // q = 9
    result_by_lpq.set(c_Key3(5, 2, 9), common_term_60);
    result_by_q.set(c_Key1(9), common_term_60);
    // q = 10
    result_by_lpq.set(c_Key3(5, 2, 10), common_term_61);
    result_by_q.set(c_Key1(10), common_term_61);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 2), result_by_q);
    result_by_q.clear();

    // l , p = (5, 3).
    // q = -10
    result_by_lpq.set(c_Key3(5, 3, -10), common_term_61);
    result_by_q.set(c_Key1(-10), common_term_61);
    // q = -9
    result_by_lpq.set(c_Key3(5, 3, -9), common_term_60);
    result_by_q.set(c_Key1(-9), common_term_60);
    // q = -8
    result_by_lpq.set(c_Key3(5, 3, -8), common_term_59);
    result_by_q.set(c_Key1(-8), common_term_59);
    // q = -7
    result_by_lpq.set(c_Key3(5, 3, -7), common_term_58);
    result_by_q.set(c_Key1(-7), common_term_58);
    // q = -6
    result_by_lpq.set(c_Key3(5, 3, -6), common_term_57);
    result_by_q.set(c_Key1(-6), common_term_57);
    // q = -5
    result_by_lpq.set(c_Key3(5, 3, -5), common_term_56);
    result_by_q.set(c_Key1(-5), common_term_56);
    // q = -4
    result_by_lpq.set(c_Key3(5, 3, -4), common_term_55);
    result_by_q.set(c_Key1(-4), common_term_55);
    // q = -3
    result_by_lpq.set(c_Key3(5, 3, -3), common_term_54);
    result_by_q.set(c_Key1(-3), common_term_54);
    // q = -2
    result_by_lpq.set(c_Key3(5, 3, -2), common_term_53);
    result_by_q.set(c_Key1(-2), common_term_53);
    // q = -1
    result_by_lpq.set(c_Key3(5, 3, -1), common_term_52);
    result_by_q.set(c_Key1(-1), common_term_52);
    // q = 0
    result_by_lpq.set(c_Key3(5, 3, 0), common_term_51);
    result_by_q.set(c_Key1(0), common_term_51);
    // q = 1
    result_by_lpq.set(c_Key3(5, 3, 1), common_term_50);
    result_by_q.set(c_Key1(1), common_term_50);
    // q = 2
    result_by_lpq.set(c_Key3(5, 3, 2), common_term_49);
    result_by_q.set(c_Key1(2), common_term_49);
    // q = 3
    result_by_lpq.set(c_Key3(5, 3, 3), common_term_48);
    result_by_q.set(c_Key1(3), common_term_48);
    // q = 4
    result_by_lpq.set(c_Key3(5, 3, 4), common_term_47);
    result_by_q.set(c_Key1(4), common_term_47);
    // q = 5
    result_by_lpq.set(c_Key3(5, 3, 5), common_term_46);
    result_by_q.set(c_Key1(5), common_term_46);
    // q = 6
    result_by_lpq.set(c_Key3(5, 3, 6), common_term_45);
    result_by_q.set(c_Key1(6), common_term_45);
    // q = 7
    result_by_lpq.set(c_Key3(5, 3, 7), common_term_44);
    result_by_q.set(c_Key1(7), common_term_44);
    // q = 8
    result_by_lpq.set(c_Key3(5, 3, 8), common_term_43);
    result_by_q.set(c_Key1(8), common_term_43);
    // q = 9
    result_by_lpq.set(c_Key3(5, 3, 9), common_term_42);
    result_by_q.set(c_Key1(9), common_term_42);
    // q = 10
    result_by_lpq.set(c_Key3(5, 3, 10), common_term_41);
    result_by_q.set(c_Key1(10), common_term_41);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 3), result_by_q);
    result_by_q.clear();

    // l , p = (5, 4).
    // q = -10
    result_by_lpq.set(c_Key3(5, 4, -10), common_term_40);
    result_by_q.set(c_Key1(-10), common_term_40);
    // q = -9
    result_by_lpq.set(c_Key3(5, 4, -9), common_term_39);
    result_by_q.set(c_Key1(-9), common_term_39);
    // q = -8
    result_by_lpq.set(c_Key3(5, 4, -8), common_term_38);
    result_by_q.set(c_Key1(-8), common_term_38);
    // q = -7
    result_by_lpq.set(c_Key3(5, 4, -7), common_term_37);
    result_by_q.set(c_Key1(-7), common_term_37);
    // q = -6
    result_by_lpq.set(c_Key3(5, 4, -6), common_term_36);
    result_by_q.set(c_Key1(-6), common_term_36);
    // q = -5
    result_by_lpq.set(c_Key3(5, 4, -5), common_term_35);
    result_by_q.set(c_Key1(-5), common_term_35);
    // q = -4
    result_by_lpq.set(c_Key3(5, 4, -4), common_term_34);
    result_by_q.set(c_Key1(-4), common_term_34);
    // q = -3
    result_by_lpq.set(c_Key3(5, 4, -3), common_term_33);
    result_by_q.set(c_Key1(-3), common_term_33);
    // q = -2
    result_by_lpq.set(c_Key3(5, 4, -2), common_term_32);
    result_by_q.set(c_Key1(-2), common_term_32);
    // q = -1
    result_by_lpq.set(c_Key3(5, 4, -1), common_term_31);
    result_by_q.set(c_Key1(-1), common_term_31);
    // q = 0
    result_by_lpq.set(c_Key3(5, 4, 0), common_term_30);
    result_by_q.set(c_Key1(0), common_term_30);
    // q = 1
    result_by_lpq.set(c_Key3(5, 4, 1), common_term_29);
    result_by_q.set(c_Key1(1), common_term_29);
    // q = 2
    result_by_lpq.set(c_Key3(5, 4, 2), common_term_28);
    result_by_q.set(c_Key1(2), common_term_28);
    // q = 3
    result_by_lpq.set(c_Key3(5, 4, 3), common_term_27);
    result_by_q.set(c_Key1(3), common_term_27);
    // q = 4
    result_by_lpq.set(c_Key3(5, 4, 4), common_term_26);
    result_by_q.set(c_Key1(4), common_term_26);
    // q = 5
    result_by_lpq.set(c_Key3(5, 4, 5), common_term_25);
    result_by_q.set(c_Key1(5), common_term_25);
    // q = 6
    result_by_lpq.set(c_Key3(5, 4, 6), common_term_24);
    result_by_q.set(c_Key1(6), common_term_24);
    // q = 7
    result_by_lpq.set(c_Key3(5, 4, 7), common_term_23);
    result_by_q.set(c_Key1(7), common_term_23);
    // q = 8
    result_by_lpq.set(c_Key3(5, 4, 8), common_term_22);
    result_by_q.set(c_Key1(8), common_term_22);
    // q = 9
    result_by_lpq.set(c_Key3(5, 4, 9), common_term_21);
    result_by_q.set(c_Key1(9), common_term_21);
    // q = 10
    result_by_lpq.set(c_Key3(5, 4, 10), common_term_20);
    result_by_q.set(c_Key1(10), common_term_20);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 4), result_by_q);
    result_by_q.clear();

    // l , p = (5, 5).
    // q = -10
    result_by_lpq.set(c_Key3(5, 5, -10), common_term_19);
    result_by_q.set(c_Key1(-10), common_term_19);
    // q = -9
    result_by_lpq.set(c_Key3(5, 5, -9), common_term_18);
    result_by_q.set(c_Key1(-9), common_term_18);
    // q = -8
    result_by_lpq.set(c_Key3(5, 5, -8), common_term_17);
    result_by_q.set(c_Key1(-8), common_term_17);
    // q = -7
    result_by_lpq.set(c_Key3(5, 5, -7), common_term_16);
    result_by_q.set(c_Key1(-7), common_term_16);
    // q = -6
    result_by_lpq.set(c_Key3(5, 5, -6), common_term_15);
    result_by_q.set(c_Key1(-6), common_term_15);
    // q = -5
    result_by_lpq.set(c_Key3(5, 5, -5), common_term_14);
    result_by_q.set(c_Key1(-5), common_term_14);
    // q = -4
    result_by_lpq.set(c_Key3(5, 5, -4), common_term_13);
    result_by_q.set(c_Key1(-4), common_term_13);
    // q = -3
    result_by_lpq.set(c_Key3(5, 5, -3), common_term_12);
    result_by_q.set(c_Key1(-3), common_term_12);
    // q = -2
    result_by_lpq.set(c_Key3(5, 5, -2), common_term_11);
    result_by_q.set(c_Key1(-2), common_term_11);
    // q = -1
    result_by_lpq.set(c_Key3(5, 5, -1), common_term_10);
    result_by_q.set(c_Key1(-1), common_term_10);
    // q = 0
    result_by_lpq.set(c_Key3(5, 5, 0), common_term_9);
    result_by_q.set(c_Key1(0), common_term_9);
    // q = 1
    result_by_lpq.set(c_Key3(5, 5, 1), common_term_8);
    result_by_q.set(c_Key1(1), common_term_8);
    // q = 2
    result_by_lpq.set(c_Key3(5, 5, 2), common_term_7);
    result_by_q.set(c_Key1(2), common_term_7);
    // q = 3
    result_by_lpq.set(c_Key3(5, 5, 3), common_term_6);
    result_by_q.set(c_Key1(3), common_term_6);
    // q = 4
    result_by_lpq.set(c_Key3(5, 5, 4), common_term_5);
    result_by_q.set(c_Key1(4), common_term_5);
    // q = 6
    result_by_lpq.set(c_Key3(5, 5, 6), common_term_4);
    result_by_q.set(c_Key1(6), common_term_4);
    // q = 7
    result_by_lpq.set(c_Key3(5, 5, 7), common_term_3);
    result_by_q.set(c_Key1(7), common_term_3);
    // q = 8
    result_by_lpq.set(c_Key3(5, 5, 8), common_term_2);
    result_by_q.set(c_Key1(8), common_term_2);
    // q = 9
    result_by_lpq.set(c_Key3(5, 5, 9), common_term_1);
    result_by_q.set(c_Key1(9), common_term_1);
    // q = 10
    result_by_lpq.set(c_Key3(5, 5, 10), common_term_0);
    result_by_q.set(c_Key1(10), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 5), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l5_e15(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 5.
    // Functions are taylor expanded and truncated at eccentricity^15.

    // Functions Calculated for l = 5.

    c_IntMap<c_Key3, double> result_by_lpq(184);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(6);
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
    double common_term_0 = 0.02333729166204778*eccentricity_15;
    double common_term_1 = 0.016016471334406767*eccentricity_14;
    double common_term_2 = 0.01655045718537782*eccentricity_15 + 0.010777041888152999*eccentricity_13;
    double common_term_3 = 0.011667415808731163*eccentricity_14 + 0.007054716535511866*eccentricity_12;
    double common_term_4 = 0.01054002637987013*eccentricity_15 + 0.0077663352272727273*eccentricity_13 + 0.0044379058441558442*eccentricity_11;
    double common_term_5 = 0.0064426308455876757*eccentricity_14 + 0.0047783104678042829*eccentricity_12 + 0.0026280707572923556*eccentricity_10;
    double common_term_6 = 0.0039962677925640889*eccentricity_15 + 0.0034696168029501363*eccentricity_13 + 0.0026102292768959436*eccentricity_11 + 0.0014109347442680776*eccentricity_9;
    double common_term_7 = 0.0016919335451993075*eccentricity_14 + 0.001507652827671596*eccentricity_12 + 0.0011653355189732143*eccentricity_10 + 0.00063563755580357143*eccentricity_8;
    double common_term_8 = 0.0004686555326746299*eccentricity_15 + 0.00046257256760728983*eccentricity_13 + 0.0004285163139329806*eccentricity_11 + 0.00034722222222222222*eccentricity_9 + 0.00019841269841269841*eccentricity_7;
    double common_term_9 = 3.8490304520801542e-5*eccentricity_14 + 3.9638763732730747e-5*eccentricity_12 + 3.8885691809275794e-5*eccentricity_10 + 3.410218253968254e-5*eccentricity_8 + 2.1701388888888889e-5*eccentricity_6;
    double common_term_10 = 0.00090825499400610366*eccentricity_14 + 0.0010550213869286593*eccentricity_12 + 0.0012444067253637566*eccentricity_10 + 0.0014946831597222222*eccentricity_8 + 0.0018229166666666667*eccentricity_6 + 0.0026041666666666667*eccentricity_4;
    double common_term_11 = 0.0012744934811630413*eccentricity_15 + 0.0012573371362433862*eccentricity_13 + 0.0011388062169312169*eccentricity_11 + 0.001099537037037037*eccentricity_9 - 0.0041666666666666667*eccentricity_7 + 0.041666666666666667*eccentricity_5 - 0.16666666666666667*eccentricity_3;
    double common_term_12 = -0.0051266195092882429*eccentricity_14 - 0.008580322265625*eccentricity_12 + 0.00580902099609375*eccentricity_10 - 0.1740234375*eccentricity_8 + 0.8525390625*eccentricity_6 - 2.25*eccentricity_4 + 1.125*eccentricity_2;
    double common_term_13 = -0.021366011392983119*eccentricity_15 - 0.068361992945326279*eccentricity_13 + 0.30132233796296296*eccentricity_11 - 2.2220486111111111*eccentricity_9 + 8.3368055555555556*eccentricity_7 - 16.833333333333333*eccentricity_5 + 11.5*eccentricity_3 - 2.0*eccentricity;
    double common_term_14 = -0.63502111012042395*eccentricity_14 + 3.6683625939451618*eccentricity_12 - 17.866980658637153*eccentricity_10 + 53.395758734809028*eccentricity_8 - 90.342881944444444*eccentricity_6 + 66.484375*eccentricity_4 - 17.5*eccentricity_2 + 1.0;
    double common_term_15 = -5.4744128667091837*eccentricity_15 + 27.873044084821429*eccentricity_13 - 105.696796875*eccentricity_11 + 262.5140625*eccentricity_9 - 388.625*eccentricity_7 + 286.5*eccentricity_5 - 88.5*eccentricity_3 + 8.0*eccentricity;
    double common_term_16 = 160.33849005591722*eccentricity_14 - 503.92212854456019*eccentricity_12 + 1075.0551734641746*eccentricity_10 - 1426.3118272569444*eccentricity_8 + 1024.1761067708333*eccentricity_6 - 338.91666666666667*eccentricity_4 + 37.375*eccentricity_2;
    double common_term_17 = 759.82512425289046*eccentricity_15 - 2049.3854050925926*eccentricity_13 + 3847.0054274140212*eccentricity_11 - 4644.6884259259259*eccentricity_9 + 3210.0166666666667*eccentricity_7 - 1090.1666666666667*eccentricity_5 + 133.16666666666667*eccentricity_3;
    double common_term_18 = -7378.0934716027124*eccentricity_14 + 12404.238924591882*eccentricity_12 - 13776.578367396763*eccentricity_10 + 9117.8461669921875*eccentricity_8 - 3105.83671875*eccentricity_6 + 400.6796875*eccentricity_4;
    double common_term_19 = -24118.051818381191*eccentricity_15 + 36805.308178380824*eccentricity_13 - 37910.616887814153*eccentricity_11 + 23983.827504960317*eccentricity_9 - 8086.7361111111111*eccentricity_7 + 1072.4583333333333*eccentricity_5;
    double common_term_20 = 102024.90101705193*eccentricity_14 - 98093.605275043216*eccentricity_12 + 59306.949686698308*eccentricity_10 - 19637.63640562996*eccentricity_8 + 2633.1491102430556*eccentricity_6;
    double common_term_21 = 267204.62059265929*eccentricity_15 - 241078.85928571429*eccentricity_13 + 139375.73253348214*eccentricity_11 - 45101.5265625*eccentricity_9 + 6049.0116071428571*eccentricity_7;
    double common_term_22 = -567140.76834504931*eccentricity_14 + 313848.34123679371*eccentricity_12 - 98957.484222180671*eccentricity_10 + 13179.87816617451*eccentricity_8;
    double common_term_23 = -1284970.2308093875*eccentricity_15 + 681474.58250646658*eccentricity_13 - 208988.53118537809*eccentricity_11 + 27504.092611882716*eccentricity_9;
    double common_term_24 = 1433984.7381258379*eccentricity_14 - 427283.52163934088*eccentricity_12 + 55372.913824353899*eccentricity_10;
    double common_term_25 = 2935945.1361715278*eccentricity_15 - 849564.73169093381*eccentricity_13 + 108149.79697137546*eccentricity_11;
    double common_term_26 = -1648688.3416402937*eccentricity_14 + 205812.7946639314*eccentricity_12;
    double common_term_27 = -3132020.7740035927*eccentricity_15 + 382953.25728189779*eccentricity_13;
    double common_term_28 = 698665.82920653284*eccentricity_14;
    double common_term_29 = 1252708.9143546075*eccentricity_15;
    double common_term_30 = 21.091149702836226*eccentricity_15;
    double common_term_31 = 15.408737763597732*eccentricity_14;
    double common_term_32 = 21.465307221747884*eccentricity_15 + 11.257742673270017*eccentricity_13;
    double common_term_33 = 17.835125904545322*eccentricity_14 + 8.2256967558179583*eccentricity_12;
    double common_term_34 = 29.942290013347956*eccentricity_15 + 14.607572256543611*eccentricity_13 + 6.0111507560726311*eccentricity_11;
    double common_term_35 = 24.964970074140129*eccentricity_14 + 11.829020482329809*eccentricity_12 + 4.3938146464029948*eccentricity_10;
    double common_term_36 = 38.251415635146104*eccentricity_15 + 20.724025720373377*eccentricity_13 + 9.49171875*eccentricity_11 + 3.2127232142857143*eccentricity_9;
    double common_term_37 = 32.425630900461131*eccentricity_14 + 17.125299496625466*eccentricity_12 + 7.5595372921575314*eccentricity_10 + 2.3502705891927083*eccentricity_8;
    double common_term_38 = 47.331421073802259*eccentricity_15 + 27.380366236772487*eccentricity_13 + 14.087541335978836*eccentricity_11 + 5.9837797619047619*eccentricity_9 + 1.7205357142857143*eccentricity_7;
    double common_term_39 = 40.606042443854468*eccentricity_14 + 23.031711687360491*eccentricity_12 + 11.538018035888672*eccentricity_10 + 4.7125279017857143*eccentricity_8 + 1.2607421875*eccentricity_6;
    double common_term_40 = 57.137932308889991*eccentricity_15 + 34.721210868606702*eccentricity_13 + 19.301029265873016*eccentricity_11 + 9.410639880952381*eccentricity_9 + 3.6958333333333333*eccentricity_7 + 0.925*eccentricity_5;
    double common_term_41 = 49.485521528195242*eccentricity_14 + 29.591459706472972*eccentricity_12 + 16.115094284784226*eccentricity_10 + 7.6452555338541667*eccentricity_8 + 2.88828125*eccentricity_6 + 0.6796875*eccentricity_4;
    double common_term_42 = 0.5*eccentricity_3*std::pow(1.0 - eccentricity_2, -4.5);
    double common_term_43 = 59.054041014140245*eccentricity_14 + 36.793448970734127*eccentricity_12 + 21.283941819932726*eccentricity_10 + 11.114388020833333*eccentricity_8 + 4.9912109375*eccentricity_6 + 1.75*eccentricity_4 + 0.375*eccentricity_2;
    double common_term_44 = 78.891815424520503*eccentricity_15 + 51.408773251488095*eccentricity_13 + 31.589105902777778*eccentricity_11 + 17.966145833333333*eccentricity_9 + 9.1875*eccentricity_7 + 4.0*eccentricity_5 + 1.5*eccentricity_3;
    double common_term_45 = 69.31039448524008*eccentricity_14 + 44.636981925964355*eccentricity_12 + 27.040718994140625*eccentricity_10 + 15.16680908203125*eccentricity_8 + 7.2265625*eccentricity_6 + 4.734375*eccentricity_4 - 1.5*eccentricity_2 + 1.0;
    double common_term_46 = 90.831207124994095*eccentricity_15 + 60.754408895502646*eccentricity_13 + 38.584227430555556*eccentricity_11 + 23.668229166666667*eccentricity_9 + 9.6770833333333333*eccentricity_7 + 15.5*eccentricity_5 - 10.5*eccentricity_3 + 6.0*eccentricity;
    double common_term_47 = 80.343051059693886*eccentricity_14 + 52.36621335953001*eccentricity_12 + 37.758710225423177*eccentricity_10 + 2.5192057291666667*eccentricity_8 + 51.7041015625*eccentricity_6 - 44.25*eccentricity_4 + 22.125*eccentricity_2;
    double common_term_48 = 104.32126953125*eccentricity_15 + 65.505837053571429*eccentricity_13 + 70.158314732142857*eccentricity_11 - 45.9796875*eccentricity_9 + 164.025*eccentricity_7 - 145.125*eccentricity_5 + 64.5*eccentricity_3;
    double common_term_49 = 63.870867509900788*eccentricity_14 + 166.38394224908617*eccentricity_12 - 231.58345404730903*eccentricity_10 + 480.89112141927083*eccentricity_8 - 407.71484375*eccentricity_6 + 163.0859375*eccentricity_4;
    double common_term_50 = -6.6140546840553351*eccentricity_15 + 469.22716910548942*eccentricity_13 - 812.08832465277778*eccentricity_11 + 1302.5502976190476*eccentricity_9 - 1029.7375*eccentricity_7 + 374.45*eccentricity_5;
    double common_term_51 = 1382.9710076344013*eccentricity_14 - 2426.4057120186942*eccentricity_12 + 3289.3376407078334*eccentricity_10 - 2404.9318359375*eccentricity_8 + 801.6439453125*eccentricity_6;
    double common_term_52 = 3964.3006371220889*eccentricity_15 - 6560.9110105268959*eccentricity_13 + 7825.6551235945767*eccentricity_11 - 5288.5431547619048*eccentricity_9 + 1627.2440476190476*eccentricity_7;
    double common_term_53 = -16504.128205581729*eccentricity_14 + 17704.166441436142*eccentricity_12 - 11086.362901814779*eccentricity_10 + 3167.5322576613653*eccentricity_8;
    double common_term_54 = -39239.608693752537*eccentricity_15 + 38382.229918831169*eccentricity_13 - 22351.217912946429*eccentricity_11 + 5960.3247767857143*eccentricity_9;
    double common_term_55 = 80243.948315284375*eccentricity_14 - 43623.675519422097*eccentricity_12 + 10905.918879855524*eccentricity_10;
    double common_term_56 = 162603.91326242611*eccentricity_15 - 82837.988055095223*eccentricity_13 + 19491.291307081229*eccentricity_11;
    double common_term_57 = -153647.68956916828*eccentricity_14 + 34143.931015370728*eccentricity_12;
    double common_term_58 = -279235.16068955808*eccentricity_15 + 58786.349618854332*eccentricity_13;
    double common_term_59 = 99699.59406914746*eccentricity_14;
    double common_term_60 = 166860.24606244119*eccentricity_15;
    double common_term_61 = 895.6770696572522*eccentricity_15;
    double common_term_62 = 615.8782870739856*eccentricity_14;
    double common_term_63 = 138.04397284286584*eccentricity_15 + 421.82299405672453*eccentricity_13;
    double common_term_64 = 178.41581534061428*eccentricity_14 + 287.62977158715766*eccentricity_12;
    double common_term_65 = 429.61432180327308*eccentricity_15 + 179.28346223165017*eccentricity_13 + 195.1342063366803*eccentricity_11;
    double common_term_66 = 341.44604407120061*eccentricity_14 + 161.079640740729*eccentricity_12 + 131.61114809308733*eccentricity_10;
    double common_term_67 = 444.23567553653961*eccentricity_15 + 272.95488521123938*eccentricity_13 + 135.5390270888448*eccentricity_11 + 88.163544422398589*eccentricity_9;
    double common_term_68 = 368.32778909907642*eccentricity_14 + 217.86436910017037*eccentricity_12 + 109.04968559359327*eccentricity_10 + 58.584195285373264*eccentricity_8;
    double common_term_69 = 487.12988477830763*eccentricity_15 + 302.38700613839286*eccentricity_13 + 172.72349330357143*eccentricity_11 + 84.790848214285714*eccentricity_9 + 38.552678571428571*eccentricity_7;
    double common_term_70 = 405.78072080117755*eccentricity_14 + 245.76078824195211*eccentricity_12 + 135.52102582416837*eccentricity_10 + 64.078047495039683*eccentricity_8 + 25.069552951388889*eccentricity_6;
    double common_term_71 = 528.39256430213385*eccentricity_15 + 335.18820862498163*eccentricity_13 + 197.59574239417989*eccentricity_11 + 104.94945436507937*eccentricity_9 + 47.195138888888889*eccentricity_7 + 16.058333333333333*eccentricity_5;
    double common_term_72 = 442.34680258410318*eccentricity_14 + 274.3350358077458*eccentricity_12 + 157.00033656529018*eccentricity_10 + 80.0341552734375*eccentricity_8 + 33.89765625*eccentricity_6 + 10.0859375*eccentricity_4;
    double common_term_73 = 569.00009836336316*eccentricity_15 + 367.3229865244709*eccentricity_13 + 222.25277364417989*eccentricity_11 + 123.10688657407407*eccentricity_9 + 59.960416666666667*eccentricity_7 + 23.708333333333333*eccentricity_5 + 6.1666666666666667*eccentricity_3;
    double common_term_74 = 478.27117265849938*eccentricity_14 + 302.30942661184482*eccentricity_12 + 178.02036545364945*eccentricity_10 + 95.092730034722222*eccentricity_8 + 44.002278645833333*eccentricity_6 + 16.083333333333333*eccentricity_4 + 3.625*eccentricity_2;
    double common_term_75 = std::pow(1.0 - eccentricity_2, -4.5)*(1.5*eccentricity_3 + 2.0*eccentricity);
    double common_term_76 = 513.54641791282898*eccentricity_14 + 329.68472556549826*eccentricity_12 + 198.49574069552951*eccentricity_10 + 109.66821967230903*eccentricity_8 + 53.677951388888889*eccentricity_6 + 21.859375*eccentricity_4 + 6.5*eccentricity_2 + 1.0;
    double common_term_77 = 648.20072550424973*eccentricity_15 + 429.72375657930996*eccentricity_13 + 269.85374421296296*eccentricity_11 + 157.90815972222222*eccentricity_9 + 83.961805555555556*eccentricity_7 + 38.916666666666667*eccentricity_5 + 14.5*eccentricity_3 + 4.0*eccentricity;
    double common_term_78 = 548.17235871945109*eccentricity_14 + 356.46074497767857*eccentricity_12 + 218.42696685791016*eccentricity_10 + 123.7505859375*eccentricity_8 + 63.0263671875*eccentricity_6 + 26.75*eccentricity_4 + 10.875*eccentricity_2;
    double common_term_79 = 686.79135353418792*eccentricity_15 + 459.98591817542989*eccentricity_13 + 292.80587590939153*eccentricity_11 + 174.4318287037037*eccentricity_9 + 95.977083333333333*eccentricity_7 + 43.041666666666667*eccentricity_5 + 24.833333333333333*eccentricity_3;
    double common_term_80 = 582.13725769954314*eccentricity_14 + 382.74131503054705*eccentricity_12 + 237.09027357959243*eccentricity_10 + 140.88095431857639*eccentricity_8 + 60.779947916666667*eccentricity_6 + 51.221354166666667*eccentricity_4;
    double common_term_81 = 724.61253892299107*eccentricity_15 + 490.30125558035714*eccentricity_13 + 311.54815848214286*eccentricity_11 + 204.35424107142857*eccentricity_9 + 72.31875*eccentricity_7 + 98.775*eccentricity_5;
    double common_term_82 = 618.84262999488003*eccentricity_14 + 393.72408189812805*eccentricity_12 + 301.28983722262912*eccentricity_10 + 60.194509548611111*eccentricity_8 + 181.54624565972222*eccentricity_6;
    double common_term_83 = 775.77064236998373*eccentricity_15 + 469.16829952968842*eccentricity_13 + 464.57999820877425*eccentricity_11 - 11.144667658730159*eccentricity_9 + 321.84662698412698*eccentricity_7;
    double common_term_84 = 498.47863376902295*eccentricity_14 + 763.95394112450736*eccentricity_12 - 209.35917576381138*eccentricity_10 + 554.7166024344308*eccentricity_8;
    double common_term_85 = 386.31665203684428*eccentricity_15 + 1341.1313508660514*eccentricity_13 - 657.54596905313051*eccentricity_11 + 934.67407958553792*eccentricity_9;
    double common_term_86 = 2473.528428908686*eccentricity_14 - 1572.4675775587832*eccentricity_12 + 1545.8537958783616*eccentricity_10;
    double common_term_87 = 4687.0080340118866*eccentricity_15 - 3326.4069472909903*eccentricity_13 + 2517.1731057224026*eccentricity_11;
    double common_term_88 = -6546.1679942372922*eccentricity_14 + 4044.9237368728242*eccentricity_12;
    double common_term_89 = -12270.047500231993*eccentricity_15 + 6426.300886272694*eccentricity_13;
    double common_term_90 = 10108.995635604223*eccentricity_14;
    double common_term_91 = 15764.316954849364*eccentricity_15;

    c_IntMap<c_Key1, double> result_by_q(31);

    // eccentricity function by mode:
    // l , p = (5, 0).
    // q = -15
    result_by_lpq.set(c_Key3(5, 0, -15), common_term_0);
    result_by_q.set(c_Key1(-15), common_term_0);
    // q = -14
    result_by_lpq.set(c_Key3(5, 0, -14), common_term_1);
    result_by_q.set(c_Key1(-14), common_term_1);
    // q = -13
    result_by_lpq.set(c_Key3(5, 0, -13), common_term_2);
    result_by_q.set(c_Key1(-13), common_term_2);
    // q = -12
    result_by_lpq.set(c_Key3(5, 0, -12), common_term_3);
    result_by_q.set(c_Key1(-12), common_term_3);
    // q = -11
    result_by_lpq.set(c_Key3(5, 0, -11), common_term_4);
    result_by_q.set(c_Key1(-11), common_term_4);
    // q = -10
    result_by_lpq.set(c_Key3(5, 0, -10), common_term_5);
    result_by_q.set(c_Key1(-10), common_term_5);
    // q = -9
    result_by_lpq.set(c_Key3(5, 0, -9), common_term_6);
    result_by_q.set(c_Key1(-9), common_term_6);
    // q = -8
    result_by_lpq.set(c_Key3(5, 0, -8), common_term_7);
    result_by_q.set(c_Key1(-8), common_term_7);
    // q = -7
    result_by_lpq.set(c_Key3(5, 0, -7), common_term_8);
    result_by_q.set(c_Key1(-7), common_term_8);
    // q = -6
    result_by_lpq.set(c_Key3(5, 0, -6), common_term_9);
    result_by_q.set(c_Key1(-6), common_term_9);
    // q = -4
    result_by_lpq.set(c_Key3(5, 0, -4), common_term_10);
    result_by_q.set(c_Key1(-4), common_term_10);
    // q = -3
    result_by_lpq.set(c_Key3(5, 0, -3), common_term_11);
    result_by_q.set(c_Key1(-3), common_term_11);
    // q = -2
    result_by_lpq.set(c_Key3(5, 0, -2), common_term_12);
    result_by_q.set(c_Key1(-2), common_term_12);
    // q = -1
    result_by_lpq.set(c_Key3(5, 0, -1), common_term_13);
    result_by_q.set(c_Key1(-1), common_term_13);
    // q = 0
    result_by_lpq.set(c_Key3(5, 0, 0), common_term_14);
    result_by_q.set(c_Key1(0), common_term_14);
    // q = 1
    result_by_lpq.set(c_Key3(5, 0, 1), common_term_15);
    result_by_q.set(c_Key1(1), common_term_15);
    // q = 2
    result_by_lpq.set(c_Key3(5, 0, 2), common_term_16);
    result_by_q.set(c_Key1(2), common_term_16);
    // q = 3
    result_by_lpq.set(c_Key3(5, 0, 3), common_term_17);
    result_by_q.set(c_Key1(3), common_term_17);
    // q = 4
    result_by_lpq.set(c_Key3(5, 0, 4), common_term_18);
    result_by_q.set(c_Key1(4), common_term_18);
    // q = 5
    result_by_lpq.set(c_Key3(5, 0, 5), common_term_19);
    result_by_q.set(c_Key1(5), common_term_19);
    // q = 6
    result_by_lpq.set(c_Key3(5, 0, 6), common_term_20);
    result_by_q.set(c_Key1(6), common_term_20);
    // q = 7
    result_by_lpq.set(c_Key3(5, 0, 7), common_term_21);
    result_by_q.set(c_Key1(7), common_term_21);
    // q = 8
    result_by_lpq.set(c_Key3(5, 0, 8), common_term_22);
    result_by_q.set(c_Key1(8), common_term_22);
    // q = 9
    result_by_lpq.set(c_Key3(5, 0, 9), common_term_23);
    result_by_q.set(c_Key1(9), common_term_23);
    // q = 10
    result_by_lpq.set(c_Key3(5, 0, 10), common_term_24);
    result_by_q.set(c_Key1(10), common_term_24);
    // q = 11
    result_by_lpq.set(c_Key3(5, 0, 11), common_term_25);
    result_by_q.set(c_Key1(11), common_term_25);
    // q = 12
    result_by_lpq.set(c_Key3(5, 0, 12), common_term_26);
    result_by_q.set(c_Key1(12), common_term_26);
    // q = 13
    result_by_lpq.set(c_Key3(5, 0, 13), common_term_27);
    result_by_q.set(c_Key1(13), common_term_27);
    // q = 14
    result_by_lpq.set(c_Key3(5, 0, 14), common_term_28);
    result_by_q.set(c_Key1(14), common_term_28);
    // q = 15
    result_by_lpq.set(c_Key3(5, 0, 15), common_term_29);
    result_by_q.set(c_Key1(15), common_term_29);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 0), result_by_q);
    result_by_q.clear();

    // l , p = (5, 1).
    // q = -15
    result_by_lpq.set(c_Key3(5, 1, -15), common_term_30);
    result_by_q.set(c_Key1(-15), common_term_30);
    // q = -14
    result_by_lpq.set(c_Key3(5, 1, -14), common_term_31);
    result_by_q.set(c_Key1(-14), common_term_31);
    // q = -13
    result_by_lpq.set(c_Key3(5, 1, -13), common_term_32);
    result_by_q.set(c_Key1(-13), common_term_32);
    // q = -12
    result_by_lpq.set(c_Key3(5, 1, -12), common_term_33);
    result_by_q.set(c_Key1(-12), common_term_33);
    // q = -11
    result_by_lpq.set(c_Key3(5, 1, -11), common_term_34);
    result_by_q.set(c_Key1(-11), common_term_34);
    // q = -10
    result_by_lpq.set(c_Key3(5, 1, -10), common_term_35);
    result_by_q.set(c_Key1(-10), common_term_35);
    // q = -9
    result_by_lpq.set(c_Key3(5, 1, -9), common_term_36);
    result_by_q.set(c_Key1(-9), common_term_36);
    // q = -8
    result_by_lpq.set(c_Key3(5, 1, -8), common_term_37);
    result_by_q.set(c_Key1(-8), common_term_37);
    // q = -7
    result_by_lpq.set(c_Key3(5, 1, -7), common_term_38);
    result_by_q.set(c_Key1(-7), common_term_38);
    // q = -6
    result_by_lpq.set(c_Key3(5, 1, -6), common_term_39);
    result_by_q.set(c_Key1(-6), common_term_39);
    // q = -5
    result_by_lpq.set(c_Key3(5, 1, -5), common_term_40);
    result_by_q.set(c_Key1(-5), common_term_40);
    // q = -4
    result_by_lpq.set(c_Key3(5, 1, -4), common_term_41);
    result_by_q.set(c_Key1(-4), common_term_41);
    // q = -3
    result_by_lpq.set(c_Key3(5, 1, -3), common_term_42);
    result_by_q.set(c_Key1(-3), common_term_42);
    // q = -2
    result_by_lpq.set(c_Key3(5, 1, -2), common_term_43);
    result_by_q.set(c_Key1(-2), common_term_43);
    // q = -1
    result_by_lpq.set(c_Key3(5, 1, -1), common_term_44);
    result_by_q.set(c_Key1(-1), common_term_44);
    // q = 0
    result_by_lpq.set(c_Key3(5, 1, 0), common_term_45);
    result_by_q.set(c_Key1(0), common_term_45);
    // q = 1
    result_by_lpq.set(c_Key3(5, 1, 1), common_term_46);
    result_by_q.set(c_Key1(1), common_term_46);
    // q = 2
    result_by_lpq.set(c_Key3(5, 1, 2), common_term_47);
    result_by_q.set(c_Key1(2), common_term_47);
    // q = 3
    result_by_lpq.set(c_Key3(5, 1, 3), common_term_48);
    result_by_q.set(c_Key1(3), common_term_48);
    // q = 4
    result_by_lpq.set(c_Key3(5, 1, 4), common_term_49);
    result_by_q.set(c_Key1(4), common_term_49);
    // q = 5
    result_by_lpq.set(c_Key3(5, 1, 5), common_term_50);
    result_by_q.set(c_Key1(5), common_term_50);
    // q = 6
    result_by_lpq.set(c_Key3(5, 1, 6), common_term_51);
    result_by_q.set(c_Key1(6), common_term_51);
    // q = 7
    result_by_lpq.set(c_Key3(5, 1, 7), common_term_52);
    result_by_q.set(c_Key1(7), common_term_52);
    // q = 8
    result_by_lpq.set(c_Key3(5, 1, 8), common_term_53);
    result_by_q.set(c_Key1(8), common_term_53);
    // q = 9
    result_by_lpq.set(c_Key3(5, 1, 9), common_term_54);
    result_by_q.set(c_Key1(9), common_term_54);
    // q = 10
    result_by_lpq.set(c_Key3(5, 1, 10), common_term_55);
    result_by_q.set(c_Key1(10), common_term_55);
    // q = 11
    result_by_lpq.set(c_Key3(5, 1, 11), common_term_56);
    result_by_q.set(c_Key1(11), common_term_56);
    // q = 12
    result_by_lpq.set(c_Key3(5, 1, 12), common_term_57);
    result_by_q.set(c_Key1(12), common_term_57);
    // q = 13
    result_by_lpq.set(c_Key3(5, 1, 13), common_term_58);
    result_by_q.set(c_Key1(13), common_term_58);
    // q = 14
    result_by_lpq.set(c_Key3(5, 1, 14), common_term_59);
    result_by_q.set(c_Key1(14), common_term_59);
    // q = 15
    result_by_lpq.set(c_Key3(5, 1, 15), common_term_60);
    result_by_q.set(c_Key1(15), common_term_60);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 1), result_by_q);
    result_by_q.clear();

    // l , p = (5, 2).
    // q = -15
    result_by_lpq.set(c_Key3(5, 2, -15), common_term_61);
    result_by_q.set(c_Key1(-15), common_term_61);
    // q = -14
    result_by_lpq.set(c_Key3(5, 2, -14), common_term_62);
    result_by_q.set(c_Key1(-14), common_term_62);
    // q = -13
    result_by_lpq.set(c_Key3(5, 2, -13), common_term_63);
    result_by_q.set(c_Key1(-13), common_term_63);
    // q = -12
    result_by_lpq.set(c_Key3(5, 2, -12), common_term_64);
    result_by_q.set(c_Key1(-12), common_term_64);
    // q = -11
    result_by_lpq.set(c_Key3(5, 2, -11), common_term_65);
    result_by_q.set(c_Key1(-11), common_term_65);
    // q = -10
    result_by_lpq.set(c_Key3(5, 2, -10), common_term_66);
    result_by_q.set(c_Key1(-10), common_term_66);
    // q = -9
    result_by_lpq.set(c_Key3(5, 2, -9), common_term_67);
    result_by_q.set(c_Key1(-9), common_term_67);
    // q = -8
    result_by_lpq.set(c_Key3(5, 2, -8), common_term_68);
    result_by_q.set(c_Key1(-8), common_term_68);
    // q = -7
    result_by_lpq.set(c_Key3(5, 2, -7), common_term_69);
    result_by_q.set(c_Key1(-7), common_term_69);
    // q = -6
    result_by_lpq.set(c_Key3(5, 2, -6), common_term_70);
    result_by_q.set(c_Key1(-6), common_term_70);
    // q = -5
    result_by_lpq.set(c_Key3(5, 2, -5), common_term_71);
    result_by_q.set(c_Key1(-5), common_term_71);
    // q = -4
    result_by_lpq.set(c_Key3(5, 2, -4), common_term_72);
    result_by_q.set(c_Key1(-4), common_term_72);
    // q = -3
    result_by_lpq.set(c_Key3(5, 2, -3), common_term_73);
    result_by_q.set(c_Key1(-3), common_term_73);
    // q = -2
    result_by_lpq.set(c_Key3(5, 2, -2), common_term_74);
    result_by_q.set(c_Key1(-2), common_term_74);
    // q = -1
    result_by_lpq.set(c_Key3(5, 2, -1), common_term_75);
    result_by_q.set(c_Key1(-1), common_term_75);
    // q = 0
    result_by_lpq.set(c_Key3(5, 2, 0), common_term_76);
    result_by_q.set(c_Key1(0), common_term_76);
    // q = 1
    result_by_lpq.set(c_Key3(5, 2, 1), common_term_77);
    result_by_q.set(c_Key1(1), common_term_77);
    // q = 2
    result_by_lpq.set(c_Key3(5, 2, 2), common_term_78);
    result_by_q.set(c_Key1(2), common_term_78);
    // q = 3
    result_by_lpq.set(c_Key3(5, 2, 3), common_term_79);
    result_by_q.set(c_Key1(3), common_term_79);
    // q = 4
    result_by_lpq.set(c_Key3(5, 2, 4), common_term_80);
    result_by_q.set(c_Key1(4), common_term_80);
    // q = 5
    result_by_lpq.set(c_Key3(5, 2, 5), common_term_81);
    result_by_q.set(c_Key1(5), common_term_81);
    // q = 6
    result_by_lpq.set(c_Key3(5, 2, 6), common_term_82);
    result_by_q.set(c_Key1(6), common_term_82);
    // q = 7
    result_by_lpq.set(c_Key3(5, 2, 7), common_term_83);
    result_by_q.set(c_Key1(7), common_term_83);
    // q = 8
    result_by_lpq.set(c_Key3(5, 2, 8), common_term_84);
    result_by_q.set(c_Key1(8), common_term_84);
    // q = 9
    result_by_lpq.set(c_Key3(5, 2, 9), common_term_85);
    result_by_q.set(c_Key1(9), common_term_85);
    // q = 10
    result_by_lpq.set(c_Key3(5, 2, 10), common_term_86);
    result_by_q.set(c_Key1(10), common_term_86);
    // q = 11
    result_by_lpq.set(c_Key3(5, 2, 11), common_term_87);
    result_by_q.set(c_Key1(11), common_term_87);
    // q = 12
    result_by_lpq.set(c_Key3(5, 2, 12), common_term_88);
    result_by_q.set(c_Key1(12), common_term_88);
    // q = 13
    result_by_lpq.set(c_Key3(5, 2, 13), common_term_89);
    result_by_q.set(c_Key1(13), common_term_89);
    // q = 14
    result_by_lpq.set(c_Key3(5, 2, 14), common_term_90);
    result_by_q.set(c_Key1(14), common_term_90);
    // q = 15
    result_by_lpq.set(c_Key3(5, 2, 15), common_term_91);
    result_by_q.set(c_Key1(15), common_term_91);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 2), result_by_q);
    result_by_q.clear();

    // l , p = (5, 3).
    // q = -15
    result_by_lpq.set(c_Key3(5, 3, -15), common_term_91);
    result_by_q.set(c_Key1(-15), common_term_91);
    // q = -14
    result_by_lpq.set(c_Key3(5, 3, -14), common_term_90);
    result_by_q.set(c_Key1(-14), common_term_90);
    // q = -13
    result_by_lpq.set(c_Key3(5, 3, -13), common_term_89);
    result_by_q.set(c_Key1(-13), common_term_89);
    // q = -12
    result_by_lpq.set(c_Key3(5, 3, -12), common_term_88);
    result_by_q.set(c_Key1(-12), common_term_88);
    // q = -11
    result_by_lpq.set(c_Key3(5, 3, -11), common_term_87);
    result_by_q.set(c_Key1(-11), common_term_87);
    // q = -10
    result_by_lpq.set(c_Key3(5, 3, -10), common_term_86);
    result_by_q.set(c_Key1(-10), common_term_86);
    // q = -9
    result_by_lpq.set(c_Key3(5, 3, -9), common_term_85);
    result_by_q.set(c_Key1(-9), common_term_85);
    // q = -8
    result_by_lpq.set(c_Key3(5, 3, -8), common_term_84);
    result_by_q.set(c_Key1(-8), common_term_84);
    // q = -7
    result_by_lpq.set(c_Key3(5, 3, -7), common_term_83);
    result_by_q.set(c_Key1(-7), common_term_83);
    // q = -6
    result_by_lpq.set(c_Key3(5, 3, -6), common_term_82);
    result_by_q.set(c_Key1(-6), common_term_82);
    // q = -5
    result_by_lpq.set(c_Key3(5, 3, -5), common_term_81);
    result_by_q.set(c_Key1(-5), common_term_81);
    // q = -4
    result_by_lpq.set(c_Key3(5, 3, -4), common_term_80);
    result_by_q.set(c_Key1(-4), common_term_80);
    // q = -3
    result_by_lpq.set(c_Key3(5, 3, -3), common_term_79);
    result_by_q.set(c_Key1(-3), common_term_79);
    // q = -2
    result_by_lpq.set(c_Key3(5, 3, -2), common_term_78);
    result_by_q.set(c_Key1(-2), common_term_78);
    // q = -1
    result_by_lpq.set(c_Key3(5, 3, -1), common_term_77);
    result_by_q.set(c_Key1(-1), common_term_77);
    // q = 0
    result_by_lpq.set(c_Key3(5, 3, 0), common_term_76);
    result_by_q.set(c_Key1(0), common_term_76);
    // q = 1
    result_by_lpq.set(c_Key3(5, 3, 1), common_term_75);
    result_by_q.set(c_Key1(1), common_term_75);
    // q = 2
    result_by_lpq.set(c_Key3(5, 3, 2), common_term_74);
    result_by_q.set(c_Key1(2), common_term_74);
    // q = 3
    result_by_lpq.set(c_Key3(5, 3, 3), common_term_73);
    result_by_q.set(c_Key1(3), common_term_73);
    // q = 4
    result_by_lpq.set(c_Key3(5, 3, 4), common_term_72);
    result_by_q.set(c_Key1(4), common_term_72);
    // q = 5
    result_by_lpq.set(c_Key3(5, 3, 5), common_term_71);
    result_by_q.set(c_Key1(5), common_term_71);
    // q = 6
    result_by_lpq.set(c_Key3(5, 3, 6), common_term_70);
    result_by_q.set(c_Key1(6), common_term_70);
    // q = 7
    result_by_lpq.set(c_Key3(5, 3, 7), common_term_69);
    result_by_q.set(c_Key1(7), common_term_69);
    // q = 8
    result_by_lpq.set(c_Key3(5, 3, 8), common_term_68);
    result_by_q.set(c_Key1(8), common_term_68);
    // q = 9
    result_by_lpq.set(c_Key3(5, 3, 9), common_term_67);
    result_by_q.set(c_Key1(9), common_term_67);
    // q = 10
    result_by_lpq.set(c_Key3(5, 3, 10), common_term_66);
    result_by_q.set(c_Key1(10), common_term_66);
    // q = 11
    result_by_lpq.set(c_Key3(5, 3, 11), common_term_65);
    result_by_q.set(c_Key1(11), common_term_65);
    // q = 12
    result_by_lpq.set(c_Key3(5, 3, 12), common_term_64);
    result_by_q.set(c_Key1(12), common_term_64);
    // q = 13
    result_by_lpq.set(c_Key3(5, 3, 13), common_term_63);
    result_by_q.set(c_Key1(13), common_term_63);
    // q = 14
    result_by_lpq.set(c_Key3(5, 3, 14), common_term_62);
    result_by_q.set(c_Key1(14), common_term_62);
    // q = 15
    result_by_lpq.set(c_Key3(5, 3, 15), common_term_61);
    result_by_q.set(c_Key1(15), common_term_61);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 3), result_by_q);
    result_by_q.clear();

    // l , p = (5, 4).
    // q = -15
    result_by_lpq.set(c_Key3(5, 4, -15), common_term_60);
    result_by_q.set(c_Key1(-15), common_term_60);
    // q = -14
    result_by_lpq.set(c_Key3(5, 4, -14), common_term_59);
    result_by_q.set(c_Key1(-14), common_term_59);
    // q = -13
    result_by_lpq.set(c_Key3(5, 4, -13), common_term_58);
    result_by_q.set(c_Key1(-13), common_term_58);
    // q = -12
    result_by_lpq.set(c_Key3(5, 4, -12), common_term_57);
    result_by_q.set(c_Key1(-12), common_term_57);
    // q = -11
    result_by_lpq.set(c_Key3(5, 4, -11), common_term_56);
    result_by_q.set(c_Key1(-11), common_term_56);
    // q = -10
    result_by_lpq.set(c_Key3(5, 4, -10), common_term_55);
    result_by_q.set(c_Key1(-10), common_term_55);
    // q = -9
    result_by_lpq.set(c_Key3(5, 4, -9), common_term_54);
    result_by_q.set(c_Key1(-9), common_term_54);
    // q = -8
    result_by_lpq.set(c_Key3(5, 4, -8), common_term_53);
    result_by_q.set(c_Key1(-8), common_term_53);
    // q = -7
    result_by_lpq.set(c_Key3(5, 4, -7), common_term_52);
    result_by_q.set(c_Key1(-7), common_term_52);
    // q = -6
    result_by_lpq.set(c_Key3(5, 4, -6), common_term_51);
    result_by_q.set(c_Key1(-6), common_term_51);
    // q = -5
    result_by_lpq.set(c_Key3(5, 4, -5), common_term_50);
    result_by_q.set(c_Key1(-5), common_term_50);
    // q = -4
    result_by_lpq.set(c_Key3(5, 4, -4), common_term_49);
    result_by_q.set(c_Key1(-4), common_term_49);
    // q = -3
    result_by_lpq.set(c_Key3(5, 4, -3), common_term_48);
    result_by_q.set(c_Key1(-3), common_term_48);
    // q = -2
    result_by_lpq.set(c_Key3(5, 4, -2), common_term_47);
    result_by_q.set(c_Key1(-2), common_term_47);
    // q = -1
    result_by_lpq.set(c_Key3(5, 4, -1), common_term_46);
    result_by_q.set(c_Key1(-1), common_term_46);
    // q = 0
    result_by_lpq.set(c_Key3(5, 4, 0), common_term_45);
    result_by_q.set(c_Key1(0), common_term_45);
    // q = 1
    result_by_lpq.set(c_Key3(5, 4, 1), common_term_44);
    result_by_q.set(c_Key1(1), common_term_44);
    // q = 2
    result_by_lpq.set(c_Key3(5, 4, 2), common_term_43);
    result_by_q.set(c_Key1(2), common_term_43);
    // q = 3
    result_by_lpq.set(c_Key3(5, 4, 3), common_term_42);
    result_by_q.set(c_Key1(3), common_term_42);
    // q = 4
    result_by_lpq.set(c_Key3(5, 4, 4), common_term_41);
    result_by_q.set(c_Key1(4), common_term_41);
    // q = 5
    result_by_lpq.set(c_Key3(5, 4, 5), common_term_40);
    result_by_q.set(c_Key1(5), common_term_40);
    // q = 6
    result_by_lpq.set(c_Key3(5, 4, 6), common_term_39);
    result_by_q.set(c_Key1(6), common_term_39);
    // q = 7
    result_by_lpq.set(c_Key3(5, 4, 7), common_term_38);
    result_by_q.set(c_Key1(7), common_term_38);
    // q = 8
    result_by_lpq.set(c_Key3(5, 4, 8), common_term_37);
    result_by_q.set(c_Key1(8), common_term_37);
    // q = 9
    result_by_lpq.set(c_Key3(5, 4, 9), common_term_36);
    result_by_q.set(c_Key1(9), common_term_36);
    // q = 10
    result_by_lpq.set(c_Key3(5, 4, 10), common_term_35);
    result_by_q.set(c_Key1(10), common_term_35);
    // q = 11
    result_by_lpq.set(c_Key3(5, 4, 11), common_term_34);
    result_by_q.set(c_Key1(11), common_term_34);
    // q = 12
    result_by_lpq.set(c_Key3(5, 4, 12), common_term_33);
    result_by_q.set(c_Key1(12), common_term_33);
    // q = 13
    result_by_lpq.set(c_Key3(5, 4, 13), common_term_32);
    result_by_q.set(c_Key1(13), common_term_32);
    // q = 14
    result_by_lpq.set(c_Key3(5, 4, 14), common_term_31);
    result_by_q.set(c_Key1(14), common_term_31);
    // q = 15
    result_by_lpq.set(c_Key3(5, 4, 15), common_term_30);
    result_by_q.set(c_Key1(15), common_term_30);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 4), result_by_q);
    result_by_q.clear();

    // l , p = (5, 5).
    // q = -15
    result_by_lpq.set(c_Key3(5, 5, -15), common_term_29);
    result_by_q.set(c_Key1(-15), common_term_29);
    // q = -14
    result_by_lpq.set(c_Key3(5, 5, -14), common_term_28);
    result_by_q.set(c_Key1(-14), common_term_28);
    // q = -13
    result_by_lpq.set(c_Key3(5, 5, -13), common_term_27);
    result_by_q.set(c_Key1(-13), common_term_27);
    // q = -12
    result_by_lpq.set(c_Key3(5, 5, -12), common_term_26);
    result_by_q.set(c_Key1(-12), common_term_26);
    // q = -11
    result_by_lpq.set(c_Key3(5, 5, -11), common_term_25);
    result_by_q.set(c_Key1(-11), common_term_25);
    // q = -10
    result_by_lpq.set(c_Key3(5, 5, -10), common_term_24);
    result_by_q.set(c_Key1(-10), common_term_24);
    // q = -9
    result_by_lpq.set(c_Key3(5, 5, -9), common_term_23);
    result_by_q.set(c_Key1(-9), common_term_23);
    // q = -8
    result_by_lpq.set(c_Key3(5, 5, -8), common_term_22);
    result_by_q.set(c_Key1(-8), common_term_22);
    // q = -7
    result_by_lpq.set(c_Key3(5, 5, -7), common_term_21);
    result_by_q.set(c_Key1(-7), common_term_21);
    // q = -6
    result_by_lpq.set(c_Key3(5, 5, -6), common_term_20);
    result_by_q.set(c_Key1(-6), common_term_20);
    // q = -5
    result_by_lpq.set(c_Key3(5, 5, -5), common_term_19);
    result_by_q.set(c_Key1(-5), common_term_19);
    // q = -4
    result_by_lpq.set(c_Key3(5, 5, -4), common_term_18);
    result_by_q.set(c_Key1(-4), common_term_18);
    // q = -3
    result_by_lpq.set(c_Key3(5, 5, -3), common_term_17);
    result_by_q.set(c_Key1(-3), common_term_17);
    // q = -2
    result_by_lpq.set(c_Key3(5, 5, -2), common_term_16);
    result_by_q.set(c_Key1(-2), common_term_16);
    // q = -1
    result_by_lpq.set(c_Key3(5, 5, -1), common_term_15);
    result_by_q.set(c_Key1(-1), common_term_15);
    // q = 0
    result_by_lpq.set(c_Key3(5, 5, 0), common_term_14);
    result_by_q.set(c_Key1(0), common_term_14);
    // q = 1
    result_by_lpq.set(c_Key3(5, 5, 1), common_term_13);
    result_by_q.set(c_Key1(1), common_term_13);
    // q = 2
    result_by_lpq.set(c_Key3(5, 5, 2), common_term_12);
    result_by_q.set(c_Key1(2), common_term_12);
    // q = 3
    result_by_lpq.set(c_Key3(5, 5, 3), common_term_11);
    result_by_q.set(c_Key1(3), common_term_11);
    // q = 4
    result_by_lpq.set(c_Key3(5, 5, 4), common_term_10);
    result_by_q.set(c_Key1(4), common_term_10);
    // q = 6
    result_by_lpq.set(c_Key3(5, 5, 6), common_term_9);
    result_by_q.set(c_Key1(6), common_term_9);
    // q = 7
    result_by_lpq.set(c_Key3(5, 5, 7), common_term_8);
    result_by_q.set(c_Key1(7), common_term_8);
    // q = 8
    result_by_lpq.set(c_Key3(5, 5, 8), common_term_7);
    result_by_q.set(c_Key1(8), common_term_7);
    // q = 9
    result_by_lpq.set(c_Key3(5, 5, 9), common_term_6);
    result_by_q.set(c_Key1(9), common_term_6);
    // q = 10
    result_by_lpq.set(c_Key3(5, 5, 10), common_term_5);
    result_by_q.set(c_Key1(10), common_term_5);
    // q = 11
    result_by_lpq.set(c_Key3(5, 5, 11), common_term_4);
    result_by_q.set(c_Key1(11), common_term_4);
    // q = 12
    result_by_lpq.set(c_Key3(5, 5, 12), common_term_3);
    result_by_q.set(c_Key1(12), common_term_3);
    // q = 13
    result_by_lpq.set(c_Key3(5, 5, 13), common_term_2);
    result_by_q.set(c_Key1(13), common_term_2);
    // q = 14
    result_by_lpq.set(c_Key3(5, 5, 14), common_term_1);
    result_by_q.set(c_Key1(14), common_term_1);
    // q = 15
    result_by_lpq.set(c_Key3(5, 5, 15), common_term_0);
    result_by_q.set(c_Key1(15), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 5), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l5_e20(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 5.
    // Functions are taylor expanded and truncated at eccentricity^20.

    // Functions Calculated for l = 5.

    c_IntMap<c_Key3, double> result_by_lpq(244);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(6);
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
    double common_term_0 = 0.13034688319858333*eccentricity_20;
    double common_term_1 = 0.093706159533454832*eccentricity_19;
    double common_term_2 = 0.04937121715117802*eccentricity_20 + 0.06700379470517017*eccentricity_18;
    double common_term_3 = 0.043622864530427556*eccentricity_19 + 0.047588579487739152*eccentricity_17;
    double common_term_4 = 0.053219446655162202*eccentricity_20 + 0.036467521288496459*eccentricity_18 + 0.033510695238077827*eccentricity_16;
    double common_term_5 = 0.040968958708190496*eccentricity_19 + 0.029171614577559725*eccentricity_17 + 0.02333729166204778*eccentricity_15;
    double common_term_6 = 0.036675216782145496*eccentricity_20 + 0.03088801522576805*eccentricity_18 + 0.022423059868169474*eccentricity_16 + 0.016016471334406767*eccentricity_14;
    double common_term_7 = 0.026926032880177501*eccentricity_19 + 0.022606128341578077*eccentricity_17 + 0.01655045718537782*eccentricity_15 + 0.010777041888152999*eccentricity_13;
    double common_term_8 = 0.020945442051497058*eccentricity_20 + 0.018896333344165057*eccentricity_18 + 0.015881591431506881*eccentricity_16 + 0.011667415808731163*eccentricity_14 + 0.007054716535511866*eccentricity_12;
    double common_term_9 = 0.013715940849747351*eccentricity_19 + 0.012474752159113208*eccentricity_17 + 0.01054002637987013*eccentricity_15 + 0.0077663352272727273*eccentricity_13 + 0.0044379058441558442*eccentricity_11;
    double common_term_10 = 0.0085268666325783791*eccentricity_20 + 0.0082019264262138981*eccentricity_18 + 0.0075470727315588646*eccentricity_16 + 0.0064426308455876757*eccentricity_14 + 0.0047783104678042829*eccentricity_12 + 0.0026280707572923556*eccentricity_10;
    double common_term_11 = 0.0043769150516835702*eccentricity_19 + 0.0042725184623332771*eccentricity_17 + 0.0039962677925640889*eccentricity_15 + 0.0034696168029501363*eccentricity_13 + 0.0026102292768959436*eccentricity_11 + 0.0014109347442680776*eccentricity_9;
    double common_term_12 = 0.0017455462329204279*eccentricity_20 + 0.0017771076895333789*eccentricity_18 + 0.0017683520226122497*eccentricity_16 + 0.0016919335451993075*eccentricity_14 + 0.001507652827671596*eccentricity_12 + 0.0011653355189732143*eccentricity_10 + 0.00063563755580357143*eccentricity_8;
    double common_term_13 = 0.00044215828704911701*eccentricity_19 + 0.00045936335772851977*eccentricity_17 + 0.0004686555326746299*eccentricity_15 + 0.00046257256760728983*eccentricity_13 + 0.0004285163139329806*eccentricity_11 + 0.00034722222222222222*eccentricity_9 + 0.00019841269841269841*eccentricity_7;
    double common_term_14 = 3.2005118802176957e-5*eccentricity_20 + 3.4277782785809149e-5*eccentricity_18 + 3.6525413383621298e-5*eccentricity_16 + 3.8490304520801542e-5*eccentricity_14 + 3.9638763732730747e-5*eccentricity_12 + 3.8885691809275794e-5*eccentricity_10 + 3.410218253968254e-5*eccentricity_8 + 2.1701388888888889e-5*eccentricity_6;
    double common_term_15 = 0.00062178430468765591*eccentricity_20 + 0.0006984841081645443*eccentricity_18 + 0.00079210976217042085*eccentricity_16 + 0.00090825499400610366*eccentricity_14 + 0.0010550213869286593*eccentricity_12 + 0.0012444067253637566*eccentricity_10 + 0.0014946831597222222*eccentricity_8 + 0.0018229166666666667*eccentricity_6 + 0.0026041666666666667*eccentricity_4;
    double common_term_16 = 0.0011968105893641711*eccentricity_19 + 0.0012464319559485317*eccentricity_17 + 0.0012744934811630413*eccentricity_15 + 0.0012573371362433862*eccentricity_13 + 0.0011388062169312169*eccentricity_11 + 0.001099537037037037*eccentricity_9 - 0.0041666666666666667*eccentricity_7 + 0.041666666666666667*eccentricity_5 - 0.16666666666666667*eccentricity_3;
    double common_term_17 = -0.001900914548669902*eccentricity_20 - 0.0026364412014087548*eccentricity_18 - 0.0036793709044553796*eccentricity_16 - 0.0051266195092882429*eccentricity_14 - 0.008580322265625*eccentricity_12 + 0.00580902099609375*eccentricity_10 - 0.1740234375*eccentricity_8 + 0.8525390625*eccentricity_6 - 2.25*eccentricity_4 + 1.125*eccentricity_2;
    double common_term_18 = -0.015330877138793653*eccentricity_19 - 0.019292902221491735*eccentricity_17 - 0.021366011392983119*eccentricity_15 - 0.068361992945326279*eccentricity_13 + 0.30132233796296296*eccentricity_11 - 2.2220486111111111*eccentricity_9 + 8.3368055555555556*eccentricity_7 - 16.833333333333333*eccentricity_5 + 11.5*eccentricity_3 - 2.0*eccentricity;
    double common_term_19 = -0.040859207705474203*eccentricity_20 - 0.055092025349983133*eccentricity_18 - 8.7820967524665584e-5*eccentricity_16 - 0.63502111012042395*eccentricity_14 + 3.6683625939451618*eccentricity_12 - 17.866980658637153*eccentricity_10 + 53.395758734809028*eccentricity_8 - 90.342881944444444*eccentricity_6 + 66.484375*eccentricity_4 - 17.5*eccentricity_2 + 1.0;
    double common_term_20 = -0.18488623545121173*eccentricity_19 + 0.63875543088329082*eccentricity_17 - 5.4744128667091837*eccentricity_15 + 27.873044084821429*eccentricity_13 - 105.696796875*eccentricity_11 + 262.5140625*eccentricity_9 - 388.625*eccentricity_7 + 286.5*eccentricity_5 - 88.5*eccentricity_3 + 8.0*eccentricity;
    double common_term_21 = -1.0884223437412307*eccentricity_20 + 6.3668344487724724*eccentricity_18 - 37.494931585554887*eccentricity_16 + 160.33849005591722*eccentricity_14 - 503.92212854456019*eccentricity_12 + 1075.0551734641746*eccentricity_10 - 1426.3118272569444*eccentricity_8 + 1024.1761067708333*eccentricity_6 - 338.91666666666667*eccentricity_4 + 37.375*eccentricity_2;
    double common_term_22 = 43.598649818944057*eccentricity_19 - 208.63642635393921*eccentricity_17 + 759.82512425289046*eccentricity_15 - 2049.3854050925926*eccentricity_13 + 3847.0054274140212*eccentricity_11 - 4644.6884259259259*eccentricity_9 + 3210.0166666666667*eccentricity_7 - 1090.1666666666667*eccentricity_5 + 133.16666666666667*eccentricity_3;
    double common_term_23 = 241.11708148646427*eccentricity_20 - 984.23438598949129*eccentricity_18 + 3114.0833738291221*eccentricity_16 - 7378.0934716027124*eccentricity_14 + 12404.238924591882*eccentricity_12 - 13776.578367396763*eccentricity_10 + 9117.8461669921875*eccentricity_8 - 3105.83671875*eccentricity_6 + 400.6796875*eccentricity_4;
    double common_term_24 = -4073.5442484584606*eccentricity_19 + 11393.061790127197*eccentricity_17 - 24118.051818381191*eccentricity_15 + 36805.308178380824*eccentricity_13 - 37910.616887814153*eccentricity_11 + 23983.827504960317*eccentricity_9 - 8086.7361111111111*eccentricity_7 + 1072.4583333333333*eccentricity_5;
    double common_term_25 = -15170.418590552063*eccentricity_20 + 38038.879129569927*eccentricity_18 - 72891.23457818038*eccentricity_16 + 102024.90101705193*eccentricity_14 - 98093.605275043216*eccentricity_12 + 59306.949686698308*eccentricity_10 - 19637.63640562996*eccentricity_8 + 2633.1491102430556*eccentricity_6;
    double common_term_26 = 117780.75740898134*eccentricity_19 - 206417.08788714045*eccentricity_17 + 267204.62059265929*eccentricity_15 - 241078.85928571429*eccentricity_13 + 139375.73253348214*eccentricity_11 - 45101.5265625*eccentricity_9 + 6049.0116071428571*eccentricity_7;
    double common_term_27 = 342324.7307661078*eccentricity_20 - 553303.6081348789*eccentricity_18 + 666911.92329358198*eccentricity_16 - 567140.76834504931*eccentricity_14 + 313848.34123679371*eccentricity_12 - 98957.484222180671*eccentricity_10 + 13179.87816617451*eccentricity_8;
    double common_term_28 = -1415021.1889027011*eccentricity_19 + 1597044.5660111935*eccentricity_17 - 1284970.2308093875*eccentricity_15 + 681474.58250646658*eccentricity_13 - 208988.53118537809*eccentricity_11 + 27504.092611882716*eccentricity_9;
    double common_term_29 = -3474344.0864772061*eccentricity_20 + 3689255.8650042882*eccentricity_18 - 2817708.9136906191*eccentricity_16 + 1433984.7381258379*eccentricity_14 - 427283.52163934088*eccentricity_12 + 55372.913824353899*eccentricity_10;
    double common_term_30 = 8257405.3851045284*eccentricity_19 - 6003920.8919581193*eccentricity_17 + 2935945.1361715278*eccentricity_15 - 849564.73169093381*eccentricity_13 + 108149.79697137546*eccentricity_11;
    double common_term_31 = 17972412.482659033*eccentricity_20 - 12472160.799540178*eccentricity_18 + 5867921.5062678656*eccentricity_16 - 1648688.3416402937*eccentricity_14 + 205812.7946639314*eccentricity_12;
    double common_term_32 = -25328686.237587594*eccentricity_19 + 11479769.985796126*eccentricity_17 - 3132020.7740035927*eccentricity_15 + 382953.25728189779*eccentricity_13;
    double common_term_33 = -50403156.156724156*eccentricity_20 + 22033475.005640188*eccentricity_18 - 5838659.8916396411*eccentricity_16 + 698665.82920653284*eccentricity_14;
    double common_term_34 = 41569091.03835248*eccentricity_19 - 10702585.837331192*eccentricity_17 + 1252708.9143546075*eccentricity_15;
    double common_term_35 = 77216615.508264751*eccentricity_20 - 19324057.864394923*eccentricity_18 + 2211698.9437125349*eccentricity_16;
    double common_term_36 = -34417472.294430327*eccentricity_19 + 3851258.1027234425*eccentricity_17;
    double common_term_37 = -60545005.512002398*eccentricity_20 + 6623379.283608481*eccentricity_18;
    double common_term_38 = 11263436.460914544*eccentricity_19;
    double common_term_39 = 18959213.164720992*eccentricity_20;
    double common_term_40 = 101.30344153893582*eccentricity_20;
    double common_term_41 = 74.022405576633975*eccentricity_19;
    double common_term_42 = 32.742752984305674*eccentricity_20 + 54.084754634960899*eccentricity_18;
    double common_term_43 = 34.173913744413025*eccentricity_19 + 39.51516483781995*eccentricity_17;
    double common_term_44 = 72.289325805799827*eccentricity_20 + 32.468413469162155*eccentricity_18 + 28.869294706593075*eccentricity_16;
    double common_term_45 = 60.382442904628525*eccentricity_19 + 29.209707000379643*eccentricity_17 + 21.091149702836226*eccentricity_15;
    double common_term_46 = 82.203082658031135*eccentricity_20 + 50.69299762308267*eccentricity_18 + 25.356777714868295*eccentricity_16 + 15.408737763597732*eccentricity_14;
    double common_term_47 = 71.160280721757899*eccentricity_19 + 42.61207083884018*eccentricity_17 + 21.465307221747884*eccentricity_15 + 11.257742673270017*eccentricity_13;
    double common_term_48 = 96.92416349492058*eccentricity_20 + 61.312850032383186*eccentricity_18 + 35.772356813847487*eccentricity_16 + 17.835125904545322*eccentricity_14 + 8.2256967558179583*eccentricity_12;
    double common_term_49 = 84.530214718504035*eccentricity_19 + 52.602874519198536*eccentricity_17 + 29.942290013347956*eccentricity_15 + 14.607572256543611*eccentricity_13 + 6.0111507560726311*eccentricity_11;
    double common_term_50 = 112.28114248699222*eccentricity_20 + 73.483653229148223*eccentricity_18 + 44.946439727013113*eccentricity_16 + 24.964970074140129*eccentricity_14 + 11.829020482329809*eccentricity_12 + 4.3938146464029948*eccentricity_10;
    double common_term_51 = 98.653415154930867*eccentricity_19 + 63.671757086907233*eccentricity_17 + 38.251415635146104*eccentricity_15 + 20.724025720373377*eccentricity_13 + 9.49171875*eccentricity_11 + 3.2127232142857143*eccentricity_9;
    double common_term_52 = 128.48301614654084*eccentricity_20 + 86.434658455066511*eccentricity_18 + 54.988036567405389*eccentricity_16 + 32.425630900461131*eccentricity_14 + 17.125299496625466*eccentricity_12 + 7.5595372921575314*eccentricity_10 + 2.3502705891927083*eccentricity_8;
    double common_term_53 = 113.59521277505208*eccentricity_19 + 75.512989852513019*eccentricity_17 + 47.331421073802259*eccentricity_15 + 27.380366236772487*eccentricity_13 + 14.087541335978836*eccentricity_11 + 5.9837797619047619*eccentricity_9 + 1.7205357142857143*eccentricity_7;
    double common_term_54 = 145.51638409014868*eccentricity_20 + 100.17931976285923*eccentricity_18 + 65.781423915345947*eccentricity_16 + 40.606042443854468*eccentricity_14 + 23.031711687360491*eccentricity_12 + 11.538018035888672*eccentricity_10 + 4.7125279017857143*eccentricity_8 + 1.2607421875*eccentricity_6;
    double common_term_55 = 129.34449168686201*eccentricity_19 + 88.12254088346208*eccentricity_17 + 57.137932308889991*eccentricity_15 + 34.721210868606702*eccentricity_13 + 19.301029265873016*eccentricity_11 + 9.410639880952381*eccentricity_9 + 3.6958333333333333*eccentricity_7 + 0.925*eccentricity_5;
    double common_term_56 = 163.37101242103358*eccentricity_20 + 114.70668577481098*eccentricity_18 + 77.316947758029956*eccentricity_16 + 49.485521528195242*eccentricity_14 + 29.591459706472972*eccentricity_12 + 16.115094284784226*eccentricity_10 + 7.6452555338541667*eccentricity_8 + 2.88828125*eccentricity_6 + 0.6796875*eccentricity_4;
    double common_term_57 = 0.5*eccentricity_3*std::pow(1.0 - eccentricity_2, -4.5);
    double common_term_58 = 182.0392914740479*eccentricity_20 + 130.00835943664295*eccentricity_18 + 89.585365880957436*eccentricity_16 + 59.054041014140245*eccentricity_14 + 36.793448970734127*eccentricity_12 + 21.283941819932726*eccentricity_10 + 11.114388020833333*eccentricity_8 + 4.9912109375*eccentricity_6 + 1.75*eccentricity_4 + 0.375*eccentricity_2;
    double common_term_59 = 163.23127634061437*eccentricity_19 + 115.60951004860012*eccentricity_17 + 78.891815424520503*eccentricity_15 + 51.408773251488095*eccentricity_13 + 31.589105902777778*eccentricity_11 + 17.966145833333333*eccentricity_9 + 9.1875*eccentricity_7 + 4.0*eccentricity_5 + 1.5*eccentricity_3;
    double common_term_60 = 201.5200179409498*eccentricity_20 + 146.08308230954051*eccentricity_18 + 102.58540450008016*eccentricity_16 + 69.31039448524008*eccentricity_14 + 44.636981925964355*eccentricity_12 + 27.040718994140625*eccentricity_10 + 15.16680908203125*eccentricity_8 + 7.2265625*eccentricity_6 + 4.734375*eccentricity_4 - 1.5*eccentricity_2 + 1.0;
    double common_term_61 = 181.3621142091418*eccentricity_19 + 130.47979079114767*eccentricity_17 + 90.831207124994095*eccentricity_15 + 60.754408895502646*eccentricity_13 + 38.584227430555556*eccentricity_11 + 23.668229166666667*eccentricity_9 + 9.6770833333333333*eccentricity_7 + 15.5*eccentricity_5 - 10.5*eccentricity_3 + 6.0*eccentricity;
    double common_term_62 = 221.80767662263434*eccentricity_20 + 162.92584278025849*eccentricity_18 + 116.30249946367143*eccentricity_16 + 80.343051059693886*eccentricity_14 + 52.36621335953001*eccentricity_12 + 37.758710225423177*eccentricity_10 + 2.5192057291666667*eccentricity_8 + 51.7041015625*eccentricity_6 - 44.25*eccentricity_4 + 22.125*eccentricity_2;
    double common_term_63 = 200.28753368868854*eccentricity_19 + 145.9886597526706*eccentricity_17 + 104.32126953125*eccentricity_15 + 65.505837053571429*eccentricity_13 + 70.158314732142857*eccentricity_11 - 45.9796875*eccentricity_9 + 164.025*eccentricity_7 - 145.125*eccentricity_5 + 64.5*eccentricity_3;
    double common_term_64 = 243.00403911477662*eccentricity_20 + 179.66055128726637*eccentricity_18 + 136.35942733855098*eccentricity_16 + 63.870867509900788*eccentricity_14 + 166.38394224908617*eccentricity_12 - 231.58345404730903*eccentricity_10 + 480.89112141927083*eccentricity_8 - 407.71484375*eccentricity_6 + 163.0859375*eccentricity_4;
    double common_term_65 = 214.45168786365633*eccentricity_19 + 191.97484058503145*eccentricity_17 - 6.6140546840553351*eccentricity_15 + 469.22716910548942*eccentricity_13 - 812.08832465277778*eccentricity_11 + 1302.5502976190476*eccentricity_9 - 1029.7375*eccentricity_7 + 374.45*eccentricity_5;
    double common_term_66 = 235.93868304345835*eccentricity_20 + 330.24015701988642*eccentricity_18 - 324.79488135284882*eccentricity_16 + 1382.9710076344013*eccentricity_14 - 2426.4057120186942*eccentricity_12 + 3289.3376407078334*eccentricity_10 - 2404.9318359375*eccentricity_8 + 801.6439453125*eccentricity_6;
    double common_term_67 = 752.35254982054338*eccentricity_19 - 1423.5265543557814*eccentricity_17 + 3964.3006371220889*eccentricity_15 - 6560.9110105268959*eccentricity_13 + 7825.6551235945767*eccentricity_11 - 5288.5431547619048*eccentricity_9 + 1627.2440476190476*eccentricity_7;
    double common_term_68 = 2082.1336630021181*eccentricity_20 - 4767.5118968235651*eccentricity_18 + 10804.857791881188*eccentricity_16 - 16504.128205581729*eccentricity_14 + 17704.166441436142*eccentricity_12 - 11086.362901814779*eccentricity_10 + 3167.5322576613653*eccentricity_8;
    double common_term_69 = -14131.050672884375*eccentricity_19 + 27940.360954393507*eccentricity_17 - 39239.608693752537*eccentricity_15 + 38382.229918831169*eccentricity_13 - 22351.217912946429*eccentricity_11 + 5960.3247767857143*eccentricity_9;
    double common_term_70 = -38736.735070583516*eccentricity_20 + 68836.881642861614*eccentricity_18 - 89099.900719965592*eccentricity_16 + 80243.948315284375*eccentricity_14 - 43623.675519422097*eccentricity_12 + 10905.918879855524*eccentricity_10;
    double common_term_71 = 162469.86779797639*eccentricity_19 - 194661.20638780322*eccentricity_17 + 162603.91326242611*eccentricity_15 - 82837.988055095223*eccentricity_13 + 19491.291307081229*eccentricity_11;
    double common_term_72 = 369285.58781710729*eccentricity_20 - 411499.41875095112*eccentricity_18 + 320689.23309976115*eccentricity_16 - 153647.68956916828*eccentricity_14 + 34143.931015370728*eccentricity_12;
    double common_term_73 = -845388.19377453576*eccentricity_19 + 617657.06957927141*eccentricity_17 - 279235.16068955808*eccentricity_15 + 58786.349618854332*eccentricity_13;
    double common_term_74 = -1693851.4245874523*eccentricity_20 + 1165040.1795688612*eccentricity_18 - 498497.9703457373*eccentricity_16 + 99699.59406914746*eccentricity_14;
    double common_term_75 = 2157170.5274973838*eccentricity_19 - 876016.29182781623*eccentricity_17 + 166860.24606244119*eccentricity_15;
    double common_term_76 = 3928595.1526072801*eccentricity_20 - 1517998.4109456487*eccentricity_18 + 275999.71108102704*eccentricity_16;
    double common_term_77 = -2597624.6277619236*eccentricity_19 + 451760.80482816062*eccentricity_17;
    double common_term_78 = -4395089.9943353128*eccentricity_20 + 732514.99905588547*eccentricity_18;
    double common_term_79 = 1177684.0382157204*eccentricity_19;
    double common_term_80 = 1878825.8558761366*eccentricity_20;
    double common_term_81 = 5564.332528881033*eccentricity_20;
    double common_term_82 = 3881.6156773759586*eccentricity_19;
    double common_term_83 = -2925.4678236390773*eccentricity_20 + 2701.4251077049553*eccentricity_18;
    double common_term_84 = -1512.3205334539453*eccentricity_19 + 1875.2342746080607*eccentricity_17;
    double common_term_85 = 2111.4442718734672*eccentricity_20 - 684.89005417649475*eccentricity_18 + 1298.0341112390616*eccentricity_16;
    double common_term_86 = 1408.0515163422989*eccentricity_19 - 220.39207753001276*eccentricity_17 + 895.6770696572522*eccentricity_15;
    double common_term_87 = 874.21169953947969*eccentricity_20 + 984.59027693740732*eccentricity_18 + 23.841435871712988*eccentricity_16 + 615.8782870739856*eccentricity_14;
    double common_term_88 = 814.83802406828662*eccentricity_19 + 721.17261793751728*eccentricity_17 + 138.04397284286584*eccentricity_15 + 421.82299405672453*eccentricity_13;
    double common_term_89 = 1101.3108660887688*eccentricity_20 + 723.43140748559176*eccentricity_18 + 549.07101660530937*eccentricity_16 + 178.41581534061428*eccentricity_14 + 287.62977158715766*eccentricity_12;
    double common_term_90 = 944.27559070679386*eccentricity_19 + 624.77040083163669*eccentricity_17 + 429.61432180327308*eccentricity_15 + 179.28346223165017*eccentricity_13 + 195.1342063366803*eccentricity_11;
    double common_term_91 = 1162.2144840169674*eccentricity_20 + 807.16550083085718*eccentricity_18 + 530.17529606347437*eccentricity_16 + 341.44604407120061*eccentricity_14 + 161.079640740729*eccentricity_12 + 131.61114809308733*eccentricity_10;
    double common_term_92 = 1004.5034637500127*eccentricity_19 + 686.51767710309666*eccentricity_17 + 444.23567553653961*eccentricity_15 + 272.95488521123938*eccentricity_13 + 135.5390270888448*eccentricity_11 + 88.163544422398589*eccentricity_9;
    double common_term_93 = 1230.4126539300074*eccentricity_20 + 863.42159013784042*eccentricity_18 + 580.31405412056381*eccentricity_16 + 368.32778909907642*eccentricity_14 + 217.86436910017037*eccentricity_12 + 109.04968559359327*eccentricity_10 + 58.584195285373264*eccentricity_8;
    double common_term_94 = 1066.5684615169986*eccentricity_19 + 737.81344454837155*eccentricity_17 + 487.12988477830763*eccentricity_15 + 302.38700613839286*eccentricity_13 + 172.72349330357143*eccentricity_11 + 84.790848214285714*eccentricity_9 + 38.552678571428571*eccentricity_7;
    double common_term_95 = 1297.7147279879589*eccentricity_20 + 919.70293081591935*eccentricity_18 + 626.51567608094325*eccentricity_16 + 405.78072080117755*eccentricity_14 + 245.76078824195211*eccentricity_12 + 135.52102582416837*eccentricity_10 + 64.078047495039683*eccentricity_8 + 25.069552951388889*eccentricity_6;
    double common_term_96 = 1127.8603763960602*eccentricity_19 + 788.59849767313756*eccentricity_17 + 528.39256430213385*eccentricity_15 + 335.18820862498163*eccentricity_13 + 197.59574239417989*eccentricity_11 + 104.94945436507937*eccentricity_9 + 47.195138888888889*eccentricity_7 + 16.058333333333333*eccentricity_5;
    double common_term_97 = 1364.2471347594452*eccentricity_20 + 975.25224154950479*eccentricity_18 + 672.07692575220977*eccentricity_16 + 442.34680258410318*eccentricity_14 + 274.3350358077458*eccentricity_12 + 157.00033656529018*eccentricity_10 + 80.0341552734375*eccentricity_8 + 33.89765625*eccentricity_6 + 10.0859375*eccentricity_4;
    double common_term_98 = 1188.3994728119306*eccentricity_19 + 838.6738190358074*eccentricity_17 + 569.00009836336316*eccentricity_15 + 367.3229865244709*eccentricity_13 + 222.25277364417989*eccentricity_11 + 123.10688657407407*eccentricity_9 + 59.960416666666667*eccentricity_7 + 23.708333333333333*eccentricity_5 + 6.1666666666666667*eccentricity_3;
    double common_term_99 = 1430.002760646794*eccentricity_20 + 1030.0666697113358*eccentricity_18 + 716.94775855110242*eccentricity_16 + 478.27117265849938*eccentricity_14 + 302.30942661184482*eccentricity_12 + 178.02036545364945*eccentricity_10 + 95.092730034722222*eccentricity_8 + 44.002278645833333*eccentricity_6 + 16.083333333333333*eccentricity_4 + 3.625*eccentricity_2;
    double common_term_100 = std::pow(1.0 - eccentricity_2, -4.5)*(1.5*eccentricity_3 + 2.0*eccentricity);
    double common_term_101 = 1494.977913219978*eccentricity_20 + 1084.141845620006*eccentricity_18 + 761.12296438316618*eccentricity_16 + 513.54641791282898*eccentricity_14 + 329.68472556549826*eccentricity_12 + 198.49574069552951*eccentricity_10 + 109.66821967230903*eccentricity_8 + 53.677951388888889*eccentricity_6 + 21.859375*eccentricity_4 + 6.5*eccentricity_2 + 1.0;
    double common_term_102 = 1307.1999153318526*eccentricity_19 + 936.67434108648805*eccentricity_17 + 648.20072550424973*eccentricity_15 + 429.72375657930996*eccentricity_13 + 269.85374421296296*eccentricity_11 + 157.90815972222222*eccentricity_9 + 83.961805555555556*eccentricity_7 + 38.916666666666667*eccentricity_5 + 14.5*eccentricity_3 + 4.0*eccentricity;
    double common_term_103 = 1559.1723844842473*eccentricity_20 + 1137.4775636000659*eccentricity_18 + 804.60234598181686*eccentricity_16 + 548.17235871945109*eccentricity_14 + 356.46074497767857*eccentricity_12 + 218.42696685791016*eccentricity_10 + 123.7505859375*eccentricity_8 + 63.0263671875*eccentricity_6 + 26.75*eccentricity_4 + 10.875*eccentricity_2;
    double common_term_104 = 1365.4592879608301*eccentricity_19 + 984.5973290939487*eccentricity_17 + 686.79135353418792*eccentricity_15 + 459.98591817542989*eccentricity_13 + 292.80587590939153*eccentricity_11 + 174.4318287037037*eccentricity_9 + 95.977083333333333*eccentricity_7 + 43.041666666666667*eccentricity_5 + 24.833333333333333*eccentricity_3;
    double common_term_105 = 1622.5853473876708*eccentricity_20 + 1190.0729790794319*eccentricity_18 + 847.386136116445*eccentricity_16 + 582.13725769954314*eccentricity_14 + 382.74131503054705*eccentricity_12 + 237.09027357959243*eccentricity_10 + 140.88095431857639*eccentricity_8 + 60.779947916666667*eccentricity_6 + 51.221354166666667*eccentricity_4;
    double common_term_106 = 1422.956059602925*eccentricity_19 + 1031.8115933473265*eccentricity_17 + 724.61253892299107*eccentricity_15 + 490.30125558035714*eccentricity_13 + 311.54815848214286*eccentricity_11 + 204.35424107142857*eccentricity_9 + 72.31875*eccentricity_7 + 98.775*eccentricity_5;
    double common_term_107 = 1685.2065946877134*eccentricity_20 + 1242.0092952469509*eccentricity_18 + 888.87784927769694*eccentricity_16 + 618.84262999488003*eccentricity_14 + 393.72408189812805*eccentricity_12 + 301.28983722262912*eccentricity_10 + 60.194509548611111*eccentricity_8 + 181.54624565972222*eccentricity_6;
    double common_term_108 = 1480.1862066148248*eccentricity_19 + 1075.3451519848438*eccentricity_17 + 775.77064236998373*eccentricity_15 + 469.16829952968842*eccentricity_13 + 464.57999820877425*eccentricity_11 - 11.144667658730159*eccentricity_9 + 321.84662698412698*eccentricity_7;
    double common_term_109 = 1749.5097877179296*eccentricity_20 + 1280.7425386831805*eccentricity_18 + 979.91205114165394*eccentricity_16 + 498.47863376902295*eccentricity_14 + 763.95394112450736*eccentricity_12 - 209.35917576381138*eccentricity_10 + 554.7166024344308*eccentricity_8;
    double common_term_110 = 1490.5958167494981*eccentricity_19 + 1279.5942505011368*eccentricity_17 + 386.31665203684428*eccentricity_15 + 1341.1313508660514*eccentricity_13 - 657.54596905313051*eccentricity_11 + 934.67407958553792*eccentricity_9;
    double common_term_111 = 1659.4989809890515*eccentricity_20 + 1794.1679169895196*eccentricity_18 - 81.449406696395013*eccentricity_16 + 2473.528428908686*eccentricity_14 - 1572.4675775587832*eccentricity_12 + 1545.8537958783616*eccentricity_10;
    double common_term_112 = 2804.5875311177122*eccentricity_19 - 1360.674104664518*eccentricity_17 + 4687.0080340118866*eccentricity_15 - 3326.4069472909903*eccentricity_13 + 2517.1731057224026*eccentricity_11;
    double common_term_113 = 4942.9542095989066*eccentricity_20 - 4382.2684311336832*eccentricity_18 + 8951.5619111832492*eccentricity_16 - 6546.1679942372922*eccentricity_14 + 4044.9237368728242*eccentricity_12;
    double common_term_114 = -10984.204108158708*eccentricity_19 + 17015.360916425528*eccentricity_17 - 12270.047500231993*eccentricity_15 + 6426.300886272694*eccentricity_13;
    double common_term_115 = -24694.653433418976*eccentricity_20 + 31966.995547568329*eccentricity_18 - 22194.685671492962*eccentricity_16 + 10108.995635604223*eccentricity_14;
    double common_term_116 = 59170.215666304854*eccentricity_19 - 39060.443628485108*eccentricity_17 + 15764.316954849364*eccentricity_15;
    double common_term_117 = 107801.26347521513*eccentricity_20 - 67249.151903960326*eccentricity_18 + 24394.695545820952*eccentricity_16;
    double common_term_118 = -113705.86390192126*eccentricity_19 + 37491.312354412894*eccentricity_17;
    double common_term_119 = -189352.74605325478*eccentricity_20 + 57264.646977309608*eccentricity_18;
    double common_term_120 = 86980.895560501239*eccentricity_19;
    double common_term_121 = 131451.80748928829*eccentricity_20;

    c_IntMap<c_Key1, double> result_by_q(41);

    // eccentricity function by mode:
    // l , p = (5, 0).
    // q = -20
    result_by_lpq.set(c_Key3(5, 0, -20), common_term_0);
    result_by_q.set(c_Key1(-20), common_term_0);
    // q = -19
    result_by_lpq.set(c_Key3(5, 0, -19), common_term_1);
    result_by_q.set(c_Key1(-19), common_term_1);
    // q = -18
    result_by_lpq.set(c_Key3(5, 0, -18), common_term_2);
    result_by_q.set(c_Key1(-18), common_term_2);
    // q = -17
    result_by_lpq.set(c_Key3(5, 0, -17), common_term_3);
    result_by_q.set(c_Key1(-17), common_term_3);
    // q = -16
    result_by_lpq.set(c_Key3(5, 0, -16), common_term_4);
    result_by_q.set(c_Key1(-16), common_term_4);
    // q = -15
    result_by_lpq.set(c_Key3(5, 0, -15), common_term_5);
    result_by_q.set(c_Key1(-15), common_term_5);
    // q = -14
    result_by_lpq.set(c_Key3(5, 0, -14), common_term_6);
    result_by_q.set(c_Key1(-14), common_term_6);
    // q = -13
    result_by_lpq.set(c_Key3(5, 0, -13), common_term_7);
    result_by_q.set(c_Key1(-13), common_term_7);
    // q = -12
    result_by_lpq.set(c_Key3(5, 0, -12), common_term_8);
    result_by_q.set(c_Key1(-12), common_term_8);
    // q = -11
    result_by_lpq.set(c_Key3(5, 0, -11), common_term_9);
    result_by_q.set(c_Key1(-11), common_term_9);
    // q = -10
    result_by_lpq.set(c_Key3(5, 0, -10), common_term_10);
    result_by_q.set(c_Key1(-10), common_term_10);
    // q = -9
    result_by_lpq.set(c_Key3(5, 0, -9), common_term_11);
    result_by_q.set(c_Key1(-9), common_term_11);
    // q = -8
    result_by_lpq.set(c_Key3(5, 0, -8), common_term_12);
    result_by_q.set(c_Key1(-8), common_term_12);
    // q = -7
    result_by_lpq.set(c_Key3(5, 0, -7), common_term_13);
    result_by_q.set(c_Key1(-7), common_term_13);
    // q = -6
    result_by_lpq.set(c_Key3(5, 0, -6), common_term_14);
    result_by_q.set(c_Key1(-6), common_term_14);
    // q = -4
    result_by_lpq.set(c_Key3(5, 0, -4), common_term_15);
    result_by_q.set(c_Key1(-4), common_term_15);
    // q = -3
    result_by_lpq.set(c_Key3(5, 0, -3), common_term_16);
    result_by_q.set(c_Key1(-3), common_term_16);
    // q = -2
    result_by_lpq.set(c_Key3(5, 0, -2), common_term_17);
    result_by_q.set(c_Key1(-2), common_term_17);
    // q = -1
    result_by_lpq.set(c_Key3(5, 0, -1), common_term_18);
    result_by_q.set(c_Key1(-1), common_term_18);
    // q = 0
    result_by_lpq.set(c_Key3(5, 0, 0), common_term_19);
    result_by_q.set(c_Key1(0), common_term_19);
    // q = 1
    result_by_lpq.set(c_Key3(5, 0, 1), common_term_20);
    result_by_q.set(c_Key1(1), common_term_20);
    // q = 2
    result_by_lpq.set(c_Key3(5, 0, 2), common_term_21);
    result_by_q.set(c_Key1(2), common_term_21);
    // q = 3
    result_by_lpq.set(c_Key3(5, 0, 3), common_term_22);
    result_by_q.set(c_Key1(3), common_term_22);
    // q = 4
    result_by_lpq.set(c_Key3(5, 0, 4), common_term_23);
    result_by_q.set(c_Key1(4), common_term_23);
    // q = 5
    result_by_lpq.set(c_Key3(5, 0, 5), common_term_24);
    result_by_q.set(c_Key1(5), common_term_24);
    // q = 6
    result_by_lpq.set(c_Key3(5, 0, 6), common_term_25);
    result_by_q.set(c_Key1(6), common_term_25);
    // q = 7
    result_by_lpq.set(c_Key3(5, 0, 7), common_term_26);
    result_by_q.set(c_Key1(7), common_term_26);
    // q = 8
    result_by_lpq.set(c_Key3(5, 0, 8), common_term_27);
    result_by_q.set(c_Key1(8), common_term_27);
    // q = 9
    result_by_lpq.set(c_Key3(5, 0, 9), common_term_28);
    result_by_q.set(c_Key1(9), common_term_28);
    // q = 10
    result_by_lpq.set(c_Key3(5, 0, 10), common_term_29);
    result_by_q.set(c_Key1(10), common_term_29);
    // q = 11
    result_by_lpq.set(c_Key3(5, 0, 11), common_term_30);
    result_by_q.set(c_Key1(11), common_term_30);
    // q = 12
    result_by_lpq.set(c_Key3(5, 0, 12), common_term_31);
    result_by_q.set(c_Key1(12), common_term_31);
    // q = 13
    result_by_lpq.set(c_Key3(5, 0, 13), common_term_32);
    result_by_q.set(c_Key1(13), common_term_32);
    // q = 14
    result_by_lpq.set(c_Key3(5, 0, 14), common_term_33);
    result_by_q.set(c_Key1(14), common_term_33);
    // q = 15
    result_by_lpq.set(c_Key3(5, 0, 15), common_term_34);
    result_by_q.set(c_Key1(15), common_term_34);
    // q = 16
    result_by_lpq.set(c_Key3(5, 0, 16), common_term_35);
    result_by_q.set(c_Key1(16), common_term_35);
    // q = 17
    result_by_lpq.set(c_Key3(5, 0, 17), common_term_36);
    result_by_q.set(c_Key1(17), common_term_36);
    // q = 18
    result_by_lpq.set(c_Key3(5, 0, 18), common_term_37);
    result_by_q.set(c_Key1(18), common_term_37);
    // q = 19
    result_by_lpq.set(c_Key3(5, 0, 19), common_term_38);
    result_by_q.set(c_Key1(19), common_term_38);
    // q = 20
    result_by_lpq.set(c_Key3(5, 0, 20), common_term_39);
    result_by_q.set(c_Key1(20), common_term_39);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 0), result_by_q);
    result_by_q.clear();

    // l , p = (5, 1).
    // q = -20
    result_by_lpq.set(c_Key3(5, 1, -20), common_term_40);
    result_by_q.set(c_Key1(-20), common_term_40);
    // q = -19
    result_by_lpq.set(c_Key3(5, 1, -19), common_term_41);
    result_by_q.set(c_Key1(-19), common_term_41);
    // q = -18
    result_by_lpq.set(c_Key3(5, 1, -18), common_term_42);
    result_by_q.set(c_Key1(-18), common_term_42);
    // q = -17
    result_by_lpq.set(c_Key3(5, 1, -17), common_term_43);
    result_by_q.set(c_Key1(-17), common_term_43);
    // q = -16
    result_by_lpq.set(c_Key3(5, 1, -16), common_term_44);
    result_by_q.set(c_Key1(-16), common_term_44);
    // q = -15
    result_by_lpq.set(c_Key3(5, 1, -15), common_term_45);
    result_by_q.set(c_Key1(-15), common_term_45);
    // q = -14
    result_by_lpq.set(c_Key3(5, 1, -14), common_term_46);
    result_by_q.set(c_Key1(-14), common_term_46);
    // q = -13
    result_by_lpq.set(c_Key3(5, 1, -13), common_term_47);
    result_by_q.set(c_Key1(-13), common_term_47);
    // q = -12
    result_by_lpq.set(c_Key3(5, 1, -12), common_term_48);
    result_by_q.set(c_Key1(-12), common_term_48);
    // q = -11
    result_by_lpq.set(c_Key3(5, 1, -11), common_term_49);
    result_by_q.set(c_Key1(-11), common_term_49);
    // q = -10
    result_by_lpq.set(c_Key3(5, 1, -10), common_term_50);
    result_by_q.set(c_Key1(-10), common_term_50);
    // q = -9
    result_by_lpq.set(c_Key3(5, 1, -9), common_term_51);
    result_by_q.set(c_Key1(-9), common_term_51);
    // q = -8
    result_by_lpq.set(c_Key3(5, 1, -8), common_term_52);
    result_by_q.set(c_Key1(-8), common_term_52);
    // q = -7
    result_by_lpq.set(c_Key3(5, 1, -7), common_term_53);
    result_by_q.set(c_Key1(-7), common_term_53);
    // q = -6
    result_by_lpq.set(c_Key3(5, 1, -6), common_term_54);
    result_by_q.set(c_Key1(-6), common_term_54);
    // q = -5
    result_by_lpq.set(c_Key3(5, 1, -5), common_term_55);
    result_by_q.set(c_Key1(-5), common_term_55);
    // q = -4
    result_by_lpq.set(c_Key3(5, 1, -4), common_term_56);
    result_by_q.set(c_Key1(-4), common_term_56);
    // q = -3
    result_by_lpq.set(c_Key3(5, 1, -3), common_term_57);
    result_by_q.set(c_Key1(-3), common_term_57);
    // q = -2
    result_by_lpq.set(c_Key3(5, 1, -2), common_term_58);
    result_by_q.set(c_Key1(-2), common_term_58);
    // q = -1
    result_by_lpq.set(c_Key3(5, 1, -1), common_term_59);
    result_by_q.set(c_Key1(-1), common_term_59);
    // q = 0
    result_by_lpq.set(c_Key3(5, 1, 0), common_term_60);
    result_by_q.set(c_Key1(0), common_term_60);
    // q = 1
    result_by_lpq.set(c_Key3(5, 1, 1), common_term_61);
    result_by_q.set(c_Key1(1), common_term_61);
    // q = 2
    result_by_lpq.set(c_Key3(5, 1, 2), common_term_62);
    result_by_q.set(c_Key1(2), common_term_62);
    // q = 3
    result_by_lpq.set(c_Key3(5, 1, 3), common_term_63);
    result_by_q.set(c_Key1(3), common_term_63);
    // q = 4
    result_by_lpq.set(c_Key3(5, 1, 4), common_term_64);
    result_by_q.set(c_Key1(4), common_term_64);
    // q = 5
    result_by_lpq.set(c_Key3(5, 1, 5), common_term_65);
    result_by_q.set(c_Key1(5), common_term_65);
    // q = 6
    result_by_lpq.set(c_Key3(5, 1, 6), common_term_66);
    result_by_q.set(c_Key1(6), common_term_66);
    // q = 7
    result_by_lpq.set(c_Key3(5, 1, 7), common_term_67);
    result_by_q.set(c_Key1(7), common_term_67);
    // q = 8
    result_by_lpq.set(c_Key3(5, 1, 8), common_term_68);
    result_by_q.set(c_Key1(8), common_term_68);
    // q = 9
    result_by_lpq.set(c_Key3(5, 1, 9), common_term_69);
    result_by_q.set(c_Key1(9), common_term_69);
    // q = 10
    result_by_lpq.set(c_Key3(5, 1, 10), common_term_70);
    result_by_q.set(c_Key1(10), common_term_70);
    // q = 11
    result_by_lpq.set(c_Key3(5, 1, 11), common_term_71);
    result_by_q.set(c_Key1(11), common_term_71);
    // q = 12
    result_by_lpq.set(c_Key3(5, 1, 12), common_term_72);
    result_by_q.set(c_Key1(12), common_term_72);
    // q = 13
    result_by_lpq.set(c_Key3(5, 1, 13), common_term_73);
    result_by_q.set(c_Key1(13), common_term_73);
    // q = 14
    result_by_lpq.set(c_Key3(5, 1, 14), common_term_74);
    result_by_q.set(c_Key1(14), common_term_74);
    // q = 15
    result_by_lpq.set(c_Key3(5, 1, 15), common_term_75);
    result_by_q.set(c_Key1(15), common_term_75);
    // q = 16
    result_by_lpq.set(c_Key3(5, 1, 16), common_term_76);
    result_by_q.set(c_Key1(16), common_term_76);
    // q = 17
    result_by_lpq.set(c_Key3(5, 1, 17), common_term_77);
    result_by_q.set(c_Key1(17), common_term_77);
    // q = 18
    result_by_lpq.set(c_Key3(5, 1, 18), common_term_78);
    result_by_q.set(c_Key1(18), common_term_78);
    // q = 19
    result_by_lpq.set(c_Key3(5, 1, 19), common_term_79);
    result_by_q.set(c_Key1(19), common_term_79);
    // q = 20
    result_by_lpq.set(c_Key3(5, 1, 20), common_term_80);
    result_by_q.set(c_Key1(20), common_term_80);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 1), result_by_q);
    result_by_q.clear();

    // l , p = (5, 2).
    // q = -20
    result_by_lpq.set(c_Key3(5, 2, -20), common_term_81);
    result_by_q.set(c_Key1(-20), common_term_81);
    // q = -19
    result_by_lpq.set(c_Key3(5, 2, -19), common_term_82);
    result_by_q.set(c_Key1(-19), common_term_82);
    // q = -18
    result_by_lpq.set(c_Key3(5, 2, -18), common_term_83);
    result_by_q.set(c_Key1(-18), common_term_83);
    // q = -17
    result_by_lpq.set(c_Key3(5, 2, -17), common_term_84);
    result_by_q.set(c_Key1(-17), common_term_84);
    // q = -16
    result_by_lpq.set(c_Key3(5, 2, -16), common_term_85);
    result_by_q.set(c_Key1(-16), common_term_85);
    // q = -15
    result_by_lpq.set(c_Key3(5, 2, -15), common_term_86);
    result_by_q.set(c_Key1(-15), common_term_86);
    // q = -14
    result_by_lpq.set(c_Key3(5, 2, -14), common_term_87);
    result_by_q.set(c_Key1(-14), common_term_87);
    // q = -13
    result_by_lpq.set(c_Key3(5, 2, -13), common_term_88);
    result_by_q.set(c_Key1(-13), common_term_88);
    // q = -12
    result_by_lpq.set(c_Key3(5, 2, -12), common_term_89);
    result_by_q.set(c_Key1(-12), common_term_89);
    // q = -11
    result_by_lpq.set(c_Key3(5, 2, -11), common_term_90);
    result_by_q.set(c_Key1(-11), common_term_90);
    // q = -10
    result_by_lpq.set(c_Key3(5, 2, -10), common_term_91);
    result_by_q.set(c_Key1(-10), common_term_91);
    // q = -9
    result_by_lpq.set(c_Key3(5, 2, -9), common_term_92);
    result_by_q.set(c_Key1(-9), common_term_92);
    // q = -8
    result_by_lpq.set(c_Key3(5, 2, -8), common_term_93);
    result_by_q.set(c_Key1(-8), common_term_93);
    // q = -7
    result_by_lpq.set(c_Key3(5, 2, -7), common_term_94);
    result_by_q.set(c_Key1(-7), common_term_94);
    // q = -6
    result_by_lpq.set(c_Key3(5, 2, -6), common_term_95);
    result_by_q.set(c_Key1(-6), common_term_95);
    // q = -5
    result_by_lpq.set(c_Key3(5, 2, -5), common_term_96);
    result_by_q.set(c_Key1(-5), common_term_96);
    // q = -4
    result_by_lpq.set(c_Key3(5, 2, -4), common_term_97);
    result_by_q.set(c_Key1(-4), common_term_97);
    // q = -3
    result_by_lpq.set(c_Key3(5, 2, -3), common_term_98);
    result_by_q.set(c_Key1(-3), common_term_98);
    // q = -2
    result_by_lpq.set(c_Key3(5, 2, -2), common_term_99);
    result_by_q.set(c_Key1(-2), common_term_99);
    // q = -1
    result_by_lpq.set(c_Key3(5, 2, -1), common_term_100);
    result_by_q.set(c_Key1(-1), common_term_100);
    // q = 0
    result_by_lpq.set(c_Key3(5, 2, 0), common_term_101);
    result_by_q.set(c_Key1(0), common_term_101);
    // q = 1
    result_by_lpq.set(c_Key3(5, 2, 1), common_term_102);
    result_by_q.set(c_Key1(1), common_term_102);
    // q = 2
    result_by_lpq.set(c_Key3(5, 2, 2), common_term_103);
    result_by_q.set(c_Key1(2), common_term_103);
    // q = 3
    result_by_lpq.set(c_Key3(5, 2, 3), common_term_104);
    result_by_q.set(c_Key1(3), common_term_104);
    // q = 4
    result_by_lpq.set(c_Key3(5, 2, 4), common_term_105);
    result_by_q.set(c_Key1(4), common_term_105);
    // q = 5
    result_by_lpq.set(c_Key3(5, 2, 5), common_term_106);
    result_by_q.set(c_Key1(5), common_term_106);
    // q = 6
    result_by_lpq.set(c_Key3(5, 2, 6), common_term_107);
    result_by_q.set(c_Key1(6), common_term_107);
    // q = 7
    result_by_lpq.set(c_Key3(5, 2, 7), common_term_108);
    result_by_q.set(c_Key1(7), common_term_108);
    // q = 8
    result_by_lpq.set(c_Key3(5, 2, 8), common_term_109);
    result_by_q.set(c_Key1(8), common_term_109);
    // q = 9
    result_by_lpq.set(c_Key3(5, 2, 9), common_term_110);
    result_by_q.set(c_Key1(9), common_term_110);
    // q = 10
    result_by_lpq.set(c_Key3(5, 2, 10), common_term_111);
    result_by_q.set(c_Key1(10), common_term_111);
    // q = 11
    result_by_lpq.set(c_Key3(5, 2, 11), common_term_112);
    result_by_q.set(c_Key1(11), common_term_112);
    // q = 12
    result_by_lpq.set(c_Key3(5, 2, 12), common_term_113);
    result_by_q.set(c_Key1(12), common_term_113);
    // q = 13
    result_by_lpq.set(c_Key3(5, 2, 13), common_term_114);
    result_by_q.set(c_Key1(13), common_term_114);
    // q = 14
    result_by_lpq.set(c_Key3(5, 2, 14), common_term_115);
    result_by_q.set(c_Key1(14), common_term_115);
    // q = 15
    result_by_lpq.set(c_Key3(5, 2, 15), common_term_116);
    result_by_q.set(c_Key1(15), common_term_116);
    // q = 16
    result_by_lpq.set(c_Key3(5, 2, 16), common_term_117);
    result_by_q.set(c_Key1(16), common_term_117);
    // q = 17
    result_by_lpq.set(c_Key3(5, 2, 17), common_term_118);
    result_by_q.set(c_Key1(17), common_term_118);
    // q = 18
    result_by_lpq.set(c_Key3(5, 2, 18), common_term_119);
    result_by_q.set(c_Key1(18), common_term_119);
    // q = 19
    result_by_lpq.set(c_Key3(5, 2, 19), common_term_120);
    result_by_q.set(c_Key1(19), common_term_120);
    // q = 20
    result_by_lpq.set(c_Key3(5, 2, 20), common_term_121);
    result_by_q.set(c_Key1(20), common_term_121);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 2), result_by_q);
    result_by_q.clear();

    // l , p = (5, 3).
    // q = -20
    result_by_lpq.set(c_Key3(5, 3, -20), common_term_121);
    result_by_q.set(c_Key1(-20), common_term_121);
    // q = -19
    result_by_lpq.set(c_Key3(5, 3, -19), common_term_120);
    result_by_q.set(c_Key1(-19), common_term_120);
    // q = -18
    result_by_lpq.set(c_Key3(5, 3, -18), common_term_119);
    result_by_q.set(c_Key1(-18), common_term_119);
    // q = -17
    result_by_lpq.set(c_Key3(5, 3, -17), common_term_118);
    result_by_q.set(c_Key1(-17), common_term_118);
    // q = -16
    result_by_lpq.set(c_Key3(5, 3, -16), common_term_117);
    result_by_q.set(c_Key1(-16), common_term_117);
    // q = -15
    result_by_lpq.set(c_Key3(5, 3, -15), common_term_116);
    result_by_q.set(c_Key1(-15), common_term_116);
    // q = -14
    result_by_lpq.set(c_Key3(5, 3, -14), common_term_115);
    result_by_q.set(c_Key1(-14), common_term_115);
    // q = -13
    result_by_lpq.set(c_Key3(5, 3, -13), common_term_114);
    result_by_q.set(c_Key1(-13), common_term_114);
    // q = -12
    result_by_lpq.set(c_Key3(5, 3, -12), common_term_113);
    result_by_q.set(c_Key1(-12), common_term_113);
    // q = -11
    result_by_lpq.set(c_Key3(5, 3, -11), common_term_112);
    result_by_q.set(c_Key1(-11), common_term_112);
    // q = -10
    result_by_lpq.set(c_Key3(5, 3, -10), common_term_111);
    result_by_q.set(c_Key1(-10), common_term_111);
    // q = -9
    result_by_lpq.set(c_Key3(5, 3, -9), common_term_110);
    result_by_q.set(c_Key1(-9), common_term_110);
    // q = -8
    result_by_lpq.set(c_Key3(5, 3, -8), common_term_109);
    result_by_q.set(c_Key1(-8), common_term_109);
    // q = -7
    result_by_lpq.set(c_Key3(5, 3, -7), common_term_108);
    result_by_q.set(c_Key1(-7), common_term_108);
    // q = -6
    result_by_lpq.set(c_Key3(5, 3, -6), common_term_107);
    result_by_q.set(c_Key1(-6), common_term_107);
    // q = -5
    result_by_lpq.set(c_Key3(5, 3, -5), common_term_106);
    result_by_q.set(c_Key1(-5), common_term_106);
    // q = -4
    result_by_lpq.set(c_Key3(5, 3, -4), common_term_105);
    result_by_q.set(c_Key1(-4), common_term_105);
    // q = -3
    result_by_lpq.set(c_Key3(5, 3, -3), common_term_104);
    result_by_q.set(c_Key1(-3), common_term_104);
    // q = -2
    result_by_lpq.set(c_Key3(5, 3, -2), common_term_103);
    result_by_q.set(c_Key1(-2), common_term_103);
    // q = -1
    result_by_lpq.set(c_Key3(5, 3, -1), common_term_102);
    result_by_q.set(c_Key1(-1), common_term_102);
    // q = 0
    result_by_lpq.set(c_Key3(5, 3, 0), common_term_101);
    result_by_q.set(c_Key1(0), common_term_101);
    // q = 1
    result_by_lpq.set(c_Key3(5, 3, 1), common_term_100);
    result_by_q.set(c_Key1(1), common_term_100);
    // q = 2
    result_by_lpq.set(c_Key3(5, 3, 2), common_term_99);
    result_by_q.set(c_Key1(2), common_term_99);
    // q = 3
    result_by_lpq.set(c_Key3(5, 3, 3), common_term_98);
    result_by_q.set(c_Key1(3), common_term_98);
    // q = 4
    result_by_lpq.set(c_Key3(5, 3, 4), common_term_97);
    result_by_q.set(c_Key1(4), common_term_97);
    // q = 5
    result_by_lpq.set(c_Key3(5, 3, 5), common_term_96);
    result_by_q.set(c_Key1(5), common_term_96);
    // q = 6
    result_by_lpq.set(c_Key3(5, 3, 6), common_term_95);
    result_by_q.set(c_Key1(6), common_term_95);
    // q = 7
    result_by_lpq.set(c_Key3(5, 3, 7), common_term_94);
    result_by_q.set(c_Key1(7), common_term_94);
    // q = 8
    result_by_lpq.set(c_Key3(5, 3, 8), common_term_93);
    result_by_q.set(c_Key1(8), common_term_93);
    // q = 9
    result_by_lpq.set(c_Key3(5, 3, 9), common_term_92);
    result_by_q.set(c_Key1(9), common_term_92);
    // q = 10
    result_by_lpq.set(c_Key3(5, 3, 10), common_term_91);
    result_by_q.set(c_Key1(10), common_term_91);
    // q = 11
    result_by_lpq.set(c_Key3(5, 3, 11), common_term_90);
    result_by_q.set(c_Key1(11), common_term_90);
    // q = 12
    result_by_lpq.set(c_Key3(5, 3, 12), common_term_89);
    result_by_q.set(c_Key1(12), common_term_89);
    // q = 13
    result_by_lpq.set(c_Key3(5, 3, 13), common_term_88);
    result_by_q.set(c_Key1(13), common_term_88);
    // q = 14
    result_by_lpq.set(c_Key3(5, 3, 14), common_term_87);
    result_by_q.set(c_Key1(14), common_term_87);
    // q = 15
    result_by_lpq.set(c_Key3(5, 3, 15), common_term_86);
    result_by_q.set(c_Key1(15), common_term_86);
    // q = 16
    result_by_lpq.set(c_Key3(5, 3, 16), common_term_85);
    result_by_q.set(c_Key1(16), common_term_85);
    // q = 17
    result_by_lpq.set(c_Key3(5, 3, 17), common_term_84);
    result_by_q.set(c_Key1(17), common_term_84);
    // q = 18
    result_by_lpq.set(c_Key3(5, 3, 18), common_term_83);
    result_by_q.set(c_Key1(18), common_term_83);
    // q = 19
    result_by_lpq.set(c_Key3(5, 3, 19), common_term_82);
    result_by_q.set(c_Key1(19), common_term_82);
    // q = 20
    result_by_lpq.set(c_Key3(5, 3, 20), common_term_81);
    result_by_q.set(c_Key1(20), common_term_81);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 3), result_by_q);
    result_by_q.clear();

    // l , p = (5, 4).
    // q = -20
    result_by_lpq.set(c_Key3(5, 4, -20), common_term_80);
    result_by_q.set(c_Key1(-20), common_term_80);
    // q = -19
    result_by_lpq.set(c_Key3(5, 4, -19), common_term_79);
    result_by_q.set(c_Key1(-19), common_term_79);
    // q = -18
    result_by_lpq.set(c_Key3(5, 4, -18), common_term_78);
    result_by_q.set(c_Key1(-18), common_term_78);
    // q = -17
    result_by_lpq.set(c_Key3(5, 4, -17), common_term_77);
    result_by_q.set(c_Key1(-17), common_term_77);
    // q = -16
    result_by_lpq.set(c_Key3(5, 4, -16), common_term_76);
    result_by_q.set(c_Key1(-16), common_term_76);
    // q = -15
    result_by_lpq.set(c_Key3(5, 4, -15), common_term_75);
    result_by_q.set(c_Key1(-15), common_term_75);
    // q = -14
    result_by_lpq.set(c_Key3(5, 4, -14), common_term_74);
    result_by_q.set(c_Key1(-14), common_term_74);
    // q = -13
    result_by_lpq.set(c_Key3(5, 4, -13), common_term_73);
    result_by_q.set(c_Key1(-13), common_term_73);
    // q = -12
    result_by_lpq.set(c_Key3(5, 4, -12), common_term_72);
    result_by_q.set(c_Key1(-12), common_term_72);
    // q = -11
    result_by_lpq.set(c_Key3(5, 4, -11), common_term_71);
    result_by_q.set(c_Key1(-11), common_term_71);
    // q = -10
    result_by_lpq.set(c_Key3(5, 4, -10), common_term_70);
    result_by_q.set(c_Key1(-10), common_term_70);
    // q = -9
    result_by_lpq.set(c_Key3(5, 4, -9), common_term_69);
    result_by_q.set(c_Key1(-9), common_term_69);
    // q = -8
    result_by_lpq.set(c_Key3(5, 4, -8), common_term_68);
    result_by_q.set(c_Key1(-8), common_term_68);
    // q = -7
    result_by_lpq.set(c_Key3(5, 4, -7), common_term_67);
    result_by_q.set(c_Key1(-7), common_term_67);
    // q = -6
    result_by_lpq.set(c_Key3(5, 4, -6), common_term_66);
    result_by_q.set(c_Key1(-6), common_term_66);
    // q = -5
    result_by_lpq.set(c_Key3(5, 4, -5), common_term_65);
    result_by_q.set(c_Key1(-5), common_term_65);
    // q = -4
    result_by_lpq.set(c_Key3(5, 4, -4), common_term_64);
    result_by_q.set(c_Key1(-4), common_term_64);
    // q = -3
    result_by_lpq.set(c_Key3(5, 4, -3), common_term_63);
    result_by_q.set(c_Key1(-3), common_term_63);
    // q = -2
    result_by_lpq.set(c_Key3(5, 4, -2), common_term_62);
    result_by_q.set(c_Key1(-2), common_term_62);
    // q = -1
    result_by_lpq.set(c_Key3(5, 4, -1), common_term_61);
    result_by_q.set(c_Key1(-1), common_term_61);
    // q = 0
    result_by_lpq.set(c_Key3(5, 4, 0), common_term_60);
    result_by_q.set(c_Key1(0), common_term_60);
    // q = 1
    result_by_lpq.set(c_Key3(5, 4, 1), common_term_59);
    result_by_q.set(c_Key1(1), common_term_59);
    // q = 2
    result_by_lpq.set(c_Key3(5, 4, 2), common_term_58);
    result_by_q.set(c_Key1(2), common_term_58);
    // q = 3
    result_by_lpq.set(c_Key3(5, 4, 3), common_term_57);
    result_by_q.set(c_Key1(3), common_term_57);
    // q = 4
    result_by_lpq.set(c_Key3(5, 4, 4), common_term_56);
    result_by_q.set(c_Key1(4), common_term_56);
    // q = 5
    result_by_lpq.set(c_Key3(5, 4, 5), common_term_55);
    result_by_q.set(c_Key1(5), common_term_55);
    // q = 6
    result_by_lpq.set(c_Key3(5, 4, 6), common_term_54);
    result_by_q.set(c_Key1(6), common_term_54);
    // q = 7
    result_by_lpq.set(c_Key3(5, 4, 7), common_term_53);
    result_by_q.set(c_Key1(7), common_term_53);
    // q = 8
    result_by_lpq.set(c_Key3(5, 4, 8), common_term_52);
    result_by_q.set(c_Key1(8), common_term_52);
    // q = 9
    result_by_lpq.set(c_Key3(5, 4, 9), common_term_51);
    result_by_q.set(c_Key1(9), common_term_51);
    // q = 10
    result_by_lpq.set(c_Key3(5, 4, 10), common_term_50);
    result_by_q.set(c_Key1(10), common_term_50);
    // q = 11
    result_by_lpq.set(c_Key3(5, 4, 11), common_term_49);
    result_by_q.set(c_Key1(11), common_term_49);
    // q = 12
    result_by_lpq.set(c_Key3(5, 4, 12), common_term_48);
    result_by_q.set(c_Key1(12), common_term_48);
    // q = 13
    result_by_lpq.set(c_Key3(5, 4, 13), common_term_47);
    result_by_q.set(c_Key1(13), common_term_47);
    // q = 14
    result_by_lpq.set(c_Key3(5, 4, 14), common_term_46);
    result_by_q.set(c_Key1(14), common_term_46);
    // q = 15
    result_by_lpq.set(c_Key3(5, 4, 15), common_term_45);
    result_by_q.set(c_Key1(15), common_term_45);
    // q = 16
    result_by_lpq.set(c_Key3(5, 4, 16), common_term_44);
    result_by_q.set(c_Key1(16), common_term_44);
    // q = 17
    result_by_lpq.set(c_Key3(5, 4, 17), common_term_43);
    result_by_q.set(c_Key1(17), common_term_43);
    // q = 18
    result_by_lpq.set(c_Key3(5, 4, 18), common_term_42);
    result_by_q.set(c_Key1(18), common_term_42);
    // q = 19
    result_by_lpq.set(c_Key3(5, 4, 19), common_term_41);
    result_by_q.set(c_Key1(19), common_term_41);
    // q = 20
    result_by_lpq.set(c_Key3(5, 4, 20), common_term_40);
    result_by_q.set(c_Key1(20), common_term_40);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 4), result_by_q);
    result_by_q.clear();

    // l , p = (5, 5).
    // q = -20
    result_by_lpq.set(c_Key3(5, 5, -20), common_term_39);
    result_by_q.set(c_Key1(-20), common_term_39);
    // q = -19
    result_by_lpq.set(c_Key3(5, 5, -19), common_term_38);
    result_by_q.set(c_Key1(-19), common_term_38);
    // q = -18
    result_by_lpq.set(c_Key3(5, 5, -18), common_term_37);
    result_by_q.set(c_Key1(-18), common_term_37);
    // q = -17
    result_by_lpq.set(c_Key3(5, 5, -17), common_term_36);
    result_by_q.set(c_Key1(-17), common_term_36);
    // q = -16
    result_by_lpq.set(c_Key3(5, 5, -16), common_term_35);
    result_by_q.set(c_Key1(-16), common_term_35);
    // q = -15
    result_by_lpq.set(c_Key3(5, 5, -15), common_term_34);
    result_by_q.set(c_Key1(-15), common_term_34);
    // q = -14
    result_by_lpq.set(c_Key3(5, 5, -14), common_term_33);
    result_by_q.set(c_Key1(-14), common_term_33);
    // q = -13
    result_by_lpq.set(c_Key3(5, 5, -13), common_term_32);
    result_by_q.set(c_Key1(-13), common_term_32);
    // q = -12
    result_by_lpq.set(c_Key3(5, 5, -12), common_term_31);
    result_by_q.set(c_Key1(-12), common_term_31);
    // q = -11
    result_by_lpq.set(c_Key3(5, 5, -11), common_term_30);
    result_by_q.set(c_Key1(-11), common_term_30);
    // q = -10
    result_by_lpq.set(c_Key3(5, 5, -10), common_term_29);
    result_by_q.set(c_Key1(-10), common_term_29);
    // q = -9
    result_by_lpq.set(c_Key3(5, 5, -9), common_term_28);
    result_by_q.set(c_Key1(-9), common_term_28);
    // q = -8
    result_by_lpq.set(c_Key3(5, 5, -8), common_term_27);
    result_by_q.set(c_Key1(-8), common_term_27);
    // q = -7
    result_by_lpq.set(c_Key3(5, 5, -7), common_term_26);
    result_by_q.set(c_Key1(-7), common_term_26);
    // q = -6
    result_by_lpq.set(c_Key3(5, 5, -6), common_term_25);
    result_by_q.set(c_Key1(-6), common_term_25);
    // q = -5
    result_by_lpq.set(c_Key3(5, 5, -5), common_term_24);
    result_by_q.set(c_Key1(-5), common_term_24);
    // q = -4
    result_by_lpq.set(c_Key3(5, 5, -4), common_term_23);
    result_by_q.set(c_Key1(-4), common_term_23);
    // q = -3
    result_by_lpq.set(c_Key3(5, 5, -3), common_term_22);
    result_by_q.set(c_Key1(-3), common_term_22);
    // q = -2
    result_by_lpq.set(c_Key3(5, 5, -2), common_term_21);
    result_by_q.set(c_Key1(-2), common_term_21);
    // q = -1
    result_by_lpq.set(c_Key3(5, 5, -1), common_term_20);
    result_by_q.set(c_Key1(-1), common_term_20);
    // q = 0
    result_by_lpq.set(c_Key3(5, 5, 0), common_term_19);
    result_by_q.set(c_Key1(0), common_term_19);
    // q = 1
    result_by_lpq.set(c_Key3(5, 5, 1), common_term_18);
    result_by_q.set(c_Key1(1), common_term_18);
    // q = 2
    result_by_lpq.set(c_Key3(5, 5, 2), common_term_17);
    result_by_q.set(c_Key1(2), common_term_17);
    // q = 3
    result_by_lpq.set(c_Key3(5, 5, 3), common_term_16);
    result_by_q.set(c_Key1(3), common_term_16);
    // q = 4
    result_by_lpq.set(c_Key3(5, 5, 4), common_term_15);
    result_by_q.set(c_Key1(4), common_term_15);
    // q = 6
    result_by_lpq.set(c_Key3(5, 5, 6), common_term_14);
    result_by_q.set(c_Key1(6), common_term_14);
    // q = 7
    result_by_lpq.set(c_Key3(5, 5, 7), common_term_13);
    result_by_q.set(c_Key1(7), common_term_13);
    // q = 8
    result_by_lpq.set(c_Key3(5, 5, 8), common_term_12);
    result_by_q.set(c_Key1(8), common_term_12);
    // q = 9
    result_by_lpq.set(c_Key3(5, 5, 9), common_term_11);
    result_by_q.set(c_Key1(9), common_term_11);
    // q = 10
    result_by_lpq.set(c_Key3(5, 5, 10), common_term_10);
    result_by_q.set(c_Key1(10), common_term_10);
    // q = 11
    result_by_lpq.set(c_Key3(5, 5, 11), common_term_9);
    result_by_q.set(c_Key1(11), common_term_9);
    // q = 12
    result_by_lpq.set(c_Key3(5, 5, 12), common_term_8);
    result_by_q.set(c_Key1(12), common_term_8);
    // q = 13
    result_by_lpq.set(c_Key3(5, 5, 13), common_term_7);
    result_by_q.set(c_Key1(13), common_term_7);
    // q = 14
    result_by_lpq.set(c_Key3(5, 5, 14), common_term_6);
    result_by_q.set(c_Key1(14), common_term_6);
    // q = 15
    result_by_lpq.set(c_Key3(5, 5, 15), common_term_5);
    result_by_q.set(c_Key1(15), common_term_5);
    // q = 16
    result_by_lpq.set(c_Key3(5, 5, 16), common_term_4);
    result_by_q.set(c_Key1(16), common_term_4);
    // q = 17
    result_by_lpq.set(c_Key3(5, 5, 17), common_term_3);
    result_by_q.set(c_Key1(17), common_term_3);
    // q = 18
    result_by_lpq.set(c_Key3(5, 5, 18), common_term_2);
    result_by_q.set(c_Key1(18), common_term_2);
    // q = 19
    result_by_lpq.set(c_Key3(5, 5, 19), common_term_1);
    result_by_q.set(c_Key1(19), common_term_1);
    // q = 20
    result_by_lpq.set(c_Key3(5, 5, 20), common_term_0);
    result_by_q.set(c_Key1(20), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(5, 5), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}
