#pragma once

#include "eccentricity_common_.hpp"

EccentricityFuncOutput c_eccentricity_function_l3_e1(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 3.
    // Functions are taylor expanded and truncated at eccentricity^1.

    // Functions Calculated for l = 3.

    c_IntMap<c_Key3, double> result_by_lpq(6);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(4);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = eccentricity*std::pow(1.0 - eccentricity_2, -2.5);

    c_IntMap<c_Key1, double> result_by_q(2);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (3, 0).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(3, 0, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 0), result_by_q);
    result_by_q.clear();

    // l , p = (3, 1).
    // q = -1
    result_by_lpq.set(c_Key3(3, 1, -1), common_term_0);
    result_by_q.set(c_Key1(-1), common_term_0);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(3, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 1), result_by_q);
    result_by_q.clear();

    // l , p = (3, 2).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(3, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(3, 2, 1), common_term_0);
    result_by_q.set(c_Key1(1), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 2), result_by_q);
    result_by_q.clear();

    // l , p = (3, 3).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(3, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 3), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l3_e2(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 3.
    // Functions are taylor expanded and truncated at eccentricity^2.

    // Functions Calculated for l = 3.

    c_IntMap<c_Key3, double> result_by_lpq(12);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(4);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = -eccentricity;
    double common_term_1 = 5.0*eccentricity;
    double common_term_2 = eccentricity*std::pow(1.0 - eccentricity_2, -2.5);
    double common_term_3 = 3.0*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(3);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (3, 0).
    // q = -1
    result_by_lpq.set(c_Key3(3, 0, -1), common_term_0);
    result_by_q.set(c_Key1(-1), common_term_0);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(3, 0, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(3, 0, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 0), result_by_q);
    result_by_q.clear();

    // l , p = (3, 1).
    // q = -1
    result_by_lpq.set(c_Key3(3, 1, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(3, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(3, 1, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 1), result_by_q);
    result_by_q.clear();

    // l , p = (3, 2).
    // q = -1
    result_by_lpq.set(c_Key3(3, 2, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(3, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(3, 2, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 2), result_by_q);
    result_by_q.clear();

    // l , p = (3, 3).
    // q = -1
    result_by_lpq.set(c_Key3(3, 3, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(3, 3, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(3, 3, 1), common_term_0);
    result_by_q.set(c_Key1(1), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 3), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l3_e3(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 3.
    // Functions are taylor expanded and truncated at eccentricity^3.

    // Functions Calculated for l = 3.

    c_IntMap<c_Key3, double> result_by_lpq(20);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(4);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = 0.125*eccentricity_2;
    double common_term_1 = -eccentricity;
    double common_term_2 = 1.0 - 6.0*eccentricity_2;
    double common_term_3 = 5.0*eccentricity;
    double common_term_4 = 15.875*eccentricity_2;
    double common_term_5 = 1.375*eccentricity_2;
    double common_term_6 = eccentricity*std::pow(1.0 - eccentricity_2, -2.5);
    double common_term_7 = 2.0*eccentricity_2 + 1.0;
    double common_term_8 = 3.0*eccentricity;
    double common_term_9 = 6.625*eccentricity_2;

    c_IntMap<c_Key1, double> result_by_q(5);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (3, 0).
    // q = -2
    result_by_lpq.set(c_Key3(3, 0, -2), common_term_0);
    result_by_q.set(c_Key1(-2), common_term_0);
    // q = -1
    result_by_lpq.set(c_Key3(3, 0, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    result_by_lpq.set(c_Key3(3, 0, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(3, 0, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // q = 2
    result_by_lpq.set(c_Key3(3, 0, 2), common_term_4);
    result_by_q.set(c_Key1(2), common_term_4);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 0), result_by_q);
    result_by_q.clear();

    // l , p = (3, 1).
    // q = -2
    result_by_lpq.set(c_Key3(3, 1, -2), common_term_5);
    result_by_q.set(c_Key1(-2), common_term_5);
    // q = -1
    result_by_lpq.set(c_Key3(3, 1, -1), common_term_6);
    result_by_q.set(c_Key1(-1), common_term_6);
    // q = 0
    result_by_lpq.set(c_Key3(3, 1, 0), common_term_7);
    result_by_q.set(c_Key1(0), common_term_7);
    // q = 1
    result_by_lpq.set(c_Key3(3, 1, 1), common_term_8);
    result_by_q.set(c_Key1(1), common_term_8);
    // q = 2
    result_by_lpq.set(c_Key3(3, 1, 2), common_term_9);
    result_by_q.set(c_Key1(2), common_term_9);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 1), result_by_q);
    result_by_q.clear();

    // l , p = (3, 2).
    // q = -2
    result_by_lpq.set(c_Key3(3, 2, -2), common_term_9);
    result_by_q.set(c_Key1(-2), common_term_9);
    // q = -1
    result_by_lpq.set(c_Key3(3, 2, -1), common_term_8);
    result_by_q.set(c_Key1(-1), common_term_8);
    // q = 0
    result_by_lpq.set(c_Key3(3, 2, 0), common_term_7);
    result_by_q.set(c_Key1(0), common_term_7);
    // q = 1
    result_by_lpq.set(c_Key3(3, 2, 1), common_term_6);
    result_by_q.set(c_Key1(1), common_term_6);
    // q = 2
    result_by_lpq.set(c_Key3(3, 2, 2), common_term_5);
    result_by_q.set(c_Key1(2), common_term_5);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 2), result_by_q);
    result_by_q.clear();

    // l , p = (3, 3).
    // q = -2
    result_by_lpq.set(c_Key3(3, 3, -2), common_term_4);
    result_by_q.set(c_Key1(-2), common_term_4);
    // q = -1
    result_by_lpq.set(c_Key3(3, 3, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    result_by_lpq.set(c_Key3(3, 3, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(3, 3, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // q = 2
    result_by_lpq.set(c_Key3(3, 3, 2), common_term_0);
    result_by_q.set(c_Key1(2), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 3), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l3_e4(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 3.
    // Functions are taylor expanded and truncated at eccentricity^4.

    // Functions Calculated for l = 3.

    c_IntMap<c_Key3, double> result_by_lpq(26);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(4);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double common_term_0 = 0.125*eccentricity_2;
    double common_term_1 = 1.25*eccentricity_3 - eccentricity;
    double common_term_2 = 1.0 - 6.0*eccentricity_2;
    double common_term_3 = -22.0*eccentricity_3 + 5.0*eccentricity;
    double common_term_4 = 15.875*eccentricity_2;
    double common_term_5 = 40.75*eccentricity_3;
    double common_term_6 = 1.9166666666666667*eccentricity_3;
    double common_term_7 = 1.375*eccentricity_2;
    double common_term_8 = eccentricity*std::pow(1.0 - eccentricity_2, -2.5);
    double common_term_9 = 2.0*eccentricity_2 + 1.0;
    double common_term_10 = 2.75*eccentricity_3 + 3.0*eccentricity;
    double common_term_11 = 6.625*eccentricity_2;
    double common_term_12 = 12.833333333333333*eccentricity_3;

    c_IntMap<c_Key1, double> result_by_q(7);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (3, 0).
    // q = -2
    result_by_lpq.set(c_Key3(3, 0, -2), common_term_0);
    result_by_q.set(c_Key1(-2), common_term_0);
    // q = -1
    result_by_lpq.set(c_Key3(3, 0, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    result_by_lpq.set(c_Key3(3, 0, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(3, 0, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // q = 2
    result_by_lpq.set(c_Key3(3, 0, 2), common_term_4);
    result_by_q.set(c_Key1(2), common_term_4);
    // q = 3
    result_by_lpq.set(c_Key3(3, 0, 3), common_term_5);
    result_by_q.set(c_Key1(3), common_term_5);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 0), result_by_q);
    result_by_q.clear();

    // l , p = (3, 1).
    // q = -3
    result_by_lpq.set(c_Key3(3, 1, -3), common_term_6);
    result_by_q.set(c_Key1(-3), common_term_6);
    // q = -2
    result_by_lpq.set(c_Key3(3, 1, -2), common_term_7);
    result_by_q.set(c_Key1(-2), common_term_7);
    // q = -1
    result_by_lpq.set(c_Key3(3, 1, -1), common_term_8);
    result_by_q.set(c_Key1(-1), common_term_8);
    // q = 0
    result_by_lpq.set(c_Key3(3, 1, 0), common_term_9);
    result_by_q.set(c_Key1(0), common_term_9);
    // q = 1
    result_by_lpq.set(c_Key3(3, 1, 1), common_term_10);
    result_by_q.set(c_Key1(1), common_term_10);
    // q = 2
    result_by_lpq.set(c_Key3(3, 1, 2), common_term_11);
    result_by_q.set(c_Key1(2), common_term_11);
    // q = 3
    result_by_lpq.set(c_Key3(3, 1, 3), common_term_12);
    result_by_q.set(c_Key1(3), common_term_12);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 1), result_by_q);
    result_by_q.clear();

    // l , p = (3, 2).
    // q = -3
    result_by_lpq.set(c_Key3(3, 2, -3), common_term_12);
    result_by_q.set(c_Key1(-3), common_term_12);
    // q = -2
    result_by_lpq.set(c_Key3(3, 2, -2), common_term_11);
    result_by_q.set(c_Key1(-2), common_term_11);
    // q = -1
    result_by_lpq.set(c_Key3(3, 2, -1), common_term_10);
    result_by_q.set(c_Key1(-1), common_term_10);
    // q = 0
    result_by_lpq.set(c_Key3(3, 2, 0), common_term_9);
    result_by_q.set(c_Key1(0), common_term_9);
    // q = 1
    result_by_lpq.set(c_Key3(3, 2, 1), common_term_8);
    result_by_q.set(c_Key1(1), common_term_8);
    // q = 2
    result_by_lpq.set(c_Key3(3, 2, 2), common_term_7);
    result_by_q.set(c_Key1(2), common_term_7);
    // q = 3
    result_by_lpq.set(c_Key3(3, 2, 3), common_term_6);
    result_by_q.set(c_Key1(3), common_term_6);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 2), result_by_q);
    result_by_q.clear();

    // l , p = (3, 3).
    // q = -3
    result_by_lpq.set(c_Key3(3, 3, -3), common_term_5);
    result_by_q.set(c_Key1(-3), common_term_5);
    // q = -2
    result_by_lpq.set(c_Key3(3, 3, -2), common_term_4);
    result_by_q.set(c_Key1(-2), common_term_4);
    // q = -1
    result_by_lpq.set(c_Key3(3, 3, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    result_by_lpq.set(c_Key3(3, 3, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(3, 3, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // q = 2
    result_by_lpq.set(c_Key3(3, 3, 2), common_term_0);
    result_by_q.set(c_Key1(2), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 3), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l3_e5(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 3.
    // Functions are taylor expanded and truncated at eccentricity^5.

    // Functions Calculated for l = 3.

    c_IntMap<c_Key3, double> result_by_lpq(34);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(4);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double common_term_0 = 0.0026041666666666667*eccentricity_4;
    double common_term_1 = 0.020833333333333333*eccentricity_4 + 0.125*eccentricity_2;
    double common_term_2 = 1.25*eccentricity_3 - eccentricity;
    double common_term_3 = 6.609375*eccentricity_4 - 6.0*eccentricity_2 + 1.0;
    double common_term_4 = -22.0*eccentricity_3 + 5.0*eccentricity;
    double common_term_5 = -63.854166666666667*eccentricity_4 + 15.875*eccentricity_2;
    double common_term_6 = 40.75*eccentricity_3;
    double common_term_7 = 92.221354166666667*eccentricity_4;
    double common_term_8 = 2.6796875*eccentricity_4;
    double common_term_9 = 1.9166666666666667*eccentricity_3;
    double common_term_10 = 3.0625*eccentricity_4 + 1.375*eccentricity_2;
    double common_term_11 = eccentricity*std::pow(1.0 - eccentricity_2, -2.5);
    double common_term_12 = 3.734375*eccentricity_4 + 2.0*eccentricity_2 + 1.0;
    double common_term_13 = 2.75*eccentricity_3 + 3.0*eccentricity;
    double common_term_14 = 2.4375*eccentricity_4 + 6.625*eccentricity_2;
    double common_term_15 = 12.833333333333333*eccentricity_3;
    double common_term_16 = 23.0859375*eccentricity_4;

    c_IntMap<c_Key1, double> result_by_q(9);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (3, 0).
    // q = -4
    result_by_lpq.set(c_Key3(3, 0, -4), common_term_0);
    result_by_q.set(c_Key1(-4), common_term_0);
    // q = -2
    result_by_lpq.set(c_Key3(3, 0, -2), common_term_1);
    result_by_q.set(c_Key1(-2), common_term_1);
    // q = -1
    result_by_lpq.set(c_Key3(3, 0, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    result_by_lpq.set(c_Key3(3, 0, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(3, 0, 1), common_term_4);
    result_by_q.set(c_Key1(1), common_term_4);
    // q = 2
    result_by_lpq.set(c_Key3(3, 0, 2), common_term_5);
    result_by_q.set(c_Key1(2), common_term_5);
    // q = 3
    result_by_lpq.set(c_Key3(3, 0, 3), common_term_6);
    result_by_q.set(c_Key1(3), common_term_6);
    // q = 4
    result_by_lpq.set(c_Key3(3, 0, 4), common_term_7);
    result_by_q.set(c_Key1(4), common_term_7);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 0), result_by_q);
    result_by_q.clear();

    // l , p = (3, 1).
    // q = -4
    result_by_lpq.set(c_Key3(3, 1, -4), common_term_8);
    result_by_q.set(c_Key1(-4), common_term_8);
    // q = -3
    result_by_lpq.set(c_Key3(3, 1, -3), common_term_9);
    result_by_q.set(c_Key1(-3), common_term_9);
    // q = -2
    result_by_lpq.set(c_Key3(3, 1, -2), common_term_10);
    result_by_q.set(c_Key1(-2), common_term_10);
    // q = -1
    result_by_lpq.set(c_Key3(3, 1, -1), common_term_11);
    result_by_q.set(c_Key1(-1), common_term_11);
    // q = 0
    result_by_lpq.set(c_Key3(3, 1, 0), common_term_12);
    result_by_q.set(c_Key1(0), common_term_12);
    // q = 1
    result_by_lpq.set(c_Key3(3, 1, 1), common_term_13);
    result_by_q.set(c_Key1(1), common_term_13);
    // q = 2
    result_by_lpq.set(c_Key3(3, 1, 2), common_term_14);
    result_by_q.set(c_Key1(2), common_term_14);
    // q = 3
    result_by_lpq.set(c_Key3(3, 1, 3), common_term_15);
    result_by_q.set(c_Key1(3), common_term_15);
    // q = 4
    result_by_lpq.set(c_Key3(3, 1, 4), common_term_16);
    result_by_q.set(c_Key1(4), common_term_16);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 1), result_by_q);
    result_by_q.clear();

    // l , p = (3, 2).
    // q = -4
    result_by_lpq.set(c_Key3(3, 2, -4), common_term_16);
    result_by_q.set(c_Key1(-4), common_term_16);
    // q = -3
    result_by_lpq.set(c_Key3(3, 2, -3), common_term_15);
    result_by_q.set(c_Key1(-3), common_term_15);
    // q = -2
    result_by_lpq.set(c_Key3(3, 2, -2), common_term_14);
    result_by_q.set(c_Key1(-2), common_term_14);
    // q = -1
    result_by_lpq.set(c_Key3(3, 2, -1), common_term_13);
    result_by_q.set(c_Key1(-1), common_term_13);
    // q = 0
    result_by_lpq.set(c_Key3(3, 2, 0), common_term_12);
    result_by_q.set(c_Key1(0), common_term_12);
    // q = 1
    result_by_lpq.set(c_Key3(3, 2, 1), common_term_11);
    result_by_q.set(c_Key1(1), common_term_11);
    // q = 2
    result_by_lpq.set(c_Key3(3, 2, 2), common_term_10);
    result_by_q.set(c_Key1(2), common_term_10);
    // q = 3
    result_by_lpq.set(c_Key3(3, 2, 3), common_term_9);
    result_by_q.set(c_Key1(3), common_term_9);
    // q = 4
    result_by_lpq.set(c_Key3(3, 2, 4), common_term_8);
    result_by_q.set(c_Key1(4), common_term_8);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 2), result_by_q);
    result_by_q.clear();

    // l , p = (3, 3).
    // q = -4
    result_by_lpq.set(c_Key3(3, 3, -4), common_term_7);
    result_by_q.set(c_Key1(-4), common_term_7);
    // q = -3
    result_by_lpq.set(c_Key3(3, 3, -3), common_term_6);
    result_by_q.set(c_Key1(-3), common_term_6);
    // q = -2
    result_by_lpq.set(c_Key3(3, 3, -2), common_term_5);
    result_by_q.set(c_Key1(-2), common_term_5);
    // q = -1
    result_by_lpq.set(c_Key3(3, 3, -1), common_term_4);
    result_by_q.set(c_Key1(-1), common_term_4);
    // q = 0
    result_by_lpq.set(c_Key3(3, 3, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(3, 3, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // q = 2
    result_by_lpq.set(c_Key3(3, 3, 2), common_term_1);
    result_by_q.set(c_Key1(2), common_term_1);
    // q = 4
    result_by_lpq.set(c_Key3(3, 3, 4), common_term_0);
    result_by_q.set(c_Key1(4), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 3), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l3_e10(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 3.
    // Functions are taylor expanded and truncated at eccentricity^10.

    // Functions Calculated for l = 3.

    c_IntMap<c_Key3, double> result_by_lpq(74);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(4);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_8 = eccentricity_4 * eccentricity_4;
    double eccentricity_9 = eccentricity_8 * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double eccentricity_6 = eccentricity_3 * eccentricity_3;
    double eccentricity_7 = eccentricity_6 * eccentricity;
    double eccentricity_5 = eccentricity_4 * eccentricity;
    double common_term_0 = 0.054241071428571429*eccentricity_9;
    double common_term_1 = 0.037844218905009921*eccentricity_8;
    double common_term_2 = 0.025396825396825397*eccentricity_9 + 0.025396825396825397*eccentricity_7;
    double common_term_3 = 0.016950334821428571*eccentricity_8 + 0.0158203125*eccentricity_6;
    double common_term_4 = 0.0083581349206349206*eccentricity_9 + 0.0090277777777777778*eccentricity_7 + 0.0083333333333333333*eccentricity_5;
    double common_term_5 = 0.0022325303819444444*eccentricity_8 + 0.0026041666666666667*eccentricity_6 + 0.0026041666666666667*eccentricity_4;
    double common_term_6 = 0.012771267361111111*eccentricity_8 + 0.017903645833333333*eccentricity_6 + 0.020833333333333333*eccentricity_4 + 0.125*eccentricity_2;
    double common_term_7 = 0.042361111111111111*eccentricity_9 + 0.079861111111111111*eccentricity_7 - 0.14583333333333333*eccentricity_5 + 1.25*eccentricity_3 - eccentricity;
    double common_term_8 = 0.45977783203125*eccentricity_8 - 1.953125*eccentricity_6 + 6.609375*eccentricity_4 - 6.0*eccentricity_2 + 1.0;
    double common_term_9 = 2.8431423611111111*eccentricity_9 - 10.888888888888889*eccentricity_7 + 25.291666666666667*eccentricity_5 - 22.0*eccentricity_3 + 5.0*eccentricity;
    double common_term_10 = -43.341200086805556*eccentricity_8 + 79.363606770833333*eccentricity_6 - 63.854166666666667*eccentricity_4 + 15.875*eccentricity_2;
    double common_term_11 = -141.94609375*eccentricity_9 + 217.8*eccentricity_7 - 161.0625*eccentricity_5 + 40.75*eccentricity_3;
    double common_term_12 = 542.69985622829861*eccentricity_8 - 369.51614583333333*eccentricity_6 + 92.221354166666667*eccentricity_4;
    double common_term_13 = 1257.2000992063492*eccentricity_9 - 791.42222222222222*eccentricity_7 + 191.90833333333333*eccentricity_5;
    double common_term_14 = -1608.5615373883929*eccentricity_8 + 376.0693359375*eccentricity_6;
    double common_term_15 = -3137.0445808531746*eccentricity_9 + 704.3968253968254*eccentricity_7;
    double common_term_16 = 1273.6733050633991*eccentricity_8;
    double common_term_17 = 2238.8170758928571*eccentricity_9;
    double common_term_18 = 14.123824680335097*eccentricity_9;
    double common_term_19 = 10.155956353081597*eccentricity_8;
    double common_term_20 = 5.8287946428571429*eccentricity_9 + 7.2933035714285714*eccentricity_7;
    double common_term_21 = 5.6296890500992064*eccentricity_8 + 5.2303602430555556*eccentricity_6;
    double common_term_22 = 7.7297123015873016*eccentricity_9 + 5.0854166666666667*eccentricity_7 + 3.7458333333333333*eccentricity_5;
    double common_term_23 = 6.7725341796875*eccentricity_8 + 4.4046875*eccentricity_6 + 2.6796875*eccentricity_4;
    double common_term_24 = 8.358275462962963*eccentricity_9 + 5.8989583333333333*eccentricity_7 + 3.7083333333333333*eccentricity_5 + 1.9166666666666667*eccentricity_3;
    double common_term_25 = 7.4273328993055556*eccentricity_8 + 5.0992838541666667*eccentricity_6 + 3.0625*eccentricity_4 + 1.375*eccentricity_2;
    double common_term_26 = eccentricity*std::pow(1.0 - eccentricity_2, -2.5);
    double common_term_27 = 8.0966050889756944*eccentricity_8 + 5.7690972222222222*eccentricity_6 + 3.734375*eccentricity_4 + 2.0*eccentricity_2 + 1.0;
    double common_term_28 = 9.6975694444444445*eccentricity_9 + 7.234375*eccentricity_7 + 5.1041666666666667*eccentricity_5 + 2.75*eccentricity_3 + 3.0*eccentricity;
    double common_term_29 = 8.72236328125*eccentricity_8 + 6.8759765625*eccentricity_6 + 2.4375*eccentricity_4 + 6.625*eccentricity_2;
    double common_term_30 = 10.008217592592593*eccentricity_9 + 9.8979166666666667*eccentricity_7 - 0.52083333333333333*eccentricity_5 + 12.833333333333333*eccentricity_3;
    double common_term_31 = 16.254814995659722*eccentricity_8 - 9.0182291666666667*eccentricity_6 + 23.0859375*eccentricity_4;
    double common_term_32 = 30.470870535714286*eccentricity_9 - 28.121875*eccentricity_7 + 39.5875*eccentricity_5;
    double common_term_33 = -66.439702690972222*eccentricity_8 + 65.638823784722222*eccentricity_6;
    double common_term_34 = -138.25908978174603*eccentricity_9 + 106.14940476190476*eccentricity_7;
    double common_term_35 = 168.38040684291295*eccentricity_8;
    double common_term_36 = 263.0183876212522*eccentricity_9;

    c_IntMap<c_Key1, double> result_by_q(19);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (3, 0).
    // q = -9
    result_by_lpq.set(c_Key3(3, 0, -9), common_term_0);
    result_by_q.set(c_Key1(-9), common_term_0);
    // q = -8
    result_by_lpq.set(c_Key3(3, 0, -8), common_term_1);
    result_by_q.set(c_Key1(-8), common_term_1);
    // q = -7
    result_by_lpq.set(c_Key3(3, 0, -7), common_term_2);
    result_by_q.set(c_Key1(-7), common_term_2);
    // q = -6
    result_by_lpq.set(c_Key3(3, 0, -6), common_term_3);
    result_by_q.set(c_Key1(-6), common_term_3);
    // q = -5
    result_by_lpq.set(c_Key3(3, 0, -5), common_term_4);
    result_by_q.set(c_Key1(-5), common_term_4);
    // q = -4
    result_by_lpq.set(c_Key3(3, 0, -4), common_term_5);
    result_by_q.set(c_Key1(-4), common_term_5);
    // q = -2
    result_by_lpq.set(c_Key3(3, 0, -2), common_term_6);
    result_by_q.set(c_Key1(-2), common_term_6);
    // q = -1
    result_by_lpq.set(c_Key3(3, 0, -1), common_term_7);
    result_by_q.set(c_Key1(-1), common_term_7);
    // q = 0
    result_by_lpq.set(c_Key3(3, 0, 0), common_term_8);
    result_by_q.set(c_Key1(0), common_term_8);
    // q = 1
    result_by_lpq.set(c_Key3(3, 0, 1), common_term_9);
    result_by_q.set(c_Key1(1), common_term_9);
    // q = 2
    result_by_lpq.set(c_Key3(3, 0, 2), common_term_10);
    result_by_q.set(c_Key1(2), common_term_10);
    // q = 3
    result_by_lpq.set(c_Key3(3, 0, 3), common_term_11);
    result_by_q.set(c_Key1(3), common_term_11);
    // q = 4
    result_by_lpq.set(c_Key3(3, 0, 4), common_term_12);
    result_by_q.set(c_Key1(4), common_term_12);
    // q = 5
    result_by_lpq.set(c_Key3(3, 0, 5), common_term_13);
    result_by_q.set(c_Key1(5), common_term_13);
    // q = 6
    result_by_lpq.set(c_Key3(3, 0, 6), common_term_14);
    result_by_q.set(c_Key1(6), common_term_14);
    // q = 7
    result_by_lpq.set(c_Key3(3, 0, 7), common_term_15);
    result_by_q.set(c_Key1(7), common_term_15);
    // q = 8
    result_by_lpq.set(c_Key3(3, 0, 8), common_term_16);
    result_by_q.set(c_Key1(8), common_term_16);
    // q = 9
    result_by_lpq.set(c_Key3(3, 0, 9), common_term_17);
    result_by_q.set(c_Key1(9), common_term_17);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 0), result_by_q);
    result_by_q.clear();

    // l , p = (3, 1).
    // q = -9
    result_by_lpq.set(c_Key3(3, 1, -9), common_term_18);
    result_by_q.set(c_Key1(-9), common_term_18);
    // q = -8
    result_by_lpq.set(c_Key3(3, 1, -8), common_term_19);
    result_by_q.set(c_Key1(-8), common_term_19);
    // q = -7
    result_by_lpq.set(c_Key3(3, 1, -7), common_term_20);
    result_by_q.set(c_Key1(-7), common_term_20);
    // q = -6
    result_by_lpq.set(c_Key3(3, 1, -6), common_term_21);
    result_by_q.set(c_Key1(-6), common_term_21);
    // q = -5
    result_by_lpq.set(c_Key3(3, 1, -5), common_term_22);
    result_by_q.set(c_Key1(-5), common_term_22);
    // q = -4
    result_by_lpq.set(c_Key3(3, 1, -4), common_term_23);
    result_by_q.set(c_Key1(-4), common_term_23);
    // q = -3
    result_by_lpq.set(c_Key3(3, 1, -3), common_term_24);
    result_by_q.set(c_Key1(-3), common_term_24);
    // q = -2
    result_by_lpq.set(c_Key3(3, 1, -2), common_term_25);
    result_by_q.set(c_Key1(-2), common_term_25);
    // q = -1
    result_by_lpq.set(c_Key3(3, 1, -1), common_term_26);
    result_by_q.set(c_Key1(-1), common_term_26);
    // q = 0
    result_by_lpq.set(c_Key3(3, 1, 0), common_term_27);
    result_by_q.set(c_Key1(0), common_term_27);
    // q = 1
    result_by_lpq.set(c_Key3(3, 1, 1), common_term_28);
    result_by_q.set(c_Key1(1), common_term_28);
    // q = 2
    result_by_lpq.set(c_Key3(3, 1, 2), common_term_29);
    result_by_q.set(c_Key1(2), common_term_29);
    // q = 3
    result_by_lpq.set(c_Key3(3, 1, 3), common_term_30);
    result_by_q.set(c_Key1(3), common_term_30);
    // q = 4
    result_by_lpq.set(c_Key3(3, 1, 4), common_term_31);
    result_by_q.set(c_Key1(4), common_term_31);
    // q = 5
    result_by_lpq.set(c_Key3(3, 1, 5), common_term_32);
    result_by_q.set(c_Key1(5), common_term_32);
    // q = 6
    result_by_lpq.set(c_Key3(3, 1, 6), common_term_33);
    result_by_q.set(c_Key1(6), common_term_33);
    // q = 7
    result_by_lpq.set(c_Key3(3, 1, 7), common_term_34);
    result_by_q.set(c_Key1(7), common_term_34);
    // q = 8
    result_by_lpq.set(c_Key3(3, 1, 8), common_term_35);
    result_by_q.set(c_Key1(8), common_term_35);
    // q = 9
    result_by_lpq.set(c_Key3(3, 1, 9), common_term_36);
    result_by_q.set(c_Key1(9), common_term_36);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 1), result_by_q);
    result_by_q.clear();

    // l , p = (3, 2).
    // q = -9
    result_by_lpq.set(c_Key3(3, 2, -9), common_term_36);
    result_by_q.set(c_Key1(-9), common_term_36);
    // q = -8
    result_by_lpq.set(c_Key3(3, 2, -8), common_term_35);
    result_by_q.set(c_Key1(-8), common_term_35);
    // q = -7
    result_by_lpq.set(c_Key3(3, 2, -7), common_term_34);
    result_by_q.set(c_Key1(-7), common_term_34);
    // q = -6
    result_by_lpq.set(c_Key3(3, 2, -6), common_term_33);
    result_by_q.set(c_Key1(-6), common_term_33);
    // q = -5
    result_by_lpq.set(c_Key3(3, 2, -5), common_term_32);
    result_by_q.set(c_Key1(-5), common_term_32);
    // q = -4
    result_by_lpq.set(c_Key3(3, 2, -4), common_term_31);
    result_by_q.set(c_Key1(-4), common_term_31);
    // q = -3
    result_by_lpq.set(c_Key3(3, 2, -3), common_term_30);
    result_by_q.set(c_Key1(-3), common_term_30);
    // q = -2
    result_by_lpq.set(c_Key3(3, 2, -2), common_term_29);
    result_by_q.set(c_Key1(-2), common_term_29);
    // q = -1
    result_by_lpq.set(c_Key3(3, 2, -1), common_term_28);
    result_by_q.set(c_Key1(-1), common_term_28);
    // q = 0
    result_by_lpq.set(c_Key3(3, 2, 0), common_term_27);
    result_by_q.set(c_Key1(0), common_term_27);
    // q = 1
    result_by_lpq.set(c_Key3(3, 2, 1), common_term_26);
    result_by_q.set(c_Key1(1), common_term_26);
    // q = 2
    result_by_lpq.set(c_Key3(3, 2, 2), common_term_25);
    result_by_q.set(c_Key1(2), common_term_25);
    // q = 3
    result_by_lpq.set(c_Key3(3, 2, 3), common_term_24);
    result_by_q.set(c_Key1(3), common_term_24);
    // q = 4
    result_by_lpq.set(c_Key3(3, 2, 4), common_term_23);
    result_by_q.set(c_Key1(4), common_term_23);
    // q = 5
    result_by_lpq.set(c_Key3(3, 2, 5), common_term_22);
    result_by_q.set(c_Key1(5), common_term_22);
    // q = 6
    result_by_lpq.set(c_Key3(3, 2, 6), common_term_21);
    result_by_q.set(c_Key1(6), common_term_21);
    // q = 7
    result_by_lpq.set(c_Key3(3, 2, 7), common_term_20);
    result_by_q.set(c_Key1(7), common_term_20);
    // q = 8
    result_by_lpq.set(c_Key3(3, 2, 8), common_term_19);
    result_by_q.set(c_Key1(8), common_term_19);
    // q = 9
    result_by_lpq.set(c_Key3(3, 2, 9), common_term_18);
    result_by_q.set(c_Key1(9), common_term_18);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 2), result_by_q);
    result_by_q.clear();

    // l , p = (3, 3).
    // q = -9
    result_by_lpq.set(c_Key3(3, 3, -9), common_term_17);
    result_by_q.set(c_Key1(-9), common_term_17);
    // q = -8
    result_by_lpq.set(c_Key3(3, 3, -8), common_term_16);
    result_by_q.set(c_Key1(-8), common_term_16);
    // q = -7
    result_by_lpq.set(c_Key3(3, 3, -7), common_term_15);
    result_by_q.set(c_Key1(-7), common_term_15);
    // q = -6
    result_by_lpq.set(c_Key3(3, 3, -6), common_term_14);
    result_by_q.set(c_Key1(-6), common_term_14);
    // q = -5
    result_by_lpq.set(c_Key3(3, 3, -5), common_term_13);
    result_by_q.set(c_Key1(-5), common_term_13);
    // q = -4
    result_by_lpq.set(c_Key3(3, 3, -4), common_term_12);
    result_by_q.set(c_Key1(-4), common_term_12);
    // q = -3
    result_by_lpq.set(c_Key3(3, 3, -3), common_term_11);
    result_by_q.set(c_Key1(-3), common_term_11);
    // q = -2
    result_by_lpq.set(c_Key3(3, 3, -2), common_term_10);
    result_by_q.set(c_Key1(-2), common_term_10);
    // q = -1
    result_by_lpq.set(c_Key3(3, 3, -1), common_term_9);
    result_by_q.set(c_Key1(-1), common_term_9);
    // q = 0
    result_by_lpq.set(c_Key3(3, 3, 0), common_term_8);
    result_by_q.set(c_Key1(0), common_term_8);
    // q = 1
    result_by_lpq.set(c_Key3(3, 3, 1), common_term_7);
    result_by_q.set(c_Key1(1), common_term_7);
    // q = 2
    result_by_lpq.set(c_Key3(3, 3, 2), common_term_6);
    result_by_q.set(c_Key1(2), common_term_6);
    // q = 4
    result_by_lpq.set(c_Key3(3, 3, 4), common_term_5);
    result_by_q.set(c_Key1(4), common_term_5);
    // q = 5
    result_by_lpq.set(c_Key3(3, 3, 5), common_term_4);
    result_by_q.set(c_Key1(5), common_term_4);
    // q = 6
    result_by_lpq.set(c_Key3(3, 3, 6), common_term_3);
    result_by_q.set(c_Key1(6), common_term_3);
    // q = 7
    result_by_lpq.set(c_Key3(3, 3, 7), common_term_2);
    result_by_q.set(c_Key1(7), common_term_2);
    // q = 8
    result_by_lpq.set(c_Key3(3, 3, 8), common_term_1);
    result_by_q.set(c_Key1(8), common_term_1);
    // q = 9
    result_by_lpq.set(c_Key3(3, 3, 9), common_term_0);
    result_by_q.set(c_Key1(9), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 3), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l3_e15(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 3.
    // Functions are taylor expanded and truncated at eccentricity^15.

    // Functions Calculated for l = 3.

    c_IntMap<c_Key3, double> result_by_lpq(114);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(4);
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
    double common_term_0 = 0.26586997874838606*eccentricity_14;
    double common_term_1 = 0.19603324996120135*eccentricity_13;
    double common_term_2 = 0.033219347952843665*eccentricity_14 + 0.14395050779565588*eccentricity_12;
    double common_term_3 = 0.043781732670621559*eccentricity_13 + 0.10507615840949174*eccentricity_11;
    double common_term_4 = 0.058003447433558021*eccentricity_14 + 0.044919827736320453*eccentricity_12 + 0.076018170015311535*eccentricity_10;
    double common_term_5 = 0.046043273133116883*eccentricity_13 + 0.040680803571428571*eccentricity_11 + 0.054241071428571429*eccentricity_9;
    double common_term_6 = 0.033450164237324899*eccentricity_14 + 0.035413253454514492*eccentricity_12 + 0.033639305693342152*eccentricity_10 + 0.037844218905009921*eccentricity_8;
    double common_term_7 = 0.02387419165196943*eccentricity_13 + 0.025573192239858907*eccentricity_11 + 0.025396825396825397*eccentricity_9 + 0.025396825396825397*eccentricity_7;
    double common_term_8 = 0.013446589708328247*eccentricity_14 + 0.015012523106166295*eccentricity_12 + 0.016447121756417411*eccentricity_10 + 0.016950334821428571*eccentricity_8 + 0.0158203125*eccentricity_6;
    double common_term_9 = 0.0064242426330099941*eccentricity_13 + 0.0073681382275132275*eccentricity_11 + 0.0083581349206349206*eccentricity_9 + 0.0090277777777777778*eccentricity_7 + 0.0083333333333333333*eccentricity_5;
    double common_term_10 = 0.001338001090873346*eccentricity_14 + 0.0015721853447969628*eccentricity_12 + 0.0018689029431216931*eccentricity_10 + 0.0022325303819444444*eccentricity_8 + 0.0026041666666666667*eccentricity_6 + 0.0026041666666666667*eccentricity_4;
    double common_term_11 = 0.0064085037893298652*eccentricity_14 + 0.0077965929142381779*eccentricity_12 + 0.0097758540400752315*eccentricity_10 + 0.012771267361111111*eccentricity_8 + 0.017903645833333333*eccentricity_6 + 0.020833333333333333*eccentricity_4 + 0.125*eccentricity_2;
    double common_term_12 = 0.027359853808421517*eccentricity_13 + 0.033930844907407407*eccentricity_11 + 0.042361111111111111*eccentricity_9 + 0.079861111111111111*eccentricity_7 - 0.14583333333333333*eccentricity_5 + 1.25*eccentricity_3 - eccentricity;
    double common_term_13 = 0.055321251226931202*eccentricity_14 + 0.070148124694824219*eccentricity_12 + 0.050625*eccentricity_10 + 0.45977783203125*eccentricity_8 - 1.953125*eccentricity_6 + 6.609375*eccentricity_4 - 6.0*eccentricity_2 + 1.0;
    double common_term_14 = 0.15072176063712522*eccentricity_13 - 0.27019675925925926*eccentricity_11 + 2.8431423611111111*eccentricity_9 - 10.888888888888889*eccentricity_7 + 25.291666666666667*eccentricity_5 - 22.0*eccentricity_3 + 5.0*eccentricity;
    double common_term_15 = 0.56136025454211698*eccentricity_14 - 2.6341507724983982*eccentricity_12 + 13.997727853280527*eccentricity_10 - 43.341200086805556*eccentricity_8 + 79.363606770833333*eccentricity_6 - 63.854166666666667*eccentricity_4 + 15.875*eccentricity_2;
    double common_term_16 = -14.141086774553571*eccentricity_13 + 55.714704241071429*eccentricity_11 - 141.94609375*eccentricity_9 + 217.8*eccentricity_7 - 161.0625*eccentricity_5 + 40.75*eccentricity_3;
    double common_term_17 = -59.040844039720763*eccentricity_14 + 189.07786709644176*eccentricity_12 - 407.87945240162037*eccentricity_10 + 542.69985622829861*eccentricity_8 - 369.51614583333333*eccentricity_6 + 92.221354166666667*eccentricity_4;
    double common_term_18 = 569.60108862893151*eccentricity_13 - 1065.7873677248677*eccentricity_11 + 1257.2000992063492*eccentricity_9 - 791.42222222222222*eccentricity_7 + 191.90833333333333*eccentricity_5;
    double common_term_19 = 1565.9589901169709*eccentricity_14 - 2589.7997371673584*eccentricity_12 + 2750.9964094979422*eccentricity_10 - 1608.5615373883929*eccentricity_8 + 376.0693359375*eccentricity_6;
    double common_term_20 = -5941.1905658952914*eccentricity_13 + 5750.0203493579145*eccentricity_11 - 3137.0445808531746*eccentricity_9 + 704.3968253968254*eccentricity_7;
    double common_term_21 = -13006.093681934419*eccentricity_14 + 11574.082490266097*eccentricity_12 - 5916.4057636994323*eccentricity_10 + 1273.6733050633991*eccentricity_8;
    double common_term_22 = 22574.062998427354*eccentricity_13 - 10853.089419642857*eccentricity_11 + 2238.8170758928571*eccentricity_9;
    double common_term_23 = 42864.585007077403*eccentricity_14 - 19449.455393032119*eccentricity_12 + 3845.1027043953817*eccentricity_10;
    double common_term_24 = -34166.353848185998*eccentricity_13 + 6477.1993013906425*eccentricity_11;
    double common_term_25 = -58993.064721979938*eccentricity_14 + 10733.41011059439*eccentricity_12;
    double common_term_26 = 17537.527175495062*eccentricity_13;
    double common_term_27 = 28306.626833220974*eccentricity_14;
    double common_term_28 = 72.297047055634958*eccentricity_14;
    double common_term_29 = 52.249927777300824*eccentricity_13;
    double common_term_30 = -20.59906285109559*eccentricity_14 + 37.730094974653501*eccentricity_12;
    double common_term_31 = -7.6452967484868527*eccentricity_13 + 27.220283345608866*eccentricity_11;
    double common_term_32 = 18.628644942045212*eccentricity_14 - 0.27764210985852526*eccentricity_12 + 19.618304683140346*eccentricity_10;
    double common_term_33 = 14.4182365476441*eccentricity_13 + 3.594399112654321*eccentricity_11 + 14.123824680335097*eccentricity_9;
    double common_term_34 = 13.666052360456176*eccentricity_14 + 11.837836126869107*eccentricity_12 + 5.3364221643518518*eccentricity_10 + 10.155956353081597*eccentricity_8;
    double common_term_35 = 12.615622209821429*eccentricity_13 + 10.107561383928571*eccentricity_11 + 5.8287946428571429*eccentricity_9 + 7.2933035714285714*eccentricity_7;
    double common_term_36 = 14.554256688472735*eccentricity_14 + 11.50170059361085*eccentricity_12 + 8.80878418210953*eccentricity_10 + 5.6296890500992064*eccentricity_8 + 5.2303602430555556*eccentricity_6;
    double common_term_37 = 13.336312365658069*eccentricity_13 + 10.403629298941799*eccentricity_11 + 7.7297123015873016*eccentricity_9 + 5.0854166666666667*eccentricity_7 + 3.7458333333333333*eccentricity_5;
    double common_term_38 = 15.210991800853184*eccentricity_14 + 12.172791692188808*eccentricity_12 + 9.35205078125*eccentricity_10 + 6.7725341796875*eccentricity_8 + 4.4046875*eccentricity_6 + 2.6796875*eccentricity_4;
    double common_term_39 = 13.996985263723545*eccentricity_13 + 11.064697007275132*eccentricity_11 + 8.358275462962963*eccentricity_9 + 5.8989583333333333*eccentricity_7 + 3.7083333333333333*eccentricity_5 + 1.9166666666666667*eccentricity_3;
    double common_term_40 = 15.875406170103706*eccentricity_14 + 12.836250222019418*eccentricity_12 + 10.014188130696615*eccentricity_10 + 7.4273328993055556*eccentricity_8 + 5.0992838541666667*eccentricity_6 + 3.0625*eccentricity_4 + 1.375*eccentricity_2;
    double common_term_41 = eccentricity*std::pow(1.0 - eccentricity_2, -2.5);
    double common_term_42 = 16.543552944493792*eccentricity_14 + 13.504675437550486*eccentricity_12 + 10.682979600694444*eccentricity_10 + 8.0966050889756944*eccentricity_8 + 5.7690972222222222*eccentricity_6 + 3.734375*eccentricity_4 + 2.0*eccentricity_2 + 1.0;
    double common_term_43 = 15.334919033151455*eccentricity_13 + 12.403249421296296*eccentricity_11 + 9.6975694444444445*eccentricity_9 + 7.234375*eccentricity_7 + 5.1041666666666667*eccentricity_5 + 2.75*eccentricity_3 + 3.0*eccentricity;
    double common_term_44 = 17.218381166287831*eccentricity_14 + 14.180451463971819*eccentricity_12 + 11.364543151855469*eccentricity_10 + 8.72236328125*eccentricity_8 + 6.8759765625*eccentricity_6 + 2.4375*eccentricity_4 + 6.625*eccentricity_2;
    double common_term_45 = 16.008694144758598*eccentricity_13 + 13.12767650462963*eccentricity_11 + 10.008217592592593*eccentricity_9 + 9.8979166666666667*eccentricity_7 - 0.52083333333333333*eccentricity_5 + 12.833333333333333*eccentricity_3;
    double common_term_46 = 17.860908650567012*eccentricity_14 + 15.153994364713235*eccentricity_12 + 10.319979228670635*eccentricity_10 + 16.254814995659722*eccentricity_8 - 9.0182291666666667*eccentricity_6 + 23.0859375*eccentricity_4;
    double common_term_47 = 18.073918805803571*eccentricity_13 + 7.4998046875*eccentricity_11 + 30.470870535714286*eccentricity_9 - 28.121875*eccentricity_7 + 39.5875*eccentricity_5;
    double common_term_48 = 23.830620920427542*eccentricity_14 - 3.7482045852598339*eccentricity_12 + 61.642893897162543*eccentricity_10 - 66.439702690972222*eccentricity_8 + 65.638823784722222*eccentricity_6;
    double common_term_49 = -35.321559227568342*eccentricity_13 + 127.12828689925044*eccentricity_11 - 138.25908978174603*eccentricity_9 + 106.14940476190476*eccentricity_7;
    double common_term_50 = -112.29285589267681*eccentricity_14 + 258.79158464431763*eccentricity_12 - 266.86307198660714*eccentricity_10 + 168.38040684291295*eccentricity_8;
    double common_term_51 = 513.40816221960428*eccentricity_13 - 489.63204606205908*eccentricity_11 + 263.0183876212522*eccentricity_9;
    double common_term_52 = 989.7490990219426*eccentricity_14 - 865.83837034585404*eccentricity_12 + 405.72292616667571*eccentricity_10;
    double common_term_53 = -1488.4829026734984*eccentricity_13 + 619.35198965097403*eccentricity_11;
    double common_term_54 = -2502.1719819251514*eccentricity_14 + 937.15483752107755*eccentricity_12;
    double common_term_55 = 1407.3438943445877*eccentricity_13;
    double common_term_56 = 2099.6289389476391*eccentricity_14;

    c_IntMap<c_Key1, double> result_by_q(29);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (3, 0).
    // q = -14
    result_by_lpq.set(c_Key3(3, 0, -14), common_term_0);
    result_by_q.set(c_Key1(-14), common_term_0);
    // q = -13
    result_by_lpq.set(c_Key3(3, 0, -13), common_term_1);
    result_by_q.set(c_Key1(-13), common_term_1);
    // q = -12
    result_by_lpq.set(c_Key3(3, 0, -12), common_term_2);
    result_by_q.set(c_Key1(-12), common_term_2);
    // q = -11
    result_by_lpq.set(c_Key3(3, 0, -11), common_term_3);
    result_by_q.set(c_Key1(-11), common_term_3);
    // q = -10
    result_by_lpq.set(c_Key3(3, 0, -10), common_term_4);
    result_by_q.set(c_Key1(-10), common_term_4);
    // q = -9
    result_by_lpq.set(c_Key3(3, 0, -9), common_term_5);
    result_by_q.set(c_Key1(-9), common_term_5);
    // q = -8
    result_by_lpq.set(c_Key3(3, 0, -8), common_term_6);
    result_by_q.set(c_Key1(-8), common_term_6);
    // q = -7
    result_by_lpq.set(c_Key3(3, 0, -7), common_term_7);
    result_by_q.set(c_Key1(-7), common_term_7);
    // q = -6
    result_by_lpq.set(c_Key3(3, 0, -6), common_term_8);
    result_by_q.set(c_Key1(-6), common_term_8);
    // q = -5
    result_by_lpq.set(c_Key3(3, 0, -5), common_term_9);
    result_by_q.set(c_Key1(-5), common_term_9);
    // q = -4
    result_by_lpq.set(c_Key3(3, 0, -4), common_term_10);
    result_by_q.set(c_Key1(-4), common_term_10);
    // q = -2
    result_by_lpq.set(c_Key3(3, 0, -2), common_term_11);
    result_by_q.set(c_Key1(-2), common_term_11);
    // q = -1
    result_by_lpq.set(c_Key3(3, 0, -1), common_term_12);
    result_by_q.set(c_Key1(-1), common_term_12);
    // q = 0
    result_by_lpq.set(c_Key3(3, 0, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(3, 0, 1), common_term_14);
    result_by_q.set(c_Key1(1), common_term_14);
    // q = 2
    result_by_lpq.set(c_Key3(3, 0, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // q = 3
    result_by_lpq.set(c_Key3(3, 0, 3), common_term_16);
    result_by_q.set(c_Key1(3), common_term_16);
    // q = 4
    result_by_lpq.set(c_Key3(3, 0, 4), common_term_17);
    result_by_q.set(c_Key1(4), common_term_17);
    // q = 5
    result_by_lpq.set(c_Key3(3, 0, 5), common_term_18);
    result_by_q.set(c_Key1(5), common_term_18);
    // q = 6
    result_by_lpq.set(c_Key3(3, 0, 6), common_term_19);
    result_by_q.set(c_Key1(6), common_term_19);
    // q = 7
    result_by_lpq.set(c_Key3(3, 0, 7), common_term_20);
    result_by_q.set(c_Key1(7), common_term_20);
    // q = 8
    result_by_lpq.set(c_Key3(3, 0, 8), common_term_21);
    result_by_q.set(c_Key1(8), common_term_21);
    // q = 9
    result_by_lpq.set(c_Key3(3, 0, 9), common_term_22);
    result_by_q.set(c_Key1(9), common_term_22);
    // q = 10
    result_by_lpq.set(c_Key3(3, 0, 10), common_term_23);
    result_by_q.set(c_Key1(10), common_term_23);
    // q = 11
    result_by_lpq.set(c_Key3(3, 0, 11), common_term_24);
    result_by_q.set(c_Key1(11), common_term_24);
    // q = 12
    result_by_lpq.set(c_Key3(3, 0, 12), common_term_25);
    result_by_q.set(c_Key1(12), common_term_25);
    // q = 13
    result_by_lpq.set(c_Key3(3, 0, 13), common_term_26);
    result_by_q.set(c_Key1(13), common_term_26);
    // q = 14
    result_by_lpq.set(c_Key3(3, 0, 14), common_term_27);
    result_by_q.set(c_Key1(14), common_term_27);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 0), result_by_q);
    result_by_q.clear();

    // l , p = (3, 1).
    // q = -14
    result_by_lpq.set(c_Key3(3, 1, -14), common_term_28);
    result_by_q.set(c_Key1(-14), common_term_28);
    // q = -13
    result_by_lpq.set(c_Key3(3, 1, -13), common_term_29);
    result_by_q.set(c_Key1(-13), common_term_29);
    // q = -12
    result_by_lpq.set(c_Key3(3, 1, -12), common_term_30);
    result_by_q.set(c_Key1(-12), common_term_30);
    // q = -11
    result_by_lpq.set(c_Key3(3, 1, -11), common_term_31);
    result_by_q.set(c_Key1(-11), common_term_31);
    // q = -10
    result_by_lpq.set(c_Key3(3, 1, -10), common_term_32);
    result_by_q.set(c_Key1(-10), common_term_32);
    // q = -9
    result_by_lpq.set(c_Key3(3, 1, -9), common_term_33);
    result_by_q.set(c_Key1(-9), common_term_33);
    // q = -8
    result_by_lpq.set(c_Key3(3, 1, -8), common_term_34);
    result_by_q.set(c_Key1(-8), common_term_34);
    // q = -7
    result_by_lpq.set(c_Key3(3, 1, -7), common_term_35);
    result_by_q.set(c_Key1(-7), common_term_35);
    // q = -6
    result_by_lpq.set(c_Key3(3, 1, -6), common_term_36);
    result_by_q.set(c_Key1(-6), common_term_36);
    // q = -5
    result_by_lpq.set(c_Key3(3, 1, -5), common_term_37);
    result_by_q.set(c_Key1(-5), common_term_37);
    // q = -4
    result_by_lpq.set(c_Key3(3, 1, -4), common_term_38);
    result_by_q.set(c_Key1(-4), common_term_38);
    // q = -3
    result_by_lpq.set(c_Key3(3, 1, -3), common_term_39);
    result_by_q.set(c_Key1(-3), common_term_39);
    // q = -2
    result_by_lpq.set(c_Key3(3, 1, -2), common_term_40);
    result_by_q.set(c_Key1(-2), common_term_40);
    // q = -1
    result_by_lpq.set(c_Key3(3, 1, -1), common_term_41);
    result_by_q.set(c_Key1(-1), common_term_41);
    // q = 0
    result_by_lpq.set(c_Key3(3, 1, 0), common_term_42);
    result_by_q.set(c_Key1(0), common_term_42);
    // q = 1
    result_by_lpq.set(c_Key3(3, 1, 1), common_term_43);
    result_by_q.set(c_Key1(1), common_term_43);
    // q = 2
    result_by_lpq.set(c_Key3(3, 1, 2), common_term_44);
    result_by_q.set(c_Key1(2), common_term_44);
    // q = 3
    result_by_lpq.set(c_Key3(3, 1, 3), common_term_45);
    result_by_q.set(c_Key1(3), common_term_45);
    // q = 4
    result_by_lpq.set(c_Key3(3, 1, 4), common_term_46);
    result_by_q.set(c_Key1(4), common_term_46);
    // q = 5
    result_by_lpq.set(c_Key3(3, 1, 5), common_term_47);
    result_by_q.set(c_Key1(5), common_term_47);
    // q = 6
    result_by_lpq.set(c_Key3(3, 1, 6), common_term_48);
    result_by_q.set(c_Key1(6), common_term_48);
    // q = 7
    result_by_lpq.set(c_Key3(3, 1, 7), common_term_49);
    result_by_q.set(c_Key1(7), common_term_49);
    // q = 8
    result_by_lpq.set(c_Key3(3, 1, 8), common_term_50);
    result_by_q.set(c_Key1(8), common_term_50);
    // q = 9
    result_by_lpq.set(c_Key3(3, 1, 9), common_term_51);
    result_by_q.set(c_Key1(9), common_term_51);
    // q = 10
    result_by_lpq.set(c_Key3(3, 1, 10), common_term_52);
    result_by_q.set(c_Key1(10), common_term_52);
    // q = 11
    result_by_lpq.set(c_Key3(3, 1, 11), common_term_53);
    result_by_q.set(c_Key1(11), common_term_53);
    // q = 12
    result_by_lpq.set(c_Key3(3, 1, 12), common_term_54);
    result_by_q.set(c_Key1(12), common_term_54);
    // q = 13
    result_by_lpq.set(c_Key3(3, 1, 13), common_term_55);
    result_by_q.set(c_Key1(13), common_term_55);
    // q = 14
    result_by_lpq.set(c_Key3(3, 1, 14), common_term_56);
    result_by_q.set(c_Key1(14), common_term_56);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 1), result_by_q);
    result_by_q.clear();

    // l , p = (3, 2).
    // q = -14
    result_by_lpq.set(c_Key3(3, 2, -14), common_term_56);
    result_by_q.set(c_Key1(-14), common_term_56);
    // q = -13
    result_by_lpq.set(c_Key3(3, 2, -13), common_term_55);
    result_by_q.set(c_Key1(-13), common_term_55);
    // q = -12
    result_by_lpq.set(c_Key3(3, 2, -12), common_term_54);
    result_by_q.set(c_Key1(-12), common_term_54);
    // q = -11
    result_by_lpq.set(c_Key3(3, 2, -11), common_term_53);
    result_by_q.set(c_Key1(-11), common_term_53);
    // q = -10
    result_by_lpq.set(c_Key3(3, 2, -10), common_term_52);
    result_by_q.set(c_Key1(-10), common_term_52);
    // q = -9
    result_by_lpq.set(c_Key3(3, 2, -9), common_term_51);
    result_by_q.set(c_Key1(-9), common_term_51);
    // q = -8
    result_by_lpq.set(c_Key3(3, 2, -8), common_term_50);
    result_by_q.set(c_Key1(-8), common_term_50);
    // q = -7
    result_by_lpq.set(c_Key3(3, 2, -7), common_term_49);
    result_by_q.set(c_Key1(-7), common_term_49);
    // q = -6
    result_by_lpq.set(c_Key3(3, 2, -6), common_term_48);
    result_by_q.set(c_Key1(-6), common_term_48);
    // q = -5
    result_by_lpq.set(c_Key3(3, 2, -5), common_term_47);
    result_by_q.set(c_Key1(-5), common_term_47);
    // q = -4
    result_by_lpq.set(c_Key3(3, 2, -4), common_term_46);
    result_by_q.set(c_Key1(-4), common_term_46);
    // q = -3
    result_by_lpq.set(c_Key3(3, 2, -3), common_term_45);
    result_by_q.set(c_Key1(-3), common_term_45);
    // q = -2
    result_by_lpq.set(c_Key3(3, 2, -2), common_term_44);
    result_by_q.set(c_Key1(-2), common_term_44);
    // q = -1
    result_by_lpq.set(c_Key3(3, 2, -1), common_term_43);
    result_by_q.set(c_Key1(-1), common_term_43);
    // q = 0
    result_by_lpq.set(c_Key3(3, 2, 0), common_term_42);
    result_by_q.set(c_Key1(0), common_term_42);
    // q = 1
    result_by_lpq.set(c_Key3(3, 2, 1), common_term_41);
    result_by_q.set(c_Key1(1), common_term_41);
    // q = 2
    result_by_lpq.set(c_Key3(3, 2, 2), common_term_40);
    result_by_q.set(c_Key1(2), common_term_40);
    // q = 3
    result_by_lpq.set(c_Key3(3, 2, 3), common_term_39);
    result_by_q.set(c_Key1(3), common_term_39);
    // q = 4
    result_by_lpq.set(c_Key3(3, 2, 4), common_term_38);
    result_by_q.set(c_Key1(4), common_term_38);
    // q = 5
    result_by_lpq.set(c_Key3(3, 2, 5), common_term_37);
    result_by_q.set(c_Key1(5), common_term_37);
    // q = 6
    result_by_lpq.set(c_Key3(3, 2, 6), common_term_36);
    result_by_q.set(c_Key1(6), common_term_36);
    // q = 7
    result_by_lpq.set(c_Key3(3, 2, 7), common_term_35);
    result_by_q.set(c_Key1(7), common_term_35);
    // q = 8
    result_by_lpq.set(c_Key3(3, 2, 8), common_term_34);
    result_by_q.set(c_Key1(8), common_term_34);
    // q = 9
    result_by_lpq.set(c_Key3(3, 2, 9), common_term_33);
    result_by_q.set(c_Key1(9), common_term_33);
    // q = 10
    result_by_lpq.set(c_Key3(3, 2, 10), common_term_32);
    result_by_q.set(c_Key1(10), common_term_32);
    // q = 11
    result_by_lpq.set(c_Key3(3, 2, 11), common_term_31);
    result_by_q.set(c_Key1(11), common_term_31);
    // q = 12
    result_by_lpq.set(c_Key3(3, 2, 12), common_term_30);
    result_by_q.set(c_Key1(12), common_term_30);
    // q = 13
    result_by_lpq.set(c_Key3(3, 2, 13), common_term_29);
    result_by_q.set(c_Key1(13), common_term_29);
    // q = 14
    result_by_lpq.set(c_Key3(3, 2, 14), common_term_28);
    result_by_q.set(c_Key1(14), common_term_28);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 2), result_by_q);
    result_by_q.clear();

    // l , p = (3, 3).
    // q = -14
    result_by_lpq.set(c_Key3(3, 3, -14), common_term_27);
    result_by_q.set(c_Key1(-14), common_term_27);
    // q = -13
    result_by_lpq.set(c_Key3(3, 3, -13), common_term_26);
    result_by_q.set(c_Key1(-13), common_term_26);
    // q = -12
    result_by_lpq.set(c_Key3(3, 3, -12), common_term_25);
    result_by_q.set(c_Key1(-12), common_term_25);
    // q = -11
    result_by_lpq.set(c_Key3(3, 3, -11), common_term_24);
    result_by_q.set(c_Key1(-11), common_term_24);
    // q = -10
    result_by_lpq.set(c_Key3(3, 3, -10), common_term_23);
    result_by_q.set(c_Key1(-10), common_term_23);
    // q = -9
    result_by_lpq.set(c_Key3(3, 3, -9), common_term_22);
    result_by_q.set(c_Key1(-9), common_term_22);
    // q = -8
    result_by_lpq.set(c_Key3(3, 3, -8), common_term_21);
    result_by_q.set(c_Key1(-8), common_term_21);
    // q = -7
    result_by_lpq.set(c_Key3(3, 3, -7), common_term_20);
    result_by_q.set(c_Key1(-7), common_term_20);
    // q = -6
    result_by_lpq.set(c_Key3(3, 3, -6), common_term_19);
    result_by_q.set(c_Key1(-6), common_term_19);
    // q = -5
    result_by_lpq.set(c_Key3(3, 3, -5), common_term_18);
    result_by_q.set(c_Key1(-5), common_term_18);
    // q = -4
    result_by_lpq.set(c_Key3(3, 3, -4), common_term_17);
    result_by_q.set(c_Key1(-4), common_term_17);
    // q = -3
    result_by_lpq.set(c_Key3(3, 3, -3), common_term_16);
    result_by_q.set(c_Key1(-3), common_term_16);
    // q = -2
    result_by_lpq.set(c_Key3(3, 3, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(3, 3, -1), common_term_14);
    result_by_q.set(c_Key1(-1), common_term_14);
    // q = 0
    result_by_lpq.set(c_Key3(3, 3, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(3, 3, 1), common_term_12);
    result_by_q.set(c_Key1(1), common_term_12);
    // q = 2
    result_by_lpq.set(c_Key3(3, 3, 2), common_term_11);
    result_by_q.set(c_Key1(2), common_term_11);
    // q = 4
    result_by_lpq.set(c_Key3(3, 3, 4), common_term_10);
    result_by_q.set(c_Key1(4), common_term_10);
    // q = 5
    result_by_lpq.set(c_Key3(3, 3, 5), common_term_9);
    result_by_q.set(c_Key1(5), common_term_9);
    // q = 6
    result_by_lpq.set(c_Key3(3, 3, 6), common_term_8);
    result_by_q.set(c_Key1(6), common_term_8);
    // q = 7
    result_by_lpq.set(c_Key3(3, 3, 7), common_term_7);
    result_by_q.set(c_Key1(7), common_term_7);
    // q = 8
    result_by_lpq.set(c_Key3(3, 3, 8), common_term_6);
    result_by_q.set(c_Key1(8), common_term_6);
    // q = 9
    result_by_lpq.set(c_Key3(3, 3, 9), common_term_5);
    result_by_q.set(c_Key1(9), common_term_5);
    // q = 10
    result_by_lpq.set(c_Key3(3, 3, 10), common_term_4);
    result_by_q.set(c_Key1(10), common_term_4);
    // q = 11
    result_by_lpq.set(c_Key3(3, 3, 11), common_term_3);
    result_by_q.set(c_Key1(11), common_term_3);
    // q = 12
    result_by_lpq.set(c_Key3(3, 3, 12), common_term_2);
    result_by_q.set(c_Key1(12), common_term_2);
    // q = 13
    result_by_lpq.set(c_Key3(3, 3, 13), common_term_1);
    result_by_q.set(c_Key1(13), common_term_1);
    // q = 14
    result_by_lpq.set(c_Key3(3, 3, 14), common_term_0);
    result_by_q.set(c_Key1(14), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 3), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l3_e20(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 3.
    // Functions are taylor expanded and truncated at eccentricity^20.

    // Functions Calculated for l = 3.

    c_IntMap<c_Key3, double> result_by_lpq(154);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(4);
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
    double common_term_0 = 1.1847183946702751*eccentricity_19;
    double common_term_1 = 0.88056561094154073*eccentricity_18;
    double common_term_2 = -0.52685810105034298*eccentricity_19 + 0.65403074613146026*eccentricity_17;
    double common_term_3 = -0.28545995377350605*eccentricity_18 + 0.48528192141496028*eccentricity_16;
    double common_term_4 = 0.24984004231063055*eccentricity_19 - 0.13483430854859426*eccentricity_17 + 0.3595581561295847*eccentricity_15;
    double common_term_5 = 0.16973444085980949*eccentricity_18 - 0.044311663124731011*eccentricity_16 + 0.26586997874838606*eccentricity_14;
    double common_term_6 = 0.075676030149952853*eccentricity_19 + 0.12193734893419965*eccentricity_17 + 0.0070011874986143339*eccentricity_15 + 0.19603324996120135*eccentricity_13;
    double common_term_7 = 0.071253524016709611*eccentricity_18 + 0.09231716116359456*eccentricity_16 + 0.033219347952843665*eccentricity_14 + 0.14395050779565588*eccentricity_12;
    double common_term_8 = 0.061730171585550774*eccentricity_19 + 0.063114655283967453*eccentricity_17 + 0.072576641465530354*eccentricity_15 + 0.043781732670621559*eccentricity_13 + 0.10507615840949174*eccentricity_11;
    double common_term_9 = 0.051100214156006684*eccentricity_18 + 0.053483731243927528*eccentricity_16 + 0.058003447433558021*eccentricity_14 + 0.044919827736320453*eccentricity_12 + 0.076018170015311535*eccentricity_10;
    double common_term_10 = 0.037907233069316286*eccentricity_19 + 0.040869272329233267*eccentricity_17 + 0.043423675933441558*eccentricity_15 + 0.046043273133116883*eccentricity_13 + 0.040680803571428571*eccentricity_11 + 0.054241071428571429*eccentricity_9;
    double common_term_11 = 0.028514431715364252*eccentricity_18 + 0.031058564370842929*eccentricity_16 + 0.033450164237324899*eccentricity_14 + 0.035413253454514492*eccentricity_12 + 0.033639305693342152*eccentricity_10 + 0.037844218905009921*eccentricity_8;
    double common_term_12 = 0.017889208963475939*eccentricity_19 + 0.019774755936792974*eccentricity_17 + 0.021821522099299877*eccentricity_15 + 0.02387419165196943*eccentricity_13 + 0.025573192239858907*eccentricity_11 + 0.025396825396825397*eccentricity_9 + 0.025396825396825397*eccentricity_7;
    double common_term_13 = 0.010700422771783038*eccentricity_18 + 0.011987028414552862*eccentricity_16 + 0.013446589708328247*eccentricity_14 + 0.015012523106166295*eccentricity_12 + 0.016447121756417411*eccentricity_10 + 0.016950334821428571*eccentricity_8 + 0.0158203125*eccentricity_6;
    double common_term_14 = 0.0043601410536988199*eccentricity_19 + 0.0049278861055300581*eccentricity_17 + 0.0056099186829438566*eccentricity_15 + 0.0064242426330099941*eccentricity_13 + 0.0073681382275132275*eccentricity_11 + 0.0083581349206349206*eccentricity_9 + 0.0090277777777777778*eccentricity_7 + 0.0083333333333333333*eccentricity_5;
    double common_term_15 = 0.0010050875644115889*eccentricity_18 + 0.0011529448717487996*eccentricity_16 + 0.001338001090873346*eccentricity_14 + 0.0015721853447969628*eccentricity_12 + 0.0018689029431216931*eccentricity_10 + 0.0022325303819444444*eccentricity_8 + 0.0026041666666666667*eccentricity_6 + 0.0026041666666666667*eccentricity_4;
    double common_term_16 = 0.0046163817770085458*eccentricity_18 + 0.0053900735212438288*eccentricity_16 + 0.0064085037893298652*eccentricity_14 + 0.0077965929142381779*eccentricity_12 + 0.0097758540400752315*eccentricity_10 + 0.012771267361111111*eccentricity_8 + 0.017903645833333333*eccentricity_6 + 0.020833333333333333*eccentricity_4 + 0.125*eccentricity_2;
    double common_term_17 = 0.016552703418200521*eccentricity_19 + 0.019218978038183509*eccentricity_17 + 0.022691984707656211*eccentricity_15 + 0.027359853808421517*eccentricity_13 + 0.033930844907407407*eccentricity_11 + 0.042361111111111111*eccentricity_9 + 0.079861111111111111*eccentricity_7 - 0.14583333333333333*eccentricity_5 + 1.25*eccentricity_3 - eccentricity;
    double common_term_18 = 0.039582024999723142*eccentricity_18 + 0.046407595599944494*eccentricity_16 + 0.055321251226931202*eccentricity_14 + 0.070148124694824219*eccentricity_12 + 0.050625*eccentricity_10 + 0.45977783203125*eccentricity_8 - 1.953125*eccentricity_6 + 6.609375*eccentricity_4 - 6.0*eccentricity_2 + 1.0;
    double common_term_19 = 0.065097676633971663*eccentricity_19 + 0.075984495105349882*eccentricity_17 + 0.08648675673973293*eccentricity_15 + 0.15072176063712522*eccentricity_13 - 0.27019675925925926*eccentricity_11 + 2.8431423611111111*eccentricity_9 - 10.888888888888889*eccentricity_7 + 25.291666666666667*eccentricity_5 - 22.0*eccentricity_3 + 5.0*eccentricity;
    double common_term_20 = 0.11316769354825195*eccentricity_18 + 0.084653751537879995*eccentricity_16 + 0.56136025454211698*eccentricity_14 - 2.6341507724983982*eccentricity_12 + 13.997727853280527*eccentricity_10 - 43.341200086805556*eccentricity_8 + 79.363606770833333*eccentricity_6 - 63.854166666666667*eccentricity_4 + 15.875*eccentricity_2;
    double common_term_21 = 0.18779875814406018*eccentricity_19 - 0.20534810467155612*eccentricity_17 + 2.8793540736607143*eccentricity_15 - 14.141086774553571*eccentricity_13 + 55.714704241071429*eccentricity_11 - 141.94609375*eccentricity_9 + 217.8*eccentricity_7 - 161.0625*eccentricity_5 + 40.75*eccentricity_3;
    double common_term_22 = -2.1532877969586101*eccentricity_18 + 13.827958417275356*eccentricity_16 - 59.040844039720763*eccentricity_14 + 189.07786709644176*eccentricity_12 - 407.87945240162037*eccentricity_10 + 542.69985622829861*eccentricity_8 - 369.51614583333333*eccentricity_6 + 92.221354166666667*eccentricity_4;
    double common_term_23 = -11.709313195699701*eccentricity_19 + 57.4814977707103*eccentricity_17 - 209.85188396531452*eccentricity_15 + 569.60108862893151*eccentricity_13 - 1065.7873677248677*eccentricity_11 + 1257.2000992063492*eccentricity_9 - 791.42222222222222*eccentricity_7 + 191.90833333333333*eccentricity_5;
    double common_term_24 = 209.83209740737523*eccentricity_18 - 663.391721942866*eccentricity_16 + 1565.9589901169709*eccentricity_14 - 2589.7997371673584*eccentricity_12 + 2750.9964094979422*eccentricity_10 - 1608.5615373883929*eccentricity_8 + 376.0693359375*eccentricity_6;
    double common_term_25 = 688.90625974122327*eccentricity_19 - 1915.7039755603983*eccentricity_17 + 4005.1746757604776*eccentricity_15 - 5941.1905658952914*eccentricity_13 + 5750.0203493579145*eccentricity_11 - 3137.0445808531746*eccentricity_9 + 704.3968253968254*eccentricity_7;
    double common_term_26 = -5145.9077733327566*eccentricity_18 + 9662.1137583953675*eccentricity_16 - 13006.093681934419*eccentricity_14 + 11574.082490266097*eccentricity_12 - 5916.4057636994323*eccentricity_10 + 1273.6733050633991*eccentricity_8;
    double common_term_27 = -13026.711263222304*eccentricity_19 + 22209.708993651002*eccentricity_17 - 27385.661229707792*eccentricity_15 + 22574.062998427354*eccentricity_13 - 10853.089419642857*eccentricity_11 + 2238.8170758928571*eccentricity_9;
    double common_term_28 = 49020.398714479228*eccentricity_18 - 55798.065830308909*eccentricity_16 + 42864.585007077403*eccentricity_14 - 19449.455393032119*eccentricity_12 + 3845.1027043953817*eccentricity_10;
    double common_term_29 = 104512.53953690455*eccentricity_19 - 110529.05207998121*eccentricity_17 + 79538.469868377929*eccentricity_15 - 34166.353848185998*eccentricity_13 + 6477.1993013906425*eccentricity_11;
    double common_term_30 = -213659.01335272217*eccentricity_18 + 144660.22113849366*eccentricity_16 - 58993.064721979938*eccentricity_14 + 10733.41011059439*eccentricity_12;
    double common_term_31 = -404269.3766396613*eccentricity_19 + 258509.75095400226*eccentricity_17 - 100337.41314053047*eccentricity_15 + 17537.527175495062*eccentricity_13;
    double common_term_32 = 454818.94950304826*eccentricity_18 - 168407.7811542272*eccentricity_16 + 28306.626833220974*eccentricity_14;
    double common_term_33 = 789165.06342623266*eccentricity_19 - 279347.37844541265*eccentricity_17 + 45201.783186148255*eccentricity_15;
    double common_term_34 = -458514.09449079658*eccentricity_18 + 71501.308430063661*eccentricity_16;
    double common_term_35 = -745498.0979588166*eccentricity_19 + 112154.53337980564*eccentricity_17;
    double common_term_36 = 174600.90380983502*eccentricity_18;
    double common_term_37 = 269978.01123975852*eccentricity_19;
    double common_term_38 = 362.92431968918894*eccentricity_19;
    double common_term_39 = 263.14619687930975*eccentricity_18;
    double common_term_40 = -353.57095853830246*eccentricity_19 + 190.69467336139785*eccentricity_17;
    double common_term_41 = -220.18646804326093*eccentricity_18 + 138.10874004193418*eccentricity_16;
    double common_term_42 = 164.71567477603604*eccentricity_19 - 133.3148816045511*eccentricity_17 + 99.959313551846628*eccentricity_15;
    double common_term_43 = 99.107628309412249*eccentricity_18 - 77.516317519644312*eccentricity_16 + 72.297047055634958*eccentricity_14;
    double common_term_44 = -2.353573999730738*eccentricity_19 + 60.699544446346065*eccentricity_17 - 42.305080914481724*eccentricity_15 + 52.249927777300824*eccentricity_13;
    double common_term_45 = 8.6230200074838971*eccentricity_18 + 38.529170036374742*eccentricity_16 - 20.59906285109559*eccentricity_14 + 37.730094974653501*eccentricity_12;
    double common_term_46 = 22.792027039469618*eccentricity_19 + 13.211580956444961*eccentricity_17 + 25.862073823940415*eccentricity_15 - 7.6452967484868527*eccentricity_13 + 27.220283345608866*eccentricity_11;
    double common_term_47 = 20.411065935038398*eccentricity_18 + 14.598081611907566*eccentricity_16 + 18.628644942045212*eccentricity_14 - 0.27764210985852526*eccentricity_12 + 19.618304683140346*eccentricity_10;
    double common_term_48 = 21.993871895110865*eccentricity_19 + 18.658622098263935*eccentricity_17 + 14.464730972888887*eccentricity_15 + 14.4182365476441*eccentricity_13 + 3.594399112654321*eccentricity_11 + 14.123824680335097*eccentricity_9;
    double common_term_49 = 20.559212238256856*eccentricity_18 + 17.180078433754987*eccentricity_16 + 13.666052360456176*eccentricity_14 + 11.837836126869107*eccentricity_12 + 5.3364221643518518*eccentricity_10 + 10.155956353081597*eccentricity_8;
    double common_term_50 = 22.676680482481386*eccentricity_19 + 19.154007210075081*eccentricity_17 + 15.829914677607549*eccentricity_15 + 12.615622209821429*eccentricity_13 + 10.107561383928571*eccentricity_11 + 5.8287946428571429*eccentricity_9 + 7.2933035714285714*eccentricity_7;
    double common_term_51 = 21.222382027929206*eccentricity_18 + 17.791540311252145*eccentricity_16 + 14.554256688472735*eccentricity_14 + 11.50170059361085*eccentricity_12 + 8.80878418210953*eccentricity_10 + 5.6296890500992064*eccentricity_8 + 5.2303602430555556*eccentricity_6;
    double common_term_52 = 23.336969753760689*eccentricity_19 + 19.813800440993656*eccentricity_17 + 16.476448718699478*eccentricity_15 + 13.336312365658069*eccentricity_13 + 10.403629298941799*eccentricity_11 + 7.7297123015873016*eccentricity_9 + 5.0854166666666667*eccentricity_7 + 3.7458333333333333*eccentricity_5;
    double common_term_53 = 21.884243104901407*eccentricity_18 + 18.452171295911074*eccentricity_16 + 15.210991800853184*eccentricity_14 + 12.172791692188808*eccentricity_12 + 9.35205078125*eccentricity_10 + 6.7725341796875*eccentricity_8 + 4.4046875*eccentricity_6 + 2.6796875*eccentricity_4;
    double common_term_54 = 24.001106127735863*eccentricity_19 + 20.477240133276819*eccentricity_17 + 17.138885871458119*eccentricity_15 + 13.996985263723545*eccentricity_13 + 11.064697007275132*eccentricity_11 + 8.358275462962963*eccentricity_9 + 5.8989583333333333*eccentricity_7 + 3.7083333333333333*eccentricity_5 + 1.9166666666666667*eccentricity_3;
    double common_term_55 = 22.549639233903079*eccentricity_18 + 19.117176184740987*eccentricity_16 + 15.875406170103706*eccentricity_14 + 12.836250222019418*eccentricity_12 + 10.014188130696615*eccentricity_10 + 7.4273328993055556*eccentricity_8 + 5.0992838541666667*eccentricity_6 + 3.0625*eccentricity_4 + 1.375*eccentricity_2;
    double common_term_56 = eccentricity*std::pow(1.0 - eccentricity_2, -2.5);
    double common_term_57 = 23.217402262204715*eccentricity_18 + 19.785108095855694*eccentricity_16 + 16.543552944493792*eccentricity_14 + 13.504675437550486*eccentricity_12 + 10.682979600694444*eccentricity_10 + 8.0966050889756944*eccentricity_8 + 5.7690972222222222*eccentricity_6 + 3.734375*eccentricity_4 + 2.0*eccentricity_2 + 1.0;
    double common_term_58 = 25.337603432757685*eccentricity_19 + 21.814138332486657*eccentricity_17 + 18.476260676615843*eccentricity_15 + 15.334919033151455*eccentricity_13 + 12.403249421296296*eccentricity_11 + 9.6975694444444445*eccentricity_9 + 7.234375*eccentricity_7 + 5.1041666666666667*eccentricity_5 + 2.75*eccentricity_3 + 3.0*eccentricity;
    double common_term_59 = 23.890366827521492*eccentricity_18 + 20.458915217810748*eccentricity_16 + 17.218381166287831*eccentricity_14 + 14.180451463971819*eccentricity_12 + 11.364543151855469*eccentricity_10 + 8.72236328125*eccentricity_8 + 6.8759765625*eccentricity_6 + 2.4375*eccentricity_4 + 6.625*eccentricity_2;
    double common_term_60 = 26.01228331556653*eccentricity_19 + 22.489724990967323*eccentricity_17 + 19.153175871002658*eccentricity_15 + 16.008694144758598*eccentricity_13 + 13.12767650462963*eccentricity_11 + 10.008217592592593*eccentricity_9 + 9.8979166666666667*eccentricity_7 - 0.52083333333333333*eccentricity_5 + 12.833333333333333*eccentricity_3;
    double common_term_61 = 24.566969745945764*eccentricity_18 + 21.140179657442457*eccentricity_16 + 17.860908650567012*eccentricity_14 + 15.153994364713235*eccentricity_12 + 10.319979228670635*eccentricity_10 + 16.254814995659722*eccentricity_8 - 9.0182291666666667*eccentricity_6 + 23.0859375*eccentricity_4;
    double common_term_62 = 26.687327121939863*eccentricity_19 + 23.196034116274351*eccentricity_17 + 19.609353585379464*eccentricity_15 + 18.073918805803571*eccentricity_13 + 7.4998046875*eccentricity_11 + 30.470870535714286*eccentricity_9 - 28.121875*eccentricity_7 + 39.5875*eccentricity_5;
    double common_term_63 = 25.410670672293916*eccentricity_18 + 20.759988900410696*eccentricity_16 + 23.830620920427542*eccentricity_14 - 3.7482045852598339*eccentricity_12 + 61.642893897162543*eccentricity_10 - 66.439702690972222*eccentricity_8 + 65.638823784722222*eccentricity_6;
    double common_term_64 = 28.152227062214286*eccentricity_19 + 19.69212737336542*eccentricity_17 + 37.672438663566051*eccentricity_15 - 35.321559227568342*eccentricity_13 + 127.12828689925044*eccentricity_11 - 138.25908978174603*eccentricity_9 + 106.14940476190476*eccentricity_7;
    double common_term_65 = 11.703295808164241*eccentricity_18 + 72.505384321380165*eccentricity_16 - 112.29285589267681*eccentricity_14 + 258.79158464431763*eccentricity_12 - 266.86307198660714*eccentricity_10 + 168.38040684291295*eccentricity_8;
    double common_term_66 = -15.639772736933843*eccentricity_19 + 157.9137975599858*eccentricity_17 - 284.9986189970886*eccentricity_15 + 513.40816221960428*eccentricity_13 - 489.63204606205908*eccentricity_11 + 263.0183876212522*eccentricity_9;
    double common_term_67 = 358.10881670547177*eccentricity_18 - 650.71611562409349*eccentricity_16 + 989.7490990219426*eccentricity_14 - 865.83837034585404*eccentricity_12 + 405.72292616667571*eccentricity_10;
    double common_term_68 = 806.49216793896847*eccentricity_19 - 1391.7311202162793*eccentricity_17 + 1856.3021100237653*eccentricity_15 - 1488.4829026734984*eccentricity_13 + 619.35198965097403*eccentricity_11;
    double common_term_69 = -2841.0286440435211*eccentricity_18 + 3395.7837828876213*eccentricity_16 - 2502.1719819251514*eccentricity_14 + 937.15483752107755*eccentricity_12;
    double common_term_70 = -5593.9623804454039*eccentricity_19 + 6075.96479767522*eccentricity_17 - 4129.9858741006514*eccentricity_15 + 1407.3438943445877*eccentricity_13;
    double common_term_71 = 10661.456712640627*eccentricity_18 - 6713.6831498565191*eccentricity_16 + 2099.6289389476391*eccentricity_14;
    double common_term_72 = 18388.884049270001*eccentricity_19 - 10773.619165930142*eccentricity_17 + 3114.5396917942013*eccentricity_15;
    double common_term_73 = -17097.71846969076*eccentricity_18 + 4596.7052660038601*eccentricity_16;
    double common_term_74 = -26873.145089983613*eccentricity_19 + 6753.7414328156405*eccentricity_17;
    double common_term_75 = 9883.0759593790787*eccentricity_18;
    double common_term_76 = 14409.997880518479*eccentricity_19;

    c_IntMap<c_Key1, double> result_by_q(39);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (3, 0).
    // q = -19
    result_by_lpq.set(c_Key3(3, 0, -19), common_term_0);
    result_by_q.set(c_Key1(-19), common_term_0);
    // q = -18
    result_by_lpq.set(c_Key3(3, 0, -18), common_term_1);
    result_by_q.set(c_Key1(-18), common_term_1);
    // q = -17
    result_by_lpq.set(c_Key3(3, 0, -17), common_term_2);
    result_by_q.set(c_Key1(-17), common_term_2);
    // q = -16
    result_by_lpq.set(c_Key3(3, 0, -16), common_term_3);
    result_by_q.set(c_Key1(-16), common_term_3);
    // q = -15
    result_by_lpq.set(c_Key3(3, 0, -15), common_term_4);
    result_by_q.set(c_Key1(-15), common_term_4);
    // q = -14
    result_by_lpq.set(c_Key3(3, 0, -14), common_term_5);
    result_by_q.set(c_Key1(-14), common_term_5);
    // q = -13
    result_by_lpq.set(c_Key3(3, 0, -13), common_term_6);
    result_by_q.set(c_Key1(-13), common_term_6);
    // q = -12
    result_by_lpq.set(c_Key3(3, 0, -12), common_term_7);
    result_by_q.set(c_Key1(-12), common_term_7);
    // q = -11
    result_by_lpq.set(c_Key3(3, 0, -11), common_term_8);
    result_by_q.set(c_Key1(-11), common_term_8);
    // q = -10
    result_by_lpq.set(c_Key3(3, 0, -10), common_term_9);
    result_by_q.set(c_Key1(-10), common_term_9);
    // q = -9
    result_by_lpq.set(c_Key3(3, 0, -9), common_term_10);
    result_by_q.set(c_Key1(-9), common_term_10);
    // q = -8
    result_by_lpq.set(c_Key3(3, 0, -8), common_term_11);
    result_by_q.set(c_Key1(-8), common_term_11);
    // q = -7
    result_by_lpq.set(c_Key3(3, 0, -7), common_term_12);
    result_by_q.set(c_Key1(-7), common_term_12);
    // q = -6
    result_by_lpq.set(c_Key3(3, 0, -6), common_term_13);
    result_by_q.set(c_Key1(-6), common_term_13);
    // q = -5
    result_by_lpq.set(c_Key3(3, 0, -5), common_term_14);
    result_by_q.set(c_Key1(-5), common_term_14);
    // q = -4
    result_by_lpq.set(c_Key3(3, 0, -4), common_term_15);
    result_by_q.set(c_Key1(-4), common_term_15);
    // q = -2
    result_by_lpq.set(c_Key3(3, 0, -2), common_term_16);
    result_by_q.set(c_Key1(-2), common_term_16);
    // q = -1
    result_by_lpq.set(c_Key3(3, 0, -1), common_term_17);
    result_by_q.set(c_Key1(-1), common_term_17);
    // q = 0
    result_by_lpq.set(c_Key3(3, 0, 0), common_term_18);
    result_by_q.set(c_Key1(0), common_term_18);
    // q = 1
    result_by_lpq.set(c_Key3(3, 0, 1), common_term_19);
    result_by_q.set(c_Key1(1), common_term_19);
    // q = 2
    result_by_lpq.set(c_Key3(3, 0, 2), common_term_20);
    result_by_q.set(c_Key1(2), common_term_20);
    // q = 3
    result_by_lpq.set(c_Key3(3, 0, 3), common_term_21);
    result_by_q.set(c_Key1(3), common_term_21);
    // q = 4
    result_by_lpq.set(c_Key3(3, 0, 4), common_term_22);
    result_by_q.set(c_Key1(4), common_term_22);
    // q = 5
    result_by_lpq.set(c_Key3(3, 0, 5), common_term_23);
    result_by_q.set(c_Key1(5), common_term_23);
    // q = 6
    result_by_lpq.set(c_Key3(3, 0, 6), common_term_24);
    result_by_q.set(c_Key1(6), common_term_24);
    // q = 7
    result_by_lpq.set(c_Key3(3, 0, 7), common_term_25);
    result_by_q.set(c_Key1(7), common_term_25);
    // q = 8
    result_by_lpq.set(c_Key3(3, 0, 8), common_term_26);
    result_by_q.set(c_Key1(8), common_term_26);
    // q = 9
    result_by_lpq.set(c_Key3(3, 0, 9), common_term_27);
    result_by_q.set(c_Key1(9), common_term_27);
    // q = 10
    result_by_lpq.set(c_Key3(3, 0, 10), common_term_28);
    result_by_q.set(c_Key1(10), common_term_28);
    // q = 11
    result_by_lpq.set(c_Key3(3, 0, 11), common_term_29);
    result_by_q.set(c_Key1(11), common_term_29);
    // q = 12
    result_by_lpq.set(c_Key3(3, 0, 12), common_term_30);
    result_by_q.set(c_Key1(12), common_term_30);
    // q = 13
    result_by_lpq.set(c_Key3(3, 0, 13), common_term_31);
    result_by_q.set(c_Key1(13), common_term_31);
    // q = 14
    result_by_lpq.set(c_Key3(3, 0, 14), common_term_32);
    result_by_q.set(c_Key1(14), common_term_32);
    // q = 15
    result_by_lpq.set(c_Key3(3, 0, 15), common_term_33);
    result_by_q.set(c_Key1(15), common_term_33);
    // q = 16
    result_by_lpq.set(c_Key3(3, 0, 16), common_term_34);
    result_by_q.set(c_Key1(16), common_term_34);
    // q = 17
    result_by_lpq.set(c_Key3(3, 0, 17), common_term_35);
    result_by_q.set(c_Key1(17), common_term_35);
    // q = 18
    result_by_lpq.set(c_Key3(3, 0, 18), common_term_36);
    result_by_q.set(c_Key1(18), common_term_36);
    // q = 19
    result_by_lpq.set(c_Key3(3, 0, 19), common_term_37);
    result_by_q.set(c_Key1(19), common_term_37);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 0), result_by_q);
    result_by_q.clear();

    // l , p = (3, 1).
    // q = -19
    result_by_lpq.set(c_Key3(3, 1, -19), common_term_38);
    result_by_q.set(c_Key1(-19), common_term_38);
    // q = -18
    result_by_lpq.set(c_Key3(3, 1, -18), common_term_39);
    result_by_q.set(c_Key1(-18), common_term_39);
    // q = -17
    result_by_lpq.set(c_Key3(3, 1, -17), common_term_40);
    result_by_q.set(c_Key1(-17), common_term_40);
    // q = -16
    result_by_lpq.set(c_Key3(3, 1, -16), common_term_41);
    result_by_q.set(c_Key1(-16), common_term_41);
    // q = -15
    result_by_lpq.set(c_Key3(3, 1, -15), common_term_42);
    result_by_q.set(c_Key1(-15), common_term_42);
    // q = -14
    result_by_lpq.set(c_Key3(3, 1, -14), common_term_43);
    result_by_q.set(c_Key1(-14), common_term_43);
    // q = -13
    result_by_lpq.set(c_Key3(3, 1, -13), common_term_44);
    result_by_q.set(c_Key1(-13), common_term_44);
    // q = -12
    result_by_lpq.set(c_Key3(3, 1, -12), common_term_45);
    result_by_q.set(c_Key1(-12), common_term_45);
    // q = -11
    result_by_lpq.set(c_Key3(3, 1, -11), common_term_46);
    result_by_q.set(c_Key1(-11), common_term_46);
    // q = -10
    result_by_lpq.set(c_Key3(3, 1, -10), common_term_47);
    result_by_q.set(c_Key1(-10), common_term_47);
    // q = -9
    result_by_lpq.set(c_Key3(3, 1, -9), common_term_48);
    result_by_q.set(c_Key1(-9), common_term_48);
    // q = -8
    result_by_lpq.set(c_Key3(3, 1, -8), common_term_49);
    result_by_q.set(c_Key1(-8), common_term_49);
    // q = -7
    result_by_lpq.set(c_Key3(3, 1, -7), common_term_50);
    result_by_q.set(c_Key1(-7), common_term_50);
    // q = -6
    result_by_lpq.set(c_Key3(3, 1, -6), common_term_51);
    result_by_q.set(c_Key1(-6), common_term_51);
    // q = -5
    result_by_lpq.set(c_Key3(3, 1, -5), common_term_52);
    result_by_q.set(c_Key1(-5), common_term_52);
    // q = -4
    result_by_lpq.set(c_Key3(3, 1, -4), common_term_53);
    result_by_q.set(c_Key1(-4), common_term_53);
    // q = -3
    result_by_lpq.set(c_Key3(3, 1, -3), common_term_54);
    result_by_q.set(c_Key1(-3), common_term_54);
    // q = -2
    result_by_lpq.set(c_Key3(3, 1, -2), common_term_55);
    result_by_q.set(c_Key1(-2), common_term_55);
    // q = -1
    result_by_lpq.set(c_Key3(3, 1, -1), common_term_56);
    result_by_q.set(c_Key1(-1), common_term_56);
    // q = 0
    result_by_lpq.set(c_Key3(3, 1, 0), common_term_57);
    result_by_q.set(c_Key1(0), common_term_57);
    // q = 1
    result_by_lpq.set(c_Key3(3, 1, 1), common_term_58);
    result_by_q.set(c_Key1(1), common_term_58);
    // q = 2
    result_by_lpq.set(c_Key3(3, 1, 2), common_term_59);
    result_by_q.set(c_Key1(2), common_term_59);
    // q = 3
    result_by_lpq.set(c_Key3(3, 1, 3), common_term_60);
    result_by_q.set(c_Key1(3), common_term_60);
    // q = 4
    result_by_lpq.set(c_Key3(3, 1, 4), common_term_61);
    result_by_q.set(c_Key1(4), common_term_61);
    // q = 5
    result_by_lpq.set(c_Key3(3, 1, 5), common_term_62);
    result_by_q.set(c_Key1(5), common_term_62);
    // q = 6
    result_by_lpq.set(c_Key3(3, 1, 6), common_term_63);
    result_by_q.set(c_Key1(6), common_term_63);
    // q = 7
    result_by_lpq.set(c_Key3(3, 1, 7), common_term_64);
    result_by_q.set(c_Key1(7), common_term_64);
    // q = 8
    result_by_lpq.set(c_Key3(3, 1, 8), common_term_65);
    result_by_q.set(c_Key1(8), common_term_65);
    // q = 9
    result_by_lpq.set(c_Key3(3, 1, 9), common_term_66);
    result_by_q.set(c_Key1(9), common_term_66);
    // q = 10
    result_by_lpq.set(c_Key3(3, 1, 10), common_term_67);
    result_by_q.set(c_Key1(10), common_term_67);
    // q = 11
    result_by_lpq.set(c_Key3(3, 1, 11), common_term_68);
    result_by_q.set(c_Key1(11), common_term_68);
    // q = 12
    result_by_lpq.set(c_Key3(3, 1, 12), common_term_69);
    result_by_q.set(c_Key1(12), common_term_69);
    // q = 13
    result_by_lpq.set(c_Key3(3, 1, 13), common_term_70);
    result_by_q.set(c_Key1(13), common_term_70);
    // q = 14
    result_by_lpq.set(c_Key3(3, 1, 14), common_term_71);
    result_by_q.set(c_Key1(14), common_term_71);
    // q = 15
    result_by_lpq.set(c_Key3(3, 1, 15), common_term_72);
    result_by_q.set(c_Key1(15), common_term_72);
    // q = 16
    result_by_lpq.set(c_Key3(3, 1, 16), common_term_73);
    result_by_q.set(c_Key1(16), common_term_73);
    // q = 17
    result_by_lpq.set(c_Key3(3, 1, 17), common_term_74);
    result_by_q.set(c_Key1(17), common_term_74);
    // q = 18
    result_by_lpq.set(c_Key3(3, 1, 18), common_term_75);
    result_by_q.set(c_Key1(18), common_term_75);
    // q = 19
    result_by_lpq.set(c_Key3(3, 1, 19), common_term_76);
    result_by_q.set(c_Key1(19), common_term_76);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 1), result_by_q);
    result_by_q.clear();

    // l , p = (3, 2).
    // q = -19
    result_by_lpq.set(c_Key3(3, 2, -19), common_term_76);
    result_by_q.set(c_Key1(-19), common_term_76);
    // q = -18
    result_by_lpq.set(c_Key3(3, 2, -18), common_term_75);
    result_by_q.set(c_Key1(-18), common_term_75);
    // q = -17
    result_by_lpq.set(c_Key3(3, 2, -17), common_term_74);
    result_by_q.set(c_Key1(-17), common_term_74);
    // q = -16
    result_by_lpq.set(c_Key3(3, 2, -16), common_term_73);
    result_by_q.set(c_Key1(-16), common_term_73);
    // q = -15
    result_by_lpq.set(c_Key3(3, 2, -15), common_term_72);
    result_by_q.set(c_Key1(-15), common_term_72);
    // q = -14
    result_by_lpq.set(c_Key3(3, 2, -14), common_term_71);
    result_by_q.set(c_Key1(-14), common_term_71);
    // q = -13
    result_by_lpq.set(c_Key3(3, 2, -13), common_term_70);
    result_by_q.set(c_Key1(-13), common_term_70);
    // q = -12
    result_by_lpq.set(c_Key3(3, 2, -12), common_term_69);
    result_by_q.set(c_Key1(-12), common_term_69);
    // q = -11
    result_by_lpq.set(c_Key3(3, 2, -11), common_term_68);
    result_by_q.set(c_Key1(-11), common_term_68);
    // q = -10
    result_by_lpq.set(c_Key3(3, 2, -10), common_term_67);
    result_by_q.set(c_Key1(-10), common_term_67);
    // q = -9
    result_by_lpq.set(c_Key3(3, 2, -9), common_term_66);
    result_by_q.set(c_Key1(-9), common_term_66);
    // q = -8
    result_by_lpq.set(c_Key3(3, 2, -8), common_term_65);
    result_by_q.set(c_Key1(-8), common_term_65);
    // q = -7
    result_by_lpq.set(c_Key3(3, 2, -7), common_term_64);
    result_by_q.set(c_Key1(-7), common_term_64);
    // q = -6
    result_by_lpq.set(c_Key3(3, 2, -6), common_term_63);
    result_by_q.set(c_Key1(-6), common_term_63);
    // q = -5
    result_by_lpq.set(c_Key3(3, 2, -5), common_term_62);
    result_by_q.set(c_Key1(-5), common_term_62);
    // q = -4
    result_by_lpq.set(c_Key3(3, 2, -4), common_term_61);
    result_by_q.set(c_Key1(-4), common_term_61);
    // q = -3
    result_by_lpq.set(c_Key3(3, 2, -3), common_term_60);
    result_by_q.set(c_Key1(-3), common_term_60);
    // q = -2
    result_by_lpq.set(c_Key3(3, 2, -2), common_term_59);
    result_by_q.set(c_Key1(-2), common_term_59);
    // q = -1
    result_by_lpq.set(c_Key3(3, 2, -1), common_term_58);
    result_by_q.set(c_Key1(-1), common_term_58);
    // q = 0
    result_by_lpq.set(c_Key3(3, 2, 0), common_term_57);
    result_by_q.set(c_Key1(0), common_term_57);
    // q = 1
    result_by_lpq.set(c_Key3(3, 2, 1), common_term_56);
    result_by_q.set(c_Key1(1), common_term_56);
    // q = 2
    result_by_lpq.set(c_Key3(3, 2, 2), common_term_55);
    result_by_q.set(c_Key1(2), common_term_55);
    // q = 3
    result_by_lpq.set(c_Key3(3, 2, 3), common_term_54);
    result_by_q.set(c_Key1(3), common_term_54);
    // q = 4
    result_by_lpq.set(c_Key3(3, 2, 4), common_term_53);
    result_by_q.set(c_Key1(4), common_term_53);
    // q = 5
    result_by_lpq.set(c_Key3(3, 2, 5), common_term_52);
    result_by_q.set(c_Key1(5), common_term_52);
    // q = 6
    result_by_lpq.set(c_Key3(3, 2, 6), common_term_51);
    result_by_q.set(c_Key1(6), common_term_51);
    // q = 7
    result_by_lpq.set(c_Key3(3, 2, 7), common_term_50);
    result_by_q.set(c_Key1(7), common_term_50);
    // q = 8
    result_by_lpq.set(c_Key3(3, 2, 8), common_term_49);
    result_by_q.set(c_Key1(8), common_term_49);
    // q = 9
    result_by_lpq.set(c_Key3(3, 2, 9), common_term_48);
    result_by_q.set(c_Key1(9), common_term_48);
    // q = 10
    result_by_lpq.set(c_Key3(3, 2, 10), common_term_47);
    result_by_q.set(c_Key1(10), common_term_47);
    // q = 11
    result_by_lpq.set(c_Key3(3, 2, 11), common_term_46);
    result_by_q.set(c_Key1(11), common_term_46);
    // q = 12
    result_by_lpq.set(c_Key3(3, 2, 12), common_term_45);
    result_by_q.set(c_Key1(12), common_term_45);
    // q = 13
    result_by_lpq.set(c_Key3(3, 2, 13), common_term_44);
    result_by_q.set(c_Key1(13), common_term_44);
    // q = 14
    result_by_lpq.set(c_Key3(3, 2, 14), common_term_43);
    result_by_q.set(c_Key1(14), common_term_43);
    // q = 15
    result_by_lpq.set(c_Key3(3, 2, 15), common_term_42);
    result_by_q.set(c_Key1(15), common_term_42);
    // q = 16
    result_by_lpq.set(c_Key3(3, 2, 16), common_term_41);
    result_by_q.set(c_Key1(16), common_term_41);
    // q = 17
    result_by_lpq.set(c_Key3(3, 2, 17), common_term_40);
    result_by_q.set(c_Key1(17), common_term_40);
    // q = 18
    result_by_lpq.set(c_Key3(3, 2, 18), common_term_39);
    result_by_q.set(c_Key1(18), common_term_39);
    // q = 19
    result_by_lpq.set(c_Key3(3, 2, 19), common_term_38);
    result_by_q.set(c_Key1(19), common_term_38);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 2), result_by_q);
    result_by_q.clear();

    // l , p = (3, 3).
    // q = -19
    result_by_lpq.set(c_Key3(3, 3, -19), common_term_37);
    result_by_q.set(c_Key1(-19), common_term_37);
    // q = -18
    result_by_lpq.set(c_Key3(3, 3, -18), common_term_36);
    result_by_q.set(c_Key1(-18), common_term_36);
    // q = -17
    result_by_lpq.set(c_Key3(3, 3, -17), common_term_35);
    result_by_q.set(c_Key1(-17), common_term_35);
    // q = -16
    result_by_lpq.set(c_Key3(3, 3, -16), common_term_34);
    result_by_q.set(c_Key1(-16), common_term_34);
    // q = -15
    result_by_lpq.set(c_Key3(3, 3, -15), common_term_33);
    result_by_q.set(c_Key1(-15), common_term_33);
    // q = -14
    result_by_lpq.set(c_Key3(3, 3, -14), common_term_32);
    result_by_q.set(c_Key1(-14), common_term_32);
    // q = -13
    result_by_lpq.set(c_Key3(3, 3, -13), common_term_31);
    result_by_q.set(c_Key1(-13), common_term_31);
    // q = -12
    result_by_lpq.set(c_Key3(3, 3, -12), common_term_30);
    result_by_q.set(c_Key1(-12), common_term_30);
    // q = -11
    result_by_lpq.set(c_Key3(3, 3, -11), common_term_29);
    result_by_q.set(c_Key1(-11), common_term_29);
    // q = -10
    result_by_lpq.set(c_Key3(3, 3, -10), common_term_28);
    result_by_q.set(c_Key1(-10), common_term_28);
    // q = -9
    result_by_lpq.set(c_Key3(3, 3, -9), common_term_27);
    result_by_q.set(c_Key1(-9), common_term_27);
    // q = -8
    result_by_lpq.set(c_Key3(3, 3, -8), common_term_26);
    result_by_q.set(c_Key1(-8), common_term_26);
    // q = -7
    result_by_lpq.set(c_Key3(3, 3, -7), common_term_25);
    result_by_q.set(c_Key1(-7), common_term_25);
    // q = -6
    result_by_lpq.set(c_Key3(3, 3, -6), common_term_24);
    result_by_q.set(c_Key1(-6), common_term_24);
    // q = -5
    result_by_lpq.set(c_Key3(3, 3, -5), common_term_23);
    result_by_q.set(c_Key1(-5), common_term_23);
    // q = -4
    result_by_lpq.set(c_Key3(3, 3, -4), common_term_22);
    result_by_q.set(c_Key1(-4), common_term_22);
    // q = -3
    result_by_lpq.set(c_Key3(3, 3, -3), common_term_21);
    result_by_q.set(c_Key1(-3), common_term_21);
    // q = -2
    result_by_lpq.set(c_Key3(3, 3, -2), common_term_20);
    result_by_q.set(c_Key1(-2), common_term_20);
    // q = -1
    result_by_lpq.set(c_Key3(3, 3, -1), common_term_19);
    result_by_q.set(c_Key1(-1), common_term_19);
    // q = 0
    result_by_lpq.set(c_Key3(3, 3, 0), common_term_18);
    result_by_q.set(c_Key1(0), common_term_18);
    // q = 1
    result_by_lpq.set(c_Key3(3, 3, 1), common_term_17);
    result_by_q.set(c_Key1(1), common_term_17);
    // q = 2
    result_by_lpq.set(c_Key3(3, 3, 2), common_term_16);
    result_by_q.set(c_Key1(2), common_term_16);
    // q = 4
    result_by_lpq.set(c_Key3(3, 3, 4), common_term_15);
    result_by_q.set(c_Key1(4), common_term_15);
    // q = 5
    result_by_lpq.set(c_Key3(3, 3, 5), common_term_14);
    result_by_q.set(c_Key1(5), common_term_14);
    // q = 6
    result_by_lpq.set(c_Key3(3, 3, 6), common_term_13);
    result_by_q.set(c_Key1(6), common_term_13);
    // q = 7
    result_by_lpq.set(c_Key3(3, 3, 7), common_term_12);
    result_by_q.set(c_Key1(7), common_term_12);
    // q = 8
    result_by_lpq.set(c_Key3(3, 3, 8), common_term_11);
    result_by_q.set(c_Key1(8), common_term_11);
    // q = 9
    result_by_lpq.set(c_Key3(3, 3, 9), common_term_10);
    result_by_q.set(c_Key1(9), common_term_10);
    // q = 10
    result_by_lpq.set(c_Key3(3, 3, 10), common_term_9);
    result_by_q.set(c_Key1(10), common_term_9);
    // q = 11
    result_by_lpq.set(c_Key3(3, 3, 11), common_term_8);
    result_by_q.set(c_Key1(11), common_term_8);
    // q = 12
    result_by_lpq.set(c_Key3(3, 3, 12), common_term_7);
    result_by_q.set(c_Key1(12), common_term_7);
    // q = 13
    result_by_lpq.set(c_Key3(3, 3, 13), common_term_6);
    result_by_q.set(c_Key1(13), common_term_6);
    // q = 14
    result_by_lpq.set(c_Key3(3, 3, 14), common_term_5);
    result_by_q.set(c_Key1(14), common_term_5);
    // q = 15
    result_by_lpq.set(c_Key3(3, 3, 15), common_term_4);
    result_by_q.set(c_Key1(15), common_term_4);
    // q = 16
    result_by_lpq.set(c_Key3(3, 3, 16), common_term_3);
    result_by_q.set(c_Key1(16), common_term_3);
    // q = 17
    result_by_lpq.set(c_Key3(3, 3, 17), common_term_2);
    result_by_q.set(c_Key1(17), common_term_2);
    // q = 18
    result_by_lpq.set(c_Key3(3, 3, 18), common_term_1);
    result_by_q.set(c_Key1(18), common_term_1);
    // q = 19
    result_by_lpq.set(c_Key3(3, 3, 19), common_term_0);
    result_by_q.set(c_Key1(19), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(3, 3), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}
