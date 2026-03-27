#pragma once

#include "eccentricity_common_.hpp"

EccentricityFuncOutput c_eccentricity_function_l2_e1(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 2.
    // Functions are taylor expanded and truncated at eccentricity^1.

    // Functions Calculated for l = 2.

    c_IntMap<c_Key3, double> result_by_lpq(3);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(3);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;

    c_IntMap<c_Key1, double> result_by_q(1);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (2, 0).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(2, 0, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 0), result_by_q);
    result_by_q.clear();

    // l , p = (2, 1).
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -1.5);
    result_by_lpq.set(c_Key3(2, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 1), result_by_q);
    result_by_q.clear();

    // l , p = (2, 2).
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(2, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 2), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l2_e2(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 2.
    // Functions are taylor expanded and truncated at eccentricity^2.

    // Functions Calculated for l = 2.

    c_IntMap<c_Key3, double> result_by_lpq(9);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(3);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = -0.5*eccentricity;
    double common_term_1 = 3.5*eccentricity;
    double common_term_2 = 1.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(3);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (2, 0).
    // q = -1
    result_by_lpq.set(c_Key3(2, 0, -1), common_term_0);
    result_by_q.set(c_Key1(-1), common_term_0);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(2, 0, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(2, 0, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 0), result_by_q);
    result_by_q.clear();

    // l , p = (2, 1).
    // q = -1
    result_by_lpq.set(c_Key3(2, 1, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -1.5);
    result_by_lpq.set(c_Key3(2, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(2, 1, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 1), result_by_q);
    result_by_q.clear();

    // l , p = (2, 2).
    // q = -1
    result_by_lpq.set(c_Key3(2, 2, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    tmp_double = 1.0;
    result_by_lpq.set(c_Key3(2, 2, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(2, 2, 1), common_term_0);
    result_by_q.set(c_Key1(1), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 2), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l2_e3(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 2.
    // Functions are taylor expanded and truncated at eccentricity^3.

    // Functions Calculated for l = 2.

    c_IntMap<c_Key3, double> result_by_lpq(13);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(3);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double common_term_0 = -0.5*eccentricity;
    double common_term_1 = 1.0 - 2.5*eccentricity_2;
    double common_term_2 = 3.5*eccentricity;
    double common_term_3 = 8.5*eccentricity_2;
    double common_term_4 = 2.25*eccentricity_2;
    double common_term_5 = 1.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(5);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (2, 0).
    // q = -1
    result_by_lpq.set(c_Key3(2, 0, -1), common_term_0);
    result_by_q.set(c_Key1(-1), common_term_0);
    // q = 0
    result_by_lpq.set(c_Key3(2, 0, 0), common_term_1);
    result_by_q.set(c_Key1(0), common_term_1);
    // q = 1
    result_by_lpq.set(c_Key3(2, 0, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // q = 2
    result_by_lpq.set(c_Key3(2, 0, 2), common_term_3);
    result_by_q.set(c_Key1(2), common_term_3);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 0), result_by_q);
    result_by_q.clear();

    // l , p = (2, 1).
    // q = -2
    result_by_lpq.set(c_Key3(2, 1, -2), common_term_4);
    result_by_q.set(c_Key1(-2), common_term_4);
    // q = -1
    result_by_lpq.set(c_Key3(2, 1, -1), common_term_5);
    result_by_q.set(c_Key1(-1), common_term_5);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -1.5);
    result_by_lpq.set(c_Key3(2, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(2, 1, 1), common_term_5);
    result_by_q.set(c_Key1(1), common_term_5);
    // q = 2
    result_by_lpq.set(c_Key3(2, 1, 2), common_term_4);
    result_by_q.set(c_Key1(2), common_term_4);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 1), result_by_q);
    result_by_q.clear();

    // l , p = (2, 2).
    // q = -2
    result_by_lpq.set(c_Key3(2, 2, -2), common_term_3);
    result_by_q.set(c_Key1(-2), common_term_3);
    // q = -1
    result_by_lpq.set(c_Key3(2, 2, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    result_by_lpq.set(c_Key3(2, 2, 0), common_term_1);
    result_by_q.set(c_Key1(0), common_term_1);
    // q = 1
    result_by_lpq.set(c_Key3(2, 2, 1), common_term_0);
    result_by_q.set(c_Key1(1), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 2), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l2_e4(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 2.
    // Functions are taylor expanded and truncated at eccentricity^4.

    // Functions Calculated for l = 2.

    c_IntMap<c_Key3, double> result_by_lpq(19);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(3);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double common_term_0 = 0.020833333333333333*eccentricity_3;
    double common_term_1 = 0.0625*eccentricity_3 - 0.5*eccentricity;
    double common_term_2 = 1.0 - 2.5*eccentricity_2;
    double common_term_3 = -7.6875*eccentricity_3 + 3.5*eccentricity;
    double common_term_4 = 8.5*eccentricity_2;
    double common_term_5 = 17.604166666666667*eccentricity_3;
    double common_term_6 = 3.3125*eccentricity_3;
    double common_term_7 = 2.25*eccentricity_2;
    double common_term_8 = 1.6875*eccentricity_3 + 1.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(7);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (2, 0).
    // q = -3
    result_by_lpq.set(c_Key3(2, 0, -3), common_term_0);
    result_by_q.set(c_Key1(-3), common_term_0);
    // q = -1
    result_by_lpq.set(c_Key3(2, 0, -1), common_term_1);
    result_by_q.set(c_Key1(-1), common_term_1);
    // q = 0
    result_by_lpq.set(c_Key3(2, 0, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(2, 0, 1), common_term_3);
    result_by_q.set(c_Key1(1), common_term_3);
    // q = 2
    result_by_lpq.set(c_Key3(2, 0, 2), common_term_4);
    result_by_q.set(c_Key1(2), common_term_4);
    // q = 3
    result_by_lpq.set(c_Key3(2, 0, 3), common_term_5);
    result_by_q.set(c_Key1(3), common_term_5);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 0), result_by_q);
    result_by_q.clear();

    // l , p = (2, 1).
    // q = -3
    result_by_lpq.set(c_Key3(2, 1, -3), common_term_6);
    result_by_q.set(c_Key1(-3), common_term_6);
    // q = -2
    result_by_lpq.set(c_Key3(2, 1, -2), common_term_7);
    result_by_q.set(c_Key1(-2), common_term_7);
    // q = -1
    result_by_lpq.set(c_Key3(2, 1, -1), common_term_8);
    result_by_q.set(c_Key1(-1), common_term_8);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -1.5);
    result_by_lpq.set(c_Key3(2, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(2, 1, 1), common_term_8);
    result_by_q.set(c_Key1(1), common_term_8);
    // q = 2
    result_by_lpq.set(c_Key3(2, 1, 2), common_term_7);
    result_by_q.set(c_Key1(2), common_term_7);
    // q = 3
    result_by_lpq.set(c_Key3(2, 1, 3), common_term_6);
    result_by_q.set(c_Key1(3), common_term_6);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 1), result_by_q);
    result_by_q.clear();

    // l , p = (2, 2).
    // q = -3
    result_by_lpq.set(c_Key3(2, 2, -3), common_term_5);
    result_by_q.set(c_Key1(-3), common_term_5);
    // q = -2
    result_by_lpq.set(c_Key3(2, 2, -2), common_term_4);
    result_by_q.set(c_Key1(-2), common_term_4);
    // q = -1
    result_by_lpq.set(c_Key3(2, 2, -1), common_term_3);
    result_by_q.set(c_Key1(-1), common_term_3);
    // q = 0
    result_by_lpq.set(c_Key3(2, 2, 0), common_term_2);
    result_by_q.set(c_Key1(0), common_term_2);
    // q = 1
    result_by_lpq.set(c_Key3(2, 2, 1), common_term_1);
    result_by_q.set(c_Key1(1), common_term_1);
    // q = 3
    result_by_lpq.set(c_Key3(2, 2, 3), common_term_0);
    result_by_q.set(c_Key1(3), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 2), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l2_e5(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 2.
    // Functions are taylor expanded and truncated at eccentricity^5.

    // Functions Calculated for l = 2.

    c_IntMap<c_Key3, double> result_by_lpq(25);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(3);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double common_term_0 = 0.041666666666666667*eccentricity_4;
    double common_term_1 = 0.020833333333333333*eccentricity_3;
    double common_term_2 = 0.0625*eccentricity_3 - 0.5*eccentricity;
    double common_term_3 = 0.8125*eccentricity_4 - 2.5*eccentricity_2 + 1.0;
    double common_term_4 = -7.6875*eccentricity_3 + 3.5*eccentricity;
    double common_term_5 = -19.166666666666667*eccentricity_4 + 8.5*eccentricity_2;
    double common_term_6 = 17.604166666666667*eccentricity_3;
    double common_term_7 = 33.3125*eccentricity_4;
    double common_term_8 = 4.8125*eccentricity_4;
    double common_term_9 = 3.3125*eccentricity_3;
    double common_term_10 = 1.75*eccentricity_4 + 2.25*eccentricity_2;
    double common_term_11 = 1.6875*eccentricity_3 + 1.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(9);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (2, 0).
    // q = -4
    result_by_lpq.set(c_Key3(2, 0, -4), common_term_0);
    result_by_q.set(c_Key1(-4), common_term_0);
    // q = -3
    result_by_lpq.set(c_Key3(2, 0, -3), common_term_1);
    result_by_q.set(c_Key1(-3), common_term_1);
    // q = -1
    result_by_lpq.set(c_Key3(2, 0, -1), common_term_2);
    result_by_q.set(c_Key1(-1), common_term_2);
    // q = 0
    result_by_lpq.set(c_Key3(2, 0, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(2, 0, 1), common_term_4);
    result_by_q.set(c_Key1(1), common_term_4);
    // q = 2
    result_by_lpq.set(c_Key3(2, 0, 2), common_term_5);
    result_by_q.set(c_Key1(2), common_term_5);
    // q = 3
    result_by_lpq.set(c_Key3(2, 0, 3), common_term_6);
    result_by_q.set(c_Key1(3), common_term_6);
    // q = 4
    result_by_lpq.set(c_Key3(2, 0, 4), common_term_7);
    result_by_q.set(c_Key1(4), common_term_7);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 0), result_by_q);
    result_by_q.clear();

    // l , p = (2, 1).
    // q = -4
    result_by_lpq.set(c_Key3(2, 1, -4), common_term_8);
    result_by_q.set(c_Key1(-4), common_term_8);
    // q = -3
    result_by_lpq.set(c_Key3(2, 1, -3), common_term_9);
    result_by_q.set(c_Key1(-3), common_term_9);
    // q = -2
    result_by_lpq.set(c_Key3(2, 1, -2), common_term_10);
    result_by_q.set(c_Key1(-2), common_term_10);
    // q = -1
    result_by_lpq.set(c_Key3(2, 1, -1), common_term_11);
    result_by_q.set(c_Key1(-1), common_term_11);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -1.5);
    result_by_lpq.set(c_Key3(2, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(2, 1, 1), common_term_11);
    result_by_q.set(c_Key1(1), common_term_11);
    // q = 2
    result_by_lpq.set(c_Key3(2, 1, 2), common_term_10);
    result_by_q.set(c_Key1(2), common_term_10);
    // q = 3
    result_by_lpq.set(c_Key3(2, 1, 3), common_term_9);
    result_by_q.set(c_Key1(3), common_term_9);
    // q = 4
    result_by_lpq.set(c_Key3(2, 1, 4), common_term_8);
    result_by_q.set(c_Key1(4), common_term_8);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 1), result_by_q);
    result_by_q.clear();

    // l , p = (2, 2).
    // q = -4
    result_by_lpq.set(c_Key3(2, 2, -4), common_term_7);
    result_by_q.set(c_Key1(-4), common_term_7);
    // q = -3
    result_by_lpq.set(c_Key3(2, 2, -3), common_term_6);
    result_by_q.set(c_Key1(-3), common_term_6);
    // q = -2
    result_by_lpq.set(c_Key3(2, 2, -2), common_term_5);
    result_by_q.set(c_Key1(-2), common_term_5);
    // q = -1
    result_by_lpq.set(c_Key3(2, 2, -1), common_term_4);
    result_by_q.set(c_Key1(-1), common_term_4);
    // q = 0
    result_by_lpq.set(c_Key3(2, 2, 0), common_term_3);
    result_by_q.set(c_Key1(0), common_term_3);
    // q = 1
    result_by_lpq.set(c_Key3(2, 2, 1), common_term_2);
    result_by_q.set(c_Key1(1), common_term_2);
    // q = 3
    result_by_lpq.set(c_Key3(2, 2, 3), common_term_1);
    result_by_q.set(c_Key1(3), common_term_1);
    // q = 4
    result_by_lpq.set(c_Key3(2, 2, 4), common_term_0);
    result_by_q.set(c_Key1(4), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 2), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l2_e10(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 2.
    // Functions are taylor expanded and truncated at eccentricity^10.

    // Functions Calculated for l = 2.

    c_IntMap<c_Key3, double> result_by_lpq(55);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(3);
    // Optimizations
    double eccentricity_2 = eccentricity * eccentricity;
    double eccentricity_4 = eccentricity_2 * eccentricity_2;
    double eccentricity_8 = eccentricity_4 * eccentricity_4;
    double eccentricity_9 = eccentricity_8 * eccentricity;
    double eccentricity_3 = eccentricity_2 * eccentricity;
    double eccentricity_6 = eccentricity_3 * eccentricity_3;
    double eccentricity_7 = eccentricity_6 * eccentricity;
    double eccentricity_5 = eccentricity_4 * eccentricity;
    double common_term_0 = 0.21719477147231867*eccentricity_9;
    double common_term_1 = 0.16272321428571429*eccentricity_8;
    double common_term_2 = 0.041628640795510913*eccentricity_9 + 0.12110150049603175*eccentricity_7;
    double common_term_3 = 0.044444444444444444*eccentricity_8 + 0.088888888888888889*eccentricity_6;
    double common_term_4 = 0.033759416852678571*eccentricity_9 + 0.03955078125*eccentricity_7 + 0.06328125*eccentricity_5;
    double common_term_5 = 0.022743055555555556*eccentricity_8 + 0.029166666666666667*eccentricity_6 + 0.041666666666666667*eccentricity_4;
    double common_term_6 = 0.0075841833043981482*eccentricity_9 + 0.010188802083333333*eccentricity_7 + 0.014322916666666667*eccentricity_5 + 0.020833333333333333*eccentricity_3;
    double common_term_7 = -0.0061692979600694444*eccentricity_9 - 0.0077582465277777778*eccentricity_7 - 0.013020833333333333*eccentricity_5 + 0.0625*eccentricity_3 - 0.5*eccentricity;
    double common_term_8 = -0.0086805555555555556*eccentricity_8 - 0.12152777777777778*eccentricity_6 + 0.8125*eccentricity_4 - 2.5*eccentricity_2 + 1.0;
    double common_term_9 = 0.082562255859375*eccentricity_9 - 0.86083984375*eccentricity_7 + 3.8203125*eccentricity_5 - 7.6875*eccentricity_3 + 3.5*eccentricity;
    double common_term_10 = -3.9527777777777778*eccentricity_8 + 12.520833333333333*eccentricity_6 - 19.166666666666667*eccentricity_4 + 8.5*eccentricity_2;
    double common_term_11 = -13.840795446325231*eccentricity_9 + 33.890787760416667*eccentricity_7 - 42.350260416666667*eccentricity_5 + 17.604166666666667*eccentricity_3;
    double common_term_12 = 81.34921875*eccentricity_8 - 86.41875*eccentricity_6 + 33.3125*eccentricity_4;
    double common_term_13 = 179.69544406467014*eccentricity_9 - 166.61648220486111*eccentricity_7 + 59.465364583333333*eccentricity_5;
    double common_term_14 = -307.8281746031746*eccentricity_8 + 101.90138888888889*eccentricity_6;
    double common_term_15 = -550.10469491141183*eccentricity_9 + 169.42345145089286*eccentricity_7;
    double common_term_16 = 275.21050347222222*eccentricity_8;
    double common_term_17 = 438.87348632274271*eccentricity_9;
    double common_term_18 = 28.063401140485491*eccentricity_9;
    double common_term_19 = 19.903013392857143*eccentricity_8;
    double common_term_20 = -9.3806243896484375*eccentricity_9 + 14.065462239583333*eccentricity_7;
    double common_term_21 = -3.9053571428571429*eccentricity_8 + 9.896875*eccentricity_6;
    double common_term_22 = 3.6474173409598214*eccentricity_9 - 0.81168619791666667*eccentricity_7 + 6.92578125*eccentricity_5;
    double common_term_23 = 2.809375*eccentricity_8 + 0.80625*eccentricity_6 + 4.8125*eccentricity_4;
    double common_term_24 = 2.58248291015625*eccentricity_9 + 2.41728515625*eccentricity_7 + 1.53515625*eccentricity_5 + 3.3125*eccentricity_3;
    double common_term_25 = 2.4625*eccentricity_8 + 2.203125*eccentricity_6 + 1.75*eccentricity_4 + 2.25*eccentricity_2;
    double common_term_26 = 2.587322998046875*eccentricity_9 + 2.3289388020833333*eccentricity_7 + 2.0390625*eccentricity_5 + 1.6875*eccentricity_3 + 1.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(19);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (2, 0).
    // q = -9
    result_by_lpq.set(c_Key3(2, 0, -9), common_term_0);
    result_by_q.set(c_Key1(-9), common_term_0);
    // q = -8
    result_by_lpq.set(c_Key3(2, 0, -8), common_term_1);
    result_by_q.set(c_Key1(-8), common_term_1);
    // q = -7
    result_by_lpq.set(c_Key3(2, 0, -7), common_term_2);
    result_by_q.set(c_Key1(-7), common_term_2);
    // q = -6
    result_by_lpq.set(c_Key3(2, 0, -6), common_term_3);
    result_by_q.set(c_Key1(-6), common_term_3);
    // q = -5
    result_by_lpq.set(c_Key3(2, 0, -5), common_term_4);
    result_by_q.set(c_Key1(-5), common_term_4);
    // q = -4
    result_by_lpq.set(c_Key3(2, 0, -4), common_term_5);
    result_by_q.set(c_Key1(-4), common_term_5);
    // q = -3
    result_by_lpq.set(c_Key3(2, 0, -3), common_term_6);
    result_by_q.set(c_Key1(-3), common_term_6);
    // q = -1
    result_by_lpq.set(c_Key3(2, 0, -1), common_term_7);
    result_by_q.set(c_Key1(-1), common_term_7);
    // q = 0
    result_by_lpq.set(c_Key3(2, 0, 0), common_term_8);
    result_by_q.set(c_Key1(0), common_term_8);
    // q = 1
    result_by_lpq.set(c_Key3(2, 0, 1), common_term_9);
    result_by_q.set(c_Key1(1), common_term_9);
    // q = 2
    result_by_lpq.set(c_Key3(2, 0, 2), common_term_10);
    result_by_q.set(c_Key1(2), common_term_10);
    // q = 3
    result_by_lpq.set(c_Key3(2, 0, 3), common_term_11);
    result_by_q.set(c_Key1(3), common_term_11);
    // q = 4
    result_by_lpq.set(c_Key3(2, 0, 4), common_term_12);
    result_by_q.set(c_Key1(4), common_term_12);
    // q = 5
    result_by_lpq.set(c_Key3(2, 0, 5), common_term_13);
    result_by_q.set(c_Key1(5), common_term_13);
    // q = 6
    result_by_lpq.set(c_Key3(2, 0, 6), common_term_14);
    result_by_q.set(c_Key1(6), common_term_14);
    // q = 7
    result_by_lpq.set(c_Key3(2, 0, 7), common_term_15);
    result_by_q.set(c_Key1(7), common_term_15);
    // q = 8
    result_by_lpq.set(c_Key3(2, 0, 8), common_term_16);
    result_by_q.set(c_Key1(8), common_term_16);
    // q = 9
    result_by_lpq.set(c_Key3(2, 0, 9), common_term_17);
    result_by_q.set(c_Key1(9), common_term_17);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 0), result_by_q);
    result_by_q.clear();

    // l , p = (2, 1).
    // q = -9
    result_by_lpq.set(c_Key3(2, 1, -9), common_term_18);
    result_by_q.set(c_Key1(-9), common_term_18);
    // q = -8
    result_by_lpq.set(c_Key3(2, 1, -8), common_term_19);
    result_by_q.set(c_Key1(-8), common_term_19);
    // q = -7
    result_by_lpq.set(c_Key3(2, 1, -7), common_term_20);
    result_by_q.set(c_Key1(-7), common_term_20);
    // q = -6
    result_by_lpq.set(c_Key3(2, 1, -6), common_term_21);
    result_by_q.set(c_Key1(-6), common_term_21);
    // q = -5
    result_by_lpq.set(c_Key3(2, 1, -5), common_term_22);
    result_by_q.set(c_Key1(-5), common_term_22);
    // q = -4
    result_by_lpq.set(c_Key3(2, 1, -4), common_term_23);
    result_by_q.set(c_Key1(-4), common_term_23);
    // q = -3
    result_by_lpq.set(c_Key3(2, 1, -3), common_term_24);
    result_by_q.set(c_Key1(-3), common_term_24);
    // q = -2
    result_by_lpq.set(c_Key3(2, 1, -2), common_term_25);
    result_by_q.set(c_Key1(-2), common_term_25);
    // q = -1
    result_by_lpq.set(c_Key3(2, 1, -1), common_term_26);
    result_by_q.set(c_Key1(-1), common_term_26);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -1.5);
    result_by_lpq.set(c_Key3(2, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(2, 1, 1), common_term_26);
    result_by_q.set(c_Key1(1), common_term_26);
    // q = 2
    result_by_lpq.set(c_Key3(2, 1, 2), common_term_25);
    result_by_q.set(c_Key1(2), common_term_25);
    // q = 3
    result_by_lpq.set(c_Key3(2, 1, 3), common_term_24);
    result_by_q.set(c_Key1(3), common_term_24);
    // q = 4
    result_by_lpq.set(c_Key3(2, 1, 4), common_term_23);
    result_by_q.set(c_Key1(4), common_term_23);
    // q = 5
    result_by_lpq.set(c_Key3(2, 1, 5), common_term_22);
    result_by_q.set(c_Key1(5), common_term_22);
    // q = 6
    result_by_lpq.set(c_Key3(2, 1, 6), common_term_21);
    result_by_q.set(c_Key1(6), common_term_21);
    // q = 7
    result_by_lpq.set(c_Key3(2, 1, 7), common_term_20);
    result_by_q.set(c_Key1(7), common_term_20);
    // q = 8
    result_by_lpq.set(c_Key3(2, 1, 8), common_term_19);
    result_by_q.set(c_Key1(8), common_term_19);
    // q = 9
    result_by_lpq.set(c_Key3(2, 1, 9), common_term_18);
    result_by_q.set(c_Key1(9), common_term_18);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 1), result_by_q);
    result_by_q.clear();

    // l , p = (2, 2).
    // q = -9
    result_by_lpq.set(c_Key3(2, 2, -9), common_term_17);
    result_by_q.set(c_Key1(-9), common_term_17);
    // q = -8
    result_by_lpq.set(c_Key3(2, 2, -8), common_term_16);
    result_by_q.set(c_Key1(-8), common_term_16);
    // q = -7
    result_by_lpq.set(c_Key3(2, 2, -7), common_term_15);
    result_by_q.set(c_Key1(-7), common_term_15);
    // q = -6
    result_by_lpq.set(c_Key3(2, 2, -6), common_term_14);
    result_by_q.set(c_Key1(-6), common_term_14);
    // q = -5
    result_by_lpq.set(c_Key3(2, 2, -5), common_term_13);
    result_by_q.set(c_Key1(-5), common_term_13);
    // q = -4
    result_by_lpq.set(c_Key3(2, 2, -4), common_term_12);
    result_by_q.set(c_Key1(-4), common_term_12);
    // q = -3
    result_by_lpq.set(c_Key3(2, 2, -3), common_term_11);
    result_by_q.set(c_Key1(-3), common_term_11);
    // q = -2
    result_by_lpq.set(c_Key3(2, 2, -2), common_term_10);
    result_by_q.set(c_Key1(-2), common_term_10);
    // q = -1
    result_by_lpq.set(c_Key3(2, 2, -1), common_term_9);
    result_by_q.set(c_Key1(-1), common_term_9);
    // q = 0
    result_by_lpq.set(c_Key3(2, 2, 0), common_term_8);
    result_by_q.set(c_Key1(0), common_term_8);
    // q = 1
    result_by_lpq.set(c_Key3(2, 2, 1), common_term_7);
    result_by_q.set(c_Key1(1), common_term_7);
    // q = 3
    result_by_lpq.set(c_Key3(2, 2, 3), common_term_6);
    result_by_q.set(c_Key1(3), common_term_6);
    // q = 4
    result_by_lpq.set(c_Key3(2, 2, 4), common_term_5);
    result_by_q.set(c_Key1(4), common_term_5);
    // q = 5
    result_by_lpq.set(c_Key3(2, 2, 5), common_term_4);
    result_by_q.set(c_Key1(5), common_term_4);
    // q = 6
    result_by_lpq.set(c_Key3(2, 2, 6), common_term_3);
    result_by_q.set(c_Key1(6), common_term_3);
    // q = 7
    result_by_lpq.set(c_Key3(2, 2, 7), common_term_2);
    result_by_q.set(c_Key1(7), common_term_2);
    // q = 8
    result_by_lpq.set(c_Key3(2, 2, 8), common_term_1);
    result_by_q.set(c_Key1(8), common_term_1);
    // q = 9
    result_by_lpq.set(c_Key3(2, 2, 9), common_term_0);
    result_by_q.set(c_Key1(9), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 2), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l2_e15(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 2.
    // Functions are taylor expanded and truncated at eccentricity^15.

    // Functions Calculated for l = 2.

    c_IntMap<c_Key3, double> result_by_lpq(85);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(3);
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
    double common_term_0 = 0.89889539032396175*eccentricity_14;
    double common_term_1 = 0.67675994590498271*eccentricity_13;
    double common_term_2 = -0.3332565249340423*eccentricity_14 + 0.50968644989912351*eccentricity_12;
    double common_term_3 = -0.16794225909493186*eccentricity_13 + 0.38386802078841569*eccentricity_11;
    double common_term_4 = 0.12532520976965421*eccentricity_14 - 0.065672599005932339*eccentricity_12 + 0.28895943562610229*eccentricity_10;
    double common_term_5 = 0.088173672850552098*eccentricity_13 - 0.0054298692868079668*eccentricity_11 + 0.21719477147231867*eccentricity_9;
    double common_term_6 = 0.047615031452922078*eccentricity_14 + 0.067123325892857143*eccentricity_12 + 0.027120535714285714*eccentricity_10 + 0.16272321428571429*eccentricity_8;
    double common_term_7 = 0.042511234558168262*eccentricity_13 + 0.05387545052449329*eccentricity_11 + 0.041628640795510913*eccentricity_9 + 0.12110150049603175*eccentricity_7;
    double common_term_8 = 0.029965461493239271*eccentricity_14 + 0.035552616108171664*eccentricity_12 + 0.043650793650793651*eccentricity_10 + 0.044444444444444444*eccentricity_8 + 0.088888888888888889*eccentricity_6;
    double common_term_9 = 0.022452848298209054*eccentricity_13 + 0.027257374354771205*eccentricity_11 + 0.033759416852678571*eccentricity_9 + 0.03955078125*eccentricity_7 + 0.06328125*eccentricity_5;
    double common_term_10 = 0.01172185593952087*eccentricity_14 + 0.014269696593915344*eccentricity_12 + 0.017786871693121693*eccentricity_10 + 0.022743055555555556*eccentricity_8 + 0.029166666666666667*eccentricity_6 + 0.041666666666666667*eccentricity_4;
    double common_term_11 = 0.0047346311649948201*eccentricity_13 + 0.0058915435952484292*eccentricity_11 + 0.0075841833043981482*eccentricity_9 + 0.010188802083333333*eccentricity_7 + 0.014322916666666667*eccentricity_5 + 0.020833333333333333*eccentricity_3;
    double common_term_12 = -0.0040929340517289634*eccentricity_13 - 0.0049673518428096065*eccentricity_11 - 0.0061692979600694444*eccentricity_9 - 0.0077582465277777778*eccentricity_7 - 0.013020833333333333*eccentricity_5 + 0.0625*eccentricity_3 - 0.5*eccentricity;
    double common_term_13 = -0.0095563234697814311*eccentricity_14 - 0.011227454668209877*eccentricity_12 - 0.013611111111111111*eccentricity_10 - 0.0086805555555555556*eccentricity_8 - 0.12152777777777778*eccentricity_6 + 0.8125*eccentricity_4 - 2.5*eccentricity_2 + 1.0;
    double common_term_14 = -0.016310808999197824*eccentricity_13 - 0.027522125244140625*eccentricity_11 + 0.082562255859375*eccentricity_9 - 0.86083984375*eccentricity_7 + 3.8203125*eccentricity_5 - 7.6875*eccentricity_3 + 3.5*eccentricity;
    double common_term_15 = -0.01352521839175485*eccentricity_14 - 0.11391658399470899*eccentricity_12 + 0.70339988425925926*eccentricity_10 - 3.9527777777777778*eccentricity_8 + 12.520833333333333*eccentricity_6 - 19.166666666666667*eccentricity_4 + 8.5*eccentricity_2;
    double common_term_16 = -0.61309564996648718*eccentricity_13 + 3.4192529304948433*eccentricity_11 - 13.840795446325231*eccentricity_9 + 33.890787760416667*eccentricity_7 - 42.350260416666667*eccentricity_5 + 17.604166666666667*eccentricity_3;
    double common_term_17 = -2.8222328404017857*eccentricity_14 + 12.780824497767857*eccentricity_12 - 40.741350446428571*eccentricity_10 + 81.34921875*eccentricity_8 - 86.41875*eccentricity_6 + 33.3125*eccentricity_4;
    double common_term_18 = 40.427591214454714*eccentricity_13 - 106.44385196544506*eccentricity_11 + 179.69544406467014*eccentricity_9 - 166.61648220486111*eccentricity_7 + 59.465364583333333*eccentricity_5;
    double common_term_19 = 113.64020786256246*eccentricity_14 - 254.98609090241035*eccentricity_12 + 373.3750496031746*eccentricity_10 - 307.8281746031746*eccentricity_8 + 101.90138888888889*eccentricity_6;
    double common_term_20 = -571.80859230995178*eccentricity_13 + 740.11936520167759*eccentricity_11 - 550.10469491141183*eccentricity_9 + 169.42345145089286*eccentricity_7;
    double common_term_21 = -1217.4109088636134*eccentricity_14 + 1413.1763491891259*eccentricity_12 - 957.14070422729277*eccentricity_10 + 275.21050347222222*eccentricity_8;
    double common_term_22 = 2617.1604006267897*eccentricity_13 - 1629.2002936322517*eccentricity_11 + 438.87348632274271*eccentricity_9;
    double common_term_23 = 4725.3268089234984*eccentricity_14 - 2722.714940137987*eccentricity_12 + 689.43969866071429*eccentricity_10;
    double common_term_24 = -4479.8494978454205*eccentricity_13 + 1069.6725544226351*eccentricity_11;
    double common_term_25 = -7272.9041147830075*eccentricity_14 + 1642.3080419026784*eccentricity_12;
    double common_term_26 = 2499.0369076334632*eccentricity_13;
    double common_term_27 = 3773.4130548187781*eccentricity_14;
    double common_term_28 = 150.7868458226344*eccentricity_14;
    double common_term_29 = 108.13325048320126*eccentricity_13;
    double common_term_30 = -153.88210445316402*eccentricity_14 + 77.418998706371753*eccentricity_12;
    double common_term_31 = -95.55095485052314*eccentricity_13 + 55.325853568183051*eccentricity_11;
    double common_term_32 = 54.18887844333401*eccentricity_14 - 57.801470339556277*eccentricity_12 + 39.452758143187831*eccentricity_10;
    double common_term_33 = 29.828531830899127*eccentricity_13 - 33.715195541381836*eccentricity_11 + 28.063401140485491*eccentricity_9;
    double common_term_34 = -2.2396417157061688*eccentricity_14 + 16.340335751488095*eccentricity_12 - 18.62028976521164*eccentricity_10 + 19.903013392857143*eccentricity_8;
    double common_term_35 = 0.93341902403183925*eccentricity_13 + 9.1335430992974175*eccentricity_11 - 9.3806243896484375*eccentricity_9 + 14.065462239583333*eccentricity_7;
    double common_term_36 = 3.2881651088169643*eccentricity_14 + 2.2009207589285714*eccentricity_12 + 5.4449497767857143*eccentricity_10 - 3.9053571428571429*eccentricity_8 + 9.896875*eccentricity_6;
    double common_term_37 = 3.0756571105758559*eccentricity_13 + 2.6084677378336589*eccentricity_11 + 3.6474173409598214*eccentricity_9 - 0.81168619791666667*eccentricity_7 + 6.92578125*eccentricity_5;
    double common_term_38 = 3.1441078404017857*eccentricity_14 + 2.9414132254464286*eccentricity_12 + 2.6606894841269841*eccentricity_10 + 2.809375*eccentricity_8 + 0.80625*eccentricity_6 + 4.8125*eccentricity_4;
    double common_term_39 = 3.0407588257108416*eccentricity_13 + 2.8249900817871094*eccentricity_11 + 2.58248291015625*eccentricity_9 + 2.41728515625*eccentricity_7 + 1.53515625*eccentricity_5 + 3.3125*eccentricity_3;
    double common_term_40 = 3.1427040318080357*eccentricity_14 + 2.9334933035714286*eccentricity_12 + 2.7083767361111111*eccentricity_10 + 2.4625*eccentricity_8 + 2.203125*eccentricity_6 + 1.75*eccentricity_4 + 2.25*eccentricity_2;
    double common_term_41 = 3.0393375761042196*eccentricity_13 + 2.822341054280599*eccentricity_11 + 2.587322998046875*eccentricity_9 + 2.3289388020833333*eccentricity_7 + 2.0390625*eccentricity_5 + 1.6875*eccentricity_3 + 1.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(29);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (2, 0).
    // q = -14
    result_by_lpq.set(c_Key3(2, 0, -14), common_term_0);
    result_by_q.set(c_Key1(-14), common_term_0);
    // q = -13
    result_by_lpq.set(c_Key3(2, 0, -13), common_term_1);
    result_by_q.set(c_Key1(-13), common_term_1);
    // q = -12
    result_by_lpq.set(c_Key3(2, 0, -12), common_term_2);
    result_by_q.set(c_Key1(-12), common_term_2);
    // q = -11
    result_by_lpq.set(c_Key3(2, 0, -11), common_term_3);
    result_by_q.set(c_Key1(-11), common_term_3);
    // q = -10
    result_by_lpq.set(c_Key3(2, 0, -10), common_term_4);
    result_by_q.set(c_Key1(-10), common_term_4);
    // q = -9
    result_by_lpq.set(c_Key3(2, 0, -9), common_term_5);
    result_by_q.set(c_Key1(-9), common_term_5);
    // q = -8
    result_by_lpq.set(c_Key3(2, 0, -8), common_term_6);
    result_by_q.set(c_Key1(-8), common_term_6);
    // q = -7
    result_by_lpq.set(c_Key3(2, 0, -7), common_term_7);
    result_by_q.set(c_Key1(-7), common_term_7);
    // q = -6
    result_by_lpq.set(c_Key3(2, 0, -6), common_term_8);
    result_by_q.set(c_Key1(-6), common_term_8);
    // q = -5
    result_by_lpq.set(c_Key3(2, 0, -5), common_term_9);
    result_by_q.set(c_Key1(-5), common_term_9);
    // q = -4
    result_by_lpq.set(c_Key3(2, 0, -4), common_term_10);
    result_by_q.set(c_Key1(-4), common_term_10);
    // q = -3
    result_by_lpq.set(c_Key3(2, 0, -3), common_term_11);
    result_by_q.set(c_Key1(-3), common_term_11);
    // q = -1
    result_by_lpq.set(c_Key3(2, 0, -1), common_term_12);
    result_by_q.set(c_Key1(-1), common_term_12);
    // q = 0
    result_by_lpq.set(c_Key3(2, 0, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(2, 0, 1), common_term_14);
    result_by_q.set(c_Key1(1), common_term_14);
    // q = 2
    result_by_lpq.set(c_Key3(2, 0, 2), common_term_15);
    result_by_q.set(c_Key1(2), common_term_15);
    // q = 3
    result_by_lpq.set(c_Key3(2, 0, 3), common_term_16);
    result_by_q.set(c_Key1(3), common_term_16);
    // q = 4
    result_by_lpq.set(c_Key3(2, 0, 4), common_term_17);
    result_by_q.set(c_Key1(4), common_term_17);
    // q = 5
    result_by_lpq.set(c_Key3(2, 0, 5), common_term_18);
    result_by_q.set(c_Key1(5), common_term_18);
    // q = 6
    result_by_lpq.set(c_Key3(2, 0, 6), common_term_19);
    result_by_q.set(c_Key1(6), common_term_19);
    // q = 7
    result_by_lpq.set(c_Key3(2, 0, 7), common_term_20);
    result_by_q.set(c_Key1(7), common_term_20);
    // q = 8
    result_by_lpq.set(c_Key3(2, 0, 8), common_term_21);
    result_by_q.set(c_Key1(8), common_term_21);
    // q = 9
    result_by_lpq.set(c_Key3(2, 0, 9), common_term_22);
    result_by_q.set(c_Key1(9), common_term_22);
    // q = 10
    result_by_lpq.set(c_Key3(2, 0, 10), common_term_23);
    result_by_q.set(c_Key1(10), common_term_23);
    // q = 11
    result_by_lpq.set(c_Key3(2, 0, 11), common_term_24);
    result_by_q.set(c_Key1(11), common_term_24);
    // q = 12
    result_by_lpq.set(c_Key3(2, 0, 12), common_term_25);
    result_by_q.set(c_Key1(12), common_term_25);
    // q = 13
    result_by_lpq.set(c_Key3(2, 0, 13), common_term_26);
    result_by_q.set(c_Key1(13), common_term_26);
    // q = 14
    result_by_lpq.set(c_Key3(2, 0, 14), common_term_27);
    result_by_q.set(c_Key1(14), common_term_27);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 0), result_by_q);
    result_by_q.clear();

    // l , p = (2, 1).
    // q = -14
    result_by_lpq.set(c_Key3(2, 1, -14), common_term_28);
    result_by_q.set(c_Key1(-14), common_term_28);
    // q = -13
    result_by_lpq.set(c_Key3(2, 1, -13), common_term_29);
    result_by_q.set(c_Key1(-13), common_term_29);
    // q = -12
    result_by_lpq.set(c_Key3(2, 1, -12), common_term_30);
    result_by_q.set(c_Key1(-12), common_term_30);
    // q = -11
    result_by_lpq.set(c_Key3(2, 1, -11), common_term_31);
    result_by_q.set(c_Key1(-11), common_term_31);
    // q = -10
    result_by_lpq.set(c_Key3(2, 1, -10), common_term_32);
    result_by_q.set(c_Key1(-10), common_term_32);
    // q = -9
    result_by_lpq.set(c_Key3(2, 1, -9), common_term_33);
    result_by_q.set(c_Key1(-9), common_term_33);
    // q = -8
    result_by_lpq.set(c_Key3(2, 1, -8), common_term_34);
    result_by_q.set(c_Key1(-8), common_term_34);
    // q = -7
    result_by_lpq.set(c_Key3(2, 1, -7), common_term_35);
    result_by_q.set(c_Key1(-7), common_term_35);
    // q = -6
    result_by_lpq.set(c_Key3(2, 1, -6), common_term_36);
    result_by_q.set(c_Key1(-6), common_term_36);
    // q = -5
    result_by_lpq.set(c_Key3(2, 1, -5), common_term_37);
    result_by_q.set(c_Key1(-5), common_term_37);
    // q = -4
    result_by_lpq.set(c_Key3(2, 1, -4), common_term_38);
    result_by_q.set(c_Key1(-4), common_term_38);
    // q = -3
    result_by_lpq.set(c_Key3(2, 1, -3), common_term_39);
    result_by_q.set(c_Key1(-3), common_term_39);
    // q = -2
    result_by_lpq.set(c_Key3(2, 1, -2), common_term_40);
    result_by_q.set(c_Key1(-2), common_term_40);
    // q = -1
    result_by_lpq.set(c_Key3(2, 1, -1), common_term_41);
    result_by_q.set(c_Key1(-1), common_term_41);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -1.5);
    result_by_lpq.set(c_Key3(2, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(2, 1, 1), common_term_41);
    result_by_q.set(c_Key1(1), common_term_41);
    // q = 2
    result_by_lpq.set(c_Key3(2, 1, 2), common_term_40);
    result_by_q.set(c_Key1(2), common_term_40);
    // q = 3
    result_by_lpq.set(c_Key3(2, 1, 3), common_term_39);
    result_by_q.set(c_Key1(3), common_term_39);
    // q = 4
    result_by_lpq.set(c_Key3(2, 1, 4), common_term_38);
    result_by_q.set(c_Key1(4), common_term_38);
    // q = 5
    result_by_lpq.set(c_Key3(2, 1, 5), common_term_37);
    result_by_q.set(c_Key1(5), common_term_37);
    // q = 6
    result_by_lpq.set(c_Key3(2, 1, 6), common_term_36);
    result_by_q.set(c_Key1(6), common_term_36);
    // q = 7
    result_by_lpq.set(c_Key3(2, 1, 7), common_term_35);
    result_by_q.set(c_Key1(7), common_term_35);
    // q = 8
    result_by_lpq.set(c_Key3(2, 1, 8), common_term_34);
    result_by_q.set(c_Key1(8), common_term_34);
    // q = 9
    result_by_lpq.set(c_Key3(2, 1, 9), common_term_33);
    result_by_q.set(c_Key1(9), common_term_33);
    // q = 10
    result_by_lpq.set(c_Key3(2, 1, 10), common_term_32);
    result_by_q.set(c_Key1(10), common_term_32);
    // q = 11
    result_by_lpq.set(c_Key3(2, 1, 11), common_term_31);
    result_by_q.set(c_Key1(11), common_term_31);
    // q = 12
    result_by_lpq.set(c_Key3(2, 1, 12), common_term_30);
    result_by_q.set(c_Key1(12), common_term_30);
    // q = 13
    result_by_lpq.set(c_Key3(2, 1, 13), common_term_29);
    result_by_q.set(c_Key1(13), common_term_29);
    // q = 14
    result_by_lpq.set(c_Key3(2, 1, 14), common_term_28);
    result_by_q.set(c_Key1(14), common_term_28);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 1), result_by_q);
    result_by_q.clear();

    // l , p = (2, 2).
    // q = -14
    result_by_lpq.set(c_Key3(2, 2, -14), common_term_27);
    result_by_q.set(c_Key1(-14), common_term_27);
    // q = -13
    result_by_lpq.set(c_Key3(2, 2, -13), common_term_26);
    result_by_q.set(c_Key1(-13), common_term_26);
    // q = -12
    result_by_lpq.set(c_Key3(2, 2, -12), common_term_25);
    result_by_q.set(c_Key1(-12), common_term_25);
    // q = -11
    result_by_lpq.set(c_Key3(2, 2, -11), common_term_24);
    result_by_q.set(c_Key1(-11), common_term_24);
    // q = -10
    result_by_lpq.set(c_Key3(2, 2, -10), common_term_23);
    result_by_q.set(c_Key1(-10), common_term_23);
    // q = -9
    result_by_lpq.set(c_Key3(2, 2, -9), common_term_22);
    result_by_q.set(c_Key1(-9), common_term_22);
    // q = -8
    result_by_lpq.set(c_Key3(2, 2, -8), common_term_21);
    result_by_q.set(c_Key1(-8), common_term_21);
    // q = -7
    result_by_lpq.set(c_Key3(2, 2, -7), common_term_20);
    result_by_q.set(c_Key1(-7), common_term_20);
    // q = -6
    result_by_lpq.set(c_Key3(2, 2, -6), common_term_19);
    result_by_q.set(c_Key1(-6), common_term_19);
    // q = -5
    result_by_lpq.set(c_Key3(2, 2, -5), common_term_18);
    result_by_q.set(c_Key1(-5), common_term_18);
    // q = -4
    result_by_lpq.set(c_Key3(2, 2, -4), common_term_17);
    result_by_q.set(c_Key1(-4), common_term_17);
    // q = -3
    result_by_lpq.set(c_Key3(2, 2, -3), common_term_16);
    result_by_q.set(c_Key1(-3), common_term_16);
    // q = -2
    result_by_lpq.set(c_Key3(2, 2, -2), common_term_15);
    result_by_q.set(c_Key1(-2), common_term_15);
    // q = -1
    result_by_lpq.set(c_Key3(2, 2, -1), common_term_14);
    result_by_q.set(c_Key1(-1), common_term_14);
    // q = 0
    result_by_lpq.set(c_Key3(2, 2, 0), common_term_13);
    result_by_q.set(c_Key1(0), common_term_13);
    // q = 1
    result_by_lpq.set(c_Key3(2, 2, 1), common_term_12);
    result_by_q.set(c_Key1(1), common_term_12);
    // q = 3
    result_by_lpq.set(c_Key3(2, 2, 3), common_term_11);
    result_by_q.set(c_Key1(3), common_term_11);
    // q = 4
    result_by_lpq.set(c_Key3(2, 2, 4), common_term_10);
    result_by_q.set(c_Key1(4), common_term_10);
    // q = 5
    result_by_lpq.set(c_Key3(2, 2, 5), common_term_9);
    result_by_q.set(c_Key1(5), common_term_9);
    // q = 6
    result_by_lpq.set(c_Key3(2, 2, 6), common_term_8);
    result_by_q.set(c_Key1(6), common_term_8);
    // q = 7
    result_by_lpq.set(c_Key3(2, 2, 7), common_term_7);
    result_by_q.set(c_Key1(7), common_term_7);
    // q = 8
    result_by_lpq.set(c_Key3(2, 2, 8), common_term_6);
    result_by_q.set(c_Key1(8), common_term_6);
    // q = 9
    result_by_lpq.set(c_Key3(2, 2, 9), common_term_5);
    result_by_q.set(c_Key1(9), common_term_5);
    // q = 10
    result_by_lpq.set(c_Key3(2, 2, 10), common_term_4);
    result_by_q.set(c_Key1(10), common_term_4);
    // q = 11
    result_by_lpq.set(c_Key3(2, 2, 11), common_term_3);
    result_by_q.set(c_Key1(11), common_term_3);
    // q = 12
    result_by_lpq.set(c_Key3(2, 2, 12), common_term_2);
    result_by_q.set(c_Key1(12), common_term_2);
    // q = 13
    result_by_lpq.set(c_Key3(2, 2, 13), common_term_1);
    result_by_q.set(c_Key1(13), common_term_1);
    // q = 14
    result_by_lpq.set(c_Key3(2, 2, 14), common_term_0);
    result_by_q.set(c_Key1(14), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 2), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}


EccentricityFuncOutput c_eccentricity_function_l2_e20(double eccentricity)
{
    // Eccentricity Functions Calculated for l = 2.
    // Functions are taylor expanded and truncated at eccentricity^20.

    // Functions Calculated for l = 2.

    c_IntMap<c_Key3, double> result_by_lpq(115);
    c_IntMap<c_Key2, c_IntMap<c_Key1, double>> result_by_lp(3);
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
    double common_term_0 = 3.7485643220704711*eccentricity_19;
    double common_term_1 = 2.8137061873419034*eccentricity_18;
    double common_term_2 = -3.7864321270486251*eccentricity_19 + 2.1133574662596977*eccentricity_17;
    double common_term_3 = -2.4759735389262424*eccentricity_18 + 1.5883603834621178*eccentricity_16;
    double common_term_4 = 1.4955082289759112*eccentricity_19 - 1.5864985892412163*eccentricity_17 + 1.1945401142522099*eccentricity_15;
    double common_term_5 = 0.88765919794491223*eccentricity_18 - 0.98878492935635793*eccentricity_16 + 0.89889539032396175*eccentricity_14;
    double common_term_6 = -0.12564674885432872*eccentricity_19 + 0.524992499702779*eccentricity_17 - 0.59216495266685987*eccentricity_15 + 0.67675994590498271*eccentricity_13;
    double common_term_7 = -0.022889993571858531*eccentricity_18 + 0.31330314056299144*eccentricity_16 - 0.3332565249340423*eccentricity_14 + 0.50968644989912351*eccentricity_12;
    double common_term_8 = 0.06735584280299031*eccentricity_19 + 0.024910391501840978*eccentricity_17 + 0.19262608014322542*eccentricity_15 - 0.16794225909493186*eccentricity_13 + 0.38386802078841569*eccentricity_11;
    double common_term_9 = 0.056412331544606677*eccentricity_18 + 0.044237792385940534*eccentricity_16 + 0.12532520976965421*eccentricity_14 - 0.065672599005932339*eccentricity_12 + 0.28895943562610229*eccentricity_10;
    double common_term_10 = 0.042337346657042359*eccentricity_19 + 0.049193112797624476*eccentricity_17 + 0.049362019568417848*eccentricity_15 + 0.088173672850552098*eccentricity_13 - 0.0054298692868079668*eccentricity_11 + 0.21719477147231867*eccentricity_9;
    double common_term_11 = 0.037283743112356394*eccentricity_18 + 0.043023031655844156*eccentricity_16 + 0.047615031452922078*eccentricity_14 + 0.067123325892857143*eccentricity_12 + 0.027120535714285714*eccentricity_10 + 0.16272321428571429*eccentricity_8;
    double common_term_12 = 0.027554583037287145*eccentricity_19 + 0.031655845142880503*eccentricity_17 + 0.036767070756508509*eccentricity_15 + 0.042511234558168262*eccentricity_13 + 0.05387545052449329*eccentricity_11 + 0.041628640795510913*eccentricity_9 + 0.12110150049603175*eccentricity_7;
    double common_term_13 = 0.021943306846759779*eccentricity_18 + 0.025485055582277805*eccentricity_16 + 0.029965461493239271*eccentricity_14 + 0.035552616108171664*eccentricity_12 + 0.043650793650793651*eccentricity_10 + 0.044444444444444444*eccentricity_8 + 0.088888888888888889*eccentricity_6;
    double common_term_14 = 0.01378362270957291*eccentricity_19 + 0.015984563621414172*eccentricity_17 + 0.018796816400119237*eccentricity_15 + 0.022452848298209054*eccentricity_13 + 0.027257374354771205*eccentricity_11 + 0.033759416852678571*eccentricity_9 + 0.03955078125*eccentricity_7 + 0.06328125*eccentricity_5;
    double common_term_15 = 0.0083823051245946728*eccentricity_18 + 0.0098275931820130315*eccentricity_16 + 0.01172185593952087*eccentricity_14 + 0.014269696593915344*eccentricity_12 + 0.017786871693121693*eccentricity_10 + 0.022743055555555556*eccentricity_8 + 0.029166666666666667*eccentricity_6 + 0.041666666666666667*eccentricity_4;
    double common_term_16 = 0.00282447998447232*eccentricity_19 + 0.0032938900851694608*eccentricity_17 + 0.0039075560377505207*eccentricity_15 + 0.0047346311649948201*eccentricity_13 + 0.0058915435952484292*eccentricity_11 + 0.0075841833043981482*eccentricity_9 + 0.010188802083333333*eccentricity_7 + 0.014322916666666667*eccentricity_5 + 0.020833333333333333*eccentricity_3;
    double common_term_17 = -0.0025515931232383317*eccentricity_19 - 0.0029421149027509295*eccentricity_17 - 0.003440899117132914*eccentricity_15 - 0.0040929340517289634*eccentricity_13 - 0.0049673518428096065*eccentricity_11 - 0.0061692979600694444*eccentricity_9 - 0.0077582465277777778*eccentricity_7 - 0.013020833333333333*eccentricity_5 + 0.0625*eccentricity_3 - 0.5*eccentricity;
    double common_term_18 = -0.007150322282650879*eccentricity_18 - 0.0082195286788036974*eccentricity_16 - 0.0095563234697814311*eccentricity_14 - 0.011227454668209877*eccentricity_12 - 0.013611111111111111*eccentricity_10 - 0.0086805555555555556*eccentricity_8 - 0.12152777777777778*eccentricity_6 + 0.8125*eccentricity_4 - 2.5*eccentricity_2 + 1.0;
    double common_term_19 = -0.011342737497068085*eccentricity_19 - 0.012844661894273393*eccentricity_17 - 0.01466690719127655*eccentricity_15 - 0.016310808999197824*eccentricity_13 - 0.027522125244140625*eccentricity_11 + 0.082562255859375*eccentricity_9 - 0.86083984375*eccentricity_7 + 3.8203125*eccentricity_5 - 7.6875*eccentricity_3 + 3.5*eccentricity;
    double common_term_20 = -0.017039937960901632*eccentricity_18 - 0.019675649259189272*eccentricity_16 - 0.01352521839175485*eccentricity_14 - 0.11391658399470899*eccentricity_12 + 0.70339988425925926*eccentricity_10 - 3.9527777777777778*eccentricity_8 + 12.520833333333333*eccentricity_6 - 19.166666666666667*eccentricity_4 + 8.5*eccentricity_2;
    double common_term_21 = -0.020351987519661285*eccentricity_19 - 0.02986415064757071*eccentricity_17 + 0.046774400199204829*eccentricity_15 - 0.61309564996648718*eccentricity_13 + 3.4192529304948433*eccentricity_11 - 13.840795446325231*eccentricity_9 + 33.890787760416667*eccentricity_7 - 42.350260416666667*eccentricity_5 + 17.604166666666667*eccentricity_3;
    double common_term_22 = -0.081789837499275278*eccentricity_18 + 0.41899421037946429*eccentricity_16 - 2.8222328404017857*eccentricity_14 + 12.780824497767857*eccentricity_12 - 40.741350446428571*eccentricity_10 + 81.34921875*eccentricity_8 - 86.41875*eccentricity_6 + 33.3125*eccentricity_4;
    double common_term_23 = -0.36123893873221676*eccentricity_19 + 2.1174997908480243*eccentricity_17 - 10.84046902853226*eccentricity_15 + 40.427591214454714*eccentricity_13 - 106.44385196544506*eccentricity_11 + 179.69544406467014*eccentricity_9 - 166.61648220486111*eccentricity_7 + 59.465364583333333*eccentricity_5;
    double common_term_24 = 8.553605499102485*eccentricity_18 - 36.074829145455882*eccentricity_16 + 113.64020786256246*eccentricity_14 - 254.98609090241035*eccentricity_12 + 373.3750496031746*eccentricity_10 - 307.8281746031746*eccentricity_8 + 101.90138888888889*eccentricity_6;
    double common_term_25 = 29.877995840631622*eccentricity_19 - 107.48655336883867*eccentricity_17 + 292.40192527141664*eccentricity_15 - 571.80859230995178*eccentricity_13 + 740.11936520167759*eccentricity_11 - 550.10469491141183*eccentricity_9 + 169.42345145089286*eccentricity_7;
    double common_term_26 = -293.67666566579011*eccentricity_18 + 702.28706134106163*eccentricity_16 - 1217.4109088636134*eccentricity_14 + 1413.1763491891259*eccentricity_12 - 957.14070422729277*eccentricity_10 + 275.21050347222222*eccentricity_8;
    double common_term_27 = -748.53426155212087*eccentricity_19 + 1596.1211457416656*eccentricity_17 - 2485.4608525285744*eccentricity_15 + 2617.1604006267897*eccentricity_13 - 1629.2002936322517*eccentricity_11 + 438.87348632274271*eccentricity_9;
    double common_term_28 = 3466.9457044305008*eccentricity_18 - 4901.6308334243881*eccentricity_16 + 4725.3268089234984*eccentricity_14 - 2722.714940137987*eccentricity_12 + 689.43969866071429*eccentricity_10;
    double common_term_29 = 7250.9703331544839*eccentricity_19 - 9389.5376448926914*eccentricity_17 + 8350.3215027223972*eccentricity_15 - 4479.8494978454205*eccentricity_13 + 1069.6725544226351*eccentricity_11;
    double common_term_30 = -17546.191875121375*eccentricity_18 + 14486.914630738745*eccentricity_16 - 7272.9041147830075*eccentricity_14 + 1642.3080419026784*eccentricity_12;
    double common_term_31 = -32094.641659543963*eccentricity_19 + 24735.01155327285*eccentricity_17 - 11670.721326916038*eccentricity_15 + 2499.0369076334632*eccentricity_13;
    double common_term_32 = 41646.006357703483*eccentricity_18 - 18537.628592581602*eccentricity_16 + 3773.4130548187781*eccentricity_14;
    double common_term_33 = 69257.765582481734*eccentricity_19 - 29180.344549111437*eccentricity_17 + 5659.3667191809489*eccentricity_15;
    double common_term_34 = -45565.400776266436*eccentricity_18 + 8437.7060119717232*eccentricity_16;
    double common_term_35 = -70639.971358844223*eccentricity_19 + 12513.987310315178*eccentricity_17;
    double common_term_36 = 18472.54310597884*eccentricity_18;
    double common_term_37 = 27153.443004169995*eccentricity_19;
    double common_term_38 = 780.24283889834832*eccentricity_19;
    double common_term_39 = 562.81178606797004*eccentricity_18;
    double common_term_40 = -1329.0451794324829*eccentricity_19 + 405.59164111798766*eccentricity_17;
    double common_term_41 = -881.91037438072872*eccentricity_18 + 291.98809610570637*eccentricity_16;
    double common_term_42 = 824.31080878400255*eccentricity_19 - 580.19774204246337*eccentricity_17 + 209.96289389476391*eccentricity_15;
    double common_term_43 = 495.22082633922241*eccentricity_18 - 377.81499005145045*eccentricity_16 + 150.7868458226344*eccentricity_14;
    double common_term_44 = -216.46744861058754*eccentricity_19 + 292.78830576776268*eccentricity_17 - 242.98749348695164*eccentricity_15 + 108.13325048320126*eccentricity_13;
    double common_term_45 = -110.91929798082777*eccentricity_18 + 170.06222147899534*eccentricity_16 - 153.88210445316402*eccentricity_14 + 77.418998706371753*eccentricity_12;
    double common_term_46 = 28.952513424707026*eccentricity_19 - 54.01649519891636*eccentricity_17 + 96.915003617084701*eccentricity_15 - 95.55095485052314*eccentricity_13 + 55.325853568183051*eccentricity_11;
    double common_term_47 = 14.260448209702564*eccentricity_18 - 24.294051109982642*eccentricity_16 + 54.18887844333401*eccentricity_14 - 57.801470339556277*eccentricity_12 + 39.452758143187831*eccentricity_10;
    double common_term_48 = 2.5353443481475767*eccentricity_19 + 7.6855999352900938*eccentricity_17 - 9.36642109506316*eccentricity_15 + 29.828531830899127*eccentricity_13 - 33.715195541381836*eccentricity_11 + 28.063401140485491*eccentricity_9;
    double common_term_49 = 3.1958957893024462*eccentricity_18 + 4.8878596534176504*eccentricity_16 - 2.2396417157061688*eccentricity_14 + 16.340335751488095*eccentricity_12 - 18.62028976521164*eccentricity_10 + 19.903013392857143*eccentricity_8;
    double common_term_50 = 3.629201282411903*eccentricity_19 + 3.3476530000734426*eccentricity_17 + 3.7496673280883718*eccentricity_15 + 0.93341902403183925*eccentricity_13 + 9.1335430992974175*eccentricity_11 - 9.3806243896484375*eccentricity_9 + 14.065462239583333*eccentricity_7;
    double common_term_51 = 3.5291701105291193*eccentricity_18 + 3.3224313996550325*eccentricity_16 + 3.2881651088169643*eccentricity_14 + 2.2009207589285714*eccentricity_12 + 5.4449497767857143*eccentricity_10 - 3.9053571428571429*eccentricity_8 + 9.896875*eccentricity_6;
    double common_term_52 = 3.6149488242532871*eccentricity_19 + 3.4351517619436414*eccentricity_17 + 3.2415165339039747*eccentricity_15 + 3.0756571105758559*eccentricity_13 + 2.6084677378336589*eccentricity_11 + 3.6474173409598214*eccentricity_9 - 0.81168619791666667*eccentricity_7 + 6.92578125*eccentricity_5;
    double common_term_53 = 3.5252937859967482*eccentricity_18 + 3.3402773052845753*eccentricity_16 + 3.1441078404017857*eccentricity_14 + 2.9414132254464286*eccentricity_12 + 2.6606894841269841*eccentricity_10 + 2.809375*eccentricity_8 + 0.80625*eccentricity_6 + 4.8125*eccentricity_4;
    double common_term_54 = 3.6137849946649657*eccentricity_19 + 3.4333311629754548*eccentricity_17 + 3.2429372967141015*eccentricity_15 + 3.0407588257108416*eccentricity_13 + 2.8249900817871094*eccentricity_11 + 2.58248291015625*eccentricity_9 + 2.41728515625*eccentricity_7 + 1.53515625*eccentricity_5 + 3.3125*eccentricity_3;
    double common_term_55 = 3.524281538813623*eccentricity_18 + 3.3389199888819339*eccentricity_16 + 3.1427040318080357*eccentricity_14 + 2.9334933035714286*eccentricity_12 + 2.7083767361111111*eccentricity_10 + 2.4625*eccentricity_8 + 2.203125*eccentricity_6 + 1.75*eccentricity_4 + 2.25*eccentricity_2;
    double common_term_56 = 3.6131874865348258*eccentricity_19 + 3.432554495346766*eccentricity_17 + 3.2418959114597498*eccentricity_15 + 3.0393375761042196*eccentricity_13 + 2.822341054280599*eccentricity_11 + 2.587322998046875*eccentricity_9 + 2.3289388020833333*eccentricity_7 + 2.0390625*eccentricity_5 + 1.6875*eccentricity_3 + 1.5*eccentricity;

    c_IntMap<c_Key1, double> result_by_q(39);
    double tmp_double;

    // eccentricity function by mode:
    // l , p = (2, 0).
    // q = -19
    result_by_lpq.set(c_Key3(2, 0, -19), common_term_0);
    result_by_q.set(c_Key1(-19), common_term_0);
    // q = -18
    result_by_lpq.set(c_Key3(2, 0, -18), common_term_1);
    result_by_q.set(c_Key1(-18), common_term_1);
    // q = -17
    result_by_lpq.set(c_Key3(2, 0, -17), common_term_2);
    result_by_q.set(c_Key1(-17), common_term_2);
    // q = -16
    result_by_lpq.set(c_Key3(2, 0, -16), common_term_3);
    result_by_q.set(c_Key1(-16), common_term_3);
    // q = -15
    result_by_lpq.set(c_Key3(2, 0, -15), common_term_4);
    result_by_q.set(c_Key1(-15), common_term_4);
    // q = -14
    result_by_lpq.set(c_Key3(2, 0, -14), common_term_5);
    result_by_q.set(c_Key1(-14), common_term_5);
    // q = -13
    result_by_lpq.set(c_Key3(2, 0, -13), common_term_6);
    result_by_q.set(c_Key1(-13), common_term_6);
    // q = -12
    result_by_lpq.set(c_Key3(2, 0, -12), common_term_7);
    result_by_q.set(c_Key1(-12), common_term_7);
    // q = -11
    result_by_lpq.set(c_Key3(2, 0, -11), common_term_8);
    result_by_q.set(c_Key1(-11), common_term_8);
    // q = -10
    result_by_lpq.set(c_Key3(2, 0, -10), common_term_9);
    result_by_q.set(c_Key1(-10), common_term_9);
    // q = -9
    result_by_lpq.set(c_Key3(2, 0, -9), common_term_10);
    result_by_q.set(c_Key1(-9), common_term_10);
    // q = -8
    result_by_lpq.set(c_Key3(2, 0, -8), common_term_11);
    result_by_q.set(c_Key1(-8), common_term_11);
    // q = -7
    result_by_lpq.set(c_Key3(2, 0, -7), common_term_12);
    result_by_q.set(c_Key1(-7), common_term_12);
    // q = -6
    result_by_lpq.set(c_Key3(2, 0, -6), common_term_13);
    result_by_q.set(c_Key1(-6), common_term_13);
    // q = -5
    result_by_lpq.set(c_Key3(2, 0, -5), common_term_14);
    result_by_q.set(c_Key1(-5), common_term_14);
    // q = -4
    result_by_lpq.set(c_Key3(2, 0, -4), common_term_15);
    result_by_q.set(c_Key1(-4), common_term_15);
    // q = -3
    result_by_lpq.set(c_Key3(2, 0, -3), common_term_16);
    result_by_q.set(c_Key1(-3), common_term_16);
    // q = -1
    result_by_lpq.set(c_Key3(2, 0, -1), common_term_17);
    result_by_q.set(c_Key1(-1), common_term_17);
    // q = 0
    result_by_lpq.set(c_Key3(2, 0, 0), common_term_18);
    result_by_q.set(c_Key1(0), common_term_18);
    // q = 1
    result_by_lpq.set(c_Key3(2, 0, 1), common_term_19);
    result_by_q.set(c_Key1(1), common_term_19);
    // q = 2
    result_by_lpq.set(c_Key3(2, 0, 2), common_term_20);
    result_by_q.set(c_Key1(2), common_term_20);
    // q = 3
    result_by_lpq.set(c_Key3(2, 0, 3), common_term_21);
    result_by_q.set(c_Key1(3), common_term_21);
    // q = 4
    result_by_lpq.set(c_Key3(2, 0, 4), common_term_22);
    result_by_q.set(c_Key1(4), common_term_22);
    // q = 5
    result_by_lpq.set(c_Key3(2, 0, 5), common_term_23);
    result_by_q.set(c_Key1(5), common_term_23);
    // q = 6
    result_by_lpq.set(c_Key3(2, 0, 6), common_term_24);
    result_by_q.set(c_Key1(6), common_term_24);
    // q = 7
    result_by_lpq.set(c_Key3(2, 0, 7), common_term_25);
    result_by_q.set(c_Key1(7), common_term_25);
    // q = 8
    result_by_lpq.set(c_Key3(2, 0, 8), common_term_26);
    result_by_q.set(c_Key1(8), common_term_26);
    // q = 9
    result_by_lpq.set(c_Key3(2, 0, 9), common_term_27);
    result_by_q.set(c_Key1(9), common_term_27);
    // q = 10
    result_by_lpq.set(c_Key3(2, 0, 10), common_term_28);
    result_by_q.set(c_Key1(10), common_term_28);
    // q = 11
    result_by_lpq.set(c_Key3(2, 0, 11), common_term_29);
    result_by_q.set(c_Key1(11), common_term_29);
    // q = 12
    result_by_lpq.set(c_Key3(2, 0, 12), common_term_30);
    result_by_q.set(c_Key1(12), common_term_30);
    // q = 13
    result_by_lpq.set(c_Key3(2, 0, 13), common_term_31);
    result_by_q.set(c_Key1(13), common_term_31);
    // q = 14
    result_by_lpq.set(c_Key3(2, 0, 14), common_term_32);
    result_by_q.set(c_Key1(14), common_term_32);
    // q = 15
    result_by_lpq.set(c_Key3(2, 0, 15), common_term_33);
    result_by_q.set(c_Key1(15), common_term_33);
    // q = 16
    result_by_lpq.set(c_Key3(2, 0, 16), common_term_34);
    result_by_q.set(c_Key1(16), common_term_34);
    // q = 17
    result_by_lpq.set(c_Key3(2, 0, 17), common_term_35);
    result_by_q.set(c_Key1(17), common_term_35);
    // q = 18
    result_by_lpq.set(c_Key3(2, 0, 18), common_term_36);
    result_by_q.set(c_Key1(18), common_term_36);
    // q = 19
    result_by_lpq.set(c_Key3(2, 0, 19), common_term_37);
    result_by_q.set(c_Key1(19), common_term_37);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 0), result_by_q);
    result_by_q.clear();

    // l , p = (2, 1).
    // q = -19
    result_by_lpq.set(c_Key3(2, 1, -19), common_term_38);
    result_by_q.set(c_Key1(-19), common_term_38);
    // q = -18
    result_by_lpq.set(c_Key3(2, 1, -18), common_term_39);
    result_by_q.set(c_Key1(-18), common_term_39);
    // q = -17
    result_by_lpq.set(c_Key3(2, 1, -17), common_term_40);
    result_by_q.set(c_Key1(-17), common_term_40);
    // q = -16
    result_by_lpq.set(c_Key3(2, 1, -16), common_term_41);
    result_by_q.set(c_Key1(-16), common_term_41);
    // q = -15
    result_by_lpq.set(c_Key3(2, 1, -15), common_term_42);
    result_by_q.set(c_Key1(-15), common_term_42);
    // q = -14
    result_by_lpq.set(c_Key3(2, 1, -14), common_term_43);
    result_by_q.set(c_Key1(-14), common_term_43);
    // q = -13
    result_by_lpq.set(c_Key3(2, 1, -13), common_term_44);
    result_by_q.set(c_Key1(-13), common_term_44);
    // q = -12
    result_by_lpq.set(c_Key3(2, 1, -12), common_term_45);
    result_by_q.set(c_Key1(-12), common_term_45);
    // q = -11
    result_by_lpq.set(c_Key3(2, 1, -11), common_term_46);
    result_by_q.set(c_Key1(-11), common_term_46);
    // q = -10
    result_by_lpq.set(c_Key3(2, 1, -10), common_term_47);
    result_by_q.set(c_Key1(-10), common_term_47);
    // q = -9
    result_by_lpq.set(c_Key3(2, 1, -9), common_term_48);
    result_by_q.set(c_Key1(-9), common_term_48);
    // q = -8
    result_by_lpq.set(c_Key3(2, 1, -8), common_term_49);
    result_by_q.set(c_Key1(-8), common_term_49);
    // q = -7
    result_by_lpq.set(c_Key3(2, 1, -7), common_term_50);
    result_by_q.set(c_Key1(-7), common_term_50);
    // q = -6
    result_by_lpq.set(c_Key3(2, 1, -6), common_term_51);
    result_by_q.set(c_Key1(-6), common_term_51);
    // q = -5
    result_by_lpq.set(c_Key3(2, 1, -5), common_term_52);
    result_by_q.set(c_Key1(-5), common_term_52);
    // q = -4
    result_by_lpq.set(c_Key3(2, 1, -4), common_term_53);
    result_by_q.set(c_Key1(-4), common_term_53);
    // q = -3
    result_by_lpq.set(c_Key3(2, 1, -3), common_term_54);
    result_by_q.set(c_Key1(-3), common_term_54);
    // q = -2
    result_by_lpq.set(c_Key3(2, 1, -2), common_term_55);
    result_by_q.set(c_Key1(-2), common_term_55);
    // q = -1
    result_by_lpq.set(c_Key3(2, 1, -1), common_term_56);
    result_by_q.set(c_Key1(-1), common_term_56);
    // q = 0
    tmp_double = std::pow(1.0 - eccentricity_2, -1.5);
    result_by_lpq.set(c_Key3(2, 1, 0), tmp_double);
    result_by_q.set(c_Key1(0), tmp_double);
    // q = 1
    result_by_lpq.set(c_Key3(2, 1, 1), common_term_56);
    result_by_q.set(c_Key1(1), common_term_56);
    // q = 2
    result_by_lpq.set(c_Key3(2, 1, 2), common_term_55);
    result_by_q.set(c_Key1(2), common_term_55);
    // q = 3
    result_by_lpq.set(c_Key3(2, 1, 3), common_term_54);
    result_by_q.set(c_Key1(3), common_term_54);
    // q = 4
    result_by_lpq.set(c_Key3(2, 1, 4), common_term_53);
    result_by_q.set(c_Key1(4), common_term_53);
    // q = 5
    result_by_lpq.set(c_Key3(2, 1, 5), common_term_52);
    result_by_q.set(c_Key1(5), common_term_52);
    // q = 6
    result_by_lpq.set(c_Key3(2, 1, 6), common_term_51);
    result_by_q.set(c_Key1(6), common_term_51);
    // q = 7
    result_by_lpq.set(c_Key3(2, 1, 7), common_term_50);
    result_by_q.set(c_Key1(7), common_term_50);
    // q = 8
    result_by_lpq.set(c_Key3(2, 1, 8), common_term_49);
    result_by_q.set(c_Key1(8), common_term_49);
    // q = 9
    result_by_lpq.set(c_Key3(2, 1, 9), common_term_48);
    result_by_q.set(c_Key1(9), common_term_48);
    // q = 10
    result_by_lpq.set(c_Key3(2, 1, 10), common_term_47);
    result_by_q.set(c_Key1(10), common_term_47);
    // q = 11
    result_by_lpq.set(c_Key3(2, 1, 11), common_term_46);
    result_by_q.set(c_Key1(11), common_term_46);
    // q = 12
    result_by_lpq.set(c_Key3(2, 1, 12), common_term_45);
    result_by_q.set(c_Key1(12), common_term_45);
    // q = 13
    result_by_lpq.set(c_Key3(2, 1, 13), common_term_44);
    result_by_q.set(c_Key1(13), common_term_44);
    // q = 14
    result_by_lpq.set(c_Key3(2, 1, 14), common_term_43);
    result_by_q.set(c_Key1(14), common_term_43);
    // q = 15
    result_by_lpq.set(c_Key3(2, 1, 15), common_term_42);
    result_by_q.set(c_Key1(15), common_term_42);
    // q = 16
    result_by_lpq.set(c_Key3(2, 1, 16), common_term_41);
    result_by_q.set(c_Key1(16), common_term_41);
    // q = 17
    result_by_lpq.set(c_Key3(2, 1, 17), common_term_40);
    result_by_q.set(c_Key1(17), common_term_40);
    // q = 18
    result_by_lpq.set(c_Key3(2, 1, 18), common_term_39);
    result_by_q.set(c_Key1(18), common_term_39);
    // q = 19
    result_by_lpq.set(c_Key3(2, 1, 19), common_term_38);
    result_by_q.set(c_Key1(19), common_term_38);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 1), result_by_q);
    result_by_q.clear();

    // l , p = (2, 2).
    // q = -19
    result_by_lpq.set(c_Key3(2, 2, -19), common_term_37);
    result_by_q.set(c_Key1(-19), common_term_37);
    // q = -18
    result_by_lpq.set(c_Key3(2, 2, -18), common_term_36);
    result_by_q.set(c_Key1(-18), common_term_36);
    // q = -17
    result_by_lpq.set(c_Key3(2, 2, -17), common_term_35);
    result_by_q.set(c_Key1(-17), common_term_35);
    // q = -16
    result_by_lpq.set(c_Key3(2, 2, -16), common_term_34);
    result_by_q.set(c_Key1(-16), common_term_34);
    // q = -15
    result_by_lpq.set(c_Key3(2, 2, -15), common_term_33);
    result_by_q.set(c_Key1(-15), common_term_33);
    // q = -14
    result_by_lpq.set(c_Key3(2, 2, -14), common_term_32);
    result_by_q.set(c_Key1(-14), common_term_32);
    // q = -13
    result_by_lpq.set(c_Key3(2, 2, -13), common_term_31);
    result_by_q.set(c_Key1(-13), common_term_31);
    // q = -12
    result_by_lpq.set(c_Key3(2, 2, -12), common_term_30);
    result_by_q.set(c_Key1(-12), common_term_30);
    // q = -11
    result_by_lpq.set(c_Key3(2, 2, -11), common_term_29);
    result_by_q.set(c_Key1(-11), common_term_29);
    // q = -10
    result_by_lpq.set(c_Key3(2, 2, -10), common_term_28);
    result_by_q.set(c_Key1(-10), common_term_28);
    // q = -9
    result_by_lpq.set(c_Key3(2, 2, -9), common_term_27);
    result_by_q.set(c_Key1(-9), common_term_27);
    // q = -8
    result_by_lpq.set(c_Key3(2, 2, -8), common_term_26);
    result_by_q.set(c_Key1(-8), common_term_26);
    // q = -7
    result_by_lpq.set(c_Key3(2, 2, -7), common_term_25);
    result_by_q.set(c_Key1(-7), common_term_25);
    // q = -6
    result_by_lpq.set(c_Key3(2, 2, -6), common_term_24);
    result_by_q.set(c_Key1(-6), common_term_24);
    // q = -5
    result_by_lpq.set(c_Key3(2, 2, -5), common_term_23);
    result_by_q.set(c_Key1(-5), common_term_23);
    // q = -4
    result_by_lpq.set(c_Key3(2, 2, -4), common_term_22);
    result_by_q.set(c_Key1(-4), common_term_22);
    // q = -3
    result_by_lpq.set(c_Key3(2, 2, -3), common_term_21);
    result_by_q.set(c_Key1(-3), common_term_21);
    // q = -2
    result_by_lpq.set(c_Key3(2, 2, -2), common_term_20);
    result_by_q.set(c_Key1(-2), common_term_20);
    // q = -1
    result_by_lpq.set(c_Key3(2, 2, -1), common_term_19);
    result_by_q.set(c_Key1(-1), common_term_19);
    // q = 0
    result_by_lpq.set(c_Key3(2, 2, 0), common_term_18);
    result_by_q.set(c_Key1(0), common_term_18);
    // q = 1
    result_by_lpq.set(c_Key3(2, 2, 1), common_term_17);
    result_by_q.set(c_Key1(1), common_term_17);
    // q = 3
    result_by_lpq.set(c_Key3(2, 2, 3), common_term_16);
    result_by_q.set(c_Key1(3), common_term_16);
    // q = 4
    result_by_lpq.set(c_Key3(2, 2, 4), common_term_15);
    result_by_q.set(c_Key1(4), common_term_15);
    // q = 5
    result_by_lpq.set(c_Key3(2, 2, 5), common_term_14);
    result_by_q.set(c_Key1(5), common_term_14);
    // q = 6
    result_by_lpq.set(c_Key3(2, 2, 6), common_term_13);
    result_by_q.set(c_Key1(6), common_term_13);
    // q = 7
    result_by_lpq.set(c_Key3(2, 2, 7), common_term_12);
    result_by_q.set(c_Key1(7), common_term_12);
    // q = 8
    result_by_lpq.set(c_Key3(2, 2, 8), common_term_11);
    result_by_q.set(c_Key1(8), common_term_11);
    // q = 9
    result_by_lpq.set(c_Key3(2, 2, 9), common_term_10);
    result_by_q.set(c_Key1(9), common_term_10);
    // q = 10
    result_by_lpq.set(c_Key3(2, 2, 10), common_term_9);
    result_by_q.set(c_Key1(10), common_term_9);
    // q = 11
    result_by_lpq.set(c_Key3(2, 2, 11), common_term_8);
    result_by_q.set(c_Key1(11), common_term_8);
    // q = 12
    result_by_lpq.set(c_Key3(2, 2, 12), common_term_7);
    result_by_q.set(c_Key1(12), common_term_7);
    // q = 13
    result_by_lpq.set(c_Key3(2, 2, 13), common_term_6);
    result_by_q.set(c_Key1(13), common_term_6);
    // q = 14
    result_by_lpq.set(c_Key3(2, 2, 14), common_term_5);
    result_by_q.set(c_Key1(14), common_term_5);
    // q = 15
    result_by_lpq.set(c_Key3(2, 2, 15), common_term_4);
    result_by_q.set(c_Key1(15), common_term_4);
    // q = 16
    result_by_lpq.set(c_Key3(2, 2, 16), common_term_3);
    result_by_q.set(c_Key1(16), common_term_3);
    // q = 17
    result_by_lpq.set(c_Key3(2, 2, 17), common_term_2);
    result_by_q.set(c_Key1(17), common_term_2);
    // q = 18
    result_by_lpq.set(c_Key3(2, 2, 18), common_term_1);
    result_by_q.set(c_Key1(18), common_term_1);
    // q = 19
    result_by_lpq.set(c_Key3(2, 2, 19), common_term_0);
    result_by_q.set(c_Key1(19), common_term_0);
    // Store the q table into the results_lm then reset it.
    
    result_by_lp.set(c_Key2(2, 2), result_by_q);
    result_by_q.clear();

    return EccentricityFuncOutput(result_by_lpq, result_by_lp);
}
