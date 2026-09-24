#pragma once
/*
 * obliquity_driver_.hpp - obliquity functions F_lmp(I) by degree and truncation level.
 *
 * c_obliquity_func returns the unsquared functions through I^N (every function starting at or below I^N), the form a
 * tidal potential uses. c_obliquity_squared_func returns each F_lmp^2 cut at I^N (every function starting at or below
 * I^(N / 2)), the form the global heating uses. See obliquity_common_.hpp for the truncation rule and the table
 * format. The truncation C_OBLIQUITY_GENERAL takes the exact half-angle form instead, with plain squares.
 */

#include "obliquity_common_.hpp"
#include "obliquity_func_l2_.hpp"
#include "obliquity_func_l3_.hpp"
#include "obliquity_func_l4_.hpp"
#include "obliquity_func_l5_.hpp"
#include "obliquity_func_l6_.hpp"
#include "obliquity_func_l7_.hpp"
#include "obliquity_func_l8_.hpp"
#include "obliquity_func_l9_.hpp"
#include "obliquity_func_l10_.hpp"

// The series table of one (degree, truncation). Sets *error_code_ptr to -1 for an untabulated truncation and -2 for an
// unsupported degree; the returned table is then invalid.
inline c_ObliquitySeriesTable c_obliquity_series_table(
        int* error_code_ptr,
        int degree_l,
        int truncation) noexcept {
    error_code_ptr[0] = 0;
    c_ObliquitySeriesTable table;
    switch (degree_l) {
        case 2: table = c_obliquity_series_l2(truncation); break;
        case 3: table = c_obliquity_series_l3(truncation); break;
        case 4: table = c_obliquity_series_l4(truncation); break;
        case 5: table = c_obliquity_series_l5(truncation); break;
        case 6: table = c_obliquity_series_l6(truncation); break;
        case 7: table = c_obliquity_series_l7(truncation); break;
        case 8: table = c_obliquity_series_l8(truncation); break;
        case 9: table = c_obliquity_series_l9(truncation); break;
        case 10: table = c_obliquity_series_l10(truncation); break;
        default:
            error_code_ptr[0] = -2;
            return table;
    }
    if (!table.valid()) {
        error_code_ptr[0] = -1;
    }
    return table;
}

// The general (half-angle) table of one degree. Sets *error_code_ptr to -2 for an unsupported degree.
inline c_ObliquityGeneralTable c_obliquity_general_table(int* error_code_ptr, int degree_l) noexcept {
    error_code_ptr[0] = 0;
    switch (degree_l) {
        case 2: return c_obliquity_general_l2();
        case 3: return c_obliquity_general_l3();
        case 4: return c_obliquity_general_l4();
        case 5: return c_obliquity_general_l5();
        case 6: return c_obliquity_general_l6();
        case 7: return c_obliquity_general_l7();
        case 8: return c_obliquity_general_l8();
        case 9: return c_obliquity_general_l9();
        case 10: return c_obliquity_general_l10();
        default:
            error_code_ptr[0] = -2;
            return c_ObliquityGeneralTable();
    }
}

namespace obliquity_detail {

// Fill the (l, m, p) and (l, m) -> (p) maps from function_value(m, p) over the functions `include(m, p)` admits,
// keeping non-zero values.
template <typename Include, typename FunctionValue>
inline ObliquityFuncOutput c_fill_obliquity_maps(
        int degree_l,
        const Include& include,
        const FunctionValue& function_value) {
    ObliquityFuncOutput output;
    output.first.reserve(static_cast<size_t>((degree_l + 1) * (degree_l + 1)));
    output.second.reserve(static_cast<size_t>(degree_l + 1));
    c_IntMap<c_Key1, double> by_p(static_cast<size_t>(degree_l + 1));
    for (int order_m = 0; order_m <= degree_l; ++order_m) {
        by_p.clear();
        for (int p = 0; p <= degree_l; ++p) {
            if (!include(order_m, p)) { continue; }
            const double value = function_value(order_m, p);
            if (value == 0.0) { continue; }
            output.first.set(c_Key3(degree_l, order_m, p), value);
            by_p.set(c_Key1(p), value);
        }
        if (by_p.size() > 0) {
            output.second.set(c_Key2(degree_l, order_m), by_p);
        }
    }
    return output;
}

inline ObliquityFuncOutput c_general_obliquity_maps(
        int* error_code_ptr, double obliquity, int degree_l, bool squared) {
    const c_ObliquityGeneralTable table = c_obliquity_general_table(error_code_ptr, degree_l);
    if (error_code_ptr[0] != 0) { return ObliquityFuncOutput(); }
    const c_ObliquityHalfAnglePowers powers(obliquity);
    return c_fill_obliquity_maps(
        degree_l,
        [&](int order_m, int p) { return table.mode(order_m, p).count > 0; },
        [&](int order_m, int p) {
            const double value = c_obliquity_general_value(table, order_m, p, powers);
            return squared ? value * value : value;
        });
}

}  // namespace obliquity_detail

// Unsquared F_lmp(I) through I^N for every non-zero function starting at or below I^N; for C_OBLIQUITY_GENERAL, the
// exact functions. Error codes as c_obliquity_series_table.
inline ObliquityFuncOutput c_obliquity_func(
        int* error_code_ptr,
        double obliquity,
        int degree_l,
        int truncation) {
    if (truncation == C_OBLIQUITY_GENERAL) {
        return obliquity_detail::c_general_obliquity_maps(error_code_ptr, obliquity, degree_l, false);
    }
    const c_ObliquitySeriesTable table = c_obliquity_series_table(error_code_ptr, degree_l, truncation);
    if (error_code_ptr[0] != 0) { return ObliquityFuncOutput(); }
    return obliquity_detail::c_fill_obliquity_maps(
        degree_l,
        [&](int order_m, int p) { return table.mode(order_m, p).count > 0; },
        [&](int order_m, int p) { return c_obliquity_mode_value(table, order_m, p, obliquity); });
}

// F_lmp(I)^2 cut at I^N for every non-zero function starting at or below I^(N / 2); the sum over modes is the Taylor
// series of the heating through I^N. For C_OBLIQUITY_GENERAL, the plain squares of the exact functions.
inline ObliquityFuncOutput c_obliquity_squared_func(
        int* error_code_ptr,
        double obliquity,
        int degree_l,
        int truncation) {
    if (truncation == C_OBLIQUITY_GENERAL) {
        return obliquity_detail::c_general_obliquity_maps(error_code_ptr, obliquity, degree_l, true);
    }
    const c_ObliquitySeriesTable table = c_obliquity_series_table(error_code_ptr, degree_l, truncation);
    if (error_code_ptr[0] != 0) { return ObliquityFuncOutput(); }
    return obliquity_detail::c_fill_obliquity_maps(
        degree_l,
        [&](int order_m, int p) {
            return (table.mode(order_m, p).count > 0)
                && (2 * c_obliquity_lead_power(degree_l, order_m, p) <= truncation);
        },
        [&](int order_m, int p) {
            return c_obliquity_cut_product(table, order_m, p, table, order_m, p, obliquity);
        });
}
