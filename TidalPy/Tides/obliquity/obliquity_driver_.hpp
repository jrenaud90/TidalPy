#pragma once
/*
 * obliquity_driver_.hpp - obliquity functions F_lmp(I) by degree and truncation level.
 *
 * c_obliquity_values returns the unsquared functions through I^N (every function starting at or below I^N), the form
 * a tidal potential uses. c_obliquity_squared_values returns each F_lmp^2 cut at I^N (every function starting at or
 * below I^(N / 2)), the form the global heating uses. Both are dense over (m, p); c_obliquity_func and
 * c_obliquity_squared_func give the same numbers as lookup maps of the non-zero functions. See obliquity_common_.hpp
 * for the truncation rule and the table format. The truncation C_OBLIQUITY_GENERAL takes the exact half-angle form
 * instead, with plain squares.
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

// function_value(m, p) over the functions `include(m, p)` admits; zero for the rest.
template <typename Include, typename FunctionValue>
inline c_ObliquityValues c_fill_obliquity_values(
        int degree_l,
        const Include& include,
        const FunctionValue& function_value) {
    c_ObliquityValues values;
    values.degree_l = degree_l;
    values.values.assign(static_cast<size_t>((degree_l + 1) * (degree_l + 1)), 0.0);
    for (int order_m = 0; order_m <= degree_l; ++order_m) {
        for (int p = 0; p <= degree_l; ++p) {
            if (!include(order_m, p)) { continue; }
            values.values[static_cast<size_t>(order_m * (degree_l + 1) + p)] = function_value(order_m, p);
        }
    }
    return values;
}

inline c_ObliquityValues c_general_obliquity_values(
        int* error_code_ptr, double obliquity, int degree_l, bool squared) {
    const c_ObliquityGeneralTable table = c_obliquity_general_table(error_code_ptr, degree_l);
    if (error_code_ptr[0] != 0) { return c_ObliquityValues(); }
    const c_ObliquityHalfAnglePowers powers(obliquity);
    return c_fill_obliquity_values(
        degree_l,
        [&](int order_m, int p) { return table.mode(order_m, p).count > 0; },
        [&](int order_m, int p) {
            const double value = c_obliquity_general_value(table, order_m, p, powers);
            return squared ? value * value : value;
        });
}

// The non-zero functions of dense values as the (l, m, p) and (l, m) -> (p) maps.
inline c_ModeFuncOutput c_obliquity_maps(const c_ObliquityValues& values) {
    return c_fill_mode_func_maps(
        values.degree_l, values.degree_l + 1, 0, values.degree_l,
        [&](int order_m, int p) { return values.value(order_m, p); });
}

}  // namespace obliquity_detail

// Unsquared F_lmp(I) through I^N for every function starting at or below I^N; for C_OBLIQUITY_GENERAL, the exact
// functions. Error codes as c_obliquity_series_table.
inline c_ObliquityValues c_obliquity_values(
        int* error_code_ptr,
        double obliquity,
        int degree_l,
        int truncation) {
    if (truncation == C_OBLIQUITY_GENERAL) {
        return obliquity_detail::c_general_obliquity_values(error_code_ptr, obliquity, degree_l, false);
    }
    const c_ObliquitySeriesTable table = c_obliquity_series_table(error_code_ptr, degree_l, truncation);
    if (error_code_ptr[0] != 0) { return c_ObliquityValues(); }
    c_ObliquityValues values = obliquity_detail::c_fill_obliquity_values(
        degree_l,
        [&](int order_m, int p) { return table.mode(order_m, p).count > 0; },
        [&](int order_m, int p) { return c_obliquity_mode_value(table, order_m, p, obliquity); });
    values.table = table;
    return values;
}

// F_lmp(I)^2 cut at I^N for every function starting at or below I^(N / 2); the sum over modes is the Taylor series of
// the heating through I^N. For C_OBLIQUITY_GENERAL, the plain squares of the exact functions.
inline c_ObliquityValues c_obliquity_squared_values(
        int* error_code_ptr,
        double obliquity,
        int degree_l,
        int truncation) {
    if (truncation == C_OBLIQUITY_GENERAL) {
        return obliquity_detail::c_general_obliquity_values(error_code_ptr, obliquity, degree_l, true);
    }
    const c_ObliquitySeriesTable table = c_obliquity_series_table(error_code_ptr, degree_l, truncation);
    if (error_code_ptr[0] != 0) { return c_ObliquityValues(); }
    c_ObliquityValues values = obliquity_detail::c_fill_obliquity_values(
        degree_l,
        [&](int order_m, int p) {
            return (table.mode(order_m, p).count > 0)
                && (2 * c_obliquity_lead_power(degree_l, order_m, p) <= truncation);
        },
        [&](int order_m, int p) {
            return c_obliquity_cut_product(table, order_m, p, table, order_m, p, obliquity);
        });
    values.table = table;
    return values;
}

// c_obliquity_values as lookup maps of the non-zero functions.
inline c_ModeFuncOutput c_obliquity_func(
        int* error_code_ptr,
        double obliquity,
        int degree_l,
        int truncation) {
    const c_ObliquityValues values = c_obliquity_values(error_code_ptr, obliquity, degree_l, truncation);
    if (error_code_ptr[0] != 0) { return c_ModeFuncOutput(); }
    return obliquity_detail::c_obliquity_maps(values);
}

// c_obliquity_squared_values as lookup maps of the non-zero functions.
inline c_ModeFuncOutput c_obliquity_squared_func(
        int* error_code_ptr,
        double obliquity,
        int degree_l,
        int truncation) {
    const c_ObliquityValues values = c_obliquity_squared_values(error_code_ptr, obliquity, degree_l, truncation);
    if (error_code_ptr[0] != 0) { return c_ModeFuncOutput(); }
    return obliquity_detail::c_obliquity_maps(values);
}
