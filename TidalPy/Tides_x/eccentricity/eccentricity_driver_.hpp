#pragma once
/*
 * eccentricity_driver_.hpp - eccentricity functions G_lpq(e) by degree and truncation level.
 *
 * c_eccentricity_func returns the unsquared functions through e^N (every mode with |q| <= N), the form a tidal
 * potential uses. c_eccentricity_squared_func returns each G_lpq^2 cut at e^N (every mode with |q| <= N / 2), the form
 * the global heating uses. See eccentricity_common_.hpp for the truncation rule and the table format. The truncation
 * C_ECCENTRICITY_EXACT takes the functions from the exact orbit instead (eccentricity_exact_.hpp), every mode within
 * the exact tolerance, with plain squares.
 */

#include <cmath>

#include "eccentricity_common_.hpp"
#include "eccentricity_exact_.hpp"
#include "eccentricity_func_l2_.hpp"
#include "eccentricity_func_l3_.hpp"
#include "eccentricity_func_l4_.hpp"
#include "eccentricity_func_l5_.hpp"
#include "eccentricity_func_l6_.hpp"
#include "eccentricity_func_l7_.hpp"
#include "eccentricity_func_l8_.hpp"
#include "eccentricity_func_l9_.hpp"
#include "eccentricity_func_l10_.hpp"

// The table of one (degree, truncation). Sets *error_code_ptr to -1 for an untabulated truncation and -2 for an
// unsupported degree; the returned table is then invalid.
inline c_EccentricitySeriesTable c_eccentricity_series_table(
        int* error_code_ptr,
        int degree_l,
        int truncation) noexcept {
    error_code_ptr[0] = 0;
    c_EccentricitySeriesTable table;
    switch (degree_l) {
        case 2:  table = c_eccentricity_series_l2(truncation);  break;
        case 3:  table = c_eccentricity_series_l3(truncation);  break;
        case 4:  table = c_eccentricity_series_l4(truncation);  break;
        case 5:  table = c_eccentricity_series_l5(truncation);  break;
        case 6:  table = c_eccentricity_series_l6(truncation);  break;
        case 7:  table = c_eccentricity_series_l7(truncation);  break;
        case 8:  table = c_eccentricity_series_l8(truncation);  break;
        case 9:  table = c_eccentricity_series_l9(truncation);  break;
        case 10: table = c_eccentricity_series_l10(truncation); break;
        default:
            error_code_ptr[0] = -2;
            return table;
    }
    if (!table.valid()) {
        error_code_ptr[0] = -1;
    }
    return table;
}

namespace eccentricity_detail {

// Fill the (l, p, q) and (l, p) -> (q) maps from mode_value(p, q) over |q| <= max_q, keeping non-zero modes.
template <typename ModeValue>
inline EccentricityFuncOutput c_fill_eccentricity_maps(
        const c_EccentricitySeriesTable& table,
        int max_q,
        const ModeValue& mode_value) {
    const int degree_l = table.degree_l;
    EccentricityFuncOutput output;
    output.first.reserve(static_cast<size_t>((degree_l + 1) * (2 * max_q + 1)));
    output.second.reserve(static_cast<size_t>(degree_l + 1));
    c_IntMap<c_Key1, double> by_q(static_cast<size_t>(2 * max_q + 1));
    for (int p = 0; p <= degree_l; ++p) {
        by_q.clear();
        for (int q = -max_q; q <= max_q; ++q) {
            if (table.mode(p, q).count == 0) { continue; }
            const double value = mode_value(p, q);
            if (value == 0.0) { continue; }
            output.first.set(c_Key3(degree_l, p, q), value);
            by_q.set(c_Key1(q), value);
        }
        if (by_q.size() > 0) {
            output.second.set(c_Key2(degree_l, p), by_q);
        }
    }
    return output;
}

// The exact functions of one degree as the lookup maps, squared or not. Modes that vanish identically (k = 0 with
// |l - 2p| >= l) are set to zero rather than left at the transform's rounding, so they add no forcing frequency.
inline EccentricityFuncOutput c_fill_exact_eccentricity_maps(
        int* error_code_ptr,
        double eccentricity,
        int degree_l,
        double exact_tolerance,
        bool squared) {
    error_code_ptr[0] = 0;
    if ((degree_l < 2) || (degree_l > 10)) {
        error_code_ptr[0] = -2;
        return EccentricityFuncOutput();
    }
    const c_ExactEccentricityModes modes = c_exact_eccentricity_modes(eccentricity, degree_l, exact_tolerance);
    const int max_q = modes.max_q;
    EccentricityFuncOutput output;
    output.first.reserve(static_cast<size_t>((degree_l + 1) * (2 * max_q + 1)));
    output.second.reserve(static_cast<size_t>(degree_l + 1));
    c_IntMap<c_Key1, double> by_q(static_cast<size_t>(2 * max_q + 1));
    for (int p = 0; p <= degree_l; ++p) {
        by_q.clear();
        const int order_m = degree_l - 2 * p;
        for (int q = -max_q; q <= max_q; ++q) {
            if ((order_m + q == 0) && (std::abs(order_m) >= degree_l)) { continue; }
            const double g = modes.value(p, q);
            const double value = squared ? g * g : g;
            if (value == 0.0) { continue; }
            output.first.set(c_Key3(degree_l, p, q), value);
            by_q.set(c_Key1(q), value);
        }
        if (by_q.size() > 0) {
            output.second.set(c_Key2(degree_l, p), by_q);
        }
    }
    return output;
}

}  // namespace eccentricity_detail

// Unsquared G_lpq(e) through e^N for every non-zero mode with |q| <= N; for C_ECCENTRICITY_EXACT, the exact functions
// of every mode within exact_tolerance (which the tabulated levels ignore). Error codes as c_eccentricity_series_table.
inline EccentricityFuncOutput c_eccentricity_func(
        int* error_code_ptr,
        double eccentricity,
        int degree_l,
        int truncation,
        double exact_tolerance) {
    if (truncation == C_ECCENTRICITY_EXACT) {
        return eccentricity_detail::c_fill_exact_eccentricity_maps(
            error_code_ptr, eccentricity, degree_l, exact_tolerance, false);
    }
    const c_EccentricitySeriesTable table = c_eccentricity_series_table(error_code_ptr, degree_l, truncation);
    if (error_code_ptr[0] != 0) { return EccentricityFuncOutput(); }
    return eccentricity_detail::c_fill_eccentricity_maps(table, truncation, [&](int p, int q) {
        return c_eccentricity_mode_value(table, p, q, eccentricity);
    });
}

// G_lpq(e)^2 cut at e^N for every non-zero mode with |q| <= N / 2 (exact for a k = 0 mode). A cut square can be
// negative for the highest-|q| modes at large e; the sum over modes is the Taylor series of the heating through e^N.
// For C_ECCENTRICITY_EXACT, the plain squares of the exact functions.
inline EccentricityFuncOutput c_eccentricity_squared_func(
        int* error_code_ptr,
        double eccentricity,
        int degree_l,
        int truncation,
        double exact_tolerance) {
    if (truncation == C_ECCENTRICITY_EXACT) {
        return eccentricity_detail::c_fill_exact_eccentricity_maps(
            error_code_ptr, eccentricity, degree_l, exact_tolerance, true);
    }
    const c_EccentricitySeriesTable table = c_eccentricity_series_table(error_code_ptr, degree_l, truncation);
    if (error_code_ptr[0] != 0) { return EccentricityFuncOutput(); }
    return eccentricity_detail::c_fill_eccentricity_maps(table, truncation / 2, [&](int p, int q) {
        return c_eccentricity_cut_product(table, p, q, table, p, q, eccentricity);
    });
}
