#pragma once
/*
 * mode_func_common_.hpp - the lookup maps the eccentricity and obliquity drivers return, filled from one degree's
 * dense function values.
 */

#include <cstddef>
#include <utility>

#include "intmap_.hpp"
#include "keys_.hpp"

// One degree's functions keyed by (l, row, column), and the same values as a (column) map for each (l, row): (l, p, q)
// and (l, p) -> (q) for the eccentricity functions, (l, m, p) and (l, m) -> (p) for the obliquity functions.
typedef std::pair<c_IntMap<c_Key3, double>, c_IntMap<c_Key2, c_IntMap<c_Key1, double>>> c_ModeFuncOutput;

// Fill the maps from value(row, column) over 0 <= row < num_rows and column_min <= column <= column_max, keeping the
// non-zero values. A row with none has no inner map.
template <typename Value>
inline c_ModeFuncOutput c_fill_mode_func_maps(
        int degree_l,
        int num_rows,
        int column_min,
        int column_max,
        const Value& value) {
    const size_t num_columns = static_cast<size_t>(column_max - column_min + 1);
    c_ModeFuncOutput output;
    output.first.reserve(static_cast<size_t>(num_rows) * num_columns);
    output.second.reserve(static_cast<size_t>(num_rows));
    for (int row = 0; row < num_rows; ++row) {
        c_IntMap<c_Key1, double> by_column(num_columns);
        for (int column = column_min; column <= column_max; ++column) {
            const double entry = value(row, column);
            if (entry == 0.0) { continue; }
            output.first.set(c_Key3(degree_l, row, column), entry);
            by_column.set(c_Key1(column), entry);
        }
        if (by_column.size() > 0) {
            // Rows arrive in key order, so each inner map is appended and moved in rather than copied.
            output.second.data.emplace_back(c_Key2(degree_l, row), std::move(by_column));
        }
    }
    return output;
}
