#pragma once
/* c_TableLookup: linear interpolation in a table with an ascending abscissa (radius or pressure), shared by every law
 * that tabulates a property. A bucket index over the abscissa seeds each search, so a read costs the same anywhere
 * in a long table. Values beyond the ends hold the end values; an empty table reads as NaN.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "interp_.hpp"
#include "../../constants_.hpp"

namespace tidalpy {

class c_TableLookup {
public:
    // Index the abscissa; it must already be validated (c_check_table).
    void build(const std::vector<double>& abscissa) {
        this->p_bucket_seed.clear();
        this->p_buckets_per_unit = 0.0;
        this->p_start = abscissa.empty() ? 0.0 : abscissa.front();
        const std::size_t num_points = abscissa.size();
        if (num_points < 3) { return; }
        const double span = abscissa.back() - abscissa.front();
        if (!(span > 0.0)) { return; }
        const std::size_t num_buckets = std::min<std::size_t>(4 * (num_points - 1), std::size_t{1} << 16);
        this->p_bucket_seed.resize(num_buckets);
        this->p_buckets_per_unit = static_cast<double>(num_buckets) / span;
        std::size_t row = 0;
        for (std::size_t bucket_i = 0; bucket_i < num_buckets; ++bucket_i) {
            const double edge = abscissa.front() + static_cast<double>(bucket_i) / this->p_buckets_per_unit;
            while ((row + 2 < num_points) && (abscissa[row + 1] <= edge)) { ++row; }
            this->p_bucket_seed[bucket_i] = static_cast<uint32_t>(row);
        }
    }

    // `values` at `x`, linear between the abscissa points; NaN for an empty table.
    double interpolate(double x, const std::vector<double>& abscissa, const std::vector<double>& values) const noexcept {
        if (values.empty() || (values.size() != abscissa.size())) { return TidalPyConstants::d_NAN; }
        return c_interp(x, abscissa.data(), values.data(), values.size(), this->p_seed(x));
    }

private:
    // The row to start the search from; only the cost of a read depends on it.
    std::size_t p_seed(double x) const noexcept {
        if (this->p_bucket_seed.empty() || !(x > this->p_start)) { return 0; }
        const double position = (x - this->p_start) * this->p_buckets_per_unit;
        const std::size_t last_bucket = this->p_bucket_seed.size() - 1;
        return this->p_bucket_seed[
            (position >= static_cast<double>(last_bucket)) ? last_bucket : static_cast<std::size_t>(position)];
    }

    std::vector<uint32_t> p_bucket_seed;
    double p_buckets_per_unit = 0.0;
    double p_start            = 0.0;
};

// Checks a table: a finite, ascending abscissa with at least one point, and each value table (an empty one means not
// provided) the same length. Throws std::invalid_argument starting with `what`.
inline void c_check_table(
        const std::string& what,
        const std::vector<double>& abscissa,
        const std::vector<const std::vector<double>*>& value_tables) {
    if (abscissa.empty()) {
        throw std::invalid_argument(what + " needs at least one table point.");
    }
    for (std::size_t point_i = 0; point_i < abscissa.size(); ++point_i) {
        if (!std::isfinite(abscissa[point_i]) || ((point_i > 0) && (abscissa[point_i] < abscissa[point_i - 1]))) {
            throw std::invalid_argument(what + " needs a finite, ascending table abscissa.");
        }
    }
    for (const std::vector<double>* table : value_tables) {
        if (!table->empty() && (table->size() != abscissa.size())) {
            throw std::invalid_argument(what + " has a table whose length does not match its abscissa.");
        }
    }
}

} // namespace tidalpy
