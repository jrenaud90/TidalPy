#pragma once
/*
 * eos_data_.hpp — c_LayerEOSData: per-layer EOS interpolation data.
 *
 * Stores sorted radius-indexed arrays for density, gravity, and pressure.
 * Populated by c_EOSHandler (Material_x/, Phase 8). Until populated,
 * all getters return NaN.
 *
 * Interpolation: linear, clamped to boundary values.
 *
 * Interface notes:
 *   - populate() is called by EOSHandler after solving the Adams-Williamson ODE.
 *   - For unit testing (Phase 1), populate() may be called directly.
 *   - All inputs and outputs are MKS.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace tidalpy {

class c_LayerEOSData {
public:
    c_LayerEOSData()  = default;
    ~c_LayerEOSData() = default;

    // -----------------------------------------------------------------------
    // State query
    // -----------------------------------------------------------------------
    bool is_populated() const noexcept { return p_populated; }

    // -----------------------------------------------------------------------
    // Getters (all const, MKS; return NaN when not populated)
    // -----------------------------------------------------------------------
    // Density at radius_m [kg/m^3].
    double get_density(double radius_m) const noexcept {
        if (!p_populated) return std::numeric_limits<double>::quiet_NaN();
        return interp_linear(p_radius_m, p_density_kgm3, radius_m);
    }

    // Gravitational acceleration at radius_m [m/s^2].
    double get_gravity(double radius_m) const noexcept {
        if (!p_populated) return std::numeric_limits<double>::quiet_NaN();
        return interp_linear(p_radius_m, p_gravity_ms2, radius_m);
    }

    // Pressure at radius_m [Pa].
    double get_pressure(double radius_m) const noexcept {
        if (!p_populated) return std::numeric_limits<double>::quiet_NaN();
        return interp_linear(p_radius_m, p_pressure_pa, radius_m);
    }

    // -----------------------------------------------------------------------
    // Populate (called by EOSHandler in Phase 8, or directly for testing)
    //
    // All four vectors must be the same length. radius_m must be sorted
    // ascending. Passing empty vectors resets to the unpopulated state.
    // -----------------------------------------------------------------------
    void populate(
        const std::vector<double>& radius_m,
        const std::vector<double>& density_kgm3,
        const std::vector<double>& gravity_ms2,
        const std::vector<double>& pressure_pa)
    {
        p_radius_m     = radius_m;
        p_density_kgm3 = density_kgm3;
        p_gravity_ms2  = gravity_ms2;
        p_pressure_pa  = pressure_pa;
        p_populated    = !p_radius_m.empty();
    }

private:
    bool p_populated = false;

    std::vector<double> p_radius_m;      // [m],      sorted ascending
    std::vector<double> p_density_kgm3;  // [kg/m^3]
    std::vector<double> p_gravity_ms2;   // [m/s^2]
    std::vector<double> p_pressure_pa;   // [Pa]

    // Linear interpolation clamped to boundary values.
    static double interp_linear(
        const std::vector<double>& x_vec,
        const std::vector<double>& y_vec,
        double x) noexcept
    {
        const std::size_t n = x_vec.size();
        if (n == 0) return std::numeric_limits<double>::quiet_NaN();
        if (x <= x_vec.front()) return y_vec.front();
        if (x >= x_vec.back())  return y_vec.back();

        // Binary search for the bounding interval [lo, hi].
        std::size_t lo = 0;
        std::size_t hi = n - 1;
        while (hi - lo > 1) {
            const std::size_t mid = (lo + hi) / 2;
            if (x_vec[mid] <= x) lo = mid;
            else hi = mid;
        }
        const double t = (x - x_vec[lo]) / (x_vec[hi] - x_vec[lo]);
        return y_vec[lo] + t * (y_vec[hi] - y_vec[lo]);
    }
};

} // namespace tidalpy
