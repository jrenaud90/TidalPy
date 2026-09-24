#pragma once
/*
 * eos_data_.hpp: c_LayerEOSData, the per-layer solved-structure store.
 *
 * Density, gravity, pressure, temperature, and heat flow are queried through CyRK dense output via a type-erased
 * callable installed by the world EOS solve. That callable co-owns the solution and is compiled in the world
 * (CyRK-owning) translation unit, so the layer extensions stay CyRK-free; without it the getters fall back to
 * linear interpolation of arrays a caller supplied through update_eos_data, which is how a layer carries a profile
 * that no solve produced. Getters return NaN until populated.
 *
 * The dense evaluator also reports the material state (static moduli, viscosities, melt fraction), which the
 * material evaluated while the structure was integrated, so reading it back calculates nothing and a layer and the
 * solve that used it cannot disagree. All MKS.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>

#include "constants_.hpp"   // TidalPyConstants::d_NAN
#include "interp_.hpp"      // Utilities_x c_interp (numpy.interp-style linear interp)
#include "eos_layout_.hpp"  // C_EOS_DY_VALUES and the evaluation-layout indices

namespace tidalpy {

class c_LayerEOSData {
public:
    // Dense EOS evaluator: fills a buffer of C_EOS_DY_VALUES doubles at a radius [m] from the CyRK dense output.
    // Type-erased so the layer stays CyRK-free.
    using DenseEval = std::function<void(double radius, double* y_out)>;

    c_LayerEOSData()  = default;
    ~c_LayerEOSData() = default;

    bool is_populated() const noexcept {
        return static_cast<bool>(this->p_dense_eval) || !this->p_radius.empty();
    }
    // True once a world EOS solve installed its dense evaluator, which is what carries the material state.
    bool has_dense_eval() const noexcept { return static_cast<bool>(this->p_dense_eval); }

    // The whole evaluation layout at a radius (Material_x/eos/eos_layout_.hpp): CyRK dense output when available,
    // else the linear fallback. Fills C_EOS_DY_VALUES doubles. Without the dense evaluator only the three
    // interpolated structure variables are filled and the rest are NaN.
    void evaluate(double radius, double* y_out) const noexcept {
        if (this->p_dense_eval) {
            this->p_dense_eval(radius, y_out);
            return;
        }
        for (std::size_t value_i = 0; value_i < C_EOS_DY_VALUES; ++value_i) {
            y_out[value_i] = TidalPyConstants::d_NAN;
        }
        if (this->p_radius.empty()) { return; }
        y_out[C_EOS_GRAVITY_INDEX]  = this->interp_fallback(radius, this->p_gravity_ms2);
        y_out[C_EOS_PRESSURE_INDEX] = this->interp_fallback(radius, this->p_pressure);
        y_out[C_EOS_DENSITY_INDEX]  = this->interp_fallback(radius, this->p_density_kgm3);
    }

    // Install the CyRK dense evaluator (set by the world EOS solve). The callable co-owns the solution, so the
    // dense data outlives this store regardless of re-solves.
    void set_dense_eval(DenseEval dense_eval) { this->p_dense_eval = std::move(dense_eval); }

    // Populate the structure slice arrays used by the linear fallback. radius must be sorted ascending.
    void populate(
        const std::vector<double>& radius,
        const std::vector<double>& density_kgm3,
        const std::vector<double>& gravity_ms2,
        const std::vector<double>& pressure)
    {
        const std::size_t num_points = radius.size();
        if ((num_points == 0) || (density_kgm3.size() != num_points) || (gravity_ms2.size() != num_points)
                || (pressure.size() != num_points)) {
            throw std::invalid_argument(
                "TidalPy: a layer EOS profile needs radius, density, gravity, and pressure arrays of one nonzero "
                "length.");
        }
        for (std::size_t point_i = 0; point_i < num_points; ++point_i) {
            if (!std::isfinite(radius[point_i]) || ((point_i > 0) && (radius[point_i] < radius[point_i - 1]))) {
                throw std::invalid_argument("TidalPy: a layer EOS profile's radii must be finite and ascending.");
            }
        }
        this->p_radius     = radius;
        this->p_density_kgm3 = density_kgm3;
        this->p_gravity_ms2 = gravity_ms2;
        this->p_pressure    = pressure;
    }

private:
    DenseEval p_dense_eval;  // CyRK dense output (type-erased, co-owns solution); empty until solved

    std::vector<double> p_radius;      // [m], sorted ascending
    std::vector<double> p_density_kgm3;  // [kg/m^3]   linear fallback
    std::vector<double> p_gravity_ms2;   // [m/s^2]    linear fallback
    std::vector<double> p_pressure;   // [Pa]       linear fallback

    double interp_fallback(double radius, const std::vector<double>& values) const noexcept {
        if (values.size() != this->p_radius.size()) { return TidalPyConstants::d_NAN; }
        return c_interp(radius, this->p_radius.data(), values.data(), this->p_radius.size());
    }
};

} // namespace tidalpy
