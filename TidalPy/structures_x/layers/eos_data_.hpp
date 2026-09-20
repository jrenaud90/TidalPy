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
 * The viscoelastic state is not stored here. Moduli and viscosities are evaluated on demand from this structure
 * through c_PhysicsLayer::calc_material_state, so a layer and the solve that used it cannot disagree about them
 * and no value is read off a slice grid. All MKS.
 */

#include <cstddef>
#include <functional>
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

    // CyRK EOS-ODE y-layout (see Material_x/eos/ode_.hpp):
    //   0 gravity, 1 pressure, 2 mass, 3 moment-of-inertia, 4 density,
    //   5/6 shear modulus re/im, 7/8 bulk modulus re/im, 9 shear visc, 10 bulk visc,
    //   11 temperature, 12 heat flow.
    static constexpr std::size_t EOS_INDEX_GRAVITY  = 0;
    static constexpr std::size_t EOS_INDEX_PRESSURE = 1;
    static constexpr std::size_t EOS_INDEX_DENSITY  = 4;

    c_LayerEOSData()  = default;
    ~c_LayerEOSData() = default;

    bool is_populated() const noexcept {
        return static_cast<bool>(this->p_dense_eval) || !this->p_radius.empty();
    }

    // Structure getters: CyRK dense output when available, else the linear fallback.
    double get_density(double radius) const noexcept {
        return this->dense_or_interp(radius, EOS_INDEX_DENSITY, this->p_density_kgm3);
    }
    double get_gravity(double radius) const noexcept {
        return this->dense_or_interp(radius, EOS_INDEX_GRAVITY, this->p_gravity_ms2);
    }
    double get_pressure(double radius) const noexcept {
        return this->dense_or_interp(radius, EOS_INDEX_PRESSURE, this->p_pressure);
    }

    // The whole evaluation layout at a radius, for a caller that needs more than one value and does not want to
    // pay for a dense call each. Fills C_EOS_DY_VALUES doubles. Without the dense evaluator only the three
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
        y_out[EOS_INDEX_GRAVITY]  = this->interp_fallback(radius, this->p_gravity_ms2);
        y_out[EOS_INDEX_PRESSURE] = this->interp_fallback(radius, this->p_pressure);
        y_out[EOS_INDEX_DENSITY]  = this->interp_fallback(radius, this->p_density_kgm3);
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

    // Prefers the CyRK dense output; otherwise linear-interpolates the slice array. NaN when neither is available.
    double dense_or_interp(
            double radius, std::size_t dense_index,
            const std::vector<double>& fallback_values) const noexcept {
        if (this->p_dense_eval) {
            double dense_output[C_EOS_DY_VALUES] = {0.0};
            this->p_dense_eval(radius, dense_output);
            return dense_output[dense_index];
        }
        if (this->p_radius.empty()) { return TidalPyConstants::d_NAN; }
        return c_interp(radius, this->p_radius.data(), fallback_values.data(), this->p_radius.size());
    }

    double interp_fallback(double radius, const std::vector<double>& values) const noexcept {
        if (values.size() != this->p_radius.size()) { return TidalPyConstants::d_NAN; }
        return c_interp(radius, this->p_radius.data(), values.data(), this->p_radius.size());
    }
};

} // namespace tidalpy
