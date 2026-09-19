#pragma once
/*
 * eos_data_.hpp: c_LayerEOSData, the per-layer frequency-independent state store.
 *
 * Density, gravity, and pressure are queried through CyRK dense output via a type-erased callable installed by the
 * world EOS solve. That callable co-owns the solution and is compiled in the world (CyRK-owning) translation unit,
 * so the layer extensions stay CyRK-free; without it the structure getters fall back to linear interpolation of the
 * stored slice arrays. The viscoelastic profiles are not part of the ODE solution (they are computed algebraically
 * at the radial slices), so they are always interpolated. Everything here depends only on the solved temperature
 * and pressure state and is cached once per EOS solve; only the complex modulus step is redone per forcing
 * frequency. Getters return NaN until populated. All MKS.
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
    //   5/6 shear modulus re/im, 7/8 bulk modulus re/im, 9 shear visc, 10 bulk visc.
    static constexpr std::size_t EOS_INDEX_GRAVITY  = 0;
    static constexpr std::size_t EOS_INDEX_PRESSURE = 1;
    static constexpr std::size_t EOS_INDEX_DENSITY  = 4;

    c_LayerEOSData()  = default;
    ~c_LayerEOSData() = default;

    bool is_populated() const noexcept {
        return static_cast<bool>(this->p_dense_eval) || !this->p_radius.empty();
    }
    bool is_viscoelastic_populated() const noexcept { return this->p_viscoelastic_populated; }

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

    // Post-melt viscoelastic getters (interpolated: these are not part of the ODE solution).
    double get_shear_modulus(double radius) const noexcept {
        return this->interp_viscoelastic(radius, this->p_postmelt_shear);
    }
    double get_bulk_modulus(double radius) const noexcept {
        return this->interp_viscoelastic(radius, this->p_postmelt_bulk);
    }
    double get_shear_viscosity(double radius) const noexcept {
        return this->interp_viscoelastic(radius, this->p_postmelt_shear_visc);
    }
    double get_bulk_viscosity(double radius) const noexcept {
        return this->interp_viscoelastic(radius, this->p_postmelt_bulk_visc);
    }

    // Pre-melt viscoelastic getters (the un-melted values).
    double get_premelt_shear_modulus(double radius) const noexcept {
        return this->interp_viscoelastic(radius, this->p_premelt_shear);
    }
    double get_premelt_bulk_modulus(double radius) const noexcept {
        return this->interp_viscoelastic(radius, this->p_premelt_bulk);
    }
    double get_premelt_shear_viscosity(double radius) const noexcept {
        return this->interp_viscoelastic(radius, this->p_premelt_shear_visc);
    }
    double get_premelt_bulk_viscosity(double radius) const noexcept {
        return this->interp_viscoelastic(radius, this->p_premelt_bulk_visc);
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

    // Populate the pre and post-melt viscoelastic profiles. All eight vectors must match the radius-grid length.
    void populate_viscoelastic(
        const std::vector<double>& premelt_shear,
        const std::vector<double>& premelt_bulk,
        const std::vector<double>& premelt_shear_visc,
        const std::vector<double>& premelt_bulk_visc,
        const std::vector<double>& postmelt_shear,
        const std::vector<double>& postmelt_bulk,
        const std::vector<double>& postmelt_shear_visc,
        const std::vector<double>& postmelt_bulk_visc)
    {
        this->p_premelt_shear = premelt_shear;
        this->p_premelt_bulk  = premelt_bulk;
        this->p_premelt_shear_visc = premelt_shear_visc;
        this->p_premelt_bulk_visc = premelt_bulk_visc;
        this->p_postmelt_shear    = postmelt_shear;
        this->p_postmelt_bulk     = postmelt_bulk;
        this->p_postmelt_shear_visc = postmelt_shear_visc;
        this->p_postmelt_bulk_visc     = postmelt_bulk_visc;
        this->p_viscoelastic_populated = !this->p_postmelt_shear.empty();
    }

private:
    bool p_viscoelastic_populated = false;

    DenseEval p_dense_eval;  // CyRK dense output (type-erased, co-owns solution); empty until solved

    std::vector<double> p_radius;      // [m], sorted ascending (also the viscoelastic grid)
    std::vector<double> p_density_kgm3;  // [kg/m^3]   linear fallback
    std::vector<double> p_gravity_ms2;   // [m/s^2]    linear fallback
    std::vector<double> p_pressure;   // [Pa]       linear fallback

    // Pre-melt (solid) static moduli + viscosities.
    std::vector<double> p_premelt_shear;        // [Pa]
    std::vector<double> p_premelt_bulk;         // [Pa]
    std::vector<double> p_premelt_shear_visc;  // [Pa·s]
    std::vector<double> p_premelt_bulk_visc;   // [Pa·s]

    // Post-melt (melt-weakened) static moduli + viscosities.
    std::vector<double> p_postmelt_shear;        // [Pa]
    std::vector<double> p_postmelt_bulk;         // [Pa]
    std::vector<double> p_postmelt_shear_visc;  // [Pa·s]
    std::vector<double> p_postmelt_bulk_visc;   // [Pa·s]

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

    // Linear interpolation over the radius grid. NaN when the viscoelastic profiles are not populated.
    double interp_viscoelastic(
            double radius, const std::vector<double>& values) const noexcept {
        if (!this->p_viscoelastic_populated) { return TidalPyConstants::d_NAN; }
        return c_interp(radius, this->p_radius.data(), values.data(), this->p_radius.size());
    }
};

} // namespace tidalpy
