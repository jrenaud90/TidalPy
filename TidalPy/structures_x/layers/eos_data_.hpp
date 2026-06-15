#pragma once
/*
 * eos_data_.hpp — c_LayerEOSData: per-layer state store (frequency-independent).
 *
 * Structure quantities (density, gravity, pressure) are queried through CyRK's
 * DENSE OUTPUT — the high-accuracy continuous interpolant produced by the ODE
 * solver, via a type-erased callable (DenseEval) installed by the world EOS
 * solve. Each Cython extension compiles its own copy of CyRK, so a CyRK CySolverResult
 * created in the world extension must only ever be touched by that same copy.
 * The lambda that performs the dense call is compiled in the world translation unit
 * (which owns CyRK and the solution) and co-owns the solution via a captured shared_ptr,
 * so there is no dangling pointer / leak and no cross-extension CyRK call.
 * The layer extensions stay CyRK-free. When no dense evaluator is installed (e.g. a manual
 * update_eos_data for testing) the structure getters fall back to linear
 * interpolation of the stored slice arrays via the shared Utilities_x `c_interp`.
 *
 * Viscoelastic quantities (pre/post-melt shear & bulk modulus + viscosity) are
 * NOT part of the ODE solution — they are computed algebraically at discrete radial
 * slices in the solve's post-pass — so they are stored as arrays and queried with
 * linear interpolation (`c_interp`).
 *
 * All of these are FREQUENCY-INDEPENDENT (they depend only on the solved
 * temperature/pressure state), so they are computed once per EOS solve and cached
 * here; only the downstream rheology (complex modulus) step is recomputed per
 * forcing frequency. Until populated, the getters return NaN. All MKS.
 */

#include <cstddef>
#include <functional>
#include <vector>

#include "constants_.hpp"   // TidalPyConstants::d_NAN
#include "interp_.hpp"      // Utilities_x c_interp (numpy.interp-style linear interp)

namespace tidalpy {

class c_LayerEOSData {
public:
    // Dense EOS evaluator: fills a buffer of EOS_DENSE_SIZE doubles at a radius [m]
    // from the CyRK dense output. Type-erased so the layer stays CyRK-free and the
    // dense call is compiled in the world (CyRK-owning) extension only.
    using DenseEval = std::function<void(double radius_m, double* y_out)>;

    // CyRK EOS-ODE y-layout (see Material_x/eos/ode_.hpp):
    //   primary y: 0 gravity, 1 pressure, 2 mass, 3 moment-of-inertia.
    //   extra:     4 density, 5-8 complex shear/bulk modulus.
    // IMPORTANT: CyRK's dense interpolation is only safe post-solve for the PRIMARY
    // y-values (interpolated from stored polynomial data). The EXTRA outputs
    // (index >= EOS_NUM_PRIMARY_Y) are recomputed at call time through a stored
    // pointer to the (now-destroyed) solver, so dense interpolation of them after
    // the solve is undefined behaviour. Density and moduli therefore use the
    // stored slice arrays (linear interp); gravity/pressure use the dense output.
    static constexpr std::size_t EOS_DENSE_SIZE     = 9;
    static constexpr std::size_t EOS_NUM_PRIMARY_Y  = 4;
    static constexpr std::size_t EOS_INDEX_GRAVITY  = 0;
    static constexpr std::size_t EOS_INDEX_PRESSURE = 1;
    static constexpr std::size_t EOS_INDEX_DENSITY  = 4;

    c_LayerEOSData()  = default;
    ~c_LayerEOSData() = default;

    // -----------------------------------------------------------------------
    // State query
    // -----------------------------------------------------------------------
    bool is_populated() const noexcept {
        return static_cast<bool>(this->p_dense_eval) || !this->p_radius_m.empty();
    }
    bool is_viscoelastic_populated() const noexcept { return this->p_viscoelastic_populated; }

    // -----------------------------------------------------------------------
    // Structure getters (CyRK dense output when available, else linear fallback).
    // -----------------------------------------------------------------------
    double get_density(double radius_m) const noexcept {
        return this->dense_or_interp(radius_m, EOS_INDEX_DENSITY, this->p_density_kgm3);
    }
    double get_gravity(double radius_m) const noexcept {
        return this->dense_or_interp(radius_m, EOS_INDEX_GRAVITY, this->p_gravity_ms2);
    }
    double get_pressure(double radius_m) const noexcept {
        return this->dense_or_interp(radius_m, EOS_INDEX_PRESSURE, this->p_pressure_pa);
    }

    // -----------------------------------------------------------------------
    // Post-melt viscoelastic getters (linear interp of the slice arrays).
    // -----------------------------------------------------------------------
    double get_shear_modulus(double radius_m) const noexcept {
        return this->interp_viscoelastic(radius_m, this->p_postmelt_shear_pa);
    }
    double get_bulk_modulus(double radius_m) const noexcept {
        return this->interp_viscoelastic(radius_m, this->p_postmelt_bulk_pa);
    }
    double get_shear_viscosity(double radius_m) const noexcept {
        return this->interp_viscoelastic(radius_m, this->p_postmelt_shear_visc_pas);
    }
    double get_bulk_viscosity(double radius_m) const noexcept {
        return this->interp_viscoelastic(radius_m, this->p_postmelt_bulk_visc_pas);
    }

    // -----------------------------------------------------------------------
    // Pre-melt viscoelastic getters (the un-melted values).
    // -----------------------------------------------------------------------
    double get_premelt_shear_modulus(double radius_m) const noexcept {
        return this->interp_viscoelastic(radius_m, this->p_premelt_shear_pa);
    }
    double get_premelt_bulk_modulus(double radius_m) const noexcept {
        return this->interp_viscoelastic(radius_m, this->p_premelt_bulk_pa);
    }
    double get_premelt_shear_viscosity(double radius_m) const noexcept {
        return this->interp_viscoelastic(radius_m, this->p_premelt_shear_visc_pas);
    }
    double get_premelt_bulk_viscosity(double radius_m) const noexcept {
        return this->interp_viscoelastic(radius_m, this->p_premelt_bulk_visc_pas);
    }

    // -----------------------------------------------------------------------
    // Install the CyRK dense evaluator (set by the world EOS solve). The callable
    // co-owns the solution (captured shared_ptr), so the dense data outlives this
    // store regardless of re-solves.
    // -----------------------------------------------------------------------
    void set_dense_eval(DenseEval dense_eval) { this->p_dense_eval = std::move(dense_eval); }

    // -----------------------------------------------------------------------
    // Populate the structure slice arrays (the radius grid + the linear-fallback
    // density/gravity/pressure). radius_m must be sorted ascending.
    // -----------------------------------------------------------------------
    void populate(
        const std::vector<double>& radius_m,
        const std::vector<double>& density_kgm3,
        const std::vector<double>& gravity_ms2,
        const std::vector<double>& pressure_pa)
    {
        this->p_radius_m     = radius_m;
        this->p_density_kgm3 = density_kgm3;
        this->p_gravity_ms2  = gravity_ms2;
        this->p_pressure_pa  = pressure_pa;
    }

    // -----------------------------------------------------------------------
    // Populate the pre/post-melt viscoelastic profiles (same radius grid as the
    // structure). All eight vectors must match the radius-grid length.
    // -----------------------------------------------------------------------
    void populate_viscoelastic(
        const std::vector<double>& premelt_shear_pa,
        const std::vector<double>& premelt_bulk_pa,
        const std::vector<double>& premelt_shear_visc_pas,
        const std::vector<double>& premelt_bulk_visc_pas,
        const std::vector<double>& postmelt_shear_pa,
        const std::vector<double>& postmelt_bulk_pa,
        const std::vector<double>& postmelt_shear_visc_pas,
        const std::vector<double>& postmelt_bulk_visc_pas)
    {
        this->p_premelt_shear_pa        = premelt_shear_pa;
        this->p_premelt_bulk_pa         = premelt_bulk_pa;
        this->p_premelt_shear_visc_pas  = premelt_shear_visc_pas;
        this->p_premelt_bulk_visc_pas   = premelt_bulk_visc_pas;
        this->p_postmelt_shear_pa       = postmelt_shear_pa;
        this->p_postmelt_bulk_pa        = postmelt_bulk_pa;
        this->p_postmelt_shear_visc_pas = postmelt_shear_visc_pas;
        this->p_postmelt_bulk_visc_pas  = postmelt_bulk_visc_pas;
        this->p_viscoelastic_populated  = !this->p_postmelt_shear_pa.empty();
    }

private:
    bool p_viscoelastic_populated = false;

    DenseEval p_dense_eval;  // CyRK dense output (type-erased, co-owns solution); empty until solved

    std::vector<double> p_radius_m;      // [m], sorted ascending (also the viscoelastic grid)
    std::vector<double> p_density_kgm3;  // [kg/m^3]   linear fallback
    std::vector<double> p_gravity_ms2;   // [m/s^2]    linear fallback
    std::vector<double> p_pressure_pa;   // [Pa]       linear fallback

    // Pre-melt (solid) static moduli + viscosities.
    std::vector<double> p_premelt_shear_pa;        // [Pa]
    std::vector<double> p_premelt_bulk_pa;         // [Pa]
    std::vector<double> p_premelt_shear_visc_pas;  // [Pa·s]
    std::vector<double> p_premelt_bulk_visc_pas;   // [Pa·s]

    // Post-melt (melt-weakened) static moduli + viscosities.
    std::vector<double> p_postmelt_shear_pa;        // [Pa]
    std::vector<double> p_postmelt_bulk_pa;         // [Pa]
    std::vector<double> p_postmelt_shear_visc_pas;  // [Pa·s]
    std::vector<double> p_postmelt_bulk_visc_pas;   // [Pa·s]

    // Structure query: for a primary y-value (gravity/pressure) use the CyRK dense
    // output (high accuracy, safe post-solve). For an extra output (density) the
    // dense interpolant would dereference the destroyed solver, so always use the
    // stored slice array. Falls back to the slice array when no dense evaluator is
    // installed. NaN when neither is available.
    double dense_or_interp(
            double radius_m, std::size_t dense_index,
            const std::vector<double>& fallback_values) const noexcept {
        if (this->p_dense_eval && dense_index < EOS_NUM_PRIMARY_Y) {
            double dense_output[EOS_DENSE_SIZE] = {0.0};
            this->p_dense_eval(radius_m, dense_output);
            return dense_output[dense_index];
        }
        if (this->p_radius_m.empty()) { return TidalPyConstants::d_NAN; }
        return c_interp(radius_m, this->p_radius_m.data(), fallback_values.data(), this->p_radius_m.size());
    }

    // Viscoelastic query: linear interpolation over the radius grid (NaN when the
    // viscoelastic profiles have not been populated).
    double interp_viscoelastic(
            double radius_m, const std::vector<double>& values) const noexcept {
        if (!this->p_viscoelastic_populated) { return TidalPyConstants::d_NAN; }
        return c_interp(radius_m, this->p_radius_m.data(), values.data(), this->p_radius_m.size());
    }
};

} // namespace tidalpy
