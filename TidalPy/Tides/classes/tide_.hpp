#pragma once
/*
 * tide_.hpp - TidalPy global (1D) tidal dissipation models: c_RheologyTide (alias "rheology"),
 * c_FixedQTide ("cpl"/"fixed_q"), c_FixedLagTide ("ctl"/"fixed_dt"), and c_CTLQTide
 * ("ctl_q"/"fixed_dt_q").
 *
 * Each returns the complex Love number k_l at a tidal frequency; the collapse (tide_collapse_.hpp)
 * uses -Im[k_l] as the per-mode dissipation multiplier. Each model declares its parameters in one table
 * (c_SpecModel, spec_model_.hpp), which gives its constructor, validation, config entries, and binary record. The
 * analytic models' parameters (k_l, Q_l, dt_l) are lists indexed from l = 2, up to l = 10, the range the
 * eccentricity and obliquity tables cover; a degree a list does not reach is 0, no contribution at that degree. All
 * quantities MKS; frequencies in rad s-1.
 */

#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "constants_.hpp"     // TidalPyConstants::d_EPS
#include "tide_base_.hpp"
#include "tide_result_.hpp"   // c_TideSolveConfig (orbital state), c_TideConfig (truncation)
#include "../../Utilities/classes/registry_.hpp"
#include "../../Utilities/classes/spec_model_.hpp"

namespace tidalpy {

// The 3D methods of c_RheologyTide are declared here but defined in Structures/worlds/world_tides_.hpp,
// compiled into the world extension alone, so this header never pulls in the world, potential-engine, or
// kernel headers. The incomplete type is legal in a declaration and complete at each definition.
class c_BaseWorld;

// Grid outputs are row-major in the order radius, colatitude, longitude, time.
struct c_Grid3DAxes {
    const double* radii           = nullptr;   // [m]
    size_t        num_radii       = 0;
    const double* colatitudes     = nullptr;   // [rad]
    size_t        num_colatitudes = 0;
    const double* longitudes      = nullptr;   // [rad]
    size_t        num_longitudes  = 0;
    const double* times           = nullptr;   // [s]
    size_t        num_times       = 0;
};

// l = 2..10, matching the eccentricity and obliquity tables.
constexpr int C_TIDE_MIN_DEGREE  = 2;
constexpr int C_TIDE_MAX_DEGREE  = 10;
constexpr int C_TIDE_NUM_DEGREES = C_TIDE_MAX_DEGREE - C_TIDE_MIN_DEGREE + 1;  // 9

// A per-degree list's value at a degree; 0 for a degree it does not reach.
inline double c_tide_degree_value(const std::vector<double>& values, int degree_l) noexcept {
    const int index = degree_l - C_TIDE_MIN_DEGREE;
    if ((index < 0) || (index >= static_cast<int>(values.size()))) {
        return 0.0;
    }
    return values[static_cast<std::size_t>(index)];
}


// k_l supplied by the radial solver (alias "rheology"); a pass-through that flags the world to run it.
class c_RheologyTide final : public c_SpecModel<c_RheologyTide, c_TideBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::RheologyTide;

    static const std::vector<c_ParamSpec<c_RheologyTide>>& parameter_specs() {
        static const std::vector<c_ParamSpec<c_RheologyTide>> specs = {};
        return specs;
    }

    c_RheologyTide() : c_RheologyTide(c_ParamMap{}) {}
    explicit c_RheologyTide(const c_ParamMap& params) : c_SpecModel("rheology") { this->p_initialize(params); }

    c_LoveNumbers calc_love_numbers(
            int /*degree_l*/, double /*frequency*/, const c_LoveNumbers& solver_love) const override {
        return solver_love;
    }

    bool needs_radial_solve() const override { return true; }

    // Secular 3D tidal volumetric heating [W m-3], averaged over longitude; the rheology model alone
    // supports the 3D path. The active modes merge into coherent waves, the radial problem is solved once
    // per (l, |omega|), and each frequency contributes (|omega|/2) Im(sigma_c : conj(eps_c)) of its summed
    // complex amplitudes, so the volume integral equals the world's 1D get_tidal_heating. NaN at the center
    // and below the solver's starting radius, 0 in liquid layers, which have no shear kernel.
    double calc_3d_tidal_heating(
            c_BaseWorld& world,
            const c_TideSolveConfig& state,
            double radius,
            double colatitude) const;

    // Batch form over paired (radii[i], colatitudes[i]) points. The wave list is built once and the radial
    // solve amortized across the points, since it depends on (l, |omega|) alone; points sharing a
    // colatitude share its angular work, and the colatitudes run on up to num_threads threads.
    void calc_3d_tidal_heating_batch(
            c_BaseWorld& world,
            const c_TideSolveConfig& state,
            const double* radii,
            const double* colatitudes,
            size_t num_points,
            double* out_heating,
            int num_threads) const;

    // Instantaneous tidal displacements [m] (radial, polar, azimuthal) on the axes' grid. Each coherent
    // wave's complex amplitude at (r, theta, phi) is (y1 U_c, y3 dU_c/dtheta, y3 dU_c/dphi / sin theta),
    // and each component at time t sums Re[amplitude e^{i |omega| t}] over the frequencies. out_disp holds
    // 3 * nr * nth * nph * nt doubles ordered (r, theta, phi, t, component), NaN where the radius has no
    // depth-resolved solution.
    void calc_3d_displacements_grid(
            c_BaseWorld& world,
            const c_TideSolveConfig& state,
            const c_Grid3DAxes& axes,
            double* out_disp,
            int num_threads) const;

    // Instantaneous stress [Pa] and strain on the axes' grid, as 6 * nr * nth * nph * nt doubles ordered
    // (r, theta, phi, t, component) with the components rr, tt, pp, rt, rp, tp. Either output may be null
    // to skip it. NaN where no wave has a shear kernel: a radius with no depth-resolved solution, or a
    // liquid layer. Modes at zero forcing frequency, the permanent tide, are excluded.
    void calc_3d_stress_strain_grid(
            c_BaseWorld& world,
            const c_TideSolveConfig& state,
            const c_Grid3DAxes& axes,
            double* out_stress,
            double* out_strain,
            int num_threads) const;

    // Collapsed secular 3D tidal heating: the density is reduced along whichever of the colatitude,
    // longitude, and radial dimensions cfg names. radii and colatitudes are the user grids for the
    // surviving axes; a summed axis uses an internal integration grid. Writes the marginal power density
    // into out_values and, when all three spatial axes are summed, the per-layer totals into
    // out_layer_totals. The radial solves run on the calling thread, the per-point evaluation on up to
    // cfg.num_threads threads.
    void calc_3d_tidal_heating_collapsed(
            c_BaseWorld& world,
            const c_TideSolveConfig& state,
            const double* radii,
            size_t num_radii,
            const double* colatitudes,
            size_t num_colatitudes,
            const double* longitudes,
            size_t num_longitudes,
            const double* times,
            size_t num_times,
            const c_Heating3DCollapseConfig& cfg,
            double* out_values,
            double* out_layer_totals) const;
};

// The shared part of the analytic tide models (Derived is the concrete model): no radial solve, h and l undefined,
// and every per-degree list no longer than the tabulated degrees.
template <class Derived>
class c_AnalyticTide : public c_SpecModel<Derived, c_TideBase> {
public:
    explicit c_AnalyticTide(const std::string& model_name) : c_SpecModel<Derived, c_TideBase>(model_name) {}

    bool needs_radial_solve() const override { return false; }

protected:
    void p_validate() const override {
        for (const c_ParamSpec<Derived>& spec : Derived::parameter_specs()) {
            const std::size_t num_values = this->get_parameter(spec.key).size();
            if (num_values > static_cast<std::size_t>(C_TIDE_NUM_DEGREES)) {
                throw std::invalid_argument(
                    this->p_describe() + " takes at most " + std::to_string(C_TIDE_NUM_DEGREES) + " values for '"
                    + spec.key + "' (degrees l = 2 to 10); got " + std::to_string(num_values) + ".");
            }
        }
    }

    // h and l are not defined for an analytic model (no radial solution).
    static c_LoveNumbers p_analytic_love(const std::complex<double>& k) {
        const std::complex<double> nan_love(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
        return c_LoveNumbers(k, nan_love, nan_love);
    }
};

// c_FixedQTide: constant phase lag / fixed Q (alias "cpl" / "fixed_q"). Frequency-independent
// dissipation per degree:
//
//   k_l(omega) = k_l * (1 - i / Q_l)            ->  -Im[k_l] = k_l / Q_l
class c_FixedQTide final : public c_AnalyticTide<c_FixedQTide> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::FixedQTide;

    static const std::vector<c_ParamSpec<c_FixedQTide>>& parameter_specs() {
        using Self = c_FixedQTide;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"fixed_k", "fixed_k", &Self::p_fixed_k, 0.0, c_ParamBounds::NonNegative,
             "Static potential Love number k_l of each degree, from l = 2."},
            {"fixed_q", "fixed_q", &Self::p_fixed_q, 0.0, c_ParamBounds::NonNegative,
             "Tidal quality factor Q_l of each degree, from l = 2; 0 is no dissipation."},
        };
        return specs;
    }

    c_FixedQTide() : c_FixedQTide(c_ParamMap{}) {}
    explicit c_FixedQTide(const c_ParamMap& params) : c_AnalyticTide("fixed_q") { this->p_initialize(params); }

    double get_fixed_k(int degree_l) const override { return c_tide_degree_value(this->p_fixed_k, degree_l); }
    double get_fixed_q(int degree_l) const override { return c_tide_degree_value(this->p_fixed_q, degree_l); }

    c_LoveNumbers calc_love_numbers(
            int degree_l, double /*frequency*/, const c_LoveNumbers& /*solver_love*/) const override {
        const double k_l = this->get_fixed_k(degree_l);
        const double q_l = this->get_fixed_q(degree_l);
        if (std::abs(q_l) <= TidalPyConstants::d_EPS) {
            // An unset or zero Q_l is no dissipation, not a divide by zero.
            return p_analytic_love(std::complex<double>(k_l, 0.0));
        }
        return p_analytic_love(k_l * c_constant_phase_lag(q_l));
    }

protected:
    std::vector<double> p_fixed_k;
    std::vector<double> p_fixed_q;
};

// c_FixedLagTide: constant time lag / CTL (alias "ctl" / "fixed_dt"):
//
//   k_l(omega) = k_l * (1 - i * |omega| * dt_l)   ->  -Im[k_l] = k_l * |omega| * dt_l
class c_FixedLagTide final : public c_AnalyticTide<c_FixedLagTide> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::FixedLagTide;

    static const std::vector<c_ParamSpec<c_FixedLagTide>>& parameter_specs() {
        using Self = c_FixedLagTide;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"fixed_k", "fixed_k", &Self::p_fixed_k, 0.0, c_ParamBounds::NonNegative,
             "Static potential Love number k_l of each degree, from l = 2."},
            {"fixed_dt", "fixed_dt_s", &Self::p_fixed_dt, 0.0, c_ParamBounds::NonNegative,
             "Tidal time lag dt_l of each degree [s], from l = 2."},
        };
        return specs;
    }

    c_FixedLagTide() : c_FixedLagTide(c_ParamMap{}) {}
    explicit c_FixedLagTide(const c_ParamMap& params) : c_AnalyticTide("fixed_dt") { this->p_initialize(params); }

    double get_fixed_k(int degree_l) const override { return c_tide_degree_value(this->p_fixed_k, degree_l); }
    double get_fixed_dt(int degree_l) const override { return c_tide_degree_value(this->p_fixed_dt, degree_l); }

    c_LoveNumbers calc_love_numbers(
            int degree_l, double frequency, const c_LoveNumbers& /*solver_love*/) const override {
        const double k_l  = this->get_fixed_k(degree_l);
        const double dt_l = this->get_fixed_dt(degree_l);
        return p_analytic_love(k_l * c_constant_time_lag(frequency, dt_l));
    }

protected:
    std::vector<double> p_fixed_k;
    std::vector<double> p_fixed_dt;
};

// c_CTLQTide: constant time lag with a quality factor (alias "ctl_q" / "fixed_dt_q"):
//
//   k_l(omega) = k_l * (1 - i * |omega| * dt_l / Q_l)  ->  -Im[k_l] = k_l * |omega| * dt_l / Q_l
class c_CTLQTide final : public c_AnalyticTide<c_CTLQTide> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::CTLQTide;

    static const std::vector<c_ParamSpec<c_CTLQTide>>& parameter_specs() {
        using Self = c_CTLQTide;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"fixed_k", "fixed_k", &Self::p_fixed_k, 0.0, c_ParamBounds::NonNegative,
             "Static potential Love number k_l of each degree, from l = 2."},
            {"fixed_dt", "fixed_dt_s", &Self::p_fixed_dt, 0.0, c_ParamBounds::NonNegative,
             "Tidal time lag dt_l of each degree [s], from l = 2."},
            {"fixed_q", "fixed_q", &Self::p_fixed_q, 0.0, c_ParamBounds::NonNegative,
             "Tidal quality factor Q_l of each degree, from l = 2; 0 is no dissipation."},
        };
        return specs;
    }

    c_CTLQTide() : c_CTLQTide(c_ParamMap{}) {}
    explicit c_CTLQTide(const c_ParamMap& params) : c_AnalyticTide("fixed_dt_q") { this->p_initialize(params); }

    double get_fixed_k(int degree_l) const override { return c_tide_degree_value(this->p_fixed_k, degree_l); }
    double get_fixed_dt(int degree_l) const override { return c_tide_degree_value(this->p_fixed_dt, degree_l); }
    double get_fixed_q(int degree_l) const override { return c_tide_degree_value(this->p_fixed_q, degree_l); }

    c_LoveNumbers calc_love_numbers(
            int degree_l, double frequency, const c_LoveNumbers& /*solver_love*/) const override {
        const double k_l  = this->get_fixed_k(degree_l);
        const double dt_l = this->get_fixed_dt(degree_l);
        const double q_l  = this->get_fixed_q(degree_l);
        if (std::abs(q_l) <= TidalPyConstants::d_EPS) {
            return p_analytic_love(std::complex<double>(k_l, 0.0));
        }
        return p_analytic_love(k_l * c_constant_time_lag(frequency, dt_l, q_l));
    }

protected:
    std::vector<double> p_fixed_k;
    std::vector<double> p_fixed_dt;
    std::vector<double> p_fixed_q;
};

inline const c_ModelRegistry<c_TideBase>& c_tide_registry() {
    static const c_ModelRegistry<c_TideBase> registry = {
        {{"rheology"},
         BinaryClassID::RheologyTide, &c_make_entry<c_TideBase, c_RheologyTide>},
        {{"fixed_q", "cpl", "constant_phase_lag"},
         BinaryClassID::FixedQTide,   &c_make_entry<c_TideBase, c_FixedQTide>},
        {{"fixed_dt", "ctl", "constant_time_lag"},
         BinaryClassID::FixedLagTide, &c_make_entry<c_TideBase, c_FixedLagTide>},
        {{"fixed_dt_q", "ctl_q", "constant_time_lag_and_q"},
         BinaryClassID::CTLQTide,     &c_make_entry<c_TideBase, c_CTLQTide>},
    };
    return registry;
}

// The family's entry points, each one line over the generic registry functions.
inline std::unique_ptr<c_TideBase> c_find_tide(const std::string& model_name, const c_ParamMap& params) {
    return c_make_model(c_tide_registry(), model_name, params);
}

inline std::unique_ptr<c_TideBase> c_tide_from_binary(std::istream& in, bool force = false) {
    return c_model_from_binary(c_tide_registry(), in, force);
}

inline std::string c_tide_canonical_name(const std::string& model_name) {
    return c_canonical_model_name(c_tide_registry(), model_name);
}

inline std::vector<std::string> c_tide_model_names() {
    return c_model_names(c_tide_registry());
}

}  // namespace tidalpy
