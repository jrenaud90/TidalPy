#pragma once
/*
 * tide_.hpp - TidalPy global (1D) tidal dissipation models: c_RheologyTide (alias "rheology"),
 * c_FixedQTide ("cpl"/"fixed_q"), c_FixedLagTide ("ctl"/"fixed_dt"), and c_CTLQTide
 * ("ctl_q"/"fixed_dt_q").
 *
 * Each returns the complex Love number k_l at a tidal frequency; the collapse (tide_collapse_.hpp)
 * uses -Im[k_l] as the per-mode dissipation multiplier. Fixed per-degree parameters (k_l, Q_l, dt_l)
 * live in fixed-size slots indexed by (degree_l - 2) for l = 2..10, the range the eccentricity and
 * obliquity tables cover; a slot left at 0 means no contribution at that degree. All quantities MKS;
 * frequencies in rad s-1.
 *
 * Binary payload: the model name followed by the per-degree slots as doubles.
 */

#include <algorithm>
#include <array>
#include <cctype>
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
#include "../../Utilities_x/classes_x/model_names_.hpp"   // c_to_lower

namespace tidalpy {

// The 3D methods of c_RheologyTide are declared here but defined in structures_x/worlds/world_tides_.hpp,
// compiled into the world extension alone, so this header never pulls in the world, potential-engine, or
// kernel headers. The incomplete type is legal in a declaration and complete at each definition.
class c_LayeredWorld;

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

// Combined construction parameters for all tide models. The per-degree vectors are indexed from l = 2 and
// may be short; missing entries default to 0.
struct c_TideModelConfig {
    std::vector<double> fixed_k;   // static potential Love numbers k_l  [dimensionless]
    std::vector<double> fixed_q;   // tidal quality factors Q_l          [dimensionless]
    std::vector<double> fixed_dt;  // tidal time lags dt_l               [s]
};

// Copy a possibly short or over-long per-degree config vector into a fixed 9-slot array.
// Throws std::invalid_argument for more values than tabulated degrees or for a value that is negative or not
// finite: every per-degree parameter (k, Q, dt) is non-negative, and a zero Q means no dissipation.
inline void c_tide_fill_degree_slots(
        const std::vector<double>& src, std::array<double, C_TIDE_NUM_DEGREES>& dst) {
    if (src.size() > static_cast<std::size_t>(C_TIDE_NUM_DEGREES)) {
        throw std::invalid_argument(
            "TidalPy: a per-degree tide parameter list holds " + std::to_string(src.size()) +
            " values, more than the " + std::to_string(C_TIDE_NUM_DEGREES) + " degrees (l = 2 to 10) it covers.");
    }
    dst.fill(0.0);
    for (std::size_t i = 0; i < src.size(); ++i) {
        if (!(std::isfinite(src[i]) && (src[i] >= 0.0))) {
            throw std::invalid_argument(
                "TidalPy: per-degree tide parameters (fixed_k, fixed_q, fixed_dt_s) must be finite and not "
                "negative; got " + std::to_string(src[i]) + " at degree " + std::to_string(i + C_TIDE_MIN_DEGREE) +
                ".");
        }
        dst[i] = src[i];
    }
}

// 0 for an out-of-range degree.
inline double c_tide_degree_value(const std::array<double, C_TIDE_NUM_DEGREES>& arr, int degree_l) {
    const int idx = degree_l - C_TIDE_MIN_DEGREE;
    if (idx < 0 || idx >= C_TIDE_NUM_DEGREES) {
        return 0.0;
    }
    return arr[idx];
}


// k_l supplied by the radial solver (alias "rheology"); a pass-through that flags the world to run it.
class c_RheologyTide : public c_TideBase {
public:
    c_RheologyTide() : c_TideBase("rheology") {}
    explicit c_RheologyTide(const c_TideModelConfig& /*cfg*/) : c_TideBase("rheology") {}
    ~c_RheologyTide() override = default;

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
            c_LayeredWorld& world,
            const c_TideSolveConfig& state,
            double radius,
            double colatitude) const;

    // Batch form over paired (radii[i], colatitudes[i]) points. The wave list is built once and the radial
    // solve amortized across the points, since it depends on (l, |omega|) alone; points sharing a
    // colatitude share its angular work, and the colatitudes run on up to num_threads threads.
    void calc_3d_tidal_heating_batch(
            c_LayeredWorld& world,
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
            c_LayeredWorld& world,
            const c_TideSolveConfig& state,
            const c_Grid3DAxes& axes,
            double* out_disp,
            int num_threads) const;

    // Instantaneous stress [Pa] and strain on the axes' grid, as 6 * nr * nth * nph * nt doubles ordered
    // (r, theta, phi, t, component) with the components rr, tt, pp, rt, rp, tp. Either output may be null
    // to skip it. NaN where no wave has a shear kernel: a radius with no depth-resolved solution, or a
    // liquid layer. Modes at zero forcing frequency, the permanent tide, are excluded.
    void calc_3d_stress_strain_grid(
            c_LayeredWorld& world,
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
            c_LayeredWorld& world,
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

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::RheologyTide));
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }
};

// The shared part of the analytic tide models: NumSlots per-degree parameter slots, the first always k_l, each with
// its config key. The slots run in one order through the config entries, the construction config, and the binary
// payload (the model name, then every slot's degrees in turn).
template <std::size_t NumSlots>
class c_AnalyticTide : public c_TideBase {
public:
    ~c_AnalyticTide() override = default;

    double get_fixed_k(int degree_l) const { return c_tide_degree_value(this->p_slots[0], degree_l); }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_TideBase::append_config_entries(out);
        // Per-degree slots for l = 2..10.
        for (std::size_t slot = 0; slot < NumSlots; ++slot) {
            out.push_back(c_config_doubles(
                this->p_config_keys[slot],
                std::vector<double>(this->p_slots[slot].begin(), this->p_slots[slot].end())));
        }
    }

    bool needs_radial_solve() const override { return false; }

    void write_binary(std::ostream& out) const override {
        std::vector<double> params;
        params.reserve(NumSlots * C_TIDE_NUM_DEGREES);
        for (const auto& slot : this->p_slots) {
            params.insert(params.end(), slot.begin(), slot.end());
        }
        this->write_physics_binary(out, static_cast<uint32_t>(this->p_class_id), params);
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, NumSlots * C_TIDE_NUM_DEGREES);
        for (std::size_t slot = 0; slot < NumSlots; ++slot) {
            for (int i = 0; i < C_TIDE_NUM_DEGREES; ++i) {
                this->p_slots[slot][i] = params[slot * C_TIDE_NUM_DEGREES + i];
            }
        }
    }

protected:
    c_AnalyticTide(
            const std::string& model_name,
            BinaryClassID class_id,
            const std::array<const char*, NumSlots>& config_keys) :
        c_TideBase(model_name),
        p_class_id(class_id),
        p_config_keys(config_keys) {}

    // Each slot from its list of the construction config, in slot order (c_tide_fill_degree_slots).
    void p_fill_slots(const std::array<const std::vector<double>*, NumSlots>& sources) {
        for (std::size_t slot = 0; slot < NumSlots; ++slot) {
            c_tide_fill_degree_slots(*sources[slot], this->p_slots[slot]);
        }
    }

    double p_slot_value(std::size_t slot, int degree_l) const {
        return c_tide_degree_value(this->p_slots[slot], degree_l);
    }

    // h and l are not defined for an analytic model (no radial solution).
    static c_LoveNumbers p_analytic_love(const std::complex<double>& k) {
        const std::complex<double> nan_love(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
        return c_LoveNumbers(k, nan_love, nan_love);
    }

    std::array<std::array<double, C_TIDE_NUM_DEGREES>, NumSlots> p_slots{};

private:
    BinaryClassID p_class_id;
    std::array<const char*, NumSlots> p_config_keys;
};

// c_FixedQTide: constant phase lag / fixed Q (alias "cpl" / "fixed_q"). Frequency-independent
// dissipation per degree:
//
//   k_l(omega) = k_l * (1 - i / Q_l)            ->  -Im[k_l] = k_l / Q_l
class c_FixedQTide : public c_AnalyticTide<2> {
public:
    c_FixedQTide() : c_AnalyticTide<2>("fixed_q", BinaryClassID::FixedQTide, {"fixed_k", "fixed_q"}) {}
    explicit c_FixedQTide(const c_TideModelConfig& cfg) : c_FixedQTide() {
        this->p_fill_slots({&cfg.fixed_k, &cfg.fixed_q});
    }
    ~c_FixedQTide() override = default;

    double get_fixed_q(int degree_l) const override { return this->p_slot_value(1, degree_l); }

    c_LoveNumbers calc_love_numbers(
            int degree_l, double /*frequency*/, const c_LoveNumbers& /*solver_love*/) const override {
        const double k_l = this->get_fixed_k(degree_l);
        const double q_l = this->p_slot_value(1, degree_l);
        if (std::abs(q_l) <= TidalPyConstants::d_EPS) {
            // An unset or zero Q_l is no dissipation, not a divide by zero.
            return p_analytic_love(std::complex<double>(k_l, 0.0));
        }
        return p_analytic_love(std::complex<double>(k_l, -k_l / q_l));
    }
};

// c_FixedLagTide: constant time lag / CTL (alias "ctl" / "fixed_dt"):
//
//   k_l(omega) = k_l * (1 - i * omega * dt_l)   ->  -Im[k_l] = k_l * omega * dt_l
class c_FixedLagTide : public c_AnalyticTide<2> {
public:
    c_FixedLagTide() : c_AnalyticTide<2>("fixed_dt", BinaryClassID::FixedLagTide, {"fixed_k", "fixed_dt_s"}) {}
    explicit c_FixedLagTide(const c_TideModelConfig& cfg) : c_FixedLagTide() {
        this->p_fill_slots({&cfg.fixed_k, &cfg.fixed_dt});
    }
    ~c_FixedLagTide() override = default;

    double get_fixed_dt(int degree_l) const override { return this->p_slot_value(1, degree_l); }

    c_LoveNumbers calc_love_numbers(
            int degree_l, double frequency, const c_LoveNumbers& /*solver_love*/) const override {
        const double k_l  = this->get_fixed_k(degree_l);
        const double dt_l = this->p_slot_value(1, degree_l);
        return p_analytic_love(std::complex<double>(k_l, -k_l * frequency * dt_l));
    }
};

// c_CTLQTide: constant time lag with a quality factor (alias "ctl_q" / "fixed_dt_q"):
//
//   k_l(omega) = k_l * (1 - i * omega * dt_l / Q_l)  ->  -Im[k_l] = k_l * omega * dt_l / Q_l
class c_CTLQTide : public c_AnalyticTide<3> {
public:
    c_CTLQTide() :
        c_AnalyticTide<3>("fixed_dt_q", BinaryClassID::CTLQTide, {"fixed_k", "fixed_dt_s", "fixed_q"}) {}
    explicit c_CTLQTide(const c_TideModelConfig& cfg) : c_CTLQTide() {
        this->p_fill_slots({&cfg.fixed_k, &cfg.fixed_dt, &cfg.fixed_q});
    }
    ~c_CTLQTide() override = default;

    double get_fixed_dt(int degree_l) const override { return this->p_slot_value(1, degree_l); }
    double get_fixed_q(int degree_l) const override { return this->p_slot_value(2, degree_l); }

    c_LoveNumbers calc_love_numbers(
            int degree_l, double frequency, const c_LoveNumbers& /*solver_love*/) const override {
        const double k_l  = this->get_fixed_k(degree_l);
        const double dt_l = this->p_slot_value(1, degree_l);
        const double q_l  = this->p_slot_value(2, degree_l);
        if (std::abs(q_l) <= TidalPyConstants::d_EPS) {
            return p_analytic_love(std::complex<double>(k_l, 0.0));
        }
        return p_analytic_love(std::complex<double>(k_l, -k_l * frequency * dt_l / q_l));
    }
};

enum class c_TideModel : uint8_t {
    Rheology = 0,
    FixedQ   = 1,
    FixedLag = 2,
    CTLQ     = 3,
};

// Model names are matched case-insensitively.
inline c_TideModel c_tide_model_from_name(const std::string& model_name) {
    const std::string name = c_to_lower(model_name);
    if (name == "rheology")                                                           { return c_TideModel::Rheology; }
    if (name == "cpl"   || name == "fixed_q"    || name == "constant_phase_lag")      { return c_TideModel::FixedQ; }
    if (name == "ctl"   || name == "fixed_dt"   || name == "constant_time_lag")       { return c_TideModel::FixedLag; }
    if (name == "ctl_q" || name == "fixed_dt_q" || name == "constant_time_lag_and_q") { return c_TideModel::CTLQ; }
    
    throw std::invalid_argument("TidalPy: unknown tide model name '" + model_name + "'");
}

inline std::unique_ptr<c_TideBase> c_find_tide(c_TideModel model, const c_TideModelConfig& cfg) {
    switch (model) {
        case c_TideModel::Rheology: return std::make_unique<c_RheologyTide>(cfg);
        case c_TideModel::FixedQ:   return std::make_unique<c_FixedQTide>(cfg);
        case c_TideModel::FixedLag: return std::make_unique<c_FixedLagTide>(cfg);
        case c_TideModel::CTLQ:     return std::make_unique<c_CTLQTide>(cfg);
    }
    throw std::invalid_argument("TidalPy: unrecognised c_TideModel enum value");
}

inline std::unique_ptr<c_TideBase> c_find_tide(const std::string& model_name, const c_TideModelConfig& cfg) {
    return c_find_tide(c_tide_model_from_name(model_name), cfg);
}

// The class id is peeked without consuming the header so the default-constructed model restores itself.
inline std::unique_ptr<c_TideBase> c_tide_from_binary(std::istream& in, bool force = false) {
    const c_BinaryHeader header = c_peek_binary_header(in);

    std::unique_ptr<c_TideBase> model;
    switch (static_cast<BinaryClassID>(header.class_id)) {
        case BinaryClassID::RheologyTide: model = std::make_unique<c_RheologyTide>(); break;
        case BinaryClassID::FixedQTide:   model = std::make_unique<c_FixedQTide>();   break;
        case BinaryClassID::FixedLagTide: model = std::make_unique<c_FixedLagTide>(); break;
        case BinaryClassID::CTLQTide:     model = std::make_unique<c_CTLQTide>();     break;
        default:
            throw std::runtime_error("TidalPy: unknown tide class id in binary stream");
    }
    model->read_binary(in, force);
    return model;
}

}  // namespace tidalpy
