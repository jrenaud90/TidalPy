#pragma once
/*
 * tide_.hpp — TidalPy global (1D) tidal dissipation models.
 *
 * Inherits c_TideBase (tide_base_.hpp) -> c_PhysicsBase. Each model returns the complex
 * Love number k_l at a tidal frequency; the collapse (tide_collapse_.hpp) uses
 * -Im[k_l] as the per-mode dissipation multiplier.
 *
 * Models (with config aliases handled by the factory):
 *   c_RheologyTide  (alias "rheology")              — k_l from the radial solver.
 *   c_FixedQTide    (alias "cpl"/"fixed_q")         — k_l*(1 - i/Q_l).
 *   c_FixedLagTide  (alias "ctl"/"fixed_dt")        — k_l*(1 - i*omega*dt_l).
 *   c_CTLQTide      (alias "ctl_q"/"fixed_dt_q")    — k_l*(1 - i*omega*dt_l/Q_l).
 *
 * Fixed per-degree parameters (k_l, Q_l, dt_l) are stored in fixed-size slots indexed by
 * (degree_l - 2) for l = 2..10 (the eccentricity/obliquity tables cover this range). A
 * value left at 0 means "no contribution at that degree". All quantities MKS; frequencies
 * in rad s-1.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::<Model> (901-904)
 *   payload: model_name length (uint32_t) | model_name bytes | model params (doubles)
 *   Rheology writes 0, FixedQ writes 18 (k[9]+q[9]), FixedLag writes 18 (k[9]+dt[9]),
 *   CTLQ writes 27 (k[9]+dt[9]+q[9]) scalars.
 */

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <complex>
#include <cstdint>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "constants_.hpp"     // TidalPyConstants::d_EPS
#include "tide_base_.hpp"

namespace tidalpy {

// Forward declarations for the on-demand 3D tidal-heating path. c_RheologyTide::calc_3d_tidal_heating is
// declared here but defined in structures_x/worlds/world_tides_.hpp (compiled only into the
// world extension), so this lightweight tide header never includes the world / potential / kernel
// headers. The method calls the world's members directly (no callbacks); references to these incomplete
// types are legal in the declaration, and the complete types are visible at the definition.
class c_LayeredWorld;
class c_TidalPotentialBase;
struct c_TidalPotentialState;

// Supported tidal degrees: l = 2..10 (matches the eccentricity/obliquity tables).
constexpr int C_TIDE_MIN_DEGREE  = 2;
constexpr int C_TIDE_MAX_DEGREE  = 10;
constexpr int C_TIDE_NUM_DEGREES = C_TIDE_MAX_DEGREE - C_TIDE_MIN_DEGREE + 1;  // 9

// -------------------------------------------------------------------------------
// c_TideModelConfig — combined construction parameters for all tide models.
//
// The per-degree vectors are indexed from l=2 (index 0 -> l=2, index 1 -> l=3, ...).
// They may be shorter than C_TIDE_NUM_DEGREES; missing entries default to 0.
// -------------------------------------------------------------------------------
struct c_TideModelConfig {
    std::vector<double> fixed_k;   // static potential Love numbers k_l  [dimensionless]
    std::vector<double> fixed_q;   // tidal quality factors Q_l          [dimensionless]
    std::vector<double> fixed_dt;  // tidal time lags dt_l               [s]
};

// Lower-case a model name for case-insensitive factory lookup.
inline std::string tide_to_lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

// Copy a (possibly short / over-long) per-degree config vector into a fixed 9-slot array.
inline void tide_fill_degree_slots(
        const std::vector<double>& src, std::array<double, C_TIDE_NUM_DEGREES>& dst) {
    dst.fill(0.0);
    const std::size_t n = std::min<std::size_t>(src.size(), C_TIDE_NUM_DEGREES);
    for (std::size_t i = 0; i < n; ++i) {
        dst[i] = src[i];
    }
}

// Fetch a per-degree slot (returns 0 for an out-of-range degree).
inline double tide_degree_value(const std::array<double, C_TIDE_NUM_DEGREES>& arr, int degree_l) {
    const int idx = degree_l - C_TIDE_MIN_DEGREE;
    if (idx < 0 || idx >= C_TIDE_NUM_DEGREES) {
        return 0.0;
    }
    return arr[idx];
}


// =====================================================================================================================
// Tide models
// =====================================================================================================================

// -------------------------------------------------------------------------------
// c_RheologyTide — k_l supplied by the radial solver (alias "rheology").
// The dissipation is whatever the viscoelastic-gravitational solution produced; this
// model is a thin pass-through that flags the world to run the radial solver.
// -------------------------------------------------------------------------------
class c_RheologyTide : public c_TideBase {
public:
    c_RheologyTide() : c_TideBase("rheology") {}
    explicit c_RheologyTide(const c_TideModelConfig& /*cfg*/) : c_TideBase("rheology") {}
    ~c_RheologyTide() override = default;

    c_LoveNumbers calc_love_numbers(
            int /*degree_l*/, double /*frequency*/, const c_LoveNumbers& solver_love) const override {
        // Pass the full radial-solver suite (k, h, l) straight through.
        return solver_love;
    }

    bool needs_radial_solve() const override { return true; }

    // On-demand 3D tidal volumetric heating [W m-3] at one point (radius [m], colatitude/longitude
    // [rad], time [s]). Only the rheology model supports the 3D path (it alone has the depth-resolved
    // radial solution). Loops the potential model's active modes, solving the world radial response once
    // per unique mode frequency, summing each mode's stress and strain tensors (each scaled by its own
    // freq_half = |omega|/2), and computing the heating once from the combined tensors (cross-mode terms
    // preserved). Defined out-of-line in structures_x/worlds/world_tides_.hpp. Returns NaN in liquid
    // layers / at the center / below the solver's starting radius.
    double calc_3d_tidal_heating(
            c_LayeredWorld& world,
            const c_TidalPotentialBase& potential,
            const c_TidalPotentialState& state,
            double radius,
            double colatitude,
            double longitude,
            double time) const;

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::RheologyTide));
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }
};

// -------------------------------------------------------------------------------
// c_FixedQTide — constant phase lag / fixed-Q (alias "cpl" / "fixed_q"):
//
//   k_l(omega) = k_l * (1 - i / Q_l)            ->  -Im[k_l] = k_l / Q_l
//
// Frequency-independent dissipation per degree.
// -------------------------------------------------------------------------------
class c_FixedQTide : public c_TideBase {
public:
    c_FixedQTide() : c_TideBase("fixed_q") {}
    explicit c_FixedQTide(const c_TideModelConfig& cfg) : c_TideBase("fixed_q") {
        tide_fill_degree_slots(cfg.fixed_k, this->p_fixed_k);
        tide_fill_degree_slots(cfg.fixed_q, this->p_fixed_q);
    }
    ~c_FixedQTide() override = default;

    double get_fixed_k(int degree_l) const { return tide_degree_value(this->p_fixed_k, degree_l); }
    double get_fixed_q(int degree_l) const { return tide_degree_value(this->p_fixed_q, degree_l); }

    c_LoveNumbers calc_love_numbers(
            int degree_l, double /*frequency*/, const c_LoveNumbers& /*solver_love*/) const override {
        const double k_l = tide_degree_value(this->p_fixed_k, degree_l);
        const double q_l = tide_degree_value(this->p_fixed_q, degree_l);
        std::complex<double> k;
        if (std::abs(q_l) <= TidalPyConstants::d_EPS) {
            // Q_l unset/zero -> treat as no dissipation (purely elastic) rather than divide by zero.
            k = std::complex<double>(k_l, 0.0);
        } else {
            k = std::complex<double>(k_l, -k_l / q_l);
        }
        // h, l are not defined for an analytic CPL model (no radial solution).
        const std::complex<double> nan_love(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
        return c_LoveNumbers(k, nan_love, nan_love);
    }

    bool needs_radial_solve() const override { return false; }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::FixedQTide),
                                   this->pack_params());
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> p = this->read_physics_binary(in, force, 2 * C_TIDE_NUM_DEGREES);
        for (int i = 0; i < C_TIDE_NUM_DEGREES; ++i) { this->p_fixed_k[i] = p[i]; }
        for (int i = 0; i < C_TIDE_NUM_DEGREES; ++i) { this->p_fixed_q[i] = p[C_TIDE_NUM_DEGREES + i]; }
    }

protected:
    std::array<double, C_TIDE_NUM_DEGREES> p_fixed_k{};
    std::array<double, C_TIDE_NUM_DEGREES> p_fixed_q{};

    std::vector<double> pack_params() const {
        std::vector<double> p;
        p.reserve(2 * C_TIDE_NUM_DEGREES);
        p.insert(p.end(), this->p_fixed_k.begin(), this->p_fixed_k.end());
        p.insert(p.end(), this->p_fixed_q.begin(), this->p_fixed_q.end());
        return p;
    }
};

// -------------------------------------------------------------------------------
// c_FixedLagTide — constant time lag / CTL (alias "ctl" / "fixed_dt"):
//
//   k_l(omega) = k_l * (1 - i * omega * dt_l)   ->  -Im[k_l] = k_l * omega * dt_l
// -------------------------------------------------------------------------------
class c_FixedLagTide : public c_TideBase {
public:
    c_FixedLagTide() : c_TideBase("fixed_dt") {}
    explicit c_FixedLagTide(const c_TideModelConfig& cfg) : c_TideBase("fixed_dt") {
        tide_fill_degree_slots(cfg.fixed_k, this->p_fixed_k);
        tide_fill_degree_slots(cfg.fixed_dt, this->p_fixed_dt);
    }
    ~c_FixedLagTide() override = default;

    double get_fixed_k(int degree_l) const  { return tide_degree_value(this->p_fixed_k, degree_l); }
    double get_fixed_dt(int degree_l) const { return tide_degree_value(this->p_fixed_dt, degree_l); }

    c_LoveNumbers calc_love_numbers(
            int degree_l, double frequency, const c_LoveNumbers& /*solver_love*/) const override {
        const double k_l  = tide_degree_value(this->p_fixed_k, degree_l);
        const double dt_l = tide_degree_value(this->p_fixed_dt, degree_l);
        const std::complex<double> k(k_l, -k_l * frequency * dt_l);
        // h, l are not defined for an analytic CTL model (no radial solution).
        const std::complex<double> nan_love(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
        return c_LoveNumbers(k, nan_love, nan_love);
    }

    bool needs_radial_solve() const override { return false; }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::FixedLagTide),
                                   this->pack_params());
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> p = this->read_physics_binary(in, force, 2 * C_TIDE_NUM_DEGREES);
        for (int i = 0; i < C_TIDE_NUM_DEGREES; ++i) { this->p_fixed_k[i] = p[i]; }
        for (int i = 0; i < C_TIDE_NUM_DEGREES; ++i) { this->p_fixed_dt[i] = p[C_TIDE_NUM_DEGREES + i]; }
    }

protected:
    std::array<double, C_TIDE_NUM_DEGREES> p_fixed_k{};
    std::array<double, C_TIDE_NUM_DEGREES> p_fixed_dt{};

    std::vector<double> pack_params() const {
        std::vector<double> p;
        p.reserve(2 * C_TIDE_NUM_DEGREES);
        p.insert(p.end(), this->p_fixed_k.begin(), this->p_fixed_k.end());
        p.insert(p.end(), this->p_fixed_dt.begin(), this->p_fixed_dt.end());
        return p;
    }
};

// -------------------------------------------------------------------------------
// c_CTLQTide — constant time lag with a quality factor (alias "ctl_q" / "fixed_dt_q"):
//
//   k_l(omega) = k_l * (1 - i * omega * dt_l / Q_l)  ->  -Im[k_l] = k_l * omega * dt_l / Q_l
// -------------------------------------------------------------------------------
class c_CTLQTide : public c_TideBase {
public:
    c_CTLQTide() : c_TideBase("fixed_dt_q") {}
    explicit c_CTLQTide(const c_TideModelConfig& cfg) : c_TideBase("fixed_dt_q") {
        tide_fill_degree_slots(cfg.fixed_k, this->p_fixed_k);
        tide_fill_degree_slots(cfg.fixed_dt, this->p_fixed_dt);
        tide_fill_degree_slots(cfg.fixed_q, this->p_fixed_q);
    }
    ~c_CTLQTide() override = default;

    double get_fixed_k(int degree_l) const  { return tide_degree_value(this->p_fixed_k, degree_l); }
    double get_fixed_dt(int degree_l) const { return tide_degree_value(this->p_fixed_dt, degree_l); }
    double get_fixed_q(int degree_l) const  { return tide_degree_value(this->p_fixed_q, degree_l); }

    c_LoveNumbers calc_love_numbers(
            int degree_l, double frequency, const c_LoveNumbers& /*solver_love*/) const override {
        const double k_l  = tide_degree_value(this->p_fixed_k, degree_l);
        const double dt_l = tide_degree_value(this->p_fixed_dt, degree_l);
        const double q_l  = tide_degree_value(this->p_fixed_q, degree_l);
        std::complex<double> k;
        if (std::abs(q_l) <= TidalPyConstants::d_EPS) {
            k = std::complex<double>(k_l, 0.0);
        } else {
            k = std::complex<double>(k_l, -k_l * frequency * dt_l / q_l);
        }
        // h, l are not defined for an analytic CTL+Q model (no radial solution).
        const std::complex<double> nan_love(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
        return c_LoveNumbers(k, nan_love, nan_love);
    }

    bool needs_radial_solve() const override { return false; }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::CTLQTide),
                                   this->pack_params());
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> p = this->read_physics_binary(in, force, 3 * C_TIDE_NUM_DEGREES);
        for (int i = 0; i < C_TIDE_NUM_DEGREES; ++i) { this->p_fixed_k[i]  = p[i]; }
        for (int i = 0; i < C_TIDE_NUM_DEGREES; ++i) { this->p_fixed_dt[i] = p[C_TIDE_NUM_DEGREES + i]; }
        for (int i = 0; i < C_TIDE_NUM_DEGREES; ++i) { this->p_fixed_q[i]  = p[2 * C_TIDE_NUM_DEGREES + i]; }
    }

protected:
    std::array<double, C_TIDE_NUM_DEGREES> p_fixed_k{};
    std::array<double, C_TIDE_NUM_DEGREES> p_fixed_dt{};
    std::array<double, C_TIDE_NUM_DEGREES> p_fixed_q{};

    std::vector<double> pack_params() const {
        std::vector<double> p;
        p.reserve(3 * C_TIDE_NUM_DEGREES);
        p.insert(p.end(), this->p_fixed_k.begin(),  this->p_fixed_k.end());
        p.insert(p.end(), this->p_fixed_dt.begin(), this->p_fixed_dt.end());
        p.insert(p.end(), this->p_fixed_q.begin(),  this->p_fixed_q.end());
        return p;
    }
};

// =====================================================================================================================
// Factory
// =====================================================================================================================

enum class c_TideModel : uint8_t {
    Rheology = 0,
    FixedQ   = 1,
    FixedLag = 2,
    CTLQ     = 3,
};

// Map a (case-insensitive) model name or alias to a c_TideModel enum value.
// Throws std::invalid_argument on an unknown name.
inline c_TideModel c_tide_model_from_name(const std::string& model_name) {
    const std::string name = tide_to_lower(model_name);
    if (name == "rheology")                                                           { return c_TideModel::Rheology; }
    if (name == "cpl"   || name == "fixed_q"    || name == "constant_phase_lag")      { return c_TideModel::FixedQ; }
    if (name == "ctl"   || name == "fixed_dt"   || name == "constant_time_lag")       { return c_TideModel::FixedLag; }
    if (name == "ctl_q" || name == "fixed_dt_q" || name == "constant_time_lag_and_q") { return c_TideModel::CTLQ; }
    
    throw std::invalid_argument("TidalPy: unknown tide model name '" + model_name + "'");
}

// Build the tide model named by the enum; returns an owning unique_ptr.
inline std::unique_ptr<c_TideBase> c_find_tide(c_TideModel model, const c_TideModelConfig& cfg) {
    switch (model) {
        case c_TideModel::Rheology: return std::make_unique<c_RheologyTide>(cfg);
        case c_TideModel::FixedQ:   return std::make_unique<c_FixedQTide>(cfg);
        case c_TideModel::FixedLag: return std::make_unique<c_FixedLagTide>(cfg);
        case c_TideModel::CTLQ:     return std::make_unique<c_CTLQTide>(cfg);
    }
    throw std::invalid_argument("TidalPy: unrecognised c_TideModel enum value");
}

// Name overload.
inline std::unique_ptr<c_TideBase> c_find_tide(const std::string& model_name, const c_TideModelConfig& cfg) {
    return c_find_tide(c_tide_model_from_name(model_name), cfg);
}

// Reconstruct a tide model from a binary stream (peek class id -> build -> read).
inline std::unique_ptr<c_TideBase> c_tide_from_binary(std::istream& in, bool force = false) {
    const std::streampos start = in.tellg();
    const c_BinaryHeader header = read_binary_header(in);
    in.seekg(start);

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
