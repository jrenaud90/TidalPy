#pragma once
/*
 * material_eos_.hpp — TidalPy material equation-of-state (EOS) models.
 *
 * A material EOS model returns a layer material's mass density [kg/m^3] given the
 * local pressure [Pa] (and, in future, temperature [K]) — the per-layer density
 * source consumed by the whole-planet radial EOS solve. Models follow the same
 * pattern as the rheology / cooling / radiogenics hierarchies: a c_PhysicsBase
 * subclass, an enum-based factory, and shared binary serialization.
 *
 * The radial solver integrates the structure ODE over radius and carries the
 * pressure as a state variable, so an analytic density(pressure) model can be
 * evaluated inline during the integration (no separate coupled iteration beyond
 * the solver's existing surface-pressure convergence loop).
 *
 * Models (with factory aliases):
 *   c_ConstantDensityEOS  (alias "constant"/"uniform")        — incompressible.
 *   c_BirchMurnaghanEOS   (alias "bm"/"birch_murnaghan")      — 3rd-order BM.
 *   c_VinetEOS            (alias "vinet")                     — Vinet/UBER.
 *   c_InterpolatedEOS     (alias "interpolate"/"interp")      — density(radius) table.
 *
 * All quantities are MKS. The analytic models are isothermal (temperature is
 * accepted for API uniformity but not yet used; thermal expansion is deferred).
 *
 * References
 * ----------
 * - Birch (1947), Phys. Rev. 71, 809 — finite-strain (Birch-Murnaghan) EOS.
 * - Vinet et al. (1987), J. Geophys. Res. 92, 9319 — universal (Vinet) EOS.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::<Model> (601-604)
 *   payload: model_name length (uint32_t) | model_name bytes | model params
 *   Constant/BM/Vinet write scalar doubles via the shared c_PhysicsBase helpers;
 *   Interpolated writes its variable-length radius/density arrays directly.
 *   The layer observer pointer (p_layer_ptr) is NOT serialized.
 */

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <istream>
#include <limits>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "physics_base_.hpp"
#include "interp_.hpp"

namespace tidalpy {

// -------------------------------------------------------------------------------
// Default settings for the analytic density-from-pressure inversion.
//
// These are NUMERICAL settings (not physical parameters): they bound the
// safeguarded Newton/bisection inversion used by the compressible models. They
// are exposed through c_MaterialEOSConfig so a caller can tune accuracy vs. cost
// per material. Convergence normally exits in well under ~10 iterations; the cap
// only guarantees termination on pathological input. The default cap is set
// generously above the worst-case pure-bisection count needed to shrink the
// [lo, hi] bracket below the relative tolerance.
// -------------------------------------------------------------------------------
inline constexpr double d_EOS_INVERT_RTOL      = 1.0e-13;
inline constexpr int    d_EOS_INVERT_MAX_ITERS = 60;

// -------------------------------------------------------------------------------
// c_MaterialEOSConfig — combined construction parameters for all EOS models.
// Each model reads only the fields it needs.
// -------------------------------------------------------------------------------
struct c_MaterialEOSConfig {
    double reference_density_kg_m3      = 3500.0;    // rho0 [kg/m^3]
    double reference_bulk_modulus_pa    = 1.0e11;    // K0   [Pa]
    double bulk_modulus_derivative      = 4.0;       // K0'  [dimensionless]

    // Numerical settings for the analytic density-from-pressure inversion.
    double invert_rtol      = d_EOS_INVERT_RTOL;       // relative convergence tol on eta
    int    invert_max_iters = d_EOS_INVERT_MAX_ITERS;  // termination-safeguard cap

    // Interpolated model: sorted-ascending radius [m] and matching density [kg/m^3].
    std::vector<double> radius_m;
    std::vector<double> density_kg_m3;
};

// =====================================================================================================================
// Analytic pressure laws and the density-from-pressure inversion
//
// All laws are written in terms of the compression ratio eta = rho / rho0 = V0 / V
// (mass conservation), and are monotonically increasing in eta, so the inversion
// (density given pressure) is well posed.
// =====================================================================================================================

// 3rd-order Birch-Murnaghan pressure [Pa] at compression eta = rho/rho0.
inline double eos_bm_pressure(double eta, double K0, double K0_prime) noexcept {
    const double f23 = std::pow(eta, 2.0 / 3.0);
    const double f53 = std::pow(eta, 5.0 / 3.0);
    const double f73 = std::pow(eta, 7.0 / 3.0);
    return 1.5 * K0 * (f73 - f53) * (1.0 + 0.75 * (K0_prime - 4.0) * (f23 - 1.0));
}

// Vinet pressure [Pa] at compression eta = rho/rho0 (x = (V/V0)^{1/3} = eta^{-1/3}).
inline double eos_vinet_pressure(double eta, double K0, double K0_prime) noexcept {
    const double x = std::pow(eta, -1.0 / 3.0);
    return 3.0 * K0 * (1.0 - x) / (x * x) * std::exp(1.5 * (K0_prime - 1.0) * (1.0 - x));
}

// Invert a monotonic pressure law for the compression eta = rho/rho0 given a
// target pressure, using a safeguarded Newton iteration (bisection fallback).
// Returns eta (caller multiplies by rho0 for density). PressureFn = double(eta, K0, K0').
//
// rtol is the relative convergence tolerance on eta; max_iters is a hard
// termination-safeguard cap (not a physical parameter). Both come from the model
// config (see c_MaterialEOSConfig). Convergence normally exits in well under ~10
// iterations.
template <typename PressureFn>
inline double eos_invert_eta(
        double pressure_target_pa, double K0, double K0_prime, PressureFn pressure_fn,
        double rtol, int max_iters) noexcept {
    if (pressure_target_pa == 0.0) { return 1.0; }

    // Bracket eta over a wide compression range; the laws are monotonic in eta.
    double lo = 1.0e-3;
    double hi = 1.0e3;
    const double p_lo = pressure_fn(lo, K0, K0_prime);
    const double p_hi = pressure_fn(hi, K0, K0_prime);
    if (pressure_target_pa <= p_lo) { return lo; }
    if (pressure_target_pa >= p_hi) { return hi; }

    double eta = 1.0;  // P(eta=1) == 0; a good starting point near typical solutions.
    for (int i = 0; i < max_iters; ++i) {
        const double p = pressure_fn(eta, K0, K0_prime);
        if (p < pressure_target_pa) { lo = eta; } else { hi = eta; }

        // Numerical derivative for the Newton step.
        const double h    = 1.0e-6 * eta;
        const double dp   = (pressure_fn(eta + h, K0, K0_prime) - p) / h;
        double next       = (dp > 0.0) ? eta + (pressure_target_pa - p) / dp : 0.5 * (lo + hi);
        if (!(next > lo && next < hi)) { next = 0.5 * (lo + hi); }  // safeguard

        // Converged once the compression stops changing to relative tolerance.
        if (std::abs(next - eta) <= rtol * eta) { 
            return next;
        }
        eta = next;
    }
    return eta;  // cap reached without full convergence; return the best estimate.
}

// -------------------------------------------------------------------------------
// Lower-case a model name for case-insensitive factory lookup.
// -------------------------------------------------------------------------------
inline std::string eos_to_lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

// =====================================================================================================================
// c_MaterialEOSBase — abstract base for all EOS models.
// =====================================================================================================================
class c_MaterialEOSBase : public c_PhysicsBase {
public:
    explicit c_MaterialEOSBase(const std::string& model_name) : c_PhysicsBase(model_name) {}
    ~c_MaterialEOSBase() override = default;

    // Mass density [kg/m^3] given local pressure [Pa], temperature [K], and radius
    // [m]. Analytic models use pressure; the interpolated model uses radius.
    virtual double calc_density(
        double pressure_pa, double temperature_k, double radius_m) const = 0;
};

// -------------------------------------------------------------------------------
// c_ConstantDensityEOS — incompressible (uniform) density (alias "constant").
// -------------------------------------------------------------------------------
class c_ConstantDensityEOS : public c_MaterialEOSBase {
public:
    c_ConstantDensityEOS() : c_MaterialEOSBase("constant") {}
    explicit c_ConstantDensityEOS(const c_MaterialEOSConfig& cfg)
        : c_MaterialEOSBase("constant"),
          p_reference_density(cfg.reference_density_kg_m3) {}
    ~c_ConstantDensityEOS() override = default;

    double get_reference_density() const noexcept { return this->p_reference_density; }

    double calc_density(
            double /*pressure_pa*/,
            double /*temperature_k*/,
            double /*radius_m*/) const override {
        return this->p_reference_density;
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::ConstantDensityEOS),
                                   {this->p_reference_density});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 1);
        this->p_reference_density = params[0];
    }

protected:
    double p_reference_density = 3500.0;
};

// -------------------------------------------------------------------------------
// c_BirchMurnaghanEOS — 3rd-order Birch-Murnaghan, density from pressure.
// -------------------------------------------------------------------------------
class c_BirchMurnaghanEOS : public c_MaterialEOSBase {
public:
    c_BirchMurnaghanEOS() : c_MaterialEOSBase("birch_murnaghan") {}
    explicit c_BirchMurnaghanEOS(const c_MaterialEOSConfig& cfg)
        : c_MaterialEOSBase("birch_murnaghan"),
          p_reference_density(cfg.reference_density_kg_m3),
          p_reference_bulk_modulus(cfg.reference_bulk_modulus_pa),
          p_bulk_modulus_derivative(cfg.bulk_modulus_derivative),
          p_invert_rtol(cfg.invert_rtol),
          p_invert_max_iters(cfg.invert_max_iters) {}
    ~c_BirchMurnaghanEOS() override = default;

    double get_reference_density()       const noexcept { return this->p_reference_density; }
    double get_reference_bulk_modulus()  const noexcept { return this->p_reference_bulk_modulus; }
    double get_bulk_modulus_derivative() const noexcept { return this->p_bulk_modulus_derivative; }
    double get_invert_rtol()             const noexcept { return this->p_invert_rtol; }
    int    get_invert_max_iters()        const noexcept { return this->p_invert_max_iters; }

    double calc_density(
            double pressure_pa,
            double /*temperature_k*/,
            double /*radius_m*/) const override {
        const double eta = eos_invert_eta(
            pressure_pa,
            this->p_reference_bulk_modulus,
            this->p_bulk_modulus_derivative,
            eos_bm_pressure,
            this->p_invert_rtol,
            this->p_invert_max_iters);
        return this->p_reference_density * eta;
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::BirchMurnaghanEOS),
            {this->p_reference_density, this->p_reference_bulk_modulus,
             this->p_bulk_modulus_derivative, this->p_invert_rtol,
             static_cast<double>(this->p_invert_max_iters)});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 5);
        this->p_reference_density        = params[0];
        this->p_reference_bulk_modulus   = params[1];
        this->p_bulk_modulus_derivative  = params[2];
        this->p_invert_rtol              = params[3];
        this->p_invert_max_iters         = static_cast<int>(params[4]);
    }

protected:
    double p_reference_density       = 3500.0;
    double p_reference_bulk_modulus  = 1.0e11;
    double p_bulk_modulus_derivative = 4.0;
    double p_invert_rtol             = d_EOS_INVERT_RTOL;
    int    p_invert_max_iters        = d_EOS_INVERT_MAX_ITERS;
};

// -------------------------------------------------------------------------------
// c_VinetEOS — Vinet (universal) EOS, density from pressure (alias "vinet").
// -------------------------------------------------------------------------------
class c_VinetEOS : public c_MaterialEOSBase {
public:
    c_VinetEOS() : c_MaterialEOSBase("vinet") {}
    explicit c_VinetEOS(const c_MaterialEOSConfig& cfg)
        : c_MaterialEOSBase("vinet"),
          p_reference_density(cfg.reference_density_kg_m3),
          p_reference_bulk_modulus(cfg.reference_bulk_modulus_pa),
          p_bulk_modulus_derivative(cfg.bulk_modulus_derivative),
          p_invert_rtol(cfg.invert_rtol),
          p_invert_max_iters(cfg.invert_max_iters) {}
    ~c_VinetEOS() override = default;

    double get_reference_density()       const noexcept { return this->p_reference_density; }
    double get_reference_bulk_modulus()  const noexcept { return this->p_reference_bulk_modulus; }
    double get_bulk_modulus_derivative() const noexcept { return this->p_bulk_modulus_derivative; }
    double get_invert_rtol()             const noexcept { return this->p_invert_rtol; }
    int    get_invert_max_iters()        const noexcept { return this->p_invert_max_iters; }

    double calc_density(
            double pressure_pa,
            double /*temperature_k*/,
            double /*radius_m*/) const override {
        const double eta = eos_invert_eta(
            pressure_pa,
            this->p_reference_bulk_modulus,
            this->p_bulk_modulus_derivative,
            eos_vinet_pressure,
            this->p_invert_rtol,
            this->p_invert_max_iters);
        return this->p_reference_density * eta;
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::VinetEOS),
            {this->p_reference_density, this->p_reference_bulk_modulus,
             this->p_bulk_modulus_derivative, this->p_invert_rtol,
             static_cast<double>(this->p_invert_max_iters)});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 5);
        this->p_reference_density        = params[0];
        this->p_reference_bulk_modulus   = params[1];
        this->p_bulk_modulus_derivative  = params[2];
        this->p_invert_rtol              = params[3];
        this->p_invert_max_iters         = static_cast<int>(params[4]);
    }

protected:
    double p_reference_density       = 3500.0;
    double p_reference_bulk_modulus  = 1.0e11;
    double p_bulk_modulus_derivative = 4.0;
    double p_invert_rtol             = d_EOS_INVERT_RTOL;
    int    p_invert_max_iters        = d_EOS_INVERT_MAX_ITERS;
};

// -------------------------------------------------------------------------------
// c_InterpolatedEOS — density(radius) lookup table (alias "interpolate").
//
// Holds sorted-ascending radius [m] and matching density [kg/m^3] arrays and
// returns density by linear interpolation in radius (clamped at the boundaries).
// This reproduces the legacy interpolation EOS (e.g. PREM Earth profiles).
// -------------------------------------------------------------------------------
class c_InterpolatedEOS : public c_MaterialEOSBase {
public:
    c_InterpolatedEOS() : c_MaterialEOSBase("interpolate") {}
    explicit c_InterpolatedEOS(const c_MaterialEOSConfig& cfg)
        : c_MaterialEOSBase("interpolate"),
          p_radius_m(cfg.radius_m),
          p_density_kg_m3(cfg.density_kg_m3) {}
    ~c_InterpolatedEOS() override = default;

    std::size_t get_num_points() const noexcept { return this->p_radius_m.size(); }

    double calc_density(
            double /*pressure_pa*/,
            double /*temperature_k*/,
            double radius_m) const override {
        // Linear interpolation in radius via the shared array utility (clamped at
        // the boundaries; NaN for an empty table).
        return c_interp(
            radius_m,
            this->p_radius_m.data(),
            this->p_density_kg_m3.data(),
            this->p_radius_m.size());
    }

    void write_binary(std::ostream& out) const override {
        const auto n = static_cast<uint64_t>(this->p_radius_m.size());
        const uint64_t payload =
            binary_string_bytes(this->p_model_name)
            + sizeof(uint64_t)                       // point count
            + n * 2 * sizeof(double);                // radius + density
        write_binary_header(out, static_cast<uint32_t>(BinaryClassID::InterpolatedEOS), payload);
        write_binary_string(out, this->p_model_name);
        out.write(reinterpret_cast<const char*>(&n), sizeof(uint64_t));
        for (uint64_t i = 0; i < n; ++i) {
            out.write(reinterpret_cast<const char*>(&this->p_radius_m[i]),      sizeof(double));
            out.write(reinterpret_cast<const char*>(&this->p_density_kg_m3[i]), sizeof(double));
        }
        if (!out) {
            throw std::runtime_error("TidalPy: failed to write interpolated EOS binary data");
        }
    }
    void read_binary(std::istream& in, bool force = false) override {
        c_TidalPyBaseClass::read_binary(in, force);
        this->p_model_name = read_binary_string(in);
        uint64_t n = 0;
        in.read(reinterpret_cast<char*>(&n), sizeof(uint64_t));
        this->p_radius_m.resize(n);
        this->p_density_kg_m3.resize(n);
        for (uint64_t i = 0; i < n; ++i) {
            in.read(reinterpret_cast<char*>(&this->p_radius_m[i]),      sizeof(double));
            in.read(reinterpret_cast<char*>(&this->p_density_kg_m3[i]), sizeof(double));
        }
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read interpolated EOS binary data");
        }
    }

protected:
    std::vector<double> p_radius_m;
    std::vector<double> p_density_kg_m3;
};

// =====================================================================================================================
// Factory
// =====================================================================================================================

enum class c_MaterialEOSModel : uint8_t {
    Constant       = 0,
    BirchMurnaghan = 1,
    Vinet          = 2,
    Interpolated   = 3,
};

// Map a (case-insensitive) model name or alias to a c_MaterialEOSModel enum value.
// Throws std::invalid_argument on an unknown name.
inline c_MaterialEOSModel c_material_eos_model_from_name(const std::string& model_name) {
    const std::string name = eos_to_lower(model_name);
    if (name == "constant" || name == "uniform" ||
        name == "constant_density")                  { return c_MaterialEOSModel::Constant; }
    if (name == "bm" || name == "birch_murnaghan" ||
        name == "birch-murnaghan")                   { return c_MaterialEOSModel::BirchMurnaghan; }
    if (name == "vinet")                             { return c_MaterialEOSModel::Vinet; }
    if (name == "interpolate" || name == "interp" ||
        name == "interpolated")                      { return c_MaterialEOSModel::Interpolated; }
    throw std::invalid_argument("TidalPy: unknown material EOS model name '" + model_name + "'");
}

// Build the EOS model named by the enum; returns an owning unique_ptr.
inline std::unique_ptr<c_MaterialEOSBase> c_find_material_eos(
        c_MaterialEOSModel model, const c_MaterialEOSConfig& cfg) {
    switch (model) {
        case c_MaterialEOSModel::Constant:       return std::make_unique<c_ConstantDensityEOS>(cfg);
        case c_MaterialEOSModel::BirchMurnaghan: return std::make_unique<c_BirchMurnaghanEOS>(cfg);
        case c_MaterialEOSModel::Vinet:          return std::make_unique<c_VinetEOS>(cfg);
        case c_MaterialEOSModel::Interpolated:   return std::make_unique<c_InterpolatedEOS>(cfg);
    }
    throw std::invalid_argument("TidalPy: unrecognised c_MaterialEOSModel enum value");
}

// Name overload.
inline std::unique_ptr<c_MaterialEOSBase> c_find_material_eos(
        const std::string& model_name, const c_MaterialEOSConfig& cfg) {
    return c_find_material_eos(c_material_eos_model_from_name(model_name), cfg);
}

// Reconstruct an EOS model from a binary stream (peek class id -> build -> read).
inline std::unique_ptr<c_MaterialEOSBase> c_material_eos_from_binary(std::istream& in, bool force = false) {
    const std::streampos start = in.tellg();
    const c_BinaryHeader header = read_binary_header(in);
    in.seekg(start);

    std::unique_ptr<c_MaterialEOSBase> model;
    switch (static_cast<BinaryClassID>(header.class_id)) {
        case BinaryClassID::ConstantDensityEOS: model = std::make_unique<c_ConstantDensityEOS>(); break;
        case BinaryClassID::BirchMurnaghanEOS:  model = std::make_unique<c_BirchMurnaghanEOS>();  break;
        case BinaryClassID::VinetEOS:           model = std::make_unique<c_VinetEOS>();           break;
        case BinaryClassID::InterpolatedEOS:    model = std::make_unique<c_InterpolatedEOS>();    break;
        default:
            throw std::runtime_error("TidalPy: unknown material EOS class id in binary stream");
    }
    model->read_binary(in, force);
    return model;
}

}  // namespace tidalpy
