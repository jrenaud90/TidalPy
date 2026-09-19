#pragma once
/*
 * material_eos_.hpp: material equation-of-state (EOS) models.
 *
 * A model returns a layer material's density [kg/m^3] from the local pressure [Pa] (analytic models) or
 * radius [m] (interpolated model); the whole-planet EOS solve evaluates it inline while integrating the
 * structure ODE.
 *
 * Every model carries a thermal expansivity alpha0 [1/K] and a reference temperature T_ref [K]. Birch-Murnaghan
 * and Vinet add the thermal pressure alpha0 K0 (T - T_ref) to their cold pressure law (alpha K_T taken constant,
 * its high-temperature limit); the constant and interpolated models scale their density by
 * exp(-alpha0 (T - T_ref)). A zero expansivity (the default) or a non-finite temperature gives the athermal EOS.
 *
 * Models (factory aliases): c_ConstantDensityEOS ("constant", "uniform"), c_BirchMurnaghanEOS ("bm",
 * "birch_murnaghan"), c_VinetEOS ("vinet"), c_InterpolatedEOS ("interpolate", "interp").
 *
 * References
 * ----------
 * Birch (1947), Phys. Rev. 71, 809. Vinet et al. (1987), J. Geophys. Res. 92, 9319. Anderson (1995), Equations of
 * State of Solids for Geophysics and Ceramic Science (thermal pressure).
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
#include "../../constants_.hpp"                    // TidalPyConstants::d_EPS
#include "../../Utilities_x/math_x/numerics_.hpp"  // c_safe_pow, c_safe_exp

namespace tidalpy {

// Defaults for the safeguarded Newton/bisection density-from-pressure inversion (c_MaterialEOSConfig). The cap
// only guarantees termination; convergence normally takes well under 10 iterations.
inline constexpr double d_EOS_INVERT_RTOL      = 1.0e-13;
inline constexpr int    d_EOS_INVERT_MAX_ITERS = 60;

// Default reference temperature of the thermal terms [K]: ambient, where mineral-physics rho0 and K0 are quoted.
inline constexpr double d_EOS_REFERENCE_TEMPERATURE = 300.0;

// Combined construction parameters for all EOS models; each model reads only the fields it needs.
struct c_MaterialEOSConfig {
    double reference_density      = 3500.0;    // rho0 [kg/m^3]
    double reference_bulk_modulus = 1.0e11;    // K0   [Pa]
    double bulk_modulus_derivative = 4.0;       // K0'  [dimensionless]

    // Thermal terms, read by every model. A zero expansivity is the athermal EOS.
    double thermal_expansion     = 0.0;                          // alpha0 [1/K]
    double reference_temperature = d_EOS_REFERENCE_TEMPERATURE;  // T_ref [K], where rho0 and K0 apply

    double invert_rtol      = d_EOS_INVERT_RTOL;       // relative convergence tol on eta
    int    invert_max_iters = d_EOS_INVERT_MAX_ITERS;  // termination-safeguard cap

    // Interpolated model: sorted-ascending radius [m] and matching density [kg/m^3].
    std::vector<double> radius;
    std::vector<double> density;
    // Interpolated model, optional radius-varying static moduli [Pa] and viscosities [Pa s]; an empty table means
    // "not provided" (the world solve falls back to the layer constant). Non-empty tables must match radius.
    std::vector<double> shear_modulus;
    std::vector<double> bulk_modulus;
    std::vector<double> shear_viscosity;
    std::vector<double> bulk_viscosity;
};

// =====================================================================================================================
// Analytic pressure laws and the density-from-pressure inversion
// =====================================================================================================================
// All laws are written in the compression ratio eta = rho / rho0 = V0 / V. They increase monotonically in eta only
// near eta = 1 (the finite-strain corrections turn them over at extreme eta), so the inversion brackets its root
// within the monotonic range.

// 3rd-order Birch-Murnaghan pressure [Pa] at compression eta = rho/rho0.
inline double eos_bm_pressure(double eta, double K0, double K0_prime) noexcept {
    const double eta_23 = c_safe_pow(eta, 2.0 / 3.0);
    const double eta_53 = c_safe_pow(eta, 5.0 / 3.0);
    const double eta_73 = c_safe_pow(eta, 7.0 / 3.0);
    return 1.5 * K0 * (eta_73 - eta_53) * (1.0 + 0.75 * (K0_prime - 4.0) * (eta_23 - 1.0));
}

// Vinet pressure [Pa] at compression eta = rho/rho0 (inv_cbrt_eta = (V/V0)^{1/3} = eta^{-1/3}).
inline double eos_vinet_pressure(double eta, double K0, double K0_prime) noexcept {
    const double inv_cbrt_eta = c_safe_pow(eta, -1.0 / 3.0);
    return 3.0 * K0 * (1.0 - inv_cbrt_eta) / (inv_cbrt_eta * inv_cbrt_eta)
        * c_safe_exp(1.5 * (K0_prime - 1.0) * (1.0 - inv_cbrt_eta));
}

// Isothermal bulk modulus K = eta dP/deta [Pa] of the 3rd-order Birch-Murnaghan law at compression eta.
inline double eos_bm_bulk_modulus(double eta, double K0, double K0_prime) noexcept {
    const double eta_23 = c_safe_pow(eta, 2.0 / 3.0);
    const double eta_53 = c_safe_pow(eta, 5.0 / 3.0);
    const double eta_73 = c_safe_pow(eta, 7.0 / 3.0);
    const double strain_coeff = 0.75 * (K0_prime - 4.0);
    return 1.5 * K0 * (
        ((7.0 / 3.0) * eta_73 - (5.0 / 3.0) * eta_53) * (1.0 + strain_coeff * (eta_23 - 1.0))
        + (eta_73 - eta_53) * strain_coeff * (2.0 / 3.0) * eta_23);
}

// Isothermal bulk modulus K = eta dP/deta [Pa] of the Vinet law at compression eta.
inline double eos_vinet_bulk_modulus(double eta, double K0, double K0_prime) noexcept {
    const double inv_cbrt_eta = c_safe_pow(eta, -1.0 / 3.0);
    const double exponent_coeff = 1.5 * (K0_prime - 1.0);
    return K0 * c_safe_exp(exponent_coeff * (1.0 - inv_cbrt_eta))
        * (2.0 - inv_cbrt_eta + exponent_coeff * inv_cbrt_eta * (1.0 - inv_cbrt_eta))
        / (inv_cbrt_eta * inv_cbrt_eta);
}

// Invert a pressure law for the compression eta = rho/rho0 at a target pressure by safeguarded Newton iteration
// with bisection fallback; PressureFn = double(eta, K0, K0'). The 3rd-order Birch-Murnaghan factor
// 1 + (3/4)(K0'-4)(eta^(2/3)-1) changes sign at large eta when K0' != 4, so P(eta) turns over there; the root is
// bracketed by expanding outward from eta = 1 (P = 0) and stopping at the turning point.
template <typename PressureFn>
inline double eos_invert_eta(
        double pressure_target,
        double K0,
        double K0_prime,
        PressureFn pressure_fn,
        double rtol,
        int max_iters) noexcept {
    if (std::abs(pressure_target) <= TidalPyConstants::d_EPS) { return 1.0; }

    double lo;
    double hi;
    if (pressure_target > 0.0) {
        // Compression: the solution has eta > 1. Grow hi while P keeps increasing.
        lo         = 1.0;
        hi         = 1.0;
        double p_hi = 0.0;  // P(eta = 1) = 0
        for (int k = 0; k < 200; ++k) {
            const double cand   = hi * 1.25;
            const double p_cand = pressure_fn(cand, K0, K0_prime);
            if (!(p_cand > p_hi)) { break; }  // reached the monotonic turning point
            hi   = cand;
            p_hi = p_cand;
            if (p_hi >= pressure_target) { break; }  // target now bracketed
        }
        if (pressure_target >= p_hi) { return hi; }  // beyond the model's valid range
    } else {
        // Tension: the solution has eta < 1. Shrink lo while P keeps decreasing.
        hi         = 1.0;
        lo         = 1.0;
        double p_lo = 0.0;
        for (int k = 0; k < 200; ++k) {
            const double cand   = lo * 0.8;
            const double p_cand = pressure_fn(cand, K0, K0_prime);
            if (!(p_cand < p_lo)) { break; }
            lo   = cand;
            p_lo = p_cand;
            if (p_lo <= pressure_target) { break; }
        }
        if (pressure_target <= p_lo) { return lo; }
    }

    double eta = 0.5 * (lo + hi);
    for (int i = 0; i < max_iters; ++i) {
        const double pressure = pressure_fn(eta, K0, K0_prime);
        if (pressure < pressure_target) { lo = eta; } else { hi = eta; }

        // Numerical derivative for the Newton step.
        const double fd_step        = 1.0e-6 * eta;
        const double pressure_slope = (pressure_fn(eta + fd_step, K0, K0_prime) - pressure) / fd_step;
        double next = (pressure_slope > TidalPyConstants::d_EPS)
            ? eta + (pressure_target - pressure) / pressure_slope
            : 0.5 * (lo + hi);
        if (!(next > lo && next < hi)) { next = 0.5 * (lo + hi); }  // safeguard

        // Converged once the compression stops changing to relative tolerance.
        if (std::abs(next - eta) <= rtol * eta) { 
            return next;
        }
        eta = next;
    }
    return eta;  // cap reached without full convergence; return the best estimate.
}

// Lower-case a model name for case-insensitive factory lookup.
inline std::string eos_to_lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

// =====================================================================================================================
// c_MaterialEOSBase: abstract base for all EOS models
// =====================================================================================================================
class c_MaterialEOSBase : public c_PhysicsBase {
public:
    explicit c_MaterialEOSBase(const std::string& model_name) : c_PhysicsBase(model_name) {}
    c_MaterialEOSBase(const std::string& model_name, const c_MaterialEOSConfig& cfg)
        : c_PhysicsBase(model_name),
          p_thermal_expansion(cfg.thermal_expansion),
          p_reference_temperature(cfg.reference_temperature) {}
    ~c_MaterialEOSBase() override = default;

    double get_thermal_expansion()     const noexcept { return this->p_thermal_expansion; }
    double get_reference_temperature() const noexcept { return this->p_reference_temperature; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_PhysicsBase::append_config_entries(out);
        out.push_back(c_config_double("thermal_expansion_1_k", this->p_thermal_expansion));
        out.push_back(c_config_double("reference_temperature_k", this->p_reference_temperature));
    }

    // Density [kg/m^3] from pressure [Pa], temperature [K], and radius [m]; analytic models use the pressure, the
    // interpolated model the radius. A non-finite temperature gives the athermal density.
    virtual double calc_density(
        double pressure, double temperature, double radius) const = 0;

    // Density and isothermal bulk modulus [Pa] together, so a model that inverts its pressure law does it once.
    // The bulk modulus is NaN ("not provided": the layer constant applies) unless the model defines one.
    virtual void calc_density_and_bulk_modulus(
            double pressure,
            double temperature,
            double radius,
            double& density,
            double& bulk_modulus) const {
        density      = this->calc_density(pressure, temperature, radius);
        bulk_modulus = this->calc_static_bulk_modulus(radius);
    }

    // Isothermal bulk modulus [Pa] at a pressure, temperature, and radius (see calc_density_and_bulk_modulus).
    double calc_bulk_modulus(double pressure, double temperature, double radius) const {
        double density      = TidalPyConstants::d_NAN;
        double bulk_modulus = TidalPyConstants::d_NAN;
        this->calc_density_and_bulk_modulus(pressure, temperature, radius, density, bulk_modulus);
        return bulk_modulus;
    }

    // Optional radius-varying static moduli [Pa] and viscosities [Pa s]. NaN means "not provided", and the world
    // solve then falls back to the layer constant; only c_InterpolatedEOS overrides these.
    virtual double calc_static_shear_modulus(double /*radius*/) const {
        return std::numeric_limits<double>::quiet_NaN();
    }
    virtual double calc_static_bulk_modulus(double /*radius*/) const {
        return std::numeric_limits<double>::quiet_NaN();
    }
    virtual double calc_shear_viscosity(double /*radius*/) const {
        return std::numeric_limits<double>::quiet_NaN();
    }
    virtual double calc_bulk_viscosity(double /*radius*/) const {
        return std::numeric_limits<double>::quiet_NaN();
    }

protected:
    // Temperature above the reference state [K]; zero (the athermal EOS) for a zero expansivity or a non-finite
    // temperature.
    double p_temperature_offset(double temperature) const noexcept {
        if (this->p_thermal_expansion == 0.0 || !std::isfinite(temperature)) { return 0.0; }
        return temperature - this->p_reference_temperature;
    }

    // Density factor exp(-alpha0 (T - T_ref)) of the models with no pressure law to carry a thermal pressure.
    double p_thermal_expansion_factor(double temperature) const noexcept {
        return c_safe_exp(-this->p_thermal_expansion * this->p_temperature_offset(temperature));
    }

    double p_thermal_expansion     = 0.0;
    double p_reference_temperature = d_EOS_REFERENCE_TEMPERATURE;
};

// Incompressible (uniform) density.
class c_ConstantDensityEOS : public c_MaterialEOSBase {
public:
    c_ConstantDensityEOS() : c_MaterialEOSBase("constant") {}
    explicit c_ConstantDensityEOS(const c_MaterialEOSConfig& cfg)
        : c_MaterialEOSBase("constant", cfg),
          p_reference_density(cfg.reference_density) {}
    ~c_ConstantDensityEOS() override = default;

    double get_reference_density() const noexcept { return this->p_reference_density; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_MaterialEOSBase::append_config_entries(out);
        out.push_back(c_config_double("reference_density_kg_m3", this->p_reference_density));
    }

    double calc_density(
            double /*pressure*/,
            double temperature,
            double /*radius*/) const override {
        return this->p_reference_density * this->p_thermal_expansion_factor(temperature);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::ConstantDensityEOS),
            {this->p_reference_density, this->p_thermal_expansion, this->p_reference_temperature});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 3);
        this->p_reference_density     = params[0];
        this->p_thermal_expansion     = params[1];
        this->p_reference_temperature = params[2];
    }

protected:
    double p_reference_density = 3500.0;
};

// 3rd-order Birch-Murnaghan, density from pressure.
class c_BirchMurnaghanEOS : public c_MaterialEOSBase {
public:
    c_BirchMurnaghanEOS() : c_MaterialEOSBase("birch_murnaghan") {}
    explicit c_BirchMurnaghanEOS(const c_MaterialEOSConfig& cfg)
        : c_MaterialEOSBase("birch_murnaghan", cfg),
          p_reference_density(cfg.reference_density),
          p_reference_bulk_modulus(cfg.reference_bulk_modulus),
          p_bulk_modulus_derivative(cfg.bulk_modulus_derivative),
          p_invert_rtol(cfg.invert_rtol),
          p_invert_max_iters(cfg.invert_max_iters) {}
    ~c_BirchMurnaghanEOS() override = default;

    double get_reference_density()       const noexcept { return this->p_reference_density; }
    double get_reference_bulk_modulus()  const noexcept { return this->p_reference_bulk_modulus; }
    double get_bulk_modulus_derivative() const noexcept { return this->p_bulk_modulus_derivative; }
    double get_invert_rtol()             const noexcept { return this->p_invert_rtol; }
    int    get_invert_max_iters()        const noexcept { return this->p_invert_max_iters; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_MaterialEOSBase::append_config_entries(out);
        out.push_back(c_config_double("reference_density_kg_m3", this->p_reference_density));
        out.push_back(c_config_double("reference_bulk_modulus_pa", this->p_reference_bulk_modulus));
        out.push_back(c_config_double("bulk_modulus_derivative", this->p_bulk_modulus_derivative));
        out.push_back(c_config_double("invert_rtol", this->p_invert_rtol));
        out.push_back(c_config_int("invert_max_iters", this->p_invert_max_iters));
    }

    double calc_density(
            double pressure,
            double temperature,
            double /*radius*/) const override {
        return this->p_reference_density * this->p_calc_compression(pressure, temperature);
    }

    void calc_density_and_bulk_modulus(
            double pressure,
            double temperature,
            double /*radius*/,
            double& density,
            double& bulk_modulus) const override {
        const double eta = this->p_calc_compression(pressure, temperature);
        density      = this->p_reference_density * eta;
        bulk_modulus = eos_bm_bulk_modulus(
            eta, this->p_reference_bulk_modulus, this->p_bulk_modulus_derivative);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::BirchMurnaghanEOS),
            {this->p_reference_density, this->p_reference_bulk_modulus,
             this->p_bulk_modulus_derivative, this->p_invert_rtol,
             static_cast<double>(this->p_invert_max_iters),
             this->p_thermal_expansion, this->p_reference_temperature});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 7);
        this->p_reference_density        = params[0];
        this->p_reference_bulk_modulus   = params[1];
        this->p_bulk_modulus_derivative  = params[2];
        this->p_invert_rtol              = params[3];
        this->p_invert_max_iters         = static_cast<int>(params[4]);
        this->p_thermal_expansion        = params[5];
        this->p_reference_temperature    = params[6];
    }

protected:
    // Compression eta = rho/rho0: the cold pressure law inverted at the pressure less the thermal pressure
    // alpha0 K0 (T - T_ref).
    double p_calc_compression(double pressure, double temperature) const noexcept {
        const double thermal_pressure = this->p_thermal_expansion * this->p_reference_bulk_modulus
            * this->p_temperature_offset(temperature);
        return eos_invert_eta(
            pressure - thermal_pressure,
            this->p_reference_bulk_modulus,
            this->p_bulk_modulus_derivative,
            eos_bm_pressure,
            this->p_invert_rtol,
            this->p_invert_max_iters);
    }

    double p_reference_density       = 3500.0;
    double p_reference_bulk_modulus  = 1.0e11;
    double p_bulk_modulus_derivative = 4.0;
    double p_invert_rtol             = d_EOS_INVERT_RTOL;
    int    p_invert_max_iters        = d_EOS_INVERT_MAX_ITERS;
};

// Vinet (universal) EOS, density from pressure.
class c_VinetEOS : public c_MaterialEOSBase {
public:
    c_VinetEOS() : c_MaterialEOSBase("vinet") {}
    explicit c_VinetEOS(const c_MaterialEOSConfig& cfg)
        : c_MaterialEOSBase("vinet", cfg),
          p_reference_density(cfg.reference_density),
          p_reference_bulk_modulus(cfg.reference_bulk_modulus),
          p_bulk_modulus_derivative(cfg.bulk_modulus_derivative),
          p_invert_rtol(cfg.invert_rtol),
          p_invert_max_iters(cfg.invert_max_iters) {}
    ~c_VinetEOS() override = default;

    double get_reference_density()       const noexcept { return this->p_reference_density; }
    double get_reference_bulk_modulus()  const noexcept { return this->p_reference_bulk_modulus; }
    double get_bulk_modulus_derivative() const noexcept { return this->p_bulk_modulus_derivative; }
    double get_invert_rtol()             const noexcept { return this->p_invert_rtol; }
    int    get_invert_max_iters()        const noexcept { return this->p_invert_max_iters; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_MaterialEOSBase::append_config_entries(out);
        out.push_back(c_config_double("reference_density_kg_m3", this->p_reference_density));
        out.push_back(c_config_double("reference_bulk_modulus_pa", this->p_reference_bulk_modulus));
        out.push_back(c_config_double("bulk_modulus_derivative", this->p_bulk_modulus_derivative));
        out.push_back(c_config_double("invert_rtol", this->p_invert_rtol));
        out.push_back(c_config_int("invert_max_iters", this->p_invert_max_iters));
    }

    double calc_density(
            double pressure,
            double temperature,
            double /*radius*/) const override {
        return this->p_reference_density * this->p_calc_compression(pressure, temperature);
    }

    void calc_density_and_bulk_modulus(
            double pressure,
            double temperature,
            double /*radius*/,
            double& density,
            double& bulk_modulus) const override {
        const double eta = this->p_calc_compression(pressure, temperature);
        density      = this->p_reference_density * eta;
        bulk_modulus = eos_vinet_bulk_modulus(
            eta, this->p_reference_bulk_modulus, this->p_bulk_modulus_derivative);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::VinetEOS),
            {this->p_reference_density, this->p_reference_bulk_modulus,
             this->p_bulk_modulus_derivative, this->p_invert_rtol,
             static_cast<double>(this->p_invert_max_iters),
             this->p_thermal_expansion, this->p_reference_temperature});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 7);
        this->p_reference_density        = params[0];
        this->p_reference_bulk_modulus   = params[1];
        this->p_bulk_modulus_derivative  = params[2];
        this->p_invert_rtol              = params[3];
        this->p_invert_max_iters         = static_cast<int>(params[4]);
        this->p_thermal_expansion        = params[5];
        this->p_reference_temperature    = params[6];
    }

protected:
    // Compression eta = rho/rho0: the cold pressure law inverted at the pressure less the thermal pressure
    // alpha0 K0 (T - T_ref).
    double p_calc_compression(double pressure, double temperature) const noexcept {
        const double thermal_pressure = this->p_thermal_expansion * this->p_reference_bulk_modulus
            * this->p_temperature_offset(temperature);
        return eos_invert_eta(
            pressure - thermal_pressure,
            this->p_reference_bulk_modulus,
            this->p_bulk_modulus_derivative,
            eos_vinet_pressure,
            this->p_invert_rtol,
            this->p_invert_max_iters);
    }

    double p_reference_density       = 3500.0;
    double p_reference_bulk_modulus  = 1.0e11;
    double p_bulk_modulus_derivative = 4.0;
    double p_invert_rtol             = d_EOS_INVERT_RTOL;
    int    p_invert_max_iters        = d_EOS_INVERT_MAX_ITERS;
};

// density(radius) lookup table (PREM-style profiles): linear interpolation in radius, clamped at the ends.
class c_InterpolatedEOS : public c_MaterialEOSBase {
public:
    c_InterpolatedEOS() : c_MaterialEOSBase("interpolate") {}
    explicit c_InterpolatedEOS(const c_MaterialEOSConfig& cfg)
        : c_MaterialEOSBase("interpolate", cfg),
          p_radius(cfg.radius),
          p_density(cfg.density),
          p_shear_modulus(cfg.shear_modulus),
          p_bulk_modulus(cfg.bulk_modulus),
          p_shear_viscosity(cfg.shear_viscosity),
          p_bulk_viscosity(cfg.bulk_viscosity)
    {
        // A table longer than the radius table would otherwise read past the end of the radius array.
        this->p_validate_tables();
    }
    ~c_InterpolatedEOS() override = default;

    std::size_t get_num_points() const noexcept { return this->p_radius.size(); }
    bool has_shear_modulus()   const noexcept { return !this->p_shear_modulus.empty(); }
    bool has_bulk_modulus()    const noexcept { return !this->p_bulk_modulus.empty(); }
    bool has_shear_viscosity() const noexcept { return !this->p_shear_viscosity.empty(); }
    bool has_bulk_viscosity()  const noexcept { return !this->p_bulk_viscosity.empty(); }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_MaterialEOSBase::append_config_entries(out);
        out.push_back(c_config_doubles("radius_m", this->p_radius));
        out.push_back(c_config_doubles("density_kg_m3", this->p_density));
        if (this->has_shear_modulus()) {
            out.push_back(c_config_doubles("shear_modulus_pa", this->p_shear_modulus));
        }
        if (this->has_bulk_modulus()) {
            out.push_back(c_config_doubles("bulk_modulus_pa", this->p_bulk_modulus));
        }
        if (this->has_shear_viscosity()) {
            out.push_back(c_config_doubles("shear_viscosity_pas", this->p_shear_viscosity));
        }
        if (this->has_bulk_viscosity()) {
            out.push_back(c_config_doubles("bulk_viscosity_pas", this->p_bulk_viscosity));
        }
    }

    double calc_density(
            double /*pressure*/,
            double temperature,
            double radius) const override {
        // NaN for an empty table.
        return this->p_thermal_expansion_factor(temperature) * c_interp(
            radius,
            this->p_radius.data(),
            this->p_density.data(),
            this->p_radius.size());
    }

    // Each returns NaN when its table is empty.
    double calc_static_shear_modulus(double radius) const override {
        return this->p_interp_optional(radius, this->p_shear_modulus);
    }
    double calc_static_bulk_modulus(double radius) const override {
        return this->p_interp_optional(radius, this->p_bulk_modulus);
    }
    double calc_shear_viscosity(double radius) const override {
        return this->p_interp_optional(radius, this->p_shear_viscosity);
    }
    double calc_bulk_viscosity(double radius) const override {
        return this->p_interp_optional(radius, this->p_bulk_viscosity);
    }

    void write_binary(std::ostream& out) const override {
        const auto n = static_cast<uint64_t>(this->p_radius.size());
        const uint64_t optional_count =
            (this->has_shear_modulus()   ? 1u : 0u) + (this->has_bulk_modulus()   ? 1u : 0u)
            + (this->has_shear_viscosity() ? 1u : 0u) + (this->has_bulk_viscosity() ? 1u : 0u);
        const uint64_t payload =
            binary_string_bytes(this->p_model_name)
            + sizeof(uint64_t)                       // point count
            + n * 2 * sizeof(double)                 // radius + density
            + 4 * sizeof(uint8_t)                    // 4 optional-array presence flags
            + optional_count * n * sizeof(double)    // present optional arrays
            + 2 * sizeof(double);                    // thermal expansivity + reference temperature
        write_binary_header(out, static_cast<uint32_t>(BinaryClassID::InterpolatedEOS), payload);
        write_binary_string(out, this->p_model_name);
        out.write(reinterpret_cast<const char*>(&n), sizeof(uint64_t));
        for (uint64_t i = 0; i < n; ++i) {
            out.write(reinterpret_cast<const char*>(&this->p_radius[i]),      sizeof(double));
            out.write(reinterpret_cast<const char*>(&this->p_density[i]), sizeof(double));
        }
        this->p_write_optional_array(out, this->p_shear_modulus);
        this->p_write_optional_array(out, this->p_bulk_modulus);
        this->p_write_optional_array(out, this->p_shear_viscosity);
        this->p_write_optional_array(out, this->p_bulk_viscosity);
        out.write(reinterpret_cast<const char*>(&this->p_thermal_expansion),     sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_reference_temperature), sizeof(double));
        if (!out) {
            throw std::runtime_error("TidalPy: failed to write interpolated EOS binary data");
        }
    }
    void read_binary(std::istream& in, bool force = false) override {
        c_TidalPyBaseClass::read_binary(in, force);
        this->p_model_name = read_binary_string(in);
        uint64_t n = 0;
        in.read(reinterpret_cast<char*>(&n), sizeof(uint64_t));
        this->p_radius.resize(n);
        this->p_density.resize(n);
        for (uint64_t i = 0; i < n; ++i) {
            in.read(reinterpret_cast<char*>(&this->p_radius[i]),      sizeof(double));
            in.read(reinterpret_cast<char*>(&this->p_density[i]), sizeof(double));
        }
        this->p_read_optional_array(in, this->p_shear_modulus, n);
        this->p_read_optional_array(in, this->p_bulk_modulus, n);
        this->p_read_optional_array(in, this->p_shear_viscosity, n);
        this->p_read_optional_array(in, this->p_bulk_viscosity, n);
        in.read(reinterpret_cast<char*>(&this->p_thermal_expansion),     sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_reference_temperature), sizeof(double));
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read interpolated EOS binary data");
        }
    }

protected:
    // Throw unless density and every non-empty optional table match the radius table in length.
    void p_validate_tables() const {
        const std::size_t num_points = this->p_radius.size();
        if (this->p_density.size() != num_points) {
            throw std::invalid_argument(
                "TidalPy: interpolated EOS density table length does not match its radius table.");
        }
        const std::vector<double>* optional_tables[4] = {
            &this->p_shear_modulus, &this->p_bulk_modulus,
            &this->p_shear_viscosity, &this->p_bulk_viscosity};
        const char* table_labels[4] = {
            "shear modulus", "bulk modulus", "shear viscosity", "bulk viscosity"};
        for (int table_i = 0; table_i < 4; ++table_i) {
            const std::vector<double>& table = *optional_tables[table_i];
            if (!table.empty() && table.size() != num_points) {
                throw std::invalid_argument(
                    std::string("TidalPy: interpolated EOS ") + table_labels[table_i] +
                    " table length does not match its radius table.");
            }
        }
    }

    // Interpolate an optional table vs radius; NaN if the table is empty.
    double p_interp_optional(double radius, const std::vector<double>& values) const {
        if (values.empty()) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return c_interp(radius, this->p_radius.data(), values.data(), values.size());
    }
    void p_write_optional_array(std::ostream& out, const std::vector<double>& values) const {
        const uint8_t present = values.empty() ? 0u : 1u;
        out.write(reinterpret_cast<const char*>(&present), sizeof(uint8_t));
        if (present) {
            for (double value : values) {
                out.write(reinterpret_cast<const char*>(&value), sizeof(double));
            }
        }
    }
    void p_read_optional_array(std::istream& in, std::vector<double>& values, uint64_t n) {
        uint8_t present = 0;
        in.read(reinterpret_cast<char*>(&present), sizeof(uint8_t));
        values.clear();
        if (present) {
            values.resize(n);
            for (uint64_t i = 0; i < n; ++i) {
                in.read(reinterpret_cast<char*>(&values[i]), sizeof(double));
            }
        }
    }

    std::vector<double> p_radius;
    std::vector<double> p_density;
    // Optional tables; empty means not provided.
    std::vector<double> p_shear_modulus;
    std::vector<double> p_bulk_modulus;
    std::vector<double> p_shear_viscosity;
    std::vector<double> p_bulk_viscosity;
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

// Map a case-insensitive model name or alias to the enum; throws std::invalid_argument on an unknown name.
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

// Build the EOS model named by the enum.
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

// Reconstruct an EOS model from a binary stream: peek the class id, build, read.
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
