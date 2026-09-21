#pragma once
/*
 * material_eos_.hpp: material equation-of-state (EOS) models.
 *
 * A model is a layer's material: it returns the density [kg/m^3] from the local pressure [Pa] (analytic models)
 * or radius [m] (interpolated model), and every other frequency-independent property with it, through
 * calc_material_state. The whole-planet EOS solve evaluates it inline while integrating the structure ODE, so
 * the solved structure is the one place those properties are read from afterwards.
 *
 * The material owns the static shear law mu = mu0 + mu'_P P + mu'_T (T - T_ref), the constant bulk modulus of
 * the models with no pressure law, the static viscosities, the optional viscosity and partial-melt models, and
 * the thermal constants (conductivity, heat capacity, and the expansivity the density law shares).
 * Nothing here depends on a forcing frequency: complex moduli are the rheology's job.
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
#include "binary_.hpp"                                   // write_optional_binary, read_optional_binary
#include "../../viscosity_x/viscosity_.hpp"              // c_ViscosityBase, c_viscosity_from_binary
#include "../../partial_melt_x/partial_melt_.hpp"        // c_PartialMeltBase, c_partial_melt_from_binary
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
    double reference_density       = 3500.0;    // rho0 [kg/m^3]
    double reference_bulk_modulus  = 1.0e11;    // K0   [Pa]
    double bulk_modulus_derivative = 4.0;       // K0'  [dimensionless]

    // Thermal terms, read by every model. A zero expansivity is the athermal EOS.
    double thermal_expansion     = 0.0;                          // alpha0 [1/K]
    double reference_temperature = d_EOS_REFERENCE_TEMPERATURE;  // T_ref [K], where rho0 and K0 apply

    double invert_rtol      = d_EOS_INVERT_RTOL;       // relative convergence tol on eta
    int    invert_max_iters = d_EOS_INVERT_MAX_ITERS;  // termination-safeguard cap

    // Static (unrelaxed) moduli [Pa] and viscosities [Pa s]. The bulk constant applies to a model with no
    // pressure law of its own; a NaN viscosity means unset (attach a viscosity model instead).
    double shear_modulus_static   = 0.0;
    double bulk_modulus_static    = 0.0;
    double shear_viscosity_static = std::numeric_limits<double>::quiet_NaN();
    double bulk_viscosity_static  = std::numeric_limits<double>::quiet_NaN();
    // Linear static shear law: mu = mu0 + mu'_P P + mu'_T (T - T_ref).
    double shear_modulus_pressure_derivative    = 0.0;                          // [dimensionless]
    double shear_modulus_temperature_derivative = 0.0;                          // [Pa/K]
    double shear_modulus_reference_temperature  = d_EOS_REFERENCE_TEMPERATURE;  // [K]
    // Thermal constants. The expansivity is thermal_expansion above: one alpha serves the density law, the
    // adiabat, and convection.
    double thermal_conductivity = 4.0;      // k   [W/(m K)]
    double heat_capacity        = 1200.0;   // c_p [J/(kg K)]

    // Interpolated model: sorted-ascending radius [m] and matching density [kg/m^3].
    std::vector<double> radius;
    std::vector<double> density;
    // Interpolated model, optional radius-varying static moduli [Pa] and viscosities [Pa s]; an empty table means
    // "not provided" (the material's law or constant applies). Non-empty tables must match radius.
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

// 3rd-order Birch-Murnaghan pressure [Pa] and isothermal bulk modulus K = eta dP/deta [Pa] at compression
// eta = rho/rho0. One cube root serves every fractional power, and the inversion wants both values at once.
inline void eos_bm_pressure_and_bulk_modulus(
        double eta,
        double K0,
        double K0_prime,
        double& pressure,
        double& bulk_modulus) noexcept {
    const double cbrt_eta     = std::cbrt(eta);
    const double eta_23       = cbrt_eta * cbrt_eta;
    const double eta_53       = eta * eta_23;
    const double eta_73       = eta_53 * eta_23;
    const double strain_coeff = 0.75 * (K0_prime - 4.0);
    const double strain_term  = 1.0 + strain_coeff * (eta_23 - 1.0);
    pressure     = 1.5 * K0 * (eta_73 - eta_53) * strain_term;
    bulk_modulus = 1.5 * K0 * (
        ((7.0 / 3.0) * eta_73 - (5.0 / 3.0) * eta_53) * strain_term
        + (eta_73 - eta_53) * strain_coeff * (2.0 / 3.0) * eta_23);
}

// Vinet pressure [Pa] and isothermal bulk modulus K = eta dP/deta [Pa] at compression eta = rho/rho0
// (inv_cbrt_eta = (V/V0)^{1/3} = eta^{-1/3}).
inline void eos_vinet_pressure_and_bulk_modulus(
        double eta,
        double K0,
        double K0_prime,
        double& pressure,
        double& bulk_modulus) noexcept {
    const double inv_cbrt_eta   = 1.0 / std::cbrt(eta);
    const double exponent_coeff = 1.5 * (K0_prime - 1.0);
    const double exponential    = c_safe_exp(exponent_coeff * (1.0 - inv_cbrt_eta));
    const double inv_square     = 1.0 / (inv_cbrt_eta * inv_cbrt_eta);
    pressure     = 3.0 * K0 * (1.0 - inv_cbrt_eta) * inv_square * exponential;
    bulk_modulus = K0 * exponential
        * (2.0 - inv_cbrt_eta + exponent_coeff * inv_cbrt_eta * (1.0 - inv_cbrt_eta)) * inv_square;
}

inline double eos_bm_pressure(double eta, double K0, double K0_prime) noexcept {
    double pressure;
    double bulk_modulus;
    eos_bm_pressure_and_bulk_modulus(eta, K0, K0_prime, pressure, bulk_modulus);
    return pressure;
}

inline double eos_vinet_pressure(double eta, double K0, double K0_prime) noexcept {
    double pressure;
    double bulk_modulus;
    eos_vinet_pressure_and_bulk_modulus(eta, K0, K0_prime, pressure, bulk_modulus);
    return pressure;
}

inline double eos_bm_bulk_modulus(double eta, double K0, double K0_prime) noexcept {
    double pressure;
    double bulk_modulus;
    eos_bm_pressure_and_bulk_modulus(eta, K0, K0_prime, pressure, bulk_modulus);
    return bulk_modulus;
}

inline double eos_vinet_bulk_modulus(double eta, double K0, double K0_prime) noexcept {
    double pressure;
    double bulk_modulus;
    eos_vinet_pressure_and_bulk_modulus(eta, K0, K0_prime, pressure, bulk_modulus);
    return bulk_modulus;
}

// The compressions over which a pressure law rises with compression, and the pressures at the two ends. Every law
// turns over in tension, and the 3rd-order Birch-Murnaghan factor 1 + (3/4)(K0'-4)(eta^(2/3)-1) changes sign at
// large eta when K0' < 4, so it turns over in compression too. An end the search never reaches stays unbounded.
// The range depends on the law's constants alone, so a model finds it once rather than at every inversion.
struct c_PressureLawRange {
    double compression_min = 0.0;
    double compression_max = TidalPyConstants::d_INF;
    double pressure_min    = -TidalPyConstants::d_INF;
    double pressure_max    = TidalPyConstants::d_INF;
};

// Find a law's monotonic range: step outward from eta = 1 until the bulk modulus K = eta dP/deta stops being
// positive, then bisect that sign change to rtol. LawFn fills the pressure and K at a compression.
template <typename LawFn>
inline c_PressureLawRange eos_find_monotonic_range(double K0, double K0_prime, LawFn law_fn, double rtol) noexcept {
    c_PressureLawRange range;
    double pressure = 0.0;
    double bulk     = 0.0;
    const auto rising = [&](double eta) {
        law_fn(eta, K0, K0_prime, pressure, bulk);
        return (bulk > 0.0) && std::isfinite(pressure);
    };
    for (int side = 0; side < 2; ++side) {
        const double growth = (side == 0) ? 0.8 : 1.25;
        double inside  = 1.0;
        double outside = 1.0;
        bool   bounded = false;
        for (int k = 0; k < 200; ++k) {
            const double candidate = inside * growth;
            if (!rising(candidate)) { outside = candidate; bounded = true; break; }
            inside = candidate;
        }
        if (!bounded) { continue; }
        for (int k = 0; (k < 200) && (std::abs(outside - inside) > rtol * inside); ++k) {
            const double middle = 0.5 * (inside + outside);
            if (rising(middle)) { inside = middle; } else { outside = middle; }
        }
        law_fn(inside, K0, K0_prime, pressure, bulk);
        if (side == 0) {
            range.compression_min = inside;
            range.pressure_min    = pressure;
        } else {
            range.compression_max = inside;
            range.pressure_max    = pressure;
        }
    }
    return range;
}

// Invert a pressure law for the compression eta = rho/rho0 at a target pressure. A target past either end of the
// monotonic range has no compression to find and takes that end, so the answer is continuous in the pressure; the
// structure solve depends on that while its central pressure is still a guess and its outer radii sit in tension.
// Inside the range this is Newton's method on the exact slope K/eta, started from the Murnaghan law (which
// inverts in closed form and tracks both laws closely over planetary compressions) and kept inside a bracket
// that every evaluation tightens; a step that leaves the bracket is replaced by its midpoint.
template <typename LawFn>
inline double eos_invert_eta(
        double pressure_target,
        double K0,
        double K0_prime,
        LawFn law_fn,
        const c_PressureLawRange& range,
        double rtol,
        int max_iters) noexcept {
    if (std::abs(pressure_target) <= TidalPyConstants::d_EPS) { return 1.0; }
    if (pressure_target <= range.pressure_min) { return range.compression_min; }
    if (pressure_target >= range.pressure_max) { return range.compression_max; }

    double lo = range.compression_min;
    double hi = range.compression_max;

    // Murnaghan: P = (K0 / K0') (eta^K0' - 1).
    const double murnaghan_base = 1.0 + K0_prime * pressure_target / K0;
    double eta = (murnaghan_base > 0.0 && K0_prime > 0.0) ? c_safe_pow(murnaghan_base, 1.0 / K0_prime) : 1.0;
    if (!(eta > lo && eta < hi)) { eta = 1.0; }

    double pressure = 0.0;
    double bulk     = 0.0;
    for (int i = 0; i < max_iters; ++i) {
        law_fn(eta, K0, K0_prime, pressure, bulk);
        if (pressure < pressure_target) { lo = eta; } else { hi = eta; }

        // A step smaller than the spacing of doubles lands on the bracket's edge, which is convergence and not
        // an escape, so the bounds are inclusive.
        double next = eta + (pressure_target - pressure) * eta / bulk;
        if (!(next >= lo && next <= hi)) { next = 0.5 * (lo + hi); }

        // Converged once the compression stops changing to relative tolerance.
        if (std::abs(next - eta) <= rtol * eta) { return next; }
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
// c_MaterialState: the frequency-independent properties of a material at one point
// =====================================================================================================================
// What c_MaterialEOSBase::calc_material_state fills and what the solved EOS reports at a radius. The moduli and
// viscosities are the values after the partial-melt model; a complex modulus is for the rheology to compute from
// them.
struct c_MaterialState {
    double density         = std::numeric_limits<double>::quiet_NaN();  // [kg/m^3]
    double melt_fraction   = 0.0;                                       // [m^3/m^3]
    double shear_modulus   = std::numeric_limits<double>::quiet_NaN();  // [Pa]
    double bulk_modulus    = std::numeric_limits<double>::quiet_NaN();  // [Pa]
    double shear_viscosity = std::numeric_limits<double>::quiet_NaN();  // [Pa s]
    double bulk_viscosity  = std::numeric_limits<double>::quiet_NaN();  // [Pa s]
};

// =====================================================================================================================
// c_MaterialEOSBase: abstract base for all EOS models
// =====================================================================================================================
class c_MaterialEOSBase : public c_PhysicsBase {
public:
    explicit c_MaterialEOSBase(const std::string& model_name) : c_PhysicsBase(model_name) {}
    c_MaterialEOSBase(const std::string& model_name, const c_MaterialEOSConfig& cfg)
        : c_PhysicsBase(model_name),
          p_thermal_expansion(cfg.thermal_expansion),
          p_reference_temperature(cfg.reference_temperature),
          p_shear_modulus_static(cfg.shear_modulus_static),
          p_bulk_modulus_static(cfg.bulk_modulus_static),
          p_shear_viscosity_static(cfg.shear_viscosity_static),
          p_bulk_viscosity_static(cfg.bulk_viscosity_static),
          p_shear_modulus_pressure_derivative(cfg.shear_modulus_pressure_derivative),
          p_shear_modulus_temperature_derivative(cfg.shear_modulus_temperature_derivative),
          p_shear_modulus_reference_temperature(cfg.shear_modulus_reference_temperature),
          p_thermal_conductivity(cfg.thermal_conductivity),
          p_heat_capacity(cfg.heat_capacity) {}
    ~c_MaterialEOSBase() override = default;

    double get_thermal_expansion()     const noexcept { return this->p_thermal_expansion; }
    double get_reference_temperature() const noexcept { return this->p_reference_temperature; }

    // Static constants and the shear law (see c_MaterialEOSConfig).
    double get_shear_modulus_static()   const noexcept { return this->p_shear_modulus_static; }
    double get_bulk_modulus_static()    const noexcept { return this->p_bulk_modulus_static; }
    double get_shear_viscosity_static() const noexcept { return this->p_shear_viscosity_static; }
    double get_bulk_viscosity_static()  const noexcept { return this->p_bulk_viscosity_static; }
    double get_shear_modulus_pressure_derivative() const noexcept {
        return this->p_shear_modulus_pressure_derivative;
    }
    double get_shear_modulus_temperature_derivative() const noexcept {
        return this->p_shear_modulus_temperature_derivative;
    }
    double get_shear_modulus_reference_temperature() const noexcept {
        return this->p_shear_modulus_reference_temperature;
    }
    // Thermal constants. The thermal diffusivity [m^2/s] is k / (rho c_p) at a density [kg/m^3]; NaN when the
    // density or heat capacity is not positive.
    double get_thermal_conductivity() const noexcept { return this->p_thermal_conductivity; }
    double get_heat_capacity()        const noexcept { return this->p_heat_capacity; }
    double calc_thermal_diffusivity(double density) const noexcept {
        if (!(density > 0.0) || !(this->p_heat_capacity > 0.0)) { return TidalPyConstants::d_NAN; }
        return this->p_thermal_conductivity / (density * this->p_heat_capacity);
    }

    void set_shear_modulus_static(double value)   noexcept { this->p_shear_modulus_static = value; }
    void set_bulk_modulus_static(double value)    noexcept { this->p_bulk_modulus_static = value; }
    void set_shear_viscosity_static(double value) noexcept { this->p_shear_viscosity_static = value; }
    void set_bulk_viscosity_static(double value)  noexcept { this->p_bulk_viscosity_static = value; }

    // Viscosity and partial-melt models (ownership transfers in; each is optional).
    void set_shear_viscosity(std::unique_ptr<c_ViscosityBase> model) {
        this->p_shear_viscosity_model = std::move(model);
    }
    void set_bulk_viscosity(std::unique_ptr<c_ViscosityBase> model) {
        this->p_bulk_viscosity_model = std::move(model);
    }
    void set_partial_melt(std::unique_ptr<c_PartialMeltBase> model) {
        this->p_partial_melt_model = std::move(model);
    }
    c_ViscosityBase*   get_shear_viscosity_model() const noexcept { return this->p_shear_viscosity_model.get(); }
    c_ViscosityBase*   get_bulk_viscosity_model()  const noexcept { return this->p_bulk_viscosity_model.get(); }
    c_PartialMeltBase* get_partial_melt_model()    const noexcept { return this->p_partial_melt_model.get(); }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_PhysicsBase::append_config_entries(out);
        out.push_back(c_config_double("thermal_expansion_1_k", this->p_thermal_expansion));
        out.push_back(c_config_double("reference_temperature_k", this->p_reference_temperature));
        out.push_back(c_config_double("shear_modulus_static_pa", this->p_shear_modulus_static));
        out.push_back(c_config_double("bulk_modulus_static_pa", this->p_bulk_modulus_static));
        // An unset (NaN) static viscosity is left out: absence means unset on the way back in, and a config
        // that holds no NaN compares equal to itself.
        if (std::isfinite(this->p_shear_viscosity_static)) {
            out.push_back(c_config_double("shear_viscosity_static_pas", this->p_shear_viscosity_static));
        }
        if (std::isfinite(this->p_bulk_viscosity_static)) {
            out.push_back(c_config_double("bulk_viscosity_static_pas", this->p_bulk_viscosity_static));
        }
        out.push_back(c_config_double("shear_modulus_pressure_derivative", this->p_shear_modulus_pressure_derivative));
        out.push_back(c_config_double(
            "shear_modulus_temperature_derivative_pa_k", this->p_shear_modulus_temperature_derivative));
        out.push_back(c_config_double(
            "shear_modulus_reference_temperature_k", this->p_shear_modulus_reference_temperature));
        out.push_back(c_config_double("thermal_conductivity_w_mk", this->p_thermal_conductivity));
        out.push_back(c_config_double("heat_capacity_j_kgk", this->p_heat_capacity));
    }

    // Density [kg/m^3] from pressure [Pa], temperature [K], and radius [m]; analytic models use the pressure, the
    // interpolated model the radius. A non-finite temperature gives the athermal density.
    virtual double calc_density(
        double pressure, double temperature, double radius) const = 0;

    // Density and isothermal bulk modulus [Pa] together, so a model that inverts its pressure law does it once.
    // The bulk modulus is NaN (the material constant applies) unless the model defines one.
    virtual void calc_density_and_bulk_modulus(
            double pressure,
            double temperature,
            double radius,
            double& density,
            double& bulk_modulus) const {
        density      = this->calc_density(pressure, temperature, radius);
        bulk_modulus = this->p_has_bulk_modulus_table
            ? this->get_tabulated_bulk_modulus(radius) : TidalPyConstants::d_NAN;
    }

    // Isothermal bulk modulus [Pa] at a pressure, temperature, and radius (see calc_density_and_bulk_modulus).
    double calc_bulk_modulus(double pressure, double temperature, double radius) const {
        double density      = TidalPyConstants::d_NAN;
        double bulk_modulus = TidalPyConstants::d_NAN;
        this->calc_density_and_bulk_modulus(pressure, temperature, radius, density, bulk_modulus);
        return bulk_modulus;
    }

    // Radius-varying static moduli [Pa] and viscosities [Pa s] from a model that stores tables of them. NaN means
    // the model holds no such table, and the law or constant of the material applies; only c_InterpolatedEOS
    // has any.
    virtual double get_tabulated_shear_modulus(double /*radius*/) const {
        return std::numeric_limits<double>::quiet_NaN();
    }
    virtual double get_tabulated_bulk_modulus(double /*radius*/) const {
        return std::numeric_limits<double>::quiet_NaN();
    }
    virtual double get_tabulated_shear_viscosity(double /*radius*/) const {
        return std::numeric_limits<double>::quiet_NaN();
    }
    virtual double get_tabulated_bulk_viscosity(double /*radius*/) const {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Every frequency-independent property of the material at a pressure [Pa], temperature [K], and radius [m]:
    // the one place a point is mapped onto the laws and models of the material, called by the EOS solve as it
    // integrates. thermal_density says whether the density law sees the temperature (the viscosity and melt
    // models always do). Each property comes from one source and only that source is evaluated: a table when
    // the model carries one (the p_has_*_table flags), else the shear law, the K of the pressure law, a viscosity
    // model, or the constant. The partial-melt model then weakens the shear pair and the bulk pair.
    void calc_material_state(
            double pressure,
            double temperature,
            bool thermal_density,
            double radius,
            c_MaterialState& out) const {
        /* Density and Bulk Modulus */
        const double density_temperature = thermal_density ? temperature : TidalPyConstants::d_NAN;
        double bulk = TidalPyConstants::d_NAN;
        this->calc_density_and_bulk_modulus(pressure, density_temperature, radius, out.density, bulk);
        if (!std::isfinite(bulk)) { bulk = this->p_bulk_modulus_static; }
        
        /* Shear Modulus */
        double shear = this->p_has_shear_modulus_table
            ? this->get_tabulated_shear_modulus(radius) : TidalPyConstants::d_NAN;
        if (!std::isfinite(shear)) {
            shear = this->p_shear_modulus_static + this->p_shear_modulus_pressure_derivative * pressure;
            if (std::isfinite(temperature)) {
                shear += this->p_shear_modulus_temperature_derivative
                    * (temperature - this->p_shear_modulus_reference_temperature);
            }
        }
        // Floored so a steep temperature derivative cannot drive the shear law negative.
        const double min_modulus = tidalpy_config_ptr->d_MIN_MODULUS;
        if (shear < min_modulus) { shear = min_modulus; }

        /* Viscosities */
        double shear_viscosity = this->p_has_shear_viscosity_table
            ? this->get_tabulated_shear_viscosity(radius) : TidalPyConstants::d_NAN;
        if (!std::isfinite(shear_viscosity)) {
            shear_viscosity = this->p_shear_viscosity_model
                ? this->p_shear_viscosity_model->calc_viscosity(temperature, pressure)
                : this->p_shear_viscosity_static;
        }
        double bulk_viscosity = this->p_has_bulk_viscosity_table
            ? this->get_tabulated_bulk_viscosity(radius) : TidalPyConstants::d_NAN;
        if (!std::isfinite(bulk_viscosity)) {
            bulk_viscosity = this->p_bulk_viscosity_model
                ? this->p_bulk_viscosity_model->calc_viscosity(temperature, pressure)
                : this->p_bulk_viscosity_static;
        }
        
        /* Partial Melting */
        out.melt_fraction = 0.0;
        if (this->p_partial_melt_model) {
            // The liquid viscosity is the pre-melt viscosity until a dedicated liquid-viscosity model exists.
            c_PartialMeltInputs inputs;
            inputs.temperature       = temperature;
            inputs.premelt_viscosity = shear_viscosity;
            inputs.premelt_shear     = shear;
            inputs.liquid_viscosity  = shear_viscosity;
            const c_PartialMeltResult shear_result = this->p_partial_melt_model->calc_partial_melt(inputs);
            inputs.premelt_viscosity = bulk_viscosity;
            inputs.premelt_shear     = bulk;
            inputs.liquid_viscosity  = bulk_viscosity;
            const c_PartialMeltResult bulk_result = this->p_partial_melt_model->calc_partial_melt(inputs);
            out.melt_fraction = shear_result.melt_fraction;
            shear             = shear_result.postmelt_shear_modulus;
            shear_viscosity   = shear_result.postmelt_viscosity;
            bulk              = bulk_result.postmelt_shear_modulus;
            bulk_viscosity    = bulk_result.postmelt_viscosity;
        }
        out.shear_modulus   = shear;
        out.bulk_modulus    = bulk;
        out.shear_viscosity = shear_viscosity;
        out.bulk_viscosity  = bulk_viscosity;
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

    // The material section every model appends after its own binary record: nine doubles (the static constants,
    // the shear law, the conductivity, and the heat capacity), then a presence flag and nested record for each
    // of the three optional models.
    void write_material_binary(std::ostream& out) const {
        const double values[9] = {
            this->p_shear_modulus_static, this->p_bulk_modulus_static,
            this->p_shear_viscosity_static, this->p_bulk_viscosity_static,
            this->p_shear_modulus_pressure_derivative, this->p_shear_modulus_temperature_derivative,
            this->p_shear_modulus_reference_temperature,
            this->p_thermal_conductivity, this->p_heat_capacity};
        out.write(reinterpret_cast<const char*>(values), sizeof(values));
        if (!out) { throw std::runtime_error("TidalPy: failed to write material EOS binary data"); }
        write_optional_binary(out, this->p_shear_viscosity_model);
        write_optional_binary(out, this->p_bulk_viscosity_model);
        write_optional_binary(out, this->p_partial_melt_model);
    }

    void read_material_binary(std::istream& in, bool force) {
        double values[9];
        in.read(reinterpret_cast<char*>(values), sizeof(values));
        if (!in) { throw std::runtime_error("TidalPy: failed to read material EOS binary data"); }
        this->p_shear_modulus_static                 = values[0];
        this->p_bulk_modulus_static                  = values[1];
        this->p_shear_viscosity_static               = values[2];
        this->p_bulk_viscosity_static                = values[3];
        this->p_shear_modulus_pressure_derivative    = values[4];
        this->p_shear_modulus_temperature_derivative = values[5];
        this->p_shear_modulus_reference_temperature  = values[6];
        this->p_thermal_conductivity                 = values[7];
        this->p_heat_capacity                        = values[8];
        this->p_shear_viscosity_model = read_optional_binary<c_ViscosityBase>(in, force, c_viscosity_from_binary);
        this->p_bulk_viscosity_model  = read_optional_binary<c_ViscosityBase>(in, force, c_viscosity_from_binary);
        this->p_partial_melt_model    = read_optional_binary<c_PartialMeltBase>(in, force, c_partial_melt_from_binary);
    }

    double p_thermal_expansion     = 0.0;
    double p_reference_temperature = d_EOS_REFERENCE_TEMPERATURE;

    double p_shear_modulus_static   = 0.0;
    double p_bulk_modulus_static    = 0.0;
    double p_shear_viscosity_static = std::numeric_limits<double>::quiet_NaN();
    double p_bulk_viscosity_static  = std::numeric_limits<double>::quiet_NaN();
    double p_shear_modulus_pressure_derivative    = 0.0;
    double p_shear_modulus_temperature_derivative = 0.0;
    double p_shear_modulus_reference_temperature  = d_EOS_REFERENCE_TEMPERATURE;
    double p_thermal_conductivity = 4.0;
    double p_heat_capacity        = 1200.0;

    std::unique_ptr<c_ViscosityBase>   p_shear_viscosity_model;
    std::unique_ptr<c_ViscosityBase>   p_bulk_viscosity_model;
    std::unique_ptr<c_PartialMeltBase> p_partial_melt_model;

    // Which properties the model tabulates by radius. Set by a model that holds tables (c_InterpolatedEOS), so
    // calc_material_state asks for a table only where there is one and evaluates a law only where there is not.
    bool p_has_shear_modulus_table   = false;
    bool p_has_bulk_modulus_table    = false;
    bool p_has_shear_viscosity_table = false;
    bool p_has_bulk_viscosity_table  = false;
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
        this->write_material_binary(out);
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 3);
        this->p_reference_density     = params[0];
        this->p_thermal_expansion     = params[1];
        this->p_reference_temperature = params[2];
        this->read_material_binary(in, force);
    }

protected:
    double p_reference_density = 3500.0;
};

// 3rd-order Birch-Murnaghan, density from pressure.
class c_BirchMurnaghanEOS : public c_MaterialEOSBase {
public:
    c_BirchMurnaghanEOS() : c_MaterialEOSBase("birch_murnaghan") { this->update_law_range(); }
    explicit c_BirchMurnaghanEOS(const c_MaterialEOSConfig& cfg)
        : c_MaterialEOSBase("birch_murnaghan", cfg),
          p_reference_density(cfg.reference_density),
          p_reference_bulk_modulus(cfg.reference_bulk_modulus),
          p_bulk_modulus_derivative(cfg.bulk_modulus_derivative),
          p_invert_rtol(cfg.invert_rtol),
          p_invert_max_iters(cfg.invert_max_iters) { this->update_law_range(); }
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
        this->write_material_binary(out);
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
        this->read_material_binary(in, force);
        this->update_law_range();
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
            eos_bm_pressure_and_bulk_modulus,
            this->p_law_range,
            this->p_invert_rtol,
            this->p_invert_max_iters);
    }

    // The law's monotonic range follows from K0 and K0' alone; found again whenever they change.
    void update_law_range() noexcept {
        this->p_law_range = eos_find_monotonic_range(
            this->p_reference_bulk_modulus,
            this->p_bulk_modulus_derivative,
            eos_bm_pressure_and_bulk_modulus,
            this->p_invert_rtol);
    }

    c_PressureLawRange p_law_range;
    double p_reference_density       = 3500.0;
    double p_reference_bulk_modulus  = 1.0e11;
    double p_bulk_modulus_derivative = 4.0;
    double p_invert_rtol             = d_EOS_INVERT_RTOL;
    int    p_invert_max_iters        = d_EOS_INVERT_MAX_ITERS;
};

// Vinet (universal) EOS, density from pressure.
class c_VinetEOS : public c_MaterialEOSBase {
public:
    c_VinetEOS() : c_MaterialEOSBase("vinet") { this->update_law_range(); }
    explicit c_VinetEOS(const c_MaterialEOSConfig& cfg)
        : c_MaterialEOSBase("vinet", cfg),
          p_reference_density(cfg.reference_density),
          p_reference_bulk_modulus(cfg.reference_bulk_modulus),
          p_bulk_modulus_derivative(cfg.bulk_modulus_derivative),
          p_invert_rtol(cfg.invert_rtol),
          p_invert_max_iters(cfg.invert_max_iters) { this->update_law_range(); }
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
        this->write_material_binary(out);
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
        this->read_material_binary(in, force);
        this->update_law_range();
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
            eos_vinet_pressure_and_bulk_modulus,
            this->p_law_range,
            this->p_invert_rtol,
            this->p_invert_max_iters);
    }

    // The law's monotonic range follows from K0 and K0' alone; found again whenever they change.
    void update_law_range() noexcept {
        this->p_law_range = eos_find_monotonic_range(
            this->p_reference_bulk_modulus,
            this->p_bulk_modulus_derivative,
            eos_vinet_pressure_and_bulk_modulus,
            this->p_invert_rtol);
    }

    c_PressureLawRange p_law_range;
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
        this->p_update_table_flags();
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
    double get_tabulated_shear_modulus(double radius) const override {
        return this->p_interp_optional(radius, this->p_shear_modulus);
    }
    double get_tabulated_bulk_modulus(double radius) const override {
        return this->p_interp_optional(radius, this->p_bulk_modulus);
    }
    double get_tabulated_shear_viscosity(double radius) const override {
        return this->p_interp_optional(radius, this->p_shear_viscosity);
    }
    double get_tabulated_bulk_viscosity(double radius) const override {
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
        this->write_material_binary(out);
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
        this->p_update_table_flags();
        in.read(reinterpret_cast<char*>(&this->p_thermal_expansion),     sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_reference_temperature), sizeof(double));
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read interpolated EOS binary data");
        }
        this->read_material_binary(in, force);
    }

protected:
    void p_update_table_flags() noexcept {
        this->p_has_shear_modulus_table   = this->has_shear_modulus();
        this->p_has_bulk_modulus_table    = this->has_bulk_modulus();
        this->p_has_shear_viscosity_table = this->has_shear_viscosity();
        this->p_has_bulk_viscosity_table  = this->has_bulk_viscosity();
    }

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
