#pragma once
/* Equation-of-state laws: a phase's density [kg m-3], isothermal and adiabatic bulk moduli [Pa], and thermal
 * expansivity [1/K] at a point (pressure, temperature, radius). All MKS.
 *
 * Every law carries the same thermal parameters, so density, thermal pressure, the adiabat, and the Rayleigh number
 * all use one expansivity:
 *   - alpha0 at the reference temperature T_ref, scaled with compression by an Anderson-Gruneisen parameter that
 *     itself falls with compression, delta_T = delta_T0 (rho0 / rho)^kappa, so
 *     alpha = alpha0 exp[(delta_T0 / kappa) ((rho0 / rho)^kappa - 1)] (kappa -> 0 gives (rho0 / rho)^delta_T0);
 *   - a Gruneisen parameter gamma, which gives the adiabatic bulk modulus K_S = K_T (1 + alpha gamma T), the one a
 *     tidal (adiabatic) deformation sees; gamma = 0 makes K_S equal K_T.
 * The pressure laws (Birch-Murnaghan, Vinet, Murnaghan) add a thermal pressure alpha0 K0 (T - T_ref), taking
 * alpha K_T constant at its high-temperature limit; the constant and tabulated laws scale their density by
 * exp(-alpha0 (T - T_ref)). A law evaluated with `thermal` off, or at a non-finite temperature, gives the athermal
 * density; its expansivity is still reported, since the adiabat needs it either way.
 *
 * References
 * ----------
 * - Birch (1947), Phys. Rev. 71, 809; Vinet et al. (1987), J. Geophys. Res. 92, 9319; Murnaghan (1944), PNAS 30,
 *   244.
 * - Anderson (1995), Equations of State of Solids for Geophysics and Ceramic Science (thermal pressure, K_S / K_T).
 * - Anderson (1967); Chopelas and Boehler (1992), GRL 19, 1983 (the Anderson-Gruneisen compression scaling).
 * - Chandrasekhar (1939), An Introduction to the Study of Stellar Structure (polytropes).
 * - Seager et al. (2007), ApJ 669, 1279 (the modified polytrope rho = rho0 + c P^n).
 */

#include <cmath>
#include <istream>
#include <memory>
#include <string>
#include <vector>

#include "broadcast_.hpp"
#include "physics_base_.hpp"
#include "registry_.hpp"
#include "spec_model_.hpp"
#include "thermo_point_.hpp"
#include "pressure_laws_.hpp"
#include "../../constants_.hpp"
#include "../../Utilities/arrays/table_lookup_.hpp"
#include "../../Utilities/math/numerics_.hpp"  // c_safe_pow, c_safe_exp

namespace tidalpy {

// What an equation-of-state law gives at a point.
struct c_EOSPoint {
    double density                = TidalPyConstants::d_NAN;   // [kg m-3]
    double bulk_modulus           = TidalPyConstants::d_NAN;   // isothermal, K_T [Pa]; NaN when the law gives none
    double adiabatic_bulk_modulus = TidalPyConstants::d_NAN;   // K_S [Pa]
    double thermal_expansion      = TidalPyConstants::d_NAN;   // alpha [1/K]
};

class c_EOSBase : public c_PhysicsBase {
public:
    static constexpr const char* C_FAMILY_NAME = "equation of state";

    explicit c_EOSBase(const std::string& model_name) : c_PhysicsBase(model_name) {}
    ~c_EOSBase() override = default;

    // The law at a point. `thermal` says whether the density sees the temperature (a layer's use_thermal_expansion).
    void calc_eos(const c_ThermoPoint& point, bool thermal, c_EOSPoint& out) const noexcept {
        this->p_calc_law(point, this->p_temperature_offset(point.temperature, thermal), out);
        out.thermal_expansion = this->calc_thermal_expansion(out.density);
        const bool adiabatic = (this->p_gruneisen_parameter != 0.0) && std::isfinite(point.temperature)
            && std::isfinite(out.thermal_expansion);
        out.adiabatic_bulk_modulus = adiabatic
            ? out.bulk_modulus * (1.0 + out.thermal_expansion * this->p_gruneisen_parameter * point.temperature)
            : out.bulk_modulus;
    }

    // The density alone [kg m-3], what the structure integration asks for.
    double calc_density(const c_ThermoPoint& point, bool thermal) const noexcept {
        c_EOSPoint law_point;
        this->p_calc_law(point, this->p_temperature_offset(point.temperature, thermal), law_point);
        return law_point.density;
    }

    // Thermal expansivity [1/K] at a density [kg m-3]: alpha0 times the Anderson-Gruneisen factor. alpha0 for a law
    // with no reference density, a zero delta_T0, or a density that is not positive.
    double calc_thermal_expansion(double density) const noexcept {
        const double reference_density = this->p_expansion_reference_density();
        const double delta_t0 = this->p_anderson_gruneisen_parameter;
        if ((delta_t0 == 0.0) || !(reference_density > 0.0) || !(density > 0.0)) {
            return this->p_thermal_expansion;
        }
        const double expansion = reference_density / density;
        const double kappa = this->p_anderson_gruneisen_exponent;
        if (kappa == 0.0) { return this->p_thermal_expansion * std::pow(expansion, delta_t0); }
        return this->p_thermal_expansion * std::exp((delta_t0 / kappa) * (std::pow(expansion, kappa) - 1.0));
    }

    // The pressures a pressure law represents (the cold pressure, less any thermal pressure); unbounded for a law
    // with no such limit.
    virtual c_PressureLawRange get_pressure_law_range() const noexcept { return c_PressureLawRange{}; }

    double get_reference_temperature() const noexcept { return this->p_reference_temperature; }

    // Element-wise over pressure, temperature, and radius; each holds one value per point or a single value used at
    // every point. Fills the four outputs.
    void calc_eos_vectorize(
            const std::vector<double>& pressure,
            const std::vector<double>& temperature,
            const std::vector<double>& radius,
            bool thermal,
            std::vector<double>& out_density,
            std::vector<double>& out_bulk_modulus,
            std::vector<double>& out_adiabatic_bulk_modulus,
            std::vector<double>& out_thermal_expansion) const {
        const std::size_t num_points = c_broadcast_length(
            {pressure.size(), temperature.size(), radius.size()}, "calc_eos_vectorize");
        const std::size_t pressure_stride    = c_broadcast_stride(pressure.size());
        const std::size_t temperature_stride = c_broadcast_stride(temperature.size());
        const std::size_t radius_stride      = c_broadcast_stride(radius.size());
        out_density.resize(num_points);
        out_bulk_modulus.resize(num_points);
        out_adiabatic_bulk_modulus.resize(num_points);
        out_thermal_expansion.resize(num_points);
        c_ThermoPoint point;
        c_EOSPoint law_point;
        for (std::size_t i = 0; i < num_points; ++i) {
            point.pressure    = pressure[i * pressure_stride];
            point.temperature = temperature[i * temperature_stride];
            point.radius      = radius[i * radius_stride];
            this->calc_eos(point, thermal, law_point);
            out_density[i]                = law_point.density;
            out_bulk_modulus[i]           = law_point.bulk_modulus;
            out_adiabatic_bulk_modulus[i] = law_point.adiabatic_bulk_modulus;
            out_thermal_expansion[i]      = law_point.thermal_expansion;
        }
    }

protected:
    // The law itself: density and isothermal bulk modulus at a point, with the temperature above the reference
    // (zero for the athermal law).
    virtual void p_calc_law(
        const c_ThermoPoint& point, double temperature_offset, c_EOSPoint& out) const noexcept = 0;

    // The density the expansivity's compression scaling is measured from; NaN for a law without one.
    virtual double p_expansion_reference_density() const noexcept { return TidalPyConstants::d_NAN; }

    double p_temperature_offset(double temperature, bool thermal) const noexcept {
        if (!thermal || (this->p_thermal_expansion == 0.0) || !std::isfinite(temperature)) { return 0.0; }
        return temperature - this->p_reference_temperature;
    }

    // The thermal parameters every law shares, as rows of a model's parameter table.
    template <class Model>
    static void p_append_thermal_specs(std::vector<c_ParamSpec<Model>>& specs) {
        specs.push_back({"thermal_expansion", "thermal_expansion_1_k", &Model::p_thermal_expansion, 0.0,
                         c_ParamBounds::Finite, "Thermal expansivity alpha0 at the reference state [1/K]."});
        specs.push_back({"reference_temperature", "reference_temperature_k", &Model::p_reference_temperature,
                         d_EOS_REFERENCE_TEMPERATURE, c_ParamBounds::Positive,
                         "Temperature where the reference density and bulk modulus apply [K]."});
        specs.push_back({"anderson_gruneisen_parameter", "anderson_gruneisen_parameter",
                         &Model::p_anderson_gruneisen_parameter, 0.0, c_ParamBounds::Finite,
                         "Anderson-Gruneisen parameter delta_T0 at the reference density; 0 keeps alpha constant."});
        specs.push_back({"anderson_gruneisen_exponent", "anderson_gruneisen_exponent",
                         &Model::p_anderson_gruneisen_exponent, 0.0, c_ParamBounds::Finite,
                         "Compression exponent kappa of delta_T = delta_T0 (rho0 / rho)^kappa."});
        specs.push_back({"gruneisen_parameter", "gruneisen_parameter", &Model::p_gruneisen_parameter, 0.0,
                         c_ParamBounds::NonNegative,
                         "Gruneisen parameter gamma in K_S = K_T (1 + alpha gamma T); 0 makes K_S equal K_T."});
    }

    double p_thermal_expansion            = 0.0;
    double p_reference_temperature        = d_EOS_REFERENCE_TEMPERATURE;
    double p_anderson_gruneisen_parameter = 0.0;
    double p_anderson_gruneisen_exponent  = 0.0;
    double p_gruneisen_parameter          = 0.0;
};

// Uniform density (aliases "uniform", "constant_density"): rho0 exp(-alpha0 (T - T_ref)), with a constant bulk
// modulus for the radial solver.
class c_ConstantEOS final : public c_SpecModel<c_ConstantEOS, c_EOSBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::ConstantEOSLaw;

    static const std::vector<c_ParamSpec<c_ConstantEOS>>& parameter_specs() {
        using Self = c_ConstantEOS;
        static const std::vector<c_ParamSpec<Self>> specs = [] {
            std::vector<c_ParamSpec<Self>> rows = {
                {"reference_density", "reference_density_kg_m3", &Self::p_reference_density, 3500.0,
                 c_ParamBounds::Positive, "Density at the reference temperature [kg m-3]."},
                {"bulk_modulus", "bulk_modulus_pa", &Self::p_bulk_modulus, 1.0e11, c_ParamBounds::NonNegative,
                 "Bulk modulus [Pa]."},
            };
            p_append_thermal_specs(rows);
            return rows;
        }();
        return specs;
    }

    c_ConstantEOS() : c_ConstantEOS(c_ParamMap{}) {}
    explicit c_ConstantEOS(const c_ParamMap& params) : c_SpecModel("constant") { this->p_initialize(params); }

protected:
    void p_calc_law(const c_ThermoPoint& /*point*/, double temperature_offset, c_EOSPoint& out) const noexcept override {
        out.density      = this->p_reference_density * c_safe_exp(-this->p_thermal_expansion * temperature_offset);
        out.bulk_modulus = this->p_bulk_modulus;
    }
    double p_expansion_reference_density() const noexcept override { return this->p_reference_density; }

    double p_reference_density = 0.0;
    double p_bulk_modulus      = 0.0;

};

// Density from pressure through an analytic pressure law (Birch-Murnaghan or Vinet), inverted at the pressure less
// the thermal pressure alpha0 K0 (T - T_ref). Beyond the range the law rises over, the density holds at that end.
template <class Derived, c_PressureLaw LAW>
class c_PressureLawEOSModel : public c_SpecModel<Derived, c_EOSBase> {
public:
    using c_SpecModel<Derived, c_EOSBase>::c_SpecModel;

    c_PressureLawRange get_pressure_law_range() const noexcept override { return this->p_law_range; }

    static const std::vector<c_ParamSpec<Derived>>& parameter_specs() {
        static const std::vector<c_ParamSpec<Derived>> specs = [] {
            std::vector<c_ParamSpec<Derived>> rows = {
                {"reference_density", "reference_density_kg_m3", &Derived::p_reference_density, 3500.0,
                 c_ParamBounds::Positive, "Density at zero pressure and the reference temperature, rho0 [kg m-3]."},
                {"reference_bulk_modulus", "reference_bulk_modulus_pa", &Derived::p_reference_bulk_modulus, 1.0e11,
                 c_ParamBounds::Positive, "Isothermal bulk modulus at the reference state, K0 [Pa]."},
                {"bulk_modulus_derivative", "bulk_modulus_derivative", &Derived::p_bulk_modulus_derivative, 4.0,
                 c_ParamBounds::Finite, "Pressure derivative of the bulk modulus, K0' [dimensionless]."},
                {"invert_rtol", "invert_rtol", &Derived::p_invert_rtol, TidalPyConstants::d_NAN, c_ParamBounds::Any,
                 "Relative tolerance of the density inversion; unset takes [numerical] eos_invert_rtol."},
                {"invert_max_iters", "invert_max_iters", &Derived::p_invert_max_iters, -1.0, c_ParamBounds::Any,
                 "Iteration cap of the density inversion; unset (-1) takes [numerical] eos_invert_max_iters."},
            };
            c_EOSBase::p_append_thermal_specs(rows);
            return rows;
        }();
        return specs;
    }

protected:
    void p_calc_law(const c_ThermoPoint& point, double temperature_offset, c_EOSPoint& out) const noexcept override {
        const double thermal_pressure = this->p_thermal_expansion * this->p_reference_bulk_modulus * temperature_offset;
        const double eta = eos_invert_eta(
            point.pressure - thermal_pressure, this->p_reference_bulk_modulus, this->p_bulk_modulus_derivative,
            p_law_function(), this->p_law_range, this->p_resolved_rtol, this->p_resolved_max_iters);
        double pressure = 0.0;
        double bulk_modulus = 0.0;
        p_law_function()(eta, this->p_reference_bulk_modulus, this->p_bulk_modulus_derivative, pressure, bulk_modulus);
        out.density      = this->p_reference_density * eta;
        out.bulk_modulus = bulk_modulus;
    }

    double p_expansion_reference_density() const noexcept override { return this->p_reference_density; }

    // The inversion settings and the law's range depend on the parameters alone, so they are found once.
    void p_update_derived() noexcept override {
        this->p_resolved_rtol      = c_resolve_eos_invert_rtol(this->p_invert_rtol);
        this->p_resolved_max_iters = c_resolve_eos_invert_max_iters(this->p_invert_max_iters);
        this->p_law_range = eos_find_monotonic_range(
            this->p_reference_bulk_modulus, this->p_bulk_modulus_derivative, p_law_function(), this->p_resolved_rtol);
    }

    static auto p_law_function() noexcept {
        if constexpr (LAW == c_PressureLaw::BirchMurnaghan) { return &eos_bm_pressure_and_bulk_modulus; }
        else { return &eos_vinet_pressure_and_bulk_modulus; }
    }

    double p_reference_density      = 0.0;
    double p_reference_bulk_modulus = 0.0;
    double p_bulk_modulus_derivative = 0.0;
    double p_invert_rtol            = TidalPyConstants::d_NAN;
    int    p_invert_max_iters       = -1;
    // Rebuilt from the parameters after every change; never serialized.
    double p_resolved_rtol          = 0.0;
    int    p_resolved_max_iters     = 0;
    c_PressureLawRange p_law_range;
};

// 3rd-order Birch-Murnaghan (aliases "bm", "birch-murnaghan").
class c_BirchMurnaghanEOS final : public c_PressureLawEOSModel<c_BirchMurnaghanEOS, c_PressureLaw::BirchMurnaghan> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::BirchMurnaghanEOSLaw;
    c_BirchMurnaghanEOS() : c_BirchMurnaghanEOS(c_ParamMap{}) {}
    explicit c_BirchMurnaghanEOS(const c_ParamMap& params) : c_PressureLawEOSModel("birch_murnaghan") {
        this->p_initialize(params);
    }
};

// Vinet (universal) law.
class c_VinetEOS final : public c_PressureLawEOSModel<c_VinetEOS, c_PressureLaw::Vinet> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::VinetEOSLaw;
    c_VinetEOS() : c_VinetEOS(c_ParamMap{}) {}
    explicit c_VinetEOS(const c_ParamMap& params) : c_PressureLawEOSModel("vinet") { this->p_initialize(params); }
};

// Murnaghan (1944) law, K = K0 + K0' P, closed form both ways: rho = rho0 (1 + K0' P / K0)^(1 / K0'), and
// rho0 exp(P / K0) for K0' = 0 or in tension, which joins it smoothly at P = 0 (the Murnaghan density has no real
// value past P = -K0 / K0'). Its rho / (d rho / dP) is its bulk modulus, so a layer of it is neutrally stratified
// under the bulk modulus the tidal equations see. Suited to melts and liquids over modest pressures.
class c_MurnaghanEOS final : public c_SpecModel<c_MurnaghanEOS, c_EOSBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::MurnaghanEOSLaw;

    static const std::vector<c_ParamSpec<c_MurnaghanEOS>>& parameter_specs() {
        using Self = c_MurnaghanEOS;
        static const std::vector<c_ParamSpec<Self>> specs = [] {
            std::vector<c_ParamSpec<Self>> rows = {
                {"reference_density", "reference_density_kg_m3", &Self::p_reference_density, 2750.0,
                 c_ParamBounds::Positive, "Density at zero pressure and the reference temperature, rho0 [kg m-3]."},
                {"reference_bulk_modulus", "reference_bulk_modulus_pa", &Self::p_reference_bulk_modulus, 2.0e10,
                 c_ParamBounds::Positive, "Isothermal bulk modulus at zero pressure, K0 [Pa]."},
                {"bulk_modulus_derivative", "bulk_modulus_derivative", &Self::p_bulk_modulus_derivative, 5.0,
                 c_ParamBounds::NonNegative, "Pressure derivative of the bulk modulus, K0' [dimensionless]."},
            };
            p_append_thermal_specs(rows);
            return rows;
        }();
        return specs;
    }

    c_MurnaghanEOS() : c_MurnaghanEOS(c_ParamMap{}) {}
    explicit c_MurnaghanEOS(const c_ParamMap& params) : c_SpecModel("murnaghan") { this->p_initialize(params); }

protected:
    void p_calc_law(const c_ThermoPoint& point, double temperature_offset, c_EOSPoint& out) const noexcept override {
        const double k0 = this->p_reference_bulk_modulus;
        const double kp = this->p_bulk_modulus_derivative;
        const double pressure = point.pressure - this->p_thermal_expansion * k0 * temperature_offset;
        if ((kp == 0.0) || (pressure <= 0.0)) {
            out.density      = this->p_reference_density * c_safe_exp(pressure / k0);
            out.bulk_modulus = k0;
            return;
        }
        out.density      = this->p_reference_density * std::pow(1.0 + kp * pressure / k0, 1.0 / kp);
        out.bulk_modulus = k0 + kp * pressure;
    }
    double p_expansion_reference_density() const noexcept override { return this->p_reference_density; }

    double p_reference_density       = 0.0;
    double p_reference_bulk_modulus  = 0.0;
    double p_bulk_modulus_derivative = 0.0;

};

// Polytrope, P = K rho^(1 + 1/n): rho = (P / K)^(n / (n + 1)) and K_T = (1 + 1/n) P. A polytrope is a barotrope, so it
// takes no temperature: the thermal parameters are carried for the adiabat only. The density is zero where the
// pressure is not positive, the surface of a gas envelope.
class c_PolytropeEOS final : public c_SpecModel<c_PolytropeEOS, c_EOSBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::PolytropeEOSLaw;

    static const std::vector<c_ParamSpec<c_PolytropeEOS>>& parameter_specs() {
        using Self = c_PolytropeEOS;
        static const std::vector<c_ParamSpec<Self>> specs = [] {
            std::vector<c_ParamSpec<Self>> rows = {
                {"polytropic_constant", "polytropic_constant", &Self::p_polytropic_constant, 2.0e5,
                 c_ParamBounds::Positive,
                 "K in P = K rho^(1 + 1/n) [Pa (m3 kg-1)^(1 + 1/n)]; about 2e5 for Jupiter at n = 1."},
                {"polytropic_index", "polytropic_index", &Self::p_polytropic_index, 1.0, c_ParamBounds::Positive,
                 "Polytropic index n [dimensionless]."},
            };
            p_append_thermal_specs(rows);
            return rows;
        }();
        return specs;
    }

    c_PolytropeEOS() : c_PolytropeEOS(c_ParamMap{}) {}
    explicit c_PolytropeEOS(const c_ParamMap& params) : c_SpecModel("polytrope") { this->p_initialize(params); }

protected:
    void p_calc_law(const c_ThermoPoint& point, double /*temperature_offset*/, c_EOSPoint& out) const noexcept override {
        if (!(point.pressure > 0.0)) {
            out.density      = 0.0;
            out.bulk_modulus = 0.0;
            return;
        }
        const double n = this->p_polytropic_index;
        out.density      = std::pow(point.pressure / this->p_polytropic_constant, n / (n + 1.0));
        out.bulk_modulus = (1.0 + 1.0 / n) * point.pressure;
    }

    double p_polytropic_constant = 0.0;
    double p_polytropic_index    = 0.0;

};

// Modified polytrope (Seager et al. 2007), rho = rho0 + c P^n for P > 0 and rho0 otherwise, a fit to the cold
// compression of planetary materials to TPa pressures. K_T = rho / (d rho / dP) = rho / (c n P^(n - 1)). Isothermal:
// the thermal parameters are carried for the adiabat only.
class c_ModifiedPolytropeEOS final : public c_SpecModel<c_ModifiedPolytropeEOS, c_EOSBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::ModifiedPolytropeEOSLaw;

    static const std::vector<c_ParamSpec<c_ModifiedPolytropeEOS>>& parameter_specs() {
        using Self = c_ModifiedPolytropeEOS;
        static const std::vector<c_ParamSpec<Self>> specs = [] {
            std::vector<c_ParamSpec<Self>> rows = {
                {"reference_density", "reference_density_kg_m3", &Self::p_reference_density, 8300.0,
                 c_ParamBounds::Positive, "Density at zero pressure, rho0 [kg m-3]."},
                {"polytrope_coefficient", "polytrope_coefficient", &Self::p_coefficient, 0.00349,
                 c_ParamBounds::Positive, "c in rho = rho0 + c P^n [kg m-3 Pa^-n]."},
                {"polytrope_exponent", "polytrope_exponent", &Self::p_exponent, 0.528, c_ParamBounds::Positive,
                 "n in rho = rho0 + c P^n [dimensionless]."},
            };
            p_append_thermal_specs(rows);
            return rows;
        }();
        return specs;
    }

    c_ModifiedPolytropeEOS() : c_ModifiedPolytropeEOS(c_ParamMap{}) {}
    explicit c_ModifiedPolytropeEOS(const c_ParamMap& params) : c_SpecModel("modified_polytrope") {
        this->p_initialize(params);
    }

protected:
    void p_calc_law(const c_ThermoPoint& point, double /*temperature_offset*/, c_EOSPoint& out) const noexcept override {
        if (!(point.pressure > 0.0)) {
            out.density      = this->p_reference_density;
            out.bulk_modulus = (this->p_exponent < 1.0) ? 0.0 : TidalPyConstants::d_INF;
            return;
        }
        const double compression_term = this->p_coefficient * std::pow(point.pressure, this->p_exponent);
        out.density      = this->p_reference_density + compression_term;
        out.bulk_modulus = out.density * point.pressure / (this->p_exponent * compression_term);
    }
    double p_expansion_reference_density() const noexcept override { return this->p_reference_density; }

    double p_reference_density = 0.0;
    double p_coefficient       = 0.0;
    double p_exponent          = 0.0;

};

// A density profile tabulated in radius (aliases "interp", "interpolated"), with an optional bulk-modulus table: a
// seismic model such as PREM. Linear in radius, held at the end values beyond the table, and scaled by
// exp(-alpha0 (T - T_ref)). Its layer must keep its volume, since the table is in radius.
class c_InterpolatedEOS final : public c_SpecModel<c_InterpolatedEOS, c_EOSBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::InterpolatedEOSLaw;

    static const std::vector<c_ParamSpec<c_InterpolatedEOS>>& parameter_specs() {
        using Self = c_InterpolatedEOS;
        static const std::vector<c_ParamSpec<Self>> specs = [] {
            std::vector<c_ParamSpec<Self>> rows = {
                // The default tables hold one point, a constant profile.
                {"radius", "radius_m", &Self::p_radius, 0.0, c_ParamBounds::Finite, "Table radii, ascending [m].",
                 {0.0}},
                {"density", "density_kg_m3", &Self::p_density, 0.0, c_ParamBounds::Positive,
                 "Density at each table radius [kg m-3].", {3500.0}},
                {"bulk_modulus", "bulk_modulus_pa", &Self::p_bulk_modulus, 0.0, c_ParamBounds::NonNegative,
                 "Bulk modulus at each table radius [Pa]; optional."},
            };
            p_append_thermal_specs(rows);
            return rows;
        }();
        return specs;
    }

    c_InterpolatedEOS() : c_InterpolatedEOS(c_ParamMap{}) {}
    explicit c_InterpolatedEOS(const c_ParamMap& params) : c_SpecModel("interpolate") { this->p_initialize(params); }

protected:
    void p_calc_law(const c_ThermoPoint& point, double temperature_offset, c_EOSPoint& out) const noexcept override {
        out.density = c_safe_exp(-this->p_thermal_expansion * temperature_offset)
            * this->p_lookup.interpolate(point.radius, this->p_radius, this->p_density);
        out.bulk_modulus = this->p_lookup.interpolate(point.radius, this->p_radius, this->p_bulk_modulus);
    }

    void p_validate() const override {
        c_check_table(this->p_describe(), this->p_radius, {&this->p_density, &this->p_bulk_modulus});
        if (this->p_density.size() != this->p_radius.size()) {
            throw std::invalid_argument(this->p_describe() + " needs a 'density_kg_m3' value at every table radius.");
        }
    }
    void p_update_derived() noexcept override { this->p_lookup.build(this->p_radius); }

    std::vector<double> p_radius;
    std::vector<double> p_density;
    std::vector<double> p_bulk_modulus;
    c_TableLookup p_lookup;

};

inline const c_ModelRegistry<c_EOSBase>& c_eos_registry() {
    static const c_ModelRegistry<c_EOSBase> registry = {
        {{"constant", "uniform", "constant_density"}, BinaryClassID::ConstantEOSLaw,
         &c_make_entry<c_EOSBase, c_ConstantEOS>},
        {{"birch_murnaghan", "bm", "birch-murnaghan"}, BinaryClassID::BirchMurnaghanEOSLaw,
         &c_make_entry<c_EOSBase, c_BirchMurnaghanEOS>},
        {{"vinet"}, BinaryClassID::VinetEOSLaw, &c_make_entry<c_EOSBase, c_VinetEOS>},
        {{"murnaghan"}, BinaryClassID::MurnaghanEOSLaw, &c_make_entry<c_EOSBase, c_MurnaghanEOS>},
        {{"polytrope"}, BinaryClassID::PolytropeEOSLaw, &c_make_entry<c_EOSBase, c_PolytropeEOS>},
        {{"modified_polytrope", "seager"}, BinaryClassID::ModifiedPolytropeEOSLaw,
         &c_make_entry<c_EOSBase, c_ModifiedPolytropeEOS>},
        {{"interpolate", "interp", "interpolated"}, BinaryClassID::InterpolatedEOSLaw,
         &c_make_entry<c_EOSBase, c_InterpolatedEOS>},
    };
    return registry;
}

inline std::unique_ptr<c_EOSBase> c_find_eos(const std::string& model_name, const c_ParamMap& params) {
    return c_make_model(c_eos_registry(), model_name, params);
}

inline std::unique_ptr<c_EOSBase> c_eos_from_binary(std::istream& in, bool force = false) {
    return c_model_from_binary(c_eos_registry(), in, force);
}

inline std::string c_eos_canonical_name(const std::string& model_name) {
    return c_canonical_model_name(c_eos_registry(), model_name);
}

inline std::vector<std::string> c_eos_model_names() {
    return c_model_names(c_eos_registry());
}

}  // namespace tidalpy
