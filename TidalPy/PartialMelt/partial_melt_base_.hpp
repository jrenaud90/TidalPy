#pragma once
/* Abstract base for TidalPy partial-melt models. Concrete models live in partial_melt_.hpp. All MKS.
 *
 * The melt fraction, the liquid limits, the melt phase's pressure law, and the optional bulk and density effects of
 * melt are model independent, so they live here rather than on the models. These quantities are frequency
 * independent: the world pipeline caches them once after the EOS solve and only repeats the downstream rheology
 * step per forcing frequency.
 *
 * References
 * ----------
 * - Fischer and Spohn (1990): temperature-based melt viscosity and shear law.
 * - Henning (2009, 2010); Renaud and Henning (2018): three-regime melt weakening.
 * - Hashin and Shtrikman (1963), J. Mech. Phys. Solids 11, 127: bounds on the moduli of a two-phase aggregate,
 *   used for the optional bulk-modulus weakening.
 * - Mavko (1980), JGR 85, 5173; Takei (2002), JGR 107, 2043: melt-bearing rock moduli, which bear on the bulk
 *   modulus far less than on the shear modulus.
 * - Murnaghan (1944), PNAS 30, 244: the melt phase's pressure law.
 * - McKenzie (1984), J. Petrol. 25, 713; Takei and Holtzman (2009), JGR 114, B06205: the compaction (bulk)
 *   viscosity of a partially molten rock, eta / phi and of order eta respectively.
 * - Kervazo et al. (2021), A&A 650, A72: bulk dissipation in Io's partially molten interior.
 */

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "physics_base_.hpp"
#include "../constants_.hpp"

namespace tidalpy {

// Combined construction parameters; each model reads only the fields it needs. Its defaults are the models' defaults.
struct c_PartialMeltConfig {
    // Shared melt envelope and liquid limits.
    double solidus             = 1600.0;  // [K]
    double liquidus            = 2000.0;  // [K]
    double liquid_shear        = 1.0e-5;  // shear modulus of the fully molten material [Pa]
    double liquid_viscosity    = 0.2;     // viscosity of the fully molten material [Pa·s]
    bool   bulk_melt_weakening = false;   // weaken the bulk modulus with melt (calc_bulk_modulus_melt)
    double liquid_bulk_modulus = 2.0e10;  // bulk modulus of the melt at zero pressure [Pa]

    // The melt phase's Murnaghan law and the switches for melt's other effects; every switch is off by default.
    double liquid_bulk_modulus_derivative  = 5.0;     // dK_l / dP of the melt [dimensionless]
    double liquid_density                  = 2750.0;  // density of the melt at zero pressure [kg m-3]
    bool   density_melt_mixing             = false;   // mix the melt into the density (calc_mixture_density)
    bool   bulk_viscosity_melt_weakening   = false;   // melt sets a bulk viscosity (calc_bulk_viscosity_melt)
    double melt_bulk_viscosity_coefficient = 1.0;     // c in zeta_melt = c eta / phi^n [dimensionless]
    double melt_bulk_viscosity_exponent    = 1.0;     // n in zeta_melt = c eta / phi^n [dimensionless]

    // Spohn (Fischer & Spohn 1990) parameters.
    double fs_visc_power_slope       = 27000.0;  // [K]
    double fs_visc_log10_at_solidus  = 15.875;   // log10 of the post-melt viscosity at the solidus [Pa s]
    double fs_shear_power_slope      = 82000.0;  // [K]
    double fs_shear_log10_at_solidus = 10.65;    // log10 of the post-melt shear modulus at the solidus [Pa]

    // Henning (2009/2010) parameters.
    double crit_melt_frac         = 0.5;      // [m^3/m^3]
    double crit_melt_frac_width   = 0.05;     // [m^3/m^3]
    double hn_visc_slope_1        = 13.5;
    double hn_visc_falloff_slope  = 370.0;
    double hn_shear_param_1       = 40000.0;  // [K], b1 in exp[b1 (1 / T - 1 / T_sol)]
    double hn_shear_falloff_slope = 700.0;
};

// Per-evaluation state. Material constants live on the model object; only what varies is passed here.
struct c_PartialMeltInputs {
    double temperature       = 0.0;   // local temperature [K]
    double premelt_viscosity = 0.0;   // solid (pre-melt) viscosity [Pa·s]
    double premelt_shear     = 0.0;   // solid (pre-melt) shear modulus [Pa]
};

struct c_PartialMeltResult {
    double melt_fraction          = 0.0;   // volumetric melt fraction φ [m^3/m^3]
    double postmelt_viscosity     = 0.0;   // post-melt viscosity [Pa·s]
    double postmelt_shear_modulus = 0.0;   // post-melt shear modulus [Pa]
};

class c_PartialMeltBase : public c_PhysicsBase {
public:
    // Number of shared parameters every model writes first in its binary payload.
    static constexpr std::size_t C_NUM_ENVELOPE_PARAMS = 12;

    c_PartialMeltBase() : c_PartialMeltBase(std::string(), c_PartialMeltConfig{}) {}

    explicit c_PartialMeltBase(const std::string& model_name)
        : c_PartialMeltBase(model_name, c_PartialMeltConfig{}) {}

    c_PartialMeltBase(const std::string& model_name, const c_PartialMeltConfig& cfg)
        : c_PhysicsBase(model_name),
          p_solidus(cfg.solidus),
          p_liquidus(cfg.liquidus),
          p_liquid_shear(cfg.liquid_shear),
          p_liquid_viscosity(cfg.liquid_viscosity),
          p_bulk_melt_weakening(cfg.bulk_melt_weakening),
          p_liquid_bulk_modulus(cfg.liquid_bulk_modulus),
          p_liquid_bulk_modulus_derivative(cfg.liquid_bulk_modulus_derivative),
          p_liquid_density(cfg.liquid_density),
          p_density_melt_mixing(cfg.density_melt_mixing),
          p_bulk_viscosity_melt_weakening(cfg.bulk_viscosity_melt_weakening),
          p_melt_bulk_viscosity_coefficient(cfg.melt_bulk_viscosity_coefficient),
          p_melt_bulk_viscosity_exponent(cfg.melt_bulk_viscosity_exponent) {}

    ~c_PartialMeltBase() override = default;

    double get_solidus()             const noexcept { return this->p_solidus; }
    double get_liquidus()            const noexcept { return this->p_liquidus; }
    double get_liquid_shear()        const noexcept { return this->p_liquid_shear; }
    double get_liquid_viscosity()    const noexcept { return this->p_liquid_viscosity; }
    bool   get_bulk_melt_weakening() const noexcept { return this->p_bulk_melt_weakening; }
    double get_liquid_bulk_modulus() const noexcept { return this->p_liquid_bulk_modulus; }
    double get_liquid_bulk_modulus_derivative()  const noexcept { return this->p_liquid_bulk_modulus_derivative; }
    double get_liquid_density()                  const noexcept { return this->p_liquid_density; }
    bool   get_density_melt_mixing()             const noexcept { return this->p_density_melt_mixing; }
    bool   get_bulk_viscosity_melt_weakening()   const noexcept { return this->p_bulk_viscosity_melt_weakening; }
    double get_melt_bulk_viscosity_coefficient() const noexcept { return this->p_melt_bulk_viscosity_coefficient; }
    double get_melt_bulk_viscosity_exponent()    const noexcept { return this->p_melt_bulk_viscosity_exponent; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_PhysicsBase::append_config_entries(out);
        out.push_back(c_config_double("solidus_k", this->p_solidus));
        out.push_back(c_config_double("liquidus_k", this->p_liquidus));
        out.push_back(c_config_double("liquid_shear_pa", this->p_liquid_shear));
        out.push_back(c_config_double("liquid_viscosity_pas", this->p_liquid_viscosity));
        out.push_back(c_config_bool("bulk_melt_weakening", this->p_bulk_melt_weakening));
        out.push_back(c_config_double("liquid_bulk_modulus_pa", this->p_liquid_bulk_modulus));
        out.push_back(c_config_double("liquid_bulk_modulus_derivative", this->p_liquid_bulk_modulus_derivative));
        out.push_back(c_config_double("liquid_density_kg_m3", this->p_liquid_density));
        out.push_back(c_config_bool("density_melt_mixing", this->p_density_melt_mixing));
        out.push_back(c_config_bool("bulk_viscosity_melt_weakening", this->p_bulk_viscosity_melt_weakening));
        out.push_back(c_config_double("melt_bulk_viscosity_coefficient", this->p_melt_bulk_viscosity_coefficient));
        out.push_back(c_config_double("melt_bulk_viscosity_exponent", this->p_melt_bulk_viscosity_exponent));
    }

    // Volumetric melt fraction: phi = clip((T - T_sol) / (T_liq - T_sol), 0, 1).
    // A non-positive envelope (solidus >= liquidus) gives 0, fully solid; a non-finite temperature gives NaN.
    double calc_melt_fraction(double temperature) const noexcept {
        if (!std::isfinite(temperature)) { return TidalPyConstants::d_NAN; }
        const double denom = this->p_liquidus - this->p_solidus;
        if (denom <= TidalPyConstants::d_EPS) { return 0.0; }
        double phi = (temperature - this->p_solidus) / denom;
        if (phi < 0.0) { phi = 0.0; }
        if (phi > 1.0) { phi = 1.0; }
        return phi;
    }

    // Post-melt viscosity and shear modulus. Models floor both at the liquid limits (liquid_viscosity,
    // liquid_shear), return the pre-melt pair below the solidus, and return NaN for a non-finite temperature.
    virtual c_PartialMeltResult calc_partial_melt(const c_PartialMeltInputs& inputs) const = 0;

    // Bulk modulus of the melt at a pressure [Pa], the Murnaghan (1944) law K_l = K_l0 + K_l' P with
    // K_l0 = liquid_bulk_modulus and K_l' = liquid_bulk_modulus_derivative. The law is isothermal. Under tension,
    // which only the structure solve's trial central pressures reach, it continues at K_l0, so the density below
    // stays finite (the Murnaghan density has no real value past P = -K_l0 / K_l').
    double calc_liquid_bulk_modulus(double pressure) const noexcept {
        if (pressure <= 0.0) { return this->p_liquid_bulk_modulus; }
        return this->p_liquid_bulk_modulus + this->p_liquid_bulk_modulus_derivative * pressure;
    }

    // Density of the melt at a pressure [kg m-3], from the same law: rho_l0 (1 + K_l' P / K_l0)^(1 / K_l'), and
    // rho_l0 exp(P / K_l0) for K_l' = 0 or under tension, which joins it smoothly at P = 0. Its rho / (d rho / dP) is
    // calc_liquid_bulk_modulus, so a layer that is wholly melt is neutrally stratified under the bulk modulus the
    // tidal equations see.
    double calc_liquid_density(double pressure) const noexcept {
        const double k0 = this->p_liquid_bulk_modulus;
        const double kp = this->p_liquid_bulk_modulus_derivative;
        if (kp == 0.0 || pressure <= 0.0) { return this->p_liquid_density * std::exp(pressure / k0); }
        return this->p_liquid_density * std::pow(1.0 + kp * pressure / k0, 1.0 / kp);
    }

    // Density of the partially molten material [kg m-3]. The solid density law's value unless density_melt_mixing
    // is on; then the two phases mix by volume at the same pressure, (1 - phi) rho_solid + phi rho_l(P). Unchanged
    // below the solidus and without a finite temperature (no melt state). Silicate melt is less dense than its
    // source, so melt lowers a rocky layer's density; water is denser than ice I, so it raises an ice shell's.
    double calc_mixture_density(double temperature, double pressure, double solid_density) const noexcept {
        if (!this->p_density_melt_mixing) { return solid_density; }
        const double phi = this->calc_melt_fraction(temperature);
        if (!(phi > 0.0)) { return solid_density; }
        return (1.0 - phi) * solid_density + phi * this->calc_liquid_density(pressure);
    }

    // Post-melt bulk modulus [Pa]. Unchanged unless bulk_melt_weakening is on. When on, the Hashin-Shtrikman
    // (1963) bound for melt of bulk modulus K_l (calc_liquid_bulk_modulus at the local pressure) in a solid
    // framework of bulk modulus K_s:
    //     K = K_s + phi / (1 / (K_l - K_s) + (1 - phi) / (K_s + (4/3) mu_frame)),
    // with the framework shear modulus mu_frame taken as the melt-weakened (post-melt) one. While the framework
    // holds (mu_frame near the solid value) this is the upper bound for isolated melt pockets, a weak reduction
    // (about 16 percent at 10 percent melt for K_s 130, mu 60, K_l 20 GPa); once the framework's shear modulus has
    // collapsed it becomes the Reuss (Wood 1955) average of a suspension, and it reaches K_l at phi = 1.
    //
    // This is the unrelaxed (undrained) modulus: the melt is held in place over a forcing cycle. Its relaxation as
    // melt moves is a bulk rheology's job (a Zener one relaxes it toward the drained value), at the rate the bulk
    // viscosity sets (calc_bulk_viscosity_melt).
    //
    // Assumptions: two-phase isotropic aggregate; the melt fraction is the model's linear solidus-to-liquidus
    // fraction.
    double calc_bulk_modulus_melt(
            double temperature,
            double pressure,
            double premelt_bulk,
            double framework_shear) const noexcept {
        if (!this->p_bulk_melt_weakening) { return premelt_bulk; }
        const double phi = this->calc_melt_fraction(temperature);
        if (!std::isfinite(phi)) { return TidalPyConstants::d_NAN; }
        if (phi <= 0.0) { return premelt_bulk; }
        const double contrast = this->calc_liquid_bulk_modulus(pressure) - premelt_bulk;
        if (std::abs(contrast) <= TidalPyConstants::d_EPS * std::abs(premelt_bulk)) { return premelt_bulk; }
        const double shear = std::max(framework_shear, 0.0);
        const double framework_term = (1.0 - phi) / (premelt_bulk + (4.0 / 3.0) * shear);
        return premelt_bulk + phi / ((1.0 / contrast) + framework_term);
    }

    // Post-melt bulk viscosity [Pa s]. Unchanged unless bulk_viscosity_melt_weakening is on. A melt-free rock has no
    // viscous compaction; melt adds one, a matrix bulk viscosity zeta_melt = c eta / phi^n with eta the post-melt
    // shear viscosity, in series with the pre-melt bulk viscosity:
    //     1 / zeta = 1 / zeta_premelt + phi^n / (c eta).
    // n = 1 with c of order 1 is McKenzie (1984); n = 0 (a bulk viscosity of order the shear one) is closer to
    // Takei and Holtzman (2009). Unchanged below the solidus; for n > 0 the series form is continuous there. A
    // non-finite pre-melt bulk viscosity counts as no pre-melt dashpot, so melt alone sets zeta.
    double calc_bulk_viscosity_melt(
            double temperature,
            double premelt_bulk_viscosity,
            double postmelt_shear_viscosity) const noexcept {
        if (!this->p_bulk_viscosity_melt_weakening) { return premelt_bulk_viscosity; }
        const double phi = this->calc_melt_fraction(temperature);
        if (!std::isfinite(phi)) { return TidalPyConstants::d_NAN; }
        if (phi <= 0.0) { return premelt_bulk_viscosity; }
        const double melt_term = std::pow(phi, this->p_melt_bulk_viscosity_exponent)
            / (this->p_melt_bulk_viscosity_coefficient * postmelt_shear_viscosity);
        const double premelt_term = (std::isfinite(premelt_bulk_viscosity) && premelt_bulk_viscosity > 0.0)
            ? 1.0 / premelt_bulk_viscosity : 0.0;
        const double inverse = premelt_term + melt_term;
        if (!(inverse > 0.0)) { return premelt_bulk_viscosity; }
        return 1.0 / inverse;
    }

    // Element-wise over temperature and the pre-melt strengths; this is the radial sweep, one entry per slice.
    void calc_partial_melt_vectorize(
            const std::vector<double>& temperature,
            const std::vector<double>& premelt_viscosity,
            const std::vector<double>& premelt_shear,
            std::vector<c_PartialMeltResult>& out_results) const {
        const std::size_t n = temperature.size();
        if (premelt_viscosity.size() != n || premelt_shear.size() != n) {
            throw std::invalid_argument(
                "TidalPy::calc_partial_melt_vectorize: temperature, premelt_viscosity, "
                "and premelt_shear vectors must have the same length");
        }
        out_results.resize(n);
        c_PartialMeltInputs inputs;
        for (std::size_t i = 0; i < n; ++i) {
            inputs.temperature       = temperature[i];
            inputs.premelt_viscosity = premelt_viscosity[i];
            inputs.premelt_shear     = premelt_shear[i];
            out_results[i] = this->calc_partial_melt(inputs);
        }
    }

protected:
    // The shared parameters, in binary order, followed by each model's own.
    std::vector<double> envelope_params() const {
        return {this->p_solidus, this->p_liquidus, this->p_liquid_shear, this->p_liquid_viscosity,
                this->p_bulk_melt_weakening ? 1.0 : 0.0, this->p_liquid_bulk_modulus,
                this->p_liquid_bulk_modulus_derivative, this->p_liquid_density,
                this->p_density_melt_mixing ? 1.0 : 0.0, this->p_bulk_viscosity_melt_weakening ? 1.0 : 0.0,
                this->p_melt_bulk_viscosity_coefficient, this->p_melt_bulk_viscosity_exponent};
    }

    void set_envelope_params(const std::vector<double>& params) {
        this->p_solidus             = params[0];
        this->p_liquidus            = params[1];
        this->p_liquid_shear        = params[2];
        this->p_liquid_viscosity    = params[3];
        this->p_bulk_melt_weakening = (params[4] != 0.0);
        this->p_liquid_bulk_modulus = params[5];
        this->p_liquid_bulk_modulus_derivative  = params[6];
        this->p_liquid_density                  = params[7];
        this->p_density_melt_mixing             = (params[8] != 0.0);
        this->p_bulk_viscosity_melt_weakening   = (params[9] != 0.0);
        this->p_melt_bulk_viscosity_coefficient = params[10];
        this->p_melt_bulk_viscosity_exponent    = params[11];
    }

    // Floor a post-melt pair at the liquid limits.
    void apply_liquid_floor(double& viscosity, double& shear) const noexcept {
        if (viscosity <= this->p_liquid_viscosity) { viscosity = this->p_liquid_viscosity; }
        if (shear     <= this->p_liquid_shear)     { shear     = this->p_liquid_shear; }
    }

    double p_solidus;              // [K]
    double p_liquidus;             // [K]
    double p_liquid_shear;         // [Pa]
    double p_liquid_viscosity;     // [Pa·s]
    bool   p_bulk_melt_weakening;
    double p_liquid_bulk_modulus;  // [Pa]
    double p_liquid_bulk_modulus_derivative;
    double p_liquid_density;       // [kg m-3]
    bool   p_density_melt_mixing;
    bool   p_bulk_viscosity_melt_weakening;
    double p_melt_bulk_viscosity_coefficient;
    double p_melt_bulk_viscosity_exponent;
};

} // namespace tidalpy
