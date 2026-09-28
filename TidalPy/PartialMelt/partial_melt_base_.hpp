#pragma once
/* Abstract base for TidalPy partial-melt models. Concrete models live in partial_melt_.hpp. All MKS.
 *
 * The melt fraction, the liquid limits, and the optional bulk-modulus weakening are model independent, so they
 * live here rather than on the models. These quantities are frequency independent: the world pipeline caches
 * them once after the EOS solve and only repeats the downstream rheology step per forcing frequency.
 *
 * References
 * ----------
 * - Fischer and Spohn (1990): temperature-based melt viscosity and shear law.
 * - Henning (2009, 2010); Renaud and Henning (2018): three-regime melt weakening.
 * - Hashin and Shtrikman (1963), J. Mech. Phys. Solids 11, 127: bounds on the moduli of a two-phase aggregate,
 *   used for the optional bulk-modulus weakening.
 * - Mavko (1980), JGR 85, 5173; Takei (2002), JGR 107, 2043: melt-bearing rock moduli, which bear on the bulk
 *   modulus far less than on the shear modulus.
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
    double liquid_bulk_modulus = 2.0e10;  // bulk modulus of the melt [Pa]

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
    static constexpr std::size_t C_NUM_ENVELOPE_PARAMS = 6;

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
          p_liquid_bulk_modulus(cfg.liquid_bulk_modulus) {}

    ~c_PartialMeltBase() override = default;

    double get_solidus()             const noexcept { return this->p_solidus; }
    double get_liquidus()            const noexcept { return this->p_liquidus; }
    double get_liquid_shear()        const noexcept { return this->p_liquid_shear; }
    double get_liquid_viscosity()    const noexcept { return this->p_liquid_viscosity; }
    bool   get_bulk_melt_weakening() const noexcept { return this->p_bulk_melt_weakening; }
    double get_liquid_bulk_modulus() const noexcept { return this->p_liquid_bulk_modulus; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_PhysicsBase::append_config_entries(out);
        out.push_back(c_config_double("solidus_k", this->p_solidus));
        out.push_back(c_config_double("liquidus_k", this->p_liquidus));
        out.push_back(c_config_double("liquid_shear_pa", this->p_liquid_shear));
        out.push_back(c_config_double("liquid_viscosity_pas", this->p_liquid_viscosity));
        out.push_back(c_config_bool("bulk_melt_weakening", this->p_bulk_melt_weakening));
        out.push_back(c_config_double("liquid_bulk_modulus_pa", this->p_liquid_bulk_modulus));
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

    // Post-melt bulk modulus [Pa]. Unchanged unless bulk_melt_weakening is on. When on, the Hashin-Shtrikman
    // (1963) bound for melt of bulk modulus K_l in a solid framework of bulk modulus K_s:
    //     K = K_s + phi / (1 / (K_l - K_s) + (1 - phi) / (K_s + (4/3) mu_frame)),
    // with the framework shear modulus mu_frame taken as the melt-weakened (post-melt) one. While the framework
    // holds (mu_frame near the solid value) this is the upper bound for isolated melt pockets, a weak reduction
    // (about 16 percent at 10 percent melt for K_s 130, mu 60, K_l 20 GPa); once the framework's shear modulus has
    // collapsed it becomes the Reuss (Wood 1955) average of a suspension, and it reaches K_l at phi = 1.
    //
    // Assumptions: two-phase isotropic aggregate; the melt fraction is the model's linear solidus-to-liquidus
    // fraction; the bulk viscosity is not changed by melt.
    double calc_bulk_modulus_melt(
            double temperature,
            double premelt_bulk,
            double framework_shear) const noexcept {
        if (!this->p_bulk_melt_weakening) { return premelt_bulk; }
        const double phi = this->calc_melt_fraction(temperature);
        if (!std::isfinite(phi)) { return TidalPyConstants::d_NAN; }
        if (phi <= 0.0) { return premelt_bulk; }
        const double contrast = this->p_liquid_bulk_modulus - premelt_bulk;
        if (std::abs(contrast) <= TidalPyConstants::d_EPS * std::abs(premelt_bulk)) { return premelt_bulk; }
        const double shear = std::max(framework_shear, 0.0);
        const double framework_term = (1.0 - phi) / (premelt_bulk + (4.0 / 3.0) * shear);
        return premelt_bulk + phi / ((1.0 / contrast) + framework_term);
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
                this->p_bulk_melt_weakening ? 1.0 : 0.0, this->p_liquid_bulk_modulus};
    }

    void set_envelope_params(const std::vector<double>& params) {
        this->p_solidus             = params[0];
        this->p_liquidus            = params[1];
        this->p_liquid_shear        = params[2];
        this->p_liquid_viscosity    = params[3];
        this->p_bulk_melt_weakening = (params[4] != 0.0);
        this->p_liquid_bulk_modulus = params[5];
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
};

} // namespace tidalpy
