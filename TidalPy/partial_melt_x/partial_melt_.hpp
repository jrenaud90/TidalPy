#pragma once
/* TidalPy's partial-melt models (melt weakening of viscosity and modulus). All MKS.
 *
 * References
 * ----------
 * - Fischer and Spohn (1990), Icarus 83, 39.
 * - Henning, O'Connell, and Sasselov (2009); Renaud and Henning (2018), ApJ 857, 98.
 *
 * Both temperature laws are anchored at the model's solidus, so they carry over to materials whose solidus is not
 * the 1600 K silicate value the published fits assume (the published forms are recovered at T_sol = 1600 K).
 *
 * Binary payload: model name then the model's doubles. Every model writes the 6 shared parameters first
 * [solidus, liquidus, liquid_shear, liquid_viscosity, bulk_melt_weakening (0 or 1), liquid_bulk_modulus]; Spohn
 * appends its 4 scalars, Henning its 6. The layer observer pointer is not serialized.
 */

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "partial_melt_base_.hpp"
#include "../Utilities_x/math_x/numerics_.hpp"  // c_safe_pow, c_safe_exp

namespace tidalpy {

inline std::string melt_to_lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

// No melt weakening (alias "none"); the melt fraction is still reported.
class c_OffPartialMelt : public c_PartialMeltBase {
public:
    c_OffPartialMelt() : c_PartialMeltBase("off") {}
    explicit c_OffPartialMelt(const c_PartialMeltConfig& cfg) : c_PartialMeltBase("off", cfg) {}
    ~c_OffPartialMelt() override = default;

    c_PartialMeltResult calc_partial_melt(const c_PartialMeltInputs& in) const override {
        c_PartialMeltResult result;
        result.melt_fraction      = this->calc_melt_fraction(in.temperature);
        result.postmelt_viscosity = in.premelt_viscosity;
        result.postmelt_shear_modulus = in.premelt_shear;
        return result;
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::OffPartialMelt), this->envelope_params());
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->set_envelope_params(this->read_physics_binary(in, force, C_NUM_ENVELOPE_PARAMS));
    }
};

// Fischer and Spohn (1990) temperature law (aliases "fischer", "fischer_spohn"). Above the solidus the post-melt
// strengths are 10^(log10_at_solidus + slope (1 / T - 1 / T_sol)), independent of the pre-melt values; at or below
// it the pre-melt pair is returned unchanged. Fischer and Spohn's absolute fits, 10^(27000 / T - 1) Pa s and
// 10^(82000 / T - 40.6) Pa, are this form at T_sol = 1600 K (log10 values 15.875 and 10.65); anchoring at the
// solidus keeps an icy solidus from giving 10^260 Pa.
class c_SpohnPartialMelt : public c_PartialMeltBase {
public:
    c_SpohnPartialMelt() : c_PartialMeltBase("spohn") {}
    explicit c_SpohnPartialMelt(const c_PartialMeltConfig& cfg)
        : c_PartialMeltBase("spohn", cfg),
          p_fs_visc_power_slope(cfg.fs_visc_power_slope),
          p_fs_visc_log10_at_solidus(cfg.fs_visc_log10_at_solidus),
          p_fs_shear_power_slope(cfg.fs_shear_power_slope),
          p_fs_shear_log10_at_solidus(cfg.fs_shear_log10_at_solidus) {}
    ~c_SpohnPartialMelt() override = default;

    double get_visc_power_slope()       const noexcept { return this->p_fs_visc_power_slope; }
    double get_visc_log10_at_solidus()  const noexcept { return this->p_fs_visc_log10_at_solidus; }
    double get_shear_power_slope()      const noexcept { return this->p_fs_shear_power_slope; }
    double get_shear_log10_at_solidus() const noexcept { return this->p_fs_shear_log10_at_solidus; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_PartialMeltBase::append_config_entries(out);
        out.push_back(c_config_double("fs_visc_power_slope_k", this->p_fs_visc_power_slope));
        out.push_back(c_config_double("fs_visc_log10_at_solidus", this->p_fs_visc_log10_at_solidus));
        out.push_back(c_config_double("fs_shear_power_slope_k", this->p_fs_shear_power_slope));
        out.push_back(c_config_double("fs_shear_log10_at_solidus", this->p_fs_shear_log10_at_solidus));
    }

    c_PartialMeltResult calc_partial_melt(const c_PartialMeltInputs& in) const override {
        c_PartialMeltResult result;
        const double phi = this->calc_melt_fraction(in.temperature);
        result.melt_fraction = phi;
        if (!std::isfinite(phi)) {
            result.postmelt_viscosity     = TidalPyConstants::d_NAN;
            result.postmelt_shear_modulus = TidalPyConstants::d_NAN;
            return result;
        }

        double post_visc  = in.premelt_viscosity;
        double post_shear = in.premelt_shear;
        if (phi > 0.0) {
            const double inv_temp_shift = (1.0 / in.temperature) - (1.0 / this->p_solidus);
            post_visc = c_safe_pow(10.0,
                this->p_fs_visc_log10_at_solidus + this->p_fs_visc_power_slope * inv_temp_shift);
            post_shear = c_safe_pow(10.0,
                this->p_fs_shear_log10_at_solidus + this->p_fs_shear_power_slope * inv_temp_shift);
        }
        this->apply_liquid_floor(post_visc, post_shear);

        result.postmelt_viscosity     = post_visc;
        result.postmelt_shear_modulus = post_shear;
        return result;
    }

    void write_binary(std::ostream& out) const override {
        std::vector<double> params = this->envelope_params();
        params.insert(params.end(), {this->p_fs_visc_power_slope, this->p_fs_visc_log10_at_solidus,
                                     this->p_fs_shear_power_slope, this->p_fs_shear_log10_at_solidus});
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::SpohnPartialMelt), params);
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, C_NUM_ENVELOPE_PARAMS + 4);
        this->set_envelope_params(params);
        const std::size_t i0 = C_NUM_ENVELOPE_PARAMS;
        this->p_fs_visc_power_slope       = params[i0];
        this->p_fs_visc_log10_at_solidus  = params[i0 + 1];
        this->p_fs_shear_power_slope      = params[i0 + 2];
        this->p_fs_shear_log10_at_solidus = params[i0 + 3];
    }

protected:
    double p_fs_visc_power_slope       = 27000.0;  // [K]
    double p_fs_visc_log10_at_solidus  = 15.875;   // log10 of the viscosity at the solidus [Pa s]
    double p_fs_shear_power_slope      = 82000.0;  // [K]
    double p_fs_shear_log10_at_solidus = 10.65;    // log10 of the shear modulus at the solidus [Pa]
};

// Henning (2009, 2010) three-regime melt weakening: exponential weakening below the critical melt
// fraction, a steeper breakdown falloff across [crit, crit + width], then liquid-like above it.
//
// The sub-critical shear law is mu_pre exp[b1 (1 / T - 1 / T_sol)], which is 1 at the solidus. Henning et al.
// (2009) Eq. 20, exp(40000 / T - 25), is this form for b1 = 40000 K at the silicate solidus of 1600 K
// (25 = 40000 / 1600); anchoring at the model's solidus keeps the law continuous for any other solidus.
class c_HenningPartialMelt : public c_PartialMeltBase {
public:
    c_HenningPartialMelt() : c_PartialMeltBase("henning") {}
    explicit c_HenningPartialMelt(const c_PartialMeltConfig& cfg)
        : c_PartialMeltBase("henning", cfg),
          p_crit_melt_frac(cfg.crit_melt_frac),
          p_crit_melt_frac_width(cfg.crit_melt_frac_width),
          p_hn_visc_slope_1(cfg.hn_visc_slope_1),
          p_hn_visc_falloff_slope(cfg.hn_visc_falloff_slope),
          p_hn_shear_param_1(cfg.hn_shear_param_1),
          p_hn_shear_falloff_slope(cfg.hn_shear_falloff_slope) {}
    ~c_HenningPartialMelt() override = default;

    double get_crit_melt_frac()        const noexcept { return this->p_crit_melt_frac; }
    double get_crit_melt_frac_width()  const noexcept { return this->p_crit_melt_frac_width; }
    double get_visc_slope_1()          const noexcept { return this->p_hn_visc_slope_1; }
    double get_visc_falloff_slope()    const noexcept { return this->p_hn_visc_falloff_slope; }
    double get_shear_param_1()         const noexcept { return this->p_hn_shear_param_1; }
    double get_shear_falloff_slope()   const noexcept { return this->p_hn_shear_falloff_slope; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_PartialMeltBase::append_config_entries(out);
        out.push_back(c_config_double("crit_melt_frac", this->p_crit_melt_frac));
        out.push_back(c_config_double("crit_melt_frac_width", this->p_crit_melt_frac_width));
        out.push_back(c_config_double("hn_visc_slope_1", this->p_hn_visc_slope_1));
        out.push_back(c_config_double("hn_visc_falloff_slope", this->p_hn_visc_falloff_slope));
        out.push_back(c_config_double("hn_shear_param_1_k", this->p_hn_shear_param_1));
        out.push_back(c_config_double("hn_shear_falloff_slope", this->p_hn_shear_falloff_slope));
    }

    c_PartialMeltResult calc_partial_melt(const c_PartialMeltInputs& in) const override {
        c_PartialMeltResult result;
        const double phi = this->calc_melt_fraction(in.temperature);
        result.melt_fraction = phi;
        if (!std::isfinite(phi)) {
            result.postmelt_viscosity     = TidalPyConstants::d_NAN;
            result.postmelt_shear_modulus = TidalPyConstants::d_NAN;
            return result;
        }

        const double crit       = this->p_crit_melt_frac;
        const double crit_plus  = crit + this->p_crit_melt_frac_width;
        const double break_temp = this->p_solidus + crit * (this->p_liquidus - this->p_solidus);

        double post_visc;
        double post_shear;
        if (phi <= 0.0) {
            post_visc  = in.premelt_viscosity;
            post_shear = in.premelt_shear;
        } else if (phi < crit) {
            // Sub-critical exponential weakening.
            post_visc  = in.premelt_viscosity * c_safe_exp(-this->p_hn_visc_slope_1 * phi);
            post_shear = in.premelt_shear
                       * c_safe_exp(this->p_hn_shear_param_1 * ((1.0 / in.temperature) - (1.0 / this->p_solidus)));
        } else if (phi <= crit_plus) {
            // Breakdown band: the full sub-critical effect, then a steep falloff.
            post_visc  = in.premelt_viscosity
                       * c_safe_exp(-this->p_hn_visc_slope_1 * crit)
                       * c_safe_exp(-this->p_hn_visc_falloff_slope * (phi - crit));
            post_shear = in.premelt_shear
                       * c_safe_exp(this->p_hn_shear_param_1 * ((1.0 / break_temp) - (1.0 / this->p_solidus)))
                       * c_safe_exp(-this->p_hn_shear_falloff_slope * (phi - crit));
        } else {
            // Past breakdown: liquid-like.
            post_visc  = this->p_liquid_viscosity;
            post_shear = this->p_liquid_shear;
        }
        this->apply_liquid_floor(post_visc, post_shear);

        result.postmelt_viscosity     = post_visc;
        result.postmelt_shear_modulus = post_shear;
        return result;
    }

    void write_binary(std::ostream& out) const override {
        std::vector<double> params = this->envelope_params();
        params.insert(params.end(), {this->p_crit_melt_frac, this->p_crit_melt_frac_width,
                                     this->p_hn_visc_slope_1, this->p_hn_visc_falloff_slope,
                                     this->p_hn_shear_param_1, this->p_hn_shear_falloff_slope});
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::HenningPartialMelt), params);
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, C_NUM_ENVELOPE_PARAMS + 6);
        this->set_envelope_params(params);
        const std::size_t i0 = C_NUM_ENVELOPE_PARAMS;
        this->p_crit_melt_frac         = params[i0];
        this->p_crit_melt_frac_width   = params[i0 + 1];
        this->p_hn_visc_slope_1        = params[i0 + 2];
        this->p_hn_visc_falloff_slope  = params[i0 + 3];
        this->p_hn_shear_param_1       = params[i0 + 4];
        this->p_hn_shear_falloff_slope = params[i0 + 5];
    }

protected:
    double p_crit_melt_frac         = 0.5;
    double p_crit_melt_frac_width   = 0.05;
    double p_hn_visc_slope_1        = 13.5;
    double p_hn_visc_falloff_slope  = 370.0;
    double p_hn_shear_param_1       = 40000.0;  // [K]
    double p_hn_shear_falloff_slope = 700.0;
};

enum class c_PartialMeltModel : uint8_t {
    Off     = 0,
    Spohn   = 1,
    Henning = 2,
};

// Model names are matched case-insensitively.
inline c_PartialMeltModel c_partial_melt_model_from_name(const std::string& model_name) {
    const std::string name = melt_to_lower(model_name);
    if (name == "off" || name == "none")            { return c_PartialMeltModel::Off; }
    if (name == "spohn" || name == "fischer" ||
        name == "fischer_spohn")                    { return c_PartialMeltModel::Spohn; }
    if (name == "henning")                          { return c_PartialMeltModel::Henning; }
    throw std::invalid_argument("TidalPy: unknown partial-melt model name '" + model_name + "'");
}

inline std::unique_ptr<c_PartialMeltBase> c_find_partial_melt(
        c_PartialMeltModel model, const c_PartialMeltConfig& cfg) {
    switch (model) {
        case c_PartialMeltModel::Off:     return std::make_unique<c_OffPartialMelt>(cfg);
        case c_PartialMeltModel::Spohn:   return std::make_unique<c_SpohnPartialMelt>(cfg);
        case c_PartialMeltModel::Henning: return std::make_unique<c_HenningPartialMelt>(cfg);
    }
    throw std::invalid_argument("TidalPy: unrecognised c_PartialMeltModel enum value");
}

inline std::unique_ptr<c_PartialMeltBase> c_find_partial_melt(
        const std::string& model_name, const c_PartialMeltConfig& cfg) {
    return c_find_partial_melt(c_partial_melt_model_from_name(model_name), cfg);
}

// The class id is peeked without consuming the header so the default-constructed model restores itself.
inline std::unique_ptr<c_PartialMeltBase> c_partial_melt_from_binary(std::istream& in, bool force = false) {
    const std::streampos start = in.tellg();
    const c_BinaryHeader header = read_binary_header(in);
    in.seekg(start);

    std::unique_ptr<c_PartialMeltBase> model;
    switch (static_cast<BinaryClassID>(header.class_id)) {
        case BinaryClassID::OffPartialMelt:     model = std::make_unique<c_OffPartialMelt>();     break;
        case BinaryClassID::SpohnPartialMelt:   model = std::make_unique<c_SpohnPartialMelt>();   break;
        case BinaryClassID::HenningPartialMelt: model = std::make_unique<c_HenningPartialMelt>(); break;
        default:
            throw std::runtime_error("TidalPy: unknown partial-melt class id in binary stream");
    }
    model->read_binary(in, force);
    return model;
}

}  // namespace tidalpy
