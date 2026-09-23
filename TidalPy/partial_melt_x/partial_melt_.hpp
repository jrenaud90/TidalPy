#pragma once
/* TidalPy's partial-melt models (melt weakening of viscosity and modulus). All MKS.
 *
 * References
 * ----------
 * - Fischer and Spohn (1990), Icarus 83, 39.
 * - Henning, O'Connell, and Sasselov (2009); Renaud and Henning (2018), ApJ 857, 98.
 *
 * Binary payload: model name then the model's doubles. Off writes [solidus, liquidus, liquid_shear],
 * Spohn appends its 4 scalars, Henning its 7. The layer observer pointer is not serialized.
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

// Combined construction parameters; each model reads only the fields it needs.
struct c_PartialMeltConfig {
    // Shared melt envelope.
    double solidus      = 1600.0;  // [K]
    double liquidus     = 2000.0;  // [K]
    double liquid_shear = 1.0e-5;  // [Pa]

    // Spohn (Fischer & Spohn 1990) parameters.
    double fs_visc_power_slope  = 27000.0;  // [K]
    double fs_visc_power_phase  = 1.0;
    double fs_shear_power_slope = 82000.0;  // [K]
    double fs_shear_power_phase = 40.6;

    // Henning (2009/2010) parameters.
    double crit_melt_frac         = 0.5;      // [m^3/m^3]
    double crit_melt_frac_width   = 0.05;     // [m^3/m^3]
    double hn_visc_slope_1        = 13.5;
    double hn_visc_falloff_slope  = 370.0;
    double hn_shear_param_1       = 40000.0;  // [K]
    double hn_shear_param_2       = 25.0;
    double hn_shear_falloff_slope = 700.0;
};

inline std::string melt_to_lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

// No melt weakening (alias "none"); the melt fraction is still reported.
class c_OffPartialMelt : public c_PartialMeltBase {
public:
    c_OffPartialMelt() : c_PartialMeltBase("off") {}
    explicit c_OffPartialMelt(const c_PartialMeltConfig& cfg)
        : c_PartialMeltBase("off", cfg.solidus, cfg.liquidus, cfg.liquid_shear) {}
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
            out, static_cast<uint32_t>(BinaryClassID::OffPartialMelt),
            {this->p_solidus, this->p_liquidus, this->p_liquid_shear});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 3);
        this->p_solidus  = params[0];
        this->p_liquidus = params[1];
        this->p_liquid_shear = params[2];
    }
};

// Fischer and Spohn (1990) temperature law (aliases "fischer", "fischer_spohn"). The post-melt
// strengths depend only on temperature, not on the pre-melt values.
class c_SpohnPartialMelt : public c_PartialMeltBase {
public:
    c_SpohnPartialMelt() : c_PartialMeltBase("spohn") {}
    explicit c_SpohnPartialMelt(const c_PartialMeltConfig& cfg)
        : c_PartialMeltBase("spohn", cfg.solidus, cfg.liquidus, cfg.liquid_shear),
          p_fs_visc_power_slope(cfg.fs_visc_power_slope),
          p_fs_visc_power_phase(cfg.fs_visc_power_phase),
          p_fs_shear_power_slope(cfg.fs_shear_power_slope),
          p_fs_shear_power_phase(cfg.fs_shear_power_phase) {}
    ~c_SpohnPartialMelt() override = default;

    double get_visc_power_slope()  const noexcept { return this->p_fs_visc_power_slope; }
    double get_visc_power_phase()  const noexcept { return this->p_fs_visc_power_phase; }
    double get_shear_power_slope() const noexcept { return this->p_fs_shear_power_slope; }
    double get_shear_power_phase() const noexcept { return this->p_fs_shear_power_phase; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_PartialMeltBase::append_config_entries(out);
        out.push_back(c_config_double("fs_visc_power_slope_k", this->p_fs_visc_power_slope));
        out.push_back(c_config_double("fs_visc_power_phase", this->p_fs_visc_power_phase));
        out.push_back(c_config_double("fs_shear_power_slope_k", this->p_fs_shear_power_slope));
        out.push_back(c_config_double("fs_shear_power_phase", this->p_fs_shear_power_phase));
    }

    c_PartialMeltResult calc_partial_melt(const c_PartialMeltInputs& in) const override {
        c_PartialMeltResult result;
        result.melt_fraction = this->calc_melt_fraction(in.temperature);

        double post_visc = c_safe_pow(10.0,
            (this->p_fs_visc_power_slope / in.temperature) - this->p_fs_visc_power_phase);
        double post_shear = c_safe_pow(10.0,
            (this->p_fs_shear_power_slope / in.temperature) - this->p_fs_shear_power_phase);

        // Floor at the liquid limits.
        if (post_visc  <= in.liquid_viscosity)     { post_visc  = in.liquid_viscosity; }
        if (post_shear <= this->p_liquid_shear) { post_shear = this->p_liquid_shear; }

        result.postmelt_viscosity     = post_visc;
        result.postmelt_shear_modulus = post_shear;
        return result;
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::SpohnPartialMelt),
            {this->p_solidus, this->p_liquidus, this->p_liquid_shear,
             this->p_fs_visc_power_slope, this->p_fs_visc_power_phase,
             this->p_fs_shear_power_slope, this->p_fs_shear_power_phase});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 7);
        this->p_solidus              = params[0];
        this->p_liquidus             = params[1];
        this->p_liquid_shear         = params[2];
        this->p_fs_visc_power_slope  = params[3];
        this->p_fs_visc_power_phase  = params[4];
        this->p_fs_shear_power_slope = params[5];
        this->p_fs_shear_power_phase = params[6];
    }

protected:
    double p_fs_visc_power_slope  = 27000.0;
    double p_fs_visc_power_phase  = 1.0;
    double p_fs_shear_power_slope = 82000.0;
    double p_fs_shear_power_phase = 40.6;
};

// Henning (2009, 2010) three-regime melt weakening: exponential weakening below the critical melt
// fraction, a steeper breakdown falloff across [crit, crit + width], then liquid-like above it.
class c_HenningPartialMelt : public c_PartialMeltBase {
public:
    c_HenningPartialMelt() : c_PartialMeltBase("henning") {}
    explicit c_HenningPartialMelt(const c_PartialMeltConfig& cfg)
        : c_PartialMeltBase("henning", cfg.solidus, cfg.liquidus, cfg.liquid_shear),
          p_crit_melt_frac(cfg.crit_melt_frac),
          p_crit_melt_frac_width(cfg.crit_melt_frac_width),
          p_hn_visc_slope_1(cfg.hn_visc_slope_1),
          p_hn_visc_falloff_slope(cfg.hn_visc_falloff_slope),
          p_hn_shear_param_1(cfg.hn_shear_param_1),
          p_hn_shear_param_2(cfg.hn_shear_param_2),
          p_hn_shear_falloff_slope(cfg.hn_shear_falloff_slope) {}
    ~c_HenningPartialMelt() override = default;

    double get_crit_melt_frac()        const noexcept { return this->p_crit_melt_frac; }
    double get_crit_melt_frac_width()  const noexcept { return this->p_crit_melt_frac_width; }
    double get_visc_slope_1()          const noexcept { return this->p_hn_visc_slope_1; }
    double get_visc_falloff_slope()    const noexcept { return this->p_hn_visc_falloff_slope; }
    double get_shear_param_1()         const noexcept { return this->p_hn_shear_param_1; }
    double get_shear_param_2()         const noexcept { return this->p_hn_shear_param_2; }
    double get_shear_falloff_slope()   const noexcept { return this->p_hn_shear_falloff_slope; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_PartialMeltBase::append_config_entries(out);
        out.push_back(c_config_double("crit_melt_frac", this->p_crit_melt_frac));
        out.push_back(c_config_double("crit_melt_frac_width", this->p_crit_melt_frac_width));
        out.push_back(c_config_double("hn_visc_slope_1", this->p_hn_visc_slope_1));
        out.push_back(c_config_double("hn_visc_falloff_slope", this->p_hn_visc_falloff_slope));
        out.push_back(c_config_double("hn_shear_param_1_k", this->p_hn_shear_param_1));
        out.push_back(c_config_double("hn_shear_param_2", this->p_hn_shear_param_2));
        out.push_back(c_config_double("hn_shear_falloff_slope", this->p_hn_shear_falloff_slope));
    }

    c_PartialMeltResult calc_partial_melt(const c_PartialMeltInputs& in) const override {
        c_PartialMeltResult result;
        const double phi = this->calc_melt_fraction(in.temperature);
        result.melt_fraction = phi;

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
                       * c_safe_exp((this->p_hn_shear_param_1 / in.temperature) - this->p_hn_shear_param_2);
        } else if (phi <= crit_plus) {
            // Breakdown band: the full sub-critical effect, then a steep falloff.
            post_visc  = in.premelt_viscosity
                       * c_safe_exp(-this->p_hn_visc_slope_1 * crit)
                       * c_safe_exp(-this->p_hn_visc_falloff_slope * (phi - crit));
            post_shear = in.premelt_shear
                       * c_safe_exp((this->p_hn_shear_param_1 / break_temp) - this->p_hn_shear_param_2)
                       * c_safe_exp(-this->p_hn_shear_falloff_slope * (phi - crit));
        } else {
            // Past breakdown: liquid-like.
            post_visc  = in.liquid_viscosity;
            post_shear = this->p_liquid_shear;
        }

        // Floor at the liquid limits.
        if (post_visc  <= in.liquid_viscosity)     { post_visc  = in.liquid_viscosity; }
        if (post_shear <= this->p_liquid_shear) { post_shear = this->p_liquid_shear; }

        result.postmelt_viscosity     = post_visc;
        result.postmelt_shear_modulus = post_shear;
        return result;
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::HenningPartialMelt),
            {this->p_solidus, this->p_liquidus, this->p_liquid_shear,
             this->p_crit_melt_frac, this->p_crit_melt_frac_width,
             this->p_hn_visc_slope_1, this->p_hn_visc_falloff_slope,
             this->p_hn_shear_param_1, this->p_hn_shear_param_2,
             this->p_hn_shear_falloff_slope});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 10);
        this->p_solidus                = params[0];
        this->p_liquidus               = params[1];
        this->p_liquid_shear           = params[2];
        this->p_crit_melt_frac         = params[3];
        this->p_crit_melt_frac_width   = params[4];
        this->p_hn_visc_slope_1        = params[5];
        this->p_hn_visc_falloff_slope  = params[6];
        this->p_hn_shear_param_1       = params[7];
        this->p_hn_shear_param_2       = params[8];
        this->p_hn_shear_falloff_slope = params[9];
    }

protected:
    double p_crit_melt_frac         = 0.5;
    double p_crit_melt_frac_width   = 0.05;
    double p_hn_visc_slope_1        = 13.5;
    double p_hn_visc_falloff_slope  = 370.0;
    double p_hn_shear_param_1       = 40000.0;
    double p_hn_shear_param_2       = 25.0;
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
