#pragma once
/*
 * partial_melt_.hpp — TidalPy partial-melt models (melt weakening of viscosity
 * and modulus).
 *
 * Inherits c_PartialMeltBase (partial_melt_base_.hpp), which itself inherits
 * c_PhysicsBase. Each model implements calc_partial_melt(c_PartialMeltInputs),
 * returning a c_PartialMeltResult (melt fraction, post-melt viscosity [Pa·s],
 * post-melt modulus [Pa]).
 *
 * Models (with config aliases handled by the factory):
 *   c_OffPartialMelt      (alias "none")        — no melt weakening (returns pre-melt).
 *   c_SpohnPartialMelt    (alias "fischer")     — Fischer & Spohn (1990) T-based law.
 *   c_HenningPartialMelt                        — Henning (2009/2010) three-regime law.
 *
 * All quantities MKS. The math mirrors the validated legacy implementation in
 * TidalPy/rheology/partial_melt/melting_models.py.
 *
 * References
 * ----------
 * - Fischer and Spohn (1990), Icarus 83, 39.
 * - Henning, O'Connell, and Sasselov (2009); Renaud and Henning (2018), ApJ 857, 98.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::<Model> (701-703)
 *   payload: model_name length (uint32_t) | model_name bytes | model params (doubles)
 *   Off writes [solidus, liquidus, liquid_shear]; Spohn appends its 4 scalars;
 *   Henning appends its 7 scalars. The layer observer pointer is NOT serialized.
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

// -------------------------------------------------------------------------------
// c_PartialMeltConfig — combined construction parameters for all melt models.
// Each model reads only the fields it needs.
// -------------------------------------------------------------------------------
struct c_PartialMeltConfig {
    // Shared melt envelope.
    double solidus_k       = 1600.0;   // [K]
    double liquidus_k      = 2000.0;   // [K]
    double liquid_shear_pa = 1.0e-5;   // [Pa]

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

// -------------------------------------------------------------------------------
// Lower-case a model name for case-insensitive factory lookup.
// -------------------------------------------------------------------------------
inline std::string melt_to_lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

// =====================================================================================================================
// Partial-melt models
// =====================================================================================================================

// -------------------------------------------------------------------------------
// c_OffPartialMelt — no melt weakening; post-melt strength equals pre-melt
// (alias "none"). The melt fraction is still reported for reference.
// -------------------------------------------------------------------------------
class c_OffPartialMelt : public c_PartialMeltBase {
public:
    c_OffPartialMelt() : c_PartialMeltBase("off") {}
    explicit c_OffPartialMelt(const c_PartialMeltConfig& cfg)
        : c_PartialMeltBase("off", cfg.solidus_k, cfg.liquidus_k, cfg.liquid_shear_pa) {}
    ~c_OffPartialMelt() override = default;

    c_PartialMeltResult calc_partial_melt(const c_PartialMeltInputs& in) const override {
        c_PartialMeltResult result;
        result.melt_fraction          = this->calc_melt_fraction(in.temperature_k);
        result.postmelt_viscosity     = in.premelt_viscosity;
        result.postmelt_shear_modulus = in.premelt_shear;
        return result;
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::OffPartialMelt),
            {this->p_solidus_k, this->p_liquidus_k, this->p_liquid_shear_pa});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> p = this->read_physics_binary(in, force, 3);
        this->p_solidus_k       = p[0];
        this->p_liquidus_k      = p[1];
        this->p_liquid_shear_pa = p[2];
    }
};

// -------------------------------------------------------------------------------
// c_SpohnPartialMelt — Fischer & Spohn (1990) temperature-based law (alias
// "fischer"/"fischer_spohn"). The post-melt viscosity and shear modulus depend
// only on temperature (the pre-melt values are not used), floored at the liquid
// limits.
// -------------------------------------------------------------------------------
class c_SpohnPartialMelt : public c_PartialMeltBase {
public:
    c_SpohnPartialMelt() : c_PartialMeltBase("spohn") {}
    explicit c_SpohnPartialMelt(const c_PartialMeltConfig& cfg)
        : c_PartialMeltBase("spohn", cfg.solidus_k, cfg.liquidus_k, cfg.liquid_shear_pa),
          p_fs_visc_power_slope(cfg.fs_visc_power_slope),
          p_fs_visc_power_phase(cfg.fs_visc_power_phase),
          p_fs_shear_power_slope(cfg.fs_shear_power_slope),
          p_fs_shear_power_phase(cfg.fs_shear_power_phase) {}
    ~c_SpohnPartialMelt() override = default;

    double get_visc_power_slope()  const noexcept { return this->p_fs_visc_power_slope; }
    double get_visc_power_phase()  const noexcept { return this->p_fs_visc_power_phase; }
    double get_shear_power_slope() const noexcept { return this->p_fs_shear_power_slope; }
    double get_shear_power_phase() const noexcept { return this->p_fs_shear_power_phase; }

    c_PartialMeltResult calc_partial_melt(const c_PartialMeltInputs& in) const override {
        c_PartialMeltResult result;
        result.melt_fraction = this->calc_melt_fraction(in.temperature_k);

        double post_visc = c_safe_pow(10.0,
            (this->p_fs_visc_power_slope / in.temperature_k) - this->p_fs_visc_power_phase);
        double post_shear = c_safe_pow(10.0,
            (this->p_fs_shear_power_slope / in.temperature_k) - this->p_fs_shear_power_phase);

        // Floor at the liquid limits (legacy sanity check).
        if (post_visc  <= in.liquid_viscosity)     { post_visc  = in.liquid_viscosity; }
        if (post_shear <= this->p_liquid_shear_pa) { post_shear = this->p_liquid_shear_pa; }

        result.postmelt_viscosity     = post_visc;
        result.postmelt_shear_modulus = post_shear;
        return result;
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::SpohnPartialMelt),
            {this->p_solidus_k, this->p_liquidus_k, this->p_liquid_shear_pa,
             this->p_fs_visc_power_slope, this->p_fs_visc_power_phase,
             this->p_fs_shear_power_slope, this->p_fs_shear_power_phase});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> p  = this->read_physics_binary(in, force, 7);
        this->p_solidus_k            = p[0];
        this->p_liquidus_k           = p[1];
        this->p_liquid_shear_pa      = p[2];
        this->p_fs_visc_power_slope  = p[3];
        this->p_fs_visc_power_phase  = p[4];
        this->p_fs_shear_power_slope = p[5];
        this->p_fs_shear_power_phase = p[6];
    }

protected:
    double p_fs_visc_power_slope  = 27000.0;
    double p_fs_visc_power_phase  = 1.0;
    double p_fs_shear_power_slope = 82000.0;
    double p_fs_shear_power_phase = 40.6;
};

// -------------------------------------------------------------------------------
// c_HenningPartialMelt — Henning (2009/2010) three-regime melt weakening.
//
// Below the critical melt fraction the strength weakens exponentially; in the
// transition band [crit, crit+width] a steeper "breakdown" falloff applies; above
// it the material is liquid-like. Floored at the liquid limits.
// -------------------------------------------------------------------------------
class c_HenningPartialMelt : public c_PartialMeltBase {
public:
    c_HenningPartialMelt() : c_PartialMeltBase("henning") {}
    explicit c_HenningPartialMelt(const c_PartialMeltConfig& cfg)
        : c_PartialMeltBase("henning", cfg.solidus_k, cfg.liquidus_k, cfg.liquid_shear_pa),
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

    c_PartialMeltResult calc_partial_melt(const c_PartialMeltInputs& in) const override {
        c_PartialMeltResult result;
        const double phi = this->calc_melt_fraction(in.temperature_k);
        result.melt_fraction = phi;

        const double crit       = this->p_crit_melt_frac;
        const double crit_plus  = crit + this->p_crit_melt_frac_width;
        const double break_temp = this->p_solidus_k + crit * (this->p_liquidus_k - this->p_solidus_k);

        double post_visc;
        double post_shear;
        if (phi <= 0.0) {
            // No melt: return pre-melt strengths.
            post_visc  = in.premelt_viscosity;
            post_shear = in.premelt_shear;
        } else if (phi < crit) {
            // Sub-critical exponential weakening.
            post_visc  = in.premelt_viscosity * c_safe_exp(-this->p_hn_visc_slope_1 * phi);
            post_shear = in.premelt_shear
                       * c_safe_exp((this->p_hn_shear_param_1 / in.temperature_k) - this->p_hn_shear_param_2);
        } else if (phi <= crit_plus) {
            // Transition / breakdown band: maximum sub-critical effect then a steep falloff.
            post_visc  = in.premelt_viscosity
                       * c_safe_exp(-this->p_hn_visc_slope_1 * crit)
                       * c_safe_exp(-this->p_hn_visc_falloff_slope * (phi - crit));
            post_shear = in.premelt_shear
                       * c_safe_exp((this->p_hn_shear_param_1 / break_temp) - this->p_hn_shear_param_2)
                       * c_safe_exp(-this->p_hn_shear_falloff_slope * (phi - crit));
        } else {
            // Past breakdown: liquid-like.
            post_visc  = in.liquid_viscosity;
            post_shear = this->p_liquid_shear_pa;
        }

        // Floor at the liquid limits (legacy sanity check).
        if (post_visc  <= in.liquid_viscosity)     { post_visc  = in.liquid_viscosity; }
        if (post_shear <= this->p_liquid_shear_pa) { post_shear = this->p_liquid_shear_pa; }

        result.postmelt_viscosity     = post_visc;
        result.postmelt_shear_modulus = post_shear;
        return result;
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::HenningPartialMelt),
            {this->p_solidus_k, this->p_liquidus_k, this->p_liquid_shear_pa,
             this->p_crit_melt_frac, this->p_crit_melt_frac_width,
             this->p_hn_visc_slope_1, this->p_hn_visc_falloff_slope,
             this->p_hn_shear_param_1, this->p_hn_shear_param_2,
             this->p_hn_shear_falloff_slope});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> p = this->read_physics_binary(in, force, 10);
        this->p_solidus_k              = p[0];
        this->p_liquidus_k             = p[1];
        this->p_liquid_shear_pa        = p[2];
        this->p_crit_melt_frac         = p[3];
        this->p_crit_melt_frac_width   = p[4];
        this->p_hn_visc_slope_1        = p[5];
        this->p_hn_visc_falloff_slope  = p[6];
        this->p_hn_shear_param_1       = p[7];
        this->p_hn_shear_param_2       = p[8];
        this->p_hn_shear_falloff_slope = p[9];
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

// =====================================================================================================================
// Factory
// =====================================================================================================================

enum class c_PartialMeltModel : uint8_t {
    Off     = 0,
    Spohn   = 1,
    Henning = 2,
};

// Map a (case-insensitive) model name or alias to a c_PartialMeltModel enum value.
// Throws std::invalid_argument on an unknown name.
inline c_PartialMeltModel c_partial_melt_model_from_name(const std::string& model_name) {
    const std::string name = melt_to_lower(model_name);
    if (name == "off" || name == "none")            { return c_PartialMeltModel::Off; }
    if (name == "spohn" || name == "fischer" ||
        name == "fischer_spohn")                    { return c_PartialMeltModel::Spohn; }
    if (name == "henning")                          { return c_PartialMeltModel::Henning; }
    throw std::invalid_argument("TidalPy: unknown partial-melt model name '" + model_name + "'");
}

// Build the partial-melt model named by the enum; returns an owning unique_ptr.
inline std::unique_ptr<c_PartialMeltBase> c_find_partial_melt(
        c_PartialMeltModel model, const c_PartialMeltConfig& cfg) {
    switch (model) {
        case c_PartialMeltModel::Off:     return std::make_unique<c_OffPartialMelt>(cfg);
        case c_PartialMeltModel::Spohn:   return std::make_unique<c_SpohnPartialMelt>(cfg);
        case c_PartialMeltModel::Henning: return std::make_unique<c_HenningPartialMelt>(cfg);
    }
    throw std::invalid_argument("TidalPy: unrecognised c_PartialMeltModel enum value");
}

// Name overload.
inline std::unique_ptr<c_PartialMeltBase> c_find_partial_melt(
        const std::string& model_name, const c_PartialMeltConfig& cfg) {
    return c_find_partial_melt(c_partial_melt_model_from_name(model_name), cfg);
}

// Reconstruct a partial-melt model from a binary stream (peek class id -> build -> read).
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
