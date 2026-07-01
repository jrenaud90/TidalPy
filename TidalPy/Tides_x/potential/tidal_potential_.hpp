#pragma once
/*
 * tidal_potential_.hpp - TidalPy tidal-potential truncation models.
 *
 * Inherits c_TidalPotentialBase (tidal_potential_base_.hpp) -> c_PhysicsBase. Each model evaluates the
 * active tidal modes (signed frequency + potential angular factor U and its theta/phi derivatives) at a
 * point. Follows the standard physics-model pattern (enum, name factory -> unique_ptr, binary dispatch):
 *
 *   c_SyncLowEPotential  (alias "sync_low_e") - synchronous rotation, low eccentricity, no obliquity.
 *                         One mode at the orbital frequency n.
 *   c_NSRModesPotential  (alias "nsr_modes")  - moderate eccentricity (e^3), non-synchronous rotation,
 *                         no obliquity. Up to nine l = 2 modes n, 2n, 3n, 2o+-kn (o = spin frequency).
 *
 * Binary format (20-byte header + payload):
 *   c_SyncLowEPotential: class_id 1001, 0 params.
 *   c_NSRModesPotential: class_id 1002, 1 param (use_static as 0.0/1.0).
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

#include "tidal_potential_base_.hpp"

namespace tidalpy {

// -------------------------------------------------------------------------------------------------------
// c_TidalPotentialConfig - construction parameters shared by the tidal-potential models.
// -------------------------------------------------------------------------------------------------------
struct c_TidalPotentialConfig {
    bool use_static = false;   // include the time-independent potential terms (NSR truncation only)
};

// Lower-case a model name for case-insensitive factory lookup.
inline std::string tidal_potential_to_lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return text;
}


// =====================================================================================================================
// Tidal potential models
// =====================================================================================================================

// -------------------------------------------------------------------------------------------------------
// c_SyncLowEPotential - synchronous rotation, low eccentricity, no obliquity (alias "sync_low_e").
// One mode at the orbital frequency n. coeff = G * host_mass * r^2 * e / a^3.
// -------------------------------------------------------------------------------------------------------
class c_SyncLowEPotential : public c_TidalPotentialBase {
public:
    c_SyncLowEPotential() : c_TidalPotentialBase("sync_low_e") {}
    explicit c_SyncLowEPotential(const c_TidalPotentialConfig& /*cfg*/) : c_TidalPotentialBase("sync_low_e") {}
    ~c_SyncLowEPotential() override = default;

    int num_modes() const override { return 1; }

    c_TidalPotentialModeSet calc_modes(
            const c_TidalPotentialState& state,
            double radius,
            double colatitude,
            double longitude,
            double time) const override {
        const double G = tidalpy_config_ptr->d_G;

        const double cos_t = std::cos(colatitude);
        const double sin_t = std::sin(colatitude);
        const double cos2t = cos_t * cos_t;
        const double sin2t = sin_t * sin_t;

        const double p_20   = 0.5 * (3.0 * cos2t - 1.0);
        const double dp_20  = -3.0 * cos_t * sin_t;
        const double d2p_20 = 3.0 * (sin2t - cos2t);
        const double p_22   = 3.0 * (1.0 - cos2t);
        const double dp_22  = 6.0 * cos_t * sin_t;
        const double d2p_22 = 6.0 * (-sin2t + cos2t);

        const double freq   = std::fabs(state.orbital_frequency);
        const double cos_nt = std::cos(freq * time);
        const double sin_nt = std::sin(freq * time);
        const double cos_2l = std::cos(2.0 * longitude);
        const double sin_2l = std::sin(2.0 * longitude);

        const double coeff = G * state.host_mass * radius * radius * state.eccentricity
                           / (state.semi_major_axis * state.semi_major_axis * state.semi_major_axis);

        const double lt_main = 3.0 * cos_nt * cos_2l + 4.0 * sin_nt * sin_2l;
        const double lt_phi  = -6.0 * cos_nt * sin_2l + 8.0 * sin_nt * cos_2l;
        const double lt_phi2 = -12.0 * cos_nt * cos_2l - 16.0 * sin_nt * sin_2l;

        c_TidalPotentialModeSet set;
        set.num_modes = 1;
        set.mode_frequency[0] = state.orbital_frequency;
        set.potential[0] = c_PotentialPoint {
            coeff * (-1.5 * p_20 * cos_nt + 0.25 * p_22 * lt_main),
            coeff * (-1.5 * dp_20 * cos_nt + 0.25 * dp_22 * lt_main),
            coeff * (0.25 * p_22 * lt_phi),
            coeff * (-1.5 * d2p_20 * cos_nt + 0.25 * d2p_22 * lt_main),
            coeff * (0.25 * p_22 * lt_phi2),
            coeff * (0.25 * dp_22 * lt_phi)
        };
        return set;
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::SyncLowEPotential));
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }
};


// -------------------------------------------------------------------------------------------------------
// c_NSRModesPotential - moderate eccentricity (e^3), non-synchronous rotation, no obliquity
// (alias "nsr_modes"). The nine l = 2 modes n, 2n, 3n, 2o+-kn; global_coefficient = 1.5 G M r^2 / a^3.
// -------------------------------------------------------------------------------------------------------
class c_NSRModesPotential : public c_TidalPotentialBase {
public:
    c_NSRModesPotential() : c_TidalPotentialBase("nsr_modes") {}
    explicit c_NSRModesPotential(const c_TidalPotentialConfig& cfg)
        : c_TidalPotentialBase("nsr_modes"), p_use_static(cfg.use_static) {}
    ~c_NSRModesPotential() override = default;

    bool get_use_static() const { return this->p_use_static; }

    int num_modes() const override { return 9; }

    c_TidalPotentialModeSet calc_modes(
            const c_TidalPotentialState& state,
            double radius,
            double colatitude,
            double longitude,
            double time) const override {
        const double G                  = tidalpy_config_ptr->d_G;
        const double min_spin_orbit_diff = tidalpy_config_ptr->d_MIN_SPIN_ORBIT_DIFF;

        const double cos_t = std::cos(colatitude);
        const double sin_t = std::sin(colatitude);
        const double cos2t = cos_t * cos_t;
        const double sin2t = sin_t * sin_t;

        const double p_20   = 0.5 * (3.0 * cos2t - 1.0);
        const double dp_20  = -3.0 * cos_t * sin_t;
        const double d2p_20 = 3.0 * (sin2t - cos2t);
        const double p_22   = 3.0 * (1.0 - cos2t);
        const double dp_22  = 6.0 * cos_t * sin_t;
        const double d2p_22 = 6.0 * (cos2t - sin2t);

        const double cos_2l = std::cos(2.0 * longitude);
        const double sin_2l = std::sin(2.0 * longitude);

        const double ecc  = state.eccentricity;
        const double ecc2 = ecc * ecc;
        const double ecc3 = ecc * ecc2;
        const double n    = state.orbital_frequency;
        const double o    = state.spin_frequency;

        const double modes[9] = {
            n,
            2.0 * n,
            3.0 * n,
            2.0 * o + n,
            2.0 * o - n,
            2.0 * o - 2.0 * n,
            2.0 * o - 3.0 * n,
            2.0 * o - 4.0 * n,
            2.0 * o - 5.0 * n
        };
        const int mode_m[9] = {
            0,
            0,
            0,
            2,
            2,
            2,
            2,
            2,
            2
        };
        const double mode_coeffs[9] = {
            -ecc - (9.0 / 8.0) * ecc3,
            -(3.0 / 2.0) * ecc2,
            -(53.0 / 24.0) * ecc3,
            (1.0 / 288.0) * ecc3,
            (-1.0 / 12.0) * ecc + (1.0 / 96.0) * ecc3,
            (1.0 / 6.0) - (5.0 / 12.0) * ecc2,
            (7.0 / 12.0) * ecc - (41.0 / 32.0) * ecc3,
            (17.0 / 12.0) * ecc2,
            (845.0 / 288.0) * ecc3
        };
        const double static_coeff = (-1.0 / 3.0) - (1.0 / 2.0) * ecc2;
        const double global_coefficient = (3.0 / 2.0) * G * state.host_mass * radius * radius
                                        / (state.semi_major_axis * state.semi_major_axis * state.semi_major_axis);

        c_TidalPotentialModeSet set;
        set.num_modes = 9;
        for (int mode_i = 0; mode_i < 9; ++mode_i) {
            const double mode = modes[mode_i];
            const double freq = std::fabs(mode);
            set.mode_frequency[mode_i] = mode;

            const double cos_mode = std::cos(mode * time);
            const double sin_mode = std::sin(mode * time);

            double legendre, legendre_dtheta, legendre_dtheta2;
            if (mode_m[mode_i] == 0) {
                legendre = p_20; legendre_dtheta = dp_20; legendre_dtheta2 = d2p_20;
            } else {
                legendre = p_22; legendre_dtheta = dp_22; legendre_dtheta2 = d2p_22;
            }

            double longitude_coeff, longitude_coeff_dphi, longitude_coeff_dphi2;
            if (mode_m[mode_i] == 0) {
                longitude_coeff       = cos_mode;
                longitude_coeff_dphi  = 0.0;
                longitude_coeff_dphi2 = 0.0;
            } else {
                longitude_coeff       = cos_2l * cos_mode - sin_2l * sin_mode;
                longitude_coeff_dphi  = -2.0 * sin_2l * cos_mode - 2.0 * cos_2l * sin_mode;
                longitude_coeff_dphi2 = -4.0 * cos_2l * cos_mode + 4.0 * sin_2l * sin_mode;
            }

            const double mode_switch  = (this->p_use_static || freq > min_spin_orbit_diff) ? 1.0 : 0.0;
            const double switch_coeff = mode_switch * mode_coeffs[mode_i];

            double potential         = switch_coeff * longitude_coeff * legendre;
            double partial_theta     = switch_coeff * longitude_coeff * legendre_dtheta;
            double partial_phi       = switch_coeff * longitude_coeff_dphi * legendre;
            double partial_theta2    = switch_coeff * longitude_coeff * legendre_dtheta2;
            double partial_phi2      = switch_coeff * longitude_coeff_dphi2 * legendre;
            double partial_theta_phi = switch_coeff * longitude_coeff_dphi * legendre_dtheta;

            if (this->p_use_static) {
                potential      += static_coeff * p_20;
                partial_theta  += static_coeff * dp_20;
                partial_theta2 += static_coeff * d2p_20;
            }

            set.potential[mode_i] = c_PotentialPoint {
                global_coefficient * potential,
                global_coefficient * partial_theta,
                global_coefficient * partial_phi,
                global_coefficient * partial_theta2,
                global_coefficient * partial_phi2,
                global_coefficient * partial_theta_phi
            };
        }
        return set;
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::NSRModesPotential),
                                   {this->p_use_static ? 1.0 : 0.0});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 1);
        this->p_use_static = (params[0] != 0.0);
    }

protected:
    bool p_use_static = false;
};


// =====================================================================================================================
// Factory
// =====================================================================================================================

enum class c_TidalPotentialModel : uint8_t {
    SyncLowE = 0,
    NSRModes = 1,
};

// Map a (case-insensitive) model name or alias to a c_TidalPotentialModel. Throws on an unknown name.
inline c_TidalPotentialModel c_tidal_potential_model_from_name(const std::string& model_name) {
    const std::string name = tidal_potential_to_lower(model_name);
    if (name == "sync_low_e" || name == "synchronous_low_e" || name == "simple") {
        return c_TidalPotentialModel::SyncLowE;
    }
    if (name == "nsr_modes" || name == "nsr" || name == "nsr_med_eccen") {
        return c_TidalPotentialModel::NSRModes;
    }
    throw std::invalid_argument("TidalPy: unknown tidal potential model name '" + model_name + "'");
}

// Build the tidal potential model named by the enum; returns an owning unique_ptr.
inline std::unique_ptr<c_TidalPotentialBase> c_find_tidal_potential(
        c_TidalPotentialModel model, const c_TidalPotentialConfig& cfg) {
    switch (model) {
        case c_TidalPotentialModel::SyncLowE: return std::make_unique<c_SyncLowEPotential>(cfg);
        case c_TidalPotentialModel::NSRModes: return std::make_unique<c_NSRModesPotential>(cfg);
    }
    throw std::invalid_argument("TidalPy: unrecognised c_TidalPotentialModel enum value");
}

// Name overload.
inline std::unique_ptr<c_TidalPotentialBase> c_find_tidal_potential(
        const std::string& model_name, const c_TidalPotentialConfig& cfg) {
    return c_find_tidal_potential(c_tidal_potential_model_from_name(model_name), cfg);
}

// Reconstruct a tidal potential model from a binary stream (peek class id -> build -> read).
inline std::unique_ptr<c_TidalPotentialBase> c_tidal_potential_from_binary(std::istream& in, bool force = false) {
    const std::streampos start = in.tellg();
    const c_BinaryHeader header = read_binary_header(in);
    in.seekg(start);

    std::unique_ptr<c_TidalPotentialBase> model;
    switch (static_cast<BinaryClassID>(header.class_id)) {
        case BinaryClassID::SyncLowEPotential: model = std::make_unique<c_SyncLowEPotential>(); break;
        case BinaryClassID::NSRModesPotential: model = std::make_unique<c_NSRModesPotential>(); break;
        default:
            throw std::runtime_error("TidalPy: unknown tidal potential class id in binary stream");
    }
    model->read_binary(in, force);
    return model;
}

}  // namespace tidalpy
