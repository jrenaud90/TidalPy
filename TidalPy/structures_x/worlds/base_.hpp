#pragma once
/*
 * base_.hpp — c_BaseWorld: base class for all TidalPy world types.
 *
 * Inherits c_StructureBase (Utilities_x/classes_x/structure_base_.hpp).
 * Stores world-level identification and orbital/thermal scalars (albedo,
 * emissivity, obliquity, spin frequency) and provides bulk geometry and
 * equilibrium-temperature calculations.  Layer ownership and whole-planet solves
 * live in the c_LayeredWorld subclass.
 *
 * All fields are MKS (radius [m], mass [kg], angles [rad], frequency [rad/s]).
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::BaseWorld (200)
 *   payload:
 *     p_radius                 (double, 8)
 *     p_mass                   (double, 8)
 *     name_len + name          (uint32_t + bytes)
 *     world_type_len + type    (uint32_t + bytes)
 *     albedo                   (double, 8)
 *     emissivity               (double, 8)
 *     obliquity_rad            (double, 8)
 *     spin_frequency_rad_s     (double, 8)
 */

#include <cmath>
#include <complex>
#include <cstdint>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

#include "structure_base_.hpp"

// Global (1D) tidal dissipation: the analytic tide pipeline (cpl/ctl/ctl_q) is common to ALL
// world types and lives here on c_BaseWorld so even a layerless star can dissipate tidally.
// These are LIGHT headers (no eccentricity/obliquity tables); the heavy global-potential engine
// is pulled only by the out-of-line calc_tides definition in world_tides_base_.hpp.
#include "../../Tides_x/classes/tide_base_.hpp"     // tidalpy::c_TideBase, c_LoveNumbers
#include "../../Tides_x/classes/tide_result_.hpp"   // c_TideConfig, c_TideSolveConfig, c_GlobalTideResult
#include "../../Utilities_x/lookups/keys_.hpp"      // c_Key4
#include "../../Utilities_x/lookups/intmap_.hpp"    // c_IntMap (per-mode solver Love-number store)

namespace tidalpy {

// -------------------------------------------------------------------------------
// c_WorldConfig — construction parameters for c_BaseWorld (and subclasses).
// -------------------------------------------------------------------------------
struct c_WorldConfig {
    std::string name;
    std::string world_type_str       = "world";  // "star", "gasgiant", "terrestrial", ...
    double      radius_m             = 0.0;      // [m]
    double      mass_kg              = 0.0;      // [kg]
    double      albedo               = 0.3;      // [dimensionless]
    double      emissivity           = 1.0;      // [dimensionless]
    double      obliquity_rad        = 0.0;      // [rad]
    double      spin_frequency_rad_s = 0.0;      // [rad/s]
};

// -------------------------------------------------------------------------------
// c_BaseWorld
// -------------------------------------------------------------------------------
class c_BaseWorld : public c_StructureBase {
public:
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------
    c_BaseWorld() = default;

    explicit c_BaseWorld(const c_WorldConfig& cfg)
        : c_StructureBase(cfg.radius_m, cfg.mass_kg),
          p_name(cfg.name),
          p_world_type(cfg.world_type_str),
          p_albedo(cfg.albedo),
          p_emissivity(cfg.emissivity),
          p_obliquity_rad(cfg.obliquity_rad),
          p_spin_frequency_rad_s(cfg.spin_frequency_rad_s)
    {}

    ~c_BaseWorld() override = default;

    // -----------------------------------------------------------------------
    // Getters (const, MKS)
    // -----------------------------------------------------------------------
    const std::string& get_name()             const noexcept { return this->p_name; }
    const std::string& get_world_type()       const noexcept { return this->p_world_type; }
    double             get_albedo()           const noexcept { return this->p_albedo; }
    double             get_emissivity()       const noexcept { return this->p_emissivity; }
    double             get_obliquity()        const noexcept { return this->p_obliquity_rad; }
    double             get_spin_frequency()   const noexcept { return this->p_spin_frequency_rad_s; }

    // -----------------------------------------------------------------------
    // Bulk geometry (const, MKS) — use the world's own stored radius/mass.
    // -----------------------------------------------------------------------
    double calc_surface_gravity() const noexcept {
        return this->c_StructureBase::calc_surface_gravity(this->p_mass, this->p_radius);
    }
    double calc_escape_velocity() const noexcept {
        return this->c_StructureBase::calc_escape_velocity(this->p_mass, this->p_radius);
    }
    double calc_mean_density() const noexcept {
        return this->c_StructureBase::calc_mean_density(
            this->p_mass, this->calc_volume_sphere(this->p_radius));
    }

    // -----------------------------------------------------------------------
    // calc_equilibrium_temperature [K]
    //
    // Fast-rotator radiative equilibrium with a uniform-temperature surface:
    //   T_eq = [ (1 - A) * F / (4 * eps * sigma) ]^(1/4)
    //
    // where F is the incident insolation flux [W/m^2], A is the bond albedo,
    // eps is the emissivity, and sigma is the Stefan-Boltzmann constant.
    // Returns 0.0 for non-positive flux or when the config pointer is null.
    // -----------------------------------------------------------------------
    double calc_equilibrium_temperature(double insolation_flux_w_m2) const noexcept {
        if (insolation_flux_w_m2 <= 0.0 || tidalpy_config_ptr == nullptr) { return 0.0; }
        const double sigma = tidalpy_config_ptr->d_SBC;
        const double eps   = (this->p_emissivity > 0.0) ? this->p_emissivity : 1.0;
        const double absorbed = (1.0 - this->p_albedo) * insolation_flux_w_m2;
        return std::pow(absorbed / (4.0 * eps * sigma), 0.25);
    }

    // -----------------------------------------------------------------------
    // Mutators (non-const)
    // -----------------------------------------------------------------------
    void set_name(const std::string& name)      { this->p_name = name; }
    void set_spin_frequency(double freq_rad_s) noexcept { this->p_spin_frequency_rad_s = freq_rad_s; }
    void set_obliquity(double obliq_rad)        noexcept { this->p_obliquity_rad = obliq_rad; }

    // -----------------------------------------------------------------------
    // Global (1D) tidal dissipation (common to all world types)
    //
    // Attach a tide model (c_TideBase) + a tide config ([tides] truncation/degree), then call
    // calc_tides(orbital state) to collapse the global tidal modes into the total heating + the
    // three orbital potential derivatives. The analytic models (cpl/ctl/ctl_q) work on any
    // world; the rheology model needs the radial solver and is only supported on c_LayeredWorld
    // (which hides this calc_tides with its own). calc_tides is defined out-of-line in
    // world_tides_base_.hpp (it needs the heavy global-potential engine).
    // -----------------------------------------------------------------------
    void set_tide_model(std::unique_ptr<c_TideBase> tide) noexcept {
        this->p_tide         = std::move(tide);
        this->p_tides_solved = false;
    }
    bool get_tide_model_set() const noexcept { return this->p_tide != nullptr; }
    const c_TideBase* get_tide_model() const noexcept { return this->p_tide.get(); }

    void set_tide_config(const c_TideConfig& cfg) noexcept {
        this->p_tide_config  = cfg;
        this->p_tides_solved = false;
    }
    const c_TideConfig& get_tide_config() const noexcept { return this->p_tide_config; }

    // Run the global tidal solve for the supplied orbital/spin state. The analytic version
    // (here) raises if the attached model needs the radial solver; c_LayeredWorld overrides it.
    void calc_tides(const c_TideSolveConfig& state);

    bool get_tides_solved() const noexcept { return this->p_tides_solved; }
    double get_tidal_heating() const noexcept {
        return this->p_tides_solved ? this->p_tide_result.tidal_heating : TidalPyConstants::d_NAN;
    }
    double get_tidal_dU_dM() const noexcept {
        return this->p_tides_solved ? this->p_tide_result.dU_dM : TidalPyConstants::d_NAN;
    }
    double get_tidal_dU_dw() const noexcept {
        return this->p_tides_solved ? this->p_tide_result.dU_dw : TidalPyConstants::d_NAN;
    }
    double get_tidal_dU_dO() const noexcept {
        return this->p_tides_solved ? this->p_tide_result.dU_dO : TidalPyConstants::d_NAN;
    }
    int get_num_tidal_modes() const noexcept {
        return this->p_tides_solved ? this->p_tide_result.num_modes : 0;
    }

    // The whole collapsed global tidal result (heating + the three potential derivatives + mode/error
    // codes) from the most recent calc_tides, as one struct. Read it after checking get_tides_solved():
    // the fields hold their unsolved defaults (zeros) when the tides have not been solved, unlike the
    // scalar getters above which return NaN.
    const c_GlobalTideResult& get_tide_result() const noexcept { return this->p_tide_result; }

    // Complex potential Love number k_l for the tidal mode (l, m, p, q) from the most recent
    // rheology calc_tides. NaN for the analytic models (no radial solution) or an inactive mode.
    std::complex<double> get_tidal_love_k(int degree_l, int m, int p, int q) const noexcept {
        bool found = false;
        c_Key4 lmpq_key(static_cast<int16_t>(degree_l), static_cast<int16_t>(m),
                        static_cast<int16_t>(p), static_cast<int16_t>(q));
        c_LoveNumbers love = this->p_tide_solver_love.get(found, lmpq_key);
        if (!found) {
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        }
        return love.k;
    }

    // -----------------------------------------------------------------------
    // Binary I/O
    // -----------------------------------------------------------------------
    void write_binary(std::ostream& out) const override {
        this->write_world_binary(out, static_cast<uint32_t>(BinaryClassID::BaseWorld));
    }

    void read_binary(std::istream& in, bool force = false) override {
        c_TidalPyBaseClass::read_binary(in, force);
        this->read_world_fields(in);
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read BaseWorld binary data");
        }
    }

protected:
    // -----------------------------------------------------------------------
    // Shared world binary helpers (reused by subclasses with their own class id).
    // -----------------------------------------------------------------------
    // Number of payload bytes the c_BaseWorld scalar/string fields occupy.
    uint64_t world_payload_bytes() const noexcept {
        return sizeof(double) * 2                       // radius, mass
             + binary_string_bytes(this->p_name)
             + binary_string_bytes(this->p_world_type)
             + sizeof(double) * 4;                      // albedo, emissivity, obliquity, spin
    }

    // Write the c_BaseWorld payload (header + fields) for the given class id.
    void write_world_binary(std::ostream& out, uint32_t class_id) const {
        write_binary_header(out, class_id, this->world_payload_bytes());
        this->write_world_fields(out);
        if (!out) {
            throw std::runtime_error("TidalPy: failed to write world binary data");
        }
    }

    // Write only the c_BaseWorld fields (no header) so subclasses can append.
    void write_world_fields(std::ostream& out) const {
        out.write(reinterpret_cast<const char*>(&this->p_radius), sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_mass),   sizeof(double));
        write_binary_string(out, this->p_name);
        write_binary_string(out, this->p_world_type);
        out.write(reinterpret_cast<const char*>(&this->p_albedo),               sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_emissivity),           sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_obliquity_rad),        sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_spin_frequency_rad_s), sizeof(double));
    }

    // Read only the c_BaseWorld fields (header already consumed by caller).
    void read_world_fields(std::istream& in) {
        in.read(reinterpret_cast<char*>(&this->p_radius), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_mass),   sizeof(double));
        this->p_name       = read_binary_string(in);
        this->p_world_type = read_binary_string(in);
        in.read(reinterpret_cast<char*>(&this->p_albedo),               sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_emissivity),           sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_obliquity_rad),        sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_spin_frequency_rad_s), sizeof(double));
    }

    std::string p_name;
    std::string p_world_type           = "world";
    double      p_albedo               = 0.3;   // [dimensionless]
    double      p_emissivity           = 1.0;   // [dimensionless]
    double      p_obliquity_rad        = 0.0;   // [rad]
    double      p_spin_frequency_rad_s = 0.0;   // [rad/s]

    // Global (1D) tidal dissipation state (results not serialized — recompute via calc_tides).
    c_TideConfig                         p_tide_config;
    std::unique_ptr<c_TideBase>          p_tide;
    c_GlobalTideResult                   p_tide_result;
    bool                                 p_tides_solved = false;
    // Per-mode radial-solver Love numbers (k, h, l) keyed by the tidal mode (l, m, p, q),
    // retained from the most recent rheology calc_tides (empty for the analytic models).
    c_IntMap<c_Key4, c_LoveNumbers>      p_tide_solver_love;
};

} // namespace tidalpy
