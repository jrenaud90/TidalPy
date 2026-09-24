#pragma once
/*
 * base_.hpp: c_BaseWorld, the base class for every TidalPy world type (extends c_StructureBase).
 *
 * Holds world-level identification and the orbital and thermal scalars (albedo, emissivity, obliquity, spin
 * frequency) and provides bulk geometry and equilibrium-temperature calculations. Layer ownership and
 * whole-planet solves live in c_LayeredWorld. All MKS: radius [m], mass [kg], angles [rad], frequency [rad/s].
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
 *     obliquity            (double, 8)
 *     spin_frequency     (double, 8)
 *   then the tide section, which every world class writes right after these fields:
 *     min_degree_l, max_degree_l, eccentricity_truncation, obliquity_truncation, love_method (int32_t, 4 each)
 *     layer_tidal_heating      (uint8_t, 1)
 *     love_fixed_q, love_fixed_dt, eccentricity_exact_tolerance (double, 8 each)
 *     tide model presence flag (uint8_t, 1), followed when set by the tide model's own complete record
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

// Global (1D) tidal dissipation: the analytic tide pipeline (cpl, ctl, ctl_q) is common to every world type and
// lives here on c_BaseWorld so even a layerless star can dissipate tidally. These headers are light (no
// eccentricity or obliquity tables); the heavy global-potential engine comes only with the out-of-line
// calc_tides definition in world_tides_base_.hpp.
#include "../../Tides_x/classes/tide_base_.hpp"     // tidalpy::c_TideBase, c_LoveNumbers
#include "../../Tides_x/classes/tide_result_.hpp"   // c_TideConfig, c_TideSolveConfig, c_GlobalTideResult
#include "../../Tides_x/classes/tide_.hpp"          // c_tide_from_binary (the tide model in the world record)
#include "../../Utilities_x/lookups/keys_.hpp"      // c_Key4
#include "../../Utilities_x/lookups/intmap_.hpp"    // c_IntMap (per-mode solver Love-number store)

namespace tidalpy {

// Construction parameters for c_BaseWorld and its subclasses.
struct c_WorldConfig {
    std::string name;
    std::string world_type_str = "world";  // "star", "gasgiant", "terrestrial", ...
    double      radius     = 0.0;      // [m]
    double      mass       = 0.0;      // [kg]
    double      albedo     = 0.3;      // [dimensionless]
    double      emissivity = 1.0;      // [dimensionless]
    double      obliquity  = 0.0;      // [rad]
    double      spin_frequency = 0.0;      // [rad/s]
};

class c_BaseWorld : public c_StructureBase {
public:
    // Construction
    c_BaseWorld() = default;

    explicit c_BaseWorld(const c_WorldConfig& cfg)
        : c_StructureBase(cfg.radius, cfg.mass),
          p_name(cfg.name),
          p_world_type(cfg.world_type_str),
          p_albedo(cfg.albedo),
          p_emissivity(cfg.emissivity),
          p_obliquity(cfg.obliquity),
          p_spin_frequency(cfg.spin_frequency)
    {}

    ~c_BaseWorld() override = default;

    // Getters (const, MKS)
    const std::string& get_name()             const noexcept { return this->p_name; }
    const std::string& get_world_type()       const noexcept { return this->p_world_type; }
    double             get_albedo()           const noexcept { return this->p_albedo; }
    double             get_emissivity()       const noexcept { return this->p_emissivity; }
    double             get_obliquity()        const noexcept { return this->p_obliquity; }
    double             get_spin_frequency()   const noexcept { return this->p_spin_frequency; }

    // Bulk geometry (const, MKS) from the world's own stored radius and mass.
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

    // Fast-rotator radiative equilibrium temperature [K] over a uniform-temperature surface:
    //   T_eq = [(1 - A) * F / (4 * eps * sigma)]^(1/4)
    // with F the incident insolation flux [W/m^2], A the bond albedo, eps the emissivity, and sigma the
    // Stefan-Boltzmann constant. Returns 0.0 for non-positive flux or when the config pointer is null.
    double calc_equilibrium_temperature(double insolation_flux) const noexcept {
        if (insolation_flux <= 0.0 || tidalpy_config_ptr == nullptr) { return 0.0; }
        const double sigma = tidalpy_config_ptr->d_SBC;
        const double eps   = (this->p_emissivity > 0.0) ? this->p_emissivity : 1.0;
        const double absorbed = (1.0 - this->p_albedo) * insolation_flux;
        return std::pow(absorbed / (4.0 * eps * sigma), 0.25);
    }

    // Mutators
    void set_name(const std::string& name)      { this->p_name = name; }

    // The system this world belongs to, as a source of the orbital state its tides are raised in (non-owning;
    // null outside a system), and this world's index there. The system sets both when the world is added and
    // clears them when it goes away.
    void set_tide_state_provider(const c_TideStateProvider* provider_ptr, std::size_t world_index) noexcept {
        this->p_tide_state_provider_ptr = provider_ptr;
        this->p_tide_state_index        = world_index;
    }
    const c_TideStateProvider* get_tide_state_provider() const noexcept { return this->p_tide_state_provider_ptr; }
    std::size_t get_tide_state_index() const noexcept { return this->p_tide_state_index; }

    // The tidal state the world's system gives it; false outside a system, or with no tidal host or usable orbit.
    bool get_tide_state(c_TideSolveConfig& state_out) const {
        if (this->p_tide_state_provider_ptr == nullptr) { return false; }
        return this->p_tide_state_provider_ptr->get_tide_state(this->p_tide_state_index, state_out);
    }
    void set_spin_frequency(double freq) noexcept { this->p_spin_frequency = freq; }
    void set_obliquity(double obliq)        noexcept { this->p_obliquity = obliq; }

    // Global (1D) tidal dissipation (common to all world types). Attach a tide model and a tide config, then
    // call calc_tides(orbital state) to collapse the global tidal modes into the total heating and the three
    // orbital potential derivatives. The analytic models (cpl, ctl, ctl_q) work on any world; the rheology
    // model needs the radial solver, so only c_LayeredWorld supports it (hiding this calc_tides with its own).
    // calc_tides is defined out-of-line in world_tides_base_.hpp, which carries the global-potential engine.
    void set_tide_model(std::unique_ptr<c_TideBase> tide) noexcept {
        this->p_tide         = std::move(tide);
        this->p_tides_solved = false;
    }
    bool get_tide_model_set() const noexcept { return this->p_tide != nullptr; }
    const c_TideBase* get_tide_model() const noexcept { return this->p_tide.get(); }

    // Throws std::invalid_argument for a degree range outside 2 <= min <= max <= 10 (the tabulated degrees), a
    // non-positive Love-method Q or negative time lag (NaN leaves either unset), or an exact eccentricity tolerance
    // outside (0, 1).
    void set_tide_config(const c_TideConfig& cfg) {
        validate_tide_config(cfg);
        this->p_tide_config  = cfg;
        this->p_tides_solved = false;
    }
    static void validate_tide_config(const c_TideConfig& cfg) {
        if (!((cfg.min_degree_l >= 2) && (cfg.min_degree_l <= cfg.max_degree_l) && (cfg.max_degree_l <= 10))) {
            throw std::invalid_argument(
                "TidalPy: the tidal degree range must satisfy 2 <= min_degree_l <= max_degree_l <= 10; got " +
                std::to_string(cfg.min_degree_l) + " to " + std::to_string(cfg.max_degree_l) + ".");
        }
        if (cfg.love_fixed_q <= 0.0) {
            throw std::invalid_argument("TidalPy: love_fixed_q must be positive (NaN leaves it unset).");
        }
        if (cfg.love_fixed_dt < 0.0) {
            throw std::invalid_argument("TidalPy: love_fixed_dt must not be negative (NaN leaves it unset).");
        }
        if (!((cfg.eccentricity_exact_tolerance > 0.0) && (cfg.eccentricity_exact_tolerance < 1.0))) {
            throw std::invalid_argument("TidalPy: eccentricity_exact_tolerance must be in (0, 1).");
        }
    }
    const c_TideConfig& get_tide_config() const noexcept { return this->p_tide_config; }

    // Run the global tidal solve for the supplied orbital and spin state. This analytic version throws when
    // the attached model needs the radial solver; c_LayeredWorld hides it.
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

    // The whole collapsed global tidal result (heating, the three potential derivatives, mode and error codes)
    // from the most recent calc_tides. Check get_tides_solved() first: these fields keep their unsolved
    // defaults of zero, where the scalar getters above return NaN.
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

    // Binary I/O
    void write_binary(std::ostream& out) const override {
        this->write_world_binary(out, static_cast<uint32_t>(BinaryClassID::BaseWorld));
    }

    void read_binary(std::istream& in, bool force = false) override {
        c_TidalPyBaseClass::read_binary(in, force);
        this->read_world_fields(in);
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read BaseWorld binary data");
        }
        this->read_tide_section(in, force);
    }

protected:
    // Shared world binary helpers (reused by subclasses with their own class id).
    // Number of payload bytes the c_BaseWorld scalar/string fields occupy.
    uint64_t world_payload_bytes() const noexcept {
        return sizeof(double) * 2                       // radius, mass
             + binary_string_bytes(this->p_name)
             + binary_string_bytes(this->p_world_type)
             + sizeof(double) * 4;                      // albedo, emissivity, obliquity, spin
    }

    // Write the c_BaseWorld record (header, fields, tide section) for the given class id.
    void write_world_binary(std::ostream& out, uint32_t class_id) const {
        write_binary_header(out, class_id, this->world_payload_bytes() + tide_section_payload_bytes());
        this->write_world_fields(out);
        if (!out) {
            throw std::runtime_error("TidalPy: failed to write world binary data");
        }
        this->write_tide_section(out);
    }

    // The tide section: the [tides] configuration, then the tide model as an optional sub-object, so a loaded
    // world dissipates as the saved one did. The payload counts the configuration and the presence flag; the
    // model's own record follows as a separate appended record. Tide results are not saved (recompute them).
    static constexpr uint64_t tide_section_payload_bytes() noexcept {
        return 5 * sizeof(int32_t) + sizeof(uint8_t) + 3 * sizeof(double) + optional_binary_flag_bytes();
    }

    void write_tide_section(std::ostream& out) const {
        const c_TideConfig& cfg = this->p_tide_config;
        const int32_t ints[5] = {
            cfg.min_degree_l, cfg.max_degree_l, cfg.eccentricity_truncation, cfg.obliquity_truncation,
            cfg.love_method};
        out.write(reinterpret_cast<const char*>(ints), sizeof(ints));
        const uint8_t layer_heating_byte = cfg.layer_tidal_heating ? 1 : 0;
        out.write(reinterpret_cast<const char*>(&layer_heating_byte), sizeof(uint8_t));
        out.write(reinterpret_cast<const char*>(&cfg.love_fixed_q),  sizeof(double));
        out.write(reinterpret_cast<const char*>(&cfg.love_fixed_dt), sizeof(double));
        out.write(reinterpret_cast<const char*>(&cfg.eccentricity_exact_tolerance), sizeof(double));
        if (!out) {
            throw std::runtime_error("TidalPy: failed to write world tide configuration binary data");
        }
        write_optional_binary(out, this->p_tide);
    }

    // Reads into locals and commits only once the whole section is read and valid, so a corrupt record throws
    // std::runtime_error without replacing the world's tide model or configuration.
    void read_tide_section(std::istream& in, bool force) {
        int32_t ints[5] = {0, 0, 0, 0, 0};
        in.read(reinterpret_cast<char*>(ints), sizeof(ints));
        uint8_t layer_heating_byte = 1;
        in.read(reinterpret_cast<char*>(&layer_heating_byte), sizeof(uint8_t));
        c_TideConfig cfg;
        cfg.min_degree_l            = ints[0];
        cfg.max_degree_l            = ints[1];
        cfg.eccentricity_truncation = ints[2];
        cfg.obliquity_truncation    = ints[3];
        cfg.love_method             = ints[4];
        cfg.layer_tidal_heating     = (layer_heating_byte != 0);
        in.read(reinterpret_cast<char*>(&cfg.love_fixed_q),  sizeof(double));
        in.read(reinterpret_cast<char*>(&cfg.love_fixed_dt), sizeof(double));
        in.read(reinterpret_cast<char*>(&cfg.eccentricity_exact_tolerance), sizeof(double));
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read world tide configuration binary data");
        }
        try {
            validate_tide_config(cfg);
        } catch (const std::invalid_argument& error) {
            throw std::runtime_error(std::string("TidalPy: corrupt world binary data: ") + error.what());
        }
        std::unique_ptr<c_TideBase> tide = read_optional_binary<c_TideBase>(in, force, c_tide_from_binary);

        this->p_tide_config  = cfg;
        this->p_tide         = std::move(tide);
        this->p_tides_solved = false;
        this->p_tide_result  = c_GlobalTideResult();
        this->p_tide_solver_love.clear();
    }

    // Write only the c_BaseWorld fields (no header) so subclasses can append.
    void write_world_fields(std::ostream& out) const {
        out.write(reinterpret_cast<const char*>(&this->p_radius), sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_mass),   sizeof(double));
        write_binary_string(out, this->p_name);
        write_binary_string(out, this->p_world_type);
        out.write(reinterpret_cast<const char*>(&this->p_albedo),               sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_emissivity),           sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_obliquity),        sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_spin_frequency), sizeof(double));
    }

    // Read only the c_BaseWorld fields (header already consumed by caller).
    void read_world_fields(std::istream& in) {
        in.read(reinterpret_cast<char*>(&this->p_radius), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_mass),   sizeof(double));
        this->p_name       = read_binary_string(in);
        this->p_world_type = read_binary_string(in);
        in.read(reinterpret_cast<char*>(&this->p_albedo),               sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_emissivity),           sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_obliquity),        sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_spin_frequency), sizeof(double));
    }

    std::string p_name;
    const c_TideStateProvider* p_tide_state_provider_ptr = nullptr;
    std::size_t p_tide_state_index = 0;
    std::string p_world_type = "world";
    double      p_albedo     = 0.3;   // [dimensionless]
    double      p_emissivity = 1.0;   // [dimensionless]
    double      p_obliquity  = 0.0;   // [rad]
    double      p_spin_frequency = 0.0;   // [rad/s]

    // Checks the orbital state a tidal solve is about to use: throws std::invalid_argument for an eccentricity
    // outside [0, 1) or a semi-major axis that is not positive, and warns once per world when an obliquity would be
    // ignored because the obliquity truncation is off, when the obliquity is past the range of the obliquity
    // truncation (c_obliquity_truncation_limit), or when the eccentricity is past the range of the eccentricity
    // truncation (c_eccentricity_truncation_limit).
    void p_check_tide_state(const c_TideSolveConfig& state) const {
        if (!((state.eccentricity >= 0.0) && (state.eccentricity < 1.0))) {
            throw std::invalid_argument(
                "TidalPy: world '" + this->get_name() + "' tides need an eccentricity in [0, 1); got " +
                std::to_string(state.eccentricity) + ".");
        }
        if (!(state.semi_major_axis > 0.0)) {
            throw std::invalid_argument(
                "TidalPy: world '" + this->get_name() + "' tides need a positive semi-major axis; got " +
                std::to_string(state.semi_major_axis) + " m.");
        }
        if ((state.obliquity != 0.0) && (this->p_tide_config.obliquity_truncation == 0)
                && !this->p_obliquity_off_warned) {
            this->p_obliquity_off_warned = true;
            TIDALPY_LOG_WARN(
                "TidalPy: world '{}' has an obliquity of {:.3e} rad but its obliquity truncation is off, so its "
                "obliquity tides are ignored. Set obliquity_trunc_lvl (2, 4, or 'gen') in its [tides] table or "
                "set_tide_config to include them. Shown once per world.",
                this->get_name(), state.obliquity);
        }
        const int obliquity_truncation = this->p_tide_config.obliquity_truncation;
        if (obliquity_truncation != C_OBLIQUITY_OFF) {
            const double obliquity_limit =
                c_obliquity_truncation_limit(obliquity_truncation, this->p_tide_config.max_degree_l);
            if ((std::abs(state.obliquity) > obliquity_limit) && !this->p_obliquity_range_warned) {
                this->p_obliquity_range_warned = true;
                TIDALPY_LOG_WARN(
                    "TidalPy: world '{}' has an obliquity of {:.3f} rad, past {:.3f}, where its obliquity truncation "
                    "(level {}) can misstate the tides by 10% or more. Raise obliquity_trunc_lvl in its [tides] table "
                    "or set_tide_config (recommend_obliquity_truncation picks a level for a tolerance), or use 'gen'. "
                    "Shown once per world.",
                    this->get_name(), state.obliquity, obliquity_limit, obliquity_truncation);
            }
        }
        const int eccentricity_truncation = this->p_tide_config.eccentricity_truncation;
        const double eccentricity_limit   =
            c_eccentricity_truncation_limit(eccentricity_truncation, this->p_tide_config.max_degree_l);
        if ((state.eccentricity > eccentricity_limit) && !this->p_eccentricity_range_warned) {
            this->p_eccentricity_range_warned = true;
            TIDALPY_LOG_WARN(
                "TidalPy: world '{}' has an eccentricity of {:.3f}, past {:.3f}, where its eccentricity truncation "
                "(level {}) can underestimate the tides by 10% or more. Raise eccentricity_trunc_lvl in its [tides] "
                "table or set_tide_config (recommend_eccentricity_truncation picks a level for a tolerance); past "
                "about {:.2f} use 'exact'. Shown once per world.",
                this->get_name(), state.eccentricity, eccentricity_limit, eccentricity_truncation,
                c_eccentricity_truncation_limit(50, this->p_tide_config.max_degree_l));
        }
    }
    mutable bool p_obliquity_off_warned = false;
    mutable bool p_obliquity_range_warned = false;
    mutable bool p_eccentricity_range_warned = false;

    // Global (1D) tidal dissipation state. The configuration and model are serialized (the tide section); the
    // results are not (recompute with calc_tides).
    c_TideConfig                         p_tide_config;
    std::unique_ptr<c_TideBase>          p_tide;
    c_GlobalTideResult                   p_tide_result;
    bool                                 p_tides_solved = false;
    // Per-mode radial-solver Love numbers (k, h, l) keyed by the tidal mode (l, m, p, q),
    // retained from the most recent rheology calc_tides (empty for the analytic models).
    c_IntMap<c_Key4, c_LoveNumbers>      p_tide_solver_love;
};

} // namespace tidalpy
