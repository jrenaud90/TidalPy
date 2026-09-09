#pragma once
/*
 * stellar_.hpp — c_StarWorld: a star (no internal layers, no EOS).
 *
 * Inherits c_BaseWorld. A star carries no layer stack and needs no equation of
 * state. It will hold a luminosity model (c_LuminosityBase), exposed
 * here as scalar effective temperature and luminosity fields so the structure is
 * usable before the luminosity hierarchy lands.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::StarWorld (203)
 *   payload: [all c_BaseWorld fields] + effective_temperature_k (double)
 *                                     + luminosity_w (double)
 */

#include <cmath>
#include <cstdint>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>

#include "base_.hpp"
#include "../../stellar_x/luminosity_base_.hpp"   // c_LuminosityBase (luminosity model attached to the star)

namespace tidalpy {

// -------------------------------------------------------------------------------
// c_StarConfig — extends c_WorldConfig with stellar scalars.
// -------------------------------------------------------------------------------
struct c_StarConfig : public c_WorldConfig {
    double effective_temperature_k = 5772.0;   // [K] (solar default)
    double luminosity_w            = 0.0;      // [W] (0 => derive from T via Stefan-Boltzmann)
};

// -------------------------------------------------------------------------------
// c_StarWorld
// -------------------------------------------------------------------------------
class c_StarWorld : public c_BaseWorld {
public:
    c_StarWorld() { this->p_world_type = "star"; }

    explicit c_StarWorld(const c_StarConfig& cfg)
        : c_BaseWorld(cfg),
          p_effective_temperature_k(cfg.effective_temperature_k),
          p_luminosity_w(cfg.luminosity_w)
    {
        if (this->p_world_type.empty() || this->p_world_type == "world") {
            this->p_world_type = "star";
        }
        // Derive luminosity from the effective temperature if not supplied.
        if (this->p_luminosity_w <= 0.0) {
            this->p_luminosity_w = this->calc_luminosity_from_temperature(this->p_effective_temperature_k);
        }
    }

    ~c_StarWorld() override = default;

    // -----------------------------------------------------------------------
    // Getters
    // -----------------------------------------------------------------------
    double get_effective_temperature() const noexcept { return this->p_effective_temperature_k; }
    double get_luminosity()            const noexcept { return this->p_luminosity_w; }

    // -----------------------------------------------------------------------
    // Stefan-Boltzmann luminosity <-> effective temperature (const, MKS)
    //   L = 4 * pi * R^2 * sigma * T^4
    // Returns 0.0 when the config pointer is null or inputs are non-positive.
    // -----------------------------------------------------------------------
    double calc_luminosity_from_temperature(double temperature_k) const noexcept {
        if (temperature_k <= 0.0 || tidalpy_config_ptr == nullptr) { return 0.0; }
        const double sigma = tidalpy_config_ptr->d_SBC;
        const double area  = this->calc_surface_area(this->p_radius);
        return area * sigma * temperature_k * temperature_k * temperature_k * temperature_k;
    }

    double calc_temperature_from_luminosity(double luminosity_w) const noexcept {
        if (luminosity_w <= 0.0 || tidalpy_config_ptr == nullptr) { return 0.0; }
        const double sigma = tidalpy_config_ptr->d_SBC;
        const double area  = this->calc_surface_area(this->p_radius);
        if (area <= 0.0 || sigma <= 0.0) { return 0.0; }
        return std::pow(luminosity_w / (area * sigma), 0.25);
    }

    // -----------------------------------------------------------------------
    // Mutators (keep T and L consistent via Stefan-Boltzmann)
    // -----------------------------------------------------------------------
    void set_effective_temperature(double temperature_k) noexcept {
        this->p_effective_temperature_k = temperature_k;
        this->p_luminosity_w = this->calc_luminosity_from_temperature(temperature_k);
    }
    void set_luminosity(double luminosity_w) noexcept {
        this->p_luminosity_w = luminosity_w;
        this->p_effective_temperature_k = this->calc_temperature_from_luminosity(luminosity_w);
    }

    // -----------------------------------------------------------------------
    // Luminosity model (c_LuminosityBase; a global-scale physics model owned by the star)
    //
    // When attached, the star can derive its luminosity (and effective temperature) from its own mass
    // via the model's mass-luminosity relation. The star's own radius drives the Stefan-Boltzmann
    // conversion. The model is optional; the star still keeps consistent scalar T/L without one.
    // -----------------------------------------------------------------------
    void set_luminosity_model(std::unique_ptr<c_LuminosityBase> model) noexcept {
        this->p_luminosity = std::move(model);
    }
    const c_LuminosityBase* get_luminosity_model() const noexcept { return this->p_luminosity.get(); }
    bool has_luminosity_model()                    const noexcept { return this->p_luminosity != nullptr; }

    // Luminosity [W] derived from the star's mass via the attached model.
    // Throws std::runtime_error if no luminosity model is attached.
    double calc_luminosity_from_mass() const {
        if (this->p_luminosity == nullptr) {
            throw std::runtime_error(
                "TidalPy: c_StarWorld::calc_luminosity_from_mass — no luminosity model attached "
                "(call set_luminosity_model first).");
        }
        return this->p_luminosity->calc_luminosity(this->get_mass());
    }

    // Effective temperature [K] derived from the star's mass (mass -> L -> T via the attached model).
    // Throws std::runtime_error if no luminosity model is attached.
    double calc_effective_temperature_from_mass() const {
        if (this->p_luminosity == nullptr) {
            throw std::runtime_error(
                "TidalPy: c_StarWorld::calc_effective_temperature_from_mass — no luminosity model "
                "attached (call set_luminosity_model first).");
        }
        return this->p_luminosity->calc_effective_temperature(this->get_mass(), this->get_radius());
    }

    // Update the star's stored luminosity and effective temperature from its mass using the attached
    // model. Throws std::runtime_error if no luminosity model is attached.
    void update_luminosity_from_mass() {
        this->p_luminosity_w = this->calc_luminosity_from_mass();
        this->p_effective_temperature_k = this->calc_temperature_from_luminosity(this->p_luminosity_w);
    }

    // -----------------------------------------------------------------------
    // Binary I/O
    // -----------------------------------------------------------------------
    void write_binary(std::ostream& out) const override {
        const uint64_t payload = this->world_payload_bytes() + sizeof(double) * 2;
        write_binary_header(out, static_cast<uint32_t>(BinaryClassID::StarWorld), payload);
        this->write_world_fields(out);
        out.write(reinterpret_cast<const char*>(&this->p_effective_temperature_k), sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_luminosity_w),            sizeof(double));
        if (!out) {
            throw std::runtime_error("TidalPy: failed to write StarWorld binary data");
        }
    }

    void read_binary(std::istream& in, bool force = false) override {
        c_TidalPyBaseClass::read_binary(in, force);
        this->read_world_fields(in);
        in.read(reinterpret_cast<char*>(&this->p_effective_temperature_k), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_luminosity_w),            sizeof(double));
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read StarWorld binary data");
        }
    }

protected:
    double p_effective_temperature_k = 5772.0;   // [K]
    double p_luminosity_w            = 0.0;       // [W]
    // Optional global-scale luminosity model (mass -> luminosity); not serialized (reattach after load).
    std::unique_ptr<c_LuminosityBase> p_luminosity {};
};

} // namespace tidalpy
