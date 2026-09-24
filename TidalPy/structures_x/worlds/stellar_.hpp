#pragma once
/*
 * stellar_.hpp: c_StarWorld, a star with no internal layers and no equation of state.
 *
 * Built on c_BaseWorld. Carries scalar effective temperature and luminosity kept consistent by the
 * Stefan-Boltzmann law, plus an optional c_LuminosityBase model that derives both from the star's mass.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::StarWorld (203)
 *   payload: [all c_BaseWorld fields and its tide section] + effective_temperature (double)
 *                                                            + luminosity (double)
 *            + luminosity model presence flag (uint8_t), followed when set by the model's own complete record
 */

#include <cmath>
#include <cstdint>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>

#include "base_.hpp"
#include "../../stellar_x/luminosity_base_.hpp"   // c_LuminosityBase (luminosity model attached to the star)
#include "../../stellar_x/luminosity_.hpp"        // c_luminosity_from_binary

namespace tidalpy {

// Construction parameters for c_StarWorld: c_WorldConfig plus the stellar scalars.
struct c_StarConfig : public c_WorldConfig {
    double effective_temperature = 5772.0;   // [K] (solar default)
    double luminosity            = 0.0;      // [W] (0 => derive from T via Stefan-Boltzmann)
};

class c_StarWorld : public c_BaseWorld {
public:
    c_StarWorld() { this->p_world_type = "star"; }

    explicit c_StarWorld(const c_StarConfig& cfg)
        : c_BaseWorld(cfg),
          p_effective_temperature(cfg.effective_temperature),
          p_luminosity(cfg.luminosity)
    {
        if (this->p_world_type.empty() || this->p_world_type == "world") {
            this->p_world_type = "star";
        }
        // Derive luminosity from the effective temperature if not supplied.
        if (this->p_luminosity <= 0.0) {
            this->p_luminosity = this->calc_luminosity_from_temperature(this->p_effective_temperature);
        }
    }

    ~c_StarWorld() override = default;

    // Getters
    double get_effective_temperature() const noexcept { return this->p_effective_temperature; }
    double get_luminosity()            const noexcept { return this->p_luminosity; }

    // Stefan-Boltzmann conversions between luminosity and effective temperature: L = 4 pi R^2 sigma T^4.
    // Both return 0.0 when the config pointer is null or the input is non-positive.
    double calc_luminosity_from_temperature(double temperature) const noexcept {
        if (temperature <= 0.0 || tidalpy_config_ptr == nullptr) { return 0.0; }
        const double sigma = tidalpy_config_ptr->d_SBC;
        const double area  = this->calc_surface_area(this->p_radius);
        return area * sigma * temperature * temperature * temperature * temperature;
    }

    double calc_temperature_from_luminosity(double luminosity) const noexcept {
        if (luminosity <= 0.0 || tidalpy_config_ptr == nullptr) { return 0.0; }
        const double sigma = tidalpy_config_ptr->d_SBC;
        const double area  = this->calc_surface_area(this->p_radius);
        if (area <= 0.0 || sigma <= 0.0) { return 0.0; }
        return std::pow(luminosity / (area * sigma), 0.25);
    }

    // Mutators (keep T and L consistent via Stefan-Boltzmann)
    void set_effective_temperature(double temperature) noexcept {
        this->p_effective_temperature = temperature;
        this->p_luminosity = this->calc_luminosity_from_temperature(temperature);
    }
    void set_luminosity(double luminosity) noexcept {
        this->p_luminosity = luminosity;
        this->p_effective_temperature = this->calc_temperature_from_luminosity(luminosity);
    }

    // Optional luminosity model owned by the star. When attached, the star derives its luminosity from its own
    // mass and its effective temperature from that luminosity and its own radius. Without one the star still
    // keeps a consistent scalar temperature and luminosity pair.
    void set_luminosity_model(std::unique_ptr<c_LuminosityBase> model) noexcept {
        this->p_luminosity_model = std::move(model);
    }
    const c_LuminosityBase* get_luminosity_model() const noexcept { return this->p_luminosity_model.get(); }
    bool has_luminosity_model()                    const noexcept { return this->p_luminosity_model != nullptr; }

    // Luminosity [W] derived from the star's mass via the attached model.
    // Throws std::runtime_error if no luminosity model is attached.
    double calc_luminosity_from_mass() const {
        if (this->p_luminosity_model == nullptr) {
            throw std::runtime_error(
                "TidalPy: c_StarWorld::calc_luminosity_from_mass: no luminosity model attached "
                "(call set_luminosity_model first).");
        }
        return this->p_luminosity_model->calc_luminosity(this->get_mass());
    }

    // Effective temperature [K] derived from the star's mass (mass -> L -> T via the attached model).
    // Throws std::runtime_error if no luminosity model is attached.
    double calc_effective_temperature_from_mass() const {
        if (this->p_luminosity_model == nullptr) {
            throw std::runtime_error(
                "TidalPy: c_StarWorld::calc_effective_temperature_from_mass: no luminosity model "
                "attached (call set_luminosity_model first).");
        }
        return this->p_luminosity_model->calc_effective_temperature(this->get_mass(), this->get_radius());
    }

    // Update the star's stored luminosity and effective temperature from its mass using the attached
    // model. Throws std::runtime_error if no luminosity model is attached.
    void update_luminosity_from_mass() {
        this->p_luminosity = this->calc_luminosity_from_mass();
        this->p_effective_temperature = this->calc_temperature_from_luminosity(this->p_luminosity);
    }

    // Binary I/O
    void write_binary(std::ostream& out) const override {
        const uint64_t payload = this->world_payload_bytes() + tide_section_payload_bytes() + sizeof(double) * 2
            + optional_binary_flag_bytes();
        write_binary_header(out, static_cast<uint32_t>(BinaryClassID::StarWorld), payload);
        this->write_world_fields(out);
        this->write_tide_section(out);
        out.write(reinterpret_cast<const char*>(&this->p_effective_temperature), sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_luminosity),            sizeof(double));
        if (!out) {
            throw std::runtime_error("TidalPy: failed to write StarWorld binary data");
        }
        write_optional_binary(out, this->p_luminosity_model);
    }

    void read_binary(std::istream& in, bool force = false) override {
        c_TidalPyBaseClass::read_binary(in, force);
        this->read_world_fields(in);
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read StarWorld binary data");
        }
        this->read_tide_section(in, force);
        in.read(reinterpret_cast<char*>(&this->p_effective_temperature), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_luminosity),            sizeof(double));
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read StarWorld binary data");
        }
        this->p_luminosity_model =
            read_optional_binary<c_LuminosityBase>(in, force, c_luminosity_from_binary);
    }

protected:
    double p_effective_temperature = 5772.0;   // [K]
    double p_luminosity            = 0.0;       // [W]
    // Optional global-scale luminosity model (mass -> luminosity); serialized as an optional sub-object.
    std::unique_ptr<c_LuminosityBase> p_luminosity_model {};
};

} // namespace tidalpy
