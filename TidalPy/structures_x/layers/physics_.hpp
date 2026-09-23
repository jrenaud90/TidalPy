#pragma once
/*
 * physics_.hpp: c_PhysicsLayer, the mechanical layer built on c_BaseLayer.
 *
 * Adds the radial-solver classification flags, the layer temperature, the three complex Love numbers held in a
 * c_LoveNumbers struct, and the shear and bulk rheology. The material itself (static moduli, the shear law,
 * viscosities, partial melt) belongs to the layer's EOS model: the setters here that take a viscosity or
 * partial-melt model hand it to that EOS, and the static getters read it back. A complex modulus is the rheology
 * applied to the static modulus and viscosity the solved EOS reports at a radius; without a rheology it is the
 * static value as a purely real number (no dissipation). All MKS.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::PhysicsLayer (101)
 *   payload:
 *     [all c_BaseLayer fields: same byte layout as the BaseLayer binary payload]
 *     love_number_k  re, im        (double x2, 16)
 *     love_number_h  re, im        (double x2, 16)
 *     love_number_l  re, im        (double x2, 16)
 *     is_solid, is_static, is_incompressible (uint8_t x3, 3)
 *     temperature                   (double, 8)
 *     use_thermal_eos               (uint8_t, 1)
 *     use_heating                   (uint8_t, 1)
 *     eos_model       presence flag (uint8_t, 1) + (if present) its binary record
 *     shear_rheology  presence flag (uint8_t, 1) + (if present) its binary record
 *     bulk_rheology   presence flag (uint8_t, 1) + (if present) its binary record
 *   The attached material EOS model (which carries its own viscosity and partial-melt models) and the two
 *   rheologies are serialized recursively: the three presence flags belong to this payload and each nested model
 *   follows as its own record. The EOS profile data is not serialized; re-run the world EOS solve after loading.
 */

#include <complex>
#include <cstdint>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

#include "base_.hpp"
#include "love_.hpp"
#include "rheology_.hpp"

namespace tidalpy {

// c_BaseLayerConfig plus the mechanical fields.
struct c_PhysicsConfig : public c_BaseLayerConfig {
    c_LoveNumbers love_numbers;                       // k, h, l [dimensionless] placeholder
    // Radial-solver layer classification flags.
    bool          is_solid          = true;   // false for liquid layers
    bool          is_static         = true;   // use static (no dynamic terms) approximation
    bool          is_incompressible = false;  // use incompressible approximation
    // The layer temperature is 0 K until set: the cold, rigid limit of the viscosity laws.
    double        temperature = 0.0;   // [K]
    bool          use_thermal_eos = false;    // the EOS density and bulk modulus see the temperature
    bool          use_heating     = false;    // the world's heat sources act inside this layer
};

class c_PhysicsLayer : public c_BaseLayer {
public:
    c_PhysicsLayer() = default;

    explicit c_PhysicsLayer(const c_PhysicsConfig& cfg)
        : c_BaseLayer(cfg),
          p_love_numbers(cfg.love_numbers),
          p_is_solid(cfg.is_solid),
          p_is_static(cfg.is_static),
          p_is_incompressible(cfg.is_incompressible),
          p_temperature(cfg.temperature),
          p_use_thermal_eos(cfg.use_thermal_eos),
          p_use_heating(cfg.use_heating)
    {}

    ~c_PhysicsLayer() override = default;

    // The unique_ptr members delete the implicit copy assignment, which Cython's stack allocation emits from
    // freshly constructed temporaries; those always have null model pointers, so resetting on copy is safe.
    c_PhysicsLayer& operator=(const c_PhysicsLayer& other) noexcept {
        if (this != &other) {
            c_BaseLayer::operator=(other);
            this->p_love_numbers      = other.p_love_numbers;
            this->p_is_solid          = other.p_is_solid;
            this->p_is_static         = other.p_is_static;
            this->p_is_incompressible = other.p_is_incompressible;
            this->p_temperature       = other.p_temperature;
            this->p_use_thermal_eos   = other.p_use_thermal_eos;
            this->p_use_heating       = other.p_use_heating;
            // Owned model pointers cannot be copied; source temporaries always hold null.
            this->p_shear_rheology.reset();
            this->p_bulk_rheology.reset();
        }
        return *this;
    }
    c_PhysicsLayer& operator=(c_PhysicsLayer&&) noexcept = default;

    uint32_t get_layer_class_id() const noexcept override {
        return static_cast<uint32_t>(BinaryClassID::PhysicsLayer);
    }

    // Read from the layer's EOS model; NaN when none is attached.
    double get_shear_modulus_static() const noexcept {
        return this->p_eos ? this->p_eos->get_shear_modulus_static() : TidalPyConstants::d_NAN;
    }
    double get_bulk_modulus_static() const noexcept {
        return this->p_eos ? this->p_eos->get_bulk_modulus_static() : TidalPyConstants::d_NAN;
    }
    double get_shear_viscosity_static() const noexcept {
        return this->p_eos ? this->p_eos->get_shear_viscosity_static() : TidalPyConstants::d_NAN;
    }
    double get_bulk_viscosity_static() const noexcept {
        return this->p_eos ? this->p_eos->get_bulk_viscosity_static() : TidalPyConstants::d_NAN;
    }

    c_LoveNumbers        get_love_numbers()   const noexcept { return this->p_love_numbers; }
    std::complex<double> get_love_number_k()  const noexcept { return this->p_love_numbers.k; }
    std::complex<double> get_love_number_h()  const noexcept { return this->p_love_numbers.h; }
    std::complex<double> get_love_number_l()  const noexcept { return this->p_love_numbers.l; }

    bool get_is_solid()          const noexcept { return this->p_is_solid; }
    bool get_is_static()         const noexcept { return this->p_is_static; }
    bool get_is_incompressible() const noexcept { return this->p_is_incompressible; }

    // These control the shooting and propagation-matrix assumptions.
    void set_is_solid(bool value)          noexcept { this->p_is_solid = value; }
    void set_is_static(bool value)         noexcept { this->p_is_static = value; }
    void set_is_incompressible(bool value) noexcept { this->p_is_incompressible = value; }

    // Layer temperature [K] and whether the EOS density law sees it.
    double get_temperature()     const noexcept { return this->p_temperature; }
    bool   get_use_thermal_eos() const noexcept { return this->p_use_thermal_eos; }
    void set_temperature(double value)   noexcept { this->p_temperature = value; }
    void set_use_thermal_eos(bool value) noexcept { this->p_use_thermal_eos = value; }

    // Whether the world's heat sources act inside this layer during a thermal EOS solve. Off, the layer
    // generates no heat whatever models it carries.
    bool get_use_heating() const noexcept { return this->p_use_heating; }
    void set_use_heating(bool value) noexcept { this->p_use_heating = value; }

    // From the material's static constants: the rheology applied to them, or the static modulus as a purely
    // real number without one. The static viscosity is NaN until set, so a viscous rheology then returns NaN.
    std::complex<double> calc_complex_shear_modulus(double frequency) const noexcept {
        return this->apply_shear_rheology(
            this->get_shear_modulus_static(), this->get_shear_viscosity_static(), frequency);
    }

    // Same rules as the shear overload above.
    std::complex<double> calc_complex_bulk_modulus(double frequency) const noexcept {
        return this->apply_bulk_rheology(
            this->get_bulk_modulus_static(), this->get_bulk_viscosity_static(), frequency);
    }

    // The only place a complex modulus comes from: the EOS supplies the two static inputs and knows nothing
    // about frequency. Purely real, with no dissipation, when no rheology is attached.
    std::complex<double> apply_shear_rheology(
            double static_modulus, double viscosity, double frequency) const noexcept {
        if (this->p_shear_rheology) {
            return this->p_shear_rheology->calc_complex_modulus(static_modulus, viscosity, frequency);
        }
        return std::complex<double>(static_modulus, 0.0);
    }
    std::complex<double> apply_bulk_rheology(
            double static_modulus, double viscosity, double frequency) const noexcept {
        if (this->p_bulk_rheology) {
            return this->p_bulk_rheology->calc_complex_modulus(static_modulus, viscosity, frequency);
        }
        return std::complex<double>(static_modulus, 0.0);
    }

    // The rheology applied to the static modulus and viscosity the solved EOS reports at that radius.
    std::complex<double> calc_complex_shear_modulus(double radius, double frequency) const noexcept {
        double state[C_EOS_DY_VALUES];
        this->p_eos_data.evaluate(radius, state);
        return this->apply_shear_rheology(
            state[C_EOS_SHEAR_MODULUS_INDEX], state[C_EOS_SHEAR_VISCOSITY_INDEX], frequency);
    }
    std::complex<double> calc_complex_bulk_modulus(double radius, double frequency) const noexcept {
        double state[C_EOS_DY_VALUES];
        this->p_eos_data.evaluate(radius, state);
        return this->apply_bulk_rheology(
            state[C_EOS_BULK_MODULUS_INDEX], state[C_EOS_BULK_VISCOSITY_INDEX], frequency);
    }

    // Ownership transfers in, and each registers this layer as the model's observer.
    void set_shear_rheology(std::unique_ptr<c_RheologyBase> shear) {
        this->p_shear_rheology = std::move(shear);
        if (this->p_shear_rheology) { this->p_shear_rheology->set_layer_ptr(this); }
    }

    void set_bulk_rheology(std::unique_ptr<c_RheologyBase> bulk) {
        this->p_bulk_rheology = std::move(bulk);
        if (this->p_bulk_rheology) { this->p_bulk_rheology->set_layer_ptr(this); }
    }

    bool get_shear_rheology_set() const noexcept { return this->p_shear_rheology != nullptr; }
    bool get_bulk_rheology_set()  const noexcept { return this->p_bulk_rheology  != nullptr; }

    // Non-owning; null when unset.
    c_RheologyBase* get_shear_rheology_model() const noexcept { return this->p_shear_rheology.get(); }
    c_RheologyBase* get_bulk_rheology_model()  const noexcept { return this->p_bulk_rheology.get(); }

    // Shared handles, for a consumer that must outlive this layer (see the member declarations).
    std::shared_ptr<const c_RheologyBase> share_shear_rheology() const noexcept { return this->p_shear_rheology; }
    std::shared_ptr<const c_RheologyBase> share_bulk_rheology()  const noexcept { return this->p_bulk_rheology; }

    // The material owns these models, so each call hands the model to the layer's EOS; they exist so a layer
    // can be configured in one place. Attach the EOS first, or there is no material to give the model to.
    void set_shear_viscosity(std::unique_ptr<c_ViscosityBase> viscosity) {
        this->p_require_eos("a shear viscosity model")->set_shear_viscosity(std::move(viscosity));
    }
    void set_bulk_viscosity(std::unique_ptr<c_ViscosityBase> viscosity) {
        this->p_require_eos("a bulk viscosity model")->set_bulk_viscosity(std::move(viscosity));
    }
    void set_partial_melt(std::unique_ptr<c_PartialMeltBase> partial_melt) {
        this->p_require_eos("a partial-melt model")->set_partial_melt(std::move(partial_melt));
    }

    bool get_shear_viscosity_set() const noexcept { return this->get_shear_viscosity_model() != nullptr; }
    bool get_bulk_viscosity_set()  const noexcept { return this->get_bulk_viscosity_model()  != nullptr; }
    bool get_partial_melt_set()    const noexcept { return this->get_partial_melt_model()    != nullptr; }

    // Non-owning; null when unset or no EOS is attached.
    c_ViscosityBase* get_shear_viscosity_model() const noexcept {
        return this->p_eos ? this->p_eos->get_shear_viscosity_model() : nullptr;
    }
    c_ViscosityBase* get_bulk_viscosity_model() const noexcept {
        return this->p_eos ? this->p_eos->get_bulk_viscosity_model() : nullptr;
    }
    c_PartialMeltBase* get_partial_melt_model() const noexcept {
        return this->p_eos ? this->p_eos->get_partial_melt_model() : nullptr;
    }

    void write_binary(std::ostream& out) const override {
        const auto     name_len = static_cast<uint32_t>(this->p_name.size());
        const auto     mat_len  = static_cast<uint32_t>(this->p_material_name.size());
        const uint64_t payload  =
            sizeof(double)   * 2 +           // p_radius, p_mass
            sizeof(uint32_t) + name_len +    // name length + bytes
            sizeof(int32_t)  +               // layer_index
            sizeof(double)   +               // radius_inner
            sizeof(uint32_t) + mat_len +     // material_name length + bytes
            sizeof(uint8_t)  * 2 +           // is_tidal, is_volume_fixed
            sizeof(double)   +               // tidal_scale
            sizeof(double)   * 6 +           // love_number k, h, l (each: re + im)
            sizeof(uint8_t)  * 3 +           // is_solid, is_static, is_incompressible
            material_law_bytes() +           // temperature, use_thermal_eos, use_heating
            optional_binary_flag_bytes() +         // material EOS model presence flag
            this->physics_models_presence_bytes(); // shear and bulk rheology presence flags

        write_binary_header(out, static_cast<uint32_t>(BinaryClassID::PhysicsLayer), payload);

        // Same layout as the c_BaseLayer::write_binary payload.
        out.write(reinterpret_cast<const char*>(&this->p_radius), sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_mass),   sizeof(double));
        out.write(reinterpret_cast<const char*>(&name_len),       sizeof(uint32_t));
        if (name_len > 0) { out.write(this->p_name.data(), name_len); }
        const int32_t idx = static_cast<int32_t>(this->p_layer_index);
        out.write(reinterpret_cast<const char*>(&idx),                  sizeof(int32_t));
        out.write(reinterpret_cast<const char*>(&this->p_radius_inner), sizeof(double));
        out.write(reinterpret_cast<const char*>(&mat_len),              sizeof(uint32_t));
        if (mat_len > 0) { out.write(this->p_material_name.data(), mat_len); }
        const uint8_t is_tidal_byte = static_cast<uint8_t>(this->p_is_tidal);
        const uint8_t is_volume_fixed_byte = static_cast<uint8_t>(this->p_is_volume_fixed);
        out.write(reinterpret_cast<const char*>(&is_tidal_byte),       sizeof(uint8_t));
        out.write(reinterpret_cast<const char*>(&is_volume_fixed_byte), sizeof(uint8_t));
        out.write(reinterpret_cast<const char*>(&this->p_tidal_scale), sizeof(double));

        // Love numbers k, h, l
        auto write_complex = [&](const std::complex<double>& c) {
            const double re = c.real(), im = c.imag();
            out.write(reinterpret_cast<const char*>(&re), sizeof(double));
            out.write(reinterpret_cast<const char*>(&im), sizeof(double));
        };
        write_complex(this->p_love_numbers.k);
        write_complex(this->p_love_numbers.h);
        write_complex(this->p_love_numbers.l);

        // Radial-solver layer classification flags.
        const uint8_t is_solid_byte          = static_cast<uint8_t>(this->p_is_solid);
        const uint8_t is_static_byte         = static_cast<uint8_t>(this->p_is_static);
        const uint8_t is_incompressible_byte = static_cast<uint8_t>(this->p_is_incompressible);
        out.write(reinterpret_cast<const char*>(&is_solid_byte),          sizeof(uint8_t));
        out.write(reinterpret_cast<const char*>(&is_static_byte),         sizeof(uint8_t));
        out.write(reinterpret_cast<const char*>(&is_incompressible_byte), sizeof(uint8_t));
        this->write_material_law_binary(out);

        if (!out) {
            throw std::runtime_error("TidalPy: failed to write PhysicsLayer binary data");
        }

        this->write_eos_model_binary(out);
        this->write_physics_models_binary(out);
    }

    void read_binary(std::istream& in, bool force = false) override {
        c_TidalPyBaseClass::read_binary(in, force);
        // A loaded layer carries no solved profile or heating until its world solves again.
        this->clear_eos_data();
        this->p_tidal_heating = TidalPyConstants::d_NAN;

        // c_BaseLayer fields
        in.read(reinterpret_cast<char*>(&this->p_radius), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_mass),   sizeof(double));

        this->p_name = read_binary_string(in);

        int32_t idx = 0;
        in.read(reinterpret_cast<char*>(&idx), sizeof(int32_t));
        this->p_layer_index = static_cast<int>(idx);

        in.read(reinterpret_cast<char*>(&this->p_radius_inner), sizeof(double));

        this->p_material_name = read_binary_string(in);

        uint8_t is_tidal_byte = 0;
        in.read(reinterpret_cast<char*>(&is_tidal_byte), sizeof(uint8_t));
        this->p_is_tidal = static_cast<bool>(is_tidal_byte);
        uint8_t is_volume_fixed_byte = 0;
        in.read(reinterpret_cast<char*>(&is_volume_fixed_byte), sizeof(uint8_t));
        this->p_is_volume_fixed = static_cast<bool>(is_volume_fixed_byte);

        in.read(reinterpret_cast<char*>(&this->p_tidal_scale), sizeof(double));


        // Love numbers k, h, l
        auto read_complex = [&](std::complex<double>& c) {
            double re = 0.0, im = 0.0;
            in.read(reinterpret_cast<char*>(&re), sizeof(double));
            in.read(reinterpret_cast<char*>(&im), sizeof(double));
            c = std::complex<double>(re, im);
        };
        read_complex(this->p_love_numbers.k);
        read_complex(this->p_love_numbers.h);
        read_complex(this->p_love_numbers.l);

        // Radial-solver layer classification flags.
        uint8_t is_solid_byte = 0;
        uint8_t is_static_byte = 0;
        uint8_t is_incompressible_byte = 0;
        in.read(reinterpret_cast<char*>(&is_solid_byte),          sizeof(uint8_t));
        in.read(reinterpret_cast<char*>(&is_static_byte),         sizeof(uint8_t));
        in.read(reinterpret_cast<char*>(&is_incompressible_byte), sizeof(uint8_t));
        this->p_is_solid          = static_cast<bool>(is_solid_byte);
        this->p_is_static         = static_cast<bool>(is_static_byte);
        this->p_is_incompressible = static_cast<bool>(is_incompressible_byte);
        this->read_material_law_binary(in);

        if (!in) {
            throw std::runtime_error("TidalPy: failed to read PhysicsLayer binary data");
        }

        this->read_eos_model_binary(in, force);
        this->read_physics_models_binary(in, force);

        this->update_physicals();
    }

protected:
    // The attached EOS model, or a clear error naming what needed it.
    c_MaterialEOSBase* p_require_eos(const char* what) const {
        if (!this->p_eos) {
            throw std::logic_error(
                std::string("TidalPy: attach an EOS model to layer '") + this->p_name + "' before giving it "
                + what + ": the material owns it.");
        }
        return this->p_eos.get();
    }

    // Recursive (de)serialization of the two optional rheology models, shared by c_PhysicsLayer and its
    // subclasses so the section keeps one byte layout: a presence flag each, followed when set by the model's
    // own record. On read the concrete model is rebuilt through the rheology binary-dispatch factory and
    // re-registered as this layer's observer.
    void write_physics_models_binary(std::ostream& out) const {
        write_optional_binary(out, this->p_shear_rheology);
        write_optional_binary(out, this->p_bulk_rheology);
    }

    void read_physics_models_binary(std::istream& in, bool force) {
        this->p_shear_rheology =
            read_optional_binary<c_RheologyBase>(in, force, c_rheology_from_binary);
        if (this->p_shear_rheology) { this->p_shear_rheology->set_layer_ptr(this); }
        this->p_bulk_rheology =
            read_optional_binary<c_RheologyBase>(in, force, c_rheology_from_binary);
        if (this->p_bulk_rheology) { this->p_bulk_rheology->set_layer_ptr(this); }
    }

    // The layer-state scalars (temperature, use_thermal_eos, use_heating), shared with the subclasses so the
    // three layer records keep one byte layout for them.
    static constexpr uint64_t material_law_bytes() { return sizeof(double) + 2 * sizeof(uint8_t); }

    void write_material_law_binary(std::ostream& out) const {
        out.write(reinterpret_cast<const char*>(&this->p_temperature), sizeof(double));
        const uint8_t use_thermal_eos_byte = static_cast<uint8_t>(this->p_use_thermal_eos);
        out.write(reinterpret_cast<const char*>(&use_thermal_eos_byte), sizeof(uint8_t));
        const uint8_t use_heating_byte = static_cast<uint8_t>(this->p_use_heating);
        out.write(reinterpret_cast<const char*>(&use_heating_byte), sizeof(uint8_t));
    }

    void read_material_law_binary(std::istream& in) {
        in.read(reinterpret_cast<char*>(&this->p_temperature), sizeof(double));
        uint8_t use_thermal_eos_byte = 0;
        in.read(reinterpret_cast<char*>(&use_thermal_eos_byte), sizeof(uint8_t));
        this->p_use_thermal_eos = static_cast<bool>(use_thermal_eos_byte);
        uint8_t use_heating_byte = 0;
        in.read(reinterpret_cast<char*>(&use_heating_byte), sizeof(uint8_t));
        this->p_use_heating = static_cast<bool>(use_heating_byte);
    }

    // Payload bytes contributed by the two rheology presence flags (the nested records follow the payload).
    static constexpr uint64_t physics_models_presence_bytes() {
        return 2 * optional_binary_flag_bytes();
    }

    c_LoveNumbers p_love_numbers;

    // Radial-solver layer classification.
    bool p_is_solid          = true;
    bool p_is_static         = true;
    bool p_is_incompressible = false;

    // Layer state (see c_PhysicsConfig).
    double p_temperature     = 0.0;
    bool   p_use_thermal_eos = false;
    bool   p_use_heating     = false;

    // Optional rheology objects (serialized recursively via write_physics_models_binary).
    // The rheology classes are shared not unique: a radial-solver solution exported to Python keeps a
    // copy of these pointers so it can reproduce the complex moduli it was solved with, at any radius, after
    // potentially after this layer is gone.
    std::shared_ptr<c_RheologyBase> p_shear_rheology;
    std::shared_ptr<c_RheologyBase> p_bulk_rheology;
};

} // namespace tidalpy
