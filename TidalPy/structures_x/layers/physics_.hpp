#pragma once
/*
 * physics_.hpp: c_PhysicsLayer, the mechanical-properties layer built on c_BaseLayer.
 *
 * Adds the static mechanical properties (shear modulus, bulk modulus, shear and bulk viscosity) and the three
 * complex Love numbers held in a c_LoveNumbers struct. Attached c_RheologyBase models give frequency-dependent
 * complex moduli; without one, calc_complex_shear/bulk_modulus return the static value as a purely real complex
 * number (no dissipation). All MKS.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::PhysicsLayer (101)
 *   payload:
 *     [all c_BaseLayer fields: same byte layout as the BaseLayer binary payload]
 *     shear_modulus_static      (double, 8)
 *     bulk_modulus_static       (double, 8)
 *     shear_viscosity_static   (double, 8)
 *     bulk_viscosity_static    (double, 8)
 *     love_number_k  re, im        (double×2, 16)
 *     love_number_h  re, im        (double×2, 16)
 *     love_number_l  re, im        (double×2, 16)
 *     is_solid, is_static, is_incompressible (uint8_t×3, 3)
 *     eos_model       presence flag (uint8_t, 1) + (if present) its binary record
 *     shear_rheology  presence flag (uint8_t, 1) + (if present) its binary record
 *     bulk_rheology   presence flag (uint8_t, 1) + (if present) its binary record
 *     shear_viscosity presence flag (uint8_t, 1) + (if present) its binary record
 *     bulk_viscosity  presence flag (uint8_t, 1) + (if present) its binary record
 *     partial_melt    presence flag (uint8_t, 1) + (if present) its binary record
 *   The attached material EOS model and physics models (rheology, viscosity, partial melt) are serialized
 *   recursively: the six presence flags belong to this payload and each nested model follows as its own record.
 *   The EOS profile data is not serialized; re-run the world EOS solve after loading.
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
#include "viscosity_.hpp"      // c_ViscosityBase (shear/bulk pre-melt viscosity)
#include "partial_melt_.hpp"   // c_PartialMeltBase (melt weakening)

namespace tidalpy {

// Construction parameters for c_PhysicsLayer: c_BaseLayerConfig plus the mechanical property fields.
struct c_PhysicsConfig : public c_BaseLayerConfig {
    double        shear_modulus_static = 0.0;   // [Pa]
    double        bulk_modulus_static  = 0.0;   // [Pa]
    double        shear_viscosity_static = TidalPyConstants::d_NAN;   // [Pa·s], NaN until set
    double        bulk_viscosity_static  = TidalPyConstants::d_NAN;   // [Pa·s], NaN until set
    c_LoveNumbers love_numbers;                       // k, h, l [dimensionless] placeholder
    // Radial-solver layer classification flags.
    bool          is_solid          = true;   // false for liquid layers
    bool          is_static         = true;   // use static (no dynamic terms) approximation
    bool          is_incompressible = false;  // use incompressible approximation
};

class c_PhysicsLayer : public c_BaseLayer {
public:
    // Construction
    c_PhysicsLayer() = default;

    explicit c_PhysicsLayer(const c_PhysicsConfig& cfg)
        : c_BaseLayer(cfg),
          p_shear_modulus_static(cfg.shear_modulus_static),
          p_bulk_modulus_static(cfg.bulk_modulus_static),
          p_shear_viscosity_static(cfg.shear_viscosity_static),
          p_bulk_viscosity_static(cfg.bulk_viscosity_static),
          p_love_numbers(cfg.love_numbers),
          p_is_solid(cfg.is_solid),
          p_is_static(cfg.is_static),
          p_is_incompressible(cfg.is_incompressible)
    {}

    ~c_PhysicsLayer() override = default;

    // The unique_ptr members delete the implicit copy assignment. Cython's stack allocation emits it from
    // freshly constructed temporaries, which always have null model pointers, so resetting on copy is safe.
    c_PhysicsLayer& operator=(const c_PhysicsLayer& other) noexcept {
        if (this != &other) {
            c_BaseLayer::operator=(other);
            this->p_shear_modulus_static = other.p_shear_modulus_static;
            this->p_bulk_modulus_static  = other.p_bulk_modulus_static;
            this->p_shear_viscosity_static = other.p_shear_viscosity_static;
            this->p_bulk_viscosity_static = other.p_bulk_viscosity_static;
            this->p_love_numbers      = other.p_love_numbers;
            this->p_is_solid          = other.p_is_solid;
            this->p_is_static         = other.p_is_static;
            this->p_is_incompressible = other.p_is_incompressible;
            // Owned model pointers cannot be copied; source temporaries always have null ptrs.
            this->p_shear_rheology.reset();
            this->p_bulk_rheology.reset();
            this->p_shear_viscosity.reset();
            this->p_bulk_viscosity.reset();
            this->p_partial_melt.reset();
        }
        return *this;
    }
    c_PhysicsLayer& operator=(c_PhysicsLayer&&) noexcept = default;

    uint32_t get_layer_class_id() const noexcept override {
        return static_cast<uint32_t>(BinaryClassID::PhysicsLayer);
    }

    // Static mechanical property getters (const, MKS)
    double        get_shear_modulus_static()   const noexcept { return this->p_shear_modulus_static; }
    double        get_bulk_modulus_static()    const noexcept { return this->p_bulk_modulus_static; }
    double        get_shear_viscosity_static() const noexcept { return this->p_shear_viscosity_static; }
    double        get_bulk_viscosity_static()  const noexcept { return this->p_bulk_viscosity_static; }

    // Love number getters: full struct or individual components
    c_LoveNumbers        get_love_numbers()   const noexcept { return this->p_love_numbers; }
    std::complex<double> get_love_number_k()  const noexcept { return this->p_love_numbers.k; }
    std::complex<double> get_love_number_h()  const noexcept { return this->p_love_numbers.h; }
    std::complex<double> get_love_number_l()  const noexcept { return this->p_love_numbers.l; }

    // Radial-solver layer classification getters.
    bool get_is_solid()          const noexcept { return this->p_is_solid; }
    bool get_is_static()         const noexcept { return this->p_is_static; }
    bool get_is_incompressible() const noexcept { return this->p_is_incompressible; }

    // Radial-solver layer classification setters (control the shooting / propagation-matrix assumptions).
    void set_is_solid(bool value)          noexcept { this->p_is_solid = value; }
    void set_is_static(bool value)         noexcept { this->p_is_static = value; }
    void set_is_incompressible(bool value) noexcept { this->p_is_incompressible = value; }

    // Complex shear modulus [Pa] at a forcing frequency, from the layer-constant static modulus and viscosity.
    // Without a rheology the static modulus is returned with no imaginary part. The static viscosity is NaN
    // until set, so a viscous rheology returns NaN; the radius-resolved overload uses the EOS profile instead.
    std::complex<double> calc_complex_shear_modulus(double frequency) const noexcept {
        if (this->p_shear_rheology) {
            return this->p_shear_rheology->calc_complex_modulus(
                this->p_shear_modulus_static, this->p_shear_viscosity_static, frequency);
        }
        return std::complex<double>(this->p_shear_modulus_static, 0.0);
    }

    // Complex bulk modulus [Pa] at a forcing frequency; same rules as the shear overload above.
    std::complex<double> calc_complex_bulk_modulus(double frequency) const noexcept {
        if (this->p_bulk_rheology) {
            return this->p_bulk_rheology->calc_complex_modulus(
                this->p_bulk_modulus_static, this->p_bulk_viscosity_static, frequency);
        }
        return std::complex<double>(this->p_bulk_modulus_static, 0.0);
    }

    // Radius-resolved complex moduli [Pa], using the post-melt static modulus and viscosity stored at radius by
    // the world EOS solve rather than the layer-constant values. Feeds the radial Love-number solve. Purely real
    // (no dissipation) when no rheology is attached; NaN until the viscoelastic state is populated.
    std::complex<double> calc_complex_shear_modulus(
            double radius, double frequency) const noexcept {
        const double static_modulus = this->get_shear_modulus(radius);    // post-melt
        const double viscosity      = this->get_shear_viscosity(radius);  // post-melt
        if (this->p_shear_rheology) {
            return this->p_shear_rheology->calc_complex_modulus(static_modulus, viscosity, frequency);
        }
        return std::complex<double>(static_modulus, 0.0);
    }

    std::complex<double> calc_complex_bulk_modulus(
            double radius, double frequency) const noexcept {
        const double static_modulus = this->get_bulk_modulus(radius);    // post-melt
        const double viscosity      = this->get_bulk_viscosity(radius);  // post-melt
        if (this->p_bulk_rheology) {
            return this->p_bulk_rheology->calc_complex_modulus(static_modulus, viscosity, frequency);
        }
        return std::complex<double>(static_modulus, 0.0);
    }

    // Rheology setters (transfer ownership; each registers this layer as the model's observer).
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

    // Non-owning observer pointers (nullptr if unset).
    c_RheologyBase* get_shear_rheology_model() const noexcept { return this->p_shear_rheology.get(); }
    c_RheologyBase* get_bulk_rheology_model()  const noexcept { return this->p_bulk_rheology.get(); }

    // Viscosity and partial-melt setters (transfer ownership; each registers this layer as the observer). The
    // viscosity models supply the pre-melt viscosities at (T, P) and the partial-melt model weakens the static
    // moduli and viscosities; both feed the frequency-independent state built by the world EOS solve.
    void set_shear_viscosity(std::unique_ptr<c_ViscosityBase> viscosity) {
        this->p_shear_viscosity = std::move(viscosity);
        if (this->p_shear_viscosity) { this->p_shear_viscosity->set_layer_ptr(this); }
    }

    void set_bulk_viscosity(std::unique_ptr<c_ViscosityBase> viscosity) {
        this->p_bulk_viscosity = std::move(viscosity);
        if (this->p_bulk_viscosity) { this->p_bulk_viscosity->set_layer_ptr(this); }
    }

    void set_partial_melt(std::unique_ptr<c_PartialMeltBase> partial_melt) {
        this->p_partial_melt = std::move(partial_melt);
        if (this->p_partial_melt) { this->p_partial_melt->set_layer_ptr(this); }
    }

    bool get_shear_viscosity_set() const noexcept { return this->p_shear_viscosity != nullptr; }
    bool get_bulk_viscosity_set()  const noexcept { return this->p_bulk_viscosity  != nullptr; }
    bool get_partial_melt_set()    const noexcept { return this->p_partial_melt    != nullptr; }

    // Non-owning observer pointers (nullptr if unset), read by the world EOS solve for the viscoelastic state.
    c_ViscosityBase*   get_shear_viscosity_model() const noexcept { return this->p_shear_viscosity.get(); }
    c_ViscosityBase*   get_bulk_viscosity_model()  const noexcept { return this->p_bulk_viscosity.get(); }
    c_PartialMeltBase* get_partial_melt_model()    const noexcept { return this->p_partial_melt.get(); }

    // Binary I/O
    void write_binary(std::ostream& out) const override {
        const auto     name_len = static_cast<uint32_t>(this->p_name.size());
        const auto     mat_len  = static_cast<uint32_t>(this->p_material_name.size());
        const uint64_t payload  =
            sizeof(double)   * 2 +           // p_radius, p_mass
            sizeof(uint32_t) + name_len +    // name length + bytes
            sizeof(int32_t)  +               // layer_index
            sizeof(double)   +               // radius_inner
            sizeof(uint32_t) + mat_len +     // material_name length + bytes
            sizeof(uint8_t)  +               // is_tidal
            sizeof(double)   +               // tidal_scale
            sizeof(uint8_t)  +               // tidal_scale_method
            sizeof(double)   * 4 +           // shear modulus, bulk modulus, shear viscosity, bulk viscosity
            sizeof(double)   * 6 +           // love_number k, h, l (each: re + im)
            sizeof(uint8_t)  * 3 +           // is_solid, is_static, is_incompressible
            optional_binary_flag_bytes() +         // material EOS model presence flag
            this->physics_models_presence_bytes(); // rheology + viscosity + partial-melt presence flags

        write_binary_header(out, static_cast<uint32_t>(BinaryClassID::PhysicsLayer), payload);

        // c_BaseLayer fields (same layout as c_BaseLayer::write_binary payload)
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
        out.write(reinterpret_cast<const char*>(&is_tidal_byte),       sizeof(uint8_t));
        out.write(reinterpret_cast<const char*>(&this->p_tidal_scale), sizeof(double));
        const uint8_t scale_method_byte = static_cast<uint8_t>(this->p_tidal_scale_method);
        out.write(reinterpret_cast<const char*>(&scale_method_byte),   sizeof(uint8_t));

        // c_PhysicsLayer scalar fields
        out.write(reinterpret_cast<const char*>(&this->p_shear_modulus_static),    sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_bulk_modulus_static),     sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_shear_viscosity_static), sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_bulk_viscosity_static),  sizeof(double));

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

        if (!out) {
            throw std::runtime_error("TidalPy: failed to write PhysicsLayer binary data");
        }

        this->write_eos_model_binary(out);
        this->write_physics_models_binary(out);
    }

    void read_binary(std::istream& in, bool force = false) override {
        // Read and validate the 20-byte TPYB header.
        c_TidalPyBaseClass::read_binary(in, force);

        // c_BaseLayer fields
        in.read(reinterpret_cast<char*>(&this->p_radius), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_mass),   sizeof(double));

        uint32_t name_len = 0;
        in.read(reinterpret_cast<char*>(&name_len), sizeof(uint32_t));
        this->p_name.resize(name_len);
        if (name_len > 0) { in.read(this->p_name.data(), name_len); }

        int32_t idx = 0;
        in.read(reinterpret_cast<char*>(&idx), sizeof(int32_t));
        this->p_layer_index = static_cast<int>(idx);

        in.read(reinterpret_cast<char*>(&this->p_radius_inner), sizeof(double));

        uint32_t mat_len = 0;
        in.read(reinterpret_cast<char*>(&mat_len), sizeof(uint32_t));
        this->p_material_name.resize(mat_len);
        if (mat_len > 0) { in.read(this->p_material_name.data(), mat_len); }

        uint8_t is_tidal_byte = 0;
        in.read(reinterpret_cast<char*>(&is_tidal_byte), sizeof(uint8_t));
        this->p_is_tidal = static_cast<bool>(is_tidal_byte);

        in.read(reinterpret_cast<char*>(&this->p_tidal_scale), sizeof(double));

        uint8_t scale_method_byte = 0;
        in.read(reinterpret_cast<char*>(&scale_method_byte), sizeof(uint8_t));
        this->p_tidal_scale_method = static_cast<c_TidalScaleMethod>(scale_method_byte);

        // c_PhysicsLayer scalar fields
        in.read(reinterpret_cast<char*>(&this->p_shear_modulus_static),    sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_bulk_modulus_static),     sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_shear_viscosity_static), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_bulk_viscosity_static),  sizeof(double));

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

        if (!in) {
            throw std::runtime_error("TidalPy: failed to read PhysicsLayer binary data");
        }

        this->read_eos_model_binary(in, force);
        this->read_physics_models_binary(in, force);

        this->update_physicals();
    }

protected:
    // Recursive (de)serialization of the optional physics models (shear and bulk rheology, shear and bulk
    // viscosity, partial melt), shared by c_PhysicsLayer and its subclasses so the section keeps one byte
    // layout: a presence flag each, followed when set by the model's own record. On read the concrete model is
    // rebuilt through that module's binary-dispatch factory and re-registered as this layer's observer.
    void write_physics_models_binary(std::ostream& out) const {
        write_optional_binary(out, this->p_shear_rheology);
        write_optional_binary(out, this->p_bulk_rheology);
        write_optional_binary(out, this->p_shear_viscosity);
        write_optional_binary(out, this->p_bulk_viscosity);
        write_optional_binary(out, this->p_partial_melt);
    }

    void read_physics_models_binary(std::istream& in, bool force) {
        this->p_shear_rheology =
            read_optional_binary<c_RheologyBase>(in, force, c_rheology_from_binary);
        if (this->p_shear_rheology) { this->p_shear_rheology->set_layer_ptr(this); }
        this->p_bulk_rheology =
            read_optional_binary<c_RheologyBase>(in, force, c_rheology_from_binary);
        if (this->p_bulk_rheology) { this->p_bulk_rheology->set_layer_ptr(this); }
        this->p_shear_viscosity =
            read_optional_binary<c_ViscosityBase>(in, force, c_viscosity_from_binary);
        if (this->p_shear_viscosity) { this->p_shear_viscosity->set_layer_ptr(this); }
        this->p_bulk_viscosity =
            read_optional_binary<c_ViscosityBase>(in, force, c_viscosity_from_binary);
        if (this->p_bulk_viscosity) { this->p_bulk_viscosity->set_layer_ptr(this); }
        this->p_partial_melt =
            read_optional_binary<c_PartialMeltBase>(in, force, c_partial_melt_from_binary);
        if (this->p_partial_melt) { this->p_partial_melt->set_layer_ptr(this); }
    }

    // Payload bytes contributed by the five model presence flags (the nested
    // model records follow as separate appended records).
    static constexpr uint64_t physics_models_presence_bytes() {
        return 5 * optional_binary_flag_bytes();
    }

    double        p_shear_modulus_static = 0.0;   // [Pa]
    double        p_bulk_modulus_static  = 0.0;   // [Pa]
    double        p_shear_viscosity_static = TidalPyConstants::d_NAN;   // [Pa·s], NaN until set
    double        p_bulk_viscosity_static  = TidalPyConstants::d_NAN;   // [Pa·s], NaN until set
    c_LoveNumbers p_love_numbers;                       // k, h, l [dimensionless] placeholder
    // Radial-solver layer classification.
    bool          p_is_solid          = true;
    bool          p_is_static         = true;
    bool          p_is_incompressible = false;

    // Optional rheology objects (serialized recursively via write_physics_models_binary).
    std::unique_ptr<c_RheologyBase> p_shear_rheology;
    std::unique_ptr<c_RheologyBase> p_bulk_rheology;

    // Optional viscosity and partial-melt objects (serialized recursively): the pre-melt viscosities and the
    // melt weakening consumed by the world EOS solve.
    std::unique_ptr<c_ViscosityBase>   p_shear_viscosity;
    std::unique_ptr<c_ViscosityBase>   p_bulk_viscosity;
    std::unique_ptr<c_PartialMeltBase> p_partial_melt;
};

} // namespace tidalpy
