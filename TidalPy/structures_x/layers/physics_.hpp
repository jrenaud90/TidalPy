#pragma once
/*
 * physics_.hpp — c_PhysicsLayer: mechanical-properties layer.
 *
 * Inherits c_BaseLayer (structures_x/layers/base_.hpp).
 * Adds static mechanical properties (shear modulus, bulk modulus, shear and bulk
 * viscosity) and three complex Love numbers (k, h, l) held in a c_LoveNumbers
 * struct.  Optional rheology objects (c_RheologyBase subclasses) enable
 * frequency-dependent complex modulus calculations.
 *
 * When no rheology is attached, calc_complex_shear/bulk_modulus() return the
 * static value as a purely real complex number (no dissipation).
 *
 * All spatial fields are in MKS units.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::PhysicsLayer (101)
 *   payload:
 *     [all c_BaseLayer fields — same byte layout as BaseLayer binary payload]
 *     shear_modulus_static_pa      (double, 8)
 *     bulk_modulus_static_pa       (double, 8)
 *     shear_viscosity_static_pas   (double, 8)
 *     bulk_viscosity_static_pas    (double, 8)
 *     love_number_k  re, im        (double×2, 16)
 *     love_number_h  re, im        (double×2, 16)
 *     love_number_l  re, im        (double×2, 16)
 *     shear_rheology presence flag (uint8_t, 1) + (if present) its binary record
 *     bulk_rheology  presence flag (uint8_t, 1) + (if present) its binary record
 *   Attached rheology objects ARE serialized recursively (presence flag + the
 *   model's own binary record); the two presence flags are part of this payload,
 *   each nested model record follows as a separate record.
 *   EOS profile data is NOT serialized (inherited rule from c_BaseLayer).
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

// -------------------------------------------------------------------------------
// c_PhysicsConfig — construction parameters for c_PhysicsLayer.
// Extends c_BaseLayerConfig with mechanical property fields.
// -------------------------------------------------------------------------------
struct c_PhysicsConfig : public c_BaseLayerConfig {
    double        shear_modulus_static_pa    = 0.0;   // [Pa]
    double        bulk_modulus_static_pa     = 0.0;   // [Pa]
    double        shear_viscosity_static_pas = 0.0;   // [Pa·s]
    double        bulk_viscosity_static_pas  = 0.0;   // [Pa·s]
    c_LoveNumbers love_numbers;                       // k, h, l [dimensionless] placeholder
};

// -------------------------------------------------------------------------------
// c_PhysicsLayer
// -------------------------------------------------------------------------------
class c_PhysicsLayer : public c_BaseLayer {
public:
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------
    c_PhysicsLayer() = default;

    explicit c_PhysicsLayer(const c_PhysicsConfig& cfg)
        : c_BaseLayer(cfg),
          p_shear_modulus_static_pa(cfg.shear_modulus_static_pa),
          p_bulk_modulus_static_pa(cfg.bulk_modulus_static_pa),
          p_shear_viscosity_static_pas(cfg.shear_viscosity_static_pas),
          p_bulk_viscosity_static_pas(cfg.bulk_viscosity_static_pas),
          p_love_numbers(cfg.love_numbers)
    {}

    ~c_PhysicsLayer() override = default;

    // unique_ptr members delete the implicit copy-assignment. Cython's stack
    // allocation emits copy-assignment from freshly-constructed temporaries that
    // always have null rheology pointers, so resetting them on copy is safe.
    c_PhysicsLayer& operator=(const c_PhysicsLayer& other) noexcept {
        if (this != &other) {
            c_BaseLayer::operator=(other);
            this->p_shear_modulus_static_pa    = other.p_shear_modulus_static_pa;
            this->p_bulk_modulus_static_pa     = other.p_bulk_modulus_static_pa;
            this->p_shear_viscosity_static_pas = other.p_shear_viscosity_static_pas;
            this->p_bulk_viscosity_static_pas  = other.p_bulk_viscosity_static_pas;
            this->p_love_numbers               = other.p_love_numbers;
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

    // -----------------------------------------------------------------------
    // Static mechanical property getters (all const, MKS)
    // -----------------------------------------------------------------------
    double        get_shear_modulus_static()   const noexcept { return this->p_shear_modulus_static_pa; }
    double        get_bulk_modulus_static()    const noexcept { return this->p_bulk_modulus_static_pa; }
    double        get_shear_viscosity_static() const noexcept { return this->p_shear_viscosity_static_pas; }
    double        get_bulk_viscosity_static()  const noexcept { return this->p_bulk_viscosity_static_pas; }

    // Love number getters — full struct or individual components
    c_LoveNumbers        get_love_numbers()   const noexcept { return this->p_love_numbers; }
    std::complex<double> get_love_number_k()  const noexcept { return this->p_love_numbers.k; }
    std::complex<double> get_love_number_h()  const noexcept { return this->p_love_numbers.h; }
    std::complex<double> get_love_number_l()  const noexcept { return this->p_love_numbers.l; }

    // -----------------------------------------------------------------------
    // Tidal susceptibility [m^3]
    //
    // Purely geometric: (3/2) * r^5 / (G * m^2).
    // Uses the TidalPy global config pointer for Newton's G.
    // Returns 0.0 when the config pointer is null or mass is zero.
    // -----------------------------------------------------------------------
    double calc_tidal_susceptibility() const noexcept {
        if (tidalpy_config_ptr == nullptr || this->p_mass == 0.0) { return 0.0; }
        const double G = tidalpy_config_ptr->d_G;
        const double r = this->p_radius;
        const double m = this->p_mass;
        return (1.5 * r * r * r * r * r) / (G * m * m);
    }

    // -----------------------------------------------------------------------
    // Complex shear modulus [Pa] at forcing frequency frequency_rad_s.
    //
    // Delegates to p_shear_rheology->calc_complex_modulus when a rheology
    // object is set; otherwise returns the static shear modulus as a purely
    // real complex value.
    // -----------------------------------------------------------------------
    std::complex<double> calc_complex_shear_modulus(double frequency_rad_s) const noexcept {
        if (this->p_shear_rheology) {
            return this->p_shear_rheology->calc_complex_modulus(
                this->p_shear_modulus_static_pa, this->p_shear_viscosity_static_pas, frequency_rad_s);
        }
        return std::complex<double>(this->p_shear_modulus_static_pa, 0.0);
    }

    // -----------------------------------------------------------------------
    // Complex bulk modulus [Pa] at forcing frequency frequency_rad_s.
    //
    // Delegates to p_bulk_rheology->calc_complex_modulus when set; otherwise
    // returns the static bulk modulus as a purely real complex value.
    // -----------------------------------------------------------------------
    std::complex<double> calc_complex_bulk_modulus(double frequency_rad_s) const noexcept {
        if (this->p_bulk_rheology) {
            return this->p_bulk_rheology->calc_complex_modulus(
                this->p_bulk_modulus_static_pa, this->p_bulk_viscosity_static_pas, frequency_rad_s);
        }
        return std::complex<double>(this->p_bulk_modulus_static_pa, 0.0);
    }

    // -----------------------------------------------------------------------
    // Rheology setters (non-const; transfer ownership via unique_ptr)
    // Each setter stores the object and registers this layer as the observer.
    // -----------------------------------------------------------------------
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

    // -----------------------------------------------------------------------
    // Viscosity + partial-melt setters (non-const; transfer ownership).
    // The shear/bulk viscosity models supply the pre-melt viscosities at (T, P);
    // the partial-melt model weakens the static moduli and viscosities. These feed
    // the frequency-independent state computed by the world EOS solve. Each setter
    // registers this layer as the model's observer.
    // -----------------------------------------------------------------------
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

    // Non-owning observer pointers (nullptr if unset) — used by the world EOS solve
    // to compute the per-layer viscoelastic state.
    c_ViscosityBase*   get_shear_viscosity_model() const noexcept { return this->p_shear_viscosity.get(); }
    c_ViscosityBase*   get_bulk_viscosity_model()  const noexcept { return this->p_bulk_viscosity.get(); }
    c_PartialMeltBase* get_partial_melt_model()    const noexcept { return this->p_partial_melt.get(); }

    // -----------------------------------------------------------------------
    // Binary I/O
    // Writes a single record with class_id = PhysicsLayer (101).
    // All c_BaseLayer fields are written first (same byte layout as BaseLayer
    // binary payload), followed by 4 scalar doubles and 3 complex love numbers.
    // Rheology and EOS data are NOT serialized.
    // -----------------------------------------------------------------------
    void write_binary(std::ostream& out) const override {
        const auto     name_len = static_cast<uint32_t>(this->p_name.size());
        const auto     mat_len  = static_cast<uint32_t>(this->p_material_name.size());
        const uint64_t payload  =
            sizeof(double)   * 2 +           // p_radius, p_mass
            sizeof(uint32_t) + name_len +    // name length + bytes
            sizeof(int32_t)  +               // layer_index
            sizeof(double)   +               // radius_inner_m
            sizeof(uint32_t) + mat_len +     // material_name length + bytes
            sizeof(uint8_t)  +               // is_tidal
            sizeof(double)   +               // tidal_scale
            sizeof(double)   * 4 +           // shear modulus, bulk modulus, shear viscosity, bulk viscosity
            sizeof(double)   * 6 +           // love_number k, h, l (each: re + im)
            this->rheology_presence_bytes(); // shear + bulk rheology presence flags

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

        // c_PhysicsLayer scalar fields
        out.write(reinterpret_cast<const char*>(&this->p_shear_modulus_static_pa),    sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_bulk_modulus_static_pa),     sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_shear_viscosity_static_pas), sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_bulk_viscosity_static_pas),  sizeof(double));

        // Love numbers k, h, l
        auto write_complex = [&](const std::complex<double>& c) {
            const double re = c.real(), im = c.imag();
            out.write(reinterpret_cast<const char*>(&re), sizeof(double));
            out.write(reinterpret_cast<const char*>(&im), sizeof(double));
        };
        write_complex(this->p_love_numbers.k);
        write_complex(this->p_love_numbers.h);
        write_complex(this->p_love_numbers.l);

        if (!out) {
            throw std::runtime_error("TidalPy: failed to write PhysicsLayer binary data");
        }

        // Attached rheology models (presence flag + recursive record each).
        this->write_rheology_binary(out);
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

        // c_PhysicsLayer scalar fields
        in.read(reinterpret_cast<char*>(&this->p_shear_modulus_static_pa),    sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_bulk_modulus_static_pa),     sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_shear_viscosity_static_pas), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_bulk_viscosity_static_pas),  sizeof(double));

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

        if (!in) {
            throw std::runtime_error("TidalPy: failed to read PhysicsLayer binary data");
        }

        // Attached rheology models (presence flag + recursive record each).
        this->read_rheology_binary(in, force);

        // Recompute derived geometry fields from loaded radii.
        this->update_physicals();
    }

protected:
    // -----------------------------------------------------------------------
    // Recursive (de)serialization of the optional shear/bulk rheology models.
    //
    // Shared by c_PhysicsLayer and its subclasses (c_SolidLiquidLayer,
    // c_GasLayer) so the rheology section has one canonical byte layout. Each
    // model is written as a presence flag followed, when set, by the model's own
    // binary record; on read the correct concrete model is rebuilt via the
    // rheology binary-dispatch factory and re-registered as this layer's observer.
    // -----------------------------------------------------------------------
    void write_rheology_binary(std::ostream& out) const {
        write_optional_binary(out, this->p_shear_rheology);
        write_optional_binary(out, this->p_bulk_rheology);
    }

    void read_rheology_binary(std::istream& in, bool force) {
        this->p_shear_rheology =
            read_optional_binary<c_RheologyBase>(in, force, c_rheology_from_binary);
        if (this->p_shear_rheology) { this->p_shear_rheology->set_layer_ptr(this); }
        this->p_bulk_rheology =
            read_optional_binary<c_RheologyBase>(in, force, c_rheology_from_binary);
        if (this->p_bulk_rheology) { this->p_bulk_rheology->set_layer_ptr(this); }
    }

    // Payload bytes contributed by the two rheology presence flags (the nested
    // model records follow as separate appended records).
    static constexpr uint64_t rheology_presence_bytes() {
        return 2 * optional_binary_flag_bytes();
    }

    double        p_shear_modulus_static_pa    = 0.0;   // [Pa]
    double        p_bulk_modulus_static_pa     = 0.0;   // [Pa]
    double        p_shear_viscosity_static_pas = 0.0;   // [Pa·s]
    double        p_bulk_viscosity_static_pas  = 0.0;   // [Pa·s]
    c_LoveNumbers p_love_numbers;                       // k, h, l [dimensionless] placeholder

    // Optional rheology objects (not serialized; set by Python layer after construction).
    std::unique_ptr<c_RheologyBase> p_shear_rheology;
    std::unique_ptr<c_RheologyBase> p_bulk_rheology;

    // Optional viscosity + partial-melt objects (not yet serialized — TODO: add to
    // the recursive binary, like rheology). Supply the pre-melt viscosities and the
    // melt weakening consumed by the world EOS solve's state computation.
    std::unique_ptr<c_ViscosityBase>   p_shear_viscosity;
    std::unique_ptr<c_ViscosityBase>   p_bulk_viscosity;
    std::unique_ptr<c_PartialMeltBase> p_partial_melt;
};

} // namespace tidalpy
