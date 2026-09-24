#pragma once
/*
 * base_.hpp: c_BaseLayer, the geometry base class for every TidalPy layer type (extends c_StructureBase).
 *
 * Holds the geometry and identification set at construction plus the EOS profile (density, gravity, pressure
 * against radius) populated by the world's EOS solve. Spatial fields in meters [m], mass in kilograms [kg].
 *
 * Binary format (20-byte header + variable payload):
 *   header: class_id = BinaryClassID::BaseLayer (100)
 *   payload layout (fixed part 47 bytes + variable string data):
 *     p_radius           (double, 8)
 *     p_mass             (double, 8)
 *     name_len           (uint32_t, 4)
 *     name               (name_len bytes, UTF-8)
 *     layer_index        (int32_t, 4)
 *     radius_inner     (double, 8)
 *     material_name_len  (uint32_t, 4)
 *     material_name      (material_name_len bytes, UTF-8)
 *     is_tidal           (uint8_t, 1)
 *     is_volume_fixed    (uint8_t, 1)
 *     tidal_scale        (double, 8; NaN when the world uses the layer's volume fraction)
 *     eos_model          presence flag (uint8_t, 1) + (if present) the model's own binary record
 *   Derived fields (thickness, volume, surface areas) are recomputed on load. The attached material EOS model is
 *   serialized; the EOS profile it produces is not and is repopulated by re-running the world EOS solve.
 */

#include <cctype>
#include <cstdint>
#include <limits>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "eos_data_.hpp"
#include "structure_base_.hpp"
#include "material_eos_.hpp"   // c_MaterialEOSBase (per-layer density source)

namespace tidalpy {

// The "class" key of a layer config table.
inline const char* c_layer_class_name(uint32_t class_id) noexcept {
    switch (class_id) {
        case static_cast<uint32_t>(BinaryClassID::PhysicsLayer):     return "physics";
        case static_cast<uint32_t>(BinaryClassID::SolidLiquidLayer): return "solidliquid";
        case static_cast<uint32_t>(BinaryClassID::GasLayer):         return "gas";
        case static_cast<uint32_t>(BinaryClassID::BaseLayer):
        default:                                                     return "base";
    }
}

// Grouped to avoid a long constructor argument list.
struct c_BaseLayerConfig {
    std::string        name;
    int                layer_index  = 0;
    double             radius_inner = 0.0;   // [m]
    double             radius_outer = 0.0;   // [m]
    double             mass         = 0.0;   // [kg]
    std::string        material_name = "Unknown";
    bool               is_tidal    = true;
    // False lets the layer grow or shrink to hold its mass while the solve redistributes the interior.
    bool               is_volume_fixed = true;
    // The layer's share of the planet's volume in the quasi-homogeneous Love methods; NaN takes the layer's volume
    // fraction when the world uses it (c_BaseLayer::calc_tidal_scale).
    double             tidal_scale = TidalPyConstants::d_NAN;   // dimensionless
};

class c_BaseLayer : public c_StructureBase {
public:
    c_BaseLayer() = default;

    explicit c_BaseLayer(const c_BaseLayerConfig& cfg)
        : c_StructureBase(cfg.radius_outer, cfg.mass),
          p_name(cfg.name),
          p_layer_index(cfg.layer_index),
          p_radius_inner(cfg.radius_inner),
          p_material_name(cfg.material_name),
          p_is_tidal(cfg.is_tidal),
          p_is_volume_fixed(cfg.is_volume_fixed),
          p_tidal_scale(cfg.tidal_scale)
    {
        // An inverted or negative shell would carry a negative volume through every mass and heating sum.
        if (!(std::isfinite(cfg.radius_inner) && std::isfinite(cfg.radius_outer)
                && (cfg.radius_inner >= 0.0) && (cfg.radius_outer >= cfg.radius_inner))) {
            throw std::invalid_argument(
                "TidalPy: layer '" + cfg.name + "' needs finite radii with 0 <= radius_inner <= radius_outer; got " +
                std::to_string(cfg.radius_inner) + " to " + std::to_string(cfg.radius_outer) + " m.");
        }
        if (cfg.mass < 0.0) {
            throw std::invalid_argument(
                "TidalPy: layer '" + cfg.name + "' has a negative mass (" + std::to_string(cfg.mass) + " kg).");
        }
        this->update_physicals();
    }

    ~c_BaseLayer() override = default;

    // p_eos is a unique_ptr, which deletes the implicit copy assignment, and the subclass operator=s call
    // this one. Source temporaries always have a null p_eos, so resetting on copy is safe.
    c_BaseLayer& operator=(const c_BaseLayer& other) noexcept {
        if (this != &other) {
            c_StructureBase::operator=(other);
            this->p_name               = other.p_name;
            this->p_layer_index        = other.p_layer_index;
            this->p_radius_inner       = other.p_radius_inner;
            this->p_radius_outer       = other.p_radius_outer;
            this->p_thickness          = other.p_thickness;
            this->p_volume             = other.p_volume;
            this->p_surface_area_inner = other.p_surface_area_inner;
            this->p_surface_area_outer = other.p_surface_area_outer;
            this->p_material_name      = other.p_material_name;
            this->p_is_tidal           = other.p_is_tidal;
            this->p_is_volume_fixed    = other.p_is_volume_fixed;
            this->p_tidal_scale        = other.p_tidal_scale;
            this->p_tidal_heating      = other.p_tidal_heating;
            this->p_eos_data           = other.p_eos_data;
            this->p_eos.reset();
        }
        return *this;
    }
    c_BaseLayer& operator=(c_BaseLayer&&) noexcept = default;

    const std::string& get_name()                const noexcept { return this->p_name; }
    int                get_layer_index()         const noexcept { return this->p_layer_index; }
    double             get_radius_inner()        const noexcept { return this->p_radius_inner; }
    double             get_radius_outer()        const noexcept { return this->p_radius_outer; }
    double             get_thickness()           const noexcept { return this->p_thickness; }
    double             get_volume()              const noexcept { return this->p_volume; }
    double             get_surface_area_inner()  const noexcept { return this->p_surface_area_inner; }
    double             get_surface_area_outer()  const noexcept { return this->p_surface_area_outer; }
    const std::string& get_material_name()       const noexcept { return this->p_material_name; }
    bool               get_is_tidal()            const noexcept { return this->p_is_tidal; }
    bool               get_is_volume_fixed()     const noexcept { return this->p_is_volume_fixed; }
    void               set_is_volume_fixed(bool value) noexcept { this->p_is_volume_fixed = value; }

    // Keeps every derived geometric quantity in step. The EOS solve calls it when a layer below has grown
    // or shrunk, or when this layer is holding its mass rather than its volume.
    void set_radii(double radius_inner, double radius_outer) noexcept {
        this->p_radius_inner = radius_inner;
        this->p_radius       = radius_outer;
        this->update_physicals();
    }
    // The configured tidal scale; NaN when the layer takes its volume fraction.
    double get_tidal_scale() const noexcept { return this->p_tidal_scale; }
    void   set_tidal_scale(double tidal_scale) noexcept { this->p_tidal_scale = tidal_scale; }

    // The share this layer carries in the quasi-homogeneous Love methods (homogeneous, cpl, ctl), which treat each
    // tidal layer as a homogeneous planet of its own averaged material and scale that planet's Im(k) by this
    // factor: the configured tidal_scale, or the layer's volume over the planet's [m3] when none is set. Zero for
    // a layer that is not tidal.
    double calc_tidal_scale(double planet_volume) const noexcept {
        if (!this->p_is_tidal) { return 0.0; }
        if (std::isfinite(this->p_tidal_scale)) { return this->p_tidal_scale; }
        return (planet_volume > TidalPyConstants::d_EPS) ? this->get_volume() / planet_volume : 0.0;
    }

    // Matches the binary class id, so a caller holding a c_BaseLayer* can build the matching wrapper.
    virtual uint32_t get_layer_class_id() const noexcept {
        return static_cast<uint32_t>(BinaryClassID::BaseLayer);
    }
    // A transient result, not serialized: the heating [W] the world's last calc_tides put in this layer, and NaN
    // until then (see c_LayeredWorld::calc_tides for how each Love method distributes it).
    double get_tidal_heating()                   const noexcept { return this->p_tidal_heating; }
    void set_tidal_heating(double heating)   noexcept { this->p_tidal_heating = heating; }

    // Each successful world EOS solve sets this to the enclosed-mass gain across the layer; before one it is
    // whatever the layer was constructed with.
    void set_mass(double mass) noexcept { this->p_mass = mass; }

    // Bulk density [kg m-3] = mass / shell volume; NaN for a zero-volume layer.
    double get_density_bulk() const noexcept {
        if (this->p_volume <= TidalPyConstants::d_EPS) { return TidalPyConstants::d_NAN; }
        return this->p_mass / this->p_volume;
    }

    bool   get_eos_data_populated()             const noexcept { return this->p_eos_data.is_populated(); }
    double get_density(double radius)         const noexcept { return this->p_eos_data.get_density(radius); }
    double get_gravity(double radius)         const noexcept { return this->p_eos_data.get_gravity(radius); }
    double get_pressure(double radius)        const noexcept { return this->p_eos_data.get_pressure(radius); }
    void   update_eos_data(const c_LayerEOSData& data) { this->p_eos_data = data; }
    // Forget the solved profile, so nothing reads a structure that no longer describes this layer.
    void   clear_eos_data() { this->p_eos_data = c_LayerEOSData(); }

    // Read from the solved EOS: the material evaluated these as the structure was integrated, so nothing is
    // calculated here and every value is the one the solve used. They are the frequency-independent moduli,
    // viscosities, and melt fraction after the partial-melt model. NaN before a solve, and for a profile
    // supplied by hand, which carries density, gravity, and pressure alone.
    bool   get_viscoelastic_populated()         const noexcept { return this->p_eos_data.has_dense_eval(); }
    double get_shear_modulus(double radius)   const noexcept {
        return this->p_eos_value(radius, C_EOS_SHEAR_MODULUS_INDEX);
    }
    double get_bulk_modulus(double radius)    const noexcept {
        return this->p_eos_value(radius, C_EOS_BULK_MODULUS_INDEX);
    }
    double get_shear_viscosity(double radius) const noexcept {
        return this->p_eos_value(radius, C_EOS_SHEAR_VISCOSITY_INDEX);
    }
    double get_bulk_viscosity(double radius)  const noexcept {
        return this->p_eos_value(radius, C_EOS_BULK_VISCOSITY_INDEX);
    }
    double get_melt_fraction(double radius)   const noexcept {
        return this->p_eos_value(radius, C_EOS_MELT_FRACTION_INDEX);
    }
    double get_temperature_at(double radius)  const noexcept {
        return this->p_eos_value(radius, C_EOS_TEMPERATURE_INDEX);
    }

    // Every solved quantity at a radius in one dense evaluation, in the layout of eos_layout_.hpp.
    void get_eos_state(double radius, double* y_out) const noexcept { this->p_eos_data.evaluate(radius, y_out); }

    // The per-layer density source of the world-level EOS solve. Ownership transfers in.
    void set_eos(std::unique_ptr<c_MaterialEOSBase> eos) {
        this->p_eos = std::move(eos);
        if (this->p_eos) { this->p_eos->set_layer_ptr(this); }
    }

    // Non-owning; null when unset.
    c_MaterialEOSBase* get_eos()      const noexcept { return this->p_eos.get(); }
    bool               get_eos_set()  const noexcept { return this->p_eos != nullptr; }

    void update_physicals() {
        this->p_radius_outer       = this->p_radius;
        this->p_thickness          = this->p_radius_outer - this->p_radius_inner;
        this->p_volume             = this->calc_volume_shell(this->p_radius_outer, this->p_radius_inner);
        this->p_surface_area_inner = this->calc_surface_area(this->p_radius_inner);
        this->p_surface_area_outer = this->calc_surface_area(this->p_radius_outer);
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
            optional_binary_flag_bytes();    // material EOS model presence flag

        write_binary_header(out, static_cast<uint32_t>(BinaryClassID::BaseLayer), payload);

        out.write(reinterpret_cast<const char*>(&this->p_radius),       sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_mass),         sizeof(double));
        out.write(reinterpret_cast<const char*>(&name_len),       sizeof(uint32_t));
        if (name_len > 0) out.write(p_name.data(), name_len);

        const int32_t idx = static_cast<int32_t>(this->p_layer_index);
        out.write(reinterpret_cast<const char*>(&idx),            sizeof(int32_t));
        out.write(reinterpret_cast<const char*>(&this->p_radius_inner), sizeof(double));
        out.write(reinterpret_cast<const char*>(&mat_len),        sizeof(uint32_t));
        if (mat_len > 0) out.write(this->p_material_name.data(), mat_len);

        const uint8_t is_tidal = static_cast<uint8_t>(this->p_is_tidal);
        out.write(reinterpret_cast<const char*>(&is_tidal),         sizeof(uint8_t));
        const uint8_t is_volume_fixed = static_cast<uint8_t>(this->p_is_volume_fixed);
        out.write(reinterpret_cast<const char*>(&is_volume_fixed),  sizeof(uint8_t));
        out.write(reinterpret_cast<const char*>(&this->p_tidal_scale),  sizeof(double));

        if (!out) {
            throw std::runtime_error("TidalPy: failed to write BaseLayer binary data");
        }

        this->write_eos_model_binary(out);
    }

    void read_binary(std::istream& in, bool force = false) override {
        c_TidalPyBaseClass::read_binary(in, force);
        // A loaded layer carries no solved profile or heating until its world solves again.
        this->clear_eos_data();
        this->p_tidal_heating = TidalPyConstants::d_NAN;

        in.read(reinterpret_cast<char*>(&this->p_radius), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_mass),   sizeof(double));

        this->p_name = read_binary_string(in);

        int32_t idx = 0;
        in.read(reinterpret_cast<char*>(&idx), sizeof(int32_t));
        this->p_layer_index = static_cast<int>(idx);

        in.read(reinterpret_cast<char*>(&this->p_radius_inner), sizeof(double));

        this->p_material_name = read_binary_string(in);

        uint8_t is_tidal = 0;
        in.read(reinterpret_cast<char*>(&is_tidal), sizeof(uint8_t));
        this->p_is_tidal = static_cast<bool>(is_tidal);
        uint8_t is_volume_fixed = 0;
        in.read(reinterpret_cast<char*>(&is_volume_fixed), sizeof(uint8_t));
        this->p_is_volume_fixed = static_cast<bool>(is_volume_fixed);

        in.read(reinterpret_cast<char*>(&this->p_tidal_scale), sizeof(double));


        if (!in) {
            throw std::runtime_error("TidalPy: failed to read BaseLayer binary data");
        }

        this->read_eos_model_binary(in, force);

        this->update_physicals();
    }

protected:
    // One entry of the evaluation layout at a radius.
    double p_eos_value(double radius, std::size_t index) const noexcept {
        double state[C_EOS_DY_VALUES];
        this->p_eos_data.evaluate(radius, state);
        return state[index];
    }

    // Shared by every layer class so the section has one byte layout: a presence flag followed, when set, by
    // the model's own binary record. On read the concrete model is rebuilt through the binary-dispatch
    // factory and re-registered as this layer's observer.
    void write_eos_model_binary(std::ostream& out) const {
        write_optional_binary(out, this->p_eos);
    }

    void read_eos_model_binary(std::istream& in, bool force) {
        this->p_eos = read_optional_binary<c_MaterialEOSBase>(in, force, c_material_eos_from_binary);
        if (this->p_eos) { this->p_eos->set_layer_ptr(this); }
    }

    // Set at construction and never modified.
    std::string p_name;
    int         p_layer_index        = 0;
    double      p_radius_inner       = 0.0;   // [m]
    double      p_radius_outer       = 0.0;   // [m]
    double      p_thickness          = 0.0;   // [m]
    double      p_volume             = 0.0;   // [m^3]
    double      p_surface_area_inner = 0.0;   // [m^2]
    double      p_surface_area_outer = 0.0;   // [m^2]
    std::string p_material_name;
    bool               p_is_tidal           = true;
    bool               p_is_volume_fixed    = true;
    double             p_tidal_scale        = TidalPyConstants::d_NAN;   // dimensionless; NaN: volume fraction
    double             p_tidal_heating      = std::numeric_limits<double>::quiet_NaN();  // [W]; set by the world tidal solve

    // Populated by the world-level EOS solve; not serialized.
    c_LayerEOSData p_eos_data;

    // Attached from Python through set_eos and serialized with the layer binary record.
    std::unique_ptr<c_MaterialEOSBase> p_eos;
};

} // namespace tidalpy
