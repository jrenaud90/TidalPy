#pragma once
/*
 * base_.hpp — c_BaseLayer: geometry base class for all TidalPy layer types.
 *
 * Inherits c_StructureBase (Utilities_x/classes_x/structure_base_.hpp).
 * Stores immutable geometry and identification data set at construction, plus
 * mutable EOS profile data (density, gravity, pressure vs. radius) populated
 * by c_EOSHandler (Material_x/, Phase 8).
 *
 * All spatial fields are in meters [m]; mass in kilograms [kg] (MKS).
 *
 * Binary format (20-byte header + variable payload):
 *   header: class_id = BinaryClassID::BaseLayer (100)
 *   payload layout (fixed part 45 bytes + variable string data):
 *     p_radius          (double, 8)
 *     p_mass            (double, 8)
 *     name_len          (uint32_t, 4)
 *     name              (name_len bytes, UTF-8)
 *     layer_index       (int32_t, 4)
 *     radius_inner_m    (double, 8)
 *     material_name_len (uint32_t, 4)
 *     material_name     (material_name_len bytes, UTF-8)
 *     is_tidal          (uint8_t, 1)
 *     tidal_scale       (double, 8)
 *   Derived fields (thickness, volume, surface areas) are recomputed on load.
 *   EOS profile data is NOT serialized; it must be repopulated after loading.
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

// ---------------------------------------------------------------------------
// c_TidalScaleMethod — how a layer's share of the world's global tidal heating
// is determined (consumed by c_LayeredWorld::calc_tides).
//   user_provided   : use the layer's tidal_scale field directly.
//   volume_fraction : scale = layer volume / planet volume.
//   tidal_timescale : Maxwell-time bell curve vs the tidal forcing period
//                     (per-layer volume-averaged eta/mu); not yet wired.
// ---------------------------------------------------------------------------
enum class c_TidalScaleMethod : uint8_t {
    user_provided   = 0,
    volume_fraction = 1,
    tidal_timescale = 2
};

// Case-insensitive string -> c_TidalScaleMethod (alias-aware). Throws
// std::invalid_argument on an unknown name (surfaced as ValueError in Cython).
inline c_TidalScaleMethod c_tidal_scale_method_from_name(const std::string& name) {
    std::string key;
    key.reserve(name.size());
    for (char ch : name) {
        key.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(ch))));
    }
    if (key == "user_provided" || key == "user_provided_scale" || key == "user") {
        return c_TidalScaleMethod::user_provided;
    }
    if (key == "volume_fraction" || key == "volume_fraction_scale" || key == "volume") {
        return c_TidalScaleMethod::volume_fraction;
    }
    if (key == "tidal_timescale" || key == "tidal_timescale_scale"
            || key == "timescale" || key == "maxwell") {
        return c_TidalScaleMethod::tidal_timescale;
    }
    throw std::invalid_argument("TidalPy: unknown tidal_scale_method '" + name + "'");
}

// Canonical name for a c_TidalScaleMethod (round-trips through the factory).
inline const char* c_tidal_scale_method_name(c_TidalScaleMethod method) noexcept {
    switch (method) {
        case c_TidalScaleMethod::volume_fraction: return "volume_fraction_scale";
        case c_TidalScaleMethod::tidal_timescale: return "tidal_timescale_scale";
        case c_TidalScaleMethod::user_provided:
        default:                                  return "user_provided_scale";
    }
}

// ---------------------------------------------------------------------------
// c_BaseLayerConfig — construction parameters for c_BaseLayer.
// Using a config struct avoids a long constructor argument list.
// ---------------------------------------------------------------------------
struct c_BaseLayerConfig {
    std::string        name;
    int                layer_index        = 0;
    double             radius_inner_m     = 0.0;   // [m]
    double             radius_outer_m     = 0.0;   // [m]
    double             mass_kg            = 0.0;   // [kg]
    std::string        material_name      = "Unknown";
    bool               is_tidal           = true;
    double             tidal_scale        = 1.0;   // dimensionless
    c_TidalScaleMethod tidal_scale_method = c_TidalScaleMethod::user_provided;
};

// ---------------------------------------------------------------------------
// c_BaseLayer
// ---------------------------------------------------------------------------
class c_BaseLayer : public c_StructureBase {
public:
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------
    c_BaseLayer() = default;

    explicit c_BaseLayer(const c_BaseLayerConfig& cfg)
        : c_StructureBase(cfg.radius_outer_m, cfg.mass_kg),
          p_name(cfg.name),
          p_layer_index(cfg.layer_index),
          p_radius_inner(cfg.radius_inner_m),
          p_material_name(cfg.material_name),
          p_is_tidal(cfg.is_tidal),
          p_tidal_scale(cfg.tidal_scale),
          p_tidal_scale_method(cfg.tidal_scale_method)
    {
        // Update derived properties.
        this->update_physicals();
    }

    ~c_BaseLayer() override = default;

    // The owned EOS model (p_eos) is a unique_ptr, which deletes the implicit
    // copy-assignment. Subclass operator=s call c_BaseLayer::operator=, so an
    // explicit one is provided here. Source temporaries always have a null
    // p_eos, so resetting on copy is safe (mirrors the rheology/cooling pattern).
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
            this->p_tidal_scale        = other.p_tidal_scale;
            this->p_tidal_scale_method = other.p_tidal_scale_method;
            this->p_tidal_heating      = other.p_tidal_heating;
            this->p_eos_data           = other.p_eos_data;
            // EOS model pointer cannot be copied; source temporaries have null ptrs.
            this->p_eos.reset();
        }
        return *this;
    }
    c_BaseLayer& operator=(c_BaseLayer&&) noexcept = default;

    // -----------------------------------------------------------------------
    // Immutable geometry getters (all const, MKS)
    // -----------------------------------------------------------------------
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
    double             get_tidal_scale()         const noexcept { return this->p_tidal_scale; }
    c_TidalScaleMethod get_tidal_scale_method()  const noexcept { return this->p_tidal_scale_method; }

    // Concrete-type discriminator (matches the binary class id) so a caller holding a
    // c_BaseLayer* can build the matching wrapper. Subclasses override.
    virtual uint32_t get_layer_class_id() const noexcept {
        return static_cast<uint32_t>(BinaryClassID::BaseLayer);
    }
    void   set_tidal_scale_method(c_TidalScaleMethod method) noexcept { this->p_tidal_scale_method = method; }

    // -----------------------------------------------------------------------
    // Tidal heating [W] deposited in this layer (transient result; NOT serialized).
    // Populated by the world's global tidal solve (c_LayeredWorld::calc_tides) as
    // the world's total tidal heating scaled by this layer's contribution. Returns
    // NaN until a tidal solve has run.
    // -----------------------------------------------------------------------
    double get_tidal_heating()                   const noexcept { return this->p_tidal_heating; }
    void   set_tidal_heating(double heating_w)   noexcept { this->p_tidal_heating = heating_w; }

    // -----------------------------------------------------------------------
    // EOS profile (mutable; populated by EOSHandler in Phase 8)
    // -----------------------------------------------------------------------
    bool   get_eos_data_populated()             const noexcept { return this->p_eos_data.is_populated(); }
    double get_density(double radius_m)         const noexcept { return this->p_eos_data.get_density(radius_m); }
    double get_gravity(double radius_m)         const noexcept { return this->p_eos_data.get_gravity(radius_m); }
    double get_pressure(double radius_m)        const noexcept { return this->p_eos_data.get_pressure(radius_m); }
    void   update_eos_data(const c_LayerEOSData& data) { this->p_eos_data = data; }

    // -----------------------------------------------------------------------
    // Viscoelastic state (post-melt by default; pre-melt via the premelt getters).
    // Populated by the world EOS solve's viscoelastic post-pass (PhysicsLayer and
    // its subclasses); NaN on a geometry-only BaseLayer or before the solve.
    // -----------------------------------------------------------------------
    bool   get_viscoelastic_populated()         const noexcept { return this->p_eos_data.is_viscoelastic_populated(); }
    double get_shear_modulus(double radius_m)   const noexcept { return this->p_eos_data.get_shear_modulus(radius_m); }
    double get_bulk_modulus(double radius_m)    const noexcept { return this->p_eos_data.get_bulk_modulus(radius_m); }
    double get_shear_viscosity(double radius_m) const noexcept { return this->p_eos_data.get_shear_viscosity(radius_m); }
    double get_bulk_viscosity(double radius_m)  const noexcept { return this->p_eos_data.get_bulk_viscosity(radius_m); }
    double get_premelt_shear_modulus(double radius_m)   const noexcept {
        return this->p_eos_data.get_premelt_shear_modulus(radius_m);
    }
    double get_premelt_bulk_modulus(double radius_m)    const noexcept {
        return this->p_eos_data.get_premelt_bulk_modulus(radius_m);
    }
    double get_premelt_shear_viscosity(double radius_m) const noexcept {
        return this->p_eos_data.get_premelt_shear_viscosity(radius_m);
    }
    double get_premelt_bulk_viscosity(double radius_m)  const noexcept {
        return this->p_eos_data.get_premelt_bulk_viscosity(radius_m);
    }

    // Store the pre/post-melt viscoelastic profiles (called by the world solve).
    void update_viscoelastic_data(
            const std::vector<double>& premelt_shear_pa,
            const std::vector<double>& premelt_bulk_pa,
            const std::vector<double>& premelt_shear_visc_pas,
            const std::vector<double>& premelt_bulk_visc_pas,
            const std::vector<double>& postmelt_shear_pa,
            const std::vector<double>& postmelt_bulk_pa,
            const std::vector<double>& postmelt_shear_visc_pas,
            const std::vector<double>& postmelt_bulk_visc_pas) {
        this->p_eos_data.populate_viscoelastic(
            premelt_shear_pa, premelt_bulk_pa, premelt_shear_visc_pas, premelt_bulk_visc_pas,
            postmelt_shear_pa, postmelt_bulk_pa, postmelt_shear_visc_pas, postmelt_bulk_visc_pas);
    }

    // -----------------------------------------------------------------------
    // Material EOS model (the per-layer density source consumed by the
    // whole-planet world-level EOS solve). Ownership is transferred in.
    // -----------------------------------------------------------------------
    void set_eos(std::unique_ptr<c_MaterialEOSBase> eos) {
        this->p_eos = std::move(eos);
        if (this->p_eos) { this->p_eos->set_layer_ptr(this); }
    }

    // Non-owning observer pointer to the attached EOS model (nullptr if unset).
    c_MaterialEOSBase* get_eos()      const noexcept { return this->p_eos.get(); }
    bool               get_eos_set()  const noexcept { return this->p_eos != nullptr; }

    // -----------------------------------------------------------------------
    // Update physical properties
    // -----------------------------------------------------------------------
    void update_physicals() {
        this->p_radius_outer       = this->p_radius;
        this->p_thickness          = this->p_radius_outer - this->p_radius_inner;
        this->p_volume             = this->calc_volume_shell(this->p_radius_outer, this->p_radius_inner);
        this->p_surface_area_inner = this->calc_surface_area(this->p_radius_inner);
        this->p_surface_area_outer = this->calc_surface_area(this->p_radius_outer);
    }

    // -----------------------------------------------------------------------
    // Binary I/O
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
            sizeof(double);                  // tidal_scale

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
        out.write(reinterpret_cast<const char*>(&this->p_tidal_scale),  sizeof(double));

        if (!out) {
            throw std::runtime_error("TidalPy: failed to write BaseLayer binary data");
        }
    }

    void read_binary(std::istream& in, bool force = false) override {
        // Reads and validates the 20-byte header; throws on version mismatch.
        c_TidalPyBaseClass::read_binary(in, force);

        in.read(reinterpret_cast<char*>(&this->p_radius), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_mass),   sizeof(double));

        uint32_t name_len = 0;
        in.read(reinterpret_cast<char*>(&name_len), sizeof(uint32_t));
        this->p_name.resize(name_len);
        if (name_len > 0) in.read(this->p_name.data(), name_len);

        int32_t idx = 0;
        in.read(reinterpret_cast<char*>(&idx), sizeof(int32_t));
        this->p_layer_index = static_cast<int>(idx);

        in.read(reinterpret_cast<char*>(&this->p_radius_inner), sizeof(double));

        uint32_t mat_len = 0;
        in.read(reinterpret_cast<char*>(&mat_len), sizeof(uint32_t));
        this->p_material_name.resize(mat_len);
        if (mat_len > 0) in.read(this->p_material_name.data(), mat_len);

        uint8_t is_tidal = 0;
        in.read(reinterpret_cast<char*>(&is_tidal), sizeof(uint8_t));
        this->p_is_tidal = static_cast<bool>(is_tidal);

        in.read(reinterpret_cast<char*>(&this->p_tidal_scale), sizeof(double));

        if (!in) {
            throw std::runtime_error("TidalPy: failed to read BaseLayer binary data");
        }

        // Restore derived fields from the loaded radii.
        this->update_physicals();
    }

protected:
    // Immutable geometry (set at construction; not modified after)
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
    double             p_tidal_scale        = 1.0;   // dimensionless
    c_TidalScaleMethod p_tidal_scale_method = c_TidalScaleMethod::user_provided;
    double             p_tidal_heating      = std::numeric_limits<double>::quiet_NaN();  // [W]; set by the world tidal solve

    // Mutable EOS profile (populated by the world-level EOS solve; not serialized)
    c_LayerEOSData p_eos_data;

    // Optional material EOS model — the per-layer density source (not serialized;
    // attached from Python via set_eos, mirroring the rheology/cooling pattern).
    std::unique_ptr<c_MaterialEOSBase> p_eos;
};

} // namespace tidalpy
