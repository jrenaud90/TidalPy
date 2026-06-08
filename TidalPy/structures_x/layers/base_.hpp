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

#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string>

#include "eos_data_.hpp"
#include "structure_base_.hpp"

namespace tidalpy {

// ---------------------------------------------------------------------------
// c_BaseLayerConfig — construction parameters for c_BaseLayer.
// Using a config struct avoids a long constructor argument list.
// ---------------------------------------------------------------------------
struct c_BaseLayerConfig {
    std::string name;
    int         layer_index    = 0;
    double      radius_inner_m = 0.0;   // [m]
    double      radius_outer_m = 0.0;   // [m]
    double      mass_kg        = 0.0;   // [kg]
    std::string material_name  = "Unknown";
    bool        is_tidal       = true;
    double      tidal_scale    = 1.0;   // dimensionless
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
          p_tidal_scale(cfg.tidal_scale)
    {
        // Update derived properties.
        this->update_physicals();
    }

    ~c_BaseLayer() override = default;

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

    // -----------------------------------------------------------------------
    // EOS profile (mutable; populated by EOSHandler in Phase 8)
    // -----------------------------------------------------------------------
    bool   get_eos_data_populated()             const noexcept { return this->p_eos_data.is_populated(); }
    double get_density(double radius_m)         const noexcept { return this->p_eos_data.get_density(radius_m); }
    double get_gravity(double radius_m)         const noexcept { return this->p_eos_data.get_gravity(radius_m); }
    double get_pressure(double radius_m)        const noexcept { return this->p_eos_data.get_pressure(radius_m); }
    void   update_eos_data(const c_LayerEOSData& data) { this->p_eos_data = data; }

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
    bool        p_is_tidal           = true;
    double      p_tidal_scale        = 1.0;   // dimensionless

    // Mutable EOS profile (populated by EOSHandler; not serialized)
    c_LayerEOSData p_eos_data;
};

} // namespace tidalpy
