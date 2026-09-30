#pragma once
/*
 * base_.hpp: c_BaseLayer, the layer class every TidalPy layer type builds on (extends c_StructureBase).
 *
 * Holds the geometry and identification set at construction, the radial-solver classification flags, the layer
 * temperature and thermal switches, the material (the layer's EOS model, which holds the static moduli, the shear
 * law, the viscosities, and partial melt), the shear and bulk rheology, and the EOS profile (density, gravity,
 * pressure against radius) populated by the world's EOS solve. A complex modulus is the rheology applied to the static
 * modulus and viscosity the solved EOS reports at a radius; without a rheology it is the static value as a purely real
 * number (no dissipation). All MKS. A layer a world owns reads and writes its profile under the world's call lock
 * (c_WorldCallLock, call_lock_.hpp), so its getters take turns with the world's solves on other threads, and tells
 * the world when a setting its EOS solve reads changes (c_LayerOwner), so the world forgets the structure it solved
 * with the old setting.
 *
 * Binary payload: radius and mass (c_StructureBase), the name, layer index, inner radius, material name, is_tidal,
 * is_volume_fixed, tidal_scale (NaN when the world uses the layer's volume fraction), the three classification flags,
 * the temperature, use_thermal_eos, and use_heating, then the material EOS model and the shear and bulk rheologies,
 * each behind a presence flag. The derived geometry is recomputed on load, and the EOS profile is not serialized:
 * re-run the world EOS solve after loading.
 */

#include <cctype>
#include <complex>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "call_lock_.hpp"      // c_WorldCallLock (the owning world's call lock, which a layer takes too)
#include "eos_data_.hpp"
#include "structure_base_.hpp"
#include "material_eos_.hpp"   // c_MaterialEOSBase (per-layer density source)
#include "rheology_.hpp"

namespace tidalpy {

// The world that owns a layer, as the layer sees it (c_BaseWorld). A layer tells its owner when a setting the
// owner's EOS solve reads changes (its material and the models the material holds, its temperature and thermal
// switches, its cooling and radiogenics, whether it holds its volume), so the owner forgets the structure it solved
// with the old setting rather than reporting it as the current one.
class c_LayerOwner {
public:
    virtual ~c_LayerOwner() = default;

    // Called by the changed layer while it holds the owner's call lock.
    virtual void update_after_layer_change() = 0;
};

// Non-owning reference to the world that owns a layer and to that world's call lock. It names the owner of the
// object it sits in, so neither a copy nor a move of a layer's contents carries it: a layer copied from a world's
// layer starts unowned, and a layer assigned into keeps its own owner.
class c_OwnerCallMutex {
public:
    c_OwnerCallMutex() noexcept = default;
    c_OwnerCallMutex(const c_OwnerCallMutex& /*other*/) noexcept {}
    c_OwnerCallMutex& operator=(const c_OwnerCallMutex& /*other*/) noexcept { return *this; }

    void set(c_LayerOwner* owner_ptr, std::recursive_mutex* mutex_ptr) noexcept {
        this->p_owner_ptr = owner_ptr;
        this->p_mutex_ptr = mutex_ptr;
    }
    std::recursive_mutex* get() const noexcept { return this->p_mutex_ptr; }
    c_LayerOwner* get_owner() const noexcept { return this->p_owner_ptr; }

private:
    c_LayerOwner*         p_owner_ptr = nullptr;
    std::recursive_mutex* p_mutex_ptr = nullptr;
};

// The "class" key of a layer config table.
inline const char* c_layer_class_name(uint32_t class_id) noexcept {
    switch (class_id) {
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
    // Radial-solver layer classification flags.
    bool               is_solid          = true;    // false for liquid layers
    bool               is_static         = true;    // use static (no dynamic terms) approximation
    bool               is_incompressible = false;   // use incompressible approximation
    // The layer temperature is 0 K until set: the cold, rigid limit of the viscosity laws.
    double             temperature     = 0.0;       // [K]
    bool               use_thermal_eos = false;     // the EOS density and bulk modulus see the temperature
    bool               use_heating     = false;     // the world's heat sources act inside this layer
};

// The world that owns layers; it reads their profile through the unlocked p_ helpers while it holds its call lock.
class c_BaseWorld;

class c_BaseLayer : public c_StructureBase {
    // The owning world reads the profile and applies the rheology through the p_ helpers while it holds its call
    // lock.
    friend class c_BaseWorld;

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
          p_tidal_scale(cfg.tidal_scale),
          p_is_solid(cfg.is_solid),
          p_is_static(cfg.is_static),
          p_is_incompressible(cfg.is_incompressible),
          p_temperature(cfg.temperature),
          p_use_thermal_eos(cfg.use_thermal_eos),
          p_use_heating(cfg.use_heating)
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

    // The owned models delete the implicit copy assignment, which Cython's stack allocation emits from freshly
    // constructed temporaries, and the subclass operator=s call this one. Source temporaries always hold null
    // models, so resetting them on copy is safe.
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
            this->p_is_solid           = other.p_is_solid;
            this->p_is_static          = other.p_is_static;
            this->p_is_incompressible  = other.p_is_incompressible;
            this->p_temperature        = other.p_temperature;
            this->p_use_thermal_eos    = other.p_use_thermal_eos;
            this->p_use_heating        = other.p_use_heating;
            this->p_eos_data           = other.p_eos_data;
            this->p_eos.reset();
            this->p_shear_rheology.reset();
            this->p_bulk_rheology.reset();
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
    // The EOS solve reads it, so the owning world forgets its solved structure (c_LayerOwner).
    void set_is_volume_fixed(bool value) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_is_volume_fixed = value;
        this->p_update_owner_after_change();
    }

    // Keeps every derived geometric quantity in step. The EOS solve calls it when a layer below has grown
    // or shrunk, or when this layer is holding its mass rather than its volume, so it leaves the owning world's
    // solved structure alone; a caller moving a world's layer by hand tells the world itself
    // (c_BaseWorld::update_after_layer_change). The owner's call lock keeps a move out of a running solve.
    void set_radii(double radius_inner, double radius_outer) noexcept {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
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

    // The binary class id, so a caller holding a c_BaseLayer* can build the matching wrapper.
    uint32_t get_layer_class_id() const noexcept { return this->get_binary_class_id(); }
    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::BaseLayer); }
    // A transient result, not serialized: the heating [W] the world's last calc_tides put in this layer, and NaN
    // until then (see c_BaseWorld::calc_tides for how each Love method distributes it).
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

    // The owning world and its call lock, set by the world as it takes the layer in; null for a layer no world owns.
    // Every read and write of the solved profile below holds the lock, so a profile read never overlaps a world call
    // that replaces the profile (solve_eos, load_binary), and so does every setter of what the world's solves read.
    // The world owns the mutex through a unique_ptr, so its address outlives any move of the world, and the world
    // outlives its layers.
    void set_owner(c_LayerOwner* owner_ptr, std::recursive_mutex* mutex_ptr) noexcept {
        this->p_owner_call_mutex.set(owner_ptr, mutex_ptr);
    }
    std::recursive_mutex* get_owner_call_mutex() const noexcept { return this->p_owner_call_mutex.get(); }

    // The profile getters are noexcept: a lock that cannot be taken (a broken process) terminates rather than
    // returning a value read without it.
    bool   get_eos_data_populated() const noexcept {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        return this->p_eos_data.is_populated();
    }
    void   update_eos_data(const c_LayerEOSData& data) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_eos_data = data;
    }
    // Forget the solved profile, so nothing reads a structure that no longer describes this layer.
    void   clear_eos_data() {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_eos_data = c_LayerEOSData();
    }

    // Read from the solved EOS: the material evaluated these as the structure was integrated, so nothing is
    // calculated here and every value is the one the solve used. They are the frequency-independent moduli,
    // viscosities, and melt fraction after the partial-melt model. NaN before a solve, and for a profile
    // supplied by hand, which carries density, gravity, and pressure alone. Any other entry of the layout is
    // read through get_eos_state or get_eos_fields.
    bool   get_viscoelastic_populated() const noexcept {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        return this->p_eos_data.has_dense_eval();
    }
    double get_shear_modulus(double radius)   const noexcept {
        return this->p_eos_value(radius, C_EOS_SHEAR_MODULUS_INDEX);
    }

    // Every solved quantity at a radius in one dense evaluation, in the layout of eos_layout_.hpp.
    void get_eos_state(double radius, double* y_out) const noexcept {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_eos_state(radius, y_out);
    }

    // Vectorized profile read: entries field_indices[0 .. num_fields) of the evaluation layout (eos_layout_.hpp) at
    // each of radii[0 .. num_radii) [m], from one dense evaluation per radius, written field-major:
    // values_out[field_i * num_radii + radius_i]. An index outside the layout gives NaN. The owner's lock is taken
    // once for the whole call, so every value comes from one solve and a loop over radii pays for one lock.
    //
    // Assumes values_out holds num_fields * num_radii doubles.
    void get_eos_fields(
            const std::size_t* field_indices,
            std::size_t num_fields,
            const double* radii,
            std::size_t num_radii,
            double* values_out) const noexcept {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        double state[C_EOS_DY_VALUES];
        for (std::size_t radius_i = 0; radius_i < num_radii; ++radius_i) {
            this->p_eos_state(radii[radius_i], state);
            for (std::size_t field_i = 0; field_i < num_fields; ++field_i) {
                const std::size_t field_index = field_indices[field_i];
                values_out[field_i * num_radii + radius_i] =
                    (field_index < C_EOS_DY_VALUES) ? state[field_index] : TidalPyConstants::d_NAN;
            }
        }
    }

    // The per-layer density source of the world-level EOS solve. Ownership transfers in. The viscosity and
    // partial-melt models the layer's previous material held carry over to a new material that has none of its own.
    // The owning world forgets the structure it solved with the old material (c_LayerOwner).
    void set_eos(std::unique_ptr<c_MaterialEOSBase> eos) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        if (eos && this->p_eos) { eos->adopt_missing_models(*this->p_eos); }
        this->p_eos = std::move(eos);
        if (this->p_eos) { this->p_eos->set_layer_ptr(this); }
        this->p_update_owner_after_change();
    }

    // Non-owning; null when unset.
    c_MaterialEOSBase* get_eos()      const noexcept { return this->p_eos.get(); }
    bool               get_eos_set()  const noexcept { return this->p_eos != nullptr; }

    bool get_is_solid()          const noexcept { return this->p_is_solid; }
    bool get_is_static()         const noexcept { return this->p_is_static; }
    bool get_is_incompressible() const noexcept { return this->p_is_incompressible; }

    // These control the shooting and propagation-matrix assumptions. Each Love solve reads them afresh, so they
    // leave the owning world's EOS solve standing; they take its call lock so they never change under a solve.
    void set_is_solid(bool value) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_is_solid = value;
    }
    void set_is_static(bool value) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_is_static = value;
    }
    void set_is_incompressible(bool value) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_is_incompressible = value;
    }

    // Layer temperature [K] and whether the EOS density law sees it. The EOS solve reads both, so setting either
    // makes the owning world forget its solved structure (c_LayerOwner).
    double get_temperature()     const noexcept { return this->p_temperature; }
    bool   get_use_thermal_eos() const noexcept { return this->p_use_thermal_eos; }
    void set_temperature(double value) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_temperature = value;
        this->p_update_owner_after_change();
    }
    void set_use_thermal_eos(bool value) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_use_thermal_eos = value;
        this->p_update_owner_after_change();
    }

    // Whether the world's heat sources act inside this layer during a thermal EOS solve. Off, the layer
    // generates no heat whatever models it carries. Setting it makes the owning world forget its solved structure.
    bool get_use_heating() const noexcept { return this->p_use_heating; }
    void set_use_heating(bool value) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_use_heating = value;
        this->p_update_owner_after_change();
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

    // The rheology applied to the static modulus and viscosity the solved EOS reports at that radius. The solved
    // profile is read under the owning world's call lock (set_owner).
    std::complex<double> calc_complex_shear_modulus(double radius, double frequency) const noexcept {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        return this->p_complex_modulus(true, radius, frequency);
    }
    std::complex<double> calc_complex_bulk_modulus(double radius, double frequency) const noexcept {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        return this->p_complex_modulus(false, radius, frequency);
    }

    // Vectorized form of the two above at one frequency [rad s-1]: the shear (is_shear) or bulk complex modulus [Pa]
    // at each of radii[0 .. num_radii) [m], taking the owner's lock once for the whole call.
    //
    // Assumes moduli_out holds num_radii values.
    void calc_complex_moduli(
            bool is_shear,
            const double* radii,
            std::size_t num_radii,
            double frequency,
            std::complex<double>* moduli_out) const noexcept {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        for (std::size_t radius_i = 0; radius_i < num_radii; ++radius_i) {
            moduli_out[radius_i] = this->p_complex_modulus(is_shear, radii[radius_i], frequency);
        }
    }

    // Ownership transfers in, and each registers this layer as the model's observer. The EOS solve does not read
    // the rheology (each Love solve applies it afresh), so the owning world's solved structure stands; the owner's
    // call lock keeps the swap out of a solve that is applying the old model.
    void set_shear_rheology(std::unique_ptr<c_RheologyBase> shear) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_shear_rheology = std::move(shear);
        if (this->p_shear_rheology) { this->p_shear_rheology->set_layer_ptr(this); }
    }

    void set_bulk_rheology(std::unique_ptr<c_RheologyBase> bulk) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
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
    // can be configured in one place. Attach the EOS first, or there is no material to give the model to. The EOS
    // solve evaluates them, so each makes the owning world forget its solved structure (c_LayerOwner).
    void set_shear_viscosity(std::unique_ptr<c_ViscosityBase> viscosity) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_require_eos("a shear viscosity model")->set_shear_viscosity(std::move(viscosity));
        this->p_update_owner_after_change();
    }
    void set_bulk_viscosity(std::unique_ptr<c_ViscosityBase> viscosity) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_require_eos("a bulk viscosity model")->set_bulk_viscosity(std::move(viscosity));
        this->p_update_owner_after_change();
    }
    void set_partial_melt(std::unique_ptr<c_PartialMeltBase> partial_melt) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_require_eos("a partial-melt model")->set_partial_melt(std::move(partial_melt));
        this->p_update_owner_after_change();
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

    void update_physicals() {
        this->p_radius_outer       = this->p_radius;
        this->p_thickness          = this->p_radius_outer - this->p_radius_inner;
        this->p_volume             = this->calc_volume_shell(this->p_radius_outer, this->p_radius_inner);
        this->p_surface_area_inner = this->calc_surface_area(this->p_radius_inner);
        this->p_surface_area_outer = this->calc_surface_area(this->p_radius_outer);
    }

    // What load_binary reads a file into first, so a bad file never reaches this layer
    // (c_TidalPyBaseClass::make_binary_scratch). Each layer class overrides it with its own.
    std::unique_ptr<c_TidalPyBaseClass> make_binary_scratch() const override {
        return std::make_unique<c_BaseLayer>();
    }

protected:
    void p_write_payload(std::ostream& out) const override {
        c_StructureBase::p_write_payload(out);
        write_binary_string(out, this->p_name);
        const int32_t layer_index = static_cast<int32_t>(this->p_layer_index);
        out.write(reinterpret_cast<const char*>(&layer_index),          sizeof(int32_t));
        out.write(reinterpret_cast<const char*>(&this->p_radius_inner), sizeof(double));
        write_binary_string(out, this->p_material_name);
        const uint8_t flag_bytes[2] = {
            static_cast<uint8_t>(this->p_is_tidal), static_cast<uint8_t>(this->p_is_volume_fixed)};
        out.write(reinterpret_cast<const char*>(flag_bytes), sizeof(flag_bytes));
        out.write(reinterpret_cast<const char*>(&this->p_tidal_scale), sizeof(double));
        const uint8_t state_bytes[5] = {
            static_cast<uint8_t>(this->p_is_solid), static_cast<uint8_t>(this->p_is_static),
            static_cast<uint8_t>(this->p_is_incompressible), static_cast<uint8_t>(this->p_use_thermal_eos),
            static_cast<uint8_t>(this->p_use_heating)};
        out.write(reinterpret_cast<const char*>(state_bytes), sizeof(state_bytes));
        out.write(reinterpret_cast<const char*>(&this->p_temperature), sizeof(double));
        write_optional_binary(out, this->p_eos);
        write_optional_binary(out, this->p_shear_rheology);
        write_optional_binary(out, this->p_bulk_rheology);
    }

    // A loaded layer carries no solved profile or heating until its world solves again.
    void p_read_payload(std::istream& in, bool force) override {
        this->clear_eos_data();
        this->p_tidal_heating = TidalPyConstants::d_NAN;
        c_StructureBase::p_read_payload(in, force);
        this->p_name = read_binary_string(in);
        int32_t layer_index = 0;
        in.read(reinterpret_cast<char*>(&layer_index),          sizeof(int32_t));
        in.read(reinterpret_cast<char*>(&this->p_radius_inner), sizeof(double));
        this->p_layer_index   = static_cast<int>(layer_index);
        this->p_material_name = read_binary_string(in);
        uint8_t flag_bytes[2] = {0, 0};
        in.read(reinterpret_cast<char*>(flag_bytes), sizeof(flag_bytes));
        this->p_is_tidal        = static_cast<bool>(flag_bytes[0]);
        this->p_is_volume_fixed = static_cast<bool>(flag_bytes[1]);
        in.read(reinterpret_cast<char*>(&this->p_tidal_scale), sizeof(double));
        uint8_t state_bytes[5] = {0, 0, 0, 0, 0};
        in.read(reinterpret_cast<char*>(state_bytes), sizeof(state_bytes));
        this->p_is_solid          = static_cast<bool>(state_bytes[0]);
        this->p_is_static         = static_cast<bool>(state_bytes[1]);
        this->p_is_incompressible = static_cast<bool>(state_bytes[2]);
        this->p_use_thermal_eos   = static_cast<bool>(state_bytes[3]);
        this->p_use_heating       = static_cast<bool>(state_bytes[4]);
        in.read(reinterpret_cast<char*>(&this->p_temperature), sizeof(double));
        this->p_eos = read_optional_binary<c_MaterialEOSBase>(in, force, c_material_eos_from_binary);
        if (this->p_eos) { this->p_eos->set_layer_ptr(this); }
        this->p_shear_rheology = read_optional_binary<c_RheologyBase>(in, force, c_rheology_from_binary);
        if (this->p_shear_rheology) { this->p_shear_rheology->set_layer_ptr(this); }
        this->p_bulk_rheology = read_optional_binary<c_RheologyBase>(in, force, c_rheology_from_binary);
        if (this->p_bulk_rheology) { this->p_bulk_rheology->set_layer_ptr(this); }
        this->update_physicals();
    }

    // The shear (is_shear) or bulk rheology applied to the solved static modulus and viscosity at a radius [m];
    // the caller holds the owner's lock.
    std::complex<double> p_complex_modulus(bool is_shear, double radius, double frequency) const noexcept {
        double state[C_EOS_DY_VALUES];
        this->p_eos_state(radius, state);
        if (is_shear) {
            return this->apply_shear_rheology(
                state[C_EOS_SHEAR_MODULUS_INDEX], state[C_EOS_SHEAR_VISCOSITY_INDEX], frequency);
        }
        return this->apply_bulk_rheology(
            state[C_EOS_BULK_MODULUS_INDEX], state[C_EOS_BULK_VISCOSITY_INDEX], frequency);
    }

    // The attached EOS model, or a clear error naming what needed it.
    c_MaterialEOSBase* p_require_eos(const char* what) const {
        if (!this->p_eos) {
            throw std::logic_error(
                std::string("TidalPy: attach an EOS model to layer '") + this->p_name + "' before giving it "
                + what + ": the material owns it.");
        }
        return this->p_eos.get();
    }

    // Every solved quantity at a radius; the caller holds the owner's lock. The owning world calls it directly
    // (friend) from inside its own locked calls, so an array read through the world takes the lock once.
    void p_eos_state(double radius, double* state_out) const noexcept {
        this->p_eos_data.evaluate(radius, state_out);
    }

    // Tells the owning world, if any, that a setting its EOS solve reads has changed; the caller holds the owner's
    // lock.
    void p_update_owner_after_change() {
        c_LayerOwner* owner_ptr = this->p_owner_call_mutex.get_owner();
        if (owner_ptr != nullptr) { owner_ptr->update_after_layer_change(); }
    }

    // One entry of the evaluation layout at a radius, under the owner's lock.
    double p_eos_value(double radius, std::size_t index) const noexcept {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        double state[C_EOS_DY_VALUES];
        this->p_eos_state(radius, state);
        return state[index];
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

    // Radial-solver layer classification.
    bool p_is_solid          = true;
    bool p_is_static         = true;
    bool p_is_incompressible = false;

    // Layer state (see c_BaseLayerConfig).
    double p_temperature     = 0.0;
    bool   p_use_thermal_eos = false;
    bool   p_use_heating     = false;

    // Populated by the world-level EOS solve; not serialized. Read and written under p_owner_call_mutex.
    c_LayerEOSData p_eos_data;

    // The owning world and its call lock (set_owner); not serialized, and never copied or moved with the layer.
    c_OwnerCallMutex p_owner_call_mutex;

    // Attached from Python through set_eos and serialized with the layer binary record.
    std::unique_ptr<c_MaterialEOSBase> p_eos;

    // Optional rheology objects. They are shared, not unique: a radial-solver solution exported to Python keeps a
    // copy of these pointers so it can reproduce the complex moduli it was solved with, at any radius, even after
    // this layer is gone.
    std::shared_ptr<c_RheologyBase> p_shear_rheology;
    std::shared_ptr<c_RheologyBase> p_bulk_rheology;
};

} // namespace tidalpy
