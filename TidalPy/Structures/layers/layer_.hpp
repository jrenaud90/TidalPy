#pragma once
/*
 * layer_.hpp: c_Layer, the one layer class of a TidalPy world (extends c_StructureBase).
 *
 * A layer holds its geometry and identification, the material it is made of (a shared, immutable c_Material), the
 * physics switches that say how much of that material it uses, its temperature, the radial-solver assumptions, optional
 * shear and bulk rheologies that override its material's defaults, an optional cooling model (how heat moves through
 * it), an optional radiogenics model, and the profile its world's EOS solve populates (density, gravity, pressure,
 * and the material state against radius). All MKS.
 *
 * A layer a world owns reads and writes its profile under the world's call lock (c_WorldCallLock, call_lock_.hpp), so
 * its getters take turns with the world's solves on other threads, and tells the world when a setting its EOS solve
 * reads changes (c_LayerOwner), so the world forgets the structure it solved with the old setting.
 *
 * Binary payload: radius and mass (c_StructureBase), the name, layer index, inner radius, ten flag bytes (use_tides,
 * is_volume_fixed, state, is_static, is_incompressible, use_thermal_expansion, use_melting, use_pressure_melting,
 * use_melt_density, use_heating), tidal_scale (NaN when the world uses the layer's volume fraction), the temperature,
 * then the material, the shear and bulk rheology overrides, the cooling model, and the radiogenics model, each behind a
 * presence flag. The derived geometry is recomputed on load, and the EOS profile is not serialized: re-run the world
 * EOS solve after loading.
 */

#include <cctype>
#include <cmath>
#include <complex>
#include <cstdint>
#include <istream>
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
#include "../../Material/material_.hpp"
#include "../../Cooling/cooling_.hpp"
#include "../../Radiogenics/radiogenics_.hpp"

namespace tidalpy {

// The world that owns a layer, as the layer sees it (c_BaseWorld). A layer tells its owner when a setting the
// owner's EOS solve reads changes (its material, temperature, switches, cooling, radiogenics, or whether it holds its
// volume), so the owner forgets the structure it solved with the old setting rather than reporting it as the current
// one.
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

// How the radial solver treats a layer. Auto takes it from the material: a liquid-only material is liquid, any other
// is solid, split into liquid zones where it melts past its rigidity threshold (c_BaseWorld::get_zones).
enum class c_LayerState : uint8_t {
    Auto   = 0,
    Solid  = 1,
    Liquid = 2,
};

inline const char* c_layer_state_name(c_LayerState state) noexcept {
    switch (state) {
        case c_LayerState::Solid:  return "solid";
        case c_LayerState::Liquid: return "liquid";
        case c_LayerState::Auto:
        default:                   return "auto";
    }
}

// Case-insensitive; throws std::invalid_argument naming the accepted values.
inline c_LayerState c_layer_state_from_name(const std::string& name) {
    std::string lowered = name;
    for (char& character : lowered) {
        character = static_cast<char>(std::tolower(static_cast<unsigned char>(character)));
    }
    if (lowered == "auto")   { return c_LayerState::Auto; }
    if (lowered == "solid")  { return c_LayerState::Solid; }
    if (lowered == "liquid") { return c_LayerState::Liquid; }
    throw std::invalid_argument(
        "TidalPy: a layer state is 'auto', 'solid', or 'liquid'; got '" + name + "'.");
}

// Grouped to avoid a long constructor argument list.
struct c_LayerConfig {
    std::string name;
    int         layer_index  = 0;
    double      radius_inner = 0.0;   // [m]
    double      radius_outer = 0.0;   // [m]
    double      mass         = 0.0;   // [kg]
    // Whether the layer takes part in the tides: its tidal scale in the quasi-homogeneous Love methods.
    bool        use_tides    = true;
    // False lets the layer grow or shrink to hold its mass while the solve redistributes the interior.
    bool        is_volume_fixed = true;
    // The layer's share of the planet's volume in the quasi-homogeneous Love methods; NaN takes the layer's volume
    // fraction (c_Layer::calc_tidal_scale).
    double      tidal_scale = TidalPyConstants::d_NAN;   // dimensionless
    // Radial-solver assumptions.
    c_LayerState state             = c_LayerState::Auto;
    bool         is_static         = true;    // static (no inertia) equations in a liquid
    bool         is_incompressible = false;   // incompressible equations
    // The layer temperature is 0 K until set: the cold, rigid limit of the viscosity laws.
    double      temperature = 0.0;   // [K]
    // How much of its material the layer uses (c_MaterialSwitches), and whether the world's heat sources act in it.
    c_MaterialSwitches switches;
    bool               use_heating = false;
};

class c_BaseWorld;

class c_Layer : public c_StructureBase {
    // The owning world reads the profile and applies the rheology through the p_ helpers while it holds its call
    // lock.
    friend class c_BaseWorld;

public:
    c_Layer() = default;

    explicit c_Layer(const c_LayerConfig& cfg)
        : c_StructureBase(cfg.radius_outer, cfg.mass),
          p_name(cfg.name),
          p_layer_index(cfg.layer_index),
          p_radius_inner(cfg.radius_inner),
          p_use_tides(cfg.use_tides),
          p_is_volume_fixed(cfg.is_volume_fixed),
          p_tidal_scale(cfg.tidal_scale),
          p_state(cfg.state),
          p_is_static(cfg.is_static),
          p_is_incompressible(cfg.is_incompressible),
          p_temperature(cfg.temperature),
          p_switches(cfg.switches),
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

    ~c_Layer() override = default;

    // The owned radiogenics model deletes the implicit copy assignment, which Cython's stack allocation emits from
    // freshly constructed temporaries. Source temporaries hold no owned model, so it is reset on copy; the shared
    // material, rheologies, and cooling model are shared.
    c_Layer& operator=(const c_Layer& other) noexcept {
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
            this->p_use_tides          = other.p_use_tides;
            this->p_is_volume_fixed    = other.p_is_volume_fixed;
            this->p_tidal_scale        = other.p_tidal_scale;
            this->p_tidal_heating      = other.p_tidal_heating;
            this->p_state              = other.p_state;
            this->p_is_static          = other.p_is_static;
            this->p_is_incompressible  = other.p_is_incompressible;
            this->p_temperature        = other.p_temperature;
            this->p_switches           = other.p_switches;
            this->p_use_heating        = other.p_use_heating;
            this->p_eos_data           = other.p_eos_data;
            this->p_material           = other.p_material;
            this->p_shear_rheology     = other.p_shear_rheology;
            this->p_bulk_rheology      = other.p_bulk_rheology;
            this->p_cooling            = other.p_cooling;
            this->p_radiogenics.reset();
        }
        return *this;
    }
    c_Layer& operator=(c_Layer&&) noexcept = default;

    // =================================================================================================================
    // Geometry and identification
    // =================================================================================================================
    const std::string& get_name()               const noexcept { return this->p_name; }
    int                get_layer_index()        const noexcept { return this->p_layer_index; }
    double             get_radius_inner()       const noexcept { return this->p_radius_inner; }
    double             get_radius_outer()       const noexcept { return this->p_radius_outer; }
    double             get_thickness()          const noexcept { return this->p_thickness; }
    double             get_volume()             const noexcept { return this->p_volume; }
    double             get_surface_area_inner() const noexcept { return this->p_surface_area_inner; }
    double             get_surface_area_outer() const noexcept { return this->p_surface_area_outer; }
    double             get_radius_mid()         const noexcept {
        return 0.5 * (this->p_radius_inner + this->p_radius_outer);
    }

    bool get_is_volume_fixed() const noexcept { return this->p_is_volume_fixed; }
    // The EOS solve reads it, so the owning world forgets its solved structure (c_LayerOwner).
    void set_is_volume_fixed(bool value) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_is_volume_fixed = value;
        this->p_update_owner_after_change();
    }

    // Keeps every derived geometric quantity in step. The EOS solve calls it when a layer below has grown or shrunk,
    // or when this layer is holding its mass rather than its volume, so it leaves the owning world's solved structure
    // alone; a caller moving a world's layer by hand tells the world itself (c_BaseWorld::update_after_layer_change).
    // The owner's call lock keeps a move out of a running solve.
    void set_radii(double radius_inner, double radius_outer) noexcept {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_radius_inner = radius_inner;
        this->p_radius       = radius_outer;
        this->update_physicals();
    }

    void update_physicals() {
        this->p_radius_outer       = this->p_radius;
        this->p_thickness          = this->p_radius_outer - this->p_radius_inner;
        this->p_volume             = this->calc_volume_shell(this->p_radius_outer, this->p_radius_inner);
        this->p_surface_area_inner = this->calc_surface_area(this->p_radius_inner);
        this->p_surface_area_outer = this->calc_surface_area(this->p_radius_outer);
    }

    // Each successful world EOS solve sets this to the enclosed-mass gain across the layer; before one it is
    // whatever the layer was constructed with.
    void set_mass(double mass) noexcept { this->p_mass = mass; }

    // Bulk density [kg m-3] = mass / shell volume; NaN for a zero-volume layer.
    double get_density_bulk() const noexcept {
        if (this->p_volume <= TidalPyConstants::d_EPS) { return TidalPyConstants::d_NAN; }
        return this->p_mass / this->p_volume;
    }

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::Layer); }

    // What load_binary reads a file into first, so a bad file never reaches this layer
    // (c_TidalPyBaseClass::make_binary_scratch).
    std::unique_ptr<c_TidalPyBaseClass> make_binary_scratch() const override {
        return std::make_unique<c_Layer>();
    }

    // =================================================================================================================
    // Tides
    // =================================================================================================================
    bool get_use_tides() const noexcept { return this->p_use_tides; }
    void set_use_tides(bool value) noexcept { this->p_use_tides = value; }

    // The configured tidal scale; NaN when the layer takes its volume fraction.
    double get_tidal_scale() const noexcept { return this->p_tidal_scale; }
    void   set_tidal_scale(double tidal_scale) noexcept { this->p_tidal_scale = tidal_scale; }

    // The share this layer carries in the quasi-homogeneous Love methods (homogeneous, cpl, ctl), which treat each
    // tidal layer as a homogeneous planet of its own averaged material and scale that planet's Im(k) by this factor:
    // the configured tidal_scale, or the layer's volume over the planet's [m3] when none is set. Zero for a layer that
    // takes no part in the tides.
    double calc_tidal_scale(double planet_volume) const noexcept {
        if (!this->p_use_tides) { return 0.0; }
        if (std::isfinite(this->p_tidal_scale)) { return this->p_tidal_scale; }
        return (planet_volume > TidalPyConstants::d_EPS) ? this->get_volume() / planet_volume : 0.0;
    }

    // A transient result, not serialized: the heating [W] the world's last calc_tides put in this layer, and NaN
    // until then (see c_BaseWorld::calc_tides for how each Love method distributes it).
    double get_tidal_heating() const noexcept { return this->p_tidal_heating; }
    void set_tidal_heating(double heating) noexcept { this->p_tidal_heating = heating; }

    // =================================================================================================================
    // Material, switches, and state
    // =================================================================================================================
    // Shared and immutable, so one material serves any number of layers and solves. The owning world forgets the
    // structure it solved with the old material (c_LayerOwner).
    void set_material(std::shared_ptr<const c_Material> material) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_material = std::move(material);
        this->p_update_owner_after_change();
    }
    // Non-owning; null when unset.
    const c_Material* get_material() const noexcept { return this->p_material.get(); }
    std::shared_ptr<const c_Material> share_material() const noexcept { return this->p_material; }
    bool get_material_set() const noexcept { return this->p_material != nullptr; }

    const c_MaterialSwitches& get_switches() const noexcept { return this->p_switches; }
    // Every switch changes what the EOS solve evaluates, so each setter makes the owning world forget its solved
    // structure (c_LayerOwner).
    void set_switches(const c_MaterialSwitches& switches) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_switches = switches;
        this->p_update_owner_after_change();
    }

    // Layer temperature [K]: the profile of a solve without a temperature contrast, and the lumped temperature of a
    // thermal solve. The EOS solve reads it, so setting it makes the owning world forget its solved structure.
    double get_temperature() const noexcept { return this->p_temperature; }
    void set_temperature(double value) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_temperature = value;
        this->p_update_owner_after_change();
    }

    // Whether the world's heat sources act inside this layer during a thermal EOS solve. Off, the layer generates no
    // heat whatever models it carries. Setting it makes the owning world forget its solved structure.
    bool get_use_heating() const noexcept { return this->p_use_heating; }
    void set_use_heating(bool value) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_use_heating = value;
        this->p_update_owner_after_change();
    }

    // The layer's material at a point (pressure [Pa], temperature [K], radius [m]) with the layer's switches. A layer
    // without a material reports NaN throughout.
    void calc_state(const c_ThermoPoint& point, c_MaterialState& out) const noexcept {
        if (!this->p_material) {
            out = c_MaterialState();
            return;
        }
        this->p_material->calc_state(point, this->p_switches, out);
    }

    // The layer's material at zero pressure, its own temperature, and its mid-radius: the layer-constant state the
    // one-argument complex moduli use.
    void calc_reference_state(c_MaterialState& out) const noexcept {
        c_ThermoPoint point;
        point.pressure    = 0.0;
        point.temperature = this->p_temperature;
        point.radius      = this->get_radius_mid();
        this->calc_state(point, out);
    }

    // How the radial solver treats the layer.
    c_LayerState get_state() const noexcept { return this->p_state; }
    // The radial solver reads it afresh at every Love solve, so the EOS solve stands, unless the change decides whether
    // the layer can change state (get_can_change_state): the EOS solve finds its solid and liquid zones only where it
    // can, so then the owning world forgets its solved structure. Takes the owner's call lock so it never changes
    // under a solve.
    void set_state(c_LayerState state) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        const bool could_change_state = this->get_can_change_state();
        this->p_state = state;
        if (could_change_state != this->get_can_change_state()) { this->p_update_owner_after_change(); }
    }

    // Whether the radial solver treats the whole layer as a liquid: a layer set liquid, or an automatic one made of a
    // liquid-only material.
    bool get_is_liquid() const noexcept {
        if (this->p_state == c_LayerState::Liquid) { return true; }
        if (this->p_state == c_LayerState::Solid)  { return false; }
        return this->p_material && this->p_material->get_is_liquid_only();
    }

    // Whether part of the layer can turn liquid during a solve: an automatic layer whose material melts, with melting
    // switched on.
    bool get_can_change_state() const noexcept {
        return (this->p_state == c_LayerState::Auto) && this->p_material && this->p_material->get_can_melt()
            && this->p_switches.use_melting;
    }

    bool get_is_static()         const noexcept { return this->p_is_static; }
    bool get_is_incompressible() const noexcept { return this->p_is_incompressible; }
    // These control the shooting and propagation-matrix assumptions. Each Love solve reads them afresh, so they leave
    // the owning world's EOS solve standing; they take its call lock so they never change under a solve.
    void set_is_static(bool value) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_is_static = value;
    }
    void set_is_incompressible(bool value) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_is_incompressible = value;
    }

    // =================================================================================================================
    // Rheology
    // =================================================================================================================
    // Overrides of the material's default rheologies; null clears an override, so the material's (if any) applies.
    // The EOS solve does not read the rheology (each Love solve applies it afresh), so the owning world's solved
    // structure stands; the owner's call lock keeps the swap out of a solve that is applying the old model.
    void set_shear_rheology(std::shared_ptr<const c_RheologyBase> rheology) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_shear_rheology = std::move(rheology);
    }
    void set_bulk_rheology(std::shared_ptr<const c_RheologyBase> rheology) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_bulk_rheology = std::move(rheology);
    }
    // The layer's own overrides; null when unset.
    const c_RheologyBase* get_shear_rheology_override() const noexcept { return this->p_shear_rheology.get(); }
    const c_RheologyBase* get_bulk_rheology_override()  const noexcept { return this->p_bulk_rheology.get(); }

    // The rheology in effect: the layer's override, else its material's default (its base phase's); null for none, in
    // which case the complex modulus is the static modulus as a purely real number. Shared, so a radial-solver
    // solution exported to Python can reproduce the complex moduli it was solved with after this layer is gone.
    std::shared_ptr<const c_RheologyBase> share_shear_rheology() const noexcept {
        if (this->p_shear_rheology || !this->p_material) { return this->p_shear_rheology; }
        return this->p_material->get_base_phase().get_components().shear_rheology;
    }
    std::shared_ptr<const c_RheologyBase> share_bulk_rheology() const noexcept {
        if (this->p_bulk_rheology || !this->p_material) { return this->p_bulk_rheology; }
        return this->p_material->get_base_phase().get_components().bulk_rheology;
    }
    const c_RheologyBase* get_shear_rheology() const noexcept { return this->share_shear_rheology().get(); }
    const c_RheologyBase* get_bulk_rheology()  const noexcept { return this->share_bulk_rheology().get(); }

    // The same through the generic model handle the Python wrappers hold; each setter checks the model's family and
    // throws std::invalid_argument for the wrong one. A null handle clears.
    void set_material_model(const std::shared_ptr<c_PhysicsBase>& model) {
        this->set_material(c_share_as<c_Material>(model, "a layer's material"));
    }
    void set_shear_rheology_model(const std::shared_ptr<c_PhysicsBase>& model) {
        this->set_shear_rheology(c_share_as<c_RheologyBase>(model, "a layer's shear rheology"));
    }
    void set_bulk_rheology_model(const std::shared_ptr<c_PhysicsBase>& model) {
        this->set_bulk_rheology(c_share_as<c_RheologyBase>(model, "a layer's bulk rheology"));
    }
    std::shared_ptr<c_PhysicsBase> share_material_model() const { return c_share_physics_of(this->p_material); }
    // The rheology in effect (is_override false) or the layer's own override (true).
    std::shared_ptr<c_PhysicsBase> share_shear_rheology_model(bool is_override) const {
        return c_share_physics_of(is_override ? this->p_shear_rheology : this->share_shear_rheology());
    }
    std::shared_ptr<c_PhysicsBase> share_bulk_rheology_model(bool is_override) const {
        return c_share_physics_of(is_override ? this->p_bulk_rheology : this->share_bulk_rheology());
    }

    // The only place a complex modulus comes from: the material supplies the two static inputs and knows nothing
    // about frequency. Purely real, with no dissipation, when no rheology is in effect.
    std::complex<double> apply_shear_rheology(
            double static_modulus,
            double viscosity,
            double frequency) const noexcept {
        const c_RheologyBase* rheology = this->get_shear_rheology();
        if (rheology) { return rheology->calc_complex_modulus(static_modulus, viscosity, frequency); }
        return std::complex<double>(static_modulus, 0.0);
    }
    std::complex<double> apply_bulk_rheology(
            double static_modulus,
            double viscosity,
            double frequency) const noexcept {
        const c_RheologyBase* rheology = this->get_bulk_rheology();
        if (rheology) { return rheology->calc_complex_modulus(static_modulus, viscosity, frequency); }
        return std::complex<double>(static_modulus, 0.0);
    }

    // The rheology applied to the layer-constant state (calc_reference_state); NaN without a material.
    std::complex<double> calc_complex_shear_modulus(double frequency) const noexcept {
        c_MaterialState state;
        this->calc_reference_state(state);
        return this->apply_shear_rheology(state.shear_modulus, state.shear_viscosity, frequency);
    }
    std::complex<double> calc_complex_bulk_modulus(double frequency) const noexcept {
        c_MaterialState state;
        this->calc_reference_state(state);
        return this->apply_bulk_rheology(state.adiabatic_bulk_modulus, state.bulk_viscosity, frequency);
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

    // =================================================================================================================
    // Cooling and radiogenics
    // =================================================================================================================
    // The cooling model is shared and immutable, like the rheologies; null clears it, which leaves the layer at one
    // temperature. The radiogenics model's ownership transfers in, and it registers this layer as its observer. A
    // thermal EOS solve reads both, so each makes the owning world forget its solved structure (c_LayerOwner).
    void set_cooling(std::shared_ptr<const c_CoolingBase> cooling) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_cooling = std::move(cooling);
        this->p_update_owner_after_change();
    }
    // Through the generic model handle the Python wrappers hold; throws std::invalid_argument for another family.
    void set_cooling_model(const std::shared_ptr<c_PhysicsBase>& model) {
        this->set_cooling(c_share_as<c_CoolingBase>(model, "a layer's cooling model"));
    }
    std::shared_ptr<c_PhysicsBase> share_cooling_model() const { return c_share_physics_of(this->p_cooling); }
    void set_radiogenics(std::unique_ptr<c_RadiogenicsBase> radiogenics) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_radiogenics = std::move(radiogenics);
        if (this->p_radiogenics) { this->p_radiogenics->set_layer_ptr(this); }
        this->p_update_owner_after_change();
    }
    // Non-owning; null when unset.
    const c_CoolingBase*     get_cooling_model()     const noexcept { return this->p_cooling.get(); }
    const c_RadiogenicsBase* get_radiogenics_model() const noexcept { return this->p_radiogenics.get(); }

    // Radiogenic heating [W] from the attached model at a time [s] for a mass [kg]; 0.0 when none is attached.
    double calc_radiogenic_heating(double time, double mass) const noexcept {
        if (!this->p_radiogenics) { return 0.0; }
        return this->p_radiogenics->calc_heating(time, mass);
    }

    // =================================================================================================================
    // Solved profile
    // =================================================================================================================
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
    bool get_eos_data_populated() const noexcept {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        return this->p_eos_data.is_populated();
    }
    void update_eos_data(const c_LayerEOSData& data) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_eos_data = data;
    }
    // Forget the solved profile, so nothing reads a structure that no longer describes this layer.
    void clear_eos_data() {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_eos_data = c_LayerEOSData();
    }

    // Read from the solved EOS: the material was evaluated there as the structure was integrated, so nothing is
    // calculated here and every value is the one the solve used. NaN before a solve, and for a profile supplied by
    // hand, which carries density, gravity, and pressure alone.
    bool get_viscoelastic_populated() const noexcept {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        return this->p_eos_data.has_dense_eval();
    }
    double get_shear_modulus(double radius) const noexcept {
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

protected:
    void p_write_payload(std::ostream& out) const override {
        c_StructureBase::p_write_payload(out);
        write_binary_string(out, this->p_name);
        const int32_t layer_index = static_cast<int32_t>(this->p_layer_index);
        out.write(reinterpret_cast<const char*>(&layer_index),          sizeof(int32_t));
        out.write(reinterpret_cast<const char*>(&this->p_radius_inner), sizeof(double));
        const uint8_t flag_bytes[10] = {
            static_cast<uint8_t>(this->p_use_tides),
            static_cast<uint8_t>(this->p_is_volume_fixed),
            static_cast<uint8_t>(this->p_state),
            static_cast<uint8_t>(this->p_is_static),
            static_cast<uint8_t>(this->p_is_incompressible),
            static_cast<uint8_t>(this->p_switches.use_thermal_expansion),
            static_cast<uint8_t>(this->p_switches.use_melting),
            static_cast<uint8_t>(this->p_switches.use_pressure_melting),
            static_cast<uint8_t>(this->p_switches.use_melt_density),
            static_cast<uint8_t>(this->p_use_heating)};
        out.write(reinterpret_cast<const char*>(flag_bytes), sizeof(flag_bytes));
        out.write(reinterpret_cast<const char*>(&this->p_tidal_scale), sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_temperature), sizeof(double));
        write_optional_binary(out, this->p_material);
        write_optional_binary(out, this->p_shear_rheology);
        write_optional_binary(out, this->p_bulk_rheology);
        write_optional_binary(out, this->p_cooling);
        write_optional_binary(out, this->p_radiogenics);
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
        this->p_layer_index = static_cast<int>(layer_index);
        uint8_t flag_bytes[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        in.read(reinterpret_cast<char*>(flag_bytes), sizeof(flag_bytes));
        if (flag_bytes[2] > static_cast<uint8_t>(c_LayerState::Liquid)) {
            throw std::runtime_error("TidalPy: corrupt layer binary data: unknown layer state.");
        }
        this->p_use_tides                      = static_cast<bool>(flag_bytes[0]);
        this->p_is_volume_fixed                = static_cast<bool>(flag_bytes[1]);
        this->p_state                          = static_cast<c_LayerState>(flag_bytes[2]);
        this->p_is_static                      = static_cast<bool>(flag_bytes[3]);
        this->p_is_incompressible              = static_cast<bool>(flag_bytes[4]);
        this->p_switches.use_thermal_expansion = static_cast<bool>(flag_bytes[5]);
        this->p_switches.use_melting           = static_cast<bool>(flag_bytes[6]);
        this->p_switches.use_pressure_melting  = static_cast<bool>(flag_bytes[7]);
        this->p_switches.use_melt_density      = static_cast<bool>(flag_bytes[8]);
        this->p_use_heating                    = static_cast<bool>(flag_bytes[9]);
        in.read(reinterpret_cast<char*>(&this->p_tidal_scale), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_temperature), sizeof(double));
        this->p_material       = read_optional_binary<c_Material>(in, force, c_material_from_binary);
        this->p_shear_rheology = read_optional_binary<c_RheologyBase>(in, force, c_rheology_from_binary);
        this->p_bulk_rheology  = read_optional_binary<c_RheologyBase>(in, force, c_rheology_from_binary);
        this->p_cooling = read_optional_binary<c_CoolingBase>(in, force, c_cooling_from_binary);
        this->p_radiogenics = read_optional_binary<c_RadiogenicsBase>(in, force, c_radiogenics_from_binary);
        if (this->p_radiogenics) { this->p_radiogenics->set_layer_ptr(this); }
        this->update_physicals();
    }

    // The shear (is_shear) or bulk rheology applied to the solved static modulus and viscosity at a radius [m]; the
    // caller holds the owner's lock.
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

    std::string p_name;
    int         p_layer_index        = 0;
    double      p_radius_inner       = 0.0;   // [m]
    double      p_radius_outer       = 0.0;   // [m]
    double      p_thickness          = 0.0;   // [m]
    double      p_volume             = 0.0;   // [m^3]
    double      p_surface_area_inner = 0.0;   // [m^2]
    double      p_surface_area_outer = 0.0;   // [m^2]
    bool        p_use_tides          = true;
    bool        p_is_volume_fixed    = true;
    double      p_tidal_scale        = TidalPyConstants::d_NAN;   // dimensionless; NaN: volume fraction
    double      p_tidal_heating      = std::numeric_limits<double>::quiet_NaN();  // [W]; set by the world tidal solve

    // Radial-solver assumptions.
    c_LayerState p_state             = c_LayerState::Auto;
    bool         p_is_static         = true;
    bool         p_is_incompressible = false;

    // Layer state (see c_LayerConfig).
    double             p_temperature = 0.0;   // [K]
    c_MaterialSwitches p_switches;
    bool               p_use_heating = false;

    // Populated by the world-level EOS solve; not serialized. Read and written under p_owner_call_mutex.
    c_LayerEOSData p_eos_data;

    // The owning world and its call lock (set_owner); not serialized, and never copied or moved with the layer.
    c_OwnerCallMutex p_owner_call_mutex;

    // Shared and immutable: a material, a rheology, or a solve's copy of either outlives any change to this layer.
    std::shared_ptr<const c_Material>     p_material;
    std::shared_ptr<const c_RheologyBase> p_shear_rheology;
    std::shared_ptr<const c_RheologyBase> p_bulk_rheology;

    // Shared and immutable, like the rheologies.
    std::shared_ptr<const c_CoolingBase> p_cooling;
    // Owned, observing this layer.
    std::unique_ptr<c_RadiogenicsBase>   p_radiogenics;
};

} // namespace tidalpy
