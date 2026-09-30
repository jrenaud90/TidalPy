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
 * Binary payload: the c_BaseLayer payload, then the Love numbers k, h, l (real then imaginary part each), the three
 * classification flags, the temperature, use_thermal_eos, and use_heating, then the shear and bulk rheologies, each
 * behind a presence flag.
 */

#include <complex>
#include <cstdint>
#include <initializer_list>
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
    // The owning world applies the rheology through p_complex_modulus while it holds its call lock.
    friend class c_LayeredWorld;

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

    uint32_t get_binary_class_id() const override {
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
    // profile is read under the owning world's call lock (c_BaseLayer::set_owner_call_mutex).
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

    // What load_binary reads a file into first (c_TidalPyBaseClass::make_binary_scratch).
    std::unique_ptr<c_TidalPyBaseClass> make_binary_scratch() const override {
        return std::make_unique<c_PhysicsLayer>();
    }

protected:
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

    void p_write_payload(std::ostream& out) const override {
        c_BaseLayer::p_write_payload(out);
        for (const std::complex<double>& love_number :
                {this->p_love_numbers.k, this->p_love_numbers.h, this->p_love_numbers.l}) {
            const double parts[2] = {love_number.real(), love_number.imag()};
            out.write(reinterpret_cast<const char*>(parts), sizeof(parts));
        }
        const uint8_t flag_bytes[3] = {
            static_cast<uint8_t>(this->p_is_solid),
            static_cast<uint8_t>(this->p_is_static),
            static_cast<uint8_t>(this->p_is_incompressible)};
        out.write(reinterpret_cast<const char*>(flag_bytes), sizeof(flag_bytes));
        out.write(reinterpret_cast<const char*>(&this->p_temperature), sizeof(double));
        const uint8_t state_bytes[2] = {
            static_cast<uint8_t>(this->p_use_thermal_eos), static_cast<uint8_t>(this->p_use_heating)};
        out.write(reinterpret_cast<const char*>(state_bytes), sizeof(state_bytes));
        write_optional_binary(out, this->p_shear_rheology);
        write_optional_binary(out, this->p_bulk_rheology);
    }

    void p_read_payload(std::istream& in, bool force) override {
        c_BaseLayer::p_read_payload(in, force);
        for (std::complex<double>* love_number :
                {&this->p_love_numbers.k, &this->p_love_numbers.h, &this->p_love_numbers.l}) {
            double parts[2] = {0.0, 0.0};
            in.read(reinterpret_cast<char*>(parts), sizeof(parts));
            *love_number = std::complex<double>(parts[0], parts[1]);
        }
        uint8_t flag_bytes[3] = {0, 0, 0};
        in.read(reinterpret_cast<char*>(flag_bytes), sizeof(flag_bytes));
        this->p_is_solid          = static_cast<bool>(flag_bytes[0]);
        this->p_is_static         = static_cast<bool>(flag_bytes[1]);
        this->p_is_incompressible = static_cast<bool>(flag_bytes[2]);
        in.read(reinterpret_cast<char*>(&this->p_temperature), sizeof(double));
        uint8_t state_bytes[2] = {0, 0};
        in.read(reinterpret_cast<char*>(state_bytes), sizeof(state_bytes));
        this->p_use_thermal_eos = static_cast<bool>(state_bytes[0]);
        this->p_use_heating     = static_cast<bool>(state_bytes[1]);
        this->p_shear_rheology = read_optional_binary<c_RheologyBase>(in, force, c_rheology_from_binary);
        if (this->p_shear_rheology) { this->p_shear_rheology->set_layer_ptr(this); }
        this->p_bulk_rheology = read_optional_binary<c_RheologyBase>(in, force, c_rheology_from_binary);
        if (this->p_bulk_rheology) { this->p_bulk_rheology->set_layer_ptr(this); }
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

    // Optional rheology objects.
    // The rheology classes are shared not unique: a radial-solver solution exported to Python keeps a
    // copy of these pointers so it can reproduce the complex moduli it was solved with, at any radius, after
    // potentially after this layer is gone.
    std::shared_ptr<c_RheologyBase> p_shear_rheology;
    std::shared_ptr<c_RheologyBase> p_bulk_rheology;
};

} // namespace tidalpy
