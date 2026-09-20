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
 *     temperature, shear-modulus pressure derivative, temperature derivative, and reference temperature
 *                                   (double×4, 32)
 *     use_thermal_eos               (uint8_t, 1)
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
    // Material state. The layer temperature is 0 K until set: the cold, rigid limit of the viscosity laws.
    double        temperature = 0.0;   // [K]
    // Static shear modulus law: mu = mu0 + mu'_P P + mu'_T (T - T_ref), with mu0 = shear_modulus_static.
    double        shear_modulus_pressure_derivative    = 0.0;                          // mu'_P [Pa/Pa]
    double        shear_modulus_temperature_derivative = 0.0;                          // mu'_T [Pa/K]
    double        shear_modulus_reference_temperature  = d_EOS_REFERENCE_TEMPERATURE;  // T_ref [K]
    bool          use_thermal_eos = false;    // the EOS density and bulk modulus see the temperature
};

// Material properties of a layer at one point: c_PhysicsLayer::calc_material_state.
struct c_MaterialState {
    double density       = TidalPyConstants::d_NAN;   // [kg/m^3]; NaN without an EOS model
    double melt_fraction = 0.0;                       // [m^3/m^3]
    // Static moduli [Pa] and viscosities [Pa·s] before the partial-melt model, then after it.
    double premelt_shear_modulus   = TidalPyConstants::d_NAN;
    double premelt_bulk_modulus    = TidalPyConstants::d_NAN;
    double premelt_shear_viscosity = TidalPyConstants::d_NAN;
    double premelt_bulk_viscosity  = TidalPyConstants::d_NAN;
    double shear_modulus   = TidalPyConstants::d_NAN;
    double bulk_modulus    = TidalPyConstants::d_NAN;
    double shear_viscosity = TidalPyConstants::d_NAN;
    double bulk_viscosity  = TidalPyConstants::d_NAN;
    // Complex moduli [Pa] at the forcing frequency; the post-melt static moduli for a non-finite frequency.
    std::complex<double> complex_shear_modulus = {TidalPyConstants::d_NAN, 0.0};
    std::complex<double> complex_bulk_modulus  = {TidalPyConstants::d_NAN, 0.0};
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
          p_is_incompressible(cfg.is_incompressible),
          p_temperature(cfg.temperature),
          p_shear_modulus_pressure_derivative(cfg.shear_modulus_pressure_derivative),
          p_shear_modulus_temperature_derivative(cfg.shear_modulus_temperature_derivative),
          p_shear_modulus_reference_temperature(cfg.shear_modulus_reference_temperature),
          p_use_thermal_eos(cfg.use_thermal_eos)
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
            this->p_temperature       = other.p_temperature;
            this->p_shear_modulus_pressure_derivative    = other.p_shear_modulus_pressure_derivative;
            this->p_shear_modulus_temperature_derivative = other.p_shear_modulus_temperature_derivative;
            this->p_shear_modulus_reference_temperature  = other.p_shear_modulus_reference_temperature;
            this->p_use_thermal_eos   = other.p_use_thermal_eos;
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
    double get_shear_modulus_static()   const noexcept { return this->p_shear_modulus_static; }
    double get_bulk_modulus_static()    const noexcept { return this->p_bulk_modulus_static; }
    double get_shear_viscosity_static() const noexcept { return this->p_shear_viscosity_static; }
    double get_bulk_viscosity_static()  const noexcept { return this->p_bulk_viscosity_static; }

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

    // Material-state parameters.
    double get_temperature()     const noexcept { return this->p_temperature; }
    bool   get_use_thermal_eos() const noexcept { return this->p_use_thermal_eos; }
    double get_shear_modulus_pressure_derivative() const noexcept {
        return this->p_shear_modulus_pressure_derivative;
    }
    double get_shear_modulus_temperature_derivative() const noexcept {
        return this->p_shear_modulus_temperature_derivative;
    }
    double get_shear_modulus_reference_temperature() const noexcept {
        return this->p_shear_modulus_reference_temperature;
    }
    void set_temperature(double value)   noexcept { this->p_temperature = value; }
    void set_use_thermal_eos(bool value) noexcept { this->p_use_thermal_eos = value; }

    // Material properties at a radius [m], pressure [Pa], and temperature [K]
    void calc_material_state(
            double radius,
            double pressure,
            double temperature,
            double frequency,
            c_MaterialState& out) const noexcept {
        const double min_modulus = tidalpy_config_ptr->d_MIN_MODULUS;

        // Shear law of the layer, floored so a steep temperature derivative cannot drive it negative.
        double shear = this->p_shear_modulus_static + this->p_shear_modulus_pressure_derivative * pressure;
        if (std::isfinite(temperature)) {
            shear += this->p_shear_modulus_temperature_derivative
                * (temperature - this->p_shear_modulus_reference_temperature);
        }
        double bulk            = this->p_bulk_modulus_static;
        double shear_viscosity = this->p_shear_viscosity
            ? this->p_shear_viscosity->calc_viscosity(temperature, pressure) : this->p_shear_viscosity_static;
        double bulk_viscosity  = this->p_bulk_viscosity
            ? this->p_bulk_viscosity->calc_viscosity(temperature, pressure) : this->p_bulk_viscosity_static;

        out.density = TidalPyConstants::d_NAN;
        if (this->p_eos) {
            const double eos_temperature = this->p_use_thermal_eos ? temperature : TidalPyConstants::d_NAN;
            double eos_bulk = TidalPyConstants::d_NAN;
            this->p_eos->calc_density_and_bulk_modulus(pressure, eos_temperature, radius, out.density, eos_bulk);
            const double eos_shear           = this->p_eos->calc_static_shear_modulus(radius);
            const double eos_shear_viscosity = this->p_eos->calc_shear_viscosity(radius);
            const double eos_bulk_viscosity  = this->p_eos->calc_bulk_viscosity(radius);
            if (std::isfinite(eos_shear))           { shear = eos_shear; }
            if (std::isfinite(eos_bulk))            { bulk = eos_bulk; }
            if (std::isfinite(eos_shear_viscosity)) { shear_viscosity = eos_shear_viscosity; }
            if (std::isfinite(eos_bulk_viscosity))  { bulk_viscosity = eos_bulk_viscosity; }
        }
        if (shear < min_modulus) { shear = min_modulus; }

        out.premelt_shear_modulus   = shear;
        out.premelt_bulk_modulus    = bulk;
        out.premelt_shear_viscosity = shear_viscosity;
        out.premelt_bulk_viscosity  = bulk_viscosity;
        out.melt_fraction           = 0.0;
        if (this->p_partial_melt) {
            // The liquid viscosity is the pre-melt viscosity until a dedicated liquid-viscosity model exists.
            c_PartialMeltInputs inputs;
            inputs.temperature       = temperature;
            inputs.premelt_viscosity = shear_viscosity;
            inputs.premelt_shear     = shear;
            inputs.liquid_viscosity  = shear_viscosity;
            const c_PartialMeltResult shear_result = this->p_partial_melt->calc_partial_melt(inputs);
            inputs.premelt_viscosity = bulk_viscosity;
            inputs.premelt_shear     = bulk;
            inputs.liquid_viscosity  = bulk_viscosity;
            const c_PartialMeltResult bulk_result = this->p_partial_melt->calc_partial_melt(inputs);
            out.melt_fraction = shear_result.melt_fraction;
            shear             = shear_result.postmelt_shear_modulus;
            shear_viscosity   = shear_result.postmelt_viscosity;
            bulk              = bulk_result.postmelt_shear_modulus;
            bulk_viscosity    = bulk_result.postmelt_viscosity;
        }
        out.shear_modulus   = shear;
        out.bulk_modulus    = bulk;
        out.shear_viscosity = shear_viscosity;
        out.bulk_viscosity  = bulk_viscosity;

        const bool use_rheology = std::isfinite(frequency);
        out.complex_shear_modulus = (use_rheology && this->p_shear_rheology)
            ? this->p_shear_rheology->calc_complex_modulus(shear, shear_viscosity, frequency)
            : std::complex<double>(shear, 0.0);
        out.complex_bulk_modulus = (use_rheology && this->p_bulk_rheology)
            ? this->p_bulk_rheology->calc_complex_modulus(bulk, bulk_viscosity, frequency)
            : std::complex<double>(bulk, 0.0);
    }

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

    // The layer's material state at a radius on the solved structure: reads the pressure and temperature from the
    // EOS profile and maps them onto the attached models. This is the one route every radius-resolved getter and
    // the radial Love solve take, so none of them can read a value off a slice grid or disagree with the solve.
    // A non-finite frequency skips the rheology and leaves the complex moduli at the static values. Returns false,
    // leaving out at its defaults, before an EOS profile is stored.
    bool calc_material_state_at(double radius, double frequency, c_MaterialState& out) const noexcept {
        if (!this->p_eos_data.is_populated()) { return false; }
        double state[C_EOS_DY_VALUES];
        this->p_eos_data.evaluate(radius, state);
        // The solve stores a temperature at every radius (its segment's uniform value when temperature was not
        // integrated); a profile supplied through update_eos_data carries none, so the layer's own is used.
        double temperature = state[C_EOS_TEMPERATURE_INDEX];
        if (!std::isfinite(temperature)) { temperature = this->p_temperature; }
        this->calc_material_state(radius, state[c_LayerEOSData::EOS_INDEX_PRESSURE], temperature, frequency, out);
        return true;
    }

    // Radius-resolved viscoelastic state, evaluated on demand through calc_material_state_at.
    bool get_viscoelastic_populated() const noexcept override { return this->p_eos_data.is_populated(); }

    double get_shear_modulus(double radius)   const noexcept override {
        return this->p_state_value(radius, &c_MaterialState::shear_modulus);
    }
    double get_bulk_modulus(double radius)    const noexcept override {
        return this->p_state_value(radius, &c_MaterialState::bulk_modulus);
    }
    double get_shear_viscosity(double radius) const noexcept override {
        return this->p_state_value(radius, &c_MaterialState::shear_viscosity);
    }
    double get_bulk_viscosity(double radius)  const noexcept override {
        return this->p_state_value(radius, &c_MaterialState::bulk_viscosity);
    }
    double get_premelt_shear_modulus(double radius)   const noexcept override {
        return this->p_state_value(radius, &c_MaterialState::premelt_shear_modulus);
    }
    double get_premelt_bulk_modulus(double radius)    const noexcept override {
        return this->p_state_value(radius, &c_MaterialState::premelt_bulk_modulus);
    }
    double get_premelt_shear_viscosity(double radius) const noexcept override {
        return this->p_state_value(radius, &c_MaterialState::premelt_shear_viscosity);
    }
    double get_premelt_bulk_viscosity(double radius)  const noexcept override {
        return this->p_state_value(radius, &c_MaterialState::premelt_bulk_viscosity);
    }

    // Melt fraction at a radius on the solved structure; 0 without a partial-melt model, NaN before a solve.
    double get_melt_fraction(double radius) const noexcept override {
        return this->p_state_value(radius, &c_MaterialState::melt_fraction);
    }

    // Radius-resolved complex moduli [Pa] at a forcing frequency, from the material state at that radius.
    // Feeds the radial Love-number solve. Purely real (no dissipation) when no rheology is attached; NaN before
    // an EOS profile is stored.
    std::complex<double> calc_complex_shear_modulus(
            double radius, double frequency) const noexcept {
        c_MaterialState state;
        if (!this->calc_material_state_at(radius, frequency, state)) {
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        }
        return state.complex_shear_modulus;
    }

    std::complex<double> calc_complex_bulk_modulus(
            double radius, double frequency) const noexcept {
        c_MaterialState state;
        if (!this->calc_material_state_at(radius, frequency, state)) {
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        }
        return state.complex_bulk_modulus;
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
            sizeof(uint8_t)  * 2 +           // is_tidal, is_volume_fixed
            sizeof(double)   +               // tidal_scale
            sizeof(uint8_t)  +               // tidal_scale_method
            sizeof(double)   * 4 +           // shear modulus, bulk modulus, shear viscosity, bulk viscosity
            sizeof(double)   * 6 +           // love_number k, h, l (each: re + im)
            sizeof(uint8_t)  * 3 +           // is_solid, is_static, is_incompressible
            material_law_bytes() +           // temperature, shear law, use_thermal_eos
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
        const uint8_t is_volume_fixed_byte = static_cast<uint8_t>(this->p_is_volume_fixed);
        out.write(reinterpret_cast<const char*>(&is_tidal_byte),       sizeof(uint8_t));
        out.write(reinterpret_cast<const char*>(&is_volume_fixed_byte), sizeof(uint8_t));
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
        this->write_material_law_binary(out);

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
        uint8_t is_volume_fixed_byte = 0;
        in.read(reinterpret_cast<char*>(&is_volume_fixed_byte), sizeof(uint8_t));
        this->p_is_volume_fixed = static_cast<bool>(is_volume_fixed_byte);

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
        this->read_material_law_binary(in);

        if (!in) {
            throw std::runtime_error("TidalPy: failed to read PhysicsLayer binary data");
        }

        this->read_eos_model_binary(in, force);
        this->read_physics_models_binary(in, force);

        this->update_physicals();
    }

protected:
    // One field of the material state at a radius, for the radius-resolved getters. The rheology is skipped
    // (NaN frequency), so this costs the EOS inversion and the static models only.
    double p_state_value(double radius, double c_MaterialState::* field) const noexcept {
        c_MaterialState state;
        if (!this->calc_material_state_at(radius, TidalPyConstants::d_NAN, state)) {
            return TidalPyConstants::d_NAN;
        }
        return state.*field;
    }

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

    // The material-state scalars (temperature, shear law, use_thermal_eos), shared with the subclasses so the
    // three layer records keep one byte layout for them.
    static constexpr uint64_t material_law_bytes() { return 4 * sizeof(double) + sizeof(uint8_t); }

    void write_material_law_binary(std::ostream& out) const {
        out.write(reinterpret_cast<const char*>(&this->p_temperature), sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_shear_modulus_pressure_derivative),    sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_shear_modulus_temperature_derivative), sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_shear_modulus_reference_temperature),  sizeof(double));
        const uint8_t use_thermal_eos_byte = static_cast<uint8_t>(this->p_use_thermal_eos);
        out.write(reinterpret_cast<const char*>(&use_thermal_eos_byte), sizeof(uint8_t));
    }

    void read_material_law_binary(std::istream& in) {
        in.read(reinterpret_cast<char*>(&this->p_temperature), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_shear_modulus_pressure_derivative),    sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_shear_modulus_temperature_derivative), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_shear_modulus_reference_temperature),  sizeof(double));
        uint8_t use_thermal_eos_byte = 0;
        in.read(reinterpret_cast<char*>(&use_thermal_eos_byte), sizeof(uint8_t));
        this->p_use_thermal_eos = static_cast<bool>(use_thermal_eos_byte);
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
    // Material state (see c_PhysicsConfig).
    double        p_temperature = 0.0;   // [K]
    double        p_shear_modulus_pressure_derivative    = 0.0;                          // [Pa/Pa]
    double        p_shear_modulus_temperature_derivative = 0.0;                          // [Pa/K]
    double        p_shear_modulus_reference_temperature  = d_EOS_REFERENCE_TEMPERATURE;  // [K]
    bool          p_use_thermal_eos = false;

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
