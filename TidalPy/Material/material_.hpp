#pragma once
/* Phases and materials: the composites that turn the property laws into the state of a material at a point.
 *
 * A c_Phase is one phase of a material (a solid, or its melt): an equation-of-state law, an optional shear-modulus
 * law, optional shear and bulk viscosity laws, optional default rheologies for the layers that use it, and its
 * thermal conductivity and heat capacity. A c_Material is a solid phase, a liquid phase, or both. With both it melts:
 * the solidus and liquidus lie between them, with a melt-weakening law and optional bulk-mixing laws. A material with
 * one phase is that phase everywhere (a water ocean is a liquid-only material).
 *
 * c_Material::calc_state is the one place a point (pressure, temperature, radius) is mapped onto a material's
 * properties: the structure and thermal integrations, the thermal network, the radial solver, and every getter read
 * it. A layer's switches (c_MaterialSwitches) say how much of the material it uses, so one material can serve
 * layers that treat it more or less simply. Phases and materials, like every law, are not changed once built:
 * several layers, and a solve, may hold the same one.
 *
 * All MKS.
 */

#include <cmath>
#include <cstdint>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "binary_.hpp"
#include "broadcast_.hpp"
#include "physics_base_.hpp"
#include "spec_model_.hpp"
#include "thermo_point_.hpp"
#include "laws/eos_law_.hpp"
#include "laws/shear_modulus_law_.hpp"
#include "../constants_.hpp"
#include "../PartialMelt/melt_mixing_.hpp"
#include "../PartialMelt/melt_weakening_.hpp"
#include "../PartialMelt/melting_curve_.hpp"
#include "../Rheology/rheology_.hpp"
#include "../Viscosity/viscosity_.hpp"

namespace tidalpy {

// How much of a material a layer uses; each switch off is the simpler, faster case.
struct c_MaterialSwitches {
    bool use_thermal_expansion = false;   // the density (and thermal pressure) follow the temperature
    bool use_melting           = false;   // the liquid phase, melt fraction, and weakening take part
    bool use_pressure_melting  = false;   // the solidus and liquidus follow the pressure (else read at zero pressure)
    bool use_melt_density      = false;   // the density mixes the phases' by melt fraction (else the solid's)
};

enum class c_MaterialPhase : uint8_t {
    Solid   = 0,
    Partial = 1,
    Liquid  = 2,
};

// One phase's properties at a point.
struct c_PhaseState {
    double density                = TidalPyConstants::d_NAN;   // [kg m-3]
    double bulk_modulus           = TidalPyConstants::d_NAN;   // isothermal [Pa]
    double adiabatic_bulk_modulus = TidalPyConstants::d_NAN;   // [Pa]
    double thermal_expansion      = TidalPyConstants::d_NAN;   // [1/K]
    double heat_capacity          = TidalPyConstants::d_NAN;   // [J kg-1 K-1]
    double thermal_conductivity   = TidalPyConstants::d_NAN;   // [W m-1 K-1]
    double shear_modulus          = 0.0;                       // static [Pa]; 0 for a phase with no shear law
    double shear_viscosity        = TidalPyConstants::d_NAN;   // [Pa s]; NaN for a phase with no viscosity law
    double bulk_viscosity         = TidalPyConstants::d_NAN;   // [Pa s]
};

// A material's properties at a point: what every consumer reads.
struct c_MaterialState {
    c_MaterialPhase phase         = c_MaterialPhase::Solid;
    double density                = TidalPyConstants::d_NAN;   // [kg m-3]
    double bulk_modulus           = TidalPyConstants::d_NAN;   // isothermal [Pa]
    double adiabatic_bulk_modulus = TidalPyConstants::d_NAN;   // [Pa], what a tidal (adiabatic) deformation sees
    double thermal_expansion      = TidalPyConstants::d_NAN;   // [1/K]
    double heat_capacity          = TidalPyConstants::d_NAN;   // effective, latent heat included [J kg-1 K-1]
    // The latent heat's share of heat_capacity [J kg-1 K-1]: L / (T_liq - T_sol) inside a melting range, else zero.
    double latent_heat_capacity = 0.0;
    // The latent heat's share of the adiabat's expansivity [1/K]: nonzero only inside a melting range whose curves
    // follow the pressure (see calc_state). An adiabat runs at
    // dT/dr = -(thermal_expansion + latent_expansion) g T / c_p.
    double latent_expansion = 0.0;
    double thermal_conductivity   = TidalPyConstants::d_NAN;   // [W m-1 K-1]
    double shear_modulus          = TidalPyConstants::d_NAN;   // static, post-melt [Pa]
    double shear_viscosity        = TidalPyConstants::d_NAN;   // post-melt [Pa s]
    double bulk_viscosity         = TidalPyConstants::d_NAN;   // post-melt [Pa s]
    double melt_fraction          = 0.0;                       // [m3 m-3]; NaN without a finite temperature
    double solidus                = TidalPyConstants::d_NAN;   // at the point [K]; NaN for a material that cannot melt
    double liquidus               = TidalPyConstants::d_NAN;   // [K]
};

// =====================================================================================================================
// Phase
// =====================================================================================================================

// The family base of c_Phase. A composite reports its sub-models as tables and has no `model` entry of its own.
class c_PhaseFamily : public c_PhysicsBase {
public:
    static constexpr const char* C_FAMILY_NAME = "phase";
    explicit c_PhaseFamily(const std::string& model_name) : c_PhysicsBase(model_name) {}
    void append_config_entries(std::vector<c_ConfigEntry>& /*out*/) const override {}
};

// A phase's sub-models, set by slot name from shared base pointers (the Python wrappers' form).
struct c_PhaseComponents {
    std::shared_ptr<const c_EOSBase>          eos;
    std::shared_ptr<const c_ShearModulusBase> shear_modulus;
    std::shared_ptr<const c_ViscosityBase>    shear_viscosity;
    std::shared_ptr<const c_ViscosityBase>    bulk_viscosity;
    std::shared_ptr<const c_RheologyBase>     shear_rheology;
    std::shared_ptr<const c_RheologyBase>     bulk_rheology;

    static const std::vector<std::string>& slots() {
        static const std::vector<std::string> names = {
            "eos", "shear_modulus", "shear_viscosity", "bulk_viscosity", "shear_rheology", "bulk_rheology"};
        return names;
    }

    // Throws std::invalid_argument for an unknown slot or a model of the wrong family.
    void set(const std::string& slot, const std::shared_ptr<c_PhysicsBase>& model) {
        const std::string what = "the phase slot '" + slot + "'";
        if (slot == "eos")                  { this->eos = c_share_as<c_EOSBase>(model, what); }
        else if (slot == "shear_modulus")   { this->shear_modulus = c_share_as<c_ShearModulusBase>(model, what); }
        else if (slot == "shear_viscosity") { this->shear_viscosity = c_share_as<c_ViscosityBase>(model, what); }
        else if (slot == "bulk_viscosity")  { this->bulk_viscosity = c_share_as<c_ViscosityBase>(model, what); }
        else if (slot == "shear_rheology")  { this->shear_rheology = c_share_as<c_RheologyBase>(model, what); }
        else if (slot == "bulk_rheology")   { this->bulk_rheology = c_share_as<c_RheologyBase>(model, what); }
        else {
            throw std::invalid_argument(
                "TidalPy: a phase has no slot '" + slot + "'" + c_did_you_mean(slot, slots()) + ".");
        }
    }

    std::shared_ptr<c_PhysicsBase> get(const std::string& slot) const {
        if (slot == "eos")             { return c_share_physics_of(this->eos); }
        if (slot == "shear_modulus")   { return c_share_physics_of(this->shear_modulus); }
        if (slot == "shear_viscosity") { return c_share_physics_of(this->shear_viscosity); }
        if (slot == "bulk_viscosity")  { return c_share_physics_of(this->bulk_viscosity); }
        if (slot == "shear_rheology")  { return c_share_physics_of(this->shear_rheology); }
        if (slot == "bulk_rheology")   { return c_share_physics_of(this->bulk_rheology); }
        throw std::invalid_argument(
            "TidalPy: a phase has no slot '" + slot + "'" + c_did_you_mean(slot, slots()) + ".");
    }
};

class c_Phase final : public c_SpecModel<c_Phase, c_PhaseFamily> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::Phase;

    static const std::vector<c_ParamSpec<c_Phase>>& parameter_specs() {
        using Self = c_Phase;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"thermal_conductivity", "thermal_conductivity_w_mk", &Self::p_thermal_conductivity, 4.0,
             c_ParamBounds::Positive, "Thermal conductivity k0 at the thermal reference temperature [W m-1 K-1]."},
            {"conductivity_temperature_exponent", "conductivity_temperature_exponent",
             &Self::p_conductivity_exponent, 0.0, c_ParamBounds::Finite,
             "n in k = k0 (T / T_ref)^n; -1 for ice [dimensionless]."},
            {"heat_capacity", "heat_capacity_j_kgk", &Self::p_heat_capacity, 1200.0, c_ParamBounds::Positive,
             "Isobaric heat capacity c_p0 at the thermal reference temperature [J kg-1 K-1]."},
            {"heat_capacity_temperature_exponent", "heat_capacity_temperature_exponent",
             &Self::p_heat_capacity_exponent, 0.0, c_ParamBounds::Finite,
             "n in c_p = c_p0 (T / T_ref)^n [dimensionless]."},
            {"thermal_reference_temperature", "thermal_reference_temperature_k", &Self::p_thermal_reference_temperature,
             300.0, c_ParamBounds::Positive, "Temperature where k0 and c_p0 apply [K]."},
        };
        return specs;
    }

    // A default phase (a default instance, which a binary read fills) has a constant-density equation of state.
    c_Phase() : c_Phase(c_ParamMap{}) {}
    explicit c_Phase(const c_ParamMap& params) : c_Phase(params, p_default_components()) {}
    c_Phase(const c_ParamMap& params, c_PhaseComponents components)
        : c_SpecModel("phase"), p_components(std::move(components)) {
        this->p_initialize(params);
    }

    const c_PhaseComponents& get_components() const noexcept { return this->p_components; }
    const c_EOSBase& get_eos() const noexcept { return *this->p_components.eos; }
    bool get_has_shear_viscosity() const noexcept { return static_cast<bool>(this->p_components.shear_viscosity); }

    // The phase at a point; `thermal` says whether its density sees the temperature. Without a law for a property the
    // phase reports none: a shear modulus of 0 (a fluid), NaN viscosities. A shear law's value is floored at
    // [numerical] minimum_modulus.
    void calc_phase_state(const c_ThermoPoint& point, bool thermal, c_PhaseState& out) const noexcept {
        c_EOSPoint law_point;
        this->p_components.eos->calc_eos(point, thermal, law_point);
        out.density                = law_point.density;
        out.bulk_modulus           = law_point.bulk_modulus;
        out.adiabatic_bulk_modulus = law_point.adiabatic_bulk_modulus;
        out.thermal_expansion      = law_point.thermal_expansion;
        this->p_calc_thermal_properties(point.temperature, out);
        out.shear_modulus = 0.0;
        if (this->p_components.shear_modulus) {
            out.shear_modulus = this->p_components.shear_modulus->calc_shear_modulus(point);
            if ((tidalpy_config_ptr != nullptr) && (out.shear_modulus < tidalpy_config_ptr->d_MIN_MODULUS)) {
                out.shear_modulus = tidalpy_config_ptr->d_MIN_MODULUS;
            }
        }
        out.shear_viscosity = this->p_components.shear_viscosity
            ? this->p_components.shear_viscosity->calc_viscosity(point) : TidalPyConstants::d_NAN;
        out.bulk_viscosity = this->p_components.bulk_viscosity
            ? this->p_components.bulk_viscosity->calc_viscosity(point) : TidalPyConstants::d_NAN;
    }

    // Only what the thermal integration needs: density, expansivity, heat capacity, conductivity.
    void calc_phase_thermal(const c_ThermoPoint& point, bool thermal, c_PhaseState& out) const noexcept {
        c_EOSPoint law_point;
        this->p_components.eos->calc_eos(point, thermal, law_point);
        out.density           = law_point.density;
        out.thermal_expansion = law_point.thermal_expansion;
        this->p_calc_thermal_properties(point.temperature, out);
    }

    double calc_density(const c_ThermoPoint& point, bool thermal) const noexcept {
        return this->p_components.eos->calc_density(point, thermal);
    }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_SpecModel::append_config_entries(out);
        for (const std::string& slot : c_PhaseComponents::slots()) {
            const std::shared_ptr<c_PhysicsBase> model = this->p_components.get(slot);
            if (model) { out.push_back(c_config_table(slot, model->get_config_entries())); }
        }
    }

protected:
    static c_PhaseComponents p_default_components() {
        c_PhaseComponents components;
        components.eos = std::make_shared<const c_ConstantEOS>();
        return components;
    }

    void p_calc_thermal_properties(double temperature, c_PhaseState& out) const noexcept {
        out.thermal_conductivity = this->p_thermal_conductivity;
        out.heat_capacity        = this->p_heat_capacity;
        if (std::isfinite(temperature) && (temperature > 0.0)) {
            const double ratio = temperature / this->p_thermal_reference_temperature;
            if (this->p_conductivity_exponent != 0.0) {
                out.thermal_conductivity *= std::pow(ratio, this->p_conductivity_exponent);
            }
            if (this->p_heat_capacity_exponent != 0.0) {
                out.heat_capacity *= std::pow(ratio, this->p_heat_capacity_exponent);
            }
        }
    }

    void p_validate() const override {
        if (!this->p_components.eos) {
            throw std::invalid_argument(this->p_describe() + " needs an equation of state ('eos').");
        }
    }

    // The parameters, then one optional record per slot, in slot order.
    void p_write_payload(std::ostream& out) const override {
        c_SpecModel::p_write_payload(out);
        write_optional_binary(out, this->p_components.eos);
        write_optional_binary(out, this->p_components.shear_modulus);
        write_optional_binary(out, this->p_components.shear_viscosity);
        write_optional_binary(out, this->p_components.bulk_viscosity);
        write_optional_binary(out, this->p_components.shear_rheology);
        write_optional_binary(out, this->p_components.bulk_rheology);
    }

    void p_read_payload(std::istream& in, bool force) override {
        c_SpecModel::p_read_payload(in, force);
        c_PhaseComponents components;
        components.eos = read_optional_binary<c_EOSBase>(in, force, c_eos_from_binary);
        components.shear_modulus =
            read_optional_binary<c_ShearModulusBase>(in, force, c_shear_modulus_from_binary);
        components.shear_viscosity = read_optional_binary<c_ViscosityBase>(in, force, c_viscosity_from_binary);
        components.bulk_viscosity  = read_optional_binary<c_ViscosityBase>(in, force, c_viscosity_from_binary);
        components.shear_rheology  = read_optional_binary<c_RheologyBase>(in, force, c_rheology_from_binary);
        components.bulk_rheology   = read_optional_binary<c_RheologyBase>(in, force, c_rheology_from_binary);
        this->p_components = std::move(components);
        try {
            this->p_validate();
        }
        catch (const std::invalid_argument& component_error) {
            throw std::runtime_error(std::string("TidalPy: corrupt binary data: ") + component_error.what());
        }
    }

    c_PhaseComponents p_components;
    double p_thermal_conductivity          = 0.0;
    double p_conductivity_exponent         = 0.0;
    double p_heat_capacity                 = 0.0;
    double p_heat_capacity_exponent        = 0.0;
    double p_thermal_reference_temperature = 0.0;
};

inline std::unique_ptr<c_Phase> c_phase_from_binary(std::istream& in, bool force = false) {
    const c_BinaryHeader header = c_peek_binary_header(in);
    if (header.class_id != static_cast<uint32_t>(BinaryClassID::Phase)) {
        throw std::runtime_error("TidalPy: expected a phase record in binary stream");
    }
    std::unique_ptr<c_Phase> phase = std::make_unique<c_Phase>();
    phase->read_binary(in, force);
    return phase;
}

// =====================================================================================================================
// Material
// =====================================================================================================================

class c_MaterialFamily : public c_PhysicsBase {
public:
    static constexpr const char* C_FAMILY_NAME = "material";
    explicit c_MaterialFamily(const std::string& model_name) : c_PhysicsBase(model_name) {}
    void append_config_entries(std::vector<c_ConfigEntry>& /*out*/) const override {}
};

// A material's sub-models, set by slot name from shared base pointers (the Python wrappers' form).
struct c_MaterialComponents {
    std::shared_ptr<const c_Phase>                   solid;
    std::shared_ptr<const c_Phase>                   liquid;
    std::shared_ptr<const c_MeltingCurveBase>        solidus;
    std::shared_ptr<const c_MeltingCurveBase>        liquidus;
    std::shared_ptr<const c_MeltWeakeningBase>       weakening;
    std::shared_ptr<const c_BulkModulusMixingBase>   bulk_modulus_mixing;
    std::shared_ptr<const c_BulkViscosityMixingBase> bulk_viscosity_mixing;

    static const std::vector<std::string>& slots() {
        static const std::vector<std::string> names = {
            "solid", "liquid", "solidus", "liquidus", "weakening", "bulk_modulus_mixing", "bulk_viscosity_mixing"};
        return names;
    }

    // The slots a config table nests under its `melting` table.
    static bool is_melting_slot(const std::string& slot) {
        return (slot != "solid") && (slot != "liquid");
    }

    void set(const std::string& slot, const std::shared_ptr<c_PhysicsBase>& model) {
        const std::string what = "the material slot '" + slot + "'";
        if (slot == "solid")          { this->solid = c_share_as<c_Phase>(model, what); }
        else if (slot == "liquid")    { this->liquid = c_share_as<c_Phase>(model, what); }
        else if (slot == "solidus")   { this->solidus = c_share_as<c_MeltingCurveBase>(model, what); }
        else if (slot == "liquidus")  { this->liquidus = c_share_as<c_MeltingCurveBase>(model, what); }
        else if (slot == "weakening") { this->weakening = c_share_as<c_MeltWeakeningBase>(model, what); }
        else if (slot == "bulk_modulus_mixing") {
            this->bulk_modulus_mixing = c_share_as<c_BulkModulusMixingBase>(model, what);
        }
        else if (slot == "bulk_viscosity_mixing") {
            this->bulk_viscosity_mixing = c_share_as<c_BulkViscosityMixingBase>(model, what);
        }
        else {
            throw std::invalid_argument(
                "TidalPy: a material has no slot '" + slot + "'" + c_did_you_mean(slot, slots()) + ".");
        }
    }

    std::shared_ptr<c_PhysicsBase> get(const std::string& slot) const {
        if (slot == "solid")                 { return c_share_physics_of(this->solid); }
        if (slot == "liquid")                { return c_share_physics_of(this->liquid); }
        if (slot == "solidus")               { return c_share_physics_of(this->solidus); }
        if (slot == "liquidus")              { return c_share_physics_of(this->liquidus); }
        if (slot == "weakening")             { return c_share_physics_of(this->weakening); }
        if (slot == "bulk_modulus_mixing")   { return c_share_physics_of(this->bulk_modulus_mixing); }
        if (slot == "bulk_viscosity_mixing") { return c_share_physics_of(this->bulk_viscosity_mixing); }
        throw std::invalid_argument(
            "TidalPy: a material has no slot '" + slot + "'" + c_did_you_mean(slot, slots()) + ".");
    }
};

class c_Material final : public c_SpecModel<c_Material, c_MaterialFamily> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::Material;

    static const std::vector<c_ParamSpec<c_Material>>& parameter_specs() {
        using Self = c_Material;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"latent_heat", "latent_heat_j_kg", &Self::p_latent_heat, 0.0, c_ParamBounds::NonNegative,
             "Latent heat of melting [J kg-1], spread over the melting range as an effective heat capacity."},
        };
        return specs;
    }

    // A default material (a default instance, which a binary read fills) is one default solid phase.
    c_Material() : c_Material(c_ParamMap{}) {}
    explicit c_Material(const c_ParamMap& params) : c_Material(params, p_default_components()) {}
    c_Material(const c_ParamMap& params, c_MaterialComponents components)
        : c_SpecModel("material"), p_components(std::move(components)) {
        this->p_initialize(params);
    }

    const c_MaterialComponents& get_components() const noexcept { return this->p_components; }
    // Null for a material without that phase.
    const c_Phase* get_solid() const noexcept { return this->p_components.solid.get(); }
    const c_Phase* get_liquid() const noexcept { return this->p_components.liquid.get(); }

    // Whether the material can melt: it has a solid and a liquid phase (and so a solidus and liquidus).
    bool get_can_melt() const noexcept {
        return static_cast<bool>(this->p_components.solid) && static_cast<bool>(this->p_components.liquid);
    }
    // Whether the material is liquid everywhere: it has a liquid phase and no solid one.
    bool get_is_liquid_only() const noexcept { return !this->p_components.solid; }
    // The latent heat of melting [J kg-1].
    double get_latent_heat() const noexcept { return this->p_latent_heat; }

    // The phase the material is when it is not melting: the solid, or the liquid of a liquid-only material. Its
    // equation of state and default rheologies are the material's own.
    const c_Phase& get_base_phase() const noexcept {
        return this->p_components.solid ? *this->p_components.solid : *this->p_components.liquid;
    }

    // The solidus and liquidus [K] at a pressure [Pa] as a layer sees them (at zero pressure unless it uses pressure
    // melting); NaN for a material that cannot melt.
    void calc_melting_range(double pressure, const c_MaterialSwitches& switches, double& solidus,
                            double& liquidus) const noexcept {
        if (!this->get_can_melt()) {
            solidus  = TidalPyConstants::d_NAN;
            liquidus = TidalPyConstants::d_NAN;
            return;
        }
        const double melting_pressure = switches.use_pressure_melting ? pressure : 0.0;
        solidus  = this->p_components.solidus->calc_melting_temperature(melting_pressure);
        liquidus = this->p_components.liquidus->calc_melting_temperature(melting_pressure);
    }

    // The material at a point, as a layer with these switches sees it.
    //
    // A liquid-only material is its liquid phase, fully molten, whatever the switches. Otherwise, below the solidus
    // (or with melting off, or for a material that cannot melt) it is the solid phase. Above it the
    // melt fraction phi runs linearly from the solidus to the liquidus (a step at a single melting temperature when
    // they coincide), and:
    //   - shear modulus and viscosity: the weakening law between the solid's and the liquid's values (none: the solid's
    //     until fully molten);
    //   - bulk modulus and bulk viscosity: the mixing laws when present, else the solid's until fully molten;
    //   - density: mixed by volume with use_melt_density, else the solid's;
    //   - expansivity, heat capacity, conductivity: linear in phi; the heat capacity adds the latent heat spread over
    //     the melting range, L / (T_liq - T_sol);
    //   - latent expansion: where the melting curves follow the pressure, the melt fraction changes with pressure too,
    //     and an isentrope, (c_p + L dphi/dT) dT = (alpha T / rho - L dphi/dP) dP, takes the latent heat of that
    //     change as an expansivity, rho L [(1 - phi) dT_sol/dP + phi dT_liq/dP] / ((T_liq - T_sol) T). Without it
    //     an adiabat inside the range is shallower than the melting curve it nears while the one outside is steeper,
    //     and the two meet on the curve, which no integrator can step along.
    // Without a finite temperature there is no melt state: the solid phase with a NaN melt fraction (a liquid-only
    // material stays liquid with a melt fraction of 1). A single melting temperature (a step) adds no latent heat:
    // there is no range to spread it over, so the boundary between the solid and liquid zones of a layer carries it
    // (the world's thermal network, c_zone_boundary_latent_capacity).
    void calc_state(
            const c_ThermoPoint& point,
            const c_MaterialSwitches& switches,
            c_MaterialState& out) const noexcept {
        this->p_evaluate(point, switches, true, out);
    }

    // Only what the thermal integration needs (density, expansivity, heat capacity, conductivity, melt fraction): no
    // moduli, viscosities, or weakening.
    void calc_thermal(
            const c_ThermoPoint& point,
            const c_MaterialSwitches& switches,
            c_MaterialState& out) const noexcept {
        this->p_evaluate(point, switches, false, out);
    }

    // The density alone [kg m-3], what the structure integration asks for.
    double calc_density(const c_ThermoPoint& point, const c_MaterialSwitches& switches) const noexcept {
        const bool thermal = switches.use_thermal_expansion;
        if (!(switches.use_melting && switches.use_melt_density && this->get_can_melt())) {
            return this->get_base_phase().calc_density(point, thermal);
        }
        c_MaterialState state;
        this->p_evaluate(point, switches, false, state);
        return state.density;
    }

    // Whether calc_density can change with pressure: through an equation of state, or through the melt fraction where
    // the melt density is mixed in and the melting range moves with pressure.
    bool get_density_depends_on_pressure(const c_MaterialSwitches& switches) const noexcept {
        if (!(switches.use_melting && switches.use_melt_density && this->get_can_melt())) {
            return this->get_base_phase().get_eos().get_density_depends_on_pressure();
        }
        return switches.use_pressure_melting
            || this->p_components.solid->get_eos().get_density_depends_on_pressure()
            || this->p_components.liquid->get_eos().get_density_depends_on_pressure();
    }

    // ln(mu / mu_min), the post-melt shear modulus against a threshold [Pa]: positive while the material is solid
    // enough for the radial solver's solid equations, negative where it is to be treated as a liquid; -inf for no
    // shear modulus.
    double calc_rigidity_margin(
            const c_ThermoPoint& point,
            const c_MaterialSwitches& switches,
            double minimum_shear_modulus) const noexcept {
        c_MaterialState state;
        this->p_evaluate(point, switches, true, state);
        if (!(state.shear_modulus > 0.0)) { return -TidalPyConstants::d_INF; }
        return std::log(state.shear_modulus / minimum_shear_modulus);
    }

    // Element-wise over pressure, temperature, and radius (each one value per point or a single value used at every
    // point), one state per point.
    void calc_state_vectorize(
            const std::vector<double>& pressure,
            const std::vector<double>& temperature,
            const std::vector<double>& radius,
            const c_MaterialSwitches& switches,
            std::vector<c_MaterialState>& out_states) const {
        const std::size_t num_points = c_broadcast_length(
            {pressure.size(), temperature.size(), radius.size()}, "calc_state_vectorize");
        const std::size_t pressure_stride    = c_broadcast_stride(pressure.size());
        const std::size_t temperature_stride = c_broadcast_stride(temperature.size());
        const std::size_t radius_stride      = c_broadcast_stride(radius.size());
        out_states.resize(num_points);
        c_ThermoPoint point;
        for (std::size_t i = 0; i < num_points; ++i) {
            point.pressure    = pressure[i * pressure_stride];
            point.temperature = temperature[i * temperature_stride];
            point.radius      = radius[i * radius_stride];
            this->calc_state(point, switches, out_states[i]);
        }
    }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_SpecModel::append_config_entries(out);
        std::vector<c_ConfigEntry> melting;
        for (const std::string& slot : c_MaterialComponents::slots()) {
            const std::shared_ptr<c_PhysicsBase> model = this->p_components.get(slot);
            if (!model) { continue; }
            if (c_MaterialComponents::is_melting_slot(slot)) {
                melting.push_back(c_config_table(slot, model->get_config_entries()));
            } else {
                out.push_back(c_config_table(slot, model->get_config_entries()));
            }
        }
        if (!melting.empty()) { out.push_back(c_config_table("melting", melting)); }
    }

protected:
    static c_MaterialComponents p_default_components() {
        c_MaterialComponents components;
        components.solid = std::make_shared<const c_Phase>();
        return components;
    }

    void p_evaluate(
            const c_ThermoPoint& point,
            const c_MaterialSwitches& switches,
            bool mechanical,
            c_MaterialState& out) const noexcept {
        const bool thermal = switches.use_thermal_expansion;
        const bool liquid_only = this->get_is_liquid_only();
        c_PhaseState solid;
        if (mechanical) { this->get_base_phase().calc_phase_state(point, thermal, solid); }
        else            { this->get_base_phase().calc_phase_thermal(point, thermal, solid); }
        p_copy_phase(solid, out);
        out.phase         = liquid_only ? c_MaterialPhase::Liquid : c_MaterialPhase::Solid;
        out.melt_fraction = liquid_only ? 1.0 : 0.0;
        out.latent_expansion = 0.0;
        out.latent_heat_capacity = 0.0;
        out.solidus       = TidalPyConstants::d_NAN;
        out.liquidus      = TidalPyConstants::d_NAN;
        if (!(switches.use_melting && this->get_can_melt())) { return; }

        this->calc_melting_range(point.pressure, switches, out.solidus, out.liquidus);
        const double temperature = point.temperature;
        if (!std::isfinite(temperature)) {
            out.melt_fraction = TidalPyConstants::d_NAN;
            return;
        }
        const double span = out.liquidus - out.solidus;
        const bool step = !(span > TidalPyConstants::d_EPS);
        double phi = step ? ((temperature > out.solidus) ? 1.0 : 0.0) : (temperature - out.solidus) / span;
        if (!(phi > 0.0)) { return; }   // solid (a NaN melting curve leaves the solid too)
        if (phi > 1.0) { phi = 1.0; }
        out.melt_fraction = phi;
        out.phase = (phi >= 1.0) ? c_MaterialPhase::Liquid : c_MaterialPhase::Partial;

        c_PhaseState liquid;
        if (mechanical) { this->p_components.liquid->calc_phase_state(point, thermal, liquid); }
        else            { this->p_components.liquid->calc_phase_thermal(point, thermal, liquid); }

        // Thermal properties: linear in the melt fraction, plus the latent heat over the melting range.
        out.thermal_expansion    = (1.0 - phi) * solid.thermal_expansion + phi * liquid.thermal_expansion;
        out.thermal_conductivity = (1.0 - phi) * solid.thermal_conductivity + phi * liquid.thermal_conductivity;
        out.heat_capacity        = (1.0 - phi) * solid.heat_capacity + phi * liquid.heat_capacity;
        if (!step && (phi < 1.0)) {
            out.latent_heat_capacity = this->p_latent_heat / span;
            out.heat_capacity += out.latent_heat_capacity;
        }
        if (switches.use_melt_density) { out.density = (1.0 - phi) * solid.density + phi * liquid.density; }
        if (!step && (phi < 1.0) && switches.use_pressure_melting) {
            const double melting_slope = (1.0 - phi) * this->p_components.solidus->calc_melting_slope(point.pressure)
                + phi * this->p_components.liquidus->calc_melting_slope(point.pressure);
            const double latent_expansion = out.density * this->p_latent_heat * melting_slope / (span * temperature);
            if (std::isfinite(latent_expansion)) { out.latent_expansion = latent_expansion; }
        }
        if (!mechanical) { return; }

        // Shear modulus and viscosity through the weakening law.
        c_MeltWeakeningInputs inputs;
        inputs.temperature      = temperature;
        inputs.solidus          = out.solidus;
        inputs.liquidus         = out.liquidus;
        inputs.melt_fraction    = phi;
        inputs.solid_shear      = solid.shear_modulus;
        inputs.solid_viscosity  = solid.shear_viscosity;
        inputs.liquid_shear     = liquid.shear_modulus;
        inputs.liquid_viscosity = liquid.shear_viscosity;
        if (this->p_components.weakening && this->p_components.weakening->get_uses_solidus_state()) {
            // A law anchored on the solid's own pair at the solidus continues it into the melting range.
            c_ThermoPoint solidus_point = point;
            solidus_point.temperature = out.solidus;
            c_PhaseState solid_at_solidus;
            this->get_base_phase().calc_phase_state(solidus_point, thermal, solid_at_solidus);
            inputs.solid_shear_at_solidus     = solid_at_solidus.shear_modulus;
            inputs.solid_viscosity_at_solidus = solid_at_solidus.shear_viscosity;
        }
        const c_MeltWeakeningResult weakened = this->p_components.weakening
            ? this->p_components.weakening->calc_weakening(inputs)
            : c_MeltWeakeningResult{(phi >= 1.0) ? liquid.shear_modulus : solid.shear_modulus,
                                    (phi >= 1.0) ? liquid.shear_viscosity : solid.shear_viscosity};
        out.shear_modulus   = weakened.shear_modulus;
        out.shear_viscosity = weakened.viscosity;

        // Bulk modulus through the mixing law. Without one it blends linearly into the liquid's across the weakening
        // law's breakdown band, as the shear modulus does, so the aggregate the radial solver treats as a liquid has
        // the liquid's moduli; with no weakening law both step at full melt.
        if (this->p_components.bulk_modulus_mixing) {
            out.bulk_modulus = this->p_components.bulk_modulus_mixing->calc_bulk_modulus(
                solid.bulk_modulus, liquid.bulk_modulus, out.shear_modulus, phi);
            out.adiabatic_bulk_modulus = this->p_components.bulk_modulus_mixing->calc_bulk_modulus(
                solid.adiabatic_bulk_modulus, liquid.adiabatic_bulk_modulus, out.shear_modulus, phi);
        } else {
            const double blend = this->p_components.weakening
                ? this->p_components.weakening->calc_band_blend(phi) : ((phi >= 1.0) ? 1.0 : 0.0);
            if (blend >= 1.0) {
                out.bulk_modulus           = liquid.bulk_modulus;
                out.adiabatic_bulk_modulus = liquid.adiabatic_bulk_modulus;
            } else if (blend > 0.0) {
                out.bulk_modulus = (1.0 - blend) * solid.bulk_modulus + blend * liquid.bulk_modulus;
                out.adiabatic_bulk_modulus =
                    (1.0 - blend) * solid.adiabatic_bulk_modulus + blend * liquid.adiabatic_bulk_modulus;
            }
        }
        // Bulk viscosity through its mixing law, else the solid's until fully molten.
        if (this->p_components.bulk_viscosity_mixing) {
            out.bulk_viscosity = this->p_components.bulk_viscosity_mixing->calc_bulk_viscosity(
                solid.bulk_viscosity, out.shear_viscosity, phi);
        } else if (phi >= 1.0) {
            out.bulk_viscosity = liquid.bulk_viscosity;
        }
    }

    static void p_copy_phase(const c_PhaseState& phase, c_MaterialState& out) noexcept {
        out.density                = phase.density;
        out.bulk_modulus           = phase.bulk_modulus;
        out.adiabatic_bulk_modulus = phase.adiabatic_bulk_modulus;
        out.thermal_expansion      = phase.thermal_expansion;
        out.heat_capacity          = phase.heat_capacity;
        out.thermal_conductivity   = phase.thermal_conductivity;
        out.shear_modulus          = phase.shear_modulus;
        out.shear_viscosity        = phase.shear_viscosity;
        out.bulk_viscosity         = phase.bulk_viscosity;
    }

    // A material needs at least one phase. With both, the liquid brings the melting curves and needs a viscosity for
    // the melt; melting laws need both phases to melt between.
    void p_validate() const override {
        const c_MaterialComponents& components = this->p_components;
        if (!components.solid && !components.liquid) {
            throw std::invalid_argument(
                this->p_describe() + " needs a solid phase ('solid'), a liquid phase ('liquid'), or both.");
        }
        if (components.solid && components.liquid) {
            if (!components.solidus || !components.liquidus) {
                throw std::invalid_argument(
                    this->p_describe() + " has a liquid phase, so it needs a 'solidus' and a 'liquidus' curve.");
            }
            if (!components.liquid->get_has_shear_viscosity()) {
                throw std::invalid_argument(
                    this->p_describe() + "'s liquid phase needs a 'shear_viscosity' law (the melt's viscosity).");
            }
            // A layer takes its default rheology from the base phase, the solid here, so a rheology on the liquid
            // phase would never be read.
            const c_PhaseComponents& liquid = components.liquid->get_components();
            if (liquid.shear_rheology || liquid.bulk_rheology) {
                throw std::invalid_argument(
                    this->p_describe() + "'s liquid phase has a 'shear_rheology' or 'bulk_rheology', which a material "
                    "with a solid phase never uses: a layer takes its rheology from the solid phase (or its own "
                    "override). Move the rheology to the solid phase, or to the layer.");
            }
            this->p_check_melting_curves_do_not_cross();
        } else if (components.solidus || components.liquidus || components.weakening
                   || components.bulk_modulus_mixing || components.bulk_viscosity_mixing) {
            throw std::invalid_argument(
                this->p_describe() + " has melting laws but not both a solid ('solid') and a liquid ('liquid') phase "
                "for them to melt between.");
        }
    }

    // The liquidus at or above the solidus at zero pressure, where every layer reads them without pressure melting.
    // Deeper, a liquidus at or below the solidus melts as a step at the solidus (c_Material::calc_state).
    void p_check_melting_curves_do_not_cross() const {
        const double solidus_temperature  = this->p_components.solidus->calc_melting_temperature(0.0);
        const double liquidus_temperature = this->p_components.liquidus->calc_melting_temperature(0.0);
        if (liquidus_temperature < solidus_temperature) {
            throw std::invalid_argument(
                this->p_describe() + " has its liquidus below its solidus at zero pressure ("
                + c_format_param_value(liquidus_temperature) + " K against " + c_format_param_value(solidus_temperature)
                + " K).");
        }
    }

    void p_write_payload(std::ostream& out) const override {
        c_SpecModel::p_write_payload(out);
        write_optional_binary(out, this->p_components.solid);
        write_optional_binary(out, this->p_components.liquid);
        write_optional_binary(out, this->p_components.solidus);
        write_optional_binary(out, this->p_components.liquidus);
        write_optional_binary(out, this->p_components.weakening);
        write_optional_binary(out, this->p_components.bulk_modulus_mixing);
        write_optional_binary(out, this->p_components.bulk_viscosity_mixing);
    }

    void p_read_payload(std::istream& in, bool force) override {
        c_SpecModel::p_read_payload(in, force);
        c_MaterialComponents components;
        components.solid    = read_optional_binary<c_Phase>(in, force, c_phase_from_binary);
        components.liquid   = read_optional_binary<c_Phase>(in, force, c_phase_from_binary);
        components.solidus  = read_optional_binary<c_MeltingCurveBase>(in, force, c_melting_curve_from_binary);
        components.liquidus = read_optional_binary<c_MeltingCurveBase>(in, force, c_melting_curve_from_binary);
        components.weakening =
            read_optional_binary<c_MeltWeakeningBase>(in, force, c_melt_weakening_from_binary);
        components.bulk_modulus_mixing =
            read_optional_binary<c_BulkModulusMixingBase>(in, force, c_bulk_modulus_mixing_from_binary);
        components.bulk_viscosity_mixing =
            read_optional_binary<c_BulkViscosityMixingBase>(in, force, c_bulk_viscosity_mixing_from_binary);
        this->p_components = std::move(components);
        try {
            this->p_validate();
        }
        catch (const std::invalid_argument& component_error) {
            throw std::runtime_error(std::string("TidalPy: corrupt binary data: ") + component_error.what());
        }
    }

    c_MaterialComponents p_components;
    double p_latent_heat = 0.0;
};

inline std::unique_ptr<c_Material> c_material_from_binary(std::istream& in, bool force = false) {
    const c_BinaryHeader header = c_peek_binary_header(in);
    if (header.class_id != static_cast<uint32_t>(BinaryClassID::Material)) {
        throw std::runtime_error("TidalPy: expected a material record in binary stream");
    }
    std::unique_ptr<c_Material> material = std::make_unique<c_Material>();
    material->read_binary(in, force);
    return material;
}

// Builders from slot components, for the Python wrappers.
inline std::unique_ptr<c_Phase> c_make_phase(const c_ParamMap& params, const c_PhaseComponents& components) {
    c_PhaseComponents filled = components;
    if (!filled.eos) { filled.eos = std::make_shared<const c_ConstantEOS>(); }
    return std::make_unique<c_Phase>(params, std::move(filled));
}

// A material given neither phase gets a default solid one.
inline std::unique_ptr<c_Material> c_make_material(const c_ParamMap& params, const c_MaterialComponents& components) {
    c_MaterialComponents filled = components;
    if (!filled.solid && !filled.liquid) { filled.solid = std::make_shared<const c_Phase>(); }
    return std::make_unique<c_Material>(params, std::move(filled));
}

}  // namespace tidalpy
