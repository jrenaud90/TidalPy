#pragma once
/*
 * cooling_.hpp - Implements TidalPy's cooling (heat-transport) models.
 *
 * Inherits c_CoolingBase (cooling_base_.hpp), which itself inherits c_PhysicsBase.
 * Each model implements calc_cooling(c_CoolingInputs), returning a c_CoolingResult
 * (heat flux [W/m^2], boundary-layer thickness [m], Rayleigh and Nusselt numbers).
 *
 * Models (with config aliases handled by the factory):
 *   c_OffCooling         (alias "none")  — cooling disabled (zero flux)
 *   c_ConvectiveCooling                  — parameterized boundary-layer convection
 *   c_ConductiveCooling                  — conduction across the layer
 *
 * All quantities are MKS. The math mirrors the validated legacy implementation in
 * TidalPy/cooling/cooling_models.py.
 *
 * References
 * ----------
 * - Turcotte and Schubert (2002), Geodynamics — Rayleigh/Nusselt convection scaling.
 * - Solomatov (1995); Schubert, Turcotte, and Olson (2001) — boundary-layer theory.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::<Model> (401-403)
 *   payload: model_name length (uint32_t) | model_name bytes | model params (doubles)
 *   Off / Conduction write zero params; Convection writes its three scalars.
 *   The layer observer pointer (p_layer_ptr) is NOT serialized.
 */

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <istream>
#include <limits>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "constants_.hpp"
#include "cooling_base_.hpp"

namespace tidalpy {

// -------------------------------------------------------------------------------
// Numerical floor used to guard denominators that may approach zero.
// -------------------------------------------------------------------------------
inline constexpr double d_COOLING_FLOOR = 1.0e-100;

// Replace a magnitude smaller than the floor with a signed floor value.
inline double cool_guard(double value) noexcept {
    if (std::abs(value) < d_COOLING_FLOOR) {
        return (value < 0.0) ? -d_COOLING_FLOOR : d_COOLING_FLOOR;
    }
    return value;
}

// -------------------------------------------------------------------------------
// c_CoolingConfig — construction parameters for the convection model.
// (The off and conduction models take no parameters.)
// -------------------------------------------------------------------------------
struct c_CoolingConfig {
    double convection_alpha  = 1.0;                  // Nu = alpha * (Ra / Ra_crit)^beta  [dimensionless]
    double convection_beta   = 0.3333333333333333;   // convection exponent (~1/3)        [dimensionless]
    double critical_rayleigh = 1100.0;               // critical Rayleigh number           [dimensionless]
};

// =====================================================================================================================
// Cooling functions
//
// Each mirrors the validated legacy implementation in
// TidalPy/cooling/cooling_models.py.
// =====================================================================================================================

// Off: no cooling. Boundary layer is half the layer thickness; flux zero.
inline c_CoolingResult cool_off(const c_CoolingInputs& in) noexcept {
    c_CoolingResult result;
    result.cooling_flux_w_m2 = 0.0;
    result.blt_m             = 0.5 * in.thickness_m;
    result.rayleigh_number   = 0.0;
    result.nusselt_number    = 1.0;
    return result;
}

// Conduction: flux = k * delta_temp / thickness; boundary layer = thickness.
inline c_CoolingResult cool_conduction(const c_CoolingInputs& in) noexcept {
    c_CoolingResult result;
    result.blt_m             = in.thickness_m;
    result.cooling_flux_w_m2 = in.thermal_conductivity_w_mk * in.delta_temp_k
                             / cool_guard(in.thickness_m);
    result.rayleigh_number   = 0.0;
    result.nusselt_number    = 1.0;
    return result;
}

// Parameterized convection via the Rayleigh number.
//
//   Ra = expansion * density * gravity * delta_temp * thickness^3 / (viscosity * diffusivity)
//   Nu = max(alpha * (Ra / Ra_crit)^beta, 2)
//   boundary layer = thickness / Nu
//   flux = k * delta_temp / boundary_layer
//
// Degenerate inputs (delta_temp <= 0, thickness below the minimum-layer-thickness
// floor) collapse to the legacy edge behaviour (Ra = 0, Nu = 2). The minimum
// thickness comes from the shared TidalPy config (tidalpy_config_ptr->d_MIN_THICKNESS).
inline c_CoolingResult cool_convection(
        const c_CoolingInputs& in, const c_CoolingConfig& cfg) noexcept {
    const double eps = std::numeric_limits<double>::epsilon();
    const double min_thickness = tidalpy_config_ptr->d_MIN_THICKNESS;
    c_CoolingResult result;

    const double rate_heat_loss   = in.thermal_diffusivity_m2_s / cool_guard(in.thickness_m);
    const double parcel_rise_rate = in.thermal_expansion_1_k * in.density_kg_m3 * in.gravity_m_s2
                                  * in.delta_temp_k * in.thickness_m * in.thickness_m
                                  / cool_guard(in.viscosity_pas);

    double rayleigh = parcel_rise_rate / cool_guard(rate_heat_loss);
    if (!(in.delta_temp_k > eps))          { rayleigh = 0.0; }
    if (!(in.thickness_m  >= min_thickness)) { rayleigh = 0.0; }

    double nusselt = cfg.convection_alpha
                   * std::pow(rayleigh / cool_guard(cfg.critical_rayleigh), cfg.convection_beta);
    if (in.delta_temp_k <= eps)            { nusselt = 2.0; }
    if (in.thickness_m  <= min_thickness)  { nusselt = 2.0; }
    if (nusselt <= 2.0)                    { nusselt = 2.0; }

    double blt = in.thickness_m / cool_guard(nusselt);
    if (in.delta_temp_k <= eps)            { blt = 1.0; }
    if (in.thickness_m  <= min_thickness)  { blt = in.thickness_m; }

    result.cooling_flux_w_m2 = in.thermal_conductivity_w_mk * in.delta_temp_k / cool_guard(blt);
    result.blt_m             = blt;
    result.rayleigh_number   = rayleigh;
    result.nusselt_number    = nusselt;
    return result;
}

// -------------------------------------------------------------------------------
// Lower-case a model name for case-insensitive factory lookup.
// -------------------------------------------------------------------------------
inline std::string cool_to_lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

// =====================================================================================================================
// Cooling models
//
// Binary serialization uses the shared c_PhysicsBase helpers
// (write_physics_binary / read_physics_binary): Off and Conduction write zero
// params; Convection writes its three scalars.
// =====================================================================================================================

// -------------------------------------------------------------------------------
// c_OffCooling — cooling disabled (alias "none").
// -------------------------------------------------------------------------------
class c_OffCooling : public c_CoolingBase {
public:
    c_OffCooling() : c_CoolingBase("off") {}
    explicit c_OffCooling(const c_CoolingConfig& /*cfg*/) : c_CoolingBase("off") {}
    ~c_OffCooling() override = default;

    c_CoolingResult calc_cooling(const c_CoolingInputs& inputs) const override {
        return cool_off(inputs);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::OffCooling));
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }
};

// -------------------------------------------------------------------------------
// c_ConductiveCooling — conduction across the layer.
// -------------------------------------------------------------------------------
class c_ConductiveCooling : public c_CoolingBase {
public:
    c_ConductiveCooling() : c_CoolingBase("conduction") {}
    explicit c_ConductiveCooling(const c_CoolingConfig& /*cfg*/) : c_CoolingBase("conduction") {}
    ~c_ConductiveCooling() override = default;

    c_CoolingResult calc_cooling(const c_CoolingInputs& inputs) const override {
        return cool_conduction(inputs);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::ConductiveCooling));
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }
};

// -------------------------------------------------------------------------------
// c_ConvectiveCooling — parameterized boundary-layer convection.
// -------------------------------------------------------------------------------
class c_ConvectiveCooling : public c_CoolingBase {
public:
    c_ConvectiveCooling() : c_CoolingBase("convection") {}
    explicit c_ConvectiveCooling(const c_CoolingConfig& cfg)
        : c_CoolingBase("convection"),
          p_convection_alpha(cfg.convection_alpha),
          p_convection_beta(cfg.convection_beta),
          p_critical_rayleigh(cfg.critical_rayleigh) {}
    ~c_ConvectiveCooling() override = default;

    double get_convection_alpha()  const noexcept { return this->p_convection_alpha; }
    double get_convection_beta()   const noexcept { return this->p_convection_beta; }
    double get_critical_rayleigh() const noexcept { return this->p_critical_rayleigh; }

    c_CoolingResult calc_cooling(const c_CoolingInputs& inputs) const override {
        c_CoolingConfig cfg;
        cfg.convection_alpha  = this->p_convection_alpha;
        cfg.convection_beta   = this->p_convection_beta;
        cfg.critical_rayleigh = this->p_critical_rayleigh;
        return cool_convection(inputs, cfg);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::ConvectiveCooling),
            {this->p_convection_alpha, this->p_convection_beta, this->p_critical_rayleigh});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 3);
        this->p_convection_alpha  = params[0];
        this->p_convection_beta   = params[1];
        this->p_critical_rayleigh = params[2];
    }

protected:
    double p_convection_alpha  = 1.0;
    double p_convection_beta   = 0.3333333333333333;
    double p_critical_rayleigh = 1100.0;
};

// =====================================================================================================================
// Factory
// =====================================================================================================================

// -------------------------------------------------------------------------------
// c_CoolingModel — one named value per cooling model.
// -------------------------------------------------------------------------------
enum class c_CoolingModel : uint8_t {
    Off        = 0,
    Convection = 1,
    Conduction = 2,
};

// -------------------------------------------------------------------------------
// c_cooling_model_from_name — map a (case-insensitive) model name or alias to a
// c_CoolingModel enum value.
//
// Recognised names and aliases:
//   "off" / "none"
//   "convection" / "convective"
//   "conduction" / "conductive"
//
// Throws std::invalid_argument if the model name is unknown.
// -------------------------------------------------------------------------------
inline c_CoolingModel c_cooling_model_from_name(const std::string& model_name) {
    const std::string name = cool_to_lower(model_name);

    if (name == "off"        || name == "none")       { return c_CoolingModel::Off; }
    if (name == "convection" || name == "convective") { return c_CoolingModel::Convection; }
    if (name == "conduction" || name == "conductive") { return c_CoolingModel::Conduction; }

    throw std::invalid_argument("TidalPy: unknown cooling model name '" + model_name + "'");
}

// -------------------------------------------------------------------------------
// c_find_cooling — build the cooling model named by a c_CoolingModel enum.
//
// Returns a unique_ptr to a newly heap-allocated concrete model constructed from
// the supplied config. Throws std::invalid_argument for an unrecognised enum value.
// -------------------------------------------------------------------------------
inline std::unique_ptr<c_CoolingBase> c_find_cooling(
        c_CoolingModel model, const c_CoolingConfig& cfg) {
    switch (model) {
        case c_CoolingModel::Off:        return std::make_unique<c_OffCooling>(cfg);
        case c_CoolingModel::Convection: return std::make_unique<c_ConvectiveCooling>(cfg);
        case c_CoolingModel::Conduction: return std::make_unique<c_ConductiveCooling>(cfg);
    }
    throw std::invalid_argument("TidalPy: unrecognised c_CoolingModel enum value");
}

// -------------------------------------------------------------------------------
// c_find_cooling (name overload) — maps a name/alias to the enum and builds the
// model. Throws std::invalid_argument on unknown names.
// -------------------------------------------------------------------------------
inline std::unique_ptr<c_CoolingBase> c_find_cooling(
        const std::string& model_name, const c_CoolingConfig& cfg) {
    return c_find_cooling(c_cooling_model_from_name(model_name), cfg);
}

// -------------------------------------------------------------------------------
// c_cooling_from_binary — reconstruct a cooling model from a binary stream.
//
// Peeks the upcoming record's BinaryClassID (without consuming the header),
// constructs the matching default-initialised concrete model, then delegates to
// its read_binary to restore the model name and parameters. Used by the layer
// recursive deserialization (see structures_x/layers). Throws std::runtime_error
// if the class id is not a known cooling model.
// -------------------------------------------------------------------------------
inline std::unique_ptr<c_CoolingBase> c_cooling_from_binary(std::istream& in, bool force = false) {
    const std::streampos start = in.tellg();
    const c_BinaryHeader header = read_binary_header(in);
    in.seekg(start);

    std::unique_ptr<c_CoolingBase> model;
    switch (static_cast<BinaryClassID>(header.class_id)) {
        case BinaryClassID::OffCooling:        model = std::make_unique<c_OffCooling>();        break;
        case BinaryClassID::ConvectiveCooling: model = std::make_unique<c_ConvectiveCooling>(); break;
        case BinaryClassID::ConductiveCooling: model = std::make_unique<c_ConductiveCooling>(); break;
        default:
            throw std::runtime_error("TidalPy: unknown cooling class id in binary stream");
    }
    model->read_binary(in, force);
    return model;
}

}  // namespace tidalpy
