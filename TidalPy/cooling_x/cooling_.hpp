#pragma once
/* TidalPy's cooling (heat-transport) models. All quantities MKS.
 *
 * References
 * ----------
 * - Turcotte and Schubert (2002), Geodynamics: Rayleigh and Nusselt convection scaling.
 * - Solomatov (1995); Schubert, Turcotte, and Olson (2001): boundary-layer theory.
 *
 * Binary payload: model name then the model's parameters as doubles. The observing layer pointer is not
 * serialized.
 */

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
#include "model_names_.hpp"

namespace tidalpy {

// Construction parameters for the convection model, and its defaults; the other models take none.
struct c_CoolingConfig {
    double convection_alpha  = 1.0;                  // Nu = alpha * (Ra / Ra_crit)^beta  [dimensionless]
    double convection_beta   = 0.3333333333333333;   // convection exponent (~1/3)        [dimensionless]
    double critical_rayleigh = 1100.0;               // critical Rayleigh number           [dimensionless]
};

// No cooling. Boundary layer is set to half the layer thickness so downstream users see a sane value.
inline c_CoolingResult cool_off(const c_CoolingInputs& in) noexcept {
    c_CoolingResult result;
    result.cooling_flux    = 0.0;
    result.blt             = 0.5 * in.thickness;
    result.rayleigh_number = 0.0;
    result.nusselt_number  = 1.0;
    return result;
}

// Conduction across the whole layer: flux = k * delta_temp / thickness.
inline c_CoolingResult cool_conduction(const c_CoolingInputs& in) noexcept {
    c_CoolingResult result;
    result.blt             = in.thickness;
    result.cooling_flux    = in.thermal_conductivity * in.delta_temp / c_guard_denominator(in.thickness);
    result.rayleigh_number = 0.0;
    result.nusselt_number  = 1.0;
    return result;
}

// Parameterized convection via the Rayleigh number.
//
//   Ra = expansion * density * gravity * delta_temp * thickness^3 / (viscosity * diffusivity)
//   Nu = max(alpha * (Ra / Ra_crit)^beta, Nu_min)
//   boundary layer = thickness / Nu
//   flux = k * delta_temp / boundary_layer
//
// Nu_min defaults to 2: Nu = 1 is conduction across the whole layer, so the floor keeps a barely convecting
// layer losing heat through a boundary layer half the layer thick. Degenerate inputs (no temperature
// contrast, or a vanishingly thin layer) collapse to Ra = 0 and Nu = Nu_min. Each test is made once so that
// Ra, Nu, and the boundary layer stay consistent with one another.
inline c_CoolingResult cool_convection(
        const c_CoolingInputs& in, const c_CoolingConfig& cfg) noexcept {
    const double eps = TidalPyConstants::d_EPS;
    // An unwired config gives NaN limits, so the result shows the missing initialization.
    const bool config_wired    = (tidalpy_config_ptr != nullptr);
    const double min_thickness = config_wired ? tidalpy_config_ptr->d_MIN_THICKNESS : TidalPyConstants::d_NAN;
    const double min_nusselt   = config_wired ? tidalpy_config_ptr->d_MIN_NUSSELT : TidalPyConstants::d_NAN;
    c_CoolingResult result;

    const double rate_heat_loss   = in.thermal_diffusivity / c_guard_denominator(in.thickness);
    const double parcel_rise_rate = in.thermal_expansion * in.density * in.gravity
                                  * in.delta_temp * in.thickness * in.thickness
                                  / c_guard_denominator(in.viscosity);

    const bool no_contrast = !(in.delta_temp > eps);
    const bool too_thin    = !(in.thickness > min_thickness);

    double rayleigh = parcel_rise_rate / c_guard_denominator(rate_heat_loss);
    if (no_contrast || too_thin) { rayleigh = 0.0; }

    double nusselt = cfg.convection_alpha
                   * std::pow(rayleigh / c_guard_denominator(cfg.critical_rayleigh), cfg.convection_beta);
    // A NaN from another input (the viscosity, say) falls through so it reaches the caller.
    if (no_contrast || too_thin || (nusselt <= min_nusselt)) { nusselt = min_nusselt; }

    // With no contrast the boundary layer still follows the floored Nusselt number: the thermal network uses it
    // as the resistance between this layer and its neighbors, whose temperatures can differ from this one's.
    double blt = in.thickness / c_guard_denominator(nusselt);
    if (too_thin) { blt = in.thickness; }

    result.cooling_flux    = in.thermal_conductivity * in.delta_temp / c_guard_denominator(blt);
    result.blt             = blt;
    result.rayleigh_number = rayleigh;
    result.nusselt_number  = nusselt;
    return result;
}

// Cooling disabled (alias "none").
class c_OffCooling final : public c_CoolingBase {
public:
    c_OffCooling() : c_OffCooling(c_CoolingConfig{}) {}
    explicit c_OffCooling(const c_CoolingConfig& /*cfg*/) : c_CoolingBase("off") {}
    ~c_OffCooling() override = default;

    c_CoolingResult calc_cooling(const c_CoolingInputs& inputs) const override {
        return cool_off(inputs);
    }
    c_CoolingModel get_model_type() const noexcept override { return c_CoolingModel::Off; }

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::OffCooling); }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, this->get_binary_class_id());
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }
};

class c_ConductiveCooling final : public c_CoolingBase {
public:
    c_ConductiveCooling() : c_ConductiveCooling(c_CoolingConfig{}) {}
    explicit c_ConductiveCooling(const c_CoolingConfig& /*cfg*/) : c_CoolingBase("conduction") {}
    ~c_ConductiveCooling() override = default;

    c_CoolingResult calc_cooling(const c_CoolingInputs& inputs) const override {
        return cool_conduction(inputs);
    }
    c_CoolingModel get_model_type() const noexcept override { return c_CoolingModel::Conduction; }

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::ConductiveCooling); }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, this->get_binary_class_id());
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }
};

// Parameterized boundary-layer convection.
class c_ConvectiveCooling final : public c_CoolingBase {
public:
    c_ConvectiveCooling() : c_ConvectiveCooling(c_CoolingConfig{}) {}
    explicit c_ConvectiveCooling(const c_CoolingConfig& cfg)
        : c_CoolingBase("convection"),
          p_config(cfg) {}
    ~c_ConvectiveCooling() override = default;

    double get_convection_alpha()  const noexcept { return this->p_config.convection_alpha; }
    double get_convection_beta()   const noexcept { return this->p_config.convection_beta; }
    double get_critical_rayleigh() const noexcept { return this->p_config.critical_rayleigh; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_CoolingBase::append_config_entries(out);
        out.push_back(c_config_double("convection_alpha", this->p_config.convection_alpha));
        out.push_back(c_config_double("convection_beta", this->p_config.convection_beta));
        out.push_back(c_config_double("critical_rayleigh", this->p_config.critical_rayleigh));
    }

    c_CoolingResult calc_cooling(const c_CoolingInputs& inputs) const override {
        return cool_convection(inputs, this->p_config);
    }
    c_CoolingModel get_model_type() const noexcept override { return c_CoolingModel::Convection; }

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::ConvectiveCooling); }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, this->get_binary_class_id(),
            {this->p_config.convection_alpha, this->p_config.convection_beta, this->p_config.critical_rayleigh});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 3);
        this->p_config.convection_alpha  = params[0];
        this->p_config.convection_beta   = params[1];
        this->p_config.critical_rayleigh = params[2];
    }

protected:
    c_CoolingConfig p_config;
};

// Model names are matched case-insensitively.
inline c_CoolingModel c_cooling_model_from_name(const std::string& model_name) {
    const std::string name = c_to_lower(model_name);

    if (name == "off"        || name == "none")       { return c_CoolingModel::Off; }
    if (name == "convection" || name == "convective") { return c_CoolingModel::Convection; }
    if (name == "conduction" || name == "conductive") { return c_CoolingModel::Conduction; }

    throw std::invalid_argument("TidalPy: unknown cooling model name '" + model_name + "'");
}

// Builds a model from its enum value and parameters; the Cython wrappers construct through it. A saved record is
// restored by c_cooling_from_binary instead.
inline std::unique_ptr<c_CoolingBase> c_find_cooling(
        c_CoolingModel model, const c_CoolingConfig& cfg) {
    switch (model) {
        case c_CoolingModel::Off:        return std::make_unique<c_OffCooling>(cfg);
        case c_CoolingModel::Convection: return std::make_unique<c_ConvectiveCooling>(cfg);
        case c_CoolingModel::Conduction: return std::make_unique<c_ConductiveCooling>(cfg);
    }
    throw std::invalid_argument("TidalPy: unrecognised c_CoolingModel enum value");
}

inline std::unique_ptr<c_CoolingBase> c_find_cooling(
        const std::string& model_name, const c_CoolingConfig& cfg) {
    return c_find_cooling(c_cooling_model_from_name(model_name), cfg);
}

// The class id is peeked without consuming the header so the matching default-constructed model can
// restore the record itself.
inline std::unique_ptr<c_CoolingBase> c_cooling_from_binary(std::istream& in, bool force = false) {
    const c_BinaryHeader header = c_peek_binary_header(in);

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
