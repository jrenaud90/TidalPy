#pragma once
/* Abstract base for TidalPy cooling models. Concrete models live in cooling_.hpp. All quantities MKS. */

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "physics_base_.hpp"

namespace tidalpy {

// Physical state passed to a cooling model.
struct c_CoolingInputs {
    double delta_temp           = 0.0;  // temperature drop across the layer [K]
    double thickness            = 0.0;  // layer (or sub-layer) thickness [m]
    double gravity              = 0.0;  // gravitational acceleration [m/s^2]
    double density              = 0.0;  // bulk density [kg/m^3]
    double viscosity            = 0.0;  // dynamic viscosity [Pa·s]
    double thermal_conductivity = 0.0;  // thermal conductivity [W/m/K]
    double thermal_diffusivity  = 0.0;  // thermal diffusivity [m^2/s]
    double thermal_expansion    = 0.0;  // thermal expansivity [1/K]
};

struct c_CoolingResult {
    double cooling_flux    = 0.0;  // heat flux leaving the layer [W/m^2]
    double blt             = 0.0;  // boundary-layer thickness [m]
    double rayleigh_number = 0.0;  // Rayleigh number [dimensionless]
    double nusselt_number  = 1.0;  // Nusselt number [dimensionless]
};

enum class c_CoolingModel : uint8_t {
    Off        = 0,
    Convection = 1,
    Conduction = 2,
};

class c_CoolingBase : public c_PhysicsBase {
public:
    c_CoolingBase() = default;

    explicit c_CoolingBase(const std::string& model_name) : c_PhysicsBase(model_name) {}

    ~c_CoolingBase() override = default;

    // Assumes steady-state boundary-layer theory.
    virtual c_CoolingResult calc_cooling(const c_CoolingInputs& inputs) const = 0;

    // Which heat transport the model describes. The thermal network builds a layer's profile from this, so a
    // model is never told apart by its name.
    virtual c_CoolingModel get_model_type() const noexcept = 0;

    // Vectorized over the temperature drop at otherwise fixed state; out_results is resized.
    void calc_cooling_vectorize_temperature(
            const std::vector<double>& delta_temp,
            const c_CoolingInputs& base_inputs,
            std::vector<c_CoolingResult>& out_results) const
    {
        const std::size_t n = delta_temp.size();
        out_results.resize(n);
        c_CoolingInputs inputs = base_inputs;
        
        for (std::size_t i = 0; i < n; ++i)
        {
            inputs.delta_temp = delta_temp[i];
            out_results[i]    = this->calc_cooling(inputs);
        }
    }

    // Vectorized over the viscosity at otherwise fixed state.
    void calc_cooling_vectorize_viscosity(
            const std::vector<double>& viscosity,
            const c_CoolingInputs& base_inputs,
            std::vector<c_CoolingResult>& out_results) const
    {
        const std::size_t n = viscosity.size();
        out_results.resize(n);
        c_CoolingInputs inputs = base_inputs;
        
        for (std::size_t i = 0; i < n; ++i)
        {
            inputs.viscosity = viscosity[i];
            out_results[i]   = this->calc_cooling(inputs);
        }
    }

    // Vectorized element-wise over temperature drop and viscosity.
    void calc_cooling_vectorize_all(
            const std::vector<double>& delta_temp,
            const std::vector<double>& viscosity,
            const c_CoolingInputs& base_inputs,
            std::vector<c_CoolingResult>& out_results) const
    {
        if (delta_temp.size() != viscosity.size())
        {
            throw std::invalid_argument(
                "TidalPy::calc_cooling_vectorize_all: delta_temp and viscosity "
                "vectors must have the same length");
        }
        const std::size_t n = delta_temp.size();
        out_results.resize(n);
        c_CoolingInputs inputs = base_inputs;
        
        for (std::size_t i = 0; i < n; ++i)
        {
            inputs.delta_temp = delta_temp[i];
            inputs.viscosity  = viscosity[i];
            out_results[i]    = this->calc_cooling(inputs);
        }
    }
};

} // namespace tidalpy
