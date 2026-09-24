#pragma once
/* Abstract base for TidalPy cooling models. Concrete models live in cooling_.hpp. All quantities MKS. */

#include <cstdint>
#include <string>
#include <vector>

#include "broadcast_.hpp"
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

    // Element-wise over the temperature drop and the viscosity at an otherwise fixed state (base_inputs). Each holds
    // one value per point or a single value used at every point; out_results is resized.
    void calc_cooling_vectorize(
            const std::vector<double>& delta_temp,
            const std::vector<double>& viscosity,
            const c_CoolingInputs& base_inputs,
            std::vector<c_CoolingResult>& out_results) const
    {
        const std::size_t num_points = c_broadcast_length(
            {delta_temp.size(), viscosity.size()}, "calc_cooling_vectorize");
        const std::size_t delta_temp_stride = c_broadcast_stride(delta_temp.size());
        const std::size_t viscosity_stride  = c_broadcast_stride(viscosity.size());
        out_results.resize(num_points);
        c_CoolingInputs inputs = base_inputs;

        for (std::size_t i = 0; i < num_points; ++i)
        {
            inputs.delta_temp = delta_temp[i * delta_temp_stride];
            inputs.viscosity  = viscosity[i * viscosity_stride];
            out_results[i]    = this->calc_cooling(inputs);
        }
    }
};

} // namespace tidalpy
