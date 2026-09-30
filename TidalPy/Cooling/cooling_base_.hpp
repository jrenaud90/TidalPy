#pragma once
/* Abstract base for TidalPy cooling models. Concrete models live in cooling_.hpp. All quantities MKS. */

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "broadcast_.hpp"
#include "physics_base_.hpp"

namespace tidalpy {

// Physical state passed to a cooling model.
struct c_CoolingInputs {
    double delta_temp           = 0.0;  // temperature drop across the layer, all its boundary layers [K]
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
    double blt             = 0.0;  // conducting thickness carrying the flux across the whole drop [m]
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
    // one value per point or a single value used at every point. The four outputs are caller-owned buffers of
    // num_points values, the broadcast length of the two inputs. Each model implements it with p_vectorize_kernel, so
    // its kernel inlines into the loop rather than costing a virtual call per point.
    virtual void calc_cooling_vectorize(
            const std::vector<double>& delta_temp,
            const std::vector<double>& viscosity,
            const c_CoolingInputs& base_inputs,
            std::size_t num_points,
            double* out_cooling_flux,
            double* out_blt,
            double* out_rayleigh,
            double* out_nusselt) const = 0;

protected:
    // The shared loop over a model's kernel, a callable taking c_CoolingInputs and returning a c_CoolingResult. A model
    // passes a lambda holding a copy of its parameters: the compiler inlines it and, unlike the model's own members,
    // the copy cannot alias the output buffers, so the parameters stay in registers across the loop.
    template <class Kernel>
    static void p_vectorize_kernel(
            const Kernel& kernel,
            const std::vector<double>& delta_temp,
            const std::vector<double>& viscosity,
            const c_CoolingInputs& base_inputs,
            std::size_t num_points,
            double* out_cooling_flux,
            double* out_blt,
            double* out_rayleigh,
            double* out_nusselt)
    {
        if (c_broadcast_length({delta_temp.size(), viscosity.size()}, "calc_cooling_vectorize") != num_points)
        {
            throw std::invalid_argument("TidalPy::calc_cooling_vectorize: the output buffers do not match the inputs");
        }
        if (num_points == 0)
        {
            return;
        }
        // Which inputs vary is fixed at compile time: a runtime broadcast stride in the index costs a sweep about an
        // eighth of its time.
        const bool delta_temp_varies = c_broadcast_stride(delta_temp.size()) == 1;
        const bool viscosity_varies  = c_broadcast_stride(viscosity.size()) == 1;
        if (delta_temp_varies && viscosity_varies)
        {
            p_sweep<true, true>(kernel, delta_temp, viscosity, base_inputs, num_points,
                                out_cooling_flux, out_blt, out_rayleigh, out_nusselt);
        }
        else if (delta_temp_varies)
        {
            p_sweep<true, false>(kernel, delta_temp, viscosity, base_inputs, num_points,
                                 out_cooling_flux, out_blt, out_rayleigh, out_nusselt);
        }
        else if (viscosity_varies)
        {
            p_sweep<false, true>(kernel, delta_temp, viscosity, base_inputs, num_points,
                                 out_cooling_flux, out_blt, out_rayleigh, out_nusselt);
        }
        else
        {
            p_sweep<false, false>(kernel, delta_temp, viscosity, base_inputs, num_points,
                                  out_cooling_flux, out_blt, out_rayleigh, out_nusselt);
        }
    }

private:
    template <bool DELTA_TEMP_VARIES, bool VISCOSITY_VARIES, class Kernel>
    static void p_sweep(
            const Kernel& kernel,
            const std::vector<double>& delta_temp,
            const std::vector<double>& viscosity,
            const c_CoolingInputs& base_inputs,
            std::size_t num_points,
            double* out_cooling_flux,
            double* out_blt,
            double* out_rayleigh,
            double* out_nusselt)
    {
        const double* const delta_temp_data = delta_temp.data();
        const double* const viscosity_data  = viscosity.data();
        c_CoolingInputs inputs = base_inputs;
        inputs.delta_temp = delta_temp_data[0];
        inputs.viscosity  = viscosity_data[0];

        for (std::size_t i = 0; i < num_points; ++i)
        {
            if constexpr (DELTA_TEMP_VARIES) { inputs.delta_temp = delta_temp_data[i]; }
            if constexpr (VISCOSITY_VARIES)  { inputs.viscosity  = viscosity_data[i]; }
            const c_CoolingResult result = kernel(inputs);
            out_cooling_flux[i] = result.cooling_flux;
            out_blt[i]          = result.blt;
            out_rayleigh[i]     = result.rayleigh_number;
            out_nusselt[i]      = result.nusselt_number;
        }
    }
};

} // namespace tidalpy
