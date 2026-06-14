#pragma once
/*
 * cooling_base_.hpp - c_CoolingBase: abstract base for all TidalPy cooling models.
 *
 * Inherits c_PhysicsBase (Utilities_x/classes_x/physics_base_.hpp).
 *
 * Defines the abstract interface that every cooling model must satisfy:
 *   calc_cooling(c_CoolingInputs) -> c_CoolingResult
 *
 * The three implemented models (Off, Convective, Conductive) live in cooling_.hpp.
 *
 * All calc_* methods are const and operate in MKS units. The physical state at a
 * cooling evaluation is bundled in c_CoolingInputs (eight quantities, so a struct
 * is used per the style guide), and the result is bundled in c_CoolingResult.
 */

#include <stdexcept>
#include <string>
#include <vector>

#include "physics_base_.hpp"

namespace tidalpy {

// -------------------------------------------------------------------------------
// c_CoolingInputs — the physical state passed to a cooling model (all MKS).
// -------------------------------------------------------------------------------
struct c_CoolingInputs {
    double delta_temp_k              = 0.0;   // temperature drop across the layer [K]
    double thickness_m               = 0.0;   // layer (or sub-layer) thickness [m]
    double gravity_m_s2              = 0.0;   // gravitational acceleration [m/s^2]
    double density_kg_m3             = 0.0;   // bulk density [kg/m^3]
    double viscosity_pas             = 0.0;   // dynamic viscosity [Pa·s]
    double thermal_conductivity_w_mk = 0.0;   // thermal conductivity [W/m/K]
    double thermal_diffusivity_m2_s  = 0.0;   // thermal diffusivity [m^2/s]
    double thermal_expansion_1_k     = 0.0;   // thermal expansivity [1/K]
};

// -------------------------------------------------------------------------------
// c_CoolingResult — the quantities every cooling model reports.
// -------------------------------------------------------------------------------
struct c_CoolingResult {
    double cooling_flux_w_m2 = 0.0;   // heat flux leaving the layer [W/m^2]
    double blt_m             = 0.0;   // boundary-layer thickness [m]
    double rayleigh_number   = 0.0;   // Rayleigh number [dimensionless]
    double nusselt_number    = 1.0;   // Nusselt number [dimensionless]
};

// -------------------------------------------------------------------------------
// c_CoolingBase
// -------------------------------------------------------------------------------
class c_CoolingBase : public c_PhysicsBase {
public:
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------
    c_CoolingBase() = default;

    explicit c_CoolingBase(const std::string& model_name) : c_PhysicsBase(model_name) {}

    ~c_CoolingBase() override = default;

    // -----------------------------------------------------------------------
    // Cooling (pure virtual)
    //
    // Each model overrides this to map the layer's physical state to a cooling
    // result (heat flux, boundary-layer thickness, Rayleigh and Nusselt numbers).
    // c_SolidLiquidLayer calls this to obtain a layer's surface heat flux.
    //
    // Assumptions
    // -----------
    // - Steady-state boundary-layer theory; all inputs/outputs MKS.
    // -----------------------------------------------------------------------
    virtual c_CoolingResult calc_cooling(const c_CoolingInputs& inputs) const = 0;

    // -----------------------------------------------------------------------
    // Vectorized cooling — vary the temperature drop at otherwise fixed state.
    //
    // Evaluates calc_cooling over each delta_temp, copying the base inputs and
    // overriding their delta_temp_k. out_results is resized to the input length.
    // -----------------------------------------------------------------------
    void calc_cooling_vectorize_temperature(
            const std::vector<double>& delta_temp_k,
            const c_CoolingInputs& base_inputs,
            std::vector<c_CoolingResult>& out_results) const {
        const std::size_t n = delta_temp_k.size();
        out_results.resize(n);
        c_CoolingInputs inputs = base_inputs;
        for (std::size_t i = 0; i < n; ++i) {
            inputs.delta_temp_k = delta_temp_k[i];
            out_results[i] = this->calc_cooling(inputs);
        }
    }

    // -----------------------------------------------------------------------
    // Vectorized cooling — vary the viscosity at otherwise fixed state.
    // -----------------------------------------------------------------------
    void calc_cooling_vectorize_viscosity(
            const std::vector<double>& viscosity_pas,
            const c_CoolingInputs& base_inputs,
            std::vector<c_CoolingResult>& out_results) const {
        const std::size_t n = viscosity_pas.size();
        out_results.resize(n);
        c_CoolingInputs inputs = base_inputs;
        for (std::size_t i = 0; i < n; ++i) {
            inputs.viscosity_pas = viscosity_pas[i];
            out_results[i] = this->calc_cooling(inputs);
        }
    }

    // -----------------------------------------------------------------------
    // Vectorized cooling — vary temperature drop and viscosity element-wise.
    //
    // The two input vectors must have the same length N; out_results is resized
    // to N. Throws std::invalid_argument if the input vectors differ in length.
    // -----------------------------------------------------------------------
    void calc_cooling_vectorize_all(
            const std::vector<double>& delta_temp_k,
            const std::vector<double>& viscosity_pas,
            const c_CoolingInputs& base_inputs,
            std::vector<c_CoolingResult>& out_results) const {
        if (delta_temp_k.size() != viscosity_pas.size()) {
            throw std::invalid_argument(
                "TidalPy: calc_cooling_vectorize_all — delta_temp and viscosity "
                "vectors must have the same length");
        }
        const std::size_t n = delta_temp_k.size();
        out_results.resize(n);
        c_CoolingInputs inputs = base_inputs;
        for (std::size_t i = 0; i < n; ++i) {
            inputs.delta_temp_k  = delta_temp_k[i];
            inputs.viscosity_pas = viscosity_pas[i];
            out_results[i] = this->calc_cooling(inputs);
        }
    }
};

} // namespace tidalpy
