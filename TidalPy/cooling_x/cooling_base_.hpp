#pragma once
// Abstract base for all TidalPy cooling models.
// Full concrete models (ConvectiveCooling, ConductiveCooling) arrive in Phase 6.
#include <string>
#include "physics_base_.hpp"

namespace tidalpy {

struct c_CoolingResult {
    double cooling_flux_w_m2      = 0.0;   // [W/m²]
    double blt_m                  = 0.0;   // boundary-layer thickness [m]
    double rayleigh_number        = 0.0;   // [dimensionless]
    double nusselt_number         = 1.0;   // [dimensionless]
};

class c_CoolingBase : public c_PhysicsBase {
public:
    c_CoolingBase() = default;
    explicit c_CoolingBase(const std::string& model_name) : c_PhysicsBase(model_name) {}
    ~c_CoolingBase() override = default;

    virtual c_CoolingResult calc_cooling(
            double delta_temp_k,
            double thickness_m,
            double gravity_m_s2,
            double density_kg_m3,
            double viscosity_pas,
            double thermal_conductivity_w_mk,
            double thermal_diffusivity_m2_s,
            double thermal_expansion_1_k) const = 0;
};

} // namespace tidalpy
