#pragma once
// Abstract base for all TidalPy radiogenics models.
// Full concrete models (IsotopeRadiogenics, FixedRadiogenics) arrive in Phase 7.
#include <string>
#include "physics_base_.hpp"

namespace tidalpy {

class c_RadiogenicsBase : public c_PhysicsBase {
public:
    c_RadiogenicsBase() = default;
    explicit c_RadiogenicsBase(const std::string& model_name) : c_PhysicsBase(model_name) {}
    ~c_RadiogenicsBase() override = default;

    // Returns total radiogenic heating [W] for the given mass and elapsed time.
    virtual double calc_heating(double time_s, double mass_kg) const = 0;
};

} // namespace tidalpy
