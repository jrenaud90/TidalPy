#pragma once
/* Abstract base for TidalPy viscosity models. Concrete models live in viscosity_.hpp.
 *
 * The result is the pre-melt (solid) viscosity that the partial-melt step weakens. Like the partial-melt outputs it
 * is frequency independent.
 *
 * References
 * ----------
 * - Moore (2006): Arrhenius (activation energy and volume) flow law.
 * - Henning (2009): reference-viscosity (relative activation) law.
 */

#include <string>
#include <vector>

#include "broadcast_.hpp"
#include "physics_base_.hpp"

namespace tidalpy {

class c_ViscosityBase : public c_PhysicsBase {
public:
    static constexpr const char* C_FAMILY_NAME = "viscosity";

    c_ViscosityBase() = default;

    explicit c_ViscosityBase(const std::string& model_name) : c_PhysicsBase(model_name) {}

    ~c_ViscosityBase() override = default;

    // Dynamic viscosity [Pa s] at a temperature [K] and pressure [Pa]. Assumes a steady-state flow law.
    virtual double calc_viscosity(double temperature, double pressure) const noexcept = 0;

    // Element-wise over temperature and pressure. Each input holds one value per point or a single value used at
    // every point.
    void calc_viscosity_vectorize(
            const std::vector<double>& temperature,
            const std::vector<double>& pressure,
            std::vector<double>& out_viscosity) const {
        const std::size_t num_points =
            c_broadcast_length({temperature.size(), pressure.size()}, "calc_viscosity_vectorize");
        const std::size_t temperature_stride = c_broadcast_stride(temperature.size());
        const std::size_t pressure_stride    = c_broadcast_stride(pressure.size());
        out_viscosity.resize(num_points);
        for (std::size_t i = 0; i < num_points; ++i) {
            out_viscosity[i] = this->calc_viscosity(
                temperature[i * temperature_stride], pressure[i * pressure_stride]);
        }
    }
};

} // namespace tidalpy
