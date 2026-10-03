#pragma once
/* Abstract base for TidalPy viscosity models. Concrete models live in viscosity_.hpp.
 *
 * The result is a phase's viscosity, which a material's melt weakening may then lower (c_Material). It is frequency
 * independent.
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
#include "thermo_point_.hpp"

namespace tidalpy {

class c_ViscosityBase : public c_PhysicsBase {
public:
    static constexpr const char* C_FAMILY_NAME = "viscosity";

    c_ViscosityBase() = default;

    explicit c_ViscosityBase(const std::string& model_name) : c_PhysicsBase(model_name) {}

    ~c_ViscosityBase() override = default;

    // Dynamic viscosity [Pa s] at a point (pressure [Pa], temperature [K], radius [m]). Assumes a steady-state flow
    // law; only a tabulated model reads the radius.
    virtual double calc_viscosity(const c_ThermoPoint& point) const noexcept = 0;

    // Element-wise over temperature, pressure, and radius. Each input holds one value per point or a single value
    // used at every point.
    void calc_viscosity_vectorize(
            const std::vector<double>& temperature,
            const std::vector<double>& pressure,
            const std::vector<double>& radius,
            std::vector<double>& out_viscosity) const {
        const std::size_t num_points = c_broadcast_length(
            {temperature.size(), pressure.size(), radius.size()}, "calc_viscosity_vectorize");
        const std::size_t temperature_stride = c_broadcast_stride(temperature.size());
        const std::size_t pressure_stride    = c_broadcast_stride(pressure.size());
        const std::size_t radius_stride      = c_broadcast_stride(radius.size());
        out_viscosity.resize(num_points);
        c_ThermoPoint point;
        for (std::size_t i = 0; i < num_points; ++i) {
            point.pressure    = pressure[i * pressure_stride];
            point.temperature = temperature[i * temperature_stride];
            point.radius      = radius[i * radius_stride];
            out_viscosity[i]  = this->calc_viscosity(point);
        }
    }
};

} // namespace tidalpy
