#pragma once
/*
 * viscosity_base_.hpp - c_ViscosityBase: abstract base for TidalPy viscosity models.
 *
 * Inherits c_PhysicsBase. A viscosity model returns a material's dynamic viscosity [Pa s] at a
 * temperature [K] and pressure [Pa]. This is the pre-melt (solid) viscosity that the partial-melt
 * step weakens; like the partial-melt outputs it is frequency independent and is cached once per
 * EOS solve. The concrete models (Arrhenius, Reference, Constant) live in viscosity_.hpp.
 *
 * References
 * ----------
 * - Moore (2006): Arrhenius (activation energy and volume) flow law.
 * - Henning (2009): reference-viscosity (relative activation) law.
 */

#include <stdexcept>
#include <string>
#include <vector>

#include "physics_base_.hpp"

namespace tidalpy {

class c_ViscosityBase : public c_PhysicsBase {
public:
    c_ViscosityBase() = default;

    explicit c_ViscosityBase(const std::string& model_name) : c_PhysicsBase(model_name) {}

    ~c_ViscosityBase() override = default;

    // Dynamic viscosity [Pa s] at a temperature [K] and pressure [Pa]. Assumes a steady-state flow law.
    virtual double calc_viscosity(double temperature, double pressure) const = 0;

    // Vectorized element-wise over temperature and pressure: the primary radial sweep, one entry per
    // slice. The two input vectors must match in length.
    void calc_viscosity_vectorize(
            const std::vector<double>& temperature,
            const std::vector<double>& pressure,
            std::vector<double>& out_viscosity) const {
        const std::size_t n = temperature.size();
        if (pressure.size() != n) {
            throw std::invalid_argument(
                "TidalPy: calc_viscosity_vectorize — temperature and pressure "
                "vectors must have the same length");
        }
        out_viscosity.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            out_viscosity[i] = this->calc_viscosity(temperature[i], pressure[i]);
        }
    }
};

} // namespace tidalpy
