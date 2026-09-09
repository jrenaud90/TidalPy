#pragma once
/*
 * viscosity_base_.hpp — c_ViscosityBase: abstract base for all TidalPy
 * (solid/liquid) viscosity models.
 *
 * Inherits c_PhysicsBase (Utilities_x/classes_x/physics_base_.hpp).
 *
 * A viscosity model returns a material's dynamic viscosity [Pa·s] as a function of
 * temperature [K] and pressure [Pa]. This is the PRE-melt ("solid") viscosity that
 * the partial-melt step (partial_melt) subsequently weakens. Like the partial-
 * melt outputs it is frequency-independent and cached once per EOS solve in the
 * whole-planet love-number pipeline.
 *
 * The three implemented models (Arrhenius, Reference, Constant) live in
 * viscosity_.hpp. All calc_* methods are const and operate in MKS units.
 *
 * References
 * ----------
 * - Moore (2006) — Arrhenius (activation energy/volume) flow law.
 * - Henning (2009) — reference-viscosity (relative activation) law.
 */

#include <stdexcept>
#include <string>
#include <vector>

#include "physics_base_.hpp"

namespace tidalpy {

// -------------------------------------------------------------------------------
// c_ViscosityBase
// -------------------------------------------------------------------------------
class c_ViscosityBase : public c_PhysicsBase {
public:
    c_ViscosityBase() = default;

    explicit c_ViscosityBase(const std::string& model_name) : c_PhysicsBase(model_name) {}

    ~c_ViscosityBase() override = default;

    // -----------------------------------------------------------------------
    // Dynamic viscosity [Pa·s] at temperature [K] and pressure [Pa] (pure virtual).
    //
    // Assumptions
    // -----------
    // - Steady-state flow law; all inputs/outputs MKS.
    // -----------------------------------------------------------------------
    virtual double calc_viscosity(double temperature_k, double pressure_pa) const = 0;

    // -----------------------------------------------------------------------
    // Vectorized viscosity — vary temperature and pressure element-wise. The two
    // input vectors must have the same length N; out_viscosity is resized to N.
    // Throws std::invalid_argument on a length mismatch. This is the primary
    // radial sweep (one entry per slice).
    // -----------------------------------------------------------------------
    void calc_viscosity_vectorize(
            const std::vector<double>& temperature_k,
            const std::vector<double>& pressure_pa,
            std::vector<double>& out_viscosity) const {
        const std::size_t n = temperature_k.size();
        if (pressure_pa.size() != n) {
            throw std::invalid_argument(
                "TidalPy: calc_viscosity_vectorize — temperature and pressure "
                "vectors must have the same length");
        }
        out_viscosity.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            out_viscosity[i] = this->calc_viscosity(temperature_k[i], pressure_pa[i]);
        }
    }
};

} // namespace tidalpy
