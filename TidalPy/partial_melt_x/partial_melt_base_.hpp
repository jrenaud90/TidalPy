#pragma once
/*
 * partial_melt_base_.hpp: c_PartialMeltBase, the abstract base for TidalPy partial-melt models, derived
 * from c_PhysicsBase.
 *
 * A partial-melt model maps a material's pre-melt (solid) viscosity and shear modulus, its temperature,
 * and its melt fraction to the post-melt viscosity and shear modulus. The melt fraction itself is
 * model-independent and lives on this base class. The three models (Off, Spohn, Henning) are in
 * partial_melt_.hpp. All calc_* methods are const and MKS.
 *
 * These quantities are frequency independent, so the world pipeline caches them once after the EOS solve
 * and only the downstream rheology (complex modulus) step repeats per forcing frequency.
 *
 * References
 * ----------
 * - Fischer and Spohn (1990): temperature-based melt viscosity and shear law.
 * - Henning (2009, 2010); Renaud and Henning (2018): three-regime melt weakening.
 */

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "physics_base_.hpp"
#include "../constants_.hpp"  // TidalPyConstants::d_EPS

namespace tidalpy {

// -------------------------------------------------------------------------------
// c_PartialMeltInputs: the per-evaluation state, all MKS. Material constants live
// on the model object; only the varying state is passed here.
// -------------------------------------------------------------------------------
struct c_PartialMeltInputs {
    double temperature       = 0.0;   // local temperature [K]
    double premelt_viscosity = 0.0;   // solid (pre-melt) viscosity [Pa·s]
    double premelt_shear    = 0.0;   // solid (pre-melt) shear modulus [Pa]
    double liquid_viscosity = 0.0;   // viscosity if fully molten at this T [Pa·s]
};

// -------------------------------------------------------------------------------
// c_PartialMeltResult: what every partial-melt model reports.
// -------------------------------------------------------------------------------
struct c_PartialMeltResult {
    double melt_fraction          = 0.0;   // volumetric melt fraction φ [m^3/m^3]
    double postmelt_viscosity     = 0.0;   // post-melt viscosity [Pa·s]
    double postmelt_shear_modulus = 0.0;   // post-melt shear modulus [Pa]
};

// -------------------------------------------------------------------------------
// c_PartialMeltBase
// -------------------------------------------------------------------------------
class c_PartialMeltBase : public c_PhysicsBase {
public:
    c_PartialMeltBase() = default;

    explicit c_PartialMeltBase(const std::string& model_name) : c_PhysicsBase(model_name) {}

    // Construct with the material's melt envelope shared by every model.
    c_PartialMeltBase(
            const std::string& model_name,
            double solidus,
            double liquidus,
            double liquid_shear)
        : c_PhysicsBase(model_name),
          p_solidus(solidus),
          p_liquidus(liquidus),
          p_liquid_shear(liquid_shear) {}

    ~c_PartialMeltBase() override = default;

    // -----------------------------------------------------------------------
    // Shared material constants (the melt envelope).
    // -----------------------------------------------------------------------
    double get_solidus()      const noexcept { return this->p_solidus; }
    double get_liquidus()     const noexcept { return this->p_liquidus; }
    double get_liquid_shear() const noexcept { return this->p_liquid_shear; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_PhysicsBase::append_config_entries(out);
        out.push_back(c_config_double("solidus_k", this->p_solidus));
        out.push_back(c_config_double("liquidus_k", this->p_liquidus));
        out.push_back(c_config_double("liquid_shear_pa", this->p_liquid_shear));
    }

    // -----------------------------------------------------------------------
    // Volumetric melt fraction, model-independent:
    //   φ = clip((T − T_solidus) / (T_liquidus − T_solidus), 0, 1)
    // A non-positive envelope (solidus >= liquidus) gives 0, fully solid.
    // -----------------------------------------------------------------------
    double calc_melt_fraction(double temperature) const noexcept {
        const double denom = this->p_liquidus - this->p_solidus;
        if (denom <= TidalPyConstants::d_EPS) { return 0.0; }
        double phi = (temperature - this->p_solidus) / denom;
        if (phi < 0.0) { phi = 0.0; }
        if (phi > 1.0) { phi = 1.0; }
        return phi;
    }

    // -----------------------------------------------------------------------
    // Partial melt (pure virtual). Maps the pre-melt state to the post-melt
    // viscosity and shear modulus, also returning the melt fraction. The shear
    // modulus and viscosity are floored at the liquid limits by the models.
    // -----------------------------------------------------------------------
    virtual c_PartialMeltResult calc_partial_melt(const c_PartialMeltInputs& inputs) const = 0;

    // -----------------------------------------------------------------------
    // Vectorized partial melt: vary temperature and the pre-melt strengths
    // element-wise at the constant liquid viscosity. All vectors must share
    // length N; out_results is resized to N and a mismatch throws
    // std::invalid_argument. This is the radial sweep, one entry per slice.
    // -----------------------------------------------------------------------
    void calc_partial_melt_vectorize(
            const std::vector<double>& temperature,
            const std::vector<double>& premelt_viscosity,
            const std::vector<double>& premelt_shear,
            double liquid_viscosity,
            std::vector<c_PartialMeltResult>& out_results) const {
        const std::size_t n = temperature.size();
        if (premelt_viscosity.size() != n || premelt_shear.size() != n) {
            throw std::invalid_argument(
                "TidalPy: calc_partial_melt_vectorize — temperature, premelt_viscosity, "
                "and premelt_shear vectors must have the same length");
        }
        out_results.resize(n);
        c_PartialMeltInputs inputs;
        inputs.liquid_viscosity = liquid_viscosity;
        for (std::size_t i = 0; i < n; ++i) {
            inputs.temperature     = temperature[i];
            inputs.premelt_viscosity = premelt_viscosity[i];
            inputs.premelt_shear     = premelt_shear[i];
            out_results[i] = this->calc_partial_melt(inputs);
        }
    }

protected:
    double p_solidus  = 1600.0;  // [K]
    double p_liquidus = 2000.0;  // [K]
    double p_liquid_shear = 1.0e-5;  // [Pa]
};

} // namespace tidalpy
