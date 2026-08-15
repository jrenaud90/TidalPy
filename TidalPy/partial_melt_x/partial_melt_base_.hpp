#pragma once
/*
 * partial_melt_base_.hpp — c_PartialMeltBase: abstract base for all TidalPy
 * partial-melt models.
 *
 * Inherits c_PhysicsBase (Utilities_x/classes_x/physics_base_.hpp).
 *
 * A partial-melt model maps a material's pre-melt (solid) viscosity and shear
 * modulus, plus its temperature and melt fraction, to the POST-melt viscosity and
 * shear modulus, i.e. it applies melt weakening (and, for some models, the bare
 * temperature dependence of the strength). The melt fraction itself is a
 * model-independent function of temperature and the material's solidus/liquidus,
 * provided here on the base class.
 *
 * These quantities are frequency-INDEPENDENT: they depend only on the (already
 * solved) temperature/pressure state, so the world-level pipeline caches them once
 * after the EOS solve and only the downstream rheology (complex modulus) step is
 * recomputed per forcing frequency. See the "Complex-Moduli Pipeline" design note
 * in the planning doc.
 *
 * The three implemented models (Off, Spohn, Henning) live in partial_melt_.hpp.
 *
 * All calc_* methods are const and operate in MKS units.
 *
 * References
 * ----------
 * - Fischer and Spohn (1990) — temperature-based melt viscosity/shear law.
 * - Henning (2009, 2010); Renaud and Henning (2018) — three-regime melt weakening.
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
// c_PartialMeltInputs — the per-evaluation state passed to a partial-melt model.
// (Material constants — solidus/liquidus, liquid shear, model parameters — live
// on the model object; only the varying state is passed here.) All MKS.
// -------------------------------------------------------------------------------
struct c_PartialMeltInputs {
    double temperature_k      = 0.0;   // local temperature [K]
    double premelt_viscosity  = 0.0;   // solid (pre-melt) viscosity [Pa·s]
    double premelt_shear      = 0.0;   // solid (pre-melt) shear modulus [Pa]
    double liquid_viscosity   = 0.0;   // viscosity if fully molten at this T [Pa·s]
};

// -------------------------------------------------------------------------------
// c_PartialMeltResult — what every partial-melt model reports.
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
            double solidus_k,
            double liquidus_k,
            double liquid_shear_pa)
        : c_PhysicsBase(model_name),
          p_solidus_k(solidus_k),
          p_liquidus_k(liquidus_k),
          p_liquid_shear_pa(liquid_shear_pa) {}

    ~c_PartialMeltBase() override = default;

    // -----------------------------------------------------------------------
    // Shared material constants (the melt envelope).
    // -----------------------------------------------------------------------
    double get_solidus()      const noexcept { return this->p_solidus_k; }
    double get_liquidus()     const noexcept { return this->p_liquidus_k; }
    double get_liquid_shear() const noexcept { return this->p_liquid_shear_pa; }

    // -----------------------------------------------------------------------
    // Volumetric melt fraction φ ∈ [0, 1] from temperature (model-independent):
    //
    //   φ = clip( (T − T_solidus) / (T_liquidus − T_solidus), 0, 1 )
    //
    // A non-positive (solidus >= liquidus) envelope yields 0 (fully solid).
    // -----------------------------------------------------------------------
    double calc_melt_fraction(double temperature_k) const noexcept {
        const double denom = this->p_liquidus_k - this->p_solidus_k;
        if (denom <= TidalPyConstants::d_EPS) { return 0.0; }
        double phi = (temperature_k - this->p_solidus_k) / denom;
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
    // Vectorized partial melt — vary temperature and the pre-melt strengths
    // element-wise at the (constant) liquid viscosity. All vectors must be the
    // same length N; out_results is resized to N. Throws std::invalid_argument on
    // a length mismatch. This is the primary radial sweep (one entry per slice).
    // -----------------------------------------------------------------------
    void calc_partial_melt_vectorize(
            const std::vector<double>& temperature_k,
            const std::vector<double>& premelt_viscosity,
            const std::vector<double>& premelt_shear,
            double liquid_viscosity,
            std::vector<c_PartialMeltResult>& out_results) const {
        const std::size_t n = temperature_k.size();
        if (premelt_viscosity.size() != n || premelt_shear.size() != n) {
            throw std::invalid_argument(
                "TidalPy: calc_partial_melt_vectorize — temperature, premelt_viscosity, "
                "and premelt_shear vectors must have the same length");
        }
        out_results.resize(n);
        c_PartialMeltInputs inputs;
        inputs.liquid_viscosity = liquid_viscosity;
        for (std::size_t i = 0; i < n; ++i) {
            inputs.temperature_k     = temperature_k[i];
            inputs.premelt_viscosity = premelt_viscosity[i];
            inputs.premelt_shear     = premelt_shear[i];
            out_results[i] = this->calc_partial_melt(inputs);
        }
    }

protected:
    double p_solidus_k       = 1600.0;  // [K]
    double p_liquidus_k      = 2000.0;  // [K]
    double p_liquid_shear_pa = 1.0e-5;  // [Pa]
};

} // namespace tidalpy
