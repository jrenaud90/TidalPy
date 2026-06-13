#pragma once
/*
 * rheology_base_.hpp — c_RheologyBase: abstract base for all TidalPy rheology models.
 *
 * Inherits c_PhysicsBase (Utilities_x/classes_x/physics_base_.hpp).
 *
 * Defines the abstract interface that every rheology model must satisfy:
 *   calc_complex_compliance(viscosity, modulus, frequency) -> c_ComplexCompliance
 *
 * The seven concrete models (Elastic, Viscous, Voigt, Maxwell, Burgers, Andrade,
 * Sundberg) are implemented in rheology_.hpp (Phase 5).  This header exists
 * separately so that c_PhysicsLayer (Phase 2) can hold unique_ptr<c_RheologyBase>
 * without depending on Phase 5.
 *
 * All calc_* methods are const.
 */

#include <complex>
#include <string>

#include "physics_base_.hpp"

namespace tidalpy {

// Complex compliance alias [Pa^-1].
using c_ComplexModulus = std::complex<double>;

// -------------------------------------------------------------------------------
// c_RheologyBase
// -------------------------------------------------------------------------------
class c_RheologyBase : public c_PhysicsBase {
public:
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------
    c_RheologyBase() = default;

    explicit c_RheologyBase(const std::string& model_name) : c_PhysicsBase(model_name) {}

    ~c_RheologyBase() override = default;

    // -----------------------------------------------------------------------
    // Complex modulus (pure virtual)
    //
    // Each concrete model overrides this to implement its constitutive law.
    //
    // Parameters
    // ----------
    // viscosity_pas   : reference dynamic (shear/bulk) viscosity [Pa·s]
    // modulus_pa      : unrelaxed (shear/bulk) modulus [Pa]
    // frequency_rad_s : tidal forcing frequency [rad/s]
    //
    // Returns
    // -------
    // c_ComplexShear [Pa]:
    //   real part  — elastic (in-phase) shear
    //   imag part  — viscous (out-of-phase) shear (positive = energy loss)
    //
    // Assumptions
    // -----------
    // - Single forcing frequency (linear viscoelastic regime).
    // - All inputs are reference (background) values at the layer mid-point.
    // -----------------------------------------------------------------------
    virtual c_ComplexModulus calc_complex_modulus(
            double viscosity_pas,
            double modulus_pa,
            double frequency_rad_s) const = 0;
};

} // namespace tidalpy
