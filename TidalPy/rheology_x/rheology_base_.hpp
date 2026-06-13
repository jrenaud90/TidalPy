#pragma once
/*
 * rheology_base_.hpp - c_RheologyBase: abstract base for all TidalPy rheology models.
 *
 * Inherits c_PhysicsBase (Utilities_x/classes_x/physics_base_.hpp).
 *
 * Defines the abstract interface that every rheology model must satisfy:
 *   calc_complex_modulus(modulus, viscosity, frequency) -> c_ComplexModulus
 *
 * The seven currently implemented models (Elastic, Viscous, Voigt, Maxwell, Burgers, Andrade,
 * Sundberg) are implemented in rheology_.hpp.
 *
 * All calc_* methods are const.
 */

#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

#include "physics_base_.hpp"

namespace tidalpy {

// Complex modulus alias [Pa] (real = storage / in-phase, imag = loss / out-of-phase).
using c_ComplexModulus = std::complex<double>;

// Complex compliance alias [Pa^-1] (inverse of complex modulus).
using c_ComplexCompliance = std::complex<double>;

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
    // Complex (shear/bulk) modulus (pure virtual) [Pa]
    //
    // Each model overrides this to implement its constitutive law and
    // return the complex modulus directly. c_PhysicsLayer calls this method to
    // obtain the frequency-dependent complex shear and bulk moduli.
    //
    // Parameters
    // ----------
    // modulus_pa      : unrelaxed (shear/bulk) modulus [Pa]
    // viscosity_pas   : reference dynamic (shear/bulk) viscosity [Pa·s]
    // frequency_rad_s : tidal forcing frequency [rad/s]
    //
    // Returns
    // -------
    // c_ComplexModulus [Pa]:
    //   real part  — storage (in-phase) modulus
    //   imag part  — loss (out-of-phase) modulus (positive = energy loss)
    //
    // Assumptions
    // -----------
    // - Single forcing frequency (linear viscoelastic regime).
    // - All inputs are reference (background) values at the layer mid-point.
    // -----------------------------------------------------------------------
    virtual c_ComplexModulus calc_complex_modulus(
            double modulus_pa,
            double viscosity_pas,
            double frequency_rad_s) const {
        // Base constructor will act like a off model.
        return c_ComplexModulus(modulus_pa, 0.0);

    };

    // -----------------------------------------------------------------------
    // Vectorized complex modulus — vary (modulus, viscosity) at one frequency.
    //
    // Evaluates calc_complex_modulus element-wise over the (viscosity, modulus)
    // pairs at a single constant frequency.  viscosity_pas and modulus_pa must
    // have the same length N; out_complex_modulus is resized to N.
    //
    // Throws std::invalid_argument if the two input vectors differ in length.
    // -----------------------------------------------------------------------
    void calc_complex_modulus_vectorize_modulus(
            const std::vector<double>& modulus_pa,
            const std::vector<double>& viscosity_pas,
            double frequency_rad_s,
            std::vector<c_ComplexModulus>& out_complex_modulus) const {
        if (viscosity_pas.size() != modulus_pa.size()) {
            throw std::invalid_argument(
                "TidalPy: calc_complex_modulus_vectorize_modulus — viscosity and "
                "modulus vectors must have the same length");
        }
        const std::size_t n = modulus_pa.size();
        out_complex_modulus.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            out_complex_modulus[i] =
                this->calc_complex_modulus(modulus_pa[i], viscosity_pas[i], frequency_rad_s);
        }
    }

    // -----------------------------------------------------------------------
    // Vectorized complex modulus — vary frequency at constant modulus/viscosity.
    //
    // Evaluates calc_complex_modulus over each frequency at the single constant
    // viscosity and modulus.  out_complex_modulus is resized to the frequency
    // vector length.
    // -----------------------------------------------------------------------
    void calc_complex_modulus_vectorize_frequency(
            double modulus_pa,    
            double viscosity_pas,
            const std::vector<double>& frequency_rad_s,
            std::vector<c_ComplexModulus>& out_complex_modulus) const {
        const std::size_t n = frequency_rad_s.size();
        out_complex_modulus.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            out_complex_modulus[i] =
                this->calc_complex_modulus(modulus_pa, viscosity_pas, frequency_rad_s[i]);
        }
    }

    // -----------------------------------------------------------------------
    // Vectorized complex modulus — vary all of (modulus, viscosity, frequency).
    //
    // Evaluates calc_complex_modulus element-wise over the three input vectors,
    // which must all have the same length N.  out_complex_modulus is resized to N.
    //
    // Throws std::invalid_argument if the input vectors differ in length.
    // -----------------------------------------------------------------------
    void calc_complex_modulus_vectorize_all(
            const std::vector<double>& modulus_pa,
            const std::vector<double>& viscosity_pas,
            const std::vector<double>& frequency_rad_s,
            std::vector<c_ComplexModulus>& out_complex_modulus) const {
        if (viscosity_pas.size() != modulus_pa.size() ||
            viscosity_pas.size() != frequency_rad_s.size()) {
            throw std::invalid_argument(
                "TidalPy: calc_complex_modulus_vectorize_all — viscosity, modulus, "
                "and frequency vectors must all have the same length");
        }
        const std::size_t n = modulus_pa.size();
        out_complex_modulus.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            out_complex_modulus[i] =
                this->calc_complex_modulus(modulus_pa[i], viscosity_pas[i], frequency_rad_s[i]);
        }
    }
};

} // namespace tidalpy
