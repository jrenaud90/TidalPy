#pragma once
/*
 * rheology_base_.hpp: c_RheologyBase, the abstract base for TidalPy rheology models, derived from
 * c_PhysicsBase. Every model implements calc_complex_modulus(modulus, viscosity, frequency). The seven
 * models (Elastic, Viscous, Voigt, Maxwell, Burgers, Andrade, Sundberg) are in rheology_.hpp. All
 * calc_* methods are const.
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
    // Complex (shear or bulk) modulus [Pa] from the unrelaxed modulus [Pa], the
    // reference dynamic viscosity [Pa s], and the forcing frequency [rad s-1].
    // The real part is the storage (in-phase) modulus and the imaginary part the
    // loss (out-of-phase) modulus, positive for energy loss. The base
    // implementation is elastic; every model overrides it. c_PhysicsLayer calls
    // this for the frequency-dependent complex shear and bulk moduli.
    //
    // Assumes a linear viscoelastic regime at a single forcing frequency, with
    // reference (background) inputs at the layer mid-point.
    // -----------------------------------------------------------------------
    virtual c_ComplexModulus calc_complex_modulus(
            double modulus,
            double viscosity,
            double frequency) const {
        return c_ComplexModulus(modulus, 0.0);

    };

    // -----------------------------------------------------------------------
    // Vectorized: vary (modulus, viscosity) element-wise at one frequency. Both
    // vectors must share length N; a mismatch throws std::invalid_argument.
    // -----------------------------------------------------------------------
    void calc_complex_modulus_vectorize_modulus(
            const std::vector<double>& modulus,
            const std::vector<double>& viscosity,
            double frequency,
            std::vector<c_ComplexModulus>& out_complex_modulus) const {
        if (viscosity.size() != modulus.size()) {
            throw std::invalid_argument(
                "TidalPy: calc_complex_modulus_vectorize_modulus — viscosity and "
                "modulus vectors must have the same length");
        }
        const std::size_t n = modulus.size();
        out_complex_modulus.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            out_complex_modulus[i] =
                this->calc_complex_modulus(modulus[i], viscosity[i], frequency);
        }
    }

    // -----------------------------------------------------------------------
    // Vectorized: vary frequency at constant modulus and viscosity.
    // out_complex_modulus is resized to the frequency vector length.
    // -----------------------------------------------------------------------
    void calc_complex_modulus_vectorize_frequency(
            double modulus,    
            double viscosity,
            const std::vector<double>& frequency,
            std::vector<c_ComplexModulus>& out_complex_modulus) const {
        const std::size_t n = frequency.size();
        out_complex_modulus.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            out_complex_modulus[i] =
                this->calc_complex_modulus(modulus, viscosity, frequency[i]);
        }
    }

    // -----------------------------------------------------------------------
    // Vectorized: vary modulus, viscosity, and frequency element-wise. All three
    // vectors must share length N; a mismatch throws std::invalid_argument.
    // -----------------------------------------------------------------------
    void calc_complex_modulus_vectorize_all(
            const std::vector<double>& modulus,
            const std::vector<double>& viscosity,
            const std::vector<double>& frequency,
            std::vector<c_ComplexModulus>& out_complex_modulus) const {
        if (viscosity.size() != modulus.size() ||
            viscosity.size() != frequency.size()) {
            throw std::invalid_argument(
                "TidalPy: calc_complex_modulus_vectorize_all — viscosity, modulus, "
                "and frequency vectors must all have the same length");
        }
        const std::size_t n = modulus.size();
        out_complex_modulus.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            out_complex_modulus[i] =
                this->calc_complex_modulus(modulus[i], viscosity[i], frequency[i]);
        }
    }
};

} // namespace tidalpy
