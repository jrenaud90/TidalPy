#pragma once
/* Abstract base for TidalPy rheology models. Concrete models live in rheology_.hpp. */

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

class c_RheologyBase : public c_PhysicsBase {
public:
    c_RheologyBase() = default;

    explicit c_RheologyBase(const std::string& model_name) : c_PhysicsBase(model_name) {}

    ~c_RheologyBase() override = default;

    // Complex (shear or bulk) modulus [Pa] from the unrelaxed modulus [Pa], the reference dynamic
    // viscosity [Pa s], and the forcing frequency [rad s-1]. Real part is storage (in-phase), imaginary
    // part is loss (out-of-phase, positive for energy loss). Assumes a linear viscoelastic regime at a
    // single forcing frequency with reference (background) inputs at the layer mid-point.
    // Pure virtual on purpose: a default body here would make a model that forgets it silently elastic.
    virtual c_ComplexModulus calc_complex_modulus(
            double modulus,
            double viscosity,
            double frequency) const = 0;

    // Element-wise over (modulus, viscosity) at one frequency.
    void calc_complex_modulus_vectorize_modulus(
            const std::vector<double>& modulus,
            const std::vector<double>& viscosity,
            double frequency,
            std::vector<c_ComplexModulus>& out_complex_modulus) const {
        if (viscosity.size() != modulus.size()) {
            throw std::invalid_argument(
                "TidalPy::calc_complex_modulus_vectorize_modulus: Viscosity and "
                "modulus vectors must have the same length");
        }
        const std::size_t n = modulus.size();
        out_complex_modulus.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            out_complex_modulus[i] =
                this->calc_complex_modulus(modulus[i], viscosity[i], frequency);
        }
    }

    // Over frequency at constant modulus and viscosity.
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

    // Element-wise over all three.
    void calc_complex_modulus_vectorize_all(
            const std::vector<double>& modulus,
            const std::vector<double>& viscosity,
            const std::vector<double>& frequency,
            std::vector<c_ComplexModulus>& out_complex_modulus) const {
        if (viscosity.size() != modulus.size() ||
            viscosity.size() != frequency.size()) {
            throw std::invalid_argument(
                "TidalPy::calc_complex_modulus_vectorize_all: viscosity, modulus, "
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
