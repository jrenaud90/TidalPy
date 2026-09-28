#pragma once
/* Abstract base for TidalPy rheology models. Concrete models live in rheology_.hpp. */

#include <complex>
#include <string>
#include <vector>

#include "broadcast_.hpp"
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

    // Element-wise over modulus, viscosity, and frequency. Each input holds one value per point or a single value used
    // at every point, so a frequency sweep at fixed material passes one modulus and one viscosity.
    void calc_complex_modulus_vectorize(
            const std::vector<double>& modulus,
            const std::vector<double>& viscosity,
            const std::vector<double>& frequency,
            std::vector<c_ComplexModulus>& out_complex_modulus) const {
        const std::size_t num_points = c_broadcast_length(
            {modulus.size(), viscosity.size(), frequency.size()}, "calc_complex_modulus_vectorize");
        const std::size_t modulus_stride   = c_broadcast_stride(modulus.size());
        const std::size_t viscosity_stride = c_broadcast_stride(viscosity.size());
        const std::size_t frequency_stride = c_broadcast_stride(frequency.size());
        out_complex_modulus.resize(num_points);
        for (std::size_t i = 0; i < num_points; ++i) {
            out_complex_modulus[i] = this->calc_complex_modulus(
                modulus[i * modulus_stride], viscosity[i * viscosity_stride], frequency[i * frequency_stride]);
        }
    }
};

} // namespace tidalpy
