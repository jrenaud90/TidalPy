#pragma once
/* Abstract base for TidalPy stellar luminosity models. Concrete models live in luminosity_.hpp.
 *
 * A luminosity model is world-level (a star), so its inherited layer observer pointer stays null. The
 * base carries the model-independent Stefan-Boltzmann conversions between effective surface temperature
 * and luminosity (the c_stefan_boltzmann_* functions); a null config pointer or a non-positive input yields NaN.
 */

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "physics_base_.hpp"
#include "constants_.hpp"

namespace tidalpy {

// 4 pi R^2 sigma [W K^-4]: the luminosity per T^4 of a sphere of radius R [m] radiating as an ideal gray body, with
// sigma the config's Stefan-Boltzmann constant. NaN when the config pointer is unwired.
inline double c_stefan_boltzmann_factor(double radius) noexcept {
    const double sigma = (tidalpy_config_ptr != nullptr) ? tidalpy_config_ptr->d_SBC : TidalPyConstants::d_NAN;
    return 4.0 * TidalPyConstants::d_PI * radius * radius * sigma;
}

// Stefan-Boltzmann luminosity L = 4 pi R^2 sigma T^4 [W] from the effective temperature [K] and radius [m]. No input
// checks: each caller chooses its own result for a non-positive input.
inline double c_stefan_boltzmann_luminosity(double temperature, double radius) noexcept {
    return c_stefan_boltzmann_factor(radius) * temperature * temperature * temperature * temperature;
}

// Effective temperature T = (L / (4 pi R^2 sigma))^(1/4) [K] from the luminosity [W] and radius [m], the inverse of
// c_stefan_boltzmann_luminosity. No input checks, as there.
inline double c_stefan_boltzmann_temperature(double luminosity, double radius) noexcept {
    return std::pow(luminosity / c_stefan_boltzmann_factor(radius), 0.25);
}

class c_LuminosityBase : public c_PhysicsBase {
public:
    c_LuminosityBase() = default;

    explicit c_LuminosityBase(const std::string& model_name) : c_PhysicsBase(model_name) {}

    ~c_LuminosityBase() override = default;

    // Stellar luminosity [W] from mass [kg]; the Fixed model ignores the mass.
    virtual double calc_luminosity(double mass) const = 0;

    // L = 4 pi R^2 sigma T^4; assumes the star radiates as an ideal gray body.
    double calc_luminosity_from_temperature(double temperature, double radius) const noexcept {
        if (temperature <= 0.0 || radius <= 0.0 || tidalpy_config_ptr == nullptr) {
            return TidalPyConstants::d_NAN;
        }
        return c_stefan_boltzmann_luminosity(temperature, radius);
    }

    // T = (L / (4 pi R^2 sigma))^(1/4).
    double calc_temperature_from_luminosity(double luminosity, double radius) const noexcept {
        if (luminosity <= 0.0 || radius <= 0.0 || tidalpy_config_ptr == nullptr) {
            return TidalPyConstants::d_NAN;
        }
        if (std::abs(c_stefan_boltzmann_factor(radius)) <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        return c_stefan_boltzmann_temperature(luminosity, radius);
    }

    // mass -> L -> T.
    double calc_effective_temperature(double mass, double radius) const noexcept {
        return this->calc_temperature_from_luminosity(this->calc_luminosity(mass), radius);
    }

    // Over mass.
    void calc_luminosity_vectorize_mass(
            const std::vector<double>& mass,
            std::vector<double>& out_luminosity) const {
        const std::size_t num_masses = mass.size();
        out_luminosity.resize(num_masses);
        for (std::size_t i = 0; i < num_masses; ++i) {
            out_luminosity[i] = this->calc_luminosity(mass[i]);
        }
    }
};

} // namespace tidalpy
