#pragma once
/*
 * luminosity_base_.hpp - c_LuminosityBase: abstract base for all TidalPy stellar luminosity models.
 *
 * Inherits c_PhysicsBase (Utilities_x/classes_x/physics_base_.hpp). A luminosity model is a
 * world-level (star) physics model, so its inherited layer observer pointer stays null.
 *
 * The base provides the model-independent Stefan-Boltzmann conversions between a star's effective
 * surface temperature and its luminosity (both need the stellar radius), plus the abstract
 * mass -> luminosity relation that each concrete model implements:
 *   calc_luminosity(mass)                       -> double [W]   (pure virtual; model specific)
 *   calc_luminosity_from_temperature(T, radius) -> double [W]   (L = 4 pi R^2 sigma T^4)
 *   calc_temperature_from_luminosity(L, radius) -> double [K]   (T = (L / (4 pi R^2 sigma))^(1/4))
 *   calc_effective_temperature(mass, radius)    -> double [K]   (mass -> L -> T)
 *
 * The concrete models (Fixed, MassToLuminosity, PowerLaw) live in luminosity_.hpp.
 *
 * All calc_* methods are const and operate in MKS units: mass [kg], radius [m], temperature [K],
 * luminosity [W]. The Stefan-Boltzmann constant is read from the shared TidalPy config
 * (tidalpy_config_ptr->d_SBC); a null config pointer or a non-positive/degenerate input yields NaN.
 */

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "physics_base_.hpp"
#include "constants_.hpp"

namespace tidalpy {

// -------------------------------------------------------------------------------
// c_LuminosityBase
// -------------------------------------------------------------------------------
class c_LuminosityBase : public c_PhysicsBase {
public:
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------
    c_LuminosityBase() = default;

    explicit c_LuminosityBase(const std::string& model_name) : c_PhysicsBase(model_name) {}

    ~c_LuminosityBase() override = default;

    // -----------------------------------------------------------------------
    // Luminosity from stellar mass (pure virtual) [W]
    //
    // Each model overrides this with its mass -> luminosity relation. The Fixed
    // model ignores the mass and returns its stored luminosity.
    //
    // Parameters
    // ----------
    // mass_kg : stellar mass [kg]
    //
    // Returns
    // -------
    // double : stellar luminosity [W]
    //
    // Assumptions
    // -----------
    // - Main-sequence mass-luminosity scaling (model specific).
    // - All inputs and outputs are MKS.
    // -----------------------------------------------------------------------
    virtual double calc_luminosity(double mass_kg) const = 0;

    // -----------------------------------------------------------------------
    // Stefan-Boltzmann luminosity [W]: L = 4 pi R^2 sigma T^4.
    // Returns NaN for a non-positive temperature/radius or a null config pointer.
    //
    // Assumptions
    // -----------
    // The star radiates as an ideal gray body at the given effective temperature.
    // -----------------------------------------------------------------------
    double calc_luminosity_from_temperature(double temperature_k, double radius_m) const noexcept {
        if (temperature_k <= 0.0 || radius_m <= 0.0 || tidalpy_config_ptr == nullptr) {
            return TidalPyConstants::d_NAN;
        }
        const double area  = 4.0 * TidalPyConstants::d_PI * radius_m * radius_m;
        return area * tidalpy_config_ptr->d_SBC * temperature_k * temperature_k * temperature_k * temperature_k;
    }

    // -----------------------------------------------------------------------
    // Effective temperature [K] from luminosity: T = (L / (4 pi R^2 sigma))^(1/4).
    // Returns NaN for a non-positive luminosity/radius or a null config pointer.
    // -----------------------------------------------------------------------
    double calc_temperature_from_luminosity(double luminosity_w, double radius_m) const noexcept {
        if (luminosity_w <= 0.0 || radius_m <= 0.0 || tidalpy_config_ptr == nullptr) {
            return TidalPyConstants::d_NAN;
        }
        const double area  = 4.0 * TidalPyConstants::d_PI * radius_m * radius_m;
        const double denom = area * tidalpy_config_ptr->d_SBC;
        if (std::abs(denom) <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        return std::pow(luminosity_w / denom, 0.25);
    }

    // -----------------------------------------------------------------------
    // Effective temperature [K] derived from the stellar mass: mass -> L -> T.
    // -----------------------------------------------------------------------
    double calc_effective_temperature(double mass_kg, double radius_m) const noexcept {
        return this->calc_temperature_from_luminosity(this->calc_luminosity(mass_kg), radius_m);
    }

    // -----------------------------------------------------------------------
    // Vectorized luminosity - vary mass.
    //
    // Evaluates calc_luminosity over each mass, writing into out_luminosity
    // (resized to the mass vector length).
    // -----------------------------------------------------------------------
    void calc_luminosity_vectorize_mass(
            const std::vector<double>& mass_kg,
            std::vector<double>& out_luminosity) const {
        const std::size_t num_masses = mass_kg.size();
        out_luminosity.resize(num_masses);
        for (std::size_t i = 0; i < num_masses; ++i) {
            out_luminosity[i] = this->calc_luminosity(mass_kg[i]);
        }
    }
};

} // namespace tidalpy
