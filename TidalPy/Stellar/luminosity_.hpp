#pragma once
/* TidalPy's stellar luminosity models. Solar anchors come from TidalPyConstants.
 *
 * Each model declares its parameters in one table (c_SpecModel, spec_model_.hpp), which gives its constructor,
 * validation, config entries, and binary record.
 *
 * References
 * ----------
 * - Cuntz and Wang (2018), doi:10.3847/2515-5172/aaaa67 - low-mass mass-luminosity polynomial exponent.
 * - Wikipedia mass-luminosity relation (piecewise main-sequence scaling) for the high/low-mass regimes.
 */

#include <cmath>
#include <cstdint>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "luminosity_base_.hpp"
#include "constants_.hpp"
#include "model_names_.hpp"
#include "registry_.hpp"
#include "spec_model_.hpp"

namespace tidalpy {

// Piecewise main-sequence relation (Cuntz and Wang 2018).
inline double lum_from_mass(double mass) noexcept {
    const double mass_solar      = TidalPyConstants::d_MASS_SOLAR;
    const double luminosity_solar = TidalPyConstants::d_LUMINOSITY_SOLAR;
    if (mass <= 0.0 || mass_solar <= 0.0) {
        return TidalPyConstants::d_NAN;
    }
    const double mass_ratio = mass / mass_solar;

    if (mass_ratio < 0.2) {
        return luminosity_solar * 0.23 * std::pow(mass_ratio, 2.3);
    }
    if (mass_ratio < 0.85) {
        // Cuntz and Wang (2018) polynomial exponent in the mass ratio.
        const double exponent =
            -141.7 * std::pow(mass_ratio, 4.0)
            + 232.4 * std::pow(mass_ratio, 3.0)
            - 129.1 * std::pow(mass_ratio, 2.0)
            + 33.29 * mass_ratio
            + 0.215;
        return luminosity_solar * std::pow(mass_ratio, exponent);
    }
    if (mass_ratio < 2.0) {
        return luminosity_solar * std::pow(mass_ratio, 4.0);
    }
    // The linear branch takes over near where it meets the 1.4 M^3.5 branch: at 55 Msun it is 1.9 percent above it.
    // The joints at 0.2 (-18.7 percent) and 2 Msun (-1.0 percent) step too; luminosity.md tabulates them.
    if (mass_ratio < 55.0) {
        return luminosity_solar * 1.4 * std::pow(mass_ratio, 3.5);
    }
    return luminosity_solar * 3.2e4 * mass_ratio;
}

inline double lum_from_power_law(double mass, double coeff, double exponent) noexcept {
    const double mass_solar       = TidalPyConstants::d_MASS_SOLAR;
    const double luminosity_solar = TidalPyConstants::d_LUMINOSITY_SOLAR;
    if (mass <= 0.0 || mass_solar <= 0.0) {
        return TidalPyConstants::d_NAN;
    }
    return luminosity_solar * coeff * std::pow(mass / mass_solar, exponent);
}

// Luminosity supplied directly, independent of mass (alias "constant").
class c_FixedLuminosity final : public c_SpecModel<c_FixedLuminosity, c_LuminosityBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::FixedLuminosity;

    static const std::vector<c_ParamSpec<c_FixedLuminosity>>& parameter_specs() {
        using Self = c_FixedLuminosity;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"luminosity", "luminosity_w", &Self::p_luminosity, 0.0, c_ParamBounds::NonNegative,
             "The luminosity reported at every mass [W]."},
        };
        return specs;
    }

    c_FixedLuminosity() : c_FixedLuminosity(c_ParamMap{}) {}
    explicit c_FixedLuminosity(const c_ParamMap& params) : c_SpecModel("fixed") { this->p_initialize(params); }

    double calc_luminosity(double /*mass*/) const override { return this->p_luminosity; }

protected:
    double p_luminosity = 0.0;
};

// Piecewise main-sequence L(M) (aliases "cuntz_wang", "cw").
class c_MassToLuminosity final : public c_SpecModel<c_MassToLuminosity, c_LuminosityBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::MassToLuminosity;

    static const std::vector<c_ParamSpec<c_MassToLuminosity>>& parameter_specs() {
        static const std::vector<c_ParamSpec<c_MassToLuminosity>> specs = {};
        return specs;
    }

    c_MassToLuminosity() : c_MassToLuminosity(c_ParamMap{}) {}
    explicit c_MassToLuminosity(const c_ParamMap& params) : c_SpecModel("mass_to_luminosity") {
        this->p_initialize(params);
    }

    double calc_luminosity(double mass) const override { return lum_from_mass(mass); }
};

// Single power law L = Lsun * coeff * (M/Msun)^exponent (alias "powerlaw").
class c_PowerLawLuminosity final : public c_SpecModel<c_PowerLawLuminosity, c_LuminosityBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::PowerLawLuminosity;

    static const std::vector<c_ParamSpec<c_PowerLawLuminosity>>& parameter_specs() {
        using Self = c_PowerLawLuminosity;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"coeff", "power_law_coeff", &Self::p_coeff, 1.0, c_ParamBounds::Positive,
             "Prefactor of the power law, in solar luminosities."},
            {"exponent", "power_law_exponent", &Self::p_exponent, 3.5, c_ParamBounds::Finite,
             "Exponent of the mass ratio; 3.5 is the classic main-sequence value."},
        };
        return specs;
    }

    c_PowerLawLuminosity() : c_PowerLawLuminosity(c_ParamMap{}) {}
    explicit c_PowerLawLuminosity(const c_ParamMap& params) : c_SpecModel("power_law") { this->p_initialize(params); }

    double calc_luminosity(double mass) const override {
        return lum_from_power_law(mass, this->p_coeff, this->p_exponent);
    }

protected:
    double p_coeff    = 1.0;
    double p_exponent = 3.5;
};

inline const c_ModelRegistry<c_LuminosityBase>& c_luminosity_registry() {
    static const c_ModelRegistry<c_LuminosityBase> registry = {
        {{"fixed", "constant"},
         BinaryClassID::FixedLuminosity,    &c_make_entry<c_LuminosityBase, c_FixedLuminosity>},
        {{"mass_to_luminosity", "cuntz_wang", "cw"},
         BinaryClassID::MassToLuminosity,   &c_make_entry<c_LuminosityBase, c_MassToLuminosity>},
        {{"power_law", "powerlaw"},
         BinaryClassID::PowerLawLuminosity, &c_make_entry<c_LuminosityBase, c_PowerLawLuminosity>},
    };
    return registry;
}

// The family's entry points, each one line over the generic registry functions.
inline std::unique_ptr<c_LuminosityBase> c_find_luminosity(const std::string& model_name, const c_ParamMap& params) {
    return c_make_model(c_luminosity_registry(), model_name, params);
}

inline std::unique_ptr<c_LuminosityBase> c_luminosity_from_binary(std::istream& in, bool force = false) {
    return c_model_from_binary(c_luminosity_registry(), in, force);
}

inline std::string c_luminosity_canonical_name(const std::string& model_name) {
    return c_canonical_model_name(c_luminosity_registry(), model_name);
}

inline std::vector<std::string> c_luminosity_model_names() {
    return c_model_names(c_luminosity_registry());
}

}  // namespace tidalpy
