#pragma once
/* TidalPy's stellar luminosity models. Solar anchors come from TidalPyConstants.
 *
 * References
 * ----------
 * - Cuntz and Wang (2018), doi:10.3847/2515-5172/aaaa67 - low-mass mass-luminosity polynomial exponent.
 * - Wikipedia mass-luminosity relation (piecewise main-sequence scaling) for the high/low-mass regimes.
 *
 * Binary payload: model name then the model's doubles, through the shared c_PhysicsBase helpers. The
 * layer observer pointer is not serialized.
 */

#include <algorithm>
#include <cctype>
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

namespace tidalpy {

// Combined construction parameters; each model reads only the fields it needs.
struct c_LuminosityConfig {
    double luminosity = 0.0;                  // [W]; Fixed model

    // Power-law model: L = Lsun * coeff * (M / Msun)^exponent.
    double power_law_coeff    = 1.0;            // dimensionless prefactor
    double power_law_exponent = 3.5;            // dimensionless exponent (classic main-sequence value)
};

inline double lum_from_fixed(double /*mass*/, double luminosity) noexcept {
    return luminosity;
}

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

inline std::string lum_to_lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
    return text;
}

// Luminosity supplied directly, independent of mass (alias "constant").
class c_FixedLuminosity : public c_LuminosityBase {
public:
    c_FixedLuminosity() : c_LuminosityBase("fixed") {}
    explicit c_FixedLuminosity(const c_LuminosityConfig& config)
        : c_LuminosityBase("fixed"),
          p_luminosity(config.luminosity) {}
    ~c_FixedLuminosity() override = default;

    double get_luminosity() const noexcept { return this->p_luminosity; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_LuminosityBase::append_config_entries(out);
        out.push_back(c_config_double("luminosity_w", this->p_luminosity));
    }

    double calc_luminosity(double mass) const override {
        return lum_from_fixed(mass, this->p_luminosity);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::FixedLuminosity), {this->p_luminosity});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 1);
        this->p_luminosity = params[0];
    }

protected:
    double p_luminosity = 0.0;
};

// Piecewise main-sequence L(M) (aliases "cuntz_wang", "cw").
class c_MassToLuminosity : public c_LuminosityBase {
public:
    c_MassToLuminosity() : c_LuminosityBase("mass_to_luminosity") {}
    explicit c_MassToLuminosity(const c_LuminosityConfig& /*config*/)
        : c_LuminosityBase("mass_to_luminosity") {}
    ~c_MassToLuminosity() override = default;

    double calc_luminosity(double mass) const override {
        return lum_from_mass(mass);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::MassToLuminosity));
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }
};

// Single power law L = Lsun * coeff * (M/Msun)^p (alias "power_law").
class c_PowerLawLuminosity : public c_LuminosityBase {
public:
    c_PowerLawLuminosity() : c_LuminosityBase("power_law") {}
    explicit c_PowerLawLuminosity(const c_LuminosityConfig& config)
        : c_LuminosityBase("power_law"),
          p_coeff(config.power_law_coeff),
          p_exponent(config.power_law_exponent) {}
    ~c_PowerLawLuminosity() override = default;

    double get_coeff()    const noexcept { return this->p_coeff; }
    double get_exponent() const noexcept { return this->p_exponent; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_LuminosityBase::append_config_entries(out);
        out.push_back(c_config_double("power_law_coeff", this->p_coeff));
        out.push_back(c_config_double("power_law_exponent", this->p_exponent));
    }

    double calc_luminosity(double mass) const override {
        return lum_from_power_law(mass, this->p_coeff, this->p_exponent);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::PowerLawLuminosity),
            {this->p_coeff, this->p_exponent});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 2);
        this->p_coeff    = params[0];
        this->p_exponent = params[1];
    }

protected:
    double p_coeff    = 1.0;
    double p_exponent = 3.5;
};

// One value per model, so c_find_luminosity dispatches without string comparisons.
enum class c_LuminosityModel : uint8_t {
    Fixed            = 0,
    MassToLuminosity = 1,
    PowerLaw         = 2,
};

// Model names are matched case-insensitively.
inline c_LuminosityModel c_luminosity_model_from_name(const std::string& model_name) {
    const std::string name = lum_to_lower(model_name);

    if (name == "fixed" || name == "constant") { return c_LuminosityModel::Fixed; }
    if (name == "mass_to_luminosity" || name == "cuntz_wang" || name == "cw") {
        return c_LuminosityModel::MassToLuminosity;
    }
    if (name == "power_law" || name == "powerlaw") { return c_LuminosityModel::PowerLaw; }

    throw std::invalid_argument("TidalPy: unknown luminosity model name '" + model_name + "'");
}

// The canonical C++ factory: worlds, binary reconstruction, and the Cython wrapper all route here.
inline std::unique_ptr<c_LuminosityBase> c_find_luminosity(
        c_LuminosityModel model, const c_LuminosityConfig& config) {
    switch (model) {
        case c_LuminosityModel::Fixed:            return std::make_unique<c_FixedLuminosity>(config);
        case c_LuminosityModel::MassToLuminosity: return std::make_unique<c_MassToLuminosity>(config);
        case c_LuminosityModel::PowerLaw:         return std::make_unique<c_PowerLawLuminosity>(config);
    }
    throw std::invalid_argument("TidalPy: unrecognised c_LuminosityModel enum value");
}

inline std::unique_ptr<c_LuminosityBase> c_find_luminosity(
        const std::string& model_name, const c_LuminosityConfig& config) {
    return c_find_luminosity(c_luminosity_model_from_name(model_name), config);
}

// The class id is peeked without consuming the header so the default-constructed model restores itself.
inline std::unique_ptr<c_LuminosityBase> c_luminosity_from_binary(std::istream& in, bool force = false) {
    const std::streampos start = in.tellg();
    const c_BinaryHeader header = read_binary_header(in);
    in.seekg(start);

    std::unique_ptr<c_LuminosityBase> model;
    switch (static_cast<BinaryClassID>(header.class_id)) {
        case BinaryClassID::FixedLuminosity:    model = std::make_unique<c_FixedLuminosity>();    break;
        case BinaryClassID::MassToLuminosity:   model = std::make_unique<c_MassToLuminosity>();   break;
        case BinaryClassID::PowerLawLuminosity: model = std::make_unique<c_PowerLawLuminosity>(); break;
        default:
            throw std::runtime_error("TidalPy: unknown luminosity class id in binary stream");
    }
    model->read_binary(in, force);
    return model;
}

}  // namespace tidalpy
