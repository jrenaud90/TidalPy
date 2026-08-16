#pragma once
/*
 * luminosity_.hpp - Implements TidalPy's stellar luminosity models.
 *
 * Inherits c_LuminosityBase (luminosity_base_.hpp), which itself inherits c_PhysicsBase. Each model
 * implements calc_luminosity(mass) [W]; the base supplies the Stefan-Boltzmann effective-temperature
 * conversions shared by every model.
 *
 * Models (with config aliases handled by the factory):
 *   c_FixedLuminosity    (alias "constant")            - luminosity set directly (mass independent).
 *   c_MassToLuminosity   (alias "cuntz_wang" / "cw")   - piecewise main-sequence L(M) relation.
 *   c_PowerLawLuminosity (alias "power_law")           - single power law L = Lsun * coeff * (M/Msun)^p.
 *
 * All quantities are MKS: mass [kg], radius [m], temperature [K], luminosity [W]. The solar mass and
 * luminosity anchors come from TidalPyConstants (d_MASS_SOLAR, d_LUMINOSITY_SOLAR).
 *
 * References
 * ----------
 * - Cuntz and Wang (2018), doi:10.3847/2515-5172/aaaa67 - low-mass mass-luminosity polynomial exponent.
 * - Wikipedia mass-luminosity relation (piecewise main-sequence scaling) for the high/low-mass regimes.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::<Model> (1001-1003)
 *   payload: model_name length (uint32_t) | model_name bytes | model params (scalar doubles)
 *   All three models serialize via the shared c_PhysicsBase helpers.
 *   The layer observer pointer (p_layer_ptr) is NOT serialized.
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

// -------------------------------------------------------------------------------
// c_LuminosityConfig - combined construction parameters for all models.
// Each model reads only the fields it needs.
// -------------------------------------------------------------------------------
struct c_LuminosityConfig {
    // Fixed model - the luminosity to report regardless of mass.
    double luminosity_w = 0.0;                  // [W]

    // Power-law model - L = Lsun * power_law_coeff * (M / Msun)^power_law_exponent.
    double power_law_coeff    = 1.0;            // dimensionless prefactor
    double power_law_exponent = 3.5;            // dimensionless exponent (classic main-sequence value)
};

// =====================================================================================================================
// Luminosity relations [W]
//
// Free helpers mirroring the validated legacy implementation in TidalPy/stellar/stellar.py
// (luminosity_from_mass), with the low-mass polynomial-exponent branch corrected to use the mass ratio
// M/Msun rather than the raw kilogram mass.
// =====================================================================================================================

// Fixed: report the stored luminosity regardless of mass.
inline double lum_from_fixed(double /*mass_kg*/, double luminosity_w) noexcept {
    return luminosity_w;
}

// Mass-to-luminosity: the piecewise main-sequence relation (Cuntz and Wang 2018).
// mass_ratio = M / Msun.
inline double lum_from_mass(double mass_kg) noexcept {
    const double mass_solar      = TidalPyConstants::d_MASS_SOLAR;
    const double luminosity_solar = TidalPyConstants::d_LUMINOSITY_SOLAR;
    if (mass_kg <= 0.0 || mass_solar <= 0.0) {
        return TidalPyConstants::d_NAN;
    }
    const double mass_ratio = mass_kg / mass_solar;

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
    // The linear branch takes over where it meets the 1.4 M^3.5 branch (~55 Msun);
    // both give ~1.75e6 Lsun there, keeping the relation continuous.
    if (mass_ratio < 55.0) {
        return luminosity_solar * 1.4 * std::pow(mass_ratio, 3.5);
    }
    return luminosity_solar * 3.2e4 * mass_ratio;
}

// Power law: L = Lsun * coeff * (M / Msun)^exponent.
inline double lum_from_power_law(double mass_kg, double coeff, double exponent) noexcept {
    const double mass_solar       = TidalPyConstants::d_MASS_SOLAR;
    const double luminosity_solar = TidalPyConstants::d_LUMINOSITY_SOLAR;
    if (mass_kg <= 0.0 || mass_solar <= 0.0) {
        return TidalPyConstants::d_NAN;
    }
    return luminosity_solar * coeff * std::pow(mass_kg / mass_solar, exponent);
}

// -------------------------------------------------------------------------------
// Lower-case a model name for case-insensitive factory lookup.
// -------------------------------------------------------------------------------
inline std::string lum_to_lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
    return text;
}

// =====================================================================================================================
// Luminosity models
// =====================================================================================================================

// -------------------------------------------------------------------------------
// c_FixedLuminosity - luminosity supplied directly (mass independent, alias "constant").
// -------------------------------------------------------------------------------
class c_FixedLuminosity : public c_LuminosityBase {
public:
    c_FixedLuminosity() : c_LuminosityBase("fixed") {}
    explicit c_FixedLuminosity(const c_LuminosityConfig& config)
        : c_LuminosityBase("fixed"),
          p_luminosity_w(config.luminosity_w) {}
    ~c_FixedLuminosity() override = default;

    double get_luminosity() const noexcept { return this->p_luminosity_w; }

    double calc_luminosity(double mass_kg) const override {
        return lum_from_fixed(mass_kg, this->p_luminosity_w);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::FixedLuminosity), {this->p_luminosity_w});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 1);
        this->p_luminosity_w = params[0];
    }

protected:
    double p_luminosity_w = 0.0;
};

// -------------------------------------------------------------------------------
// c_MassToLuminosity - piecewise main-sequence L(M) (aliases "cuntz_wang", "cw").
// -------------------------------------------------------------------------------
class c_MassToLuminosity : public c_LuminosityBase {
public:
    c_MassToLuminosity() : c_LuminosityBase("mass_to_luminosity") {}
    explicit c_MassToLuminosity(const c_LuminosityConfig& /*config*/)
        : c_LuminosityBase("mass_to_luminosity") {}
    ~c_MassToLuminosity() override = default;

    double calc_luminosity(double mass_kg) const override {
        return lum_from_mass(mass_kg);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::MassToLuminosity));
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }
};

// -------------------------------------------------------------------------------
// c_PowerLawLuminosity - single power law L = Lsun * coeff * (M/Msun)^p (alias "power_law").
// -------------------------------------------------------------------------------
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

    double calc_luminosity(double mass_kg) const override {
        return lum_from_power_law(mass_kg, this->p_coeff, this->p_exponent);
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

// =====================================================================================================================
// Factory
// =====================================================================================================================

// -------------------------------------------------------------------------------
// c_LuminosityModel - one named value per luminosity model.
// -------------------------------------------------------------------------------
enum class c_LuminosityModel : uint8_t {
    Fixed            = 0,
    MassToLuminosity = 1,
    PowerLaw         = 2,
};

// -------------------------------------------------------------------------------
// c_luminosity_model_from_name - map a (case-insensitive) model name or alias to a
// c_LuminosityModel enum value.
//
// Recognized names and aliases:
//   "fixed" / "constant"
//   "mass_to_luminosity" / "cuntz_wang" / "cw"
//   "power_law" / "powerlaw"
//
// Throws std::invalid_argument if the model name is unknown.
// -------------------------------------------------------------------------------
inline c_LuminosityModel c_luminosity_model_from_name(const std::string& model_name) {
    const std::string name = lum_to_lower(model_name);

    if (name == "fixed" || name == "constant") { return c_LuminosityModel::Fixed; }
    if (name == "mass_to_luminosity" || name == "cuntz_wang" || name == "cw") {
        return c_LuminosityModel::MassToLuminosity;
    }
    if (name == "power_law" || name == "powerlaw") { return c_LuminosityModel::PowerLaw; }

    throw std::invalid_argument("TidalPy: unknown luminosity model name '" + model_name + "'");
}

// -------------------------------------------------------------------------------
// c_find_luminosity - build the luminosity model named by a c_LuminosityModel.
//
// Returns a unique_ptr to a newly heap-allocated concrete model constructed from the supplied config.
// This is the canonical C++ factory. Throws std::invalid_argument for an unrecognised enum value.
// -------------------------------------------------------------------------------
inline std::unique_ptr<c_LuminosityBase> c_find_luminosity(
        c_LuminosityModel model, const c_LuminosityConfig& config) {
    switch (model) {
        case c_LuminosityModel::Fixed:            return std::make_unique<c_FixedLuminosity>(config);
        case c_LuminosityModel::MassToLuminosity: return std::make_unique<c_MassToLuminosity>(config);
        case c_LuminosityModel::PowerLaw:         return std::make_unique<c_PowerLawLuminosity>(config);
    }
    throw std::invalid_argument("TidalPy: unrecognised c_LuminosityModel enum value");
}

// -------------------------------------------------------------------------------
// c_find_luminosity (name overload) - maps a name/alias to the enum and builds the model.
// -------------------------------------------------------------------------------
inline std::unique_ptr<c_LuminosityBase> c_find_luminosity(
        const std::string& model_name, const c_LuminosityConfig& config) {
    return c_find_luminosity(c_luminosity_model_from_name(model_name), config);
}

// -------------------------------------------------------------------------------
// c_luminosity_from_binary - reconstruct a luminosity model from a binary stream.
//
// Peeks the upcoming record's BinaryClassID (without consuming the header), constructs the matching
// default-initialized concrete model, then delegates to its read_binary to restore the model name and
// parameters. Throws std::runtime_error if the class id is not a known luminosity model.
// -------------------------------------------------------------------------------
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
