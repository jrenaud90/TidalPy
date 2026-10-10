#pragma once
/* Bulk mixing: how melt changes a partially molten aggregate's bulk modulus [Pa] and bulk viscosity [Pa s]. Each is
 * an optional law a material holds; without one the solid phase's value stands below the liquidus. All MKS.
 *
 * References
 * ----------
 * - Hashin and Shtrikman (1963), J. Mech. Phys. Solids 11, 127: bounds on the moduli of a two-phase aggregate.
 * - Mavko (1980), JGR 85, 5173; Takei (2002), JGR 107, 2043: melt bears on the bulk modulus far less than on the
 *   shear modulus.
 * - McKenzie (1984), J. Petrol. 25, 713; Takei and Holtzman (2009), JGR 114, B06205: the compaction (bulk) viscosity
 *   of a partially molten rock, eta / phi and of order eta respectively.
 * - Kervazo et al. (2021), A&A 650, A72: bulk dissipation in Io's partially molten interior.
 */

#include <algorithm>
#include <cmath>
#include <istream>
#include <memory>
#include <string>
#include <vector>

#include "physics_base_.hpp"
#include "registry_.hpp"
#include "spec_model_.hpp"
#include "../constants_.hpp"

namespace tidalpy {

class c_BulkModulusMixingBase : public c_PhysicsBase {
public:
    static constexpr const char* C_FAMILY_NAME = "bulk-modulus mixing";

    explicit c_BulkModulusMixingBase(const std::string& model_name) : c_PhysicsBase(model_name) {}
    ~c_BulkModulusMixingBase() override = default;

    // The aggregate's bulk modulus [Pa] from the solid's and the liquid's at the local pressure, the framework's
    // (post-melt) shear modulus, and the melt fraction. The solid's value with no melt; NaN for a non-finite one.
    virtual double calc_bulk_modulus(
        double solid_bulk_modulus,
        double liquid_bulk_modulus,
        double framework_shear_modulus,
        double melt_fraction) const noexcept = 0;
};

// The Hashin-Shtrikman (1963) bound for melt of bulk modulus K_l in a solid framework of bulk modulus K_s and shear
// modulus mu (aliases "hs", "hashin-shtrikman"):
//     K = K_s + phi / (1 / (K_l - K_s) + (1 - phi) / (K_s + (4/3) mu)).
// While the framework holds this is the upper bound for isolated melt pockets, a weak reduction (about 16 percent at
// 10 percent melt for K_s 130, mu 60, K_l 20 GPa); once the framework's shear modulus has collapsed it becomes the
// Reuss (Wood 1955) average of a suspension, and it reaches K_l at phi = 1. This is the unrelaxed (undrained) modulus:
// its relaxation as melt moves is a bulk rheology's job, at the rate the bulk viscosity sets.
class c_HashinShtrikmanMixing final : public c_SpecModel<c_HashinShtrikmanMixing, c_BulkModulusMixingBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::HashinShtrikmanMixing;
    static const std::vector<c_ParamSpec<c_HashinShtrikmanMixing>>& parameter_specs() {
        static const std::vector<c_ParamSpec<c_HashinShtrikmanMixing>> specs;
        return specs;
    }

    c_HashinShtrikmanMixing() : c_HashinShtrikmanMixing(c_ParamMap{}) {}
    explicit c_HashinShtrikmanMixing(const c_ParamMap& params) : c_SpecModel("hashin_shtrikman") {
        this->p_initialize(params);
    }

    double calc_bulk_modulus(
            double solid_bulk_modulus,
            double liquid_bulk_modulus,
            double framework_shear_modulus,
            double melt_fraction) const noexcept override {
        if (!std::isfinite(melt_fraction)) { return TidalPyConstants::d_NAN; }
        if (melt_fraction <= 0.0) { return solid_bulk_modulus; }
        const double contrast = liquid_bulk_modulus - solid_bulk_modulus;
        if (std::abs(contrast) <= TidalPyConstants::d_EPS * std::abs(solid_bulk_modulus)) { return solid_bulk_modulus; }
        const double framework_term =
            (1.0 - melt_fraction) / (solid_bulk_modulus + (4.0 / 3.0) * std::max(framework_shear_modulus, 0.0));
        return solid_bulk_modulus + melt_fraction / ((1.0 / contrast) + framework_term);
    }
};

class c_BulkViscosityMixingBase : public c_PhysicsBase {
public:
    static constexpr const char* C_FAMILY_NAME = "bulk-viscosity mixing";

    explicit c_BulkViscosityMixingBase(const std::string& model_name) : c_PhysicsBase(model_name) {}
    ~c_BulkViscosityMixingBase() override = default;

    // The aggregate's bulk viscosity [Pa s] from the solid's, the post-melt shear viscosity, and the melt fraction.
    // The solid's value with no melt; NaN for a non-finite melt fraction.
    virtual double calc_bulk_viscosity(
        double solid_bulk_viscosity,
        double postmelt_shear_viscosity,
        double melt_fraction) const noexcept = 0;
};

// Compaction viscosity (alias "mckenzie"): melt adds a matrix bulk viscosity zeta_melt = c eta / phi^n, with eta the
// post-melt shear viscosity, in series with the solid's: 1 / zeta = 1 / zeta_solid + phi^n / (c eta). n = 1 with c of
// order 1 is McKenzie (1984); n = 0 (a bulk viscosity of order the shear one) is closer to Takei and Holtzman (2009).
// A non-finite or non-positive solid bulk viscosity counts as no solid dashpot, so melt alone sets zeta.
class c_CompactionViscosity final : public c_SpecModel<c_CompactionViscosity, c_BulkViscosityMixingBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::CompactionViscosity;

    static const std::vector<c_ParamSpec<c_CompactionViscosity>>& parameter_specs() {
        using Self = c_CompactionViscosity;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"coefficient", "coefficient", &Self::p_coefficient, 1.0, c_ParamBounds::Positive,
             "c in zeta_melt = c eta / phi^n [dimensionless]."},
            {"exponent", "exponent", &Self::p_exponent, 1.0, c_ParamBounds::NonNegative,
             "n in zeta_melt = c eta / phi^n [dimensionless]."},
        };
        return specs;
    }

    c_CompactionViscosity() : c_CompactionViscosity(c_ParamMap{}) {}
    explicit c_CompactionViscosity(const c_ParamMap& params) : c_SpecModel("compaction") {
        this->p_initialize(params);
    }

    double calc_bulk_viscosity(
            double solid_bulk_viscosity,
            double postmelt_shear_viscosity,
            double melt_fraction) const noexcept override {
        if (!std::isfinite(melt_fraction)) { return TidalPyConstants::d_NAN; }
        if (melt_fraction <= 0.0) { return solid_bulk_viscosity; }
        const double melt_term = std::pow(melt_fraction, this->p_exponent)
            / (this->p_coefficient * postmelt_shear_viscosity);
        const double solid_term = (std::isfinite(solid_bulk_viscosity) && (solid_bulk_viscosity > 0.0))
            ? 1.0 / solid_bulk_viscosity : 0.0;
        const double inverse = solid_term + melt_term;
        if (!(inverse > 0.0)) { return solid_bulk_viscosity; }
        return 1.0 / inverse;
    }

protected:
    double p_coefficient = 0.0;
    double p_exponent    = 0.0;
};

inline const c_ModelRegistry<c_BulkModulusMixingBase>& c_bulk_modulus_mixing_registry() {
    static const c_ModelRegistry<c_BulkModulusMixingBase> registry = {
        {{"hashin_shtrikman", "hs", "hashin-shtrikman"}, BinaryClassID::HashinShtrikmanMixing,
         &c_make_entry<c_BulkModulusMixingBase, c_HashinShtrikmanMixing>},
    };
    return registry;
}

inline const c_ModelRegistry<c_BulkViscosityMixingBase>& c_bulk_viscosity_mixing_registry() {
    static const c_ModelRegistry<c_BulkViscosityMixingBase> registry = {
        {{"compaction", "mckenzie"}, BinaryClassID::CompactionViscosity,
         &c_make_entry<c_BulkViscosityMixingBase, c_CompactionViscosity>},
    };
    return registry;
}

inline std::unique_ptr<c_BulkModulusMixingBase> c_find_bulk_modulus_mixing(
        const std::string& model_name, const c_ParamMap& params) {
    return c_make_model(c_bulk_modulus_mixing_registry(), model_name, params);
}

inline std::unique_ptr<c_BulkModulusMixingBase> c_bulk_modulus_mixing_from_binary(std::istream& in, bool force = false) {
    return c_model_from_binary(c_bulk_modulus_mixing_registry(), in, force);
}

inline std::string c_bulk_modulus_mixing_canonical_name(const std::string& model_name) {
    return c_canonical_model_name(c_bulk_modulus_mixing_registry(), model_name);
}

inline std::unique_ptr<c_BulkViscosityMixingBase> c_find_bulk_viscosity_mixing(
        const std::string& model_name, const c_ParamMap& params) {
    return c_make_model(c_bulk_viscosity_mixing_registry(), model_name, params);
}

inline std::unique_ptr<c_BulkViscosityMixingBase> c_bulk_viscosity_mixing_from_binary(
        std::istream& in, bool force = false) {
    return c_model_from_binary(c_bulk_viscosity_mixing_registry(), in, force);
}

inline std::string c_bulk_viscosity_mixing_canonical_name(const std::string& model_name) {
    return c_canonical_model_name(c_bulk_viscosity_mixing_registry(), model_name);
}

}  // namespace tidalpy
