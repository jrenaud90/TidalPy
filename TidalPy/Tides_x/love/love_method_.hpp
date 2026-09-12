#pragma once
/*
 * love_method_.hpp - Love-number solution methods and the homogeneous-sphere formulas.
 *
 * c_LoveMethod names how a world's Love numbers are obtained:
 *   RadialSolver           ("radial_solver", "shooting", "rs")                 - shooting-method radial solve (default)
 *   PropagationMatrix      ("propagation_matrix", "prop_matrix", "pm", "prop") - propagation-matrix radial solve
 *   Homogeneous            ("homogeneous", "homogen")                          - analytic homogeneous-sphere formulas with
 *                                                                                the volume-averaged complex shear modulus
 *   HomogeneousCPL         ("cpl")                                             - homogeneous formulas on the static modulus,
 *                                                                                then k_l (1 - i / Q_l)
 *   HomogeneousCTL         ("ctl")                                             - homogeneous formulas on the static modulus,
 *                                                                                then k_l (1 - i omega dt_l)
 *   LaterallyInhomogeneous ("laterally_inhomogeneous", "3d", "lat_inhom")      - reserved for the 3D Love solver
 *
 * Homogeneous incompressible sphere (Love 1911; Munk & MacDonald 1960):
 *   mu_eff_l = (2 l^2 + 4 l + 3) / l * mu / (rho g R)
 *   k_l = 3 / (2 (l - 1)) / (1 + mu_eff_l)
 *   h_l = (2 l + 1) / (2 (l - 1)) / (1 + mu_eff_l)
 *   l_l = 3 / (2 l (l - 1)) / (1 + mu_eff_l)
 * With a complex (viscoelastic) shear modulus these give complex Love numbers; with the real static
 * modulus they give the static (elastic) Love numbers that the CPL/CTL structures multiply.
 */

#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>

#include "love_.hpp"

namespace tidalpy {

enum class c_LoveMethod : int {
    RadialSolver           = 0,
    PropagationMatrix      = 1,
    Homogeneous            = 2,
    HomogeneousCPL         = 3,
    HomogeneousCTL         = 4,
    LaterallyInhomogeneous = 5
};

// Canonical name of a method (the first alias listed above).
inline const char* c_love_method_name(c_LoveMethod method) noexcept
{
    switch (method) {
        case c_LoveMethod::RadialSolver:           return "radial_solver";
        case c_LoveMethod::PropagationMatrix:      return "propagation_matrix";
        case c_LoveMethod::Homogeneous:            return "homogeneous";
        case c_LoveMethod::HomogeneousCPL:         return "cpl";
        case c_LoveMethod::HomogeneousCTL:         return "ctl";
        case c_LoveMethod::LaterallyInhomogeneous: return "laterally_inhomogeneous";
    }
    return "unknown";
}

// Parse a method name or alias (case-insensitive). Throws std::invalid_argument for unknown names.
inline c_LoveMethod c_parse_love_method(const std::string& name)
{
    std::string key = name;
    std::transform(key.begin(), key.end(), key.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (key == "radial_solver" || key == "shooting" || key == "rs") {
        return c_LoveMethod::RadialSolver;
    }
    if (key == "propagation_matrix" || key == "prop_matrix" || key == "pm" || key == "prop") {
        return c_LoveMethod::PropagationMatrix;
    }
    if (key == "homogeneous" || key == "homogen") {
        return c_LoveMethod::Homogeneous;
    }
    if (key == "cpl") {
        return c_LoveMethod::HomogeneousCPL;
    }
    if (key == "ctl") {
        return c_LoveMethod::HomogeneousCTL;
    }
    if (key == "laterally_inhomogeneous" || key == "3d" || key == "lat_inhom") {
        return c_LoveMethod::LaterallyInhomogeneous;
    }
    throw std::invalid_argument(
        "TidalPy: unknown Love-number method '" + name + "'. Choose from radial_solver (shooting, rs), "
        "propagation_matrix (prop_matrix, pm, prop), homogeneous (homogen), cpl, ctl, "
        "laterally_inhomogeneous (3d, lat_inhom).");
}

// Convert a stored integer (e.g. from a config struct) back to the enum; out-of-range values throw.
inline c_LoveMethod c_love_method_from_int(int value)
{
    if (value < 0 || value > static_cast<int>(c_LoveMethod::LaterallyInhomogeneous)) {
        throw std::invalid_argument("TidalPy: Love-number method index " + std::to_string(value) + " is out of range.");
    }
    return static_cast<c_LoveMethod>(value);
}

// Integer-typed conveniences for the Cython layer (the enum is stored as an int in the config structs).
inline int c_parse_love_method_int(const std::string& name)
{
    return static_cast<int>(c_parse_love_method(name));
}

inline std::string c_love_method_name_int(int value)
{
    return std::string(c_love_method_name(c_love_method_from_int(value)));
}

// True for the methods that run the radial solver (and therefore provide radial y-functions).
inline bool c_love_method_uses_radial_solver(c_LoveMethod method) noexcept
{
    return method == c_LoveMethod::RadialSolver || method == c_LoveMethod::PropagationMatrix;
}

// True for the analytic homogeneous-sphere methods (homogeneous, cpl, ctl).
inline bool c_love_method_is_homogeneous(c_LoveMethod method) noexcept
{
    return method == c_LoveMethod::Homogeneous || method == c_LoveMethod::HomogeneousCPL
        || method == c_LoveMethod::HomogeneousCTL;
}

// ---------------------------------------------------------------------------------------------------------------------
// Homogeneous incompressible sphere
// ---------------------------------------------------------------------------------------------------------------------

// Degree-l effective rigidity: (2 l^2 + 4 l + 3) / l * mu / (rho g R). Templated so the shear modulus may be real
// (static rigidity) or complex (viscoelastic).
template <typename ModulusT>
inline ModulusT c_calc_effective_rigidity(
        ModulusT shear_modulus,
        double density,
        double gravity,
        double radius,
        int degree_l)
{
    if (degree_l < 2) {
        throw std::invalid_argument("TidalPy: the homogeneous Love-number formulas need degree_l >= 2.");
    }
    const double l = static_cast<double>(degree_l);
    const double structure = density * gravity * radius;
    if (!(structure > 0.0) || !std::isfinite(structure)) {
        throw std::invalid_argument(
            "TidalPy: the homogeneous Love-number formulas need positive, finite density, gravity, and radius.");
    }
    return ((2.0 * l * l + 4.0 * l + 3.0) / l) * shear_modulus / structure;
}

// Non-template spellings of the effective rigidity for the Cython layer.
inline double c_calc_effective_rigidity_real(
        double shear_modulus, double density, double gravity, double radius, int degree_l)
{
    return c_calc_effective_rigidity(shear_modulus, density, gravity, radius, degree_l);
}

inline std::complex<double> c_calc_effective_rigidity_complex(
        std::complex<double> shear_modulus, double density, double gravity, double radius, int degree_l)
{
    return c_calc_effective_rigidity(shear_modulus, density, gravity, radius, degree_l);
}

// Love numbers k_l, h_l, l_l of a homogeneous incompressible sphere from its (complex) shear modulus.
inline c_LoveNumbers c_calc_homogeneous_love_numbers(
        std::complex<double> complex_shear_modulus,
        double density,
        double gravity,
        double radius,
        int degree_l)
{
    const std::complex<double> mu_eff = c_calc_effective_rigidity(
        complex_shear_modulus, density, gravity, radius, degree_l);
    const double l = static_cast<double>(degree_l);
    const std::complex<double> response = 1.0 / (std::complex<double>(1.0, 0.0) + mu_eff);
    c_LoveNumbers love;
    love.k = (3.0 / (2.0 * (l - 1.0))) * response;
    love.h = ((2.0 * l + 1.0) / (2.0 * (l - 1.0))) * response;
    love.l = (3.0 / (2.0 * l * (l - 1.0))) * response;
    return love;
}

// Constant phase lag: every Love number is multiplied by (1 - i / Q), so -Im[k] = Re[k] / Q.
inline c_LoveNumbers c_apply_fixed_q(const c_LoveNumbers& love, double fixed_q)
{
    if (!(fixed_q > 0.0) || !std::isfinite(fixed_q)) {
        throw std::invalid_argument("TidalPy: the cpl Love-number method needs a positive, finite fixed_q.");
    }
    const std::complex<double> lag(1.0, -1.0 / fixed_q);
    return c_LoveNumbers(love.k * lag, love.h * lag, love.l * lag);
}

// Constant time lag: every Love number is multiplied by (1 - i omega dt), so -Im[k] = Re[k] omega dt.
inline c_LoveNumbers c_apply_fixed_dt(const c_LoveNumbers& love, double frequency, double fixed_dt)
{
    if (!(fixed_dt >= 0.0) || !std::isfinite(fixed_dt)) {
        throw std::invalid_argument("TidalPy: the ctl Love-number method needs a non-negative, finite fixed_dt.");
    }
    const std::complex<double> lag(1.0, -std::abs(frequency) * fixed_dt);
    return c_LoveNumbers(love.k * lag, love.h * lag, love.l * lag);
}

} // namespace tidalpy
