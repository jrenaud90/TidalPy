#pragma once
/* The warnings a global or 3D tidal solve gives when its truncations misstate the tides, shared by a world's solves and
 * the standalone tide functions (collapse_global_tides, global_potential, tidal_potential_3d_modes).
 *
 * Three cases, each shown once per owner: an obliquity that the 'off' obliquity truncation ignores, an obliquity past
 * the range of its truncation (c_obliquity_truncation_limit), and an eccentricity past the range of its truncation
 * (c_eccentricity_truncation_limit), where a low level can underestimate the heating by 10% or more, or even make it
 * negative near e = 1.
 */

#include <cmath>
#include <mutex>
#include <string>

#include "../eccentricity/eccentricity_accuracy_.hpp"   // c_eccentricity_truncation_limit
#include "../obliquity/obliquity_accuracy_.hpp"         // c_obliquity_truncation_limit, C_OBLIQUITY_OFF
#include "../../Utilities/logging/logger_.hpp"

// The highest tabulated eccentricity truncation, whose range the warning quotes as the place to switch to 'exact'.
inline constexpr int C_HIGHEST_ECCENTRICITY_TRUNCATION = 50;

// Which of the three warnings an owner has already shown.
struct c_TruncationWarningsShown {
    bool obliquity_off        = false;
    bool obliquity_range      = false;
    bool eccentricity_range   = false;
};

// Warn, once per owner and case, when the truncations misstate the tides at this orbital state. `owner` names who is
// warned about ("world 'io'", or the standalone function's name) and `remedy` says where the truncation levels are set
// ("its [tides] table or set_tide_config", or "the eccentricity_truncation and obliquity_truncation arguments").
inline void c_warn_tide_truncations(
        const std::string& owner,
        const std::string& remedy,
        double eccentricity,
        double obliquity,
        int eccentricity_truncation,
        int obliquity_truncation,
        int max_degree_l,
        c_TruncationWarningsShown& shown) {
    if ((obliquity != 0.0) && (obliquity_truncation == C_OBLIQUITY_OFF) && !shown.obliquity_off) {
        shown.obliquity_off = true;
        TIDALPY_LOG_WARN(
            "TidalPy: {} has an obliquity of {:.3e} rad but its obliquity truncation is off, so its obliquity tides "
            "are ignored. Set the obliquity truncation (2, 4, or 'gen') in {} to include them. Shown once.",
            owner, obliquity, remedy);
    }
    if (obliquity_truncation != C_OBLIQUITY_OFF) {
        const double obliquity_limit = c_obliquity_truncation_limit(obliquity_truncation, max_degree_l);
        if ((std::abs(obliquity) > obliquity_limit) && !shown.obliquity_range) {
            shown.obliquity_range = true;
            TIDALPY_LOG_WARN(
                "TidalPy: {} has an obliquity of {:.3f} rad, past {:.3f}, where its obliquity truncation (level {}) "
                "can misstate the tides by 10% or more. Raise the obliquity truncation in {} "
                "(recommend_obliquity_truncation picks a level for a tolerance), or use 'gen'. Shown once.",
                owner, obliquity, obliquity_limit, obliquity_truncation, remedy);
        }
    }
    const double eccentricity_limit = c_eccentricity_truncation_limit(eccentricity_truncation, max_degree_l);
    if ((eccentricity > eccentricity_limit) && !shown.eccentricity_range) {
        shown.eccentricity_range = true;
        TIDALPY_LOG_WARN(
            "TidalPy: {} has an eccentricity of {:.3f}, past {:.3f}, where its eccentricity truncation (level {}) can "
            "underestimate the tides by 10% or more (far past it, the heating can come out negative). Raise the "
            "eccentricity truncation in {} (recommend_eccentricity_truncation picks a level for a tolerance); past "
            "about {:.2f} use 'exact'. Shown once.",
            owner, eccentricity, eccentricity_limit, eccentricity_truncation, remedy,
            c_eccentricity_truncation_limit(C_HIGHEST_ECCENTRICITY_TRUNCATION, max_degree_l));
    }
}

// The same for the standalone tide functions, which have no world to remember the warnings: each case is shown once
// per process by each extension module that calls this (the flags are this header's, compiled into each module).
inline void c_warn_standalone_tide_truncations(
        const std::string& function_name,
        double eccentricity,
        double obliquity,
        int eccentricity_truncation,
        int obliquity_truncation,
        int max_degree_l) {
    static std::mutex shown_mutex;
    static c_TruncationWarningsShown shown;
    const std::lock_guard<std::mutex> guard(shown_mutex);
    c_warn_tide_truncations(
        function_name, "the eccentricity_truncation and obliquity_truncation arguments", eccentricity, obliquity,
        eccentricity_truncation, obliquity_truncation, max_degree_l, shown);
}
