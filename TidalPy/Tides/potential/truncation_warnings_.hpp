#pragma once
/* The warnings a global or 3D tidal solve gives when its truncations misstate the tides, shared by a world's solves and
 * the standalone tide functions (collapse_global_tides, global_potential, tidal_potential_3d_modes).
 *
 * Three cases: an obliquity that the 'off' obliquity truncation ignores (shown once per owner), an obliquity past the
 * range of its truncation (c_obliquity_truncation_limit), and an eccentricity past the range of its truncation
 * (c_eccentricity_truncation_limit), where a low level can misstate the heating by 10% or more, or even make it
 * negative near e = 1. The two range warnings are shown once per owner and truncation level, so an owner that moves to
 * another level is warned again when that level is out of its range.
 */

#include <cmath>
#include <cstdio>
#include <limits>
#include <mutex>
#include <set>
#include <string>

#include "../eccentricity/eccentricity_accuracy_.hpp"   // c_eccentricity_truncation_limit
#include "../obliquity/obliquity_accuracy_.hpp"         // c_obliquity_truncation_limit, C_OBLIQUITY_OFF
#include "../../Utilities/logging/logger_.hpp"

// The highest tabulated eccentricity truncation, whose range the warning quotes as the place to switch to 'exact'.
inline constexpr int C_HIGHEST_ECCENTRICITY_TRUNCATION = 50;

// The fewest significant digits a range warning prints its value and limit with.
inline constexpr int C_RANGE_WARNING_MIN_DIGITS = 4;

// Which warnings an owner has already shown: the ignored obliquity once, and each range warning once per level.
struct c_TruncationWarningsShown {
    bool obliquity_off = false;
    std::set<int> obliquity_range_levels;
    std::set<int> eccentricity_range_levels;
};

// The fewest significant digits, at least C_RANGE_WARNING_MIN_DIGITS, at which `value` and `limit` print differently
// in the "g" format, so a warning never quotes a value past its limit as equal to it.
inline int c_distinct_digits(double value, double limit) noexcept {
    constexpr int max_digits = std::numeric_limits<double>::max_digits10;
    // Room for a sign, max_digits digits, a decimal point, and an exponent.
    char value_text[max_digits + 16];
    char limit_text[max_digits + 16];
    for (int digits = C_RANGE_WARNING_MIN_DIGITS; digits < max_digits; ++digits) {
        std::snprintf(value_text, sizeof(value_text), "%.*g", digits, value);
        std::snprintf(limit_text, sizeof(limit_text), "%.*g", digits, limit);
        if (std::string(value_text) != std::string(limit_text)) {
            return digits;
        }
    }
    return max_digits;
}

// Warn when the truncations misstate the tides at this orbital state: an ignored obliquity once per owner, and each
// range case once per owner and truncation level. `owner` names who is warned about ("world 'io'", or the standalone
// function's name) and `remedy` says where the truncation levels are set ("its [tides] table or set_tide_config", or
// "the eccentricity_truncation and obliquity_truncation arguments").
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
        // insert reports whether the level is new, so the warning is shown on its first occurrence only.
        const bool past_limit = std::abs(obliquity) > obliquity_limit;
        if (past_limit && shown.obliquity_range_levels.insert(obliquity_truncation).second) {
            const int digits = c_distinct_digits(std::abs(obliquity), obliquity_limit);
            TIDALPY_LOG_WARN(
                "TidalPy: {} has an obliquity of {:.{}g} rad, past {:.{}g}, where its obliquity truncation (level {}) "
                "can misstate the tides by 10% or more. Raise the obliquity truncation in {} "
                "(recommend_obliquity_truncation picks a level for a tolerance), or use 'gen'. Shown once per level.",
                owner, obliquity, digits, obliquity_limit, digits, obliquity_truncation, remedy);
        }
    }
    const double eccentricity_limit = c_eccentricity_truncation_limit(eccentricity_truncation, max_degree_l);
    if ((eccentricity > eccentricity_limit) && shown.eccentricity_range_levels.insert(eccentricity_truncation).second) {
        const int digits = c_distinct_digits(eccentricity, eccentricity_limit);
        TIDALPY_LOG_WARN(
            "TidalPy: {} has an eccentricity of {:.{}g}, past {:.{}g}, where its eccentricity truncation (level {}) "
            "can misstate the tides by 10% or more (far past it, the heating can come out negative). Raise the "
            "eccentricity truncation in {} (recommend_eccentricity_truncation picks a level for a tolerance); past "
            "about {:.2f} use 'exact'. Shown once per level.",
            owner, eccentricity, digits, eccentricity_limit, digits, eccentricity_truncation, remedy,
            c_eccentricity_truncation_limit(C_HIGHEST_ECCENTRICITY_TRUNCATION, max_degree_l));
    }
}

// The same for the standalone tide functions, which have no world to remember the warnings: each case is shown once
// (each range case once per level) per process by each extension module that calls this (the flags are this
// header's, compiled into each module).
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
