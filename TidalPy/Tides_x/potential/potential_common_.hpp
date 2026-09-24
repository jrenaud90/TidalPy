#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "intmap_.hpp"
#include "keys_.hpp"
#include "tolerance_index_.hpp"
#include "../../Utilities_x/math_x/numerics_.hpp"
#include "constants_.hpp"

struct c_FrequencyStorage
{
    size_t num_instances;
    double frequency;
    c_FrequencyStorage() : 
        num_instances(0),
        frequency(TidalPyConstants::d_NAN)
    {
    }
    c_FrequencyStorage(double freq) :
        num_instances(1),
        frequency(freq)
    {
    }
};

struct c_ModeStorage
{
    double mode;
    double mode_strength;
    int n_coeff;
    int o_coeff;
    c_ModeStorage() : 
        mode(TidalPyConstants::d_NAN),
        mode_strength(TidalPyConstants::d_NAN),
        n_coeff(0),
        o_coeff(0)
    {
    }

    c_ModeStorage(double mode_, double mode_strength_, int n_coeff_, int o_coeff_) :
        mode(mode_),
        mode_strength(mode_strength_),
        n_coeff(n_coeff_),
        o_coeff(o_coeff_)
    {
    }
};

typedef c_IntMap<c_Key4, c_ModeStorage> c_ModeMap;
typedef c_IntMap<c_Key4, size_t> c_UniqueFreqIndexMap;
typedef std::vector<c_FrequencyStorage> c_UniqueFreqMap;


// Modes are collapsed onto a shared frequency when they agree to this much; it also decides when a frequency
// counts as zero. The fallback covers a call made before the config is loaded.
inline constexpr double d_FREQUENCY_MATCH_RTOL_FALLBACK = 1.0e-9;

inline double c_frequency_match_rtol() noexcept
{
    if (tidalpy_config_ptr != nullptr && std::isfinite(tidalpy_config_ptr->d_FREQUENCY_MATCH_RTOL))
    {
        return tidalpy_config_ptr->d_FREQUENCY_MATCH_RTOL;
    }
    return d_FREQUENCY_MATCH_RTOL_FALLBACK;
}

// The frequency tolerances of one engine call, read from the config once: modes whose frequencies agree to match_rtol
// share a frequency, and min_frequency is the floor at which a mode counts as static (zero before the config is
// loaded).
struct c_FrequencyTolerance
{
    double match_rtol;
    double min_frequency;
};

inline c_FrequencyTolerance c_read_frequency_tolerance() noexcept
{
    return c_FrequencyTolerance{
        c_frequency_match_rtol(),
        (tidalpy_config_ptr != nullptr) ? tidalpy_config_ptr->d_MIN_FREQUENCY : 0.0};
}

// Whether two frequencies agree to a relative tolerance, NaN-safe (c_isclose with no absolute floor).
struct c_IsCloseMatch
{
    double rtol;
    bool operator()(double stored, double query) const noexcept
    {
        return c_isclose(query, stored, this->rtol, 0.0);
    }
};

// The distinct nonzero mode frequencies of one engine call, indexed for c_record_unique_frequencies.
struct c_UniqueFrequencyTracker
{
    c_FrequencyTolerance tolerance;
    c_ToleranceIndex<c_IsCloseMatch> index;

    explicit c_UniqueFrequencyTracker(const c_FrequencyTolerance& tolerance_in) :
        tolerance(tolerance_in),
        index(tolerance_in.match_rtol, c_IsCloseMatch{tolerance_in.match_rtol})
    {
    }
};

// Records a mode's frequency against the first recorded frequency it matches, or as a new one, and returns that
// frequency's index. A zero-frequency mode is a static deformation: it dissipates nothing, so it is not recorded and
// the result is -1.
inline std::ptrdiff_t c_record_unique_frequencies(
        const c_Key4& lmpq_key,
        double frequency,
        c_UniqueFrequencyTracker& tracker,
        c_UniqueFreqIndexMap& frequency_index_map,
        c_UniqueFreqMap& frequency_map)
{
    const c_FrequencyTolerance& tolerance = tracker.tolerance;
    if (c_isclose(frequency, 0.0, tolerance.match_rtol, tolerance.min_frequency))
    {
        return -1;
    }

    std::ptrdiff_t frequency_index = tracker.index.find(0, frequency);
    if (frequency_index >= 0)
    {
        frequency_map[static_cast<size_t>(frequency_index)].num_instances += 1;
    }
    else
    {
        frequency_index = static_cast<std::ptrdiff_t>(tracker.index.insert(0, frequency));
        frequency_map.emplace_back(frequency);
    }
    frequency_index_map.set(lmpq_key, static_cast<size_t>(frequency_index));
    return frequency_index;
}

// (l - m)! / (l + m)! (2 - delta_m0) for l = 2 to 10 and 0 <= m <= l, degree l's row starting at entry
// l (l + 1) / 2 - 3. The (2 - delta_m0) of the potential's real form is already in each entry, so callers must not
// apply it again.
inline constexpr double C_LM_COEFFS[63] = {
    // l = 2
    1.000000000000000e+0, 3.333333333333333e-1, 8.333333333333333e-2,
    // l = 3
    1.000000000000000e+0, 1.666666666666667e-1, 1.666666666666667e-2, 2.777777777777778e-3,
    // l = 4
    1.000000000000000e+0, 1.000000000000000e-1, 5.555555555555556e-3, 3.968253968253968e-4, 4.960317460317460e-5,
    // l = 5
    1.000000000000000e+0, 6.666666666666667e-2, 2.380952380952381e-3, 9.920634920634921e-5, 5.511463844797178e-6,
    5.511463844797178e-7,
    // l = 6
    1.000000000000000e+0, 4.761904761904762e-2, 1.190476190476190e-3, 3.306878306878307e-5, 1.102292768959436e-6,
    5.010421677088344e-8, 4.175351397573620e-9,
    // l = 7
    1.000000000000000e+0, 3.571428571428571e-2, 6.613756613756614e-4, 1.322751322751323e-5, 3.006253006253006e-7,
    8.350702795147240e-9, 3.211808767364323e-10, 2.294149119545945e-11,
    // l = 8
    1.000000000000000e+0, 2.777777777777778e-2, 3.968253968253968e-4, 6.012506012506013e-6, 1.002084335417669e-7,
    1.927085260418594e-9, 4.588298239091890e-11, 1.529432746363963e-12, 9.558954664774771e-14,
    // l = 9
    1.000000000000000e+0, 2.222222222222222e-2, 2.525252525252525e-4, 3.006253006253006e-6, 3.854170520837188e-8,
    5.505957886910268e-10, 9.176596478183780e-12, 1.911790932954954e-13, 5.622914508691042e-15, 3.123841393717245e-16,
    // l = 10
    1.000000000000000e+0, 1.818181818181818e-2, 1.683501683501684e-4, 1.618751618751619e-6, 1.651787366073080e-8,
    1.835319295636756e-10, 2.294149119545945e-12, 3.373748705214625e-14, 6.247682787434491e-16, 1.644127049324866e-17,
    8.220635246624330e-19,
};

// The (l, m) coefficient of C_LM_COEFFS; NaN outside l = 2 to 10, 0 <= m <= l.
inline double c_lm_coeff(int degree_l, int order_m) noexcept
{
    if ((degree_l < 2) || (degree_l > 10) || (order_m < 0) || (order_m > degree_l))
    {
        return TidalPyConstants::d_NAN;
    }
    return C_LM_COEFFS[degree_l * (degree_l + 1) / 2 - 3 + order_m];
}
