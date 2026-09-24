#pragma once

#include <cmath>
#include <cstddef>

#include "obliquity_driver_.hpp"
#include "eccentricity_driver_.hpp"
#include "potential_common_.hpp"
#include "constants_.hpp"


struct c_GlobalPotentialResultAtMode
{
    double dU_dM;
    double dU_dw;
    double dU_dO;
    double E_dot;
    // dU_dM - dU_dw of this mode, which is q times the common coefficient. Summed on its own it carries none of the
    // cancellation of the two separate sums, whose q = 0 modes agree and dominate at small eccentricity.
    double dU_dM_minus_dw;
    // The mode's entry in unique_freq_map, the same one unique_freq_index_map holds.
    size_t frequency_index;
    c_GlobalPotentialResultAtMode(
            double dU_dM_,
            double dU_dw_,
            double dU_dO_,
            double E_dot_,
            double dU_dM_minus_dw_,
            size_t frequency_index_) :
        dU_dM(dU_dM_),
        dU_dw(dU_dw_),
        dU_dO(dU_dO_),
        E_dot(E_dot_),
        dU_dM_minus_dw(dU_dM_minus_dw_),
        frequency_index(frequency_index_)
    {
    }
};


// mode_map, unique_freq_index_map, and potential_map hold the same (l, m, p, q) keys in the same order: every mode
// with a nonzero frequency.
struct c_GlobalPotentialStorage
{
    c_ModeMap mode_map;
    c_UniqueFreqIndexMap unique_freq_index_map;
    c_UniqueFreqMap unique_freq_map;
    c_IntMap<c_Key4, c_GlobalPotentialResultAtMode> potential_map;
    int error_code = 0;
    int working_on_l = -1;
};

// The tidal potential's mode decomposition for degrees min_degree_l to max_degree_l: each (l, m, p, q) mode's
// frequency and its contributions to the potential derivatives and the heating, before the Love numbers. The
// result's error_code is 0 on success, -1 for a truncation level the obliquity or eccentricity functions do not
// tabulate at a degree, and -2 for a degree they do not tabulate; working_on_l then names the degree.
inline c_GlobalPotentialStorage c_global_potential(
        double planet_radius,
        double semi_major_axis,
        double orbital_frequency,
        double spin_frequency,
        double obliquity,
        double eccentricity,
        double host_mass,
        double G_to_use,
        int min_degree_l,
        int max_degree_l,
        int obliquity_truncation,
        int eccentricity_truncation,
        double eccentricity_exact_tolerance
    )
{
    c_GlobalPotentialStorage result;
    result.error_code = 0;

    // Only degrees 2 to 10 are tabulated; checked before anything is sized from them.
    if ((min_degree_l < 2) || (max_degree_l > 10) || (min_degree_l > max_degree_l)
        || ((eccentricity_truncation < 0) && (eccentricity_truncation != C_ECCENTRICITY_EXACT)))
    {
        result.error_code   = -2;
        result.working_on_l = (min_degree_l < 2) ? min_degree_l : max_degree_l;
        return result;
    }

    // Upper bound on the number of modes; an overestimate, since some modes are skipped.
    size_t target_size = 0;
    for (int degree_l = min_degree_l; degree_l <= max_degree_l; degree_l++)
    {
        target_size += static_cast<size_t>((degree_l + 1) * (degree_l + 1));
    }
    // Modes with |q| <= N / 2 enter the heating; the exact functions' range depends on e, so this is a first guess.
    const int reserve_q = (eccentricity_truncation == C_ECCENTRICITY_EXACT) ? 10 : (eccentricity_truncation / 2);
    target_size *= static_cast<size_t>(2 * reserve_q + 1);
    result.mode_map.reserve(target_size);
    result.unique_freq_index_map.reserve(target_size);
    result.unique_freq_map.reserve(target_size);
    result.potential_map.reserve(target_size);
    // The four maps take about 96 bytes on the stack, and their heap storage grows with the number of modes
    // reserved above.

    // The config's frequency tolerances, read once for the whole call.
    c_UniqueFrequencyTracker frequency_tracker(c_read_frequency_tolerance());

    // For later calculation of the maximum relative mode.
    double max_mode_strength = 0;

    c_Key4 lmpq_key = c_Key4();

    double R_a = planet_radius / semi_major_axis;
    double R_a_2 = R_a * R_a;
    double ra_l_coeff = 0;
    switch (min_degree_l)
    {
    case 2:
        ra_l_coeff = R_a_2 * R_a_2 * R_a;  // 2(l=2) + 1 == 5
        break;
    case 3:
        ra_l_coeff = R_a_2 * R_a_2 * R_a_2 * R_a;  // 2(l=3) + 1 == 7
        break;
    default:
        ra_l_coeff = std::pow(R_a, 2.0 * static_cast<double>(min_degree_l) + 1);
        break;
    }
    // Fold the outermost coefficient in here so nothing outside this function has to apply it; outside, the
    // focus is only on the Love-number multiplier.
    ra_l_coeff *= G_to_use * host_mass / semi_major_axis;

    for (int degree_l = min_degree_l; degree_l < (max_degree_l + 1); degree_l++)
    {
        // Set the degree l we are working on for error reporting.
        result.working_on_l = degree_l;

        // The heating goes as F^2 too: each square cut at the obliquity truncation's power (plain squares for the
        // general functions).
        const c_ObliquityValues obliquity_squared = c_obliquity_squared_values(
            &result.error_code,
            obliquity,
            degree_l,
            obliquity_truncation
        );
        if (result.error_code != 0)
        {
            return result;
        }

        // The heating goes as G^2: each square cut at the truncation's power, so the mode sum is the heating's
        // Taylor series through e^N.
        const c_EccentricityValues eccentricity_squared = c_eccentricity_squared_values(
            &result.error_code,
            eccentricity,
            degree_l,
            eccentricity_truncation,
            eccentricity_exact_tolerance
        );
        if (result.error_code != 0)
        {
            return result;
        }
        if (degree_l > min_degree_l)
        {
            // Every sequential degree l grows the coeff by (R/a)^2
            ra_l_coeff *= R_a_2;
        }

        lmpq_key.a = degree_l;
        const int max_q = eccentricity_squared.max_q;

        // A zero F^2 or G^2 (a mode the truncation leaves out, or one that vanishes) contributes nothing.
        for (int order_m = 0; order_m <= degree_l; ++order_m)
        {
            // (l - m)! / (l + m)! (2 - d_m0) coefficient
            const double lm_coeff = c_lm_coeff(degree_l, order_m);
            lmpq_key.b = order_m;

            for (int p = 0; p <= degree_l; ++p)
            {
                const double F_lmp_squared = obliquity_squared.value(order_m, p);
                if (F_lmp_squared == 0.0)
                {
                    continue;
                }
                lmpq_key.c = p;

                // The global potential goes as F^2 (the 3D path uses F); fold in (l - m)!/(l + m)!(2 - d_m0).
                double lmp_coeff = F_lmp_squared * ra_l_coeff * lm_coeff;

                for (int q = -max_q; q <= max_q; ++q)
                {
                    const double G_lpq_squared = eccentricity_squared.value(p, q);
                    if (G_lpq_squared == 0.0)
                    {
                        continue;
                    }

                    lmpq_key.d = q;
                    lmpq_key.rebuild_reference();

                    // The full tidal mode is
                    //   omega_lmpq = (l - 2p) periastron_dot + (l - 2p + q) n + m (node_dot - spin),
                    // which reduces to the form below once periastron_dot and node_dot are taken as zero.
                    // Periapse and node precession are not modeled, so both rates are zero here.
                    const int n_coeff = lmpq_key.a - 2 * lmpq_key.c + lmpq_key.d;
                    const int o_coeff = -lmpq_key.b;
                    const double d_n_coeff = static_cast<double>(n_coeff);
                    const double d_o_coeff = static_cast<double>(o_coeff);
                    const double mode =
                        d_n_coeff * orbital_frequency +
                        d_o_coeff * spin_frequency;
                    double mode_sign = 1.0;
                    if (mode < 0.0)
                    {
                        mode_sign = -1.0;
                    }

                    // Records the mode's frequency; a negative index marks a zero frequency.
                    const std::ptrdiff_t frequency_index = c_record_unique_frequencies(
                        lmpq_key,
                        std::abs(mode),
                        frequency_tracker,
                        result.unique_freq_index_map,
                        result.unique_freq_map
                    );

                    if (frequency_index >= 0)
                    {
                        // mode_strength doubles as the common coefficient G^2 times the lmp coeff (which
                        // carries F^2), so a user can see which modes matter and lower a truncation level
                        // that is higher than the problem needs.
                        double common_coeff = G_lpq_squared * lmp_coeff;

                        // The per-mode strength keeps the sign of the tidal mode.
                        const c_ModeStorage mode_storage(mode, mode_sign * common_coeff, n_coeff, o_coeff);
                        result.mode_map.set(lmpq_key, mode_storage);

                        max_mode_strength = std::max(max_mode_strength, std::abs(mode_storage.mode_strength));

                        // Each potential component has a different coefficient but all share the common one.
                        // See Eq. 7 in Renaud+ (2021; PSJ).
                        result.potential_map.set(
                            lmpq_key,
                            c_GlobalPotentialResultAtMode(
                                // dU_dM (mean anomaly)
                                d_n_coeff * mode_sign * common_coeff,
                                // dU_dw (periapsis)
                                (d_n_coeff - static_cast<double>(lmpq_key.d)) * mode_sign * common_coeff,
                                // dU_dSig (node)
                                -d_o_coeff * mode_sign * common_coeff,
                                // Heating
                                std::abs(mode) * host_mass * common_coeff,
                                // dU_dM - dU_dw, exactly: the (l - 2p) parts cancel before any rounding
                                static_cast<double>(lmpq_key.d) * mode_sign * common_coeff,
                                static_cast<size_t>(frequency_index)
                            )
                        );
                    }
                }
            }
        }
    }

    // Normalize the mode strengths against the strongest mode.
    if (max_mode_strength > 0.0)
    {
        for (auto& [lmpq_key, mode_data] : result.mode_map)
        {
            mode_data.mode_strength /= max_mode_strength;
        }
    }

    return result;
}
