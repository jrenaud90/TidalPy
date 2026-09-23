#pragma once

#include <cmath>

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
    c_GlobalPotentialResultAtMode(
            double dU_dM_,
            double dU_dw_,
            double dU_dO_,
            double E_dot_,
            double dU_dM_minus_dw_) :
        dU_dM(dU_dM_),
        dU_dw(dU_dw_),
        dU_dO(dU_dO_),
        E_dot(E_dot_),
        dU_dM_minus_dw(dU_dM_minus_dw_)
    {
    }
};


struct c_GlobalPotentialStorage
{
    c_ModeMap mode_map;
    c_UniqueFreqIndexMap unique_freq_index_map;
    c_UniqueFreqMap unique_freq_map;
    c_IntMap<c_Key4, c_GlobalPotentialResultAtMode> potential_map;
    int error_code = 0;
    int working_on_l = -1;
};

struct c_GlobalPotentialResult
{
    c_ModeMap mode_map;
    c_UniqueFreqIndexMap unique_freq_index_map;
    c_UniqueFreqMap unique_freq_map;
    c_IntMap<c_Key4, c_GlobalPotentialResultAtMode> potential_map;
    int error_code = 0;
};

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
        int eccentricity_truncation
    )
{
    c_GlobalPotentialStorage result;
    result.error_code = 0;

    // Upper bound on the number of modes; an overestimate, since some modes are skipped.
    int target_size = 0;
    for (int degree_l = min_degree_l; degree_l <= max_degree_l; degree_l++)
    {
        target_size += (degree_l + 1) * (degree_l + 1);
    }
    target_size *= (2 * eccentricity_truncation + 1);
    result.mode_map.reserve(target_size);
    result.unique_freq_index_map.reserve(target_size);
    result.unique_freq_map.reserve(target_size);
    result.potential_map.reserve(target_size);
    // The four maps take about 96 bytes on the stack; the heap cost at eccentricity truncation 6 is roughly
    // 11 kB at l = 2, 31 kB at l = 3, and 62 kB at l = 4.

    // For later calculation of the maximum relative mode.
    double max_mode_strength = 0;

    c_Key2 lm_key   = c_Key2();
    c_Key4 lmpq_key = c_Key4();
    auto& lm_coeff_map = c_get_lm_coeff_map();

    double R_a = planet_radius / semi_major_axis;
    double R_a_2 = R_a * R_a;
    double ra_l_coeff = 0;
    switch (min_degree_l)
    {
    case 1:
        ra_l_coeff = R_a_2 * R_a;
        break;
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

        ObliquityFuncOutput obliquity_funcs = c_obliquity_func(
            &result.error_code,
            obliquity,
            degree_l,
            obliquity_truncation
        );
        if (result.error_code != 0)
        {
            return result;
        }

        EccentricityFuncOutput eccentricity_funcs = c_eccentricity_func(
            &result.error_code,
            eccentricity,
            degree_l,
            eccentricity_truncation
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


        lm_key.a = degree_l;
        // Unphysical m sentinel, so the first pass through the loop registers as a new m.
        lm_key.b = -1;
        lmpq_key.a = degree_l;
        
        double lm_coeff = TidalPyConstants::d_NAN;

        for (const auto& [lmp_key, F_lmp] : obliquity_funcs.first) {

            if (F_lmp == 0.0)
            {
                continue;
            }
            
            bool found = false;
            if (lmp_key.b != lm_key.b)
            {
                lm_key.b = lmp_key.b;
                lm_key.rebuild_reference();
                
                // Get (l - m)! / (l + m)! (2 - d_m0) coefficient
                lm_coeff = lm_coeff_map.get(found, lm_key);
                if (!found)
                {
                    // No (l, m) coefficient; likely an unsupported degree l.
                    result.error_code = -20;
                    return result;
                }

                lmpq_key.b = lmp_key.b;
            }
            lmpq_key.c = lmp_key.c;

            // The global potential goes as F^2 (the 3D path uses F); fold in (l - m)!/(l + m)!(2 - d_m0).
            double lmp_coeff = F_lmp * F_lmp * ra_l_coeff * lm_coeff;

            found = false;
            const c_IntMap<c_Key1, double>* eccentricity_by_q_ptr = 
                eccentricity_funcs.second.get_ptr(found, c_Key2(lmp_key.a, lmp_key.c));  // a == l; b == m; c == p
            
            // No entry for this (l, p) means G_lpq = 0 for every q.
            if (found)
            {
                for (const auto& [q_key, G_lpq] : *eccentricity_by_q_ptr)
                {
                    if (G_lpq == 0.0)
                    {
                        continue;
                    }

                    lmpq_key.d = q_key.a;
                    lmpq_key.rebuild_reference();

                    // The full tidal mode is
                    //   omega_lmpq = (l - 2p) periastron_dot + (l - 2p + q) n + m (node_dot - spin),
                    // which reduces to the form below once periastron_dot and node_dot are taken as zero.
                    // TODO: support nonzero periapse and node precession.
                    c_ModeStorage mode_storage = c_ModeStorage(
                        lmpq_key.a - 2 * lmpq_key.c + lmpq_key.d,  // n coeff 
                        -lmpq_key.b                                // o coeff
                    );
                    const double d_n_coeff = static_cast<double>(mode_storage.n_coeff);
                    const double d_o_coeff = static_cast<double>(mode_storage.o_coeff);
                    mode_storage.mode = 
                        d_n_coeff * orbital_frequency + 
                        d_o_coeff * spin_frequency;
                    double mode_sign = 1.0;
                    if (mode_storage.mode < 0.0)
                    {
                        mode_sign = -1.0;
                    }

                    // Records the mode's frequency and reports whether it is nonzero.
                    bool nonzero_freq = record_unique_frequencies(                
                        lmpq_key,
                        std::abs(mode_storage.mode),
                        result.unique_freq_index_map,
                        result.unique_freq_map
                    );

                    if (nonzero_freq)
                    {
                        // mode_strength doubles as the common coefficient G^2 times the lmp coeff (which
                        // carries F^2), so a user can see which modes matter and lower a truncation level
                        // that is higher than the problem needs.
                        double common_coeff = G_lpq * G_lpq * lmp_coeff;

                        // The per-mode strength keeps the sign of the tidal mode.
                        mode_storage.mode_strength = mode_sign * common_coeff;
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
                                std::abs(mode_storage.mode) * host_mass * common_coeff,
                                // dU_dM - dU_dw, exactly: the (l - 2p) parts cancel before any rounding
                                static_cast<double>(lmpq_key.d) * mode_sign * common_coeff
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
