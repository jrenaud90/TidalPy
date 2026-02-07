#pragma once

#include <cmath>

#include "obliquity_driver_.hpp"
#include "eccentricity_driver_.hpp"
#include "potential_common_.hpp"
#include "constants_.hpp"


struct PotentialResultAtMode
{
    double dU_dM;
    double dU_dw;
    double dU_dO;
    double E_dot;
    PotentialResultAtMode(
            double dU_dM_,
            double dU_dw_,
            double dU_dO_,
            double E_dot_) :
        dU_dM(dU_dM_),
        dU_dw(dU_dw_),
        dU_dO(dU_dO_),
        E_dot(E_dot_)
    {
    }
};


struct GlobalPotentialResult
{
    ModeMap mode_map;
    UniqueFreqIndexMap unique_freq_index_map;
    UniqueFreqMap unique_freq_map;
    c_IntMap<c_Key4, PotentialResultAtMode> potential_map;
    int error_code = 0;
    int working_on_l = -1;
};


GlobalPotentialResult global_potential(
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
    // Setup output
    GlobalPotentialResult result;
    result.error_code = 0;

    // We can determine an upper bound on the size of our arrays. This will be overkill because some modes will be 
    //  skipped.
    int target_size = 0;
    for (size_t degree_l = min_degree_l; degree_l < (max_degree_l + 1); degree_l++)
    {
        target_size += (degree_l + 1) * (degree_l + 1);
    }
    target_size *= (2 * eccentricity_truncation + 1);
    result.mode_map.reserve(target_size);
    result.unique_freq_index_map.reserve(target_size);
    result.unique_freq_map.reserve(target_size);
    result.potential_map.reserve(target_size);
    // Currently the size of the structure is a bit larger than:
    // IntMaps are  vector(3*8) + (key(8) + sizeof(value)) * size
    //   mode_map = 3*8 + (8 + 8 + 8) * size
    // + unique_freq_index_map = 3*8 + (8 + 8 + 8) * size
    // + unique_freq_map = 3*8 + 8 * size
    // + potential_map = 3*8 + (8 + 8 * 4) * size
    // So stack memory is going to be around 96 bytes
    // Heap memory for l=2 at eccentricity truncation of 6: > 11.23 KB; l=3 > 31.2 KB; l=4 > 62.4 KB

    // For later calculation of the maximum relative mode.
    double max_mode_strength = 0;

    // Optimizations
    c_Key2 lm_key   = c_Key2();
    c_Key4 lmpq_key = c_Key4();
    auto& lm_coeff_map = get_lm_coeff_map();

    // Setup R/a coeff
    double R_a = planet_radius / semi_major_axis;
    double R_a_2 = R_a * R_a;
    double ra_l_coeff = 0;
    switch (min_degree_l)
    {
    case 1:
        ra_l_coeff = R_a_2 * R_a;
    case 2:
        ra_l_coeff = R_a_2 * R_a_2 * R_a;  // 2(l=2) + 1 == 5
        break;
    case 3:
        ra_l_coeff = R_a_2 * R_a_2 * R_a_2 * R_a;  // 2(l=3) + 1 == 7
        break;
    default:
        ra_l_coeff *= std::pow(R_a, 2.0 * static_cast<double>(min_degree_l) + 1);
        break;
    }
    // Multiple the ra_l coeff by the outer most coefficient now. Slightly inefficient doing it now but it will 
    //  allow us having to multiple it outside of this function (outside of this function the focus should only
    //  be on the Love number multiplier).
    ra_l_coeff *= G_to_use * host_mass / semi_major_axis;

    // Step through each degree l and find the potential.
    for (int degree_l = min_degree_l; degree_l < (max_degree_l + 1); degree_l++)
    {
        // Set the degree l we are working on for error reporting.
        result.working_on_l = degree_l;

        // Determine Obliquity functions
        ObliquityFuncOutput obliquity_funcs = c_obliquity_func(
            &result.error_code,
            obliquity,
            degree_l,
            obliquity_truncation
        );
        if (result.error_code != 0)
        {
            // Error, return early.
            return result;
        }

        // Determine Eccentricity functions
        EccentricityFuncOutput eccentricity_funcs = c_eccentricity_func(
            &result.error_code,
            eccentricity,
            degree_l,
            eccentricity_truncation
        );
        if (result.error_code != 0)
        {
            // Error, return early.
            return result;
        }
        // Correct the R_a coeff if we are not at the min degree_l (it has already been initialized for that).
        if (degree_l > min_degree_l)
        {
            // Every sequential degree l grows the coeff by (R/a)^2
            ra_l_coeff *= R_a_2;
        }

        // The combination of obliquity and eccentricity functions determines the number of unique modes required 
        //  for the tidal potential.

        // Set up l in various keys
        lm_key.a = degree_l;
        // Set m to a negative (which is not physical) so we can check if it changes later
        lm_key.b = -1;
        lmpq_key.a = degree_l;
        
        // Prepare the lm_coeff
        double lm_coeff = TidalPyConstants::d_NAN;

        // Step through the outer lmp loop defined by the obliquity function.
        for (const auto& [lmp_key, F_lmp] : obliquity_funcs.first) {
            bool found = false;
            if (lmp_key.b != lm_key.b)
            {
                // New m.
                lm_key.b = lmp_key.b;
                lm_key.rebuild_reference();
                
                // Get (l - m)! / (l + m)! (2 - d_m0) coefficient
                lm_coeff = lm_coeff_map.get(found, lm_key);
                if (!found)
                {
                    // Can not find l,m combo. Perhaps unsupported degree l.
                    result.error_code = -20;
                    return result;
                }

                // Also update the m for our lmpq key
                lmpq_key.b = lmp_key.b;
            }
            lmpq_key.c = lmp_key.c;

            // We will need F^2 for the global potential (local potential only used F).
            // We also need a common coefficient of (l - m)! / (l + m)!
            // Might as well multiple this by the F^2 term.
            double lmp_coeff = F_lmp * F_lmp * ra_l_coeff * lm_coeff;

            // Now use the current l and p to find the eccentricity vector of G_lpq results sorted by q.
            found = false;
            const c_IntMap<c_Key1, double>* eccentricity_by_q_ptr = 
                eccentricity_funcs.second.get_ptr(found, c_Key2(lmp_key.a, lmp_key.c));  // a == l; b == m; c == p
            
            // It is possible that the eccentricity_by_q has no data for this l,p; meaning that G_lp(q) = 0 for all q.
            //   Continue to the next mode if that is the case. Otherwise, carry on.
            if (found)
            {
                for (const auto& [q_key, G_lpq] : *eccentricity_by_q_ptr)
                {
                    // Build key for this unique lmpq
                    lmpq_key.d = q_key.a;
                    lmpq_key.rebuild_reference();

                    // Find the tidal mode at this lmpq. The full equation for tidal mode is:
                    //  $\omega_{lmpq} = (l - 2p) * periastron_dot + (l - 2p + q) * orbital_motion + m * (node_dot - spin_freq)
                    // If we assume node_dot ~ 0 and periastron_dot ~ 0; this reduces to the following.
                    // TODO: Perhaps put a function here to switch between cases. But then we'd need to track these other
                    // parameters. Perhaps better to make whole new potential functions.
                    ModeStorage mode_storage = ModeStorage(
                        lmpq_key.a - 2 * lmpq_key.c + lmpq_key.d,
                        -lmpq_key.b
                    );
                    mode_storage.mode = 
                        static_cast<double>(mode_storage.n_coeff) * orbital_frequency + 
                        static_cast<double>(mode_storage.o_coeff) * spin_frequency;
                    double mode_sign = 1.0;
                    if (mode_storage.mode < 0)
                    {
                        mode_sign = -1.0;
                    }

                    // Use the mode to build up our frequency arrays and also check if we should skip this one.
                    bool nonzero_freq = record_unique_frequencies(                
                        lmpq_key,
                        std::abs(mode_storage.mode),
                        result.unique_freq_index_map,
                        result.unique_freq_map
                    );

                    if (nonzero_freq)
                    {
                        // Actually start calculating some potentials! First let's save this mode as an active one.

                        // We want to track which modes are actually important so the user can adjust truncation levels
                        //  if they are too high for the problem at hand.
                        // Instead of creating its own parameter, just use mode strength as our common coefficient.
                        // Like the obliquity function, we will need G^2 and lets take this opportunity to multiple by F^2
                        double common_coeff = G_lpq * G_lpq * lmp_coeff;
                        mode_storage.mode_strength = mode_sign * common_coeff;
                        max_mode_strength = std::max(max_mode_strength, std::abs(mode_storage.mode_strength));
                        result.mode_map.set(lmpq_key, mode_storage);

                        // Each potential component has a different coefficient but all share the common one.
                        // See Eq. 7 in Renaud+ (2021; PSJ).
                        result.potential_map.set(lmpq_key,
                            PotentialResultAtMode(
                                static_cast<double>(mode_storage.n_coeff) * mode_sign * common_coeff, // dU_dM (mean anomaly)
                                static_cast<double>(mode_storage.n_coeff - lmpq_key.d) * mode_sign * common_coeff,  // dU_dw (periapsis)
                                -static_cast<double>(mode_storage.o_coeff) * mode_sign * common_coeff, // dU_dSig (node)
                                std::abs(mode_storage.mode) * host_mass * common_coeff  // Heating 
                            )
                        );
                    }
                }
            }
        }
    }

    // Step through our modes to calculate the relative mode strength
    if (max_mode_strength > 0.0)
    {
        for (auto& [lmpq_key, mode_data] : result.mode_map)
        {
            mode_data.mode_strength /= max_mode_strength;
        }
    }

    return result;
}
