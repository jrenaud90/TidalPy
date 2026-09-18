#pragma once
/*
 * world_tides_.hpp: out-of-line definition of c_LayeredWorld::calc_tides and the 3D tidal paths.
 *
 * c_LayeredWorld extends the common analytic tide path (c_BaseWorld::calc_tides in world_tides_base_.hpp) with
 * two layered-world capabilities: the rheology model, whose -Im[k_l(omega)] comes from the world radial solver
 * run at each unique tidal frequency (the EOS must be solved first), and the distribution of the heating to the
 * layers by each layer's tidal_scale_method. The tide-model holder, config, and result state stay on
 * c_BaseWorld.
 *
 * This header pulls in the heavy global-potential tables; force-include it in the layered and gas-giant world
 * extension only.
 */

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <vector>

#include "layered_.hpp"
#include "../../Tides_x/classes/tide_collapse_.hpp"      // c_global_potential, c_collapse_global_tides
#include "../../Tides_x/classes/tide_.hpp"               // c_RheologyTide (3D orchestration target)
#include "../../Tides_x/potential/potential_3d_.hpp"     // c_tidal_potential_3d_modes (dynamic Kaula engine)
#include "../../Tides_x/multilayer/kernel_.hpp"          // strain/stress/heating kernel + c_StrainRadialCoeffs
#include "../../Tides_x/multilayer/angular_collapse_.hpp" // c_theta_integrated_heating (analytic colatitude collapse)

namespace tidalpy {

inline void c_LayeredWorld::calc_tides(const c_TideSolveConfig& state) {
    if (!this->p_tide) {
        throw std::runtime_error(
            "TidalPy: no tide model attached to the world — call set_tide_model() first");
    }

    const double planet_radius = this->get_radius();
    const double G_to_use = c_get_G();
    const c_TideConfig& tcfg = this->p_tide_config;

    // Model-independent per-mode terms + the unique-frequency maps.
    c_GlobalPotentialStorage potential = c_global_potential(
        planet_radius,
        state.semi_major_axis,
        state.orbital_frequency,
        state.spin_frequency,
        state.obliquity,
        state.eccentricity,
        state.host_mass,
        G_to_use,
        tcfg.min_degree_l,
        tcfg.max_degree_l,
        tcfg.obliquity_truncation,
        tcfg.eccentricity_truncation
    );

    if (potential.error_code != 0) {
        this->p_tides_solved           = false;
        this->p_tide_result            = c_GlobalTideResult();
        this->p_tide_result.error_code = potential.error_code;
        throw std::runtime_error("TidalPy: global potential failed during calc_tides");
    }

    // Collapse the per-mode potential terms with the tide model's dissipation multiplier.
    this->p_tide_solver_love.clear();
    if (this->p_tide->needs_radial_solve()) {
        // Rheology model: the per-mode -Im[k_l(omega)] comes from the world radial solver,
        // which needs the EOS solved first.
        if (!this->p_eos_solved || !this->p_eos_solution) {
            this->p_tides_solved = false;
            throw std::runtime_error(
                "TidalPy: the rheology tide model needs the EOS solved first — call "
                "solve_eos() before calc_tides()");
        }

        // Solve the Love numbers once per unique (degree_l, frequency) pair, caching by that
        // pair so modes that share a degree and frequency reuse one radial solve. Then record
        // each active mode's Love numbers keyed by its (l, m, p, q).
        c_IntMap<c_Key2, tidalpy::c_LoveNumbers> love_by_l_freq;
        c_LoveSolveConfig love_cfg = this->make_love_solve_config();

        // The homogeneous Love methods reuse their node values and per-frequency averages across this call's solves.
        c_HomogeneousLoveCache homogeneous_cache;
        for (const auto& mode_entry : potential.potential_map) {
            const c_Key4& lmpq_key = mode_entry.first;
            const int degree_l     = static_cast<int>(lmpq_key.a);

            bool found = false;
            const std::size_t freq_index = potential.unique_freq_index_map.get(found, lmpq_key);
            if (!found) {
                // Inactive (zero-frequency) mode; contributes nothing.
                continue;
            }
            const double frequency = potential.unique_freq_map[freq_index].frequency;

            c_Key2 lf_key(static_cast<int16_t>(degree_l), static_cast<int16_t>(freq_index));
            bool cached = false;
            love_by_l_freq.get(cached, lf_key);
            if (!cached) {
                love_cfg.degree_l        = degree_l;
                love_cfg.frequency = frequency;
                this->solve_love_numbers(love_cfg, &homogeneous_cache);
                if (!this->get_love_success()) {
                    this->p_tides_solved = false;
                    throw std::runtime_error(
                        "TidalPy: radial-solver Love-number solve failed during calc_tides: "
                        + this->get_love_message());
                }
                tidalpy::c_LoveNumbers solved_love;
                solved_love.k = this->get_love_number_k(0);
                solved_love.h = this->get_love_number_h(0);
                solved_love.l = this->get_love_number_l(0);
                love_by_l_freq.set(lf_key, solved_love);
            }

            bool have = false;
            tidalpy::c_LoveNumbers mode_love = love_by_l_freq.get(have, lf_key);
            this->p_tide_solver_love.set(lmpq_key, mode_love);
        }

        this->p_tide_result = c_collapse_global_tides(potential, *this->p_tide, &this->p_tide_solver_love);
    } else {
        // Analytic models: no radial-solver Love numbers needed.
        this->p_tide_result = c_collapse_global_tides(potential, *this->p_tide, nullptr);
    }
    this->p_tides_solved = true;

    // Distribute heating to the layers by their effective tidal scale (0 for non-tidal
    // layers), and store the result on each layer so layer.get_tidal_heating() reports it.
    const double planet_volume =
        (4.0 / 3.0) * TidalPyConstants::d_PI * planet_radius * planet_radius * planet_radius;
    const std::size_t n_layers = this->p_layers.size();
    this->p_layer_tidal_heating.assign(n_layers, 0.0);
    for (std::size_t i = 0; i < n_layers; ++i) {
        c_BaseLayer* layer = this->p_layers[i].get();
        const double scale = this->effective_tidal_scale(layer, planet_volume, state);
        const double heat  = this->p_tide_result.tidal_heating * scale;
        this->p_layer_tidal_heating[i] = heat;
        layer->set_tidal_heating(heat);
    }
}

// Effective per-layer tidal-heating scale for the layer's tidal_scale_method (0 for a non-tidal layer).
//   user_provided   : the layer's tidal_scale field.
//   volume_fraction : layer volume / planet volume.
//   tidal_timescale : a log-Gaussian bell in the layer's Maxwell time tau = eta/mu (from its static shear
//                     modulus and viscosity) about the tidal forcing period 2*pi/|orbital_frequency|, with the
//                     width [decades] from the tide config. Returns 0 for a geometry-only layer or when mu,
//                     eta, or the forcing are unusable.
inline double c_LayeredWorld::effective_tidal_scale(
        const c_BaseLayer* layer, double planet_volume, const c_TideSolveConfig& state) const {
    if (!layer->get_is_tidal()) {
        return 0.0;
    }
    switch (layer->get_tidal_scale_method()) {
        case c_TidalScaleMethod::user_provided:
            return layer->get_tidal_scale();
        case c_TidalScaleMethod::volume_fraction:
            return (planet_volume > TidalPyConstants::d_EPS)
                 ? layer->get_volume() / planet_volume : 0.0;
        case c_TidalScaleMethod::tidal_timescale: {
            const auto* phys = dynamic_cast<const c_PhysicsLayer*>(layer);
            if (phys == nullptr) {
                return 0.0;   // geometry-only layer has no Maxwell time
            }
            const double shear_modulus  = phys->get_shear_modulus_static();
            const double shear_viscosity = phys->get_shear_viscosity_static();
            const double orbital_freq   = std::abs(state.orbital_frequency);
            if (!std::isfinite(shear_modulus) || shear_modulus <= TidalPyConstants::d_EPS
             || !std::isfinite(shear_viscosity) || shear_viscosity <= TidalPyConstants::d_EPS
             || orbital_freq <= TidalPyConstants::d_EPS) {
                return 0.0;
            }
            const double maxwell_time   = shear_viscosity / shear_modulus;          // [s]
            const double forcing_period = 2.0 * TidalPyConstants::d_PI / orbital_freq;  // [s]
            double width = this->p_tide_config.tidal_timescale_width_decades;
            if (width <= TidalPyConstants::d_EPS) {
                width = 1.0;
            }
            const double z = std::log10(maxwell_time / forcing_period) / width;
            return std::exp(-0.5 * z * z);
        }
        default:
            return 0.0;
    }
}


// =====================================================================================================================
// On-demand 3D tidal heating
// =====================================================================================================================
//
// Every 3D path works from the coherent wave list (c_coherent_tidal_waves_3d, potential_3d_.hpp): each active
// (l, m, p, q) mode is mapped onto its non-negative frequency and merged with the modes that share its real spatial
// function. The radial (Love-number) solve depends on (degree l, |omega|) only, so it runs once per unique pair and
// its strain radial coefficients are reused across waves, points, longitudes, and times.
//
// Secular (cycle and orbit-averaged) heating density:
//     h_bar(r, theta, phi) = sum over |omega| of (|omega|/2) Im( sigma_c : conj(eps_c) )
// with sigma_c, eps_c the total complex stress and strain amplitude at that frequency, every wave at that |omega|
// summed before the bilinear form. Cross terms between different frequencies average to zero over the orbit and are
// dropped; those between waves at the same frequency survive the average and are kept. They are what the m = 0 pairs
// contribute (each pair is one real sinusoid) and what makes the heating of a synchronously rotating body, whose
// active modes all sit at multiples of n, depend on longitude. The scalar and batch paths take no longitude and
// return the longitudinal mean of h_bar: cross terms between waves with different azimuthal structure e^{i mu phi}
// integrate to zero over phi, so the mean is the sum over (|omega|, mu) groups evaluated at phi = 0. The volume
// integral of h_bar is the 1D global heating (get_tidal_heating). At nonzero obliquity the 3D value is for the
// geometry with zero argument of periapse and node (the engine drops precession): same-frequency modes of the same
// (l, m) then combine coherently, a cross term the precession-averaged 1D formula does not carry, so the two agree
// only to the size of those terms there.

namespace tides3d {

struct c_RadialGroup3D {
    int degree_l = 0;
    double frequency = 0.0;   // |omega| [rad s-1]
};

// The coherent waves plus the two groupings the 3D paths need.
struct c_WaveSet3D {
    std::vector<c_TidalWave3D> waves;
    std::vector<c_RadialGroup3D> radial_groups;   // unique (degree l, |omega|): one radial solve each
    std::vector<int> wave_radial_group;           // per wave, index into radial_groups
    std::vector<double> frequencies;              // unique |omega|: waves sharing one combine coherently
    std::vector<int> wave_frequency_group;        // per wave, index into frequencies
    std::vector<std::array<int, 2>> angular_pairs; // unique (degree l, order m): one Legendre evaluation each
    std::vector<int> wave_angular_pair;           // per wave, index into angular_pairs
    std::vector<int> azimuthal_orders;            // unique signed azimuthal wavenumber mu = azimuthal_sign * m
    std::vector<int> wave_azimuthal_order;        // per wave, index into azimuthal_orders
};

inline c_WaveSet3D c_build_wave_set_3d(
        const std::vector<c_TidalPotential3DModeCoeff>& modes,
        double min_frequency) {
    c_WaveSet3D set;
    set.waves = c_coherent_tidal_waves_3d(modes, min_frequency);
    const size_t num_waves = set.waves.size();
    set.wave_radial_group.assign(num_waves, -1);
    set.wave_frequency_group.assign(num_waves, -1);
    set.wave_angular_pair.assign(num_waves, -1);
    set.wave_azimuthal_order.assign(num_waves, -1);
    for (size_t w = 0; w < num_waves; ++w) {
        const c_TidalWave3D& wave = set.waves[w];
        int radial_group = -1;
        for (size_t g = 0; g < set.radial_groups.size(); ++g) {
            if (set.radial_groups[g].degree_l == wave.degree_l
                && c_tidal_wave_same_frequency(set.radial_groups[g].frequency, wave.frequency)) {
                radial_group = static_cast<int>(g);
                break;
            }
        }
        if (radial_group < 0) {
            radial_group = static_cast<int>(set.radial_groups.size());
            set.radial_groups.push_back(c_RadialGroup3D{wave.degree_l, wave.frequency});
        }
        set.wave_radial_group[w] = radial_group;

        int frequency_group = -1;
        for (size_t f = 0; f < set.frequencies.size(); ++f) {
            if (c_tidal_wave_same_frequency(set.frequencies[f], wave.frequency)) {
                frequency_group = static_cast<int>(f);
                break;
            }
        }
        if (frequency_group < 0) {
            frequency_group = static_cast<int>(set.frequencies.size());
            set.frequencies.push_back(wave.frequency);
        }
        set.wave_frequency_group[w] = frequency_group;

        // Waves of one (l, m) share their Legendre values at a colatitude, and waves of one mu their e^{i mu phi}.
        const std::array<int, 2> pair{wave.degree_l, wave.order_m};
        const auto pair_found = std::find(set.angular_pairs.begin(), set.angular_pairs.end(), pair);
        set.wave_angular_pair[w] = static_cast<int>(pair_found - set.angular_pairs.begin());
        if (pair_found == set.angular_pairs.end()) { set.angular_pairs.push_back(pair); }
        const int mu = wave.azimuthal_sign * wave.order_m;
        const auto mu_found = std::find(set.azimuthal_orders.begin(), set.azimuthal_orders.end(), mu);
        set.wave_azimuthal_order[w] = static_cast<int>(mu_found - set.azimuthal_orders.begin());
        if (mu_found == set.azimuthal_orders.end()) { set.azimuthal_orders.push_back(mu); }
    }
    return set;
}

// The wave set for the world's tide config and an orbital state.
inline c_WaveSet3D c_world_wave_set_3d(
        c_LayeredWorld& world,
        const c_TideSolveConfig& state,
        const char* what) {
    const c_TideConfig& tide_cfg = world.get_tide_config();
    int engine_error = 0;
    const std::vector<c_TidalPotential3DModeCoeff> modes = c_tidal_potential_3d_mode_coeffs(
        world.get_radius(),
        state.semi_major_axis,
        state.orbital_frequency,
        state.spin_frequency,
        state.obliquity,
        state.eccentricity,
        state.host_mass,
        c_get_G(),
        tide_cfg.min_degree_l,
        tide_cfg.max_degree_l,
        tide_cfg.obliquity_truncation,
        tide_cfg.eccentricity_truncation,
        &engine_error);
    if (engine_error != 0) {
        throw std::runtime_error(
            std::string("TidalPy: tidal potential engine failed during ") + what + " (error "
            + std::to_string(engine_error) + "); check degree/truncation levels");
    }
    const double min_freq =
        (tidalpy_config_ptr != nullptr) ? tidalpy_config_ptr->d_MIN_SPIN_ORBIT_DIFF : 1.0e-9;
    return c_build_wave_set_3d(modes, min_freq);
}

// Run the world radial solve for one radial group and return its solution storage.
inline const ::c_RadialSolutionStorage* c_solve_radial_group_3d(
        c_LayeredWorld& world,
        c_LoveSolveConfig& love_cfg,
        const c_RadialGroup3D& group,
        const char* what) {
    love_cfg.degree_l = group.degree_l;
    love_cfg.frequency = group.frequency;
    world.solve_love_numbers(love_cfg);
    if (!world.get_love_success()) {
        throw std::runtime_error(
            std::string("TidalPy: radial solve failed during ") + what + ": " + world.get_love_message());
    }
    const ::c_RadialSolutionStorage* storage = world.get_love_storage();
    if (storage == nullptr) {
        throw std::runtime_error(std::string("TidalPy: missing radial solution during ") + what);
    }
    return storage;
}

// Strain radial coefficients of one radial group at one radius from the world's current radial solution.
// Returns false when the radius has no depth-resolved solution (center / below the solver start); a liquid
// layer returns true with out.valid == false (no shear kernel there: contributes nothing).
inline bool c_strain_coeffs_at_radius_3d(
        c_LayeredWorld& world,
        const ::c_RadialSolutionStorage* storage,
        double radius,
        const c_RadialGroup3D& group,
        tides::c_StrainRadialCoeffs& out) {
    std::complex<double> y_at_r[C_MAX_NUM_Y];
    if (!storage->get_radial_solution(radius, 0, y_at_r)
        || !std::isfinite(y_at_r[0].real()) || !std::isfinite(y_at_r[1].real())
        || !std::isfinite(y_at_r[2].real()) || !std::isfinite(y_at_r[3].real())) {
        out.valid = false;
        return false;
    }
    bool is_solid = true;
    bool is_incompressible = false;
    const auto* physics_layer = dynamic_cast<const c_PhysicsLayer*>(world.find_layer_for_radius(radius));
    if (physics_layer != nullptr) {
        is_solid = physics_layer->get_is_solid();
        is_incompressible = physics_layer->get_is_incompressible();
    }
    const std::complex<double> shear = world.calc_complex_shear_modulus(radius, group.frequency);
    const std::complex<double> bulk  = world.calc_complex_bulk_modulus(radius, group.frequency);
    out = tides::c_compute_strain_radial_coeffs(
        y_at_r[0],
        y_at_r[1],
        y_at_r[2],
        y_at_r[3],
        shear,
        bulk,
        radius,
        static_cast<double>(group.degree_l),
        is_solid,
        is_incompressible);
    return true;
}

// Strain radial coefficients at a set of radii for every radial group of a wave set, and which radii are unusable.
struct c_RadialCoefficients3D {
    std::vector<std::vector<tides::c_StrainRadialCoeffs>> by_radius;   // [radius][radial group]
    std::vector<unsigned char> radius_failed;                          // 1 where no radial group has a solution
};

// Solve the radial problem once per radial group (l, |omega|) and evaluate each group's strain radial coefficients
// at every radius. A radius is unusable only when no group has a depth-resolved solution there: the solver's start
// radius grows with degree l, so a higher-degree group whose solution starts further out simply contributes nothing
// below it.
inline c_RadialCoefficients3D c_radial_coefficients_3d(
        c_LayeredWorld& world,
        const c_WaveSet3D& set,
        const double* radii,
        size_t num_radii,
        const char* what) {
    const size_t num_groups = set.radial_groups.size();
    c_RadialCoefficients3D out;
    out.by_radius.assign(num_radii, std::vector<tides::c_StrainRadialCoeffs>(num_groups));
    out.radius_failed.assign(num_radii, 0);
    std::vector<size_t> radius_missing(num_radii, 0);
    c_LoveSolveConfig love_cfg = world.make_radial_love_solve_config();
    for (size_t g = 0; g < num_groups; ++g) {
        const ::c_RadialSolutionStorage* storage =
            c_solve_radial_group_3d(world, love_cfg, set.radial_groups[g], what);
        for (size_t ir = 0; ir < num_radii; ++ir) {
            if (!c_strain_coeffs_at_radius_3d(
                world,
                storage,
                radii[ir],
                set.radial_groups[g],
                out.by_radius[ir][g])) {
                radius_missing[ir] += 1;
            }
        }
    }
    for (size_t ir = 0; ir < num_radii; ++ir) {
        out.radius_failed[ir] = (num_groups > 0 && radius_missing[ir] == num_groups) ? 1 : 0;
    }
    return out;
}

// cos(|omega_f| t) and sin(|omega_f| t) for every frequency f of a wave set and every time, row-major
// [f * num_times + it]. They depend on neither position nor wave, so a grid tabulates them once.
struct c_PhaseTable3D {
    std::vector<double> cos_phase;
    std::vector<double> sin_phase;
};

inline c_PhaseTable3D c_phase_table_3d(
        const std::vector<double>& frequencies,
        const double* times,
        size_t num_times) {
    c_PhaseTable3D table;
    table.cos_phase.resize(frequencies.size() * num_times);
    table.sin_phase.resize(frequencies.size() * num_times);
    for (size_t f = 0; f < frequencies.size(); ++f) {
        for (size_t it = 0; it < num_times; ++it) {
            const double phase = frequencies[f] * times[it];
            table.cos_phase[f * num_times + it] = std::cos(phase);
            table.sin_phase[f * num_times + it] = std::sin(phase);
        }
    }
    return table;
}

// Total complex stress and strain amplitude of each frequency of a wave set at one point. active[f] is 1 where
// frequency f received at least one wave.
struct c_FrequencyAmplitudes3D {
    std::vector<tides::c_Tensor6> stress;
    std::vector<tides::c_Tensor6> strain;
    std::vector<unsigned char> active;

    explicit c_FrequencyAmplitudes3D(size_t num_frequencies) :
        stress(num_frequencies),
        strain(num_frequencies),
        active(num_frequencies, 0) {}
};

// The radius-independent part of every wave of a set at one (colatitude, longitude), which a grid forms once per point
// and reuses at every radius. c_wave_angular_colatitude_3d sets the colatitude: the Legendre values of each (l, m) and
// the sine and cotangent factors. c_wave_angular_longitude_3d then sets the longitude: the phasor e^{i mu phi} of each
// mu, every wave's potential point, and, when asked, every wave's angular strain factors.
struct c_WaveAngular3D {
    std::vector<c_LegendreValue> legendre;                  // [angular pair] at the colatitude
    tides::c_ColatitudeTrig trig;                           // sine and cotangent factors of the colatitude
    std::vector<std::complex<double>> phasors;              // [azimuthal order] at the longitude
    std::vector<c_PotentialPointC> potentials;              // [wave] at the point
    std::vector<tides::c_AngularStrainFactors> factors;     // [wave] at the point, when strain factors were asked for

    explicit c_WaveAngular3D(const c_WaveSet3D& set) :
        legendre(set.angular_pairs.size()),
        phasors(set.azimuthal_orders.size()),
        potentials(set.waves.size()),
        factors(set.waves.size()) {}
};

inline void c_wave_angular_colatitude_3d(const c_WaveSet3D& set, double colatitude, c_WaveAngular3D& angular) {
    angular.trig = tides::c_colatitude_trig(colatitude);
    for (size_t pair = 0; pair < set.angular_pairs.size(); ++pair) {
        angular.legendre[pair] = c_legendre(set.angular_pairs[pair][0], set.angular_pairs[pair][1], colatitude);
    }
}

inline void c_wave_angular_longitude_3d(
        const c_WaveSet3D& set,
        double longitude,
        bool strain_factors,
        c_WaveAngular3D& angular) {
    for (size_t order = 0; order < set.azimuthal_orders.size(); ++order) {
        angular.phasors[order] = c_azimuthal_phasor(static_cast<double>(set.azimuthal_orders[order]), longitude);
    }
    for (size_t w = 0; w < set.waves.size(); ++w) {
        angular.potentials[w] = c_wave_point_from_parts(
            set.waves[w],
            angular.legendre[static_cast<size_t>(set.wave_angular_pair[w])],
            angular.phasors[static_cast<size_t>(set.wave_azimuthal_order[w])]);
        if (strain_factors) {
            angular.factors[w] = tides::c_angular_strain_factors(angular.potentials[w], angular.trig);
        }
    }
}

// Fill the amplitudes at one point from its wave angular state (with strain factors) and the strain radial
// coefficients of every radial group at one radius. Every wave with a shear kernel there is added into the total of
// its frequency: waves at one frequency superpose, so every cross term between them is kept. Returns false when no
// wave has a shear kernel (a liquid layer).
inline bool c_frequency_amplitudes_3d(
        const c_WaveSet3D& set,
        const std::vector<tides::c_StrainRadialCoeffs>& coeffs_at_radius,
        const c_WaveAngular3D& angular,
        c_FrequencyAmplitudes3D& amplitudes) {
    std::fill(amplitudes.active.begin(), amplitudes.active.end(), 0);
    bool any_kernel = false;
    for (size_t w = 0; w < set.waves.size(); ++w) {
        const tides::c_StrainRadialCoeffs& radial = coeffs_at_radius[set.wave_radial_group[w]];
        if (!radial.valid) { continue; }
        tides::c_Tensor6 strain, stress;
        tides::c_compute_strain_stress_from_factors(
            radial,
            angular.factors[w],
            strain,
            stress);
        const size_t f = static_cast<size_t>(set.wave_frequency_group[w]);
        if (!amplitudes.active[f]) {
            amplitudes.stress[f] = stress;
            amplitudes.strain[f] = strain;
            amplitudes.active[f] = 1;
        } else {
            for (size_t k = 0; k < 6; ++k) {
                amplitudes.stress[f].c[k] += stress.c[k];
                amplitudes.strain[f].c[k] += strain.c[k];
            }
        }
        any_kernel = true;
    }
    return any_kernel;
}

// Signed azimuthal wavenumber of a wave (its longitude structure is e^{i mu phi}).
inline int c_wave_mu(const c_TidalWave3D& wave) {
    return wave.azimuthal_sign * wave.order_m;
}

// The waves the secular density sums coherently: one group per frequency, or per (frequency, mu) for the longitude
// mean, over which waves of different mu average out. Groups are listed in order of first appearance and their members
// in wave order, the order the sums run in.
struct c_SecularGroup3D {
    int frequency_group = 0;
    std::vector<size_t> waves;
};

inline std::vector<c_SecularGroup3D> c_secular_groups_3d(const c_WaveSet3D& set, bool split_by_mu) {
    std::vector<c_SecularGroup3D> groups;
    const size_t num_waves = set.waves.size();
    std::vector<unsigned char> done(num_waves, 0);
    for (size_t w0 = 0; w0 < num_waves; ++w0) {
        if (done[w0]) { continue; }
        c_SecularGroup3D group;
        group.frequency_group = set.wave_frequency_group[w0];
        const int mu0 = c_wave_mu(set.waves[w0]);
        for (size_t w = w0; w < num_waves; ++w) {
            if (done[w] || set.wave_frequency_group[w] != group.frequency_group) { continue; }
            if (split_by_mu && c_wave_mu(set.waves[w]) != mu0) { continue; }
            done[w] = 1;
            group.waves.push_back(w);
        }
        groups.push_back(std::move(group));
    }
    return groups;
}

// Secular density at one point: each group's waves summed into a total complex stress and strain, then
// (|omega|/2) Im(sigma_c : conj(eps_c)) summed over the groups. coeffs_by_group[g] is radial group g's strain radial
// coefficient set at this radius (valid == false where the group has no shear kernel there), and the wave angular
// state carries strain factors at the point (the longitude mean evaluates it at phi = 0).
inline double c_secular_density_3d(
        const c_WaveSet3D& set,
        const std::vector<c_SecularGroup3D>& groups,
        const std::vector<tides::c_StrainRadialCoeffs>& coeffs_by_group,
        const c_WaveAngular3D& angular) {
    double heating = 0.0;
    for (const c_SecularGroup3D& group : groups) {
        tides::c_Tensor6 stress_total;
        tides::c_Tensor6 strain_total;
        bool any = false;
        for (const size_t w : group.waves) {
            const tides::c_StrainRadialCoeffs& radial = coeffs_by_group[set.wave_radial_group[w]];
            if (!radial.valid) { continue; }
            tides::c_Tensor6 strain;
            tides::c_Tensor6 stress;
            tides::c_compute_strain_stress_from_factors(radial, angular.factors[w], strain, stress);
            for (size_t k = 0; k < 6; ++k) {
                stress_total.c[k] += stress.c[k];
                strain_total.c[k] += strain.c[k];
            }
            any = true;
        }
        if (any) {
            heating += 0.5 * set.frequencies[group.frequency_group]
                     * tides::c_volumetric_heating_signed(stress_total, strain_total);
        }
    }
    return heating;
}

// Run body(task) for every task in 0..num_tasks-1 on up to num_threads threads, the calling thread among them. Tasks
// start in index order but may finish in any order, so a task must write only outputs no other task writes. The first
// exception a task throws stops further tasks from starting and is rethrown on the calling thread once every thread
// has finished. If the system refuses a thread, the threads already running finish the work.
template <typename Body>
inline void c_parallel_tasks_3d(size_t num_tasks, int num_threads, const Body& body) {
    const size_t workers = std::min(num_tasks, static_cast<size_t>(std::max(num_threads, 1)));
    if (workers <= 1) {
        for (size_t task = 0; task < num_tasks; ++task) { body(task); }
        return;
    }
    std::atomic<size_t> next_task(0);
    std::atomic<bool> failed(false);
    std::exception_ptr first_error;
    std::mutex error_mutex;
    const auto work = [&]() {
        for (;;) {
            const size_t task = next_task.fetch_add(1);
            if (task >= num_tasks || failed.load()) { return; }
            try {
                body(task);
            } catch (...) {
                const std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error) { first_error = std::current_exception(); }
                failed.store(true);
                return;
            }
        }
    };
    std::vector<std::thread> threads;
    threads.reserve(workers - 1);
    try {
        for (size_t worker = 1; worker < workers; ++worker) { threads.emplace_back(work); }
    } catch (const std::system_error&) {
        // No more threads available: the ones already started and this one run the remaining tasks.
    }
    work();
    for (std::thread& thread : threads) { thread.join(); }
    if (first_error) { std::rethrow_exception(first_error); }
}

// The axis grids and integration weights of one collapsed 3D heating call: the user axis for each surviving
// dimension and an internal quadrature for each summed one.
struct c_CollapseGrids3D {
    bool instantaneous = false;
    bool any_summed = false;
    bool latitude_full_sphere = true;
    std::vector<double> r_grid, r_wsum;     // r_wsum carries the r^2 Jacobian
    std::vector<size_t> r_layer;            // layer index of each summed radial node
    std::vector<double> th_grid, th_wsum;   // th_wsum absorbs sin theta
    std::vector<double> ph_grid, ph_wsum;
    std::vector<double> t_grid;
};

inline c_CollapseGrids3D c_collapse_grids_3d(
        c_LayeredWorld& world,
        const double* radii,
        size_t num_radii,
        const double* colatitudes,
        size_t num_colatitudes,
        const double* longitudes,
        size_t num_longitudes,
        const double* times,
        size_t num_times,
        const c_Heating3DCollapseConfig& cfg) {
    const double two_pi = 2.0 * TidalPyConstants::d_PI;
    c_CollapseGrids3D grids;
    grids.instantaneous = !cfg.orbit_averaged;
    grids.any_summed = cfg.latitude_summed || cfg.longitude_summed || cfg.radial_summed;

    // Radius: user array unless summed. Summed: Gauss-Legendre nodes inside each layer, with r_wsum carrying the r^2
    // Jacobian. No node sits on a layer boundary, where the modulus and radial-solution lookups take the layer below,
    // so no node can weigh the lower layer's heating into the upper layer's integral.
    if (cfg.radial_summed) {
        const size_t num_layers = world.get_num_layers();
        const int nodes_per_layer = (cfg.radial_slices > 0) ? cfg.radial_slices : 16;
        std::vector<double> gl_r, gl_w;
        tides::c_gauss_legendre_nodes(nodes_per_layer, gl_r, gl_w);
        for (size_t layer_i = 0; layer_i < num_layers; ++layer_i) {
            const c_BaseLayer* layer = world.get_layer(layer_i);
            const double r_mid  = 0.5 * (layer->get_radius_outer() + layer->get_radius_inner());
            const double r_half = 0.5 * (layer->get_radius_outer() - layer->get_radius_inner());
            for (int node = 0; node < nodes_per_layer; ++node) {
                const double rr = r_mid + r_half * gl_r[node];
                grids.r_grid.push_back(rr);
                grids.r_wsum.push_back(r_half * gl_w[node] * rr * rr);  // Gauss-Legendre weight x r^2
                grids.r_layer.push_back(layer_i);
            }
        }
    } else {
        grids.r_grid.assign(radii, radii + num_radii);
    }
    // Colatitude: user array unless summed (Gauss-Legendre in cos theta; the weight absorbs sin theta).
    // A colatitude band [min, max] narrower than [0, pi] maps the nodes onto [cos(max), cos(min)]
    // (affine substitution; the weights scale by the half-width) so the integral covers only the band.
    grids.latitude_full_sphere =
        (cfg.colatitude_min <= TidalPyConstants::d_EPS)
        && (cfg.colatitude_max >= TidalPyConstants::d_PI - TidalPyConstants::d_EPS);
    if (cfg.latitude_summed) {
        if (cfg.colatitude_min < 0.0 || cfg.colatitude_max > TidalPyConstants::d_PI
                || cfg.colatitude_min >= cfg.colatitude_max) {
            throw std::invalid_argument(
                "TidalPy: colatitude band must satisfy 0 <= colatitude_min < colatitude_max <= pi");
        }
        const int num_nodes = (cfg.latitude_nodes > 1) ? cfg.latitude_nodes : 48;
        std::vector<double> gl_x;
        tides::c_gauss_legendre_nodes(num_nodes, gl_x, grids.th_wsum);
        const double x_lo = std::cos(cfg.colatitude_max);
        const double x_hi = std::cos(cfg.colatitude_min);
        const double x_mid  = 0.5 * (x_hi + x_lo);
        const double x_half = 0.5 * (x_hi - x_lo);
        grids.th_grid.resize(num_nodes);
        for (int i = 0; i < num_nodes; ++i) {
            grids.th_grid[i] = std::acos(x_mid + x_half * gl_x[i]);
            grids.th_wsum[i] *= x_half;
        }
    } else {
        grids.th_grid.assign(colatitudes, colatitudes + num_colatitudes);
    }
    // Longitude: user array unless summed. Secular -> analytic 2*pi on the longitude-mean density (single
    // point). Instantaneous -> periodic trapezoid over [0, 2*pi) (uniform weight, no endpoint halving since the
    // field is periodic).
    if (cfg.longitude_summed) {
        if (grids.instantaneous) {
            const int num_nodes = (cfg.longitude_nodes > 1) ? cfg.longitude_nodes : 64;
            const double dph = two_pi / static_cast<double>(num_nodes);
            for (int k = 0; k < num_nodes; ++k) { grids.ph_grid.push_back(k * dph); grids.ph_wsum.push_back(dph); }
        } else {
            grids.ph_grid.push_back(0.0);
            grids.ph_wsum.push_back(two_pi);
        }
    } else {
        grids.ph_grid.assign(longitudes, longitudes + num_longitudes);
    }
    // Time: user array (instantaneous only); a single dummy sample when orbit-averaged.
    if (grids.instantaneous) { grids.t_grid.assign(times, times + num_times); }
    else { grids.t_grid.push_back(0.0); }
    return grids;
}

// The reported axes, output shape, and layer and time counts of a collapsed call (values and layer_totals empty).
inline c_Heating3DCollapsed c_collapse_layout_3d(
        const c_CollapseGrids3D& grids,
        const c_Heating3DCollapseConfig& cfg,
        size_t num_layers) {
    c_Heating3DCollapsed result;
    const bool surv_r = !cfg.radial_summed;
    const bool surv_th = !cfg.latitude_summed;
    const bool surv_ph = !cfg.longitude_summed;
    const bool surv_t = grids.instantaneous;
    // Reported axes (surviving axes carry their grid; summed axes report empty).
    if (surv_r)  { result.radii = grids.r_grid; }
    if (surv_th) { result.colatitudes = grids.th_grid; }
    if (surv_ph) { result.longitudes = grids.ph_grid; }
    if (surv_t)  { result.times = grids.t_grid; }
    result.all_spatial_summed = cfg.radial_summed && cfg.latitude_summed && cfg.longitude_summed;
    result.n_layers = num_layers;
    result.n_times = surv_t ? grids.t_grid.size() : 1;
    if (surv_r)  { result.shape.push_back(grids.r_grid.size()); }
    if (surv_th) { result.shape.push_back(grids.th_grid.size()); }
    if (surv_ph) { result.shape.push_back(grids.ph_grid.size()); }
    if (surv_t)  { result.shape.push_back(grids.t_grid.size()); }
    return result;
}

// Number of output values a layout holds (the product of its shape; 1 for a fully summed secular total).
inline size_t c_collapse_size_3d(const c_Heating3DCollapsed& layout) {
    size_t size = 1;
    for (const size_t axis_length : layout.shape) { size *= axis_length; }
    return size;
}

typedef std::map<std::array<int, 3>, std::array<double, 36>> c_GramCache3D;

// Analytic colatitude integral of the longitude-mean secular density at one radius: the sum over (|omega|, mu)
// groups and ordered wave pairs (a, b) within a group of (|omega|/2) int Im(c_a conj(c_b) sigma~_a : conj(eps~_b))
// sin(theta) dtheta, the angular integral coming from the (cross-)degree Gram matrices (angular_collapse_.hpp).
// Multiplying by 2*pi*r^2 gives the radial power density dP/dr.
inline double c_secular_theta_integral_3d(
        const c_WaveSet3D& set,
        const std::vector<tides::c_StrainRadialCoeffs>& coeffs_by_group,
        c_GramCache3D& gram_cache) {
    double total = 0.0;
    const size_t num_waves = set.waves.size();
    std::vector<unsigned char> done(num_waves, 0);
    std::vector<size_t> members;
    for (size_t a0 = 0; a0 < num_waves; ++a0) {
        if (done[a0]) {
            continue;
        }

        const int frequency_group = set.wave_frequency_group[a0];
        const int mu0 = c_wave_mu(set.waves[a0]);
        members.clear();
        for (size_t w = a0; w < num_waves; ++w) {
            if (done[w] || set.wave_frequency_group[w] != frequency_group || c_wave_mu(set.waves[w]) != mu0) {
                continue;
            }
            done[w] = 1;
            members.push_back(w);
        }
        const double frequency = set.frequencies[frequency_group];
        for (const size_t a : members) {
            const tides::c_StrainRadialCoeffs& radial_a = coeffs_by_group[set.wave_radial_group[a]];
            if (!radial_a.valid) { continue; }
            const c_TidalWave3D& wave_a = set.waves[a];
            for (const size_t b : members) {
                const tides::c_StrainRadialCoeffs& radial_b = coeffs_by_group[set.wave_radial_group[b]];
                if (!radial_b.valid) { continue; }
                const c_TidalWave3D& wave_b = set.waves[b];
                const std::array<int, 3> key{wave_a.degree_l, wave_b.degree_l, wave_a.order_m};
                auto found = gram_cache.find(key);
                if (found == gram_cache.end()) {
                    double gram[6][6];
                    if (!tides::c_angular_gram_pair(wave_a.degree_l, wave_b.degree_l, wave_a.order_m, gram)) {
                        throw std::runtime_error(
                            "TidalPy: angular Gram matrix out of range during the analytic 3D heating collapse");
                    }
                    std::array<double, 36> flat;
                    for (int i = 0; i < 6; ++i) {
                        for (int j = 0; j < 6; ++j) { flat[static_cast<size_t>(i * 6 + j)] = gram[i][j]; }
                    }
                    found = gram_cache.emplace(key, flat).first;
                }
                double gram[6][6];
                for (int i = 0; i < 6; ++i) {
                    for (int j = 0; j < 6; ++j) { gram[i][j] = found->second[static_cast<size_t>(i * 6 + j)]; }
                }
                const std::complex<double> c_pair = wave_a.amplitude * std::conj(wave_b.amplitude);
                total += 0.5 * frequency * tides::c_theta_integrated_heating_pair(
                    radial_a,
                    radial_b,
                    wave_a.order_m,
                    wave_a.azimuthal_sign,
                    c_pair,
                    gram);
            }
        }
    }
    return total;
}

}  // namespace tides3d

// The 3D orchestration lives on the rheology tide model (the only TideBase with a depth-resolved solution) and
// calls the world's members directly, with no callbacks. It is defined here, in the world extension, where
// c_LayeredWorld and the kernel and potential headers are complete and CyRK lives, so every radial solve and dense
// call stays in its owning extension.

// Scalar form: the batch path with one point.
inline double c_RheologyTide::calc_3d_tidal_heating(
        c_LayeredWorld& world,
        const c_TideSolveConfig& state,
        double radius,
        double colatitude) const {
    double heating = 0.0;
    this->calc_3d_tidal_heating_batch(world, state, &radius, &colatitude, 1, &heating, 1);
    return heating;
}

// World delegation: validate preconditions, map the solve state into the potential's state struct, and
// hand off to the rheology tide model's 3D orchestration.
inline double c_LayeredWorld::get_3d_tidal_heating(
        const c_TideSolveConfig& state,
        double radius,
        double colatitude) {
    if (!this->p_tide) {
        throw std::runtime_error(
            "TidalPy: no tide model attached to the world — call set_tide_model() first");
    }
    auto* rheology = dynamic_cast<c_RheologyTide*>(this->p_tide.get());
    if (rheology == nullptr) {
        throw std::runtime_error(
            "TidalPy: 3D tidal heating requires the rheology tide model (the analytic cpl/ctl/ctl_q "
            "models have no depth-resolved radial solution)");
    }
    if (!this->p_eos_solved || !this->p_eos_solution) {
        throw std::runtime_error(
            "TidalPy: 3D tidal heating needs the EOS solved first — call solve_eos()");
    }
    return rheology->calc_3d_tidal_heating(*this, state, radius, colatitude);
}

// World delegation for the batch path: same preconditions as the scalar get_3d_tidal_heating.
inline void c_LayeredWorld::get_3d_tidal_heating_array(
        const c_TideSolveConfig& state,
        const double* radii,
        const double* colatitudes,
        size_t num_points,
        double* out_heating,
        int num_threads) {
    if (!this->p_tide) {
        throw std::runtime_error(
            "TidalPy: no tide model attached to the world — call set_tide_model() first");
    }
    auto* rheology = dynamic_cast<c_RheologyTide*>(this->p_tide.get());
    if (rheology == nullptr) {
        throw std::runtime_error(
            "TidalPy: 3D tidal heating requires the rheology tide model (the analytic cpl/ctl/ctl_q "
            "models have no depth-resolved radial solution)");
    }
    if (!this->p_eos_solved || !this->p_eos_solution) {
        throw std::runtime_error(
            "TidalPy: 3D tidal heating needs the EOS solved first — call solve_eos()");
    }
    rheology->calc_3d_tidal_heating_batch(
        *this,
        state,
        radii,
        colatitudes,
        num_points,
        out_heating,
        num_threads);
}

// World delegation for the displacement grid: same preconditions as the 3D heating paths.
inline void c_LayeredWorld::get_3d_displacements_grid(
        const c_TideSolveConfig& state,
        const c_Grid3DAxes& axes,
        double* out_disp,
        int num_threads) {
    if (!this->p_tide) {
        throw std::runtime_error(
            "TidalPy: no tide model attached to the world — call set_tide_model() first");
    }
    auto* rheology = dynamic_cast<c_RheologyTide*>(this->p_tide.get());
    if (rheology == nullptr) {
        throw std::runtime_error(
            "TidalPy: 3D tidal displacements require the rheology tide model (the analytic cpl/ctl/ctl_q "
            "models have no depth-resolved radial solution)");
    }
    if (!this->p_eos_solved || !this->p_eos_solution) {
        throw std::runtime_error(
            "TidalPy: 3D tidal displacements need the EOS solved first — call solve_eos()");
    }
    rheology->calc_3d_displacements_grid(
        *this,
        state,
        axes,
        out_disp,
        num_threads);
}

// Instantaneous displacement grid. The coherent wave list is built once, and the radial solve and the y1/y3 samples at
// every radius once per radial group (l, |omega|). At each (r, theta, phi) every wave's complex displacement amplitude
// is added into the total of its frequency, and each component at time t is the sum over frequencies of
// Re[amplitude e^{i |omega| t}], with the phase factors tabulated once. The colatitude rows run on up to num_threads
// threads, each writing only its own cells.
inline void c_RheologyTide::calc_3d_displacements_grid(
        c_LayeredWorld& world,
        const c_TideSolveConfig& state,
        const c_Grid3DAxes& axes,
        double* out_disp,
        int num_threads) const {
    const tides3d::c_WaveSet3D set = tides3d::c_world_wave_set_3d(world, state, "3D tidal displacements");
    const size_t nr = axes.num_radii;
    const size_t nth = axes.num_colatitudes;
    const size_t nph = axes.num_longitudes;
    const size_t nt = axes.num_times;
    if (nr == 0 || nth == 0 || nph == 0 || nt == 0) {
        return;
    }
    const size_t num_groups = set.radial_groups.size();
    const size_t num_waves = set.waves.size();
    const size_t num_frequencies = set.frequencies.size();

    // y1 and y3 of every radial group at every radius, [radius * num_groups + group]. A group without a depth-resolved
    // solution at a radius (the solver's start radius grows with degree l) contributes nothing there, and a radius
    // where no group has one is NaN.
    std::vector<std::complex<double>> y1_at(nr * num_groups);
    std::vector<std::complex<double>> y3_at(nr * num_groups);
    std::vector<unsigned char> group_missing(nr * num_groups, 0);
    std::vector<size_t> radius_missing(nr, 0);
    c_LoveSolveConfig love_cfg = world.make_radial_love_solve_config();
    for (size_t g = 0; g < num_groups; ++g) {
        const ::c_RadialSolutionStorage* storage =
            tides3d::c_solve_radial_group_3d(world, love_cfg, set.radial_groups[g], "3D tidal displacements");
        for (size_t ir = 0; ir < nr; ++ir) {
            std::complex<double> y_at_r[C_MAX_NUM_Y];
            const size_t slot = ir * num_groups + g;
            if (!storage->get_radial_solution(axes.radii[ir], 0, y_at_r)
                || !std::isfinite(y_at_r[0].real()) || !std::isfinite(y_at_r[2].real())) {
                group_missing[slot] = 1;
                radius_missing[ir] += 1;
                continue;
            }
            y1_at[slot] = y_at_r[0];
            y3_at[slot] = y_at_r[2];
        }
    }
    const tides3d::c_PhaseTable3D phase = tides3d::c_phase_table_3d(set.frequencies, axes.times, nt);
    const double nan_v = TidalPyConstants::d_NAN;

    tides3d::c_parallel_tasks_3d(nth, num_threads, [&](size_t ith) {
        const double colatitude = axes.colatitudes[ith];
        tides3d::c_WaveAngular3D angular(set);
        tides3d::c_wave_angular_colatitude_3d(set, colatitude, angular);
        std::vector<tides::c_Vector3> amplitudes(num_frequencies);
        std::vector<unsigned char> active(num_frequencies, 0);
        for (size_t iph = 0; iph < nph; ++iph) {
            tides3d::c_wave_angular_longitude_3d(set, axes.longitudes[iph], false, angular);
            for (size_t ir = 0; ir < nr; ++ir) {
                double* out = out_disp + 3 * (((ir * nth + ith) * nph + iph) * nt);
                if (num_groups > 0 && radius_missing[ir] == num_groups) {
                    std::fill(out, out + 3 * nt, nan_v);
                    continue;
                }
                std::fill(active.begin(), active.end(), 0);
                for (size_t w = 0; w < num_waves; ++w) {
                    const size_t slot = ir * num_groups + static_cast<size_t>(set.wave_radial_group[w]);
                    if (group_missing[slot]) { continue; }
                    tides::c_Vector3 amplitude;
                    tides::c_compute_displacements(
                        y1_at[slot],
                        y3_at[slot],
                        angular.potentials[w],
                        colatitude,
                        amplitude);
                    const size_t f = static_cast<size_t>(set.wave_frequency_group[w]);
                    if (!active[f]) {
                        amplitudes[f] = amplitude;
                        active[f] = 1;
                    } else {
                        for (size_t k = 0; k < 3; ++k) { amplitudes[f].c[k] += amplitude.c[k]; }
                    }
                }
                for (size_t it = 0; it < nt; ++it) {
                    for (size_t k = 0; k < 3; ++k) {
                        double value = 0.0;
                        for (size_t f = 0; f < num_frequencies; ++f) {
                            if (!active[f]) { continue; }
                            const std::complex<double>& a = amplitudes[f].c[k];
                            value += a.real() * phase.cos_phase[f * nt + it] - a.imag() * phase.sin_phase[f * nt + it];
                        }
                        out[3 * it + k] = value;
                    }
                }
            }
        }
    });
}

// Instantaneous stress and strain grid. The coherent wave list, the radial solves with their strain radial
// coefficients, and the phase tables are built once; at each (r, theta, phi) every wave's complex amplitude is added
// into its frequency's total, and each component at time t is the sum over frequencies of
// Re[amplitude e^{i |omega| t}]. The colatitude rows run on up to num_threads threads, each writing only its own cells.
inline void c_RheologyTide::calc_3d_stress_strain_grid(
        c_LayeredWorld& world,
        const c_TideSolveConfig& state,
        const c_Grid3DAxes& axes,
        double* out_stress,
        double* out_strain,
        int num_threads) const {
    if (out_stress == nullptr && out_strain == nullptr) {
        throw std::invalid_argument("TidalPy: the 3D stress and strain grid needs at least one output buffer");
    }
    const size_t nr = axes.num_radii;
    const size_t nth = axes.num_colatitudes;
    const size_t nph = axes.num_longitudes;
    const size_t nt = axes.num_times;
    if (nr == 0 || nth == 0 || nph == 0 || nt == 0) {
        return;
    }

    const char* what = "3D tidal stress and strain";
    const tides3d::c_WaveSet3D set = tides3d::c_world_wave_set_3d(world, state, what);
    const tides3d::c_RadialCoefficients3D radial_coefficients =
        tides3d::c_radial_coefficients_3d(world, set, axes.radii, nr, what);
    const tides3d::c_PhaseTable3D phase = tides3d::c_phase_table_3d(set.frequencies, axes.times, nt);
    const size_t num_frequencies = set.frequencies.size();
    const double nan_v = TidalPyConstants::d_NAN;

    tides3d::c_parallel_tasks_3d(nth, num_threads, [&](size_t ith) {
        tides3d::c_WaveAngular3D angular(set);
        tides3d::c_wave_angular_colatitude_3d(set, axes.colatitudes[ith], angular);
        tides3d::c_FrequencyAmplitudes3D amplitudes(num_frequencies);
        for (size_t iph = 0; iph < nph; ++iph) {
            tides3d::c_wave_angular_longitude_3d(set, axes.longitudes[iph], true, angular);
            for (size_t ir = 0; ir < nr; ++ir) {
                // Undefined at a radius without a solution, and where waves exist but none has a shear kernel; with
                // no active waves at all the tensors are zero.
                const bool defined = !radial_coefficients.radius_failed[ir]
                    && (tides3d::c_frequency_amplitudes_3d(
                            set,
                            radial_coefficients.by_radius[ir],
                            angular,
                            amplitudes)
                        || set.waves.empty());
                const size_t point_offset = ((ir * nth + ith) * nph + iph) * nt;
                for (size_t it = 0; it < nt; ++it) {
                    const size_t offset = 6 * (point_offset + it);
                    for (size_t k = 0; k < 6; ++k) {
                        double stress = nan_v;
                        double strain = nan_v;
                        if (defined) {
                            stress = 0.0;
                            strain = 0.0;
                            for (size_t f = 0; f < num_frequencies; ++f) {
                                if (!amplitudes.active[f]) { continue; }
                                const double cos_wt = phase.cos_phase[f * nt + it];
                                const double sin_wt = phase.sin_phase[f * nt + it];
                                const std::complex<double>& sc = amplitudes.stress[f].c[k];
                                const std::complex<double>& ec = amplitudes.strain[f].c[k];
                                stress += sc.real() * cos_wt - sc.imag() * sin_wt;
                                strain += ec.real() * cos_wt - ec.imag() * sin_wt;
                            }
                        }
                        if (out_stress != nullptr) { out_stress[offset + k] = stress; }
                        if (out_strain != nullptr) { out_strain[offset + k] = strain; }
                    }
                }
            }
        }
    });
}

// World delegation for the stress and strain grid: same preconditions as the 3D heating paths.
inline void c_LayeredWorld::get_3d_stress_strain_grid(
        const c_TideSolveConfig& state,
        const c_Grid3DAxes& axes,
        double* out_stress,
        double* out_strain,
        int num_threads) {
    if (!this->p_tide) {
        throw std::runtime_error("TidalPy: no tide model attached to the world: call set_tide_model() first");
    }
    auto* rheology = dynamic_cast<c_RheologyTide*>(this->p_tide.get());
    if (rheology == nullptr) {
        throw std::runtime_error(
            "TidalPy: 3D tidal stress and strain require the rheology tide model (the analytic cpl/ctl/ctl_q "
            "models have no depth-resolved radial solution)");
    }
    if (!this->p_eos_solved || !this->p_eos_solution) {
        throw std::runtime_error("TidalPy: 3D tidal stress and strain need the EOS solved first: call solve_eos()");
    }
    rheology->calc_3d_stress_strain_grid(*this, state, axes, out_stress, out_strain, num_threads);
}

// Batch form of the secular 3D heating: the longitude-mean secular density at num_points paired (radius,
// colatitude) points. The coherent wave list is built once, the radial solve runs once per radial group
// (l, |omega|), and its strain radial coefficients are evaluated once per unique radius (points on a map share
// radii). Points that share a colatitude share its angular work, and the colatitudes run on up to num_threads
// threads. A point whose radius has no depth-resolved solution (the center, below the solver start) is NaN.
inline void c_RheologyTide::calc_3d_tidal_heating_batch(
        c_LayeredWorld& world,
        const c_TideSolveConfig& state,
        const double* radii,
        const double* colatitudes,
        size_t num_points,
        double* out_heating,
        int num_threads) const {
    const tides3d::c_WaveSet3D set = tides3d::c_world_wave_set_3d(world, state, "secular 3D tidal heating");
    for (size_t i = 0; i < num_points; ++i) {
        out_heating[i] = 0.0;
    }
    if (num_points == 0) {
        return;
    }

    // Unique radii.
    std::vector<double> unique_radii;
    std::vector<size_t> point_radius(num_points, 0);
    std::vector<unsigned char> point_invalid(num_points, 0);
    std::map<double, size_t> radius_index;
    for (size_t i = 0; i < num_points; ++i) {
        if (!std::isfinite(radii[i])) {
            point_invalid[i] = 1;
            continue;
        }
        auto found = radius_index.find(radii[i]);
        if (found == radius_index.end()) {
            found = radius_index.emplace(radii[i], unique_radii.size()).first;
            unique_radii.push_back(radii[i]);
        }
        point_radius[i] = found->second;
    }
    // Strain radial coefficients [unique radius][radial group]
    const tides3d::c_RadialCoefficients3D radial_coefficients = tides3d::c_radial_coefficients_3d(
        world,
        set,
        unique_radii.data(),
        unique_radii.size(),
        "secular 3D tidal heating");

    // One row per distinct finite colatitude, and one row for each point whose colatitude is not finite.
    std::vector<std::vector<size_t>> rows;
    std::map<double, size_t> colatitude_row;
    for (size_t i = 0; i < num_points; ++i) {
        if (!std::isfinite(colatitudes[i])) {
            rows.push_back(std::vector<size_t>{i});
            continue;
        }
        auto found = colatitude_row.find(colatitudes[i]);
        if (found == colatitude_row.end()) {
            found = colatitude_row.emplace(colatitudes[i], rows.size()).first;
            rows.emplace_back();
        }
        rows[found->second].push_back(i);
    }
    const std::vector<tides3d::c_SecularGroup3D> groups = tides3d::c_secular_groups_3d(set, true);

    tides3d::c_parallel_tasks_3d(rows.size(), num_threads, [&](size_t row) {
        const std::vector<size_t>& points = rows[row];
        tides3d::c_WaveAngular3D angular(set);
        tides3d::c_wave_angular_colatitude_3d(set, colatitudes[points[0]], angular);
        tides3d::c_wave_angular_longitude_3d(set, 0.0, true, angular);
        for (const size_t i : points) {
            if (point_invalid[i] || radial_coefficients.radius_failed[point_radius[i]]) {
                out_heating[i] = TidalPyConstants::d_NAN;
                continue;
            }
            out_heating[i] = tides3d::c_secular_density_3d(
                set,
                groups,
                radial_coefficients.by_radius[point_radius[i]],
                angular);
        }
    });
}

// =====================================================================================================================
// Collapsed (summed / averaged) 3D tidal heating
// =====================================================================================================================

// Produce the 3D tidal heating as a full grid over (radius, colatitude, longitude[, time]) or integrated along any
// spatial dimension, written into caller buffers sized by c_LayeredWorld::calc_3d_tides_layout.
// orbit_averaged=true gives the secular density h_bar (the pointwise time average, or its longitude mean when
// longitude is summed); orbit_averaged=false gives the instantaneous power sigma_ij(t) eps_dot_ij(t) at each user
// time. See c_Heating3DCollapseConfig and c_Heating3DCollapsed for the flags and output conventions. The radial
// solves run on the calling thread and the per-point evaluation on up to cfg.num_threads threads over colatitude
// rows.
inline void c_RheologyTide::calc_3d_tidal_heating_collapsed(
        c_LayeredWorld& world,
        const c_TideSolveConfig& state,
        const double* radii,
        size_t num_radii,
        const double* colatitudes,
        size_t num_colatitudes,
        const double* longitudes,
        size_t num_longitudes,
        const double* times,
        size_t num_times,
        const c_Heating3DCollapseConfig& cfg,
        double* out_values,
        double* out_layer_totals) const {
    // Coherent wave list (built once, reused across every grid point).
    const tides3d::c_WaveSet3D set = tides3d::c_world_wave_set_3d(world, state, "3D tidal heating");

    // Axis grids, output layout, and zeroed outputs.
    const tides3d::c_CollapseGrids3D grids = tides3d::c_collapse_grids_3d(
        world,
        radii,
        num_radii,
        colatitudes,
        num_colatitudes,
        longitudes,
        num_longitudes,
        times,
        num_times,
        cfg);
    const size_t num_layers = world.get_num_layers();
    const c_Heating3DCollapsed layout = tides3d::c_collapse_layout_3d(grids, cfg, num_layers);
    const size_t num_values = tides3d::c_collapse_size_3d(layout);
    const bool totals = layout.all_spatial_summed;
    const size_t num_layer_totals = totals ? num_layers * layout.n_times : 0;
    if (num_values > 0) { std::fill(out_values, out_values + num_values, 0.0); }
    if (num_layer_totals > 0) { std::fill(out_layer_totals, out_layer_totals + num_layer_totals, 0.0); }

    const size_t nr = grids.r_grid.size();
    const size_t nth = grids.th_grid.size();
    const size_t nph = grids.ph_grid.size();
    const size_t nt = grids.t_grid.size();
    if (nr == 0 || nth == 0 || nph == 0 || nt == 0) {
        return;
    }
    const bool instantaneous = grids.instantaneous;
    const bool any_summed = grids.any_summed;
    const bool surv_r = !cfg.radial_summed;
    const bool surv_th = !cfg.latitude_summed;
    const bool surv_ph = !cfg.longitude_summed;
    const bool surv_t = instantaneous;

    // Row-major flat index over the surviving axes, in the fixed order [radius, colatitude, longitude, time].
    auto surv_index = [&](size_t ir, size_t ith, size_t iph, size_t it) -> size_t {
        size_t idx = 0;
        if (surv_r)  { idx = idx * nr + ir; }
        if (surv_th) { idx = idx * nth + ith; }
        if (surv_ph) { idx = idx * nph + iph; }
        if (surv_t)  { idx = idx * nt + it; }
        return idx;
    };
    // Combined per-point weight: summed axes -> integration weight (with Jacobian); surviving spatial axes
    // -> their Jacobian when any axis is summed, else 1 (raw density). Longitude/time Jacobian is 1.
    auto combined_weight = [&](size_t ir, size_t ith, size_t iph) -> double {
        const double wr = cfg.radial_summed ? grids.r_wsum[ir]
                                            : (any_summed ? grids.r_grid[ir] * grids.r_grid[ir] : 1.0);
        const double wth = cfg.latitude_summed ? grids.th_wsum[ith]
                                               : (any_summed ? std::sin(grids.th_grid[ith]) : 1.0);
        const double wph = cfg.longitude_summed ? grids.ph_wsum[iph] : 1.0;
        return wr * wth * wph;
    };

    // Radial solve once per radial group (l, |omega|) and strain radial coefficients [radius][group]
    const tides3d::c_RadialCoefficients3D radial_coefficients = tides3d::c_radial_coefficients_3d(
        world,
        set,
        grids.r_grid.data(),
        nr,
        "3D tidal heating");
    const std::vector<std::vector<tides::c_StrainRadialCoeffs>>& coeffs = radial_coefficients.by_radius;
    const std::vector<unsigned char>& radius_solve_failed = radial_coefficients.radius_failed;
    const double nan_v = TidalPyConstants::d_NAN;

    if (!instantaneous && cfg.latitude_summed && cfg.latitude_analytic && grids.latitude_full_sphere) {
        // Analytic colatitude collapse: integrate the longitude-mean secular density over theta with the
        // (cross-)degree Gram matrices (exact, no theta grid). theta is summed away, so scatter over (radius, phi).
        // There is no per-point grid to spread over threads, so this runs on the calling thread.
        tides3d::c_GramCache3D gram_cache;
        for (size_t ir = 0; ir < nr; ++ir) {
            double theta_integral = 0.0;
            if (!radius_solve_failed[ir]) {
                theta_integral = tides3d::c_secular_theta_integral_3d(set, coeffs[ir], gram_cache);
            }
            // theta already integrated (Gram absorbs sin theta); apply only the radial + longitude factors.
            const double radial_factor = cfg.radial_summed ? grids.r_wsum[ir] : (grids.r_grid[ir] * grids.r_grid[ir]);
            for (size_t iph = 0; iph < nph; ++iph) {
                const double longitude_factor = cfg.longitude_summed ? grids.ph_wsum[iph] : 1.0;
                const double contrib = theta_integral * radial_factor * longitude_factor;
                out_values[surv_index(ir, 0, iph, 0)] += contrib;
                if (totals) {
                    out_layer_totals[grids.r_layer[ir] * nt + 0] += contrib;
                }
            }
        }
        return;
    }

    // Evaluate over colatitude rows. A row owns its cells when colatitude survives and writes them directly. When
    // colatitude is summed every row adds into the same cells, so each row fills a buffer of its own and the rows are
    // merged in row order below; either way the result is identical for any thread count.
    const bool rows_share_cells = !surv_th;
    std::vector<std::vector<double>> row_values(rows_share_cells ? nth : 0);
    std::vector<std::vector<double>> row_layer_totals(totals ? nth : 0);
    const size_t num_frequencies = set.frequencies.size();
    // Secular: each frequency's waves (split by mu for the longitude mean, which the single phi node then carries
    // exactly) summed coherently, then (|omega|/2) Im(sigma_c : conj(eps_c)). Instantaneous: every wave's complex
    // stress and strain amplitude added into the total of its frequency, each frequency evolved in time as
    // Re[. e^{i |omega| t}] with the phase factors tabulated once, and the real fields summed.
    const bool longitude_averaged = cfg.longitude_summed;
    const std::vector<tides3d::c_SecularGroup3D> groups = instantaneous
        ? std::vector<tides3d::c_SecularGroup3D>()
        : tides3d::c_secular_groups_3d(set, longitude_averaged);
    const tides3d::c_PhaseTable3D phase =
        tides3d::c_phase_table_3d(set.frequencies, grids.t_grid.data(), instantaneous ? nt : 0);

    tides3d::c_parallel_tasks_3d(nth, cfg.num_threads, [&](size_t ith) {
        double* values = out_values;
        if (rows_share_cells) {
            row_values[ith].assign(num_values, 0.0);
            values = row_values[ith].data();
        }
        double* layer_totals = out_layer_totals;
        if (totals) {
            row_layer_totals[ith].assign(num_layer_totals, 0.0);
            layer_totals = row_layer_totals[ith].data();
        }
        tides3d::c_WaveAngular3D angular(set);
        tides3d::c_wave_angular_colatitude_3d(set, grids.th_grid[ith], angular);

        if (!instantaneous) {
            for (size_t iph = 0; iph < nph; ++iph) {
                tides3d::c_wave_angular_longitude_3d(
                    set,
                    longitude_averaged ? 0.0 : grids.ph_grid[iph],
                    true,
                    angular);
                for (size_t ir = 0; ir < nr; ++ir) {
                    const size_t index = surv_index(ir, ith, iph, 0);
                    if (radius_solve_failed[ir]) {
                        if (!any_summed) { values[index] = nan_v; }
                        continue;
                    }
                    const double density = tides3d::c_secular_density_3d(set, groups, coeffs[ir], angular);
                    if (!any_summed) {
                        values[index] = density;
                    } else {
                        const double contrib = density * combined_weight(ir, ith, iph);
                        values[index] += contrib;
                        if (totals) {
                            layer_totals[grids.r_layer[ir] * nt + 0] += contrib;
                        }
                    }
                }
            }
            return;
        }

        tides3d::c_FrequencyAmplitudes3D amplitudes(num_frequencies);
        for (size_t iph = 0; iph < nph; ++iph) {
            tides3d::c_wave_angular_longitude_3d(set, grids.ph_grid[iph], true, angular);
            for (size_t ir = 0; ir < nr; ++ir) {
                if (radius_solve_failed[ir]) {
                    if (!any_summed) {
                        for (size_t it = 0; it < nt; ++it) {
                            values[surv_index(ir, ith, iph, it)] = nan_v;
                        }
                    }
                    continue;
                }
                tides3d::c_frequency_amplitudes_3d(
                    set,
                    coeffs[ir],
                    angular,
                    amplitudes);
                const double point_weight = any_summed ? combined_weight(ir, ith, iph) : 1.0;
                for (size_t it = 0; it < nt; ++it) {
                    double power = 0.0;
                    for (size_t k = 0; k < 6; ++k) {
                        double sigma = 0.0;
                        double eps_dot = 0.0;
                        for (size_t f = 0; f < num_frequencies; ++f) {
                            if (!amplitudes.active[f]) { continue; }
                            const double omega = set.frequencies[f];
                            const double cos_wt = phase.cos_phase[f * nt + it];
                            const double sin_wt = phase.sin_phase[f * nt + it];
                            const std::complex<double>& sc = amplitudes.stress[f].c[k];
                            const std::complex<double>& ec = amplitudes.strain[f].c[k];
                            sigma   += sc.real() * cos_wt - sc.imag() * sin_wt;
                            eps_dot += -omega * (ec.real() * sin_wt + ec.imag() * cos_wt);
                        }
                        power += ((k < 3) ? 1.0 : 2.0) * sigma * eps_dot;
                    }
                    if (!any_summed) {
                        values[surv_index(ir, ith, iph, it)] = power;
                    } else {
                        const double contrib = power * point_weight;
                        values[surv_index(ir, ith, iph, it)] += contrib;
                        if (totals) {
                            layer_totals[grids.r_layer[ir] * nt + it] += contrib;
                        }
                    }
                }
            }
        }
    });

    // Merge the rows that share cells, in row order.
    for (size_t ith = 0; ith < nth; ++ith) {
        if (rows_share_cells) {
            for (size_t k = 0; k < num_values; ++k) { out_values[k] += row_values[ith][k]; }
        }
        if (totals) {
            for (size_t k = 0; k < num_layer_totals; ++k) { out_layer_totals[k] += row_layer_totals[ith][k]; }
        }
    }
}

// World delegation for the collapse layout: the axes and output shape calc_3d_tides produces, from the layer geometry
// alone.
inline c_Heating3DCollapsed c_LayeredWorld::calc_3d_tides_layout(
        const double* radii,
        size_t num_radii,
        const double* colatitudes,
        size_t num_colatitudes,
        const double* longitudes,
        size_t num_longitudes,
        const double* times,
        size_t num_times,
        const c_Heating3DCollapseConfig& cfg) {
    return tides3d::c_collapse_layout_3d(
        tides3d::c_collapse_grids_3d(
            *this,
            radii,
            num_radii,
            colatitudes,
            num_colatitudes,
            longitudes,
            num_longitudes,
            times,
            num_times,
            cfg),
        cfg,
        this->get_num_layers());
}

// World delegation for the collapse path into caller buffers: same preconditions as the scalar get_3d_tidal_heating.
inline void c_LayeredWorld::calc_3d_tides_into(
        const c_TideSolveConfig& state,
        const double* radii,
        size_t num_radii,
        const double* colatitudes,
        size_t num_colatitudes,
        const double* longitudes,
        size_t num_longitudes,
        const double* times,
        size_t num_times,
        const c_Heating3DCollapseConfig& cfg,
        double* out_values,
        double* out_layer_totals) {
    if (!this->p_tide) {
        throw std::runtime_error(
            "TidalPy: no tide model attached to the world — call set_tide_model() first");
    }
    auto* rheology = dynamic_cast<c_RheologyTide*>(this->p_tide.get());
    if (rheology == nullptr) {
        throw std::runtime_error(
            "TidalPy: 3D tidal heating requires the rheology tide model (the analytic cpl/ctl/ctl_q "
            "models have no depth-resolved radial solution)");
    }
    if (!this->p_eos_solved || !this->p_eos_solution) {
        throw std::runtime_error(
            "TidalPy: 3D tidal heating needs the EOS solved first — call solve_eos()");
    }
    rheology->calc_3d_tidal_heating_collapsed(
        *this,
        state,
        radii,
        num_radii,
        colatitudes,
        num_colatitudes,
        longitudes,
        num_longitudes,
        times,
        num_times,
        cfg,
        out_values,
        out_layer_totals);
}

// World delegation for the collapse path returning vectors: the layout, then calc_3d_tides_into.
inline c_Heating3DCollapsed c_LayeredWorld::calc_3d_tides(
        const c_TideSolveConfig& state,
        const double* radii,
        size_t num_radii,
        const double* colatitudes,
        size_t num_colatitudes,
        const double* longitudes,
        size_t num_longitudes,
        const double* times,
        size_t num_times,
        const c_Heating3DCollapseConfig& cfg) {
    c_Heating3DCollapsed result = this->calc_3d_tides_layout(
        radii,
        num_radii,
        colatitudes,
        num_colatitudes,
        longitudes,
        num_longitudes,
        times,
        num_times,
        cfg);
    result.values.assign(tides3d::c_collapse_size_3d(result), 0.0);
    if (result.all_spatial_summed) {
        result.layer_totals.assign(result.n_layers * result.n_times, 0.0);
    }
    this->calc_3d_tides_into(
        state,
        radii,
        num_radii,
        colatitudes,
        num_colatitudes,
        longitudes,
        num_longitudes,
        times,
        num_times,
        cfg,
        result.values.data(),
        result.all_spatial_summed ? result.layer_totals.data() : nullptr);
    return result;
}

} // namespace tidalpy
