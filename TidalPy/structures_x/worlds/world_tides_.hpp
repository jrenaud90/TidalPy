#pragma once
/*
 * world_tides_.hpp — out-of-line definition of c_LayeredWorld::calc_tides.
 *
 * c_LayeredWorld extends the common analytic tide path (c_BaseWorld::calc_tides,
 * world_tides_base_.hpp) with two layered-world-only capabilities:
 *   - the rheology model: -Im[k_l(omega)] from the world radial solver run at each unique
 *     tidal frequency (the EOS must be solved first); and
 *   - per-layer heating distribution by each layer's tidal_scale_method.
 *
 * It runs the global-potential engine, collapses (analytic or rheology), stores the result,
 * then distributes the heat to the layers. The tide-model holder / config / result state is
 * inherited from c_BaseWorld.
 *
 * This header pulls in the heavy global-potential tables; force-include it in the layered
 * (and gas-giant) world extension only.
 */

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
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
                this->solve_love_numbers(love_cfg);
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

// Effective per-layer tidal-heating scale for the layer's tidal_scale_method (0 for a
// non-tidal layer).
//   user_provided   : the layer's tidal_scale field.
//   volume_fraction : layer volume / planet volume.
//   tidal_timescale : a log-Gaussian bell curve in the layer's Maxwell time tau = eta/mu
//                     (from its static shear modulus + viscosity) about the tidal forcing
//                     period 2*pi/|orbital_frequency|; width [decades] from the tide config.
//                     Returns 0 for a geometry-only layer or when mu/eta/forcing are unusable.
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
// Secular (cycle/orbit-averaged) heating density:
//     h_bar(r, theta, phi) = sum over |omega| of (|omega|/2) Im( sigma_c : conj(eps_c) )
// with sigma_c, eps_c the total complex stress/strain amplitude at that frequency, i.e. every wave at that |omega|
// summed before the bilinear form. Cross terms between different frequencies average to zero over the orbit and are
// dropped; cross terms between waves at the same frequency survive the average and are kept. They are what the
// m = 0 pairs contribute (each pair is one real sinusoid), and what makes the heating of a synchronously rotating
// body, whose active modes all sit at multiples of n, depend on longitude. The scalar / batch paths take no
// longitude and return the longitudinal mean of h_bar: cross terms between waves with different azimuthal structure
// e^{i mu phi} integrate to zero over phi, so the mean is the sum over (|omega|, mu) groups evaluated at phi = 0.
// The volume integral of h_bar is the 1D global heating (get_tidal_heating). At nonzero obliquity the 3D value is
// for the geometry with zero argument of periapse and node (the engine drops precession): same-frequency modes of
// the same (l, m) then combine coherently, a cross term the precession-averaged 1D formula does not carry, so the
// two agree only to the size of those terms there.

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
};

inline c_WaveSet3D c_build_wave_set_3d(
        const std::vector<c_TidalPotential3DModeCoeff>& modes,
        double min_frequency) {
    c_WaveSet3D set;
    set.waves = c_coherent_tidal_waves_3d(modes, min_frequency);
    const size_t num_waves = set.waves.size();
    set.wave_radial_group.assign(num_waves, -1);
    set.wave_frequency_group.assign(num_waves, -1);
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
        world.get_radius(), state.semi_major_axis, state.orbital_frequency, state.spin_frequency,
        state.obliquity, state.eccentricity, state.host_mass, c_get_G(),
        tide_cfg.min_degree_l, tide_cfg.max_degree_l,
        tide_cfg.obliquity_truncation, tide_cfg.eccentricity_truncation, &engine_error);
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
        y_at_r[0], y_at_r[1], y_at_r[2], y_at_r[3], shear, bulk, radius,
        static_cast<double>(group.degree_l), is_solid, is_incompressible);
    return true;
}

// Signed azimuthal wavenumber of a wave (its longitude structure is e^{i mu phi}).
inline int c_wave_mu(const c_TidalWave3D& wave) {
    return wave.azimuthal_sign * wave.order_m;
}

// Secular density at one point. coeffs_by_group[g] is radial group g's strain radial coefficient set at this
// radius (valid == false where the group has no shear kernel there). longitude_averaged -> the longitude mean
// (waves split by mu, evaluated at phi = 0); otherwise the pointwise value at `longitude`.
inline double c_secular_density_3d(
        const c_WaveSet3D& set,
        const std::vector<tides::c_StrainRadialCoeffs>& coeffs_by_group,
        double colatitude,
        double longitude,
        bool longitude_averaged) {
    double heating = 0.0;
    const size_t num_waves = set.waves.size();
    std::vector<unsigned char> done(num_waves, 0);
    for (size_t w0 = 0; w0 < num_waves; ++w0) {
        if (done[w0]) { continue; }
        const int frequency_group = set.wave_frequency_group[w0];
        const int mu0 = c_wave_mu(set.waves[w0]);
        tides::c_Tensor6 stress_total;
        tides::c_Tensor6 strain_total;
        bool any = false;
        for (size_t w = w0; w < num_waves; ++w) {
            if (done[w] || set.wave_frequency_group[w] != frequency_group) { continue; }
            if (longitude_averaged && c_wave_mu(set.waves[w]) != mu0) { continue; }
            done[w] = 1;
            const tides::c_StrainRadialCoeffs& radial = coeffs_by_group[set.wave_radial_group[w]];
            if (!radial.valid) { continue; }
            const c_PotentialPointC potential =
                c_eval_wave_point_3d(set.waves[w], colatitude, longitude_averaged ? 0.0 : longitude);
            tides::c_Tensor6 strain;
            tides::c_Tensor6 stress;
            tides::c_compute_strain_stress(radial, potential, colatitude, strain, stress);
            for (size_t k = 0; k < 6; ++k) {
                stress_total.c[k] += stress.c[k];
                strain_total.c[k] += strain.c[k];
            }
            any = true;
        }
        if (any) {
            heating += 0.5 * set.frequencies[frequency_group]
                     * tides::c_volumetric_heating_signed(stress_total, strain_total);
        }
    }
    return heating;
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
                    radial_a, radial_b, wave_a.order_m, wave_a.azimuthal_sign, c_pair, gram);
            }
        }
    }
    return total;
}

}  // namespace tides3d

// The 3D orchestration lives on the rheology tide model (the only TideBase with a depth-resolved solution); it
// calls the world's members directly (no callbacks). Defined here, in the world extension, where c_LayeredWorld
// + the kernel/potential headers are complete and CyRK lives, so every radial solve + dense call stays in its
// owning extension.

// Scalar form: the batch path with one point.
inline double c_RheologyTide::calc_3d_tidal_heating(
        c_LayeredWorld& world,
        const c_TideSolveConfig& state,
        double radius,
        double colatitude) const {
    double heating = 0.0;
    this->calc_3d_tidal_heating_batch(world, state, &radius, &colatitude, 1, &heating);
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
        double* out_heating) {
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
    rheology->calc_3d_tidal_heating_batch(*this, state, radii, colatitudes, num_points, out_heating);
}

// World delegation for the displacement grid: same preconditions as the 3D heating paths.
inline void c_LayeredWorld::get_3d_displacements_grid(
        const c_TideSolveConfig& state,
        const double* radii,
        size_t num_radii,
        const double* colatitudes,
        size_t num_colatitudes,
        const double* longitudes,
        size_t num_longitudes,
        const double* times,
        size_t num_times,
        double* out_disp) {
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
        *this, state, radii, num_radii, colatitudes, num_colatitudes, longitudes, num_longitudes,
        times, num_times, out_disp);
}

// Instantaneous displacement grid. The coherent wave list is built once; the radial solve and the y1/y3 samples
// at every radius are computed once per radial group (l, |omega|); each wave's complex displacement amplitude at
// (r, theta, phi) is then evolved over the time grid as Re[u_c e^{i |omega| t}] and the waves are summed.
inline void c_RheologyTide::calc_3d_displacements_grid(
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
        double* out_disp) const {
    const tides3d::c_WaveSet3D set = tides3d::c_world_wave_set_3d(world, state, "3D tidal displacements");
    const size_t nr = num_radii, nth = num_colatitudes, nph = num_longitudes, nt = num_times;
    const size_t total = nr * nth * nph * nt;
    for (size_t i = 0; i < 3 * total; ++i) {
        out_disp[i] = 0.0;
    }
    // A radius is unusable only if NO radial group has a depth-resolved solution there (the solver's start
    // radius grows with degree l; a group without a solution at a radius contributes nothing there).
    const size_t num_groups = set.radial_groups.size();
    std::vector<size_t> radius_missing(nr, 0);
    std::vector<unsigned char> group_missing(nr, 0);
    std::vector<std::complex<double>> y1_at_r(nr), y3_at_r(nr);
    c_LoveSolveConfig love_cfg = world.make_radial_love_solve_config();
    for (size_t g = 0; g < num_groups; ++g) {
        const ::c_RadialSolutionStorage* storage =
            tides3d::c_solve_radial_group_3d(world, love_cfg, set.radial_groups[g], "3D tidal displacements");
        for (size_t ir = 0; ir < nr; ++ir) {
            std::complex<double> y_at_r[C_MAX_NUM_Y];
            if (!storage->get_radial_solution(radii[ir], 0, y_at_r)
                || !std::isfinite(y_at_r[0].real()) || !std::isfinite(y_at_r[2].real())) {
                group_missing[ir] = 1;
                radius_missing[ir] += 1;
                continue;
            }
            group_missing[ir] = 0;
            y1_at_r[ir] = y_at_r[0];
            y3_at_r[ir] = y_at_r[2];
        }
        for (size_t w = 0; w < set.waves.size(); ++w) {
            if (set.wave_radial_group[w] != static_cast<int>(g)) { continue; }
            const c_TidalWave3D& wave = set.waves[w];
            for (size_t ith = 0; ith < nth; ++ith) {
                for (size_t iph = 0; iph < nph; ++iph) {
                    const c_PotentialPointC potential =
                        c_eval_wave_point_3d(wave, colatitudes[ith], longitudes[iph]);
                    for (size_t ir = 0; ir < nr; ++ir) {
                        if (group_missing[ir]) {
                            continue;
                        }
                        tides::c_Vector3 amplitude;
                        tides::c_compute_displacements(
                            y1_at_r[ir], y3_at_r[ir], potential, colatitudes[ith], amplitude);
                        for (size_t it = 0; it < nt; ++it) {
                            const double phase = wave.frequency * times[it];
                            const double cos_wt = std::cos(phase);
                            const double sin_wt = std::sin(phase);
                            double* out = out_disp + 3 * (((ir * nth + ith) * nph + iph) * nt + it);
                            for (size_t k = 0; k < 3; ++k) {
                                out[k] += amplitude.c[k].real() * cos_wt - amplitude.c[k].imag() * sin_wt;
                            }
                        }
                    }
                }
            }
        }
    }
    for (size_t ir = 0; ir < nr; ++ir) {
        if (!(num_groups > 0 && radius_missing[ir] == num_groups)) {
            continue;
        }
        double* out = out_disp + 3 * (ir * nth * nph * nt);
        for (size_t i = 0; i < 3 * nth * nph * nt; ++i) {
            out[i] = TidalPyConstants::d_NAN;
        }
    }
}

// Batch form of the secular 3D heating: the longitude-mean secular density at num_points paired (radius,
// colatitude) points. The coherent wave list is built once, the radial solve runs once per radial group
// (l, |omega|), and its strain radial coefficients are evaluated once per UNIQUE radius (points on a map share
// radii), then every point sums its frequency groups coherently. A point whose radius has no depth-resolved
// solution (below the solver start / center) is NaN.
inline void c_RheologyTide::calc_3d_tidal_heating_batch(
        c_LayeredWorld& world,
        const c_TideSolveConfig& state,
        const double* radii,
        const double* colatitudes,
        size_t num_points,
        double* out_heating) const {
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
    const size_t num_radii = unique_radii.size();
    const size_t num_groups = set.radial_groups.size();

    // Strain radial coefficients [unique radius][radial group]. A radius is unusable only if NO group has a
    // depth-resolved solution there: the solver's start radius grows with degree l, so a higher-degree group
    // whose solution starts further out simply contributes nothing below it.
    std::vector<std::vector<tides::c_StrainRadialCoeffs>> coeffs(
        num_radii, std::vector<tides::c_StrainRadialCoeffs>(num_groups));
    std::vector<size_t> radius_missing(num_radii, 0);
    c_LoveSolveConfig love_cfg = world.make_radial_love_solve_config();
    for (size_t g = 0; g < num_groups; ++g) {
        const ::c_RadialSolutionStorage* storage =
            tides3d::c_solve_radial_group_3d(world, love_cfg, set.radial_groups[g], "secular 3D tidal heating");
        for (size_t ur = 0; ur < num_radii; ++ur) {
            if (!tides3d::c_strain_coeffs_at_radius_3d(
                    world, storage, unique_radii[ur], set.radial_groups[g], coeffs[ur][g])) {
                radius_missing[ur] += 1;
            }
        }
    }

    for (size_t i = 0; i < num_points; ++i) {
        const bool unusable = point_invalid[i]
            || (num_groups > 0 && radius_missing[point_radius[i]] == num_groups);
        if (unusable) {
            out_heating[i] = TidalPyConstants::d_NAN;
            continue;
        }
        out_heating[i] = tides3d::c_secular_density_3d(set, coeffs[point_radius[i]], colatitudes[i], 0.0, true);
    }
}

// =====================================================================================================================
// Collapsed (summed / averaged) 3D tidal heating
// =====================================================================================================================

// Produce the 3D tidal heating as a full grid over (radius, colatitude, longitude[, time]) or reduced (integrated)
// along any spatial dimension. orbit_averaged=true gives the secular density h_bar(r, theta, phi) (the pointwise
// time average; the longitude mean when longitude is summed); orbit_averaged=false gives the instantaneous power
// sigma_ij(t) eps_dot_ij(t) at each user time. See c_Heating3DCollapseConfig / c_Heating3DCollapsed for the
// flags + output conventions.
inline c_Heating3DCollapsed c_RheologyTide::calc_3d_tidal_heating_collapsed(
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
        const c_Heating3DCollapseConfig& cfg) const {
    c_Heating3DCollapsed result;
    const double two_pi = 2.0 * TidalPyConstants::d_PI;
    const bool instantaneous = !cfg.orbit_averaged;
    const bool any_summed = cfg.latitude_summed || cfg.longitude_summed || cfg.radial_summed;

    // Coherent wave list (built once, reused across every grid point).
    const tides3d::c_WaveSet3D set = tides3d::c_world_wave_set_3d(world, state, "3D tidal heating");
    const size_t num_waves = set.waves.size();
    const size_t num_groups = set.radial_groups.size();

    // Axis Grids
    // Radius: user array unless summed (then a per-layer trapezoid grid; r_wsum carries the r^2 Jacobian).
    std::vector<double> r_grid, r_wsum;
    std::vector<size_t> r_layer;
    const size_t num_layers = world.get_num_layers();
    if (cfg.radial_summed) {
        const int slices = (cfg.radial_slices > 1) ? cfg.radial_slices : 80;
        for (size_t layer_i = 0; layer_i < num_layers; ++layer_i) {
            const c_BaseLayer* layer = world.get_layer(layer_i);
            const double r_inner = layer->get_radius_inner();
            const double r_outer = layer->get_radius_outer();
            const double dr = (r_outer - r_inner) / static_cast<double>(slices);
            for (int s = 0; s <= slices; ++s) {
                const double rr = r_inner + dr * static_cast<double>(s);
                r_grid.push_back(rr);
                r_wsum.push_back(((s == 0 || s == slices) ? 0.5 * dr : dr) * rr * rr);  // trapezoid x r^2
                r_layer.push_back(layer_i);
            }
        }
    } else {
        r_grid.assign(radii, radii + num_radii);
    }
    // Colatitude: user array unless summed (Gauss-Legendre in cos theta; the weight absorbs sin theta).
    // A colatitude band [min, max] narrower than [0, pi] maps the nodes onto [cos(max), cos(min)]
    // (affine substitution; the weights scale by the half-width) so the integral covers only the band.
    const bool latitude_full_sphere =
        (cfg.colatitude_min <= TidalPyConstants::d_EPS)
        && (cfg.colatitude_max >= TidalPyConstants::d_PI - TidalPyConstants::d_EPS);
    std::vector<double> th_grid, th_wsum;
    if (cfg.latitude_summed) {
        if (cfg.colatitude_min < 0.0 || cfg.colatitude_max > TidalPyConstants::d_PI
                || cfg.colatitude_min >= cfg.colatitude_max) {
            throw std::invalid_argument(
                "TidalPy: colatitude band must satisfy 0 <= colatitude_min < colatitude_max <= pi");
        }
        const int num_nodes = (cfg.latitude_nodes > 1) ? cfg.latitude_nodes : 48;
        std::vector<double> gl_x;
        tides::c_gauss_legendre_nodes(num_nodes, gl_x, th_wsum);
        const double x_lo = std::cos(cfg.colatitude_max);
        const double x_hi = std::cos(cfg.colatitude_min);
        const double x_mid  = 0.5 * (x_hi + x_lo);
        const double x_half = 0.5 * (x_hi - x_lo);
        th_grid.resize(num_nodes);
        for (int i = 0; i < num_nodes; ++i) {
            th_grid[i] = std::acos(x_mid + x_half * gl_x[i]);
            th_wsum[i] *= x_half;
        }
    } else {
        th_grid.assign(colatitudes, colatitudes + num_colatitudes);
    }
    // Longitude: user array unless summed. Secular -> analytic 2*pi on the longitude-mean density (single
    // point). Instantaneous -> periodic trapezoid over [0, 2*pi) (uniform weight, no endpoint halving since the
    // field is periodic).
    std::vector<double> ph_grid, ph_wsum;
    if (cfg.longitude_summed) {
        if (instantaneous) {
            const int num_nodes = (cfg.longitude_nodes > 1) ? cfg.longitude_nodes : 64;
            const double dph = two_pi / static_cast<double>(num_nodes);
            for (int k = 0; k < num_nodes; ++k) { ph_grid.push_back(k * dph); ph_wsum.push_back(dph); }
        } else {
            ph_grid.push_back(0.0);
            ph_wsum.push_back(two_pi);
        }
    } else {
        ph_grid.assign(longitudes, longitudes + num_longitudes);
    }
    // Time: user array (instantaneous only); a single dummy sample when orbit-averaged.
    std::vector<double> t_grid;
    if (instantaneous) { t_grid.assign(times, times + num_times); }
    else { t_grid.push_back(0.0); }

    const size_t nr = r_grid.size();
    const size_t nth = th_grid.size();
    const size_t nph = ph_grid.size();
    const size_t nt = t_grid.size();

    const bool surv_r = !cfg.radial_summed;
    const bool surv_th = !cfg.latitude_summed;
    const bool surv_ph = !cfg.longitude_summed;
    const bool surv_t = instantaneous;

    // Reported axes (surviving axes carry their grid; summed axes report empty).
    if (surv_r)  { result.radii = r_grid; }
    if (surv_th) { result.colatitudes = th_grid; }
    if (surv_ph) { result.longitudes = ph_grid; }
    if (surv_t)  { result.times = t_grid; }
    result.all_spatial_summed = cfg.radial_summed && cfg.latitude_summed && cfg.longitude_summed;
    result.n_layers = num_layers;
    result.n_times = surv_t ? nt : 1;
    if (surv_r)  { result.shape.push_back(nr); }
    if (surv_th) { result.shape.push_back(nth); }
    if (surv_ph) { result.shape.push_back(nph); }
    if (surv_t)  { result.shape.push_back(nt); }
    if (nr == 0 || nth == 0 || nph == 0 || nt == 0) {
        return result;
    }

    size_t out_size = 1;
    for (size_t s : result.shape) { out_size *= s; }
    result.values.assign(out_size, 0.0);
    if (result.all_spatial_summed) {
        result.layer_totals.assign(num_layers * nt, 0.0);
    }

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
        const double wr = cfg.radial_summed ? r_wsum[ir]
                                            : (any_summed ? r_grid[ir] * r_grid[ir] : 1.0);
        const double wth = cfg.latitude_summed ? th_wsum[ith]
                                               : (any_summed ? std::sin(th_grid[ith]) : 1.0);
        const double wph = cfg.longitude_summed ? ph_wsum[iph] : 1.0;
        return wr * wth * wph;
    };

    // Radial solve once per radial group (l, |omega|); strain radial coefficients [radius][group]. A radius has
    // no depth-resolved solution (radius_solve_failed) only if NO group has one there: the solver's start radius
    // grows with degree l, so a higher-degree group whose solution starts further out simply contributes
    // nothing below it.
    std::vector<std::vector<tides::c_StrainRadialCoeffs>> coeffs(
        nr, std::vector<tides::c_StrainRadialCoeffs>(num_groups));
    std::vector<size_t> radius_missing(nr, 0);
    c_LoveSolveConfig love_cfg = world.make_radial_love_solve_config();
    for (size_t g = 0; g < num_groups; ++g) {
        const ::c_RadialSolutionStorage* storage =
            tides3d::c_solve_radial_group_3d(world, love_cfg, set.radial_groups[g], "3D tidal heating");
        for (size_t ir = 0; ir < nr; ++ir) {
            if (!tides3d::c_strain_coeffs_at_radius_3d(world, storage, r_grid[ir], set.radial_groups[g], coeffs[ir][g])) {
                radius_missing[ir] += 1;
            }
        }
    }
    std::vector<unsigned char> radius_solve_failed(nr, 0);
    for (size_t ir = 0; ir < nr; ++ir) {
        radius_solve_failed[ir] = (num_groups > 0 && radius_missing[ir] == num_groups) ? 1 : 0;
    }

    // Evaluate and Reduce
    const double nan_v = TidalPyConstants::d_NAN;
    if (!instantaneous && cfg.latitude_summed && cfg.latitude_analytic && latitude_full_sphere) {
        // Analytic colatitude collapse: integrate the longitude-mean secular density over theta with the
        // (cross-)degree Gram matrices (exact, no theta grid). theta is summed away, so scatter over (radius, phi).
        tides3d::c_GramCache3D gram_cache;
        for (size_t ir = 0; ir < nr; ++ir) {
            double theta_integral = 0.0;
            if (!radius_solve_failed[ir]) {
                theta_integral = tides3d::c_secular_theta_integral_3d(set, coeffs[ir], gram_cache);
            }
            // theta already integrated (Gram absorbs sin theta); apply only the radial + longitude factors.
            const double radial_factor = cfg.radial_summed ? r_wsum[ir] : (r_grid[ir] * r_grid[ir]);
            for (size_t iph = 0; iph < nph; ++iph) {
                const double longitude_factor = cfg.longitude_summed ? ph_wsum[iph] : 1.0;
                const double contrib = theta_integral * radial_factor * longitude_factor;
                result.values[surv_index(ir, 0, iph, 0)] += contrib;
                if (result.all_spatial_summed) {
                    result.layer_totals[r_layer[ir] * nt + 0] += contrib;
                }
            }
        }
    } else if (!instantaneous) {
        // Secular density h_bar(r, theta, phi): each frequency's waves summed coherently, then
        // (|omega|/2) Im(sigma_c : conj(eps_c)). With longitude summed the single phi node carries the
        // longitude MEAN (exact); otherwise the pointwise value at each user longitude.
        const bool longitude_averaged = cfg.longitude_summed;
        for (size_t ir = 0; ir < nr; ++ir) {
            for (size_t ith = 0; ith < nth; ++ith) {
                for (size_t iph = 0; iph < nph; ++iph) {
                    double density = 0.0;
                    if (!radius_solve_failed[ir]) {
                        density = tides3d::c_secular_density_3d(
                            set, coeffs[ir], th_grid[ith], ph_grid[iph], longitude_averaged);
                    }
                    if (!any_summed) {
                        result.values[surv_index(ir, ith, iph, 0)] =
                            radius_solve_failed[ir] ? nan_v : density;
                    } else if (!radius_solve_failed[ir]) {
                        const double contrib = density * combined_weight(ir, ith, iph);
                        result.values[surv_index(ir, ith, iph, 0)] += contrib;
                        if (result.all_spatial_summed) {
                            result.layer_totals[r_layer[ir] * nt + 0] += contrib;
                        }
                    }
                }
            }
        }
    } else {
        // Instantaneous power sigma_ij(t) eps_dot_ij(t). Build each wave's complex stress/strain amplitude at
        // (r, theta, phi), then evolve in time as Re[. e^{i |omega| t}] and sum the real fields (every cross
        // term present).
        std::vector<tides::c_Tensor6> wave_stress, wave_strain;
        std::vector<double> wave_freq;
        for (size_t ir = 0; ir < nr; ++ir) {
            for (size_t ith = 0; ith < nth; ++ith) {
                for (size_t iph = 0; iph < nph; ++iph) {
                    if (radius_solve_failed[ir]) {
                        if (!any_summed) {
                            for (size_t it = 0; it < nt; ++it) {
                                result.values[surv_index(ir, ith, iph, it)] = nan_v;
                            }
                        }
                        continue;
                    }
                    wave_stress.clear();
                    wave_strain.clear();
                    wave_freq.clear();
                    for (size_t w = 0; w < num_waves; ++w) {
                        const tides::c_StrainRadialCoeffs& radial = coeffs[ir][set.wave_radial_group[w]];
                        if (!radial.valid) { continue; }
                        const c_PotentialPointC potential =
                            c_eval_wave_point_3d(set.waves[w], th_grid[ith], ph_grid[iph]);
                        tides::c_Tensor6 strain, stress;
                        tides::c_compute_strain_stress(radial, potential, th_grid[ith], strain, stress);
                        wave_stress.push_back(stress);
                        wave_strain.push_back(strain);
                        wave_freq.push_back(set.waves[w].frequency);
                    }
                    for (size_t it = 0; it < nt; ++it) {
                        const double time = t_grid[it];
                        double power = 0.0;
                        for (size_t k = 0; k < 6; ++k) {
                            double sigma = 0.0;
                            double eps_dot = 0.0;
                            for (size_t ww = 0; ww < wave_freq.size(); ++ww) {
                                const double omega = wave_freq[ww];
                                const double cos_wt = std::cos(omega * time);
                                const double sin_wt = std::sin(omega * time);
                                const std::complex<double>& sc = wave_stress[ww].c[k];
                                const std::complex<double>& ec = wave_strain[ww].c[k];
                                sigma   += sc.real() * cos_wt - sc.imag() * sin_wt;
                                eps_dot += -omega * (ec.real() * sin_wt + ec.imag() * cos_wt);
                            }
                            power += ((k < 3) ? 1.0 : 2.0) * sigma * eps_dot;
                        }
                        if (!any_summed) {
                            result.values[surv_index(ir, ith, iph, it)] = power;
                        } else {
                            const double contrib = power * combined_weight(ir, ith, iph);
                            result.values[surv_index(ir, ith, iph, it)] += contrib;
                            if (result.all_spatial_summed) {
                                result.layer_totals[r_layer[ir] * nt + it] += contrib;
                            }
                        }
                    }
                }
            }
        }
    }
    return result;
}

// World delegation for the collapse path: same preconditions as the scalar get_3d_tidal_heating.
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
    return rheology->calc_3d_tidal_heating_collapsed(
        *this, state, radii, num_radii, colatitudes, num_colatitudes,
        longitudes, num_longitudes, times, num_times, cfg);
}

} // namespace tidalpy
