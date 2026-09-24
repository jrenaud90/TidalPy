#pragma once
/*
 * world_tides_.hpp: out-of-line definition of c_LayeredWorld::calc_tides and the 3D tidal paths.
 *
 * c_LayeredWorld extends the common analytic tide path (c_BaseWorld::calc_tides in world_tides_base_.hpp) with
 * two layered-world capabilities: the rheology model, whose -Im[k_l(omega)] comes from the world's Love solve at
 * each unique tidal frequency (the EOS must be solved first), and each layer's share of the heating. The Love
 * solves of these paths go into workspaces of their own, so they leave the world's last solve_love_numbers result
 * alone. The tide-model holder, config, and result state stay on c_BaseWorld.
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
#include "world_tides_base_.hpp"                         // c_world_global_potential
#include "../../Tides_x/classes/tide_collapse_.hpp"      // c_global_potential, c_collapse_global_tides
#include "../../Tides_x/classes/tide_.hpp"               // c_RheologyTide (3D orchestration target)
#include "../../Tides_x/potential/potential_3d_.hpp"     // c_tidal_potential_3d_modes (dynamic Kaula engine)
#include "../../Tides_x/multilayer/kernel_.hpp"          // strain/stress/heating kernel + c_StrainRadialCoeffs
#include "../../Tides_x/multilayer/angular_collapse_.hpp" // c_theta_integrated_heating_pair (analytic collapse)

namespace tidalpy {

// The heating [W] each layer takes depends on where the Love numbers come from:
//   radial_solver, propagation_matrix : the volume integral of the radial solution's orbit-averaged heating density
//                                       over the layer (calc_layer_tidal_heating_radial), when the tides config's
//                                       layer_tidal_heating is on; NaN otherwise.
//   homogeneous, cpl, ctl             : each tidal layer is a homogeneous planet of its own averaged material whose
//                                       Im(k) is scaled by its tidal scale; the layer takes the heating of those
//                                       scaled Love numbers, and the layers sum to the total.
//   an analytic tide model            : the whole-body heating times the layer's tidal scale.
inline void c_LayeredWorld::calc_tides(const c_TideSolveConfig& state) {
    const c_WorldCallLock call_lock(this->p_call_mutex.get());
    if (!this->p_tide) {
        throw std::runtime_error(
            "TidalPy: no tide model attached to the world: call set_tide_model() first");
    }
    this->p_check_tide_state(state);

    const double planet_radius = this->get_radius();
    const double planet_volume =
        (4.0 / 3.0) * TidalPyConstants::d_PI * planet_radius * planet_radius * planet_radius;
    const c_TideConfig& tcfg = this->p_tide_config;
    const std::size_t n_layers = this->p_layers.size();

    // Model-independent per-mode terms, plus the unique-frequency maps.
    const c_GlobalPotentialStorage potential =
        c_world_global_potential(*this, state, this->p_tides_solved, this->p_tide_result);

    // Everything is gathered here and committed to the world at the end.
    c_GlobalTideResult tide_result;
    c_IntMap<c_Key4, tidalpy::c_LoveNumbers> tide_love;
    std::vector<double> layer_heating(n_layers, TidalPyConstants::d_NAN);
    if (this->p_tide->needs_radial_solve()) {
        // The per-mode -Im[k_l(omega)] comes from the world's Love solve, which needs a solved EOS.
        if (!this->p_eos_solved || !this->p_eos_solution) {
            this->p_tides_solved = false;
            throw std::runtime_error(
                "TidalPy: the rheology tide model needs the EOS solved first. Call "
                "solve_eos() before calc_tides()");
        }

        // Once per unique (degree_l, frequency) pair, so modes sharing a degree and frequency reuse one solve; each
        // active mode's Love numbers are then recorded by its (l, m, p, q). The quasi-homogeneous methods also keep
        // each tidal layer's scaled Love numbers per solve, in the fixed layer order of the averages cache. When the
        // per-layer heating integral will follow, each solve gets a workspace of its own and is kept for it: the
        // integral needs the radial solution of exactly these (degree, |omega|) groups.
        c_LoveWorkspace workspace;
        c_LoveSolveConfig love_cfg = this->make_love_solve_config();
        const bool quasi_homogeneous = c_love_method_is_homogeneous(c_love_method_from_int(love_cfg.love_method));
        const bool retain_radial_solves = !quasi_homogeneous && tcfg.layer_tidal_heating;
        std::vector<std::unique_ptr<c_LoveWorkspace>> retained_workspaces;
        std::vector<c_RetainedRadialSolve> retained_solves;
        c_HomogeneousLoveCache homogeneous_cache;
        c_IntMap<c_Key2, std::size_t> solve_by_l_freq;
        std::vector<tidalpy::c_LoveNumbers> world_love_by_solve;
        std::vector<std::vector<tidalpy::c_LoveNumbers>> layer_love_by_solve;
        std::vector<std::size_t> part_layer_index;
        c_IntMap<c_Key4, std::size_t> solve_by_mode;
        for (const auto& mode_entry : potential.potential_map) {
            const c_Key4& lmpq_key = mode_entry.first;
            const int degree_l     = static_cast<int>(lmpq_key.a);
            // Every mode of the potential map has a nonzero frequency.
            const std::size_t freq_index = mode_entry.second.frequency_index;
            const double frequency = potential.unique_freq_map[freq_index].frequency;

            c_Key2 lf_key(static_cast<int16_t>(degree_l), static_cast<int16_t>(freq_index));
            bool cached = false;
            std::size_t solve_index = solve_by_l_freq.get(cached, lf_key);
            if (!cached) {
                love_cfg.degree_l  = degree_l;
                love_cfg.frequency = frequency;
                c_LoveWorkspace* solve_workspace = &workspace;
                if (retain_radial_solves) {
                    retained_workspaces.push_back(std::make_unique<c_LoveWorkspace>());
                    solve_workspace = retained_workspaces.back().get();
                }
                this->solve_love_numbers(love_cfg, &homogeneous_cache, *solve_workspace);
                if (!solve_workspace->get_success()) {
                    this->p_tides_solved = false;
                    throw std::runtime_error(
                        "TidalPy: Love-number solve failed during calc_tides: " + solve_workspace->get_message());
                }
                solve_index = world_love_by_solve.size();
                world_love_by_solve.push_back(solve_workspace->get_love(0));
                if (retain_radial_solves) {
                    retained_solves.push_back(
                        c_RetainedRadialSolve{degree_l, frequency, solve_workspace->get_storage()});
                }
                if (quasi_homogeneous) {
                    std::vector<tidalpy::c_LoveNumbers> scaled_parts;
                    scaled_parts.reserve(workspace.analytic_layers.size());
                    part_layer_index.clear();
                    for (const c_LayerLove& part : workspace.analytic_layers) {
                        scaled_parts.emplace_back(
                            part.tidal_scale * part.love.k,
                            part.tidal_scale * part.love.h,
                            part.tidal_scale * part.love.l);
                        part_layer_index.push_back(part.layer_index);
                    }
                    layer_love_by_solve.push_back(std::move(scaled_parts));
                }
                solve_by_l_freq.set(lf_key, solve_index);
            }
            tide_love.set(lmpq_key, world_love_by_solve[solve_index]);
            solve_by_mode.set(lmpq_key, solve_index);
        }

        tide_result = c_collapse_global_tides(potential, *this->p_tide, &tide_love);

        if (quasi_homogeneous) {
            // The collapse is linear in each mode's -Im[k], so the layers' own collapses sum to the total.
            std::fill(layer_heating.begin(), layer_heating.end(), 0.0);
            for (std::size_t part_i = 0; part_i < part_layer_index.size(); ++part_i) {
                c_IntMap<c_Key4, tidalpy::c_LoveNumbers> part_love;
                for (const auto& mode_entry : solve_by_mode.data) {
                    part_love.set(mode_entry.first, layer_love_by_solve[mode_entry.second][part_i]);
                }
                layer_heating[part_layer_index[part_i]] =
                    c_collapse_global_tides(potential, *this->p_tide, &part_love).tidal_heating;
            }
        } else if (tcfg.layer_tidal_heating) {
            this->calc_layer_tidal_heating_radial(state, tide_result.tidal_heating, layer_heating, &retained_solves);
        }
    } else {
        // The analytic models need no Love solve. They fix the whole body's heating, so each tidal layer takes the
        // fraction its tidal scale is of all the tidal layers' scales, and the layers sum to the total. With no
        // tidal layer the heating has nowhere to go and every layer reports NaN.
        tide_result = c_collapse_global_tides(potential, *this->p_tide, nullptr);
        double scale_sum = 0.0;
        for (std::size_t i = 0; i < n_layers; ++i) {
            layer_heating[i] = this->p_layers[i]->calc_tidal_scale(planet_volume);
            scale_sum += layer_heating[i];
        }
        for (std::size_t i = 0; i < n_layers; ++i) {
            layer_heating[i] = (scale_sum > 0.0)
                ? tide_result.tidal_heating * layer_heating[i] / scale_sum
                : TidalPyConstants::d_NAN;
        }
    }

    // Commit: layer.get_tidal_heating() reports each layer's share.
    this->p_tide_result        = tide_result;
    this->p_tide_solver_love   = tide_love;
    this->p_layer_tidal_heating = layer_heating;
    for (std::size_t i = 0; i < n_layers; ++i) {
        this->p_layers[i]->set_tidal_heating(layer_heating[i]);
    }
    this->p_tides_solved = true;
}

// Each layer's orbit-averaged heating [W] from the radial solution: the secular heating density integrated over
// the layer's volume, with the analytic colatitude integral, the 2 pi longitude integral, and Gauss-Legendre nodes
// inside each layer (the same integral as calc_3d_tides with every axis summed). The layers are scaled so they sum
// to `total_heating`, the 1D global result, which removes the small radial-quadrature residual between the two. A
// liquid layer carries no shear dissipation and takes 0; with no usable integral every layer is NaN.
// `retained_solves`, when given, are the radial solves the 1D pass already ran; the integral reuses them instead of
// solving its (degree, |omega|) groups again, which would double the cost of calc_tides.
inline void c_LayeredWorld::calc_layer_tidal_heating_radial(
        const c_TideSolveConfig& state,
        double total_heating,
        std::vector<double>& out,
        const std::vector<c_RetainedRadialSolve>* retained_solves) {
    const std::size_t n_layers = this->p_layers.size();
    out.assign(n_layers, TidalPyConstants::d_NAN);
    c_Heating3DCollapseConfig cfg;
    cfg.orbit_averaged   = true;
    cfg.latitude_summed  = true;
    cfg.longitude_summed = true;
    cfg.radial_summed    = true;
    // Lent for this call only; the call lock is held throughout, so no other call sees them.
    struct c_RetainedSolvesLoan {
        const std::vector<c_RetainedRadialSolve>*& slot;
        c_RetainedSolvesLoan(const std::vector<c_RetainedRadialSolve>*& slot_in,
                             const std::vector<c_RetainedRadialSolve>* loan) : slot(slot_in) { slot = loan; }
        ~c_RetainedSolvesLoan() { slot = nullptr; }
    } loan(this->p_retained_radial_solves, retained_solves);
    const c_Heating3DCollapsed integrated = this->calc_3d_tides(
        state,
        nullptr,
        0,
        nullptr,
        0,
        nullptr,
        0,
        nullptr,
        0,
        cfg);
    if (integrated.layer_totals.size() < n_layers) { return; }
    double integral_sum = 0.0;
    for (std::size_t i = 0; i < n_layers; ++i) {
        const double layer_total = integrated.layer_totals[i];
        if (std::isfinite(layer_total)) { integral_sum += layer_total; }
    }
    if (!std::isfinite(total_heating)) { return; }
    if (!(std::abs(integral_sum) > 0.0)) {
        // Nothing dissipates (a body with no viscosity, say): every layer takes the zero total.
        if (total_heating == 0.0) { std::fill(out.begin(), out.end(), 0.0); }
        return;
    }
    for (std::size_t i = 0; i < n_layers; ++i) {
        const double layer_total = integrated.layer_totals[i];
        out[i] = std::isfinite(layer_total) ? total_heating * layer_total / integral_sum : 0.0;
    }
}

namespace tides3d {

struct c_RadialGroup3D {
    int degree_l = 0;
    double frequency = 0.0;   // |omega| [rad s-1]
};

struct c_WaveSet3D {
    std::vector<c_TidalWave3D> waves;
    std::vector<c_RadialGroup3D> radial_groups;   // unique (degree l, |omega|): one radial solve each
    std::vector<int> wave_radial_group;           // per wave, index into radial_groups
    std::vector<double> frequencies;              // unique |omega|: waves sharing one combine coherently
    std::vector<int> wave_frequency_group;        // per wave, index into frequencies
    std::vector<int> wave_frequency_slot;         // per wave, its place among its frequency's waves, in wave order
    std::vector<size_t> frequency_num_waves;      // per frequency, its number of waves
    std::vector<std::array<int, 2>> angular_pairs; // unique (degree l, order m): one Legendre evaluation each
    std::vector<int> wave_angular_pair;           // per wave, index into angular_pairs
    std::vector<int> azimuthal_orders;            // unique signed azimuthal wavenumber mu = azimuthal_sign * m
    std::vector<int> wave_azimuthal_order;        // per wave, index into azimuthal_orders
    // The frequency match tolerance the set was built with, which later matches of its frequencies use too.
    double frequency_match_rtol = 0.0;
    // A secular set's cut amplitude products (c_wave_pair_power_3d) of the ordered pairs of waves of one frequency,
    // which the secular heating takes in place of amplitude_a * conj(amplitude_b): frequency f's block starts at
    // pair_power_offset[f] and holds its waves' pairs row-major by wave_frequency_slot. Empty for an instantaneous
    // set, whose fields are linear in the unsquared amplitudes.
    bool secular = false;
    std::vector<std::complex<double>> pair_power;
    std::vector<size_t> pair_power_offset;

    // The cut product of waves a and b of one frequency.
    const std::complex<double>& pair_power_of(size_t a, size_t b) const noexcept {
        const size_t f = static_cast<size_t>(this->wave_frequency_group[a]);
        return this->pair_power[this->pair_power_offset[f]
            + static_cast<size_t>(this->wave_frequency_slot[a]) * this->frequency_num_waves[f]
            + static_cast<size_t>(this->wave_frequency_slot[b])];
    }
};

// The nonzero cut products of a secular set's waves before the filter: row a holds (b, product) for each wave b >= a
// whose frequency a's matches, sorted by b.
typedef std::vector<std::vector<std::pair<size_t, std::complex<double>>>> c_PairPowerRows3D;

// The product of waves a and b from their rows: row a's entry for a < b, and otherwise the conjugate of row b's entry
// for a (a wave's product with itself is conjugated too); zero for a pair with no entry.
inline std::complex<double> c_pair_power_from_rows_3d(const c_PairPowerRows3D& rows, size_t a, size_t b) {
    const bool direct = a < b;
    const std::vector<std::pair<size_t, std::complex<double>>>& row = rows[direct ? a : b];
    const size_t partner = direct ? b : a;
    const auto found = std::lower_bound(
        row.begin(), row.end(), partner,
        [](const std::pair<size_t, std::complex<double>>& entry, size_t index) { return entry.first < index; });
    if ((found == row.end()) || (found->first != partner)) { return std::complex<double>(0.0, 0.0); }
    return direct ? found->second : std::conj(found->second);
}

// A secular set keeps only the waves with a non-zero cut product with some wave of their frequency (at zero obliquity,
// the modes with |q| <= N / 2), so the radial solves cover exactly what the heating uses.
inline c_WaveSet3D c_build_wave_set_3d(
        const std::vector<c_TidalPotential3DModeCoeff>& modes,
        const c_FrequencyTolerance& tolerance,
        double eccentricity,
        double obliquity,
        bool secular) {
    c_WaveSet3D set;
    set.secular = secular;
    set.frequency_match_rtol = tolerance.match_rtol;
    set.waves = c_coherent_tidal_waves_3d(modes, tolerance);
    c_PairPowerRows3D pair_rows;
    std::vector<size_t> kept;   // per kept wave, its index before the filter
    if (secular) {
        const size_t num_all = set.waves.size();
        c_WaveFrequencyIndex all_frequencies = c_wave_frequency_index(tolerance.match_rtol);
        for (size_t a = 0; a < num_all; ++a) { all_frequencies.insert(0, set.waves[a].frequency); }
        pair_rows.resize(num_all);
        std::vector<unsigned char> keep(num_all, 0);
        std::vector<size_t> partners;
        for (size_t a = 0; a < num_all; ++a) {
            partners.clear();
            all_frequencies.for_each_match(0, set.waves[a].frequency, [&](size_t b) {
                if (b >= a) { partners.push_back(b); }
            });
            std::sort(partners.begin(), partners.end());
            for (const size_t b : partners) {
                const std::complex<double> power =
                    c_wave_pair_power_3d(set.waves[a], set.waves[b], eccentricity, obliquity);
                if (std::abs(power) == 0.0) { continue; }
                pair_rows[a].emplace_back(b, power);
                keep[a] = 1;
                keep[b] = 1;
            }
        }
        std::vector<c_TidalWave3D> kept_waves;
        for (size_t a = 0; a < num_all; ++a) {
            if (!keep[a]) { continue; }
            kept.push_back(a);
            kept_waves.push_back(std::move(set.waves[a]));
        }
        set.waves = std::move(kept_waves);
    }
    const size_t num_waves = set.waves.size();
    set.wave_radial_group.assign(num_waves, -1);
    set.wave_frequency_group.assign(num_waves, -1);
    set.wave_frequency_slot.assign(num_waves, -1);
    set.wave_angular_pair.assign(num_waves, -1);
    set.wave_azimuthal_order.assign(num_waves, -1);
    // Each wave joins the first group whose frequency it matches, as a scan of the groups in order would find.
    c_WaveFrequencyIndex radial_index = c_wave_frequency_index(tolerance.match_rtol);
    c_WaveFrequencyIndex frequency_index = c_wave_frequency_index(tolerance.match_rtol);
    for (size_t w = 0; w < num_waves; ++w) {
        const c_TidalWave3D& wave = set.waves[w];
        std::ptrdiff_t radial_group = radial_index.find(wave.degree_l, wave.frequency);
        if (radial_group < 0) {
            radial_group = static_cast<std::ptrdiff_t>(radial_index.insert(wave.degree_l, wave.frequency));
            set.radial_groups.push_back(c_RadialGroup3D{wave.degree_l, wave.frequency});
        }
        set.wave_radial_group[w] = static_cast<int>(radial_group);

        std::ptrdiff_t frequency_group = frequency_index.find(0, wave.frequency);
        if (frequency_group < 0) {
            frequency_group = static_cast<std::ptrdiff_t>(frequency_index.insert(0, wave.frequency));
            set.frequencies.push_back(wave.frequency);
            set.frequency_num_waves.push_back(0);
        }
        set.wave_frequency_group[w] = static_cast<int>(frequency_group);
        set.wave_frequency_slot[w] = static_cast<int>(set.frequency_num_waves[static_cast<size_t>(frequency_group)]++);

        // Waves of one (l, m) share their Legendre values, and waves of one mu their e^{i mu phi}.
        const std::array<int, 2> pair{wave.degree_l, wave.order_m};
        const auto pair_found = std::find(set.angular_pairs.begin(), set.angular_pairs.end(), pair);
        set.wave_angular_pair[w] = static_cast<int>(pair_found - set.angular_pairs.begin());
        if (pair_found == set.angular_pairs.end()) { set.angular_pairs.push_back(pair); }
        const int mu = wave.azimuthal_sign * wave.order_m;
        const auto mu_found = std::find(set.azimuthal_orders.begin(), set.azimuthal_orders.end(), mu);
        set.wave_azimuthal_order[w] = static_cast<int>(mu_found - set.azimuthal_orders.begin());
        if (mu_found == set.azimuthal_orders.end()) { set.azimuthal_orders.push_back(mu); }
    }
    if (secular) {
        // Only pairs within one frequency are ever read, so each frequency keeps a block of its own.
        const size_t num_frequencies = set.frequencies.size();
        std::vector<std::vector<size_t>> frequency_waves(num_frequencies);
        for (size_t w = 0; w < num_waves; ++w) {
            frequency_waves[static_cast<size_t>(set.wave_frequency_group[w])].push_back(w);
        }
        set.pair_power_offset.assign(num_frequencies, 0);
        size_t num_pairs = 0;
        for (size_t f = 0; f < num_frequencies; ++f) {
            set.pair_power_offset[f] = num_pairs;
            num_pairs += set.frequency_num_waves[f] * set.frequency_num_waves[f];
        }
        set.pair_power.assign(num_pairs, std::complex<double>(0.0, 0.0));
        for (size_t f = 0; f < num_frequencies; ++f) {
            const std::vector<size_t>& members = frequency_waves[f];
            const size_t block = members.size();
            for (size_t i = 0; i < block; ++i) {
                for (size_t j = 0; j < block; ++j) {
                    set.pair_power[set.pair_power_offset[f] + i * block + j] =
                        c_pair_power_from_rows_3d(pair_rows, kept[members[i]], kept[members[j]]);
                }
            }
        }
    }
    return set;
}

// The wave set for the world's tide config and an orbital state.
inline c_WaveSet3D c_world_wave_set_3d(
        c_LayeredWorld& world,
        const c_TideSolveConfig& state,
        const char* what,
        bool secular) {
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
        tide_cfg.eccentricity_exact_tolerance,
        &engine_error);
    if (engine_error != 0) {
        throw std::runtime_error(
            std::string("TidalPy: tidal potential engine failed during ") + what + " (error "
            + std::to_string(engine_error) + "); check degree/truncation levels");
    }
    // The floor below which a mode is inactive is the 1D path's (c_record_unique_frequencies), so both paths keep the
    // same modes; a slow mode near a Maxwell peak otherwise went missing from the 3D heating alone.
    return c_build_wave_set_3d(modes, c_read_frequency_tolerance(), state.eccentricity, state.obliquity, secular);
}

// Run the world radial solve for one radial group into `workspace` and return its solution storage, which lives
// until the workspace's next solve. The world's own last Love solve is left alone.
inline const ::c_RadialSolutionStorage* c_solve_radial_group_3d(
        const c_LayeredWorld& world,
        c_LoveSolveConfig& love_cfg,
        const c_RadialGroup3D& group,
        const char* what,
        c_LoveWorkspace& workspace) {
    love_cfg.degree_l = group.degree_l;
    love_cfg.frequency = group.frequency;
    world.solve_love_numbers(love_cfg, nullptr, workspace);
    if (!workspace.get_success()) {
        throw std::runtime_error(
            std::string("TidalPy: radial solve failed during ") + what + ": " + workspace.get_message());
    }
    const ::c_RadialSolutionStorage* storage = workspace.get_storage();
    if (storage == nullptr) {
        throw std::runtime_error(std::string("TidalPy: missing radial solution during ") + what);
    }
    return storage;
}

// The layer at a radius as the strain there needs it, looked up once per radius.
struct c_RadiusLayer3D {
    const c_PhysicsLayer* physics_layer = nullptr;   // null for a geometry-only layer
    // A liquid layer, or a molten stretch the radial solver treats as a static liquid.
    bool liquid = false;
};

inline c_RadiusLayer3D c_radius_layer_3d(const c_LayeredWorld& world, double radius) {
    c_RadiusLayer3D layer;
    layer.physics_layer = dynamic_cast<const c_PhysicsLayer*>(world.find_layer_for_radius(radius));
    layer.liquid = ((layer.physics_layer != nullptr) && !layer.physics_layer->get_is_solid())
        || world.get_is_molten_at(radius);
    return layer;
}

// Strain radial coefficients of one radial group at one radius. False where there is no depth-resolved strain
// solution (the center, below the solver start), where a point-wise quantity is NaN and a radial sum takes it as
// zero. A liquid point is not missing: it has no shear kernel, so its coefficients are invalid and it contributes no
// heating (0) while its stress and strain are NaN. A geometry-only layer has NaN moduli.
inline bool c_strain_coeffs_at_radius_3d(
        const c_RadiusLayer3D& layer,
        const ::c_RadialSolutionStorage* storage,
        double radius,
        const c_RadialGroup3D& group,
        tides::c_StrainRadialCoeffs& out) {
    if (layer.liquid) {
        out = tides::c_StrainRadialCoeffs();
        out.valid = false;
        return true;
    }
    std::complex<double> y_at_r[C_MAX_NUM_Y];
    if (!storage->get_radial_solution(radius, 0, y_at_r)
        || !std::isfinite(y_at_r[0].real()) || !std::isfinite(y_at_r[1].real())
        || !std::isfinite(y_at_r[2].real()) || !std::isfinite(y_at_r[3].real())) {
        out.valid = false;
        return false;
    }
    bool is_solid = true;
    bool is_incompressible = false;
    std::complex<double> shear(TidalPyConstants::d_NAN, 0.0);
    std::complex<double> bulk(TidalPyConstants::d_NAN, 0.0);
    if (layer.physics_layer != nullptr) {
        is_solid = layer.physics_layer->get_is_solid();
        is_incompressible = layer.physics_layer->get_is_incompressible();
        shear = layer.physics_layer->calc_complex_shear_modulus(radius, group.frequency);
        bulk  = layer.physics_layer->calc_complex_bulk_modulus(radius, group.frequency);
    }
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

// Solve once per radial group (l, |omega|) and evaluate each group's strain radial coefficients at every
// radius. A radius is unusable only when no group has a solution there: the solver's start radius grows
// with l, so a higher-degree group starting further out just contributes nothing below it.
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
    c_LoveWorkspace workspace;
    // A solve calc_tides already ran for a group (same degree and |omega|) is used as it is: the first such solve.
    c_WaveFrequencyIndex retained_index = c_wave_frequency_index(set.frequency_match_rtol);
    std::vector<const ::c_RadialSolutionStorage*> retained_storage;
    const std::vector<c_RetainedRadialSolve>* retained = world.get_retained_radial_solves();
    if (retained != nullptr) {
        for (const c_RetainedRadialSolve& solve : *retained) {
            if (solve.storage == nullptr) { continue; }
            retained_index.insert(solve.degree_l, solve.frequency);
            retained_storage.push_back(solve.storage);
        }
    }
    std::vector<c_RadiusLayer3D> radius_layers(num_radii);
    for (size_t ir = 0; ir < num_radii; ++ir) {
        radius_layers[ir] = c_radius_layer_3d(world, radii[ir]);
    }
    for (size_t g = 0; g < num_groups; ++g) {
        const std::ptrdiff_t retained_solve =
            retained_index.find(set.radial_groups[g].degree_l, set.radial_groups[g].frequency);
        const ::c_RadialSolutionStorage* storage = (retained_solve >= 0)
            ? retained_storage[static_cast<size_t>(retained_solve)]
            : c_solve_radial_group_3d(world, love_cfg, set.radial_groups[g], what, workspace);
        for (size_t ir = 0; ir < num_radii; ++ir) {
            if (!c_strain_coeffs_at_radius_3d(
                radius_layers[ir],
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

// cos(|omega_f| t) and sin(|omega_f| t) per frequency and time, row-major [f * num_times + it]. They depend
// on neither position nor wave, so a grid tabulates them once.
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

// Total complex stress and strain amplitude per frequency at one point; active[f] is 1 where frequency f
// received at least one wave.
struct c_FrequencyAmplitudes3D {
    std::vector<tides::c_Tensor6> stress;
    std::vector<tides::c_Tensor6> strain;
    std::vector<unsigned char> active;

    explicit c_FrequencyAmplitudes3D(size_t num_frequencies) :
        stress(num_frequencies),
        strain(num_frequencies),
        active(num_frequencies, 0) {}
};

// The radius-independent part of every wave at one (colatitude, longitude), which a grid forms once per
// point and reuses at every radius. c_wave_angular_colatitude_3d sets the colatitude: the Legendre values
// of each (l, m) and the sine and cotangent factors. c_wave_angular_longitude_3d then sets the longitude:
// the phasor e^{i mu phi} of each mu, every wave's potential point, and its angular strain factors. For a secular set
// the potentials are per unit amplitude: its heating weighs pairs of waves by their cut amplitude products.
struct c_WaveAngular3D {
    std::vector<c_LegendreValue> legendre;                  // [angular pair] at the colatitude
    tides::c_ColatitudeTrig trig;                           // sine and cotangent factors of the colatitude
    std::vector<std::complex<double>> phasors;              // [azimuthal order] at the longitude
    std::vector<c_PotentialPointC> potentials;              // [wave] at the point
    std::vector<tides::c_AngularStrainFactors> factors;     // [wave] at the point, when strain factors were asked for
    // Per-wave unit stress and strain at the current radius, a workspace for c_secular_density_3d so the per-point
    // call allocates nothing; an angular state belongs to one thread.
    mutable std::vector<tides::c_Tensor6> unit_stress;
    mutable std::vector<tides::c_Tensor6> unit_strain;
    mutable std::vector<unsigned char> unit_valid;

    explicit c_WaveAngular3D(const c_WaveSet3D& set) :
        legendre(set.angular_pairs.size()),
        phasors(set.azimuthal_orders.size()),
        potentials(set.waves.size()),
        factors(set.waves.size()),
        unit_stress(set.waves.size()),
        unit_strain(set.waves.size()),
        unit_valid(set.waves.size(), 0) {}
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
            angular.phasors[static_cast<size_t>(set.wave_azimuthal_order[w])],
            set.secular ? std::complex<double>(1.0, 0.0) : set.waves[w].amplitude);
        if (strain_factors) {
            angular.factors[w] = tides::c_angular_strain_factors(angular.potentials[w], angular.trig);
        }
    }
}

// Fill the amplitudes at one point from its wave angular state and the strain radial coefficients of every
// radial group at one radius. Every wave with a shear kernel there is added into its frequency's total:
// waves at one frequency superpose, so every cross term between them is kept. False in a liquid layer,
// where no wave has a shear kernel.
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

// Signed azimuthal wavenumber; a wave's longitude structure is e^{i mu phi}.
inline int c_wave_mu(const c_TidalWave3D& wave) {
    return wave.azimuthal_sign * wave.order_m;
}

// The waves the secular density sums coherently: one group per frequency, or per (frequency, mu) for the
// longitude mean, over which waves of different mu average out. Groups are listed in order of first
// appearance and their members in wave order, which is the order the sums run in.
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

// Secular density at one point of a secular set: for each group, (|omega|/2) Im(sigma_c : conj(eps_c)) of the group's
// total stress and strain with every product of two amplitudes cut at e^N, i.e. the sum over wave pairs (a, b) of
// Im(P_ab sigma~_a : conj(eps~_b)), sigma~ and eps~ per unit amplitude and P_ab the cut product. coeffs_by_group[g] is
// radial group g's strain radial coefficient set at this radius (valid == false where the group has no shear kernel
// there), and the wave angular state carries unit-amplitude strain factors at the point (the longitude mean evaluates
// it at phi = 0).
inline double c_secular_density_3d(
        const c_WaveSet3D& set,
        const std::vector<c_SecularGroup3D>& groups,
        const std::vector<tides::c_StrainRadialCoeffs>& coeffs_by_group,
        const c_WaveAngular3D& angular) {
    if (!set.secular) {
        throw std::logic_error("TidalPy: the secular 3D heating needs a secular wave set (cut pair products)");
    }
    std::vector<tides::c_Tensor6>& stress = angular.unit_stress;
    std::vector<tides::c_Tensor6>& strain = angular.unit_strain;
    std::vector<unsigned char>& valid = angular.unit_valid;
    double heating = 0.0;
    for (const c_SecularGroup3D& group : groups) {
        for (const size_t w : group.waves) {
            const tides::c_StrainRadialCoeffs& radial = coeffs_by_group[set.wave_radial_group[w]];
            valid[w] = radial.valid ? 1 : 0;
            if (valid[w]) {
                tides::c_compute_strain_stress_from_factors(radial, angular.factors[w], strain[w], stress[w]);
            }
        }
        double group_heating = 0.0;
        for (const size_t a : group.waves) {
            if (!valid[a]) { continue; }
            for (const size_t b : group.waves) {
                if (!valid[b]) { continue; }
                const std::complex<double>& power = set.pair_power_of(a, b);
                if (std::abs(power) == 0.0) { continue; }
                double pair_heating = 0.0;
                for (size_t k = 0; k < 6; ++k) {
                    const double term = std::imag(power * stress[a].c[k] * std::conj(strain[b].c[k]));
                    pair_heating += (k < 3) ? term : 2.0 * term;
                }
                group_heating += pair_heating;
            }
        }
        heating += 0.5 * set.frequencies[group.frequency_group] * group_heating;
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

    // Radius: the user array, or Gauss-Legendre nodes inside each layer with r_wsum carrying the r^2
    // Jacobian. No node sits on a layer boundary, where the modulus and radial-solution lookups take the
    // layer below, so no node can weigh the lower layer's heating into the upper layer's integral.
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
    // Colatitude: the user array, or Gauss-Legendre in cos theta with the weight absorbing sin theta. A band
    // narrower than [0, pi] maps the nodes onto [cos(max), cos(min)] by affine substitution, the weights
    // scaling by the half-width, so the integral covers only that band.
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
    // Longitude: the user array, an analytic 2*pi on the longitude-mean secular density at a single point,
    // or a periodic trapezoid over [0, 2*pi) with no endpoint halving, the field being periodic.
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

// The angular Gram matrices of (l_a, l_b, m), computed once each.
struct c_GramMatrix3D {
    double values[6][6];
};
typedef std::map<std::array<int, 3>, c_GramMatrix3D> c_GramCache3D;

// Analytic colatitude integral of the longitude-mean secular density at one radius: over the (|omega|, mu) groups
// (c_secular_groups_3d split by mu) and the ordered wave pairs (a, b) within a group, (|omega|/2) int Im(c_a conj(c_b)
// sigma~_a : conj(eps~_b)) sin(theta) dtheta, the angular part coming from the Gram matrices. Times 2*pi*r^2 this is
// the radial power density dP/dr.
inline double c_secular_theta_integral_3d(
        const c_WaveSet3D& set,
        const std::vector<c_SecularGroup3D>& groups,
        const std::vector<tides::c_StrainRadialCoeffs>& coeffs_by_group,
        c_GramCache3D& gram_cache) {
    if (!set.secular) {
        throw std::logic_error("TidalPy: the secular 3D heating needs a secular wave set (cut pair products)");
    }
    double total = 0.0;
    for (const c_SecularGroup3D& group : groups) {
        const double frequency = set.frequencies[group.frequency_group];
        for (const size_t a : group.waves) {
            const tides::c_StrainRadialCoeffs& radial_a = coeffs_by_group[set.wave_radial_group[a]];
            if (!radial_a.valid) { continue; }
            const c_TidalWave3D& wave_a = set.waves[a];
            for (const size_t b : group.waves) {
                const tides::c_StrainRadialCoeffs& radial_b = coeffs_by_group[set.wave_radial_group[b]];
                if (!radial_b.valid) { continue; }
                const std::complex<double>& c_pair = set.pair_power_of(a, b);
                if (std::abs(c_pair) == 0.0) { continue; }
                const c_TidalWave3D& wave_b = set.waves[b];
                const std::array<int, 3> key{wave_a.degree_l, wave_b.degree_l, wave_a.order_m};
                auto found = gram_cache.find(key);
                if (found == gram_cache.end()) {
                    c_GramMatrix3D gram;
                    if (!tides::c_angular_gram_pair(wave_a.degree_l, wave_b.degree_l, wave_a.order_m, gram.values)) {
                        throw std::runtime_error(
                            "TidalPy: angular Gram matrix out of range during the analytic 3D heating collapse");
                    }
                    found = gram_cache.emplace(key, gram).first;
                }
                total += 0.5 * frequency * tides::c_theta_integrated_heating_pair(
                    radial_a,
                    radial_b,
                    wave_a.order_m,
                    wave_a.azimuthal_sign,
                    c_pair,
                    found->second.values);
            }
        }
    }
    return total;
}

}  // namespace tides3d

// The 3D orchestration lives on the rheology tide model, the only TideBase with a depth-resolved solution,
// and calls the world's members directly. It is defined here, in the world extension, where c_LayeredWorld
// and the kernel and potential headers are complete and CyRK lives, so every radial solve and dense call
// stays in its owning extension.

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

namespace tides3d {

// The rheology tide model a world's 3D path `what` runs on, once the world is checked to have one attached and its
// EOS solved; `plural` says whether `what` takes a plural verb.
inline const c_RheologyTide& c_require_3d_rheology(const c_LayeredWorld& world, const char* what, bool plural) {
    if (!world.get_tide_model_set()) {
        throw std::runtime_error("TidalPy: no tide model attached to the world. Call set_tide_model() first");
    }
    const auto* rheology = dynamic_cast<const c_RheologyTide*>(world.get_tide_model());
    if (rheology == nullptr) {
        throw std::runtime_error(
            std::string("TidalPy: ") + what + (plural ? " require" : " requires")
            + " the rheology tide model (the analytic cpl/ctl/ctl_q models have no depth-resolved radial solution)");
    }
    if (!world.get_eos_solved() || (world.get_eos_solution() == nullptr)) {
        throw std::runtime_error(
            std::string("TidalPy: ") + what + (plural ? " need" : " needs")
            + " the EOS solved first. Call solve_eos()");
    }
    return *rheology;
}

}  // namespace tides3d

// World delegation: validate preconditions, map the solve state into the potential's state struct, and
// hand off to the rheology tide model's 3D orchestration.
inline double c_LayeredWorld::get_3d_tidal_heating(
        const c_TideSolveConfig& state,
        double radius,
        double colatitude) {
    const c_WorldCallLock call_lock(this->p_call_mutex.get());
    const c_RheologyTide& rheology = tides3d::c_require_3d_rheology(*this, "3D tidal heating", false);
    return rheology.calc_3d_tidal_heating(*this, state, radius, colatitude);
}

// World delegation for the batch path: same preconditions as the scalar get_3d_tidal_heating.
inline void c_LayeredWorld::get_3d_tidal_heating_array(
        const c_TideSolveConfig& state,
        const double* radii,
        const double* colatitudes,
        size_t num_points,
        double* out_heating,
        int num_threads) {
    const c_WorldCallLock call_lock(this->p_call_mutex.get());
    const c_RheologyTide& rheology = tides3d::c_require_3d_rheology(*this, "3D tidal heating", false);
    rheology.calc_3d_tidal_heating_batch(
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
    const c_WorldCallLock call_lock(this->p_call_mutex.get());
    const c_RheologyTide& rheology = tides3d::c_require_3d_rheology(*this, "3D tidal displacements", true);
    rheology.calc_3d_displacements_grid(
        *this,
        state,
        axes,
        out_disp,
        num_threads);
}

// Instantaneous displacement grid. The coherent wave list is built once, the radial solve and the y1/y3
// samples once per radial group (l, |omega|). At each (r, theta, phi) every wave's complex displacement
// amplitude is added into its frequency's total, and each component at time t sums
// Re[amplitude e^{i |omega| t}] over the frequencies, the phase factors tabulated once. The colatitude rows
// run on up to num_threads threads, each writing only its own cells.
inline void c_RheologyTide::calc_3d_displacements_grid(
        c_LayeredWorld& world,
        const c_TideSolveConfig& state,
        const c_Grid3DAxes& axes,
        double* out_disp,
        int num_threads) const {
    const tides3d::c_WaveSet3D set = tides3d::c_world_wave_set_3d(world, state, "3D tidal displacements", false);
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

    // y1 and y3 per radial group and radius, [radius * num_groups + group]. A group with no depth-resolved
    // solution at a radius contributes nothing there, and a radius where no group has one is NaN.
    std::vector<std::complex<double>> y1_at(nr * num_groups);
    std::vector<std::complex<double>> y3_at(nr * num_groups);
    std::vector<unsigned char> group_missing(nr * num_groups, 0);
    std::vector<size_t> radius_missing(nr, 0);
    c_LoveSolveConfig love_cfg = world.make_radial_love_solve_config();
    c_LoveWorkspace workspace;
    for (size_t g = 0; g < num_groups; ++g) {
        const ::c_RadialSolutionStorage* storage = tides3d::c_solve_radial_group_3d(
            world, love_cfg, set.radial_groups[g], "3D tidal displacements", workspace);
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
// coefficients, and the phase tables are built once; at each (r, theta, phi) every wave's complex amplitude
// is added into its frequency's total, and each component at time t sums Re[amplitude e^{i |omega| t}] over
// the frequencies. The colatitude rows run on up to num_threads threads.
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
    const tides3d::c_WaveSet3D set = tides3d::c_world_wave_set_3d(world, state, what, false);
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
                // Undefined at a radius with no solution, and where waves exist but none has a shear kernel;
                // with no active waves at all the tensors are zero.
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
    const c_WorldCallLock call_lock(this->p_call_mutex.get());
    const c_RheologyTide& rheology = tides3d::c_require_3d_rheology(*this, "3D tidal stress and strain", true);
    rheology.calc_3d_stress_strain_grid(*this, state, axes, out_stress, out_strain, num_threads);
}

// Batch form of the secular 3D heating: the longitude-mean density at paired (radius, colatitude) points.
// The coherent wave list is built once, the radial solve runs once per radial group, and its strain radial
// coefficients are evaluated once per unique radius, since points on a map share radii. Points sharing a
// colatitude share its angular work, and the colatitudes run on up to num_threads threads.
inline void c_RheologyTide::calc_3d_tidal_heating_batch(
        c_LayeredWorld& world,
        const c_TideSolveConfig& state,
        const double* radii,
        const double* colatitudes,
        size_t num_points,
        double* out_heating,
        int num_threads) const {
    const tides3d::c_WaveSet3D set = tides3d::c_world_wave_set_3d(world, state, "secular 3D tidal heating", true);
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

// Collapsed (summed or averaged) 3D tidal heating: a full grid over (radius, colatitude, longitude[, time])
// or an integral along any spatial dimension, written into caller buffers. orbit_averaged gives the secular
// density h_bar, the pointwise time average or its longitude mean when longitude is summed; otherwise the
// instantaneous power sigma_ij(t) eps_dot_ij(t) at each user time. The radial solves run on the calling
// thread and the per-point evaluation on up to cfg.num_threads threads over colatitude rows.
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
    const tides3d::c_WaveSet3D set =
        tides3d::c_world_wave_set_3d(world, state, "3D tidal heating", cfg.orbit_averaged);

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
    // A summed axis contributes its integration weight with the Jacobian; a surviving spatial axis its
    // Jacobian when any axis is summed, else 1 for the raw density. Longitude and time have Jacobian 1.
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
    // Secular: each frequency's waves summed coherently, split by mu for the longitude mean, which the
    // single phi node then carries exactly, then (|omega|/2) Im(sigma_c : conj(eps_c)).
    const bool longitude_averaged = cfg.longitude_summed;
    const std::vector<tides3d::c_SecularGroup3D> groups = instantaneous
        ? std::vector<tides3d::c_SecularGroup3D>()
        : tides3d::c_secular_groups_3d(set, longitude_averaged);

    // The Gram path integrates the longitude-mean density, so it serves only a call that also sums longitude; a
    // call that keeps longitudes takes the Gauss-Legendre colatitude quadrature below.
    if (!instantaneous && cfg.latitude_summed && cfg.longitude_summed && cfg.latitude_analytic
            && grids.latitude_full_sphere) {
        // Integrate the longitude-mean secular density over theta with the Gram matrices, exactly and with
        // no theta grid. theta is summed away, so scatter over (radius, phi). With no per-point grid to
        // spread over threads, this runs on the calling thread.
        tides3d::c_GramCache3D gram_cache;
        for (size_t ir = 0; ir < nr; ++ir) {
            double theta_integral = 0.0;
            if (!radius_solve_failed[ir]) {
                theta_integral = tides3d::c_secular_theta_integral_3d(set, groups, coeffs[ir], gram_cache);
            }
            // theta is already integrated, the Gram absorbing sin theta, so only radius and longitude remain.
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

    // Evaluate over colatitude rows. A row owns its cells when colatitude survives and writes them directly.
    // When colatitude is summed every row adds into the same cells, so each fills a buffer of its own and
    // the rows are merged in row order below; either way the result is identical for any thread count.
    const bool rows_share_cells = !surv_th;
    std::vector<std::vector<double>> row_values(rows_share_cells ? nth : 0);
    std::vector<std::vector<double>> row_layer_totals(totals ? nth : 0);
    const size_t num_frequencies = set.frequencies.size();
    // Instantaneous: every wave's complex amplitude added into its frequency's total, each frequency evolved as
    // Re[. e^{i |omega| t}] with the phase factors tabulated once, and the real fields summed.
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

// World delegation for the collapse layout, from the layer geometry alone.
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
    const c_WorldCallLock call_lock(this->p_call_mutex.get());
    const c_RheologyTide& rheology = tides3d::c_require_3d_rheology(*this, "3D tidal heating", false);
    rheology.calc_3d_tidal_heating_collapsed(
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
