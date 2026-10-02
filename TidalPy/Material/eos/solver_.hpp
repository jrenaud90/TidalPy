#pragma once

#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "c_common.hpp"    // CyRK: DiffeqFuncType, PreEvalFunc, ODEMethod
#include "cysolution.hpp"  // CyRK: CySolverResult
#include "cysolve.hpp"     // CyRK: baseline_cysolve_ivp_noreturn
#include "c_events.hpp"    // CyRK: Event, EventFunc

#include "constants_.hpp"  // TidalPy: TidalPyConstants, d_PI, d_INF, d_EPS_100

#include "ode_.hpp"            // c_eos_diffeq, c_EOS_ODEInput, c_EOSSegment, c_eos_mass_event, C_EOS_Y_VALUES
#include "eos_solution_.hpp"   // c_EOSSolution, c_EOSPiece


/// How the structure integration bounds one layer.
///
/// A layer holds either its volume or its mass. One that holds its volume ends at a radius: the top it was given, or,
/// once a layer below has moved its base, the radius that keeps its volume. One that holds its mass ends where the
/// enclosed mass has gained that mass over its base, a terminal event, and every layer above it moves with it. A layer
/// whose material can change state carries the rigidity margin as a second terminal event: where it changes sign the
/// integration stops and restarts in the other state, so the solid and liquid zones of the layer are found by the
/// integration itself.
struct c_EOSLayerBounds
{
    // The geometry the layer's segments were laid out on [solve units], and the same in SI [m], from which the radius
    // grid of a layer the solve did not move is built exactly as the caller built it.
    double radius_inner = 0.0;
    double radius_outer = 0.0;
    double radius_inner_si = 0.0;
    double radius_outer_si = 0.0;
    // r_outer^3 - r_inner^3 of the layer's own volume [solve units^3], kept when a layer below moves its base.
    double volume_term = 0.0;
    // A layer holding its mass ends where it encloses target_mass [solve units] over its base.
    bool holds_mass = false;
    double target_mass = TidalPyConstants::d_NAN;
    // The rigidity margin of a layer whose state can change (c_preeval_rigidity_margin), null for one that cannot,
    // with the shear modulus [Pa] it is measured against; and the state of a layer that cannot change state.
    EventFunc state_event = nullptr;
    double minimum_shear_modulus = 0.0;
    bool liquid = false;
};


/// The scalar settings of a structure solve, in the units it runs in except where marked.
struct c_EOSSolverSettings
{
    double planet_bulk_density = 0.0;   // seeds the central-pressure guess
    double surface_pressure = 0.0;      // target surface pressure
    double G_to_use = -1.0;             // negative selects the shared config's SI value
    ODEMethod integration_method = ODEMethod::DOP853;
    double rtol = 1.0e-10;
    double atol = 1.0e-14;
    // Surface-pressure tolerance relative to the central-pressure scale (2/3) pi G rho_bulk^2 R^2 + surface_pressure.
    // Must exceed rtol, the integrator's own noise, to converge.
    double pressure_tol = 1.0e-8;
    size_t max_iters = 100;
    bool verbose = false;
    // Carry temperature and heat flow as two extra states, each segment's gradient form taken from its kind.
    bool integrate_temperature = false;
    // First central pressure to try: the last converged value when the caller has one. NaN or a non-positive value
    // starts from a uniform sphere at the bulk density.
    double central_pressure_guess = TidalPyConstants::d_NAN;
    // Radial samples per layer of the reported profile (>= 2), and the length [m] of one solve unit.
    size_t slices_per_layer = 100;
    double length_scale = 1.0;
};


/// What one pass of the structure integration produced.
struct c_EOSPassResult
{
    bool failed = false;
    std::string failure_message;
    // The pressure at the top of the outermost layer.
    double surface_pressure = TidalPyConstants::d_INF;
    // Each layer's top, and whether a layer below it that holds its mass moved it off the geometry it was given.
    std::vector<double> layer_tops;
    std::vector<char> layer_moved;
    std::vector<c_EOSPiece> pieces;
    // One integration per piece, kept only by a pass that captured its dense output.
    std::vector<std::unique_ptr<CySolverResult>> results;
};


/// One integration the pass asks for: a stretch of a layer's segment in one state, toward a stop radius.
struct c_EOSPieceRequest
{
    size_t layer_index = 0;
    size_t segment_index = 0;
    const c_EOSSegment* segment_ptr = nullptr;
    const c_EOSLayerBounds* bounds_ptr = nullptr;
    double radius_stop = 0.0;
    bool liquid = false;
};


// The most times the integration toward a layer's mass may be extended past its estimated top before the layer is
// taken to be unable to hold it, and the most pieces one segment may be split into.
inline constexpr size_t C_EOS_MAX_MASS_EXTENSIONS = 64;
inline constexpr size_t C_EOS_MAX_SEGMENT_PIECES = 4096;
// How far past its estimated top the first integration toward a layer's mass runs, as a fraction of the estimated
// stretch, so a layer near its last geometry ends in one piece.
inline constexpr double C_EOS_MASS_STOP_MARGIN = 0.25;
// A state root this close to the end of its stretch [multiples of machine epsilon, relative] leaves nothing to
// integrate in the new state there.
inline constexpr double C_EOS_ROOT_AT_STOP_EPS = 16.0;


/// One pass of the structure integration at a given central pressure, from the center through every layer and its
/// segments, each split into pieces by its events (c_EOSLayerBounds).
class c_EOSPassIntegrator
{
public:
    c_EOSPassIntegrator(
            c_EOSSolution* eos_solution_ptr,
            std::vector<PreEvalFunc>& eos_function_bylayer,
            std::vector<c_EOS_ODEInput>& eos_input_bylayer,
            const std::vector<c_EOSLayerBounds>& layer_bounds,
            const std::vector<c_EOSSegment>& segment_vec,
            const c_EOSSolverSettings& settings) :
                p_solution_ptr(eos_solution_ptr),
                p_functions(eos_function_bylayer),
                p_inputs(eos_input_bylayer),
                p_bounds(layer_bounds),
                p_segments(segment_vec),
                p_settings(settings),
                p_num_y(settings.integrate_temperature ? C_EOS_THERMAL_Y_VALUES : C_EOS_Y_VALUES),
                p_diffeq(settings.integrate_temperature ? c_eos_diffeq_thermal : c_eos_diffeq),
                p_args_vec(sizeof(c_EOS_ODEInput))
    {
    }

    size_t get_num_y() const noexcept { return this->p_num_y; }

    /// Integrate from the center state `y0` (num_y values) into `out`. A pass that captures its dense output keeps
    /// every piece's integration.
    void integrate(const double* y0, bool capture_dense, c_EOSPassResult& out)
    {
        const size_t num_layers   = this->p_bounds.size();
        const size_t num_segments = this->p_segments.size();
        out = c_EOSPassResult();
        out.layer_tops.assign(num_layers, TidalPyConstants::d_NAN);
        out.layer_moved.assign(num_layers, 0);

        std::vector<double> y_start(y0, y0 + this->p_num_y);
        double radius = 0.0;
        // Whether a layer holding its mass lies below, which moves every layer above it off its given geometry.
        bool moved = false;
        size_t segment_i = 0;
        for (size_t layer_i = 0; layer_i < num_layers; ++layer_i)
        {
            const c_EOSLayerBounds& bounds = this->p_bounds[layer_i];
            out.layer_moved[layer_i] = moved ? 1 : 0;
            const double radius_base = radius;
            const double mass_base   = y_start[C_EOS_MASS_INDEX];
            bool liquid              = bounds.liquid;
            bool layer_started       = false;

            while ((segment_i < num_segments) && (this->p_segments[segment_i].layer_index == layer_i))
            {
                const c_EOSSegment& segment = this->p_segments[segment_i];
                if ((segment_i > 0) && this->p_settings.integrate_temperature)
                {
                    // A segment that sets its own base temperature breaks the profile there, an isothermal
                    // layer against its neighbor; the rest continue from the segment below.
                    if (std::isfinite(segment.start_temperature))
                    {
                        y_start[4] = segment.start_temperature;
                    }
                    y_start[5] = segment.start_heat_flow;
                }
                const bool last_of_layer = (segment_i + 1 == num_segments)
                    || (this->p_segments[segment_i + 1].layer_index != layer_i);

                // The input the integration reads, copied into the integrator's arguments for every piece.
                c_EOS_ODEInput& layer_input = this->p_inputs[layer_i];
                layer_input.temperature_kind      = segment.temperature_kind;
                // The integration needs only the density, and a thermal one the thermal properties, so the moduli,
                // viscosities, and melt weakening are skipped.
                layer_input.full_state            = false;
                layer_input.thermal_state         = this->p_settings.integrate_temperature;
                layer_input.minimum_shear_modulus = bounds.minimum_shear_modulus;
                layer_input.target_mass           = bounds.holds_mass
                    ? mass_base + segment.upper_mass_fraction * bounds.target_mass : TidalPyConstants::d_NAN;
                std::memcpy(this->p_args_vec.data(), &layer_input, sizeof(c_EOS_ODEInput));

                // The layer's state where it starts; inside the layer it changes only at the roots of its margin.
                if (!layer_started && (bounds.state_event != nullptr))
                {
                    liquid = !(bounds.state_event(radius, y_start.data(), this->p_args_vec.data()) > 0.0);
                }
                layer_started = true;

                // Where the segment ends: its given radius, the radius that keeps its volume over a moved base, or,
                // for a layer holding its mass, an estimate past which the integration is extended until the mass
                // event fires.
                const double base_cubed = radius_base * radius_base * radius_base;
                const double volume_to_top = last_of_layer
                    ? bounds.volume_term
                    : (segment.upper_radius * segment.upper_radius * segment.upper_radius
                       - bounds.radius_inner * bounds.radius_inner * bounds.radius_inner);
                double radius_stop = segment.upper_radius;
                if (moved || bounds.holds_mass)
                {
                    radius_stop = std::cbrt(base_cubed + volume_to_top);
                }
                if (bounds.holds_mass)
                {
                    // A margin past the estimate, so a layer near its last geometry ends in one piece.
                    radius_stop += C_EOS_MASS_STOP_MARGIN * (radius_stop - radius)
                        + TidalPyConstants::d_EPS_100 * radius_stop;
                }

                c_EOSPieceRequest request;
                request.layer_index   = layer_i;
                request.segment_index = segment_i;
                request.segment_ptr   = &segment;
                request.bounds_ptr    = &bounds;
                size_t extensions = 0;
                size_t pieces_in_segment = 0;
                while (true)
                {
                    if (++pieces_in_segment > C_EOS_MAX_SEGMENT_PIECES)
                    {
                        out.failed = true;
                        out.failure_message = std::string("layer ") + std::to_string(layer_i)
                            + std::string(" changed state more than ") + std::to_string(C_EOS_MAX_SEGMENT_PIECES)
                            + std::string(" times inside one temperature segment");
                        return;
                    }
                    request.radius_stop = radius_stop;
                    request.liquid      = liquid;
                    const bool reached_end = this->p_integrate_piece(request, capture_dense, radius, y_start, out);
                    if (out.failed) { return; }
                    if (this->p_last_state_root)
                    {
                        liquid = !liquid;
                        // A root at the end of the stretch leaves nothing to integrate in the new state there.
                        const bool at_stop = (radius_stop - radius)
                            <= C_EOS_ROOT_AT_STOP_EPS * TidalPyConstants::d_EPS * std::fabs(radius_stop);
                        if (!at_stop) { continue; }
                        if (!bounds.holds_mass) { break; }
                    }
                    else
                    {
                        if (!reached_end) { break; }   // the mass event ended the segment
                        if (!bounds.holds_mass) { break; }
                    }
                    // The segment's mass lies further out: integrate on, a growing stretch at a time.
                    if (++extensions > C_EOS_MAX_MASS_EXTENSIONS)
                    {
                        out.failed = true;
                        out.failure_message = std::string("layer ") + std::to_string(layer_i)
                            + std::string(" holds its mass, but its integration reached ")
                            + std::to_string(radius) + std::string(" (solve units) without enclosing it; its "
                              "material's density may be near zero there");
                        return;
                    }
                    radius_stop = radius + 2.0 * std::max(radius - radius_base, TidalPyConstants::d_EPS_100);
                }
                ++segment_i;
            }
            out.layer_tops[layer_i] = radius;
            if (bounds.holds_mass) { moved = true; }
        }
        out.surface_pressure = y_start[C_EOS_PRESSURE_INDEX];
    }

private:
    /// One CyRK integration from `radius` toward the request's stop radius, which a terminal event may end early.
    /// Records the piece, moves `radius` and `y_start` to its end, and returns whether the integration reached the stop
    /// radius; whether a state root ended it is left in p_last_state_root.
    bool p_integrate_piece(
            const c_EOSPieceRequest& request,
            bool capture_dense,
            double& radius,
            std::vector<double>& y_start,
            c_EOSPassResult& out)
    {
        const c_EOSLayerBounds& bounds = *request.bounds_ptr;
        const bool liquid              = request.liquid;
        const double radius_stop       = request.radius_stop;
        this->p_last_state_root = false;
        // The state event looks only for the crossing out of the current state, so a restart on a root does not find
        // that root again. Only a pass that keeps its integrations reports zones, so a pass that only steers the
        // central-pressure search integrates through a state change without stopping there.
        std::vector<Event> events_vec;
        size_t state_event_index = static_cast<size_t>(-1);
        if ((bounds.state_event != nullptr) && capture_dense)
        {
            state_event_index = events_vec.size();
            events_vec.emplace_back(bounds.state_event, 1, liquid ? 1 : -1);
        }
        if (bounds.holds_mass)
        {
            events_vec.emplace_back(c_eos_mass_event, 1, 1);
        }

        if (!this->p_result_uptr)
        {
            this->p_result_uptr = std::make_unique<CySolverResult>(this->p_settings.integration_method);
        }
        CySolverResult* result_ptr = this->p_result_uptr.get();
        // Maximum step of one third of the stretch.
        const double max_step = 0.33 * (radius_stop - radius);
        baseline_cysolve_ivp_noreturn(
            result_ptr,
            this->p_diffeq,                         // Differential equation [DiffeqFuncType]
            radius,                                 // Start radius for this piece
            radius_stop,                            // Stop radius for this piece
            y_start,                                // y0 array vector<double>
            this->p_expected_size,                  // Expected final integration size [size_t]
            0,                                      // Number of extra outputs tracked [size_t]
            this->p_args_vec,                       // Extra input args to diffeq vector[char]
            this->p_max_num_steps,                  // Max number of steps [size_t]
            this->p_max_ram_MB,                     // Max amount of RAM allowed [size_t]
            capture_dense,                          // Use dense output [bool]
            this->p_t_eval_vec,                     // Interpolate at radius array vector[double]
            this->p_functions[request.layer_index], // Pre-eval function used in diffeq [PreEvalFunc]
            events_vec,                             // Events vector [vector<Event>]
            this->p_rtols_vec,                      // Relative Tolerance vector[double]
            this->p_atols_vec,                      // Absolute Tolerance vector[double]
            max_step,                               // Maximum step size [double]
            0.0,                                    // Initial step size [double]
            true,                                   // Force retain solver [bool]
            nullptr                                 // Analytic jacobian (null = numerical) [JacobianFuncType]
        );
        this->p_solution_ptr->save_steps_taken(result_ptr->steps_taken);
        if (!result_ptr->success)
        {
            out.failed = true;
            // Where it failed: the layer, its segment, and the stretch (solve units) the piece covered.
            out.failure_message = result_ptr->message + std::string(" (layer ") + std::to_string(request.layer_index)
                + std::string(", segment ") + std::to_string(request.segment_index) + std::string(", from r = ")
                + std::to_string(radius) + std::string(" toward ") + std::to_string(radius_stop)
                + std::string(" in solve units)");
            return false;
        }

        // A terminal event stopped the integration at its root, the last stored step.
        const size_t last_index      = result_ptr->size - 1;
        const bool   terminated      = result_ptr->event_terminated;
        const double radius_end      = terminated ? result_ptr->time_domain_vec[last_index] : radius_stop;
        this->p_last_state_root      = terminated && (result_ptr->event_terminate_index == state_event_index);
        std::memcpy(y_start.data(), &result_ptr->solution[this->p_num_y * last_index], sizeof(double) * this->p_num_y);

        c_EOSPiece piece;
        piece.layer_index   = request.layer_index;
        piece.segment_index = request.segment_index;
        piece.radius_lower  = radius;
        piece.radius_upper  = radius_end;
        piece.liquid        = liquid;
        piece.temperature   = request.segment_ptr->start_temperature;
        out.pieces.push_back(piece);
        if (capture_dense)
        {
            // Held until the pass is judged; a converged pass hands these to the solution.
            out.results.push_back(std::move(this->p_result_uptr));
        }
        radius = radius_end;
        return !terminated;
    }

    c_EOSSolution*                      p_solution_ptr;
    std::vector<PreEvalFunc>&           p_functions;
    std::vector<c_EOS_ODEInput>&        p_inputs;
    const std::vector<c_EOSLayerBounds>& p_bounds;
    const std::vector<c_EOSSegment>&    p_segments;
    const c_EOSSolverSettings&          p_settings;
    size_t                              p_num_y;
    DiffeqFuncType                      p_diffeq;
    // CyRK reads its arguments from a byte buffer and takes the number of states from y0, so both are kept sized.
    std::vector<char>                   p_args_vec;
    // One tolerance for every y; CyRK still takes vectors.
    std::vector<double>                 p_rtols_vec  = {this->p_settings.rtol};
    std::vector<double>                 p_atols_vec  = {this->p_settings.atol};
    std::vector<double>                 p_t_eval_vec = std::vector<double>(0);
    size_t                              p_expected_size = 64;
    size_t                              p_max_num_steps = 10000;
    size_t                              p_max_ram_MB    = 512;
    // Reused between pieces of a pass that keeps no dense output.
    std::unique_ptr<CySolverResult>     p_result_uptr;
    bool                                p_last_state_root = false;
};


/// Hand a converged pass to the solution: the layer tops, the reported radius grid, the pieces, and their integrations.
/// A layer the solve did not move samples the caller's SI geometry exactly as the caller would; one that moved, or that
/// holds its mass, samples its solved radii.
inline void c_commit_eos_pass(
        c_EOSSolution* eos_solution_ptr,
        const std::vector<c_EOSLayerBounds>& layer_bounds,
        const c_EOSSolverSettings& settings,
        c_EOSPassResult& pass)
{
    const size_t num_layers = layer_bounds.size();
    const size_t slices     = settings.slices_per_layer;
    std::vector<double> grid(num_layers * slices);
    for (size_t layer_i = 0; layer_i < num_layers; ++layer_i)
    {
        const c_EOSLayerBounds& bounds = layer_bounds[layer_i];
        const bool kept = !pass.layer_moved[layer_i] && !bounds.holds_mass;
        const double radius_inner = kept ? bounds.radius_inner_si
            : ((layer_i == 0) ? 0.0 : pass.layer_tops[layer_i - 1]);
        const double radius_outer = kept ? bounds.radius_outer_si : pass.layer_tops[layer_i];
        const double scale        = kept ? settings.length_scale : 1.0;
        for (size_t slice_i = 0; slice_i < slices; ++slice_i)
        {
            const double fraction = static_cast<double>(slice_i) / static_cast<double>(slices - 1);
            grid[layer_i * slices + slice_i] = (radius_inner + fraction * (radius_outer - radius_inner)) / scale;
        }
    }
    eos_solution_ptr->upper_radius_bylayer_vec = pass.layer_tops;
    eos_solution_ptr->change_radius_array(grid.data(), grid.size());
    eos_solution_ptr->set_pieces(pass.pieces);
    for (std::unique_ptr<CySolverResult>& result_uptr : pass.results)
    {
        eos_solution_ptr->save_cyresult(std::move(result_uptr));
    }
    pass.results.clear();
}


/// Solve the equation of state for a layered planet.
///
/// Integrates gravity, pressure, mass, and moment of inertia from center to surface in the caller's units, layer by
/// layer and segment by segment, each split into pieces where an event ends it (c_EOSLayerBounds): the enclosed mass of
/// a layer that holds its mass, and the rigidity margin of a layer that can change state. The central pressure comes
/// from a safeguarded secant iteration on the surface-pressure mismatch at the top of the outermost layer, wherever
/// the layers holding their mass put it. The first update assumes a unit slope, exact for an incompressible planet,
/// and later ones use the measured slope. Where that slope is not positive (the surface pressure of a compressible
/// planet can fall as the central pressure first rises), the step grows geometrically in the unit-slope direction
/// until the mismatch changes sign, and until then no secant step may grow more than fourfold over the one before.
/// Once it has, the root is bracketed: a secant step that would leave the bracket is replaced by an Illinois
/// false-position step inside it, and any step not under half the one two passes back by a bisection (Brent's
/// safeguard), so a mismatch with a near-jump cannot stall the iteration.
///
/// Only the pass that converges needs its dense output, and capturing it roughly doubles the cost of a
/// pass, so it is switched on for the first pass of a warm start (which usually converges there) and for
/// any pass the secant's own error model expects to converge. A pass that converges without it is repeated. The same
/// passes alone stop at state changes: the zones are reported from the kept pass, and the others only steer the
/// central pressure.
///
/// The converged pass sets the solution's layer tops, its pieces, and its reported radius grid (slices_per_layer per
/// layer); a layer the solve did not move keeps the grid the caller's SI geometry gives.
///
/// Parameters
/// ----------
/// eos_solution_ptr : c_EOSSolution*
///     Output; must be constructed with the number of layers alone (its tops and grid are outputs).
/// eos_function_bylayer_ptr_vec, eos_input_bylayer_vec
///     EOS evaluation function and arguments per layer; used during the integration and kept for later.
/// layer_bounds : vector of c_EOSLayerBounds
///     How each layer ends, inner to outer.
/// segment_vec : vector of c_EOSSegment
///     The temperature segments, ascending, each layer's a contiguous run with at least one per layer.
/// settings : c_EOSSolverSettings
///
/// Assumptions
/// -----------
/// - Spherical symmetry and hydrostatic equilibrium.
/// - The surface pressure rises with the central pressure near the root, and the unit-slope direction leads
///   toward it: too high a surface pressure lowers the central pressure. A structure still off its target surface
///   pressure at `max_iters` is reported as a failure.
/// - A state change inside one integration step that starts and ends in the same state is not found: such a band is
///   thinner than the step's own error, because the temperature and pressure are smooth inside a segment.
inline void c_solve_eos(
        c_EOSSolution* eos_solution_ptr,
        std::vector<PreEvalFunc>& eos_function_bylayer_ptr_vec,
        std::vector<c_EOS_ODEInput>& eos_input_bylayer_vec,
        const std::vector<c_EOSLayerBounds>& layer_bounds,
        const std::vector<c_EOSSegment>& segment_vec,
        const c_EOSSolverSettings& settings)
{
    eos_solution_ptr->message = std::string("Equation of state solver finished without issue.");

    const double G_to_use             = (settings.G_to_use < 0.0) ? c_get_G() : settings.G_to_use;
    const double surface_pressure     = settings.surface_pressure;
    const double planet_bulk_density  = settings.planet_bulk_density;
    const double pressure_tol         = settings.pressure_tol;
    const size_t max_iters            = settings.max_iters;
    const bool   verbose              = settings.verbose;
    const bool   integrate_temperature = settings.integrate_temperature && !segment_vec.empty();
    const size_t num_layers           = layer_bounds.size();

    // The outermost layer's given top, as the caller's grid reaches it.
    const c_EOSLayerBounds& outermost = layer_bounds.back();
    const double planet_radius = (outermost.radius_inner_si + (outermost.radius_outer_si - outermost.radius_inner_si))
        / settings.length_scale;
    const double r0_gravity        = 0.0;
    const double r0_pressure_guess = (
        (2.0 / 3.0) * TidalPyConstants::d_PI * G_to_use * planet_radius *
        planet_radius * planet_bulk_density * planet_bulk_density)
        + surface_pressure;

    const double r0_mass = 0.0;
    const double r0_moi  = 0.0;

    // One isothermal segment per layer unless the caller supplied a layout.
    std::vector<c_EOSSegment> default_segments;
    if (segment_vec.empty())
    {
        default_segments.resize(num_layers);
        for (size_t layer_i = 0; layer_i < num_layers; ++layer_i)
        {
            default_segments[layer_i].upper_radius = layer_bounds[layer_i].radius_outer;
            default_segments[layer_i].layer_index  = layer_i;
        }
    }
    const std::vector<c_EOSSegment>& segments = segment_vec.empty() ? default_segments : segment_vec;
    c_EOSSolverSettings pass_settings = settings;
    pass_settings.integrate_temperature = integrate_temperature;
    c_EOSPassIntegrator integrator(
        eos_solution_ptr, eos_function_bylayer_ptr_vec, eos_input_bylayer_vec, layer_bounds, segments, pass_settings);

    // The caller's central pressure, or that of a uniform sphere at the bulk density. A thermal solve
    // starts at the innermost segment's temperature with whatever heat flow it carries, zero at a regular
    // center.
    const double guess       = settings.central_pressure_guess;
    const bool warm_start    = std::isfinite(guess) && (guess > 0.0);
    const double r0_pressure = warm_start ? guess : r0_pressure_guess;
    double y0[C_EOS_THERMAL_Y_VALUES] = {r0_gravity, r0_pressure, r0_mass, r0_moi, 0.0, 0.0};
    if (integrate_temperature)
    {
        y0[4] = segments[0].start_temperature;
        y0[5] = segments[0].start_heat_flow;
    }
    const size_t num_y = integrator.get_num_y();
    eos_solution_ptr->num_y_solved = num_y;

    // Secant state on f(P_c) = P_surface(P_c) - surface_pressure. The convergence test is relative to the
    // central-pressure scale, the size of the integrator's own noise on the surface pressure.
    const double pressure_scale   = (r0_pressure_guess > 0.0) ? r0_pressure_guess : 1.0;
    const double pressure_tol_abs = pressure_tol * pressure_scale;
    double previous_central       = TidalPyConstants::d_NAN;
    double previous_diff          = TidalPyConstants::d_NAN;
    double oldest_diff            = TidalPyConstants::d_NAN;
    // The latest central pressure on each side of the root, with its mismatch; the root lies between them once
    // both are known. `last_replaced` is the side the last pass replaced (-1 below, +1 above, 0 none), and
    // `last_bracket_step` whether the step to it replaced a secant step inside the bracket (false position or
    // bisection), for the Illinois halving of an end kept twice.
    double below_central          = TidalPyConstants::d_NAN;   // mismatch < 0 there
    double below_diff             = TidalPyConstants::d_NAN;
    double above_central          = TidalPyConstants::d_NAN;   // mismatch >= 0 there
    double above_diff             = TidalPyConstants::d_NAN;
    int    last_replaced          = 0;
    bool   last_bracket_step      = false;
    // Sizes of the last two steps: a run of unusable slopes doubles the last, and the one before it bounds a
    // bracketed step (Brent's safeguard).
    double last_step_abs          = 0.0;
    double step_two_ago_abs       = TidalPyConstants::d_INF;

    double calculated_surf_pressure = TidalPyConstants::d_INF;
    double pressure_diff            = TidalPyConstants::d_INF;
    double pressure_diff_abs        = TidalPyConstants::d_INF;
    int iterations                  = 0;
    bool failed                     = false;
    bool max_iters_hit              = false;
    // Whether this pass keeps its dense output, and whether it is kept whatever it finds: the repeat of a
    // converged pass, or the pass after the iteration cap. Capturing roughly doubles the cost of a pass,
    // so the first pass does it only from a warm start, which usually converges there.
    bool capture_dense              = warm_start;
    bool final_pass                 = false;
    c_EOSPassResult pass;
    std::string integrator_failure_message;

    while (true)
    {
        calculated_surf_pressure = TidalPyConstants::d_INF;

        if (!final_pass)
        {
            iterations++;
        }

        integrator.integrate(y0, capture_dense, pass);
        if (pass.failed)
        {
            failed = true;
            integrator_failure_message = pass.failure_message;
            break;
        }
        calculated_surf_pressure = pass.surface_pressure;

        pressure_diff     = calculated_surf_pressure - surface_pressure;
        pressure_diff_abs = std::fabs(pressure_diff);
        const bool converged = (pressure_diff_abs <= pressure_tol_abs);

        if (final_pass || (converged && capture_dense))
        {
            // The pass's layout and the whole result of every piece are kept for later dense calls.
            c_commit_eos_pass(eos_solution_ptr, layer_bounds, settings, pass);
            break;
        }
        else
        {
            if (converged)
            {
                // Converged without its dense output, so run the same central pressure again and keep it.
                final_pass    = true;
                capture_dense = true;
            }
            else
            {
                // Record which side of the root this pass landed on. Two false-position passes in a row on one side
                // leave the other end stale, so its mismatch is halved (Illinois), which keeps the step moving.
                const double current_central = y0[1];
                if (pressure_diff < 0.0)
                {
                    if (last_bracket_step && (last_replaced == -1)) { above_diff *= 0.5; }
                    below_central = current_central;
                    below_diff    = pressure_diff;
                    last_replaced = -1;
                }
                else
                {
                    if (last_bracket_step && (last_replaced == +1)) { below_diff *= 0.5; }
                    above_central = current_central;
                    above_diff    = pressure_diff;
                    last_replaced = +1;
                }
                const bool bracketed = std::isfinite(below_central) && std::isfinite(above_central);

                // The secant step on the measured slope when it is positive (and, once bracketed, stays inside
                // the bracket). Otherwise a false-position step inside the bracket, or before there is one, a
                // unit-slope step that doubles on every pass the slope stays unusable, so a mismatch that first
                // moves the wrong way costs a few passes rather than one per residual's worth of pressure.
                double step      = TidalPyConstants::d_NAN;
                bool   secant_ok = false;
                last_bracket_step = false;
                if (std::isfinite(previous_central))
                {
                    const double slope = (pressure_diff - previous_diff) / (current_central - previous_central);
                    if (std::isfinite(slope) && (slope > 0.0))
                    {
                        step      = -pressure_diff / slope;
                        secant_ok = std::isfinite(step);
                    }
                }
                if (bracketed)
                {
                    const double bracket_lo = std::min(below_central, above_central);
                    const double bracket_hi = std::max(below_central, above_central);
                    double next_central = current_central + step;
                    if (!(secant_ok && (next_central > bracket_lo) && (next_central < bracket_hi)))
                    {
                        next_central = above_central
                            - above_diff * (above_central - below_central) / (above_diff - below_diff);
                        last_bracket_step = true;
                    }
                    // A converging step is under half the one two passes back. One that is not (a mismatch with a
                    // near-jump, from an EOS at the end of its range, stalls both the secant and false position)
                    // is replaced by a bisection, which shrinks the bracket steadily.
                    const bool stalled = std::fabs(next_central - current_central) > 0.5 * step_two_ago_abs;
                    if (stalled || !(std::isfinite(next_central) && (next_central > bracket_lo)
                                     && (next_central < bracket_hi)))
                    {
                        next_central      = 0.5 * (bracket_lo + bracket_hi);
                        last_bracket_step = true;
                    }
                    step = next_central - current_central;
                }
                else if (!secant_ok)
                {
                    const double direction = (pressure_diff > 0.0) ? -1.0 : 1.0;
                    step = direction * std::max(std::fabs(pressure_diff), 2.0 * last_step_abs);
                }
                else if ((last_step_abs > 0.0) && (std::fabs(step) > 4.0 * last_step_abs))
                {
                    // Before the root is bracketed a nearly flat mismatch can send the secant orders of magnitude
                    // past it; growing at most fourfold per pass still reaches a distant root geometrically.
                    step = std::copysign(4.0 * last_step_abs, step);
                }
                // The secant error model e[k+1] = M e[k] e[k-1], with M measured from the residuals in
                // hand, says whether the next pass should converge and so whether to capture its dense
                // output. With one residual there is no model yet, and a wrong guess costs more than a
                // repeated pass saves, so that pass runs without it.
                const double reference_diff = std::isfinite(oldest_diff) ? oldest_diff : previous_diff;
                const double predicted_diff = std::isfinite(reference_diff)
                    ? pressure_diff * pressure_diff / std::fabs(reference_diff) : TidalPyConstants::d_INF;
                capture_dense = (predicted_diff <= pressure_tol_abs);

                oldest_diff      = previous_diff;
                previous_central = y0[1];
                previous_diff    = pressure_diff;

                // A non-finite step has nowhere to go, and this pass kept no dense output to report.
                if (!std::isfinite(step))
                {
                    failed = true;
                    integrator_failure_message =
                        std::string("the secant step on the central pressure is not finite (surface pressure "
                                    "residual ") + std::to_string(pressure_diff) + std::string(" Pa)");
                    break;
                }
                // Keep the central pressure positive by halving an overshooting step.
                double next_central = y0[1] + step;
                while ((next_central <= 0.0) && (std::fabs(step) > TidalPyConstants::d_EPS * pressure_scale))
                {
                    step        *= 0.5;
                    next_central = y0[1] + step;
                }
                // Before a second step there is no step two passes back to bound the next one.
                step_two_ago_abs = (last_step_abs > 0.0) ? last_step_abs : TidalPyConstants::d_INF;
                last_step_abs    = std::fabs(step);
                y0[1] = next_central;
            }
        }

        if (!final_pass && (iterations >= static_cast<int>(max_iters)))
        {
            max_iters_hit = true;
            eos_solution_ptr->max_iters_hit = true;
            // Still produce output: one more pass, kept whatever it finds.
            final_pass    = true;
            capture_dense = true;
        }
    }
    eos_solution_ptr->iterations = iterations;

    // The pass the cap forces is kept for its diagnostics, but a structure whose surface pressure misses the
    // target is not hydrostatic and must not be reported as solved.
    const bool unconverged = max_iters_hit && !failed && !(pressure_diff_abs <= pressure_tol_abs);

    if (failed)
    {
        eos_solution_ptr->success = false;
        if (iterations > 1)
        {
            // A later pass fails where the search for the central pressure has taken it, usually far past any
            // structure the layers can hold.
            eos_solution_ptr->message =
                std::string("`c_solve_eos` found no hydrostatic structure: the structure integration failed at "
                            "iteration ") + std::to_string(iterations) + std::string(", at a central pressure of ") +
                std::to_string(y0[1]) + std::string(" (") + std::to_string(y0[1] / pressure_scale) +
                std::string(" times the uniform-sphere estimate), while searching for the central pressure that "
                            "meets the target surface pressure. The layers' equations of state may have no "
                            "hydrostatic solution at this radius and mass");
        }
        else
        {
            eos_solution_ptr->message =
                std::string("Warning in `c_solve_eos`: Integrator failed at iteration ") + std::to_string(iterations);
        }
        if (!integrator_failure_message.empty())
        {
            eos_solution_ptr->message += std::string(". Message: ") + integrator_failure_message;
        }
        if (verbose)
        {
            std::printf("%s", eos_solution_ptr->message.c_str());
        }
    }
    else
    {
        eos_solution_ptr->success = !unconverged;
        eos_solution_ptr->pressure_error = pressure_diff_abs;
        if (unconverged)
        {
            eos_solution_ptr->message =
                std::string("`c_solve_eos` found no hydrostatic structure: after ") + std::to_string(iterations) +
                std::string(" iterations the surface pressure misses its target by ") +
                std::to_string(pressure_diff_abs) + std::string(" Pa (tolerance ") +
                std::to_string(pressure_tol_abs) + std::string(" Pa). The layers' equations of state may have no "
                "hydrostatic solution at this radius and mass; otherwise raise `max_iters`.");
            if (verbose)
            {
                std::printf("%s\n", eos_solution_ptr->message.c_str());
            }
        }

        // Keep the layer EOS functions for on-demand evaluation, then sample onto the radius array.
        eos_solution_ptr->save_eos_functions(eos_function_bylayer_ptr_vec, eos_input_bylayer_vec);
        eos_solution_ptr->interpolate_full_planet();
    }
}
