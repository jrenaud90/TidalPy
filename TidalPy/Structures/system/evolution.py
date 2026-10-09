"""Coupled orbital, spin, and thermal evolution of a world about its tidal host, its spin tracking its spin-orbit
equilibria (:meth:`System.evolve <TidalPy.Structures.system.System.evolve>`).

Near a stable spin-orbit equilibrium a world's spin relaxes onto it within years to kyr, while its orbit and interior
change over Myr to Gyr. Integrated directly, the spin forces an implicit integrator to steps of its relaxation time, and
an equilibrium narrower than the integrator's tolerance on the spin (about 1e-6 in spin ratio in a cold mantle) is
stepped over. Here the spin is *free* (a state variable) between equilibria and *tracked* on one: while a stable root
s*(x) of the spin balance b(s; x) = ds/dt exists for the slow state x (the orbit, the host's spin, the layer
temperatures), the spin sits on it and only x is integrated (a quasi-static spin; Walterova and Behounkova 2020). Each
tracked right-hand side finds s*(x) from the last root and returns the Filippov combination of the rates on the root's
two sides, so the spin and the orbit take one torque.

Transitions are CyRK terminal events, each a function of the state alone and positive at its segment's start:

* a free spin entering the band |s - k/2| < ``capture_band`` of the next commensurability, or coming within
  ``capture_margin`` of a stable root found ahead of it, runs a capture test: the first sign change of b in the
  direction the spin moves brackets the stable root it reaches;
* a tracked spin's window (an interval about its root in which sample points spaced by factors of 4 show no other sign
  change) has an edge whose balance loses its restoring sign: the root is solved again there and the window rebuilt
  about it, or, with no root left nearby (the equilibrium vanished with its unstable neighbor), the spin is released.

Assumptions
-----------
* The spin sits exactly on a tracked equilibrium; the neglected lag is its relaxation time over the evolution time.
* The world raises tides on the orbit about its tidal host and the host on the world (``calc_pair_evolution``); the
  host's spin is always free.
"""
import math
import time as clock
from dataclasses import dataclass, field

import numpy as np
from CyRK import pysolve_ivp
from CyRK.cy.pysolver import PySolver

import TidalPy.constants as tidalpy_constants
from TidalPy.Utilities.logging import log_debug, log_info

# The driver integrates in Myr, so the balances it compares are of order one.
TIME_UNIT = 1.0e6 * tidalpy_constants.year   # [s]
# Half the spacing of the half-integer commensurabilities k/2 [spin ratio].
HALF_SPACING = 0.25
# The finest search offset keeps the resonant mode's frequency 2 n |s - k/2| this many times above [numerical]
# minimum_frequency, below which the mode is static.
STATIC_MARGIN = 30.0
# Spin-ratio spacing of the four probes that sample the balance's noise at a root.
NOISE_STEP = 1.0e-13
# Consecutive segments that end where they began, or that end on a failed evaluation, before the run stops.
MAX_STALLED_SEGMENTS = 5
MAX_FAILED_SEGMENTS = 10


class ProbeFailed(RuntimeError):
    """A tide solve at one spin ratio raised or returned non-finite rates."""


@dataclass
class SpinProbe:
    """The rates of a fixed slow state at one spin ratio.

    Attributes
    ----------
    spin_ratio : float
        The world's spin frequency over its orbital motion.
    balance : float
        ds/dt [Myr-1].
    rates : np.ndarray
        The slow state's scaled rates [Myr-1].
    heating : float
        The world's tidal heating [W].
    weakest : bool
        Set when no root exists and this is the weakest balance found in its place.
    """
    spin_ratio: float
    balance: float
    rates: np.ndarray
    heating: float
    weakest: bool = False


@dataclass
class SpinWindow:
    """An interval [lower, upper] in spin ratio about a stable root: the balance is positive at `lower`, negative at
    `upper`, and at its sample points between them it changes sign only at the root. `noise` is the balance's noise at
    the root [Myr-1]."""
    lower: float
    upper: float
    noise: float


@dataclass
class EvolutionResult:
    """The evolution of a world about its tidal host (:meth:`System.evolve`), at every step the integrator took.

    Attributes
    ----------
    time : np.ndarray
        [s]
    semi_major_axis : np.ndarray
        [m]
    eccentricity : np.ndarray
    spin_ratio : np.ndarray
        The world's spin frequency over its orbital motion.
    spin_frequency, host_spin_frequency : np.ndarray
        [rad s-1]
    temperature : np.ndarray
        Each layer's temperature [K], shape (layers, steps); ``None`` when the thermal state was not evolved.
    tidal_heating : np.ndarray
        The world's tidal heating [W].
    tracked : np.ndarray
        True where the spin sat on an equilibrium.
    segments : list of dict
        One entry per integration segment: its start ``time`` [s], ``mode`` ("free", "approach", or "tracked"), the
        ``spin_ratio`` it started from, and for a tracked segment its ``root`` and ``window``, for an approach its
        ``target``, and for a free one its ``reason``; ``ended`` says how it ended when it did not reach the end.
    success : bool
        The run reached the end of its span and left the system at its final state.
    message : str
    num_tide_solves, num_eos_solves : int
        Pair evolutions (each a tide solve of the world and of its host) and EOS solves.
    elapsed : float
        Wall-clock time [s].
    """
    time: np.ndarray
    semi_major_axis: np.ndarray
    eccentricity: np.ndarray
    spin_ratio: np.ndarray
    spin_frequency: np.ndarray
    host_spin_frequency: np.ndarray
    temperature: object
    tidal_heating: np.ndarray
    tracked: np.ndarray
    segments: list = field(default_factory=list)
    success: bool = False
    message: str = ""
    num_tide_solves: int = 0
    num_eos_solves: int = 0
    elapsed: float = 0.0


# ======================================================================================================================
# The rates at a fixed slow state
# ======================================================================================================================
class PairRates:
    """The rates of a world and its tidal host at one slow state, as functions of the world's spin ratio.

    The slow state x is [a / a0, e, host spin / its initial value] and, with the thermal state evolved, each layer's
    temperature over its initial value. With the thermal state evolved, a new slow state solves the world's EOS with its
    temperature profile, its surface at `surface_temperature` (None: the system's insolation temperature when it has a
    star, else no flow through the surface), and its heat sources at the integration time.

    Assumptions
    -----------
    * The eccentricity stays in [0, 1); a trial state outside it is clipped to the nearest valid value.
    """

    def __init__(
            self,
            system,
            world,
            host,
            evolve_thermal,
            surface_temperature,
            deadline):
        self.system, self.world, self.host = system, world, host
        self.evolve_thermal = evolve_thermal
        self.surface_temperature = surface_temperature
        self.deadline = deadline
        self.semi_major_axis0 = system.get_semi_major_axis(world)
        self.host_spin0 = host.spin_frequency if host.spin_frequency > 0.0 else system.calc_orbital_frequency(world)
        self.temperature0 = np.array([layer.temperature for layer in world])
        self.num_layers = len(world)
        self.orbital_frequency = math.nan
        self.num_probes = 0
        self.num_eos_solves = 0
        self.p_key = None

    def initial_state(self):
        """The slow state of the system as it is [scaled]."""
        state = [1.0, self.system.get_eccentricity(self.world), self.host.spin_frequency / self.host_spin0]
        if self.evolve_thermal:
            state.extend([1.0] * self.num_layers)
        return np.array(state)

    def set_slow_state(self, time, x):
        """Puts the system at slow state `x` at `time` [Myr]: the orbit, the host's spin, and with the thermal state
        evolved the layer temperatures and the world's EOS. A repeated (time, x) is not solved again; a failed EOS solve
        raises (SolutionFailedError, a RuntimeError)."""
        key = (time, x.tobytes())
        if key == self.p_key:
            return
        self.p_key = None   # Until the solve below succeeds
        self.system.set_semi_major_axis(self.world, x[0] * self.semi_major_axis0)
        self.system.set_eccentricity(self.world, min(max(x[1], 0.0), 1.0 - 1.0e-12))
        self.host.set_spin_frequency(x[2] * self.host_spin0)
        self.orbital_frequency = self.system.calc_orbital_frequency(self.world)
        if self.evolve_thermal:
            for layer, temperature in zip(self.world, x[3:] * self.temperature0):
                layer.temperature = temperature
            self.world.clear_tidal_heating()
            surface_temperature = self.surface_temperature
            if (surface_temperature is None) and self.system.has_star:
                surface_temperature = self.system.calc_equilibrium_temperature(self.world)
            self.world.solve_eos(
                solve_temperature=True,
                surface_temperature=surface_temperature,
                time=time * TIME_UNIT,
                raise_on_fail=True)
            self.num_eos_solves += 1
        self.p_key = key

    def finest_offset(self, finest):
        """The finest spin-ratio offset searched about a commensurability: `finest`, or more where the resonant mode's
        frequency 2 n |s - k/2| would come within STATIC_MARGIN of [numerical] minimum_frequency."""
        return max(finest, STATIC_MARGIN * tidalpy_constants.min_frequency / (2.0 * self.orbital_frequency))

    def probe(self, spin_ratio):
        """The rates at `spin_ratio` (a SpinProbe). The world keeps that spin and its tide solve. Raises ProbeFailed
        when the tide solve raises or returns non-finite rates, and TimeoutError past the deadline."""
        if (self.deadline is not None) and (clock.perf_counter() > self.deadline):
            raise TimeoutError("TidalPy: System.evolve reached its max_wall_time.")
        n = self.orbital_frequency
        self.world.set_spin_frequency(spin_ratio * n)
        try:
            pair = self.system.calc_pair_evolution(self.world)
        except (RuntimeError, ValueError) as error:
            lines = str(error).strip().splitlines()
            reason = lines[-1].strip() if lines else type(error).__name__
            raise ProbeFailed(f"TidalPy: the tide solve failed at spin ratio {spin_ratio:.15f}: {reason}") from error
        self.num_probes += 1
        world_part, host_part = (pair["worlds"][name] for name in pair["world_names"])
        balance = (world_part["dspin_dt"] - spin_ratio * pair["dn_dt"]) / n * TIME_UNIT
        rates = [pair["da_dt"] / self.semi_major_axis0, pair["de_dt"], host_part["dspin_dt"] / self.host_spin0]
        if self.evolve_thermal:
            rates.extend(self.world.calc_layer_temperature_rate(layer_i) / self.temperature0[layer_i]
                         for layer_i in range(self.num_layers))
        rates = np.array(rates) * TIME_UNIT
        if not (math.isfinite(balance) and np.all(np.isfinite(rates))):
            raise ProbeFailed(f"TidalPy: the rates at spin ratio {spin_ratio:.15f} are not finite.")
        return SpinProbe(spin_ratio, balance, rates, world_part["tidal_heating"])


# ======================================================================================================================
# Roots of the spin balance
# ======================================================================================================================
def refine_root(rates, lower, upper, tolerance):
    """Shrinks a bracket below `tolerance` in spin ratio by the Illinois method, bisecting whenever two trials in a row
    fail to halve it; a trial whose tide solve fails is replaced by the midpoint once.

    Parameters
    ----------
    rates : PairRates
    lower, upper : SpinProbe
        The bracket: lower.spin_ratio < upper.spin_ratio and lower.balance > 0 > upper.balance.
    tolerance : float
        Final width [spin ratio].

    Returns
    -------
    tuple of SpinProbe
        The final (lower, upper), or one probe twice where the balance is exactly zero.
    """
    lower_value, upper_value = lower.balance, upper.balance
    retained, slow_trials = 0, 0
    while upper.spin_ratio - lower.spin_ratio > tolerance:
        width = upper.spin_ratio - lower.spin_ratio
        if slow_trials >= 2:
            trial = lower.spin_ratio + 0.5 * width
            slow_trials = 0
        else:
            trial = upper.spin_ratio - upper_value * width / (upper_value - lower_value)
            margin = 0.25 * tolerance
            trial = min(max(trial, lower.spin_ratio + margin), upper.spin_ratio - margin)
        try:
            probe = rates.probe(trial)
        except ProbeFailed:
            probe = rates.probe(lower.spin_ratio + 0.5 * width)
        if probe.balance > 0.0:
            lower, lower_value = probe, probe.balance
            if retained == +1:
                upper_value *= 0.5
            retained = +1
        elif probe.balance < 0.0:
            upper, upper_value = probe, probe.balance
            if retained == -1:
                lower_value *= 0.5
            retained = -1
        else:
            return probe, probe
        slow_trials = slow_trials + 1 if (upper.spin_ratio - lower.spin_ratio) > 0.5 * width else 0
    return lower, upper


def combine_sides(lower, upper):
    """The Filippov combination of a final bracket's two sides: the weights that zero the balance give the spin ratio,
    the rates, and the heating (a tuple)."""
    if lower is upper:
        return lower.spin_ratio, lower.rates, lower.heating
    weight = upper.balance / (upper.balance - lower.balance)
    return (weight * lower.spin_ratio + (1.0 - weight) * upper.spin_ratio,
            weight * lower.rates + (1.0 - weight) * upper.rates,
            weight * lower.heating + (1.0 - weight) * upper.heating)


def search_points(start, center, moving_up, finest):
    """Points from `start` toward the far side of `center` (up to HALF_SPACING past it), in the order a spin moving that
    way meets them, spaced by factors of 2 in their distance from `center` from `finest` out."""
    distances = finest * 2.0 ** np.arange(int(math.log2(HALF_SPACING / finest)) + 1)
    points = np.concatenate([
        center - distances,
        center + distances,
        [center + (HALF_SPACING if moving_up else -HALF_SPACING)]])
    if moving_up:
        return np.sort(points[points > start])
    return np.sort(points[points < start])[::-1]


def locate_root(rates, spin_ratio, tolerance, finest):
    """The stable root a spin at `spin_ratio` reaches at the current slow state: the first sign change of the balance in
    the direction it points, searched up to HALF_SPACING past the nearest commensurability on `search_points`. Returns
    the refined bracket, or None. A search point whose tide solve fails is skipped.

    Assumptions
    -----------
    * A pair of roots closer together than the search spacing near their location is not resolved.
    """
    start = rates.probe(spin_ratio)
    if start.balance == 0.0:
        return start, start
    moving_up = start.balance > 0.0
    previous = start
    for point in search_points(spin_ratio, 0.5 * round(2.0 * spin_ratio), moving_up, rates.finest_offset(finest)):
        try:
            probe = rates.probe(point)
        except ProbeFailed:
            continue
        if (probe.balance > 0.0) != moving_up:
            lower, upper = (previous, probe) if moving_up else (probe, previous)
            return refine_root(rates, lower, upper, tolerance)
        previous = probe
    return None


def build_window(rates, root_ratio, resolution, kept=(None, None)):
    """The window about a stable root (a SpinWindow), or None when either side has no reliably restoring point.

    The net torque near a root is a small difference of large mode torques, so the balance carries the radial solve's
    error amplified by about the quality factor. Its noise is measured from four probes NOISE_STEP apart at the root,
    and a balance within three times their spread has no reliable sign. Points step out from `resolution` past the
    root by factors of 4 (up to HALF_SPACING) until one has a reliably wrong sign; each edge is the restoring point one
    step inside the last one, keeping a margin from the sign change beyond. A rebuilt window passes the edge that did
    not fire in `kept` (lower, upper), reused while it still restores reliably. Points whose tide solve fails are
    skipped.

    Assumptions
    -----------
    * A root pair between two sample points is not detected.
    """
    samples = []
    for j in range(4):
        try:
            samples.append(rates.probe(root_ratio + j * NOISE_STEP).balance)
        except ProbeFailed:
            continue
    noise = 3.0 * (max(samples) - min(samples)) if len(samples) > 1 else 0.0
    edges = []
    for direction, kept_edge in zip((-1.0, +1.0), kept):
        if (kept_edge is not None) and (direction * (kept_edge - root_ratio) > resolution):
            try:
                balance = rates.probe(kept_edge).balance
                if (balance if direction < 0.0 else -balance) > noise:
                    edges.append(kept_edge)
                    continue
            except ProbeFailed:
                pass
        distance, restoring = resolution, []
        while distance <= HALF_SPACING:
            try:
                probe = rates.probe(root_ratio + direction * distance)
            except ProbeFailed:
                distance *= 4.0
                continue
            signed = probe.balance if direction < 0.0 else -probe.balance   # Positive when restoring
            if signed > noise:
                restoring.append(probe.spin_ratio)
            elif signed < -noise:
                break
            distance *= 4.0
        if not restoring:
            return None
        edges.append(restoring[-2] if len(restoring) > 1 else restoring[-1])
    return SpinWindow(edges[0], edges[1], noise)


# ======================================================================================================================
# The driver
# ======================================================================================================================
class TerminalEvent:
    """A CyRK terminal event, triggered when `function(t, y)` crosses zero downward."""

    def __init__(self, function):
        self.function = function
        self.terminal = 1
        self.direction = -1

    def __call__(self, t, y, *args):
        return self.function(t, np.asarray(y))


class SpinTracker:
    """Integrates the slow state, the spin free or tracked, over a span [Myr] in LSODA segments (see the module).

    Assumptions
    -----------
    * At most one stable root is reachable inside the band about a commensurability that a capture test searches.
    """

    def __init__(
            self,
            rates,
            spin_ratio,
            rtol,
            spin_rtol,
            atol,
            capture_band,
            capture_margin,
            root_tolerance,
            resolution,
            finest):
        self.rates = rates
        self.spin_ratio = spin_ratio
        self.rtol, self.spin_rtol, self.atol = rtol, spin_rtol, atol
        self.capture_band, self.capture_margin = capture_band, capture_margin
        self.root_tolerance, self.resolution, self.finest = root_tolerance, resolution, finest
        self.segments = []
        self.times, self.states, self.extras, self.tracked = [], [], [], []
        self.p_root = math.nan
        self.p_step = 4.0 * root_tolerance
        self.p_window = None
        # The spin ratio and the world's tidal heating of the last right-hand side call at each time of the current
        # segment. LSODA's corrector iterations of a step all run at the step's time, so the last call there is at its
        # accepted state to within the tolerances.
        self.p_recorded = {}

    # ==================================================================================================================
    # Right-hand sides and events
    # ==================================================================================================================
    def free_rates(self, t, y):
        """The scaled rates of a free segment's state y = [a / a0, e, s, host spin, temperatures]."""
        self.rates.set_slow_state(t, np.delete(y, 2))
        probe = self.rates.probe(y[2])
        self.p_recorded[t] = (y[2], probe.heating)
        return np.insert(probe.rates, 2, probe.balance)

    def tracked_root(self, t, x):
        """The tracked root's final bracket at slow state `x` (a tuple of SpinProbe).

        The search steps out from the last root by factors of 8, toward the side the balance points, as far as the
        window's edge, and refines the first sign change. A trial state past an edge event (an implicit step can reach
        well beyond one) can have an edge with the wrong sign; the search then follows the root past the window on a
        grid refined by factors of 2. With no root at all (past the state where the equilibrium vanishes) it returns
        the weakest balance found, flagged, which keeps the right-hand side defined; the edge event ends the step
        before that state.
        """
        rates = self.rates
        rates.set_slow_state(t, x)
        center = rates.probe(self.p_root)
        if center.balance == 0.0:
            return center, center
        moving_up = center.balance > 0.0
        limit = self.p_window.upper if moving_up else self.p_window.lower
        step, previous = self.p_step, center
        while True:
            point = self.p_root + step if moving_up else self.p_root - step
            at_edge = (point >= limit) if moving_up else (point <= limit)
            probe = rates.probe(limit if at_edge else point)
            if (probe.balance > 0.0) != moving_up:
                lower, upper = (previous, probe) if moving_up else (probe, previous)
                return refine_root(rates, lower, upper, self.root_tolerance)
            if at_edge:
                break
            previous = probe
            step *= 8.0
        previous, weakest = center, center
        for point in search_points(self.p_root, self.p_root, moving_up, rates.finest_offset(self.finest)):
            try:
                probe = rates.probe(point)
            except ProbeFailed:
                continue
            if (probe.balance > 0.0) != moving_up:
                lower, upper = (previous, probe) if moving_up else (probe, previous)
                return refine_root(rates, lower, upper, self.root_tolerance)
            if abs(probe.balance) < abs(weakest.balance):
                weakest = probe
            previous = probe
        weakest.weakest = True
        return weakest, weakest

    def tracked_rates(self, t, x):
        """The scaled rates of a tracked segment's slow state, at the Filippov combination of its root's sides."""
        lower, upper = self.tracked_root(t, x)
        spin_ratio, rates, heating = combine_sides(lower, upper)
        # Only a root inside the window becomes the next warm start, which keeps every search starting inside it.
        if (not lower.weakest) and (self.p_window.lower < spin_ratio < self.p_window.upper):
            # The next warm start brackets about as far as this root moved, within the tolerance and a hundredth.
            self.p_step = min(max(2.0 * abs(spin_ratio - self.p_root), 4.0 * self.root_tolerance), 1.0e-2)
            self.p_root = spin_ratio
        self.p_recorded[t] = (spin_ratio, heating)
        return rates

    def p_edge_balance(self, t, x, edge):
        self.rates.set_slow_state(t, x)
        return self.rates.probe(edge).balance

    def tracked_events(self):
        """The window's two edges losing their restoring signs."""
        window = self.p_window
        return [TerminalEvent(lambda t, x: self.p_edge_balance(t, x, window.lower)),
                TerminalEvent(lambda t, x: -self.p_edge_balance(t, x, window.upper))]

    def band_events(self, spin_ratio):
        """The spin entering the band of the next commensurability below or above, a band it starts in (or on the edge
        of) excluded, so both events are positive at the segment's start."""
        width = self.capture_band
        nearest = 0.5 * round(2.0 * spin_ratio)
        if abs(spin_ratio - nearest) <= width * (1.0 + 1.0e-6):
            below, above = nearest - 0.5, nearest + 0.5
        else:
            below, above = 0.5 * math.floor(2.0 * spin_ratio), 0.5 * math.ceil(2.0 * spin_ratio)
        return [TerminalEvent(lambda t, y: y[2] - (below + width)),
                TerminalEvent(lambda t, y: (above - width) - y[2])]

    def approach_events(self, spin_ratio, root_ratio):
        """The spin coming within capture_margin of `root_ratio`, or entering the next band either way (the root can
        move or vanish while the spin approaches it)."""
        margin, side = self.capture_margin, math.copysign(1.0, spin_ratio - root_ratio)
        return [TerminalEvent(lambda t, y: side * (y[2] - root_ratio) - margin)] + self.band_events(spin_ratio)

    # ==================================================================================================================
    # Transitions
    # ==================================================================================================================
    def p_decide(self, t, x, spin_ratio):
        """The next segment from a state: ("tracked", root bracket), ("approach", root ratio), or ("free", reason)."""
        self.rates.set_slow_state(t, x)
        # A band event stops the spin on the band's edge, which counts as inside.
        if abs(spin_ratio - 0.5 * round(2.0 * spin_ratio)) > self.capture_band * (1.0 + 1.0e-6):
            return "free", "outside the bands"
        bracket = locate_root(self.rates, spin_ratio, self.root_tolerance, self.finest)
        if bracket is None:
            return "free", "no stable root ahead"
        root_ratio = combine_sides(*bracket)[0]
        if abs(root_ratio - spin_ratio) <= self.capture_margin * (1.0 + 1.0e-6):   # An approach event ends on it
            return "tracked", bracket
        return "approach", root_ratio

    def p_run_segment(
            self,
            function,
            t_start,
            t_end,
            y0,
            events,
            first_step,
            free):
        """One LSODA segment. The wall-clock cap or a failed tide or EOS solve (a trial state the solvers reject) ends
        it, its stored steps kept. Returns (solver, error)."""
        solver = PySolver()
        rtol = np.asarray(self.rtol, dtype=float)
        if free:
            rtol = np.insert(rtol, 2, self.spin_rtol)
        self.p_recorded = {}
        try:
            pysolve_ivp(
                function,
                (t_start, t_end),
                y0,
                method="LSODA",
                rtol=rtol,
                atol=self.atol,
                events=events,
                first_step=first_step,
                solution_reuse=solver)
            error = None
        except (TimeoutError, RuntimeError, ValueError) as caught:
            error = caught
        return solver, error

    def p_end_state(
            self,
            mode,
            t,
            x,
            states,
            extras,
            error,
            solver):
        """The spin ratio at a segment's end, with the end's spin and heating written into `extras`, and for an edge
        event whose root is still near the old window the (root bracket, kept edges) of the window to rebuild (else
        None). A root farther than the window's width from it is left to the capture test."""
        if mode != "tracked":
            spin_ratio = float(states[2, -1])
            if math.isnan(extras[1, -1]):   # An event's end state is interpolated, never evaluated
                self.rates.set_slow_state(t, x)
                extras[1, -1] = self.rates.probe(spin_ratio).heating
            return spin_ratio, None
        # The last right-hand side call can be a trial state past the segment's end: the root is solved again at the end
        # state, from the last real root (kept when the end state has none).
        lower, upper = self.tracked_root(t, x)
        recenter = None
        if not lower.weakest:
            root_ratio = combine_sides(lower, upper)[0]
            window = self.p_window
            width = window.upper - window.lower
            if (window.lower - width) < root_ratio < (window.upper + width):
                self.p_root = root_ratio
                if (error is None) and solver.event_terminated:
                    fired_lower = solver.event_terminate_index == 0
                    recenter = ((lower, upper), (None, window.upper) if fired_lower else (window.lower, None))
        extras[:, -1] = (self.p_root, combine_sides(lower, upper)[2])
        return self.p_root, recenter

    def run(self, t_start, t_end):
        """Integrates from `t_start` to `t_end` [Myr]; returns (success, message)."""
        rates = self.rates
        t, x, spin_ratio = t_start, rates.initial_state(), self.spin_ratio
        first_step = 0.0   # LSODA's own choice, except after a failed evaluation
        recenter = None    # (root bracket, kept window edges) after an edge event with the root still there
        stalled, failed = 0, 0
        while t < t_end:
            if stalled >= MAX_STALLED_SEGMENTS:
                return False, f"TidalPy: System.evolve stalled: {stalled} segments in a row ended at t = {t:.10g} Myr."
            if failed >= MAX_FAILED_SEGMENTS:
                return False, (f"TidalPy: System.evolve stopped at t = {t:.10g} Myr: {failed} segments in a row ended "
                               f"on a failed evaluation.")
            try:
                if recenter is not None:
                    (mode, detail), kept, recenter = ("tracked", recenter[0]), recenter[1], None
                else:
                    (mode, detail), kept = self.p_decide(t, x, spin_ratio), (None, None)
                window = None
                if mode == "tracked":
                    window = build_window(rates, combine_sides(*detail)[0], self.resolution, kept)
                    if window is None:
                        mode, detail = "free", "no window about the root"
            except (TimeoutError, RuntimeError, ValueError) as error:   # The cap, or a solve failing at a reached state
                return False, str(error)
            entry = dict(time=t * TIME_UNIT, mode=mode, spin_ratio=float(spin_ratio))
            if mode == "tracked":
                self.p_root, self.p_window, self.p_step = combine_sides(*detail)[0], window, 4.0 * self.root_tolerance
                entry.update(root=self.p_root, window=(window.lower, window.upper))
                solver, error = self.p_run_segment(
                    self.tracked_rates,
                    t,
                    t_end,
                    x,
                    self.tracked_events(),
                    first_step,
                    free=False)
            else:
                if mode == "approach":
                    entry.update(target=detail)
                    events = self.approach_events(spin_ratio, detail)
                else:
                    entry.update(reason=detail)
                    events = self.band_events(spin_ratio)
                solver, error = self.p_run_segment(
                    self.free_rates,
                    t,
                    t_end,
                    np.insert(x, 2, spin_ratio),
                    events,
                    first_step,
                    free=True)
            self.segments.append(entry)
            log_debug(f"TidalPy: System.evolve segment at {t:.6g} Myr: {mode}, spin ratio {spin_ratio:.10f}, "
                      + ", ".join(f"{key} {value}" for key, value in entry.items() if key not in ("time", "mode")))

            # A segment that took no step failed at its own start, a state the run already reached: the run ends there.
            # Its arrays hold at most that start, so they are not read.
            if (solver.steps_taken == 0) or (solver.size == 0):
                return False, f"TidalPy: System.evolve took no step at t = {t:.10g} Myr ({error or solver.message})."
            times, states = np.asarray(solver.t), np.asarray(solver.y)
            slow = states if mode == "tracked" else np.delete(states, 2, axis=0)
            extras = np.array([self.p_recorded.get(time, (math.nan, math.nan)) for time in times]).T
            if mode != "tracked":
                extras[0] = states[2]   # A free spin is a state variable
            self.times.append(times)
            self.states.append(slow)
            self.extras.append(extras)
            self.tracked.append(np.full(times.size, mode == "tracked"))
            stalled = stalled + 1 if float(times[-1]) <= t else 0
            t, x = float(times[-1]), slow[:, -1].copy()
            try:
                spin_ratio, recenter = self.p_end_state(
                    mode,
                    t,
                    x,
                    states,
                    extras,
                    error,
                    solver)
            except (TimeoutError, RuntimeError, ValueError) as failure:
                return False, str(failure)
            if isinstance(error, TimeoutError):
                return False, str(error)
            if error is not None:   # Retest from the last stored step, LSODA restarting with a short first step
                entry["ended"] = str(error)
                log_info(f"TidalPy: System.evolve segment ended at {t:.6g} Myr: {error}")
                first_step = 1.0e-6 * max(t, 1.0e-3)
                failed += 1
                continue
            first_step, failed = 0.0, 0
            if solver.event_terminated:
                entry["ended"] = "event"
                continue
            if not solver.success:
                return False, f"TidalPy: System.evolve failed at t = {t:.10g} Myr: {solver.message}"
        return True, "TidalPy: System.evolve reached the end of its span."


# ======================================================================================================================
# Public entry point
# ======================================================================================================================
def join_segments(arrays, num_rows):
    """Segment arrays (1D when `num_rows` is 1, else rows by steps) joined along the steps, without the first step of
    every segment after the first, which repeats the step that ended the segment before."""
    if not arrays:
        return np.empty(0) if num_rows == 1 else np.empty((num_rows, 0))
    return np.concatenate([array[..., (i > 0):] for i, array in enumerate(arrays)], axis=-1)


def evolve_pair(
        system,
        world,
        time_span,
        *,
        evolve_thermal=True,
        surface_temperature=None,
        orbit_rtol=1.0e-9,
        thermal_rtol=1.0e-4,
        spin_rtol=1.0e-4,
        atol=1.0e-10,
        capture_band=1.0e-2,
        capture_margin=1.0e-3,
        root_tolerance=1.0e-10,
        resolution=1.0e-8,
        max_wall_time=None) -> EvolutionResult:
    """Evolves a world and its tidal host, the world's spin tracking its spin-orbit equilibria (see the module).

    The implementation behind :meth:`System.evolve <TidalPy.Structures.system.System.evolve>`, which documents the
    parameters. Returns an :class:`EvolutionResult` and leaves the system at the final state.

    Raises
    ------
    ValueError
        For the argument and state checks System.evolve lists.
    """
    if system.get_tidal_host_index(world) < 0:
        raise ValueError("TidalPy: System.evolve needs a world with a tidal host.")
    host = system.get_tidal_host(world)
    t_start, t_end = (float(value) for value in time_span)
    if not (math.isfinite(t_start) and math.isfinite(t_end) and (t_end > t_start)):
        raise ValueError("TidalPy: System.evolve needs a finite time span (t_start, t_end) with t_end > t_start [s].")
    orbital_frequency = system.calc_orbital_frequency(world)
    spin_ratio = world.spin_frequency / orbital_frequency
    if not (spin_ratio > 0.0):
        raise ValueError("TidalPy: System.evolve needs a world spinning prograde (spin_frequency > 0).")
    if evolve_thermal:
        temperatures = [layer.temperature for layer in world]
        if not temperatures:
            raise ValueError("TidalPy: System.evolve with evolve_thermal needs a world with layers.")
        if not all(math.isfinite(value) and (value > 0.0) for value in temperatures):
            raise ValueError("TidalPy: System.evolve with evolve_thermal needs every layer's temperature finite and "
                             "positive [K].")
    for name, value in (("capture_band", capture_band), ("capture_margin", capture_margin),
                        ("root_tolerance", root_tolerance), ("resolution", resolution)):
        if not (0.0 < value < HALF_SPACING):
            raise ValueError(f"TidalPy: System.evolve's {name} must be in (0, {HALF_SPACING}).")
    if not (root_tolerance < resolution < capture_margin <= capture_band):
        raise ValueError("TidalPy: System.evolve needs root_tolerance < resolution < capture_margin <= capture_band.")

    started = clock.perf_counter()
    deadline = None if max_wall_time is None else started + float(max_wall_time)
    rates = PairRates(
        system,
        world,
        host,
        evolve_thermal,
        surface_temperature,
        deadline)
    rtol = [orbit_rtol, orbit_rtol, orbit_rtol] + ([thermal_rtol] * rates.num_layers if evolve_thermal else [])
    tracker = SpinTracker(
        rates,
        spin_ratio,
        rtol,
        spin_rtol,
        atol,
        capture_band,
        capture_margin,
        root_tolerance,
        resolution,
        finest=root_tolerance)
    success, message = tracker.run(t_start / TIME_UNIT, t_end / TIME_UNIT)

    time = join_segments(tracker.times, 1) * TIME_UNIT
    slow = join_segments(tracker.states, len(rtol))
    extras = join_segments(tracker.extras, 2)
    tracked = join_segments(tracker.tracked, 1).astype(bool)
    semi_major_axis = slow[0] * rates.semi_major_axis0
    # The spin frequency follows the orbital motion of each step (Kepler's third law about the total mass).
    orbital_motion = orbital_frequency * (semi_major_axis / rates.semi_major_axis0) ** -1.5
    result = EvolutionResult(
        time=time,
        semi_major_axis=semi_major_axis,
        eccentricity=slow[1],
        spin_ratio=extras[0],
        spin_frequency=extras[0] * orbital_motion,
        host_spin_frequency=slow[2] * rates.host_spin0,
        temperature=(slow[3:] * rates.temperature0[:, np.newaxis]) if evolve_thermal else None,
        tidal_heating=extras[1],
        tracked=tracked,
        segments=tracker.segments,
        success=success,
        message=message,
        num_tide_solves=rates.num_probes)

    # Leave the system at the final state, the world's spin on its last ratio.
    if time.size > 0:
        rates.deadline = None
        try:
            rates.set_slow_state(time[-1] / TIME_UNIT, slow[:, -1].copy())
            world.set_spin_frequency(result.spin_ratio[-1] * rates.orbital_frequency)
        except (RuntimeError, ValueError) as error:
            result.success = False
            result.message += f" The system could not be set to the final state: {error}"
    result.num_eos_solves = rates.num_eos_solves
    result.elapsed = clock.perf_counter() - started
    log_info(f"TidalPy: System.evolve: {result.message} {time.size} steps, {len(tracker.segments)} segments, "
             f"{rates.num_probes} tide solves, {result.elapsed:.1f} s.")
    return result
