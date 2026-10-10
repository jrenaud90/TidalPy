"""Long forcing periods: the dynamic liquids' pressure form (derivatives/odes_.hpp) and the re-orthonormalized
shooting solutions (orthonormalize_.hpp).

The original-form values are TS72's y2 form at rtol 1e-12 / atol 1e-14, recorded before the dynamic liquids moved to
P = y2 - rho g y1 + rho y5 and before a layer's solutions were integrated together. The long-period references come from
an independent model of the same equations (Chebyshev fits of each world's solved EOS and moduli, integrated with QR
re-orthonormalization after every step, matching these worlds' short-period k2 to 5e-10).
"""
import math

import numpy as np
import pytest

import TidalPy
from TidalPy.constants import G, update_constants
from TidalPy.Material import Material
from TidalPy.RadialSolver import build_rs_input_homogeneous_layers
from TidalPy.RadialSolver.solver import radial_solver
from TidalPy.Structures import build_world
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.base import BaseWorld

DAY = 86400.0
TIGHT = dict(rtol=1.0e-12, atol=1.0e-14)

# label: (world, period [days], liquid (is_static, is_incompressible) or None for the world's own, (k, h, l)).
ORIGINAL_FORM = {
    "europa_dynamic": ("europa_dynamic", 3.551, None, (
        2.627784064331e-01 - 2.418464849288e-03j, 1.214064335183e+00 - 1.094510992851e-02j,
        3.044668197581e-01 - 6.981840824166e-03j)),
    "europa_incompressible_ocean": ("europa_dynamic", 3.551, (False, True), (
        2.622863168331e-01 - 2.440617835474e-03j, 1.219995557816e+00 - 1.098519390893e-02j,
        3.056668616102e-01 - 6.964699737671e-03j)),
    "luna_dynamic": ("luna_dynamic", 1.0, None, (
        2.390590250415e-02 - 5.706667422348e-05j, 4.149691028257e-02 - 9.852015719327e-05j,
        1.099480900030e-02 - 1.878176403927e-05j)),
    "mercury_dynamic": ("mercury_dynamic", 87.9691, None, (
        5.689999631269e-01 - 4.725872605770e-03j, 1.018746557624e+00 - 8.240828701902e-03j,
        2.019458305361e-01 - 2.085373865180e-03j)),
    "pluto_dynamic": ("pluto_dynamic", 6.3872, None, (
        3.202609308434e-01 - 7.520191778652e-02j, 9.710086137672e-01 - 2.280788057424e-01j,
        2.325519517818e-01 - 6.133394945598e-02j)),
    "earth_simple_incompressible_core": ("earth_simple", 0.5, (False, True), (
        3.004764482866e-01 - 1.209089609990e-03j, 6.312910423482e-01 - 2.220886494122e-03j,
        7.540986967100e-02 - 5.294090649346e-04j)),
    # Solids and static liquids: no pressure form, only the solutions integrated together.
    "io": ("io", 1.769, None, (
        3.571305182280e-02 - 1.499696390563e-02j, 6.298322955852e-02 - 2.720676614890e-02j,
        5.709801946622e-02 - 1.263155588023e-01j)),
    "earth_simple": ("earth_simple", 0.5, None, (
        2.999390132507e-01 - 1.205405399491e-03j, 6.301492455065e-01 - 2.213961542832e-03j,
        7.542227914550e-02 - 5.293780003045e-04j)),
    "earth_prem": ("earth_prem", 0.5, None, (
        2.983639333614e-01 + 0.0j, 6.044228970077e-01 + 0.0j, 8.385741288555e-02 + 0.0j)),
}
SOLID_AND_STATIC = ("io", "earth_simple", "earth_prem")

# (world, period [days]): k2 of the independent model.
INDEPENDENT_LONG_PERIOD = {
    ("europa_dynamic", 1.0e4): 0.2758323838 - 0.0009730946j,
    ("luna_dynamic", 1.0e3): 0.0262045293 - 0.0004664797j,
}


def _world(name, liquid=None, **eos_kwargs):
    world = build_world(name)
    if liquid is not None:
        for layer in world.layers:
            if layer.is_liquid:
                layer.is_static, layer.is_incompressible = liquid
    world.solve_eos(**eos_kwargs)
    return world


def _frequency(period_days):
    return 2.0 * math.pi / (period_days * DAY)


def _love(world):
    return np.asarray([complex(world.love_number_k), complex(world.love_number_h), complex(world.love_number_l)])


def _relative(value, reference):
    return np.max(np.abs(np.asarray(value) - np.asarray(reference)) / np.abs(np.asarray(reference)))


@pytest.mark.parametrize("label", tuple(ORIGINAL_FORM))
def test_short_periods_match_the_original_form(label):
    """At short periods, where TS72's y2 form is well conditioned, the pressure form and the joint integration give
    its Love numbers to the integration tolerance."""
    name, period, liquid, reference = ORIGINAL_FORM[label]
    world = _world(name, liquid)
    world.solve_love_numbers(_frequency(period), 2, warnings=False, **TIGHT)
    assert _relative(_love(world), reference) < 1.0e-8


@pytest.mark.parametrize("label", SOLID_AND_STATIC)
def test_solid_and_static_liquid_worlds_are_unchanged_at_default_tolerances(label):
    name, period, liquid, reference = ORIGINAL_FORM[label]
    world = _world(name, liquid)
    world.solve_love_numbers(_frequency(period), 2, warnings=False)
    assert _relative(_love(world), reference) < 1.0e-6


@pytest.mark.parametrize("name, period", tuple(INDEPENDENT_LONG_PERIOD))
def test_long_periods_match_an_independent_model(name, period, spdlog_text):
    """Default tolerances reproduce the independent model at periods where TS72's y2 form is off by 3e-5 (Europa,
    10^4 days) and 1.5e-2 (Luna, 1000 days), and the surface solve stays well conditioned."""
    world = _world(name)
    world.solve_love_numbers(_frequency(period), 2)
    assert _relative(complex(world.love_number_k), INDEPENDENT_LONG_PERIOD[(name, period)]) < 1.0e-8
    assert world.love_surface_rcond > 1.0e-3
    assert "poorly conditioned" not in spdlog_text()


@pytest.mark.parametrize("name, period", [
    ("europa_dynamic", 1.0e3), ("europa_dynamic", 1.0e4), ("luna_dynamic", 1.0e2), ("luna_dynamic", 1.0e3),
    ("mercury_dynamic", 1.0e3), ("pluto_dynamic", 1.0e4)])
def test_long_periods_converge_in_few_steps(name, period):
    """Default and tight tolerances agree in a few hundred steps at most, where TS72's y2 form takes 6000 to 44000
    per Europa ocean solve at 10^3 to 10^4 days."""
    world = _world(name)
    world.solve_love_numbers(_frequency(period), 2, warnings=False)
    default = _love(world)
    steps = world.release_radial_solution().steps_taken.max(axis=1)
    world.solve_love_numbers(_frequency(period), 2, warnings=False, **TIGHT)
    assert _relative(default, _love(world)) < 1.0e-7
    assert int(np.sum(steps)) < 500


@pytest.mark.parametrize("frequency", (_frequency(1.0e4), 1.0e-10, 1.0e-14))
def test_dense_solution_is_continuous_across_restarts(frequency, restore_config):
    """The radial functions through Europa's ocean, rebuilt segment by segment from the carried basis changes, agree
    between a solve forced to re-orthonormalize often and a tight one, with no step at a restart. y3, rebuilt as
    -P / (rho omega^2 r), holds its digits at long periods too."""
    world = _world("europa_dynamic")
    ocean_i = [layer.name for layer in world.layers].index("ocean")
    ocean = world.layers[ocean_i]
    # Inside the ocean: y3 steps across its solid-liquid interfaces, which belong to the layers below them.
    span = ocean.radius_outer - ocean.radius_inner
    radii = np.linspace(ocean.radius_inner + 1.0e-6 * span, ocean.radius_outer - 1.0e-6 * span, 1501)
    world.solve_love_numbers(frequency, 2, warnings=False, **TIGHT)
    tight = np.asarray([world.get_love_radial_y(radii, 0, y_i) for y_i in (0, 1, 2, 4, 5)])
    TidalPy.config["numerical"]["minimum_solution_independence"] = 0.9
    update_constants()
    world.solve_love_numbers(frequency, 2, warnings=False, rtol=1.0e-10, atol=1.0e-10)
    assert int(world.release_radial_solution().orthonormalizations[ocean_i]) > 1
    world.solve_love_numbers(frequency, 2, warnings=False, rtol=1.0e-10, atol=1.0e-10)
    restarted = np.asarray([world.get_love_radial_y(radii, 0, y_i) for y_i in (0, 1, 2, 4, 5)])
    assert np.all(np.isfinite(restarted))
    scale = np.max(np.abs(tight), axis=1, keepdims=True)
    assert np.max(np.abs(restarted - tight) / scale) < 1.0e-8
    # No step between neighboring radii larger than the smooth profile allows.
    assert np.max(np.abs(np.diff(restarted, axis=1)) / scale) < 1.0e-2


@pytest.mark.parametrize("degree_l, solve_for", [(1, ("loading",)), (2, ("tidal", "loading")), (3, ("tidal",))])
def test_other_degrees_and_boundary_conditions_share_the_path(degree_l, solve_for):
    """Degree-1 loading (whose frame fixes k' = 0), loading, and higher degrees through a dynamic ocean at 1000 days."""
    world = _world("europa_dynamic")
    results = []
    for kwargs in ({}, TIGHT):
        world.solve_love_numbers(_frequency(1.0e3), degree_l, solve_for=solve_for, warnings=False, **kwargs)
        results.append([np.atleast_1d(np.asarray(getattr(world, f"love_number_{x}"), dtype=complex)) for x in "khl"])
    (k, h, l), (k_tight, h_tight, l_tight) = results
    assert _relative(h, h_tight) < 1.0e-7
    assert _relative(l, l_tight) < 1.0e-7
    if degree_l == 1:
        assert np.max(np.abs(k)) < 1.0e-10
    else:
        assert _relative(k, k_tight) < 1.0e-7


def test_an_unstably_stratified_liquid_solves():
    """A constant-density compressible liquid core (N^2 = -rho g^2 / K < 0) solved dynamically at 100 days agrees with
    a tight solve and lies close to the static core. TS72's y2 form gives k2 = -0.40 for it at half a day."""
    world = _world("earth_simple", liquid=(False, False))
    world.solve_love_numbers(_frequency(100.0), 2, warnings=False)
    dynamic = complex(world.love_number_k)
    world.solve_love_numbers(_frequency(100.0), 2, warnings=False, **TIGHT)
    assert _relative(dynamic, complex(world.love_number_k)) < 1.0e-7
    static = _world("earth_simple")
    static.solve_love_numbers(_frequency(100.0), 2, warnings=False)
    assert _relative(dynamic, complex(static.love_number_k)) < 1.0e-4


def test_a_thermal_dynamic_core_keeps_a_positive_q():
    """Luna-Dynamic with its temperature solved. TS72's y2 form gives it a negative monthly Q (k2 = 0.0211 + 3.4e-4 i)
    at the default tolerances while reporting success."""
    world = _world("luna_dynamic", solve_temperature=True)
    world.solve_love_numbers(_frequency(27.3217), 2, warnings=False)
    k2 = complex(world.love_number_k)
    world.solve_love_numbers(_frequency(27.3217), 2, warnings=False, **TIGHT)
    assert _relative(k2, complex(world.love_number_k)) < 1.0e-7
    assert -k2.imag > 0.0


def test_near_static_solids_stay_conditioned():
    """Earth-Simple at 1e-16 rad/s, where its viscous mantle is nearly fluid. Without re-orthonormalization the surface
    rcond is 1e-6 and k2 is off by 1e-2 with the wrong sign of its imaginary part."""
    world = _world("earth_simple")
    world.solve_love_numbers(1.0e-16, 2, warnings=False)
    k2 = complex(world.love_number_k)
    assert world.love_surface_rcond > 1.0e-3
    world.solve_love_numbers(1.0e-16, 2, warnings=False, **TIGHT)
    reference = complex(world.love_number_k)
    assert abs(k2.real - reference.real) < 1.0e-8 * abs(reference.real)
    assert abs(k2.imag - reference.imag) < 1.0e-6 * abs(reference.imag)


def test_si_and_non_dimensional_solves_restart_alike():
    """The independence measure scales each radial function by its characteristic size, so Earth-Simple near zero
    frequency re-orthonormalizes about as often in SI units as in non-dimensional ones, to the same k2."""
    world = _world("earth_simple")
    counts, love = [], []
    for nondimensionalize in (True, False):
        world.solve_love_numbers(1.0e-16, 2, warnings=False, nondimensionalize=nondimensionalize)
        love.append(complex(world.love_number_k))
        counts.append(int(np.sum(world.release_radial_solution().orthonormalizations)))
    assert 0 < counts[1] <= 2 * counts[0]
    assert _relative(love[1], love[0]) < 1.0e-7


def test_minimum_solution_independence_is_wired_through(restore_config):
    """The configured value reaches the constants module, which C++ reads."""
    from TidalPy import constants
    assert constants.minimum_solution_independence == TidalPy.config["numerical"]["minimum_solution_independence"]
    TidalPy.reinit(provided_config={"numerical": {"minimum_solution_independence": 1.0e-3}})
    assert constants.minimum_solution_independence == 1.0e-3


@pytest.mark.parametrize("value", (-1.0, 1.0, 1.5))
def test_minimum_solution_independence_is_a_fraction(value, restore_config):
    """A value outside [0, 1) is refused by the configuration and by update_constants."""
    message = "minimum_solution_independence must be at least 0 and below 1"
    with pytest.raises(ValueError, match=message):
        TidalPy.reinit(provided_config={"numerical": {"minimum_solution_independence": value}})
    TidalPy.config["numerical"]["minimum_solution_independence"] = value
    with pytest.raises(ValueError, match=message):
        update_constants()


def _unstable_ocean(frequency):
    """A solid planet under a 190 km constant-density compressible ocean (N^2 < 0), every layer dynamic."""
    return build_rs_input_homogeneous_layers(
        6.371e6, frequency, (5500.0, 3300.0, 1000.0), (1.5e11, 1.2e11, 2.2e9), (6.0e10, 7.0e10, 0.0),
        (1.0e30,) * 3, (1.0e30,) * 3, ("solid", "solid", "liquid"), (False, False, False), (False, False, False),
        "elastic", "elastic", radius_fraction_tuple=(0.5, 0.97, 1.0), slice_per_layer=40)


def test_max_num_steps_holds_over_a_layer():
    """The step cap holds for a layer's re-orthonormalized segments together, and a layer that reaches it fails with
    its steps recorded."""
    full = radial_solver(*_unstable_ocean(1.0e-7), warnings=False)
    assert full.success, full.message
    ocean_steps = int(np.max(full.steps_taken[2]))
    assert (ocean_steps > 1000) and (int(full.orthonormalizations[2]) > 10)
    capped = radial_solver(*_unstable_ocean(1.0e-7), warnings=False, max_num_steps=1000)
    assert not capped.success
    assert 0 < int(np.max(capped.steps_taken[2])) <= 1000


def test_a_degree_one_load_on_an_unstable_ocean_reaches_the_static_limit():
    """At 1e-9 rad/s an unstably stratified dynamic surface ocean grows by about exp(1e5) through its 190 km, and its
    degree-1 load h' approaches the static ocean's 1 - rho_bar / rho_s at every tolerance."""
    static_limit = 1.0 - (5500.0 * 0.125 + 3300.0 * (0.97**3 - 0.125) + 1000.0 * (1.0 - 0.97**3)) / 1000.0
    for rtol in (3.0e-8, 1.0e-9):
        # Its 200,000 to 400,000 steps over 30,000 re-orthonormalized segments keep about 500 MB of dense solution.
        solution = radial_solver(*_unstable_ocean(1.0e-9), degree_l=1, solve_for=("loading",), warnings=False,
                                 integration_rtol=rtol, integration_atol=rtol, max_ram_MB=2000,
                                 max_num_steps=2_000_000)
        assert solution.success, solution.message
        assert abs(complex(solution.h).real - static_limit) < 1.0e-4 * abs(static_limit)


def test_calc_tides_at_a_long_period_dynamic_ocean():
    """Eccentricity tides on Europa-Dynamic in a synchronous 10^4-day orbit take their Love number from the same
    long-period radial solve, so the heating is -(21/2) Im(k2) (n R)^5 e^2 / G."""
    world = _world("europa_dynamic")
    mean_motion, eccentricity, host_mass = _frequency(1.0e4), 1.0e-3, 1.898e27
    semi_major_axis = (G * host_mass / mean_motion**2) ** (1.0 / 3.0)
    world.set_spin_frequency(mean_motion)
    world.calc_tides(mean_motion, mean_motion, eccentricity, 0.0, semi_major_axis, host_mass)
    heating = world.get_tidal_heating()
    world.solve_love_numbers(mean_motion, 2, warnings=False)
    k2 = complex(world.love_number_k)
    expected = -10.5 * k2.imag * (mean_motion * world.radius) ** 5 * eccentricity**2 / G
    assert heating > 0.0
    assert heating == pytest.approx(expected, rel=1.0e-4)


def test_pluto_dynamic_solves_at_very_low_frequencies():
    """pluto_dynamic's ocean solves at 1e-13 rad/s, where default and tight tolerances agree."""
    world = _world("pluto_dynamic")
    world.solve_love_numbers(1.0e-13, 2, warnings=False)
    default = _love(world)
    world.solve_love_numbers(1.0e-13, 2, warnings=False, **TIGHT)
    assert _relative(default, _love(world)) < 1.0e-7


# k2 of Europa-Dynamic at a 3.551 day period with the ocean's supplied bulk modulus scaled, from TS72's y2 form at
# rtol 1e-13 (which reads no density gradient, so it needs no buoyancy frequency against the supplied modulus).
SUPPLIED_OCEAN_BULK = {0.5: 2.638151844339e-01 - 2.419553473276e-03j, 2.0: 2.610968649382e-01 - 2.372567339167e-03j}


@pytest.mark.parametrize("factor", tuple(SUPPLIED_OCEAN_BULK))
def test_a_supplied_liquid_bulk_modulus_sets_the_stratification(factor):
    """A dynamic liquid given a bulk modulus other than its EOS's is stratified against the supplied one: the buoyancy
    frequency the solve reads moves with the modulus."""
    frequency = _frequency(3.551)
    slices = 200
    world = build_world("europa_dynamic")
    radius = np.ascontiguousarray(world.solve_eos(slices_per_layer=slices)["radius"], dtype=np.float64)
    shear = np.empty(radius.size, dtype=np.complex128)
    bulk = np.empty(radius.size, dtype=np.complex128)
    for layer_i, layer in enumerate(world.layers):
        part = slice(layer_i * slices, (layer_i + 1) * slices)
        shear[part] = layer.calc_complex_shear_modulus(radius[part], frequency)
        bulk[part] = layer.calc_complex_bulk_modulus(radius[part], frequency) * (factor if layer.is_liquid else 1.0)
    result = world.solve_love_numbers_supplied(
        shear, bulk, radius.copy(), frequency=frequency, rtol=1.0e-10, atol=1.0e-12)
    assert result["success"], result["message"]
    assert _relative(complex(result["love_number_k"]), SUPPLIED_OCEAN_BULK[factor]) < 2.0e-8


@pytest.mark.parametrize("name", ("luna", "mercury"))
@pytest.mark.parametrize("frequency", (1.0e-12, 1.0e-16))
@pytest.mark.parametrize("nondimensionalize", (True, False))
def test_constant_density_incompressible_liquids_reach_the_static_limit(name, frequency, nondimensionalize):
    """A constant-density incompressible liquid is neutral, and its density gradient is exactly 0, so solved dynamically
    it reaches the static liquid's Love numbers at the lowest frequencies, in either units."""
    dynamic = _world(name, liquid=(False, True))
    dynamic.solve_love_numbers(frequency, 2, warnings=False, nondimensionalize=nondimensionalize)
    static = _world(name)
    static.solve_love_numbers(frequency, 2, warnings=False, **TIGHT)
    assert _relative(_love(dynamic)[:2], _love(static)[:2]) < 1.0e-8


@pytest.mark.parametrize("nondimensionalize", (True, False))
@pytest.mark.parametrize("frequency", (1.0e-3, 1.0e-5))
def test_a_tabulated_dynamic_core(frequency, nondimensionalize):
    """PREM's outer core, a density tabulated in radius, solved dynamically: its density gradient holds at the core's
    own end radii, so default and tight tolerances agree in either units."""
    world = build_world("earth_prem")
    for layer in world.layers:
        if layer.is_liquid:
            layer.is_static = False
    world.solve_eos()
    world.solve_love_numbers(frequency, 2, warnings=False, nondimensionalize=nondimensionalize)
    default = complex(world.love_number_k)
    world.solve_love_numbers(frequency, 2, warnings=False, nondimensionalize=nondimensionalize, **TIGHT)
    assert _relative(default, complex(world.love_number_k)) < 1.0e-6


def _tabulated_ocean(frequency, profile, slices=60):
    """_unstable_ocean's three layers with the ocean's density tabulated to vary linearly through it."""
    inputs = list(build_rs_input_homogeneous_layers(
        6.371e6, frequency, (5500.0, 3300.0, 1000.0), (1.5e11, 1.2e11, 2.2e9), (6.0e10, 7.0e10, 0.0),
        (1.0e30,) * 3, (1.0e30,) * 3, ("solid", "solid", "liquid"), (False, False, False), (False, False, False),
        "elastic", "elastic", radius_fraction_tuple=(0.5, 0.97, 1.0), slice_per_layer=slices))
    radius, density = inputs[0], inputs[1].copy()
    fraction = (radius[-slices:] - radius[-slices]) / (radius[-1] - radius[-slices])
    density[-slices:] = profile(fraction)
    inputs[1] = np.ascontiguousarray(density)
    return inputs


@pytest.mark.parametrize("profile", (lambda x: 1100.0 - 100.0 * x, lambda x: 1000.0 + 50.0 * x),
                         ids=("stable", "unstable"))
def test_a_supplied_tabulated_liquid_solves_in_either_units(profile):
    """A dynamic ocean whose supplied density varies through it solves from its own base radius on, to the same Love
    number in SI and non-dimensional units."""
    love = []
    for nondimensionalize in (True, False):
        solution = radial_solver(
            *_tabulated_ocean(1.0e-3, profile), warnings=False, nondimensionalize=nondimensionalize)
        assert solution.success, solution.message
        love.append(complex(solution.k))
    assert _relative(love[1], love[0]) < 1.0e-7


def test_turning_re_orthonormalization_off(restore_config):
    """[numerical] minimum_solution_independence = 0 integrates the solutions together without restarts, and the
    surface solve of a near-static solid loses its conditioning again."""
    world = _world("earth_simple")
    world.solve_love_numbers(1.0e-16, 2, warnings=False)
    assert world.love_surface_rcond > 1.0e-3
    assert int(np.sum(world.release_radial_solution().orthonormalizations)) > 0
    TidalPy.config["numerical"]["minimum_solution_independence"] = 0.0
    update_constants()
    world.solve_love_numbers(1.0e-16, 2, warnings=False)
    assert world.love_surface_rcond < 1.0e-5
    assert int(np.sum(world.release_radial_solution().orthonormalizations)) == 0


@pytest.mark.parametrize("name", ("europa_dynamic", "luna_dynamic", "mercury_dynamic", "pluto_dynamic"))
@pytest.mark.parametrize("frequency", (1.0e-13, 1.0e-15, 1.0e-16))
@pytest.mark.parametrize("nondimensionalize", (True, False))
def test_neutral_liquids_reach_the_static_limit(name, frequency, nondimensionalize):
    """Far below every frequency of the body, a neutral dynamic liquid's k and h are a static liquid's, to the
    integration tolerance, in non-dimensional and SI units alike. (Pluto's l, set by its nearly fluid ice shell, moves
    by 1e-6 between tolerances in either form.)"""
    dynamic = _world(name)
    dynamic.solve_love_numbers(frequency, 2, warnings=False, nondimensionalize=nondimensionalize)
    static = build_world(name)
    for layer in static.layers:
        if layer.is_liquid or layer.can_change_state:
            layer.is_static = True
    static.solve_eos()
    static.solve_love_numbers(frequency, 2, warnings=False, **TIGHT)
    assert _relative(_love(dynamic)[:2], _love(static)[:2]) < 1.0e-8


def _buoyancy_difference(solution, radii, step):
    """The EOS's N^2 at the radii, N^2 from a centered difference of its density, and the size of g rho' / rho. The
    EOS's own density gradient matches the difference too."""
    state = solution.eos_call(radii)
    above, below = solution.eos_call(radii + step)["density"], solution.eos_call(radii - step)["density"]
    density, gravity = state["density"], state["gravity"]
    gradient = (above - below) / (2.0 * step)
    assert np.max(np.abs(state["density_gradient"] - gradient)) <= 1.0e-6 * np.max(np.abs(gradient)) + 1.0e-12
    gradient_term = gravity * gradient / density
    difference = -(gradient_term + gravity * gravity * density / state["bulk_modulus"])
    return state["buoyancy_frequency_squared"], difference, np.max(np.abs(gradient_term))


@pytest.mark.parametrize("name, eos_kwargs, layer_name", [
    ("europa_dynamic", {}, "ocean"), ("luna_dynamic", dict(solve_temperature=True), "outer_core"),
    ("earth_prem", {}, "layer_1")])
def test_the_buoyancy_frequency_follows_the_density(name, eos_kwargs, layer_name):
    """The N^2 = -g (rho' / rho + rho g / K) the pressure form reads matches one from a centered difference of the
    reported density: a Birch-Murnaghan ocean, a thermal Birch-Murnaghan core, and PREM's tabulated outer core."""
    world = _world(name, **eos_kwargs)
    world.solve_love_numbers(1.0e-5, 2, warnings=False)
    solution = world.release_radial_solution()
    layer = world.layers[[layer.name for layer in world.layers].index(layer_name)]
    span = layer.radius_outer - layer.radius_inner
    radii = np.linspace(layer.radius_inner + 0.05 * span, layer.radius_outer - 0.05 * span, 21)
    buoyancy, difference, scale = _buoyancy_difference(solution, radii, 1.0e-4 * span)
    assert np.max(np.abs(buoyancy - difference)) < 1.0e-6 * scale


def test_a_neutral_liquid_has_no_buoyancy():
    """Europa-Dynamic's ocean has a density that follows its bulk modulus, so its N^2 is exactly 0, not the roundoff
    of rho' against rho^2 g / K that would set its response at very long periods."""
    world = _world("europa_dynamic")
    world.solve_love_numbers(1.0e-5, 2, warnings=False)
    solution = world.release_radial_solution()
    ocean = world.layers[[layer.name for layer in world.layers].index("ocean")]
    radii = np.linspace(ocean.radius_inner, ocean.radius_outer, 41)[1:-1]
    assert np.all(solution.eos_call(radii)["buoyancy_frequency_squared"] == 0.0)


@pytest.mark.parametrize("solve_temperature", [False, True])
def test_the_buoyancy_frequency_through_a_melting_range(solve_temperature):
    """Inside a melting range that follows the pressure, with a lighter liquid's density mixed in, N^2 carries the melt
    fraction's change with pressure and temperature, under a uniform or a solved (adiabatic) temperature."""
    radius, density = 2.0e6, 4000.0
    thermal = {"thermal_conductivity_w_mk": 3.0, "heat_capacity_j_kgk": 1000.0}

    def curve(temperature):
        return {"model": "simon_glatzel", "temperature_k": temperature, "simon_a_pa": 2.0e10, "simon_c": 1.0}

    material = Material(config={
        "solid": {**thermal, "shear_modulus": {"model": "constant", "shear_modulus_pa": 5.0e10},
                  "eos": {"model": "birch_murnaghan", "reference_density_kg_m3": density,
                          "reference_bulk_modulus_pa": 1.0e11, "bulk_modulus_derivative": 4.0},
                  "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e20}},
        "liquid": {**thermal, "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0},
                   "eos": {"model": "birch_murnaghan", "reference_density_kg_m3": 3700.0,
                           "reference_bulk_modulus_pa": 3.0e10, "bulk_modulus_derivative": 5.0}},
        "melting": {"solidus": curve(1500.0), "liquidus": curve(1700.0)},
        "latent_heat_j_kg": 4.0e5})
    world = BaseWorld("melting", radius, 4.0 / 3.0 * math.pi * density * radius**3)
    world.add_layer(Layer("body", 0, 0.0, radius, material=material, temperature=1690.0, use_melting=True,
                          use_pressure_melting=True, use_melt_density=True))
    assert world.solve_eos(G_to_use=G, solve_temperature=solve_temperature)["success"]
    world.solve_love_numbers(1.0e-5, 2, warnings=False)
    solution = world.release_radial_solution()
    radii = np.linspace(0.05 * radius, 0.95 * radius, 41)
    buoyancy, difference, scale = _buoyancy_difference(solution, radii, 1.0e-5 * radius)
    assert np.nanmax(solution.eos_call(radii)["melt_fraction"]) > 0.3
    assert np.max(np.abs(buoyancy - difference)) < 1.0e-6 * scale
