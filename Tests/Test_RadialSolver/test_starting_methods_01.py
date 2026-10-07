"""Starting methods beyond the published closed forms: the power series (Martens 2016) spans the closed-form solutions
and reproduces Martens' homogeneous-Earth Love numbers; its compressible-solid solutions are the series of the Takeuchi
and Saito wavenumber solutions; the derived Takeuchi incompressible and Kamata static incompressible forms span the
series' solutions; every method starts a static incompressible solid; a weak solid start layer gives the closed-form
Love numbers or a refusal; and the configuration's starting_method reaches the solver."""
import numpy as np
import pytest

import TidalPy
from TidalPy.exceptions import SolutionFailedError
from TidalPy.RadialSolver.derivatives.odes import find_num_shooting_solutions
from TidalPy.RadialSolver.solver import radial_solver
from TidalPy.RadialSolver.starting.kamata import (
    kamata_liquid_dynamic_incompressible,
    kamata_solid_dynamic_incompressible,
    kamata_solid_static_incompressible,
)
from TidalPy.RadialSolver.starting.power_series import (
    power_series_liquid_dynamic_compressible,
    power_series_liquid_dynamic_incompressible,
    power_series_solid_dynamic_compressible,
    power_series_solid_dynamic_incompressible,
    power_series_solid_static_compressible,
    power_series_solid_static_incompressible,
)
from TidalPy.RadialSolver.starting.takeuchi import (
    takeuchi_liquid_dynamic_compressible,
    takeuchi_liquid_dynamic_incompressible,
    takeuchi_solid_dynamic_compressible,
    takeuchi_solid_dynamic_incompressible,
    takeuchi_solid_static_compressible,
    takeuchi_solid_static_incompressible,
)
from TidalPy.Rheology import Maxwell

from starting_methods import STARTING_METHODS

# Martens (2016) Table C.4: a homogeneous Earth with VP = 10 km/s, VS = 5 km/s, and rho = 5 g/cc, at the M2 period.
# n: (h, n l, n k, h', n l', n k'). The table used G = 6.672e-11 (LoadDef's default); TidalPy's SciPy G is 3.5e-4
# larger, which moves these values by up to about 2e-4.
TABLE_C4 = {
    2: (0.5079, 0.2729, 0.5845, -0.6000, 0.0026, -0.4313),
    3: (0.3159, 0.1157, 0.3881, -0.6886, 0.0542, -0.5597),
    5: (0.1935, 0.0432, 0.2537, -0.7846, 0.1133, -0.7141),
    10: (0.1025, 0.0116, 0.1427, -0.8799, 0.1730, -0.8827),
}
TABLE_C4_ATOL = 5.0e-4
EARTH_RADIUS = 6371.0e3
M2_FREQUENCY = 2.0 * np.pi / (12.42 * 3600.0)
HOMOGENEOUS_DENSITY = 5000.0
HOMOGENEOUS_SHEAR = HOMOGENEOUS_DENSITY * 5.0e3**2
HOMOGENEOUS_BULK = HOMOGENEOUS_DENSITY * 10.0e3**2 - 4.0 * HOMOGENEOUS_SHEAR / 3.0
NUM_SLICES = 100


def _subspace_sine(first, second):
    """Sine of the largest principal angle between the row spaces of two (num_solutions, num_ys) arrays, after each
    component is scaled by its largest magnitude in either (a change of units, which leaves equal spaces equal)."""
    scale = np.maximum(np.max(np.abs(first), axis=0), np.max(np.abs(second), axis=0))
    scale[scale == 0.0] = 1.0
    q_first, _ = np.linalg.qr((first / scale).T)
    q_second, _ = np.linalg.qr((second / scale).T)
    return np.linalg.norm(q_second - q_first @ (q_first.conj().T @ q_second), 2)


def _empty(layer_type, is_static, is_incompressible):
    num_sols = find_num_shooting_solutions(layer_type, is_static, is_incompressible)
    return np.zeros((num_sols, 2 * num_sols), dtype=np.complex128, order="C")


# (series wrapper, closed-form wrapper, arguments before the view, (layer_type, is_static, is_incompressible)).
SOLID = (5000.0, 3.3e11 + 1.0e9j, 1.25e11 + 2.0e9j)
LIQUID = (10000.0, 1.0e11 + 0.0j)
DAY = 2.0 * np.pi / 86400.0
SPAN_CASES = (
    (power_series_solid_dynamic_compressible, takeuchi_solid_dynamic_compressible,
     lambda r, l: (M2_FREQUENCY, r, *SOLID, l), (0, False, False)),
    (power_series_solid_static_compressible, takeuchi_solid_static_compressible,
     lambda r, l: (r, *SOLID, l), (0, True, False)),
    (power_series_solid_dynamic_incompressible, kamata_solid_dynamic_incompressible,
     lambda r, l: (M2_FREQUENCY, r, SOLID[0], SOLID[2], l), (0, False, True)),
    (power_series_liquid_dynamic_compressible, takeuchi_liquid_dynamic_compressible,
     lambda r, l: (DAY, r, *LIQUID, l), (1, False, False)),
    (power_series_liquid_dynamic_incompressible, kamata_liquid_dynamic_incompressible,
     lambda r, l: (DAY, r, LIQUID[0], l), (1, False, True)),
    # The derived forms.
    (power_series_solid_dynamic_incompressible, takeuchi_solid_dynamic_incompressible,
     lambda r, l: (M2_FREQUENCY, r, SOLID[0], SOLID[2], l), (0, False, True)),
    (power_series_solid_static_incompressible, takeuchi_solid_static_incompressible,
     lambda r, l: (r, SOLID[0], SOLID[2], l), (0, True, True)),
    (power_series_solid_static_incompressible, kamata_solid_static_incompressible,
     lambda r, l: (r, SOLID[0], SOLID[2], l), (0, True, True)),
    (power_series_liquid_dynamic_incompressible, takeuchi_liquid_dynamic_incompressible,
     lambda r, l: (DAY, r, LIQUID[0], l), (1, False, True)),
)


@pytest.mark.parametrize("series, closed_form, arguments, layer", SPAN_CASES,
                         ids=[f"{case[1].__name__}" for case in SPAN_CASES])
@pytest.mark.parametrize("degree_l", (1, 2, 3, 10))
@pytest.mark.parametrize("radius", (1.0e4, 3.0e5, 1.0e6))
def test_series_and_closed_forms_span_the_same_solutions(series, closed_form, arguments, layer, degree_l, radius):
    """Both starts describe the same regular solutions at a homogeneous center. A dynamic compressible liquid at the
    largest radius is past the series' reach (its solutions grow exponentially), where the series may refuse."""
    from_series = _empty(*layer)
    try:
        series(*arguments(radius, degree_l), None, from_series)
    except RuntimeError:
        assert series is power_series_liquid_dynamic_compressible and radius == 1.0e6
        return
    from_closed_form = _empty(*layer)
    closed_form(*arguments(radius, degree_l), None, from_closed_form)
    assert _subspace_sine(from_series, from_closed_form) < 1.0e-9


@pytest.mark.parametrize("shear", (1.25e11 + 2.0e9j, 1.0e6 + 1.0e6j))
@pytest.mark.parametrize("frequency", (0.0, M2_FREQUENCY))
@pytest.mark.parametrize("degree_l", (1, 2, 10))
@pytest.mark.parametrize("radius", (1.0e4, 3.0e5))
def test_series_solid_solutions_are_the_takeuchi_wavenumber_solutions(shear, frequency, degree_l, radius):
    """Each non-polynomial series solution is one Takeuchi and Saito wavenumber solution (index 1 against the k2_neg
    solution, index 2 against k2_pos), not a mix of the two, for a stiff and a weak solid."""
    from_series = _empty(0, False, False)
    from_takeuchi = _empty(0, False, False)
    try:
        power_series_solid_dynamic_compressible(frequency, radius, SOLID[0], SOLID[1], shear, degree_l, None,
                                                from_series)
    except RuntimeError:
        # Only a weak solid's fast wave grows steeply enough for the series to refuse.
        assert abs(shear) < 1.0e9
        return
    takeuchi_solid_dynamic_compressible(frequency, radius, SOLID[0], SOLID[1], shear, degree_l, None, from_takeuchi)
    for series_row, takeuchi_row in ((1, 0), (2, 1)):
        assert _subspace_sine(from_series[series_row:series_row + 1], from_takeuchi[takeuchi_row:takeuchi_row + 1]) \
            < 1.0e-7


def _homogeneous_solve(degree_l, starting_method, is_static=False, is_incompressible=False, shear=HOMOGENEOUS_SHEAR,
                       solve_for=("tidal", "loading")):
    radius = np.linspace(0.0, EARTH_RADIUS, NUM_SLICES)
    return radial_solver(
        radius,
        HOMOGENEOUS_DENSITY * np.ones_like(radius),
        HOMOGENEOUS_BULK * np.ones(NUM_SLICES, dtype=np.complex128),
        shear * np.ones(NUM_SLICES, dtype=np.complex128),
        M2_FREQUENCY,
        HOMOGENEOUS_DENSITY,
        ("solid",),
        (is_static,),
        (is_incompressible,),
        np.asarray((EARTH_RADIUS,)),
        degree_l=degree_l,
        solve_for=solve_for,
        starting_method=starting_method,
        raise_on_fail=True)


@pytest.mark.parametrize("degree_l", tuple(TABLE_C4))
def test_power_series_reproduces_martens_homogeneous_earth(degree_l):
    """Martens' Table C.4 potential and load Love numbers, and the same values as the Takeuchi start."""
    solution = _homogeneous_solve(degree_l, "power_series")
    values = []
    for ytype in (0, 1):
        values += [solution.h[ytype].real, degree_l * solution.l[ytype].real, degree_l * solution.k[ytype].real]
    np.testing.assert_allclose(values, TABLE_C4[degree_l], rtol=0.0, atol=TABLE_C4_ATOL)

    reference = _homogeneous_solve(degree_l, "takeuchi")
    np.testing.assert_allclose(solution.k, reference.k, rtol=1.0e-6)
    np.testing.assert_allclose(solution.h, reference.h, rtol=1.0e-6)


@pytest.mark.parametrize("starting_method", ("takeuchi", "kamata", "power_series"))
@pytest.mark.parametrize("shear", (1.0e10, 1.0e11))
def test_static_incompressible_solid_matches_love_1911(starting_method, shear):
    """Every regular start covers a static incompressible solid; k2 matches Love's (1911) closed form."""
    solution = _homogeneous_solve(2, starting_method, is_static=True, is_incompressible=True, shear=shear,
                                  solve_for=("tidal",))
    gravity = 4.0 * np.pi * TidalPy.constants.G * HOMOGENEOUS_DENSITY * EARTH_RADIUS / 3.0
    expected = 1.5 / (1.0 + 19.0 * shear / (2.0 * HOMOGENEOUS_DENSITY * gravity * EARTH_RADIUS))
    assert complex(solution.k).real == pytest.approx(expected, rel=1.0e-6)


@pytest.mark.parametrize("degree_l", (2, 3, 10))
def test_unity_matches_the_closed_form_start(degree_l):
    """Unit vectors leave singular content that decays outward as (r0 / r)^(2l - 1) in a solid; at the automatic start
    it is far below the solve's tolerance on a homogeneous body."""
    unity = _homogeneous_solve(degree_l, "unity")
    reference = _homogeneous_solve(degree_l, "takeuchi")
    np.testing.assert_allclose(unity.k, reference.k, rtol=1.0e-5)


def _weak_start_solve(starting_method, viscosity, degree_l, weak_radius):
    """A Maxwell solid (60 GPa, the given viscosity) out to weak_radius under an elastic 60 GPa lid, at a 1.8 day
    period: the start layer is a weak solid (|mu*| of a few MPa)."""
    frequency = 4.11e-5
    weak = np.linspace(0.0, weak_radius, 60)
    lid = np.linspace(weak_radius, EARTH_RADIUS, 60)
    weak_shear = complex(Maxwell().calc_complex_modulus(6.0e10, viscosity, frequency))
    return radial_solver(
        np.concatenate((weak, lid)),
        HOMOGENEOUS_DENSITY * np.ones(120),
        2.0e11 * np.ones(120, dtype=np.complex128),
        np.concatenate((weak_shear * np.ones(60), 6.0e10 * np.ones(60))).astype(np.complex128),
        frequency,
        HOMOGENEOUS_DENSITY,
        ("solid", "solid"),
        (False, False),
        (False, False),
        np.asarray((weak_radius, EARTH_RADIUS)),
        degree_l=degree_l,
        solve_for=("tidal",),
        starting_method=starting_method,
        integration_rtol=1.0e-10,
        integration_atol=1.0e-10,
        warnings=False,
        raise_on_fail=True)


@pytest.mark.parametrize("viscosity, degree_l, weak_radius, series_starts", (
    (1.0e11, 5, 6000.0e3, True), (1.0e12, 5, 6200.0e3, True), (1.0e13, 10, 6200.0e3, True),
    (1.0e11, 10, 6200.0e3, False), (1.0e12, 10, 6200.0e3, False)))
def test_weak_solid_start_layer(viscosity, degree_l, weak_radius, series_starts):
    """With a weak solid starting layer, the series' Love number agrees with Kamata's at a tight tolerance, or, where
    the fast wavenumber times the starting radius is past the series' bound (the degree 10 starts), the series
    refuses."""
    if not series_starts:
        with pytest.raises(SolutionFailedError, match="power series starting conditions refused"):
            _weak_start_solve("power_series", viscosity, degree_l, weak_radius)
        return
    reference = complex(_weak_start_solve("kamata", viscosity, degree_l, weak_radius).k)
    series = complex(_weak_start_solve("power_series", viscosity, degree_l, weak_radius).k)
    assert abs(series - reference) <= 1.0e-6 * abs(reference)


@pytest.mark.parametrize("shear, radius, degree_l", ((1.0e-8, 0.15, 2), (1.0e-7, 0.2, 3), (1.0e-5, 0.5, 5)))
def test_series_refuses_past_its_wavenumber_bound(shear, radius, degree_l):
    """In a very weak static solid (non-dimensional rho = gamma = R = 1, |k| r of 16 to 18), roundoff from the fast
    wavenumber solution takes the slow one's digits; the series refuses."""
    with pytest.raises(RuntimeError, match="refused to start"):
        power_series_solid_static_compressible(radius, 1.0, 3.0 + 0.0j, shear + 0.0j, degree_l, 3.0 / (4.0 * np.pi),
                                               _empty(0, True, False))


def test_series_starts_it_accepts_match_takeuchi_in_weak_solids():
    """Across weak solids (|mu| = 1e-6 to 1e-3 rho g R, static and dynamic, l = 2 to 10, r = 0.01 to 0.5 R), every start
    the series accepts matches the Takeuchi solutions one by one, and it accepts a good share of them. Degree 1 is left
    out: a rigid translation solves the static degree-1 equations, so single solutions are not defined there. Below
    about 1e-6 rho g R an accepted static start can differ by up to 1e-5 (starting_conditions.md)."""
    G = 3.0 / (4.0 * np.pi)
    accepted = 0
    total = 0
    for shear in (1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3):
        for radius in (0.01, 0.05, 0.15, 0.3, 0.5):
            for degree_l in (2, 3, 5, 10):
                for frequency in (0.0, 0.3):
                    total += 1
                    from_series = _empty(0, False, False)
                    try:
                        power_series_solid_dynamic_compressible(frequency, radius, 1.0, 3.0 + 0.0j, shear + 0.0j,
                                                                degree_l, G, from_series)
                    except RuntimeError:
                        continue
                    accepted += 1
                    from_takeuchi = _empty(0, False, False)
                    takeuchi_solid_dynamic_compressible(frequency, radius, 1.0, 3.0 + 0.0j, shear + 0.0j, degree_l,
                                                        G, from_takeuchi)
                    for series_row, takeuchi_row in ((1, 0), (2, 1)):
                        sine = _subspace_sine(from_series[series_row:series_row + 1],
                                              from_takeuchi[takeuchi_row:takeuchi_row + 1])
                        assert sine < 1.0e-6, (shear, radius, degree_l, frequency, sine)
    assert accepted >= total // 3


def test_series_refuses_where_it_cannot_converge():
    """A dynamic compressible liquid at a long period and a large radius: the wrapper raises rather than return an
    inaccurate start."""
    with pytest.raises(RuntimeError, match="refused to start"):
        power_series_liquid_dynamic_compressible(2.0 * np.pi / (30.0 * 86400.0), 3.0e6, 10000.0, 1.0e11, 10, None,
                                                 _empty(1, False, False))


def test_configured_starting_method_reaches_the_solver(restore_config):
    """[radial_solver] starting_method is the default of every solve that does not pass its own; a call's own
    argument still wins."""
    by_name = {method: complex(_homogeneous_solve(3, method, solve_for=("tidal",)).k) for method in STARTING_METHODS}
    # The methods differ at the integration-error level, which is what tells them apart here.
    assert len(set(by_name.values())) == len(STARTING_METHODS)
    for method in STARTING_METHODS:
        TidalPy.reinit(provided_config={"radial_solver": {"starting_method": method}})
        assert complex(_homogeneous_solve(3, None, solve_for=("tidal",)).k) == by_name[method]
        assert complex(_homogeneous_solve(3, "kamata", solve_for=("tidal",)).k) == by_name["kamata"]
