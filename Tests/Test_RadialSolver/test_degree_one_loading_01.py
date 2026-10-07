"""Degree-1 load Love numbers. A rigid translation meets every degree-1 surface condition of a load, so the surface
solve fixes the frame of the body's own center of mass (CE: y5 = 1, k' = 0) in place of its last condition (Guo et al.
2004; Martens 2016), and degree1_frame shifts the result to the frames of Blewitt (2003). Checked against the exact
homogeneous-sphere value and Martens' Table C.4, Guo et al.'s PREM, Blewitt's frame relations, a static ocean's
hydrostatic surface, and the weak-lid limit of that ocean."""
from pathlib import Path

import numpy as np
import pytest

import TidalPy
from TidalPy.RadialSolver import build_rs_input_from_data, build_rs_input_homogeneous_layers, radial_solver
from TidalPy.RadialSolver.rs_solution import DEGREE1_FRAME_RESIDUAL_WARNING, check_degree1_frame_residual
from TidalPy.Rheology import Elastic
from TidalPy.Structures import build_world

from starting_methods import STARTING_METHODS

_RADIUS = 6.371e6                       # [m]
_ONE_DAY = 2.0 * np.pi / 86400.0        # [rad s-1]
_GUO_DATA = Path(__file__).resolve().parents[2] / "Benchmarks" / "RadialSolver" / "Guo+2004.npy"

# Martens (2016) Table C.4 body: Vp = 10 km/s, Vs = 5 km/s, rho = 5000 kg m-3.
_C4_DENSITY = 5000.0
_C4_SHEAR = _C4_DENSITY * 5.0e3**2
_C4_BULK = _C4_DENSITY * 10.0e3**2 - (4.0 / 3.0) * _C4_SHEAR
# The exact CE value of that body, from the Takeuchi and Saito solutions at the surface (2026-10-02 research).
_C4_H, _C4_L = -0.206731, 0.161380
# k' = y5 - 1 at the surface is held at 0 by the frame row; the dense output there differs from the integration's last
# step by the integration error.
_K_ZERO = 1.0e-8


def _homogeneous(frequency=1.0e-5, is_static=True, shear=_C4_SHEAR, bulk=_C4_BULK, num_slices=40):
    radius = np.linspace(0.0, _RADIUS, num_slices)
    return (radius, np.full(num_slices, _C4_DENSITY), np.full(num_slices, bulk + 0.0j),
            np.full(num_slices, shear + 0.0j), frequency, _C4_DENSITY, ("solid",), (is_static,), (False,),
            np.asarray((_RADIUS,)))


def _load(inputs, **kwargs):
    kwargs.setdefault("warnings", False)
    return radial_solver(*inputs, degree_l=1, solve_for=("loading",), **kwargs)


def _love(solution):
    return complex(solution.k), complex(solution.h), complex(solution.l)


# ======================================================================================================================
# The center-of-mass frame (CE)
# ======================================================================================================================
@pytest.mark.parametrize("starting_method", STARTING_METHODS)
def test_martens_table_c4_body(starting_method):
    """Every start gives the exact CE value of the C.4 body, within 0.3 percent of Martens' table (-0.2069, 0.1617)."""
    solution = _load(_homogeneous(), starting_method=starting_method)
    assert solution.success, solution.message
    k, h, l = _love(solution)
    assert abs(k) < _K_ZERO
    tolerance = 1.0e-3 if starting_method == "unity" else 2.0e-6
    assert h.real == pytest.approx(_C4_H, abs=tolerance)
    assert l.real == pytest.approx(_C4_L, abs=tolerance)
    assert h.real == pytest.approx(-0.2069, rel=3.0e-3)
    assert l.real == pytest.approx(0.1617, rel=3.0e-3)
    # Unit vectors carry singular content, which the regular solutions of a static body do not.
    assert solution.surface_frame_residual < (1.0e-3 if starting_method == "unity" else 1.0e-6)


@pytest.mark.parametrize("frequency", (1.0e-4, 1.0e-6, 1.0e-8, 1.0e-9, 1.0e-11))
def test_a_dynamic_body_is_frame_fixed_at_every_frequency(frequency):
    """With inertia the exact solution is already CE (Saito 1974 App. 2), but the system it was found from turned
    singular as omega^2: k' grew from 8e-6 at 1e-6 rad/s to 7.9 at 1e-9. The frame row holds k' at 0 and h', l' at
    the static values."""
    solution = _load(_homogeneous(frequency, is_static=False))
    assert solution.success, solution.message
    k, h, l = _love(solution)
    assert abs(k) < _K_ZERO
    # At 1e-4 rad/s (a 17 hour period) inertia moves h' and l' by about 2e-4.
    tolerance = 5.0e-4 if frequency >= 1.0e-4 else 2.0e-6
    assert h.real == pytest.approx(_C4_H, abs=tolerance)
    assert l.real == pytest.approx(_C4_L, abs=tolerance)


def test_the_answer_does_not_depend_on_tolerance_or_units():
    love = [_love(_load(_homogeneous(), integration_rtol=rtol, integration_atol=rtol))
            for rtol in (1.0e-6, 1.0e-8, 1.0e-12)]
    love += [_love(_load(_homogeneous(), nondimensionalize=False))]
    np.testing.assert_allclose(np.array(love), np.array([love[2]] * len(love)), atol=2.0e-7)


def test_a_viscoelastic_body():
    """A complex shear modulus gives complex h' and l' and still k' = 0 (so Q and lag of k' are undefined)."""
    solution = _load(_homogeneous(shear=_C4_SHEAR * (1.0 + 0.125j)))
    assert solution.success, solution.message
    k, h, l = _love(solution)
    assert abs(k) < _K_ZERO
    assert h.imag != 0.0 and l.imag != 0.0
    assert h.real == pytest.approx(_C4_H, abs=2.0e-2)
    assert solution.surface_frame_residual < 1.0e-6


def test_the_near_incompressible_limit_is_rigid():
    """In CE a homogeneous incompressible sphere does not deform under a degree-1 load: h', l' -> 0 as K/mu grows."""
    solution = _load(_homogeneous(bulk=1.0e6 * _C4_SHEAR))
    assert solution.success, solution.message
    _, h, l = _love(solution)
    assert abs(h) < 1.0e-6 and abs(l) < 1.0e-6


def _guo_earth_inputs(is_static):
    data = np.load(_GUO_DATA)
    radius = data[:, 0] * 1.0e3
    density = data[:, 1] * 1.0e3
    shear = (data[:, 2] * 1.0e3)**2 * density
    bulk = (data[:, 3] * 1.0e3)**2 * density - (4.0 / 3.0) * shear
    viscosity = np.full_like(radius, 1.0e30)
    return build_rs_input_from_data(
        _ONE_DAY, radius, density, bulk, shear, viscosity, viscosity, (1.2225e6, 3.4810e6, 6.3710e6),
        ("solid", "liquid", "solid"), is_static, (False, False, False), Elastic(), Elastic(), perform_checks=False,
        warnings=False)


@pytest.mark.skipif(not _GUO_DATA.exists(), reason="Guo+2004.npy benchmark data not found")
@pytest.mark.parametrize("is_static", ((True, True, True), (False, True, False), (True, True, False)))
def test_guo_prem(is_static):
    """Guo et al. (2004) Table 1, PREM in CE: h' = -0.285694, l' = 0.103633 (within 0.05 and 0.15 percent)."""
    solution = _load(_guo_earth_inputs(is_static), integration_method="RK45", integration_rtol=1.0e-8,
                     integration_atol=1.0e-12)
    assert solution.success, solution.message
    k, h, l = _love(solution)
    assert abs(k) < _K_ZERO
    assert h.real == pytest.approx(-0.285694, rel=5.0e-4)
    assert l.real == pytest.approx(0.103633, rel=1.5e-3)


# ======================================================================================================================
# Other frames (Blewitt 2003)
# ======================================================================================================================
@pytest.mark.parametrize("frame", ("CE", "CM", "CF", "CL", "CH", "cm"))
def test_frames_follow_blewitt(frame):
    """h', l', and 1 + k' drop by the same alpha, and the radial functions move with them: y1 and y3 by a constant c,
    y5 by c g(r), the stresses not at all."""
    reference = _load(_homogeneous())
    _, h_ce, l_ce = _love(reference)
    alpha = {"CE": 0.0, "CM": 1.0, "CF": (h_ce + 2.0 * l_ce) / 3.0, "CL": l_ce, "CH": h_ce}[frame.upper()]
    solution = _load(_homogeneous(), degree1_frame=frame)
    assert solution.success, solution.message
    k, h, l = _love(solution)
    assert abs(k - (0.0 - alpha)) < 1.0e-14
    assert abs(h - (h_ce - alpha)) < 1.0e-14
    assert abs(l - (l_ce - alpha)) < 1.0e-14

    # Blewitt's frame definitions.
    if frame.upper() == "CM":
        assert abs(1.0 + k) < 1.0e-14
    elif frame == "CF":
        assert abs(h + 2.0 * l) < 1.0e-14
    elif frame == "CL":
        assert abs(l) < 1.0e-14
    elif frame == "CH":
        assert abs(h) < 1.0e-14

    y, y_ce = solution.result, reference.result
    gravity = solution.surface_gravity
    shift = -alpha / gravity
    np.testing.assert_allclose(y[4, -1] - 1.0, k, atol=1.0e-14)
    np.testing.assert_allclose(y[0, -1] * gravity, h, atol=1.0e-14)
    np.testing.assert_allclose(y[2, -1] * gravity, l, atol=1.0e-14)
    solved = np.isfinite(y_ce[0])   # radii below the starting radius are NaN
    scale = np.max(np.abs(y_ce[0, solved]))
    np.testing.assert_allclose(y[0, solved] - y_ce[0, solved], shift, rtol=0.0, atol=1.0e-12 * scale)
    np.testing.assert_allclose(y[2, solved] - y_ce[2, solved], shift, rtol=0.0, atol=1.0e-12 * scale)
    for stress_row in (1, 3, 5):
        np.testing.assert_allclose(y[stress_row, solved], y_ce[stress_row, solved], rtol=1.0e-12, atol=0.0)


def test_an_unknown_frame_is_refused():
    with pytest.raises(ValueError, match="degree-1 frame"):
        _load(_homogeneous(), degree1_frame="CX")


# ======================================================================================================================
# A liquid surface layer
# ======================================================================================================================
_CORE = ("solid", 5500.0, 1.5e11, 6.0e10)
_MANTLE = ("solid", 3300.0, 1.2e11, 7.0e10)
_OCEAN = ("liquid", 1000.0, 2.2e9, 0.0)
_FRACTIONS = (0.5, 0.97, 1.0)


def _layered(layers, is_static, is_incompressible=None, frequency=_ONE_DAY):
    count = len(layers)
    return build_rs_input_homogeneous_layers(
        _RADIUS, frequency, tuple(layer[1] for layer in layers), tuple(layer[2] for layer in layers),
        tuple(layer[3] for layer in layers), (1.0e30,) * count, (1.0e30,) * count,
        tuple(layer[0] for layer in layers), is_static, is_incompressible or (False,) * count, "elastic", "elastic",
        radius_fraction_tuple=_FRACTIONS, slice_per_layer=40)


def _bulk_density(layers):
    shells = np.diff(np.concatenate(([0.0], np.asarray(_FRACTIONS)**3)))
    return float(np.sum(shells * np.array([layer[1] for layer in layers])))


@pytest.mark.parametrize("starting_method", STARTING_METHODS)
def test_a_static_ocean_surface_is_hydrostatic(starting_method):
    """A static liquid surface has no y3 and a degree-1 load leaves its y7 condition zero; in CE (y5 = 1) its
    hydrostatic surface, y2 = rho (g y1 - y5), gives h' = 1 - rho_bar / rho_surface exactly, whatever lies below."""
    layers = [_CORE, _MANTLE, _OCEAN]
    solution = _load(_layered(layers, (True, True, True)), starting_method=starting_method)
    assert solution.success, solution.message
    k, h, l = _love(solution)
    assert abs(k) < _K_ZERO
    assert h.real == pytest.approx(1.0 - _bulk_density(layers) / _OCEAN[1], rel=1.0e-12)
    assert np.isnan(l)


def test_a_static_ocean_is_the_limit_of_a_weakening_lid():
    """An incompressible solid lid of the ocean's density tends to the ocean's h' as its rigidity falls."""
    layers = [_CORE, _MANTLE, _OCEAN]
    ocean_h = complex(_load(_layered(layers, (True, True, True))).h).real
    errors = []
    for shear in (1.0e8, 1.0e6, 1.0e4):
        lid = ("solid", _OCEAN[1], _OCEAN[2], shear)
        solution = _load(_layered([_CORE, _MANTLE, lid], (True, True, True), (False, False, True)),
                         integration_rtol=1.0e-10, integration_atol=1.0e-12)
        assert solution.success, solution.message
        errors.append(abs(complex(solution.h).real - ocean_h))
    assert errors[0] > errors[1] > errors[2]
    assert errors[2] < 0.03 * abs(ocean_h)


@pytest.mark.parametrize("frame, offset", (("CM", 1.0), ("CH", None)))
def test_a_static_ocean_surface_takes_the_frames_that_need_no_l(frame, offset):
    layers = [_CORE, _MANTLE, _OCEAN]
    ce = complex(_load(_layered(layers, (True, True, True))).h)
    solution = _load(_layered(layers, (True, True, True)), degree1_frame=frame)
    assert solution.success, solution.message
    alpha = ce if offset is None else offset
    k, h, _ = _love(solution)
    assert abs(k + alpha) < 1.0e-12
    assert abs(h - (ce - alpha)) < 1.0e-12


@pytest.mark.parametrize("frame", ("CF", "CL"))
def test_a_static_ocean_surface_refuses_the_frames_that_need_l(frame):
    solution = _load(_layered([_CORE, _MANTLE, _OCEAN], (True, True, True)), degree1_frame=frame)
    assert not solution.success
    assert solution.error_code == -16
    assert "static liquid surface" in solution.message


def test_a_dynamic_ocean_surface():
    """A dynamic liquid surface takes the frame row in place of its y6 condition; with every layer dynamic the
    replaced condition holds to roundoff, and h' is close to that of a weak dynamic lid of the same density."""
    layers = [_CORE, _MANTLE, _OCEAN]
    solution = _load(_layered(layers, (False, False, False)))
    assert solution.success, solution.message
    k, h, l = _love(solution)
    assert abs(k) < _K_ZERO
    assert np.isfinite(h) and np.isfinite(l)
    assert solution.surface_frame_residual < 1.0e-8
    lid = ("solid", _OCEAN[1], _OCEAN[2], 1.0e4)
    weak_lid = _load(_layered([_CORE, _MANTLE, lid], (False, False, False)))
    assert complex(weak_lid.h).real == pytest.approx(h.real, rel=1.0e-2)


# ======================================================================================================================
# Mixed static and dynamic layers
# ======================================================================================================================
def test_mixed_inertia_leaves_a_residual_that_is_warned_about(spdlog_text):
    """Dynamic solids under a static ocean have no exact degree-1 frame: the replaced condition holds only to about
    omega^2 R / g. It is reported, and warned about once it is large."""
    layers = [_CORE, _MANTLE, _OCEAN]
    at_a_day = _load(_layered(layers, (False, False, True)))
    assert at_a_day.success, at_a_day.message
    assert 1.0e-4 < at_a_day.surface_frame_residual < DEGREE1_FRAME_RESIDUAL_WARNING
    fast = _load(_layered(layers, (False, False, True), frequency=1.0e-3), warnings=True)
    assert fast.success, fast.message
    assert fast.surface_frame_residual > DEGREE1_FRAME_RESIDUAL_WARNING
    assert "reference frame replaced" in spdlog_text()
    assert check_degree1_frame_residual(fast.surface_frame_residual)
    assert not check_degree1_frame_residual(float("nan"))


def test_other_solves_report_no_residual():
    solution = radial_solver(*_homogeneous(), degree_l=2, solve_for=("loading",), warnings=False)
    assert np.isnan(solution.surface_frame_residual)


# ======================================================================================================================
# World path and configuration
# ======================================================================================================================
def test_a_bundled_world_gives_degree_one_load_numbers():
    """earth_prem (every layer static) solves in CE, and a pinned frame or the call's own argument moves it."""
    world = build_world("earth_prem")
    world.solve_eos()
    result = world.solve_love_numbers(frequency=_ONE_DAY, degree_l=1, solve_for="loading")
    assert result["success"], result["message"]
    h_ce = complex(result["love_number_h"])
    assert abs(complex(result["love_number_k"])) < _K_ZERO
    assert h_ce.real == pytest.approx(-0.2857, rel=1.0e-2)
    assert world.love_surface_frame_residual < 1.0e-6

    world.set_solver_defaults(radial_solver={"degree1_frame": "CM"})
    assert world.get_solver_defaults()["radial_solver"] == {"degree1_frame": "CM"}
    pinned = world.solve_love_numbers(frequency=_ONE_DAY, degree_l=1, solve_for="loading")
    assert complex(pinned["love_number_k"]) == pytest.approx(-1.0, abs=1.0e-14)
    assert complex(pinned["love_number_h"]) == pytest.approx(h_ce - 1.0, abs=1.0e-12)
    argument = world.solve_love_numbers(frequency=_ONE_DAY, degree_l=1, solve_for="loading", degree1_frame="CE")
    assert abs(complex(argument["love_number_k"])) < _K_ZERO


def test_the_configured_frame_reaches_the_solver(restore_config):
    TidalPy.reinit(provided_config={"radial_solver": {"degree1_frame": "CM"}})
    assert complex(_load(_homogeneous()).k) == pytest.approx(-1.0, abs=1.0e-14)
    assert abs(complex(_load(_homogeneous(), degree1_frame="CE").k)) < _K_ZERO


def test_degree_two_ignores_the_frame():
    plain = radial_solver(*_homogeneous(), degree_l=2, solve_for=("loading",), warnings=False)
    framed = radial_solver(*_homogeneous(), degree_l=2, solve_for=("loading",), degree1_frame="CM", warnings=False)
    assert _love(plain) == _love(framed)
