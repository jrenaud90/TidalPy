"""Love solves reject undefined inputs (degree, frequency, start tolerance) and keep settings unit-independent."""
import numpy as np
import pytest

from TidalPy.RadialSolver import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.Rheology import Maxwell, Elastic
from TidalPy.Structures import build_world

_IO_FREQUENCY = 4.11e-5   # [rad s-1]


@pytest.mark.parametrize("kwargs", [
    dict(degree_l=1),
    dict(degree_l=0),
    dict(frequency=0.0),
    dict(frequency=-_IO_FREQUENCY),
    dict(frequency=float("nan")),
    dict(start_radius_tol=2.0),
    dict(start_radius_tol=-1.0e-5),
])
def test_the_world_path_refuses_invalid_input(io, kwargs):
    """The world Love solve raises ValueError for an undefined degree, frequency, or start tolerance."""
    call = dict(frequency=_IO_FREQUENCY, degree_l=2)
    call.update(kwargs)
    with pytest.raises(ValueError):
        io.solve_love_numbers(**call)


def _one_layer_inputs():
    return build_rs_input_homogeneous_layers(
        1.0e6, _IO_FREQUENCY, (3000.0,), (1.0e11,), (5.0e10,), (1.0e30,), (1.0e20,),
        ("solid",), (True,), (False,), Maxwell(), Elastic(),
        radius_fraction_tuple=(1.0,), slice_per_layer=50)


@pytest.mark.parametrize("degree_l, solve_for", [
    (1, ("tidal",)), (1, ("free",)), (1, ("tidal", "loading")), (0, ("loading",)), (-1, ("tidal",))])
def test_the_standalone_solver_refuses_an_undefined_degree(degree_l, solve_for):
    """The standalone solver raises ValueError for a degree undefined for the requested solve."""
    with pytest.raises(ValueError):
        radial_solver(*_one_layer_inputs(), degree_l=degree_l, solve_for=solve_for)


def test_a_degree_one_load_fails_as_singular_without_a_frame():
    """A degree-1 load solve fails cleanly as singular."""
    # Without a reference-frame condition (Farrell 1972, not implemented) a rigid translation solves the system.
    solution = radial_solver(*_one_layer_inputs(), degree_l=1, solve_for=("loading",))
    assert not solution.success
    assert solution.error_code == -13


def test_homogeneous_methods_use_the_solved_mass(io):
    """Homogeneous k2 uses the EOS mass, not the declared mass."""
    config = io.get_config_dict()
    k2 = io.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2, love_method="homogeneous")["love_number_k"]
    lighter = dict(config, mass_kg=0.8 * config["mass_kg"])
    other = build_world(lighter)
    other.solve_eos()
    assert other.planet_mass_eos == pytest.approx(io.planet_mass_eos, rel=1e-10)
    k2_other = other.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2, love_method="homogeneous")["love_number_k"]
    assert k2_other == pytest.approx(k2, rel=1e-10)


def test_a_step_limit_in_meters_is_the_same_in_either_unit_system(io):
    """max_step in meters gives the same step count with and without nondimensionalization."""
    steps = {}
    for nondimensionalize in (True, False):
        io.solve_love_numbers(
            frequency=_IO_FREQUENCY,
            degree_l=2,
            nondimensionalize=nondimensionalize,
            max_step=5.0e3)
        steps[nondimensionalize] = int(np.sum(io.release_radial_solution().steps_taken))
    # Io is 1.8e6 m across, so a 5e3 m limit forces about 360 steps per integrated solution, well above what the
    # tolerances alone take in either unit system (about 100 and 230 in all).
    assert steps[True] > 900
    assert steps[False] == pytest.approx(steps[True], rel=0.2)
