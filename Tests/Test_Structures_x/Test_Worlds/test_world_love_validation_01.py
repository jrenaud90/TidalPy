"""A Love-number solve refuses an input no method can answer, and its settings mean the same in any units.

Degree 1 (for a tidal or free-surface solve; a loading solve may ask for it), a non-positive or out-of-range frequency,
and a starting-radius tolerance outside (0, 1) raise ValueError on the world path and the standalone
``radial_solver`` alike, rather than failing mid-integration or starting deep in the planet and reporting success. The homogeneous methods pair the EOS surface gravity with the EOS mass, and a
step limit given in meters is honored whether the solve runs non-dimensionally or not.
"""
import numpy as np
import pytest

from TidalPy.RadialSolver_x import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.rheology_x import Maxwell, Elastic
from TidalPy.structures_x import build_world, build_world_from_dict

_IO_FREQUENCY = 4.11e-5   # [rad s-1]


@pytest.fixture(scope="module")
def io():
    world = build_world("io")
    world.solve_eos()
    return world


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
    with pytest.raises(ValueError):
        radial_solver(*_one_layer_inputs(), degree_l=degree_l, solve_for=solve_for)


def test_a_degree_one_load_fails_as_singular_without_a_frame():
    """A degree-1 load is allowed, but with no reference frame imposed a rigid translation satisfies the surface
    conditions, so the surface system is singular and the solve fails cleanly instead of returning arbitrary Love
    numbers (Farrell 1972; the frame condition is not implemented yet)."""
    solution = radial_solver(*_one_layer_inputs(), degree_l=1, solve_for=("loading",))
    assert not solution.success
    assert solution.error_code == -13


def test_homogeneous_methods_use_the_solved_mass(io):
    """The declared mass does not enter: two worlds with the same structure give the same homogeneous k2."""
    config = io.get_config_dict()
    k2 = io.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2, love_method="homogeneous")["love_number_k"]
    lighter = dict(config, mass_kg=0.8 * config["mass_kg"])
    other = build_world_from_dict(lighter)
    other.solve_eos()
    assert other.planet_mass_eos == pytest.approx(io.planet_mass_eos, rel=1e-10)
    k2_other = other.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2, love_method="homogeneous")["love_number_k"]
    assert k2_other == pytest.approx(k2, rel=1e-10)


def test_a_step_limit_in_meters_is_the_same_in_either_unit_system(io):
    """max_step is converted into the units the integration runs in, so both give the same step count."""
    steps = {}
    for nondimensionalize in (True, False):
        io.solve_love_numbers(
            frequency=_IO_FREQUENCY,
            degree_l=2,
            nondimensionalize=nondimensionalize,
            max_step=2.0e4)
        steps[nondimensionalize] = int(np.sum(io.release_radial_solution().steps_taken))
    # Io is 1.8e6 m across, so a 2e4 m limit forces about 90 steps per integrated solution.
    assert steps[True] > 200
    assert steps[False] == pytest.approx(steps[True], rel=0.2)
