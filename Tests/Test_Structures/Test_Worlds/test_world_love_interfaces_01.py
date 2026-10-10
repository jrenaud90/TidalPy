"""Love numbers across layer interfaces: each duplicated interface radius must carry its own layer's moduli.

A stiff interior under a thin soft top layer makes a wrong interface value visible (1.6 percent at 10 slices).
"""

import cmath
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Material import Material, Phase
from TidalPy.RadialSolver.solver import radial_solver
from TidalPy.Rheology.rheology import Maxwell
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.base import BaseWorld


_RADIUS    = 1.0e6
_FREQUENCY = 2.0 * np.pi / 86400.0
_BULK      = 1.0e11

# name, inner radius [m], outer radius [m], density [kg m-3], static shear modulus [Pa], shear viscosity [Pa s]
_LAYERS = (
    ("core",       0.0,           0.5 * _RADIUS,  8000.0, 5.0e10, 1.0e22),
    ("mantle",     0.5 * _RADIUS, 0.95 * _RADIUS, 3500.0, 5.0e10, 1.0e22),
    ("soft_layer", 0.95 * _RADIUS, _RADIUS,       3000.0, 1.0e9,  1.0e13),
)

_CORE_STATES = [pytest.param("solid", id="solid_core"), pytest.param("liquid", id="liquid_core"),
                pytest.param("liquid_zero_shear", id="liquid_core_zero_shear")]


def _material(density, shear, viscosity):
    """Constant-density solid; a zero shear modulus leaves out the shear-modulus law (a fluid)."""
    return Material(solid=Phase(
        eos={"model": "constant", "reference_density_kg_m3": density, "bulk_modulus_pa": _BULK},
        shear_modulus=None if shear == 0.0 else {"model": "constant", "shear_modulus_pa": shear},
        shear_viscosity={"model": "constant", "reference_viscosity_pas": viscosity},
        bulk_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e30}))


def _build_world(core_state):
    layer_masses = [(4.0 / 3.0) * math.pi * (r_outer ** 3 - r_inner ** 3) * density
                    for _, r_inner, r_outer, density, _, _ in _LAYERS]
    world = BaseWorld("interfaces", _RADIUS, sum(layer_masses))
    for index, (layer_data, mass) in enumerate(zip(_LAYERS, layer_masses)):
        name, r_inner, r_outer, density, shear, viscosity = layer_data
        if name == "core" and core_state == "liquid_zero_shear":
            shear = 0.0
        world.add_layer(Layer(
            name,
            index,
            r_inner,
            r_outer,
            mass,
            _material(density, shear, viscosity),
            state="liquid" if (name == "core" and core_state != "solid") else "solid",
            shear_rheology=Maxwell(),
            bulk_rheology=Maxwell()))
    return world


def _solve_love_number_k(world, slices_per_layer):
    """Solve the EOS at the given slice count, then k2."""
    eos = world.solve_eos(slices_per_layer=slices_per_layer, G_to_use=G)
    assert eos["success"], eos["message"]
    world.solve_love_numbers(frequency=_FREQUENCY)
    assert world.love_success, world.love_message
    return world.love_number_k, eos


def _per_layer_moduli(world, eos, slices_per_layer):
    """Each layer's complex moduli over its own slices, including both copies of each interface radius."""
    radius = np.ascontiguousarray(eos["radius"], dtype=np.float64)
    shear = np.empty(radius.size, dtype=np.complex128)
    bulk = np.empty(radius.size, dtype=np.complex128)
    for index, layer in enumerate(world):
        layer_slices = slice(index * slices_per_layer, (index + 1) * slices_per_layer)
        shear[layer_slices] = layer.calc_complex_shear_modulus(radius[layer_slices], _FREQUENCY)
        bulk[layer_slices] = layer.calc_complex_bulk_modulus(radius[layer_slices], _FREQUENCY)
    return radius, shear, bulk


@pytest.mark.parametrize("core_state", _CORE_STATES)
@pytest.mark.parametrize("slices_per_layer", [10, 25])
def test_love_number_does_not_depend_on_the_slice_count(core_state, slices_per_layer):
    """With layer-constant properties, a coarse slice grid already gives the refined k2."""
    world = _build_world(core_state)
    refined_k, _ = _solve_love_number_k(world, 200)
    coarse_k, _ = _solve_love_number_k(world, slices_per_layer)
    assert cmath.isclose(coarse_k, refined_k, rel_tol=1.0e-5), (coarse_k, refined_k)


@pytest.mark.parametrize("core_state", _CORE_STATES)
def test_love_number_matches_the_standalone_solver_given_each_layers_moduli(core_state):
    """World k2 matches the standalone radial solver fed each layer's own moduli."""
    slices_per_layer = 20
    world = _build_world(core_state)
    world_k, eos = _solve_love_number_k(world, slices_per_layer)
    radius, shear, bulk = _per_layer_moduli(world, eos, slices_per_layer)
    layers = list(world)
    solution = radial_solver(
        radius,
        np.ascontiguousarray(eos["density"], dtype=np.float64),
        bulk,
        shear,
        _FREQUENCY,
        world.mass / ((4.0 / 3.0) * math.pi * _RADIUS ** 3),
        tuple("liquid" if layer.is_liquid else "solid" for layer in layers),
        tuple(bool(layer.is_static) for layer in layers),
        tuple(bool(layer.is_incompressible) for layer in layers),
        np.array([layer.radius_outer for layer in layers], dtype=np.float64),
        degree_l=2,
        starting_method="kamata",
        integration_rtol=1.0e-6,
        integration_atol=1.0e-10,
        scale_rtols_bylayer_type=True)
    assert solution.success, solution.message
    standalone_k = complex(np.ravel(solution.k)[0])
    assert cmath.isclose(world_k, standalone_k, rel_tol=1.0e-5), (world_k, standalone_k)


@pytest.mark.parametrize("core_state", _CORE_STATES)
def test_supplied_moduli_keep_each_layers_interface_values(core_state):
    """Supplying each layer's own moduli reproduces the rheology-driven solve."""
    slices_per_layer = 10
    world = _build_world(core_state)
    world_k, eos = _solve_love_number_k(world, slices_per_layer)
    radius, shear, bulk = _per_layer_moduli(world, eos, slices_per_layer)
    result = world.solve_love_numbers_supplied(shear, bulk, radius, frequency=_FREQUENCY)
    assert result["success"] is True
    assert cmath.isclose(result["love_number_k"], world_k, rel_tol=1.0e-6, abs_tol=1.0e-9), \
        (result["love_number_k"], world_k)


def test_dynamic_liquid_layer_reports_a_finite_y3_on_the_world():
    """y3 in a dynamic liquid layer is finite on the world and matches the exported solution."""
    world = _build_world("liquid_zero_shear")
    core = list(world)[0]
    core.is_static = False
    _solve_love_number_k(world, 20)
    radius = 0.5 * (core.radius_inner + core.radius_outer)
    world_y3 = world.get_love_radial_y(radius, 0, 2)
    assert np.isfinite(world_y3.real) and np.isfinite(world_y3.imag), world_y3

    exported_y3 = complex(world.release_radial_solution().get_radial_solution(radius, 0)[2])
    assert cmath.isclose(world_y3, exported_y3, rel_tol=1.0e-12), (world_y3, exported_y3)


def test_exported_solution_reports_the_moduli_of_the_solve_that_made_it():
    """An exported solution reports the moduli of the solve that produced it."""
    slices_per_layer = 10
    world = _build_world("solid")
    _, eos = _solve_love_number_k(world, slices_per_layer)
    radius, shear, bulk = _per_layer_moduli(world, eos, slices_per_layer)
    mantle = list(world)[1]
    probe = 0.5 * (mantle.radius_inner + mantle.radius_outer)

    rheology_solution = world.release_radial_solution()
    assert cmath.isclose(
        rheology_solution.get_complex_shear_modulus(probe),
        mantle.calc_complex_shear_modulus(probe, _FREQUENCY), rel_tol=1.0e-12)
    assert rheology_solution.love_frequency == _FREQUENCY

    # Twice the stiffness, so the two sources cannot be confused.
    world.solve_love_numbers(frequency=_FREQUENCY)
    supplied_frequency = 3.0 * _FREQUENCY
    result = world.solve_love_numbers_supplied(2.0 * shear, bulk, radius, frequency=supplied_frequency)
    assert result["success"] is True
    supplied_solution = world.release_radial_solution()
    assert cmath.isclose(
        supplied_solution.get_complex_shear_modulus(probe),
        2.0 * mantle.calc_complex_shear_modulus(probe, _FREQUENCY), rel_tol=1.0e-9)
    assert supplied_solution.love_frequency == supplied_frequency
