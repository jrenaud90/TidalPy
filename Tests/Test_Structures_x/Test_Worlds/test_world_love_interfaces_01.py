"""
Love numbers of a world whose layers change stiffness across an interface.

Every interface radius appears twice in the world's slice grid: once as the top of the lower layer and once as the
base of the upper one. Each copy must carry its own layer's moduli, and every lookup during the radial integration
must read the calling layer's copy. If both carried the lower layer's, the upper layer's first slice interval would
ramp between the two moduli and k2 would converge only at first order in the slice count; if a lookup ignored the
layer, a solid layer ending on a zero-shear liquid would read the liquid's zero modulus at its own top and divide by
it. These tests use a stiff interior under a thin, soft top layer, where the ramp error is 1.6 percent at 10 slices
per layer, with a solid core, a liquid core that keeps a shear modulus, or a zero-shear liquid core.

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""

import cmath
import math

import numpy as np
import pytest

from TidalPy.constants import G


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


def _build_world(core_state):
    from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
    from TidalPy.rheology_x.rheology import Maxwell
    from TidalPy.structures_x.layers.physics import PhysicsLayer
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    from TidalPy.viscosity_x import make_viscosity

    layer_masses = [(4.0 / 3.0) * math.pi * (r_outer ** 3 - r_inner ** 3) * density
                    for _, r_inner, r_outer, density, _, _ in _LAYERS]
    world = LayeredWorld("interfaces", _RADIUS, sum(layer_masses))
    for index, (layer_data, mass) in enumerate(zip(_LAYERS, layer_masses)):
        name, r_inner, r_outer, density, shear, viscosity = layer_data
        if name == "core" and core_state == "liquid_zero_shear":
            shear = 0.0
        layer = PhysicsLayer(name, index, r_inner, r_outer, mass,
                             shear_modulus_static=shear,
                             bulk_modulus_static=_BULK)
        layer.set_eos(ConstantDensityEOS(reference_density=density))
        layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": viscosity}))
        layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e30}))
        layer.set_shear_rheology(Maxwell())
        layer.set_bulk_rheology(Maxwell())
        layer.is_solid = not (name == "core" and core_state != "solid")
        world.add_layer(layer)
    return world


def _solve_love_number_k(world, slices_per_layer):
    """Solve the EOS at the given slice count, then k2 from the layer rheology."""
    eos = world.solve_eos(slices_per_layer=slices_per_layer, G_to_use=G)
    assert eos["success"], eos["message"]
    world.solve_love_numbers(frequency=_FREQUENCY)
    assert world.love_success, world.love_message
    return world.love_number_k, eos


def _per_layer_moduli(world, eos, slices_per_layer):
    """Each layer's own complex moduli over its own slices, both copies of every interface radius included."""
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
    """Density and moduli are constant within each layer, so a coarse grid must already give the refined k2."""
    world = _build_world(core_state)
    refined_k, _ = _solve_love_number_k(world, 200)
    coarse_k, _ = _solve_love_number_k(world, slices_per_layer)
    assert cmath.isclose(coarse_k, refined_k, rel_tol=1.0e-5), (coarse_k, refined_k)


@pytest.mark.parametrize("core_state", _CORE_STATES)
def test_love_number_matches_the_standalone_solver_given_each_layers_moduli(core_state):
    """The standalone solver takes the arrays as given, so it is the reference for what the world should fill."""
    from TidalPy.RadialSolver_x.solver import radial_solver

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
        tuple("solid" if layer.is_solid else "liquid" for layer in layers),
        tuple(bool(layer.is_static) for layer in layers),
        tuple(bool(layer.is_incompressible) for layer in layers),
        np.array([layer.radius_outer for layer in layers], dtype=np.float64),
        degree_l=2,
        use_kamata=True,
        integration_rtol=1.0e-6,
        integration_atol=1.0e-10,
        scale_rtols_bylayer_type=True)
    assert solution.success, solution.message
    standalone_k = complex(np.ravel(solution.k)[0])
    assert cmath.isclose(world_k, standalone_k, rel_tol=1.0e-5), (world_k, standalone_k)


@pytest.mark.parametrize("core_state", _CORE_STATES)
def test_supplied_moduli_keep_each_layers_interface_values(core_state):
    """Supplying each layer's own moduli, interface copies included, must reproduce the rheology-driven solve."""
    slices_per_layer = 10
    world = _build_world(core_state)
    world_k, eos = _solve_love_number_k(world, slices_per_layer)
    radius, shear, bulk = _per_layer_moduli(world, eos, slices_per_layer)
    result = world.solve_love_numbers_supplied(shear, bulk, radius, frequency=_FREQUENCY)
    assert result["success"] is True
    assert cmath.isclose(result["love_number_k"], world_k, rel_tol=1.0e-6, abs_tol=1.0e-9), \
        (result["love_number_k"], world_k)
