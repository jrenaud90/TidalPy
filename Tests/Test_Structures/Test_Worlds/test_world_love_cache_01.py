"""The cached radial Love solver on LayeredWorld reproduces uncached results and rebuilds after an EOS re-solve."""

import cmath
import math

import pytest

from TidalPy.constants import G


_PLANET_RADIUS = 6.0e6     # [m]
_DENSITY       = 4000.0    # [kg/m^3]
_STATIC_SHEAR  = 6.0e10    # [Pa]
_STATIC_BULK   = 1.3e11    # [Pa]
_SHEAR_VISC    = 1.0e21    # [Pa s]


def _solid_world():
    """Single solid uniform Maxwell sphere."""
    from TidalPy.Structures.worlds.layered import LayeredWorld
    from TidalPy.Structures.layers.physics import PhysicsLayer
    from TidalPy.Material.eos.material_eos import ConstantDensityEOS
    from TidalPy.Viscosity import make_viscosity
    from TidalPy.Rheology.rheology import Maxwell

    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = LayeredWorld("solid_planet", _PLANET_RADIUS, mass)
    layer = PhysicsLayer(
        "mantle",
        0,
        0.0,
        _PLANET_RADIUS,
        mass,
    )
    layer.set_eos(ConstantDensityEOS(
        reference_density=_DENSITY, shear_modulus_static=_STATIC_SHEAR, bulk_modulus_static=_STATIC_BULK))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _SHEAR_VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e30}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Maxwell())
    world.add_layer(layer)
    return world


def _solved_world():
    world = _solid_world()
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    return world


def _klh_after(world, frequencies):
    """Solve at each frequency in turn and return (k, h, l) from the last solve."""
    for frequency in frequencies:
        world.solve_love_numbers(frequency=frequency, verbose=False)
    return (world.love_number_k, world.love_number_h, world.love_number_l)


def _assert_close(a, b, rel=1.0e-12):
    assert cmath.isclose(a, b, rel_tol=rel, abs_tol=1.0e-300), f"{a} != {b}"


@pytest.mark.parametrize(
    "first_sweep, second_sweep, same_world",
    [
        ((1.0e-5,), (1.0e-5,), True),
        ((1.0e-5,), (3.0e-6, 1.0e-5), True),
        ((2.0e-6, 1.0e-5), (1.0e-5,), False),
    ],
    ids=["repeated_same_frequency", "interleaved_sweep", "swept_vs_fresh_world"],
)
def test_cache_does_not_change_results(first_sweep, second_sweep, same_world):
    """Solving at 1e-5 after any prior solves matches the first result exactly."""
    world = _solved_world()
    first = _klh_after(world, first_sweep)
    second_world = world if same_world else _solved_world()
    second = _klh_after(second_world, second_sweep)
    for a, b in zip(first, second):
        _assert_close(a, b)


def test_frequency_changes_dissipation():
    """Im(k2) is larger near the Maxwell resonance than in the near-elastic regime."""
    world = _solved_world()
    # Maxwell tau = eta / mu ~ 1.7e10 s: 1e-5 is near-elastic, 6e-11 is near omega * tau ~ 1.
    world.solve_love_numbers(frequency=1.0e-5, verbose=False)
    im_elastic = abs(world.love_number_k.imag)
    world.solve_love_numbers(frequency=6.0e-11, verbose=False)
    im_resonant = abs(world.love_number_k.imag)
    assert im_resonant > im_elastic


def test_eos_resolve_invalidates_cache():
    """A Love solve after re-solving the EOS with the same inputs gives the same k2."""
    world = _solved_world()
    world.solve_love_numbers(frequency=1.0e-5, verbose=False)
    k_before = world.love_number_k

    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    world.solve_love_numbers(frequency=1.0e-5, verbose=False)
    assert world.love_solved is True
    k_after = world.love_number_k
    _assert_close(k_before, k_after, rel=1.0e-9)
