"""Love numbers of a gas-giant world, whose gas layer is solved as a static liquid."""

import math

from TidalPy.Material import Material, Phase
from TidalPy.Structures import build_world
from TidalPy.Structures.layers import Layer


def _gas_material():
    """A gas envelope: a liquid-only constant-density material."""
    return Material(liquid=Phase(eos={"model": "constant", "reference_density_kg_m3": 1300.0}))


def test_gas_layer_defaults_to_a_static_liquid():
    """A layer of a liquid-only (gas) material defaults to a static liquid."""
    layer = Layer(
        "envelope",
        0,
        0.0,
        1.0e7,
        1.0e27,
        _gas_material(),
    )
    assert layer.state == "auto"
    assert layer.is_liquid is True
    assert layer.is_static is True


def test_bundled_jupiter_envelope_is_a_static_liquid():
    """The bundled Jupiter's single gas layer is a liquid."""
    world = build_world("jupiter_simple")
    assert [layer.is_liquid for layer in world] == [True]


def test_jupiter_simple_love_number_is_the_uniform_fluid_value():
    """A single constant-density fluid layer has k2 = 3/2."""
    world = build_world("jupiter_simple")
    eos = world.solve_eos()
    assert eos["success"]
    result = world.solve_love_numbers(frequency=2.0 * math.pi / (0.41 * 86400.0), degree_l=2)
    assert result["success"], result["message"]
    k2 = world.love_number_k
    assert math.isclose(k2.real, 1.5, rel_tol=1.0e-6)
    assert abs(k2.imag) < 1.0e-9
