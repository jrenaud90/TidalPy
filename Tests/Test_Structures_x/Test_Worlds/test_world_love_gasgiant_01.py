"""
Love numbers of a gas-giant world.

A gas layer carries no shear stress, so it is flagged liquid (and static) and the radial solver integrates it as a
static liquid layer. A single constant-density fluid layer has the closed-form degree-2 potential Love number 3/2.

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""

import math

from TidalPy.structures_x import build_world
from TidalPy.structures_x.layers.gas import GasLayer


def test_gas_layer_defaults_to_a_static_liquid():
    layer = GasLayer("envelope", 0, 0.0, 1.0e7, 1.0e27)
    assert layer.is_solid is False
    assert layer.is_static is True
    world = build_world("jupiter_simple")
    assert [layer.is_solid for layer in world] == [False]


def test_jupiter_simple_love_number_is_the_uniform_fluid_value():
    world = build_world("jupiter_simple")
    eos = world.solve_eos()
    assert eos["success"]
    result = world.solve_love_numbers(frequency=2.0 * math.pi / (0.41 * 86400.0), degree_l=2)
    assert result["success"], result["message"]
    k2 = world.love_number_k
    assert math.isclose(k2.real, 1.5, rel_tol=1.0e-6)
    assert abs(k2.imag) < 1.0e-9
