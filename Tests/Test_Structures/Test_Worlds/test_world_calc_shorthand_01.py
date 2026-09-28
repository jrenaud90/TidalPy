"""World shorthand bundles (get_static_viscoelastics, get_state) and the calc_* variants that solve the EOS first."""

import math

import numpy as np
import pytest

from TidalPy.constants import G


_PLANET_RADIUS = 6.0e6
_DENSITY       = 4000.0
_STATIC_SHEAR  = 6.0e10
_STATIC_BULK   = 1.3e11
_SHEAR_VISC    = 1.0e21


def _world(solve=False):
    from TidalPy.Structures.worlds.layered import LayeredWorld
    from TidalPy.Structures.layers.physics import PhysicsLayer
    from TidalPy.Material.eos.material_eos import ConstantDensityEOS
    from TidalPy.Viscosity import make_viscosity

    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = LayeredWorld("rocky", _PLANET_RADIUS, mass)
    layer = PhysicsLayer("mantle", 0, 0.0, _PLANET_RADIUS, mass)
    layer.set_eos(ConstantDensityEOS(
        reference_density=_DENSITY, shear_modulus_static=_STATIC_SHEAR, bulk_modulus_static=_STATIC_BULK))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _SHEAR_VISC}))
    world.add_layer(layer)
    if solve:
        world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    return world


def test_get_static_viscoelastics_tuple():
    world = _world(solve=True)
    radius = _PLANET_RADIUS * 0.5
    shear_mod, shear_visc, bulk_mod, bulk_visc = world.get_static_viscoelastics(radius)
    assert shear_mod == pytest.approx(world.get_shear_modulus(radius))
    assert shear_visc == pytest.approx(world.get_shear_viscosity(radius))
    assert bulk_mod == pytest.approx(world.get_bulk_modulus(radius))
    assert math.isnan(bulk_visc)  # No bulk viscosity model is attached.


def test_get_static_viscoelastics_array():
    world = _world(solve=True)
    radii = np.linspace(0.0, _PLANET_RADIUS, 9)
    bundle = world.get_static_viscoelastics(radii)
    assert len(bundle) == 4
    for arr in bundle:
        assert isinstance(arr, np.ndarray)
        assert arr.shape == radii.shape


def test_get_state_keys_and_values():
    world = _world(solve=True)
    radius = _PLANET_RADIUS * 0.5
    state = world.get_state(radius)
    for key in ("density", "gravity", "pressure", "shear_modulus",
                "shear_viscosity", "bulk_modulus", "bulk_viscosity"):
        assert key in state
    assert state["density"] == pytest.approx(world.get_density(radius))
    assert state["shear_modulus"] == pytest.approx(world.get_shear_modulus(radius))


def test_get_state_array():
    world = _world(solve=True)
    radii = np.linspace(0.0, _PLANET_RADIUS, 7)
    state = world.get_state(radii)
    assert state["density"].shape == radii.shape
    assert state["pressure"].shape == radii.shape


# =====================================================================================================================
# calc_* solve-triggering
# =====================================================================================================================
@pytest.mark.parametrize("call, expected, rel", [
    pytest.param(lambda world, radius: world.calc_density(radius), _DENSITY, 0.05, id="calc_density"),
    # No rheology is attached, so the complex modulus is the real static one.
    pytest.param(lambda world, radius: world.calc_complex_shear_modulus(radius, 1.0e-5).real, _STATIC_SHEAR, 1e-6,
                 id="calc_complex_shear_modulus"),
    pytest.param(lambda world, radius: world.calc_state(radius)["density"], _DENSITY, 0.05, id="calc_state"),
])
def test_calc_solves_when_unsolved(call, expected, rel):
    world = _world(solve=False)
    assert world.eos_solved is False
    value = call(world, _PLANET_RADIUS * 0.5)
    assert world.eos_solved is True
    assert value == pytest.approx(expected, rel=rel)


def test_force_recalc_reruns():
    world = _world(solve=True)
    value = world.calc_density(_PLANET_RADIUS * 0.5, force_recalc=True)
    assert value == pytest.approx(_DENSITY, rel=0.05)
