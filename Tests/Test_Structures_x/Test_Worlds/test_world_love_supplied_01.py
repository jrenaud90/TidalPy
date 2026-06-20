"""
Tests for LayeredWorld.solve_love_numbers_supplied: the externally-supplied-moduli
radial-solve path used by the standalone RadialSolver_x.radial_solver API.

The supplied-moduli path must reproduce the rheology-driven path when fed the exact
complex moduli the layer rheology would produce on the EOS grid.

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""

import cmath
import math

import numpy as np
import pytest

from TidalPy.constants import G


_PLANET_RADIUS = 6.0e6
_DENSITY       = 4000.0
_STATIC_SHEAR  = 6.0e10
_STATIC_BULK   = 1.3e11
_SHEAR_VISC    = 1.0e21
_FREQ          = 1.0e-5


def _maxwell_world():
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    from TidalPy.structures_x.layers.physics import PhysicsLayer
    from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
    from TidalPy.viscosity_x import make_viscosity
    from TidalPy.rheology_x.rheology import Maxwell

    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = LayeredWorld("supplied_planet", _PLANET_RADIUS, mass)
    layer = PhysicsLayer("mantle", 0, 0.0, _PLANET_RADIUS, mass,
                         shear_modulus_static_pa=_STATIC_SHEAR,
                         bulk_modulus_static_pa=_STATIC_BULK)
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _SHEAR_VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": 1.0e30}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Maxwell())
    world.add_layer(layer)
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    return world


def test_supplied_matches_rheology():
    """Feeding the rheology-computed moduli back through the supplied path must agree."""
    world = _maxwell_world()

    # Rheology-driven reference.
    world.solve_love_numbers(frequency_rad_s=_FREQ, verbose=False)
    k_ref = world.love_number_k
    h_ref = world.love_number_h
    l_ref = world.love_number_l

    # Sample the moduli the rheology produces on the EOS grid, then feed them back.
    eos = world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    radius = np.ascontiguousarray(eos["radius"], dtype=np.float64)
    shear = np.ascontiguousarray(
        world.calc_complex_shear_modulus(radius, _FREQ), dtype=np.complex128)
    bulk = np.ascontiguousarray(
        world.calc_complex_bulk_modulus(radius, _FREQ), dtype=np.complex128)

    res = world.solve_love_numbers_supplied(shear, bulk, radius, frequency_rad_s=_FREQ)
    assert res["success"] is True
    assert cmath.isclose(res["love_number_k"], k_ref, rel_tol=1e-6, abs_tol=1e-9)
    assert cmath.isclose(res["love_number_h"], h_ref, rel_tol=1e-6, abs_tol=1e-9)
    assert cmath.isclose(res["love_number_l"], l_ref, rel_tol=1e-6, abs_tol=1e-9)


def test_supplied_k2_reasonable():
    world = _maxwell_world()
    eos = world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    radius = np.ascontiguousarray(eos["radius"], dtype=np.float64)
    shear = np.ascontiguousarray(
        world.calc_complex_shear_modulus(radius, _FREQ), dtype=np.complex128)
    bulk = np.ascontiguousarray(
        world.calc_complex_bulk_modulus(radius, _FREQ), dtype=np.complex128)
    res = world.solve_love_numbers_supplied(shear, bulk, radius, frequency_rad_s=_FREQ)
    k2 = res["love_number_k"]
    assert not cmath.isnan(k2)
    assert 0.0 < k2.real < 1.5


def test_supplied_length_mismatch_raises():
    world = _maxwell_world()
    radius = np.linspace(0.0, _PLANET_RADIUS, 50, dtype=np.float64)
    shear = np.ones(50, dtype=np.complex128) * _STATIC_SHEAR
    bulk = np.ones(49, dtype=np.complex128) * _STATIC_BULK   # wrong length
    with pytest.raises(ValueError):
        world.solve_love_numbers_supplied(shear, bulk, radius, frequency_rad_s=_FREQ)
