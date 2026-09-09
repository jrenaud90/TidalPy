"""LayeredWorld Love solves with non-tidal surface boundary conditions (``solve_for``).

``solve_love_numbers`` / ``solve_love_numbers_supplied`` accept ``solve_for`` = ``'tidal'``
(default), ``'loading'`` (load Love numbers k', h', l'), or ``'free'``, using the same names as the
standalone ``radial_solver``. The loading path must reproduce the standalone solver's load Love
numbers on the same homogeneous structure.
"""

import cmath
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
from TidalPy.viscosity_x import make_viscosity
from TidalPy.rheology_x.rheology import Maxwell
from TidalPy.RadialSolver_x import homogeneous_love_numbers

_PLANET_RADIUS = 6.0e6
_DENSITY = 4000.0
_STATIC_SHEAR = 6.0e10
_STATIC_BULK = 1.3e11
_SHEAR_VISC = 1.0e21
_FREQ = 1.0e-5


def _maxwell_world():
    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = LayeredWorld("loading_planet", _PLANET_RADIUS, mass)
    layer = PhysicsLayer("mantle", 0, 0.0, _PLANET_RADIUS, mass,
                         shear_modulus_static_pa=_STATIC_SHEAR,
                         bulk_modulus_static_pa=_STATIC_BULK)
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _SHEAR_VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": 1.0e30}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Maxwell())
    world.add_layer(layer)
    world.solve_eos(G_to_use=G, verbose=False)
    return world


def test_loading_matches_standalone_solver():
    """World-level load Love numbers agree with the standalone solver on the same structure."""
    world = _maxwell_world()
    layer = world.mantle
    result = world.solve_love_numbers(frequency_rad_s=_FREQ, solve_for='loading', verbose=False)
    assert result["success"] is True

    mu_complex = layer.calc_complex_shear_modulus(_FREQ)
    bulk_complex = layer.calc_complex_bulk_modulus(_FREQ)
    solution = homogeneous_love_numbers(
        _PLANET_RADIUS, _DENSITY, mu_complex, _FREQ,
        complex_bulk_modulus=bulk_complex,
        num_slices=100,
        layer_is_static=layer.is_static,
        layer_is_incompressible=layer.is_incompressible,
        solve_for=('loading',))
    assert solution.success

    assert cmath.isclose(result["love_number_k"], solution.k[0], rel_tol=1e-4, abs_tol=1e-8)
    assert cmath.isclose(result["love_number_h"], solution.h[0], rel_tol=1e-4, abs_tol=1e-8)
    assert cmath.isclose(result["love_number_l"], solution.l[0], rel_tol=1e-4, abs_tol=1e-8)


def test_loading_differs_from_tidal():
    """Load Love numbers are a different quantity from tidal Love numbers (k' is negative)."""
    world = _maxwell_world()
    tidal = world.solve_love_numbers(frequency_rad_s=_FREQ, solve_for='tidal', verbose=False)
    loading = world.solve_love_numbers(frequency_rad_s=_FREQ, solve_for='loading', verbose=False)
    assert tidal["success"] and loading["success"]
    assert tidal["love_number_k"].real > 0.0
    assert loading["love_number_k"].real < 0.0
    assert not cmath.isclose(tidal["love_number_k"], loading["love_number_k"], rel_tol=1e-3)


def test_free_surface_runs():
    world = _maxwell_world()
    result = world.solve_love_numbers(frequency_rad_s=_FREQ, solve_for='free', verbose=False)
    assert result["success"] is True


def test_supplied_path_accepts_solve_for():
    """The supplied-moduli path takes the same solve_for names and matches the rheology path."""
    world = _maxwell_world()
    reference = world.solve_love_numbers(frequency_rad_s=_FREQ, solve_for='loading', verbose=False)

    eos = world.solve_eos(G_to_use=G, verbose=False)
    radius = np.ascontiguousarray(eos["radius"], dtype=np.float64)
    shear = np.ascontiguousarray(world.calc_complex_shear_modulus(radius, _FREQ), dtype=np.complex128)
    bulk = np.ascontiguousarray(world.calc_complex_bulk_modulus(radius, _FREQ), dtype=np.complex128)
    supplied = world.solve_love_numbers_supplied(
        shear, bulk, radius, frequency_rad_s=_FREQ, solve_for='loading')
    assert supplied["success"] is True
    assert cmath.isclose(supplied["love_number_k"], reference["love_number_k"],
                         rel_tol=1e-6, abs_tol=1e-9)


def test_unknown_solve_for_raises():
    world = _maxwell_world()
    with pytest.raises(ValueError, match="solve_for"):
        world.solve_love_numbers(frequency_rad_s=_FREQ, solve_for='pressure', verbose=False)
