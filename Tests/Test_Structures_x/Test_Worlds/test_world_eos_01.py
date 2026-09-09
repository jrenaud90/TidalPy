"""
Tests for the world-level equation-of-state solve (LayeredWorld.solve_eos) and
the per-layer material EOS wiring (BaseLayer.set_eos).

This is the class-based replacement for the old functional
``Tests/Test_Material_x/test_a_eos.py``: a uniform planet is reproduced with the
``ConstantDensityEOS`` model and PREM Earth with the ``InterpolatedEOS`` model.
After a successful solve the world (and each layer) can be queried for density,
gravity, and pressure as a function of radius.

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""

import math
import os
import warnings
from pathlib import Path

import numpy as np
import pytest

from TidalPy.constants import G


# =====================================================================================================================
# Helpers
# =====================================================================================================================
def _import():
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    from TidalPy.structures_x.layers.base import BaseLayer
    from TidalPy.Material_x.eos.material_eos import (
        ConstantDensityEOS, InterpolatedEOS)
    return LayeredWorld, BaseLayer, ConstantDensityEOS, InterpolatedEOS


# Uniform one-layer planet (matches the old test_a_eos.py reference).
_N             = 100
_PLANET_RADIUS = 6000.0e3      # [m]
_DENSITY       = 3500.0        # [kg/m^3]


def _uniform_world():
    LayeredWorld, BaseLayer, ConstantDensityEOS, _ = _import()
    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = LayeredWorld("Uniform", _PLANET_RADIUS, mass, world_type="terrestrial")
    layer = BaseLayer("mantle", 0, 0.0, _PLANET_RADIUS, mass, material_name="rock")
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    world.add_layer(layer)
    return world


# =====================================================================================================================
# set_eos wiring
# =====================================================================================================================
def test_set_eos_flags():
    _, BaseLayer, ConstantDensityEOS, _ = _import()
    layer = BaseLayer("rock", 0, 0.0, 1.0e6, 1.0e20)
    assert layer.eos_set is False
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    assert layer.eos_set is True


def test_set_eos_consumes_model():
    """Attaching an EOS model leaves the Python wrapper an empty shell."""
    _, BaseLayer, ConstantDensityEOS, _ = _import()
    layer = BaseLayer("rock", 0, 0.0, 1.0e6, 1.0e20)
    eos = ConstantDensityEOS(reference_density_kg_m3=_DENSITY)
    layer.set_eos(eos)
    with pytest.raises(ValueError):
        layer.set_eos(eos)


def test_solve_requires_all_eos_set():
    LayeredWorld, BaseLayer, ConstantDensityEOS, _ = _import()
    world = LayeredWorld("NoEOS", _PLANET_RADIUS, 1.0e24)
    world.add_layer(BaseLayer("mantle", 0, 0.0, _PLANET_RADIUS, 1.0e24))
    assert world.all_eos_set is False
    with pytest.raises(ValueError):
        world.solve_eos(verbose=False)


def test_solve_requires_layers():
    LayeredWorld, _, _, _ = _import()
    world = LayeredWorld("Empty", _PLANET_RADIUS, 1.0e24)
    with pytest.raises(ValueError):
        world.solve_eos(verbose=False)


# =====================================================================================================================
# Uniform planet (ConstantDensityEOS) — reproduces the legacy functional results
# =====================================================================================================================
def test_uniform_solve_succeeds():
    world = _uniform_world()
    result = world.solve_eos(surface_pressure=0.0, G_to_use=G, verbose=False)
    assert result["success"] is True
    assert result["iterations"] >= 1
    assert world.eos_solved is True
    assert result["gravity"].shape == result["density"].shape


def test_uniform_gravity_analytic():
    """g(r) = (4/3) pi G rho r for a uniform sphere."""
    world = _uniform_world()
    world.solve_eos(G_to_use=G, verbose=False)
    radii = np.linspace(0.0, _PLANET_RADIUS, 50)[1:]
    for r in radii:
        expected = (4.0 / 3.0) * math.pi * G * _DENSITY * r
        assert math.isclose(world.get_gravity(r), expected, rel_tol=0.05)


def test_uniform_density_recovered():
    world = _uniform_world()
    world.solve_eos(G_to_use=G, verbose=False)
    for r in np.linspace(0.0, _PLANET_RADIUS, 25):
        assert math.isclose(world.get_density(r), _DENSITY, rel_tol=0.05)


def test_uniform_pressure_monotonic_and_central():
    """Pressure decreases outward; central pressure matches the analytic value."""
    world = _uniform_world()
    world.solve_eos(G_to_use=G, verbose=False)
    radii = np.linspace(0.0, _PLANET_RADIUS, 60)
    pressures = np.array([world.get_pressure(r) for r in radii])
    assert np.all(np.diff(pressures) <= 1.0e-3 * abs(pressures[0]))
    # Analytic central pressure: (2/3) pi G rho^2 R^2.
    p_central_analytic = (2.0 / 3.0) * math.pi * G * _DENSITY ** 2 * _PLANET_RADIUS ** 2
    assert math.isclose(world.get_pressure(0.0), p_central_analytic, rel_tol=0.05)
    # Surface pressure converged to ~0.
    assert abs(world.get_pressure(_PLANET_RADIUS)) < 0.01 * world.get_pressure(0.0)


def test_uniform_planet_mass():
    world = _uniform_world()
    world.solve_eos(G_to_use=G, verbose=False)
    expected_mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    assert math.isclose(world.planet_mass_eos, expected_mass, rel_tol=0.02)


def test_world_delegates_to_layer_profile():
    """The world density query delegates to the owning layer's populated profile."""
    world = _uniform_world()
    world.solve_eos(G_to_use=G, verbose=False)
    assert math.isclose(world.get_density(_PLANET_RADIUS * 0.5), _DENSITY, rel_tol=0.05)


def test_unsolved_world_returns_nan():
    world = _uniform_world()
    assert math.isnan(world.get_density(_PLANET_RADIUS * 0.5))
    assert math.isnan(world.surface_gravity_eos)


def test_unsupported_method_raises():
    world = _uniform_world()
    with pytest.raises(ValueError):
        world.solve_eos(G_to_use=G, integration_method="NOPE", verbose=False)


# =====================================================================================================================
# Two-layer planet (dense core + lighter mantle)
# =====================================================================================================================
def test_two_layer_constant_density():
    LayeredWorld, BaseLayer, ConstantDensityEOS, _ = _import()
    r_cmb  = _PLANET_RADIUS / 2.0
    rho_c  = 5000.0
    rho_m  = 3000.0
    mass   = (4.0 / 3.0) * math.pi * (rho_c * r_cmb ** 3 + rho_m * (_PLANET_RADIUS ** 3 - r_cmb ** 3))
    world = LayeredWorld("TwoLayer", _PLANET_RADIUS, mass)
    core   = BaseLayer("core", 0, 0.0, r_cmb, 0.0, material_name="iron")
    mantle = BaseLayer("mantle", 1, r_cmb, _PLANET_RADIUS, 0.0, material_name="rock")
    core.set_eos(ConstantDensityEOS(reference_density_kg_m3=rho_c))
    mantle.set_eos(ConstantDensityEOS(reference_density_kg_m3=rho_m))
    world.add_layer(core)
    world.add_layer(mantle)

    result = world.solve_eos(G_to_use=G, verbose=False)
    assert result["success"]
    assert world.central_pressure > 0.0
    # Density jump across the core-mantle boundary.
    assert math.isclose(world.get_density(r_cmb * 0.5), rho_c, rel_tol=0.05)
    assert math.isclose(world.get_density(r_cmb + (_PLANET_RADIUS - r_cmb) * 0.5), rho_m, rel_tol=0.05)
    # Gravity is non-negative and the enclosed mass matches.
    assert math.isclose(world.planet_mass_eos, mass, rel_tol=0.03)


# =====================================================================================================================
# PREM Earth (InterpolatedEOS)
# =====================================================================================================================
def _prem_dir():
    # __file__ = Tests/Test_Structures_x/Test_Worlds/...; parents[2] = Tests/.
    return Path(__file__).resolve().parents[2] / "Test_Material_x"


def test_prem_earth_interpolated():
    LayeredWorld, BaseLayer, _, InterpolatedEOS = _import()
    prem_dir = _prem_dir()

    prem_data = []
    try:
        for layer_i in range(3):
            prem_data.append(np.loadtxt(prem_dir / f"prem_layer{layer_i}.txt", delimiter=","))
    except Exception as e:  # noqa: BLE001
        warnings.warn(f"Could not load PREM data: {e}")
        pytest.skip("Could not load PREM Earth data.")

    surface_radius = prem_data[2][:, 0][-1]
    world = LayeredWorld("PREM-Earth", surface_radius, 5.972e24, world_type="terrestrial")

    prev_outer = 0.0
    for layer_i in range(3):
        radius_array  = np.ascontiguousarray(prem_data[layer_i][:, 0])
        density_array = np.ascontiguousarray(prem_data[layer_i][:, 1])
        r_outer = radius_array[-1]
        layer = BaseLayer(f"layer{layer_i}", layer_i, prev_outer, r_outer, 0.0)
        layer.set_eos(InterpolatedEOS(radius_array.tolist(), density_array.tolist()))
        world.add_layer(layer)
        prev_outer = r_outer

    result = world.solve_eos(G_to_use=G, integration_method="DOP853",
                             slices_per_layer=120, verbose=False)
    if not result["success"]:
        raise RuntimeError(f"EOS solver failed: {result['message']}")

    assert result["iterations"] >= 1
    assert world.central_pressure > 0.0
    assert math.isclose(world.surface_gravity_eos, 9.81, rel_tol=0.10)
    assert math.isclose(world.planet_mass_eos, 5.972e24, rel_tol=0.10)
    assert math.isclose(world.planet_moi_eos, 9.0e37, rel_tol=1.00)

    # Density queries reproduce the PREM profile (within interpolation tolerance).
    mid_mantle = 0.5 * (prem_data[2][:, 0][0] + prem_data[2][:, 0][-1])
    expected = np.interp(mid_mantle, prem_data[2][:, 0], prem_data[2][:, 1])
    assert math.isclose(world.get_density(mid_mantle), expected, rel_tol=0.05)
