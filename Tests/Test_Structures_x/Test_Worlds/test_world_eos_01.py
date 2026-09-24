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
LAYER_MASS_RTOL = 1e-12  # constant-density layer mass versus rho V after the solve


def _uniform_world():
    LayeredWorld, BaseLayer, ConstantDensityEOS, _ = _import()
    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = LayeredWorld("Uniform", _PLANET_RADIUS, mass, world_type="terrestrial")
    layer = BaseLayer("mantle", 0, 0.0, _PLANET_RADIUS, mass, material_name="rock")
    layer.set_eos(ConstantDensityEOS(reference_density=_DENSITY))
    world.add_layer(layer)
    return world


# =====================================================================================================================
# set_eos wiring
# =====================================================================================================================
def test_set_eos_flags():
    _, BaseLayer, ConstantDensityEOS, _ = _import()
    layer = BaseLayer("rock", 0, 0.0, 1.0e6, 1.0e20)
    assert layer.eos_set is False
    layer.set_eos(ConstantDensityEOS(reference_density=_DENSITY))
    assert layer.eos_set is True


def test_set_eos_consumes_model():
    """Attaching an EOS model leaves the Python wrapper an empty shell."""
    _, BaseLayer, ConstantDensityEOS, _ = _import()
    layer = BaseLayer("rock", 0, 0.0, 1.0e6, 1.0e20)
    eos = ConstantDensityEOS(reference_density=_DENSITY)
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
    core.set_eos(ConstantDensityEOS(reference_density=rho_c))
    mantle.set_eos(ConstantDensityEOS(reference_density=rho_m))
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
# Layer mass and bulk density (set by every successful solve)
# =====================================================================================================================
_R_CMB = _PLANET_RADIUS / 2.0


def _two_layer_world(rho_core=5000.0, rho_mantle=3000.0):
    """Two constant-density layers built without layer masses (as the TOML builder does)."""
    LayeredWorld, BaseLayer, ConstantDensityEOS, _ = _import()
    mass = (4.0 / 3.0) * math.pi * (rho_core * _R_CMB ** 3 + rho_mantle * (_PLANET_RADIUS ** 3 - _R_CMB ** 3))
    world = LayeredWorld("TwoLayer", _PLANET_RADIUS, mass)
    core = BaseLayer("core", 0, 0.0, _R_CMB, 0.0, material_name="iron")
    mantle = BaseLayer("mantle", 1, _R_CMB, _PLANET_RADIUS, 0.0, material_name="rock")
    core.set_eos(ConstantDensityEOS(reference_density=rho_core))
    mantle.set_eos(ConstantDensityEOS(reference_density=rho_mantle))
    world.add_layer(core)
    world.add_layer(mantle)
    return world


def test_solve_sets_layer_mass_and_bulk_density():
    """Constant-density layers get mass = rho V and bulk density = rho from the solve."""
    world = _two_layer_world()
    core, mantle = world.layers
    assert core.mass == 0.0 and mantle.mass == 0.0
    assert world.solve_eos(G_to_use=G, verbose=False)["success"]
    core_volume = (4.0 / 3.0) * math.pi * _R_CMB ** 3
    mantle_volume = (4.0 / 3.0) * math.pi * (_PLANET_RADIUS ** 3 - _R_CMB ** 3)
    assert math.isclose(core.mass, 5000.0 * core_volume, rel_tol=LAYER_MASS_RTOL)
    assert math.isclose(mantle.mass, 3000.0 * mantle_volume, rel_tol=LAYER_MASS_RTOL)
    assert math.isclose(core.density_bulk, 5000.0, rel_tol=LAYER_MASS_RTOL)
    assert math.isclose(mantle.density_bulk, 3000.0, rel_tol=LAYER_MASS_RTOL)


def test_layer_masses_sum_to_planet_mass():
    """Adjacent layers share the interface slice, so the layer masses telescope to the planet mass."""
    world = _two_layer_world()
    world.solve_eos(G_to_use=G, verbose=False)
    assert math.isclose(world.calc_total_mass(), world.planet_mass_eos, rel_tol=1e-14)


def test_resolve_overwrites_layer_mass():
    """Each successful solve sets the layer masses again, so a changed EOS changes them."""
    _, _, ConstantDensityEOS, _ = _import()
    world = _two_layer_world()
    world.solve_eos(G_to_use=G, verbose=False)
    mantle = world.layers[1]
    mass_before = mantle.mass
    mantle.set_eos(ConstantDensityEOS(reference_density=4500.0))
    world.solve_eos(G_to_use=G, verbose=False)
    assert math.isclose(mantle.mass / mass_before, 4500.0 / 3000.0, rel_tol=LAYER_MASS_RTOL)


def test_bundled_world_internal_heating_after_solve():
    """A bundled world has no per-layer masses; its radiogenic heating becomes nonzero once the EOS is solved."""
    from TidalPy.structures_x import build_world
    world = build_world("earth_simple")
    assert world.calc_internal_heating(0.0) == 0.0
    assert world.solve_eos(verbose=False)["success"]
    heating = world.calc_internal_heating(0.0)
    assert heating > 0.0
    expected = sum(layer.calc_radiogenic_heating(0.0, layer.mass)
                   for layer in world.layers if hasattr(layer, "calc_radiogenic_heating"))
    assert math.isclose(heating, expected, rel_tol=1e-12)


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


# =====================================================================================================================
# Binary round trip
# =====================================================================================================================
def test_loaded_world_solves_eos_without_reattaching(tmp_path):
    """A world reloaded from binary keeps its material EOS models, so solve_eos reproduces the original."""
    from TidalPy.structures_x import build_world
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    world = build_world("earth_prem")
    reference = world.solve_eos(verbose=False)
    assert reference["success"]
    path = str(tmp_path / "earth_prem.tpyb")
    world.save_binary(path)

    loaded = LayeredWorld("placeholder", 1.0, 1.0)
    loaded.load_binary(path)
    assert loaded.all_eos_set
    # The binary record is the structure alone: like the tide configuration, the solver settings a data-file world
    # pins on itself (RK45 for its EOS solve) are not in it, so they are restored here to compare like with like.
    assert loaded.get_solver_defaults() == {}
    loaded.set_solver_defaults(**world.get_solver_defaults())
    result = loaded.solve_eos(verbose=False)
    assert result["success"]
    assert math.isclose(result["planet_mass"], reference["planet_mass"], rel_tol=1e-12)
    assert math.isclose(result["planet_moi"], reference["planet_moi"], rel_tol=1e-12)


# =====================================================================================================================
# Non-dimensional solve, central-pressure iteration, and configuration defaults
# =====================================================================================================================
def _compressible_world(radius=6.371e6):
    """Two Birch-Murnaghan layers: the central-pressure iteration has to converge on a compressible planet."""
    from TidalPy.structures_x import build_world
    return build_world({
        "schema_version": "0.2.0", "name": "bm", "type": "terrestrial", "radius_m": radius, "mass_kg": 6.0e24,
        "layers": {
            "core": {"class": "physics", "type": "iron", "layer_index": 0, "radius_fraction": 0.55,
                     "material": {"model": "birch_murnaghan", "reference_density_kg_m3": 8300.0,
                                  "reference_bulk_modulus_pa": 1.6e11, "bulk_modulus_derivative": 5.0}},
            "mantle": {"class": "physics", "type": "mantle_rock", "layer_index": 1, "radius_fraction": 1.0,
                       "material": {"model": "birch_murnaghan", "reference_density_kg_m3": 3300.0,
                                    "reference_bulk_modulus_pa": 1.3e11, "bulk_modulus_derivative": 4.0}}}})


@pytest.fixture
def restore_config_x():
    """Restore ``TidalPy.config_x`` and the C++ solver defaults after a test changes them."""
    import copy
    import TidalPy
    from TidalPy.constants import update_constants_x
    original = copy.deepcopy(TidalPy.config_x)
    yield
    TidalPy.config_x = original
    update_constants_x()


@pytest.mark.parametrize("build", [_two_layer_world, _compressible_world], ids=["constant", "birch_murnaghan"])
def test_nondimensional_and_si_solves_agree(build):
    """The default non-dimensional solve and an SI solve give the same structure at tight tolerances."""
    world = build()
    nondim = world.solve_eos(G_to_use=G, rtol=1.0e-10, atol=1.0e-14, pressure_tol=1.0e-9, nondimensionalize=True)
    radii = np.linspace(0.05, 0.99, 7) * world.radius
    density_nd = np.array([world.get_density(r) for r in radii])
    gravity_nd = np.array([world.get_gravity(r) for r in radii])
    pressure_nd = np.array([world.get_pressure(r) for r in radii])
    si = world.solve_eos(G_to_use=G, rtol=1.0e-10, atol=1.0e-14, pressure_tol=1.0e-9, nondimensionalize=False)
    assert nondim["success"] and si["success"]
    for key in ("planet_mass", "planet_moi", "surface_gravity", "central_pressure"):
        assert math.isclose(nondim[key], si[key], rel_tol=1.0e-8), key
    np.testing.assert_allclose(density_nd, [world.get_density(r) for r in radii], rtol=1.0e-8)
    np.testing.assert_allclose(gravity_nd, [world.get_gravity(r) for r in radii], rtol=1.0e-8)
    np.testing.assert_allclose(pressure_nd, [world.get_pressure(r) for r in radii], rtol=1.0e-7)


def test_secant_iteration_converges_on_a_compressible_planet():
    """The central-pressure iteration converges in a few steps where the unit-slope update crawled."""
    world = _compressible_world()
    result = world.solve_eos(G_to_use=G, rtol=1.0e-10, atol=1.0e-14, pressure_tol=1.0e-9)
    assert result["success"] is True
    assert result["max_iters_hit"] is False
    assert result["iterations"] <= 12
    # The surface-pressure mismatch is below pressure_tol times the central-pressure scale.
    assert result["pressure_error"] < 1.0e-8 * world.central_pressure
    assert abs(result["surface_pressure"]) < 1.0e-8 * world.central_pressure


def _one_layer_bm_world(radius, reference_density, bulk_modulus):
    from TidalPy.structures_x import build_world
    return build_world({
        "schema_version": "0.2.0", "name": "bm1", "type": "terrestrial", "radius_m": radius, "mass_kg": 6.0e24,
        "layers": {"mantle": {"class": "physics", "type": "mantle_rock", "layer_index": 0, "radius_fraction": 1.0,
                               "material": {"model": "birch_murnaghan",
                                            "reference_density_kg_m3": reference_density,
                                            "reference_bulk_modulus_pa": bulk_modulus,
                                            "bulk_modulus_derivative": 4.0}}}})


@pytest.mark.parametrize("radius, reference_density, bulk_modulus, max_passes", [
    (2.0e7, 3300.0, 1.3e11, 16),   # 21 passes when the step crawled by the residual
    (6.4e6, 5000.0, 1.0e10, 30),   # 274 passes
    (1.2e7, 3000.0, 5.0e9, 60),    # stopped at the 300-pass cap
])
def test_secant_iteration_does_not_crawl_where_the_surface_pressure_first_falls(
        radius, reference_density, bulk_modulus, max_passes):
    """Where the surface pressure first falls as the central pressure rises, the step grows geometrically until the
    mismatch changes sign, then the root is bracketed, instead of stepping one residual's worth of pressure per pass."""
    world = _one_layer_bm_world(radius, reference_density, bulk_modulus)
    result = world.solve_eos(G_to_use=G, max_iters=300)
    assert result["success"] is True, result["message"]
    assert result["iterations"] <= max_passes
    assert abs(result["surface_pressure"]) < 1.0e-7 * world.central_pressure


def test_max_iters_hit_is_reported():
    """Stopping at the iteration cap off the target surface pressure is a failure: the structure is not
    hydrostatic, so the world stays unsolved and the message says why."""
    world = _compressible_world()
    result = world.solve_eos(G_to_use=G, pressure_tol=1.0e-12, max_iters=1)
    assert result["success"] is False
    assert result["max_iters_hit"] is True
    assert result["iterations"] == 1
    assert "no hydrostatic structure" in result["message"]
    assert world.eos_solved is False


def test_pressure_tolerance_is_relative_to_the_central_pressure():
    """A target surface pressure is met to pressure_tol of the central-pressure scale."""
    world = _two_layer_world()
    target = 1.0e5
    result = world.solve_eos(G_to_use=G, surface_pressure=target, rtol=1.0e-10, atol=1.0e-14, pressure_tol=1.0e-9)
    assert result["success"] is True and result["max_iters_hit"] is False
    assert abs(result["surface_pressure"] - target) < 1.0e-8 * world.central_pressure


def test_eos_solver_defaults_come_from_the_config(restore_config_x):
    """Arguments left as None take the [eos_solver] values; an explicit argument still wins."""
    import TidalPy
    world = _compressible_world()
    default = world.solve_eos(G_to_use=G)
    assert default["max_iters_hit"] is False
    TidalPy.reinit(provided_config_x={"eos_solver": {"max_iters": 1, "pressure_tol": 1.0e-12}})
    capped = world.solve_eos(G_to_use=G)
    assert capped["max_iters_hit"] is True and capped["iterations"] == 1
    explicit = world.solve_eos(G_to_use=G, max_iters=100, pressure_tol=1.0e-5)
    assert explicit["max_iters_hit"] is False
