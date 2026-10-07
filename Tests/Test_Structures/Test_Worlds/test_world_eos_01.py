"""The world equation-of-state solve (``BaseWorld.solve_eos``) and the per-layer material wiring (``material``)."""
import copy
import math
import re
import warnings
from pathlib import Path

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Material import Material, Phase
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.base import BaseWorld


_PLANET_RADIUS = 6000.0e3      # [m]
_DENSITY       = 3500.0        # [kg/m^3]
_R_CMB         = _PLANET_RADIUS / 2.0
LAYER_MASS_RTOL = 1e-12  # Constant-density layer mass versus rho V after the solve.


def _material(density):
    """A constant-density material with only an equation of state."""
    return Material(solid=Phase(eos={"model": "constant", "reference_density_kg_m3": density}))


def _uniform_world():
    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = BaseWorld("Uniform", _PLANET_RADIUS, mass, world_type="terrestrial")
    world.add_layer(Layer("mantle", 0, 0.0, _PLANET_RADIUS, mass, _material(_DENSITY)))
    return world


def _two_layer_world(rho_core=5000.0, rho_mantle=3000.0):
    """Two constant-density layers built without layer masses (as the TOML builder does)."""
    mass = (4.0 / 3.0) * math.pi * (rho_core * _R_CMB ** 3 + rho_mantle * (_PLANET_RADIUS ** 3 - _R_CMB ** 3))
    world = BaseWorld("TwoLayer", _PLANET_RADIUS, mass)
    core = Layer("core", 0, 0.0, _R_CMB, 0.0, _material(rho_core))
    mantle = Layer("mantle", 1, _R_CMB, _PLANET_RADIUS, 0.0, _material(rho_mantle))
    world.add_layer(core)
    world.add_layer(mantle)
    return world


@pytest.fixture(scope="module")
def solved_uniform():
    world = _uniform_world()
    world.solve_eos(G_to_use=G, verbose=False)
    return world


# =====================================================================================================================
# Material wiring
# =====================================================================================================================
def test_material_set_flags():
    layer = Layer("rock", 0, 0.0, 1.0e6, 1.0e20)
    assert layer.material_set is False
    layer.material = _material(_DENSITY)
    assert layer.material_set is True


def test_solve_requires_all_materials_set():
    world = BaseWorld("NoMaterial", _PLANET_RADIUS, 1.0e24)
    world.add_layer(Layer("mantle", 0, 0.0, _PLANET_RADIUS, 1.0e24))
    assert world.all_materials_set is False
    with pytest.raises(ValueError):
        world.solve_eos(verbose=False)


def test_solve_requires_layers():
    world = BaseWorld("Empty", _PLANET_RADIUS, 1.0e24)
    with pytest.raises(ValueError):
        world.solve_eos(verbose=False)


# =====================================================================================================================
# Uniform planet (constant-density material)
# =====================================================================================================================
def test_uniform_solve_succeeds():
    world = _uniform_world()
    result = world.solve_eos(surface_pressure=0.0, G_to_use=G, verbose=False)
    assert result["success"] is True
    assert result["iterations"] >= 1
    assert world.eos_solved is True
    assert result["gravity"].shape == result["density"].shape


def test_uniform_gravity_analytic(solved_uniform):
    """g(r) = (4/3) pi G rho r for a uniform sphere."""
    for r in np.linspace(0.0, _PLANET_RADIUS, 50)[1:]:
        expected = (4.0 / 3.0) * math.pi * G * _DENSITY * r
        assert math.isclose(solved_uniform.get_gravity(r), expected, rel_tol=0.05)


def test_uniform_density_recovered(solved_uniform):
    """The world density query, delegated to the owning layer's profile, recovers the layer density."""
    for r in np.linspace(0.0, _PLANET_RADIUS, 25):
        assert math.isclose(solved_uniform.get_density(r), _DENSITY, rel_tol=0.05)


def test_uniform_pressure_monotonic_and_central(solved_uniform):
    """Pressure decreases outward from the analytic central value (2/3) pi G rho^2 R^2 to about 0."""
    world = solved_uniform
    radii = np.linspace(0.0, _PLANET_RADIUS, 60)
    pressures = np.array([world.get_pressure(r) for r in radii])
    assert np.all(np.diff(pressures) <= 1.0e-3 * abs(pressures[0]))
    p_central_analytic = (2.0 / 3.0) * math.pi * G * _DENSITY ** 2 * _PLANET_RADIUS ** 2
    assert math.isclose(world.get_pressure(0.0), p_central_analytic, rel_tol=0.05)
    assert abs(world.get_pressure(_PLANET_RADIUS)) < 0.01 * world.get_pressure(0.0)


def test_uniform_planet_mass(solved_uniform):
    expected_mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    assert math.isclose(solved_uniform.planet_mass_eos, expected_mass, rel_tol=0.02)


def test_unsolved_world_returns_nan():
    world = _uniform_world()
    assert math.isnan(world.get_density(_PLANET_RADIUS * 0.5))
    assert math.isnan(world.surface_gravity_eos)


def test_unsupported_method_raises():
    world = _uniform_world()
    with pytest.raises(ValueError):
        world.solve_eos(G_to_use=G, integration_method="NOPE", verbose=False)


# =====================================================================================================================
# Two-layer planet: densities, layer masses, and bulk densities
# =====================================================================================================================
def test_two_layer_constant_density():
    world = _two_layer_world()
    result = world.solve_eos(G_to_use=G, verbose=False)
    assert result["success"]
    assert world.central_pressure > 0.0
    assert math.isclose(world.get_density(_R_CMB * 0.5), 5000.0, rel_tol=0.05)
    assert math.isclose(world.get_density(_R_CMB + (_PLANET_RADIUS - _R_CMB) * 0.5), 3000.0, rel_tol=0.05)
    expected_mass = (4.0 / 3.0) * math.pi * (5000.0 * _R_CMB ** 3 + 3000.0 * (_PLANET_RADIUS ** 3 - _R_CMB ** 3))
    assert math.isclose(world.planet_mass_eos, expected_mass, rel_tol=0.03)


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
    """Each successful solve sets the layer masses again, so a changed material changes them."""
    world = _two_layer_world()
    world.solve_eos(G_to_use=G, verbose=False)
    mantle = world.layers[1]
    mass_before = mantle.mass
    mantle.material = _material(4500.0)
    world.solve_eos(G_to_use=G, verbose=False)
    assert math.isclose(mantle.mass / mass_before, 4500.0 / 3000.0, rel_tol=LAYER_MASS_RTOL)


def test_bundled_world_internal_heating_after_solve():
    """A bundled world's radiogenic heating is zero until the EOS solve sets its layer masses."""
    from TidalPy.Structures import build_world
    # The bundled Earth names no radiogenic source, so its mantle is given the chondritic isotopes.
    config = copy.deepcopy(build_world("earth_simple").source_config)
    config["layers"]["mantle"]["radiogenics"] = {"model": "isotope", "isotopes": "modern_day_chondritic"}
    world = build_world(config)
    assert world.calc_internal_heating(0.0) == 0.0
    assert world.solve_eos(verbose=False)["success"]
    heating = world.calc_internal_heating(0.0)
    assert heating > 0.0
    expected = sum(layer.calc_radiogenic_heating(0.0, layer.mass) for layer in world.layers)
    assert math.isclose(heating, expected, rel_tol=1e-12)


# =====================================================================================================================
# PREM Earth (interpolated material)
# =====================================================================================================================
def test_prem_earth_interpolated():
    prem_dir = Path(__file__).resolve().parents[2] / "Test_Material"

    prem_data = []
    try:
        for layer_i in range(3):
            prem_data.append(np.loadtxt(prem_dir / f"prem_layer{layer_i}.txt", delimiter=","))
    except Exception as e:  # noqa: BLE001
        warnings.warn(f"Could not load PREM data: {e}")
        pytest.skip("Could not load PREM Earth data.")

    surface_radius = prem_data[2][:, 0][-1]
    world = BaseWorld("PREM-Earth", surface_radius, 5.972e24, world_type="terrestrial")

    prev_outer = 0.0
    for layer_i in range(3):
        radius_array  = np.ascontiguousarray(prem_data[layer_i][:, 0])
        density_array = np.ascontiguousarray(prem_data[layer_i][:, 1])
        r_outer = radius_array[-1]
        material = Material(solid=Phase(eos={
            "model": "interpolate", "radius_m": radius_array.tolist(), "density_kg_m3": density_array.tolist()}))
        world.add_layer(Layer(f"layer{layer_i}", layer_i, prev_outer, r_outer, 0.0, material))
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

    mid_mantle = 0.5 * (prem_data[2][:, 0][0] + prem_data[2][:, 0][-1])
    expected = np.interp(mid_mantle, prem_data[2][:, 0], prem_data[2][:, 1])
    assert math.isclose(world.get_density(mid_mantle), expected, rel_tol=0.05)


def test_loaded_world_solves_eos_without_reattaching(tmp_path):
    """A world reloaded from binary keeps its materials and pinned solver settings, and reproduces its solve."""
    from TidalPy.Structures import build_world
    world = build_world("earth_prem")
    reference = world.solve_eos(verbose=False)
    assert reference["success"]
    path = str(tmp_path / "earth_prem.tpyb")
    world.save_binary(path)

    loaded = type(world)("placeholder", 1.0, 1.0)
    loaded.load_binary(path)
    assert loaded.all_materials_set
    assert loaded.get_solver_defaults() == world.get_solver_defaults()
    result = loaded.solve_eos(verbose=False)
    assert result["success"]
    assert math.isclose(result["planet_mass"], reference["planet_mass"], rel_tol=1e-12)
    assert math.isclose(result["planet_moi"], reference["planet_moi"], rel_tol=1e-12)


# =====================================================================================================================
# Non-dimensional solve, central-pressure iteration, and configuration defaults
# =====================================================================================================================
def _bm_solid(reference_density, bulk_modulus, bulk_modulus_derivative, shear_modulus, shear_viscosity):
    """A solid phase table on a Birch-Murnaghan law, with a constant shear modulus and the given viscosity law."""
    return {"eos": {"model": "birch_murnaghan", "reference_density_kg_m3": reference_density,
                    "reference_bulk_modulus_pa": bulk_modulus, "bulk_modulus_derivative": bulk_modulus_derivative},
            "shear_modulus": {"model": "constant", "shear_modulus_pa": shear_modulus},
            "shear_viscosity": shear_viscosity,
            "bulk_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e22}}


# The viscosity law of the silicate mantle layers: 1e22 Pa s at 1000 K, with an activation energy of 300 kJ/mol.
_ROCK_VISCOSITY = {"model": "reference", "reference_viscosity_pas": 1.0e22, "reference_temperature_k": 1000.0,
                   "molar_activation_energy_j_mol": 3.0e5, "molar_activation_volume_m3_mol": 0.0}


def _compressible_world(radius=6.371e6):
    """Two Birch-Murnaghan layers, so the central-pressure iteration has real work to do."""
    from TidalPy.Structures import build_world
    return build_world({
        "schema_version": "0.2.0", "name": "bm", "type": "terrestrial", "radius_m": radius, "mass_kg": 6.0e24,
        "layers": {
            "core": {"layer_index": 0, "radius_fraction": 0.55,
                     "material": {"solid": _bm_solid(8300.0, 1.6e11, 5.0, 5.25e10,
                                                     {"model": "constant", "reference_viscosity_pas": 1.0e20})},
                     "shear_rheology": {"model": "maxwell"}},
            "mantle": {"layer_index": 1, "radius_fraction": 1.0,
                       "material": {"solid": _bm_solid(3300.0, 1.3e11, 4.0, 6.0e10, _ROCK_VISCOSITY)},
                       "shear_rheology": {"model": "andrade", "alpha": 0.3, "zeta": 1.0}}}})


def _one_layer_bm_world(radius, reference_density, bulk_modulus):
    from TidalPy.Structures import build_world
    return build_world({
        "schema_version": "0.2.0", "name": "bm1", "type": "terrestrial", "radius_m": radius, "mass_kg": 6.0e24,
        "layers": {"mantle": {"layer_index": 0, "radius_fraction": 1.0,
                              "material": {"solid": _bm_solid(reference_density, bulk_modulus, 4.0, 6.0e10,
                                                              _ROCK_VISCOSITY)},
                              "shear_rheology": {"model": "andrade", "alpha": 0.3, "zeta": 1.0}}}})


@pytest.fixture
def without_mass_check():
    """Lift `[numerical] maximum_eos_mass_ratio` for a test whose world carries a placeholder mass."""
    import TidalPy
    import TidalPy.constants
    numerical = TidalPy.config["numerical"]
    original = numerical["maximum_eos_mass_ratio"]
    numerical["maximum_eos_mass_ratio"] = math.inf
    TidalPy.constants.update_constants()
    yield
    numerical["maximum_eos_mass_ratio"] = original
    TidalPy.constants.update_constants()


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
    """The central-pressure iteration converges in a few steps to within pressure_tol of the central pressure."""
    world = _compressible_world()
    result = world.solve_eos(G_to_use=G, rtol=1.0e-10, atol=1.0e-14, pressure_tol=1.0e-9)
    assert result["success"] is True
    assert result["max_iters_hit"] is False
    assert result["iterations"] <= 12
    assert result["pressure_error"] < 1.0e-8 * world.central_pressure
    assert abs(result["surface_pressure"]) < 1.0e-8 * world.central_pressure


@pytest.mark.parametrize("radius, reference_density, bulk_modulus, max_passes", [
    (2.0e7, 3300.0, 1.3e11, 16),   # 21 passes when the step crawled by the residual
    (6.4e6, 5000.0, 1.0e10, 30),   # 274 passes
    (1.2e7, 3000.0, 5.0e9, 60),    # stopped at the 300-pass cap
])
def test_secant_iteration_does_not_crawl_where_the_surface_pressure_first_falls(
        radius,
        reference_density,
        bulk_modulus,
        max_passes,
        without_mass_check):
    """Where the surface pressure first falls with central pressure, the step grows until the root is bracketed."""
    # These bodies carry a placeholder mass (solved masses are 30 to 1400 times it), so the mass check is lifted.
    world = _one_layer_bm_world(radius, reference_density, bulk_modulus)
    result = world.solve_eos(G_to_use=G, max_iters=300)
    assert result["success"] is True, result["message"]
    assert result["iterations"] <= max_passes
    assert abs(result["surface_pressure"]) < 1.0e-7 * world.central_pressure


def test_max_iters_hit_is_reported():
    """Stopping at the iteration cap off the target surface pressure is a reported failure."""
    world = _compressible_world()
    result = world.solve_eos(G_to_use=G, pressure_tol=1.0e-12, max_iters=1)
    assert result["success"] is False
    assert result["max_iters_hit"] is True
    assert result["iterations"] == 1
    assert "no hydrostatic structure" in result["message"]
    # The mismatch is reported in Pa and against the central-pressure scale, never rounded to zero.
    match = re.search(r"misses its target by (\S+) Pa, (\S+) of the central-pressure scale \((\S+) Pa\)",
                      result["message"])
    assert match is not None, result["message"]
    mismatch, relative, scale = (float(match.group(index)) for index in (1, 2, 3))
    assert mismatch > 0.0
    assert relative > 1.0e-12
    assert mismatch == pytest.approx(relative * scale, rel=1.0e-2)
    assert 1.0e10 < scale < 1.0e13   # [Pa] an Earth-sized body's central-pressure scale
    assert world.eos_solved is False


def test_pressure_tolerance_is_relative_to_the_central_pressure():
    """A target surface pressure is met to pressure_tol of the central-pressure scale."""
    world = _two_layer_world()
    target = 1.0e5
    result = world.solve_eos(G_to_use=G, surface_pressure=target, rtol=1.0e-10, atol=1.0e-14, pressure_tol=1.0e-9)
    assert result["success"] is True and result["max_iters_hit"] is False
    assert abs(result["surface_pressure"] - target) < 1.0e-8 * world.central_pressure


def test_eos_solver_defaults_come_from_the_config(restore_config):
    """Arguments left as None take the [eos_solver] values; an explicit argument still wins."""
    import TidalPy
    world = _compressible_world()
    default = world.solve_eos(G_to_use=G)
    assert default["max_iters_hit"] is False
    TidalPy.reinit(provided_config={"eos_solver": {"max_iters": 1, "pressure_tol": 1.0e-12}})
    capped = world.solve_eos(G_to_use=G)
    assert capped["max_iters_hit"] is True and capped["iterations"] == 1
    explicit = world.solve_eos(G_to_use=G, max_iters=100, pressure_tol=1.0e-5)
    assert explicit["max_iters_hit"] is False


# =====================================================================================================================
# Structure integrations, and densities that do not depend on pressure
# =====================================================================================================================
def _melting_material(density):
    """A constant-density solid and liquid melting between fixed temperatures."""
    return Material(
        solid=Phase(eos={"model": "constant", "reference_density_kg_m3": density},
                    shear_modulus={"model": "constant", "shear_modulus_pa": 6.0e10},
                    shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e21}),
        liquid=Phase(eos={"model": "constant", "reference_density_kg_m3": 0.85 * density},
                     shear_viscosity={"model": "constant", "reference_viscosity_pas": 0.2}),
        solidus={"model": "constant", "temperature_k": 1600.0},
        liquidus={"model": "constant", "temperature_k": 2000.0})


@pytest.mark.parametrize("build", [_uniform_world, _two_layer_world], ids=["uniform", "two_layer"])
def test_a_density_independent_of_pressure_solves_in_two_integrations(build):
    """With every density set by radius alone, the surface pressure falls one for one with the central pressure, so
    the first unit-slope step lands on the root to rounding and the integration after it keeps its output."""
    result = build().solve_eos(G_to_use=G)
    assert result["success"] is True
    assert result["structure_integrations"] == 2
    assert result["pressure_error"] < 1.0e-12 * result["central_pressure"]


def test_a_changed_world_with_densities_independent_of_pressure_resolves_in_two_integrations():
    """A re-solve starts from the last central pressure and, when nothing changed, keeps its first integration; after a
    change its unit-slope step is exact as from scratch."""
    world = _two_layer_world()
    assert world.solve_eos(G_to_use=G)["success"] is True
    assert world.solve_eos(G_to_use=G)["structure_integrations"] == 1
    world.mantle.material = _material(3100.0)
    result = world.solve_eos(G_to_use=G)
    assert result["success"] is True
    assert result["structure_integrations"] == 2


def test_a_compressible_world_repeats_a_converged_integration_that_kept_no_output():
    """Where the density follows the pressure, the iteration measures the slope, and an integration that converges
    without the dense output it was not expected to need is run again to keep it."""
    result = _compressible_world().solve_eos(G_to_use=G)
    assert result["success"] is True
    assert result["structure_integrations"] >= 3
    assert result["structure_integrations"] - result["iterations"] in (0, 1)


@pytest.mark.parametrize("use_melt_density, use_pressure_melting, integrations", [
    (False, False, 2),
    (True, False, 2),
    (True, True, 3),
])
def test_a_layer_mixing_its_melt_density_under_pressure_melting_takes_the_general_start(
        use_melt_density,
        use_pressure_melting,
        integrations):
    """A layer that mixes its melt's density into its own is treated as depending on the pressure when its melting
    range moves with pressure, since its melt fraction then can; at a fixed melting range it is not."""
    mass = (4.0 / 3.0) * math.pi * (5000.0 * _R_CMB ** 3 + 3000.0 * (_PLANET_RADIUS ** 3 - _R_CMB ** 3))
    world = BaseWorld("Melting", _PLANET_RADIUS, mass)
    world.add_layer(Layer("core", 0, 0.0, _R_CMB, 0.0, _material(5000.0)))
    world.add_layer(Layer("mantle", 1, _R_CMB, _PLANET_RADIUS, 0.0, _melting_material(3000.0), temperature=1800.0,
                          use_melting=True, use_melt_density=use_melt_density,
                          use_pressure_melting=use_pressure_melting))
    result = world.solve_eos(G_to_use=G, solve_temperature=False)
    assert result["success"] is True, result["message"]
    assert result["structure_integrations"] == integrations


def test_a_thermal_solve_counts_the_integrations_of_every_pass():
    """structure_integrations sums every thermal pass's integrations; iterations is the last pass's alone."""
    from TidalPy.Structures import build_world
    result = build_world("earth_thermal").solve_eos(solve_temperature=True, surface_temperature=300.0)
    assert result["success"] is True
    assert result["thermal_passes"] >= 1
    assert result["structure_integrations"] >= result["iterations"] + result["thermal_passes"]
