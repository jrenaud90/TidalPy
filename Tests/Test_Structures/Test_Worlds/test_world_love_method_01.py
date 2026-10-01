"""Love-number method dispatch on BaseWorld (radial solver, propagation matrix, homogeneous, cpl, ctl)."""
import cmath
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Material import Material, Phase
from TidalPy.Rheology import Elastic, Maxwell
from TidalPy.Structures import build_world
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Tides.classes import make_tide
from TidalPy.Tides.love import calc_homogeneous_love_numbers

RADIUS = 6.0e6
DENSITY = 4000.0
MASS = (4.0 / 3.0) * math.pi * RADIUS**3 * DENSITY
SHEAR = 6.0e10
BULK = 1.3e11
VISCOSITY = 1.0e19
FREQ = 1.0e-5


def _material(shear, viscosity, bulk):
    return Material(solid=Phase(
        eos={"model": "constant", "reference_density_kg_m3": DENSITY, "bulk_modulus_pa": bulk},
        shear_modulus={"model": "constant", "shear_modulus_pa": shear},
        shear_viscosity={"model": "constant", "reference_viscosity_pas": viscosity},
        bulk_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e30}))


def _layer(name, index, r_inner, r_outer, mass, shear=SHEAR, use_tides=True, rheology=None, viscosity=VISCOSITY,
           incompressible=True, bulk=BULK):
    return Layer(
        name,
        index,
        r_inner,
        r_outer,
        mass,
        _material(shear, viscosity, bulk),
        use_tides=use_tides,
        is_incompressible=incompressible,
        shear_rheology=rheology if rheology is not None else Maxwell(),
        bulk_rheology=Elastic())


def _uniform_world(rheology=None, viscosity=VISCOSITY, incompressible=True, bulk=BULK):
    """Single solid static uniform sphere (shooting cannot start a static incompressible solid: use _stiff_world)."""
    world = BaseWorld("uniform", RADIUS, MASS)
    world.add_layer(_layer("mantle", 0, 0.0, RADIUS, MASS, rheology=rheology, viscosity=viscosity,
                           incompressible=incompressible, bulk=bulk))
    world.solve_eos()
    return world


def _stiff_world(rheology=None, viscosity=VISCOSITY):
    """Compressible uniform sphere, incompressible to about 1e-4, usable by the shooting solver."""
    return _uniform_world(rheology=rheology, viscosity=viscosity, incompressible=False, bulk=1.0e15)


def _two_layer_world(core_tidal):
    """Uniform density, two layers with different shear moduli; the core is optionally non-tidal."""
    r_core = 0.5 * RADIUS
    core_mass = MASS * (r_core / RADIUS) ** 3
    world = BaseWorld("two_layer", RADIUS, MASS)
    world.add_layer(_layer("core", 0, 0.0, r_core, core_mass, shear=3.0 * SHEAR, use_tides=core_tidal,
                           rheology=Elastic()))
    world.add_layer(_layer("mantle", 1, r_core, RADIUS, MASS - core_mass, rheology=Elastic()))
    world.solve_eos()
    return world


def _reference_love(world, shear, degree_l=2):
    """Closed-form homogeneous Love numbers at the world's bulk density and EOS surface gravity."""
    density_bulk = world.mass / ((4.0 / 3.0) * math.pi * world.radius**3)
    gravity = float(world.get_gravity(world.radius))
    return calc_homogeneous_love_numbers(shear, density_bulk, gravity, world.radius, degree_l)


@pytest.mark.parametrize("alias, canonical", (
    ("radial_solver", "radial_solver"), ("shooting", "radial_solver"), ("pm", "propagation_matrix"),
    ("homogen", "homogeneous"), ("CPL", "cpl"), ("ctl", "ctl"),
))
def test_method_aliases_and_result_key(alias, canonical):
    """Method aliases resolve to the canonical name reported by the result and the world."""
    world = _uniform_world() if canonical == "propagation_matrix" else _stiff_world()
    result = world.solve_love_numbers(frequency=FREQ, love_method=alias, fixed_q=50.0, fixed_dt=100.0)
    assert result["success"], result["message"]
    assert result["love_method"] == canonical
    assert world.love_method == canonical
    assert world.love_success and world.love_solved


def test_default_method_is_radial_solver():
    """The default method is the radial solver, which leaves the analytic-only getters NaN."""
    world = _stiff_world()
    assert world.love_method == "radial_solver"
    result = world.solve_love_numbers(frequency=FREQ)
    assert result["success"] and result["love_method"] == "radial_solver"
    assert math.isnan(world.love_effective_shear_modulus.real)
    assert math.isnan(world.love_tidal_volume)


def test_unknown_and_reserved_methods():
    """Unknown, reserved, and unsupported-for-supplied methods raise."""
    world = _uniform_world()
    with pytest.raises(ValueError, match="unknown Love-number method"):
        world.solve_love_numbers(frequency=FREQ, love_method="shoot")
    with pytest.raises(NotImplementedError, match="laterally_inhomogeneous"):
        world.solve_love_numbers(frequency=FREQ, love_method="3d")
    with pytest.raises(ValueError, match="only the radial_solver"):
        world.solve_love_numbers_supplied(
            np.full(5, SHEAR + 0j), np.full(5, BULK + 0j), np.linspace(0.0, RADIUS, 5), frequency=FREQ,
            love_method="homogeneous")


@pytest.mark.parametrize("degree_l", (2, 3))
def test_homogeneous_matches_radial_solvers_for_uniform_sphere(degree_l):
    """Homogeneous, propagation matrix, and shooting agree with the closed form on a uniform Maxwell sphere."""
    world = _uniform_world()
    analytic = world.solve_love_numbers(frequency=FREQ, degree_l=degree_l, love_method="homogeneous")
    matrix = world.solve_love_numbers(frequency=FREQ, degree_l=degree_l, love_method="propagation_matrix")
    assert matrix["success"], matrix["message"]
    stiff = _stiff_world()
    shooting = stiff.solve_love_numbers(frequency=FREQ, degree_l=degree_l, love_method="radial_solver",
                                        rtol=1e-9, atol=1e-12)
    assert shooting["success"], shooting["message"]
    for key in ("love_number_k", "love_number_h", "love_number_l"):
        np.testing.assert_allclose(analytic[key], matrix[key], rtol=1e-8)
        np.testing.assert_allclose(analytic[key], shooting[key], rtol=1e-3)   # bulk 1e15 Pa: ~1e-4 compressibility
    mu = Maxwell().calc_complex_modulus(SHEAR, VISCOSITY, FREQ)
    world.solve_love_numbers(frequency=FREQ, degree_l=degree_l, love_method="homogeneous")
    assert world.love_effective_shear_modulus == pytest.approx(mu, rel=1e-12)
    assert world.love_tidal_volume == pytest.approx((4.0 / 3.0) * math.pi * RADIUS**3, rel=1e-12)
    reference = _reference_love(world, mu, degree_l)
    assert analytic["love_number_k"] == pytest.approx(reference.k, rel=1e-12)
    assert analytic["love_number_h"] == pytest.approx(reference.h, rel=1e-12)
    assert analytic["love_number_l"] == pytest.approx(reference.l, rel=1e-12)


def test_analytic_getters_do_not_leak_radial_results():
    """After an analytic solve the radial-only getters do not report the previous radial solution."""
    world = _stiff_world()
    assert world.solve_love_numbers(frequency=FREQ, love_method="radial_solver")["success"]
    assert world.love_surface_amplification > 0.0
    world.solve_love_numbers(frequency=FREQ, love_method="homogeneous")
    assert world.love_surface_amplification == 0.0
    assert world.love_num_ytypes == 1
    assert math.isnan(world.get_love_surface_y(0, 0).real)
    world.solve_love_numbers(frequency=FREQ, love_method="radial_solver")
    assert not math.isnan(world.get_love_surface_y(0, 0).real)


def test_non_tidal_layers_take_no_part():
    """Only use_tides layers enter the homogeneous solve, each weighted by its tidal scale."""
    included = _two_layer_world(core_tidal=True)
    excluded = _two_layer_world(core_tidal=False)
    result_in = included.solve_love_numbers(frequency=FREQ, love_method="homogeneous")
    result_ex = excluded.solve_love_numbers(frequency=FREQ, love_method="homogeneous")
    core_volume = (4.0 / 3.0) * math.pi * (0.5 * RADIUS) ** 3
    total_volume = (4.0 / 3.0) * math.pi * RADIUS**3
    core_fraction = core_volume / total_volume
    mu_scaled = 3.0 * SHEAR * core_fraction + SHEAR * (1.0 - core_fraction)
    assert included.love_effective_shear_modulus == pytest.approx(mu_scaled, rel=1e-10)
    assert included.love_tidal_volume == pytest.approx(total_volume, rel=1e-12)
    assert excluded.love_effective_shear_modulus == pytest.approx(SHEAR, rel=1e-12)
    assert excluded.love_tidal_volume == pytest.approx(total_volume - core_volume, rel=1e-12)
    k_mantle = _reference_love(excluded, SHEAR).k
    k_core = _reference_love(included, 3.0 * SHEAR).k
    assert result_ex["love_number_k"] == pytest.approx((1.0 - core_fraction) * k_mantle, rel=1e-12)
    assert result_in["love_number_k"] == pytest.approx(
        core_fraction * k_core + (1.0 - core_fraction) * k_mantle, rel=1e-12)


def test_cpl_and_ctl_structure():
    """cpl and ctl apply their constant Q and time-lag phase lags to the elastic Love numbers."""
    world = _uniform_world(rheology=Elastic())
    static = world.solve_love_numbers(frequency=FREQ, love_method="homogeneous")
    assert static["love_number_k"].imag == 0.0
    cpl = world.solve_love_numbers(frequency=FREQ, love_method="cpl", fixed_q=50.0)
    assert cpl["love_number_k"] == pytest.approx(static["love_number_k"] * (1.0 - 1.0j / 50.0), rel=1e-12)
    assert cpl["love_number_h"] == pytest.approx(static["love_number_h"] * (1.0 - 1.0j / 50.0), rel=1e-12)
    assert world.love_effective_shear_modulus == pytest.approx(SHEAR)
    ctl = world.solve_love_numbers(frequency=FREQ, love_method="ctl", fixed_dt=600.0)
    assert ctl["love_number_k"] == pytest.approx(static["love_number_k"] * (1.0 - 1.0j * FREQ * 600.0), rel=1e-12)


def test_cpl_uses_static_modulus():
    """cpl and ctl average the unrelaxed modulus even for a viscous rheology."""
    world = _uniform_world(viscosity=1.0e15)   # strongly viscoelastic at FREQ
    viscoelastic = world.solve_love_numbers(frequency=FREQ, love_method="homogeneous")
    assert abs(viscoelastic["love_number_k"].imag) > 1e-3
    cpl = world.solve_love_numbers(frequency=FREQ, love_method="cpl", fixed_q=100.0)
    elastic = _reference_love(world, SHEAR)
    assert cpl["love_number_k"] == pytest.approx(elastic.k * (1.0 - 1.0j / 100.0), rel=1e-12)


def test_cpl_ctl_parameters_from_tide_model():
    """cpl and ctl take Q and the time lag from the tide model when not given explicitly."""
    world = _uniform_world(rheology=Elastic())
    with pytest.raises(ValueError, match="fixed_q"):
        world.solve_love_numbers(frequency=FREQ, love_method="cpl")
    with pytest.raises(ValueError, match="fixed_dt"):
        world.solve_love_numbers(frequency=FREQ, love_method="ctl")
    world.set_tide_model(make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [40.0]}))
    from_model = world.solve_love_numbers(frequency=FREQ, love_method="cpl")
    explicit = world.solve_love_numbers(frequency=FREQ, love_method="cpl", fixed_q=40.0)
    assert from_model["love_number_k"] == explicit["love_number_k"]
    with pytest.raises(ValueError, match="degree 3"):
        world.solve_love_numbers(frequency=FREQ, degree_l=3, love_method="cpl")   # no Q_3 on the model
    world.set_tide_model(make_tide("ctl", {"fixed_k": [0.3], "fixed_dt_s": [300.0]}))
    from_model = world.solve_love_numbers(frequency=FREQ, love_method="ctl")
    assert from_model["love_number_k"] == world.solve_love_numbers(
        frequency=FREQ, love_method="ctl", fixed_dt=300.0)["love_number_k"]


def test_tide_config_love_method():
    """set_tide_config stores the default Love method and its parameters."""
    world = _uniform_world(rheology=Elastic())
    assert world.get_tide_config()["love_method"] == "radial_solver"
    assert "love_fixed_q" not in world.get_tide_config()
    world.set_tide_config(love_method="cpl", love_fixed_q=30.0)
    cfg = world.get_tide_config()
    assert cfg["love_method"] == "cpl" and cfg["love_fixed_q"] == 30.0 and "love_fixed_dt_s" not in cfg
    with pytest.raises(ValueError, match="unknown Love-number method"):
        world.set_tide_config(love_method="nope")


def test_calc_tides_uses_configured_method():
    """With the rheology tide model, calc_tides takes k_l from the configured Love method."""
    world = _stiff_world()
    world.set_tide_model(make_tide("rheology"))
    orbit = dict(orbital_frequency=2.0e-5, spin_frequency=1.0e-5, eccentricity=0.01, obliquity=0.0,
                 semi_major_axis=4.0e8, host_mass=1.9e27)
    world.set_tide_config(max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
    world.calc_tides(**orbit)
    heating_radial = world.get_tidal_heating()
    world.set_tide_config(max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0, love_method="homogeneous")
    world.calc_tides(**orbit)
    heating_analytic = world.get_tidal_heating()
    assert world.get_tide_config()["love_method"] == "homogeneous"
    assert heating_analytic > 0.0
    np.testing.assert_allclose(heating_analytic, heating_radial, rtol=1e-3)   # same physics; bulk = 1e15 Pa
    # The analytic methods have no depth-resolved solution for the 3D path.
    with pytest.raises(RuntimeError, match="homogeneous"):
        world.get_3d_tidal_heating(radius=0.9 * RADIUS, colatitude=1.0, **orbit)


def test_build_world_round_trips_love_method():
    """The Love method and its Q survive a config-dict round trip through build_world."""
    world = build_world("earth_simple")
    cfg = world.get_config_dict()
    cfg["tides"]["love_method"] = "cpl"
    cfg["tides"]["love_fixed_q"] = 25.0
    rebuilt = build_world(cfg)
    tide_cfg = rebuilt.get_tide_config()
    assert tide_cfg["love_method"] == "cpl" and tide_cfg["love_fixed_q"] == 25.0
    assert rebuilt.get_config_dict()["tides"]["love_method"] == "cpl"
    rebuilt.solve_eos()
    result = rebuilt.solve_love_numbers(frequency=FREQ, love_method="cpl")
    assert result["success"] and cmath.isfinite(result["love_number_k"])
    assert -result["love_number_k"].imag == pytest.approx(result["love_number_k"].real / 25.0)
