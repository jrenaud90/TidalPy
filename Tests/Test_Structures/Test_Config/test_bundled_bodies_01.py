"""Bundled solar-system, TRAPPIST-1, and gas-giant worlds reproduce their published mass, C/MR2, and tidal observables.

Expected values come from the TOML comments, which record where each number came from.
"""
import dataclasses
import math

import pytest

from TidalPy.Structures.configs import build_world
from TidalPy.Structures.worlds.gasgiant import GasGiantWorld
from TidalPy.Structures.worlds.base import BaseWorld


_G = 6.674e-11
_MASS_JUPITER = 1.898e27


@dataclasses.dataclass
class Body:
    """Published bulk properties and expected interior layout for one bundled world."""
    name: str
    radius: float
    mass: float
    gravity: float
    layers: list
    tidal_layers: list
    spin_period_days: float
    moi_factor: float
    moi_tolerance: float
    # False where no C/MR2 is measured: the expected value is then the file's comment, which pins the model only.
    moi_measured: bool = True
    liquid_layers: list = dataclasses.field(default_factory=list)
    # True for the *_dynamic worlds, whose liquid layers are dynamic and compressible; the rest are static.
    dynamic_liquids: bool = False
    semi_major_axis: float = None
    eccentricity: float = None
    love_k: float = None
    love_period_days: float = None
    love_tolerance: float = None


_BODIES = [
    Body("io", 1821490.0, 8.9298e22, 1.796,
         ["core", "mantle", "asthenosphere"], ["core", "mantle", "asthenosphere"],
         1.769, moi_factor=0.37685, moi_tolerance=2.0e-2,
         semi_major_axis=4.217e8, eccentricity=0.0041),
    Body("europa", 1561000.0, 4.7998e22, 1.315,
         ["core", "mantle", "ice_shell"], ["core", "mantle", "ice_shell"],
         3.551181, moi_factor=0.346, moi_tolerance=1.0e-2),
    Body("luna", 1737400.0, 7.34579e22, 1.6242,
         ["inner_core", "outer_core", "lower_mantle", "mantle", "crust"],
         ["inner_core", "lower_mantle", "mantle", "crust"],
         27.322, moi_factor=0.3931, moi_tolerance=1.0e-3, liquid_layers=["outer_core"],
         love_k=0.02422, love_period_days=27.3217, love_tolerance=1.0e-2),
    Body("mercury", 2440000.0, 3.30103e23, 3.7007,
         ["inner_core", "outer_core", "mantle"], ["inner_core", "mantle"],
         58.6462, moi_factor=0.346, moi_tolerance=1.0e-2, liquid_layers=["outer_core"],
         love_k=0.569, love_period_days=87.9691, love_tolerance=1.0e-2),
    # Love number is the solid-Earth k2 at the M2 period (12.4206 h); the spin is the file's sidereal rate.
    Body("earth_simple", 6371000.0, 5.972e24, 9.820,
         ["inner_core", "outer_core", "mantle"], ["inner_core", "mantle"],
         2.0 * math.pi / 7.292e-5 / 86400.0, moi_factor=0.3307, moi_tolerance=1.0e-3, liquid_layers=["outer_core"],
         love_k=0.30, love_period_days=12.4206 / 24.0, love_tolerance=1.0e-2),
    # earth_simple's structure with a pressure-dependent mantle viscosity, solved with its temperature profile.
    Body("earth_thermal", 6371000.0, 5.972e24, 9.820,
         ["inner_core", "outer_core", "mantle"], ["inner_core", "mantle"],
         2.0 * math.pi / 7.292e-5 / 86400.0, moi_factor=0.3307, moi_tolerance=1.0e-3, liquid_layers=["outer_core"],
         love_k=0.30, love_period_days=12.4206 / 24.0, love_tolerance=1.0e-2),
    # Pluto and Charon are mutually synchronous; Triton is synchronous and retrograde about Neptune.
    Body("pluto", 1188300.0, 1.3024587e22, 0.615627,
         ["core", "ocean", "ice_shell"], ["core", "ice_shell"],
         6.3872304, moi_factor=0.318550, moi_tolerance=1.0e-3, moi_measured=False,
         liquid_layers=["ocean"]),
    # The *_dynamic worlds: their liquid layers are dynamic, compressible, and carry a pressure-dependent EOS.
    Body("luna_dynamic", 1737400.0, 7.34579e22, 1.6242,
         ["inner_core", "outer_core", "lower_mantle", "mantle", "crust"],
         ["inner_core", "lower_mantle", "mantle", "crust"],
         27.322, moi_factor=0.3931, moi_tolerance=1.0e-3, liquid_layers=["outer_core"], dynamic_liquids=True,
         love_k=0.02422, love_period_days=27.3217, love_tolerance=1.0e-2),
    Body("mercury_dynamic", 2440000.0, 3.30103e23, 3.7006,
         ["inner_core", "outer_core", "mantle"], ["inner_core", "mantle"],
         58.6462, moi_factor=0.346, moi_tolerance=1.0e-2, liquid_layers=["outer_core"], dynamic_liquids=True,
         love_k=0.569, love_period_days=87.9691, love_tolerance=1.0e-2),
    Body("pluto_dynamic", 1188300.0, 1.3024587e22, 0.615627,
         ["core", "ocean", "ice_shell"], ["core", "ice_shell"],
         6.3872304, moi_factor=0.31979, moi_tolerance=1.0e-3, moi_measured=False,
         liquid_layers=["ocean"], dynamic_liquids=True),
    Body("europa_dynamic", 1561000.0, 4.7998e22, 1.3147,
         ["core", "mantle", "ocean", "ice_shell"], ["core", "mantle", "ice_shell"],
         3.551181, moi_factor=0.346, moi_tolerance=1.0e-2, liquid_layers=["ocean"], dynamic_liquids=True),
    Body("charon", 606000.0, 1.5896798e21, 0.288915,
         ["core", "ice_shell"], ["core", "ice_shell"],
         6.3872304, moi_factor=0.311564, moi_tolerance=1.0e-3, moi_measured=False),
    Body("triton", 1352600.0, 2.1402926e22, 0.780801,
         ["core", "ice_shell"], ["core", "ice_shell"],
         5.876854, moi_factor=0.315542, moi_tolerance=1.0e-3, moi_measured=False),
    # The TRAPPIST-1 planets are assumed tidally locked, so the spin period is the orbital period.
    Body("trappist1b", 7110044.9, 8.2058009e24, 10.833830,
         ["core", "mantle"], ["core", "mantle"],
         1.510826, moi_factor=0.338021, moi_tolerance=1.0e-3, moi_measured=False),
    Body("trappist1c", 6988995.8, 7.8116358e24, 10.673778,
         ["core", "mantle"], ["core", "mantle"],
         2.421937, moi_factor=0.336954, moi_tolerance=1.0e-3, moi_measured=False),
    Body("trappist1d", 5020354.3, 2.3172131e24, 6.136249,
         ["core", "mantle"], ["core", "mantle"],
         4.049219, moi_factor=0.361525, moi_tolerance=1.0e-3, moi_measured=False),
    Body("trappist1e", 5861327.4, 4.1327614e24, 8.028864,
         ["core", "mantle"], ["core", "mantle"],
         6.101013, moi_factor=0.344955, moi_tolerance=1.0e-3, moi_measured=False),
    Body("trappist1f", 6657703.4, 6.2051143e24, 9.343436,
         ["core", "mantle"], ["core", "mantle"],
         9.207540, moi_factor=0.346947, moi_tolerance=1.0e-3, moi_measured=False),
    Body("trappist1g", 7192868.0, 7.8892744e24, 10.177441,
         ["core", "mantle"], ["core", "mantle"],
         12.352446, moi_factor=0.350022, moi_tolerance=1.0e-3, moi_measured=False),
    Body("trappist1h", 4810111.0, 1.9469367e24, 5.616262,
         ["core", "mantle"], ["core", "mantle"],
         18.772866, moi_factor=0.371397, moi_tolerance=1.0e-3, moi_measured=False),
]

# k2 with the fluid outer core made solid, as each file's comment quotes it.
_SOLID_CORE_LOVE_K = {"luna": 0.02306, "mercury": 0.1158, "earth_simple": 0.2395}

_BODY_CASES = [pytest.param(body, id=body.name) for body in _BODIES]
_LOVE_CASES = [pytest.param(body, id=body.name) for body in _BODIES if body.love_k is not None]
_SOLID_CORE_CASES = [pytest.param(body, id=body.name) for body in _BODIES if body.name in _SOLID_CORE_LOVE_K]


# The solid laws an outer core takes when the counterfactual below freezes it to iron, with the shear modulus and
# viscosity law the quoted counterfactual k2 was solved with.
_IRON_SOLID_LAWS = {
    "shear_modulus":   {"model": "constant", "shear_modulus_pa": 5.25e10},
    "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e20},
}


def _freeze(layer, solid_laws):
    """Make a liquid layer solid: its liquid's equation of state with the given shear laws and a Maxwell rheology."""
    liquid_eos = layer.material.get_config_dict()["liquid"]["eos"]
    layer.material = {"solid": {"eos": liquid_eos, **solid_laws}}
    layer.shear_rheology = "maxwell"
    assert not layer.is_liquid


def _synchronous_tides(world, period_days, host_mass, eccentricity):
    """Each layer's share of the heating on a synchronous orbit of the given period about the given host [kg]."""
    world.solve_eos()
    mean_motion = 2.0 * math.pi / (period_days * 86400.0)
    semi_major_axis = (_G * host_mass / mean_motion ** 2) ** (1.0 / 3.0)
    world.set_spin_frequency(mean_motion)
    world.calc_tides(
        mean_motion,
        mean_motion,
        eccentricity,
        0.0,
        semi_major_axis,
        host_mass,
    )
    total = world.get_tidal_heating()
    return {layer.name: world.get_layer_tidal_heating(index) / total for index, layer in enumerate(world)}


# =====================================================================================================================
# Every bundled body
# =====================================================================================================================
@pytest.mark.parametrize("body", _BODY_CASES)
def test_bundled_body_builds_with_its_stated_bulk_properties(body):
    world = build_world(body.name)
    assert isinstance(world, BaseWorld)
    assert world.radius == pytest.approx(body.radius, rel=1e-12)
    assert world.mass == pytest.approx(body.mass, rel=1e-12)
    assert [layer.name for layer in world] == body.layers


@pytest.mark.parametrize("body", _BODY_CASES)
def test_bundled_body_interior_reproduces_its_mass(body):
    """The solved EOS mass matches the stated mass the densities were fitted to; gravity checks the radius."""
    world = build_world(body.name)
    result = world.solve_eos()
    assert result["success"], result["message"]
    assert world.eos_solved
    assert world.planet_mass_eos == pytest.approx(body.mass, rel=1e-6)
    assert world.surface_gravity_eos == pytest.approx(body.gravity, rel=1e-3)


@pytest.mark.parametrize("body", _BODY_CASES)
def test_bundled_body_moment_of_inertia_factor_matches_its_expected_value(body):
    world = build_world(body.name)
    world.solve_eos()
    moi_factor = world.planet_moi_eos / (world.planet_mass_eos * body.radius ** 2)
    assert moi_factor == pytest.approx(body.moi_factor, rel=body.moi_tolerance)
    # A uniform sphere gives exactly 0.4; every one of these bodies is centrally condensed.
    assert moi_factor < 0.4


@pytest.mark.parametrize("body", _BODY_CASES)
def test_bundled_body_density_profile_is_physical(body):
    """Density is positive and falls outward; pressure falls outward and vanishes at the surface."""
    world = build_world(body.name)
    world.solve_eos()
    radii = [body.radius * fraction for fraction in (0.02, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)]
    densities = [world.get_density(r) for r in radii]
    assert all(d > 0.0 for d in densities), densities
    assert densities[0] >= densities[-1], densities
    pressures = [world.get_pressure(r) for r in radii]
    assert all(earlier >= later for earlier, later in zip(pressures, pressures[1:])), pressures
    assert world.get_pressure(body.radius) == pytest.approx(0.0, abs=1.0e5)


@pytest.mark.parametrize("body", _BODY_CASES)
def test_bundled_body_dissipates_in_the_intended_layers(body):
    world = build_world(body.name)
    assert [layer.name for layer in world if layer.use_tides] == body.tidal_layers


@pytest.mark.parametrize("body", _BODY_CASES)
def test_bundled_body_declares_its_liquid_layers(body):
    """Every liquid layer is a static liquid, except in the *_dynamic worlds, where it is dynamic; none is
    incompressible, and every solid layer is static."""
    world = build_world(body.name)
    assert [layer.name for layer in world if layer.is_liquid] == body.liquid_layers
    for layer in world:
        assert layer.is_static == (not layer.is_liquid or not body.dynamic_liquids), layer.name
    assert not any(layer.is_incompressible for layer in world)


@pytest.mark.parametrize("body", _BODY_CASES)
def test_bundled_body_layer_heating_partitions_the_total(body):
    """Layer shares of the heating sum to one and a liquid layer takes none."""
    world = build_world(body.name)
    # Only the quasi-homogeneous Love methods read tidal scales, so the bundled worlds leave them unset.
    assert all(layer.tidal_scale is None for layer in world)
    shares = _synchronous_tides(world, body.spin_period_days, _MASS_JUPITER, 0.01)
    assert sum(shares.values()) == pytest.approx(1.0, rel=1e-10), shares
    for name in body.liquid_layers:
        assert shares[name] == 0.0, name
    assert all(share >= 0.0 for share in shares.values()), shares


@pytest.mark.parametrize("body", _BODY_CASES)
def test_bundled_body_states_its_temperatures_and_every_solid_layer_dissipates(body):
    """Every layer has a temperature, and every solid layer is tidal with a finite viscosity."""
    world = build_world(body.name)
    world.solve_eos()
    for layer in world:
        assert layer.temperature > 0.0, layer.name
        assert layer.use_tides == (not layer.is_liquid), layer.name
        if not layer.is_liquid:
            radius = 0.5 * (layer.radius_inner + layer.radius_outer)
            assert 0.0 < world.get_shear_viscosity(radius) < math.inf, layer.name


@pytest.mark.parametrize("body", _BODY_CASES)
def test_bundled_body_spins_at_its_stated_rotation_period(body):
    world = build_world(body.name)
    assert world.spin_frequency == pytest.approx(2.0 * math.pi / (body.spin_period_days * 86400.0), rel=1e-6)


@pytest.mark.parametrize("body", _LOVE_CASES)
def test_bundled_body_love_number_matches_its_measured_value(body):
    """The mantle rigidity of each liquid-core world is fitted to its measured k2."""
    world = build_world(body.name)
    world.solve_eos()
    world.solve_love_numbers(2.0 * math.pi / (body.love_period_days * 86400.0))
    assert world.love_success, world.love_message
    assert world.love_number_k.real == pytest.approx(body.love_k, rel=body.love_tolerance)


@pytest.mark.parametrize("body", _SOLID_CORE_CASES)
def test_bundled_body_needs_its_fluid_outer_core_for_its_measured_love_number(body):
    """With the outer core made solid, k2 drops to the quoted counterfactual."""
    frequency = 2.0 * math.pi / (body.love_period_days * 86400.0)
    world = build_world(body.name)
    world.solve_eos()
    world.solve_love_numbers(frequency)
    fluid_core_k = world.love_number_k.real

    _freeze(world.outer_core, _IRON_SOLID_LAWS)
    world.solve_eos()
    world.solve_love_numbers(frequency)
    solid_core_k = world.love_number_k.real
    assert solid_core_k == pytest.approx(_SOLID_CORE_LOVE_K[body.name], rel=1e-2)
    assert solid_core_k < fluid_core_k


def test_mercury_fluid_core_is_worth_a_factor_of_five_in_its_love_number():
    body = next(body for body in _BODIES if body.name == "mercury")
    assert _SOLID_CORE_LOVE_K["mercury"] * 4.0 < body.love_k


# =====================================================================================================================
# Io's heat budget
# =====================================================================================================================
# Lainey et al. (2009), from the astrometric orbital acceleration: (9.33 +/- 1.83)e13 W.
_IO_HEATING = 9.33e13
_IO_HEATING_UNCERTAINTY = 1.83e13


@pytest.fixture(scope="module")
def io_tides():
    """Io on its orbit about Jupiter after calc_tides; the tests only read from it."""
    body = _BODIES[0]
    world = build_world(body.name)
    world.solve_eos()
    orbital_frequency = math.sqrt(_G * _MASS_JUPITER / body.semi_major_axis ** 3)
    world.set_spin_frequency(orbital_frequency)
    world.calc_tides(
        orbital_frequency,
        orbital_frequency,
        body.eccentricity,
        0.0,
        body.semi_major_axis,
        _MASS_JUPITER,
    )
    return world, body


def test_io_total_heating_matches_the_astrometric_value(io_tides):
    """The asthenosphere viscosity was fitted to this."""
    world, _ = io_tides
    assert world.get_tidal_heating() == pytest.approx(_IO_HEATING, rel=1e-3)


def test_io_heating_is_within_the_observational_uncertainty(io_tides):
    world, _ = io_tides
    assert abs(world.get_tidal_heating() - _IO_HEATING) < _IO_HEATING_UNCERTAINTY


def test_io_surface_heat_flux_is_about_two_watts_per_square_metre(io_tides):
    world, body = io_tides
    flux = world.get_tidal_heating() / (4.0 * math.pi * body.radius ** 2)
    assert flux == pytest.approx(2.24, abs=0.05)


def test_io_per_layer_heating_sums_to_the_total(io_tides):
    """With several dissipating layers the scales must partition the heat rather than repeat it."""
    world, _ = io_tides
    per_layer = sum(world.get_layer_tidal_heating(index) for index in range(world.num_layers))
    assert per_layer == pytest.approx(world.get_tidal_heating(), rel=1e-9)


def test_io_asthenosphere_dominates_the_dissipation(io_tides):
    world, _ = io_tides
    heating = {layer.name: world.get_layer_tidal_heating(index) for index, layer in enumerate(world)}
    # The solid iron core dissipates, as every solid layer does, but only a few parts per million.
    assert 0.0 < heating["core"] < 1.0e-5 * world.get_tidal_heating()
    assert heating["asthenosphere"] > 20.0 * heating["mantle"]
    assert heating["asthenosphere"] / world.get_tidal_heating() > 0.9


# =====================================================================================================================
# Luna's dissipation at two periods
# =====================================================================================================================
# Williams and Boggs (2015): Q = 38 +/- 4 at the month and 41 +/- 9 at the year, with a peak near 100 to 120 days.
_LUNA_Q_MONTH, _LUNA_Q_MONTH_UNCERTAINTY = 38.0, 4.0
_LUNA_Q_YEAR, _LUNA_Q_YEAR_UNCERTAINTY = 41.0, 9.0


def _luna_quality(world, period_days):
    world.solve_love_numbers(2.0 * math.pi / (period_days * 86400.0))
    assert world.love_success, world.love_message
    return world.love_number_k.real / (-world.love_number_k.imag)


def test_luna_quality_factor_at_the_month_and_the_year():
    """The low-viscosity zone's viscosity and top are fitted to these two."""
    world = build_world("luna")
    world.solve_eos()
    q_month = _luna_quality(world, 27.3217)
    q_year = _luna_quality(world, 365.25)
    assert q_month == pytest.approx(_LUNA_Q_MONTH, abs=0.5)
    assert q_year == pytest.approx(_LUNA_Q_YEAR, abs=0.5)
    assert abs(q_month - _LUNA_Q_MONTH) < _LUNA_Q_MONTH_UNCERTAINTY
    assert abs(q_year - _LUNA_Q_YEAR) < _LUNA_Q_YEAR_UNCERTAINTY


def test_luna_relaxation_peak_sits_between_the_month_and_the_year():
    """The whole-body Q reaches its minimum near 100 days."""
    world = build_world("luna")
    world.solve_eos()
    q_by_period = {period: _luna_quality(world, period) for period in (10.0, 27.3217, 96.0, 110.0, 365.25, 1000.0)}
    assert q_by_period[96.0] < q_by_period[27.3217] and q_by_period[96.0] < q_by_period[365.25]
    assert q_by_period[110.0] == pytest.approx(q_by_period[96.0], rel=0.05)
    assert q_by_period[10.0] > q_by_period[27.3217] and q_by_period[1000.0] > q_by_period[365.25]


def test_luna_deep_zone_does_most_of_the_dissipating():
    """About 85 percent of the heat is in the deep zone and 15 in the mantle above it."""
    shares = _synchronous_tides(build_world("luna"), 27.321661, 5.972e24, 0.0549)
    assert shares["lower_mantle"] == pytest.approx(0.852, abs=0.005)
    assert shares["mantle"] == pytest.approx(0.148, abs=0.005)
    assert shares["crust"] < 1.0e-3 and shares["inner_core"] < 1.0e-6


# =====================================================================================================================
# Earth's dissipation with its temperature profile solved
# =====================================================================================================================
# The solid Earth's M2 quality factor is about 280 (Ray et al. 2001); earth_thermal's header quotes 315.
_M2_PERIOD_DAYS = 12.4206 / 24.0
_EARTH_SURFACE_TEMPERATURE = 288.0


def _m2_quality(world_name):
    world = build_world(world_name)
    world.solve_eos(solve_temperature=True, surface_temperature=_EARTH_SURFACE_TEMPERATURE)
    world.solve_love_numbers(2.0 * math.pi / (_M2_PERIOD_DAYS * 86400.0))
    assert world.love_success, world.love_message
    return world, world.love_number_k.real / (-world.love_number_k.imag)


def test_earth_thermal_quality_factor_at_m2_is_the_solid_earths():
    """With the profile solved, the pressure-dependent viscosity keeps the deep mantle stiff and Q near Earth's;
    earth_simple's pressure-independent viscosity, fitted to an isothermal mantle, softens there and dissipates ten
    times more."""
    _, quality_thermal = _m2_quality("earth_thermal")
    _, quality_simple = _m2_quality("earth_simple")
    assert quality_thermal == pytest.approx(315.0, rel=0.05)
    assert quality_simple < 0.2 * quality_thermal


def test_earth_thermal_mantle_viscosity_rises_with_depth():
    """Below the conducting lid the adiabat warms the mantle, yet its viscosity rises with depth."""
    world, _ = _m2_quality("earth_thermal")
    depths = (400.0e3, 1000.0e3, 2000.0e3)
    viscosities = [world.get_shear_viscosity(world.radius - depth) for depth in depths]
    temperatures = [world.get_temperature(world.radius - depth) for depth in depths]
    assert temperatures[0] < temperatures[1] < temperatures[2]
    assert viscosities[0] < viscosities[1] < viscosities[2]


# =====================================================================================================================
# The TRAPPIST-1 system
# =====================================================================================================================
# Agol et al. (2021) stellar mass, luminosity, and orbits; the radius is Delrez et al. (2018) from constants_.hpp.
_TRAPPIST1_MASS = 1.785615e29
_TRAPPIST1_RADIUS = 82927440.0
_TRAPPIST1_LUMINOSITY = 2.095447e23
_TRAPPIST1_EFFECTIVE_TEMPERATURE = 2566.0

# Planet: (semi-major axis [AU], orbital period [days], core mass fraction the file's comment quotes).
_TRAPPIST1_PLANETS = {
    "trappist1b": (0.01154, 1.510826, 0.289980),
    "trappist1c": (0.01580, 2.421937, 0.301773),
    "trappist1d": (0.02227, 4.049219, 0.125329),
    "trappist1e": (0.02925, 6.101013, 0.226176),
    "trappist1f": (0.03849, 9.207540, 0.204083),
    "trappist1g": (0.04683, 12.352446, 0.173900),
    "trappist1h": (0.06189, 18.772866, 0.075819),
}

_AU = 1.495978707e11

# Earth's core mass fraction through the same two-layer recipe; larger than the accepted 0.325 because a uniform
# incompressible core is crude, which is why the planets are compared against this and not the accepted value.
_EARTH_RECIPE_CORE_MASS_FRACTION = 0.354


def test_trappist1_star_states_its_published_properties():
    star = build_world("trappist1")
    assert star.radius == pytest.approx(_TRAPPIST1_RADIUS, rel=1e-12)
    assert star.mass == pytest.approx(_TRAPPIST1_MASS, rel=1e-12)
    assert star.effective_temperature == pytest.approx(_TRAPPIST1_EFFECTIVE_TEMPERATURE, rel=1e-12)


def test_trappist1_luminosity_is_the_measured_one_not_the_derived_one():
    """Stefan-Boltzmann on the published radius and temperature misses the published luminosity by 1.4 percent."""
    star = build_world("trappist1")
    assert star.luminosity == pytest.approx(_TRAPPIST1_LUMINOSITY, rel=1e-9)

    stefan_boltzmann = 5.670374419e-8
    derived = 4.0 * math.pi * _TRAPPIST1_RADIUS ** 2 * stefan_boltzmann * _TRAPPIST1_EFFECTIVE_TEMPERATURE ** 4
    assert derived == pytest.approx(_TRAPPIST1_LUMINOSITY, rel=2.0e-2)
    assert derived != pytest.approx(_TRAPPIST1_LUMINOSITY, rel=1.0e-2)


@pytest.mark.parametrize("name", sorted(_TRAPPIST1_PLANETS))
def test_trappist1_planet_period_follows_from_the_stellar_mass_and_its_semi_major_axis(name):
    """Kepler's third law ties the star file to each planet's orbital period, which is also its spin."""
    semi_major_axis_au, period_days, _ = _TRAPPIST1_PLANETS[name]
    mean_motion = math.sqrt(_G * _TRAPPIST1_MASS / (semi_major_axis_au * _AU) ** 3)
    # The published semi-major axes carry four significant digits.
    assert mean_motion == pytest.approx(2.0 * math.pi / (period_days * 86400.0), rel=1.0e-3)
    assert build_world(name).spin_frequency == pytest.approx(2.0 * math.pi / (period_days * 86400.0), rel=1e-6)


@pytest.mark.parametrize("name", sorted(_TRAPPIST1_PLANETS))
def test_trappist1_planet_is_iron_depleted_relative_to_earth(name):
    """The core mass fraction follows from the fitted core radius and lands below Earth's by the same recipe."""
    _, _, expected_core_mass_fraction = _TRAPPIST1_PLANETS[name]
    world = build_world(name)
    world.solve_eos()
    core_mass_fraction = world.core.mass / world.planet_mass_eos
    assert core_mass_fraction == pytest.approx(expected_core_mass_fraction, rel=1e-3)
    assert core_mass_fraction < _EARTH_RECIPE_CORE_MASS_FRACTION


# =====================================================================================================================
# The gas giants
# =====================================================================================================================
# Two outer densities of each are fitted to the mass and C/MR2 at once, so C/MR2 is reproduced by construction.
_GAS_GIANTS = {
    # name: (radius, mass, gravity, C/MR2, spin period [hours], k2, Q)
    "jupiter": (69911000.0, 1.898125e27, 25.920269, 0.2640, 9.925, 0.590, 53500.0),
    "neptune": (24622000.0, 1.02409e26, 11.274498, 0.2400, 16.11, 0.410, 9000.0),
}

# Layer masses in Earth masses, as each file's comment quotes them.
_GAS_GIANT_LAYER_MASSES = {
    "jupiter": {"core": 16.177, "mantle": 268.439, "envelope": 33.211},
    "neptune": {"core": 0.503, "mantle": 14.119, "envelope": 2.526},
}

# k2 from solving the fitted interior; nothing is fitted to it.
_GAS_GIANT_SOLVED_LOVE_K = {"jupiter": 0.53438, "neptune": 0.42667}

_MASS_EARTH = 5.9721986e24

# Satellites used to raise a tide on each gas giant: Io on Jupiter, Triton on Neptune.
_TRITON_SEMI_MAJOR_AXIS = 3.547590e8
_TRITON_MASS = 2.1402926e22
_IO_SEMI_MAJOR_AXIS = 4.217e8
_IO_MASS = 8.9298e22


@pytest.mark.parametrize("name", sorted(_GAS_GIANTS))
def test_gas_giant_builds_with_its_stated_bulk_properties(name):
    radius, mass, _, _, spin_hours, _, _ = _GAS_GIANTS[name]
    world = build_world(name)
    assert isinstance(world, GasGiantWorld)
    assert world.radius == pytest.approx(radius, rel=1e-12)
    assert world.mass == pytest.approx(mass, rel=1e-12)
    assert [layer.name for layer in world] == ["core", "mantle", "envelope"]
    assert world.spin_frequency == pytest.approx(2.0 * math.pi / (spin_hours * 3600.0), rel=1e-6)


@pytest.mark.parametrize("name", sorted(_GAS_GIANTS))
def test_gas_giant_interior_reproduces_its_mass_and_moment_of_inertia(name):
    radius, mass, gravity, moi_factor, _, _, _ = _GAS_GIANTS[name]
    world = build_world(name)
    result = world.solve_eos()
    assert result["success"], result["message"]
    assert world.planet_mass_eos == pytest.approx(mass, rel=1e-6)
    assert world.surface_gravity_eos == pytest.approx(gravity, rel=1e-3)
    assert world.planet_moi_eos / (world.planet_mass_eos * radius ** 2) == pytest.approx(moi_factor, rel=1e-4)


@pytest.mark.parametrize("name", sorted(_GAS_GIANTS))
def test_gas_giant_layers_are_all_fluid(name):
    """A rigid rock core would also fail numerically: its shear modulus is negligible against rho g R."""
    world = build_world(name)
    assert [layer.name for layer in world if not layer.is_liquid] == []
    assert all(layer.is_static for layer in world)


@pytest.mark.parametrize("name", sorted(_GAS_GIANTS))
def test_gas_giant_solved_love_number_is_near_its_published_value(name):
    """The fitted interior lands within 15 percent of the published k2, where one uniform layer gives 1.5."""
    _, _, _, _, _, published_k, _ = _GAS_GIANTS[name]
    world = build_world(name)
    world.solve_eos()
    result = world.solve_love_numbers()
    assert result["success"], result["message"]
    assert world.love_number_k.real == pytest.approx(_GAS_GIANT_SOLVED_LOVE_K[name], rel=1e-3)
    assert abs(world.love_number_k.real - published_k) / published_k < 0.15


@pytest.mark.parametrize("name", sorted(_GAS_GIANTS))
def test_gas_giant_layer_masses_match_its_file(name):
    """The layer masses are a result of the fit and tie it to published interior models."""
    world = build_world(name)
    world.solve_eos()
    for layer in world:
        expected = _GAS_GIANT_LAYER_MASSES[name][layer.name]
        assert layer.mass / _MASS_EARTH == pytest.approx(expected, rel=1e-3), layer.name
    # The envelope is the only tidal layer, and it carries the whole scale.
    assert [layer.name for layer in world if layer.use_tides] == ["envelope"]
    assert sum(layer.tidal_scale for layer in world if layer.use_tides) == pytest.approx(1.0, abs=1e-4)


@pytest.mark.parametrize("name, semi_major_axis, satellite_mass, spin_multiple", [
    ("neptune", _TRITON_SEMI_MAJOR_AXIS, _TRITON_MASS, None),
    ("neptune", _TRITON_SEMI_MAJOR_AXIS, _TRITON_MASS, 3.0),
    ("neptune", 6.0e8, _TRITON_MASS, None),
    ("neptune", _TRITON_SEMI_MAJOR_AXIS, 1.0e23, None),
    ("jupiter", _IO_SEMI_MAJOR_AXIS, _IO_MASS, None),
    ("jupiter", _IO_SEMI_MAJOR_AXIS, _IO_MASS, 2.0),
])
def test_gas_giant_dissipation_matches_the_constant_phase_lag_closed_form(
        name, semi_major_axis, satellite_mass, spin_multiple):
    """With one active mode the heating is (3/4)(k2/Q)(G M_s^2 R^5 / a^6)|2(spin - n)| exactly."""
    radius, mass, _, _, _, love_k, quality = _GAS_GIANTS[name]
    world = build_world(name)
    mean_motion = math.sqrt(_G * mass / semi_major_axis ** 3)
    spin = world.spin_frequency if spin_multiple is None else spin_multiple * mean_motion

    world.calc_tides(
        mean_motion,
        spin,
        0.0,
        0.0,
        semi_major_axis,
        satellite_mass,
    )
    closed_form = (0.75 * (love_k / quality)
                   * (_G * satellite_mass ** 2 * radius ** 5 / semi_major_axis ** 6)
                   * abs(2.0 * (spin - mean_motion)))
    assert world.get_tidal_heating() == pytest.approx(closed_form, rel=1e-3)


def test_jupiter_simple_reproduces_its_mass_but_not_its_moment_of_inertia_or_love_number():
    """One uniform layer matches the mass only: C/MR2 is 0.4 and k2 is the fluid-sphere 1.5."""
    world = build_world("jupiter_simple")
    result = world.solve_eos()
    assert result["success"], result["message"]
    assert world.planet_mass_eos == pytest.approx(world.mass, rel=1e-6)
    assert world.planet_moi_eos / (world.planet_mass_eos * world.radius ** 2) == pytest.approx(0.4, rel=1e-6)
    assert world.solve_love_numbers()["success"]
    assert world.love_number_k.real == pytest.approx(1.5, rel=1e-3)


# =====================================================================================================================
# Pluto's ocean
# =====================================================================================================================
# The 165 km shell is the Nimmo et al. (2016) depth to the ocean; the ocean is what is left above the 0.7 core.
_PLUTO_SHELL_THICKNESS = 165.0e3
_PLUTO_OCEAN_THICKNESS = 112.4e3
_PLUTO_OCEAN_K = 0.145781
_PLUTO_FROZEN_K = 0.004574


def test_pluto_ocean_thickness_follows_from_the_observed_shell():
    """The ocean thickness is not fitted, and lands inside the 50 to 150 km of thermal evolution models."""
    world = build_world("pluto")
    shell = world.radius - world.ocean.radius_outer
    ocean = world.ocean.radius_outer - world.core.radius_outer
    assert shell == pytest.approx(_PLUTO_SHELL_THICKNESS, rel=1e-6)
    assert ocean == pytest.approx(_PLUTO_OCEAN_THICKNESS, rel=1e-3)
    assert 50.0e3 < ocean < 150.0e3


def test_pluto_ocean_is_worth_a_factor_of_thirty_in_its_love_number():
    """The static liquid ocean decouples the shell from the core."""
    world = build_world("pluto")
    frequency = world.spin_frequency
    world.solve_eos()
    world.solve_love_numbers(frequency)
    assert world.love_success, world.love_message
    ocean_k = world.love_number_k.real
    assert ocean_k == pytest.approx(_PLUTO_OCEAN_K, rel=1e-3)

    # Frozen into the shell, as demo P09 freezes it: the shell's ice, rheology, and temperature, viscoelastic like the
    # shell rather than held at its static moduli.
    world.ocean.material = world.ice_shell.material
    world.ocean.shear_rheology = world.ice_shell.shear_rheology
    world.ocean.temperature = world.ice_shell.temperature
    world.ocean.use_tides = True
    world.solve_eos()
    assert not world.ocean.is_liquid
    world.solve_love_numbers(frequency)
    frozen_k = world.love_number_k.real
    assert frozen_k == pytest.approx(_PLUTO_FROZEN_K, rel=1e-3)
    assert ocean_k > 25.0 * frozen_k


def test_pluto_ocean_takes_none_of_the_tidal_heating():
    world = build_world("pluto")
    assert not world.ocean.use_tides
    # Charon's mass on the mutual orbit.
    shares = _synchronous_tides(world, 6.3872, 1.586e21, 0.005)
    assert shares["ocean"] == 0.0
    # The decoupled shell takes nearly all of it; without the ocean the core would take 4 percent.
    assert shares["ice_shell"] > 0.99


# =====================================================================================================================
# The *_dynamic worlds
# =====================================================================================================================
_DYNAMIC_CASES = [pytest.param(body, id=body.name) for body in _BODIES if body.dynamic_liquids]


@pytest.mark.parametrize("body", _DYNAMIC_CASES)
def test_dynamic_world_liquid_follows_its_equation_of_state(body):
    """The liquid is denser with depth, near d rho / dr = -rho^2 g / K, which keeps its dynamic solve neutral."""
    world = build_world(body.name)
    world.solve_eos()
    for name in body.liquid_layers:
        layer = getattr(world, name)
        radius = 0.5 * (layer.radius_inner + layer.radius_outer)
        step = 1.0e-3 * (layer.radius_outer - layer.radius_inner)
        density = world.get_density(radius)
        slope = (world.get_density(radius + step) - world.get_density(radius - step)) / (2.0 * step)
        neutral_slope = -density ** 2 * world.get_gravity(radius) / world.get_bulk_modulus(radius)
        assert slope < 0.0, name
        assert slope == pytest.approx(neutral_slope, rel=1.0e-2), name


# k2 at the orbital or mutual period as each file's comment quotes it, and the most the same liquid solved as a
# static one moves it.
_DYNAMIC_LOVE_K = {"pluto_dynamic": (6.3872304, 0.1556, 2.0e-4), "europa_dynamic": (3.551181, 0.2594, 2.0e-3)}


@pytest.mark.parametrize("name", list(_DYNAMIC_LOVE_K))
def test_dynamic_ocean_world_love_number_and_its_static_counterpart(name):
    period_days, love_k, static_difference = _DYNAMIC_LOVE_K[name]
    frequency = 2.0 * math.pi / (period_days * 86400.0)
    world = build_world(name)
    world.solve_eos()
    world.solve_love_numbers(frequency)
    dynamic_k = world.love_number_k
    assert dynamic_k.real == pytest.approx(love_k, rel=1.0e-3)
    world.ocean.is_static = True
    world.solve_love_numbers(frequency)
    assert abs(world.love_number_k - dynamic_k) < static_difference * abs(dynamic_k)


def test_europa_dynamic_ocean_decouples_the_shell():
    """With its ocean Europa's k2 is an ocean world's, about 17 times the ocean-free europa.toml's 0.0153."""
    frequency = 2.0 * math.pi / (3.551181 * 86400.0)
    ocean_world = build_world("europa_dynamic")
    ocean_world.solve_eos()
    ocean_world.solve_love_numbers(frequency)
    solid_world = build_world("europa")
    solid_world.solve_eos()
    solid_world.solve_love_numbers(frequency)
    assert solid_world.love_number_k.real == pytest.approx(0.0153, rel=1.0e-2)
    assert ocean_world.love_number_k.real > 15.0 * solid_world.love_number_k.real


def test_luna_dynamic_keeps_the_quality_factor_fit_at_the_default_tolerances():
    """The default radial tolerances keep luna_dynamic's yearly Q on luna.toml's fit (its dynamic core once needed a
    pinned rtol of 1e-8 for that); the month is unaffected."""
    world = build_world("luna_dynamic")
    assert "radial_solver" not in world.source_config
    world.solve_eos()
    assert _luna_quality(world, 27.3217) == pytest.approx(_LUNA_Q_MONTH, abs=0.5)
    assert _luna_quality(world, 365.25) == pytest.approx(_LUNA_Q_YEAR, abs=0.5)
