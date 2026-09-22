"""Bundled solar-system bodies rebuilt from the classic BurnMan WorldPack.

Each of these worlds states a published mass and radius and then describes an interior that has to
reproduce that mass when the equation of state is solved. A layer density copied from the wrong source, a
mistyped boundary radius, or an EOS model that compresses differently than intended all show up as a mass
mismatch, so the mass check is the substance of these tests rather than a formality. Surface gravity is
checked against its published value as an independent handle: it follows from mass and radius, so it
catches a radius that disagrees with the one the interior was fitted to.

The moment-of-inertia factor is the stronger check, because only Luna's crust thickness and Earth-Simple's
densities are fitted to it. C/MR2 measures how centrally concentrated the mass is, so it constrains the core
size rather than the total mass, and each body carries its own tolerance recording how well its model does
against the measured value.

Every layer states its temperature and every solid layer carries a reference-law viscosity at it, so every
solid layer dissipates, and each world's per-layer tidal_scale values come from integrating the depth-resolved
3D heating (author decisions 2026-09-19 and 2026-09-21). Four bodies pin a tidal observable as well. Io has
its asthenosphere viscosity fitted so that the total dissipation matches the astrometric value of Lainey et
al. (2009). Mercury, Luna, and Earth-Simple declare their fluid outer cores as static liquids and have their
mantle rigidity fitted to the measured k2; for each the fluid core is worth a large factor in k2, so the
solid-core counterfactual is pinned too. Luna's low-viscosity zone at the base of the mantle is fitted to
the lunar Q at the month and at the year (Williams and Boggs 2015).

Expected values come from the TOML comments, which record where each number came from.
"""
import math

import pytest

from TidalPy.structures_x.configs import build_world
from TidalPy.structures_x.worlds.layered import LayeredWorld


# Jupiter and the Galilean orbits, for the tidal solves below.
_G = 6.674e-11
_MASS_JUPITER = 1.898e27


class Body:
    """Published bulk properties and expected interior layout for one bundled world."""

    def __init__(self,
                 name,
                 radius,
                 mass,
                 gravity,
                 layers,
                 tidal_layers,
                 spin_period_days,
                 moi_factor,
                 moi_tolerance,
                 liquid_layers=(),
                 semi_major_axis=None,
                 eccentricity=None,
                 love_k=None,
                 love_period_days=None,
                 love_tolerance=None):
        self.name = name
        self.radius = radius
        self.mass = mass
        self.gravity = gravity
        self.layers = layers
        self.tidal_layers = tidal_layers
        self.liquid_layers = list(liquid_layers)
        self.spin_period_days = spin_period_days
        self.moi_factor = moi_factor
        self.moi_tolerance = moi_tolerance
        self.semi_major_axis = semi_major_axis
        self.eccentricity = eccentricity
        self.love_k = love_k
        self.love_period_days = love_period_days
        self.love_tolerance = love_tolerance


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
    # Williams (1994) for C/MR2; the solid-Earth k2 at the M2 period (12.4206 h); the spin is the file's sidereal rate.
    Body("earth_simple", 6371000.0, 5.972e24, 9.820,
         ["inner_core", "outer_core", "mantle"], ["inner_core", "mantle"],
         2.0 * math.pi / 7.292e-5 / 86400.0, moi_factor=0.3307, moi_tolerance=1.0e-3, liquid_layers=["outer_core"],
         love_k=0.30, love_period_days=12.4206 / 24.0, love_tolerance=1.0e-2),
]

# The solid-core counterfactual each liquid-core world's comment quotes: k2 with the outer core solid.
_SOLID_CORE_LOVE_K = {"luna": 0.02306, "mercury": 0.1158, "earth_simple": 0.2395}


def _cases():
    return [pytest.param(body, id=body.name) for body in _BODIES]


def _love_cases():
    return [pytest.param(body, id=body.name) for body in _BODIES if body.love_k is not None]


@pytest.mark.parametrize("body", _cases())
def test_bundled_body_builds_with_its_stated_bulk_properties(body):
    world = build_world(body.name)
    assert isinstance(world, LayeredWorld)
    assert world.radius == pytest.approx(body.radius, rel=1e-12)
    assert world.mass == pytest.approx(body.mass, rel=1e-12)
    assert [layer.name for layer in world] == body.layers


@pytest.mark.parametrize("body", _cases())
def test_bundled_body_interior_reproduces_its_mass(body):
    """The solved EOS mass must match the stated mass: that is what the layer densities were fitted to."""
    world = build_world(body.name)
    result = world.solve_eos()
    assert result["success"], result["message"]
    assert world.eos_solved
    assert world.planet_mass_eos == pytest.approx(body.mass, rel=1e-6)
    assert world.surface_gravity_eos == pytest.approx(body.gravity, rel=1e-3)


@pytest.mark.parametrize("body", _cases())
def test_bundled_body_moment_of_inertia_factor_matches_its_published_value(body):
    """Nothing is fitted to C/MR2, so it is an independent test of where the mass sits."""
    world = build_world(body.name)
    world.solve_eos()
    moi_factor = world.planet_moi_eos / (world.planet_mass_eos * body.radius ** 2)
    assert moi_factor == pytest.approx(body.moi_factor, rel=body.moi_tolerance)
    # A uniform sphere gives exactly 0.4; every one of these bodies is centrally condensed.
    assert moi_factor < 0.4


@pytest.mark.parametrize("body", _cases())
def test_bundled_body_density_profile_is_physical(body):
    """Density must not increase outward across the body, and every slice must be positive."""
    world = build_world(body.name)
    world.solve_eos()
    radii = [body.radius * fraction for fraction in (0.02, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)]
    densities = [world.get_density(r) for r in radii]
    assert all(d > 0.0 for d in densities), densities
    assert densities[0] >= densities[-1], densities
    # Pressure falls monotonically outward and vanishes at the surface.
    pressures = [world.get_pressure(r) for r in radii]
    assert all(earlier >= later for earlier, later in zip(pressures, pressures[1:])), pressures
    assert world.get_pressure(body.radius) == pytest.approx(0.0, abs=1.0e5)


@pytest.mark.parametrize("body", _cases())
def test_bundled_body_dissipates_in_the_intended_layers(body):
    world = build_world(body.name)
    assert [layer.name for layer in world if layer.is_tidal] == body.tidal_layers


@pytest.mark.parametrize("body", _cases())
def test_bundled_body_declares_its_liquid_layers(body):
    """Fluid cores are static liquids; no icy moon carries a liquid layer (author decision 2026-09-19)."""
    world = build_world(body.name)
    assert [layer.name for layer in world if not layer.is_solid] == body.liquid_layers
    assert all(layer.is_static for layer in world)
    assert not any(layer.is_incompressible for layer in world)


@pytest.mark.parametrize("body", _cases())
def test_bundled_body_tidal_scales_sum_to_one(body):
    """Per-layer scales must partition the world's heating, not repeat it."""
    world = build_world(body.name)
    scales = [layer.tidal_scale for layer in world if layer.is_tidal]
    assert sum(scales) == pytest.approx(1.0, abs=1e-4), scales


@pytest.mark.parametrize("body", _cases())
def test_bundled_body_states_its_temperatures_and_every_solid_layer_dissipates(body):
    """Every layer has a temperature, every solid layer a finite viscosity at it, and every solid layer is tidal."""
    world = build_world(body.name)
    world.solve_eos()
    for layer in world:
        assert layer.temperature > 0.0, layer.name
        assert layer.is_tidal == layer.is_solid, layer.name
        if layer.is_solid:
            radius = 0.5 * (layer.radius_inner + layer.radius_outer)
            assert 0.0 < world.get_shear_viscosity(radius) < math.inf, layer.name


@pytest.mark.parametrize("body", _cases())
def test_bundled_body_spins_at_its_classic_period(body):
    """Spin came from the classic file's spin_period, synchronous for every body here except Mercury and Earth."""
    world = build_world(body.name)
    assert world.spin_frequency == pytest.approx(
        2.0 * math.pi / (body.spin_period_days * 86400.0), rel=1e-6)


@pytest.mark.parametrize("body", _love_cases())
def test_bundled_body_love_number_matches_its_measured_value(body):
    """The mantle rigidity of each liquid-core world is fitted to this; the file's comment quotes the result."""
    world = build_world(body.name)
    world.solve_eos()
    world.solve_love_numbers(2.0 * math.pi / (body.love_period_days * 86400.0))
    assert world.love_success, world.love_message
    assert world.love_number_k.real == pytest.approx(body.love_k, rel=body.love_tolerance)


@pytest.mark.parametrize("body", _love_cases())
def test_bundled_body_needs_its_fluid_outer_core_for_its_measured_love_number(body):
    """The measured k2 is reached only with the outer core liquid; solid, it drops to the quoted counterfactual."""
    frequency = 2.0 * math.pi / (body.love_period_days * 86400.0)
    world = build_world(body.name)
    world.solve_eos()
    world.solve_love_numbers(frequency)
    fluid_core_k = world.love_number_k.real

    world.outer_core.is_solid = True
    world.solve_eos()
    world.solve_love_numbers(frequency)
    solid_core_k = world.love_number_k.real
    assert solid_core_k == pytest.approx(_SOLID_CORE_LOVE_K[body.name], rel=1e-2)
    assert solid_core_k < fluid_core_k


# =====================================================================================================================
# Io's heat budget
# =====================================================================================================================
# Lainey et al. (2009), from the astrometric orbital acceleration: (9.33 +/- 1.83)e13 W.
_IO_HEATING = 9.33e13
_IO_HEATING_UNCERTAINTY = 1.83e13


def _io_with_tides():
    body = _BODIES[0]
    world = build_world(body.name)
    world.solve_eos()
    orbital_frequency = math.sqrt(_G * _MASS_JUPITER / body.semi_major_axis ** 3)
    world.set_spin_frequency(orbital_frequency)
    world.calc_tides(orbital_frequency, orbital_frequency, body.eccentricity, 0.0,
                     body.semi_major_axis, _MASS_JUPITER)
    return world, body


def test_io_total_heating_matches_the_astrometric_value():
    """The asthenosphere viscosity was fitted to this; it is the file's one observational claim."""
    world, _ = _io_with_tides()
    assert world.get_tidal_heating() == pytest.approx(_IO_HEATING, rel=1e-3)


def test_io_surface_heat_flux_is_about_two_watts_per_square_metre():
    world, body = _io_with_tides()
    flux = world.get_tidal_heating() / (4.0 * math.pi * body.radius ** 2)
    assert flux == pytest.approx(2.24, abs=0.05)


def test_io_heating_is_within_the_observational_uncertainty():
    world, _ = _io_with_tides()
    assert abs(world.get_tidal_heating() - _IO_HEATING) < _IO_HEATING_UNCERTAINTY


def test_io_per_layer_heating_sums_to_the_total():
    """With two dissipating layers the scales must partition the heat; equal scales would double it."""
    world, _ = _io_with_tides()
    total = world.get_tidal_heating()
    per_layer = sum(world.get_layer_tidal_heating(index) for index in range(world.num_layers))
    assert per_layer == pytest.approx(total, rel=1e-9)


def test_io_asthenosphere_dominates_the_dissipation():
    """The point of the asthenosphere end-member: the thin soft layer, not the mantle, does the work."""
    world, body = _io_with_tides()
    names = [layer.name for layer in world]
    heating = {name: world.get_layer_tidal_heating(index) for index, name in enumerate(names)}
    # The solid iron core dissipates, as every solid layer does, but its share is a few parts per million.
    assert 0.0 < heating["core"] < 1.0e-5 * world.get_tidal_heating()
    assert heating["asthenosphere"] > 20.0 * heating["mantle"]
    assert heating["asthenosphere"] / world.get_tidal_heating() > 0.9


# =====================================================================================================================
# Mercury's fluid outer core
# =====================================================================================================================
def test_mercury_fluid_core_is_worth_a_factor_of_five_in_its_love_number():
    """Mercury's large k2 is the evidence for its fluid outer core: with that layer solid, k2 drops fivefold."""
    body = next(body for body in _BODIES if body.name == "mercury")
    assert _SOLID_CORE_LOVE_K["mercury"] * 4.0 < body.love_k


# =====================================================================================================================
# Luna's dissipation at two periods
# =====================================================================================================================
# Williams and Boggs (2015): Q = 38 +/- 4 at the monthly period and 41 +/- 9 at the yearly one, with a relaxation
# peak near 100 to 120 days.
_LUNA_Q_MONTH, _LUNA_Q_MONTH_UNCERTAINTY = 38.0, 4.0
_LUNA_Q_YEAR, _LUNA_Q_YEAR_UNCERTAINTY = 41.0, 9.0


def _luna_quality(world, period_days):
    world.solve_love_numbers(2.0 * math.pi / (period_days * 86400.0))
    assert world.love_success, world.love_message
    return world.love_number_k.real / (-world.love_number_k.imag)


def test_luna_quality_factor_at_the_month_and_the_year():
    """The zone viscosity and its top are fitted to these two; the file's comment quotes the result."""
    world = build_world("luna")
    world.solve_eos()
    q_month = _luna_quality(world, 27.3217)
    q_year = _luna_quality(world, 365.25)
    assert q_month == pytest.approx(_LUNA_Q_MONTH, abs=0.5)
    assert q_year == pytest.approx(_LUNA_Q_YEAR, abs=0.5)
    assert abs(q_month - _LUNA_Q_MONTH) < _LUNA_Q_MONTH_UNCERTAINTY
    assert abs(q_year - _LUNA_Q_YEAR) < _LUNA_Q_YEAR_UNCERTAINTY


def test_luna_relaxation_peak_sits_between_the_month_and_the_year():
    """The whole-body Q reaches its minimum near 100 days, where Williams and Boggs (2015) place the peak."""
    world = build_world("luna")
    world.solve_eos()
    q_by_period = {period: _luna_quality(world, period) for period in (10.0, 27.3217, 96.0, 110.0, 365.25, 1000.0)}
    assert q_by_period[96.0] < q_by_period[27.3217] and q_by_period[96.0] < q_by_period[365.25]
    assert q_by_period[110.0] == pytest.approx(q_by_period[96.0], rel=0.05)
    assert q_by_period[10.0] > q_by_period[27.3217] and q_by_period[1000.0] > q_by_period[365.25]


def test_luna_deep_zone_does_most_of_the_dissipating():
    """The 3D split the tidal_scale values record: about 85 percent in the zone, 15 in the mantle above it."""
    world = build_world("luna")
    scales = {layer.name: layer.tidal_scale for layer in world}
    assert scales["lower_mantle"] == pytest.approx(0.852, abs=0.005)
    assert scales["mantle"] == pytest.approx(0.148, abs=0.005)
    assert scales["crust"] < 1.0e-3 and scales["inner_core"] < 1.0e-6
