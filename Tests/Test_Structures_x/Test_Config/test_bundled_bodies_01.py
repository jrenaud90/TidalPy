"""Bundled worlds rebuilt from the classic BurnMan WorldPack: the solar-system bodies and the TRAPPIST-1 system.

Each of these worlds states a published mass and radius and then describes an interior that has to
reproduce that mass when the equation of state is solved. A layer density copied from the wrong source, a
mistyped boundary radius, or an EOS model that compresses differently than intended all show up as a mass
mismatch, so the mass check is the substance of these tests rather than a formality. Surface gravity is
checked against its published value as an independent handle: it follows from mass and radius, so it
catches a radius that disagrees with the one the interior was fitted to.

The moment-of-inertia factor is the stronger check, because only Luna's crust thickness and Earth-Simple's
densities are fitted to it. C/MR2 measures how centrally concentrated the mass is, so it constrains the core
size rather than the total mass, and each body carries its own tolerance recording how well its model does
against the measured value. Pluto, Charon, Triton, and the seven TRAPPIST-1 planets have no measured moment of
inertia, so for them the expected C/MR2 is the value the file's comment records: that pins the model and keeps
the comment from rotting, but it is not a check on the interior. ``Body.moi_measured`` says which is which.

Every layer states its temperature and every solid layer carries a reference-law viscosity at it, so every
solid layer dissipates, and each world's per-layer tidal_scale values come from integrating the depth-resolved
3D heating (author decisions 2026-09-19 and 2026-09-21). Four bodies pin a tidal observable as well. Io has
its asthenosphere viscosity fitted so that the total dissipation matches the astrometric value of Lainey et
al. (2009). Mercury, Luna, and Earth-Simple declare their fluid outer cores as static liquids and have their
mantle rigidity fitted to the measured k2; for each the fluid core is worth a large factor in k2, so the
solid-core counterfactual is pinned too. Luna's low-viscosity zone at the base of the mantle is fitted to
the lunar Q at the month and at the year (Williams and Boggs 2015).

Pluto carries a subsurface ocean whose thickness follows from the observed ice-shell thickness rather than
from a fit, and the ocean is worth a factor of 31 in its degree-2 Love number, which is pinned below.

The gas giants and the star TRAPPIST-1 are checked at the end of the file rather than through the table,
because neither fits the table's assumptions. Jupiter and Neptune have their two outer densities fitted to
the mass and the published moment-of-inertia factor at once, and their tidal response is analytic, so it is
checked against its closed form. TRAPPIST-1 is checked against the orbital periods of its seven planets.

Expected values come from the TOML comments, which record where each number came from.
"""
import math

import pytest

from TidalPy.structures_x.configs import build_world
from TidalPy.structures_x.worlds.gasgiant import GasGiantWorld
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
                 moi_measured=True,
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
        self.moi_measured = moi_measured
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
    # Pluto and Charon are mutually synchronous on the 6.3872304 day mutual orbit; Triton is synchronous and
    # retrograde about Neptune. None of the three has a measured moment of inertia. Masses come from the GM
    # values of Brozovic and Jacobson (2024) for the Pluto system and JPL NEP097 for Triton.
    Body("pluto", 1188300.0, 1.3024587e22, 0.615627,
         ["core", "ocean", "ice_shell"], ["core", "ice_shell"],
         6.3872304, moi_factor=0.318550, moi_tolerance=1.0e-3, moi_measured=False,
         liquid_layers=["ocean"]),
    Body("charon", 606000.0, 1.5896798e21, 0.288915,
         ["core", "ice_shell"], ["core", "ice_shell"],
         6.3872304, moi_factor=0.311564, moi_tolerance=1.0e-3, moi_measured=False),
    Body("triton", 1352600.0, 2.1402926e22, 0.780801,
         ["core", "ice_shell"], ["core", "ice_shell"],
         5.876854, moi_factor=0.315542, moi_tolerance=1.0e-3, moi_measured=False),
    # The seven TRAPPIST-1 planets, all assumed tidally locked, so the spin period is the orbital period.
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
def test_bundled_body_moment_of_inertia_factor_matches_its_expected_value(body):
    """Nothing is fitted to C/MR2, so where one is measured it is an independent test of where the mass sits.

    Where none is measured (``moi_measured=False``) the expected value is the one the file's comment quotes,
    which pins the model rather than testing it.
    """
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
    """Every liquid layer is a static liquid: the fluid outer cores, and Pluto's subsurface ocean."""
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
def test_bundled_body_spins_at_its_stated_rotation_period(body):
    """Spin is the synchronous rate of each body's rotation period, except for Mercury and Earth.

    The solar-system bodies take that period from the classic file's spin_period; the TRAPPIST-1 planets are
    assumed tidally locked, so theirs is the measured orbital period.
    """
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


# =====================================================================================================================
# The TRAPPIST-1 system
# =====================================================================================================================
# Agol et al. (2021): the stellar mass and luminosity, and each planet's semi-major axis and orbital period. The
# stellar radius is the Delrez et al. (2018) value carried in TidalPy's constants_.hpp.
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

# Earth's core mass fraction through the identical two-layer recipe the seven planets use. It is larger than
# Earth's accepted 0.325 because a uniform incompressible core is a crude stand-in for a real one, which is
# exactly why the comparison below is made against this number and not against the accepted one.
_EARTH_RECIPE_CORE_MASS_FRACTION = 0.354


def test_trappist1_star_states_its_published_properties():
    star = build_world("trappist1")
    assert star.radius == pytest.approx(_TRAPPIST1_RADIUS, rel=1e-12)
    assert star.mass == pytest.approx(_TRAPPIST1_MASS, rel=1e-12)
    assert star.effective_temperature == pytest.approx(_TRAPPIST1_EFFECTIVE_TEMPERATURE, rel=1e-12)


def test_trappist1_luminosity_is_the_measured_one_not_the_derived_one():
    """The file states the luminosity, because Stefan-Boltzmann on the published radius and temperature misses
    the published luminosity by 1.4 percent."""
    star = build_world("trappist1")
    assert star.luminosity == pytest.approx(_TRAPPIST1_LUMINOSITY, rel=1e-9)

    stefan_boltzmann = 5.670374419e-8
    derived = (4.0 * math.pi * _TRAPPIST1_RADIUS ** 2 * stefan_boltzmann
               * _TRAPPIST1_EFFECTIVE_TEMPERATURE ** 4)
    assert derived == pytest.approx(_TRAPPIST1_LUMINOSITY, rel=2.0e-2)
    assert derived != pytest.approx(_TRAPPIST1_LUMINOSITY, rel=1.0e-2)


@pytest.mark.parametrize("name", sorted(_TRAPPIST1_PLANETS))
def test_trappist1_planet_period_follows_from_the_stellar_mass_and_its_semi_major_axis(name):
    """Kepler's third law ties the star file to the seven planet files: the mean motion from the stellar mass
    and the published semi-major axis must reproduce the measured orbital period, which is also each planet's
    spin. The published semi-major axes carry four significant digits, which sets the tolerance."""
    semi_major_axis_au, period_days, _ = _TRAPPIST1_PLANETS[name]
    mean_motion = math.sqrt(_G * _TRAPPIST1_MASS / (semi_major_axis_au * _AU) ** 3)
    assert mean_motion == pytest.approx(2.0 * math.pi / (period_days * 86400.0), rel=1.0e-3)

    planet = build_world(name)
    assert planet.spin_frequency == pytest.approx(2.0 * math.pi / (period_days * 86400.0), rel=1e-6)


@pytest.mark.parametrize("name", sorted(_TRAPPIST1_PLANETS))
def test_trappist1_planet_is_iron_depleted_relative_to_earth(name):
    """Nothing is fitted to the core mass fraction: the core radius is fitted to the mass and the fraction
    follows. Every one of the seven lands below the Earth the same recipe produces, which is how these files
    reproduce the iron depletion Agol et al. (2021) infer from the densities."""
    _, _, expected_core_mass_fraction = _TRAPPIST1_PLANETS[name]
    world = build_world(name)
    world.solve_eos()
    core_mass_fraction = world.core.mass / world.planet_mass_eos
    assert core_mass_fraction == pytest.approx(expected_core_mass_fraction, rel=1e-3)
    assert core_mass_fraction < _EARTH_RECIPE_CORE_MASS_FRACTION


# =====================================================================================================================
# The gas giants
# =====================================================================================================================
# Jupiter and Neptune each carry a three-layer interior whose two outer densities are fitted to the mass and to the
# published moment-of-inertia factor at once, so C/MR2 is reproduced by construction rather than tested. What the
# tests check is that the fit still holds, that the layer masses it implies are the ones the files quote, and that
# the analytic tidal response matches its closed form. jupiter_simple keeps one uniform layer, so it reproduces the
# mass but not C/MR2, which is the difference between the two files.
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
    """Two fitted densities against two observables, so both must come back."""
    radius, mass, gravity, moi_factor, _, _, _ = _GAS_GIANTS[name]
    world = build_world(name)
    result = world.solve_eos()
    assert result["success"], result["message"]
    assert world.planet_mass_eos == pytest.approx(mass, rel=1e-6)
    assert world.surface_gravity_eos == pytest.approx(gravity, rel=1e-3)
    assert world.planet_moi_eos / (world.planet_mass_eos * radius ** 2) == pytest.approx(
        moi_factor, rel=1e-4)


# Love number from solving the fitted interior, against the k2 each file states from the literature. Nothing is
# fitted to this: the densities come from the mass and C/MR2, and the Love number falls out of them.
_GAS_GIANT_SOLVED_LOVE_K = {"jupiter": 0.53438, "neptune": 0.42667}


@pytest.mark.parametrize("name", sorted(_GAS_GIANTS))
def test_gas_giant_layers_are_all_fluid(name):
    """A gas giant is a fluid body. A rigid rock core is wrong physically at these temperatures and fails
    numerically as well: a rock shear modulus is negligible against rho g R, so the core integrates as a
    near-fluid solid and the shooting solve dies on its step size."""
    world = build_world(name)
    assert [layer.name for layer in world if layer.is_solid] == []
    assert all(layer.is_static for layer in world)


@pytest.mark.parametrize("name", sorted(_GAS_GIANTS))
def test_gas_giant_solved_love_number_is_near_its_published_value(name):
    """The interior is fitted to the mass and C/MR2 only, so its Love number is an independent result. Both
    land near the published value, where a single uniform layer gives the fluid-sphere 1.5."""
    _, _, _, _, _, published_k, _ = _GAS_GIANTS[name]
    world = build_world(name)
    world.solve_eos()
    result = world.solve_love_numbers()
    assert result["success"], result["message"]
    assert world.love_number_k.real == pytest.approx(_GAS_GIANT_SOLVED_LOVE_K[name], rel=1e-3)
    # Within 15 percent of the published Love number, with nothing tuned to it.
    assert abs(world.love_number_k.real - published_k) / published_k < 0.15


@pytest.mark.parametrize("name", sorted(_GAS_GIANTS))
def test_gas_giant_layer_masses_match_its_file(name):
    """The layer masses are a result of the fit, not an input, and they are what ties it to published
    interior models: Neptune's 14.6 Earth masses of heavy elements, Jupiter's 16 Earth mass core."""
    world = build_world(name)
    world.solve_eos()
    for layer in world:
        expected = _GAS_GIANT_LAYER_MASSES[name][layer.name]
        assert layer.mass / _MASS_EARTH == pytest.approx(expected, rel=1e-3), layer.name
    # The envelope is the only tidal layer, and it carries the whole scale.
    assert [layer.name for layer in world if layer.is_tidal] == ["envelope"]
    assert sum(layer.tidal_scale for layer in world if layer.is_tidal) == pytest.approx(1.0, abs=1e-4)


@pytest.mark.parametrize(
    "name, semi_major_axis, satellite_mass, spin_multiple",
    [("neptune", _TRITON_SEMI_MAJOR_AXIS, _TRITON_MASS, None),
     ("neptune", _TRITON_SEMI_MAJOR_AXIS, _TRITON_MASS, 3.0),
     ("neptune", 6.0e8, _TRITON_MASS, None),
     ("neptune", _TRITON_SEMI_MAJOR_AXIS, 1.0e23, None),
     ("jupiter", _IO_SEMI_MAJOR_AXIS, _IO_MASS, None),
     ("jupiter", _IO_SEMI_MAJOR_AXIS, _IO_MASS, 2.0)])
def test_gas_giant_dissipation_matches_the_constant_phase_lag_closed_form(
        name, semi_major_axis, satellite_mass, spin_multiple):
    """With one active mode (circular orbit, no obliquity, degree 2) the constant-phase-lag heating is
    (3/4)(k2/Q)(G M_s^2 R^5 / a^6)|2(spin - n)|. Neither gas giant carries a viscoelastic interior, so this
    is the whole of its tidal response and the check is exact rather than a model comparison."""
    radius, mass, _, _, _, love_k, quality = _GAS_GIANTS[name]
    world = build_world(name)
    mean_motion = math.sqrt(_G * mass / semi_major_axis ** 3)
    spin = world.spin_frequency if spin_multiple is None else spin_multiple * mean_motion

    world.calc_tides(mean_motion, spin, 0.0, 0.0, semi_major_axis, satellite_mass)
    closed_form = (0.75 * (love_k / quality)
                   * (_G * satellite_mass ** 2 * radius ** 5 / semi_major_axis ** 6)
                   * abs(2.0 * (spin - mean_motion)))
    assert world.get_tidal_heating() == pytest.approx(closed_form, rel=1e-3)


def test_jupiter_simple_reproduces_its_mass_but_not_its_moment_of_inertia_or_love_number():
    """One uniform layer can match the mass and nothing else: C/MR2 is the uniform-sphere 0.4 against
    Jupiter's 0.264 and k2 is the fluid-sphere 1.5 against 0.590, which is what jupiter.toml exists to fix."""
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
# The ice shell is 165 km thick, the mean global depth to the top of the ocean that Nimmo et al. (2016) derive from
# the reorientation of Sputnik Planitia. The ocean is what is left above the core, and its thickness is therefore a
# consequence of that shell and the classic core radius rather than a fitted number.
_PLUTO_SHELL_THICKNESS = 165.0e3
_PLUTO_OCEAN_THICKNESS = 112.4e3
_PLUTO_OCEAN_K = 0.145781
_PLUTO_FROZEN_K = 0.004682


def test_pluto_ocean_thickness_follows_from_the_observed_shell():
    """The shell thickness is the observational input; the ocean is the rest of the water layer, and it lands
    inside the 50 to 150 km that thermal evolution models give without being fitted to."""
    world = build_world("pluto")
    shell = world.radius - world.ocean.radius_outer
    ocean = world.ocean.radius_outer - world.core.radius_outer
    assert shell == pytest.approx(_PLUTO_SHELL_THICKNESS, rel=1e-6)
    assert ocean == pytest.approx(_PLUTO_OCEAN_THICKNESS, rel=1e-3)
    assert 50.0e3 < ocean < 150.0e3


def test_pluto_ocean_is_worth_a_factor_of_thirty_in_its_love_number():
    """A static liquid layer decouples the shell from the core. Nothing is fitted to this, since Pluto has no
    measured k2, but it is the largest single consequence of the ocean."""
    world = build_world("pluto")
    frequency = world.spin_frequency
    world.solve_eos()
    world.solve_love_numbers(frequency)
    assert world.love_success, world.love_message
    ocean_k = world.love_number_k.real
    assert ocean_k == pytest.approx(_PLUTO_OCEAN_K, rel=1e-3)

    world.ocean.is_solid = True
    world.solve_eos()
    world.solve_love_numbers(frequency)
    frozen_k = world.love_number_k.real
    assert frozen_k == pytest.approx(_PLUTO_FROZEN_K, rel=1e-3)
    assert ocean_k > 25.0 * frozen_k


def test_pluto_ocean_takes_none_of_the_tidal_heating():
    """A liquid layer contributes no shear dissipation, so the ocean carries no tidal_scale and the solid
    layers' scales still partition the whole."""
    world = build_world("pluto")
    assert not world.ocean.is_tidal
    scales = [layer.tidal_scale for layer in world if layer.is_tidal]
    assert sum(scales) == pytest.approx(1.0, abs=1e-4)
    # The decoupled shell takes nearly all of it; without the ocean the core would take 4 percent.
    assert world.ice_shell.tidal_scale > 0.99
