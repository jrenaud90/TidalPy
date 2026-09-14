"""Bundled solar-system bodies rebuilt from the classic BurnMan WorldPack.

Each of these worlds states a published mass and radius and then describes an interior that has to
reproduce that mass when the equation of state is solved. A layer density copied from the wrong source, a
mistyped boundary radius, or an EOS model that compresses differently than intended all show up as a mass
mismatch, so the mass check is the substance of these tests rather than a formality. Surface gravity is
checked against its published value as an independent handle: it follows from mass and radius, so it
catches a radius that disagrees with the one the interior was fitted to.

The moment-of-inertia factor is the stronger check, because no world file is fitted to it. C/MR2 measures
how centrally concentrated the mass is, so it constrains the core size rather than the total mass, and
each body carries its own tolerance recording how well its model does against the measured value.

Three bodies pin a tidal observable as well. Io has its asthenosphere viscosity fitted so that the total
dissipation matches the astrometric value of Lainey et al. (2009), and its per-layer tidal_scale values
come from integrating the depth-resolved 3D heating; both are claims its comments make and either could
rot silently. Luna's k2 is not fitted to anything and comes out within half a percent of the lunar laser
ranging value, which is the check that its density and rigidity profile is right. Mercury's k2 only works
with a fluid outer core, which the schema cannot declare, so the size of that gap is pinned too.

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
         ["core", "mantle", "asthenosphere"], ["mantle", "asthenosphere"],
         1.769, moi_factor=0.37685, moi_tolerance=2.0e-2,
         semi_major_axis=4.217e8, eccentricity=0.0041),
    Body("europa", 1561000.0, 4.7998e22, 1.315,
         ["core", "mantle", "ice_shell"], ["ice_shell"],
         3.551181, moi_factor=0.346, moi_tolerance=1.0e-2),
    Body("luna", 1737400.0, 7.34579e22, 1.6242,
         ["inner_core", "outer_core", "mantle", "crust"], ["mantle"],
         27.322, moi_factor=0.3931, moi_tolerance=1.0e-3,
         love_k=0.02422, love_period_days=27.3217, love_tolerance=1.0e-2),
    Body("mercury", 2440000.0, 3.30103e23, 3.7007,
         ["inner_core", "outer_core", "mantle"], ["mantle"],
         58.6462, moi_factor=0.346, moi_tolerance=1.0e-2),
]


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
def test_bundled_body_tidal_scales_sum_to_one(body):
    """Per-layer scales must partition the world's heating, not repeat it."""
    world = build_world(body.name)
    scales = [layer.tidal_scale for layer in world if layer.is_tidal]
    assert sum(scales) == pytest.approx(1.0, abs=1e-4), scales


@pytest.mark.parametrize("body", _cases())
def test_bundled_body_spins_at_its_classic_period(body):
    """Spin came from the classic file's spin_period, synchronous for every body here except Mercury."""
    world = build_world(body.name)
    assert world.spin_frequency == pytest.approx(
        2.0 * math.pi / (body.spin_period_days * 86400.0), rel=1e-6)


@pytest.mark.parametrize("body", _love_cases())
def test_bundled_body_love_number_matches_its_measured_value(body):
    """An untuned k2 that lands on the measured one is the best evidence the interior is right."""
    world = build_world(body.name)
    world.solve_eos()
    world.solve_love_numbers(2.0 * math.pi / (body.love_period_days * 86400.0))
    assert world.love_success, world.love_message
    assert world.love_number_k.real == pytest.approx(body.love_k, rel=body.love_tolerance)


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
    assert heating["core"] == 0.0
    assert heating["asthenosphere"] > 20.0 * heating["mantle"]
    assert heating["asthenosphere"] / world.get_tidal_heating() > 0.9


# =====================================================================================================================
# Mercury's fluid outer core
# =====================================================================================================================
# Genova et al. (2019), from the MESSENGER gravity field: k2 = 0.569 +/- 0.025.
_MERCURY_LOVE_K = 0.569
_MERCURY_TIDAL_PERIOD_DAYS = 87.9691


def test_mercury_needs_a_fluid_outer_core_for_its_measured_love_number():
    """The schema cannot mark a layer fluid, so the built world understates k2 by a factor of four.

    Setting is_solid on the built layer is the documented way to recover it. This pins both numbers, so
    that the file's comment stays true and so that a future schema flag has a value to reproduce.
    """
    frequency = 2.0 * math.pi / (_MERCURY_TIDAL_PERIOD_DAYS * 86400.0)
    world = build_world("mercury")
    world.solve_eos()
    world.solve_love_numbers(frequency)
    solid_core_k = world.love_number_k.real
    assert solid_core_k == pytest.approx(0.1106, rel=1e-2)
    assert solid_core_k < 0.5 * _MERCURY_LOVE_K

    world.outer_core.is_solid = False
    world.solve_eos()
    world.solve_love_numbers(frequency)
    fluid_core_k = world.love_number_k.real
    assert fluid_core_k == pytest.approx(0.4784, rel=1e-2)
    assert fluid_core_k > 4.0 * solid_core_k
    # What is left is the uncalibrated mantle rigidity, not the core state.
    assert abs(fluid_core_k - _MERCURY_LOVE_K) / _MERCURY_LOVE_K < 0.2
