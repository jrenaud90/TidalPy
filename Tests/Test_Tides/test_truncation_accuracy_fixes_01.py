"""Truncation accuracy limits by degree: a solve through degree L takes the tightest limit of degrees 2 to L."""
import numpy as np
import pytest

from TidalPy.Structures import build_world
from TidalPy.Structures.system import System
from TidalPy.Tides.classes.collapse import collapse_global_tides
from TidalPy.Tides.eccentricity import (
    ECCENTRICITY_TRUNCATIONS, eccentricity_accuracy_limit, recommend_eccentricity_truncation)
from TidalPy.Tides.obliquity import obliquity_accuracy_limit, recommend_obliquity_truncation

DEGREES = (2, 3, 4, 5, 6, 7, 8, 9, 10)
TOLERANCES = (1.0e-8, 1.0e-6, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1)
OBLIQUITY_LEVELS = (2, 4)
# 1.5 sets several of the constant-time-lag limits between the earlier spin rates.
SPIN_RATIOS = (0.5, 1.0, 1.5, 2.3)

_BODY = dict(planet_radius=1.8215e6, orbital_frequency=4.11e-5, semi_major_axis=4.217e8, host_mass=1.898e27,
             G_to_use=6.674e-11)
_TIDES = (("cpl", {"fixed_k": [0.3] * 9, "fixed_q": [100.0] * 9}),
          ("ctl", {"fixed_k": [0.3] * 9, "fixed_dt_s": [100.0] * 9}))

# Demo 17's homogeneous Maxwell Io on Io's orbit about Jupiter. Its rigidity stretches its relaxation time to about 50
# orbital periods, so -Im k falls as 1 / frequency across the tidal modes: the case that sets the obliquity limits at
# synchronous rotation and the high eccentricity levels' limits for a fast rotator.
_IO_RADIUS, _IO_MASS = 1.8216e6, 8.9319e22
_IO_SHEAR, _IO_VISCOSITY = 6.0e10, 1.0e16
_IO_SEMI_MAJOR_AXIS = 4.217e8
_IO_LOW_ECCENTRICITY = 1.0e-3
_IO_REFERENCE_TOLERANCE = 1.0e-12


def _single_degree_heating(
        degree_l,
        spin_ratio,
        tide_model,
        tide_config,
        eccentricity,
        obliquity,
        eccentricity_truncation,
        obliquity_truncation):
    """Heating of one degree alone."""
    return collapse_global_tides(
        **_BODY,
        spin_frequency=spin_ratio * _BODY["orbital_frequency"],
        eccentricity=eccentricity,
        obliquity=obliquity,
        tide_model=tide_model,
        tide_config=tide_config,
        min_degree_l=degree_l,
        max_degree_l=degree_l,
        eccentricity_truncation=eccentricity_truncation,
        eccentricity_exact_tolerance=1.0e-10,
        obliquity_truncation=obliquity_truncation)["tidal_heating"]


def _worst_eccentricity_error(level, degree_l, eccentricity):
    """Largest relative heating error of one degree at a level against the exact functions."""
    worst = 0.0
    for tide_model, tide_config in _TIDES:
        for spin_ratio in SPIN_RATIOS:
            args = (degree_l, spin_ratio, tide_model, tide_config, eccentricity, 0.0)
            truncated = _single_degree_heating(*args, level, 0)
            exact = _single_degree_heating(*args, "exact", 0)
            worst = max(worst, abs(truncated / exact - 1.0))
    return worst


def _worst_obliquity_error(level, degree_l, obliquity):
    """Largest relative heating error of one degree at a level against the general functions, at e = 0."""
    worst = 0.0
    for tide_model, tide_config in _TIDES:
        for spin_ratio in SPIN_RATIOS:
            args = (degree_l, spin_ratio, tide_model, tide_config, 0.0, obliquity, 2)
            worst = max(worst, abs(_single_degree_heating(*args, level) / _single_degree_heating(*args, "gen") - 1.0))
    return worst


@pytest.mark.parametrize("level, degree_two, degree_three", [
    (2, 0.075, 0.06),
    (10, 0.395, 0.335),
    (20, 0.5, 0.435),
    (50, 0.57, 0.525),
])
def test_eccentricity_limits_through_degree_three_are_the_measured_values(level, degree_two, degree_three):
    """The 10% limits of a solve through degree 2 or 3 are the measured degree-2 and degree-3 values."""
    assert eccentricity_accuracy_limit(level, 1.0e-1, 2) == degree_two
    assert eccentricity_accuracy_limit(level, 1.0e-1, 3) == degree_three


@pytest.mark.parametrize("level, degree_two, degree_three", [(2, 0.35, 0.24), (4, 0.695, 0.48)])
def test_obliquity_limits_through_degree_three_are_the_measured_values(level, degree_two, degree_three):
    """The 10% limits of a solve through degree 2 or 3 are the measured degree-2 and degree-3 values."""
    assert obliquity_accuracy_limit(level, 1.0e-1, 2) == degree_two
    assert obliquity_accuracy_limit(level, 1.0e-1, 3) == degree_three


@pytest.mark.parametrize("tolerance", TOLERANCES)
@pytest.mark.parametrize("level", ECCENTRICITY_TRUNCATIONS)
def test_eccentricity_limits_never_grow_with_degree(level, tolerance):
    """Including a higher degree never raises the limit; a degree past the table uses degree 10's."""
    limits = [eccentricity_accuracy_limit(level, tolerance, degree_l) for degree_l in DEGREES]
    assert all(higher <= lower for lower, higher in zip(limits, limits[1:]))
    assert eccentricity_accuracy_limit(level, tolerance, 11) == limits[-1]


@pytest.mark.parametrize("tolerance", TOLERANCES)
@pytest.mark.parametrize("level", OBLIQUITY_LEVELS)
def test_obliquity_limits_never_grow_with_degree(level, tolerance):
    """Including a higher degree never raises the limit; a degree past the table uses degree 10's."""
    limits = [obliquity_accuracy_limit(level, tolerance, degree_l) for degree_l in DEGREES]
    assert all(higher <= lower for lower, higher in zip(limits, limits[1:]))
    assert obliquity_accuracy_limit(level, tolerance, 11) == limits[-1]


def test_high_degrees_tighten_the_limits():
    """Degrees 4 to 10 lose accuracy sooner than degree 3, so a solve through degree 10 has tighter limits."""
    for level in ECCENTRICITY_TRUNCATIONS:
        assert eccentricity_accuracy_limit(level, 1.0e-1, 10) < eccentricity_accuracy_limit(level, 1.0e-1, 3)
    for level in OBLIQUITY_LEVELS:
        assert obliquity_accuracy_limit(level, 1.0e-1, 10) < obliquity_accuracy_limit(level, 1.0e-1, 3)


@pytest.mark.parametrize("tolerance", (1.0e-2, 1.0e-1))
@pytest.mark.parametrize("level", ECCENTRICITY_TRUNCATIONS)
def test_eccentricity_limit_through_degree_ten_holds_at_every_degree(level, tolerance):
    """At a level's limit for a solve through degree 10, every degree's heating is within the tolerance."""
    eccentricity = eccentricity_accuracy_limit(level, tolerance, 10)
    assert eccentricity > 0.0
    for degree_l in DEGREES:
        assert _worst_eccentricity_error(level, degree_l, eccentricity) <= tolerance, degree_l


@pytest.mark.parametrize("tolerance", (1.0e-2, 1.0e-1))
@pytest.mark.parametrize("level", OBLIQUITY_LEVELS)
def test_obliquity_limit_through_degree_ten_holds_at_every_degree(level, tolerance):
    """At a level's limit for a solve through degree 10, every degree's heating is within the tolerance."""
    obliquity = obliquity_accuracy_limit(level, tolerance, 10)
    assert obliquity > 0.0
    for degree_l in DEGREES:
        assert _worst_obliquity_error(level, degree_l, obliquity) <= tolerance, degree_l


@pytest.mark.parametrize("eccentricity, expected", [(0.05, 4), (0.2, 20), (0.3, 50)])
def test_recommended_eccentricity_level_holds_at_every_degree(eccentricity, expected):
    """The level recommended for a solve through degree 10 holds 10% at every degree."""
    level = recommend_eccentricity_truncation(eccentricity, tolerance=1.0e-1, max_degree_l=10)
    assert level == expected
    for degree_l in DEGREES:
        assert _worst_eccentricity_error(level, degree_l, eccentricity) <= 1.0e-1, degree_l


@pytest.mark.parametrize("obliquity", (0.05, 0.15))
def test_recommended_obliquity_level_holds_at_every_degree(obliquity):
    """The level recommended for a solve through degree 10 holds 10% at every degree."""
    level = recommend_obliquity_truncation(obliquity, tolerance=1.0e-1, max_degree_l=10)
    assert level in OBLIQUITY_LEVELS
    for degree_l in DEGREES:
        assert _worst_obliquity_error(level, degree_l, obliquity) <= 1.0e-1, degree_l


def test_level_two_is_recommended_for_tight_tolerances_at_small_eccentricity():
    """Measured down to e = 1e-5, level 2's 1e-4 limit is no longer 0.0, so it is recommended at small e."""
    assert eccentricity_accuracy_limit(2, 1.0e-4) > 0.0
    assert recommend_eccentricity_truncation(5.0e-4, tolerance=1.0e-4) == 2


@pytest.fixture(scope="module")
def maxwell_io():
    """Demo 17's Maxwell Io and its Jupiter system; each test sets the spin rate it needs."""
    density = _IO_MASS / (4.0 / 3.0 * np.pi * _IO_RADIUS**3)
    io = build_world({
        "schema_version": "0.2.0", "name": "Io", "type": "terrestrial", "radius_m": _IO_RADIUS, "mass_kg": _IO_MASS,
        "tides": {"global_tidal_model": "rheology", "love_method": "homogeneous"},
        "layers": {"body": {
            "layer_index": 0, "radius_fraction": 1.0, "use_tides": True, "is_static": False,
            "material": {"solid": {
                "eos": {"model": "constant", "reference_density_kg_m3": density, "bulk_modulus_pa": 2.0 * _IO_SHEAR},
                "shear_modulus": {"model": "constant", "shear_modulus_pa": _IO_SHEAR},
                "shear_viscosity": {"model": "constant", "reference_viscosity_pas": _IO_VISCOSITY}}},
            "shear_rheology": {"model": "maxwell"}, "bulk_rheology": {"model": "elastic"}}}})
    io.solve_eos()
    jupiter = build_world({"schema_version": "0.2.0", "name": "Jupiter", "type": "gasgiant", "radius_m": 6.9911e7,
                           "mass_kg": 1.898e27,
                           "layers": {"envelope": {"material": "simple_gas", "radius_fraction": 1.0}}})
    system = System("Jupiter-Io")
    system.add_world(jupiter)
    system.add_world(io, tidal_host=jupiter, semi_major_axis=_IO_SEMI_MAJOR_AXIS, eccentricity=_IO_LOW_ECCENTRICITY)
    return io, system


def _io_heating(maxwell_io, spin_ratio, eccentricity, eccentricity_truncation, obliquity=0.0,
                obliquity_truncation=0):
    """Degree-2 heating of the Maxwell Io spinning at spin_ratio times its mean motion."""
    io, system = maxwell_io
    system.set_eccentricity(io, eccentricity)
    io.set_spin_frequency(spin_ratio * system.calc_orbital_frequency(io))
    io.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=eccentricity_truncation,
                       eccentricity_exact_tolerance=_IO_REFERENCE_TOLERANCE, obliquity_truncation=obliquity_truncation)
    io.set_obliquity(obliquity)
    return system.calc_world_evolution(io)["tidal_heating"]


@pytest.mark.parametrize("degrees, expected, lower", [(8.0, 4, 2), (23.4, "gen", 4)])
def test_recommended_obliquity_level_holds_for_the_synchronous_maxwell_io(maxwell_io, degrees, expected, lower):
    """Demo 17's synchronous Io at e = 0.001 needs the recommended level for 1%; the level below it misses.

    The earlier table recommended level 2 at 8 degrees (1.5% off here) and level 4 at 23.4 degrees (1.02% off).
    """
    obliquity = np.radians(degrees)
    level = recommend_obliquity_truncation(obliquity, tolerance=1.0e-2)
    assert level == expected

    general = _io_heating(maxwell_io, 1.0, _IO_LOW_ECCENTRICITY, 10, obliquity, "gen")

    def error(truncation):
        return abs(_io_heating(maxwell_io, 1.0, _IO_LOW_ECCENTRICITY, 10, obliquity, truncation) / general - 1.0)

    assert error(level) < 1.0e-2
    assert error(lower) > 1.0e-2


def test_level_fifty_limit_holds_for_a_fast_rotating_maxwell_io(maxwell_io):
    """Spinning at 11.5 times its mean motion, the Maxwell Io holds level 50 within 10% at the tabulated limit.

    Past it the error grows fast: about 70% at e = 0.6, well inside the 0.78 measured earlier at slower spins.
    """
    fast_spin_ratio = 11.5

    def error(eccentricity):
        exact = _io_heating(maxwell_io, fast_spin_ratio, eccentricity, "exact")
        return abs(_io_heating(maxwell_io, fast_spin_ratio, eccentricity, 50) / exact - 1.0)

    assert error(eccentricity_accuracy_limit(50, 1.0e-1)) < 1.0e-1
    assert error(0.6) > 1.0e-1
