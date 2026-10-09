"""A Love solve whose dynamic liquid layer loses accuracy at the forcing period logs one warning per world.

A dynamic liquid whose density does not follow its bulk modulus (the bundled worlds' constant-density liquids) is
unstably stratified, and its solutions grow as exp(E) with E proportional to the forcing period. The same layer
solved as static, or as dynamic but incompressible, is neutral and is not flagged.
"""
import math

import numpy as np
import pytest

from TidalPy.RadialSolver.solver import radial_solver
from TidalPy.Rheology import Maxwell
from TidalPy.Structures import build_world


_WARNING_TEXT = "grows unstable at long periods"
_DAY = 86400.0


def _world(name, is_static=True, is_incompressible=False):
    world = build_world(name)
    world.solve_eos()
    for layer in world.layers:
        if layer.is_liquid:
            layer.is_static = is_static
            layer.is_incompressible = is_incompressible
    return world


def _solve(world, period_days, **kwargs):
    """Solve; the warning is logged before the solve runs, and an unstable solve reports its failure rather than
    raising."""
    world.solve_love_numbers(frequency=2.0 * math.pi / (period_days * _DAY), degree_l=2, **kwargs)


def test_dynamic_constant_density_liquid_warns_once(spdlog_text):
    world = _world("earth_simple", is_static=False)
    _solve(world, 10.0)
    _solve(world, 10.0)
    text = spdlog_text()
    assert text.count(_WARNING_TEXT) == 1
    assert f"'{world.name}'" in text and "'outer_core'" in text


@pytest.mark.parametrize("is_static, is_incompressible", [(True, False), (False, True)])
def test_neutral_liquid_is_not_flagged(spdlog_text, is_static, is_incompressible):
    """A static liquid, or a dynamic incompressible one of constant density, has no unstable stratification."""
    _solve(_world("earth_simple", is_static=is_static, is_incompressible=is_incompressible), 100.0)
    assert _WARNING_TEXT not in spdlog_text()


def test_flag_follows_the_forcing_period(spdlog_text):
    """Luna's small outer core stays accurate at 3.55 days and is flagged at 100 days."""
    world = _world("luna", is_static=False)
    _solve(world, 3.55)
    assert _WARNING_TEXT not in spdlog_text()
    _solve(world, 100.0)
    assert spdlog_text().count(_WARNING_TEXT) == 1


def test_warnings_false_silences_it(spdlog_text):
    _solve(_world("earth_simple", is_static=False), 10.0, warnings=False)
    assert _WARNING_TEXT not in spdlog_text()


def test_calc_tides_checks_before_its_solves(spdlog_text):
    world = _world("earth_simple", is_static=False)
    orbital_frequency = 2.0 * math.pi / (10.0 * _DAY)
    world.calc_tides(orbital_frequency, orbital_frequency, 0.01, 0.0, 1.0e9, 7.3e22)
    assert spdlog_text().count(_WARNING_TEXT) == 1


def _three_layer(frequency, liquid_is_static, warnings):
    num_per_layer = 30
    radius = 6.0e6
    icb, cmb = radius / 3.0, 2.0 * radius / 3.0
    radii = np.concatenate((np.linspace(0.0, icb, num_per_layer), np.linspace(icb, cmb, num_per_layer),
                            np.linspace(cmb, radius, num_per_layer)))
    density = np.repeat((8500.0, 7000.0, 3500.0), num_per_layer)
    shear = Maxwell().calc_complex_modulus_vectorize_modulus(
        np.repeat((1.0e11, 0.0, 5.0e10), num_per_layer), np.repeat((1.0e26, 1.0e6, 1.0e20), num_per_layer), frequency)
    radial_solver(radii, density, np.full(radii.size, 1.0e11 + 0.0j), shear, frequency, float(np.mean(density)),
                  ("solid", "liquid", "solid"), (False, liquid_is_static, False), (False, False, False),
                  np.asarray((icb, cmb, radius)), degree_l=2, raise_on_fail=False, warnings=warnings)


def test_standalone_radial_solver_uses_the_same_check(spdlog_text):
    frequency = 2.0 * math.pi / (10.0 * _DAY)
    _three_layer(frequency, liquid_is_static=True, warnings=True)
    _three_layer(frequency, liquid_is_static=False, warnings=False)
    assert _WARNING_TEXT not in spdlog_text()
    _three_layer(frequency, liquid_is_static=False, warnings=True)
    assert spdlog_text().count(_WARNING_TEXT) == 1
